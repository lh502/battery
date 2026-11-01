#!/usr/bin/env python3
"""Fetch earliest official CALCE workbook members using verified ZIP ranges.

Selection reproduces run_multicell.read_arbin's candidate construction and stops
after the first chronological workbook that takes the count to at least 101.
No model fitting is performed. Whole archive hashes are deliberately omitted.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import hashlib
import io
import json
import re
import struct
import time
import urllib.request
import zlib
from pathlib import Path

import pandas as pd


def date_key(name):
    match = re.search(r'_(\d+)_(\d+)_(\d+)\.xlsx$', name)
    if not match:
        raise ValueError(f'Unexpected CALCE workbook name: {name}')
    return datetime.datetime.strptime('_'.join(match.groups()), '%m_%d_%y')


def sha(data):
    return hashlib.sha256(data).hexdigest()


class RemoteZip:
    def __init__(self, url, workers=8, block_bytes=131072):
        self.url = url
        self.workers = workers
        self.block_bytes = block_bytes
        self.total = None
        self.etag = None
        self.last_modified = None
        tail, start = self._request(None, None, suffix=65536)
        offset = tail.rfind(b'PK\x05\x06')
        if offset < 0:
            raise ValueError(f'{url}: no ZIP end-of-central-directory signature')
        end = struct.unpack_from('<4s4H2LH', tail, offset)
        _, disk, central_disk, disk_count, total_count, central_size, central_start, comment_size = end
        if disk or central_disk or disk_count != total_count:
            raise ValueError('Multi-disk ZIP archives are unsupported')
        if offset + 22 + comment_size != len(tail):
            raise ValueError('ZIP tail length does not match end-of-central-directory record')
        if central_start >= start and central_start + central_size <= start + len(tail):
            central = tail[central_start-start:central_start-start+central_size]
        else:
            central = self.range(central_start, central_start + central_size - 1)
        self.central_sha256 = sha(central)
        self.entries = []
        pos = 0
        while pos < len(central):
            h = struct.unpack_from('<4s6H3L5H2L', central, pos)
            if h[0] != b'PK\x01\x02':
                raise ValueError('Invalid ZIP central-directory entry')
            name_length, extra_length, comment_length = h[10:13]
            name_bytes = central[pos+46:pos+46+name_length]
            name = name_bytes.decode('utf-8' if h[3] & 0x800 else 'cp437')
            self.entries.append({'member': name, 'flags': h[3], 'method': h[4],
                                 'crc32': h[7], 'compressed_bytes': h[8],
                                 'size_bytes': h[9], 'local_offset': h[16]})
            pos += 46 + name_length + extra_length + comment_length
        if len(self.entries) != total_count or pos != len(central):
            raise ValueError('ZIP central-directory entry count/length mismatch')

    def _request(self, start, end, suffix=None):
        range_value = f'bytes=-{suffix}' if suffix else f'bytes={start}-{end}'
        last_error = None
        for attempt in range(6):
            try:
                req = urllib.request.Request(self.url, headers={
                    'Range': range_value, 'User-Agent': 'Mozilla/5.0',
                    'Accept-Encoding': 'identity'})
                with urllib.request.urlopen(req, timeout=30) as response:
                    content_range = response.headers.get('Content-Range', '')
                    match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', content_range)
                    if response.status != 206 or not match:
                        raise ValueError(f'Expected HTTP 206 + Content-Range, got {response.status} {content_range}')
                    actual_start, actual_end, total = map(int, match.groups())
                    if suffix:
                        if actual_end != total - 1 or actual_start != max(0, total-suffix):
                            raise ValueError(f'Unexpected suffix range: {content_range}')
                    elif actual_start != start or actual_end != end:
                        raise ValueError(f'Unexpected returned range: {content_range}')
                    if self.total is not None and total != self.total:
                        raise ValueError('Archive length changed between requests')
                    self.total = total
                    self.etag = response.headers.get('ETag')
                    self.last_modified = response.headers.get('Last-Modified')
                    data = response.read()
                    if len(data) != actual_end - actual_start + 1:
                        raise ValueError('Incomplete HTTP range response')
                    return data, actual_start
            except Exception as error:
                last_error = error
                print(f'  Retry {range_value} {attempt+1}: {type(error).__name__}: {error}', flush=True)
                time.sleep(min(1+attempt, 4))
        raise RuntimeError(f'{self.url}: failed {range_value}') from last_error

    def range(self, start, end):
        chunks = [(a, min(end, a+self.block_bytes-1))
                  for a in range(start, end+1, self.block_bytes)]
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as pool:
            data = list(pool.map(lambda pair: self._request(*pair)[0], chunks))
        return b''.join(data)

    def member(self, entry):
        offset = entry['local_offset']
        hraw = self.range(offset, offset+29)
        h = struct.unpack('<4s5H3L2H', hraw)
        if h[0] != b'PK\x03\x04' or h[3] != entry['method']:
            raise ValueError('ZIP local header signature/method mismatch')
        if h[2] & 1:
            raise ValueError('Encrypted ZIP member unsupported')
        start = offset + 30 + h[-2] + h[-1]
        payload = self.range(start, start+entry['compressed_bytes']-1)
        if entry['method'] == 8:
            raw = zlib.decompress(payload, -15)
        elif entry['method'] == 0:
            raw = payload
        else:
            raise ValueError(f"Unsupported compression method {entry['method']}")
        if len(raw) != entry['size_bytes'] or (zlib.crc32(raw) & 0xffffffff) != entry['crc32']:
            raise ValueError(f"{entry['member']}: size or CRC32 mismatch")
        return raw


def count_candidates(data, seen):
    ex = pd.ExcelFile(io.BytesIO(data))
    sheet = next(s for s in ex.sheet_names if s.startswith('Channel'))
    frame = ex.parse(sheet)
    required = ['Test_Time(s)', 'Date_Time', 'Cycle_Index', 'Current(A)',
                'Voltage(V)', 'Discharge_Capacity(Ah)']
    missing = set(required) - set(frame.columns)
    if missing:
        raise RuntimeError(f'Missing columns {missing}')
    for col in required:
        if col != 'Date_Time':
            frame[col] = pd.to_numeric(frame[col], errors='coerce')
    frame['Date_Time'] = pd.to_datetime(frame['Date_Time'], errors='coerce')
    bad_rows = int(frame[required].isna().any(axis=1).sum())
    frame = frame.dropna(subset=required)
    discharge = frame[frame['Current(A)'] < -0.1]
    count, duplicates = 0, 0
    dates = []
    for local_cycle, group in discharge.groupby('Cycle_Index', sort=True):
        group = group.sort_values('Test_Time(s)').drop_duplicates('Test_Time(s)')
        if len(group) < 2:
            continue
        key = (str(group['Date_Time'].iloc[0]), str(group['Date_Time'].iloc[-1]))
        if key in seen:
            duplicates += 1
            continue
        seen.add(key)
        count += 1
        dates.append({'local_cycle': int(local_cycle), 'date_start': key[0], 'date_end': key[1]})
    return {'sheet': sheet, 'valid_rows': int(len(frame)),
            'required_rows_dropped': bad_rows, 'candidate_cycles_in_file': count,
            'duplicate_candidates_in_file': duplicates, 'candidate_dates': dates}


def fetch_cell(number, root, workers):
    cell = f'CS2_{number}'
    url = f'https://web.calce.umd.edu/batteries/data/{cell}.zip'
    target = root / 'data' / cell
    target.mkdir(parents=True, exist_ok=True)
    print(f'{cell}: reading official remote ZIP central directory', flush=True)
    remote = RemoteZip(url, workers=workers)
    entries = sorted((e for e in remote.entries if e['member'].endswith('.xlsx')),
                     key=lambda e: date_key(e['member']))
    manifest = {'url': url, 'source_url': url,
                'retrieval_mode': 'partial ZIP entry retrieval using HTTP byte ranges',
                'whole_archive_downloaded': False,
                'whole_archive_sha256': None,
                'remote_archive_size_bytes': remote.total,
                'remote_archive_entry_count': len(remote.entries),
                'central_directory_sha256': remote.central_sha256,
                'etag': remote.etag, 'last_modified': remote.last_modified,
                'selection': 'Chronological date_key order; exact read_arbin candidate rules; stop at first workbook giving >=101 candidates',
                'retrieved_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'members': []}
    seen, total = set(), 0
    for entry in entries:
        path = target / Path(entry['member']).name
        if path.exists():
            raw = path.read_bytes()
            if len(raw) != entry['size_bytes'] or zlib.crc32(raw) & 0xffffffff != entry['crc32']:
                raise ValueError(f'{path}: existing file does not match official ZIP member CRC')
            print(f'{cell}: verified existing {path.name}', flush=True)
        else:
            print(f"{cell}: retrieving {entry['member']} ({entry['compressed_bytes']:,} ZIP bytes)", flush=True)
            raw = remote.member(entry)
            temporary = path.with_suffix('.xlsx.part')
            temporary.write_bytes(raw)
            temporary.replace(path)
        counts = count_candidates(raw, seen)
        total += counts['candidate_cycles_in_file']
        manifest['members'].append({'member': entry['member'], 'local_path': str(path),
                                    'size_bytes': len(raw), 'compressed_bytes': entry['compressed_bytes'],
                                    'sha256': sha(raw), 'crc32': f"{entry['crc32']:08x}",
                                    'crc_validation': 'passed against official ZIP central directory',
                                    **counts})
        manifest['candidate_cycles_total'] = total
        manifest['selection_complete'] = total >= 101
        (target / 'archive_source.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
        print(f'{cell}: {path.name} contributes {counts["candidate_cycles_in_file"]}; total {total}', flush=True)
        if total >= 101:
            break
    if total < 101:
        raise RuntimeError(f'{cell}: exhausted workbook entries with only {total} candidates')
    print(f'{cell}: READY, {len(manifest["members"])} workbooks, {total} candidates', flush=True)
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cells', nargs='+', type=int, required=True)
    parser.add_argument('--root', type=Path, default=Path('/Users/aaa/Documents/Codex/2026-10-12/pao/work/code2_run'))
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(args.cells)) as pool:
        results = list(pool.map(lambda n: fetch_cell(n, args.root, args.workers), args.cells))
    print(json.dumps([{'cell': n, 'members': len(r['members']), 'candidates': r['candidate_cycles_total']}
                      for n, r in zip(args.cells, results)], indent=2))


if __name__ == '__main__':
    main()
