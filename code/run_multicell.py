#!/usr/bin/env python3
"""Fixed-horizon replication of the supplied CALCE analysis on additional cells.

Raw workbook values are read only. Model implementation is in locked_core.py.
No cell, source file, cycle window, or parameter bound is selected on test error.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import math
import re
import shutil
import sys
import urllib.request
import zipfile
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
import locked_core as core

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(core.json_ready(obj), indent=2, ensure_ascii=False), encoding='utf-8')

def ci(values):
    return core.moving_block_bootstrap_mean_ci(np.asarray(values), samples=20000,
                                              block_length=5, seed=20260806)

def date_key(name):
    m = re.search(r'_(\d+)_(\d+)_(\d+)\.xlsx$', name)
    return datetime.strptime('_'.join(m.groups()), '%m_%d_%y')

def get_zip(cell, downloads):
    path = downloads / f'{cell}.zip'
    url = f'https://web.calce.umd.edu/batteries/data/{cell}.zip'
    if path.is_file() and zipfile.is_zipfile(path):
        return path, url
    downloads.mkdir(parents=True, exist_ok=True)
    for attempt in range(2):
        temporary = path.with_suffix('.part')
        with urllib.request.urlopen(url, timeout=60) as r, temporary.open('wb') as f:
            expected = r.headers.get('Content-Length')
            shutil.copyfileobj(r, f)
        if (expected is None or temporary.stat().st_size == int(expected)) and zipfile.is_zipfile(temporary):
            temporary.replace(path)
            return path, url
    raise RuntimeError(f'Incomplete download: {url}')

def read_arbin(number, downloads):
    cell = f'CS2_{number}'
    target = ROOT / 'data' / cell
    target.mkdir(parents=True, exist_ok=True)
    existing = sorted(target.glob('*.xlsx'), key=lambda p: date_key(p.name))
    candidates, sources, seen, duplicate_count = [], [], set(), 0
    archive = None
    if existing:
        entries = [(p.name, p.read_bytes()) for p in existing]
    else:
        archive_path, url = get_zip(cell, downloads)
        archive = zipfile.ZipFile(archive_path)
        names = sorted([n for n in archive.namelist() if n.endswith('.xlsx')], key=date_key)
        entries = ((n, archive.read(n)) for n in names)
        dump(target / 'archive_source.json', {'url': url, 'sha256': sha(archive_path),
                                             'size_bytes': archive_path.stat().st_size})
    for name, data in entries:
        local = target / Path(name).name
        if not local.exists():
            local.write_bytes(data)
        ex = pd.ExcelFile(io.BytesIO(data))
        sheet = next(s for s in ex.sheet_names if s.startswith('Channel'))
        frame = ex.parse(sheet)
        required = ['Test_Time(s)', 'Date_Time', 'Cycle_Index', 'Current(A)',
                    'Voltage(V)', 'Discharge_Capacity(Ah)']
        missing = set(required) - set(frame.columns)
        if missing:
            raise RuntimeError(f'{name}: missing columns {missing}')
        for col in required:
            if col != 'Date_Time':
                frame[col] = pd.to_numeric(frame[col], errors='coerce')
        frame['Date_Time'] = pd.to_datetime(frame['Date_Time'], errors='coerce')
        bad_rows = int(frame[required].isna().any(axis=1).sum())
        frame = frame.dropna(subset=required)
        discharge = frame[frame['Current(A)'] < -0.1]
        before = len(candidates)
        for local_cycle, group in discharge.groupby('Cycle_Index', sort=True):
            group = group.sort_values('Test_Time(s)').drop_duplicates('Test_Time(s)')
            if len(group) < 2:
                continue
            key = (str(group['Date_Time'].iloc[0]), str(group['Date_Time'].iloc[-1]))
            if key in seen:
                duplicate_count += 1
                continue
            seen.add(key)
            candidates.append({'source': local.name, 'local_cycle': int(local_cycle),
                               'date_start': key[0], 'date_end': key[1], 'group': group})
        sources.append({'file': local.name, 'sha256': sha(local), 'sheet': sheet,
                        'rows': int(len(frame)), 'required_rows_dropped': bad_rows,
                        'candidate_cycles_in_file': len(candidates) - before})
        print(f'  {cell}: {local.name}, {len(candidates)} candidate cycles', flush=True)
        if len(candidates) >= 101:
            break
    if archive:
        archive.close()
    candidates.sort(key=lambda c: pd.Timestamp(c['date_start']))
    profiles, rows = {}, []
    nominal_current = 0.55 if number in (33, 34) else 1.1
    rate = nominal_current / 1.1
    prior_end = None
    for cycle, candidate in enumerate(candidates[:101], 1):
        g = candidate.pop('group')
        t = g['Test_Time(s)'].to_numpy(float)
        t = t - t[0]
        current = -g['Current(A)'].to_numpy(float)
        voltage = g['Voltage(V)'].to_numpy(float)
        cutoff = core.first_crossing_time(t, voltage, 3.0)
        current_med = float(np.median(current))
        current_cv = float(np.std(current) / np.mean(current))
        reasons = []
        if len(t) < 80:
            reasons.append('fewer than 80 discharge samples')
        if t[-1] / 3600 < 0.75 / rate:
            reasons.append('duration below 0.75 / nominal_C_rate hours')
        counter_start = float(g['Discharge_Capacity(Ah)'].iloc[0])
        counter_end = float(g['Discharge_Capacity(Ah)'].iloc[-1])
        recorded_delta = counter_end - counter_start
        integrated_full = core.capacity_integral_ah(t, current)
        if recorded_delta < 0.88:
            reasons.append('recorded discharged capacity below 80% of 1.1 Ah')
        if abs(recorded_delta - integrated_full) > 0.002:
            reasons.append('cycle charge-counter increment and current integral differ by >2 mAh')
        if current_cv > 0.05:
            reasons.append('discharge-current coefficient of variation exceeds 5%')
        if abs(current_med / nominal_current - 1) > 0.1:
            reasons.append('median discharge current differs from nominal by more than 10%')
        if not math.isfinite(cutoff):
            reasons.append('no measured 3.0-V first crossing')
        start = pd.Timestamp(candidate['date_start'])
        if prior_end is not None and start < prior_end:
            raise RuntimeError(f'{cell}: overlapping discharge episodes, inspect source records')
        prior_end = pd.Timestamp(candidate['date_end'])
        capacity = math.nan
        if math.isfinite(cutoff):
            tc = np.append(t[t < cutoff - 1e-9], cutoff)
            ic = np.interp(tc, t, current)
            capacity = core.capacity_integral_ah(tc, ic)
        if not reasons:
            profiles[cycle] = core.Profile(cycle, t, current, voltage, cutoff, capacity, len(t))
        rows.append({**candidate, 'cycle': cycle, 'accepted': not reasons,
                     'reasons': '; '.join(reasons), 'samples': len(t),
                     'median_current_A': current_med, 'current_CV': current_cv,
                     'duration_h': t[-1] / 3600,
                     'recorded_discharge_counter_start_Ah': counter_start,
                     'recorded_discharge_counter_end_Ah': counter_end,
                     'recorded_discharge_Ah': recorded_delta,
                     'integrated_full_discharge_Ah': integrated_full,
                     'counter_minus_integral_mAh': (recorded_delta-integrated_full)*1000,
                     'capacity_to_3V_mAh': capacity * 1000,
                     'time_to_3V_min': cutoff / 60})
    audit = {'cell': cell.replace('_','-'), 'source_files': sources,
             'candidate_cycles_read': len(candidates), 'candidate_horizon': 101,
             'duplicates_removed': duplicate_count, 'accepted_cycles': len(profiles),
             'cycle_numbering': 'chronological discharge ordinal before QC, first 101 candidates',
             'nominal_current_A': nominal_current, 'nominal_rate_C': rate,
             'rows': rows}
    return profiles, audit

def evaluate(profile, parameters, relation, window, setting):
    q = float(core.predict_capacity_ah(profile.cycle, parameters))
    td = np.arange(0., profile.time_s[-1], 1.)
    if len(td) == 0 or td[-1] < profile.time_s[-1]:
        td = np.append(td, profile.time_s[-1])
    zd = 1 - cumulative_trapezoid(np.interp(td, profile.time_s, profile.current_a), td, initial=0) / (3600*q)
    ts = core.first_crossing_time(td, zd, 0.)
    tv = core.first_crossing_time(td, relation(zd), 3.)
    endpoint, event = core.classify_event(tv, ts)
    end = min(profile.measured_cutoff_s, endpoint) if math.isfinite(endpoint) else profile.measured_cutoff_s
    tc, ic, vc = core.arrays_to_time(profile, end)
    discharged = cumulative_trapezoid(ic, tc, initial=0)
    zc = 1 - discharged / (3600*q)
    err = (relation(zc) - vc) * 1000
    observed_soc = 1 - discharged / (3600*profile.capacity_to_cutoff_ah)
    diagnostic = (relation(observed_soc) - vc) * 1000
    t, current, measured = core.arrays_to_cutoff(profile)
    z = 1 - cumulative_trapezoid(current, t, initial=0)/(3600*q)
    prediction = relation(z)
    ef = (prediction - measured)*1000
    row = {'cycle':profile.cycle, 'window':window,'setting':setting,
           'voltage_MAE_mV':float(np.abs(err).mean()),
           'voltage_RMSE_mV':float(np.sqrt(np.mean(err**2))),
           'voltage_bias_mV':float(err.mean()), 'voltage_points':len(err),
           'sum_absolute_error_mV':float(np.abs(err).sum()),
           'sum_squared_error_mV2':float(np.square(err).sum()), 'sum_error_mV':float(err.sum()),
           'full_support_MAE_mV':float(np.abs(ef).mean()),
           'mapping_only_MAE_mV_diagnostic':float(np.abs(diagnostic).mean()),
           'sum_mapping_only_absolute_error_mV':float(np.abs(diagnostic).sum()),
           'common_support_fraction':end/profile.measured_cutoff_s,
           'observed_endpoint_min':profile.measured_cutoff_s/60,
           'predicted_endpoint_min':endpoint/60,
           'endpoint_error_min':(endpoint-profile.measured_cutoff_s)/60,
           'endpoint_evaluable':math.isfinite(endpoint), 'endpoint_event':event}
    long = pd.DataFrame({'cycle':profile.cycle,'window':window,'setting':setting,
                         'time_s':t,'current_A':current,'measured_voltage_V':measured,
                         'predicted_voltage_V':prediction,'predicted_SOC':z,
                         'within_common_support':t <= end+1e-9})
    return row,long

def voltage_summary(frame):
    rows=[]
    for (window,setting),g in frame.groupby(['window','setting'],sort=False):
        count = g['voltage_points'].sum()
        c=ci(g['voltage_MAE_mV'])
        ep=g['endpoint_error_min'].dropna()
        rows.append({'window':window,'setting':setting,'cycles':len(g),
                     'pooled_voltage_MAE_mV':g['sum_absolute_error_mV'].sum()/count,
                     'pooled_voltage_RMSE_mV':np.sqrt(g['sum_squared_error_mV2'].sum()/count),
                     'pooled_voltage_bias_mV':g['sum_error_mV'].sum()/count,
                     'cycle_mean_voltage_MAE_mV':c['mean'],
                     'cycle_mean_MAE_CI95_low_mV':c['ci95_low'],
                     'cycle_mean_MAE_CI95_high_mV':c['ci95_high'],
                     'mapping_only_MAE_mV_diagnostic':g['sum_mapping_only_absolute_error_mV'].sum()/count,
                     'mean_common_support_fraction':g['common_support_fraction'].mean(),
                     'endpoint_evaluable_count':len(ep),
                     'endpoint_MAE_min_diagnostic':np.abs(ep).mean()})
    return pd.DataFrame(rows)

def run_cell(number, downloads, source_relation):
    label=f'CS2-{number}'
    out=ROOT/'results'/label
    out.mkdir(parents=True,exist_ok=True)
    if number in (8,21):
        raw=ROOT/'data'/f'CS2_{number}_6_8_10.txt'
        profiles,audit=core.read_profiles(raw)
        profiles={n:p for n,p in profiles.items() if n<=101}
        audit['source']={'file':raw.name,'sha256':sha(raw)}
        audit['cycle_numbering']='original CADEX discharge count within fixed June 8 record'
    else:
        profiles,audit=read_arbin(number,downloads)
        pd.DataFrame(audit['rows']).to_csv(out/'cycle_quality_control.csv',index=False)
    dump(out/'extraction_audit.json',audit)
    trajectory_rows=[]
    for n,p in profiles.items():
        trajectory_rows.append(pd.DataFrame({'cycle':n,'time_s':p.time_s,'current_A':p.current_a,'voltage_V':p.voltage_v}))
    pd.concat(trajectory_rows).to_csv(out/'accepted_discharge_profiles.csv.gz',index=False,compression='gzip')
    inventory=pd.DataFrame([{'cycle':n,'capacity_mAh':p.capacity_to_cutoff_ah*1000,
                             'observed_endpoint_min':p.measured_cutoff_s/60,
                             'median_current_A':np.median(p.current_a)} for n,p in profiles.items()])
    inventory.to_csv(out/'observed_capacity.csv',index=False)
    windows={'fit':[p for n,p in profiles.items() if n<=40],
             'dev':[p for n,p in profiles.items() if n<=60],
             'validation_41_60':[p for n,p in profiles.items() if 40<n<=60],
             'robustness_61_80':[p for n,p in profiles.items() if 60<n<=80],
             'confirmation_81_101':[p for n,p in profiles.items() if 80<n<=101]}
    if any(len(v)==0 for v in windows.values()):
        raise RuntimeError('One or more required chronological windows is empty')
    params={};relations={};failures=[]
    for key in ['fit','dev']:
        relations[key]=core.build_empirical_relation(windows[key])
        try:params[key]=core.fit_saturating_capacity(windows[key])
        except core.ExperimentError as e:failures.append(f'{key}: {e}')
    dump(out/'model_parameters.json',{'parameters':params,'failures':failures,
                                     'cycle_windows':{k:[p.cycle for p in v] for k,v in windows.items()}})
    pd.DataFrame({'SOC':relations['dev'].soc_grid,'voltage_V':relations['dev'].voltage_grid}).to_csv(out/'frozen_voltage_relation.csv',index=False)
    if 'dev' not in params:
        return {'cell':label,'status':'frozen_capacity_fit_failed','accepted_cycles':len(profiles),'error':'; '.join(failures)}
    frozen=windows['dev'];test=windows['robustness_61_80']+windows['confirmation_81_101']
    recent=float(np.median([p.capacity_to_cutoff_ah*1000 for p in frozen[-10:]]))
    linear=np.polyfit([p.cycle for p in frozen],[p.capacity_to_cutoff_ah*1000 for p in frozen],1)
    capacity=[]
    for p in test:
        obs=p.capacity_to_cutoff_ah*1000
        capacity.append({'cell':label,'cycle':p.cycle,'observed_mAh':obs,
                         'recent_state_mAh':recent,'saturating_mAh':float(core.predict_capacity_ah(p.cycle,params['dev']))*1000,
                         'linear_mAh':float(np.polyval(linear,p.cycle))})
    cap=pd.DataFrame(capacity)
    cap.to_csv(out/'capacity_predictions.csv',index=False)
    caprows=[]
    for model in ['recent_state','linear','saturating']:
        error=cap[f'{model}_mAh']-cap['observed_mAh']
        caprows.append({'model':model,'test_cycles':len(cap),'MAE_mAh':np.abs(error).mean(),
                        'RMSE_mAh':np.sqrt(np.mean(error**2)),'bias_mAh':error.mean()})
    pd.DataFrame(caprows).to_csv(out/'capacity_summary.csv',index=False)
    gain=np.abs(cap['recent_state_mAh']-cap['observed_mAh'])-np.abs(cap['saturating_mAh']-cap['observed_mAh'])
    gainci=ci(gain)
    dump(out/'capacity_paired_gain_CI.json',gainci)
    all_rows=[];all_long=[]
    for window in ['validation_41_60','robustness_61_80','confirmation_81_101']:
        stage='fit' if window.startswith('validation') else 'dev'
        if stage not in params:continue
        settings=[('cell_specific',relations[stage])]
        if window.startswith('confirmation'):
            settings.append(('CS2_21_voltage_map_transfer',source_relation))
        for setting,relation in settings:
            for p in windows[window]:
                row,long=evaluate(p,params[stage],relation,window,setting)
                all_rows.append(row);all_long.append(long)
    bycycle=pd.DataFrame(all_rows)
    bycycle.to_csv(out/'voltage_by_cycle.csv',index=False)
    pd.concat(all_long).to_csv(out/'voltage_predictions.csv.gz',index=False,compression='gzip')
    summary=voltage_summary(bycycle)
    summary.to_csv(out/'voltage_summary.csv',index=False)
    final=summary[(summary['window']=='confirmation_81_101')&(summary['setting']=='cell_specific')].iloc[0]
    transfer=summary[(summary['window']=='confirmation_81_101')&(summary['setting']=='CS2_21_voltage_map_transfer')].iloc[0]
    result={'cell':label,'status':'complete' if not failures else 'development_complete_validation_fit_failed',
            'accepted_cycles':len(profiles),'development_cycles':len(frozen),
            'capacity_test_cycles':len(test),'confirmation_cycles':int(final['cycles']),
            'measured_median_current_A':inventory['median_current_A'].median(),
            'recent_state_capacity_MAE_mAh':caprows[0]['MAE_mAh'],
            'linear_capacity_MAE_mAh':caprows[1]['MAE_mAh'],
            'saturating_capacity_MAE_mAh':caprows[2]['MAE_mAh'],
            'capacity_MAE_reduction_percent':100*(1-caprows[2]['MAE_mAh']/caprows[0]['MAE_mAh']),
            'capacity_paired_gain_CI95_low_mAh':gainci['ci95_low'],
            'capacity_paired_gain_CI95_high_mAh':gainci['ci95_high'],
            **{k:final[k] for k in final.index if k not in ['window','setting','cycles']},
            'CS2_21_transfer_MAE_mV':transfer['pooled_voltage_MAE_mV'],
            'error':'; '.join(failures)}
    print(label,json.dumps(core.json_ready(result)),flush=True)
    return result

def main():
    a=argparse.ArgumentParser()
    a.add_argument('--downloads',type=Path,default=ROOT/'downloads')
    a.add_argument('--cells',nargs='+',type=int,default=[21,8,33,34,35,36,37,38])
    args=a.parse_args()
    p,_=core.read_profiles(ROOT/'data/CS2_21_6_8_10.txt')
    source=core.build_empirical_relation([v for n,v in p.items() if n<=60])
    summary_path=ROOT/'results'/'multicell_summary.csv'
    selected_labels={f'CS2-{n}' for n in args.cells}
    results=[]
    if summary_path.is_file():
        previous=pd.read_csv(summary_path)
        results=previous[~previous['cell'].isin(selected_labels)].to_dict('records')
    for number in args.cells:
        try:results.append(run_cell(number,args.downloads,source))
        except Exception as e:
            import traceback
            traceback.print_exc()
            results.append({'cell':f'CS2-{number}','status':'failed','error':str(e)})
        pd.DataFrame(results).to_csv(ROOT/'results'/'multicell_summary.csv',index=False)
    dump(ROOT/'results'/'runtime_versions.json',{'python':sys.version,'numpy':np.__version__,
         'pandas':pd.__version__,'core_sha256':sha(Path(core.__file__)),
         'protocol_sha256':sha(ROOT/'protocol_before_run.json'),
         'execution_utc':datetime.now().isoformat()})

if __name__=='__main__':main()
