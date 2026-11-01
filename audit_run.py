import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE / 'code2_run' if (HERE / 'code2_run').is_dir() else HERE
sys.path.insert(0, str(ROOT / 'code'))
import locked_core as core


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a, b, label, atol=1e-7):
    np.testing.assert_allclose(a, b, rtol=1e-9, atol=atol, err_msg=label, equal_nan=True)


def main():
    summary = pd.read_csv(ROOT / 'results/multicell_summary.csv')
    assert summary.cell.tolist() == ['CS2-21', 'CS2-8', 'CS2-33', 'CS2-34', 'CS2-35', 'CS2-36', 'CS2-37', 'CS2-38']
    assert summary.status.eq('complete').all(), summary[['cell', 'status', 'error']]
    supplied = Path('/Users/aaa/Library/Containers/com.tencent.xinWeChat/Data/Documents/xwechat_files/wxid_5b3vqyr9aeq222_20fd/msg/file/2026-10/code 2')
    hashes = {p.name: sha(p) for p in (ROOT / 'code').glob('*.py')}
    assert len(hashes) == 4
    if supplied.is_dir():
        expected_hashes = {name: sha(supplied / name) for name in hashes}
    else:
        expected_hashes = json.loads((ROOT / 'results/rerun_validation.json').read_text())['source_sha256']
    for name, digest in hashes.items():
        assert digest == expected_hashes[name], name
    checks = []
    source_profiles, _ = core.read_profiles(ROOT / 'data/CS2_21_6_8_10.txt')
    source_relation = core.build_empirical_relation([p for n, p in source_profiles.items() if n <= 60])
    for _, row in summary.iterrows():
        cell = row.cell
        out = ROOT / 'results' / cell
        cap = pd.read_csv(out / 'capacity_predictions.csv')
        for model in ['recent_state', 'linear', 'saturating']:
            value = np.abs(cap[model + '_mAh'] - cap.observed_mAh).mean()
            close(value, row[model + '_capacity_MAE_mAh'], cell + ' capacity ' + model)
        bycycle = pd.read_csv(out / 'voltage_by_cycle.csv')
        final = bycycle[(bycycle.window == 'confirmation_81_101') & (bycycle.setting == 'cell_specific')]
        close(final.sum_absolute_error_mV.sum() / final.voltage_points.sum(), row.pooled_voltage_MAE_mV, cell + ' pooled MAE')
        close(np.sqrt(final.sum_squared_error_mV2.sum() / final.voltage_points.sum()), row.pooled_voltage_RMSE_mV, cell + ' pooled RMSE')
        close(final.voltage_MAE_mV.mean(), row.cycle_mean_voltage_MAE_mV, cell + ' cycle-mean MAE')
        assert len(final) == row.confirmation_cycles
        if cell in ['CS2-21', 'CS2-8']:
            core_summary = pd.read_csv(ROOT / 'core_results' / cell / 'Table3_terminal_voltage_reconstruction_summary.csv')
            core_final = core_summary[core_summary.window == 'final_confirmation_81_101'].iloc[0]
            close(core_final.common_support_voltage_MAE_mV, row.pooled_voltage_MAE_mV, cell + ' standalone core MAE')
            profiles, _ = core.read_profiles(ROOT / 'data' / (cell.replace('-', '_') + '_6_8_10.txt'))
        else:
            qc = pd.read_csv(out / 'cycle_quality_control.csv')
            assert len(qc) == 101
            assert qc.cycle.tolist() == list(range(1, 102))
            assert int(qc.accepted.sum()) == row.accepted_cycles
            raw_profiles = pd.read_csv(out / 'accepted_discharge_profiles.csv.gz')
            profiles = {}
            for n, g in raw_profiles.groupby('cycle'):
                t = g.time_s.to_numpy()
                i = g.current_A.to_numpy()
                v = g.voltage_V.to_numpy()
                endpoint = core.first_crossing_time(t, v, 3.0)
                keep = t < endpoint - 1e-9
                tc = np.append(t[keep], endpoint)
                ic = np.append(i[keep], np.interp(endpoint, t, i))
                profiles[int(n)] = core.Profile(int(n), t, i, v, endpoint, core.capacity_integral_ah(tc, ic), len(t))
        params = json.loads((out / 'model_parameters.json').read_text())['parameters']['dev']
        relation_df = pd.read_csv(out / 'frozen_voltage_relation.csv')
        relation = core.EmpiricalRelation(relation_df.SOC.to_numpy(), relation_df.voltage_V.to_numpy(), 'saved frozen relation')
        chosen = [profiles[int(n)] for n in final.cycle]
        # Compare the batch implementation with the original core on complete final windows.
        for setting, rel in [('cell_specific', relation), ('CS2_21_voltage_map_transfer', source_relation)]:
            metrics, _, _, _ = core.evaluate_window('confirmation_81_101', chosen, params, rel)
            expected = bycycle[(bycycle.window == 'confirmation_81_101') & (bycycle.setting == setting)].set_index('cycle')
            for _, actual in metrics.iterrows():
                saved = expected.loc[int(actual.cycle)]
                for key in ['voltage_MAE_mV', 'voltage_RMSE_mV', 'voltage_bias_mV']:
                    close(actual[key], saved[key], cell + ' original core ' + setting + ' ' + key)
                close(actual.endpoint_time_error_min, saved.endpoint_error_min, cell + ' endpoint error')
        checks.append({'cell': cell, 'final_cycles_checked': len(chosen), 'settings_checked': 2, 'capacity_and_voltage_summary_recomputed': True, 'original_core_comparison_passed': True})
    figures = sorted((ROOT / 'figures').glob('*'))
    assert len(figures) == 12
    assert all(p.stat().st_size > 1000 for p in figures)
    for p in figures:
        data = p.read_bytes()
        if p.suffix == '.png':
            assert data.startswith(b'\x89PNG\r\n\x1a\n')
        elif p.suffix == '.pdf':
            assert data.startswith(b'%PDF-') and b'%%EOF' in data[-1024:]
        elif p.suffix == '.svg':
            assert b'<svg' in data[:1000]
    audit = {'passed': True, 'supplied_source_files_unchanged': True, 'source_sha256': hashes, 'cell_checks': checks, 'figure_files': len(figures)}
    (ROOT / 'results/rerun_validation.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2))
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
