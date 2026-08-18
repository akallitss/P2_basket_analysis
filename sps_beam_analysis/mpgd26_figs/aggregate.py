#!/usr/bin/env python3
"""Aggregate every per-run product fetched from EOS into tidy report datasets.

Inputs (under scratchpad):
  products/analysis/<det>/<run>/scan*/22_tag_probe_efficiency/tag_probe_efficiency_spark_vetoed.csv
  products/analysis/<det>/<run>/scan*/26_hv_spark_qa/spark_qa_<det>.csv
  products/analysis/<det>/<run>/scan*/30_raw_stream_efficiency/raw_stream_efficiency_<det>.csv
  products/analysis/telescope/<run>/scan*/29_waveform_timing/timing_vs_{mesh,drift}.csv
  products/analysis/telescope/<run>/<sub>/29_waveform_timing/waveform_timing_summary.json
  products/analysis/vmm/<run>/<sub>/<capture>/scalars.json
  eos_inventory/hv_setpoints.csv        (per run/sub_run/det mesh+drift setpoints)
  products/hv/<run>/<sub_run>/hv_monitor.csv

Outputs in report_data/.
"""
import os, re, csv, json, glob, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import S, A, RD as OUT  # noqa: E402

DETS = ['P2_IN', 'P2_MID', 'P2_OUT']

# Gas assignment: verified changeover window = run_60 overnight (P2 HV off,
# Aug 1 21:19 -> Aug 2 morning). Named July runs + run_21..58 = old gas.
def gas_of(run):
    m = re.match(r'run_(\d+)$', run)
    if m:
        n = int(m.group(1))
        if n == 1:
            return 'Ar/CO2/Iso 93/5/2'
        return 'Ar/CO2/Iso 93/5/2' if n <= 58 else 'Ar/CF4/Iso 88/10/2'
    return 'Ar/CO2/Iso 93/5/2'

# ---------------------------------------------------------------- tag-probe --
rows = []
for f in glob.glob(f'{A}/*/*/scan*/22_tag_probe_efficiency/tag_probe_efficiency*.csv'):
    parts = f[len(A) + 1:].split('/')
    det, run, prod = parts[0], parts[1], parts[2]
    if det == 'telescope':
        continue
    try:
        df = pd.read_csv(f)
    except Exception:
        continue
    df['det'], df['run'], df['prod_sub'] = det, run, prod
    df['vetoed'] = 'spark_vetoed' in os.path.basename(f)
    rows.append(df)
if rows:
    tp = pd.concat(rows, ignore_index=True)
    tp['gas'] = tp['run'].map(gas_of)
    tp.to_csv(f'{OUT}/dream_tag_probe.csv', index=False)
    print(f'dream_tag_probe.csv  {len(tp)} rows, runs={tp.run.nunique()}')

# --------------------------------------------------------------- raw stream --
rows = []
for f in glob.glob(f'{A}/*/*/scan*/30_raw_stream_efficiency/raw_stream_efficiency_*.csv'):
    parts = f[len(A) + 1:].split('/')
    det, run = parts[0], parts[1]
    try:
        df = pd.read_csv(f)
    except Exception:
        continue
    df['det'], df['run'] = det, run
    rows.append(df)
if rows:
    rs = pd.concat(rows, ignore_index=True)
    rs.to_csv(f'{OUT}/dream_raw_stream.csv', index=False)
    print(f'dream_raw_stream.csv {len(rs)} rows, runs={rs.run.nunique()}')

# ------------------------------------------------------------------- sparks --
rows = []
for f in glob.glob(f'{A}/*/*/scan*/26_hv_spark_qa/spark_qa_*.csv'):
    parts = f[len(A) + 1:].split('/')
    det, run = parts[0], parts[1]
    try:
        df = pd.read_csv(f)
    except Exception:
        continue
    df['det'], df['run'] = det, run
    rows.append(df)
if rows:
    sp = pd.concat(rows, ignore_index=True)
    sp['gas'] = sp['run'].map(gas_of)
    sp.to_csv(f'{OUT}/dream_sparks.csv', index=False)
    print(f'dream_sparks.csv     {len(sp)} rows, runs={sp.run.nunique()}')

# ------------------------------------------------------------ timing scans --
rows = []
for f in glob.glob(f'{A}/telescope/*/scan*/29_waveform_timing/timing_vs_*.csv'):
    parts = f[len(A) + 1:].split('/')
    run = parts[1]
    axis = 'mesh' if 'mesh' in os.path.basename(f) else 'drift'
    try:
        df = pd.read_csv(f)
    except Exception:
        continue
    df['run'], df['axis'] = run, axis
    rows.append(df)
if rows:
    tm = pd.concat(rows, ignore_index=True)
    tm.to_csv(f'{OUT}/dream_timing_scans.csv', index=False)
    print(f'dream_timing_scans.csv {len(tm)} rows, runs={tm.run.nunique()}')

# --------------------------------------------- per-sub_run timing summaries --
rows = []
for f in glob.glob(f'{A}/telescope/*/*/29_waveform_timing/waveform_timing_summary.json'):
    parts = f[len(A) + 1:].split('/')
    run, sub = parts[1], parts[2]
    if sub.startswith('scan'):
        continue
    try:
        d = json.load(open(f))
    except Exception:
        continue
    for st in (d.get('stations') or []):
        det = st.get('detector')
        best = st.get('best_algorithm')
        alg = (st.get('algorithms') or {}).get(best, {})
        rows.append(dict(run=run, sub_run=sub, det=det, kind='station',
                         algo=best, n=st.get('n_benchmark'),
                         sigma_ns=alg.get('sigma_walk_ns',
                                          alg.get('sigma_ftst_ns'))))
    for pr in (d.get('pairs') or []):
        rows.append(dict(run=run, sub_run=sub, det=pr.get('pair'),
                         kind='pair', algo='', n=pr.get('n_events'),
                         sigma_ns=pr.get('sigma_single_ns')))
if rows:
    ts = pd.DataFrame(rows)
    ts.to_csv(f'{OUT}/dream_timing_persubrun.csv', index=False)
    print(f'dream_timing_persubrun.csv {len(ts)} rows, runs={ts.run.nunique()}')
else:
    print('dream_timing_persubrun: no rows -- check summary JSON layout')

# ---------------------------------------------------------------- VMM scans --
# HV setpoints from the VMM DAQ's own run_config.json (the DREAM and VMM DAQs
# number their runs independently after run_58, so the DREAM configs must not
# be used to label VMM data).
P2_HV = {'P2_IN': ('0', '1'), 'P2_MID': ('2', '3'), 'P2_OUT': ('4', '5')}
hvmap = {}
for f in glob.glob(f'{S}/eos_inventory/vmm_meta/runs/run_*/run_config.json'):
    m = re.search(r'/(run_\d+)/run_config\.json$', f)
    if not m:
        continue
    run = m.group(1)
    try:
        d = json.load(open(f))
    except Exception:
        continue
    for sr in d.get('sub_runs', []):
        hvs = (sr.get('hvs') or {}).get('8', {})
        for det, (dch, mch) in P2_HV.items():
            try:
                hvmap[(run, sr['sub_run_name'], det)] = (
                    float(hvs[mch]), float(hvs[dch]))
            except (KeyError, TypeError, ValueError):
                pass
rows = []
for f in glob.glob(f'{A}/vmm/*/*/*/scalars.json'):
    parts = f[len(A) + 1:].split('/')
    run, sub, cap = parts[1], parts[2], parts[3]
    try:
        d = json.load(open(f))
    except Exception:
        continue
    eff = d.get('efficiency') or {}
    for st, v in eff.items():
        if not isinstance(v, dict):
            continue
        mv, dv = hvmap.get((run, sub, st), ('', ''))
        rows.append(dict(run=run, sub_run=sub, capture=cap, det=st,
                         mesh_v=mv, drift_v=dv, gas=gas_of(run),
                         n_hits=d.get('n_hits'),
                         n_triggers=v.get('n_triggers'),
                         n_signal=v.get('n_signal'),
                         n_sideband=v.get('n_sideband'),
                         raw_eff=v.get('raw_efficiency'),
                         acc_eff=v.get('accidental_efficiency'),
                         eff=v.get('efficiency'),
                         eff_corr=v.get('efficiency_corrected'),
                         eff_lo=v.get('efficiency_lo'),
                         eff_hi=v.get('efficiency_hi'),
                         mu_ns=v.get('mu_ns'), sigma_ns=v.get('sigma_ns'),
                         contrast=v.get('contrast')))
vm = pd.DataFrame(rows)
if len(vm):
    vm.to_csv(f'{OUT}/vmm_captures.csv', index=False)
    print(f'vmm_captures.csv     {len(vm)} rows, runs={vm.run.nunique()}')

    # per-sub_run aggregate: trigger-weighted, sideband re-subtracted
    def agg(g):
        n_trig = g['n_triggers'].sum()
        n_sig = g['n_signal'].sum()
        # per-capture sideband scale = acc_eff * n_trig / n_sideband
        with np.errstate(divide='ignore', invalid='ignore'):
            scale = (g['acc_eff'] * g['n_triggers'] / g['n_sideband']).replace(
                [np.inf, -np.inf], np.nan).median()
        n_acc = (g['n_sideband'] * (scale if np.isfinite(scale) else 0)).sum()
        eff = (n_sig - n_acc) / n_trig if n_trig else np.nan
        err = np.sqrt(max(n_sig, 1)) / n_trig if n_trig else np.nan
        w = g['n_triggers'].clip(lower=1)
        return pd.Series(dict(
            n_captures=len(g), n_triggers=n_trig,
            eff=eff, eff_err=err,
            raw_eff=n_sig / n_trig if n_trig else np.nan,
            mu_ns=np.average(g['mu_ns'], weights=w),
            sigma_ns=np.average(g['sigma_ns'], weights=w),
            mesh_v=g['mesh_v'].iloc[0], drift_v=g['drift_v'].iloc[0],
            gas=g['gas'].iloc[0]))
    vs = (vm.dropna(subset=['n_triggers'])
            .groupby(['run', 'sub_run', 'det'], as_index=False)
            .apply(agg, include_groups=False))
    vs.to_csv(f'{OUT}/vmm_subruns.csv', index=False)
    print(f'vmm_subruns.csv      {len(vs)} rows')

# ------------------------------------------------------------- HV campaign --
HVD = os.path.join(S, 'products', 'hv')
CHAN = {'8:0': ('P2_IN', 'drift'), '8:1': ('P2_IN', 'mesh'),
        '8:2': ('P2_MID', 'drift'), '8:3': ('P2_MID', 'mesh'),
        '8:4': ('P2_OUT', 'drift'), '8:5': ('P2_OUT', 'mesh'),
        '8:6': ('uRW_front', 'drift'), '8:7': ('uRW_back', 'drift'),
        '12:0': ('uRW_front', 'resist'), '12:1': ('uRW_back', 'resist')}
rows = []
for f in glob.glob(f'{HVD}/*/*/*/hv_monitor.csv') + glob.glob(f'{HVD}/*/*/hv_monitor.csv'):
    rel = os.path.relpath(f, HVD).split('/')
    rel = [p for p in rel if p != '.']
    run, sub = rel[-3], rel[-2]
    try:
        df = pd.read_csv(f)
    except Exception:
        continue
    if 'timestamp' not in df.columns or not len(df):
        continue
    t = pd.to_datetime(df['timestamp'], errors='coerce')
    df = df[t.notna()]; t = t[t.notna()]
    # downsample to 30 s medians + per-window imon max (sparks live in the max)
    df.index = t
    for ch, (det, electrode) in CHAN.items():
        vcol, icol = f'{ch} vmon', f'{ch} imon'
        if vcol not in df.columns:
            continue
        g = df[[vcol, icol]].resample('30s')
        med = g.median()
        imax = df[icol].resample('30s').max()
        for tt, r in med.iterrows():
            rows.append((run, sub, det, electrode, tt.isoformat(),
                         r[vcol], r[icol], imax.loc[tt]))
hvc = pd.DataFrame(rows, columns=['run', 'sub_run', 'det', 'electrode', 'time',
                                  'vmon', 'imon_med', 'imon_max'])
hvc.to_csv(f'{OUT}/hv_campaign_30s.csv', index=False)
print(f'hv_campaign_30s.csv  {len(hvc)} rows, runs={hvc.run.nunique()}')
print('done')
