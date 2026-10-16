#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
32_low_field_timing.py

Why does the P2 timing degrade far more than Magboltz predicts below
~375 V/cm drift field? (Test B of HANDOFF_gap_timing_2026-10-15.md.)

Two hypotheses to tell apart with data:
  charging   a surface charging up with ions screens the low drift fields ->
             the result depends on the scan DIRECTION (up vs down in field)
             and drifts with TIME inside a sub_run after an HV step
  static     a field/geometry effect -> same result for the same HV,
             whatever came before, and flat in time

Two modes:

  extract  (on lxplus, reads the nTOF EOS tree directly). For every trigger
           of every sub_run listed in RUNS, the leading hit of each P2
           station (FEU 3/4/5 = P2_IN/MID/OUT): the largest hit with
           amplitude in [AMP_FLOOR, AMP_MAX], 0 <= time <= window, not
           saturated or truncated. `time` is the 30 % leading edge,
           ftst-corrected [ns]. One compressed table per sub_run, triggers
           with a good hit in at least two stations:
             <out>/<run>/<sub_run>.npz   eventId, t_s (s since the first
             trigger), per station time / amp (NaN = no good hit), and the
             sub_run's n_triggers / span_s

  analyse  (local, on the extracted tables). Per sub_run and station:
           median time and 68 % half-width; single-chamber resolution from
           (MID - OUT)/sqrt(2); hit fraction (triggers with a good hit).
           The same in time slices inside each sub_run. Then the up- vs
           down-scan comparison at matching drift field and mesh voltage.

The noise pitfall of 28_timing_qa (mean_ns without an amplitude cut pulled
toward the window centre at low field) is avoided by the amplitude window
and by medians/quantiles throughout.

Usage:
  python3 32_low_field_timing.py extract --out DIR [--run RUN ...] [--jobs N]
  python3 32_low_field_timing.py analyse --data DIR --out DIR
"""

import argparse
import glob
import os
import sys

import numpy as np

EOS_RUNS = '/eos/experiment/ntof/data/x17/p2_sps_july/runs'
STATIONS = {3: 'P2_IN', 4: 'P2_MID', 5: 'P2_OUT'}
AMP_MAX = 3500.0
# extraction floor; analyse applies the physics cut. 200 ADC (the sim session's
# cut) keeps only the high tail at low drift field, where the median amplitude
# falls to ~150 ADC; the MID-OUT time coincidence is what removes the noise.
AMP_FLOOR = 60.0
T_WINDOW = 16 * 60.0                   # 16 samples x 60 ns
DRIFT_GAP_CM = 0.4

# Runs with drift-field points, Ar/CO2/iso 93/5/2 (all <= run_58).
#   up   = field increasing through the run, down = decreasing
RUNS = {
    'drift_scan_1':         'up',      # 07-24 12:28, mesh 450, 0 -> 750 V/cm
    'drift_scan_2':         'up',      # 07-24 14:37, 875 -> 1125
    'drift_scan_final':     'up',      # 07-24 23:43, drift_450 only (0 V/cm)
    'drift_mesh_scan_1':    'up',      # 07-25, mesh 450, 125 -> 1125 (+ mesh scan)
    'drift_mesh_2d_1':      'up',      # 07-25 eve, 25 / 150 / 275 V/cm rows
    'mesh_drift_scan_up_1': 'down',    # 07-26 eve, drift2_*: 1062 -> 25 V/cm
    'drift_mesh_2d_2':      'mixed',   # 07-27, rows alternate up / down
    'p2_mesh_drift_eff_1':  'up',      # 07-28, drift_*: 375 -> 1125
}


# --------------------------------------------------------------------------- #
# extract (lxplus)
# --------------------------------------------------------------------------- #
def leading_hits(files, amp_min=AMP_FLOOR):
    """Per trigger: eventId, trigger time, and per station the largest good hit."""
    import uproot
    br = ['eventId', 'trigger_timestamp_ns', 'feu', 'amplitude', 'time',
          'saturated', 'trunc_left', 'trunc_right']
    ev_all, ts_all, good = [], [], []
    for fn in files:
        for a in uproot.open(fn)['hits'].iterate(br, step_size='200 MB', library='np'):
            # every trigger, whatever it holds (denominator of the hit fraction)
            u, i = np.unique(a['eventId'], return_index=True)
            ev_all.append(u)
            ts_all.append(a['trigger_timestamp_ns'][i].astype(np.int64))
            ok = (np.isin(a['feu'], list(STATIONS))
                  & (a['amplitude'] >= amp_min) & (a['amplitude'] <= AMP_MAX)
                  & (a['time'] >= 0) & (a['time'] <= T_WINDOW)
                  & ~a['saturated'] & ~a['trunc_left'] & ~a['trunc_right'])
            good.append(np.stack([a['eventId'][ok].astype(np.float64),
                                  a['feu'][ok].astype(np.float64),
                                  a['time'][ok], a['amplitude'][ok]], 1))
    ev = np.concatenate(ev_all)
    ts = np.concatenate(ts_all)
    ev, i = np.unique(ev, return_index=True)       # an event can span chunks
    ts = ts[i]
    g = np.concatenate(good) if good else np.zeros((0, 4))
    out = {'eventId': ev.astype(np.int64),
           't_s': ((ts - ts.min()) / 1e9).astype(np.float32)}
    for feu, name in STATIONS.items():
        x = g[g[:, 1] == feu]
        x = x[np.lexsort((-x[:, 3], x[:, 0]))]      # largest amplitude first
        x = x[np.r_[True, x[1:, 0] != x[:-1, 0]]] if len(x) else x
        t = np.full(len(ev), np.nan, np.float32)
        amp = np.full(len(ev), np.nan, np.float32)
        j = np.searchsorted(ev, x[:, 0].astype(np.int64))
        t[j], amp[j] = x[:, 2], x[:, 3]
        out[f'{name}_time'], out[f'{name}_amp'] = t, amp
    # keep only triggers with good hits in >= 2 stations (a track candidate;
    # most of the 3-4 kHz of SPS triggers miss the fan, and single-station
    # hits above the floor are mostly noise); trigger count and span kept
    nst = sum(np.isfinite(out[f'{name}_time']).astype(int)
              for name in STATIONS.values())
    keep = nst >= 2
    out = {k: v[keep] for k, v in out.items()}
    out['n_triggers'] = np.int64(len(ev))
    out['span_s'] = np.float32((ts.max() - ts.min()) / 1e9)
    return out


def extract_one(job):
    run, sub, out = job
    dest = os.path.join(out, run, f'{sub}.npz')
    if os.path.exists(dest):
        return f'{run}/{sub}: exists'
    files = sorted(f for f in glob.glob(os.path.join(EOS_RUNS, run, sub,
                                                     'combined_hits_root', '*.root'))
                   if '.sys.v#' not in f)
    if not files:
        return f'{run}/{sub}: no combined hits'
    try:
        res = leading_hits(files)
    except Exception as e:                          # one bad file must not stop a sweep
        return f'{run}/{sub}: FAILED {e!r}'
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    np.savez_compressed(dest + '.part.npz', **res)
    os.replace(dest + '.part.npz', dest)
    n = {k[:-5]: int(np.isfinite(v).sum()) for k, v in res.items() if k.endswith('_time')}
    return f'{run}/{sub}: {int(res["n_triggers"]):,} triggers, good hits {n}'


def cmd_extract(args):
    from concurrent.futures import ProcessPoolExecutor
    runs = args.run or list(RUNS)
    jobs = []
    for run in runs:
        for d in sorted(glob.glob(os.path.join(EOS_RUNS, run, '*', ''))):
            jobs.append((run, os.path.basename(d.rstrip('/')), args.out))
    print(f'{len(jobs)} sub_runs in {len(runs)} runs -> {args.out}', flush=True)
    with ProcessPoolExecutor(args.jobs) as ex:
        for msg in ex.map(extract_one, jobs):
            print(msg, flush=True)


# --------------------------------------------------------------------------- #
# analyse (local)
# --------------------------------------------------------------------------- #
# Sub_runs to leave out of every curve (see the 2026-07-26 logbook entry).
EXCLUDE = {('mesh_drift_scan_up_1', 'drift2_11_m450_d525'):
           'decode hung: P2_IN and P2_OUT missing',
           ('drift_mesh_scan_1', 'drift_700'):
           'MID-OUT coincidences only in the last ~1/6 of the sub_run'}
COINC_WIN = 150.0          # |dt - median| for the MID-OUT time coincidence [ns]
MIN_COINC = 150            # sub_runs with fewer coincidences are not fitted
N_SLICES = 6               # time slices inside a sub_run
# The simulation session's SPS numbers (HANDOFF_gap_timing_2026-10-15.md, sec. 1.5)
SIM_SPS = {125: 48, 250: 20, 375: 15, 500: 14, 625: 15}


def hw68(x):
    q = np.percentile(x, [16, 84])
    return (q[1] - q[0]) / 2


def boot_err(x, f, n=200, seed=0):
    rng = np.random.default_rng(seed)
    return float(np.std([f(rng.choice(x, len(x))) for _ in range(n)]))


def subrun_table(path):
    """MID-OUT coincidences of one sub_run -> dict of arrays, or None."""
    d = np.load(path)
    tm, to = d['P2_MID_time'], d['P2_OUT_time']
    c = np.isfinite(tm) & np.isfinite(to)
    if c.sum() < 20:
        return None
    dt = tm[c] - to[c]
    w = np.abs(dt - np.median(dt)) < COINC_WIN
    return dict(t_s=d['t_s'][c][w], tm=tm[c][w], to=to[c][w], dt=dt[w],
                am=d['P2_MID_amp'][c][w], ao=d['P2_OUT_amp'][c][w],
                n_trig=int(d['n_triggers']), span=float(d['span_s']))


def summarise(x):
    """Width, medians, amplitudes and yield of one coincidence sample."""
    n = len(x['dt'])
    r = dict(n_coinc=n, yield_per_trig=n / x['n_trig'])
    if n < MIN_COINC:
        return r
    r.update(sigma1=hw68(x['dt']) / np.sqrt(2),
             sigma1_err=boot_err(x['dt'], hw68) / np.sqrt(2),
             med_mid=float(np.median(x['tm'])),
             med_mid_err=boot_err(x['tm'], np.median),
             med_out=float(np.median(x['to'])),
             amp_mid=float(np.median(x['am'])), amp_out=float(np.median(x['ao'])))
    return r


def cmd_analyse(args):
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import sps_config  # noqa: F401  (puts cosmic_bench_analysis on sys.path)
    import p2_provenance as prov

    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                         'font.size': 11, 'figure.facecolor': 'white',
                         'axes.facecolor': 'white', 'savefig.facecolor': 'white',
                         'axes.spines.top': False, 'axes.spines.right': False})
    os.makedirs(args.out, exist_ok=True)
    hv = pd.read_csv(args.hv)
    hv = hv[hv['det'] == 'P2_MID'].set_index(['run', 'sub_run'])
    starts = {}
    if args.run_times and os.path.isfile(args.run_times):
        rt = pd.read_csv(args.run_times)
        starts = {(r.run, r.sub_run): r.start for r in rt.itertuples()}

    rows, slices = [], []
    for run, direction in RUNS.items():
        for f in sorted(glob.glob(os.path.join(args.data, run, '*.npz'))):
            sub = os.path.basename(f)[:-4]
            if (run, sub) in EXCLUDE or (run, sub) not in hv.index:
                continue
            x = subrun_table(f)
            if x is None:
                continue
            drift, mesh = hv.loc[(run, sub), ['drift', 'mesh_or_resist']]
            base = dict(run=run, sub_run=sub, direction=direction,
                        drift_v=drift, mesh_v=mesh,
                        field_v_cm=(drift - mesh) / DRIFT_GAP_CM,
                        n_triggers=x['n_trig'], span_s=x['span'],
                        start_unix=starts.get((run, sub), np.nan))
            rows.append({**base, **summarise(x)})
            edges = np.linspace(0, x['span'], N_SLICES + 1)
            for i in range(N_SLICES):
                m = (x['t_s'] >= edges[i]) & (x['t_s'] < edges[i + 1])
                xs = {k: (v[m] if isinstance(v, np.ndarray) else v) for k, v in x.items()}
                xs['n_trig'] = x['n_trig'] / N_SLICES
                slices.append({**base, 'slice': i, 't_mid_s': (edges[i] + edges[i + 1]) / 2,
                               **summarise(xs)})
    t = pd.DataFrame(rows)
    sl = pd.DataFrame(slices)

    # relative yield: coincidences per trigger over the run's own high-field
    # plateau (mesh 450, E >= 750 V/cm); sqrt -> one chamber, both at the same
    # field. A proxy only: it assumes the muon fraction of the triggers is
    # stable within a run.
    ref = (t[(t.mesh_v == 450) & (t.field_v_cm >= 750)]
           .groupby('run')['yield_per_trig'].median())
    t['rel_yield_1ch'] = np.sqrt(t['yield_per_trig'] / t['run'].map(ref))
    t.to_csv(os.path.join(args.out, 'low_field_timing_subruns.csv'), index=False)
    sl.to_csv(os.path.join(args.out, 'low_field_timing_slices.csv'), index=False)

    lines = []

    def log(s=''):
        print(s)
        lines.append(s)

    # ---- time dependence inside sub_runs: slope over the sub_run ----
    log('time dependence inside each sub_run (mesh 450 V; linear fit over '
        f'{N_SLICES} slices; change first -> last slice):')
    log(f'  {"run":22s} {"sub_run":22s} {"E":>6s} {"n":>7s} '
        f'{"d sigma1 [ns]":>14s} {"d median [ns]":>14s} {"d yield [%]":>12s}')
    trend = []
    for (run, sub), g in sl[sl.mesh_v == 450].groupby(['run', 'sub_run'], sort=False):
        g = g.dropna(subset=['sigma1'])
        if len(g) < 4:
            continue
        span = g.t_mid_s.max() - g.t_mid_s.min()

        def slope(col, err=None):
            w = None if err is None else 1 / np.maximum(g[err], 1e-3) ** 2
            p, cov = np.polyfit(g.t_mid_s, g[col], 1, w=None if w is None else np.sqrt(w), cov='unscaled' if w is not None else True)
            return p[0] * span, np.sqrt(cov[0, 0]) * span
        ds, dse = slope('sigma1', 'sigma1_err')
        dm, dme = slope('med_mid', 'med_mid_err')
        y = g.yield_per_trig / g.yield_per_trig.mean()
        py = np.polyfit(g.t_mid_s, y, 1)[0] * span * 100
        e = t.loc[(t.run == run) & (t.sub_run == sub), 'field_v_cm'].iloc[0]
        n = int(t.loc[(t.run == run) & (t.sub_run == sub), 'n_coinc'].iloc[0])
        trend.append(dict(run=run, sub_run=sub, field_v_cm=e, n_coinc=n,
                          d_sigma1=ds, d_sigma1_err=dse, d_median=dm,
                          d_median_err=dme, d_yield_pct=py))
        log(f'  {run:22s} {sub:22s} {e:6.0f} {n:7d} {ds:+7.1f} ± {dse:4.1f} '
            f'{dm:+7.1f} ± {dme:4.1f} {py:+10.1f}')
    tr = pd.DataFrame(trend)
    tr.to_csv(os.path.join(args.out, 'low_field_timing_trends.csv'), index=False)
    low = tr[tr.field_v_cm <= 400]
    if len(low):
        log(f'  fields <= 400 V/cm, {len(low)} sub_runs: mean change of sigma1 '
            f'{np.average(low.d_sigma1, weights=1/low.d_sigma1_err**2):+.2f} ns, '
            f'of the median {np.average(low.d_median, weights=1/low.d_median_err**2):+.2f} ns, '
            f'chi2/ndf vs no change: sigma1 '
            f'{np.sum((low.d_sigma1/low.d_sigma1_err)**2)/len(low):.2f}, median '
            f'{np.sum((low.d_median/low.d_median_err)**2)/len(low):.2f}')

    # ---- up vs down at mesh 450 ----
    m450 = t[(t.mesh_v == 450) & t.sigma1.notna()]
    up = m450[m450.run == 'drift_mesh_scan_1'].sort_values('field_v_cm')
    dn = m450[m450.run == 'mesh_drift_scan_up_1'].sort_values('field_v_cm')
    log('\nup (drift_mesh_scan_1, 07-25) vs down (mesh_drift_scan_up_1, 07-26), '
        'mesh 450 V; the down scan interpolated to the up-scan fields:')
    log(f'  {"E":>6s} | {"sigma1 up":>11s} {"down":>11s} | {"median MID up":>14s} {"down":>7s} '
        f'| {"rel. yield up":>13s} {"down":>6s} | {"amp MID up":>10s} {"down":>6s}')
    for r in up.itertuples():
        if r.field_v_cm > dn.field_v_cm.max() or r.field_v_cm < dn.field_v_cm.min():
            continue
        ip = {c: np.interp(r.field_v_cm, dn.field_v_cm, dn[c])
              for c in ['sigma1', 'sigma1_err', 'med_mid', 'rel_yield_1ch', 'amp_mid']}
        exact = r.field_v_cm in set(dn.field_v_cm)
        log(f'  {r.field_v_cm:6.0f} | {r.sigma1:5.1f}±{r.sigma1_err:3.1f} '
            f'{ip["sigma1"]:5.1f}±{ip["sigma1_err"]:3.1f} | {r.med_mid:14.1f} {ip["med_mid"]:7.1f} '
            f'| {r.rel_yield_1ch:13.3f} {ip["rel_yield_1ch"]:6.3f} | {r.amp_mid:10.0f} '
            f'{ip["amp_mid"]:6.0f}{"  (same HV)" if exact else ""}')

    # ---- the same mesh 450 V points, day by day ----
    # (the up/down pair above is also a 07-25 vs 07-26 pair; the later up
    # scans tell direction from day)
    day = {'drift_scan_1': '07-24 up', 'drift_scan_2': '07-24 up',
           'drift_mesh_scan_1': '07-25 up', 'mesh_drift_scan_up_1': '07-26 down',
           'drift_mesh_2d_2': '07-27 up', 'p2_mesh_drift_eff_1': '07-28 up'}
    dd = m450[m450.run.isin(day) & ~m450.sub_run.str.startswith('eff_')].copy()
    dd['day'] = dd.run.map(day)
    for q, name in [('sigma1', 'sigma1 [ns]'), ('amp_mid', 'median MID amplitude [ADC]')]:
        piv = dd.pivot_table(index='field_v_cm', columns='day', values=q)
        piv = piv[piv.drop(columns='07-26 down', errors='ignore').notna().any(axis=1)]
        log(f'\n{name} at mesh 450 V, by day (fields of the up scans; the 07-26 down '
            'scan interpolated):')
        if '07-26 down' in piv:
            d6 = dd[dd.day == '07-26 down'].sort_values('field_v_cm')
            piv['07-26 down'] = np.interp(piv.index, d6.field_v_cm, d6[q],
                                          left=np.nan, right=np.nan)
        log('  ' + f'{"E":>6s} ' + ' '.join(f'{c:>11s}' for c in piv.columns))
        for e, r in piv.iterrows():
            log(f'  {e:6.0f} ' + ' '.join(f'{v:11.1f}' for v in r.values))

    # ---- drift_mesh_2d_2: alternating row direction ----
    g2 = t[(t.run == 'drift_mesh_2d_2') & t.sigma1.notna()].copy()
    if len(g2):
        g2['row'] = g2.sub_run.str.slice(4, 6).astype(int)
        g2['row_dir'] = np.where(g2.row % 2 == 0, 'up', 'down')
        log('\ndrift_mesh_2d_2 (rows alternate: even rows up in field, odd rows '
            'down): sigma1 [ns] by mesh row')
        piv = g2.pivot_table(index='field_v_cm', columns='mesh_v', values='sigma1')
        log('  mesh  ' + ' '.join(f'{m:6.0f}{"u" if (i % 2 == 0) else "d"}'
                                  for i, m in enumerate(piv.columns)))
        for e, r in piv.iterrows():
            log(f'  {e:5.0f} ' + ' '.join(f'{v:7.1f}' for v in r.values))
        # each odd (down) row against the mean of its two even (up) neighbours
        log('  down row minus the mean of its two up neighbours, per field:')
        for odd in sorted(g2.row[g2.row % 2 == 1].unique()):
            a = piv.columns[odd - 1], piv.columns[odd], piv.columns[odd + 1]
            dif = piv[a[1]] - (piv[a[0]] + piv[a[2]]) / 2
            log(f'    mesh {a[1]:.0f}: ' + ', '.join(f'{e:.0f}: {v:+.1f}' for e, v in dif.items()))

    with open(os.path.join(args.out, 'low_field_timing_summary.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')

    # ---------------- figures ----------------
    col = {'drift_mesh_scan_1': '#2a78d6', 'mesh_drift_scan_up_1': '#eb6834',
           'drift_scan_1': '#1baf7a', 'drift_scan_2': '#1baf7a'}
    lab = {'drift_mesh_scan_1': 'up: drift_mesh_scan_1 (07-25)',
           'mesh_drift_scan_up_1': 'down: mesh_drift_scan_up_1 (07-26)',
           'drift_scan_1': 'up: drift_scan_1/2 (07-24)', 'drift_scan_2': None}
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    for run in col:
        g = m450[m450.run == run].sort_values('field_v_cm')
        if not len(g):
            continue
        kw = dict(color=col[run], marker='o', ms=5, lw=1.5, label=lab[run])
        axs[0].errorbar(g.field_v_cm, g.sigma1, g.sigma1_err, **kw)
        axs[1].errorbar(g.field_v_cm, g.med_mid - g[g.field_v_cm >= 750].med_mid.median(),
                        g.med_mid_err, **kw)
        if run in ref.index:
            axs[2].plot(g.field_v_cm, g.rel_yield_1ch, **kw)
    se = sorted(SIM_SPS)
    axs[0].plot(se, [SIM_SPS[e] for e in se], color='black', ls='--', lw=1.5,
                label='simulation (handoff 10-15)')
    axs[0].set_ylabel('single-chamber width σ(MID−OUT)/√2 [ns]')
    axs[0].set_title('time resolution')
    axs[1].set_ylabel('median P2_MID time − plateau [ns]')
    axs[1].set_title('median arrival time')
    axs[2].set_ylabel('coincidences per trigger, √(rel. to plateau)')
    axs[2].set_title('relative yield (one chamber)')
    for ax in axs:
        ax.set_xlabel('drift field [V/cm]')
    axs[0].legend(loc='upper center', bbox_to_anchor=(1.7, -0.18), ncol=4, frameon=False)
    fig.suptitle('SPS, Ar/CO2/iso 93/5/2, mesh 450 V: P2_MID–P2_OUT coincidences vs drift field')
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'low_field_timing_vs_field.png'), dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    # time slices inside the low-field sub_runs
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    lowsub = sl[(sl.mesh_v == 450) & (sl.field_v_cm <= 400) & sl.sigma1.notna()
                & sl.run.isin(['drift_mesh_scan_1', 'mesh_drift_scan_up_1'])]
    for (run, sub), g in lowsub.groupby(['run', 'sub_run']):
        kw = dict(color=col[run], alpha=0.6, lw=1, marker='.', ms=4)
        axs[0].plot(g.t_mid_s / 60, g.sigma1 - g.sigma1.mean(), **kw)
        axs[1].plot(g.t_mid_s / 60, g.med_mid - g.med_mid.mean(), **kw)
    for ax, yl in zip(axs, ['σ1 − sub_run mean [ns]', 'median MID time − sub_run mean [ns]']):
        ax.axhline(0, color='#999999', lw=0.8)
        ax.set_xlabel('time since the first trigger of the sub_run [min]')
        ax.set_ylabel(yl)
    axs[0].set_title('width inside each sub_run (E ≤ 400 V/cm)')
    axs[1].set_title('median inside each sub_run (E ≤ 400 V/cm)')
    from matplotlib.lines import Line2D
    axs[0].legend([Line2D([], [], color=col[r]) for r in ['drift_mesh_scan_1', 'mesh_drift_scan_up_1']],
                  ['up scan', 'down scan'], loc='upper center',
                  bbox_to_anchor=(1.1, -0.18), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'low_field_timing_within_subrun.png'), dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    prov.write(args.out, '', __file__,
               inputs={'extracted tables': args.data, 'hv setpoints': args.hv,
                       'runs': ', '.join(f'{r} ({d})' for r, d in RUNS.items())},
               settings={'amp floor [ADC]': AMP_FLOOR, 'coincidence window [ns]': COINC_WIN,
                         'min coincidences': MIN_COINC, 'time slices': N_SLICES,
                         'excluded': EXCLUDE})
    print(f'wrote {args.out}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    sp = ap.add_subparsers(dest='cmd', required=True)
    e = sp.add_parser('extract', help='per-trigger leading hits from EOS (lxplus)')
    e.add_argument('--out', required=True)
    e.add_argument('--run', action='append', help='only these runs (repeatable)')
    e.add_argument('--jobs', type=int, default=4)
    e.set_defaults(func=cmd_extract)
    a = sp.add_parser('analyse', help='timing vs field, direction and time')
    a.add_argument('--data', required=True, help='extract output dir')
    a.add_argument('--out', required=True)
    a.add_argument('--hv', required=True, help='hv_setpoints.csv (eos_inventory)')
    a.add_argument('--run-times', help='optional CSV run,sub_run,start (unix s)')
    a.set_defaults(func=cmd_analyse)
    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    sys.exit(main())
