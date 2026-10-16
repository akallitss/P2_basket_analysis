#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
25_drift_scan_time_slices.py

Does the timing of a drift-scan point change WITHIN the sub_run? (Cosmic part
of test B in HANDOFF_gap_timing_2026-10-15.md.) A dielectric surface charging
up with ions after a drift-HV step would show as a trend over the sub_run,
strongest at low field; a static field or gas effect stays flat.

Per drift-scan sub_run, cut into --slices equal time slices:
  sigma_ns     (p84 - p16)/2 of the per-event leading-pad time_of_max, the
               same estimator as stage 16's time_res_ns (leading pad =
               highest amplitude >= --sig-amp; hot pads dropped; spark veto)
  median_ns    median of that time
  sig_frac     events with a signal pad / triggers in the slice (yield proxy)
  drift_vmon   mean drift-HV monitor voltage in the slice (hv_monitor.csv),
               so an HV still ramping at the start of a sub_run is visible
Then the first -> last slice change of each, with bootstrap errors.

Products (<Analysis>/<det>/<run>/<sub_run>/25_drift_scan_time_slices/):
  time_slices<sfx>.csv     one row per (sub_run, slice)
  time_trends<sfx>.csv     first -> last change per sub_run
  time_slices<sfx>.png     sigma / median / yield / drift vmon vs time
  summary<sfx>.txt, inputs<sfx>.txt

Usage:
  python3 25_drift_scan_time_slices.py [run_key] [--slices 6] [--sig-amp 300]
"""

import argparse
import glob
import importlib
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import p2_qa_config as qa
import p2_io as p2io
import p2_mapping as pmap
import p2_sparks as ps
import p2_provenance as prov

stage16 = importlib.import_module('16_drift_scan_efficiency')

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 11, 'figure.facecolor': 'white',
                     'axes.facecolor': 'white', 'savefig.facecolor': 'white',
                     'axes.spines.top': False, 'axes.spines.right': False})


def sigma68(x):
    q = np.percentile(x, [15.9, 84.1])
    return (q[1] - q[0]) / 2


def boot(x, f, n=300, seed=0):
    rng = np.random.default_rng(seed)
    return float(np.std([f(rng.choice(x, len(x))) for _ in range(n)]))


def leading_times(hits_dir, ct, sig_amp, drop, bad):
    """eventId -> (trigger time [s since first trigger], leading-pad time_of_max);
    plus every trigger's time (denominator of the yield proxy)."""
    ev_t, lead = {}, {}
    for df in p2io.iter_hits(hits_dir, ['eventId', 'trigger_timestamp_ns', 'channel',
                                        'amplitude', 'feu', 'time_of_max'],
                             ct.attrs['feus'], progress=False):
        e = df.drop_duplicates('eventId')
        ev_t.update(zip(e['eventId'], e['trigger_timestamp_ns'].astype(np.int64)))
        h = pmap.attach_pads_to_hits(df, ct)
        h = h[h['mapped'] & h['pad_cx'].notna() & (h['amplitude'] >= sig_amp)]
        if drop:
            h = h[~h['channel_id'].isin(set(drop))]
        if len(h):
            h = h.loc[h.groupby('eventId')['amplitude'].idxmax()]
            for ev, a, t in h[['eventId', 'amplitude', 'time_of_max']].itertuples(index=False):
                if ev not in lead or a > lead[ev][0]:
                    lead[ev] = (float(a), float(t))
    trig = pd.Series(ev_t, dtype='int64')
    if bad:
        trig = trig[~trig.index.isin(bad)]
    t0 = trig.min()
    trig_s = (trig - t0) / 1e9
    ld = pd.DataFrame([(ev, v[1]) for ev, v in lead.items() if ev in trig_s.index],
                      columns=['eventId', 'time_of_max'])
    ld['t_s'] = trig_s.loc[ld['eventId']].to_numpy()
    return trig_s, ld


def drift_vmon(hv_csv, chan, span_s, edges):
    """Mean drift-HV vmon per slice, slice times taken from the HV file's own
    start (the hit clock starts at 0 each sub_run)."""
    if not os.path.isfile(hv_csv):
        return [np.nan] * (len(edges) - 1), np.nan
    hv = pd.read_csv(hv_csv)
    col = f'{chan} vmon'
    if col not in hv:
        return [np.nan] * (len(edges) - 1), np.nan
    t = pd.to_datetime(hv['timestamp'])
    s = (t - t.iloc[0]).dt.total_seconds().to_numpy()
    v = hv[col].to_numpy(float)
    set_v = float(hv[f'{chan} v0'].iloc[-1])
    return [float(np.mean(v[(s >= a) & (s < b)])) if ((s >= a) & (s < b)).any() else np.nan
            for a, b in zip(edges[:-1], edges[1:])], set_v


def main():
    ap = argparse.ArgumentParser(description='Timing inside drift-scan sub_runs.')
    ap.add_argument('run_key', nargs='?', default='det1_driftscan2')
    ap.add_argument('--slices', type=int, default=6)
    ap.add_argument('--sig-amp', type=float, default=300.0)
    ap.add_argument('--drift-chan', default='1:3',
                    help='CAEN board:channel of this detector\'s drift HV')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction, default=True)
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    out = cfg.out_dir('25_drift_scan_time_slices')
    sfx = cfg.product_suffix(args.veto_sparks)
    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy='reverse',
                                  drop_connectors=cfg.DEAD_CONNECTORS,
                                  strategy_overrides=cfg.STRATEGY_OVERRIDES)
    subruns = stage16.find_subruns(cfg, 'drift')

    rows, trends, lines = [], [], []

    def log(s=''):
        print(s)
        lines.append(s)

    for name, mesh, drift in subruns:
        sub = cfg.subrun_dir(name)
        hits_dir = os.path.join(sub, 'combined_hits_root')
        hv_csv = os.path.join(sub, 'hv_monitor.csv')
        drop = p2io.drop_pads_for(cfg, ct, hits_dir=hits_dir)
        bad = None
        if args.veto_sparks and os.path.isfile(hv_csv):
            bad = ps.SparkVeto.from_csv(hv_csv, cfg).vetoed_ids_from_hits(
                hits_dir, ct.attrs['feus'], min_amp=cfg.MIN_AMP)
        trig, ld = leading_times(hits_dir, ct, args.sig_amp, drop, bad)
        span = float(trig.max())
        edges = np.linspace(0, span, args.slices + 1)
        vmon, set_v = drift_vmon(hv_csv, args.drift_chan, span, edges)
        field = (drift - mesh) / qa.DRIFT_GAP_CM if hasattr(qa, 'DRIFT_GAP_CM') \
            else (drift - mesh) / 0.4
        per = []
        for i in range(args.slices):
            m = (ld['t_s'] >= edges[i]) & (ld['t_s'] < edges[i + 1])
            nt = int(((trig >= edges[i]) & (trig < edges[i + 1])).sum())
            x = ld.loc[m, 'time_of_max'].to_numpy()
            r = dict(sub_run=name, mesh_v=mesh, drift_v=drift, field_v_cm=field,
                     slice=i, t_mid_min=(edges[i] + edges[i + 1]) / 120,
                     n_trig=nt, n_sig=len(x), sig_frac=len(x) / max(nt, 1),
                     drift_vmon=vmon[i], drift_set=set_v)
            if len(x) >= 40:
                r.update(sigma_ns=sigma68(x), sigma_err=boot(x, sigma68),
                         median_ns=float(np.median(x)), median_err=boot(x, np.median))
            per.append(r)
        rows += per
        p = pd.DataFrame(per)
        first, last = p.iloc[0], p.iloc[-1]
        tr = dict(sub_run=name, field_v_cm=field, n_sig=int(p.n_sig.sum()),
                  d_sigma=last.get('sigma_ns', np.nan) - first.get('sigma_ns', np.nan),
                  d_sigma_err=np.hypot(last.get('sigma_err', np.nan), first.get('sigma_err', np.nan)),
                  d_median=last.get('median_ns', np.nan) - first.get('median_ns', np.nan),
                  d_median_err=np.hypot(last.get('median_err', np.nan), first.get('median_err', np.nan)),
                  d_sig_frac_pct=100 * (last.sig_frac / first.sig_frac - 1) if first.sig_frac else np.nan,
                  vmon_first=first.drift_vmon, vmon_rest=p.drift_vmon.iloc[1:].mean(),
                  drift_set=set_v)
        trends.append(tr)
        log(f'{name}: E {field:5.0f} V/cm, {tr["n_sig"]:5d} signal events | first -> last '
            f'slice: sigma {tr["d_sigma"]:+5.1f} ± {tr["d_sigma_err"]:4.1f} ns, median '
            f'{tr["d_median"]:+5.1f} ± {tr["d_median_err"]:4.1f} ns, yield '
            f'{tr["d_sig_frac_pct"]:+5.1f} % | drift vmon first slice '
            f'{tr["vmon_first"]:.1f} V, rest {tr["vmon_rest"]:.1f} V (set {set_v:.0f})')

    sl = pd.DataFrame(rows)
    tt = pd.DataFrame(trends)
    sl.to_csv(os.path.join(out, f'time_slices{sfx}.csv'), index=False)
    tt.to_csv(os.path.join(out, f'time_trends{sfx}.csv'), index=False)
    ok = tt.dropna(subset=['d_sigma', 'd_median'])
    low = ok[ok.field_v_cm <= 400]
    if len(low):
        log(f'\nfields <= 400 V/cm ({len(low)} sub_runs): chi2/ndf of the first -> last '
            f'change vs zero: sigma {np.mean((low.d_sigma / low.d_sigma_err) ** 2):.2f}, '
            f'median {np.mean((low.d_median / low.d_median_err) ** 2):.2f}')
    log(f'all sub_runs ({len(ok)}): chi2/ndf sigma '
        f'{np.mean((ok.d_sigma / ok.d_sigma_err) ** 2):.2f}, median '
        f'{np.mean((ok.d_median / ok.d_median_err) ** 2):.2f}')
    with open(os.path.join(out, f'summary{sfx}.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')

    fig, axs = plt.subplots(1, 4, figsize=(18, 4.2))
    cmap = plt.get_cmap('viridis')
    fields = sorted(sl.field_v_cm.unique())
    for j, e in enumerate(fields):
        g = sl[sl.field_v_cm == e]
        kw = dict(color=cmap(j / max(len(fields) - 1, 1)), marker='o', ms=3, lw=1.2,
                  label=f'{e:.0f}')
        axs[0].plot(g.t_mid_min, g.sigma_ns - g.sigma_ns.mean(), **kw)
        axs[1].plot(g.t_mid_min, g.median_ns - g.median_ns.mean(), **kw)
        axs[2].plot(g.t_mid_min, g.sig_frac / g.sig_frac.mean(), **kw)
        axs[3].plot(g.t_mid_min, g.drift_vmon - g.drift_set, **kw)
    for ax, (yl, ti) in zip(axs, [('σ − sub_run mean [ns]', 'time resolution'),
                                  ('median − sub_run mean [ns]', 'median peak time'),
                                  ('signal events per trigger / mean', 'yield proxy'),
                                  ('drift vmon − setpoint [V]', 'drift HV')]):
        ax.set_xlabel('time in the sub_run [min]')
        ax.set_ylabel(yl)
        ax.set_title(ti)
    axs[3].legend(title='E [V/cm]', loc='center left', bbox_to_anchor=(1.02, 0.5),
                  frameon=False, fontsize=8)
    fig.suptitle(f'{cfg.DET_TAG} cosmic drift scan (mesh {subruns[0][1]} V): '
                 'inside each 30 min sub_run')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'time_slices{sfx}.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)

    prov.write(out, sfx, __file__,
               inputs={'run dir': cfg.run_dir,
                       'sub_runs': ', '.join(n for n, _, _ in subruns)},
               settings={'slices': args.slices, 'sig_amp [ADC]': args.sig_amp,
                         'drift HV channel': args.drift_chan,
                         'spark veto': args.veto_sparks})
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
