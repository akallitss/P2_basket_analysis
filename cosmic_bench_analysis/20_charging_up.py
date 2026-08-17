#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
20_charging_up.py

Charging-up: does the efficiency quoted for a long run describe the detector,
or the average of a detector that was still charging?

An MPGD's insulating surfaces charge under irradiation and the gain creeps up
over hours. A run-integrated efficiency then UNDERSTATES the plateau value, and
by an amount nobody can guess from the integrated number alone. `det4_long2`
went 31 % -> 92 % over 10 h; quoting its 88.3 % run average as "the detector"
is wrong by 4 points.

This stage separates the two by construction:

  * efficiency, rate, arrival time and amplitude vs time, from the per-ray list
    written by stage 06 (so the efficiency definition is identical to it)
  * a plateau test: split the run into an early window and a late window and
    compare. If the late window agrees with the last-half window, the detector
    has plateaued and the late number is the detector's efficiency.

The distinction matters beyond bookkeeping: a rising rate is charging-up, while
an ON/OFF rate with a matching ARRIVAL-TIME shift is an intermittent drift
field (a broken HV contact to the drift foil lets charge amplify but not
drift). The timeline plot shows both, so they can be told apart by eye.

Products (<Analysis>/<detN>/<run>/<sub_run>/20_charging_up/):
  charging_timeline<sfx>.png   efficiency / rate / arrival time / amplitude vs t
  charging_summary<sfx>.txt    early vs late windows, plateau verdict
  charging_vs_time<sfx>.csv    the binned series

Usage:
  python3 20_charging_up.py [run_key] [--bin-min 20] [--late-h 2] [--early-h 1]
"""

import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import p2_qa_config as qa
qa.setup_paths()
import p2_io as p2io
import p2_mapping as pmap
import p2_align as pa
import p2_sparks as ps


def leading_pad_times(cfg, ct, drop, bad, sig_amp=300.0):
    """Per-event (leading-pad time_of_max, amplitude) -- same estimator as the
    quoted time resolution, so the timeline is comparable with stage 16/19."""
    best_t, best_a = {}, {}
    for df in p2io.iter_hits(cfg.combined_hits_dir,
                             ['eventId', 'channel', 'feu', 'amplitude',
                              'time_of_max'],
                             ct.attrs['feus'], progress=False,
                             min_amp=cfg.MIN_AMP,
                             t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H):
        h = pmap.attach_pads_to_hits(df, ct)
        h = h[h['mapped'] & h['pad_cx'].notna() & (h['amplitude'] >= sig_amp)]
        if drop:
            h = h[~h['channel_id'].isin(drop)]
        if bad:
            h = h[~h['eventId'].isin(bad)]
        if not len(h):
            continue
        sel = h.loc[h.groupby('eventId')['amplitude'].idxmax(),
                    ['eventId', 'amplitude', 'time_of_max']]
        for ev, a, t in sel.itertuples(index=False):
            if ev not in best_a or a > best_a[ev]:
                best_a[ev] = float(a)
                best_t[ev] = float(t)
    return pd.DataFrame({'eventId': list(best_t),
                         't_pad': list(best_t.values()),
                         'amp': [best_a[e] for e in best_t]})


def window(d, lo_h, hi_h):
    m = d[(d['h'] >= lo_h) & (d['h'] < hi_h)]
    if not len(m):
        return dict(n=0, eff=np.nan, any=np.nan)
    return dict(n=len(m), eff=100 * m['within'].mean(),
                any=100 * m['has_any'].mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--bin-min', type=float, default=20.0, help='bin [minutes]')
    ap.add_argument('--early-h', type=float, default=1.0,
                    help='length of the early window [h]')
    ap.add_argument('--late-h', type=float, default=2.0,
                    help='length of the late (plateau) window [h]')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True)
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    out = cfg.out_dir('20_charging_up')
    sfx = cfg.product_suffix(args.veto_sparks)

    rays_csv = os.path.join(cfg.OUT_BASE, '06_efficiency',
                            f'ray_hit_miss_list{sfx}.csv')
    if not os.path.isfile(rays_csv):
        print(f'!! need stage 06 first: {rays_csv}')
        return
    rays = pd.read_csv(rays_csv)
    rays = rays[rays['in_active']] if 'in_active' in rays.columns else rays
    if 't_sec' not in rays.columns:
        print('!! ray list has no t_sec column -- rerun stage 06')
        return
    rays = rays.dropna(subset=['t_sec'])
    rays['h'] = rays['t_sec'] / 3600.0
    dur = float(rays['h'].max())
    print(f'  {len(rays):,} active-area rays over {dur:.2f} h')

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    drop = set(p2io.drop_pads_for(cfg, ct))
    bad = set()
    if args.veto_sparks:
        sv = ps.SparkVeto.from_cfg(cfg)
        bad = set(int(b) for b in sv.vetoed_ids_from_hits(
            cfg.combined_hits_dir, ct.attrs['feus'], min_amp=cfg.MIN_AMP))
    tt = leading_pad_times(cfg, ct, drop, bad)
    rays = rays.merge(tt, on='eventId', how='left')

    b = args.bin_min / 60.0
    rays['bin'] = (rays['h'] // b).astype(int)
    g = rays.groupby('bin').agg(
        h=('h', 'mean'), n=('within', 'size'),
        eff=('within', lambda v: 100 * v.mean()),
        any=('has_any', lambda v: 100 * v.mean()),
        t_pad=('t_pad', 'median'), amp=('amp', 'median')).reset_index()
    g = g[g['n'] >= 20]
    g['eff_err'] = np.sqrt(g['eff'] * (100 - g['eff']) / g['n'])
    g.to_csv(os.path.join(out, f'charging_vs_time{sfx}.csv'), index=False)

    # The windows must not overlap or swallow the run: on a 2 h run a 2 h
    # "late window" IS the whole run, and comparing it with the run tells you
    # nothing. Clamp both to a third of the duration each.
    early_h = min(args.early_h, dur / 3.0)
    late_h = min(args.late_h, dur / 3.0)
    short = (early_h < args.early_h) or (late_h < args.late_h)
    early = window(rays, 0.0, early_h)
    late = window(rays, max(0.0, dur - late_h), dur + 1)
    half = window(rays, dur / 2.0, dur + 1)
    full = window(rays, -1, dur + 1)
    gain = late['eff'] - early['eff']

    # Classify by SIGN as well as magnitude. Charging-up rises; a degrading
    # detector falls; calling a fall "still rising" (as an earlier version did
    # for det3, which drops 89.6 -> 82.5 % inside 2 h) inverts the physics.
    err = np.hypot(np.sqrt(max(early['eff'] * (100 - early['eff']), 1) / max(early['n'], 1)),
                   np.sqrt(max(late['eff'] * (100 - late['eff']), 1) / max(late['n'], 1)))
    if gain < -3 * err and gain < -2.0:
        verdict = ('DEGRADING — efficiency FALLS during the run; this is not '
                   'charging-up and no window is "the" efficiency')
    elif abs(late['eff'] - half['eff']) < 2.0 and gain > -2.0:
        verdict = 'PLATEAUED — quote the late-window value'
    else:
        verdict = 'STILL RISING — the late window is a lower limit'
    if short:
        verdict += (f'  [run only {dur:.1f} h: windows clamped to '
                    f'{early_h:.2f}/{late_h:.2f} h — weak test]')

    lines = [
        f'CHARGING-UP — {cfg.DET_TAG} {cfg.RUN}/{cfg.SUB_RUN}',
        f'  run duration            : {dur:.2f} h   ({len(rays):,} active rays)',
        f'  early  (0–{early_h:.2f} h)     : eff {early["eff"]:.1f} %   '
        f'any {early["any"]:.1f} %   n={early["n"]:,}',
        f'  second half (>{dur/2:.1f} h)   : eff {half["eff"]:.1f} %   '
        f'any {half["any"]:.1f} %   n={half["n"]:,}',
        f'  late   (last {late_h:.2f} h)   : eff {late["eff"]:.1f} %   '
        f'any {late["any"]:.1f} %   n={late["n"]:,}',
        f'  FULL RUN (as quoted)    : eff {full["eff"]:.1f} %   '
        f'any {full["any"]:.1f} %',
        '',
        f'  charging-up gain        : {gain:+.1f} points (early -> late)',
        f'  run average understates : {late["eff"] - full["eff"]:+.1f} points',
        f'  VERDICT                 : {verdict}',
    ]
    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(out, f'charging_summary{sfx}.txt'), 'w') as fh:
        fh.write(txt + '\n')

    fig, axs = plt.subplots(4, 1, figsize=(12, 11), sharex=True)
    axs[0].errorbar(g['h'], g['eff'], yerr=g['eff_err'], fmt='o-', ms=3,
                    color='steelblue', label='eff (reco within R)')
    axs[0].plot(g['h'], g['any'], 's--', ms=3, color='darkorange',
                alpha=.8, label='has_any')
    axs[0].axhline(full['eff'], color='grey', ls=':',
                   label=f'run average {full["eff"]:.1f} %')
    axs[0].axhline(late['eff'], color='crimson', ls='--',
                   label=f'late window {late["eff"]:.1f} %')
    axs[0].set_ylabel('efficiency [%]'); axs[0].legend(fontsize=7, ncol=2)
    axs[1].plot(g['h'], g['n'], 'o-', ms=3, color='seagreen')
    axs[1].set_ylabel(f'rays / {args.bin_min:g} min')
    axs[2].plot(g['h'], g['t_pad'], 'o-', ms=3, color='purple')
    axs[2].set_ylabel('median peak time [ns]')
    axs[3].plot(g['h'], g['amp'], 'o-', ms=3, color='darkorange')
    axs[3].set_ylabel('median amplitude [ADC]')
    axs[3].set_xlabel('time since run start [h]')
    for a in axs:
        a.grid(alpha=.3)
        a.axvspan(0, early_h, color='crimson', alpha=.07)
        a.axvspan(max(0, dur - late_h), dur, color='seagreen', alpha=.07)
    fig.suptitle(f'{cfg.DET_NAME} charging-up — {cfg.RUN}/{cfg.SUB_RUN}\n'
                 f'early {early["eff"]:.1f} % -> late {late["eff"]:.1f} % '
                 f'({gain:+.1f} pts); run average {full["eff"]:.1f} %  —  {verdict}',
                 fontsize=10)
    fig.tight_layout()
    f = os.path.join(out, f'charging_timeline{sfx}.png')
    fig.savefig(f, dpi=150, bbox_inches='tight'); plt.close(fig)
    print('  saved', f)


if __name__ == '__main__':
    main()
