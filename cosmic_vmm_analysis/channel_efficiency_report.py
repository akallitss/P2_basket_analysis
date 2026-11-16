#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verification figures of a channel_efficiency.py measurement (for CHANNEL_EFFICIENCY.md and the CERN-site note).

Figures (PNG, light background) in FIG_DIR, numbers in FIG_DIR/numbers.json:
  v1_scint_timing     tag time minus scintillator time: the muon peak and the --scint window
  v2_prediction_grid  where the predicted points land relative to the nearest pad centre (the lattice steps)
  v3_counting_curves  the counting curve per chamber, the eps + b R^2 line, and the R_lo stability scan
  v4_eff_vs_distance  absolute efficiency against threshold - baseline, per chamber
  v5_runs             one point per run: the three repeated runs at K = 0, and K across the scan
  v6_methods          the same pads by four methods: naive own pad, all tags, scintillator-confirmed counting
  v7_toy_closure      simulated muons: true vs measured efficiency for several resolutions and effects

Usage: channel_efficiency_report.py MEASURE_DIR TAGS_DIR FIG_DIR [--alltags MEASURE_DIR_NO_SCINT] [--scint 140:310]
"""
import argparse
import glob
import json
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import channel_efficiency as CE  # noqa: E402
from event_display import style, INK, MUTED, COL, MRK  # noqa: E402

DETS = (('P2_IN', 3), ('P2_MID', 2), ('P2_OUT', 1))


def save(fig, d, name):
    fig.tight_layout()
    fig.savefig(f'{d}/{name}.png', dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('measure_dir')
    ap.add_argument('tags_dir')
    ap.add_argument('fig_dir')
    ap.add_argument('--alltags', default='')
    ap.add_argument('--scint', default='140:310')
    ap.add_argument('--toy-muons', type=int, default=20000)
    a = ap.parse_args()
    os.makedirs(a.fig_dir, exist_ok=True)
    lo, hi = map(float, a.scint.split(':'))
    R = pd.read_csv(f'{a.measure_dir}/efficiency.csv')
    Z = np.load(f'{a.measure_dir}/curves.npz', allow_pickle=True)
    num = {}

    # ---- tags: timing, prediction offsets, naive own-pad efficiency
    dts, offs, naive = {d: [] for d, _ in DETS}, {d: [] for d, _ in DETS}, {d: [0.0, 0.0, 0.0, 0.0] for d, _ in DETS}
    for fz in sorted(glob.glob(f'{a.tags_dir}/tags_run*.npz')):
        for det, p in DETS:
            t = CE.load_tags(fz, p)
            dts[det].append(t['dt_scint'][t['on']])
            sc = t['on'] & (t['dt_scint'] >= lo) & (t['dt_scint'] <= hi)
            offs[det].append(np.c_[t['px'][sc] - t['cx'][sc], t['py'][sc] - t['cy'][sc]])
            w = np.where(t['on'], 1.0, -float(t['r_scale']))
            m = (t['dt_scint'] >= lo) & (t['dt_scint'] <= hi)
            naive[det][0] += (w * m).sum(); naive[det][1] += (w * m * t['own']).sum()
            naive[det][2] += w.sum(); naive[det][3] += (w * t['own']).sum()

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    e = np.arange(-600, 801, 25)
    for det, p in DETS:
        d = np.concatenate(dts[det])
        ax.hist(d, e, histtype='step', color=COL[p], lw=1.6, label=det)
    ax.axvspan(lo, hi, color=MUTED, alpha=0.15, lw=0, label=f'window {lo:g}-{hi:g} ns')
    ax.set_yscale('log')
    ax.set_xlabel('tag time - nearest scintillator coincidence (ns)', fontsize=8, color=INK)
    ax.set_ylabel('on-time tags / 25 ns', fontsize=8, color=INK)
    ax.legend(frameon=False, fontsize=7)
    style(ax)
    save(fig, a.fig_dir, 'v1_scint_timing')
    alld = np.concatenate([np.concatenate(dts[d]) for d, _ in DETS])
    num['scint'] = dict(window=[lo, hi], in_window=int(((alld >= lo) & (alld <= hi)).sum()),
                        on_tags=int(len(alld)), outside_peak=int(((alld > hi + 100) | (alld < lo - 100)).sum()))

    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.4))
    for k, (det, p) in enumerate(DETS):
        o = np.concatenate(offs[det])
        r = np.hypot(o[:, 0], o[:, 1])
        ax[k].hist(r, np.arange(0, 9.01, 0.25), color=COL[p], alpha=0.85)
        ax[k].set_title(f'{det}: {np.mean(r < 0.5) * 100:.0f} % exactly on a pad centre', fontsize=8.5, color=INK,
                        loc='left')
        ax[k].set_xlabel('predicted point - nearest pad centre (mm)', fontsize=8, color=INK)
        style(ax[k])
        num.setdefault('on_centre_frac', {})[det] = float(np.mean(r < 0.5))
    ax[0].set_ylabel('scintillator-confirmed tags', fontsize=8, color=INK)
    save(fig, a.fig_dir, 'v2_prediction_grid')

    # ---- counting curves and R_lo scan
    pl = R[R.kind == 'plateau'].set_index('det')
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8), gridspec_kw=dict(width_ratios=[1.6, 1]))
    radii = Z['radii']
    rlo_cols = [c for c in R.columns if c.startswith('eps_rlo')]
    for det, p in DETS:
        m = (Z['det'] == det) & Z['clean'] & (Z['dist_mV'] >= 18) & (Z['dist_mV'] < 40)
        g = CE.fit_group(Z['n'][m], Z['k'][m], Z['s'][m], Z['v'][m], Z['area'][m], radii)
        ax[0].errorbar(radii, g['K'] / g['E'], g['Kerr'] / g['E'], marker=MRK[p], color=COL[p], ms=4, lw=1.2,
                       label=f'{det}')
        rr = np.linspace(0, radii[-1], 50)
        ax[0].plot(rr, g['eps'] + g['b'] * rr ** 2, color=COL[p], lw=0.9, ls='--')
        xs = [int(c[7:]) for c in rlo_cols]
        ax[1].plot(xs, [pl.loc[det, c] for c in rlo_cols], marker=MRK[p], color=COL[p], lw=1.5, label=det)
        num.setdefault('curve', {})[det] = dict(R=radii.tolist(), y=(g['K'] / g['E']).round(4).tolist(),
                                                ye=(g['Kerr'] / g['E']).round(4).tolist(), eps=g['eps'], b=g['b'])
    ax[0].axvline(40, color=MUTED, lw=0.8, ls=':')
    ax[0].set_xlabel('R (mm): tags counted within R of the pad centre', fontsize=8, color=INK)
    ax[0].set_ylabel('pad firings (random removed) / expected crossings', fontsize=8, color=INK)
    ax[0].set_title('counting curve, pads 18-40 mV above baseline', fontsize=9, color=INK, loc='left')
    ax[1].set_xlabel('R_lo (mm)', fontsize=8, color=INK)
    ax[1].set_ylabel('efficiency', fontsize=8, color=INK)
    ax[1].set_title('stability: result vs R_lo', fontsize=9, color=INK, loc='left')
    for x in ax:
        x.legend(frameon=False, fontsize=7)
        style(x)
    save(fig, a.fig_dir, 'v3_counting_curves')
    num['plateau'] = {det: dict(eps=float(pl.loc[det, 'eps']), err=float(pl.loc[det, 'err']),
                                chi2ndf=float(pl.loc[det, 'chi2ndf']), pads=int(pl.loc[det, 'pads']),
                                **{c: float(pl.loc[det, c]) for c in rlo_cols}) for det, _ in DETS}
    bl = R[R.kind == 'blind'].iloc[0]
    num['blind'] = dict(eps=float(bl.eps), err=float(bl.err), pads=int(bl.pads))

    # ---- efficiency vs distance
    dd = R[R.kind == 'distance']
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    for det, p in DETS:
        d = dd[(dd.det == det) & dd.eps.notna() & (dd.err < 0.2)]
        ax.errorbar(d.dist_mV, d.eps, d.err, marker=MRK[p], color=COL[p], ms=4, lw=1.5, label=det)
        num.setdefault('vs_distance', {})[det] = d[['bin', 'dist_mV', 'eps', 'err', 'pads']].round(4).to_dict('records')
    ax.set_ylim(0, 1.1)
    ax.set_xlabel('threshold - baseline of the channel (mV)', fontsize=8, color=INK)
    ax.set_ylabel('single-channel efficiency', fontsize=8, color=INK)
    ax.legend(frameon=False, fontsize=7)
    style(ax)
    save(fig, a.fig_dir, 'v4_eff_vs_distance')

    # ---- per run
    rr = R[R.kind == 'run']
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    for k, (det, p) in enumerate(DETS):
        d = rr[(rr.det == det) & (rr.err < 0.15)]
        ax.errorbar(d.K + (k - 1) * 1.2, d.eps, d.err, marker=MRK[p], color=COL[p], ms=4, lw=0, elinewidth=1, label=det)
    ax.set_xlabel('K (sdt counts from the chip base; repeated runs at K = 0 drawn side by side)', fontsize=8, color=INK)
    ax.set_ylabel('efficiency, all live pads of the run', fontsize=8, color=INK)
    ax.set_ylim(0, 1.2)
    ax.legend(frameon=False, fontsize=7)
    style(ax)
    save(fig, a.fig_dir, 'v5_runs')
    k0 = rr[rr.K == 0]
    num['k0_repeats'] = k0[['run', 'det', 'eps', 'err']].round(4).to_dict('records')

    # ---- methods
    meth = {det: {'naive own pad, all tags': naive[det][3] / naive[det][2],
                  'naive own pad, scint tags': naive[det][1] / naive[det][0]} for det, _ in DETS}
    if a.alltags:
        A = pd.read_csv(f'{a.alltags}/efficiency.csv')
        for det, _ in DETS:
            meth[det]['counting, all tags'] = float(A[(A.kind == 'plateau') & (A.det == det)].eps.iloc[0])
    for det, _ in DETS:
        meth[det]['counting, scint tags'] = float(pl.loc[det, 'eps'])
    num['methods'] = meth
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    labels = list(meth['P2_IN'].keys())
    shades = ['#c9c8c1', '#a6a59d', '#6b6a64', '#1f1f1e']
    for i, lab in enumerate(labels):
        for k, (det, p) in enumerate(DETS):
            ax.bar(k + (i - 1.5) * 0.2, meth[det][lab], 0.18, color=shades[i] if i < 3 else COL[p],
                   label=lab if k == 0 else None)
    ax.set_xticks(range(3), [d for d, _ in DETS])
    ax.set_ylabel('efficiency (plateau pads)', fontsize=8, color=INK)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    style(ax)
    save(fig, a.fig_dir, 'v6_methods')

    # ---- toy closure
    cases = [('Gaussian', {}), ('on the pad lattice', dict(quantize=True)), ('random hits', dict(noise=0.0005)),
             ('15 % tail at 20 mm', dict(tail=(0.15, 20.0))), ('30 % fake tags', dict(impure=0.3))]
    tt = []
    for nm, kw in cases:
        for sg in (3.0, 8.0, 12.0):
            r = CE.toy(sg, 0.85, a.toy_muons, seed=11, **kw)
            tt.append(dict(case=nm, sigma=sg, est=r['eps_est'], err=r['err'], naive=r['naive'],
                           rlo56=r.get('eps_rlo56')))
    T = pd.DataFrame(tt)
    num['toy'] = T.round(4).to_dict('records')
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    for i, (nm, _) in enumerate(cases):
        d = T[T.case == nm]
        ax.errorbar(i + (d.sigma.to_numpy() - 8) / 14, d.est, d.err, fmt='o', color=INK, ms=4)
        ax.plot(i + (d.sigma.to_numpy() - 8) / 14, d.naive, 'x', color=MUTED, ms=6)
    ax.axhline(0.85, color=COL[2], lw=1, ls='--', label='true 0.85')
    ax.plot([], [], 'o', color=INK, label='counting (sigma 3 / 8 / 12 mm, left to right)')
    ax.plot([], [], 'x', color=MUTED, label='naive "did the predicted pad fire"')
    ax.set_xticks(range(len(cases)), [c for c, _ in cases], fontsize=7.5)
    ax.set_ylim(0, 1)
    ax.set_ylabel('measured efficiency', fontsize=8, color=INK)
    ax.legend(frameon=False, fontsize=7, loc='lower left')
    style(ax)
    save(fig, a.fig_dir, 'v7_toy_closure')
    json.dump(num, open(f'{a.fig_dir}/numbers.json', 'w'), indent=1, default=float)
    print(json.dumps({k: num[k] for k in ('scint', 'on_centre_frac', 'plateau', 'blind', 'methods')}, indent=1,
                     default=float))


if __name__ == '__main__':
    main()
