#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pad efficiency of a global-threshold scan against each channel's (threshold - baseline) in mV, channels grouped.

The wk2 adaptive scan (3-5 Nov 2026) moved ONE threshold per chip (sdt = chip base + K). The distance between a
channel's threshold and its own baseline is therefore different for every channel of a chip. Here:

  * threshold and baseline per channel in mV: the monitor map of the K=-18 file (vmm_monitor_tool, 13 Nov, sd=0,
    measured unmasked); for any other run thr = thr_map + (sdt_run - sdt_map) * mV_per_sdt(chip) (sdtcal.json);
    sdt_run read from the run's own ext file (run_sdt.json: {cfg: {'sdt': {chip: sdt}}}).
  * efficiency per pad: the per-pad on/off tag counts of efficiency_tp.py (eff_runN.json), off-time subtracted,
    runs at the same K summed. A pad counts as hit when ANY live pad within probe_r (15 mm) of the prediction fired,
    so a pad's efficiency also contains its neighbours (blind pads, baseline above threshold, still read 0.06-0.23).
  * relative efficiency: k / (n * channel mean over the whole scan) - removes the position/chamber dependence and
    keeps only how the channel changes with its threshold.

Groupings: (a) all channels in bins of the distance, per chamber; (b) channels grouped by their baseline relative to
the chip median (low / middle / high) against K; (c) channels grouped by absolute baseline (mixes chips/chambers).

Usage: eff_vs_threshold_distance.py ANALYSIS_DIR MAP_DIR OUT_DIR
       (MAP_DIR holds map_Km18_i5_unmasked.json, sdtcal.json, run_sdt.json)
"""
import glob
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from event_display import style, INK, MUTED, COL, MRK  # noqa: E402

PROBES = (('probe_det3', 3, 'P2_IN'), ('probe_det2', 2, 'P2_MID'), ('probe_det1', 1, 'P2_OUT'))
DBINS = np.r_[0, 14, 18, 22, 26, 30, 34, 40, 46, 52, 58, 64, 70, 100]


def load(adir, mdir):
    mp = json.load(open(f'{mdir}/map_Km18_i5_unmasked.json'))
    cal = json.load(open(f'{mdir}/sdtcal.json'))['fit']
    rsdt = json.load(open(f'{mdir}/run_sdt.json'))
    rows = []
    for f in sorted(glob.glob(f'{adir}/eff_run*.json')):
        E = json.load(open(f))
        run, tag = E['title'].split(' ', 1)
        cfg = f'p2b-config-cern-ext_gain4.5_peaktime200_opt_cosmic_{tag}.txt'
        K = int(tag.split('_adwk2_K')[1].split('_')[0].replace('p', '+').replace('m', '-'))
        sdt = rsdt[cfg]['sdt']
        for key, det, _ in PROBES:
            r = E[key]['r_scale']
            for v, c, x, y, n_on, n_off, k_on, k_off in E[key]['per_pad']:
                rows.append((run, K, det, v, c, x, y, n_on - r * n_off, k_on - r * k_off, sdt[str(v)]))
    t = pd.DataFrame(rows, columns=['run', 'K', 'det', 'vmm', 'ch', 'x', 'y', 'n', 'k', 'sdt'])
    t['baseline'] = [mp['baseline'][str(v)][str(c)] for v, c in zip(t.vmm, t.ch)]
    thr = np.array([mp['threshold'][str(v)][str(c)] for v, c in zip(t.vmm, t.ch)], float)
    t['thr'] = thr + (t.sdt - t.vmm.map(lambda v: mp['sdt'][str(v)])) * t.vmm.map(lambda v: cal[str(v)]['mV_per_sdt'])
    t['dist'] = t.thr - t.baseline
    return t


def rel(d):
    e = d.exp.sum()
    return pd.Series({'rel': d.k.sum() / e, 'err': np.sqrt(e) / e, 'n': d.n.sum(), 'eff': d.k.sum() / d.n.sum(),
                      'dist': np.average(d.dist, weights=d.n.clip(1e-9)),
                      'channels': d[['vmm', 'ch']].drop_duplicates().shape[0]})


def main():
    adir, mdir, odir = sys.argv[1:4]
    os.makedirs(odir, exist_ok=True)
    t = load(adir, mdir)
    t.to_csv(f'{odir}/rows.csv', index=False)
    g = t.groupby(['K', 'det', 'vmm', 'ch'], as_index=False).agg(n=('n', 'sum'), k=('k', 'sum'), dist=('dist', 'mean'),
                                                                baseline=('baseline', 'first'))
    blind = g[g.dist < 0]
    g = g[g.dist >= 0].copy()
    c = g.groupby(['vmm', 'ch']).agg(N=('n', 'sum'), Ksum=('k', 'sum'))
    c['m'] = c.Ksum / c.N
    g = g.join(c[['m', 'N']], on=['vmm', 'ch'])
    g = g[g.N >= 100].copy()
    g['exp'] = g.n * g.m
    g['dbin'] = pd.cut(g.dist, DBINS)
    g['rb'] = g.baseline - g.groupby('vmm').baseline.transform('median')
    g['rbg'] = pd.cut(g.rb, [-99, -2, 2, 99], labels=['low (< -2 mV)', 'middle', 'high (> +2 mV)'])
    g['abg'] = pd.cut(g.baseline, [0, 163, 166, 168, 170, 172, 175, 1000], right=False)

    a = g.groupby(['det', 'dbin'], observed=True).apply(rel, include_groups=False).reset_index()
    a_all = g.groupby('dbin', observed=True).apply(rel, include_groups=False).reset_index()
    b = g.groupby(['rbg', 'K'], observed=True).apply(rel, include_groups=False).reset_index()
    cabs = g.groupby(['abg', 'K'], observed=True).apply(rel, include_groups=False).reset_index()
    perk = g.groupby(['det', 'K']).apply(rel, include_groups=False).reset_index()
    bl = blind.groupby('vmm').apply(lambda d: pd.Series({'n': d.n.sum(), 'eff': d.k.sum() / d.n.sum(),
                                                         'channels': d.ch.nunique()}),
                                    include_groups=False)
    for nm, df in (('by_distance_det', a), ('by_distance_all', a_all), ('by_rel_baseline_K', b),
                   ('by_abs_baseline_K', cabs), ('per_K_det', perk)):
        df.astype({col: str for col in df.columns if df[col].dtype.name == 'category'}).to_csv(f'{odir}/{nm}.csv', index=False)
    bl.to_csv(f'{odir}/blind_pads.csv')

    pd.set_option('display.width', 200)
    print(f'{g[["vmm", "ch"]].drop_duplicates().shape[0]} channels (>= 100 tags over the scan), '
          f'{g.n.sum():.0f} tags; tags per channel and K: median {g.n.median():.1f}')
    print('\n(a) relative efficiency vs threshold - baseline, all chambers:')
    print(a_all[['dbin', 'rel', 'err', 'n', 'channels']].round(3).to_string(index=False))
    print('\n(b) relative efficiency vs K, channels grouped by baseline relative to the chip median:')
    print(b.pivot(index='K', columns='rbg', values='rel').round(3))
    print(b.pivot(index='K', columns='rbg', values='dist').round(1))
    print('\npads whose baseline is above the threshold (blind blocks), absolute efficiency:')
    print(bl.round(3))

    # figure
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for _, det, name in PROBES:
        s = a[a.det == det]
        ax[0].errorbar(s.dist, s.rel, s.err, color=COL[det], marker=MRK[det], ms=4, lw=1.5, capsize=0, label=name)
    ax[0].axhline(1, color=MUTED, lw=0.8, ls=':')
    ax[0].set_xlabel('threshold - baseline of the channel (mV)', color=INK, fontsize=8)
    ax[0].set_ylabel('efficiency / channel mean over the scan', color=INK, fontsize=8)
    ax[0].set_title('all channels, binned by their distance', color=INK, fontsize=9, loc='left')
    ax[0].set_ylim(0.55, 1.2)
    ax[0].legend(frameon=False, fontsize=7)
    blues = {'low (< -2 mV)': '#9ec5f0', 'middle': '#4f94e0', 'high (> +2 mV)': '#1a4f94'}
    mk = {'low (< -2 mV)': 'v', 'middle': 'o', 'high (> +2 mV)': '^'}
    for lab in blues:
        s = b[b.rbg == lab].sort_values('K')
        ax[1].errorbar(s.K, s.rel, s.err, color=blues[lab], marker=mk[lab], ms=4, lw=1.5, capsize=0,
                       label=f'baseline {lab} vs chip median')
    ax[1].axhline(1, color=MUTED, lw=0.8, ls=':')
    ax[1].set_xlabel('K (sdt counts from the chip base)', color=INK, fontsize=8)
    ax[1].set_title('channels grouped by baseline within their chip', color=INK, fontsize=9, loc='left')
    ax[1].set_ylim(0.55, 1.2)
    ax[1].legend(frameon=False, fontsize=7)
    for x in ax:
        style(x)
    fig.suptitle('wk2 global-threshold scan (runs 181-225): pad efficiency vs threshold - baseline', color=INK,
                 fontsize=10, x=0.01, ha='left')
    fig.tight_layout()
    fig.savefig(f'{odir}/eff_vs_threshold_distance.png', dpi=150)
    print(f'\nwrote {odir}/')


if __name__ == '__main__':
    main()
