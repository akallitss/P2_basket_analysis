#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Does a channel's ADC baseline predict its efficiency? (the per-channel threshold spread of untrimmed chips)

A chip compares all its channels with ONE threshold voltage, so a channel sitting on a higher baseline needs less
signal to cross it: with no trims, efficiency should rise with the channel's baseline relative to its chip.
Efficiency per channel: the per-pad counts of efficiency_tp.py outputs (on/off-time tags, with/without probe hit),
summed over the runs given (same map and thresholds only), off-time subtracted. Baselines: channel_baselines csv of
threshold_table.py (HV-on value, else HV-off). Per chip, the channel baseline minus the chip median; then the
efficiency in bins of that offset, all chips together, and the slope per chip.

Usage: eff_vs_baseline.py baselines.csv OUT_STEM eff_runA.json [eff_runB.json ...] [--min-tags 5] [--title ...]
"""
import argparse
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from event_display import style, INK, MUTED, COL  # noqa: E402

PROBES = (('probe_det3', 3), ('probe_det2', 2), ('probe_det1', 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('baselines')
    ap.add_argument('out_stem')
    ap.add_argument('effs', nargs='+')
    ap.add_argument('--min-tags', type=float, default=5)
    ap.add_argument('--title', default='')
    a = ap.parse_args()

    rows = []
    for f in a.effs:
        E = json.load(open(f))
        for key, d in PROBES:
            r = E[key]['r_scale']
            for v, c, _, _, n_on, n_off, k_on, k_off in E[key]['per_pad']:
                rows.append((d, v, c, n_on, n_off, k_on, k_off, r))
    t = pd.DataFrame(rows, columns=['det', 'vmm', 'ch', 'n_on', 'n_off', 'k_on', 'k_off', 'r'])
    t['n'] = t['n_on'] - t['r'] * t['n_off']
    t['k'] = t['k_on'] - t['r'] * t['k_off']
    g = t.groupby(['det', 'vmm', 'ch'], as_index=False)[['n', 'k']].sum()
    g = g[g['n'] >= a.min_tags]
    g['eff'] = (g['k'] / g['n']).clip(0, 1)

    b = pd.read_csv(a.baselines)
    b['baseline'] = b['baseline_adc_HVon'].fillna(b['baseline_adc_HVoff'])
    g = g.merge(b[['vmm', 'ch', 'baseline', 'high_block_flag']], on=['vmm', 'ch'], how='left').dropna(subset=['baseline'])
    g['rel'] = g['baseline'] - g.groupby('vmm')['baseline'].transform('median')
    print(f'{len(g)} channels with >= {a.min_tags:.0f} signal tags and a baseline; tags per channel median {g["n"].median():.0f}')

    edges = np.array([-50, -6, -3, -1, 1, 3, 6, 50])
    g['bin'] = pd.cut(g['rel'], edges)
    tab = g.groupby('bin', observed=True).agg(channels=('eff', 'size'), n=('n', 'sum'), k=('k', 'sum'))
    tab['eff'] = tab['k'] / tab['n']
    tab['err'] = np.sqrt(tab['eff'] * (1 - tab['eff']) / tab['n'])
    print('efficiency by channel baseline minus chip median (ADC), all chips pooled (tags summed):')
    print(tab.round(3).to_string())

    slopes = {}
    for v, gv in g.groupby('vmm'):
        if len(gv) >= 15 and gv['rel'].std() > 0.5:
            w = gv['n']
            k = np.polyfit(gv['rel'], gv['eff'], 1, w=np.sqrt(w))
            r = np.corrcoef(gv['rel'], gv['eff'])[0, 1]
            slopes[int(v)] = (float(k[0]), float(r), int(len(gv)), float(gv['rel'].std()))
    print('per chip: efficiency slope per ADC of baseline, correlation r, channels, baseline spread (ADC):')
    for v, (s, r, n, sd) in slopes.items():
        print(f'   VMM {v:2d}: slope {s:+.4f} /ADC, r {r:+.2f}, {n} ch, spread {sd:.1f} ADC')
    w = np.array([s[2] for s in slopes.values()])
    print(f'mean slope {np.average([s[0] for s in slopes.values()], weights=w):+.4f} /ADC; chips with positive r: '
          f'{sum(s[1] > 0 for s in slopes.values())} of {len(slopes)}')
    g.drop(columns='bin').to_csv(a.out_stem + '_eff_vs_baseline.csv', index=False)

    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
    ax = axs[0]
    style(ax)
    for d in (3, 2, 1):
        gd = g[g['det'] == d]
        ax.scatter(gd['rel'], gd['eff'], s=np.clip(gd['n'], 5, 60), color=COL[d], alpha=0.5, lw=0,
                   label=f'det{d} channels')
    mid = [iv.mid for iv in tab.index]
    ax.errorbar(np.clip(mid, -10, 10), tab['eff'], yerr=tab['err'], fmt='s', color=INK, ms=7, lw=2, label='pooled per bin')
    ax.set_xlim(-12, 12)
    ax.set_xlabel('channel baseline - chip median [ADC]', fontsize=8)
    ax.set_ylabel('efficiency (chamber-only tags, lower bound)', fontsize=8)
    ax.set_title('per channel (marker area ~ tags)', fontsize=9, color=INK)
    ax.legend(fontsize=7, frameon=False)
    ax = axs[1]
    style(ax)
    vs = sorted(slopes)
    ax.bar([str(v) for v in vs], [slopes[v][0] * 100 for v in vs], color=[COL[3] if v < 8 else COL[2] if v < 14 else COL[1] for v in vs])
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xlabel('VMM', fontsize=8)
    ax.set_ylabel('efficiency change per ADC of baseline [% / ADC]', fontsize=8)
    ax.set_title('slope per chip', fontsize=9, color=INK)
    fig.suptitle(f'{a.title}: does the channel baseline predict the efficiency?', fontsize=11)
    fig.savefig(a.out_stem + '_eff_vs_baseline.png', dpi=120)


if __name__ == '__main__':
    sys.exit(main())
