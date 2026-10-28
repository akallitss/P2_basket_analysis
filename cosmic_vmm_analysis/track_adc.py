#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pulse height of the hits on straight tracks against the threshold edge: is the efficiency threshold-limited?

Straight tracks as in event_display.py (best triplet within --tol mm of the line); their three hits are the muon hits
(about 15 % of the groups are chance, see the residual plot). For each chamber: ADC of the muon hits, against the ADC
of all hits of that chamber in the candidate file (mostly noise), whose lower edge marks where the threshold cuts.
A muon spectrum (Landau) that starts right at that edge and has no rising side below its peak means part of it is cut
away by the threshold: the efficiency is then set by signal against threshold.

Usage: track_adc.py cand.npz OUT_STEM [--conn-shift -1] [--order ...] [--veto ...] [--tol 15] [--title run_124]
"""
import argparse
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import orientation_check as O  # noqa: E402
from straight_track_rate import groups, parse_veto  # noqa: E402
from event_display import best_triplets, shift_of, best_per_group, style, INK, MUTED, COL  # noqa: E402

NAME = {3: 'det3 P2_IN', 2: 'det2 P2_MID', 1: 'det1 P2_OUT'}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('cand')
    ap.add_argument('out_stem')
    ap.add_argument('--conn-shift', type=int, default=-1)
    ap.add_argument('--order', action='append')
    ap.add_argument('--veto', default='')
    ap.add_argument('--tol', type=float, default=15)
    ap.add_argument('--title', default='')
    a = ap.parse_args()
    X, Y = O.Mapper(a.conn_shift).xy_table(O.parse_orders(a.order))
    h, _ = groups(a.cand.split(','), 12, X, Y, parse_veto(a.veto))
    t, _ = best_triplets(h)
    b = best_per_group(t, shift_of(t))
    b = b[b['r'] < a.tol]
    z = np.load(a.cand.split(',')[0])
    allv, alla = z['vmm'].astype(int), z['adc'].astype(int)
    det_all = z['det'].astype(int)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    bins = np.arange(0, 1040, 16)
    print(f'{len(b)} straight tracks (det2 within {a.tol:.0f} mm)')
    for col, d in enumerate((3, 2, 1)):
        adc = h.loc[b[f'i{d}'], 'adc'].to_numpy()
        noise = alla[det_all == d]
        edge = np.percentile(noise, 1)
        ax = axs[col]
        style(ax)
        ax.hist(noise, bins=bins, histtype='step', color=MUTED, lw=1.5, density=True, label=f'all hits ({len(noise)})')
        ax.hist(adc, bins=bins, histtype='step', color=COL[d], lw=2.2, density=True, label=f'hits on tracks ({len(adc)})')
        ax.axvline(edge, color=INK, ls=':', lw=1, label=f'threshold edge (1st pct of all hits) {edge:.0f}')
        pk = bins[np.argmax(np.histogram(adc, bins=bins)[0])]
        frac_low = float(np.mean(adc < edge + 32))
        ax.set_title(f'{NAME[d]}: track hits median {np.median(adc):.0f}, peak bin {pk}, '
                     f'{100 * frac_low:.0f} % within 32 ADC of the edge', fontsize=8.5, color=INK)
        ax.set_xlabel('ADC (largest pad of the hit)', fontsize=8)
        ax.set_ylabel('normalised', fontsize=8)
        ax.legend(fontsize=7, frameon=False)
        print(f'{NAME[d]}: track hits n {len(adc)}, ADC p10/50/90 {np.percentile(adc, [10, 50, 90]).round(0)}, '
              f'saturated {np.mean(adc >= 1023):.2f}; all hits edge {edge:.0f}, p50 {np.median(noise):.0f}; '
              f'per chip edge: ' + ', '.join(f'{v}:{np.percentile(alla[allv == v], 1):.0f}'
                                              for v in sorted(set(allv[det_all == d]))))
    fig.suptitle(f'{a.title}: pulse height of muon hits against the threshold edge', fontsize=11)
    fig.savefig(a.out_stem + '_track_adc.png', dpi=120)


if __name__ == '__main__':
    sys.exit(main())
