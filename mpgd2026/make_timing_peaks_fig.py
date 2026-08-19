#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_timing_peaks_fig.py -- the dt distributions themselves, with their
Gaussian fits, at the best-timing working point of each detector and gas.

Rows are detectors, columns are the npz files given on the command line (one
per gas), each produced by vmm_dream_matching/vmm_timing_peaks.py.  The flat
pedestal under every peak is the ~18-fold BCID ambiguity inside one SRS marker
interval, which the campaign measurement subtracts by sideband; it is drawn
rather than hidden because its height is the honest measure of how clean the
coincidence is.

Usage
-----
  python3 make_timing_peaks_fig.py gasA.npz gasB.npz \
      --col-labels "gas A  Ar/CO2/iC4H10 93/5/2" "gas B  Ar/CF4/iC4H10 88/10/2" \
      --stations P2_MID P2_OUT -o figs
"""

import argparse
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

C = {'P2_IN': '#1f77b4', 'P2_MID': '#ff7f0e', 'P2_OUT': '#2ca02c'}
NL = chr(10)

plt.rcParams.update({
    'font.size': 11.5, 'axes.titlesize': 12, 'axes.labelsize': 11.5,
    'legend.fontsize': 9.5, 'figure.facecolor': 'white',
    'axes.grid': True, 'grid.alpha': 0.3, 'axes.axisbelow': True,
})


def gauss_bg(x, a, mu, s, bg):
    return a * np.exp(-0.5 * ((x - mu) / s) ** 2) + bg


def draw(ax, z, meta, station):
    key = f'{station}_counts'
    if key not in z:
        ax.text(0.5, 0.5, 'no data', transform=ax.transAxes, ha='center')
        return
    h = z[key].astype(float)
    e = z[f'{station}_edges']
    ctr = 0.5 * (e[:-1] + e[1:])
    f = meta.get(station, {})

    ax.step(ctr, h, where='mid', color=C[station], lw=1.6,
            label=f'{int(h.sum()):,} pairs'.replace(',', ' '))
    if f.get('ok'):
        xs = np.linspace(ctr[0], ctr[-1], 800)
        ax.plot(xs, gauss_bg(xs, f['amp'], f['mu'], f['sigma'], f['bg']),
                color='#111111', lw=2.0,
                label=f'Gaussian + flat bg')
        ax.axhline(f['bg'], color='#888', ls=':', lw=1.4)
        lo, hi = f['mu'] - f['sigma'], f['mu'] + f['sigma']
        ax.axvspan(lo, hi, color=C[station], alpha=0.16, zorder=0)
        ax.text(0.03, 0.95,
                f'$\\sigma$ = {f["sigma"]:.1f} ns' + NL
                + f'$\\mu$ = {f["mu"]:.0f} ns' + NL
                + f'peak/bg = {f["peak_over_bg"]:.0f}',
                transform=ax.transAxes, va='top', fontsize=11,
                fontweight='bold', color=C[station],
                bbox=dict(fc='white', ec=C[station], alpha=0.93,
                          boxstyle='round,pad=0.4'))
        ax.set_xlim(f['mu'] - 260, f['mu'] + 260)
    ax.legend(loc='upper right', framealpha=0.95)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('npz', nargs='+')
    ap.add_argument('--col-labels', nargs='*', default=None)
    ap.add_argument('--stations', nargs='*', default=['P2_MID', 'P2_OUT'])
    ap.add_argument('-o', '--out', default='figs')
    ap.add_argument('--name', default='vmm_timing_peaks')
    ap.add_argument('--title', default='VMM3a coincidence time, at each '
                                       'detector\'s best-timing working point')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    files = [(np.load(p, allow_pickle=True)) for p in a.npz]
    metas = [json.loads(str(z['meta'])) for z in files]
    cols = a.col_labels or [m.get('label') or m.get('tag') for m in metas]

    nr, nc = len(a.stations), len(files)
    fig, axes = plt.subplots(nr, nc, figsize=(6.9 * nc, 3.9 * nr),
                             squeeze=False, layout='constrained')
    for i, st in enumerate(a.stations):
        for j, (z, m) in enumerate(zip(files, metas)):
            ax = axes[i][j]
            draw(ax, z, m, st)
            if i == 0:
                ax.set_title(cols[j], fontsize=11.5)
            if i == nr - 1:
                ax.set_xlabel('t(station) - t(trigger), BCID phase [ns]')
            if j == 0:
                ax.set_ylabel(f'{st}' + NL + 'pairs per bin')
    for j, m in enumerate(metas):
        axes[-1][j].text(0.02, 0.03, m.get('label', ''),
                         transform=axes[-1][j].transAxes, fontsize=8.2,
                         color='#666', style='italic', va='bottom')
    fig.suptitle(a.title, fontsize=13.5)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(a.out, f'{a.name}.{ext}'), dpi=170,
                    bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  [{a.name}]')


if __name__ == '__main__':
    main()
