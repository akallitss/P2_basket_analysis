#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Channel by channel, chip by chip: the baseline from every neighbour read-out group, and the noise of each channel.

For every VMM one figure, channel 0-63 on the horizontal axis:
  top     the baseline (median ADC of the neighbour read-outs) in each group, with the error of that median
          (1.25 * spread / sqrt(read-outs));
  bottom  the noise of the channel in each group: the spread (std) of its neighbour read-outs, and how often it fired
          by itself (over-threshold hits per neighbour read-out, a rate that does not depend on the run length).
Then one summary: per chip, the 10-90 % range of the baselines of group --ref with all channels and without the
noisy ones (neighbour spread above --noisy-std ADC), and the correlation between a channel's noise and how far its
baseline sits from the chip median.

Usage: baseline_per_vmm.py OUT_DIR --groups A:dirA,C:dirC,... [--ref A] [--min-nb 50] [--noisy-std 10]
"""
import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from event_display import style, INK, MUTED, COL  # noqa: E402

GCOL = {'A': MUTED, 'B': '#7a3e9d', 'C': COL[2], 'D': COL[1]}
GLAB = {'A': 'A: 23 Oct, HV off', 'B': 'B: 23 Oct, HV off, VMM 4/7/9', 'C': 'C: 23 Oct, HV on', 'D': 'D: 31 Oct, new connectors'}


def load(d):
    out = {}
    for r in csv.DictReader(open(os.path.join(d, 'neighbor_baseline.csv'))):
        if float(r['adc_nb_median'] or 0) < 20:          # masked channel read out as a neighbour (ADC 0): no baseline
            continue
        out[(int(r['vmm']), int(r['ch']))] = {'n_nb': int(r['n_nb']), 'n_over': int(r['n_over']),
                                              'med': float(r['adc_nb_median']) if r['adc_nb_median'] else np.nan,
                                              'std': float(r['adc_nb_std']) if r['adc_nb_std'] else np.nan}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out_dir')
    ap.add_argument('--groups', required=True)
    ap.add_argument('--ref', default='A')
    ap.add_argument('--min-nb', type=int, default=50)
    ap.add_argument('--noisy-std', type=float, default=10)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    G = {g.split(':')[0]: load(g.split(':')[1]) for g in a.groups.split(',')}
    ch = np.arange(64)
    summary = []
    for v in range(2, 20):
        fig, (ax, bx) = plt.subplots(2, 1, figsize=(13, 7), sharex=True, constrained_layout=True,
                                     gridspec_kw={'height_ratios': [1.3, 1]})
        style(ax)
        style(bx)
        bx2 = bx.twinx()
        for i, (g, D) in enumerate(G.items()):
            rows = [D.get((v, c)) for c in ch]
            ok = np.array([r is not None and r['n_nb'] >= a.min_nb for r in rows])
            if not ok.any():
                continue
            med = np.array([r['med'] if r else np.nan for r in rows])
            sd = np.array([r['std'] if r else np.nan for r in rows])
            nnb = np.array([r['n_nb'] if r else 0 for r in rows], float)
            nov = np.array([r['n_over'] if r else 0 for r in rows], float)
            x = ch + (i - (len(G) - 1) / 2) * 0.18
            err = 1.25 * sd / np.sqrt(np.maximum(nnb, 1))
            ax.errorbar(x[ok], med[ok], yerr=err[ok], fmt='o', ms=3.5, color=GCOL[g], lw=1, label=f'{GLAB[g]} ({ok.sum()} ch)')
            ax.axhline(np.median(med[ok]), color=GCOL[g], lw=0.8, ls='--')
            bx.plot(x[ok], sd[ok], 'o', ms=3.5, color=GCOL[g])
            bx2.plot(x[ok], (nov / np.maximum(nnb, 1))[ok], 'x', ms=4, color=GCOL[g], alpha=0.7)
        ax.axhline(38, color=INK, lw=0.8, ls=':')
        ax.set_ylabel('baseline [ADC] (median of the neighbour read-outs)', fontsize=8)
        top = max([np.nanmax([G[g][(v, c)]['med'] for c in ch if (v, c) in G[g] and G[g][(v, c)]['n_nb'] >= a.min_nb] or [60])
                   for g in G])
        ax.set_ylim(25, min(max(top + 5, 80), 260))
        ax.legend(fontsize=7.5, frameon=False, ncol=4, loc='upper left')
        ax.set_title(f'VMM {v}: baseline of every channel (dashed: group median; dotted: 38)', fontsize=10, color=INK)
        bx.set_ylabel('spread of the neighbour read-outs [ADC] (dots)', fontsize=8)
        bx2.set_ylabel('own hits per neighbour read-out (x)', fontsize=8)
        bx2.set_yscale('symlog', linthresh=0.01)
        bx2.set_ylim(0, None)
        bx.axhline(a.noisy_std, color=INK, lw=0.8, ls=':')
        bx.set_xlabel('channel', fontsize=8)
        bx.set_xticks(range(0, 64, 4))
        fig.savefig(os.path.join(a.out_dir, f'baseline_vmm{v:02d}.png'), dpi=110)
        plt.close(fig)

        R = G[a.ref]
        rows = [(c, R[(v, c)]) for c in ch if (v, c) in R and R[(v, c)]['n_nb'] >= a.min_nb and R[(v, c)]['med'] < 90]
        if len(rows) >= 10:
            med = np.array([r['med'] for _, r in rows])
            sd = np.array([r['std'] for _, r in rows])
            quiet = sd <= a.noisy_std
            dev = np.abs(med - np.median(med))
            q_all = np.percentile(med, [10, 90])
            q_quiet = np.percentile(med[quiet], [10, 90]) if quiet.sum() >= 5 else (np.nan, np.nan)
            rr = float(np.corrcoef(sd, dev)[0, 1]) if len(rows) > 3 and sd.std() > 0 else np.nan
            summary.append((v, len(rows), int((~quiet).sum()), q_all, q_quiet, rr))
    print(f'reference group {a.ref}: per chip, baseline 10-90 % range with all channels and without the noisy ones '
          f'(neighbour spread > {a.noisy_std} ADC); correlation of noise with |baseline - chip median|')
    for v, n, nn, qa, qq, rr in summary:
        print(f'VMM {v:2d}: {n:2d} ch, {nn:2d} noisy; range all {qa[1] - qa[0]:5.1f} ADC ({qa[0]:.0f}-{qa[1]:.0f}), '
              f'quiet only {qq[1] - qq[0]:5.1f} ADC; r(noise, deviation) {rr:+.2f}')

    fig, ax = plt.subplots(figsize=(12, 4.3), constrained_layout=True)
    style(ax)
    for v, n, nn, qa, qq, rr in summary:
        ax.plot([v - 0.12, v - 0.12], qa, color=MUTED, lw=6, solid_capstyle='round', label='all channels' if v == summary[0][0] else None)
        ax.plot([v + 0.12, v + 0.12], qq, color=COL[2], lw=6, solid_capstyle='round', label=f'without the noisy ones (spread > {a.noisy_std:g} ADC)' if v == summary[0][0] else None)
        ax.text(v, max(qa[1], qq[1] if np.isfinite(qq[1]) else 0) + 1.5, f'{nn}', fontsize=8, ha='center', color=INK)
    ax.set_xticks(range(2, 20))
    ax.set_xlabel('VMM (number above: noisy channels removed)', fontsize=8)
    ax.set_ylabel(f'baseline, 10-90 % of the channels [ADC], group {a.ref}', fontsize=8)
    ax.set_title('Does the noise of a few channels make the baselines of a chip look spread?', fontsize=10, color=INK)
    ax.legend(fontsize=8, frameon=False)
    fig.savefig(os.path.join(a.out_dir, 'baseline_spread_summary.png'), dpi=120)


if __name__ == '__main__':
    sys.exit(main())
