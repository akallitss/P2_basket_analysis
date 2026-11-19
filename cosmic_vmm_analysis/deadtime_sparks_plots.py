#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plots and tables from the spark_gap_chunk.py outputs (CHUNKS/sg_run<pid>.npz).

Every fraction is on muons: ON-time minus r x OFF-time tags (weights 1 / -r), scintillator window 140:310 ns, as in
silent_muons.py. silent = no live hit anywhere in the probe chamber within +-win of tag + latency.

  fig1_deadtime.png    silent fraction vs time since the previous hit on the predicted chip, against the fraction
                       expected for that chip in that run (removes busy/bad chips and configurations)
  fig2_sparks.png      silent fraction vs time since the previous spark event of the probe chamber (k >= 1, 3, 10
                       saturated hits within 1 us), and vs the time to the NEXT one (control: a spark after the muon
                       cannot blind the chamber), both against the chip-and-run expectation
  fig3_dark.png        length of the hit-free interval containing the muon (window removed), for silent and detected
                       muons and for random times: chamber and hybrid scope
  summary.txt, *.csv   the numbers

Usage: deadtime_sparks_plots.py CHUNKS_DIR PSEUDO_DIR OUT_DIR [--exclude-runs 270]
"""
import argparse
import glob
import json
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
COL = {'P2_IN': '#2a78d6', 'P2_MID': '#eb6834', 'P2_OUT': '#1baf7a'}
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                     'lines.linewidth': 2, 'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb'})


def load(a):
    chunks = {int(k): v for k, v in json.load(open(f'{a.pseudo}/chunks.json')).items()}
    excl = {int(x) for x in a.exclude_runs.split(',') if x}
    T, R, SP = [], [], []
    for f in sorted(glob.glob(f'{a.chunks}/sg_run*.npz')):
        pid = int(re.search(r'run(\d+)\.npz', f)[1])
        run = chunks[pid]['run']
        Z = np.load(f)
        span = float(Z['span_s'])
        for p in (3, 2, 1):
            ev_n = Z[f'ev_p{p}_n']
            SP.append(dict(run=run, pid=pid, det=NAME[p], span_s=span, **{f'ev_k{k}': int((ev_n >= k).sum())
                                                                        for k in (1, 3, 10)}))
            if run in excl:
                continue
            g = lambda k: Z[f'p{p}_{k}']
            r = float(g('r_scale'))
            on = g('on')
            d = dict(run=run, pid=pid, det=NAME[p], on=on, w=np.where(on, 1.0, -r), vmm=g('vmm'),
                     silent=(g('n_any') == 0).astype(float), prev_chip=g('prev_chip'))
            for k in ('spark_prev_1', 'spark_prev_3', 'spark_prev_10', 'spark_next_1', 'spark_next_3',
                      'spark_next_10', 'spark_any_prev_3', 'dark_hyb', 'dark_cham', 'dark_all3'):
                d[k] = g(k)
            T.append(pd.DataFrame(d))
            if f'p{p}_rnd_cham' in Z.files:
                R.append(pd.DataFrame(dict(run=run, det=NAME[p], kind='cham', dark=g('rnd_cham'))))
                for kk in Z.files:
                    if kk.startswith(f'p{p}_rnd_hyb'):
                        R.append(pd.DataFrame(dict(run=run, det=NAME[p], kind='hyb', dark=Z[kk])))
    t = pd.concat(T, ignore_index=True)
    # expected silent fraction of the muon's chip in its run (weighted)
    t['ws'] = t.w * t.silent
    grp = t.groupby(['run', 'vmm'])
    t['base'] = grp.ws.transform('sum') / grp.w.transform('sum')
    return t, pd.concat(R, ignore_index=True), pd.DataFrame(SP)


def binned(t, col, edges):
    """weighted silent fraction, chip-and-run expectation and their difference per bin of col (binomial error)."""
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (t[col] >= lo) & (t[col] < hi)
        W = t.w[m].sum()
        n_on = int((m & t.on).sum())
        if W <= 3:
            rows.append(dict(lo=lo, hi=hi, muons=W, n_on=n_on, silent=np.nan, expected=np.nan, err=np.nan))
            continue
        f = (t.ws[m]).sum() / W
        e = (t.w[m] * t.base[m]).sum() / W
        rows.append(dict(lo=lo, hi=hi, muons=round(W, 1), n_on=n_on, silent=f, expected=e,
                         err=np.sqrt(max(f * (1 - f), 1e-4) / W)))
    df = pd.DataFrame(rows)
    df['excess'] = df.silent - df.expected
    df['lost_muons'] = df.excess * df.muons
    return df


def centers(df):
    return np.sqrt(np.maximum(df.lo, 1) * np.where(np.isinf(df.hi), df.lo * 10, df.hi))


def tick_time(ax):
    ticks = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11]
    labs = ['10 ns', '100 ns', '1 µs', '10 µs', '100 µs', '1 ms', '10 ms', '100 ms', '1 s', '10 s', '100 s']
    lo, hi = ax.get_xlim()
    k = [i for i, x in enumerate(ticks) if lo <= x <= hi]
    ax.set_xticks([ticks[i] for i in k]); ax.set_xticklabels([labs[i] for i in k])


def panel(ax, t, col, edges, title, note=None, min_mu=8):
    """all chambers together (ink, with error bars) and per chamber (thin, coloured); dotted = expectation."""
    df = binned(t, col, edges)
    x = centers(df); ok = df.muons > min_mu
    for det in NAME.values():
        d = binned(t[t.det == det], col, edges)
        k = d.muons > min_mu
        ax.plot(x[k], d.silent[k], color=COL[det], marker='o', ms=3.5, lw=1.1, alpha=0.9, label=det)
    ax.errorbar(x[ok], df.silent[ok], df.err[ok], color=INK, marker='o', ms=6, lw=2, capsize=0,
                label='all chambers')
    ax.plot(x[ok], df.expected[ok], color=INK, ls=':', lw=1.6, label='expected for that chip and run')
    for xx, n, y in zip(x[ok], df.muons[ok], df.silent[ok]):
        ax.text(xx, 0.015, f'{n:.0f}', ha='center', fontsize=7, color=INK2)
    ax.set_xscale('log'); tick_time(ax)
    ax.set_ylim(0, 0.8); ax.set_ylabel('silent fraction (muons)')
    ax.set_title(title, loc='left', fontsize=10, color=INK)
    if note:
        ax.text(0.02, 0.96, note, transform=ax.transAxes, ha='left', va='top', fontsize=9, color=INK)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('chunks'); ap.add_argument('pseudo'); ap.add_argument('out')
    ap.add_argument('--exclude-runs', default='270')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t, R, SP = load(a)
    log = open(f'{a.out}/summary.txt', 'w')

    def say(*x):
        s = ' '.join(str(v) for v in x); print(s, flush=True); log.write(s + '\n')

    W = t.w.sum()
    say(f'runs {sorted(t.run.unique())} (excluded {a.exclude_runs}); chunks {t.pid.nunique()}; tags {len(t)} '
        f'({int(t.on.sum())} on); muons {W:.0f}; silent fraction {t.ws.sum() / W:.4f}')
    for det in NAME.values():
        s = t[t.det == det]
        say(f'  {det}: muons {s.w.sum():.0f}, silent {s.ws.sum() / s.w.sum():.4f}')

    # ---------------------------------------------------------------- fig 1: dead time
    E1 = [0, 100, 300, 1e3, 3e3, 1e4, 1e5, 1e6, 1e7, 1e8, np.inf]
    d1 = binned(t, 'prev_chip', E1)
    d1.to_csv(f'{a.out}/deadtime_prev_chip.csv', index=False)
    lost = d1[d1.hi <= 1e4].lost_muons.sum()
    say('\n== silent vs time since the previous hit on the predicted chip (ns), all chambers')
    say(d1.round(4).to_string(index=False))
    say(f'muons lost to the < 10 us excess: {lost:.1f} of {W:.0f} ({100 * lost / W:.2f} %); '
        f'all silent muons {t.ws.sum():.0f} ({100 * lost / t.ws.sum():.1f} % of them)')
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    panel(ax, t, 'prev_chip', [10, 100, 300, 1e3, 3e3, 1e4, 3e4, 1e5, 1e6, 1e7, 1e8, 1e9], '',
          note=f'loss below 10 µs: {lost:.0f} of {W:.0f} muons ({100 * lost / W:.2f} %)')
    ax.set_xlim(10, 3e9); tick_time(ax)
    ax.set_xlabel('time since the previous hit on the predicted chip (before the window opens); grey numbers: muons per bin')
    ax.set_title('Dead time: is the chamber silent because the chip had just fired?', loc='left', fontsize=11)
    ax.legend(loc='upper right', fontsize=8.5, ncol=1)
    fig.tight_layout(); fig.savefig(f'{a.out}/fig1_deadtime.png', dpi=150); plt.close(fig)

    # ---------------------------------------------------------------- fig 2: sparks
    E2 = [2e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, np.inf]
    say('\n== spark events (saturated hits of a chamber grouped with a 1 us gap) per second, per run and chamber')
    sp = SP.groupby(['run', 'det'])[['ev_k1', 'ev_k3', 'ev_k10', 'span_s']].sum()
    for k in (1, 3, 10):
        sp[f'rate_k{k}_Hz'] = (sp[f'ev_k{k}'] / sp.span_s).round(3)
    say(sp[[c for c in sp.columns if c.startswith('rate')]].unstack('det').to_string())
    sp.to_csv(f'{a.out}/spark_rates.csv')
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.6), sharey=True)
    for ax, (col, ttl) in zip(axs.flat, [
            ('spark_prev_3', 'after a spark (≥ 3 saturated hits) in the same chamber'),
            ('spark_next_3', 'CONTROL: before the next spark (≥ 3) in the same chamber'),
            ('spark_prev_10', 'after a big spark (≥ 10 saturated hits) in the same chamber'),
            ('spark_prev_1', 'after any saturated hit (≥ 1) in the same chamber')]):
        panel(ax, t, col, [2e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12], ttl)
        ax.set_xlim(2e3, 1e12); tick_time(ax)
        ax.set_xlabel('time since the spark' if 'prev' in col else 'time until the spark')
    axs[0, 0].legend(loc='upper right', fontsize=8)
    fig.suptitle('Sparks: is the chamber blind after a discharge? (muon\'s own signal ±2 µs excluded)',
                 x=0.01, ha='left', fontsize=11)
    fig.tight_layout(); fig.savefig(f'{a.out}/fig2_sparks.png', dpi=150); plt.close(fig)
    for col in ('spark_prev_1', 'spark_prev_3', 'spark_prev_10', 'spark_next_3', 'spark_any_prev_3'):
        d = binned(t, col, E2)
        d.to_csv(f'{a.out}/{col}.csv', index=False)
        say(f'\n== silent vs {col} (ns), all chambers')
        say(d.round(4).to_string(index=False))

    # ---------------------------------------------------------------- fig 3: dark intervals
    EB = np.logspace(2, 10, 41)
    fig, axs = plt.subplots(2, 3, figsize=(12, 6.6), sharex=True)
    say('\n== hit-free interval containing the muon (window removed): median and fraction > 100 us / > 1 ms')
    rows = []
    for j, det in enumerate(NAME.values()):
        for i, (scope, rk) in enumerate((('dark_cham', 'cham'), ('dark_hyb', 'hyb'))):
            ax = axs[i, j]
            s = t[(t.det == det) & t.on]
            rr = R[(R.det == det) & (R.kind == rk)].dark
            for lab, x, c, ls in (('silent muons', s[s.silent == 1][scope], INK, '-'),
                                  ('detected muons', s[s.silent == 0][scope], COL[det], '-'),
                                  ('random times', rr, INK2, ':')):
                x = np.clip(x.to_numpy(), 101, 9.9e9)
                h, _ = np.histogram(x, EB)
                ax.stairs(h / max(h.sum(), 1), EB, color=c, ls=ls, lw=1.8, label=f'{lab} ({len(x)})')
                rows.append(dict(det=det, scope=rk, sample=lab, n=len(x), median_us=np.median(x) / 1e3,
                                 gt100us=(x > 1e5).mean(), gt1ms=(x > 1e6).mean()))
            ax.set_xscale('log')
            ax.set_xticks([1e3, 1e5, 1e7, 1e9]); ax.set_xticklabels(['1 µs', '100 µs', '10 ms', '1 s'])
            ax.set_title(f'{det}: {"whole chamber" if rk == "cham" else "hybrid of the predicted pad"}',
                         loc='left', fontsize=10)
            if j == 0:
                ax.set_ylabel('fraction per bin')
            if i == 1:
                ax.set_xlabel('length of the hit-free interval around the muon')
            ax.legend(fontsize=7.5, loc='upper right')
    fig.suptitle('Read-out going dark? Hit-free interval that contains the muon (ON-time tags)', x=0.01, ha='left',
                 fontsize=11)
    fig.tight_layout(); fig.savefig(f'{a.out}/fig3_dark.png', dpi=150); plt.close(fig)
    dk = pd.DataFrame(rows).round(4)
    dk.to_csv(f'{a.out}/dark_intervals.csv', index=False)
    say(dk.to_string(index=False))
    # all three chambers dark at once around a silent muon?
    s = t[t.on & (t.silent == 1)]
    say(f'\nsilent ON muons: {len(s)}; the three chambers together hit-free over > 100 us around them: '
        f'{(s.dark_all3 > 1e5).mean():.3f}; detected ON muons: {(t[t.on & (t.silent == 0)].dark_all3 > 1e5).mean():.3f}')
    say(f'\noutputs in {a.out}')


if __name__ == '__main__':
    main()
