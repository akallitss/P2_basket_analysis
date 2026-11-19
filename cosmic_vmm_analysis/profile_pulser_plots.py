#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig4_spark_profile.png: hit rate after a spark (and before it) relative to the mean rate, per chamber, from the
spark_profile_chunk.py outputs; fig5_pulser.png + pulser_summary.csv: test pulses received against the real duration
(pcapng packet clock) per chip and file, and whether the pulsed chips were alive around the silent muons.

Usage: profile_pulser_plots.py PROFILE_DIR PULSER_DIR OUT_DIR
"""
import glob
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from deadtime_sparks_plots import COL, INK, INK2, tick_time  # noqa: F401  (also sets the rcParams)

NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}


def profiles(pdir, pulser_only):
    acc = {}
    edges = None
    for f in sorted(glob.glob(f'{pdir}/sp_run*.npz')):
        Z = np.load(f)
        if pulser_only is not None and bool(Z['pulser']) != pulser_only:
            continue
        edges = Z['edges']
        span_w = None
        for k in Z.files:
            if k.startswith(('after_', 'before_', 'expo_')):
                acc[k] = acc.get(k, 0) + Z[k]
            elif k.startswith('n_ev_'):
                acc[k] = acc.get(k, 0) + int(Z[k])
            elif k.startswith('mean_rate_'):
                # exposure-weighted mean rate: weight by the chunk's after-exposure of the widest bin (~ duration)
                acc.setdefault(k, []).append(float(Z[k]))
    for k in list(acc):
        if k.startswith('mean_rate_'):
            acc[k] = float(np.mean(acc[k]))
    return acc, edges


def main():
    pdir, qdir, out = sys.argv[1:4]
    os.makedirs(out, exist_ok=True)
    log = open(f'{out}/summary_profile_pulser.txt', 'w')

    def say(*x):
        s = ' '.join(str(v) for v in x); print(s, flush=True); log.write(s + '\n')

    # ------------------------------------------------------------------ fig 4
    rows = []
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.4), sharex=True)
    for col_i, (label, pulser_only, scopes) in enumerate((
            ('cosmic runs (no pulser)', False, ('cham', 'other')),
            ('pulser run 270', True, ('cham', 'other', 'pulse')))):
        acc, edges = profiles(pdir, pulser_only)
        if not acc:
            continue
        x = np.sqrt(np.maximum(edges[:-1], 100) * edges[1:])
        for K_i, K in enumerate((3, 10)):
            for j, p in enumerate((3, 2, 1)):
                ax = axs[K_i, j]
                n_ev = acc.get(f'n_ev_K{K}_p{p}', 0)
                for sc, ls, mk in (('cham', '-', 'o'), ('other', '--', 's'), ('pulse', ':', 'D')):
                    if sc not in scopes or f'after_K{K}_{sc}_p{p}' not in acc:
                        continue
                    mr = acc[f'mean_rate_{sc}_p{p}']
                    ra = acc[f'after_K{K}_{sc}_p{p}'] / np.maximum(acc[f'expo_after_K{K}_{sc}_p{p}'], 1) / mr
                    rb = acc[f'before_K{K}_{sc}_p{p}'] / np.maximum(acc[f'expo_before_K{K}_{sc}_p{p}'], 1) / mr
                    ea = np.sqrt(np.maximum(acc[f'after_K{K}_{sc}_p{p}'], 1)) / np.maximum(
                        acc[f'expo_after_K{K}_{sc}_p{p}'], 1) / mr
                    for i in range(len(x)):
                        rows.append(dict(sample=label, K=K, det=NAME[p], scope=sc, lo_ns=edges[i], hi_ns=edges[i + 1],
                                         sparks=n_ev, counts_after=acc[f'after_K{K}_{sc}_p{p}'][i],
                                         after_rel=ra[i], after_err=ea[i], before_rel=rb[i]))
                    if col_i == 0 and sc == 'pulse':
                        continue
                    c = COL[NAME[p]] if col_i == 0 else INK
                    ok = acc[f'after_K{K}_{sc}_p{p}'] >= 5
                    lab = {'cham': 'same chamber', 'other': 'other two chambers', 'pulse': 'test pulses, same chamber'}[sc]
                    ax.errorbar(x[ok], ra[ok], ea[ok], color=c, ls=ls, marker=mk, ms=4, lw=1.5,
                                label=f'{lab} ({"cosmic" if col_i == 0 else "pulser run"})')
                ax.axhline(1, color=INK2, lw=1)
                ax.set_xscale('log'); ax.set_yscale('log')
                ax.set_ylim(3e-3, 300)
                ax.set_title(f'{NAME[p]}: after a spark with ≥ {K} saturated hits', loc='left', fontsize=10)
                if j == 0:
                    ax.set_ylabel('hit rate after the spark / mean rate')
                if K_i == 1:
                    ax.set_xlabel('time after the end of the spark')
    for ax in axs.flat:
        ax.set_xlim(150, 1.5e9); tick_time(ax)
        ax.set_xticks([1e3, 1e5, 1e7, 1e9]); ax.set_xticklabels(['1 µs', '100 µs', '10 ms', '1 s'])
    axs[0, 0].legend(fontsize=7.5, loc='lower left')
    axs[0, 1].legend(fontsize=7.5, loc='lower left')
    fig.suptitle('After a spark: does the read-out go blind (rate → 0) or recover? (1 = the normal rate)', x=0.01,
                 ha='left', fontsize=11)
    fig.tight_layout(); fig.savefig(f'{out}/fig4_spark_profile.png', dpi=150); plt.close(fig)
    pr = pd.DataFrame(rows)
    pr.to_csv(f'{out}/spark_profile.csv', index=False)
    say('== hit rate after a spark / mean rate (K = 3), per chamber and scope')
    if len(pr):
        q = pr[pr.K == 3].copy()
        q['bin'] = q.lo_ns.map(lambda v: f'{v:g}') + '-' + q.hi_ns.map(lambda v: f'{v:g}')
        say(q.pivot_table(index=['sample', 'scope', 'det'], columns='lo_ns', values='after_rel').round(2).to_string())
        say(pr.groupby(['sample', 'K', 'det']).sparks.first().to_string())

    # ------------------------------------------------------------------ fig 5
    F, M = [], []
    for f in sorted(glob.glob(f'{qdir}/ps_run*.json')):
        j = json.load(open(f))
        F += j['files']; M += j['muons']
    if F:
        F = pd.DataFrame(F); M = pd.DataFrame(M)
        F['eff'] = F.pulses / F.expected
        per_chip = F.groupby('chip').agg(files=('file', 'size'), pulses=('pulses', 'sum'), expected=('expected', 'sum'),
                                         real_s=('real_s', 'sum'), decoded_s=('decoded_s', 'sum'))
        per_chip['received'] = per_chip.pulses / per_chip.expected
        per_chip['decoded_over_real'] = per_chip.decoded_s / per_chip.real_s
        say('\n== pulser run 270: test pulses received against the real duration (pcapng packet clock)')
        say(per_chip.round(5).to_string())
        per_chip.to_csv(f'{out}/pulser_summary.csv')
        fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw=dict(width_ratios=[1.4, 1]))
        ax = axs[0]
        for v, g in F.groupby('chip'):
            ax.plot(np.arange(len(g)), g.eff.to_numpy(), marker='o', ms=3, lw=1, label=f'chip {v}')
        ax.set_ylim(0.95, 1.01); ax.set_xlabel('capture file of run 270 (60 s each); chip 8 > 1: extra noise hits on its pulsed channel')
        ax.set_ylabel('test pulses received / expected')
        ax.set_title('Read-out efficiency from the test pulses (real time)', loc='left', fontsize=10)
        ax.legend(ncol=3, fontsize=8, loc='lower left')
        ax = axs[1]
        M['w'] = np.where(M.on, 1.0, -M.r)
        rows = []
        for (p, s), g in M.groupby(['det', 'silent']):
            alive = ((g.pulsedA > 0) & (g.pulsedB > 0)).astype(float)
            W = g.w.sum()
            rows.append(dict(det=NAME[p], muons='silent' if s else 'detected', n=round(W, 1),
                             both_pulsed_chips_alive=(g.w * alive).sum() / W if W > 0 else np.nan))
        al = pd.DataFrame(rows)
        say('\n== muons of run 270: both pulsed chips of the chamber had hits within +-100 us')
        say(al.round(4).to_string(index=False))
        al.to_csv(f'{out}/pulser_alive_at_muons.csv', index=False)
        xs = np.arange(3)
        for k, (mu, c) in enumerate((('silent', INK), ('detected', INK2))):
            s = al[al.muons == mu].set_index('det').reindex(['P2_IN', 'P2_MID', 'P2_OUT'])
            ax.bar(xs + (k - 0.5) * 0.36, s.both_pulsed_chips_alive, 0.34, color=c, label=f'{mu} muons')
            for xx, (val, n) in zip(xs + (k - 0.5) * 0.36, zip(s.both_pulsed_chips_alive, s.n)):
                ax.text(xx, val + 0.01, f'{val:.2f}\n({n:.0f})', ha='center', va='bottom', fontsize=8, color=INK)
        ax.set_xticks(xs); ax.set_xticklabels(['P2_IN', 'P2_MID', 'P2_OUT']); ax.set_ylim(0, 1.25)
        ax.set_ylabel('both pulsed chips read out ±100 µs')
        ax.set_title('Was the read-out alive when the chamber missed the muon?', loc='left', fontsize=10)
        ax.legend(fontsize=8, loc='upper right')
        fig.tight_layout(); fig.savefig(f'{out}/fig5_pulser.png', dpi=150); plt.close(fig)
    say(f'\noutputs in {out}')


if __name__ == '__main__':
    main()
