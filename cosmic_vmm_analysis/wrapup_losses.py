#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Wrap-up of the 19 Nov 2026 analysis: pedagogical figures of the three "is it the read-out?" tests and the
side-by-side table of the efficiency definitions. Reads only summary products already on EOS (no decoding).

  fig_w1_readout.png    test pulses received per chip (pulser runs 270 at gain 4.5 and 283 at gain 9),
                        and how the missing pulses are grouped (isolated misses vs dark periods)
  fig_w2_deadtime.png   fraction of chamber-silent muons vs time since the previous hit on the predicted chip,
                        against what that chip's own rate predicts; number of muons per bin below
  fig_w3_sparks.png     the same vs time since the previous spark (>= 3 and >= 10 saturated hits in 1 us)
  fig_w4_budget.png     the three losses next to the gap to the Dream efficiency (same chambers, same connectors)
  eff_definitions.csv / .md   per configuration and chamber: counting eps, A20, beam-style HIT, ANY, reco-free
                        scintillator any-hit, with the Dream M3 bench and SPS beam reference rows

Usage (lxplus, from $W = .../weekend_20261119, with pipeline/env.sh sourced):
  wrapup_losses.py $W OUT_DIR
"""
import argparse
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
CH_COL = {'P2_IN': BLUE, 'P2_MID': ORANGE, 'P2_OUT': AQUA}
plt.rcParams.update({'font.size': 12, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8, 'legend.frameon': False,
                     'figure.facecolor': 'white', 'savefig.dpi': 150, 'savefig.bbox': 'tight'})

# Dream references on the same chambers and connectors (handoff 20261119d, finding 6):
# M3 bench from cosmic_bench_analysis stage 06 ray_hit_miss_list*.csv (has_any, active area);
# beam from sps_beam_analysis/dream_vmm_comparison/data/compare_dream_vmm_P2_OUT.csv
DREAM_BENCH = {'P2_OUT': (0.984, 0.968, 0.943), 'P2_MID': (0.923, 0.938, 0.942), 'P2_IN': (0.891, 0.890, 0.902)}
BEAM = {'P2_OUT c5': {'Dream': 0.960, 'VMM': 0.872}}
SCINT_RUNS = {268: 'trimmed_Km18', 269: 'trimmed_Km18', 275: 'bscan_hvon10', 278: 'gain6.0_gscan_hvon10',
              281: 'gain9.0_gscan_hvon10'}


def note(ax, txt, xy=(0.02, 0.96), **kw):
    ax.text(*xy, txt, transform=ax.transAxes, va='top', ha='left', fontsize=11, color=INK, **kw)


# ---------------------------------------------------------------- 1. read-out (test pulses)
GAP_BINS = ((1, 1), (2, 4), (5, 20), (21, 10 ** 9))


def fig_readout(W, out):
    P = os.path.join(W, 'analysis/pulser_tick225')
    runs = [(270, 'results_run270', 'run 270, gain 4.5 (REF)', BLUE), (283, 'results_run283b', 'run 283, gain 9', ORANGE)]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={'width_ratios': [1.1, 1]})
    chips = [2, 5, 8, 11, 14, 17]
    det = {2: 'IN', 5: 'IN', 8: 'MID', 11: 'MID', 14: 'OUT', 17: 'OUT'}
    rec = {}
    for k, (run, d, lab, col) in enumerate(runs):
        C = pd.read_csv(os.path.join(P, d, 'pulser_channels.csv'))
        g = C.groupby('chip')[['pulses', 'expected']].sum()
        f = (g.pulses / g.expected).reindex(chips)
        rec[run] = f
        x = np.arange(len(chips)) + (k - 0.5) * 0.28
        a1.bar(x, 100 * (1 - f.values), width=0.26, color=col, label=lab)
        G = pd.read_csv(os.path.join(P, d, 'pulser_gaps.csv'))
        G = G[G.gap_frames > 1].assign(lost=lambda q: (q.gap_frames - 1) * q.n)
        for j, chip in enumerate((2, 8)):
            g8 = G[G.chip == chip]
            tot = g8.lost.sum()
            fr = [g8.lost[(g8.gap_frames - 1 >= lo) & (g8.gap_frames - 1 <= hi)].sum() / tot for lo, hi in GAP_BINS]
            xx = np.arange(len(GAP_BINS)) + (2 * j + k - 1.5) * 0.2
            a2.bar(xx, fr, width=0.19, color=col, alpha=1.0 if chip == 8 else 0.45,
                   label=f'{lab}, chip {chip}')
    a1.set_xticks(np.arange(len(chips)), [f'chip {c}\n{det[c]}' for c in chips])
    a1.set_ylabel('test pulses NOT received (%)')
    a1.set_title('Read-out: test pulses lost per chip (ch 31)', loc='left', color=INK)
    a1.axhline(10, color=INK2, ls='--', lw=1)
    a1.text(-0.4, 10.3, 'needed to explain the ~10-point gap', ha='left', va='bottom', fontsize=10, color=INK2)
    a1.set_ylim(0, 14)
    for k, run in enumerate(rec):
        for i, v in enumerate(rec[run].values):
            a1.text(i + (k - 0.5) * 0.28, 100 * (1 - v) + 0.2, f'{100 * (1 - v):.2f}', ha='center', fontsize=8,
                    color=INK2, rotation=90, va='bottom')
    a1.legend(loc='upper right')
    a2.set_xticks(np.arange(len(GAP_BINS)), ['1 (isolated)', '2-4', '5-20\n(<= 2 ms)', '> 20'])
    a2.set_xlabel('consecutive pulses missing (92 us apart)')
    a2.set_ylabel('share of the missing pulses')
    a2.set_ylim(0, 1)
    a2.set_title('How the missing pulses are grouped', loc='left', color=INK)
    a2.legend(fontsize=9, loc='upper right')
    fig.text(0.01, -0.06, 'Gain 4.5: 99.82-99.87 % received on all six chips - the read-out chain does not lose data.\n'
             'Gain 9: chips 5, 8, 11 lose 1-12 %, the same in every one of the 49 files, as isolated pulses and '
             'short bursts (<= 2 ms), no long dark periods.\nPossibly the test-pulse amplitude near the gain-9 '
             'threshold (higher in charge) - not yet checked; it does not change the gain-4.5 conclusion.',
             fontsize=10.5, color=INK2, va='top')
    fig.savefig(os.path.join(out, 'fig_w1_readout.png'))
    plt.close(fig)
    return float(1 - rec[270].mean())


# ---------------------------------------------------------------- 2./3. silent fraction vs a time variable
def silent_vs_time(ax, axn, D, col, lab, unit_us=True):
    D = D[D.n_on >= 20].copy()
    x = np.sqrt(D.lo.clip(lower=10) * np.where(np.isfinite(D.hi), D.hi, D.lo * 10)) / 1e3   # us, geometric centre
    ok = D.silent.notna()
    ax.errorbar(x[ok], D.silent[ok], D.err[ok], fmt='o', ms=7, color=col, lw=2, capsize=0, label=lab)
    ax.plot(x[ok], D.expected[ok], color=INK2, lw=2, ls='--', drawstyle='default')
    axn.bar(x, D.muons, width=x * 0.9, color=col, alpha=0.85)
    return x


def fig_deadtime(R, out, n_total):
    D = pd.read_csv(os.path.join(R, 'deadtime_prev_chip.csv'))
    fig, (ax, axn) = plt.subplots(2, 1, figsize=(9, 7), sharex=True, gridspec_kw={'height_ratios': [2.2, 1]})
    silent_vs_time(ax, axn, D, BLUE, 'measured (all three chambers)')
    ax.plot([], [], color=INK2, ls='--', lw=2, label="expected from the chip's own rate (no dead time)")
    ax.axvspan(0.03, 10, color=YELLOW, alpha=0.12, lw=0)
    ax.set_xscale('log'); axn.set_yscale('log')
    ax.set_ylabel('fraction of muons with\nthe chamber silent')
    ax.set_ylim(0, 0.8)
    ax.set_title('Dead time: is a muon missed because the chip was still busy?', loc='left', color=INK)
    lost = D.lost_muons[D.hi <= 1e4].sum()
    note(ax, f'excess only below 10 us: {lost:.0f} of {n_total} muons = {100 * lost / n_total:.2f} %\n'
             'most muons (lower panel) come > 10 us after the previous hit', xy=(0.33, 0.95))
    ax.legend(loc='center right', fontsize=10)
    axn.set_xlabel('time since the previous hit on the predicted chip (us)')
    axn.set_ylabel('muons per bin')
    fig.savefig(os.path.join(out, 'fig_w2_deadtime.png'))
    plt.close(fig)
    return float(lost / n_total)


def fig_sparks(R, out, n_total):
    fig, (ax, axn) = plt.subplots(2, 1, figsize=(9, 7), sharex=True, gridspec_kw={'height_ratios': [2.2, 1]})
    lost = 0.0
    for k, col, lab in ((3, BLUE, 'after a spark of >= 3 saturated hits in 1 us'),
                        (10, ORANGE, 'after a big spark (>= 10 saturated hits)')):
        D = pd.read_csv(os.path.join(R, f'spark_prev_{k}.csv'))
        silent_vs_time(ax, axn, D, col, lab)
        lost = max(lost, D.lost_muons[(D.hi <= 1e8) & D.lost_muons.notna()].sum())
    ax.plot([], [], color=INK2, ls='--', lw=2, label='expected (no effect of the spark)')
    ax.set_xscale('log'); axn.set_yscale('log')
    ax.set_ylim(0, 0.5)
    ax.set_ylabel('fraction of muons with\nthe chamber silent')
    ax.set_title('Sparks: does a discharge blind the chamber?', loc='left', color=INK)
    note(ax, f'small excess 1-100 ms after a spark:\n<= {lost:.0f} of {n_total} muons = {100 * lost / n_total:.2f} %\n'
             'sparks are 0.2-1 /s per chamber,\nso few muons come close to one', xy=(0.04, 0.95))
    ax.legend(loc='center right', fontsize=10)
    axn.set_xlabel('time since the previous spark in the same chamber (us)')
    axn.set_ylabel('muons per bin')
    fig.savefig(os.path.join(out, 'fig_w3_sparks.png'))
    plt.close(fig)
    return float(lost / n_total)


# ---------------------------------------------------------------- 4. the budget
def fig_budget(out, l_ro, l_dt, l_sp, conn_any):
    fig, ax = plt.subplots(figsize=(10, 4.6))
    items = [('read-out (test pulses lost)', 100 * l_ro, BLUE), ('dead time (< 10 us after a hit)', 100 * l_dt, BLUE),
             ('sparks (1-100 ms after)', 100 * l_sp, BLUE)]
    gaps = []
    for d in ('P2_OUT', 'P2_MID', 'P2_IN'):
        v = conn_any.get(d)
        if v is None:
            continue
        gap = 100 * (np.mean(DREAM_BENCH[d]) - np.mean(v))
        gaps.append((f'gap to Dream, {d} (connectors 3-5)', gap, ORANGE))
    rows = items + [('all three together', 100 * (l_ro + l_dt + l_sp), INK2)] + gaps
    y = np.arange(len(rows))[::-1]
    for yy, (lab, v, col) in zip(y, rows):
        ax.barh(yy, v, height=0.62, color=col)
        ax.text(v + 0.15, yy, f'{v:.2f}' if v < 1 else f'{v:.1f}', va='center', fontsize=11, color=INK)
    ax.set_yticks(y, [r[0] for r in rows])
    ax.set_xlabel('efficiency points')
    ax.grid(axis='y', visible=False)
    ax.set_title('What the three tests can explain vs the gap to Dream', loc='left', color=INK)
    fig.text(0.01, -0.04, 'Gap = mean over connectors 3-5 of (Dream M3 bench any-hit - VMM cosmic any-hit), REF '
             'configuration. The losses are measured on the VMM\nitself and are < 0.6 points together: '
             'the 3-9 point gap is not lost data (consistent with signal below the VMM threshold, beam autopsy run_46).', fontsize=10.5,
             color=INK2, va='top')
    fig.savefig(os.path.join(out, 'fig_w4_budget.png'))
    plt.close(fig)


# ---------------------------------------------------------------- the definitions table
def conn_key(k):
    g, name = k.split(' ', 1)
    if not name.endswith('_tpA'):
        name = re.sub(r'(_v1)?_hvm$', '', name)
    return f'gain{g[1:]}_{name}' if name.startswith('gscan') else name


def read_conn(path):
    """chamber-level HIT and ANY (muon-weighted over connectors) and per-connector ANY, per configuration."""
    rows, per = [], {}
    for line in open(path):
        m = re.match(r'^(g\S+ \S+)\s+(IN .*)$', line.strip())
        if not m:
            continue
        cfg = conn_key(m[1])
        for cell in m[2].split('|'):
            q = re.match(r'\s*(IN|MID|OUT) c(\d) ([\d.]+)/([\d.]+) \((\d+)\)', cell)
            if q:
                rows.append(dict(config=cfg, det=f'P2_{q[1]}', conn=int(q[2]), hit=float(q[3]), anyhit=float(q[4]),
                                 n=int(q[5])))
    C = pd.DataFrame(rows)
    agg = C.assign(wh=C.hit * C.n, wa=C.anyhit * C.n).groupby(['config', 'det'])[['wh', 'wa', 'n']].sum()
    agg['HIT15'] = agg.wh / agg.n
    agg['ANY'] = agg.wa / agg.n
    return agg[['HIT15', 'ANY', 'n']].reset_index(), C


def read_scint(path):
    rows, run = [], None
    for line in open(path):
        m = re.match(r'^run_(\d+)', line)
        if m:
            run = int(m[1])
        m = re.match(r'\s*140:310 \(now\)\s*\|(.*)\(n probe: (\d+)\)', line)
        if m and run in SCINT_RUNS:
            parts = [p.split() for p in m[1].split('|')]
            for d, p in zip(('P2_IN', 'P2_MID', 'P2_OUT'), parts):
                rows.append(dict(run=run, config=SCINT_RUNS[run], det=d, scint_probe=float(p[1]), n=int(m[2])))
    S = pd.DataFrame(rows)
    S['w'] = S.scint_probe * S.n
    g = S.groupby(['config', 'det'])[['w', 'n']].sum()
    return (g.w / g.n).rename('SCINT_ANY').reset_index()


def table(W, out):
    cur = os.path.join(W, 'results/current')
    B = pd.read_csv(os.path.join(cur, 'by_config.csv'))
    B = B[B.config != 'config'][['config', 'det', 'eps', 'err']].astype({'eps': float, 'err': float})
    A = pd.read_csv(os.path.join(cur, 'chamber_level.csv'))[['config', 'det', 'A20', 'A20_err', 'tags']]
    H, C = read_conn(os.path.join(W, 'analysis/pulser_tick225/conn_eff.log'))
    S = read_scint(os.path.join(W, 'analysis/pulser_tick225/scint_any_hit.log'))
    T = B.merge(A, on=['config', 'det'], how='outer').merge(H, on=['config', 'det'], how='outer') \
         .merge(S, on=['config', 'det'], how='left')
    order = ['trimmed_Km18', 'trimmed_Km18_v1_hvm_tpA', 'bscan_hvon2', 'bscan_hvon10', 'bscan_hvon50',
             'gain3.0_gscan_hvon2', 'gain3.0_gscan_hvon10', 'gain3.0_gscan_hvon50', 'gain6.0_gscan_hvon2',
             'gain6.0_gscan_hvon10', 'gain6.0_gscan_hvon50', 'gain9.0_gscan_hvon2', 'gain9.0_gscan_hvon10',
             'gain9.0_gscan_hvon10_hvm_tpA', 'gain9.0_gscan_hvon50']
    T['config'] = pd.Categorical(T.config, [c for c in order if c in set(T.config)] +
                                 sorted(set(T.config) - set(order)))
    T['det'] = pd.Categorical(T.det, ['P2_IN', 'P2_MID', 'P2_OUT'])
    T = T.sort_values(['config', 'det']).rename(columns={'eps': 'COUNT', 'err': 'COUNT_err', 'n': 'conn_muons'})
    T.to_csv(os.path.join(out, 'eff_definitions.csv'), index=False, float_format='%.4f')

    def f(v, e=None):
        if pd.isna(v):
            return '-'
        return f'{v:.3f}' + (f' ± {e:.3f}' if e is not None and not pd.isna(e) else '')
    L = ['| config | chamber | COUNT (tag-and-probe, density model) | A20 (scint tags, pad ≤ 20 mm) | '
         'HIT15 (beam def.: ≤ 15 mm, ±150 ns) | ANY (any chamber hit, ±100 ns) | SCINT_ANY (no reco, MID bracketed) |',
         '|---|---|---|---|---|---|---|']
    for _, r in T.iterrows():
        L.append(f'| {r.config} | {r.det} | {f(r.COUNT, r.COUNT_err)} | {f(r.A20, r.A20_err)} | {f(r.HIT15)} | '
                 f'{f(r.ANY)} | {f(r.SCINT_ANY)} |')
    L += ['', '**Dream references on the same pads (any-hit; connectors 3 / 4 / 5):**', '',
          '| chamber | Dream M3 bench | VMM cosmic ANY, REF (c3 / c4 / c5) | SPS beam, P2_OUT c5 |', '|---|---|---|---|']
    ref = C[C.config == 'trimmed_Km18'].set_index(['det', 'conn']).anyhit
    conn_any = {}
    for d in ('P2_IN', 'P2_MID', 'P2_OUT'):
        v = [ref.get((d, c), np.nan) for c in (3, 4, 5)]
        conn_any[d] = v
        beam = 'Dream 0.960 / VMM 0.872' if d == 'P2_OUT' else '-'
        L.append(f'| {d} | {" / ".join(f"{x:.3f}" for x in DREAM_BENCH[d])} | '
                 f'{" / ".join(f"{x:.3f}" for x in v)} | {beam} |')
    open(os.path.join(out, 'eff_definitions.md'), 'w').write('\n'.join(L) + '\n')
    return conn_any


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('W'); ap.add_argument('out')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    R = os.path.join(a.W, 'analysis/deadtime_sparks/results')
    n_total = int(re.search(r'muons (\d+);', open(os.path.join(R, 'summary.txt')).readline())[1])
    l_ro = fig_readout(a.W, a.out)
    l_dt = fig_deadtime(R, a.out, n_total)
    l_sp = fig_sparks(R, a.out, n_total)
    conn_any = table(a.W, a.out)
    fig_budget(a.out, l_ro, l_dt, l_sp, conn_any)
    print(f'losses: read-out {100 * l_ro:.3f} %, dead time {100 * l_dt:.3f} %, sparks {100 * l_sp:.3f} %')
    print(open(os.path.join(a.out, 'eff_definitions.md')).read())


if __name__ == '__main__':
    main()
