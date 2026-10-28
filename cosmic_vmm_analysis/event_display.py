#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Event displays of three-chamber groups, and the det2 residual against the chance level -- to look at the tracks.

Input: .npz of dump_candidates.py, groups as in straight_track_rate.py (touching all three chambers, at most
--max-hits hits, hits on mapped pads). For every group the best triplet (one hit per chamber) is the one whose det2 pad
is closest to the midpoint of the det3 and det1 pads, after removing the median det2 offset of the run.

Figures (OUT_*.png):
  residuals  best-triplet residual of every group, for the given map and the old c4-c6 map, against the same with
             the det2 hits taken from another group (chance)
  straight   --n groups drawn at RANDOM among those with residual < --tol: top view x-y and side views x-z, y-z
  random     --n groups drawn at RANDOM among all groups (no selection)
Marker area ~ ADC; det3 circle, det2 square, det1 triangle; grey = pads read out; dashed = det3-det1 line of the
best triplet, x = its midpoint (where det2 should be).

Usage: event_display.py cand.npz OUT_STEM [--conn-shift -1] [--tol 30] [--max-hits 12] [--n 8] [--seed 1]
       [--order 'P2_IN:3=bot:rev|top:lin' ...] [--veto 2:27,6:44]
"""
import argparse
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import orientation_check as O  # noqa: E402
from straight_track_rate import groups, parse_veto  # noqa: E402

Z = {3: 320.0, 2: 630.0, 1: 940.0}
NAME = {3: 'det3 P2_IN', 2: 'det2 P2_MID', 1: 'det1 P2_OUT'}
COL = {3: '#2a78d6', 2: '#eb6834', 1: '#1baf7a'}
MRK = {3: 'o', 2: 's', 1: '^'}
INK, MUTED, GRID = '#1f1f1e', '#6b6a64', '#d8d7d2'


def best_triplets(h, perm=None):
    """per group: residual (rx, ry) of the best triplet and the indices of its three hits; perm as in
    straight_track_rate.count_straight (det2 hits of another group)."""
    evs = np.unique(h['ev'].to_numpy())
    idx = pd.Series(np.arange(len(evs)), index=evs)
    parts = {}
    for d in (1, 2, 3):
        p = h[h['det'] == d]
        g = idx[p['ev']].to_numpy()
        if d == 2 and perm is not None:
            g = perm[g]
        parts[d] = pd.DataFrame({'g': g, f'i{d}': p.index.to_numpy(), f'x{d}': p['x'].to_numpy(),
                                 f'y{d}': p['y'].to_numpy()})
    t = parts[3].merge(parts[1], on='g').merge(parts[2], on='g')
    t['rx'] = t['x2'] - 0.5 * (t['x3'] + t['x1'])
    t['ry'] = t['y2'] - 0.5 * (t['y3'] + t['y1'])
    return t, evs


def shift_of(t):
    near = np.hypot(t['rx'], t['ry']) < 80
    return (float(np.median(t.loc[near, 'rx'])), float(np.median(t.loc[near, 'ry']))) if near.sum() > 20 else (0., 0.)


def best_per_group(t, s):
    t = t.assign(r=np.hypot(t['rx'] - s[0], t['ry'] - s[1]))
    return t.loc[t.groupby('g')['r'].idxmin()].set_index('g')


def style(ax):
    for k in ('top', 'right'):
        ax.spines[k].set_visible(False)
    for k in ('left', 'bottom'):
        ax.spines[k].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=7)
    ax.grid(color=GRID, lw=0.5)
    ax.set_axisbelow(True)


def residual_fig(cand, a, out):
    fig, axs = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    rng = np.random.default_rng(a.seed)
    bins = np.arange(0, 302, 6)
    for col, (shift, label) in enumerate(((a.conn_shift, f'map used (connector shift {a.conn_shift:+d})'),
                                          (0, 'old map c4-c6 (shift 0)'))):
        X, Y = O.Mapper(shift).xy_table(a.choice if shift == a.conn_shift else None)
        h, _ = groups([cand], a.max_hits, X, Y, a.veto_set)
        t, evs = best_triplets(h)
        s = shift_of(t)
        b = best_per_group(t, s)
        ax = axs[col]
        style(ax)
        ax.hist(b['r'], bins=bins, histtype='step', color=INK, lw=2, label=f'data ({len(b)} groups)')
        ch = []
        for _ in range(a.n_shuffle):
            tc, _ = best_triplets(h, rng.permutation(len(evs)))
            ch.append(np.histogram(best_per_group(tc, s)['r'], bins=bins)[0])
        ch = np.mean(ch, 0)
        ax.stairs(ch, bins, color=MUTED, lw=2, ls='--', label='chance (det2 from another group)')
        n_in, c_in = int((b['r'] < a.tol).sum()), float(ch[bins[:-1] < a.tol].sum())
        ax.axvline(a.tol, color=MUTED, lw=1, ls=':')
        ax.set_title(f'{label}\n{n_in} groups < {a.tol:.0f} mm, chance {c_in:.0f}', fontsize=9, color=INK)
        ax.set_xlabel('det2 distance from the det3-det1 midpoint [mm]', fontsize=8, color=INK)
        ax.set_ylabel('groups / 6 mm', fontsize=8, color=INK)
        ax.legend(fontsize=7, frameon=False)
        if col == 0:
            ax2 = axs[2]
            style(ax2)
            ax2.hist2d(b['rx'], b['ry'], bins=[np.arange(-150, 151, 12)] * 2, cmap='Blues')
            ax2.plot(*s, 'x', color=INK, ms=8)
            ax2.set_aspect('equal')
            ax2.set_title(f'data, map used: best-triplet residual\nmedian offset ({s[0]:.0f}, {s[1]:.0f}) mm',
                          fontsize=9, color=INK)
            ax2.set_xlabel('x(det2) - midpoint x [mm]', fontsize=8, color=INK)
            ax2.set_ylabel('y(det2) - midpoint y [mm]', fontsize=8, color=INK)
    fig.suptitle(f'{a.title}: is det2 on the line through det3 and det1?', fontsize=11, color=INK)
    fig.savefig(out + '_residuals.png', dpi=130)


def display(h, b, gsel, X, Y, s, title, out):
    n = len(gsel)
    rows = (n + 1) // 2
    fig, axs = plt.subplots(rows, 6, figsize=(17, 2.9 * rows), constrained_layout=True,
                            gridspec_kw={'width_ratios': [1.4, 1, 1] * 2})
    axs = np.atleast_2d(axs)
    px, py = X[np.isfinite(X)], Y[np.isfinite(Y)]
    evs = np.unique(h['ev'].to_numpy())
    for k, g in enumerate(gsel):
        ax_xy, ax_xz, ax_yz = axs[k // 2, 3 * (k % 2): 3 * (k % 2) + 3]
        hg = h[h['ev'] == evs[g]]
        bt = b.loc[g]
        for ax in (ax_xy, ax_xz, ax_yz):
            style(ax)
        ax_xy.scatter(px, py, s=2, color=GRID, lw=0)
        for d in (3, 2, 1):
            hd = hg[hg['det'] == d]
            sz = 10 + 60 * np.clip(hd['adc'] / 1023, 0, 1)
            kw = dict(s=sz, marker=MRK[d], facecolor=COL[d], edgecolor='white', lw=0.6, zorder=3)
            ax_xy.scatter(hd['x'], hd['y'], label=NAME[d] if k == 0 else None, **kw)
            ax_xz.scatter(hd['x'], np.full(len(hd), Z[d]), **kw)
            ax_yz.scatter(hd['y'], np.full(len(hd), Z[d]), **kw)
        x3, y3, x1, y1 = bt['x3'], bt['y3'], bt['x1'], bt['y1']
        mx, my = 0.5 * (x3 + x1) + s[0], 0.5 * (y3 + y1) + s[1]
        ax_xy.plot([x3, x1], [y3, y1], ls='--', color=INK, lw=1)
        ax_xy.plot(mx, my, 'x', color=INK, ms=7, mew=1.5)
        ax_xz.plot([x3, x1], [Z[3], Z[1]], ls='--', color=INK, lw=1)
        ax_yz.plot([y3, y1], [Z[3], Z[1]], ls='--', color=INK, lw=1)
        ax_xz.plot(mx, Z[2], 'x', color=INK, ms=7, mew=1.5)
        ax_yz.plot(my, Z[2], 'x', color=INK, ms=7, mew=1.5)
        ax_xy.set_aspect('equal')
        ax_xy.set_xlim(px.min() - 15, px.max() + 15)
        ax_xy.set_ylim(py.min() - 15, py.max() + 15)
        dt = hg['t'].max() - hg['t'].min()
        ax_xy.set_title(f'group {int(evs[g])}: {len(hg)} hits, span {dt:.0f} ns, det2 off line {bt["r"]:.0f} mm',
                        fontsize=8, color=INK, loc='left')
        for ax, lab in ((ax_xz, 'x [mm]'), (ax_yz, 'y [mm]')):
            ax.set_ylim(1000, 260)
            ax.set_yticks([Z[3], Z[2], Z[1]])
            ax.set_yticklabels(['IN', 'MID', 'OUT'])
            ax.set_xlabel(lab, fontsize=7, color=INK)
        ax_xy.set_xlabel('x [mm]', fontsize=7, color=INK)
        ax_xy.set_ylabel('y [mm]', fontsize=7, color=INK)
        ax_xz.set_xlim(hg['x'].min() - 40, hg['x'].max() + 40)
        ax_yz.set_xlim(hg['y'].min() - 40, hg['y'].max() + 40)
    for k in range(n, rows * 2):
        for ax in axs[k // 2, 3 * (k % 2): 3 * (k % 2) + 3]:
            ax.set_visible(False)
    fig.legend(loc='outside lower center', ncol=3, fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=11, color=INK)
    fig.savefig(out, dpi=120)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('cand')
    ap.add_argument('out_stem')
    ap.add_argument('--conn-shift', type=int, default=-1)
    ap.add_argument('--tol', type=float, default=30)
    ap.add_argument('--max-hits', type=int, default=12)
    ap.add_argument('--n', type=int, default=8)
    ap.add_argument('--n-shuffle', type=int, default=10)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--title', default='run_121')
    ap.add_argument('--order', action='append', help="channel order of one connector, e.g. 'P2_IN:3=bot:rev|top:lin'")
    ap.add_argument('--veto', default='', help='vmm:ch,... whose groups are dropped, e.g. 2:27,6:44')
    a = ap.parse_args()
    a.choice = O.parse_orders(a.order)
    a.veto_set = parse_veto(a.veto)

    residual_fig(a.cand, a, a.out_stem)
    X, Y = O.Mapper(a.conn_shift).xy_table(a.choice)
    h, _ = groups([a.cand], a.max_hits, X, Y, a.veto_set)
    t, evs = best_triplets(h)
    s = shift_of(t)
    b = best_per_group(t, s)
    rng = np.random.default_rng(a.seed)
    straight = b.index[b['r'] < a.tol].to_numpy()
    pick = np.sort(rng.choice(straight, min(a.n, len(straight)), replace=False))
    display(h, b, pick, X, Y, s, f'{a.title}: {len(pick)} groups drawn at random among the {len(straight)} with det2 '
            f'< {a.tol:.0f} mm off the line (connector shift {a.conn_shift:+d})', a.out_stem + '_straight.png')
    pick = np.sort(rng.choice(b.index.to_numpy(), min(a.n, len(b)), replace=False))
    display(h, b, pick, X, Y, s, f'{a.title}: {len(pick)} groups drawn at random among all {len(b)} three-chamber '
            f'groups (no selection)', a.out_stem + '_random.png')
    print('groups', len(b), 'straight', len(straight), 'det2 offset', s)


if __name__ == '__main__':
    sys.exit(main())
