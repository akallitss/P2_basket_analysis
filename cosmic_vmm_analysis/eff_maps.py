#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Maps of the tag-and-probe efficiency of ONE run (output of efficiency_tp.py): is the inefficiency uniform?

Per probe chamber, the per-pad counts (on- and off-time tags, with and without a probe hit) are summed in square cells
of --cell mm and the efficiency is computed per cell as in efficiency_tp.py (off-time subtraction; the small
random-hit correction is left out). Top row: efficiency per cell (cells with fewer than --min-tags signal tags blank);
bottom row: signal tags per cell. Pads of the probe that were not live (masked, dead, vetoed) are drawn as x.
A second figure gives the efficiency per chip and, inside each chip, per group of 8 channels.

Usage: eff_maps.py eff.json OUT_STEM [--cell 40] [--min-tags 8]
"""
import argparse
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import orientation_check as O  # noqa: E402
from event_display import style, INK, MUTED, GRID, COL  # noqa: E402

PROBES = (('probe_det3', 3, 'det3 P2_IN'), ('probe_det2', 2, 'det2 P2_MID'), ('probe_det1', 1, 'det1 P2_OUT'))


def cell_eff(pads, r, cell):
    p = np.array(pads, float).reshape(-1, 8)
    ix, iy = np.floor(p[:, 2] / cell).astype(int), np.floor(p[:, 3] / cell).astype(int)
    out = {}
    for key in set(zip(ix, iy)):
        m = (ix == key[0]) & (iy == key[1])
        n_on, n_off, k_on, k_off = p[m, 4].sum(), p[m, 5].sum(), p[m, 6].sum(), p[m, 7].sum()
        n = n_on - r * n_off
        out[key] = ((k_on - r * k_off) / n if n > 0 else np.nan, n)
    return out


def chip_groups(pads, r, group=8):
    p = np.array(pads, float).reshape(-1, 8)
    res = {}
    for v in np.unique(p[:, 0]).astype(int):
        for g0 in range(0, 64, group):
            m = (p[:, 0] == v) & (p[:, 1] >= g0) & (p[:, 1] < g0 + group)
            n = p[m, 4].sum() - r * p[m, 5].sum()
            k = p[m, 6].sum() - r * p[m, 7].sum()
            res[(v, g0)] = (k / n if n > 0 else np.nan, n)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('eff_json')
    ap.add_argument('out_stem')
    ap.add_argument('--cell', type=float, default=40)
    ap.add_argument('--min-tags', type=float, default=8)
    a = ap.parse_args()
    E = json.load(open(a.eff_json))
    prm = E['params']
    X, Y = O.Mapper(prm['conn_shift']).xy_table(O.parse_orders(prm['order']))
    title = E.get('title') or a.eff_json

    fig, axs = plt.subplots(2, 3, figsize=(16, 7.5), constrained_layout=True)
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('white')
    for col, (key, d, name) in enumerate(PROBES):
        P = E[key]
        ce = cell_eff(P['per_pad'], P['r_scale'], a.cell)
        vm = range(14, 20) if d == 1 else range(8, 14) if d == 2 else range(2, 8)
        px, py = X[list(vm)].ravel(), Y[list(vm)].ravel()
        dead = {tuple(int(x) for x in s.split(':')) for s in E['dead_or_masked'][name]}
        dv = np.array([[X[v, c], Y[v, c]] for v, c in dead if np.isfinite(X[v, c])]).reshape(-1, 2)
        for row in (0, 1):
            ax = axs[row, col]
            style(ax)
            ax.grid(False)
            ax.scatter(px, py, s=2, color=GRID, lw=0)
            if len(dv):
                ax.scatter(dv[:, 0], dv[:, 1], s=10, marker='x', color=MUTED, lw=0.8,
                           label='pad not live' if col == 0 and row == 0 else None)
            vals = [(k, v) for k, v in ce.items() if (v[1] >= a.min_tags or row == 1)]
            norm = plt.Normalize(0, 1) if row == 0 else plt.Normalize(0, max(v[1] for _, v in ce.items()))
            cm = cmap if row == 0 else plt.get_cmap('Blues')
            for (ix, iy), (e, n) in vals:
                val = e if row == 0 else n
                if not np.isfinite(val):
                    continue
                ax.add_patch(plt.Rectangle((ix * a.cell, iy * a.cell), a.cell, a.cell, color=cm(norm(val)), alpha=0.85,
                                           lw=0.5, ec='white'))
                ax.text((ix + 0.5) * a.cell, (iy + 0.5) * a.cell, f'{val:.2f}' if row == 0 else f'{val:.0f}',
                        ha='center', va='center', fontsize=6, color='white' if (row == 0 and val < 0.6) or
                        (row == 1 and norm(val) > 0.5) else INK)
            ax.set_aspect('equal')
            ax.set_xlim(np.nanmin(px) - 20, np.nanmax(px) + 20)
            ax.set_ylim(np.nanmin(py) - 20, np.nanmax(py) + 20)
            ax.set_xlabel('x [mm]', fontsize=8)
            ax.set_ylabel('y [mm]', fontsize=8)
            s = P['all']
            ax.set_title(f'{name}: efficiency per {a.cell:.0f} mm cell (all {s["eff"]:.2f})' if row == 0 else
                         f'{name}: signal tags per cell ({s["n_sig"]:.0f})', fontsize=9, color=INK)
            fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cm), ax=ax, shrink=0.7)
    axs[0, 0].legend(fontsize=7, frameon=False, loc='upper left')
    fig.suptitle(f'{title}: where are the muons missed? (cells with < {a.min_tags:.0f} tags blank)', fontsize=11)
    fig.savefig(a.out_stem + '_effmap.png', dpi=120)

    fig, axs = plt.subplots(1, 3, figsize=(16, 4), constrained_layout=True)
    for col, (key, d, name) in enumerate(PROBES):
        P = E[key]
        cg = chip_groups(P['per_pad'], P['r_scale'])
        ax = axs[col]
        style(ax)
        for (v, g0), (e, n) in sorted(cg.items()):
            if n < 3 or not np.isfinite(e):
                continue
            err = np.sqrt(max(e * (1 - e), 0.02) / n)
            ax.errorbar(v + (g0 + 4) / 64 - 0.5, e, yerr=err, fmt='o', ms=3, color=COL[d], lw=1)
        for v, c in P['per_chip'].items():
            if c['eff'] is not None:
                ax.hlines(c['eff'], int(v) - 0.5, int(v) + 0.5, color=INK, lw=2)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel('chip (VMM); dots = groups of 8 channels inside it, bar = chip', fontsize=8)
        ax.set_ylabel('efficiency', fontsize=8)
        ax.set_title(f'{name}', fontsize=9, color=INK)
    fig.suptitle(f'{title}: efficiency per chip and per 8-channel group', fontsize=11)
    fig.savefig(a.out_stem + '_effchips.png', dpi=120)


if __name__ == '__main__':
    sys.exit(main())
