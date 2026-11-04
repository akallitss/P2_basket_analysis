#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Surface maps of the scintillator-confirmed efficiency on the three pad planes (scint_efficiency.py output).

Per probe chamber, the tags (predicted point from the clean clusters of the other two chambers, scintillator
coincidence) are binned in square cells of --cell mm; per cell the efficiency is the fraction of tags with a probe hit
(--hit any: anywhere in the chamber, the headline number; w50: within 50 mm of the prediction). Top row: efficiency
per cell (cells with fewer than --min-tags tags blank, the number of tags under the value); bottom row: the same on
the live area only (tags whose predicted point has every nearby pad live). Pads are grey dots, pads not live
(masked, silent) are red crosses.

Several json of the same settings (e.g. the two K0 runs) can be given comma-separated: their tags are summed.

Usage: scint_eff_maps.py scintor_run187.json[,scintor_run196.json] OUT_STEM [--cell 40] [--min-tags 10] [--hit any]
"""
import argparse
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from event_display import style, INK, GRID  # noqa: E402

DETS = ((3, 'det3 P2_IN'), (2, 'det2 P2_MID'), (1, 'det1 P2_OUT'))
DEAD = '#c0392b'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('jsons')
    ap.add_argument('out_stem')
    ap.add_argument('--cell', type=float, default=40)
    ap.add_argument('--min-tags', type=int, default=10)
    ap.add_argument('--hit', default='any', choices=('any', 'w50'))
    ap.add_argument('--title', default='')
    a = ap.parse_args()
    Es = [json.load(open(f)) for f in a.jsons.split(',')]
    title = a.title or ' + '.join(E['title'].split()[0] for E in Es)

    fig, axs = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True)
    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(0, 1)
    for col, (d, nm) in enumerate(DETS):
        Ps = [E[f'probe_{nm}'] for E in Es]
        t = {k: np.concatenate([np.array(P['tags_xy'][k], float) for P in Ps]) for k in ('x', 'y', 'live', 'any', 'w50')}
        hit = t[a.hit].astype(bool)
        live_t = t['live'].astype(bool)
        pads = Ps[-1]['pads']
        gx, gy = np.array(pads['x']), np.array(pads['y'])
        # a pad is drawn as not live if it was not live in any of the runs summed
        glive = np.all([np.array(P['pads']['live'], bool) for P in Ps], axis=0)
        x0, x1 = gx.min() - 15, gx.max() + 15
        y0, y1 = gy.min() - 15, gy.max() + 15
        ex = np.arange(np.floor(x0 / a.cell) * a.cell, x1 + a.cell, a.cell)
        ey = np.arange(np.floor(y0 / a.cell) * a.cell, y1 + a.cell, a.cell)
        for row, m, lab in ((0, np.ones(len(hit), bool), 'all tags'), (1, live_t, 'live area only')):
            ax = axs[row, col]
            style(ax)
            ax.grid(False)
            N, _, _ = np.histogram2d(t['x'][m], t['y'][m], bins=(ex, ey))
            K, _, _ = np.histogram2d(t['x'][m & hit], t['y'][m & hit], bins=(ex, ey))
            for i in range(len(ex) - 1):
                for j in range(len(ey) - 1):
                    if N[i, j] < a.min_tags:
                        continue
                    e = K[i, j] / N[i, j]
                    ax.add_patch(plt.Rectangle((ex[i], ey[j]), a.cell, a.cell, color=cmap(norm(e)), lw=0.5, ec='white'))
                    c = 'white' if e < 0.6 else INK
                    ax.text(ex[i] + a.cell / 2, ey[j] + a.cell * 0.58, f'{e:.2f}', ha='center', va='center', fontsize=6.5, color=c)
                    ax.text(ex[i] + a.cell / 2, ey[j] + a.cell * 0.28, f'{N[i, j]:.0f}', ha='center', va='center', fontsize=5, color=c)
            ax.scatter(gx[glive], gy[glive], s=1.5, color=GRID, lw=0, zorder=3)
            ax.scatter(gx[~glive], gy[~glive], s=9, marker='x', color=DEAD, lw=0.7, zorder=4,
                       label=f'pad not live ({(~glive).sum()} of {len(glive)})')
            n, k = int(m.sum()), int((m & hit).sum())
            e = k / n if n else np.nan
            ax.set_title(f'{nm}, {lab}: {e:.3f} +- {np.sqrt(max(e * (1 - e), 0) / max(n, 1)):.3f} ({n} tags)',
                         fontsize=9, color=INK)
            ax.set_aspect('equal')
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            ax.set_xlabel('x [mm]', fontsize=8)
            ax.set_ylabel('y [mm]', fontsize=8)
            ax.legend(fontsize=7, frameon=False, loc='upper right')
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs, shrink=0.5,
                 label='efficiency per cell (small number: tags)')
    hit_lab = 'probe hit anywhere' if a.hit == 'any' else 'probe hit within 50 mm'
    fig.suptitle(f'{title}: scintillator-confirmed efficiency per {a.cell:.0f} mm cell at the predicted point '
                 f'({hit_lab}; cells with < {a.min_tags} tags blank)', fontsize=11)
    fig.savefig(a.out_stem + '.png', dpi=120)
    print('wrote', a.out_stem + '.png')


if __name__ == '__main__':
    sys.exit(main())
