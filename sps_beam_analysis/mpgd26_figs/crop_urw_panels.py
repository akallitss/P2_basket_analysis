#!/usr/bin/env python3
"""Split the 8-panel uRWELL-reference QA sheet into slide-sized figures.

`urw_panels_<station>.png` is a diagnostic contact sheet: eight panels, each
legible at a desk and none of them legible from the back of a room.  For the
talk it has to become two or three figures, each making one argument.

This crops rather than re-renders, and that is a deliberate limitation: the
per-track residual arrays and the uRWELL-frame track density behind those
panels are NOT in the local workspace (only the summary statistics are), so
they cannot be redrawn without re-running the uRWELL-referenced stage on
lxplus with a histogram dump.  Cropping still roughly doubles the effective
font size on a slide, and the composed titles carry the chamber identity that
the original panels never did.

  urw_residuals_<station>.png    1D residual + 2D residual
  urw_efficiency_<station>.png   efficiency in the uRWELL frame + in the pad
                                 frame
  urw_acceptance_<station>.png   track density + matching distance

Usage:  python3 crop_urw_panels.py [station ...]
"""
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2style as st  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import chamber_history as ch  # noqa: E402

from paths import OUT  # noqa: E402

RUN = 'highstat_eff_1'      # the sheet is made from the working-point run
# the sheet is a 2x4 grid under a title band; the four columns are even
# quarters of the plot area and the row gutter is a clean blank band
TITLE_H = 83
ROW_SPLIT = (500, 508)
PAD = 18


def _grid(im):
    h, w = im.shape[:2]
    xs = np.linspace(PAD, w - PAD, 5).astype(int)
    return xs, (TITLE_H, ROW_SPLIT[0]), (ROW_SPLIT[1], h - 14)


def _trim(t, thr=246, pad=6):
    """Drop the blank margin the grid crop leaves around each panel."""
    g = t.mean(axis=2) if t.ndim == 3 else t
    rows = np.where((g < thr).any(axis=1))[0]
    cols = np.where((g < thr).any(axis=0))[0]
    if not len(rows) or not len(cols):
        return t
    r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad + 1, g.shape[0])
    c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad + 1, g.shape[1])
    return t[r0:r1, c0:c1]


def _panel(im, xs, ys, c0, c1, row):
    return _trim(im[ys[row][0]:ys[row][1], xs[c0]:xs[c1]])


def compose(tiles, title, sub, path, width=15.0):
    n = len(tiles)
    hs = [t.shape[0] for t in tiles]
    ws = [t.shape[1] for t in tiles]
    fig, axes = plt.subplots(1, n, figsize=(width, width * max(hs) / sum(ws)
                                            * 1.06),
                             gridspec_kw=dict(width_ratios=ws, wspace=0.02))
    axes = np.atleast_1d(axes)
    for ax, t in zip(axes, tiles):
        ax.imshow(t)
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
    fig.suptitle(f'{title}\n{sub}', fontweight='bold', fontsize=16, y=1.0)
    fig.savefig(path, bbox_inches='tight', dpi=170)
    fig.savefig(os.path.splitext(path)[0] + '.pdf', bbox_inches='tight')
    plt.close(fig)
    print('wrote', path, '(+pdf)')


def main(stations):
    for station in stations:
        src = f'{OUT}/urw_panels_{station}.png'
        if not os.path.exists(src):
            print('! missing', src)
            continue
        im = np.asarray(Image.open(src).convert('RGB'))
        xs, r0, r1 = _grid(im)
        ys = (r0, r1)
        lab = ch.label(station, RUN)

        compose([_panel(im, xs, ys, 0, 1, 0), _panel(im, xs, ys, 1, 2, 0)],
                f'{lab} — how well the reference track and the P2 cluster agree',
                'left: the 1D residual, both projections.  right: the same in '
                '2D — the square is the 12 mm pad, the fringe is the '
                'accidental floor',
                f'{OUT}/urw_residuals_{station}.png')

        compose([_panel(im, xs, ys, 1, 2, 1), _panel(im, xs, ys, 2, 3, 1)],
                f'{lab} — the absolute efficiency, mapped two ways',
                'left: in the uRWELL frame, so any structure here belongs to '
                'the REFERENCE.  right: in the P2 pad frame, where the big '
                'mesh pillar shows up',
                f'{OUT}/urw_efficiency_{station}.png')

        compose([_panel(im, xs, ys, 0, 1, 1), _panel(im, xs, ys, 3, 4, 1)],
                f'{lab} — what went into the number',
                'left: where the reference tracks landed.  right: distance '
                'from track to nearest P2 cluster, and the 15 mm probe radius',
                f'{OUT}/urw_acceptance_{station}.png')


if __name__ == '__main__':
    main(sys.argv[1:] or ['P2_MID'])
