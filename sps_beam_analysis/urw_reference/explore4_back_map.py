#!/usr/bin/env python3
"""Fourth look: is the BACK uRWELL's channel -> strip permutation wrong?

explore3 measured the back uRWELL pointing at 4.44 mm against the P2 planes,
versus 0.77-1.10 mm for the front, and found the best front->back
interpolation fraction to be 0.000 at every plane: the back currently adds
nothing to a track.  4.4 mm is far too large for a 1 mm-pitch strip detector,
so this is a mapping error, not a resolution.

ORDERING.md left exactly one degree of freedom unresolved.  A view is read out
on two 64-channel Dream connectors, and every test available in July separated
"the two connectors are interchanged" from "they are not" but could NOT see a
reversal of channel order WITHIN a connector, because for a uniform map that
reversal is degenerate with an overall axis reflection.  There are therefore
four candidate wirings per view:

    AB      map connector 0 -> low FEU connector,  channel order as-is
    BA      the two connectors interchanged        (what the code applies today)
    AB_rev  as AB, but channel order reversed inside each connector
    BA_rev  as BA, but channel order reversed inside each connector

The front is now an external reference good to <1 mm (explore3), and the beam
is parallel to well under a mrad (best interpolation fraction 0.0), so the back
must reproduce the front position up to a constant.  That breaks the
degeneracy: score each candidate by the width of (back - front).
"""
import os
import sys
import itertools

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import urw_lib as U  # noqa: E402

RUN_DIR = '/local/home/banco/P2_data/TB_July2026_H4/runs/highstat_eff_1'
SUB = 'beam_commissioning_00'
MAX_HITS = 4e6
MODES = ('AB', 'BA', 'AB_rev', 'BA_rev')


def view_connectors(geo, view):
    conns = [c for c in geo.feu_connectors
             if (geo.view[(c - 1) * 64:c * 64] == view).any()]
    if len(conns) != 2:
        raise RuntimeError(f'{geo.name} view {view}: connectors {conns}')
    return sorted(conns)


def apply_mode(geo, view, mode):
    """Re-lay this view's 128 map entries onto its 128 channels per `mode`.

    Operates on a geometry built with connector_swap disabled, so 'AB' is the
    map's own connector column taken at face value.
    """
    lo, hi = view_connectors(geo, view)
    a = slice((lo - 1) * 64, lo * 64)
    b = slice((hi - 1) * 64, hi * 64)
    for arr in (geo.pos, geo.pitch, geo.interpitch, geo.zone):
        blk = [arr[a].copy(), arr[b].copy()]
        if mode.endswith('_rev'):
            blk = [x[::-1] for x in blk]
        if mode.startswith('BA'):
            blk = blk[::-1]
        arr[a], arr[b] = blk[0], blk[1]


def build(name, run_json, x_mode, y_mode):
    geo = U.UrwGeometry(name, run_json, sub_run_name=SUB,
                        view_mode={'x': 'AB', 'y': 'AB'},
                        axis_flip={'x': False, 'y': False})
    apply_mode(geo, 'x', x_mode)
    apply_mode(geo, 'y', y_mode)
    return geo


def points(geo, sub_dir, max_hits=MAX_HITS):
    parts = []
    for chunk in U.iter_hits(U.feu_hit_files(sub_dir, feu=geo.feu_num),
                             max_hits=max_hits, progress=False, feu=geo.feu_num):
        cl = U.cluster_hits(chunk, geo)
        if len(cl):
            parts.append(U.leading_points(cl))
    return pd.concat(parts, ignore_index=True)


def score(a, b):
    """Robust width of b vs a after a linear fit, plus the core fraction."""
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    s, o = np.polyfit(a, b, 1)
    for _ in range(5):
        r = b - (s * a + o)
        k = np.abs(r) < 2.5 * max(r.std(), 1e-6)
        if k.sum() < 100:
            break
        s, o = np.polyfit(a[k], b[k], 1)
    r = b - (s * a + o)
    q = np.percentile(r, [25, 75])
    return dict(slope=s, offset=o, rms=r.std(),
                sigma_iqr=(q[1] - q[0]) / 1.349,
                core=float((np.abs(r - np.median(r)) < 3.0).mean()))


def main():
    run_json = os.path.join(RUN_DIR, 'run_config.json')
    sub_dir = os.path.join(RUN_DIR, SUB)

    # reference: the front with the mapping in use today (connector_swap
    # x=True, y=False -> modes BA / AB in this script's language)
    front = build('EIC_uRWELL_front', run_json, 'BA', 'AB')
    pf = points(front, sub_dir).rename(columns={'x': 'xf', 'y': 'yf'})
    print(f'front reference (BA/AB): {len(pf)} events')

    print('\n--- FRONT control: does the front itself prefer BA/AB? ---')
    _scan('EIC_uRWELL_front', run_json, sub_dir, pf, base=('BA', 'AB'))

    print('\n--- BACK: which permutation reproduces the front? ---')
    _scan('EIC_uRWELL_back', run_json, sub_dir, pf, base=('BA', 'BA'))


def _scan(name, run_json, sub_dir, pf, base):
    for view, other_view in (('x', 'y'), ('y', 'x')):
        print(f'  {name} view {view} (other view held at {base[0 if other_view == "x" else 1]}):')
        rows = []
        for mode in MODES:
            xm = mode if view == 'x' else base[0]
            ym = mode if view == 'y' else base[1]
            geo = build(name, run_json, xm, ym)
            p = points(geo, sub_dir)
            j = pf.merge(p[['eventId', view]], on='eventId').dropna(
                subset=[f'{view}f', view])
            s = score(j[f'{view}f'].to_numpy(), j[view].to_numpy())
            rows.append((mode, len(j), s))
            print(f'    {mode:7s} n={len(j):7d}  slope {s["slope"]:+.4f}  '
                  f'rms {s["rms"]:6.2f}  sigma(IQR) {s["sigma_iqr"]:6.2f} mm  '
                  f'core(|r|<3mm) {s["core"]:.1%}')
        best = min(rows, key=lambda r: r[2]['sigma_iqr'])
        print(f'    -> best: {best[0]}  (sigma {best[2]["sigma_iqr"]:.2f} mm)')


if __name__ == '__main__':
    main()
