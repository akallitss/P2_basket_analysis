#!/usr/bin/env python3
"""Fifth look: is the BACK uRWELL's position scale piecewise-wrong?

With the mapping fixed (urw_lib VIEW_MODE_DEFAULT), the back reproduces the
front to a ~0.8 mm core, but the fitted front->back straight line still comes
out with a slope of +1.04 in y (and the fit is unstable: different clipping
gives 1.029 or 1.042), while x sits at +0.9995.  A single global scale error
would be stable; an unstable one says the relation is not a straight line.

The back's map is a three-zone patchwork per view - strips 0-63 at 1.0 mm
pitch, 64-95 at 1.5 mm, 96-127 at 0.5 mm - so if one zone's pitch or offset is
wrong the position-vs-strip relation is piecewise linear and a global fit
splits the difference.  This plots (back - front) against the front position
and against the back strip index, split by the back's zone, so a bad zone shows
up as a step or a wrong local slope.

The front, whose map is uniform 1.0 mm, is the reference throughout.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import urw_lib as U  # noqa: E402

RUN_DIR = '/local/home/banco/P2_data/TB_July2026_H4/runs/highstat_eff_1'
SUB = 'beam_commissioning_00'
MAX_HITS = 4e6


def leading_with_zone(geo, sub_dir, max_hits):
    """Leading cluster per view per event, keeping the cluster's zone."""
    parts = []
    for chunk in U.iter_hits(U.feu_hit_files(sub_dir, feu=geo.feu_num),
                             max_hits=max_hits, progress=False, feu=geo.feu_num):
        cl = U.cluster_hits(chunk, geo)
        if not len(cl):
            continue
        cl = cl.sort_values('charge', ascending=False)
        lead = cl.groupby(['eventId', 'view'], sort=False).first().reset_index()
        parts.append(lead[['eventId', 'view', 'pos', 'zone', 'size', 'charge']])
    return pd.concat(parts, ignore_index=True)


def main():
    run_json = os.path.join(RUN_DIR, 'run_config.json')
    sub_dir = os.path.join(RUN_DIR, SUB)
    f = leading_with_zone(U.UrwGeometry('EIC_uRWELL_front', run_json,
                                        sub_run_name=SUB), sub_dir, MAX_HITS)
    b = leading_with_zone(U.UrwGeometry('EIC_uRWELL_back', run_json,
                                        sub_run_name=SUB), sub_dir, MAX_HITS)
    geo_b = U.UrwGeometry('EIC_uRWELL_back', run_json, sub_run_name=SUB)
    print('back zone table:')
    print(geo_b.zone_table().to_string(index=False))

    fig, ax = plt.subplots(2, 3, figsize=(18, 9))
    for row, view in enumerate(('x', 'y')):
        fv = f[f['view'] == view][['eventId', 'pos']].rename(columns={'pos': 'pf'})
        bv = b[b['view'] == view][['eventId', 'pos', 'zone', 'size']].rename(
            columns={'pos': 'pb'})
        j = fv.merge(bv, on='eventId')
        j['d'] = j['pb'] - j['pf']
        # kill the bulk offset so the shape is visible
        j['d'] -= j['d'].median()
        core = j[j['d'].abs() < 8]
        print(f'\nview {view}: {len(j)} pairs, {len(core)} within 8 mm')
        for z, g in core.groupby('zone'):
            if len(g) < 500:
                continue
            s, o = np.polyfit(g['pf'], g['pb'], 1)
            print(f'   back zone {z}: n={len(g):7d}  pb range '
                  f'{g["pb"].min():6.1f}..{g["pb"].max():6.1f}  '
                  f'local slope {s:+.4f}  offset {o:+.2f} mm  '
                  f'median (b-f) {g["d"].median():+.2f}')

        a = ax[row, 0]
        hh = a.hist2d(core['pf'], core['d'], bins=[130, 80], cmin=1)
        fig.colorbar(hh[3], ax=a)
        a.set(xlabel=f'front {view} [mm]', ylabel='back - front [mm]',
              title=f'{view}: residual vs front position')

        a = ax[row, 1]
        hh = a.hist2d(core['pb'], core['d'], bins=[130, 80], cmin=1)
        for z, g in core.groupby('zone'):
            if len(g) >= 500:
                a.axvline(g['pb'].min(), color='r', lw=0.8, ls='--')
        fig.colorbar(hh[3], ax=a)
        a.set(xlabel=f'back {view} [mm]', ylabel='back - front [mm]',
              title=f'{view}: residual vs back position (zone edges dashed)')

        a = ax[row, 2]
        g = core.groupby(pd.cut(core['pb'], 128), observed=True)['d']
        m = g.median()
        a.plot([iv.mid for iv in m.index], m.values, '.-', ms=4)
        a.axhline(0, color='k', lw=0.8)
        for z, gg in core.groupby('zone'):
            if len(gg) >= 500:
                a.axvline(gg['pb'].min(), color='r', lw=0.8, ls='--')
        a.set(xlabel=f'back {view} [mm]', ylabel='median (back - front) [mm]',
              title=f'{view}: profile')
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig('explore5_back_zones.png', dpi=95)
    print('\nwrote explore5_back_zones.png')


if __name__ == '__main__':
    main()
