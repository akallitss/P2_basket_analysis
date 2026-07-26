#!/usr/bin/env python3
"""Sixth look: the leftover few-percent scale is a DIVERGENT BEAM, not a bad map.

Two loose ends pointed the same way:
  * the front->back fit needs a slope of ~+1.04 in y but +0.9995 in x, and
    explore5 showed that y slope is a smooth ramp across the whole plane, not a
    step at a pitch-zone edge - so it is not the back's three-zone map;
  * the free affine uRWELL -> P2 comes out orthogonal but with its larger
    singular value growing with z: 1.0071 (z=320), 1.0133 (z=630), 1.0230
    (z=940), while the other stays at 1.000.

An anisotropic magnification growing linearly with z is what a beam diverging
from a virtual source upstream looks like, focused in one plane and not the
other - normal for an SPS extraction line.  If that is the cause then fitting
    A(z) = R . diag(1 + z/Lx, 1 + z/Ly)
with the front alone as the source must give a consistent Lx, Ly across all
three stations, and 1 + 1370/Ly must reproduce the measured front->back y
slope.  A bad strip map cannot do that: it would give the same wrong scale at
every z.
"""
import os
import sys

import numpy as np
import pandas as pd

# this package, then sps_beam_analysis for sps_config, whose setup_paths()
# adds the shared core in cosmic_bench_analysis (p2_io, p2_mapping, ...)
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
import sps_config as _sc  # noqa: E402
_sc.setup_paths()

import explore_urw_p2_frame as E   # noqa: E402
from explore2_resid import affine_window  # noqa: E402
import sps_config as sc            # noqa: E402
import p2_io as p2io               # noqa: E402

N_CHUNK = 1


def stretch(A):
    """Polar decomposition A = R . S with S symmetric, expressed in the SOURCE
    (uRWELL) basis.  Returns (R, S)."""
    u, s, vt = np.linalg.svd(A)
    return u @ vt, vt.T @ np.diag(s) @ vt


def main():
    run_json = os.path.join(E.RUN_DIR, 'run_config.json')
    sub_dir = os.path.join(E.RUN_DIR, E.SUB)
    ev, dz = E.urw_tracks(run_json, sub_dir, N_CHUNK)

    print('front -> back straight-line fit (no scale removed):')
    fb = {}
    for ax in ('x', 'y'):
        s, o = np.polyfit(ev[f'{ax}f'], ev[f'{ax}b'], 1)
        for _ in range(5):
            r = ev[f'{ax}b'] - (s * ev[f'{ax}f'] + o)
            m = np.abs(r) < 2.5 * r.std()
            s, o = np.polyfit(ev.loc[m, f'{ax}f'], ev.loc[m, f'{ax}b'], 1)
        fb[ax] = s
        print(f'   {ax}: slope {s:+.5f}  offset {o:+.2f} mm')

    cfg = sc.RunConfig(key='x', run='highstat_eff_1',
                       data_root='/local/home/banco/P2_data/TB_July2026_H4/runs')
    allf = p2io.hit_files(cfg.combined_hits_dir(E.SUB))[:N_CHUNK]

    print('\nfront-only -> P2, stretch tensor in the uRWELL basis:')
    rows = []
    for det in cfg.mappable_detectors():
        ct = cfg.channel_table(det)
        p2 = E._clusters_from_files(allf, ct, cfg.CLUSTER_R)
        j = ev.merge(p2[['eventId', 'x', 'y', 'single']], on='eventId')
        j = j[j['single']]
        A, t, r, keep = affine_window(j[['xf', 'yf']].to_numpy(),
                                      j[['x', 'y']].to_numpy())
        R, S = stretch(A)
        ang = np.degrees(np.arctan2(R[1, 0], R[0, 0]))
        print(f'   {det.name:7s} z={det.z:4.0f}  rot {ang:+7.3f} deg   '
              f'S = [[{S[0,0]:.5f} {S[0,1]:+.5f}], [{S[1,0]:+.5f} {S[1,1]:.5f}]]')
        rows.append((det.z, S[0, 0], S[1, 1], S[0, 1]))

    z = np.array([r[0] for r in rows])
    print('\nfit  S_ii(z) = 1 + z / L_i  :')
    out = {}
    for i, ax in ((1, 'x'), (2, 'y')):
        v = np.array([r[i] for r in rows])
        a, b = np.polyfit(z, v, 1)          # v = a z + b
        L = 1.0 / a
        out[ax] = (a, b, L)
        print(f'   {ax}: d(scale)/dz = {a:+.3e} /mm   intercept {b:.5f}   '
              f'virtual source L = {L / 1000:+.1f} m upstream')
        pred = b + a * 1370.0
        print(f'      -> predicted front->back slope at dz=1370 mm: {pred:.5f}'
              f'   measured {fb[ax]:.5f}   diff {pred - fb[ax]:+.5f}')
    off = np.array([r[3] for r in rows])
    print(f'\n   off-diagonal of S: {off.round(5)}  (shear; should be ~0)')


if __name__ == '__main__':
    main()
