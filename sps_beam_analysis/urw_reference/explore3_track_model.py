#!/usr/bin/env python3
"""Third look: which uRWELL track model actually points best at the P2 planes?

explore2 found the P2 residual growing with z (3.78 / 4.20 / 4.80 mm at
z = 320 / 630 / 940).  Decomposing that as  sigma_pred(z)^2 = sf^2 (1-t)^2 +
sb^2 t^2  with t = z/1370 and a 12 mm-pad term of 3.46 mm gives sf ~ 1.3 mm and
sb ~ 4.8 mm: the BACK uRWELL would be four times worse than the front, and the
whole front-back spread (4.5-4.8 mm, previously read as beam divergence) would
be back-plane noise rather than a real angle.

If that is right, a front-only "parallel beam" model must point BETTER at
P2_MID and P2_OUT than the two-point track does.  That is the test here:
the same P2 clusters, three track models.

  front  : (x, y) = front cluster, slope 0            (parallel beam)
  back   : (x, y) = back cluster,  slope 0
  2point : (x, y) = front + t * (back_corrected - front)
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
PAD_SIGMA = 12.0 / np.sqrt(12)     # 3.46 mm, a uniformly-lit 12 mm pad


def fit_and_report(src, dst, label):
    A, t, r, keep = affine_window(src, dst)
    sx, sy = r[keep, 0].std(), r[keep, 1].std()
    tot = np.sqrt(0.5 * (sx ** 2 + sy ** 2))
    implied = np.sqrt(max(tot ** 2 - PAD_SIGMA ** 2, 0.0))
    u, s, vt = np.linalg.svd(A)
    rot = u @ vt
    print(f'   {label:8s} rms dx {sx:5.2f}  dy {sy:5.2f}  -> pointing '
          f'{implied:5.2f} mm   (rot {np.degrees(np.arctan2(rot[1,0], rot[0,0])):+7.3f} '
          f'deg, sing {s[0]:.4f}/{s[1]:.4f}, {keep.mean():.1%} in window)')
    return tot


def main():
    run_json = os.path.join(E.RUN_DIR, 'run_config.json')
    sub_dir = os.path.join(E.RUN_DIR, E.SUB)
    ev, dz = E.urw_tracks(run_json, sub_dir, N_CHUNK)
    for ax in ('x', 'y'):
        s, o = np.polyfit(ev[f'{ax}f'], ev[f'{ax}b'], 1)
        for _ in range(4):
            r = ev[f'{ax}b'] - (s * ev[f'{ax}f'] + o)
            m = np.abs(r) < 2.5 * r.std()
            s, o = np.polyfit(ev.loc[m, f'{ax}f'], ev.loc[m, f'{ax}b'], 1)
        ev[f'{ax}b_c'] = (ev[f'{ax}b'] - o) / s
    print(f'front-back spread: x {(ev.xb_c - ev.xf).std():.2f}  '
          f'y {(ev.yb_c - ev.yf).std():.2f} mm')

    cfg = sc.RunConfig(key='x', run='highstat_eff_1',
                       data_root='/local/home/banco/P2_data/TB_July2026_H4/runs')
    allf = p2io.hit_files(cfg.combined_hits_dir(E.SUB))[:N_CHUNK]

    for det in cfg.mappable_detectors():
        ct = cfg.channel_table(det)
        p2 = E._clusters_from_files(allf, ct, cfg.CLUSTER_R)
        j = ev.merge(p2[['eventId', 'x', 'y', 'single', 'n_clus']],
                     on='eventId', how='inner')
        j = j[j['single']].reset_index(drop=True)
        t = det.z / dz
        dst = j[['x', 'y']].to_numpy()
        print(f'\n=== {det.name}  z={det.z:.0f}  t={t:.3f}  n={len(j)}')
        fit_and_report(j[['xf', 'yf']].to_numpy(), dst, 'front')
        fit_and_report(j[['xb_c', 'yb_c']].to_numpy(), dst, 'back')
        fit_and_report(np.column_stack([
            j['xf'] + t * (j['xb_c'] - j['xf']),
            j['yf'] + t * (j['yb_c'] - j['yf'])]), dst, '2point')
        # best-fit lever arm: scan the interpolation fraction
        best = (1e9, None)
        for tt in np.linspace(-0.4, 1.2, 33):
            A, o, r, keep = affine_window(np.column_stack([
                j['xf'] + tt * (j['xb_c'] - j['xf']),
                j['yf'] + tt * (j['yb_c'] - j['yf'])]), dst)
            w = np.sqrt(0.5 * (r[keep, 0].var() + r[keep, 1].var()))
            if w < best[0]:
                best = (w, tt)
        print(f'   best interpolation fraction t = {best[1]:+.3f} '
              f'(survey {t:.3f}) with rms {best[0]:.2f} mm')


if __name__ == '__main__':
    main()
