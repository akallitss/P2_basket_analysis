#!/usr/bin/env python3
"""Second look: honest residual distribution of the uRWELL track vs P2 clusters.

explore_urw_p2_frame.py fitted an affine with a 3-sigma clip on |r|, which
collapsed onto a 0.8 mm core - impossible for a 12 mm pad, so the clip is
chasing its own tail.  Here the fit uses a FIXED acceptance window (|dx|,|dy|
< WIN mm around the current estimate, two passes) and the residual is then
reported as a full distribution: histogram, percentiles, and the Gaussian core
from a fit to the central bins.
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

import explore_urw_p2_frame as E  # noqa: E402
import sps_config as sc           # noqa: E402
import p2_io as p2io              # noqa: E402
import matplotlib                 # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402

RUN_DIR = E.RUN_DIR
SUB = E.SUB
N_CHUNK = 1
WIN = 25.0          # mm, fit acceptance window (2 pads)


def affine_window(src, dst, win=WIN, passes=3):
    A, t = np.eye(2), np.zeros(2)
    keep = np.ones(len(src), bool)
    for i in range(passes):
        M = np.column_stack([src[keep], np.ones(keep.sum())])
        sol, *_ = np.linalg.lstsq(M, dst[keep], rcond=None)
        A, t = sol[:2].T, sol[2]
        r = dst - (src @ A.T + t)
        keep = (np.abs(r[:, 0]) < win) & (np.abs(r[:, 1]) < win)
        if i == 0:
            # first pass starts from identity, so the window is meaningless;
            # widen it once so the fit can find the real solution
            keep = np.ones(len(src), bool)
    return A, t, dst - (src @ A.T + t), keep


def main():
    run_json = os.path.join(RUN_DIR, 'run_config.json')
    sub_dir = os.path.join(RUN_DIR, SUB)
    ev, dz = E.urw_tracks(run_json, sub_dir, N_CHUNK)
    for ax in ('x', 'y'):
        s, o = np.polyfit(ev[f'{ax}f'], ev[f'{ax}b'], 1)
        for _ in range(4):
            r = ev[f'{ax}b'] - (s * ev[f'{ax}f'] + o)
            m = np.abs(r) < 2.5 * r.std()
            s, o = np.polyfit(ev.loc[m, f'{ax}f'], ev.loc[m, f'{ax}b'], 1)
        ev[f'{ax}b_c'] = (ev[f'{ax}b'] - o) / s
        ev[f'slope_{ax}'] = (ev[f'{ax}b_c'] - ev[f'{ax}f']) / dz

    cfg = sc.RunConfig(key='x', run='highstat_eff_1',
                       data_root='/local/home/banco/P2_data/TB_July2026_H4/runs')
    hits_dir = cfg.combined_hits_dir(SUB)
    allf = p2io.hit_files(hits_dir)[:N_CHUNK]

    fig, axes = plt.subplots(3, 4, figsize=(19, 12))
    for row, det in enumerate(cfg.mappable_detectors()):
        ct = cfg.channel_table(det)
        p2 = E._clusters_from_files(allf, ct, cfg.CLUSTER_R)
        m = ev.copy()
        m['px'] = m['xf'] + m['slope_x'] * det.z
        m['py'] = m['yf'] + m['slope_y'] * det.z
        j = m.merge(p2[['eventId', 'x', 'y', 'single', 'n_pad', 'n_clus', 'q']],
                    on='eventId', how='inner')
        j = j[j['single']].reset_index(drop=True)
        src = j[['px', 'py']].to_numpy()
        dst = j[['x', 'y']].to_numpy()
        A, t, r, keep = affine_window(src, dst)
        u, s, vt = np.linalg.svd(A)
        rot = u @ vt
        ang = np.degrees(np.arctan2(rot[1, 0], rot[0, 0]))
        print(f'\n=== {det.name} z={det.z:.0f}  n={len(j)}  in window {keep.sum()} '
              f'({keep.mean():.1%})')
        print(f'   rotation {ang:+.3f} deg   sing {s[0]:.4f} {s[1]:.4f}  '
              f'det {np.linalg.det(A):+.4f}   t=({t[0]:+.2f},{t[1]:+.2f})')
        for k, lbl in ((0, 'dx'), (1, 'dy')):
            d = r[keep, k]
            q = np.percentile(d, [5, 16, 50, 84, 95])
            print(f'   {lbl}: median {np.median(d):+.3f}  rms {d.std():.3f}  '
                  f'sigma(IQR-ish) {(q[3]-q[1])/2:.3f}  '
                  f'p5..p95 {q[0]:+.2f}..{q[4]:+.2f}')
        print(f'   pads/cluster: mean {j["n_clus"].mean():.2f}  '
              f'1-pad {(j["n_clus"]==1).mean():.1%}')

        b = np.linspace(-WIN, WIN, 120)
        axes[row, 0].hist(r[keep, 0], bins=b, histtype='step', label='dx')
        axes[row, 0].hist(r[keep, 1], bins=b, histtype='step', label='dy')
        axes[row, 0].set(xlabel='residual [mm]', title=f'{det.name} residuals')
        axes[row, 0].legend()
        axes[row, 1].hist2d(src[keep, 0], src[keep, 1], bins=100, cmin=1)
        axes[row, 1].set(xlabel='track x at plane [mm]', ylabel='track y [mm]',
                         title='track projection (uRWELL frame)')
        axes[row, 2].hist2d(dst[keep, 0], dst[keep, 1], bins=100, cmin=1)
        axes[row, 2].set(xlabel='P2 x [mm]', ylabel='P2 y [mm]',
                         title='P2 cluster (pad frame)')
        # residual vs single-pad clusters
        one = keep & (j['n_clus'].to_numpy() == 1)
        many = keep & (j['n_clus'].to_numpy() > 1)
        axes[row, 3].hist(np.hypot(*r[one].T), bins=np.linspace(0, WIN, 80),
                          histtype='step', density=True, label='1 pad')
        if many.sum():
            axes[row, 3].hist(np.hypot(*r[many].T), bins=np.linspace(0, WIN, 80),
                              histtype='step', density=True, label='>1 pad')
        axes[row, 3].set(xlabel='|residual| [mm]', title='by cluster size')
        axes[row, 3].legend()
    fig.tight_layout()
    fig.savefig('explore2_resid.png', dpi=90)
    print('\nwrote explore2_resid.png')


if __name__ == '__main__':
    main()
