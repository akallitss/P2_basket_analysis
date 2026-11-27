#!/usr/bin/env python3
"""38_nohit_vs_angle_pillar_distance.py -- is the angle dependence of the no-hit
rate a pillar shadow? Model-free check: no-hit vs zenith angle for all interior
rays and for rays >= FAR mm from the nearest small pillar (and within NEAR mm).
Needs 34 (corrected alignment, rays_aligned.pkl). A pillar shadow would make the
trend vanish far from pillars.
Usage: python3 38_nohit_vs_angle_pillar_distance.py [run_key]
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_study_cache as sc, p2_study_common as cm, p2_mapping as pmap  # noqa: E401,E402

FAR, NEAR = 1.4, 0.5


def main():
    run_key = sys.argv[1] if len(sys.argv) > 1 else 'det2_long_repos'
    hits, rays, g, cfg = sc.load(run_key)
    d34 = cfg.out_dir('34_clustering_matching'); out = cfg.out_dir('38_nohit_vs_angle')
    A = np.load(f'{d34}/alignment.npz')
    T = cm.make_T(g, R=A['R'], s=A['s'], pad_mean=A['pad_mean'], m3_mean=A['m3_mean'])
    ra = pd.read_pickle(f'{d34}/rays_aligned.pkl')
    pil = pmap.load_pillars(cfg.MASK_GBR_PATH)
    px, py = T.apply(pil.x.to_numpy(), pil.y.to_numpy())
    big = np.c_[px, py][pil.big.to_numpy()]; sm = np.c_[px, py][~pil.big.to_numpy()]
    xy = ra[['x', 'y']].to_numpy()
    dpil, _ = cKDTree(sm).query(xy); dbig, _ = cKDTree(big).query(xy)
    nohit = ~ra.eventId.isin(hits.eventId).to_numpy()
    inter = (ra.sd.to_numpy() >= 4) & (dbig > 10)
    th = ra.theta.to_numpy()
    e = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 23, 26])
    rows = []
    for nm, sel in [('all interior', inter), (f'>= {FAR} mm from any small pillar', inter & (dpil >= FAR)),
                    (f'< {NEAR} mm from a small pillar', inter & (dpil < NEAR))]:
        for lo, hi in zip(e[:-1], e[1:]):
            k = sel & (th >= lo) & (th < hi)
            if k.sum() >= 30:
                p = nohit[k].mean()
                rows.append(dict(sel=nm, th=0.5 * (lo + hi), n=int(k.sum()), nohit=100 * p,
                                 err=100 * np.sqrt(max(p * (1 - p), 1e-4) / k.sum())))
    t = pd.DataFrame(rows); t.to_csv(f'{out}/nohit_vs_angle.csv', index=False)
    fig, ax = plt.subplots(figsize=(8, 4.6))
    for (nm, col, mk, dx) in [('all interior', '#257de0', 'o', 0), (f'>= {FAR} mm from any small pillar', '#5c5c5c', 's', .3),
                              (f'< {NEAR} mm from a small pillar', '#ec4647', 'D', -.3)]:
        s = t[t.sel == nm]
        ax.errorbar(s.th + dx, s.nohit, s.err, fmt=mk + '-', color=col, lw=2, ms=5,
                    label=f'{nm} (n={s.n.sum():,})')
    ax.set(xlabel='track zenith angle [deg]', ylabel='no P2 hit at all [%]', ylim=(0, None),
           title='No-hit rate vs angle: the trend is the same far from the pillars')
    ax.grid(color='#e6e6e6'); ax.legend(frameon=False, fontsize=8.5)
    for s_ in ('top', 'right'):
        ax.spines[s_].set_visible(False)
    fig.tight_layout(); fig.savefig(f'{out}/nohit_vs_angle_pillar_distance.png', dpi=150); plt.close(fig)
    print(t.pivot(index='th', columns='sel', values='nohit').round(2).to_string())


if __name__ == '__main__':
    main()
