#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Chamber-level efficiency with a fiducial cut on the MUONS (19 Nov 2026).

Why: a fiducial on the pads (channel_efficiency.py measure --fiducial) is biased - the counting estimator takes the
expected crossings of a pad from the density of the PREDICTED positions around it, smeared by the tracking, so it
over-estimates central pads and under-estimates edge pads; only the whole chamber balances (eps > 1 at 15-20 mm).
Here the cut is on the prediction of each muon, and the efficiency is counted per muon, with no density model:

  A(R) = P(any pad of the chamber fires within R mm of the prediction) in the muon window
         - the same in the sideband window (random hits),
  over scintillator-confirmed tags (dt_scint in --scint), on-time minus r_scale x off-time (as checks.py section 3),
  whose prediction lies within 6 mm of a clean live pad (no dead/masked pad within 40 mm), and
  d_edge >= N mm, d_edge = distance of the prediction to the outline (convex hull) of the chamber's pad centres.

Outputs (OUT_DIR): fiducial_cut.csv  per configuration, chamber, N: tags, A10/A20/A30, A20_err
                   fiducial_profile.csv  the same in bands of d_edge (where the losses are)
                   fiducial_muons.png  A20 vs N per configuration (left column) and A20 vs d_edge band (right)
Usage: fiducial_muons.py TAGS_DIR PSEUDO_DIR OUT_DIR [--cuts 0 5 10 15 20 30] [--scint 140:310]
       (PSEUDO_DIR/chunks.json maps the pseudo-run ids to the configurations; label() of the pipeline's common.py
       when it is importable from --scripts, else the cfg name)"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import channel_efficiency as CE  # noqa: E402

RADII = (10, 20, 30)
BANDS = [0, 5, 10, 15, 20, 30, 40, 60, 1000]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('tags_dir'); ap.add_argument('pseudo'); ap.add_argument('out')
    ap.add_argument('--cuts', type=float, nargs='+', default=[0, 5, 10, 15, 20, 30])
    ap.add_argument('--scint', default='140:310')
    ap.add_argument('--scripts', default='', help='pipeline scripts dir with common.py (configuration labels)')
    a = ap.parse_args()
    from scipy.spatial import cKDTree
    os.makedirs(a.out, exist_ok=True)
    lo, hi = map(float, a.scint.split(':'))
    chunks = {int(k): v for k, v in json.load(open(f'{a.pseudo}/chunks.json')).items()}
    label = lambda c: c  # noqa: E731
    if a.scripts:
        sys.path.insert(0, a.scripts)
        from common import label  # noqa: F811
    X, Y, _ = CE.pad_geometry()
    allxy = {q: np.array([(X[v, c], Y[v, c]) for v in CE.DETV[q] for c in range(64) if np.isfinite(X[v, c])])
             for q in (3, 2, 1)}
    recs = []                     # one row per (chunk, chamber, tag) kept: cfg, det, w, d_edge, fired/side per R
    for fz in sorted(glob.glob(f'{a.tags_dir}/tags_run*.npz')):
        p = int(re.search(r'run(\d+)', fz)[1])
        if p not in chunks:
            continue
        info = CE.run_info(a.tags_dir, p)
        cfg = label(chunks[p]['cfg'])
        for q in (3, 2, 1):
            T = CE.load_tags(fz, q)
            allp = [(v, c) for v in CE.DETV[q] for c in range(64) if np.isfinite(X[v, c])]
            live = np.array([(X[v, c], Y[v, c]) for v, c in allp if (v, c) not in info['dead']])
            dead = np.array([(X[v, c], Y[v, c]) for v, c in allp if (v, c) in info['dead']]).reshape(-1, 2)
            clean = cKDTree(dead).query(live)[0] > 40 if len(dead) else np.ones(len(live), bool)
            dt = T['dt_scint']
            sel = (dt >= lo) & (dt <= hi)
            P = np.c_[T['px'], T['py']]
            dd, ii = cKDTree(live).query(P)
            sel &= (dd < 6) & clean[ii]
            if not sel.any():
                continue
            w = np.where(T['on'], 1.0, -float(T['r_scale']))
            got = {}
            for pre in ('fired', 'side'):
                ti = T[f'{pre}_i']
                r = np.hypot(T[f'{pre}_x'] - T['px'][ti], T[f'{pre}_y'] - T['py'][ti])
                for R in RADII:
                    f = np.zeros(len(w), bool)
                    f[ti[r < R]] = True
                    got[(pre, R)] = f
            idx = np.nonzero(sel)[0]
            de = CE.hull_edge_distance(allxy[q], P[idx])
            d = dict(cfg=np.full(len(idx), cfg, object), det=np.full(len(idx), CE.NAME[q], object),
                     w=w[idx], d_edge=de)
            for R in RADII:
                d[f'h{R}'] = got[('fired', R)][idx].astype(float) - got[('side', R)][idx]
            recs.append(pd.DataFrame(d))
        print(f'run{p} done', flush=True)
    if not recs:
        raise SystemExit('no tags')
    D = pd.concat(recs, ignore_index=True)

    def agg(g):
        W = g.w.sum()
        out = dict(tags=len(g), muons=W)
        for R in RADII:
            out[f'A{R}'] = (g.w * g[f'h{R}']).sum() / W if W > 0 else np.nan
        A = out['A20']
        out['A20_err'] = np.sqrt(max(A * (1 - A), 1e-4) / W) if W > 0 else np.nan
        return pd.Series(out)

    cut = []
    for N in a.cuts:
        g = D[D.d_edge >= N].groupby(['cfg', 'det']).apply(agg, include_groups=False).reset_index()
        g.insert(2, 'cut_mm', N)
        cut.append(g)
    C = pd.concat(cut, ignore_index=True)
    C.to_csv(f'{a.out}/fiducial_cut.csv', index=False)
    BL = [f'{l:g}-{h:g}' if h < 1000 else f'>={l:g}' for l, h in zip(BANDS[:-1], BANDS[1:])]
    D['band'] = pd.Categorical(pd.cut(D.d_edge, BANDS, right=False, labels=BL), BL, ordered=True)
    B = D.groupby(['cfg', 'det', 'band'], observed=True).apply(agg, include_groups=False).reset_index()
    B.to_csv(f'{a.out}/fiducial_profile.csv', index=False)
    Ball = D.groupby(['det', 'band'], observed=True).apply(agg, include_groups=False).reset_index()

    pd.set_option('display.width', 220)
    print('\n== A20 (any pad within 20 mm), muons predicted >= N mm inside the pad outline; all configurations together')
    rows = []
    for N in a.cuts:
        for det, g in D[D.d_edge >= N].groupby('det'):
            s = agg(g); rows.append(dict(det=det, cut_mm=N, muons=round(s.muons), A20=s.A20, A30=s.A30, err=s.A20_err))
    print(pd.DataFrame(rows).pivot(index='cut_mm', columns='det', values=['A20', 'muons']).round(3).to_string())
    print('\n== A20 per band of distance to the outline (mm), all configurations together')
    print(Ball.pivot(index='band', columns='det', values='A20').round(3).to_string())
    print('\n== A20 per configuration and cut')
    print(C.pivot_table(index=['cfg', 'det'], columns='cut_mm', values='A20').round(3).to_string())

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    dets = ['P2_IN', 'P2_MID', 'P2_OUT']
    fig, ax = plt.subplots(3, 2, figsize=(15, 13))
    for r, det in enumerate(dets):
        for cfg, g in C[C.det == det].groupby('cfg'):
            ax[r, 0].errorbar(g.cut_mm, g.A20, g.A20_err, marker='o', ms=3, lw=1, capsize=2, label=cfg)
        ax[r, 0].set_title(f'{det}: chamber efficiency A20, muons >= N mm inside the pad outline')
        ax[r, 0].set_xlabel('fiducial cut N (mm)'); ax[r, 0].set_ylabel('A20')
        b = Ball[Ball.det == det]
        b = b.set_index('band').reindex(BL).dropna(subset=['A20'])
        order = list(b.index)
        ax[r, 1].errorbar(range(len(b)), b.A20, b.A20_err, marker='s', color='k')
        ax[r, 1].set_xticks(range(len(b))); ax[r, 1].set_xticklabels(order, rotation=30, fontsize=8)
        ax[r, 1].set_title(f'{det}: A20 per band of distance to the outline (all configurations)')
        ax[r, 1].set_xlabel('distance of the predicted point to the pad outline (mm)')
        for x_ in ax[r]:
            x_.set_ylim(0.4, 1.02); x_.grid(alpha=0.3)
    ax[0, 0].legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(f'{a.out}/fiducial_muons.png', dpi=90)
    print(f'\nwrote {a.out}/fiducial_cut.csv, fiducial_profile.csv, fiducial_muons.png')


if __name__ == '__main__':
    main()
