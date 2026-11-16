#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ABSOLUTE single-channel efficiency WITHOUT a resolution model (counting / conservation argument).

The prediction error moves each muon by an unknown amount, but it does not create or remove muons. For a pad j:
the number of times it fired = eps_j x (muons that really crossed it) + random hits, and the muons that really crossed
it = (local muon density) x (pad area). The local density is counted from the predicted points within a radius R of
the pad centre: rho = n_R / (pi R^2). With k_R = number of those tags for which pad j fired,

    k_R pi R^2 / (n_R area_j) = eps_j + b' R^2      once R >> resolution (all the pad's own hits are inside R),

the R^2 term being hits that are not the pad's own muon (noise, other particles). eps = intercept of a straight line
in R^2 over R in [r_lo, r_hi] (r_lo >= ~3 sigma). Only assumption: the muon density is smooth over R (tags near the
edge of the live area / masked pads bias it: pads with a dead pad within R are dropped, --clean).

Pooled per run and chamber, and per (chamber, threshold - baseline bin) over all runs; off-time tags subtracted;
optionally only scintillator-confirmed tags (--scint LO:HI ns).

Usage: eff_integral.py TAGS_DIR MAP_DIR OUT_DIR [--scint 140:310] [--rlo 20 --rhi 40]
"""
import argparse
import glob
import json
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

import eff_unfold as EU  # noqa: E402
from event_display import style, INK, MUTED, COL, MRK  # noqa: E402

RADII = np.arange(6, 41, 2.0)


def per_pad_counts(tags, pads, radii):
    """n[j, r], k[j, r]: off-time-subtracted tags predicted within radius r of pad j, and those with pad j fired."""
    tree = cKDTree(np.c_[tags['px'], tags['py']])
    w = np.where(tags['on'], 1.0, -float(tags['r_scale']))
    key = lambda x, y: np.round(x * 10).astype(np.int64) * 100000 + np.round(y * 10).astype(np.int64)  # noqa: E731
    fired = {}
    for i, kk in zip(tags['fired_i'].tolist(), key(tags['fired_x'], tags['fired_y']).tolist()):
        fired.setdefault(kk, set()).add(i)
    pk = key(pads.x.to_numpy(), pads.y.to_numpy())
    n = np.zeros((len(pads), len(radii)))
    k = np.zeros_like(n)
    v = np.zeros_like(n)
    for j, (x, y) in enumerate(zip(pads.x, pads.y)):
        idx = np.asarray(tree.query_ball_point([x, y], radii[-1]), int)
        if not len(idx):
            continue
        d = np.hypot(tags['px'][idx] - x, tags['py'][idx] - y)
        f = np.isin(idx, list(fired.get(pk[j], ())))
        for ir, r in enumerate(radii):
            s = d < r
            n[j, ir] = w[idx[s]].sum()
            k[j, ir] = (w[idx[s]] * f[s]).sum()
            v[j, ir] = (w[idx[s]] ** 2 * f[s]).sum()
    return n, k, v


def line_fit(R, y, e, rlo, rhi):
    """eps + b R^2 over [rlo, rhi]: intercept, its error (from the scatter of the cumulative points), slope."""
    s = (R >= rlo) & (R <= rhi) & np.isfinite(y)
    if s.sum() < 3:
        return np.nan, np.nan, np.nan
    A = np.c_[np.ones(s.sum()), R[s] ** 2]
    coef, *_ = np.linalg.lstsq(A, y[s], rcond=None)
    # the cumulative points are correlated: quote the error of the point at rlo (conservative)
    return coef[0], e[s][0], coef[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('tags_dir')
    ap.add_argument('map_dir')
    ap.add_argument('out_dir')
    ap.add_argument('--scint', default='')
    ap.add_argument('--rlo', type=float, default=20)
    ap.add_argument('--rhi', type=float, default=40)
    ap.add_argument('--clean', type=float, default=40, help='drop pads with a dead/masked/unmapped-in-area pad within '
                    'this distance (0 = keep all)')
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    X, Y, G = EU.pad_geometry()
    mp = json.load(open(f'{a.map_dir}/map_Km18_i5_unmasked.json'))
    cal = json.load(open(f'{a.map_dir}/sdtcal.json'))['fit']
    rsdt = json.load(open(f'{a.map_dir}/run_sdt.json'))
    area = np.full((32, 64), np.nan)
    pcsv = pd.read_csv(EU.VS.MAP_CSV)
    tree = cKDTree(np.c_[pcsv.pad_cx, pcsv.pad_cy])
    ok = np.isfinite(X)
    _, ii = tree.query(np.c_[X[ok], Y[ok]])
    area[ok] = pcsv.pad_area.to_numpy()[ii]

    rows = []
    for fz in sorted(glob.glob(f'{a.tags_dir}/tags_run*.npz')):
        run = int(re.search(r'run(\d+)', fz)[1])
        E = json.load(open(f'{a.tags_dir}/eff_run{run}.json'))
        tag = E['title'].split(' ', 1)[1]
        K = int(tag.split('_adwk2_K')[1].split('_')[0].replace('p', '+').replace('m', '-'))
        sdt = rsdt[f'p2b-config-cern-ext_gain4.5_peaktime200_opt_cosmic_{tag}.txt']['sdt']
        Z = np.load(fz)
        dead = {tuple(map(int, s.split(':'))) for lst in E['dead_or_masked'].values() for s in lst}
        for p in (3, 2, 1):
            tags = {k.split('_', 1)[1]: Z[k] for k in Z.files if k.startswith(f'p{p}_')}
            if a.scint:
                lo, hi = map(float, a.scint.split(':'))
                tags = EU.select_tags(tags, (tags['dt_scint'] >= lo) & (tags['dt_scint'] <= hi))
            allp = [(v, c) for v in EU.DETV[p] for c in range(64) if np.isfinite(X[v, c])]
            live = [(v, c) for v, c in allp if (v, c) not in dead]
            deadxy = np.array([(X[v, c], Y[v, c]) for v, c in allp if (v, c) in dead]).reshape(-1, 2)
            pads = pd.DataFrame([(v, c, X[v, c], Y[v, c], area[v, c],
                                  mp['threshold'][str(v)][str(c)] + (sdt[str(v)] - mp['sdt'][str(v)]) *
                                  cal[str(v)]['mV_per_sdt'] - mp['baseline'][str(v)][str(c)]) for v, c in live],
                                columns=['vmm', 'ch', 'x', 'y', 'area', 'dist'])
            if a.clean and len(deadxy):
                dd = cKDTree(deadxy).query(np.c_[pads.x, pads.y])[0]
                pads['clean'] = dd > a.clean
            else:
                pads['clean'] = True
            n, k, v = per_pad_counts(tags, pads, RADII)
            for j, r in pads.iterrows():
                rows.append(dict(run=run, K=K, det=EU.NAME[p], vmm=r.vmm, ch=r.ch, dist=r.dist, area=r.area,
                                 clean=r.clean, **{f'n{int(R)}': n[j, i] for i, R in enumerate(RADII)},
                                 **{f'k{int(R)}': k[j, i] for i, R in enumerate(RADII)},
                                 **{f'v{int(R)}': v[j, i] for i, R in enumerate(RADII)}))
        print(f'run_{run} K={K:+d} done', flush=True)
    T = pd.DataFrame(rows)
    T.to_csv(f'{a.out_dir}/pads.csv', index=False)

    def pooled(d):
        # sum over pads: k_R pi R^2 / (n_R area) with n_R, k_R summed using per-pad weights (area-normalised)
        out = {}
        y, e = [], []
        for R in RADII:
            nn = (d[f'n{int(R)}'] * d.area / (np.pi * R ** 2)).sum()   # expected crossings (all pads summed)
            kk = d[f'k{int(R)}'].sum()
            vv = d[f'v{int(R)}'].sum()
            y.append(kk / nn if nn > 0 else np.nan)
            e.append(np.sqrt(max(vv, 1)) / nn if nn > 0 else np.nan)
        y, e = np.array(y), np.array(e)
        eps, err, b = line_fit(RADII, y, e, a.rlo, a.rhi)
        out.update(eps=eps, err=err, slope=b, pads=len(d), **{f'y{int(R)}': y[i] for i, R in enumerate(RADII)})
        return pd.Series(out)

    use = T[T.clean & (T.dist >= 0)]
    blind = T[T.clean & (T.dist < 0)]
    perrun = use.groupby(['run', 'K', 'det']).apply(pooled, include_groups=False).reset_index()
    perrun.to_csv(f'{a.out_dir}/per_run.csv', index=False)
    bins = np.r_[0, 14, 18, 22, 26, 30, 34, 40, 46, 52, 58, 64, 70, 100]
    use = use.assign(dbin=pd.cut(use.dist, bins))
    perd = use.groupby(['det', 'dbin'], observed=True).apply(pooled, include_groups=False).reset_index()
    perd['dist'] = use.groupby(['det', 'dbin'], observed=True).dist.mean().to_numpy()
    perd.astype({'dbin': str}).to_csv(f'{a.out_dir}/per_distance.csv', index=False)
    pd.set_option('display.width', 220)
    print(f'pads used (clean, normal) per run/chamber: median {perrun.pads.median():.0f}')
    print('per run: eps (intercept) per chamber')
    print(perrun.pivot(index=['K', 'run'], columns='det', values='eps').round(3))
    print('ratio y(R) / y(20) at R = 10, 14, 20, 30, 40 (stability in R), runs at K <= 12 pooled:')
    pl = use[use.K <= 12]
    for det, d in pl.groupby('det'):
        s = pooled(d)
        print(f'  {det}: ' + ' '.join(f'R{R}:{s[f"y{R}"]:.3f}' for R in (8, 10, 14, 20, 26, 30, 36, 40)) +
              f'  -> eps {s.eps:.3f} +- {s.err:.3f}, slope {s.slope:.2e}/mm2')
    if len(blind):
        sb = pooled(blind)
        print(f'blind pads (baseline above threshold): eps {sb.eps:.3f} +- {sb.err:.3f} ({len(blind)} pad-runs)')
    print(perd.pivot(index='dbin', columns='det', values='eps').round(3))
    print(perd.pivot(index='dbin', columns='det', values='err').round(3))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for det, d in pl.groupby('det'):
        p = {'P2_IN': 3, 'P2_MID': 2, 'P2_OUT': 1}[det]
        s = pooled(d)
        ax[0].plot(RADII ** 2, [s[f'y{int(R)}'] for R in RADII], marker=MRK[p], color=COL[p], ms=4, lw=1.5, label=det)
        ax[0].plot([0, RADII[-1] ** 2], [s.eps, s.eps + s.slope * RADII[-1] ** 2], color=COL[p], lw=0.8, ls='--')
    ax[0].axvline(a.rlo ** 2, color=MUTED, lw=0.8, ls=':')
    ax[0].set_xlabel('R$^2$ (mm$^2$): tags counted within R of the pad centre', fontsize=8, color=INK)
    ax[0].set_ylabel('fired / expected crossings', fontsize=8, color=INK)
    ax[0].set_title('counting estimator vs R (runs K <= 12); dashed: eps + b R$^2$', fontsize=9, color=INK, loc='left')
    ax[0].legend(frameon=False, fontsize=7)
    for det, d in perd.groupby('det'):
        p = {'P2_IN': 3, 'P2_MID': 2, 'P2_OUT': 1}[det]
        ax[1].errorbar(d.dist, d.eps, d.err, marker=MRK[p], color=COL[p], ms=4, lw=1.5, label=det)
    ax[1].set_ylim(0, 1.1)
    ax[1].set_xlabel('threshold - baseline of the channel (mV)', fontsize=8, color=INK)
    ax[1].set_ylabel('single-channel efficiency (absolute)', fontsize=8, color=INK)
    ax[1].set_title('absolute efficiency, no resolution model' + (' (scintillator-confirmed tags)' if a.scint else ''),
                    fontsize=9, color=INK, loc='left')
    ax[1].legend(frameon=False, fontsize=7)
    for x in ax:
        style(x)
    fig.tight_layout()
    fig.savefig(f'{a.out_dir}/eff_integral.png', dpi=150)


if __name__ == '__main__':
    main()
