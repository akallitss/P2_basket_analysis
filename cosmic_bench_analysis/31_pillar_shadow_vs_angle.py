#!/usr/bin/env python3
"""31_pillar_shadow_vs_angle.py -- are near-vertical tracks lost on the small
insulation pillars?

Usage:  python3 31_pillar_shadow_vs_angle.py <run_key> [--interior 6]

Physics: the ionisation column of a track crossing the drift gap g projects
onto the mesh as a segment of length L = g * tan(theta). If that whole segment
lies on a pillar (disc radius r on an a x a grid) every electron is lost -> no
hit. Probability (segment midpoint uniform in a cell):
    P(theta) = pi * (r^2 - L^2/4) / a^2      for L < 2r, else 0
The pointing resolution (~3 mm) is far larger than a pillar, so a distance-to-
pillar test cannot see this; the ANGLE dependence can. We fit
    nohit(theta) = b + P(theta; r_eff)
to interior rays (far from dead pads / edges), b = angle-independent rest.
Output: 31_pillar_shadow/nohit_vs_theta<sfx>.png + .csv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import awkward as ak
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa      # noqa: E402
import p2_mapping as pmap      # noqa: E402
import p2_align as pa          # noqa: E402,F401  (puts M3RefTracking on the path)
from M3RefTracking import M3RefTracking  # noqa: E402
from scipy.spatial import cKDTree      # noqa: E402

INK, INK2, GRID = '#1f1f1f', '#5c5c5c', '#e6e6e6'
C_DATA, C_FIT, C_GEO = '#ec4647', '#257de0', '#9c9c9c'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument('--interior', type=float, default=6.0,
                    help='use rays at least this far [mm] from dead pads / edge (needs stage 30)')
    args = ap.parse_args()
    cfg = qa.get_config(args.run_key)
    out_dir = cfg.out_dir('31_pillar_shadow')
    sfx = cfg.product_suffix(args.veto_sparks)
    g = cfg.DRIFT_GAP_MM

    pil = pmap.load_pillars(cfg.MASK_GBR_PATH)
    sm = pil[~pil['big']]
    r_geo = float(sm['r'].median())
    a = float(np.median(cKDTree(sm[['x', 'y']].to_numpy()).query(sm[['x', 'y']].to_numpy(), k=2)[0][:, 1]))

    rays = pd.read_csv(f'{cfg.out_dir("06_efficiency")}/ray_hit_miss_list{sfx}.csv')
    for c in ('in_active', 'has_any', 'within'):
        rays[c] = rays[c].astype(str).str.lower().isin(('true', '1'))
    ra = rays[rays['in_active']].copy()
    # interior selection from stage 30's distance map if available (recompute cheaply otherwise: all)
    m30 = os.path.join(cfg.out_dir('30_miss_anatomy'), f'miss_vs_boundary{sfx}.csv')
    use_interior = os.path.isfile(m30)
    m3 = M3RefTracking(os.path.join(cfg.m3_tracking_dir, ''), single_track=True,
                       chi2_cut=qa.M3_CHI2_CUT, min_nclus=qa.M3_MIN_NCLUS)
    rd = m3.ray_data

    def col(n):
        return np.asarray(ak.to_numpy(ak.ravel(rd[n])), dtype=float)
    ang = pd.DataFrame({'eventId': np.asarray(ak.to_numpy(ak.ravel(rd['evn'])), dtype=np.int64),
                        'theta': np.degrees(np.arctan2(np.hypot(col('X_Up') - col('X_Down'),
                                                                col('Y_Up') - col('Y_Down')),
                                                       np.abs(col('Z_Up') - col('Z_Down'))))})
    ra = ra.merge(ang.drop_duplicates('eventId'), on='eventId', how='inner')
    if use_interior:
        # stage 30 stores per-ray distances only implicitly; recompute the interior flag the
        # same way: rays whose nearest LIVE pad is far from any dropped pad (pad-centre proxy)
        ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
                                      det_name=cfg.DET_NAME, drop_connectors=cfg.DEAD_CONNECTORS,
                                      drop_channels=cfg.DEAD_CHANNELS,
                                      strategy_overrides=cfg.STRATEGY_OVERRIDES)
        padc = ct[ct['mapped']].drop_duplicates('channel_id').reset_index(drop=True)
        fp = pd.read_csv(f'{cfg.out_dir("06_efficiency")}/pad_footprint{sfx}.csv')
        P = np.c_[padc.pad_cx, padc.pad_cy, np.ones(len(padc))]
        coef, *_ = np.linalg.lstsq(P, fp[['x', 'y']].to_numpy(), rcond=None)
        dead = ct.attrs['dead_pad_xy']
        dq = np.c_[dead, np.ones(len(dead))] @ coef
        dd, _ = cKDTree(dq).query(ra[['x', 'y']].to_numpy())
        half_pitch = 5.9
        ra = ra[dd - half_pitch >= args.interior]
    n_all = len(ra)
    edges = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 23, 26, 30])
    ra['tb'] = pd.cut(ra['theta'], edges)
    t = ra.groupby('tb', observed=False).agg(n=('has_any', 'size'), nohit=('has_any', lambda s: 1 - s.mean()),
                                             theta=('theta', 'mean'),
                                             miss=('within', lambda s: 1 - s.mean())).reset_index()
    t = t[t.n >= 50].copy()
    t['err'] = np.sqrt(t.nohit * (1 - t.nohit) / t.n)
    t['tb'] = t['tb'].astype(str)

    def model(theta_deg, r, b):
        L = g * np.tan(np.radians(theta_deg))
        return b + np.where(L < 2 * r, np.pi * np.clip(r ** 2 - L ** 2 / 4, 0, None) / a ** 2, 0.0)

    # integrate the model over each bin's actual theta distribution
    def model_binned(_, r, b):
        out = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            th = ra.theta[(ra.theta > lo) & (ra.theta <= hi)].to_numpy()
            if len(th) >= 50:
                out.append(model(th, r, b).mean())
        return np.array(out)
    popt, pcov = curve_fit(model_binned, None, t.nohit.to_numpy(), p0=(0.5, 0.01),
                           sigma=t.err.to_numpy(), absolute_sigma=True, bounds=([0.05, 0], [2.0, 0.2]))
    r_eff, b = popt
    er, eb = np.sqrt(np.diag(pcov))
    t['fit'] = model_binned(None, r_eff, b)
    t['geo'] = model_binned(None, r_geo, b)
    chi2 = (((t.nohit - t.fit) / t.err) ** 2).sum()
    tot_shadow = (model(ra.theta.to_numpy(), r_eff, 0)).mean()
    tot_geo = (model(ra.theta.to_numpy(), r_geo, 0)).mean()
    t.to_csv(f'{out_dir}/nohit_vs_theta{sfx}.csv', index=False)
    print(t[['tb', 'n', 'theta', 'nohit', 'err', 'fit', 'miss']].round(4).to_string(index=False))
    print(f'\nfit: r_eff = {r_eff:.3f} +- {er:.3f} mm (Gerber pillar r = {r_geo:.2f} mm, grid {a:.1f} mm, gap {g:g} mm), '
          f'b = {100*b:.2f} +- {100*eb:.2f} %, chi2/ndf = {chi2:.1f}/{len(t)-2}')
    print(f'pillar shadow averaged over this run\'s angles: {100*tot_shadow:.2f}% of rays '
          f'(geometric r: {100*tot_geo:.2f}%); angle-independent rest b = {100*b:.2f}%  (n={n_all})')

    plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 5.2), gridspec_kw={'width_ratios': [1.6, 1]})
    ax.errorbar(t.theta, 100 * t.nohit, 100 * t.err, fmt='o', color=C_DATA, ms=6, lw=1.5,
                label=f'data: interior no-hit (n={n_all:,})')
    th = np.linspace(0, 30, 300)
    ax.plot(th, 100 * model(th, r_eff, b), color=C_FIT, lw=2.2,
            label=f'shadow model fit: r_eff = {r_eff:.2f} ± {er:.2f} mm, b = {100*b:.2f}%')
    ax.plot(th, 100 * model(th, r_geo, b), color=C_GEO, lw=1.6, ls='--',
            label=f'geometric pillar r = {r_geo:.2f} mm (same b)')
    ax.axhline(100 * b, color=INK2, lw=0.8, ls=':')
    ax.set_xlabel('M3 track zenith angle θ [deg]'); ax.set_ylabel('no hit at all [%]')
    ax.set_title('Near-vertical tracks lose all their electrons on a pillar', loc='left',
                 fontweight='bold', color=INK)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(True, color=GRID); ax.set_axisbelow(True)
    ax.text(0.99, 0.62, f'projected ionisation streak L = {g:g} mm × tanθ\n'
            f'lost if the whole streak lies on a pillar\n(disc r on a {a:g} mm grid)',
            transform=ax.transAxes, ha='right', fontsize=9, color=INK2)
    ax2.hist(ra.theta, bins=60, color='#9c9c9c')
    ax2.set_xlabel('M3 track zenith angle θ [deg]'); ax2.set_ylabel('tracks')
    ax2.set_title('track-angle distribution (M3 acceptance)', loc='left', color=INK)
    for a_ in (ax, ax2):
        for s_ in ('top', 'right'):
            a_.spines[s_].set_visible(False)
    fig.suptitle(f'{cfg.DET_NAME} pillar shadowing vs track angle — {cfg.RUN}/{cfg.SUB_RUN}  '
                 f'(run-averaged shadow {100*tot_shadow:.2f}%)', x=0.01, ha='left', color=INK, fontsize=12)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/nohit_vs_theta{sfx}.png', dpi=150)
    plt.close(fig)
    print(f'  saved {out_dir}/nohit_vs_theta{sfx}.png')


if __name__ == '__main__':
    main()
