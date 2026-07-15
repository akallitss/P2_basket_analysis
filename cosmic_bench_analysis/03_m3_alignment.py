#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_m3_alignment.py

Align the P2 pad detector to the M3 reference telescope by scanning the
rotation of the P2 pad plane about the beam (z) axis and looking for the angle
at which the per-event P2 hit position starts to correlate with the M3-predicted
track impact point.

Why a rotation scan
-------------------
The M3 track impact is reconstructed in the *bench* frame (origin ~ detector
centre, mm). The P2 pad positions come from the Gerber map in the *PCB* frame (a
fan sector sitting in the first quadrant). Between the two there is an unknown
rigid transform: a translation, a rotation about z, and possibly a reflection
(the connectors are mounted 'rotated_inverted'). We do not know the rotation a
priori, so we scan it.

The figure of merit is built from Pearson correlation coefficients, which are
invariant to translation and scale -- so we can find the rotation angle without
first solving the translation:

    rotate the pad centroid by theta:
        x' =  x_pad*cos - y_pad*sin
        y' =  x_pad*sin + y_pad*cos
    FOM(theta) = corr(x_m3, x') + corr(y_m3, y')

At the correct angle both correlations go strongly positive. Pure rotation spans
only SO(2); a reflection cannot be reached by rotating, so we run the whole scan
twice: once with the pads as-is (mirror=+1) and once with y_pad -> -y_pad
(mirror=-1), and report the best over both.

A Procrustes (best-fit rigid transform) is also computed as an independent
cross-check -- its implied rotation angle should coincide with the scan peak.

Products (written to <Analysis>/<detN>/<run>/<sub_run>/03_m3_alignment/):
  rotation_scan.png     FOM (and r_xx, r_yy) vs theta, for both mirror states
  best_angle_corr.png   x_m3 vs x' and y_m3 vs y' scatter at the best angle
  matched_positions.png M3 impact vs (rotated) P2 centroid density
  alignment_summary.txt  best angle, correlations, Procrustes transform, N

Usage
  python3 03_m3_alignment.py [run_key] [--z 232] [--dtheta 0.5]
                             [--min-amp 0] [--m3-fiducial 1500] [--leading-pad]
"""

import matplotlib
matplotlib.use('Agg')  # headless

import os
import sys
import glob
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot

import p2_qa_config as qa
import p2_mapping as pmap
import p2_align as pa
import p2_sparks as ps

# M3RefTracking lives in the nTof_x17 analysis repo (same reader m3_comparison uses).
_M3_PKG = os.path.expanduser(
    '~/Documents/PostDocSaclay/nTof_x17/cosmic_bench_analysis')
if _M3_PKG not in sys.path:
    sys.path.insert(0, _M3_PKG)
from M3RefTracking import M3RefTracking  # noqa: E402

_HIT_BRANCHES = ['eventId', 'trigger_timestamp_ns', 'channel', 'amplitude', 'feu']


# --------------------------------------------------------------------------- #
# Loaders
# --------------------------------------------------------------------------- #
def load_m3_positions(m3_dir, z, chi2_cut=qa.M3_CHI2_CUT,
                      min_nclus=qa.M3_MIN_NCLUS):
    """Single-track M3 impact positions at plane z, keyed by event number.

    Good-track recipe: Chi2X,Chi2Y < chi2_cut AND (tracking-v2 rays) NClusX,NClusY
    >= min_nclus (p2_qa_config.M3_CHI2_CUT / M3_MIN_NCLUS)."""
    # M3RefTracking builds paths as ray_dir + file_name, so it needs a trailing sep.
    m3 = M3RefTracking(os.path.join(m3_dir, ''), single_track=True,
                       chi2_cut=chi2_cut, min_nclus=min_nclus)
    x, y, evn = m3.get_xy_positions(z)
    df = pd.DataFrame({'eventId': np.asarray(evn, dtype=np.int64),
                       'x_m3': np.asarray(x, dtype=float),
                       'y_m3': np.asarray(y, dtype=float)})
    df = df[np.isfinite(df['x_m3']) & np.isfinite(df['y_m3'])]
    # a few evn repeat across files; keep the first
    df = df.drop_duplicates('eventId')
    print(f'M3: {len(df):,} single-track events projected to z={z}.')
    return df


def load_p2_centroids(hits_dir, channel_table, min_amp=0.0, leading_pad=False,
                      spark_veto=None):
    """Per-event charge-weighted P2 pad centroid (or leading-pad position)."""
    files = sorted(glob.glob(os.path.join(hits_dir, '*.root')))
    if not files:
        raise FileNotFoundError(f'No combined-hits ROOT files in {hits_dir}')
    feu_set = set(channel_table.attrs['feus'])
    parts = []
    for fp in files:
        arr = uproot.open(f'{fp}:hits').arrays(_HIT_BRANCHES, library='pd')
        parts.append(arr[arr['feu'].isin(feu_set)].copy())
    hits = pd.concat(parts, ignore_index=True)
    if spark_veto is not None:
        hits, n_rm = spark_veto.apply(hits)
        print(f'Spark veto: dropped {n_rm:,} hits in {len(spark_veto.sparks)} '
              f'sparks ({100*(1-spark_veto.live_fraction()):.2f}% deadtime) + '
              f'{spark_veto.last_burst_events} burst events '
              f'(>= {spark_veto.burst_npads} pads).')
    if min_amp > 0:
        hits = hits[hits['amplitude'] >= min_amp]

    hits = pmap.attach_pads_to_hits(hits, channel_table)
    hits = hits[hits['mapped'] & hits['pad_cx'].notna()]

    if leading_pad:
        idx = hits.groupby('eventId')['amplitude'].idxmax()
        cen = hits.loc[idx, ['eventId', 'pad_cx', 'pad_cy']].copy()
        cen = cen.rename(columns={'pad_cx': 'x_pad', 'pad_cy': 'y_pad'})
        cen['n_pad'] = 1
    else:
        w = hits['amplitude'].clip(lower=0)
        hits = hits.assign(_wx=w * hits['pad_cx'], _wy=w * hits['pad_cy'], _w=w)
        g = hits.groupby('eventId')
        cen = pd.DataFrame({
            'x_pad': g['_wx'].sum() / g['_w'].sum(),
            'y_pad': g['_wy'].sum() / g['_w'].sum(),
            'n_pad': g.size(),
        }).reset_index()
    cen = cen[np.isfinite(cen['x_pad']) & np.isfinite(cen['y_pad'])]
    print(f'P2: {len(cen):,} events with a pad centroid '
          f'({"leading-pad" if leading_pad else "charge-weighted"}).')
    return cen


# --------------------------------------------------------------------------- #
# Rotation scan
# --------------------------------------------------------------------------- #
def _corr(a, b):
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def rotation_scan(x_m3, y_m3, x_pad, y_pad, dtheta=0.5):
    """Scan theta over [0,360) for mirror=+1 and mirror=-1.

    Returns dict {mirror: (thetas, r_xx, r_yy, fom)} and the best (mirror, theta,
    r_xx, r_yy, fom).
    """
    thetas = np.arange(0.0, 360.0, dtheta)
    rad = np.radians(thetas)
    cos, sin = np.cos(rad), np.sin(rad)
    out = {}
    best = None
    # centre the pad cloud so rotation is about its own centroid (numerically
    # cleaner; correlation is translation-invariant so this does not bias FOM).
    xp0 = x_pad - x_pad.mean()
    yp0 = y_pad - y_pad.mean()
    for mirror in (+1, -1):
        yp = mirror * yp0
        r_xx = np.empty_like(thetas)
        r_yy = np.empty_like(thetas)
        for i in range(len(thetas)):
            xr = xp0 * cos[i] - yp * sin[i]
            yr = xp0 * sin[i] + yp * cos[i]
            r_xx[i] = _corr(x_m3, xr)
            r_yy[i] = _corr(y_m3, yr)
        fom = r_xx + r_yy
        out[mirror] = (thetas, r_xx, r_yy, fom)
        j = int(np.argmax(fom))
        cand = (mirror, thetas[j], r_xx[j], r_yy[j], fom[j])
        if best is None or cand[4] > best[4]:
            best = cand
    return out, best


def procrustes(x_m3, y_m3, x_pad, y_pad):
    """Best-fit rigid transform (rotation + optional reflection + translation +
    isotropic scale) mapping pad -> m3, via SVD. Returns (theta_deg, reflection,
    scale, rmse)."""
    P = np.column_stack([x_pad - x_pad.mean(), y_pad - y_pad.mean()])
    Q = np.column_stack([x_m3 - x_m3.mean(), y_m3 - y_m3.mean()])
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    reflection = np.linalg.det(R) < 0
    scale = S.sum() / (P ** 2).sum()
    theta = np.degrees(np.arctan2(R[1, 0], R[0, 0]))
    fit = scale * (P @ R.T)
    rmse = float(np.sqrt(np.mean(np.sum((fit - Q) ** 2, axis=1))))
    return theta % 360.0, reflection, float(scale), rmse


# --------------------------------------------------------------------------- #
# z-height alignment (scan z, re-optimise the angle at each z, iterate)
# --------------------------------------------------------------------------- #
def _fom_over_theta(x_m3, y_m3, xp0, yp0, thetas):
    """Vectorised: Pearson r_xx, r_yy for every theta at once.
    xp0,yp0 must already be mean-centred; x_m3,y_m3 are the (fixed) targets."""
    rad = np.radians(thetas)
    cos, sin = np.cos(rad)[:, None], np.sin(rad)[:, None]
    xr = xp0[None, :] * cos - yp0[None, :] * sin        # (T, N)
    yr = xp0[None, :] * sin + yp0[None, :] * cos
    def rows_corr(A, b):
        Ac = A - A.mean(axis=1, keepdims=True)
        bc = b - b.mean()
        num = (Ac * bc).sum(axis=1)
        den = np.sqrt((Ac ** 2).sum(axis=1) * (bc ** 2).sum())
        return np.where(den > 0, num / den, 0.0)
    return rows_corr(xr, x_m3), rows_corr(yr, y_m3)


def zscan(ep_matched, x_pad, y_pad, z_grid, thetas, mirror):
    """For each z: reproject M3, return FOM(z, theta) and best (theta, fom) per z.
    ep_matched is the endpoint frame aligned row-for-row with x_pad/y_pad."""
    xp0 = x_pad - x_pad.mean()
    yp0 = mirror * (y_pad - y_pad.mean())
    fom2d = np.empty((len(z_grid), len(thetas)))
    best_theta = np.empty(len(z_grid))
    best_fom = np.empty(len(z_grid))
    for i, z in enumerate(z_grid):
        xm, ym = pa.project_to_z(ep_matched, z)
        r_xx, r_yy = _fom_over_theta(xm, ym, xp0, yp0, thetas)
        fom = r_xx + r_yy
        fom2d[i] = fom
        j = int(np.argmax(fom))
        best_theta[i] = thetas[j]
        best_fom[i] = fom[j]
    return fom2d, best_theta, best_fom


def joint_zt_align(ep_matched, x_pad, y_pad, mirror,
                   z0=232.0, z_half=250.0, iters=3):
    """Coarse-to-fine coordinate descent on (z, theta). Returns best (z, theta,
    fom) plus the coarse 2D grid (z_grid, thetas, fom2d) for the heatmap."""
    # coarse 2D grid for the picture + starting point
    z_grid = np.arange(z0 - z_half, z0 + z_half + 1, 10.0)
    thetas = np.arange(0.0, 360.0, 1.0)
    fom2d, bt, bf = zscan(ep_matched, x_pad, y_pad, z_grid, thetas, mirror)
    k = int(np.argmax(bf))
    z_best, t_best = float(z_grid[k]), float(bt[k])

    # refine: alternate fine z-scan (fixed theta) and fine theta-scan (fixed z)
    span_z, span_t = 40.0, 10.0
    for _ in range(iters):
        zg = np.arange(z_best - span_z, z_best + span_z + 0.1, 1.0)
        f2, _, bfz = zscan(ep_matched, x_pad, y_pad, zg,
                           np.array([t_best]), mirror)
        z_best = float(zg[int(np.argmax(bfz))])
        tg = np.arange(t_best - span_t, t_best + span_t + 0.01, 0.1)
        f2, _, bft = zscan(ep_matched, x_pad, y_pad, np.array([z_best]),
                           tg, mirror)
        t_best = float(tg[int(np.argmax(f2[0]))])
        span_z, span_t = span_z / 2, span_t / 2
    # final fom at the refined point
    xm, ym = pa.project_to_z(ep_matched, z_best)
    r_xx, r_yy = _fom_over_theta(xm, ym, x_pad - x_pad.mean(),
                                 mirror * (y_pad - y_pad.mean()),
                                 np.array([t_best]))
    return (z_best, t_best, float(r_xx[0] + r_yy[0]), float(r_xx[0]),
            float(r_yy[0])), (z_grid, thetas, fom2d)


def plot_zscan(grid, z_best, t_best, out_dir, tag, suffix=''):
    z_grid, thetas, fom2d = grid
    fig, (ax, axl) = plt.subplots(1, 2, figsize=(15, 5.5),
                                  gridspec_kw={'width_ratios': [1.5, 1]})
    im = ax.imshow(fom2d, origin='lower', aspect='auto', cmap='magma',
                   extent=[thetas[0], thetas[-1], z_grid[0], z_grid[-1]])
    cb = fig.colorbar(im, ax=ax, label='FOM = corr(x)+corr(y)')
    ax.scatter([t_best], [z_best], marker='*', s=180, c='cyan',
               edgecolors='k', label=f'best ({t_best:.1f}°, {z_best:.0f} mm)')
    ax.set_xlabel('rotation θ about z  [deg]')
    ax.set_ylabel('detector plane z  [mm]')
    ax.set_title('Joint z × θ alignment FOM'); ax.legend(loc='upper right', fontsize=8)

    # FOM vs z at the best theta (a 1-D slice through the map)
    ti = int(np.argmin(np.abs(thetas - t_best)))
    axl.plot(z_grid, fom2d[:, ti], 'o-', ms=3, color='purple')
    axl.axvline(z_best, color='crimson', ls='--', label=f'z* = {z_best:.0f} mm')
    axl.set_xlabel('detector plane z  [mm]'); axl.set_ylabel('FOM at best θ')
    axl.set_title('z sensitivity (slice at best θ)')
    axl.grid(True, alpha=0.3); axl.legend(fontsize=8)
    fig.suptitle(f'P2↔M3 z-height alignment — {tag}', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = os.path.join(out_dir, f'z_theta_scan{suffix}.png')
    fig.savefig(p, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {p}')


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
def plot_scan(scan, best, out_dir, tag, suffix=''):
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), sharey=True)
    for ax, mirror in zip(axes, (+1, -1)):
        thetas, r_xx, r_yy, fom = scan[mirror]
        ax.plot(thetas, r_xx, lw=1, alpha=0.7, label='corr(x_m3, x\')')
        ax.plot(thetas, r_yy, lw=1, alpha=0.7, label='corr(y_m3, y\')')
        ax.plot(thetas, fom, lw=1.8, color='k', label='FOM = sum')
        ax.axhline(0, color='grey', lw=0.6)
        if best[0] == mirror:
            ax.axvline(best[1], color='crimson', ls='--',
                       label=f'peak {best[1]:.1f}°')
        ax.set_xlabel('rotation θ about z  [deg]')
        ax.set_title(f'mirror = {"+1 (as-is)" if mirror>0 else "-1 (flip y)"}',
                     fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc='lower right')
    axes[0].set_ylabel('Pearson correlation')
    fig.suptitle(f'P2↔M3 rotation scan — {tag}', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = os.path.join(out_dir, f'rotation_scan{suffix}.png')
    fig.savefig(p, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {p}')


def plot_best_corr(x_m3, y_m3, x_pad, y_pad, best, out_dir, tag, suffix=''):
    mirror, theta, r_xx, r_yy, _ = best
    rad = np.radians(theta)
    xp0 = x_pad - x_pad.mean()
    yp0 = mirror * (y_pad - y_pad.mean())
    xr = xp0 * np.cos(rad) - yp0 * np.sin(rad)
    yr = xp0 * np.sin(rad) + yp0 * np.cos(rad)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 6))
    hb1 = a1.hexbin(xr, x_m3, gridsize=60, cmap='viridis', mincnt=1)
    fig.colorbar(hb1, ax=a1, label='matched events / bin')
    a1.set_xlabel("P2 x'  (rotated pad centroid) [mm]")
    a1.set_ylabel('M3 x_m3 [mm]')
    a1.set_title(f'x correlation  r={r_xx:.3f}', fontsize=10)
    hb2 = a2.hexbin(yr, y_m3, gridsize=60, cmap='viridis', mincnt=1)
    fig.colorbar(hb2, ax=a2, label='matched events / bin')
    a2.set_xlabel("P2 y'  (rotated pad centroid) [mm]")
    a2.set_ylabel('M3 y_m3 [mm]')
    a2.set_title(f'y correlation  r={r_yy:.3f}', fontsize=10)
    fig.suptitle(f'P2↔M3 at best angle θ={theta:.1f}°, mirror={mirror} — {tag}',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = os.path.join(out_dir, f'best_angle_corr{suffix}.png')
    fig.savefig(p, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {p}')


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description='P2<->M3 rotation-scan alignment.')
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'],
                    help='within-connector ordering (default reverse, validated vs M3)')
    ap.add_argument('--z', type=float, default=None,
                    help='nominal P2 plane z [mm]; scan starts here. Default: '
                         'the detector det_center z from run_config.json '
                         '(232 if unavailable), so P2_2 at p2_z works too.')
    ap.add_argument('--scan-z', action=argparse.BooleanOptionalAction, default=True,
                    help='fit the detector z height jointly with the angle '
                         '(default on; --no-scan-z to fix z at --z).')
    ap.add_argument('--z-half', type=float, default=250.0,
                    help='half-width of the coarse z scan window [mm] (default 250).')
    ap.add_argument('--dtheta', type=float, default=0.5,
                    help='rotation scan step [deg] (default 0.5).')
    ap.add_argument('--chi2-cut', type=float, default=qa.M3_CHI2_CUT)
    ap.add_argument('--min-amp', type=float, default=0.0,
                    help='minimum P2 hit amplitude [ADC].')
    ap.add_argument('--m3-fiducial', type=float, default=1500.0,
                    help='keep |x_m3|,|y_m3| < this [mm] to drop wild tracks '
                         '(default 1500; 0 disables).')
    ap.add_argument('--leading-pad', action='store_true',
                    help='use the highest-amplitude pad instead of the '
                         'charge-weighted centroid.')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True,
                    help='drop events taken during an HV spark (default on).')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    if args.z is None:
        args.z = cfg.det_plane_z()
        print(f'Nominal plane z from run_config: {args.z:.0f} mm')
    out_dir = cfg.out_dir('03_m3_alignment')

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy=args.strategy,
                                  drop_connectors=cfg.DEAD_CONNECTORS)

    # matched sample carries the FULL M3 endpoints so tracks can be reprojected
    # to any z during the height scan.
    sv = ps.SparkVeto.from_cfg(cfg) if args.veto_sparks else None
    ep = pa.load_m3_endpoints(cfg.m3_tracking_dir, args.chi2_cut)
    p2 = load_p2_centroids(cfg.combined_hits_dir, ct, args.min_amp,
                           args.leading_pad, spark_veto=sv)
    m = ep.merge(p2, on='eventId', how='inner')
    x0, y0 = pa.project_to_z(m, args.z)
    if args.m3_fiducial > 0:
        keep = (np.abs(x0) < args.m3_fiducial) & (np.abs(y0) < args.m3_fiducial)
        m = m[keep].reset_index(drop=True)
    print(f'Matched sample: {len(m):,} events (M3 single-track ∩ P2 hit).')
    if len(m) < 50:
        print('Too few matched events to scan reliably.')
        return

    x_pad = m['x_pad'].to_numpy()
    y_pad = m['y_pad'].to_numpy()
    tag = f'{cfg.DET_TAG} {cfg.RUN}/{cfg.SUB_RUN} [{args.strategy}]'
    sfx = cfg.product_suffix(args.veto_sparks) + qa.chi2_tag(args.chi2_cut)

    # -- joint z x theta alignment (scan z, re-optimise angle at each z, iterate) --
    z_used = args.z
    if args.scan_z:
        best_joint, best_grid = None, None
        for mir in (+1, -1):
            res, grid = joint_zt_align(m, x_pad, y_pad, mir, z0=args.z,
                                       z_half=args.z_half)
            if best_joint is None or res[2] > best_joint[0][2]:
                best_joint, best_grid, best_mir = (res, mir), grid, mir
        (z_used, t_z, fom_z, rxx_z, ryy_z), _ = best_joint
        print(f'z-height fit: z* = {z_used:.1f} mm, θ* = {t_z:.2f}°, '
              f'mirror {best_mir:+d}, FOM {fom_z:.4f} (r_x {rxx_z:.3f}, r_y {ryy_z:.3f})')
        plot_zscan(best_grid, z_used, t_z, out_dir, tag, sfx)

    # project the matched tracks to the fitted plane for the angle products
    xm, ym = pa.project_to_z(m, z_used)
    x_m3, y_m3 = np.asarray(xm), np.asarray(ym)

    scan, best = rotation_scan(x_m3, y_m3, x_pad, y_pad, args.dtheta)
    mirror, theta, r_xx, r_yy, fom = best
    theta_pc, refl_pc, scale_pc, rmse_pc = procrustes(x_m3, y_m3, x_pad, y_pad)

    plot_scan(scan, best, out_dir, tag, sfx)
    plot_best_corr(x_m3, y_m3, x_pad, y_pad, best, out_dir, tag, sfx)

    summary = [
        f'P2<->M3 alignment  {tag}',
        f'  matched events        : {len(m)}',
        f'  z-height fit          : {"ON" if args.scan_z else "OFF"}  '
        f'(nominal {args.z:.0f} mm)',
        f'  M3 projection z       : {z_used:.1f} mm',
        f'  centroid              : {"leading-pad" if args.leading_pad else "charge-weighted"}',
        '',
        'Rotation scan (best):',
        f'  mirror                : {mirror:+d}',
        f'  best angle theta      : {theta:.2f} deg',
        f'  corr(x_m3, x\')        : {r_xx:.4f}',
        f'  corr(y_m3, y\')        : {r_yy:.4f}',
        f'  FOM (sum)             : {fom:.4f}',
        '',
        'Procrustes cross-check (best-fit rigid transform pad->m3):',
        f'  implied rotation      : {theta_pc:.2f} deg',
        f'  reflection            : {refl_pc}',
        f'  isotropic scale       : {scale_pc:.4f}',
        f'  residual RMSE         : {rmse_pc:.1f} mm',
    ]
    print('\n'.join(summary))
    sp = os.path.join(out_dir, f'alignment_summary{sfx}.txt')
    with open(sp, 'w') as fh:
        fh.write('\n'.join(summary) + '\n')
    print(f'  saved {sp}')
    print('Done.')


if __name__ == '__main__':
    main()
