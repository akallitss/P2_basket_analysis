#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
10_efficiency_map_sliding.py

Smooth sliding-window efficiency map for the P2 pad detector, adapted from
nTof_x17/mx_june_cosmic_qa/12_efficiency_map_sliding.py.

Consumes the per-ray hit/miss list written by 06_efficiency_maps.py
(<06_efficiency>/ray_hit_miss_list<suffix>.csv):

    columns: eventId, x, y, in_active, has_any, within, resid
      x, y     = clean M3 single-track projection at the aligned P2 plane
                 (reference frame, mm)
      within   = a P2 pad centroid within R mm of the projection (R baked into 06,
                 default 20 mm) -> efficiency numerator
      has_any  = P2 fired any pad in that event

Instead of hard 2D bins, every point on a regular grid averages `within` /
`has_any` over the rays inside a kernel, so overlapping kernels give a smooth
map. Two kernels:

  * FIXED    : all rays within --kernel mm of the grid point (needs --min rays).
  * ADAPTIVE : (--adaptive) k-NN kernel -- the smallest radius capturing --target
               rays, so the resolution follows the local ray density and every
               point carries the same binomial error. The kernel-radius panel is
               then the local-resolution map (capped at --max-kernel).

The suffix (which CSV is read and how the outputs are named) follows the same
dead-connector / spark-veto convention as the other stages, so the sliding map
matches the efficiency map it is derived from.

Products (written next to the efficiency maps, <06_efficiency>/):
  efficiency_map_sliding<suffix>.png / .json      (fixed kernel)
  efficiency_map_adaptive<suffix>.png / .json     (--adaptive)

Usage:
  python3 10_efficiency_map_sliding.py [run_key] [--kernel 25] [--grid 120]
         [--min 30] [--r 20] [--adaptive --target 60 --max-kernel 30]
         [--no-veto-sparks]
"""

import os
import json
import argparse
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import p2_qa_config as qa
import p2_mapping as pmap


def sliding_map(x, y, val, x_grid, y_grid, kernel, min_n):
    """Mean of `val` over rays within `kernel` mm of each grid point."""
    r2 = kernel ** 2
    eff = np.full((len(x_grid), len(y_grid)), np.nan)
    cnt = np.zeros_like(eff, dtype=int)
    for i, xg in enumerate(x_grid):
        dx2 = (x - xg) ** 2
        for j, yg in enumerate(y_grid):
            mask = (dx2 + (y - yg) ** 2) <= r2
            n = int(mask.sum())
            cnt[i, j] = n
            if n >= min_n:
                eff[i, j] = float(val[mask].mean())
    return eff, cnt


def adaptive_map(x, y, within, has_any, x_grid, y_grid, target, max_r):
    """k-NN adaptive kernel. Each grid point takes the `target` nearest rays; the
    kernel radius is the distance to the target-th ray (= local resolution).
    Returns (efficiency, has_any, radius); points with radius > max_r are NaN."""
    from scipy.spatial import cKDTree
    tree = cKDTree(np.column_stack([x, y]))
    gx, gy = np.meshgrid(x_grid, y_grid, indexing='ij')
    gpts = np.column_stack([gx.ravel(), gy.ravel()])
    k = min(target, len(x))
    dist, idx = tree.query(gpts, k=k)
    if k == 1:
        dist = dist[:, None]; idx = idx[:, None]
    radius = dist[:, -1]
    eff = within[idx].mean(axis=1)
    anyh = has_any[idx].mean(axis=1)
    bad = radius > max_r
    eff[bad] = np.nan; anyh[bad] = np.nan; radius[bad] = np.nan
    shp = (len(x_grid), len(y_grid))
    return eff.reshape(shp), anyh.reshape(shp), radius.reshape(shp)


def _outside_mask(x_grid, y_grid, ftree, mask_r):
    """Boolean grid: True where the nearest transformed pad is farther than
    mask_r mm -> the grid point is outside the detector footprint."""
    from scipy.spatial import cKDTree  # noqa: F401 (ftree already a cKDTree)
    gx, gy = np.meshgrid(x_grid, y_grid, indexing='ij')
    nn, _ = ftree.query(np.column_stack([gx.ravel(), gy.ravel()]))
    return (nn > mask_r).reshape(len(x_grid), len(y_grid))


def _draw_footprint(ax, fpx, fpy):
    """Outline the transformed pad footprint (convex hull) as the detector edge."""
    try:
        from scipy.spatial import ConvexHull
        h = ConvexHull(np.column_stack([fpx, fpy]))
        v = np.append(h.vertices, h.vertices[0])
        ax.plot(fpx[v], fpy[v], color='red', lw=1.2, alpha=0.8)
    except Exception:
        pass


def main():
    ap = argparse.ArgumentParser(description='P2 sliding-window efficiency map.')
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--kernel', type=float, default=25.0,
                    help='fixed smoothing kernel radius [mm] (default 25).')
    ap.add_argument('--grid', type=int, default=120,
                    help='grid points per axis (default 120).')
    ap.add_argument('--min', dest='min_rays', type=int, default=30,
                    help='min rays in a fixed kernel to colour a grid point (default 30).')
    ap.add_argument('--r', type=float, default=None,
                    help='match radius label [mm]; the match is baked into "within" '
                         'by stage 06 (default = the run config MATCH_R).')
    ap.add_argument('--adaptive', action='store_true',
                    help='use an adaptive k-NN kernel instead of a fixed radius.')
    ap.add_argument('--target', type=int, default=60,
                    help='rays per adaptive kernel (sets the binomial error, default 60).')
    ap.add_argument('--max-kernel', type=float, default=30.0,
                    help='cap on the adaptive kernel radius [mm] (default 30).')
    ap.add_argument('--mask-r', type=float, default=15.0,
                    help='grid points farther than this [mm] from any transformed '
                         'pad are masked (outside the detector) (default 15).')
    ap.add_argument('--margin', type=float, default=15.0,
                    help='grid margin beyond the pad footprint [mm] so the detector '
                         'borders are not clipped (default 15).')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True,
                    help='read the spark-vetoed ray list (default on), matching 06.')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    if args.r is None:
        args.r = cfg.MATCH_R
    eff_dir = cfg.out_dir('06_efficiency')
    suffix = cfg.product_suffix(args.veto_sparks)
    csv = os.path.join(eff_dir, f'ray_hit_miss_list{suffix}.csv')
    if not os.path.isfile(csv):
        raise FileNotFoundError(
            f'{csv} not found — run 06_efficiency_maps.py '
            f'{"" if args.veto_sparks else "--no-veto-sparks "}first.')

    d = pd.read_csv(csv)
    d = d[np.isfinite(d['x']) & np.isfinite(d['y'])].copy()
    for c in ('within', 'has_any'):
        d[c] = d[c].astype(str).str.lower().isin(('true', '1'))
    x = d['x'].to_numpy(float)
    y = d['y'].to_numpy(float)
    within = d['within'].to_numpy(float)
    has_any = d['has_any'].to_numpy(float)

    # zone = the real transformed pad footprint (from stage 06), not a ray box.
    from scipy.spatial import cKDTree
    fp_csv = os.path.join(eff_dir, f'pad_footprint{suffix}.csv')
    if not os.path.isfile(fp_csv):
        raise FileNotFoundError(
            f'{fp_csv} not found — re-run 06_efficiency_maps.py to write the pad '
            'footprint used to define the detector zone.')
    fp = pd.read_csv(fp_csv)
    fpx, fpy = fp['x'].to_numpy(float), fp['y'].to_numpy(float)
    ftree = cKDTree(np.column_stack([fpx, fpy]))

    # insulation-mask pillars (M3 frame, written by stage 06): pillar shadows
    # vs real dead zones on the smooth map.
    pil_csv = os.path.join(eff_dir, f'pillars_m3{suffix}.csv')
    pil = pd.read_csv(pil_csv) if os.path.isfile(pil_csv) else None
    ax0, ax1 = float(fpx.min()), float(fpx.max())
    ay0, ay1 = float(fpy.min()), float(fpy.max())

    # integrated numbers over rays that actually cross the detector footprint
    nn_ray, _ = ftree.query(np.column_stack([x, y]))
    inact = nn_ray <= args.mask_r
    integ_within = float(within[inact].mean()) if inact.any() else float('nan')
    integ_hasany = float(has_any[inact].mean()) if inact.any() else float('nan')
    tag = f'{cfg.DET_TAG} {cfg.RUN}/{cfg.SUB_RUN}'
    cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('lightgrey')

    if args.adaptive:
        grid_n = max(args.grid, 200)
        x_grid = np.linspace(ax0 - args.margin, ax1 + args.margin, grid_n)
        y_grid = np.linspace(ay0 - args.margin, ay1 + args.margin, grid_n)
        print(f'{tag}: ADAPTIVE k-NN, target={args.target} rays/kernel, '
              f'maxR={args.max_kernel:.0f} mm, {grid_n}x{grid_n} grid, {len(d)} rays')
        eff_w, eff_a, radius = adaptive_map(x, y, within, has_any, x_grid, y_grid,
                                            args.target, args.max_kernel)
        outside = _outside_mask(x_grid, y_grid, ftree, args.mask_r)
        eff_w[outside] = np.nan; eff_a[outside] = np.nan; radius[outside] = np.nan
        rad = radius[np.isfinite(radius)]
        med = float(np.median(rad)) if len(rad) else float('nan')
        p90 = float(np.percentile(rad, 90)) if len(rad) else float('nan')
        if len(rad):
            print(f'  kernel radius: median {med:.1f}, p90 {p90:.1f}, max '
                  f'{rad.max():.1f} mm  (fixed-kernel equiv: --kernel {p90:.0f} '
                  f'--min {args.target})')
        extent = [x_grid[0], x_grid[-1], y_grid[0], y_grid[-1]]
        cmap_r = plt.get_cmap('turbo').copy(); cmap_r.set_bad('lightgrey')
        fig, axes = plt.subplots(1, 3, figsize=(19, 6))
        for ax, data, label in [
                (axes[0], eff_w, f'efficiency (reco within {args.r:g} mm)'),
                (axes[1], eff_a, 'has_any (fired any pad)')]:
            im = ax.imshow(data.T, origin='lower', extent=extent, aspect='equal',
                           cmap=cmap, vmin=0, vmax=1)
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=label)
            _draw_footprint(ax, fpx, fpy)
            if pil is not None:
                pmap.draw_pillars(ax, pil, small=False)
                ax.legend(loc='upper right', fontsize=7, framealpha=0.9)
            ax.set_xlabel('reference X [mm]'); ax.set_ylabel('reference Y [mm]')
            ax.set_title(f'{cfg.DET_NAME}  {label}\nadaptive k-NN, {args.target} rays/kernel')
        im3 = axes[2].imshow(radius.T, origin='lower', extent=extent, aspect='equal',
                             cmap=cmap_r, vmin=0, vmax=args.max_kernel)
        plt.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04,
                     label='kernel radius [mm] (= local resolution)')
        _draw_footprint(axes[2], fpx, fpy)
        axes[2].set_xlabel('reference X [mm]'); axes[2].set_ylabel('reference Y [mm]')
        axes[2].set_title(f'local kernel radius\n(smaller = finer; median {med:.1f} mm)')
        fig.suptitle(f'{cfg.DET_NAME} ADAPTIVE-kernel efficiency — {tag}',
                     y=1.02, fontsize=13)
        fig.tight_layout()
        out_png = os.path.join(eff_dir, f'efficiency_map_adaptive{suffix}.png')
        fig.savefig(out_png, dpi=150, bbox_inches='tight'); plt.close(fig)
        summary = dict(det=cfg.DET_NAME, run=cfg.RUN, sub_run=cfg.SUB_RUN,
                       mode='adaptive', target_rays=args.target,
                       max_kernel_mm=args.max_kernel, grid=grid_n, n_rays=int(len(d)),
                       n_rays_active=int(inact.sum()), integrated_within=integ_within,
                       integrated_has_any=integ_hasany, kernel_radius_median_mm=med,
                       kernel_radius_p90_mm=p90, spark_vetoed=bool(args.veto_sparks))
        with open(os.path.join(eff_dir, f'efficiency_map_adaptive{suffix}.json'), 'w') as f:
            json.dump(summary, f, indent=2)
        print(f'  written: {out_png}')
        return

    x_grid = np.linspace(ax0 - args.margin, ax1 + args.margin, args.grid)
    y_grid = np.linspace(ay0 - args.margin, ay1 + args.margin, args.grid)
    print(f'{tag}: sliding kernel r={args.kernel:.0f} mm, {args.grid}x{args.grid} '
          f'grid, {len(d)} rays')
    eff_w, cnt = sliding_map(x, y, within, x_grid, y_grid, args.kernel, args.min_rays)
    eff_a, _ = sliding_map(x, y, has_any, x_grid, y_grid, args.kernel, args.min_rays)
    outside = _outside_mask(x_grid, y_grid, ftree, args.mask_r)
    eff_w[outside] = np.nan; eff_a[outside] = np.nan
    cnt = np.where(outside, 0, cnt)
    n_fit = int(np.sum(~np.isnan(eff_w)))
    print(f'  fitted {n_fit}/{args.grid**2} grid points; integrated '
          f'within{args.r:g}mm={integ_within*100:.1f}%  has_any={integ_hasany*100:.1f}%')

    extent = [x_grid[0], x_grid[-1], y_grid[0], y_grid[-1]]
    cmap_c = plt.get_cmap('plasma').copy(); cmap_c.set_bad('lightgrey')
    fig, axes = plt.subplots(1, 3, figsize=(19, 6))
    for ax, data, label in [
            (axes[0], eff_w, f'efficiency (reco within {args.r:g} mm)'),
            (axes[1], eff_a, 'has_any (fired any pad)')]:
        im = ax.imshow(data.T, origin='lower', extent=extent, aspect='equal',
                       cmap=cmap, vmin=0, vmax=1)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=label)
        _draw_footprint(ax, fpx, fpy)
        if pil is not None:
            pmap.draw_pillars(ax, pil, small=False)
            ax.legend(loc='upper right', fontsize=7, framealpha=0.9)
        ax.set_xlabel('reference X [mm]'); ax.set_ylabel('reference Y [mm]')
        ax.set_title(f'{cfg.DET_NAME}  {label}\nsliding kernel r={args.kernel:.0f} mm')
    cnt_m = np.where((cnt >= args.min_rays) & ~outside, cnt, np.nan)
    im3 = axes[2].imshow(cnt_m.T, origin='lower', extent=extent, aspect='equal', cmap=cmap_c)
    plt.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04, label='rays in kernel')
    _draw_footprint(axes[2], fpx, fpy)
    axes[2].set_xlabel('reference X [mm]'); axes[2].set_ylabel('reference Y [mm]')
    axes[2].set_title(f'rays per kernel\n(grey < {args.min_rays})')
    fig.suptitle(f'{cfg.DET_NAME} sliding-window efficiency — {tag}', y=1.02, fontsize=13)
    fig.tight_layout()
    out_png = os.path.join(eff_dir, f'efficiency_map_sliding{suffix}.png')
    fig.savefig(out_png, dpi=150, bbox_inches='tight'); plt.close(fig)

    summary = dict(det=cfg.DET_NAME, run=cfg.RUN, sub_run=cfg.SUB_RUN,
                   r_mm=args.r, kernel_mm=args.kernel, grid=args.grid,
                   n_rays=int(len(d)), n_rays_active=int(inact.sum()),
                   integrated_within=integ_within, integrated_has_any=integ_hasany,
                   spark_vetoed=bool(args.veto_sparks),
                   active_box=dict(x0=ax0, x1=ax1, y0=ay0, y1=ay1))
    with open(os.path.join(eff_dir, f'efficiency_map_sliding{suffix}.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(f'  written: {out_png}')


if __name__ == '__main__':
    main()
