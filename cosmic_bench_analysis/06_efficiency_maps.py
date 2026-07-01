#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_efficiency_maps.py

Reference-track based efficiency for the P2 pad detector, following the method
of nTof_x17/mx_june_cosmic_qa/08_efficiency_maps.py + 09_efficiency_breakdown.py,
adapted for pads:

For every clean M3 single track, project it to the P2 plane (z from run_config)
and, using the validated pad->M3 transform (p2_align, reverse ordering), record
  within : P2 has a reconstructed pad centroid within R mm of the projection
  has_any: P2 fired ANY pad in that event (detector responded at all)
A ray with no P2 event (DAQ only saves events with a valid hit) is a genuine
MISS (both False) and stays in the denominator.

The P2 "reconstructed position" is the per-event charge-weighted pad centroid
mapped into the M3 frame -- the pad analogue of the micro-TPC hit the reference
uses. Active area is the transformed pad footprint (fan-aware, via nearest-pad
distance), not a bounding box.

Products (written to <Analysis>/<detN>/<run>/<sub_run>/06_efficiency/):
  scatter_within_Rmm.png / scatter_has_any.png   green hit / red miss at projections
  map_within_Rmm.png / map_has_any.png           binned efficiency map (>=5 rays/bin)
  radial_residual.png                            |r| residual (core + tail)
  reco_positions_detector.png                    reconstructed positions in pad frame
  nonreco_ray_positions.png                      projections of muons P2 did not see
  efficiency_breakdown.png                       where do crossing muons go?
  efficiency_summary.txt / ray_hit_miss_list.csv

Usage: python3 06_efficiency_maps.py [run_key] [--r 20] [--active-r 30]
"""

import os
import json
import argparse
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

import p2_qa_config as qa
import p2_mapping as pmap
import p2_align as pa

BINS = 40


def _det_plane_z(cfg):
    with open(cfg.run_config_path) as f:
        for d in json.load(f)['detectors']:
            if d.get('name') == cfg.DET_NAME:
                return float(d['det_center_coords']['z'])
    return 232.0


def main():
    ap = argparse.ArgumentParser(description='P2 reference-track efficiency maps.')
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'])
    ap.add_argument('--r', type=float, default=20.0,
                    help='match radius [mm] for "within" (default 20 ~ 1.7 pad pitch).')
    ap.add_argument('--active-r', type=float, default=30.0,
                    help='a ray is "in active area" if within this of a transformed '
                         'pad centre [mm] (default 30).')
    ap.add_argument('--chi2-cut', type=float, default=1.5)
    ap.add_argument('--fit-fiducial', type=float, default=300.0,
                    help='|x_m3|,|y_m3| window used only to FIT the transform.')
    ap.add_argument('--z', type=float, default=None,
                    help='override the M3 projection plane z [mm]; default uses '
                         'run_config det_center z. Use the z fitted by '
                         '03_m3_alignment (e.g. 246) for best alignment.')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    out_dir = cfg.out_dir('06_efficiency')
    det_z = args.z if args.z is not None else _det_plane_z(cfg)
    R = args.r

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy=args.strategy)

    # --- inputs ---
    m3 = pa.load_m3_positions(cfg.m3_tracking_dir, det_z, args.chi2_cut)
    p2, hit_events = pa.load_p2_centroids(cfg.combined_hits_dir, ct)
    print(f'M3 single-track rays: {len(m3):,} | P2 events with centroid: {len(p2):,}')

    # --- fit the pad->M3 transform on matched events (fiducial only for the fit) ---
    matched = m3.merge(p2, on='eventId', how='inner')
    fitm = matched
    if args.fit_fiducial > 0:
        fitm = matched[(matched['x_m3'].abs() < args.fit_fiducial) &
                       (matched['y_m3'].abs() < args.fit_fiducial)]
    T = pa.fit_transform(fitm['x_m3'], fitm['y_m3'], fitm['x_pad'], fitm['y_pad'])
    print(f'transform: rotation {T.rotation_deg:.2f} deg, reflection {T.reflection}, '
          f'scale {T.s:.3f}, fit RMSE {T.rmse:.1f} mm (N={len(fitm)})')

    # --- P2 reco position (pad centroid) in the M3 frame, per event ---
    rx, ry = T.apply(p2['x_pad'].to_numpy(), p2['y_pad'].to_numpy())
    reco = {int(e): (xx, yy) for e, xx, yy in zip(p2['eventId'], rx, ry)}

    # --- transformed pad footprint -> active-area KD-tree ---
    padc = ct.drop_duplicates('channel_id')
    pcx, pcy = T.apply(padc['pad_cx'].to_numpy(), padc['pad_cy'].to_numpy())
    pad_xy = np.column_stack([pcx, pcy])
    tree = cKDTree(pad_xy)

    # --- ray list (every clean M3 single track) with hit/miss flags ---
    d = m3.rename(columns={'x_m3': 'x', 'y_m3': 'y'}).copy()
    nn_dist, _ = tree.query(np.column_stack([d['x'], d['y']]))
    d['in_active'] = nn_dist <= args.active_r
    d['has_any'] = d['eventId'].isin(hit_events)
    within = np.zeros(len(d), dtype=bool)
    resid = np.full(len(d), np.nan)
    evids = d['eventId'].to_numpy()
    xs, ys = d['x'].to_numpy(), d['y'].to_numpy()
    for i in range(len(d)):
        rc = reco.get(int(evids[i]))
        if rc is not None:
            resid[i] = np.hypot(xs[i] - rc[0], ys[i] - rc[1])
            within[i] = resid[i] <= R
    d['within'] = within
    d['resid'] = resid

    da = d[d['in_active']]
    n_act = len(da)
    eff_within = da['within'].mean() * 100 if n_act else float('nan')
    eff_any = da['has_any'].mean() * 100 if n_act else float('nan')
    print(f'active-area rays: {n_act:,}  |  within {R:g}mm = {eff_within:.1f}%  '
          f'has_any = {eff_any:.1f}%')

    # ---------------- plots ----------------
    xb = [pad_xy[:, 0].min() - 40, pad_xy[:, 0].max() + 40]
    yb = [pad_xy[:, 1].min() - 40, pad_xy[:, 1].max() + 40]

    # scatter hit/miss
    for col, name, ttl in [('within', f'within_{R:g}mm', f'hit within {R:g} mm'),
                           ('has_any', 'has_any', 'any hit on detector')]:
        fig, ax = plt.subplots(figsize=(7.5, 7))
        sub = d[d['in_active']]
        miss = sub[~sub[col]]; hit = sub[sub[col]]
        ax.scatter(miss['x'], miss['y'], s=6, c='red', alpha=0.15, linewidths=0,
                   label=f'miss ({len(miss)})')
        ax.scatter(hit['x'], hit['y'], s=6, c='green', alpha=0.15, linewidths=0,
                   label=f'hit ({len(hit)})')
        ax.scatter(pad_xy[:, 0], pad_xy[:, 1], s=2, c='k', alpha=0.15,
                   label='pad footprint')
        ax.set_xlabel('M3 X [mm]'); ax.set_ylabel('M3 Y [mm]'); ax.set_aspect('equal')
        ax.set_title(f'{cfg.DET_NAME} efficiency scatter — {ttl}\n{cfg.RUN}/{cfg.SUB_RUN}')
        lg = ax.legend(loc='upper right', framealpha=0.9)
        for h in lg.legend_handles:
            h.set_alpha(1)
        fig.tight_layout()
        fig.savefig(f'{out_dir}/scatter_{name}.png', dpi=150, bbox_inches='tight')
        plt.close(fig)

    # binned efficiency maps (over active-area rays)
    rng = [xb, yb]
    den, xe, ye = np.histogram2d(da['x'], da['y'], bins=BINS, range=rng)
    for col, name, ttl in [('within', f'within_{R:g}mm', f'hit within {R:g} mm'),
                           ('has_any', 'has_any', 'any hit')]:
        num, _, _ = np.histogram2d(da.loc[da[col], 'x'], da.loc[da[col], 'y'],
                                   bins=[xe, ye])
        with np.errstate(invalid='ignore', divide='ignore'):
            eff = np.where(den >= 5, num / den, np.nan)
        fig, ax = plt.subplots(figsize=(7.5, 6))
        cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('lightgrey')
        im = ax.imshow(eff.T, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]],
                       vmin=0, vmax=1, cmap=cmap, aspect='equal')
        plt.colorbar(im, ax=ax, label='efficiency')
        ax.set_xlabel('M3 X [mm]'); ax.set_ylabel('M3 Y [mm]')
        ax.set_title(f'{cfg.DET_NAME} efficiency map — {ttl}  (>=5 rays/bin)\n'
                     f'{cfg.RUN}/{cfg.SUB_RUN}')
        fig.tight_layout()
        fig.savefig(f'{out_dir}/map_{name}.png', dpi=150, bbox_inches='tight')
        plt.close(fig)

    # radial residual (core + tail) for active-area rays that had a P2 event
    rr = da['resid'].to_numpy()
    rr = rr[np.isfinite(rr)]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    axes[0].hist(rr[rr <= 60], bins=60, color='steelblue')
    axes[0].axvline(R, color='crimson', ls='--', label=f'R = {R:g} mm')
    axes[0].set_xlabel('|r| residual [mm]'); axes[0].set_ylabel('events')
    axes[0].set_title('|r| (0-60 mm)'); axes[0].legend(); axes[0].grid(True, alpha=0.3)
    axes[1].hist(rr, bins=80, color='indianred'); axes[1].set_yscale('log')
    axes[1].axvline(R, color='k', ls='--')
    axes[1].set_xlabel('|r| residual [mm]'); axes[1].set_ylabel('events (log)')
    axes[1].set_title('|r| full range (core + tail)'); axes[1].grid(True, alpha=0.3)
    med = np.median(rr) if len(rr) else float('nan')
    fig.suptitle(f'{cfg.DET_NAME} radial residual (active area, median {med:.1f} mm) '
                 f'— {cfg.RUN}/{cfg.SUB_RUN}')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/radial_residual.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # reconstructed positions in the pad (detector) frame, for well-reco rays
    good_ev = set(da.loc[da['within'], 'eventId'].astype(int))
    recop = p2[p2['eventId'].isin(good_ev)]
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(recop['x_pad'], recop['y_pad'], s=6, c='teal', alpha=0.2, linewidths=0)
    ax.set_aspect('equal'); ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
    ax.set_title(f'{cfg.DET_NAME} pad positions of well-reco tracks (|r|<={R:g}mm)\n'
                 f'{cfg.RUN}/{cfg.SUB_RUN}  ({len(recop):,} events)')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/reco_positions_detector.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # projections of non-reconstructed muons (in active area, no within-R hit)
    nonreco = da[~da['within']]
    fig, ax = plt.subplots(figsize=(7.5, 7))
    hb = ax.hexbin(nonreco['x'], nonreco['y'], gridsize=45, cmap='magma', mincnt=1)
    plt.colorbar(hb, ax=ax, label='non-reco muons')
    ax.scatter(pad_xy[:, 0], pad_xy[:, 1], s=2, c='cyan', alpha=0.2)
    ax.set_aspect('equal'); ax.set_xlabel('M3 X [mm]'); ax.set_ylabel('M3 Y [mm]')
    ax.set_title(f'{cfg.DET_NAME} projections of NON-reconstructed muons\n'
                 f'{cfg.RUN}/{cfg.SUB_RUN}  ({len(nonreco):,} rays)')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/nonreco_ray_positions.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # breakdown: where do active-area crossing muons go?
    n = n_act
    n_within = int(da['within'].sum())
    n_fired_not_within = int((da['has_any'] & ~da['within']).sum())
    n_nohit = int((~da['has_any']).sum())
    cats = ['reconstructed\n(|r|<=%g mm)' % R, 'fired but\nmis-reco', 'no hit\n(dead miss)']
    vals = [100 * n_within / n, 100 * n_fired_not_within / n, 100 * n_nohit / n]
    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar(cats, vals, color=['seagreen', 'goldenrod', 'firebrick'])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.5, f'{v:.1f}%', ha='center')
    ax.set_ylabel('% of active-area crossing muons')
    ax.set_title(f'{cfg.DET_NAME} where do the crossing muons go?\n'
                 f'{cfg.RUN}/{cfg.SUB_RUN}  (N={n:,})')
    ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/efficiency_breakdown.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # summary
    summary = [
        f'P2 efficiency — {cfg.DET_TAG} {cfg.RUN}/{cfg.SUB_RUN} [{args.strategy}]',
        f'  match radius R            : {R:g} mm',
        f'  active-area radius        : {args.active_r:g} mm',
        f'  transform rotation/scale  : {T.rotation_deg:.2f} deg / {T.s:.3f} '
        f'(reflection {T.reflection}, fit RMSE {T.rmse:.1f} mm)',
        f'  clean M3 rays (total)     : {len(d):,}',
        f'  rays in active area       : {n_act:,}',
        f'  INTEGRATED efficiency (active area):',
        f'     within {R:g} mm          : {eff_within:.1f}%',
        f'     has_any (fired)        : {eff_any:.1f}%',
        f'  breakdown: reco {vals[0]:.1f}% | fired-not-reco {vals[1]:.1f}% | '
        f'no-hit {vals[2]:.1f}%',
        f'  median |r| residual       : {med:.1f} mm',
    ]
    print('\n'.join(summary))
    with open(f'{out_dir}/efficiency_summary.txt', 'w') as f:
        f.write('\n'.join(summary) + '\n')
    d.to_csv(f'{out_dir}/ray_hit_miss_list.csv', index=False)
    print(f'\nWritten to: {out_dir}')


if __name__ == '__main__':
    main()
