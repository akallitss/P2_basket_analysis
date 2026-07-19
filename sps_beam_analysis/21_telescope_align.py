#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
21_telescope_align.py

Mutual alignment of the P2 telescope planes from cluster-position correlations
(replaces the bench's M3-referenced 03_m3_alignment). There is no external
tracker at the SPS, but the mu/pi beam at 80-150 GeV/c is essentially parallel
to z with negligible multiple scattering, so a beam particle hits every plane at
the SAME transverse (x, y) up to each plane's mechanical mounting offset. The
map between two planes (identical PCB geometry) is therefore a RIGID 2D
transform -- translation + rotation about z, no scale, no track fit needed.

Method (per sub_run)
--------------------
  1. For every station, stream the combined hits into one clean single cluster
     per event (sps_cluster.stream_event_clusters: leading pad + pads within
     --cluster-r, `single` = no pad outside that radius). HV-settle cut + HV
     spark veto applied.
  2. For each plane paired with the configured reference plane, select events
     with EXACTLY ONE clean cluster in BOTH planes (shared eventId, both
     `single`), giving 2D residuals dx, dy = probe - ref.
  3. Fit the rigid transform plane -> reference (scale-free Kabsch,
     sps_cluster.fit_rigid): translation (dx, dy) + rotation theta about z.
  4. Report residual widths before and after; write the alignment JSON that 22
     (tag-and-probe) consumes.

Products (<Analysis>/telescope/<run>/<sub_run>/21_telescope_align/):
  residuals_<PLANE>_to_<REF><suffix>.png   dx,dy hists + residual map, pre/post
  alignment<suffix>.png                    all planes' fitted offsets summary
  alignment.json                           {ref, planes:{name:{dx,dy,theta_deg,
                                            rmse_pre,rmse_post,n}}}

The reference plane gets the identity transform; a plane's transform maps its
cluster (x, y) into the reference frame. 22 composes these to carry a tag from
any plane into any probe plane.

Usage:
  python3 21_telescope_align.py [run_key] [--sub-run NAME] [--cluster-r 15]
        [--min-pairs 100] [--no-veto-sparks]
"""

import os
import json
import argparse

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

import sps_config as sc
import p2_io as p2io
import p2_sparks as ps
import sps_cluster as scl


def load_clusters(cfg, det, sub_run, cluster_r, veto_sparks):
    """Clean single clusters for one station in one sub_run (spark-vetoed,
    HV-settle cut). Returns a DataFrame keyed on eventId."""
    hits_dir = cfg.combined_hits_dir(sub_run)
    chunk0 = p2io.hit_files(hits_dir)[0]
    hv_csv = cfg.hv_monitor_csv(sub_run)
    t_min = scl.settle_t_min(hv_csv, det.spark_channel, chunk0)
    veto = None
    if veto_sparks and os.path.isfile(hv_csv) and det.spark_channel:
        veto = _spark_veto(cfg, det, hv_csv)
    ct = cfg.channel_table(det)
    ev = scl.stream_event_clusters(hits_dir, ct, cluster_r,
                                   min_amp=cfg.MIN_AMP, veto=veto, t_min=t_min)
    return ev


def _spark_veto(cfg, det, hv_csv):
    """Build a SparkVeto for one station (its own mesh channel/knobs)."""
    class _Shim:
        SPARK_CHANNEL = det.spark_channel
        SPARK_IMON_THR = cfg.SPARK_IMON_THR
        SPARK_GUARD_BEFORE = cfg.SPARK_GUARD_BEFORE
        SPARK_GUARD_AFTER = cfg.SPARK_GUARD_AFTER
        BURST_NPADS = cfg.BURST_NPADS
    return ps.SparkVeto.from_csv(hv_csv, _Shim)


def plot_residuals(dx, dy, dx2, dy2, xr, yr, tf, rmse_pre, rmse_post,
                   plane, ref, out_png):
    """dx,dy before / dx2,dy2 after alignment; residual map coloured pre-fit."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

    ax = axes[0]
    rng = np.nanpercentile(np.abs(np.concatenate([dx, dy])), 99) if len(dx) else 1
    b = np.linspace(-rng, rng, 60)
    ax.hist(dx, bins=b, histtype='step', color='#1f77b4', lw=1.6,
            label=f'dx pre (rms {np.std(dx):.1f})')
    ax.hist(dy, bins=b, histtype='step', color='#d62728', lw=1.6,
            label=f'dy pre (rms {np.std(dy):.1f})')
    ax.hist(dx2, bins=b, histtype='step', color='#1f77b4', lw=1.6, ls='--',
            label=f'dx post (rms {np.std(dx2):.1f})')
    ax.hist(dy2, bins=b, histtype='step', color='#d62728', lw=1.6, ls='--',
            label=f'dy post (rms {np.std(dy2):.1f})')
    ax.set_xlabel('residual [mm]'); ax.set_ylabel('event pairs')
    ax.set_title(f'{plane} - {ref} residuals, pre/post rigid fit')
    ax.grid(True, alpha=0.3); ax.legend(fontsize=7.5)

    ax = axes[1]
    sccol = np.hypot(dx, dy)
    sc_ = ax.scatter(xr, yr, c=sccol, s=8, cmap='viridis',
                     vmax=np.nanpercentile(sccol, 95) if len(sccol) else 1)
    fig.colorbar(sc_, ax=ax, label='|residual| pre [mm]')
    ax.set_xlabel(f'{ref} x [mm]'); ax.set_ylabel(f'{ref} y [mm]')
    ax.set_title('residual magnitude vs position (pre-fit)')
    ax.set_aspect('equal')

    ax = axes[2]
    ax.axis('off')
    txt = (f'PLANE  {plane}\nREF    {ref}\n\n'
           f'pairs (both clean single)  {len(dx)}\n\n'
           f'fitted rigid transform plane -> ref:\n'
           f'  dx     = {tf.dx:8.2f} mm\n'
           f'  dy     = {tf.dy:8.2f} mm\n'
           f'  theta  = {tf.theta_deg:8.3f} deg\n\n'
           f'residual RMSE\n'
           f'  pre-fit  = {rmse_pre:6.2f} mm\n'
           f'  post-fit = {rmse_post:6.2f} mm')
    ax.text(0.02, 0.98, txt, va='top', ha='left', family='monospace',
            fontsize=11, transform=ax.transAxes)

    fig.suptitle(f'Telescope alignment {plane} -> {ref}', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='P2 telescope mutual alignment.')
    ap.add_argument('run_key', nargs='?', default=sc.DEFAULT_RUN)
    ap.add_argument('--sub-run', default=None,
                    help='sub_run (default: the highest-stats discovered one).')
    ap.add_argument('--cluster-r', type=float, default=15.0,
                    help='pads within this of the leading pad join a cluster '
                         '[mm].')
    ap.add_argument('--min-pairs', type=int, default=100,
                    help='minimum clean-cluster pairs to trust a plane fit.')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True)
    args = ap.parse_args()

    cfg = sc.get_config(args.run_key)
    print(cfg)
    dets = cfg.detectors()
    if len(dets) < 2:
        print('Need >=2 stations to align. Only found:',
              [d.name for d in dets])
        return
    ref = cfg.ref_det()
    print(f'Stations: {[d.name for d in dets]}   reference = {ref.name}')

    sub_run = args.sub_run
    if sub_run is None:
        subs = cfg.find_subruns()
        if not subs:
            print('No sub_runs with combined hits.')
            return
        sub_run = subs[0]
    print(f'sub_run: {sub_run}')
    suffix = cfg.product_suffix(args.veto_sparks)
    out_dir = cfg.out_dir(sc.TELESCOPE_TAG, sub_run, '21_telescope_align')

    # cluster each plane once
    clusters = {}
    for det in dets:
        print(f'  clustering {det.name} ...', flush=True)
        clusters[det.name] = load_clusters(cfg, det, sub_run, args.cluster_r,
                                           args.veto_sparks)
        n_single = int(clusters[det.name]['single'].sum()) if \
            len(clusters[det.name]) else 0
        print(f'    {len(clusters[det.name])} events, {n_single} clean single '
              f'clusters')

    cref = clusters[ref.name]
    cref_s = cref[cref['single']].set_index('eventId')
    planes_json = {ref.name: {'det_tag': ref.det_tag, 'dx': 0.0, 'dy': 0.0,
                              'theta_deg': 0.0, 'rmse_pre': 0.0,
                              'rmse_post': 0.0, 'n': int(len(cref_s)),
                              'is_reference': True}}

    for det in dets:
        if det.name == ref.name:
            continue
        cs = clusters[det.name]
        cs_s = cs[cs['single']].set_index('eventId')
        common = cref_s.index.intersection(cs_s.index)
        n = len(common)
        print(f'  {det.name} -> {ref.name}: {n} clean-single pairs')
        if n < 3:
            print('    too few pairs to fit, skipping plane')
            planes_json[det.name] = {'det_tag': det.det_tag, 'dx': 0.0,
                                     'dy': 0.0, 'theta_deg': 0.0,
                                     'n': int(n), 'fit_ok': False}
            continue
        xr = cref_s.loc[common, 'x'].to_numpy()
        yr = cref_s.loc[common, 'y'].to_numpy()
        xp = cs_s.loc[common, 'x'].to_numpy()
        yp = cs_s.loc[common, 'y'].to_numpy()
        dx, dy = xp - xr, yp - yr
        rmse_pre = float(np.sqrt(np.mean(dx ** 2 + dy ** 2)))
        tf, rmse_post = scl.fit_rigid(xp, yp, xr, yr)
        fx, fy = tf.apply(xp, yp)
        dx2, dy2 = fx - xr, fy - yr
        print(f'    dx={tf.dx:.2f} dy={tf.dy:.2f} theta={tf.theta_deg:.3f} deg '
              f'| RMSE {rmse_pre:.2f} -> {rmse_post:.2f} mm')
        if n < args.min_pairs:
            print(f'    WARNING: only {n} pairs (< --min-pairs '
                  f'{args.min_pairs}); fit is statistics-limited.')
        plot_residuals(dx, dy, dx2, dy2, xr, yr, tf, rmse_pre, rmse_post,
                       det.name, ref.name,
                       os.path.join(out_dir,
                                    f'residuals_{det.det_tag}_to_'
                                    f'{ref.det_tag}{suffix}.png'))
        planes_json[det.name] = {'det_tag': det.det_tag, **tf.to_dict(),
                                 'rmse_pre': rmse_pre, 'rmse_post': rmse_post,
                                 'n': int(n), 'fit_ok': True,
                                 'low_stats': bool(n < args.min_pairs)}

    align = {'run': cfg.RUN, 'sub_run': sub_run, 'ref': ref.name,
             'cluster_r': args.cluster_r, 'veto_sparks': bool(args.veto_sparks),
             'planes': planes_json}
    with open(cfg.alignment_json(sub_run), 'w') as f:
        json.dump(align, f, indent=2)

    # summary scatter of fitted offsets
    fig, ax = plt.subplots(figsize=(6.5, 6))
    for name, p in planes_json.items():
        ax.annotate(name, (p['dx'], p['dy']), fontsize=10,
                    xytext=(4, 4), textcoords='offset points')
        ax.plot(p['dx'], p['dy'], 'o', ms=9,
                color='crimson' if p.get('is_reference') else 'steelblue')
    ax.axhline(0, color='grey', lw=0.6); ax.axvline(0, color='grey', lw=0.6)
    ax.set_xlabel('dx to reference [mm]'); ax.set_ylabel('dy to reference [mm]')
    ax.set_title(f'Fitted plane offsets -> {ref.name}  ({sub_run})')
    ax.grid(True, alpha=0.3); ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'alignment{suffix}.png'), dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    print(f'\nAlignment JSON -> {cfg.alignment_json(sub_run)}')
    print(f'Products       -> {out_dir}')


if __name__ == '__main__':
    main()
