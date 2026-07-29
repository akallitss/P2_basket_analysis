#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
26_hv_spark_qa.py

HV spark QA for the P2 telescope at the SPS, the beam-arm counterpart of the
cosmic bench 07_hv_spark_qa.py. Sparks are micro-discharges across the
amplification mesh: brief imon spikes on a station's mesh HV channel. This stage
reads only hv_monitor.csv (already partial-fetched with every sub_run), so it
runs with no extra data and on a live run.

It works per P2 station (the uRWELL references have no mesh channel and are
skipped) and produces, for each station:

  Per sub_run  (<Analysis>/<det_tag>/<run>/<sub_run>/26_hv_spark_qa/):
    hv_timeline_<sub_run>.png   vmon + imon vs time with the spark windows shaded

  Scan level   (<Analysis>/<det_tag>/<run>/scan/26_hv_spark_qa/):
    spark_vs_hv_<det>.png       spark rate, peak imon, total charge and live
                                fraction vs mesh HV -- for a mesh/gain scan this
                                is the sparking turn-on, the HV-headroom curve
                                that says how hard you can push the gain.
    spark_qa_<det>.csv          the per-sub_run numbers

Usage:
  SPS_DATA_ROOT=.../runs SPS_ANALYSIS_ROOT=.../analysis SPS_RUN=<run> \
      python3 26_hv_spark_qa.py live [--sub-run NAME] [--i-thr uA]
"""

import os
import csv
import argparse

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

import sps_config as sc
import p2_sparks as ps


def spark_shim(cfg, det, i_thr=None):
    """The cfg surface SparkVeto.from_csv expects, for one station's mesh."""
    class _Shim:
        SPARK_CHANNEL = det.spark_channel
        SPARK_IMON_THR = i_thr if i_thr is not None else cfg.SPARK_IMON_THR
        SPARK_GUARD_BEFORE = cfg.SPARK_GUARD_BEFORE
        SPARK_GUARD_AFTER = cfg.SPARK_GUARD_AFTER
        BURST_NPADS = cfg.BURST_NPADS
    return _Shim


def shade(ax, intervals, t_scale=1 / 60.0, label='spark window'):
    for k, (lo, hi) in enumerate(intervals):
        ax.axvspan(lo * t_scale, hi * t_scale, color='crimson', alpha=0.18,
                   lw=0, label=label if k == 0 else None)


def plot_timeline(sv, det, sub_run, mesh_v, out_png):
    hv = sv.hv
    tmin = hv['t'].to_numpy() / 60.0
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)
    a0.plot(tmin, hv['vmon'], lw=0.8, color='navy')
    a0.set_ylabel('vmon [V]')
    a0.set_title(f'{det.name} mesh HV (ch {sv.channel}, set {mesh_v} V) — '
                 f'{sub_run}')
    shade(a0, sv.intervals)
    if a0.get_legend_handles_labels()[0]:
        a0.legend(loc='lower left', fontsize=8)
    a0.grid(True, alpha=0.3)

    a1.plot(tmin, hv['imon'], lw=0.7, color='darkred')
    a1.axhline(sv.i_thr, color='k', ls='--', lw=1,
               label=f'spark threshold {sv.i_thr:g} µA')
    shade(a1, sv.intervals)
    a1.set_ylabel('imon [µA]'); a1.set_xlabel('time since sub_run start [min]')
    a1.legend(loc='upper right', fontsize=8); a1.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140, bbox_inches='tight'); plt.close(fig)


def subrun_metrics(sv, mesh_v):
    n = len(sv.sparks)
    dt_min = sv.t_run / 60.0 if sv.t_run > 0 else np.nan
    return dict(
        mesh_v=mesh_v,
        duration_s=round(sv.t_run, 1),
        n_sparks=int(n),
        spark_rate_per_min=float(n / dt_min) if dt_min and dt_min > 0 else None,
        peak_imon_uA=(float(sv.sparks['peak_imon'].max()) if n else
                      float(sv.hv['imon'].max()) if len(sv.hv) else None),
        mean_imon_uA=float(sv.hv['imon'].mean()) if len(sv.hv) else None,
        total_charge_uC=float(sv.sparks['charge'].sum()) if n else 0.0,
        live_fraction=float(sv.live_fraction()),
    )


def plot_scan(rows, det, out_png, xlabel='mesh HV [V]'):
    r = [x for x in rows if x['mesh_v'] is not None]
    if len(r) < 2:
        return False
    r.sort(key=lambda x: x['mesh_v'])
    v = np.array([x['mesh_v'] for x in r], float)
    fig, ax = plt.subplots(2, 2, figsize=(13, 8))
    panels = [
        ('spark_rate_per_min', 'spark rate [/min]', 'steelblue'),
        ('peak_imon_uA', 'peak imon [µA]', 'crimson'),
        ('total_charge_uC', 'total spark charge [µC]', 'seagreen'),
        ('live_fraction', 'live fraction (1 = no veto)', 'darkorange'),
    ]
    for a, (key, ylab, c) in zip(ax.ravel(), panels):
        y = np.array([(x[key] if x[key] is not None else np.nan) for x in r],
                     float)
        a.plot(v, y, 'o-', lw=1.8, ms=7, color=c)
        a.set_xlabel(xlabel); a.set_ylabel(ylab)
        a.grid(True, alpha=0.3)
        if key == 'live_fraction':
            a.set_ylim(0, 1.02)
        else:
            a.set_ylim(bottom=0)
    fig.suptitle(f'{det.name} — sparking vs {xlabel} (the HV headroom)',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=150, bbox_inches='tight'); plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_key', nargs='?', default=None)
    ap.add_argument('--sub-run', default=None,
                    help='only this sub_run (default: all on disk).')
    ap.add_argument('--i-thr', type=float, default=None,
                    help='override spark imon threshold [µA].')
    ap.add_argument('--scan-only', action='store_true',
                    help='build ONLY the scan-level curve, from the '
                         'scan_row.json files earlier per-sub_run runs left '
                         'behind (the merge step after an HTCondor sweep).')
    args = ap.parse_args()

    cfg = sc.get_config(args.run_key)
    print(cfg)
    dets = [d for d in cfg.detectors() if d.spark_channel]
    skipped = [d.name for d in cfg.detectors() if not d.spark_channel]
    if skipped:
        print(f'  (no mesh channel, skipped: {", ".join(skipped)})')
    if not dets:
        raise SystemExit('No station with a mesh HV channel.')

    subs = [args.sub_run] if args.sub_run else cfg.find_subruns()
    if args.scan_only:
        # The merge runs where the DATA is not: take the sub_run list from the
        # analysis tree instead of from hits files on disk.
        subs = cfg.scan_row_subruns('26_hv_spark_qa') or subs
        if not subs:
            raise SystemExit('No persisted scan rows found — run the '
                             'per-sub_run pass first.')
    elif not subs:
        raise SystemExit('No sub_runs with combined hits on disk.')

    for det in dets:
        scan_ax, scan_label = cfg.scan_axis(subs, det)
        print(f'\n== {det.name}  (mesh channel {det.spark_channel}, '
              f'scan axis: {scan_label})')
        rows = []
        for sub in (() if args.scan_only else subs):
            hv_csv = cfg.hv_monitor_csv(sub)
            if not os.path.isfile(hv_csv):
                print(f'  {sub}: no hv_monitor.csv, skipping')
                continue
            sv = ps.SparkVeto.from_csv(hv_csv, spark_shim(cfg, det, args.i_thr))
            mesh_v = (cfg.subrun_scan_hv(sub, det, scan_ax) if scan_ax
                      else cfg.subrun_mesh_hv(sub, det))
            out = cfg.out_dir(det.det_tag, sub, '26_hv_spark_qa')
            plot_timeline(sv, det, sub, mesh_v,
                          os.path.join(out, f'hv_timeline_{sub}.png'))
            m = subrun_metrics(sv, mesh_v)
            m['sub_run'] = sub
            rows.append(m)
            # persist for a later --scan-only merge (HTCondor sweeps)
            cfg.save_scan_row(det.det_tag, sub, '26_hv_spark_qa', m)
            print(f'  {sub}: mesh {mesh_v} V  {m["n_sparks"]} sparks  '
                  f'peak imon {m["peak_imon_uA"]:.2f} µA  '
                  f'live {m["live_fraction"]:.3f}')

        if args.scan_only:
            rows = cfg.load_scan_rows(det.det_tag, '26_hv_spark_qa')
            print(f'  merged {len(rows)} persisted sub_run row(s)')

        if rows and (args.scan_only or len(subs) > 1):
            scan_out = cfg.out_dir(det.det_tag, 'scan', '26_hv_spark_qa')
            if plot_scan(rows, det,
                         os.path.join(scan_out,
                                      f'spark_vs_hv_{det.det_tag}.png'),
                         xlabel=scan_label):
                print(f'  -> spark-vs-HV scan: {scan_out}')
            keys = ['sub_run', 'mesh_v', 'duration_s', 'n_sparks',
                    'spark_rate_per_min', 'peak_imon_uA', 'mean_imon_uA',
                    'total_charge_uC', 'live_fraction']
            with open(os.path.join(scan_out, f'spark_qa_{det.det_tag}.csv'),
                      'w', newline='') as fh:
                w = csv.DictWriter(fh, fieldnames=keys, extrasaction='ignore')
                w.writeheader()
                w.writerows(sorted(rows, key=lambda x: (x['mesh_v'] or 0)))


if __name__ == '__main__':
    main()
