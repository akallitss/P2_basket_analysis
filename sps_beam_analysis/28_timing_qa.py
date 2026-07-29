#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
28_timing_qa.py

Per-station timing QA for the P2 telescope at the SPS, from the analysed
combined-hits (the `time`, `time_of_max`, `max_sample`, `time_over_threshold`
branches that analyze_waveforms already wrote). It needs no decoded waveforms,
so it runs on the partial-fetched combined hits and on a live run.

This is the commissioning-level timing check: is the signal timed into the
sampling window, how wide is the per-hit time distribution, and is there
amplitude slewing (time-walk)? The full waveform-level TOA-algorithm study
(cosmic bench 13_timing_waveforms) needs the decoded_root waveforms and is a
separate, heavier stage.

Per P2 station (uRWELL refs have no pad map and are skipped):
  Per sub_run  (<Analysis>/<det_tag>/<run>/<sub_run>/28_timing_qa/):
    timing_<sub_run>.png   hit-time & max-sample distributions, time-walk
                           (time vs amplitude) with a profile, time vs channel
  Scan level   (<Analysis>/<det_tag>/<run>/scan/28_timing_qa/):
    timing_vs_hv_<det>.png mean/sigma of the hit time, and the walk slope, vs
                           mesh HV -- timing sharpens as the gain rises
    timing_qa_<det>.csv    the per-sub_run numbers

Usage:
  SPS_DATA_ROOT=.../runs SPS_ANALYSIS_ROOT=.../analysis SPS_RUN=<run> \
      python3 28_timing_qa.py live [--sub-run NAME] [--min-amp ADC]
"""

import os
import csv
import argparse

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sps_config as sc
import p2_io as p2io

BR = ['eventId', 'feu', 'channel', 'amplitude', 'time', 'max_sample',
      'time_over_threshold']


def reduce_timing(cfg, det, sub_run, n_samples, sample_period, min_amp):
    """Streamed per-hit timing accumulators for one station in one sub_run."""
    hits_dir = cfg.combined_hits_dir(sub_run)
    feus = set(det.feus)
    t_edges = np.linspace(0, n_samples * sample_period, 121)
    ms_edges = np.linspace(0, n_samples, 4 * n_samples + 1)
    # 2D time-walk accumulate as a hexbin-like sum via sampled arrays (cap size)
    t_samp, a_samp, ch_samp, tot_samp = [], [], [], []
    t_hist = np.zeros(len(t_edges) - 1)
    ms_hist = np.zeros(len(ms_edges) - 1)
    n_hits = 0
    tsum = tsum2 = 0.0
    rng_cap = 400000                       # cap points kept for scatter/profile
    for df in p2io.iter_hits(hits_dir, BR, feus=list(feus), progress=False,
                             min_amp=min_amp):
        if not len(df):
            continue
        t = df['time'].to_numpy(np.float64)
        a = df['amplitude'].to_numpy(np.float64)
        ok = np.isfinite(t) & (t >= t_edges[0]) & (t <= t_edges[-1])
        t, a = t[ok], a[ok]
        ch = df['channel'].to_numpy()[ok]
        tot = df['time_over_threshold'].to_numpy(np.float64)[ok]
        ms = df['max_sample'].to_numpy(np.float64)[ok]
        n_hits += len(t)
        tsum += t.sum(); tsum2 += (t * t).sum()
        t_hist += np.histogram(t, bins=t_edges)[0]
        ms_hist += np.histogram(np.clip(ms, ms_edges[0], ms_edges[-1]),
                                bins=ms_edges)[0]
        if sum(len(x) for x in t_samp) < rng_cap:
            t_samp.append(t); a_samp.append(a); ch_samp.append(ch)
            tot_samp.append(tot)
    if not n_hits:
        return None
    samp = dict(
        t=np.concatenate(t_samp), a=np.concatenate(a_samp),
        ch=np.concatenate(ch_samp), tot=np.concatenate(tot_samp))
    mean = tsum / n_hits
    sigma = float(np.sqrt(max(tsum2 / n_hits - mean ** 2, 0.0)))
    return dict(n_hits=n_hits, t_edges=t_edges, t_hist=t_hist,
                ms_edges=ms_edges, ms_hist=ms_hist, mean=float(mean),
                sigma=sigma, samp=samp)


def walk_profile(a, t, nbins=25):
    """Median time in amplitude bins (log-spaced) -> slewing curve + slope."""
    good = (a > 0) & np.isfinite(t)
    a, t = a[good], t[good]
    if len(a) < 50:
        return None
    edges = np.logspace(np.log10(max(a.min(), 1)), np.log10(a.max()), nbins + 1)
    idx = np.clip(np.digitize(a, edges) - 1, 0, nbins - 1)
    cx, cy = [], []
    for b in range(nbins):
        m = idx == b
        if m.sum() > 20:
            cx.append(np.sqrt(edges[b] * edges[b + 1]))
            cy.append(np.median(t[m]))
    if len(cx) < 3:
        return None
    cx, cy = np.array(cx), np.array(cy)
    slope = np.polyfit(np.log10(cx), cy, 1)[0]     # ns per decade of amplitude
    return cx, cy, float(slope)


def plot_subrun(res, det, sub_run, mesh_v, n_samples, sample_period, out_png):
    s = res['samp']
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    a0 = ax[0, 0]
    c = 0.5 * (res['t_edges'][:-1] + res['t_edges'][1:])
    a0.step(c, res['t_hist'], where='mid', color='navy')
    a0.axvline(res['mean'], color='crimson', ls='--',
               label=f'mean {res["mean"]:.0f} ns, sigma {res["sigma"]:.0f} ns')
    a0.set_xlabel('hit time in window [ns]'); a0.set_ylabel('hits')
    a0.legend(fontsize=8); a0.grid(alpha=0.3)
    a0.set_title('hit-time distribution')

    a1 = ax[0, 1]
    mc = 0.5 * (res['ms_edges'][:-1] + res['ms_edges'][1:])
    a1.step(mc, res['ms_hist'], where='mid', color='teal')
    a1.axvspan(0, 2, color='r', alpha=0.08); a1.axvspan(n_samples - 2, n_samples,
                                                        color='r', alpha=0.08)
    a1.set_xlim(0, n_samples)
    a1.set_xlabel('sample of waveform max'); a1.set_ylabel('hits')
    a1.grid(alpha=0.3); a1.set_title('peak sample (latency placement)')

    a2 = ax[1, 0]
    a2.scatter(s['a'], s['t'], s=3, alpha=0.15, c='slategrey', linewidths=0)
    wp = walk_profile(s['a'], s['t'])
    if wp:
        cx, cy, slope = wp
        a2.plot(cx, cy, 'o-', color='crimson', ms=4,
                label=f'median profile ({slope:+.1f} ns/decade)')
        a2.legend(fontsize=8)
    a2.set_xscale('log'); a2.set_xlabel('hit amplitude [ADC] (log)')
    a2.set_ylabel('hit time [ns]'); a2.grid(alpha=0.3)
    a2.set_title('time-walk (slewing)')

    a3 = ax[1, 1]
    a3.scatter(s['ch'], s['t'], s=3, alpha=0.15, c='slateblue', linewidths=0)
    a3.set_xlabel('FEU channel'); a3.set_ylabel('hit time [ns]')
    a3.grid(alpha=0.3); a3.set_title('time vs channel (per-channel offsets)')

    fig.suptitle(f'{det.name} timing — {sub_run}'
                 + (f'  (mesh {mesh_v} V)' if mesh_v is not None else ''),
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=140, bbox_inches='tight'); plt.close(fig)


def plot_scan(rows, det, out_png):
    r = [x for x in rows if x['mesh_v'] is not None]
    if len(r) < 2:
        return False
    r.sort(key=lambda x: x['mesh_v'])
    v = np.array([x['mesh_v'] for x in r], float)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for a, key, ylab in ((ax[0], 'mean_ns', 'mean hit time [ns]'),
                         (ax[1], 'sigma_ns', 'hit-time sigma [ns]'),
                         (ax[2], 'walk_slope', 'time-walk [ns/decade]')):
        y = np.array([(x[key] if x[key] is not None else np.nan) for x in r],
                     float)
        a.plot(v, y, 'o-', lw=1.8, ms=7)
        a.set_xlabel('mesh HV [V]'); a.set_ylabel(ylab); a.grid(alpha=0.3)
    fig.suptitle(f'{det.name} — timing vs mesh HV', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_png, dpi=150, bbox_inches='tight'); plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_key', nargs='?', default=None)
    ap.add_argument('--sub-run', default=None)
    ap.add_argument('--min-amp', type=float, default=0.0)
    ap.add_argument('--scan-only', action='store_true',
                    help='skip the per-sub_run pass and build ONLY the '
                         'scan-level curve, from the scan_row.json files that '
                         'earlier per-sub_run runs left behind. This is the '
                         'merge step after an HTCondor sweep, where each '
                         'sub_run was processed on a different machine.')
    args = ap.parse_args()

    cfg = sc.get_config(args.run_key)
    print(cfg)
    daq = cfg.daq_info()
    n_samples = int(daq.get('n_samples_per_waveform') or 16)
    sample_period = float(daq.get('sample_period') or 60)
    dets = cfg.mappable_detectors()
    skipped = [d.name for d in cfg.detectors() if d not in dets]
    if skipped:
        print(f'  (no pad map, skipped: {", ".join(skipped)})')
    subs = [args.sub_run] if args.sub_run else cfg.find_subruns()
    if args.scan_only:
        # The merge runs where the DATA is not: take the sub_run list from the
        # analysis tree instead of from hits files on disk.
        subs = cfg.scan_row_subruns('28_timing_qa') or subs
        if not subs:
            raise SystemExit('No persisted scan rows found — run the '
                             'per-sub_run pass first.')
    elif not subs:
        raise SystemExit('No sub_runs with combined hits on disk.')

    for det in dets:
        print(f'\n== {det.name}  FEUs {det.feus}')
        rows = []
        for sub in (() if args.scan_only else subs):
            res = reduce_timing(cfg, det, sub, n_samples, sample_period,
                                args.min_amp)
            if res is None:
                print(f'  {sub}: no hits'); continue
            mesh_v = cfg.subrun_mesh_hv(sub, det)
            out = cfg.out_dir(det.det_tag, sub, '28_timing_qa')
            plot_subrun(res, det, sub, mesh_v, n_samples, sample_period,
                        os.path.join(out, f'timing_{sub}.png'))
            wp = walk_profile(res['samp']['a'], res['samp']['t'])
            row = dict(sub_run=sub, mesh_v=mesh_v, n_hits=res['n_hits'],
                       mean_ns=round(res['mean'], 2),
                       sigma_ns=round(res['sigma'], 2),
                       walk_slope=(round(wp[2], 2) if wp else None))
            rows.append(row)
            # persist for a later --scan-only merge (HTCondor sweeps)
            cfg.save_scan_row(det.det_tag, sub, '28_timing_qa', row)
            print(f'  {sub}: mesh {mesh_v} V  {res["n_hits"]} hits  '
                  f'time {res["mean"]:.0f}+-{res["sigma"]:.0f} ns'
                  + (f'  walk {wp[2]:+.1f} ns/dec' if wp else ''))

        if args.scan_only:
            rows = cfg.load_scan_rows(det.det_tag, '28_timing_qa')
            print(f'  merged {len(rows)} persisted sub_run row(s)')

        if rows and (args.scan_only or len(subs) > 1):
            scan_out = cfg.out_dir(det.det_tag, 'scan', '28_timing_qa')
            if plot_scan(rows, det,
                         os.path.join(scan_out, f'timing_vs_hv_{det.det_tag}.png')):
                print(f'  -> timing-vs-HV scan: {scan_out}')
            with open(os.path.join(scan_out, f'timing_qa_{det.det_tag}.csv'),
                      'w', newline='') as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(sorted(rows, key=lambda x: (x['mesh_v'] or 0)))


if __name__ == '__main__':
    main()
