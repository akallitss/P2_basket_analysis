#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07_hv_spark_qa.py

HV spark QA for the P2 mesh: characterise the sparking behaviour over the run
from hv_monitor.csv and validate the veto (p2_sparks) that the other stages use.

Sparks are micro-discharges across the amplification mesh; they appear as brief
imon spikes on the mesh HV channel (cfg.SPARK_CHANNEL). This stage detects them
(p2_sparks.SparkVeto), shows how the sparking rate evolves during the run, and
cross-checks the HV-based spark times against DAQ high-multiplicity bursts (a
spark dumps charge onto many pads at once), i.e. the two independent spark
signatures should line up in time.

Products (written to <Analysis>/<detN>/<run>/<sub_run>/07_hv_spark_qa/):
  hv_current_vmon_timeline.png  vmon & imon vs time, spark windows shaded
  spark_rate_evolution.png      sparks/hour (sliding) + cumulative count
  spark_amplitude_dist.png      peak-current and per-spark charge distributions
  spark_daq_crosscheck.png      DAQ pad-multiplicity vs time vs HV spark times
  spark_qa_summary.txt          spark list + summary (rate, charge, deadtime)

Usage: python3 07_hv_spark_qa.py [run_key] [--i-thr 2.0] [--rate-bin-min 30]
"""

import os
import glob
import argparse
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot

import p2_qa_config as qa
import p2_mapping as pmap
import p2_sparks as ps

_BRANCHES = ['eventId', 'trigger_timestamp_ns', 'channel', 'feu']


# --------------------------------------------------------------------------- #
def _shade_sparks(ax, intervals, t_scale=1 / 3600.0, label='spark veto'):
    for k, (lo, hi) in enumerate(intervals):
        ax.axvspan(lo * t_scale, hi * t_scale, color='crimson', alpha=0.18,
                   lw=0, label=label if k == 0 else None)


def plot_timeline(sv, cfg, out_dir, suffix):
    hv = sv.hv
    th = hv['t'].to_numpy() / 3600.0
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(13, 7), sharex=True)
    a0.plot(th, hv['vmon'], lw=0.8, color='navy')
    a0.set_ylabel('vmon [V]')
    a0.set_title(f'{cfg.DET_NAME} mesh HV (channel {sv.channel}) over the run '
                 f'— {cfg.RUN}/{cfg.SUB_RUN}')
    _shade_sparks(a0, sv.intervals)
    a0.legend(loc='lower left', fontsize=8); a0.grid(True, alpha=0.3)

    a1.plot(th, hv['imon'], lw=0.7, color='darkred')
    a1.axhline(sv.i_thr, color='k', ls='--', lw=1,
               label=f'spark threshold {sv.i_thr:g} µA')
    _shade_sparks(a1, sv.intervals)
    a1.set_ylabel('imon [µA]'); a1.set_xlabel('time since run start [h]')
    a1.legend(loc='upper right', fontsize=8); a1.grid(True, alpha=0.3)
    fig.tight_layout()
    p = f'{out_dir}/hv_current_vmon_timeline{suffix}.png'
    fig.savefig(p, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f'  saved {p}')


def plot_rate(sv, cfg, out_dir, suffix, bin_min=30.0):
    t_run_h = sv.t_run / 3600.0
    edges = np.arange(0, t_run_h + bin_min / 60.0, bin_min / 60.0)
    tp = sv.sparks['t_peak'].to_numpy() / 3600.0 if len(sv.sparks) else np.array([])
    counts, _ = np.histogram(tp, bins=edges)
    rate = counts / (bin_min / 60.0)              # sparks per hour
    ctr = 0.5 * (edges[:-1] + edges[1:])

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(ctr, rate, width=bin_min / 60.0 * 0.9, color='steelblue',
           label=f'spark rate ({bin_min:g}-min bins)')
    ax.set_xlabel('time since run start [h]'); ax.set_ylabel('sparks / hour')
    ax.set_title(f'{cfg.DET_NAME} spark-rate evolution — {cfg.RUN}/{cfg.SUB_RUN}')
    ax.grid(True, alpha=0.3, axis='y')

    axc = ax.twinx()
    tps = np.sort(tp)
    axc.plot(tps, np.arange(1, len(tps) + 1), color='crimson', lw=1.8,
             label='cumulative sparks')
    axc.set_ylabel('cumulative # sparks', color='crimson')
    axc.tick_params(axis='y', labelcolor='crimson')
    l0, la0 = ax.get_legend_handles_labels()
    l1, la1 = axc.get_legend_handles_labels()
    ax.legend(l0 + l1, la0 + la1, loc='upper left', fontsize=8)
    fig.tight_layout()
    p = f'{out_dir}/spark_rate_evolution{suffix}.png'
    fig.savefig(p, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f'  saved {p}')


def plot_amplitude(sv, cfg, out_dir, suffix):
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 4.5))
    if len(sv.sparks):
        a0.hist(sv.sparks['peak_imon'], bins=25, color='indianred')
        a1.hist(sv.sparks['charge'], bins=np.logspace(
            np.log10(max(sv.sparks['charge'].min(), 1e-2)),
            np.log10(sv.sparks['charge'].max() + 1), 25), color='seagreen')
        a1.set_xscale('log')
    a0.set_xlabel('peak imon per spark [µA]'); a0.set_ylabel('sparks')
    a0.set_title('Spark peak current'); a0.grid(True, alpha=0.3)
    a1.set_xlabel('charge per spark [µC = µA·s]'); a1.set_ylabel('sparks')
    a1.set_title('Spark charge (log x)'); a1.grid(True, alpha=0.3)
    fig.suptitle(f'{cfg.DET_NAME} spark amplitude — {cfg.RUN}/{cfg.SUB_RUN}')
    fig.tight_layout()
    p = f'{out_dir}/spark_amplitude_dist{suffix}.png'
    fig.savefig(p, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f'  saved {p}')


def load_event_multiplicity(cfg, channel_table):
    """Per-event pad multiplicity and event time [s] from the combined hits.
    Streamed chunk by chunk (p2_io) — the old concat of every chunk was the
    OOM that froze the machine on the det4 run."""
    import p2_io as p2io
    parts = []
    for df in p2io.iter_hits(cfg.combined_hits_dir, _BRANCHES,
                             channel_table.attrs['feus'],
                             t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H, min_amp=cfg.MIN_AMP):
        df = pmap.attach_pads_to_hits(df, channel_table)
        df = df[df['mapped'] & df['pad_cx'].notna()]
        if not len(df):
            continue
        g = df.groupby('eventId')
        parts.append(pd.DataFrame({
            't': g['trigger_timestamp_ns'].first().to_numpy() / 1e9,
            'npad': g['channel_id'].nunique().to_numpy(),
        }))
    return pd.concat(parts, ignore_index=True)


def plot_crosscheck(sv, ev, cfg, out_dir, suffix):
    fig, ax = plt.subplots(figsize=(13, 5))
    th = ev['t'].to_numpy() / 3600.0
    ax.scatter(th, ev['npad'], s=4, alpha=0.25, c='slategrey', linewidths=0,
               label='event pad multiplicity')
    _shade_sparks(ax, sv.intervals, label='HV spark window')
    ax.set_yscale('log')
    ax.set_xlabel('time since run start [h]')
    ax.set_ylabel('pads firing per event (log)')
    ax.set_title(f'{cfg.DET_NAME} DAQ multiplicity vs HV sparks '
                 f'(bursts should sit under red bands) — {cfg.RUN}/{cfg.SUB_RUN}')
    ax.legend(loc='upper right', fontsize=8); ax.grid(True, alpha=0.3)

    # how much of the high-multiplicity tail is inside a spark window?
    hi = ev[ev['npad'] >= 10]
    if len(hi):
        in_win = ~sv.event_mask((hi['t'].to_numpy() * 1e9))
        frac = in_win.mean() * 100
        ax.text(0.01, 0.02,
                f'{in_win.sum()}/{len(hi)} events with ≥10 pads fall in a spark '
                f'window ({frac:.1f}%)', transform=ax.transAxes, fontsize=9,
                bbox=dict(boxstyle='round', fc='white', alpha=0.8))
    fig.tight_layout()
    p = f'{out_dir}/spark_daq_crosscheck{suffix}.png'
    fig.savefig(p, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f'  saved {p}')
    return hi, sv


def main():
    ap = argparse.ArgumentParser(description='P2 HV spark QA.')
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'])
    ap.add_argument('--i-thr', type=float, default=None,
                    help='override spark imon threshold [µA] (default cfg value).')
    ap.add_argument('--rate-bin-min', type=float, default=30.0,
                    help='spark-rate time bin [min] (default 30).')
    ap.add_argument('--no-daq-crosscheck', action='store_true',
                    help='skip loading the hits for the DAQ multiplicity check.')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    if args.i_thr is not None:
        cfg.SPARK_IMON_THR = args.i_thr
    out_dir = cfg.out_dir('07_hv_spark_qa')
    suffix = cfg.dead_suffix   # this stage is about sparks; note dead connectors only

    sv = ps.SparkVeto.from_cfg(cfg)
    print(sv.summary())

    plot_timeline(sv, cfg, out_dir, suffix)
    plot_rate(sv, cfg, out_dir, suffix, args.rate_bin_min)
    plot_amplitude(sv, cfg, out_dir, suffix)

    summary = [sv.summary(), '']
    if not args.no_daq_crosscheck:
        ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                      det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                      strategy=args.strategy,
                                      drop_connectors=cfg.DEAD_CONNECTORS)
        ev = load_event_multiplicity(cfg, ct)
        hi, _ = plot_crosscheck(sv, ev, cfg, out_dir, suffix)
        n_hi = len(hi)
        n_hi_spark = int((~sv.event_mask(hi['t'].to_numpy() * 1e9)).sum()) if n_hi else 0
        summary.append('DAQ cross-check:')
        summary.append(f'  events with >=10 pads         : {n_hi}')
        summary.append(f'  of those inside a spark window : {n_hi_spark} '
                       f'({(100*n_hi_spark/n_hi if n_hi else 0):.1f}%)')
        summary.append('')

    # full spark list
    summary.append(f'Spark list (channel {sv.channel}, imon>={sv.i_thr:g} µA):')
    summary.append('   t_peak[s]  dur[s]  peak[µA]  charge[µC]  n')
    for _, r in sv.sparks.sort_values('t_peak').iterrows():
        summary.append(f'  {r.t_peak:9.0f}  {r.dur:6.0f}  {r.peak_imon:8.2f}  '
                       f'{r.charge:9.1f}  {int(r.n_samples)}')

    txt = '\n'.join(summary)
    with open(f'{out_dir}/spark_qa_summary{suffix}.txt', 'w') as f:
        f.write(txt + '\n')
    print(f'  saved {out_dir}/spark_qa_summary{suffix}.txt')
    print('Done.')


if __name__ == '__main__':
    main()
