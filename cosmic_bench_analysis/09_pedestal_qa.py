#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09_pedestal_qa.py

Per-channel pedestal QA for the P2 pad detector, adapted from
nTof_x17/mx_june_cosmic_qa/11_pedestal_qa.py.

The DREAM processing keeps a pulse only if amplitude >= THRESHOLD_SIGMA *
pedestalRMS_channel (THRESHOLD_SIGMA = 5). So each channel's hit threshold is set
by its own pedestal RMS. This matters two ways:

  * a channel with an anomalously LOW RMS gets a LOW threshold -> can over-trigger
    on baseline noise -> excess (low-amplitude) hits: a genuine noisy channel;
  * a channel with an inflated RMS (e.g. discharge pickup during the pedestal
    run) gets a HIGH threshold -> real signal is suppressed -> a cold pad.

This stage reads the per-channel pedestal mean / RMS that the processing actually
applied (the `pedestals` TTree in the per-FEU hits_root files -- which the
combined-hits files drop), maps it onto the pad plane, and checks whether the
pedestal explains any of the per-pad hit-rate structure. Crucially it plots
pedestal RMS vs firing rate: if noisy channels were inflating the rate this would
correlate; on P2 det1 it does not (r~0.05), which is how we know the hot pads are
discharges, not noise.

Products (written to <Analysis>/<detN>/<run>/<sub_run>/09_pedestal_qa/):
  pedestal_rms_map.png       pedestal RMS on the pad plane
  pedestal_mean_map.png      pedestal baseline (mean) on the pad plane
  pedestal_rms_dist.png      RMS distribution + nominal threshold + outliers
  pedestal_rms_vs_rate.png   RMS (threshold) vs per-pad firing rate
  pedestal_qa_summary.txt / pedestal_outliers.csv

Usage: python3 09_pedestal_qa.py [run_key] [--nsigma 5] [--veto-sparks]
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
import p2_pad_sparks as pps
import p2_io as p2io

THRESHOLD_SIGMA = 5.0   # DREAM WaveformAnalyzer nominal hit threshold = 5*RMS


def load_pedestals(hits_root_dir, channel_table):
    """Per-channel pedestal (mean, rms) mapped onto connected P2 pads.

    The pedestal is per-FEU and (per nTof) identical across file indices, so we
    read the first hits_root file for each FEU. Returns the channel table
    restricted to mapped pads with ped_mean / ped_rms columns.
    """
    feus = channel_table.attrs['feus']
    peds = []
    for feu in feus:
        files = sorted(glob.glob(os.path.join(hits_root_dir, f'*_{feu:02d}_hits.root')))
        if not files:
            raise FileNotFoundError(
                f'No per-FEU hits_root file for FEU {feu} in {hits_root_dir}. '
                'Fetch the hits_root/ directory from the DAQ host.')
        a = uproot.open(files[0])['pedestals'].arrays(
            ['channel', 'mean', 'rms'], library='np')
        peds.append(pd.DataFrame({'feu': feu, 'channel': a['channel'],
                                  'ped_mean': a['mean'], 'ped_rms': a['rms']}))
    ped = pd.concat(peds, ignore_index=True)
    m = channel_table.merge(ped, on=['feu', 'channel'], how='left')
    return m[m['mapped'] & m['ped_rms'].notna()].copy()


def load_firing(cfg, ct, veto_sparks=True):
    """Per-pad firing fraction, streamed chunk by chunk (p2_io) so the full
    hit table never sits in memory."""
    br = ['eventId', 'channel', 'feu'] + (
        ['trigger_timestamp_ns'] if veto_sparks else [])
    sv = ps.SparkVeto.from_cfg(cfg) if veto_sparks else None
    fire = None
    ev_parts = []
    for a in p2io.iter_hits(cfg.combined_hits_dir, br, ct.attrs['feus'],
                            t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H, min_amp=cfg.MIN_AMP):
        if sv is not None:
            a, _ = sv.apply(a)
        h = pmap.attach_pads_to_hits(a, ct)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        if not len(h):
            continue
        ev_parts.append(h['eventId'].unique())
        c = h.groupby('channel_id')['eventId'].nunique()
        fire = c if fire is None else fire.add(c, fill_value=0)
    nev = len(np.unique(np.concatenate(ev_parts))) if ev_parts else 0
    return (fire / max(nev, 1)).rename('fire'), nev


def _pad_map(ax, df, val, label, cmap='viridis'):
    sc = ax.scatter(df['pad_cx'], df['pad_cy'], c=df[val], s=15, marker='s',
                    cmap=cmap, linewidths=0)
    ax.set_aspect('equal'); ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
    plt.colorbar(sc, ax=ax, label=label, fraction=0.046, pad=0.04)


def main():
    ap = argparse.ArgumentParser(description='P2 per-channel pedestal QA.')
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'])
    ap.add_argument('--nsigma', type=float, default=THRESHOLD_SIGMA,
                    help='threshold = nsigma * pedestal RMS (default 5).')
    ap.add_argument('--z-rms', type=float, default=5.0,
                    help='robust-z on pedestal RMS above which a channel is '
                         'flagged pedestal-anomalous (default 5).')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True)
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    out_dir = cfg.out_dir('09_pedestal_qa')
    suffix = cfg.dead_suffix

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy=args.strategy,
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    ped = load_pedestals(cfg.hits_root_dir, ct).set_index('channel_id')
    fire, nev = load_firing(cfg, ct, args.veto_sparks)
    ped = ped.join(fire); ped['fire'] = ped['fire'].fillna(0.0)
    ped['thr'] = args.nsigma * ped['ped_rms']
    ped['rms_z'] = pps.robust_z(ped['ped_rms'])
    ped['ped_noisy'] = ped['rms_z'] < -3    # anomalously low RMS -> low threshold
    ped['ped_hot_rms'] = ped['rms_z'] > args.z_rms   # inflated RMS

    r = float(ped[['ped_rms', 'fire']].corr().iloc[0, 1])
    print(f'connected pads with pedestal: {len(ped)}')
    print(f'ped_rms median {ped.ped_rms.median():.2f}, range '
          f'[{ped.ped_rms.min():.2f}, {ped.ped_rms.max():.2f}] ADC')
    print(f'nominal threshold ({args.nsigma:g}*RMS): median {ped.thr.median():.1f}, '
          f'max {ped.thr.max():.1f} ADC')
    print(f'corr(ped_rms, firing rate) = {r:.3f}')

    # -- plots --
    fig, ax = plt.subplots(figsize=(9, 8))
    _pad_map(ax, ped, 'ped_rms', 'pedestal RMS [ADC]', cmap='plasma')
    ax.set_title(f'{cfg.DET_NAME} pedestal RMS on the pad plane — {cfg.RUN}/{cfg.SUB_RUN}')
    fig.tight_layout(); fig.savefig(f'{out_dir}/pedestal_rms_map{suffix}.png',
                                    dpi=150, bbox_inches='tight'); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 8))
    _pad_map(ax, ped, 'ped_mean', 'pedestal baseline (mean) [ADC]', cmap='viridis')
    ax.set_title(f'{cfg.DET_NAME} pedestal baseline on the pad plane — {cfg.RUN}/{cfg.SUB_RUN}')
    fig.tight_layout(); fig.savefig(f'{out_dir}/pedestal_mean_map{suffix}.png',
                                    dpi=150, bbox_inches='tight'); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(ped['ped_rms'], bins=60, color='teal')
    ax.axvline(ped.ped_rms.median(), color='k', ls='--',
               label=f'median {ped.ped_rms.median():.2f}')
    ax.set_yscale('log'); ax.set_xlabel('pedestal RMS [ADC]'); ax.set_ylabel('channels')
    ax.set_title(f'{cfg.DET_NAME} pedestal RMS distribution — sets the {args.nsigma:g}σ '
                 f'hit threshold ({args.nsigma*ped.ped_rms.median():.0f} ADC median)')
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(f'{out_dir}/pedestal_rms_dist{suffix}.png',
                                    dpi=150, bbox_inches='tight'); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.scatter(ped['ped_rms'], ped['fire'], s=10, c='steelblue', alpha=0.5,
               linewidths=0, label='pads')
    hot = ped[ped['ped_hot_rms']]
    if len(hot):
        ax.scatter(hot['ped_rms'], hot['fire'], s=40, c='crimson', edgecolors='k',
                   linewidths=0.5, label=f'inflated RMS (z>{args.z_rms:g})')
    ax.set_xlabel('pedestal RMS [ADC]  (∝ hit threshold)')
    ax.set_ylabel('per-pad firing fraction')
    ax.set_title(f'{cfg.DET_NAME} pedestal RMS vs hit rate  (r={r:.3f}: noise does '
                 f'NOT drive rate)\n{cfg.RUN}/{cfg.SUB_RUN}')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(f'{out_dir}/pedestal_rms_vs_rate{suffix}.png',
                                    dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f'  saved 4 plots to {out_dir}')

    # outliers + summary
    out = ped[ped['ped_noisy'] | ped['ped_hot_rms']].sort_values('ped_rms')
    out.reset_index()[['channel_id', 'connector_N', 'radius', 'pad_cx', 'pad_cy',
                       'ped_mean', 'ped_rms', 'thr', 'rms_z', 'fire',
                       'ped_noisy', 'ped_hot_rms']].to_csv(
        f'{out_dir}/pedestal_outliers{suffix}.csv', index=False)
    summary = [
        f'Pedestal QA — {cfg.DET_NAME}  {cfg.RUN}/{cfg.SUB_RUN}',
        f'  connected pads w/ pedestal : {len(ped)}',
        f'  pedestal RMS               : median {ped.ped_rms.median():.2f}, '
        f'90pct {ped.ped_rms.quantile(.9):.2f}, max {ped.ped_rms.max():.2f} ADC',
        f'  nominal threshold {args.nsigma:g}*RMS   : median {ped.thr.median():.1f}, '
        f'max {ped.thr.max():.1f} ADC (signal MPV ~386)',
        f'  corr(ped_rms, firing rate) : {r:.3f}  '
        f'(=> pedestal noise does {"" if abs(r)>0.3 else "NOT "}drive the hit rate)',
        f'  low-RMS (noisy) channels   : {int(ped["ped_noisy"].sum())} (z<-3)',
        f'  inflated-RMS channels      : {int(ped["ped_hot_rms"].sum())} (z>{args.z_rms:g})',
        f'  spark-veto applied to rate : {"yes" if args.veto_sparks else "no"}',
    ]
    txt = '\n'.join(summary)
    print(txt)
    with open(f'{out_dir}/pedestal_qa_summary{suffix}.txt', 'w') as f:
        f.write(txt + '\n')
    print(f'  saved {out_dir}/pedestal_outliers{suffix}.csv')
    print('Done.')


if __name__ == '__main__':
    main()
