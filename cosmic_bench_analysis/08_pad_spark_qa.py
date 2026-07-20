#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_pad_spark_qa.py

Review the hit-level (per-pad) spark flag (p2_pad_sparks) BEFORE deciding whether
to mask these pads in the analysis. Individual pads can micro-spark in a way the
HV-monitor veto cannot see (too small to move the total mesh current, shorter
than the 2 s HV log); they survive as "hot" pads with an abnormal saturated /
high-amplitude fraction. This stage flags them by that discharge signature and
shows what gets flagged and what the pad plane looks like once they are removed.

By default it operates on the HV-veto-cleaned hits, so what it flags is the
*residual* single-pad sparking left after the whole-detector veto.

Products (written to <Analysis>/<detN>/<run>/<sub_run>/08_pad_spark_qa/):
  pad_spark_map.png          pad plane: flagged (sparking) pads highlighted
  spark_pad_scatter.png      sat_rate vs high-amp fraction, thresholds + flag
  firing_before_after.png    per-pad firing map before vs after masking flags
  flagged_per_connector.png  flagged-pad count per connector
  flagged_pads.csv           the flagged pad list + metrics
  pad_spark_qa_summary.txt

Usage: python3 08_pad_spark_qa.py [run_key] [--z-sat 4] [--z-hi 3]
                                  [--hi-adc 1000] [--min-hits 20] [--no-veto-sparks]
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

_BRANCHES = ['eventId', 'channel', 'amplitude', 'saturated', 'feu']


def load_hits(cfg, ct, veto_sparks=True):
    """Compact mapped-hit table (eventId, channel_id, amplitude, saturated).

    Streams the combined-hits chunks (p2_io) and keeps only the downcast
    columns the pad-spark metrics need, so the full run costs ~1.5 GB instead
    of OOM-ing when every branch of every chunk is concatenated."""
    br = _BRANCHES + (['trigger_timestamp_ns'] if veto_sparks else [])
    sv = ps.SparkVeto.from_cfg(cfg) if veto_sparks else None
    parts = []
    n_rm = n_burst = 0
    for a in p2io.iter_hits(cfg.combined_hits_dir, br, ct.attrs['feus'],
                            t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H, min_amp=cfg.MIN_AMP):
        if sv is not None:
            a, rm = sv.apply(a)
            n_rm += rm
            n_burst += sv.last_burst_events
        h = pmap.attach_pads_to_hits(a, ct)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        del a
        parts.append(pd.DataFrame({
            'eventId': h['eventId'].astype(np.int32),
            'channel_id': h['channel_id'].astype(np.int32),
            'amplitude': h['amplitude'].astype(np.float32),
            'saturated': h['saturated'].astype(bool)}))
    h = pd.concat(parts, ignore_index=True)
    if sv is not None:
        print(f'HV spark veto: dropped {n_rm:,} hits '
              f'({100*(1-sv.live_fraction()):.2f}% deadtime) + '
              f'{n_burst} burst events (>= {sv.burst_npads} pads).')
    return h


def plot_pad_map(mask, out_dir, cfg, suffix):
    m = mask.metrics
    fl = mask.flagged
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.scatter(m['pad_cx'], m['pad_cy'], s=12, marker='s', c='lightgrey',
               linewidths=0, label=f'pads ({len(m)})')
    sc = ax.scatter(fl['pad_cx'], fl['pad_cy'], s=60, marker='s',
                    c=fl['sat_rate'], cmap='autumn_r', edgecolors='k',
                    linewidths=0.6, label=f'flagged sparking ({len(fl)})')
    plt.colorbar(sc, ax=ax, label='saturated-hit fraction', fraction=0.046, pad=0.04)
    if mask.noise_ids:
        nz = mask.noise
        ax.scatter(nz['pad_cx'], nz['pad_cy'], s=70, marker='D', c='purple',
                   edgecolors='k', linewidths=0.6, label=f'noise ({len(nz)})')
    ax.set_aspect('equal'); ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
    ax.set_title(f'{cfg.DET_NAME} per-pad spark flag on the pad plane\n'
                 f'(sat_z>{mask.z_sat:g} & hi_z>{mask.z_hi:g}) — {cfg.RUN}/{cfg.SUB_RUN}')
    ax.legend(loc='upper right', fontsize=8)
    fig.tight_layout()
    p = f'{out_dir}/pad_spark_map{suffix}.png'
    fig.savefig(p, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f'  saved {p}')


def plot_scatter(mask, out_dir, cfg, suffix):
    m = mask.metrics[mask.metrics['n'] >= mask.min_hits]
    fl = m['spark_pad']
    # threshold values from robust median/MAD
    def thr(col, z):
        med = np.median(m[col]); mad = np.median(np.abs(m[col] - med)) * 1.4826
        return med + z * mad
    fig, ax = plt.subplots(figsize=(8, 6.5))
    ax.scatter(m.loc[~fl, 'sat_rate'], m.loc[~fl, 'hi_frac'], s=10,
               c='steelblue', alpha=0.5, linewidths=0, label='normal pads')
    ax.scatter(m.loc[fl, 'sat_rate'], m.loc[fl, 'hi_frac'], s=40, c='crimson',
               edgecolors='k', linewidths=0.5, label='flagged sparking')
    ax.axvline(thr('sat_rate', mask.z_sat), color='crimson', ls='--', lw=1,
               label=f'sat threshold (z={mask.z_sat:g})')
    ax.axhline(thr('hi_frac', mask.z_hi), color='darkorange', ls='--', lw=1,
               label=f'hi-amp threshold (z={mask.z_hi:g})')
    ax.set_xlabel('saturated-hit fraction per pad')
    ax.set_ylabel(f'fraction of hits > {mask.hi_adc:g} ADC')
    ax.set_title(f'{cfg.DET_NAME} pad discharge signature — {cfg.RUN}/{cfg.SUB_RUN}')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = f'{out_dir}/spark_pad_scatter{suffix}.png'
    fig.savefig(p, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f'  saved {p}')


def plot_rate_vs_amp(mask, out_dir, cfg, suffix):
    """Per-pad firing rate vs median amplitude: separates the NOISE population
    (high rate, low amplitude) from the DISCHARGE population (high amplitude)."""
    m = mask.metrics[mask.metrics['n'] >= mask.min_hits]
    normal = m[~m['noise_pad'] & ~m['spark_pad']]
    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.scatter(normal['med_amp'], normal['fire_frac'], s=10, c='steelblue',
               alpha=0.5, linewidths=0, label=f'normal ({len(normal)})')
    if mask.flagged_ids:
        f = mask.flagged
        ax.scatter(f['med_amp'], f['fire_frac'], s=45, c='crimson',
                   edgecolors='k', linewidths=0.5,
                   label=f'discharge/spark ({len(f)})')
    if mask.noise_ids:
        nz = mask.noise
        ax.scatter(nz['med_amp'], nz['fire_frac'], s=45, c='purple',
                   marker='D', edgecolors='k', linewidths=0.5,
                   label=f'noise ({len(nz)})')
    ax.axvline(mask.noise_amp, color='purple', ls='--', lw=1,
               label=f'noise band (<{mask.noise_amp:g} ADC)')
    ax.axvline(m['med_amp'].median(), color='green', ls=':', lw=1,
               label=f'median amp ({m["med_amp"].median():.0f})')
    ax.set_xlabel('per-pad median amplitude [ADC]')
    ax.set_ylabel('per-pad firing fraction')
    ax.set_title(f'{cfg.DET_NAME} pad rate vs amplitude — noise vs discharge '
                 f'separation\n{cfg.RUN}/{cfg.SUB_RUN}')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = f'{out_dir}/rate_vs_amplitude{suffix}.png'
    fig.savefig(p, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f'  saved {p}')


def plot_before_after(mask, hits, out_dir, cfg, suffix):
    m = mask.metrics
    ids = mask.flagged_ids
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
    for ax, keep_all, ttl in [(axes[0], True, 'before (all pads)'),
                              (axes[1], False, 'after masking flagged pads')]:
        d = m if keep_all else m[~m.index.isin(ids)]
        vmax = m['fire_frac'].max()
        sc = ax.scatter(d['pad_cx'], d['pad_cy'], s=14, marker='s',
                        c=d['fire_frac'], cmap='viridis',
                        norm=matplotlib.colors.LogNorm(
                            vmin=max(m['fire_frac'][m['fire_frac'] > 0].min(), 1e-5),
                            vmax=vmax), linewidths=0)
        ax.set_aspect('equal'); ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
        ax.set_title(ttl, fontsize=10)
        plt.colorbar(sc, ax=ax, label='firing fraction', fraction=0.046, pad=0.04)
    fig.suptitle(f'{cfg.DET_NAME} per-pad firing fraction, spark-pad mask effect '
                 f'— {cfg.RUN}/{cfg.SUB_RUN}', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = f'{out_dir}/firing_before_after{suffix}.png'
    fig.savefig(p, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f'  saved {p}')


def plot_per_connector(mask, out_dir, cfg, suffix):
    m = mask.metrics
    act = m[m['n'] >= mask.min_hits].groupby('connector_N').size()
    fla = mask.flagged.groupby('connector_N').size().reindex(act.index, fill_value=0)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(act.index, act.values, color='lightgrey', label='active pads')
    ax.bar(fla.index, fla.values, color='crimson', label='flagged sparking')
    ax.set_xlabel('physical connector N'); ax.set_ylabel('pads')
    ax.set_title(f'{cfg.DET_NAME} flagged sparking pads per connector '
                 f'— {cfg.RUN}/{cfg.SUB_RUN}')
    for x, v in zip(fla.index, fla.values):
        if v:
            ax.text(x, v + 1, str(int(v)), ha='center', color='crimson', fontsize=9)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    p = f'{out_dir}/flagged_per_connector{suffix}.png'
    fig.savefig(p, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f'  saved {p}')


def main():
    ap = argparse.ArgumentParser(description='Review the per-pad spark flag.')
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'])
    ap.add_argument('--z-sat', type=float, default=4.0,
                    help='saturation-rate robust-z threshold (default 4).')
    ap.add_argument('--z-hi', type=float, default=3.0,
                    help='high-amplitude-fraction robust-z threshold (default 3).')
    ap.add_argument('--hi-adc', type=float, default=pps.HI_ADC,
                    help='amplitude [ADC] above which a hit is "high" (default 1000).')
    ap.add_argument('--min-hits', type=int, default=20,
                    help='minimum hits for a pad to be flag-eligible (default 20).')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True,
                    help='flag on HV-veto-cleaned hits (default on).')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    out_dir = cfg.out_dir('08_pad_spark_qa')
    suffix = cfg.product_suffix(args.veto_sparks)

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy=args.strategy,
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    hits = load_hits(cfg, ct, args.veto_sparks)
    n_events = hits['eventId'].nunique()

    mask = pps.PadSparkMask.from_hits(hits, ct, n_events, args.z_sat, args.z_hi,
                                      args.min_hits, args.hi_adc)
    print(mask.summary())

    _, n_rm = mask.apply_hits(hits)
    frac = 100 * n_rm / len(hits) if len(hits) else 0.0

    plot_rate_vs_amp(mask, out_dir, cfg, suffix)
    plot_pad_map(mask, out_dir, cfg, suffix)
    plot_scatter(mask, out_dir, cfg, suffix)
    plot_before_after(mask, hits, out_dir, cfg, suffix)
    plot_per_connector(mask, out_dir, cfg, suffix)

    # masked-pad list (noise + discharge), tagged by flag type
    cols = ['connector_N', 'radius', 'pad_cx', 'pad_cy', 'n', 'fire_frac',
            'sat_rate', 'hi_frac', 'lo_frac', 'mean_amp', 'med_amp',
            'fire_z', 'sat_z', 'hi_z']
    out = pd.concat([
        mask.noise.assign(flag_type='noise'),
        mask.flagged.assign(flag_type='spark')])
    out.reset_index()[['channel_id', 'flag_type'] + cols].to_csv(
        f'{out_dir}/flagged_pads{suffix}.csv', index=False)

    summary = [mask.summary(), '',
               f'  hits on masked pads  : {n_rm:,} / {len(hits):,} ({frac:.2f}%)',
               f'  HV spark veto applied: {"yes" if args.veto_sparks else "no"}',
               f'  n_events             : {n_events:,}']
    txt = '\n'.join(summary)
    print(txt)
    with open(f'{out_dir}/pad_spark_qa_summary{suffix}.txt', 'w') as f:
        f.write(txt + '\n')
    print(f'  saved {out_dir}/flagged_pads{suffix}.csv')
    print('Done.')


if __name__ == '__main__':
    main()
