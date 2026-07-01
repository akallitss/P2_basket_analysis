#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_map_validation.py

Validate the P2 (BASKET) pad mapping against real cosmic-bench data.

The point of this stage is to confirm that the Gerber-derived channel_id -> pad
map (Detector_Mapping/P2_BASKET/P2_BASKET_mapping.csv), combined with the DREAM
FEU wiring from run_config.json, places DAQ hits on the *physically correct*
pads. We do this by projecting the measured hit occupancy onto the readout PCB
and checking that it looks like a real detector illumination (contiguous active
region, dead connectors show up as coherent blocks, no scrambled speckle) rather
than noise sprayed at random.

Because the one uncertain link is the order of the 64 channels inside each
connector half (see p2_mapping.STRATEGY), --compare-strategies renders the
pad-plane occupancy under each candidate ordering side by side so the correct
one can be picked by eye.

Products (written to <Analysis>/<detN>/<run>/<sub_run>/02_map_validation/):
  pad_occupancy.png        hits-per-pad on the physical pad plane (x,y) +
                           radial and phi profiles + per-connector occupancy
  pad_amplitude.png        mean pulse height per pad on the pad plane
  dead_pad_map.png         active vs silent pads (dead-connector diagnostic)
  strategy_compare.png     (with --compare-strategies) pad occupancy per ordering
  map_validation_stats.csv per-pad hit counts + geometry

Usage
  python3 02_map_validation.py [run_key] [--strategy linear] [--compare-strategies]
"""

import matplotlib
matplotlib.use('Agg')  # headless

import os
import glob
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot

import p2_qa_config as qa
import p2_mapping as pmap

_HIT_BRANCHES = ['eventId', 'channel', 'amplitude', 'saturated', 'feu', 'time']


# --------------------------------------------------------------------------- #
# Data loading
# --------------------------------------------------------------------------- #
def load_p2_hits(hits_dir, feus):
    """Load only the P2 FEU hits from every combined-hits ROOT file."""
    files = sorted(glob.glob(os.path.join(hits_dir, '*.root')))
    if not files:
        raise FileNotFoundError(f'No combined-hits ROOT files in {hits_dir}')
    feu_set = set(feus)
    parts = []
    for fp in files:
        arr = uproot.open(f'{fp}:hits').arrays(_HIT_BRANCHES, library='pd')
        parts.append(arr[arr['feu'].isin(feu_set)].copy())
    df = pd.concat(parts, ignore_index=True)
    print(f'Loaded {len(df):,} P2 hits (FEUs {sorted(feu_set)}) from '
          f'{len(files)} files over {df["eventId"].nunique():,} events.')
    return df


# --------------------------------------------------------------------------- #
# Per-pad aggregation
# --------------------------------------------------------------------------- #
def per_pad_stats(hits, channel_table):
    """Aggregate hits onto pads. Returns the full pad table (all 1280 pads,
    even silent ones) with n_hits and mean amplitude."""
    h = pmap.attach_pads_to_hits(hits, channel_table)
    agg = (h.groupby('channel_id')
           .agg(n_hits=('amplitude', 'size'),
                mean_amp=('amplitude', 'mean'),
                n_sat=('saturated', 'sum'))
           .reset_index())
    # start from the full channel table so silent pads are kept with n_hits=0
    base = channel_table.drop_duplicates('channel_id').copy()
    full = base.merge(agg, on='channel_id', how='left')
    full['n_hits'] = full['n_hits'].fillna(0).astype(int)
    full['n_sat'] = full['n_sat'].fillna(0).astype(int)
    return full


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
def _pad_scatter(ax, full, value, title, cbar_label, cmap='viridis',
                 log=False, mask_zero=False):
    d = full
    if mask_zero:
        d = full[full['n_hits'] > 0]
    v = d[value].to_numpy(dtype=float)
    norm = matplotlib.colors.LogNorm(vmin=max(v[v > 0].min(), 1), vmax=v.max()) \
        if (log and (v > 0).any()) else None
    sc = ax.scatter(d['pad_cx'], d['pad_cy'], c=v, s=14, cmap=cmap,
                    norm=norm, marker='s', linewidths=0)
    ax.set_aspect('equal')
    ax.set_xlabel('pad_cx [mm]')
    ax.set_ylabel('pad_cy [mm]')
    ax.set_title(title, fontsize=9)
    plt.colorbar(sc, ax=ax, label=cbar_label, fraction=0.046, pad=0.04)


def plot_pad_occupancy(full, out_dir, tag):
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.4, 1])

    ax0 = fig.add_subplot(gs[0, 0])
    _pad_scatter(ax0, full, 'n_hits', 'Hits per pad on the readout PCB',
                 'hits', log=True, mask_zero=True)

    # per-connector (=sector) occupancy
    ax1 = fig.add_subplot(gs[0, 1])
    conn = full.groupby('connector_N')['n_hits'].sum()
    ax1.bar(conn.index, conn.values, color='steelblue')
    ax1.set_xlabel('physical connector N (1..10)')
    ax1.set_ylabel('total hits')
    ax1.set_title('Occupancy per connector (silent connector = mapping/DAQ flag)',
                  fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')

    # radial profile
    ax2 = fig.add_subplot(gs[1, 0])
    hit = full[full['n_hits'] > 0]
    ax2.hist(hit['radius'], bins=50, weights=hit['n_hits'], color='darkorange')
    ax2.set_xlabel('radius [mm]')
    ax2.set_ylabel('hits')
    ax2.set_title('Radial hit profile', fontsize=9)
    ax2.grid(True, alpha=0.3)

    # phi profile
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.hist(hit['phi'], bins=50, weights=hit['n_hits'], color='seagreen')
    ax3.set_xlabel('phi [rad]')
    ax3.set_ylabel('hits')
    ax3.set_title('Azimuthal (phi) hit profile', fontsize=9)
    ax3.grid(True, alpha=0.3)

    fig.suptitle(f'P2 pad-map validation — {tag}', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p = os.path.join(out_dir, 'pad_occupancy.png')
    fig.savefig(p, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {p}')


def plot_pad_amplitude(full, out_dir, tag):
    fig, ax = plt.subplots(figsize=(9, 8))
    _pad_scatter(ax, full, 'mean_amp', 'Mean pulse height per pad',
                 'mean amplitude [ADC]', cmap='plasma', mask_zero=True)
    fig.suptitle(f'P2 pad-map validation — mean amplitude — {tag}', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = os.path.join(out_dir, 'pad_amplitude.png')
    fig.savefig(p, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {p}')


def plot_dead_map(full, out_dir, tag):
    fig, ax = plt.subplots(figsize=(9, 8))
    silent = full[full['n_hits'] == 0]
    active = full[full['n_hits'] > 0]
    ax.scatter(active['pad_cx'], active['pad_cy'], s=14, marker='s',
               c='steelblue', label=f'active ({len(active)})', linewidths=0)
    ax.scatter(silent['pad_cx'], silent['pad_cy'], s=18, marker='s',
               c='crimson', label=f'silent ({len(silent)})', linewidths=0)
    ax.set_aspect('equal')
    ax.set_xlabel('pad_cx [mm]')
    ax.set_ylabel('pad_cy [mm]')
    ax.set_title('Active vs silent pads (coherent silent block = dead connector)',
                 fontsize=9)
    ax.legend(loc='upper right', fontsize=8)
    fig.suptitle(f'P2 pad-map validation — dead-pad map — {tag}', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = os.path.join(out_dir, 'dead_pad_map.png')
    fig.savefig(p, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {p}')


def plot_strategy_compare(hits, cfg, out_dir, tag,
                          strategies=('linear', 'reverse', 'pairswap')):
    fig, axes = plt.subplots(1, len(strategies),
                             figsize=(6 * len(strategies), 6.5))
    if len(strategies) == 1:
        axes = [axes]
    for ax, strat in zip(axes, strategies):
        ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                      det_type=cfg.DET_TYPE,
                                      det_name=cfg.DET_NAME, strategy=strat)
        full = per_pad_stats(hits, ct)
        _pad_scatter(ax, full, 'n_hits', f'strategy = {strat}',
                     'hits', log=True, mask_zero=True)
    fig.suptitle(f'P2 within-connector ordering comparison — {tag}\n'
                 '(the correct ordering gives a smooth, physical illumination)',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    p = os.path.join(out_dir, 'strategy_compare.png')
    fig.savefig(p, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {p}')


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description='Validate the P2 pad mapping.')
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN,
                    help='run key in p2_qa_config.RUNS')
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'],
                    help='within-connector channel ordering to use for the '
                         'main products (default reverse, validated vs M3)')
    ap.add_argument('--compare-strategies', action='store_true',
                    help='also render pad occupancy for every ordering side by side')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    out_dir = cfg.out_dir('02_map_validation')

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy=args.strategy)
    print(f'Mapping: {cfg.DET_NAME}  FEUs {ct.attrs["feus"]}  '
          f'strategy={args.strategy}  '
          f'{int(ct["mapped"].sum())}/{len(ct)} channels resolved to pads')

    hits = load_p2_hits(cfg.combined_hits_dir, ct.attrs['feus'])
    full = per_pad_stats(hits, ct)

    n_active = int((full['n_hits'] > 0).sum())
    print(f'Pads: {len(full)} total, {n_active} active, '
          f'{len(full) - n_active} silent.')
    print('Total hits per connector:')
    print(full.groupby('connector_N')['n_hits'].sum().to_string())

    tag = f'{cfg.DET_TAG} {cfg.RUN}/{cfg.SUB_RUN} [{args.strategy}]'
    plot_pad_occupancy(full, out_dir, tag)
    plot_pad_amplitude(full, out_dir, tag)
    plot_dead_map(full, out_dir, tag)
    if args.compare_strategies:
        plot_strategy_compare(hits, cfg, out_dir, tag)

    stats_csv = os.path.join(out_dir, 'map_validation_stats.csv')
    full.sort_values('channel_id').to_csv(stats_csv, index=False)
    print(f'  saved {stats_csv}')
    print('Done.')


if __name__ == '__main__':
    main()
