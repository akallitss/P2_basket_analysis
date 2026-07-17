#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_detector_deep_qa.py

Deep noise / pathology QA for the P2 pad detector. Pad-adapted version of
nTof_x17/mx_june_cosmic_qa/04_detector_deep_qa.py: the reference splits hits into
X-strip and Y-strip FEUs, but a P2 pad already carries a full (x,y) position, so
there is no plane split -- everything is per pad.

Products (written to <Analysis>/<detN>/<run>/<sub_run>/05_detector_deep_qa/):
  surface_hitmap.png        2D hitmap on the pad plane (linear + log)
  centroid_hitmap.png       per-event charge-weighted centroid map (all events +
                            muon-like n_pad<=3), live pad centres overlaid
  pad_firing_fraction.png   fraction of events each pad fires, on the pad plane +
                            distribution -- flags always-firing / hot pads
  event_multiplicity.png    # pads firing per event (total + per FEU), log y
  multiplicity_vs_time.png  mean pads/event over the run (spark clustering?)
  deep_qa_summary.txt       hottest pads, spark fraction, etc.

Usage: python3 05_detector_deep_qa.py [run_key] [--strategy reverse]
"""

import os
import argparse
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import p2_qa_config as qa
import p2_mapping as pmap
import p2_sparks as ps
import p2_io as p2io

_BRANCHES = ['eventId', 'trigger_timestamp_ns', 'channel', 'amplitude', 'feu']


def load_aggregates(cfg, channel_table, spark_veto=None, drop_pads=()):
    """Stream the combined-hits chunks and reduce each one to the aggregates
    the deep-QA plots need, so the full per-hit table never sits in memory
    (the det4 long run has ~18M hits per chunk -> OOM if concatenated).

    Returns (padc, events, n_hits):
      padc    per-channel_id: n (hits), nev (events the pad fired in)
      events  per-event: ts, npad (distinct pads), n_hit, charge-weighted
              centroid x/y, npad_feu<F> per FEU
      n_hits  total mapped hits
    """
    feus = sorted(channel_table.attrs['feus'])
    padc = None
    ev_parts = []
    n_hits = 0
    n_rm = n_burst = 0
    for df in p2io.iter_hits(cfg.combined_hits_dir, _BRANCHES, feus,
                             t_max_h=cfg.T_MAX_H, min_amp=cfg.MIN_AMP):
        if spark_veto is not None:
            df, rm = spark_veto.apply(df)
            n_rm += rm
            n_burst += spark_veto.last_burst_events
        h = pmap.attach_pads_to_hits(df, channel_table)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        if drop_pads:   # known-noisy channels (cfg.NOISY_PADS), on request
            h = h[~h['channel_id'].isin(set(drop_pads))]
        del df
        if not len(h):
            continue
        n_hits += len(h)
        c = h.groupby('channel_id').agg(n=('eventId', 'size'),
                                        nev=('eventId', 'nunique'))
        padc = c if padc is None else padc.add(c, fill_value=0)
        w = h['amplitude'].clip(lower=0).astype(np.float64)
        t = h.assign(_wx=w * h['pad_cx'], _wy=w * h['pad_cy'], _w=w)
        g = t.groupby('eventId')
        ev = pd.DataFrame({'ts': g['trigger_timestamp_ns'].first(),
                           'npad': g['channel_id'].nunique(),
                           'n_hit': g.size(),
                           'x': g['_wx'].sum() / g['_w'].sum(),
                           'y': g['_wy'].sum() / g['_w'].sum()})
        for f in feus:
            nf = h[h['feu'] == f].groupby('eventId')['channel_id'].nunique()
            ev[f'npad_feu{f}'] = nf.reindex(ev.index).fillna(0).astype(int)
        ev_parts.append(ev.reset_index())
    events = pd.concat(ev_parts, ignore_index=True)
    padc = padc.astype({'n': int, 'nev': int})
    if spark_veto is not None:
        print(f'Spark veto: dropped {n_rm:,} hits in {len(spark_veto.sparks)} '
              f'sparks ({100*(1-spark_veto.live_fraction()):.2f}% deadtime) + '
              f'{n_burst} burst events (>= {spark_veto.burst_npads} pads).')
    return padc, events, n_hits


def _pad_range(df, pad=15):
    return [[df['pad_cx'].min() - pad, df['pad_cx'].max() + pad],
            [df['pad_cy'].min() - pad, df['pad_cy'].max() + pad]]


def plot_surface_hitmap(padc, channel_table, out_dir, cfg, n_hits, n_events,
                        suffix='', pillars=None):
    """Hit counts drawn on the real pad tiles (not a uniform-grid histogram, so
    the fan geometry tiles contiguously). Never-fired live pads are grey."""
    from matplotlib.collections import PolyCollection

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    if pmap.has_tile_geometry(channel_table):
        pads, verts = pmap.pad_tiles(channel_table)
        counts = padc['n'].reindex(pads['channel_id']).fillna(0).to_numpy()
        fired = counts > 0
        for ax, norm, tag in [(axes[0], None, 'linear'),
                              (axes[1], matplotlib.colors.LogNorm(vmin=1), 'log')]:
            dead = PolyCollection(verts[~fired], facecolors='0.92',
                                  edgecolors='0.7', linewidths=0.3)
            ax.add_collection(dead)
            pc = PolyCollection(verts[fired], array=counts[fired], cmap='inferno',
                                norm=norm, edgecolors='face', linewidths=0.2)
            ax.add_collection(pc)
            fig.colorbar(pc, ax=ax, label=f'pad hits ({tag})')
            ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
            ax.autoscale_view(); ax.set_aspect('equal')
            ax.set_title(f'Surface hitmap ({tag}) — grey = never fired')
            if pillars is not None and len(pillars):
                pmap.draw_pillars(ax, pillars, small=False)
                ax.legend(loc='upper right', fontsize=7, framealpha=0.9)
    else:  # old map CSV without tile geometry: fall back to the grid histogram
        gm = (channel_table.drop_duplicates('channel_id')
              .merge(padc.reset_index(), on='channel_id', how='left'))
        gm['n'] = gm['n'].fillna(0)
        rng = _pad_range(gm)
        for ax, norm, tag in [(axes[0], None, 'linear'),
                              (axes[1], matplotlib.colors.LogNorm(), 'log')]:
            h = ax.hist2d(gm['pad_cx'], gm['pad_cy'], bins=120, range=rng,
                          weights=gm['n'], cmap='inferno', norm=norm)
            fig.colorbar(h[3], ax=ax, label=f'pad hits ({tag})')
            ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
            ax.set_aspect('equal'); ax.set_title(f'Surface hitmap ({tag})')
    fig.suptitle(f'{cfg.DET_NAME} pad surface hitmap — {cfg.RUN}/{cfg.SUB_RUN}\n'
                 f'{n_hits:,} hits from {n_events:,} events')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/surface_hitmap{suffix}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_centroid_hitmap(events, channel_table, out_dir, cfg, suffix='',
                         pillars=None):
    """Per-EVENT charge-weighted centroid hitmap (complements the per-hit
    surface_hitmap): left all events, right muon-like events only (n_pad <= 3),
    with the live pad centres overlaid. Interior white holes at this statistics
    are genuinely dead zones; the striping is pad-row quantisation."""
    cen = events[np.isfinite(events['x']) & np.isfinite(events['y'])]
    sel = cen[cen['n_hit'] <= 3]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, c, tag in [(axes[0], cen, f'all events ({len(cen):,})'),
                       (axes[1], sel, f'n_pad <= 3 ({len(sel):,})')]:
        hb = ax.hexbin(c['x'], c['y'], gridsize=90, cmap='viridis', mincnt=1)
        ax.scatter(channel_table['pad_cx'], channel_table['pad_cy'],
                   s=1.0, c='red', alpha=0.2, linewidths=0)
        fig.colorbar(hb, ax=ax, label='events / bin')
        ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
        ax.set_aspect('equal'); ax.set_title(tag)
        if pillars is not None and len(pillars):
            pmap.draw_pillars(ax, pillars, small=False)
            ax.legend(loc='upper right', fontsize=7, framealpha=0.9)
    fig.suptitle(f'{cfg.DET_NAME} event centroid hitmap — {cfg.RUN}/{cfg.SUB_RUN}\n'
                 'charge-weighted pad centroid per event; red = live pad centres')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/centroid_hitmap{suffix}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    return len(cen), len(sel)


def pad_firing(padc, n_events, channel_table):
    fire = (padc['nev'] / max(n_events, 1)).rename('fire_frac')
    base = channel_table.drop_duplicates('channel_id').set_index('channel_id')
    out = base[['pad_cx', 'pad_cy', 'connector_N']].join(fire)
    out['fire_frac'] = out['fire_frac'].fillna(0.0)
    return out.reset_index()


def plot_pad_firing(fdf, out_dir, cfg, summary, suffix=''):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    live = fdf[fdf['fire_frac'] > 0]
    sc = axes[0].scatter(live['pad_cx'], live['pad_cy'], c=live['fire_frac'],
                         s=16, marker='s', cmap='viridis',
                         norm=matplotlib.colors.LogNorm(
                             vmin=max(live['fire_frac'][live['fire_frac'] > 0].min(), 1e-5),
                             vmax=live['fire_frac'].max()))
    axes[0].set_aspect('equal'); axes[0].set_xlabel('pad_cx [mm]')
    axes[0].set_ylabel('pad_cy [mm]'); axes[0].set_title('Per-pad firing fraction')
    plt.colorbar(sc, ax=axes[0], label='fraction of events pad fires',
                 fraction=0.046, pad=0.04)

    axes[1].hist(live['fire_frac'], bins=60, color='steelblue')
    med = live['fire_frac'].median()
    axes[1].axvline(med, color='green', ls='--', label=f'median {med:.4f}')
    axes[1].set_yscale('log'); axes[1].set_xlabel('fire fraction')
    axes[1].set_ylabel('pads'); axes[1].legend(); axes[1].grid(True, alpha=0.3)
    axes[1].set_title('Firing-fraction distribution')
    fig.suptitle(f'{cfg.DET_NAME} per-pad firing fraction — {cfg.RUN}/{cfg.SUB_RUN}')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/pad_firing_fraction{suffix}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    top = fdf.sort_values('fire_frac', ascending=False).head(10)
    summary.append('  hottest pads (channel_id, connector, pos_mm, fire_frac):')
    for _, r in top.iterrows():
        summary.append(f'    ch_id {int(r.channel_id):4d}  conn {int(r.connector_N):2d}  '
                       f'({r.pad_cx:6.1f},{r.pad_cy:6.1f})  {r.fire_frac:.4f}')


def plot_multiplicity(events, feus, mult, out_dir, cfg, summary, suffix=''):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    bins = np.arange(0, max(int(mult['npad'].max()), 10) + 2)
    # per-FEU multiplicity (events where that FEU fired at all)
    for feu in sorted(feus):
        nf = events.loc[events[f'npad_feu{feu}'] > 0, f'npad_feu{feu}']
        axes[0].hist(nf, bins=bins, histtype='step', lw=1.4, label=f'FEU {feu}')
    axes[0].set_yscale('log'); axes[0].set_xlabel('pads firing per event (per FEU)')
    axes[0].set_ylabel('events'); axes[0].legend(); axes[0].grid(True, alpha=0.3)
    axes[0].set_title('Per-FEU pad multiplicity')

    axes[1].hist(mult['npad'], bins=bins, color='purple', alpha=0.8)
    axes[1].set_yscale('log'); axes[1].set_xlabel('total pads firing per event')
    axes[1].set_ylabel('events'); axes[1].grid(True, alpha=0.3)
    axes[1].set_title('Total pad multiplicity')
    fig.suptitle(f'{cfg.DET_NAME} event multiplicity — {cfg.RUN}/{cfg.SUB_RUN}')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/event_multiplicity{suffix}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    for thr in (5, 10, 20):
        frac = (mult['npad'] > thr).mean()
        summary.append(f'  events with >{thr} pads: {frac*100:.2f}%')
    summary.append(f'  median pads/event: {mult["npad"].median():.0f}, '
                   f'mean: {mult["npad"].mean():.2f}, max: {int(mult["npad"].max())}')


def plot_mult_vs_time(mult, out_dir, cfg, suffix=''):
    from scipy.stats import binned_statistic
    t = (mult['ts'] - mult['ts'].min()) / 1e9
    mean_m, edges, _ = binned_statistic(t, mult['npad'], statistic='mean', bins=80)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(ctr, mean_m, 'o-', ms=3, color='purple')
    ax.set_xlabel('Time since run start [s]'); ax.set_ylabel('mean pads/event')
    ax.set_title(f'{cfg.DET_NAME} mean multiplicity vs time — {cfg.RUN}/{cfg.SUB_RUN}')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/multiplicity_vs_time{suffix}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'])
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True,
                    help='drop events taken during an HV spark (default on).')
    ap.add_argument('--mask-noisy-pads', action='store_true',
                    help='ALSO drop the known-noisy pads (cfg.NOISY_PADS) so '
                         'the deep-QA maps show the clean detector; products '
                         'get a _noisy_masked suffix. Default off: the QA '
                         'should show the pathology.')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    out_dir = cfg.out_dir('05_detector_deep_qa')
    sfx = cfg.product_suffix(args.veto_sparks)
    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy=args.strategy,
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    if cfg.DEAD_CONNECTORS:
        print(f'  dropped dead connectors: {list(cfg.DEAD_CONNECTORS)}')
    n_pads = ct['channel_id'].nunique()
    drop = cfg.NOISY_PADS if args.mask_noisy_pads else ()
    if drop:
        sfx += '_noisy_masked'
        print(f'  masking noisy pads: {list(drop)}')
    sv = ps.SparkVeto.from_cfg(cfg) if args.veto_sparks else None
    padc, events, n_hits = load_aggregates(cfg, ct, spark_veto=sv,
                                           drop_pads=drop)
    n_events = len(events)
    summary = [f'Deep QA — {cfg.DET_NAME}  {cfg.RUN}/{cfg.SUB_RUN}',
               f'  spark veto: {"ON" if args.veto_sparks else "OFF"}'
               + (f' ({100*(1-sv.live_fraction()):.2f}% deadtime removed)' if sv else ''),
               f'  events with any pad hit: {n_events:,}',
               f'  total pad hits: {n_hits:,}  ({n_hits/max(n_events, 1):.2f} hits/event)',
               f'  distinct pads fired: {len(padc)} / {n_pads}']

    pil = pmap.load_pillars(cfg.MASK_GBR_PATH)
    if len(pil):
        print(f'  insulation-mask pillars: {int(pil["big"].sum())} big + '
              f'{int((~pil["big"]).sum()):,} small (overlaid on hitmaps)')
    plot_surface_hitmap(padc, ct, out_dir, cfg, n_hits, n_events, sfx, pillars=pil)
    n_cen, n_muonlike = plot_centroid_hitmap(events, ct, out_dir, cfg, sfx,
                                             pillars=pil)
    summary.append(f'  event centroids: {n_cen:,}  (muon-like n_pad<=3: {n_muonlike:,})')
    fdf = pad_firing(padc, n_events, ct)
    plot_pad_firing(fdf, out_dir, cfg, summary, sfx)
    mult = events[['npad', 'ts']]
    plot_multiplicity(events, ct.attrs['feus'], mult, out_dir, cfg, summary, sfx)
    plot_mult_vs_time(mult, out_dir, cfg, sfx)

    txt = '\n'.join(summary)
    print(txt)
    with open(f'{out_dir}/deep_qa_summary{sfx}.txt', 'w') as f:
        f.write(txt + '\n')
    print(f'\nDeep QA written to: {out_dir}')


if __name__ == '__main__':
    main()
