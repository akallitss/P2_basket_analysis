#!/usr/bin/env python3
"""28_eff_vs_m3_chi2.py -- P2 efficiency vs the M3 reference-track quality.

Usage:  python3 28_eff_vs_m3_chi2.py <run_key> [--no-veto-sparks]

Loads EVERY single-track M3 event (no chi2 cut), applies exactly the stage-06
selection otherwise (time window, spark veto, active area incl. dead-pad cells,
clustered P2 position, 06's pad->M3 transform), and measures the P2 efficiency
(reco within MATCH_R, and any-pad) against chi2_max = max(Chi2X, Chi2Y):
  * differential: efficiency in chi2 bins  -> does it fall for poorer tracks?
  * cumulative:   efficiency for a cut chi2 < c, with the track count kept
A reference-limited efficiency rises as the cut tightens and only plateaus once
the bad/fake M3 tracks are gone; the plateau height is the detector's own.
Outputs: 28_eff_vs_m3_chi2/eff_vs_chi2<sfx>.png + .csv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import awkward as ak
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa      # noqa: E402
import p2_mapping as pmap      # noqa: E402
import p2_align as pa          # noqa: E402
import p2_io as p2io           # noqa: E402
import p2_sparks as ps         # noqa: E402
from M3RefTracking import M3RefTracking  # noqa: E402

C_RECO, C_ANY = '#257de0', '#f47230'
INK, INK2, GRID = '#1f1f1f', '#5c5c5c', '#e6e6e6'


def binom_err(k, n):
    p = k / np.maximum(n, 1)
    return np.sqrt(np.clip(p * (1 - p), 0, None) / np.maximum(n, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction, default=True)
    args = ap.parse_args()
    cfg = qa.get_config(args.run_key)
    out_dir = cfg.out_dir('28_eff_vs_m3_chi2')
    sfx = cfg.product_suffix(args.veto_sparks)
    R, AR = cfg.MATCH_R, cfg.ACTIVE_R

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
                                  det_name=cfg.DET_NAME, drop_connectors=cfg.DEAD_CONNECTORS,
                                  drop_channels=cfg.DEAD_CHANNELS,
                                  cluster_link=cfg.CLUSTER_LINK_MM,
                                  strategy_overrides=cfg.STRATEGY_OVERRIDES)
    padc = ct[ct['mapped']].drop_duplicates('channel_id').reset_index(drop=True)
    fp = pd.read_csv(f'{cfg.out_dir("06_efficiency")}/pad_footprint{sfx}.csv')
    if len(fp) != len(padc):
        raise SystemExit('06 pad footprint does not match the channel table -- rerun 06 first')
    P = np.c_[padc.pad_cx, padc.pad_cy, np.ones(len(padc))]
    coef, *_ = np.linalg.lstsq(P, fp[['x', 'y']].to_numpy(), rcond=None)

    def to_m3(x, y):
        q = np.c_[x, y, np.ones(len(x))] @ coef
        return q[:, 0], q[:, 1]

    # ---- M3: all single-track events, chi2 kept per event ----
    z = cfg.det_plane_z()
    m3 = M3RefTracking(os.path.join(cfg.m3_tracking_dir, ''), single_track=True,
                       chi2_cut=1e12, min_nclus=qa.M3_MIN_NCLUS)
    rd = m3.ray_data

    def col(n):
        return np.asarray(ak.to_numpy(ak.ravel(rd[n])), dtype=float)
    ev = np.asarray(ak.to_numpy(ak.ravel(rd['evn'])), dtype=np.int64)
    t = (z - col('Z_Down')) / (col('Z_Up') - col('Z_Down'))
    x = col('X_Down') + t * (col('X_Up') - col('X_Down'))
    y = col('Y_Down') + t * (col('Y_Up') - col('Y_Down'))
    d = pd.DataFrame({'eventId': ev, 'x': x, 'y': y,
                      'chi2': np.maximum(col('Chi2X'), col('Chi2Y'))})
    d = d[np.isfinite(d.x) & np.isfinite(d.y)].drop_duplicates('eventId')
    d = pa.filter_events_by_time(d, cfg.m3_tracking_dir, cfg.T_MAX_H, t_min_h=cfg.T_MIN_H)

    p2, hit_events = pa.load_p2_centroids(cfg.combined_hits_dir, ct, min_amp=cfg.MIN_AMP,
                                          t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H,
                                          drop_pads=p2io.drop_pads_for(cfg, ct))
    if args.veto_sparks:
        sv = ps.SparkVeto.from_cfg(cfg)
        bad = sv.vetoed_ids_from_hits(cfg.combined_hits_dir, ct.attrs['feus'], min_amp=cfg.MIN_AMP)
        d = d[~d['eventId'].isin(bad)]
        p2 = p2[~p2['eventId'].isin(bad)]
        hit_events = hit_events - bad

    # ---- active area (live pad within ACTIVE_R, not in a dead-pad cell) ----
    tree = cKDTree(fp[['x', 'y']].to_numpy())
    nn, _ = tree.query(d[['x', 'y']].to_numpy())
    dead = ct.attrs.get('dead_pad_xy')
    if dead is not None and len(dead):
        dx, dy = to_m3(dead[:, 0], dead[:, 1])
        dd, _ = cKDTree(np.c_[dx, dy]).query(d[['x', 'y']].to_numpy())
        nn = np.where(dd < nn, np.inf, nn)
    d = d[nn <= AR].copy()
    rx, ry = to_m3(p2['x_pad'].to_numpy(), p2['y_pad'].to_numpy())
    reco = pd.DataFrame({'eventId': p2['eventId'].to_numpy(), 'rx': rx, 'ry': ry})
    d = d.merge(reco, on='eventId', how='left')
    d['within'] = np.hypot(d.x - d.rx, d.y - d.ry) <= R
    d['has_any'] = d['eventId'].isin(hit_events)
    N = len(d)
    print(f'{cfg.DET_TAG} {cfg.RUN}/{cfg.SUB_RUN}: active-area single-track rays (no chi2 cut) {N:,}')

    edges = np.array([0, 0.5, 1, 2, 3, 5, 10, 20, 50, 200, np.inf])
    d['bin'] = pd.cut(d['chi2'], edges, right=False)
    g = d.groupby('bin', observed=False).agg(n=('within', 'size'), k=('within', 'sum'), a=('has_any', 'sum'))
    g['eff'] = g.k / g.n.clip(lower=1); g['eff_err'] = binom_err(g.k, g.n)
    g['any'] = g.a / g.n.clip(lower=1); g['any_err'] = binom_err(g.a, g.n)
    cuts = np.array([0.25, 0.5, 1, 1.5, 2, 3, 5, 7.5, 10, 20, 50, 200, 1e12])
    cum = []
    for c in cuts:
        s = d[d.chi2 < c]
        cum.append((c, len(s), s.within.mean() if len(s) else np.nan, s.has_any.mean() if len(s) else np.nan,
                    binom_err(s.within.sum(), len(s)), binom_err(s.has_any.sum(), len(s))))
    cum = pd.DataFrame(cum, columns=['chi2_cut', 'n', 'eff', 'any', 'eff_err', 'any_err'])
    g.reset_index().assign(bin=lambda t: t['bin'].astype(str)).to_csv(f'{out_dir}/eff_vs_chi2_bins{sfx}.csv', index=False)
    cum.to_csv(f'{out_dir}/eff_vs_chi2_cut{sfx}.csv', index=False)
    print('differential:'); print(g[['n', 'eff', 'eff_err', 'any']].round(4).to_string())
    print('cumulative:'); print(cum.round(4).to_string(index=False))

    # ---- figure: 2 columns (differential | cumulative), count panel under each ----
    plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'xtick.color': INK2,
                         'ytick.color': INK2, 'axes.labelcolor': INK})
    fig, axs = plt.subplots(2, 2, figsize=(13, 7.5), sharex='col',
                            gridspec_kw={'height_ratios': [3, 1.1], 'hspace': 0.08, 'wspace': 0.18})
    labels = ['0-0.5', '0.5-1', '1-2', '2-3', '3-5', '5-10', '10-20', '20-50', '50-200', '>200']
    xb = np.arange(len(g))
    ax = axs[0, 0]
    ax.errorbar(xb, 100 * g['eff'], 100 * g['eff_err'], fmt='o-', color=C_RECO, lw=2, ms=6,
                label=f'reco within {R:g} mm')
    ax.errorbar(xb, 100 * g['any'], 100 * g['any_err'], fmt='s--', color=C_ANY, lw=2, ms=6,
                label='any pad fired')
    ax.axvspan(-0.5, 4.5, color='#cfe3fa', alpha=0.35, lw=0)
    ax.text(0.0, 0.03, f'shipped cut chi2 < {qa.M3_CHI2_CUT:g}', transform=ax.transAxes, color=INK2, fontsize=9)
    ax.set_ylabel('P2 efficiency [%]'); ax.set_title('differential: efficiency in M3 chi2 bins',
                                                      loc='left', color=INK, fontweight='bold')
    ax.legend(frameon=False, loc='lower left', bbox_to_anchor=(0, 0.08))
    axs[1, 0].bar(xb, g['n'], color='#9c9c9c', width=0.62)
    axs[1, 0].set_xticks(xb); axs[1, 0].set_xticklabels(labels, rotation=0, fontsize=9)
    axs[1, 0].set_xlabel('M3 track chi2_max = max(Chi2X, Chi2Y)'); axs[1, 0].set_ylabel('tracks')

    ax = axs[0, 1]
    xc = np.arange(len(cum))
    ax.errorbar(xc, 100 * cum['eff'], 100 * cum['eff_err'], fmt='o-', color=C_RECO, lw=2, ms=6,
                label=f'reco within {R:g} mm')
    ax.errorbar(xc, 100 * cum['any'], 100 * cum['any_err'], fmt='s--', color=C_ANY, lw=2, ms=6,
                label='any pad fired')
    i5 = int(np.argmin(np.abs(cum['chi2_cut'] - qa.M3_CHI2_CUT)))
    ax.axvline(i5, color=INK2, lw=1, ls=':')
    ax.annotate(f'{100*cum.eff[i5]:.1f}%', (i5, 100 * cum.eff[i5]), textcoords='offset points',
                xytext=(6, -14), color=C_RECO, fontsize=9)
    ax.set_title('cumulative: efficiency for tracks with chi2 < cut', loc='left', color=INK, fontweight='bold')
    ax.legend(frameon=False, loc='lower left')
    axs[1, 1].bar(xc, cum['n'], color='#9c9c9c', width=0.62)
    axs[1, 1].set_xticks(xc)
    axs[1, 1].set_xticklabels(['all' if c > 1e9 else f'{c:g}' for c in cum['chi2_cut']], fontsize=9)
    axs[1, 1].set_xlabel('chi2 cut'); axs[1, 1].set_ylabel('tracks kept')
    for a in axs.ravel():
        a.grid(True, axis='y', color=GRID, lw=0.8); a.set_axisbelow(True)
        for s in ('top', 'right'):
            a.spines[s].set_visible(False)
    lo = np.nanmin([100 * g['eff'].min(), 100 * cum['eff'].min()])
    for a in (axs[0, 0], axs[0, 1]):
        a.set_ylim(max(lo - 5, 0), 100)
    fig.suptitle(f'{cfg.DET_NAME} efficiency vs M3 reference-track quality — {cfg.RUN}/{cfg.SUB_RUN}  '
                 f'(N={N:,} active-area single tracks, no chi2 cut)', x=0.01, ha='left', color=INK, fontsize=12)
    fig.savefig(f'{out_dir}/eff_vs_chi2{sfx}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {out_dir}/eff_vs_chi2{sfx}.png')


if __name__ == '__main__':
    main()
