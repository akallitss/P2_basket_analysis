#!/usr/bin/env python3
"""33_reference_vs_p2_blame.py -- of the 'only distant hits' (and the resolution
tail), how much would remain with a PERFECT reference track?

Usage:  python3 33_reference_vs_p2_blame.py <run_key> [--no-veto-sparks]

Only-distant = a P2 cluster fired but none within MATCH_R and none within 2R
of the track. Two populations:
  M3 fault : P2 saw the muon, the reference track points elsewhere
  P2 fault : P2 missed the muon AND an unrelated in-time particle fired P2
             elsewhere (it would have been a no-hit otherwise)
Three independent estimates of the P2-fault share:
  1. chi2 extrapolation: rate in M3 chi2 bins, linear fit -> chi2 = 0
  2. coincidence model : P(P2 miss) x P(extra in-time particle), with
     P(extra) measured on tracks that miss the detector (> 40 mm outside)
  3. distance shape    : the random-coincidence distance distribution is
     built by event mixing (each track paired with the clusters of events
     whose track missed P2), normalised in a sideband (> SIDE mm, where
     mis-pointing does not reach) and integrated
Outputs: 33_reference_blame/blame<sfx>.png + csv
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

INK, INK2, GRID = '#1f1f1f', '#5c5c5c', '#e6e6e6'
C_M3, C_P2, C_DATA = '#257de0', '#ec4647', '#1f1f1f'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument('--side', type=float, default=150.0, help='sideband start [mm] for method 3')
    ap.add_argument('--mix', type=int, default=20, help='random partner events per track (method 3)')
    args = ap.parse_args()
    cfg = qa.get_config(args.run_key)
    out_dir = cfg.out_dir('33_reference_blame')
    sfx = cfg.product_suffix(args.veto_sparks)
    R = cfg.MATCH_R

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
                                  det_name=cfg.DET_NAME, drop_connectors=cfg.DEAD_CONNECTORS,
                                  drop_channels=cfg.DEAD_CHANNELS, cluster_link=cfg.CLUSTER_LINK_MM,
                                  dead_guard=cfg.DEAD_GUARD_MM, strategy_overrides=cfg.STRATEGY_OVERRIDES)
    padc = ct[ct['mapped']].drop_duplicates('channel_id').reset_index(drop=True)
    ed = cfg.out_dir('06_efficiency')
    fp = pd.read_csv(f'{ed}/pad_footprint{sfx}.csv')
    coef, *_ = np.linalg.lstsq(np.c_[padc.pad_cx, padc.pad_cy, np.ones(len(padc))],
                               fp[['x', 'y']].to_numpy(), rcond=None)

    class _T:
        @staticmethod
        def apply(x, y):
            q = np.c_[np.asarray(x), np.asarray(y), np.ones(len(x))] @ coef
            return q[:, 0], q[:, 1]

    rays = pd.read_csv(f'{ed}/ray_hit_miss_list{sfx}.csv')
    for c in ('in_active', 'has_any', 'within'):
        rays[c] = rays[c].astype(str).str.lower().isin(('true', '1'))
    act = rays[rays['in_active']].copy()
    N = len(act)
    act['only_distant'] = act['has_any'] & ~act['within'] & (act['resid'] > 2 * R)
    act['rtail'] = act['has_any'] & ~act['within'] & (act['resid'] <= 2 * R)
    act['nohit'] = ~act['has_any']

    m3 = M3RefTracking(os.path.join(cfg.m3_tracking_dir, ''), single_track=True,
                       chi2_cut=qa.M3_CHI2_CUT, min_nclus=qa.M3_MIN_NCLUS)
    rd = m3.ray_data
    c2 = pd.Series(np.maximum(ak.to_numpy(ak.ravel(rd['Chi2X'])), ak.to_numpy(ak.ravel(rd['Chi2Y']))),
                   index=ak.to_numpy(ak.ravel(rd['evn'])).astype(np.int64))
    act['chi2'] = c2[~c2.index.duplicated()].reindex(act['eventId']).to_numpy()

    od_rate, tail_rate, nohit_rate = act.only_distant.mean(), act['rtail'].mean(), act.nohit.mean()
    print(f'{cfg.DET_TAG}: N={N:,}  only-distant {100*od_rate:.3f}%  tail {100*tail_rate:.3f}%  no-hit {100*nohit_rate:.2f}%')

    # ---------------- method 1: chi2 extrapolation ----------------
    edges = np.array([0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5])
    act['cb'] = pd.cut(act['chi2'], edges, right=False)
    g = act.groupby('cb', observed=False).agg(n=('only_distant', 'size'), od=('only_distant', 'mean'),
                                             tl=('rtail', 'mean'), c=('chi2', 'mean')).dropna()
    g = g[g.n >= 200]
    res1 = {}
    for colname in ('od', 'tl'):
        w = g.n / np.clip(g[colname] * (1 - g[colname]), 1e-6, None)
        A = np.c_[np.ones(len(g)), g.c]
        W = np.diag(w)
        cov = np.linalg.inv(A.T @ W @ A)
        beta = cov @ A.T @ W @ g[colname].to_numpy()
        res1[colname] = (beta[0], np.sqrt(cov[0, 0]), beta[1])
    print(f'method 1 (chi2 -> 0): only-distant {100*res1["od"][0]:.3f} +- {100*res1["od"][1]:.3f}%, '
          f'tail {100*res1["tl"][0]:.3f} +- {100*res1["tl"][1]:.3f}%')

    # ---------------- clusters (for methods 2 and 3) ----------------
    p2, hit_events, clusters = pa.load_p2_centroids(cfg.combined_hits_dir, ct, min_amp=cfg.MIN_AMP,
                                                    t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H,
                                                    drop_pads=p2io.drop_pads_for(cfg, ct), with_clusters=True)
    if args.veto_sparks:
        bad = ps.SparkVeto.from_cfg(cfg).vetoed_ids_from_hits(cfg.combined_hits_dir, ct.attrs['feus'],
                                                              min_amp=cfg.MIN_AMP)
        clusters = clusters[~clusters['eventId'].isin(bad)]
    cx, cy = _T.apply(clusters['x_pad'].to_numpy(), clusters['y_pad'].to_numpy())
    clusters = clusters.assign(mx=cx, my=cy)

    # ---------------- method 2: coincidence model ----------------
    dfp, _ = cKDTree(fp[['x', 'y']].to_numpy()).query(rays[['x', 'y']].to_numpy())
    outside = rays[dfp > 40]
    out_ev = set(outside['eventId'])
    p_extra = outside['eventId'].isin(set(clusters['eventId'])).mean()
    # P(muon missed) = nohit + p2part, p2part = P(miss) * p_extra  ->  p2part = nohit * p/(1-p)
    p2part_2 = nohit_rate * p_extra / (1 - p_extra)
    print(f'method 2 (coincidence): P(extra in-time particle) = {100*p_extra:.2f}% (n={len(outside):,} tracks '
          f'outside det), P2-fault only-distant = {100*p2part_2:.3f}%')

    # ---------------- method 3: distance shape, event mixing ----------------
    od = act[act.only_distant]
    extra = clusters[clusters['eventId'].isin(out_ev)]
    ex_by_ev = [g_[['mx', 'my']].to_numpy() for _, g_ in extra.groupby('eventId')]
    rng = np.random.default_rng(3)
    bins = np.arange(2 * R, 460, 20.0)
    h_data, _ = np.histogram(od['resid'].clip(upper=bins[-1] - 1), bins)
    # random template: track of each only-distant event + clusters of a random 'extra' event
    mix = []
    for x, y in zip(od['x'], od['y']):
        for k in rng.integers(0, len(ex_by_ev), args.mix):
            c = ex_by_ev[k]
            mix.append(np.hypot(c[:, 0] - x, c[:, 1] - y).min())
    mix = np.asarray(mix)
    mix = mix[mix > 2 * R]                    # same class definition (nothing within 2R)
    h_mix, _ = np.histogram(np.minimum(mix, bins[-1] - 1), bins)
    side = bins[:-1] >= args.side
    scale = h_data[side].sum() / max(h_mix[side].sum(), 1)
    n_random = scale * h_mix.sum()
    n_random_err = np.sqrt(h_data[side].sum()) / max(h_data[side].sum(), 1) * n_random
    p2part_3 = n_random / N
    print(f'method 3 (distance shape): sideband > {args.side:g} mm holds {h_data[side].sum()} of {len(od)} events; '
          f'random-coincidence component {n_random:.0f} +- {n_random_err:.0f} events = {100*p2part_3:.3f}% of tracks')

    est = pd.DataFrame([
        dict(method='1 chi2 -> 0', p2=res1['od'][0], err=res1['od'][1]),
        dict(method='2 coincidence model', p2=p2part_2, err=p2part_2 * 0.2),
        dict(method='3 distance shape (mixing)', p2=p2part_3, err=n_random_err / N)])
    est['m3'] = od_rate - est['p2']
    est.to_csv(f'{out_dir}/blame_estimates{sfx}.csv', index=False)
    print(est.assign(p2=lambda t: 100 * t.p2, err=lambda t: 100 * t.err, m3=lambda t: 100 * t.m3).round(3).to_string(index=False))
    # the three methods have different systematics (method 2's error is only a
    # nominal 20%), so quote their mean and their spread, not a weighted average
    comb = est.p2.mean(); comb_e = (est.p2.max() - est.p2.min()) / 2
    eff_now = act['within'].mean()
    tail_m3 = tail_rate - max(res1['tl'][0], 0)
    print(f'\nP2-fault share of only-distant: {100*comb:.3f} +- {100*comb_e:.3f}% of tracks (mean +- half-spread of the 3 methods) '
          f'(M3 fault {100*(od_rate-comb):.3f}%); tail: M3 ~{100*tail_m3:.3f}%')
    print(f'efficiency now {100*eff_now:.2f}% -> with a perfect reference ~{100*(eff_now + od_rate - comb + tail_m3):.2f}%')

    # ---------------- figure ----------------
    plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.4), gridspec_kw={'width_ratios': [1.1, 1.3, 1]})
    ax = axs[0]
    e_od = np.sqrt(g.od * (1 - g.od) / g.n); e_tl = np.sqrt(g.tl * (1 - g.tl) / g.n)
    ax.errorbar(g.c, 100 * g.od, 100 * e_od, fmt='o', color=C_DATA, ms=6, label='only distant hits')
    ax.errorbar(g.c, 100 * g.tl, 100 * e_tl, fmt='s', color='#9c9c9c', ms=5, label=f'resolution tail ({R:g}-{2*R:g} mm)')
    xx = np.linspace(0, 5, 50)
    for colname, c_ in (('od', C_DATA), ('tl', '#9c9c9c')):
        b0, eb0, b1 = res1[colname]
        ax.plot(xx, 100 * (b0 + b1 * xx), color=c_, lw=1.5, ls='--')
        ax.errorbar([0], [100 * b0], [100 * eb0], fmt='*', color=C_P2, ms=13)
    ax.annotate(f'chi2 -> 0: {100*res1["od"][0]:.2f}%', (0, 100 * res1['od'][0]), xytext=(18, 22),
                textcoords='offset points', color=C_P2, fontsize=9, arrowprops=dict(arrowstyle='-', color=C_P2))
    ax.set_xlabel('M3 track chi2_max'); ax.set_ylabel('% of tracks')
    ax.set_title('1  Extrapolate to a perfect track', loc='left', fontweight='bold', color=INK)
    ax.legend(frameon=False, fontsize=9)
    ax = axs[1]
    cen = 0.5 * (bins[1:] + bins[:-1])
    ax.step(cen, h_data, where='mid', color=C_DATA, lw=2, label=f'only-distant events (n={len(od)})')
    ax.fill_between(cen, scale * h_mix, step='mid', color=C_P2, alpha=0.35,
                    label=f'random coincidence (event mixing), {n_random:.0f} events')
    ax.axvline(args.side, color=INK2, ls=':', lw=1)
    ax.text(args.side, ax.get_ylim()[1] * 0.9, ' sideband (normalisation)', fontsize=8.5, color=INK2)
    ax.set_xlabel('distance from the track to the nearest P2 cluster [mm]'); ax.set_ylabel('events / 20 mm')
    ax.set_title('3  Near the track = mis-pointing; flat = coincidence', loc='left', fontweight='bold', color=INK)
    ax.legend(frameon=False, fontsize=9)
    ax = axs[2]
    yy = np.arange(len(est))[::-1]
    ax.barh(yy, 100 * est.m3, color=C_M3, height=0.55, label='M3 (reference) fault')
    ax.barh(yy, 100 * est.p2, left=100 * est.m3, color=C_P2, height=0.55, label='P2 fault')
    for yv, m_, p_ in zip(yy, est.m3, est.p2):
        ax.text(100 * (m_ + p_) + 0.01, yv, f'P2 {100*p_:.2f}%', va='center', fontsize=9, color=INK)
    ax.set_yticks(yy); ax.set_yticklabels(est.method)
    ax.set_xlabel(f'% of tracks (only-distant total {100*od_rate:.2f}%)')
    ax.set_xlim(0, 100 * od_rate * 1.35)
    ax.set_title('Who is to blame?', loc='left', fontweight='bold', color=INK)
    ax.legend(frameon=False, fontsize=9, loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=2)
    for a_ in axs:
        a_.grid(True, color=GRID, lw=0.8); a_.set_axisbelow(True)
        for s_ in ('top', 'right'):
            a_.spines[s_].set_visible(False)
    fig.suptitle(f'{cfg.DET_NAME} only-distant hits: reference (M3) vs P2 — {cfg.RUN}/{cfg.SUB_RUN}   '
                 f'P2 share {100*est.p2.min():.2f}-{100*est.p2.max():.2f}% of tracks (3 methods)', x=0.01, ha='left',
                 color=INK, fontsize=12)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/blame{sfx}.png', dpi=140)
    plt.close(fig)
    print(f'  saved {out_dir}/blame{sfx}.png')


if __name__ == '__main__':
    main()
