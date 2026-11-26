#!/usr/bin/env python3
"""27_loss_breakdown.py -- what makes up each efficiency-loss block?

Usage:  python3 27_loss_breakdown.py <run_key> [--no-veto-sparks]

Starts from stage 06's ray list (same active area, same R, same spark veto,
so the three top-level numbers are identical to 06's efficiency_breakdown) and
breaks the two loss blocks into their causes.

NO HIT (no P2 pad fired in the event) -- causes are attributed as an EXCESS
over the clean-area no-hit rate b, so a zone only claims the misses it adds:
    excess(zone) = n_zone * (rate_zone - b) / N
  * big insulation pillars (6.15 mm)            -- Saclay/CERN mask Gerber
  * small insulation pillars (core, r <= pillar radius)
  * pad-edge band (ray within EDGE mm of a pad-tile edge or in an inter-pad
    gap: charge sharing pushes both halves below threshold) -- Gerber tiles
  * low-response pads (core has_any < 50 %: weak gain, not dead)
  * signal below threshold: flat extrapolation of the leading-pad amplitude
    density of reconstructed events from [THR, PLATEAU] down to 0
  * unexplained remainder (M3 fake/mismatched tracks, gain fluctuations not
    caught above, ...)
FIRED BUT MIS-RECO (a pad fired, reconstructed point > R from the track):
  * resolution tail: reco within (R, 2R]
  * wrong cluster won: a pad within R of the track fired, but a larger
    cluster elsewhere (crosstalk copy, delta ray, second particle) was kept
  * only distant hits: nothing fired within R of the track -- the hits are
    elsewhere (noise, crosstalk, M3 track mismatch / another particle)

Outputs: 27_loss_breakdown/loss_breakdown<sfx>.png and .csv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa      # noqa: E402
import p2_mapping as pmap      # noqa: E402
import p2_io as p2io           # noqa: E402

# palette (dataviz reference instance, light): one hue per block
C_RECO, C_MIS, C_NOHIT = '#257de0', '#f47230', '#ec4647'
INK, INK2, GRID = '#1f1f1f', '#5c5c5c', '#e6e6e6'


def seg_dist(p, a, b):
    ab = b - a
    t = np.clip(((p - a) @ ab) / (ab @ ab), 0, 1)
    return np.hypot(*(p - (a + t[:, None] * ab)).T)


def signed_edge_distance(q, verts, tree_c):
    """Signed distance [mm] of pad-frame points q to the nearest pad-tile edge
    (>0 inside a tile, <0 in a gap/outside)."""
    _, cand = tree_c.query(q, k=4)
    ins = np.zeros(len(q), bool)
    din = np.zeros(len(q))
    dmin = np.full(len(q), np.inf)
    for k in range(cand.shape[1]):
        for i_ in np.unique(cand[:, k]):
            sel = np.where(cand[:, k] == i_)[0]
            v = verts[i_]
            cross, dist = [], np.full(len(sel), np.inf)
            for e in range(4):
                a, b = v[e], v[(e + 1) % 4]
                w = q[sel] - a
                cross.append((b - a)[0] * w[:, 1] - (b - a)[1] * w[:, 0])
                dist = np.minimum(dist, seg_dist(q[sel], a, b))
            cross = np.array(cross)
            inside = (cross >= 0).all(0) | (cross <= 0).all(0)
            din[sel[inside]] = dist[inside]
            ins[sel[inside]] = True
            dmin[sel] = np.minimum(dmin[sel], dist)
    return np.where(ins, din, -dmin)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument('--edge', type=float, default=1.0,
                    help='pad-edge band width [mm] (default 1)')
    ap.add_argument('--thr', type=float, default=40.0,
                    help='effective hit threshold [ADC] for the below-threshold extrapolation')
    ap.add_argument('--plateau', type=float, default=150.0,
                    help='upper end [ADC] of the flat density window above --thr')
    args = ap.parse_args()
    cfg = qa.get_config(args.run_key)
    out_dir = cfg.out_dir('27_loss_breakdown')
    sfx = cfg.product_suffix(args.veto_sparks)
    R = cfg.MATCH_R

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
                                  det_name=cfg.DET_NAME, drop_connectors=cfg.DEAD_CONNECTORS,
                                  drop_channels=cfg.DEAD_CHANNELS,
                                  cluster_link=cfg.CLUSTER_LINK_MM,
                                  strategy_overrides=cfg.STRATEGY_OVERRIDES)
    padc = ct[ct['mapped']].drop_duplicates('channel_id').reset_index(drop=True)
    ed = cfg.out_dir('06_efficiency')
    fp = pd.read_csv(f'{ed}/pad_footprint{sfx}.csv')
    if len(fp) != len(padc):
        raise SystemExit('06 pad footprint does not match the channel table -- rerun 06 first')
    # pad -> M3 affine (exact on the footprint rows) and its inverse
    P = np.c_[padc.pad_cx, padc.pad_cy, np.ones(len(padc))]
    coef, *_ = np.linalg.lstsq(P, fp[['x', 'y']].to_numpy(), rcond=None)
    M, b0 = coef[:2].T, coef[2]
    Minv = np.linalg.inv(M)
    scale = np.sqrt(abs(np.linalg.det(M)))

    rays = pd.read_csv(f'{ed}/ray_hit_miss_list{sfx}.csv')
    for c in ('in_active', 'has_any', 'within'):
        rays[c] = rays[c].astype(str).str.lower().isin(('true', '1'))
    ra = rays[rays['in_active']].copy().reset_index(drop=True)
    N = len(ra)
    q_pad = (ra[['x', 'y']].to_numpy() - b0) @ Minv.T          # rays in the pad frame

    # ---------------- zones for the NO-HIT attribution ---------------------- #
    zones = {}
    pil = pmap.load_pillars(cfg.MASK_GBR_PATH)
    if len(pil):
        big, small = pil[pil['big']], pil[~pil['big']]
        if len(big):
            d_big = np.min(np.hypot(q_pad[:, 0, None] - big['x'].to_numpy()[None],
                                    q_pad[:, 1, None] - big['y'].to_numpy()[None]), axis=1)
            zones['big insulation pillars'] = d_big <= big['r'].median()
        if len(small):
            d_s, _ = cKDTree(small[['x', 'y']].to_numpy()).query(q_pad)
            zones['small insulation pillars'] = d_s <= small['r'].median()
    if pmap.has_tile_geometry(ct):
        tpads, verts = pmap.pad_tiles(ct)
        verts = np.asarray(verts, float)
        sd = signed_edge_distance(q_pad, verts, cKDTree(verts.mean(1)))
        zones[f'pad edges / gaps (< {args.edge:g} mm)'] = sd < args.edge
    # low-response pads: core rays (<= 4 mm from a pad centre) with has_any < 50 %
    tree = cKDTree(fp[['x', 'y']].to_numpy())
    dpad, ipad = tree.query(ra[['x', 'y']].to_numpy())
    core = pd.DataFrame({'pad': ipad, 'any': ra['has_any']})[dpad <= 4.0]
    g = core.groupby('pad')['any'].agg(['size', 'mean'])
    lowpads = set(g.index[(g['size'] >= 8) & (g['mean'] < 0.5)])
    zones['low-response pads'] = np.isin(ipad, list(lowpads))

    anyzone = np.zeros(N, bool)
    for z in zones.values():
        anyzone |= z
    nohit = ~ra['has_any'].to_numpy()
    b = nohit[~anyzone].mean()
    parts_nohit = []
    claimed = np.zeros(N, bool)
    for name, z in zones.items():
        zz = z & ~claimed                      # first zone wins, no double counting
        claimed |= z
        exc = zz.sum() * (nohit[zz].mean() - b) / N if zz.any() else 0.0
        parts_nohit.append((name, max(exc, 0.0), int(zz.sum())))

    # ---------------- hits (needed for threshold + mis-reco split) ---------- #
    mpos = ct[ct['mapped']].set_index(['feu', 'channel'])[['pad_cx', 'pad_cy']]
    want = set(ra['eventId'])
    hp = []
    for df in p2io.iter_hits(cfg.combined_hits_dir, ['eventId', 'feu', 'channel', 'amplitude'],
                             ct.attrs['feus'], progress=False, t_max_h=cfg.T_MAX_H,
                             t_min_h=cfg.T_MIN_H, min_amp=cfg.MIN_AMP):
        df = pd.DataFrame(df) if not isinstance(df, pd.DataFrame) else df
        df = df[df['eventId'].isin(want)]
        idx = pd.MultiIndex.from_arrays([df['feu'].astype(int), df['channel'].astype(int)])
        ok = idx.isin(mpos.index)
        xy = mpos.reindex(idx[ok]).to_numpy()
        hp.append(pd.DataFrame({'eventId': df['eventId'].to_numpy()[ok],
                                'a': df['amplitude'].to_numpy()[ok],
                                'px': xy[:, 0], 'py': xy[:, 1]}))
    h = pd.concat(hp, ignore_index=True)
    hm = np.c_[h['px'], h['py']] @ M.T + b0
    h['mx'], h['my'] = hm[:, 0], hm[:, 1]

    # below threshold: leading-pad amplitude of reconstructed events
    good = set(ra.loc[ra['within'], 'eventId'])
    lead = h[h['eventId'].isin(good)].groupby('eventId')['a'].max()
    rho = ((lead >= args.thr) & (lead < args.plateau)).sum() / (args.plateau - args.thr)
    n_below_obs = (lead < args.thr).sum()
    missing = max(rho * args.thr - n_below_obs, 0.0)
    # those muons would have been reconstructed -> scale by the reco share
    f_bt = missing / N
    parts_nohit.append((f'signal below threshold (< {args.thr:g} ADC)', f_bt, int(round(missing))))
    nohit_tot = nohit.mean()
    explained = sum(p[1] for p in parts_nohit)
    parts_nohit.append(('unexplained (flat vs M3 chi2 -> not reference quality; see 28)', max(nohit_tot - explained, 0.0), 0))

    # ---------------- FIRED BUT MIS-RECO split ------------------------------ #
    fm = ra[ra['has_any'] & ~ra['within']]
    rx = ra.set_index('eventId')
    hh = h[h['eventId'].isin(fm['eventId'])].copy()
    hh = hh.join(rx[['x', 'y']], on='eventId')
    hh['d'] = np.hypot(hh['mx'] - hh['x'], hh['my'] - hh['y'])
    near = hh.groupby('eventId')['d'].min() <= R
    near = near.reindex(fm['eventId']).fillna(False).to_numpy()
    resid = fm['resid'].to_numpy()
    tail = resid <= 2 * R
    p_tail = (tail).sum() / N
    p_wrong = (~tail & near).sum() / N
    p_far = (~tail & ~near).sum() / N
    parts_mis = [(f'resolution tail ({R:g}-{2*R:g} mm)', p_tail, int(tail.sum())),
                 ('wrong cluster won (signal at track, larger hit elsewhere)', p_wrong, int((~tail & near).sum())),
                 ('only distant hits (nothing fired at the track)', p_far, int((~tail & ~near).sum()))]

    within_tot = ra['within'].mean()
    mis_tot = len(fm) / N

    # ---------------- table ------------------------------------------------- #
    rows = [('reconstructed', '', within_tot)]
    rows += [('fired but mis-reco', n, v) for n, v, _ in parts_mis]
    rows += [('no hit', n, v) for n, v, _ in parts_nohit]
    tab = pd.DataFrame(rows, columns=['block', 'cause', 'frac_of_active_rays'])
    tab['percent'] = (100 * tab['frac_of_active_rays']).round(3)
    tab.to_csv(f'{out_dir}/loss_breakdown{sfx}.csv', index=False)
    print(f'{cfg.DET_TAG} {cfg.RUN}/{cfg.SUB_RUN}: N={N:,}  reco {100*within_tot:.2f}%  '
          f'mis-reco {100*mis_tot:.2f}%  no-hit {100*nohit_tot:.2f}%  (clean-area no-hit {100*b:.2f}%)')
    print(tab.to_string(index=False))

    # ---------------- figure ------------------------------------------------ #
    plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                         'xtick.color': INK2, 'ytick.color': INK2})
    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.9], height_ratios=[len(parts_mis) + 1.2, len(parts_nohit) + 1.2],
                          wspace=0.12, hspace=0.35)
    ax0 = fig.add_subplot(gs[:, 0])
    axm = fig.add_subplot(gs[0, 1])
    axn = fig.add_subplot(gs[1, 1])

    cats = [f'reconstructed\n(|r| ≤ {R:g} mm)', 'fired but\nmis-reco', 'no hit\n(dead miss)']
    vals = [100 * within_tot, 100 * mis_tot, 100 * nohit_tot]
    cols = [C_RECO, C_MIS, C_NOHIT]
    bars = ax0.bar(cats, vals, color=cols, width=0.62, edgecolor='white', linewidth=2)
    for bb, v in zip(bars, vals):
        ax0.text(bb.get_x() + bb.get_width() / 2, v + 1.0, f'{v:.1f}%', ha='center', va='bottom',
                 color=INK, fontsize=11, fontweight='bold')
    ax0.set_ylabel('% of active-area crossing muons')
    ax0.set_ylim(0, 105)
    ax0.grid(True, axis='y', color=GRID, lw=0.8); ax0.set_axisbelow(True)
    for s in ('top', 'right'):
        ax0.spines[s].set_visible(False)

    def breakout(ax, parts, total, color, title):
        parts = sorted(parts, key=lambda p: -p[1])
        names = [p[0] for p in parts]
        v = np.array([100 * p[1] for p in parts])
        y = np.arange(len(parts))[::-1] * 1.6           # room for a name above each bar
        ax.barh(y, v, color=color, height=0.62, edgecolor='white', linewidth=2)
        xmax = max(v.max() * 1.35, 0.5)
        for yi, vi, nm in zip(y, v, names):
            share = 100 * vi / (100 * total) if total else 0
            ax.text(0, yi + 0.42, nm, va='bottom', ha='left', color=INK, fontsize=9.5)
            ax.text(vi + xmax * 0.012, yi, f'{vi:.2f}%  ({share:.0f}% of block)', va='center',
                    color=INK2, fontsize=9)
        ax.set_yticks([])
        ax.set_ylim(-0.6, y.max() + 1.25)
        ax.set_xlim(0, xmax)
        ax.set_xlabel('% of active-area crossing muons', color=INK2)
        ax.set_title(title, loc='left', color=INK, fontsize=11, fontweight='bold')
        ax.grid(True, axis='x', color=GRID, lw=0.8); ax.set_axisbelow(True)
        for s_ in ('top', 'right', 'left'):
            ax.spines[s_].set_visible(False)

    breakout(axm, parts_mis, mis_tot, C_MIS, f'fired but mis-reco  ({100*mis_tot:.2f}%)')
    breakout(axn, parts_nohit, nohit_tot, C_NOHIT, f'no hit  ({100*nohit_tot:.2f}%)')

    # connector lines from the loss bars to their breakouts
    for bb, axb, color in ((bars[1], axm, C_MIS), (bars[2], axn, C_NOHIT)):
        x_r = bb.get_x() + bb.get_width()
        for yb, ya in ((bb.get_height(), 1.0), (0.0, 0.0)):
            con = ConnectionPatch(xyA=(x_r, yb), coordsA=ax0.transData,
                                  xyB=(0, ya), coordsB=axb.transAxes,
                                  color=color, lw=1.0, alpha=0.55, ls='--')
            fig.add_artist(con)

    fig.suptitle(f'{cfg.DET_NAME} — where do the crossing muons go, and why are some lost?   '
                 f'{cfg.RUN}/{cfg.SUB_RUN}  (N={N:,})', color=INK, fontsize=12, x=0.01, ha='left')
    fig.text(0.01, 0.005,
             f'No-hit causes are EXCESS over the clean-area no-hit rate ({100*b:.2f}%); first matching zone wins. '
             f'Below-threshold: flat leading-pad amplitude density [{args.thr:g},{args.plateau:g}] ADC '
             f'extrapolated to 0. Pillars/tiles from {os.path.basename(cfg.MASK_GBR_PATH)} + pad-map Gerber.',
             color=INK2, fontsize=8)
    fig.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.10)
    fig.savefig(f'{out_dir}/loss_breakdown{sfx}.png', dpi=150)
    plt.close(fig)
    print(f'  saved {out_dir}/loss_breakdown{sfx}.png')


if __name__ == '__main__':
    main()
