#!/usr/bin/env python3
"""34_clustering_matching.py -- pointing resolution, pad-aware clustering and
track matching with accidental subtraction.

Answers: (1) what is the real pointing resolution at the P2 plane once the
alignment is not biased (35: the plain Procrustes fit shrank the scale by 2.3 %
because 4.7 % of the fit pairs are not the muon)? (2) does a clustering/matching
built for a pad detector do better than "connected components at 20 mm +
centroid within 20 mm"?

Alignment: robust (trimmed) 2D transform at the slope-fitted plane z (see 35).
Denominator: rays at least FID_MM inside the live-pad area (signed distance to
the union of live tiles; dead/shorted/unplugged tiles and the outer edge are
both boundaries). 35 shows the miss rate is flat from 2 mm inside.

Pointing resolution, two independent ways:
  a) single-pad clusters: the track-minus-pad-centre residual along each pad
     axis is box(pad width) (x) Gauss(sigma); fit sigma (per axis, vs angle).
  b) the boundary step: miss(sd) = b + (1-b) * Phi(-(sd - s0)/sigma).

Clustering:  CC20  connected components, centres < 20 mm (stage 06 today)
             ADJ   pads whose TILES touch (gap < 1 mm) -- the pad-detector rule
Seed threshold S: a cluster counts only if its largest pad has amp >= S.
Matching:    CEN   cluster centroid within 20 mm of the track (06 today)
             TILE  track within k*sigma_point of any pad tile of the cluster
Accidentals: each ray is also matched against the clusters of a DIFFERENT
trigger (event mixing, shift of 1000 events); eff_corr = (e - a) / (1 - a).

Usage: python3 34_clustering_matching.py [run_key]
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_study_cache as sc     # noqa: E402
import p2_study_common as cm    # noqa: E402
import p2_align as pa           # noqa: E402

INK2, GRID = '#5c5c5c', '#e6e6e6'
C1, C2, C3, C4 = '#9c9c9c', '#257de0', '#ec4647', '#1baf7a'
FID_MM = 2.0
SEEDS = [0, 50, 100, 150, 200, 300, 500, 800]


# ------------------------------------------------------------------ geometry
def tile_adjacency(tiles, gap=1.0):
    """Pairs of pads (row indices) whose tiles are closer than `gap` mm."""
    cen = tiles.mean(1)
    pairs = cKDTree(cen).query_pairs(25.0, output_type='ndarray')
    a, b = pairs[:, 0], pairs[:, 1]
    d = np.full(len(pairs), np.inf)
    for k in range(4):   # vertex of one tile to the other tile (both ways)
        d = np.minimum(d, pt_quad_dist(tiles[a, k], tiles[b]))
        d = np.minimum(d, pt_quad_dist(tiles[b, k], tiles[a]))
    return pairs[d < gap]


def pt_quad_dist(p, q):
    """Distance of points p (n,2) to convex quads q (n,4,2), row-wise; 0 inside."""
    dmin = np.full(len(p), np.inf); sgn = []
    for e in range(4):
        a, b = q[:, e], q[:, (e + 1) % 4]
        ab = b - a; ap = p - a
        t = np.clip(np.einsum('ij,ij->i', ap, ab) / np.einsum('ij,ij->i', ab, ab), 0, 1)
        dmin = np.minimum(dmin, np.hypot(*(p - (a + t[:, None] * ab)).T))
        sgn.append(ab[:, 0] * ap[:, 1] - ab[:, 1] * ap[:, 0])
    sgn = np.array(sgn)
    inside = (sgn >= 0).all(0) | (sgn <= 0).all(0)
    return np.where(inside, 0.0, dmin)


def signed_distance(xy, live_q, dead_q, res=0.2):
    from PIL import Image, ImageDraw
    from scipy.ndimage import distance_transform_edt
    lo = live_q.reshape(-1, 2).min(0) - 30; hi = live_q.reshape(-1, 2).max(0) + 30
    W, H = ((hi - lo) / res).astype(int) + 1
    im = Image.new('1', (W, H), 0); dr = ImageDraw.Draw(im)
    for q in live_q:
        dr.polygon([tuple(p) for p in (q - lo) / res], fill=1)
    for q in dead_q:
        dr.polygon([tuple(p) for p in (q - lo) / res], fill=0)
    m = np.array(im, bool)
    sd = np.where(m, distance_transform_edt(m) * res, -distance_transform_edt(~m) * res)
    ij = ((xy - lo) / res).astype(int)
    ok = (ij[:, 0] >= 0) & (ij[:, 0] < W) & (ij[:, 1] >= 0) & (ij[:, 1] < H)
    out = np.full(len(xy), -99.0); out[ok] = sd[ij[ok, 1], ij[ok, 0]]
    return out


# ------------------------------------------------------------------ alignment
def robust_alignment(rays, clu, g, z_start, n_iter=6, trim=15.0):
    """Alternate (robust 2D Procrustes on trimmed pairs) <-> (plane z from the
    residual-vs-slope regression). Returns (xy, T, z)."""
    ev = rays['eventId'].to_numpy()
    tx, ty, _ = cm.ray_dirs(rays)
    T, z = cm.make_T(g), z_start
    cxy_pad = clu[['x_pad', 'y_pad']].to_numpy()
    for it in range(n_iter):
        x, y = pa.project_to_z(rays, z); xy = np.column_stack([x, y])
        d, j = cm.nearest(ev, xy, cm.apply_T(T, cxy_pad), clu['eventId'].to_numpy())
        k = d < (trim if it else 40.0)
        T = pa.fit_transform(xy[k, 0], xy[k, 1], cxy_pad[j[k], 0], cxy_pad[j[k], 1])
        r = cm.apply_T(T, cxy_pad[j[k]]) - xy[k]
        tt = np.concatenate([tx[k], ty[k]]); rr = np.concatenate([r[:, 0], r[:, 1]])
        m = np.abs(rr) < trim
        z += float((tt[m] * rr[m]).sum() / (tt[m] ** 2).sum())
    x, y = pa.project_to_z(rays, z)
    return np.column_stack([x, y]), T, z


# ------------------------------------------------------------------ clustering
def cluster_labels(hits, mode, adj_pairs=None, pad_row=None, link=20.0):
    """Per-hit global cluster label. hits must be sorted by eventId."""
    lab = np.empty(len(hits), np.int64)
    ev = hits['eventId'].to_numpy(); xy = hits[['x_pad', 'y_pad']].to_numpy()
    rows = hits['row'].to_numpy()
    adj = None
    if mode == 'ADJ':
        adj = set(map(tuple, np.sort(adj_pairs, 1)))
    starts = np.r_[0, np.flatnonzero(np.diff(ev)) + 1, len(ev)]
    nxt = 0
    for s, e in zip(starts[:-1], starts[1:]):
        n = e - s
        if n == 1:
            lab[s] = nxt; nxt += 1; continue
        if mode == 'CC20':
            pr = cKDTree(xy[s:e]).query_pairs(link, output_type='ndarray')
        else:
            rr = rows[s:e]
            i, j = np.triu_indices(n, 1)
            keep = [(min(rr[a], rr[b]), max(rr[a], rr[b])) in adj for a, b in zip(i, j)]
            pr = np.column_stack([i, j])[np.asarray(keep, bool)] if n > 1 else np.empty((0, 2), int)
        G = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(n, n)) if len(pr) \
            else coo_matrix((n, n))
        nc, l_ = connected_components(G, directed=False)
        lab[s:e] = l_ + nxt; nxt += nc
    return lab


def cluster_table(hits, lab):
    h = hits.assign(lab=lab, wx=hits.amp * hits.x_pad, wy=hits.amp * hits.y_pad)
    c = h.groupby('lab').agg(eventId=('eventId', 'first'), q=('amp', 'sum'), qmax=('amp', 'max'),
                             n_pad=('amp', 'size'), wx=('wx', 'sum'), wy=('wy', 'sum'),
                             xm=('x_pad', 'mean'), ym=('y_pad', 'mean'))
    qs = c.q.where(c.q > 0)
    c['x_pad'] = (c.wx / qs).fillna(c.xm); c['y_pad'] = (c.wy / qs).fillna(c.ym)
    return c.reset_index()[['lab', 'eventId', 'x_pad', 'y_pad', 'q', 'qmax', 'n_pad']]


# ------------------------------------------------------------------ matching
def match_distances(ray_ev, xy, hits, lab, clu, T, tiles_m3, ev_shift=0):
    """Per ray and per seed-threshold-independent cluster: returns a DataFrame of
    (ray i, cluster lab, centroid distance, tile distance, qmax). ev_shift != 0
    pairs ray k with the hits of the event at index k+shift (event mixing)."""
    ue = np.unique(ray_ev)
    if ev_shift:
        mp = dict(zip(ue, np.roll(ue, -ev_shift)))
        key = np.array([mp[e] for e in ray_ev])
    else:
        key = ray_ev
    t = pd.DataFrame({'eventId': key, 'i': np.arange(len(ray_ev))})
    hh = pd.DataFrame({'eventId': hits['eventId'].to_numpy(), 'lab': lab, 'row': hits['row'].to_numpy()})
    m = t.merge(hh, on='eventId')
    dt = pt_quad_dist(xy[m.i.to_numpy()], tiles_m3[m.row.to_numpy()])
    m['dtile'] = dt
    per = m.groupby(['i', 'lab'], as_index=False)['dtile'].min()
    cc = clu.set_index('lab')
    cxy = cm.apply_T(T, cc.loc[per.lab, ['x_pad', 'y_pad']].to_numpy())
    per['dcen'] = np.hypot(cxy[:, 0] - xy[per.i, 0], cxy[:, 1] - xy[per.i, 1])
    per['qmax'] = cc.loc[per.lab, 'qmax'].to_numpy()
    per['n_pad'] = cc.loc[per.lab, 'n_pad'].to_numpy()
    return per


def eff_from(per, n_rays, fid, seed, kind, win):
    p = per[per.qmax >= seed]
    col = 'dcen' if kind == 'CEN' else 'dtile'
    ok = np.zeros(n_rays, bool)
    ok[p.loc[p[col] <= win, 'i'].to_numpy()] = True
    return ok[fid]


# ------------------------------------------------------------------ resolution fits
def box_gauss(x, A, mu, w, s, b):
    return A * 0.5 * (erf((x - mu + w / 2) / (np.sqrt(2) * s)) - erf((x - mu - w / 2) / (np.sqrt(2) * s))) + b


def fit_box(res, w0):
    h, e = np.histogram(res, bins=np.arange(-20, 20.01, 0.5))
    x = 0.5 * (e[1:] + e[:-1])
    p, cv = curve_fit(box_gauss, x, h, p0=[h.max(), 0, w0, 1.0, 0], sigma=np.sqrt(h + 1),
                      bounds=([0, -5, 5, 0.05, 0], [np.inf, 5, 20, 8, np.inf]))
    return p, np.sqrt(np.diag(cv)), (x, h)


def step(sd, b, s0, sig):
    return b + (1 - b) * 0.5 * (1 - erf((sd - s0) / (np.sqrt(2) * sig)))


def main():
    run_key = sys.argv[1] if len(sys.argv) > 1 else 'det2_long_repos'
    hits, rays, g, cfg = sc.load(run_key)
    out = cfg.out_dir('34_clustering_matching')
    L = []
    pad_row = {int(c): k for k, c in enumerate(g['pad_id'])}
    hits = hits[hits['channel_id'].isin(pad_row)].copy()
    hits['row'] = hits['channel_id'].map(pad_row).astype(int)
    hits = hits.sort_values(['eventId', 'row']).reset_index(drop=True)
    ev = rays['eventId'].to_numpy()
    tx, ty, theta = cm.ray_dirs(rays)

    # --- clustering (pad frame; alignment-independent) ---
    adj = tile_adjacency(g['tile_verts'])
    L.append(f'tile adjacency: {len(adj)} touching pad pairs among {len(g["pad_id"])} live pads')
    labs, clus = {}, {}
    for mode in ('CC20', 'ADJ'):
        labs[mode] = cluster_labels(hits, mode, adj)
        clus[mode] = cluster_table(hits, labs[mode])
        L.append(f'{mode}: {len(clus[mode]):,} clusters, {clus[mode].n_pad.gt(1).mean()*100:.1f}% multi-pad')

    # --- robust alignment on CC20 clusters, seed 0 (as 06) ---
    xy, T, z = robust_alignment(rays, clus['CC20'], g, float(g['plane_z']))
    L.append(f'robust alignment: z = {z:.2f} mm (config {float(g["plane_z"]):g}), scale {T.s:.4f}, '
             f'rot {T.rotation_deg:.3f} deg (stage-06 fit: scale {float(g["s"]):.4f})')
    tiles_m3 = cm.apply_T(T, g['tile_verts'])
    sd = signed_distance(xy, tiles_m3, cm.apply_T(T, g['dead_verts']))
    fid = sd >= FID_MM
    L.append(f'fiducial: {fid.sum():,} rays >= {FID_MM:g} mm inside the live area')

    # --- pointing resolution (a): single isolated pads, residual in the pad frame ---
    per0 = match_distances(ev, xy, hits, labs['ADJ'], clus['ADJ'], T, tiles_m3)
    cA = clus['ADJ'].set_index('lab')
    best = per0.loc[per0.groupby('i')['dcen'].idxmin()]
    best = best[(best.n_pad == 1) & (best.dcen < 20) & (best.qmax >= 300)]
    rowof = pd.Series(np.arange(len(hits))).groupby(labs['ADJ']).first()  # lab -> first hit index
    hrow = hits['row'].to_numpy()[rowof.loc[best.lab].to_numpy()]
    q = tiles_m3[hrow]                                   # (n,4,2) tile of the hit pad
    cen = q.mean(1)
    e1 = q[:, 1] - q[:, 0]; e2 = q[:, 3] - q[:, 0]
    l1 = np.hypot(*e1.T); l2 = np.hypot(*e2.T)
    u1 = e1 / l1[:, None]; u2 = e2 / l2[:, None]
    dxy = xy[best.i.to_numpy()] - cen
    r1 = np.einsum('ij,ij->i', dxy, u1); r2 = np.einsum('ij,ij->i', dxy, u2)
    th_b = theta[best.i.to_numpy()]
    res_rows = []
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.4))
    for ax, (r, l, nm) in zip(axs[:2], [(r1, l1, 'pad axis 1'), (r2, l2, 'pad axis 2')]):
        p, e, (x, h) = fit_box(r, float(np.median(l)))
        ax.step(x, h, where='mid', color=C1, lw=2, label='data')
        xx = np.linspace(-20, 20, 400)
        ax.plot(xx, box_gauss(xx, *p), color=C2, lw=2,
                label=f'box {p[2]:.2f} mm (x) Gauss {p[3]:.2f}+-{e[3]:.2f} mm')
        ax.set(xlabel=f'track - pad centre along {nm} [mm]', ylabel='rays / 0.5 mm',
               title=f'single-pad hits ({len(r):,}); pad size {np.median(l):.2f} mm')
        ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
        L.append(f'pointing (single pads, {nm}): sigma = {p[3]:.3f} +- {e[3]:.3f} mm, '
                 f'fitted width {p[2]:.2f} (tile {np.median(l):.2f}) mm')
        for lo, hi in [(0, 10), (10, 20), (20, 30), (30, 45)]:
            k = (th_b >= lo) & (th_b < hi)
            if k.sum() > 300:
                pp, ee, _ = fit_box(r[k], float(np.median(l)))
                res_rows.append(dict(axis=nm, th_lo=lo, th_hi=hi, n=int(k.sum()), sigma=pp[3], err=ee[3]))
    rs = pd.DataFrame(res_rows); rs.to_csv(f'{out}/pointing_sigma_vs_theta.csv', index=False)
    for nm, col in [('pad axis 1', C2), ('pad axis 2', C3)]:
        s_ = rs[rs.axis == nm]
        axs[2].errorbar(0.5 * (s_.th_lo + s_.th_hi), s_.sigma, s_.err, fmt='-o', color=col, lw=2,
                        label=nm)
    axs[2].set(xlabel='track zenith angle [deg]', ylabel='pointing sigma [mm]',
               title='pointing resolution at the P2 plane vs angle', ylim=(0, None))
    axs[2].grid(color=GRID); axs[2].legend(frameon=False)
    fig.tight_layout(); fig.savefig(f'{out}/pointing_resolution.png', dpi=140); plt.close(fig)
    sig_pt = float(np.mean([r_ for r_ in rs.groupby('axis').apply(
        lambda d: np.average(d.sigma, weights=1 / d.err ** 2))]))

    # --- pointing resolution (b): boundary step (CC20 / CEN 20 hit definition) ---
    hit06 = np.zeros(len(ev), bool)
    hit06[per0.loc[per0.dcen <= 20, 'i'].to_numpy()] = True
    bb = np.arange(-5, 8, 0.5)
    c_ = np.digitize(sd, bb)
    xs = 0.5 * (bb[1:] + bb[:-1])
    ms = np.array([1 - hit06[c_ == k].mean() if (c_ == k).sum() > 100 else np.nan for k in range(1, len(bb))])
    ok = np.isfinite(ms)
    ps, pe = curve_fit(step, xs[ok], ms[ok], p0=[0.04, 0, 1])
    L.append(f'pointing (boundary step): sigma = {ps[2]:.3f} +- {np.sqrt(pe[2, 2]):.3f} mm, '
             f'edge at {ps[1]:+.2f} mm, plateau miss {100*ps[0]:.2f} %')
    L.append(f'=> pointing sigma used for TILE matching: {sig_pt:.2f} mm')

    # --- efficiency matrix: clustering x matching x seed, with accidentals ---
    rows = []
    n = len(ev)
    for mode in ('CC20', 'ADJ'):
        per = per0 if mode == 'ADJ' else match_distances(ev, xy, hits, labs[mode], clus[mode], T, tiles_m3)
        mix = match_distances(ev, xy, hits, labs[mode], clus[mode], T, tiles_m3, ev_shift=1000)
        for kind, win in [('CEN', 20.0), ('CEN', 12.0), ('TILE', 3 * sig_pt), ('TILE', 5 * sig_pt)]:
            for s in SEEDS:
                e_ = eff_from(per, n, fid, s, kind, win)
                a_ = eff_from(mix, n, fid, s, kind, win)
                E, A = e_.mean(), a_.mean()
                rows.append(dict(cluster=mode, match=kind, window=round(win, 2), seed=s, n=int(fid.sum()),
                                 eff=100 * E, eff_err=cm.binom_err(e_.sum(), len(e_)), acc=100 * A,
                                 eff_corr=100 * (E - A) / (1 - A)))
    tab = pd.DataFrame(rows); tab.to_csv(f'{out}/efficiency_matrix.csv', index=False)
    L.append('\nefficiency matrix (fiducial rays; acc = accidental rate from event mixing):')
    L.append(tab.to_string(index=False, float_format=lambda v: f'{v:.2f}'))

    # --- plot: eff and accidentals vs seed threshold ---
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
    styles = {('CC20', 'CEN', 20.0): (C1, 'CC20 + centroid<20 mm (today)'),
              ('ADJ', 'CEN', 20.0): (C4, 'ADJ + centroid<20 mm'),
              ('ADJ', 'TILE', round(3 * sig_pt, 2)): (C2, f'ADJ + tile<3 sigma ({3*sig_pt:.1f} mm)'),
              ('ADJ', 'TILE', round(5 * sig_pt, 2)): (C3, f'ADJ + tile<5 sigma ({5*sig_pt:.1f} mm)')}
    for (cl, mt, w), (col, lab) in styles.items():
        s_ = tab[(tab.cluster == cl) & (tab.match == mt) & (tab.window == w)]
        axs[0].plot(s_.seed, s_.eff, '-o', color=col, lw=2, ms=4, label=lab + ' raw')
        axs[0].plot(s_.seed, s_.eff_corr, '--', color=col, lw=1.5)
        axs[1].plot(s_.seed, s_.acc, '-o', color=col, lw=2, ms=4, label=lab)
    axs[0].set(xlabel='seed threshold: max pad amplitude of the cluster [ADC]', ylabel='efficiency [%]',
               title='efficiency (solid raw, dashed accidental-corrected)')
    axs[1].set(xlabel='seed threshold [ADC]', ylabel='accidental match rate [%] (event mixing)',
               title='random-coincidence contribution', yscale='log')
    for ax in axs:
        ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f'{out}/eff_vs_seed.png', dpi=140); plt.close(fig)

    # boundary-step plot
    fig, ax = plt.subplots(figsize=(6.5, 4.3))
    ax.plot(xs, 100 * ms, 'o', color=C1, label='data (CC20, centroid < 20 mm)')
    xx = np.linspace(-5, 8, 300)
    ax.plot(xx, 100 * step(xx, *ps), color=C2, lw=2, label=f'step (x) Gauss, sigma {ps[2]:.2f} mm')
    ax.set(xlabel='signed distance to live-area boundary [mm] (+ inside)', ylabel='miss rate [%]',
           title='Boundary step (robust alignment)', yscale='log')
    ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f'{out}/boundary_step.png', dpi=140); plt.close(fig)

    pd.DataFrame({'eventId': ev, 'x': xy[:, 0], 'y': xy[:, 1], 'sd': sd, 'fid': fid,
                  'theta': theta, 'tx': tx, 'ty': ty}).to_pickle(f'{out}/rays_aligned.pkl')
    np.savez(f'{out}/alignment.npz', R=T.R, s=T.s, pad_mean=T.pad_mean, m3_mean=T.m3_mean, z=z,
             sigma_point=sig_pt)
    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
