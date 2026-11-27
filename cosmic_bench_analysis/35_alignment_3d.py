#!/usr/bin/env python3
"""35_alignment_3d.py -- is the P2 plane TILTED (out-of-plane rotation) w.r.t. M3?

Stage 03/06 align P2 to M3 in 2D only: in-plane rotation + translation + scale
at one fitted plane z. A tilt of the detector (rotation about an in-plane axis,
any direction in XY) makes the true crossing height position-dependent,
    z(x, y) = z0 + a*x + b*y (+ c*r^2 for a bow),
so the track is projected to the wrong height away from the tilt axis. The
in-plane error that produces is (tx, ty) * dz(x, y): it grows with the track
slope and with the distance from the axis, i.e. it is LARGEST at the detector
edges -- exactly where edge / dead-border misses would appear.

Signatures used (all on matched track-cluster pairs):
  1. local z map: in each 40 mm cell fit  r = d + t * dz_cell  (d = in-plane
     offset, t = slope); a tilt is a plane in dz_cell(x, y), a bow is a bowl.
  2. global fit of the full linear model
        r_x = dx - w*y + e*x + tx*(dz + a*x + b*y + c*rho^2)
        r_y = dy + w*x + e*y + ty*(dz + a*x + b*y + c*rho^2)
     then the exact refit: alternate (track to tilted plane) <-> (2D Procrustes).
  3. what it does to the efficiency, the edge and dead-border miss profiles and
     the guard band the dead pads need.

Usage: python3 35_alignment_3d.py [run_key]
Out:   <out>/35_alignment_3d/{local_z_map,residual_maps,miss_profiles}.png, summary.txt, fit.csv
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_study_cache as sc     # noqa: E402
import p2_study_common as cm    # noqa: E402
import p2_align as pa           # noqa: E402

INK, INK2, GRID = '#1f1f1f', '#5c5c5c', '#e6e6e6'
C_OLD, C_NEW, C_GEO = '#9c9c9c', '#257de0', '#ec4647'


def project_plane(rays, z0, a=0.0, b=0.0, c=0.0, xc=0.0, yc=0.0, n_iter=4):
    """Track intersection with z = z0 + a(x-xc) + b(y-yc) + c((x-xc)^2+(y-yc)^2)."""
    x, y = pa.project_to_z(rays, z0)
    for _ in range(n_iter):
        dx, dy = x - xc, y - yc
        x, y = pa.project_to_z(rays, z0 + a * dx + b * dy + c * (dx * dx + dy * dy))
    return np.column_stack([x, y])


def matched_pairs(ev, xy, clu, T, r=20.0):
    cxy = cm.apply_T(T, clu[['x_pad', 'y_pad']].to_numpy())
    d, j = cm.nearest(ev, xy, cxy, clu['eventId'].to_numpy())
    ok = d <= r
    return ok, j


def robust_lsq(A, y, n_iter=5, k=3.0):
    keep = np.ones(len(y), bool)
    for _ in range(n_iter):
        p, *_ = np.linalg.lstsq(A[keep], y[keep], rcond=None)
        res = y - A @ p
        s = 1.4826 * np.median(np.abs(res[keep]))
        keep = np.abs(res) < k * s
    cov = np.linalg.inv(A[keep].T @ A[keep]) * np.var(res[keep])
    return p, np.sqrt(np.diag(cov)), keep


def main():
    run_key = sys.argv[1] if len(sys.argv) > 1 else 'det2_long_repos'
    hits, rays, g, cfg = sc.load(run_key)
    out = cfg.out_dir('35_alignment_3d')
    cpath = cfg.out_dir('34_study_cache') + '/clusters_cc20.pkl'
    clu = pd.read_pickle(cpath) if os.path.isfile(cpath) else cm.cc_clusters(hits, 20.0)
    z0 = float(g['plane_z'])
    ev = rays['eventId'].to_numpy()
    tx, ty, theta = cm.ray_dirs(rays)
    T0 = cm.make_T(g)
    xy0 = project_plane(rays, z0)
    act0, _, _ = cm.active_mask(xy0, g, T0)
    ok, j = matched_pairs(ev, xy0, clu, T0)
    sel = ok & act0
    pxy = cm.apply_T(T0, clu[['x_pad', 'y_pad']].to_numpy())[j[sel]]
    X, Y = xy0[sel, 0], xy0[sel, 1]
    xc, yc = float(np.median(X)), float(np.median(Y))
    rx, ry = pxy[:, 0] - X, pxy[:, 1] - Y          # cluster - track (M3 frame)
    TX, TY = tx[sel], ty[sel]
    lines = [f'run {run_key}: {sel.sum():,} matched active pairs, plane z0 = {z0:g} mm, '
             f'centre ({xc:.0f}, {yc:.0f}) mm']

    # ---- 1. local z map: per 40 mm cell, r = d + t*dz  (both axes jointly) ----
    cell = 40.0
    ix = np.floor((X - xc) / cell).astype(int); iy = np.floor((Y - yc) / cell).astype(int)
    rows = []
    for (a_, b_), idx in pd.DataFrame({'a': ix, 'b': iy}).groupby(['a', 'b']).groups.items():
        idx = np.asarray(idx)
        if len(idx) < 150:
            continue
        A = np.zeros((2 * len(idx), 3)); yv = np.concatenate([rx[idx], ry[idx]])
        A[:len(idx), 0] = 1; A[len(idx):, 1] = 1
        A[:len(idx), 2] = TX[idx]; A[len(idx):, 2] = TY[idx]
        p, e, _ = robust_lsq(A, yv)
        rows.append(dict(xc=xc + (a_ + .5) * cell, yc=yc + (b_ + .5) * cell, n=len(idx),
                         dx=p[0], dy=p[1], dz=p[2], dz_err=e[2]))
    zm = pd.DataFrame(rows)
    zm.to_csv(f'{out}/local_z_map.csv', index=False)
    # plane fit to the local dz values (weighted)
    w = 1 / zm.dz_err ** 2
    Ap = np.column_stack([np.ones(len(zm)), zm.xc - xc, zm.yc - yc])
    pp = np.linalg.solve(Ap.T @ (Ap * w.to_numpy()[:, None]), Ap.T @ (w * zm.dz))
    chi_flat = float((w * (zm.dz - np.average(zm.dz, weights=w)) ** 2).sum())
    chi_plane = float((w * (zm.dz - Ap @ pp) ** 2).sum())
    lines.append(f'local-z cells: {len(zm)} (>=150 pairs, {cell:g} mm); dz spread '
                 f'{zm.dz.min():+.1f} .. {zm.dz.max():+.1f} mm, median err {zm.dz_err.median():.2f} mm')
    lines.append(f'  plane fit to cells: dz0 {pp[0]:+.2f} mm, slope_x {pp[1]*1e3:+.2f} mm/m, '
                 f'slope_y {pp[2]*1e3:+.2f} mm/m;  chi2 flat {chi_flat:.0f} -> plane {chi_plane:.0f} '
                 f'(ndf {len(zm)-3})')

    # ---- 2. global linear model (+bow) ---------------------------------------
    u, v = X - xc, Y - yc
    rho2 = (u * u + v * v) / 1e4                      # per (100 mm)^2
    n = len(u)
    names = ['dx', 'dy', 'w_rad', 'e_scale', 'dz', 'a_x', 'b_y', 'c_bow']
    A = np.zeros((2 * n, 8))
    A[:n, 0] = 1; A[n:, 1] = 1
    A[:n, 2] = -v; A[n:, 2] = u
    A[:n, 3] = u; A[n:, 3] = v
    for k, f in enumerate([np.ones(n), u, v, rho2]):
        A[:n, 4 + k] = TX * f; A[n:, 4 + k] = TY * f
    yv = np.concatenate([rx, ry])
    p, e, keep = robust_lsq(A, yv)
    fit = pd.DataFrame({'par': names, 'val': p, 'err': e})
    fit.to_csv(f'{out}/fit_linear.csv', index=False)
    tilt = np.degrees(np.arctan(np.hypot(p[5], p[6])))
    tilt_dir = np.degrees(np.arctan2(p[6], p[5]))
    lines.append('global linear model (cluster - track):')
    for nm, val, er in zip(names, p, e):
        lines.append(f'  {nm:8s} {val:+.5f} +- {er:.5f}')
    lines.append(f'  => tilt {tilt:.3f} deg (gradient direction {tilt_dir:.0f} deg in M3 XY), '
                 f'dz over +-200 mm = +-{200*np.hypot(p[5], p[6]):.2f} mm; '
                 f'bow c = {p[7]:+.2f} mm per (100 mm)^2')
    s0 = 1.4826 * np.median(np.abs(yv - A[:, :2] @ p[:2]))
    s1 = 1.4826 * np.median(np.abs(yv - A @ p))
    lines.append(f'  robust residual sigma/axis: offsets only {s0:.3f} -> full model {s1:.3f} mm')

    # ---- 2b. exact refit: tilted plane <-> 2D Procrustes, iterate -------------
    zz, a, b, c = z0, 0.0, 0.0, 0.0
    T = T0
    for it in range(6):
        xy = project_plane(rays, zz, a, b, c, xc, yc)
        okm, jm = matched_pairs(ev, xy, clu, T)
        act, _, _ = cm.active_mask(xy, g, T)
        s_ = okm & act
        # refit 2D transform on these pairs (pad -> M3 frame)
        cp = clu[['x_pad', 'y_pad']].to_numpy()[jm[s_]]
        T = pa.fit_transform(xy[s_, 0], xy[s_, 1], cp[:, 0], cp[:, 1])
        pxy_ = cm.apply_T(T, cp)
        r_x, r_y = pxy_[:, 0] - xy[s_, 0], pxy_[:, 1] - xy[s_, 1]
        uu, vv = xy[s_, 0] - xc, xy[s_, 1] - yc
        Tx_, Ty_ = tx[s_], ty[s_]
        nn = len(uu)
        B = np.zeros((2 * nn, 4))
        for k, f in enumerate([np.ones(nn), uu, vv, (uu * uu + vv * vv)]):
            B[:nn, k] = Tx_ * f; B[nn:, k] = Ty_ * f
        q, _, _ = robust_lsq(B, np.concatenate([r_x, r_y]))
        zz += q[0]; a += q[1]; b += q[2]; c += q[3]
        lines.append(f'  iter {it}: z0 {zz:.2f}  a {a*1e3:+.2f}  b {b*1e3:+.2f} mm/m  '
                     f'c {c*1e4:+.3f} mm/(100mm)^2  | rot {T.rotation_deg:.3f} s {T.s:.4f}  '
                     f'step |dz| {abs(q[0]):.3f}')
    # tilt-only variant (no bow) for comparison
    variants = {}
    for name, (A_, B_, C_) in {'2D (stage 06)': (0, 0, 0), 'tilt': (a, b, 0), 'tilt+bow': (a, b, c)}.items():
        if name == '2D (stage 06)':
            xy_, T_ = xy0, T0
        else:
            zt = zz
            T_ = T
            for _ in range(3):
                xy_ = project_plane(rays, zt, A_, B_, C_, xc, yc)
                okm, jm = matched_pairs(ev, xy_, clu, T_)
                act_, _, _ = cm.active_mask(xy_, g, T_)
                s_ = okm & act_
                cp = clu[['x_pad', 'y_pad']].to_numpy()[jm[s_]]
                T_ = pa.fit_transform(xy_[s_, 0], xy_[s_, 1], cp[:, 0], cp[:, 1])
                # re-centre z for this variant (slope-weighted mean residual)
                pr = cm.apply_T(T_, cp) - xy_[s_]
                tt = np.concatenate([tx[s_], ty[s_]]); rr = np.concatenate([pr[:, 0], pr[:, 1]])
                m_ = np.abs(rr) < 15
                zt += float((tt[m_] * rr[m_]).sum() / (tt[m_] ** 2).sum())
            xy_ = project_plane(rays, zt, A_, B_, C_, xc, yc)
        variants[name] = (xy_, T_)

    # ---- 3. efficiency + miss profiles per variant ---------------------------
    bins = np.arange(-10, 31, 2.0)
    prof = {}
    for name, (xy_, T_) in variants.items():
        act_, nn_, _ = cm.active_mask(xy_, g, T_)
        d_, _ = cm.nearest(ev, xy_, cm.apply_T(T_, clu[['x_pad', 'y_pad']].to_numpy()),
                           clu['eventId'].to_numpy())
        hit = d_ <= float(g['match_r'])
        # residual width on matched
        okm, jm = matched_pairs(ev, xy_, clu, T_)
        mm = okm & act_
        pr = cm.apply_T(T_, clu[['x_pad', 'y_pad']].to_numpy())[jm[mm]] - xy_[mm]
        sig = 1.4826 * np.median(np.abs(pr), axis=0)
        # distance to dead tiles and to the outer array boundary (no guard, no cell mask)
        act_ng, nn_ng, _ = cm.active_mask(xy_, g, T_, guard=0.0)
        dd = pa.dist_to_quads(xy_, cm.apply_T(T_, g['dead_verts']))
        live_q = cm.apply_T(T_, g['tile_verts'])
        dl = pa.dist_to_quads(xy_, live_q)
        inside_live = dl == 0
        # signed distance to the edge of the live area: + inside (approx by nn to pad
        # tiles from the outside); here: rays inside a live tile -> distance to the
        # nearest point not inside a live tile is approximated by (half-pitch - nn)
        k_ = (nn_ng <= 8.4)
        e_ = 100 * (1 - hit[act_].mean())
        prof[name] = dict(eff=100 * hit[act_].mean(), n=int(act_.sum()), sig=sig,
                          dd=dd, hit=hit, act_ng=act_ng, inside_live=inside_live, dl=dl, xy=xy_)
        lines.append(f'{name:14s}: active {act_.sum():,}  eff {100*hit[act_].mean():.2f} '
                     f'+- {cm.binom_err(hit[act_].sum(), act_.sum()):.2f} %  (miss {e_:.2f} %)  '
                     f'residual sigma x/y {sig[0]:.2f}/{sig[1]:.2f} mm')

    # miss rate vs distance to the nearest dead tile (no guard) and to the array edge
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    for name, col in zip(variants, (C_OLD, C_NEW, C_GEO)):
        P = prof[name]
        m = P['act_ng'] & np.isfinite(P['dd'])
        c_ = np.digitize(P['dd'][m], bins)
        mr = [100 * (1 - P['hit'][m][c_ == k].mean()) if (c_ == k).sum() > 50 else np.nan
              for k in range(1, len(bins))]
        axs[0].plot(bins[:-1] + 1, mr, '-o', color=col, lw=2, ms=4, label=name)
        # edge: rays in the outermost ring of pads, vs distance OUTSIDE the live tiles
        # (dl > 0 means outside every live tile) and inside -> use -distance to edge
        ok = np.isfinite(P['dl'])
        eb = np.arange(0, 12.1, 1.0)
        out_ = ok & (P['dl'] > 0) & (P['dd'] > 12)   # outside live tiles, not next to a dead pad
        c2 = np.digitize(P['dl'][out_], eb)
        mr2 = [100 * (1 - P['hit'][out_][c2 == k].mean()) if (c2 == k).sum() > 30 else np.nan
               for k in range(1, len(eb))]
        axs[1].plot(eb[:-1] + .5, mr2, '-o', color=col, lw=2, ms=4, label=name)
    axs[0].set(xlabel='distance from the track to the nearest dead/shorted/unplugged pad tile [mm]',
               ylabel='miss rate [%]', title='Dead-border profile (live Voronoi cells, no guard)')
    axs[0].axvline(float(g['dead_guard']), color=INK2, ls=':', lw=1)
    axs[0].text(float(g['dead_guard']) + .3, axs[0].get_ylim()[1] * .9, 'current 6 mm guard',
                color=INK2, fontsize=8)
    axs[1].set(xlabel='distance OUTSIDE the outermost live pad tile [mm] (far from dead pads)',
               ylabel='miss rate [%]', title='Outer-edge profile')
    for ax in axs:
        ax.grid(color=GRID); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(f'{out}/miss_profiles.png', dpi=140); plt.close(fig)

    # local z map + residual maps (before / after)
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    sc_ = axs[0].scatter(zm.xc, zm.yc, c=zm.dz, cmap='RdBu_r', s=260, marker='s',
                         vmin=-np.abs(zm.dz).max(), vmax=np.abs(zm.dz).max())
    plt.colorbar(sc_, ax=axs[0], label='local dz [mm] (fit r = d + slope*dz per cell)')
    axs[0].set(title=f'Local crossing-height offset, 40 mm cells (z0 = {z0:g})', xlabel='x M3 [mm]',
               ylabel='y M3 [mm]', aspect='equal')
    for ax, name in zip(axs[1:], ['2D (stage 06)', 'tilt+bow']):
        xy_, T_ = variants[name]
        act_, _, _ = cm.active_mask(xy_, g, T_)
        okm, jm = matched_pairs(ev, xy_, clu, T_)
        mm = okm & act_
        pr = cm.apply_T(T_, clu[['x_pad', 'y_pad']].to_numpy())[jm[mm]] - xy_[mm]
        df = pd.DataFrame({'a': np.floor(xy_[mm, 0] / 40), 'b': np.floor(xy_[mm, 1] / 40),
                           'rx': pr[:, 0], 'ry': pr[:, 1]})
        gb = df.groupby(['a', 'b']).agg(rx=('rx', 'median'), ry=('ry', 'median'), n=('rx', 'size'))
        gb = gb[gb.n > 150].reset_index()
        mag = np.hypot(gb.rx, gb.ry)
        q = ax.quiver((gb.a + .5) * 40, (gb.b + .5) * 40, gb.rx, gb.ry, mag, cmap='viridis',
                      clim=(0, 2.0), angles='xy', scale_units='xy', scale=0.06)
        plt.colorbar(q, ax=ax, label='|median residual| [mm]')
        ax.set(title=f'{name}: median cluster-track residual (x16.7)', xlabel='x M3 [mm]',
               aspect='equal')
    fig.tight_layout(); fig.savefig(f'{out}/local_z_and_residual_maps.png', dpi=140); plt.close(fig)

    # residual vs slope, in edge thirds (the direct tilt signature)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    for ax, (coord, slope, rr, lab) in zip(axs, [(u, TX, rx, 'x'), (v, TY, ry, 'y')]):
        for lo, hi, col in [(-1e9, -80, C_OLD), (-80, 80, INK2), (80, 1e9, C_NEW)]:
            m = (coord > lo) & (coord < hi) & (np.abs(rr) < 15)
            sb = np.linspace(-0.6, 0.6, 13)
            c_ = np.digitize(slope[m], sb)
            med = [np.median(rr[m][c_ == k]) if (c_ == k).sum() > 80 else np.nan for k in range(1, len(sb))]
            ax.plot(.5 * (sb[1:] + sb[:-1]), med, '-o', color=col, lw=2, ms=4,
                    label=f'{lab}-{lab}c in ({max(lo,-999):.0f}, {min(hi,999):.0f}) mm')
        ax.set(xlabel=f'track slope d{lab}/dz', ylabel=f'median residual {lab} [mm] (2D alignment)',
               title=f'{lab}: residual vs slope by position (slope of line = local dz)')
        ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f'{out}/residual_vs_slope_by_position.png', dpi=140); plt.close(fig)

    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
