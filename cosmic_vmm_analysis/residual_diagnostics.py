#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Why is the det2 residual wider than the pads allow? Tests on the best triplets of the straight groups.

Residual (rx, ry) = det2 pad - midpoint of the det3 and det1 pads (best triplet per group, as in event_display.py),
groups with |residual - median| < --tol (wide on purpose, so the cut does not set the width). Robust fits (3 passes of
2.5-sigma clipping) of
  1. z spacing:  rx = a + k * (x1 - x3),  ry = b + k * (y1 - y3)   -> det2 is k * 620 mm away from half way
  2. rotation:   rx = a - phi * y2,       ry = b + phi * x2        -> det2 turned by phi relative to det3/det1
                 (the same with the det3 and det1 pads: a rotation of an outer chamber shows up there)
  3. per connector half: median residual of the tracks through each det2 / det3 / det1 VMM (a mis-mapped half
     sits off zero)
and the expected width from the pads alone: positions drawn uniformly inside each pad of the same triplets
(pad_w x pad_h, rotated by pad_angle), sigma_r = core width of rx, ry. Prints everything; OUT.png shows the fits.

Usage: residual_diagnostics.py cand.npz OUT_STEM [--conn-shift -1] [--tol 60] [--max-hits 12]
       [--order 'P2_IN:3=bot:rev|top:lin' ...]  (channel order of a connector, else the SPS one)
"""
import argparse
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import orientation_check as O  # noqa: E402
from straight_track_rate import groups, parse_veto  # noqa: E402
from event_display import best_triplets, shift_of, best_per_group, style, INK, MUTED, COL  # noqa: E402


def robust_fit(x, y, n_pass=3, k=2.5):
    """straight line y = a + b x with clipping; returns a, b, err_b, sigma, mask."""
    m = np.ones(len(x), bool)
    for _ in range(n_pass):
        A = np.vstack([np.ones(m.sum()), x[m]]).T
        coef, *_ = np.linalg.lstsq(A, y[m], rcond=None)
        r = y - (coef[0] + coef[1] * x)
        sig = 1.4826 * np.median(np.abs(r[m]))
        m = np.abs(r) < k * sig
    A = np.vstack([np.ones(m.sum()), x[m]]).T
    coef, *_ = np.linalg.lstsq(A, y[m], rcond=None)
    r = y[m] - A @ coef
    cov = np.linalg.inv(A.T @ A) * np.var(r)
    return coef[0], coef[1], float(np.sqrt(cov[1, 1])), float(1.4826 * np.median(np.abs(r))), m


def core(v):
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def pad_geometry(shift, choice=None):
    """(32, 64) arrays of pad width, height, angle (deg) for the same map as Mapper.xy_table."""
    mp = O.Mapper(shift)
    ch = {**mp.choice, **(choice or {})}
    W, H, A = (np.full((32, 64), np.nan) for _ in range(3))
    for (st, n), vm in mp.conn.items():
        sb, stp = O.OPTIONS[ch[(st, n)]]
        for v, strips in ((vm['bot'], sb), (vm['top'], stp)):
            r = mp.pad.reindex((n - 1) * 128 + strips)
            W[v] = r['pad_w'].to_numpy() if 'pad_w' in r else 12.07
            H[v] = r['pad_h'].to_numpy() if 'pad_h' in r else 11.61
            A[v] = r['pad_angle'].to_numpy() if 'pad_angle' in r else 0.0
    return W, H, A


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('cand')
    ap.add_argument('out_stem')
    ap.add_argument('--conn-shift', type=int, default=-1)
    ap.add_argument('--tol', type=float, default=60)
    ap.add_argument('--max-hits', type=int, default=12)
    ap.add_argument('--title', default='run_121')
    ap.add_argument('--order', action='append', help="channel order of one connector, e.g. 'P2_IN:3=bot:rev|top:lin'")
    ap.add_argument('--veto', default='', help='vmm:ch,... whose groups are dropped, e.g. 2:27,6:44')
    a = ap.parse_args()
    choice = O.parse_orders(a.order)

    X, Y = O.Mapper(a.conn_shift).xy_table(choice)
    h, _ = groups([a.cand], a.max_hits, X, Y, parse_veto(a.veto))
    t, _ = best_triplets(h)
    s = shift_of(t)
    b = best_per_group(t, s)
    b = b[b['r'] < a.tol].copy()
    for d in (3, 2, 1):
        b[f'v{d}'] = h.loc[b[f'i{d}'], 'vmm'].to_numpy()
        b[f'c{d}'] = h.loc[b[f'i{d}'], 'ch'].to_numpy()
    rx, ry = b['rx'].to_numpy(), b['ry'].to_numpy()
    print(f'{len(b)} groups with det2 within {a.tol:.0f} mm of the line (offset {s[0]:.1f}, {s[1]:.1f} mm)')
    print(f'measured core width: rx {core(rx):.1f} mm, ry {core(ry):.1f} mm')

    # expected from the pads alone
    W, H, A = pad_geometry(a.conn_shift, choice)
    rng = np.random.default_rng(3)
    sim = []
    for _ in range(50):
        pos = {}
        for d in (3, 2, 1):
            v, c = b[f'v{d}'].to_numpy(), b[f'c{d}'].to_numpy()
            u = (rng.random(len(b)) - 0.5) * W[v, c]
            w = (rng.random(len(b)) - 0.5) * H[v, c]
            ang = np.radians(np.nan_to_num(A[v, c]))
            pos[d] = (u * np.cos(ang) - w * np.sin(ang), u * np.sin(ang) + w * np.cos(ang))
        sim.append((pos[2][0] - 0.5 * (pos[3][0] + pos[1][0]), pos[2][1] - 0.5 * (pos[3][1] + pos[1][1])))
    sx = np.concatenate([q[0] for q in sim])
    sy = np.concatenate([q[1] for q in sim])
    print(f'expected from pad size alone: rx {np.std(sx):.1f} mm, ry {np.std(sy):.1f} mm (rms)')

    fig, axs = plt.subplots(2, 3, figsize=(14, 8.5), constrained_layout=True)
    res = {}
    # 1. z spacing
    for j, (lab, r, slope) in enumerate((('x', rx, (b['x1'] - b['x3']).to_numpy()),
                                          ('y', ry, (b['y1'] - b['y3']).to_numpy()))):
        a0, k, ek, sg, m = robust_fit(slope, r)
        res[f'z_{lab}'] = (k, ek, sg)
        print(f'1. z spacing from r{lab} vs ({lab}1 - {lab}3): k = {k:.3f} +- {ek:.3f} -> det2 is {k * 620:.0f} +- '
              f'{ek * 620:.0f} mm from half way; width after {sg:.1f} mm')
        ax = axs[0, j]
        style(ax)
        ax.scatter(slope[~m], r[~m], s=8, color=MUTED, alpha=0.4, lw=0, label='clipped')
        ax.scatter(slope[m], r[m], s=10, color=COL[2], lw=0, label='used')
        xx = np.linspace(slope.min(), slope.max(), 2)
        ax.plot(xx, a0 + k * xx, color=INK, lw=1.5, label=f'slope {k:.3f} +- {ek:.3f}')
        ax.set_xlabel(f'{lab}(det1) - {lab}(det3)  [mm]  (track inclination x 620 mm)', fontsize=8)
        ax.set_ylabel(f'r{lab} [mm]', fontsize=8)
        ax.set_title(f'1. z spacing: r{lab} vs track inclination', fontsize=9, color=INK)
        ax.legend(fontsize=7, frameon=False)
    # 2. rotation of det2 (and of the outer chambers)
    for d in (2, 3, 1):
        xd, yd = b[f'x{d}'].to_numpy(), b[f'y{d}'].to_numpy()
        _, kx, ekx, _, _ = robust_fit(yd - yd.mean(), rx)
        _, ky, eky, _, _ = robust_fit(xd - xd.mean(), ry)
        phi = 0.5 * (-kx + ky)
        ephi = 0.5 * np.hypot(ekx, eky)
        res[f'rot_det{d}'] = (phi, ephi, kx, ky)
        print(f'2. rotation vs det{d} position: rx vs y slope {kx:+.4f} +- {ekx:.4f}, ry vs x slope {ky:+.4f} +- '
              f'{eky:.4f} -> phi {np.degrees(phi):+.2f} +- {np.degrees(ephi):.2f} deg '
              f'(a rotation gives opposite signs; equal signs = scale / z)')
    ax = axs[0, 2]
    style(ax)
    q = ax.quiver(b['x2'], b['y2'], rx, ry, color=COL[2], angles='xy', scale_units='xy', scale=0.5, width=0.004)
    ax.set_aspect('equal')
    ax.set_title('2. residual arrows at the det2 pad (x2 length)', fontsize=9, color=INK)
    ax.set_xlabel('x [mm]', fontsize=8)
    ax.set_ylabel('y [mm]', fontsize=8)
    del q
    # 3. per VMM
    print('3. median residual per VMM of the track (n, rx, ry):')
    rows = []
    for d in (3, 2, 1):
        for v, g in b.groupby(f'v{d}'):
            rows.append((d, int(v), len(g), float(np.median(g['rx'])), float(np.median(g['ry'])),
                         core(g['rx'].to_numpy()) if len(g) > 4 else np.nan))
    per = pd.DataFrame(rows, columns=['det', 'vmm', 'n', 'rx_med', 'ry_med', 'rx_core'])
    print(per.round(1).to_string(index=False))
    ax = axs[1, 0]
    style(ax)
    for d in (3, 2, 1):
        p = per[per['det'] == d]
        ax.errorbar(p['vmm'], p['rx_med'], yerr=12 / np.sqrt(p['n']), fmt='o', color=COL[d], ms=5, label=f'det{d} rx')
        ax.errorbar(p['vmm'] + 0.25, p['ry_med'], yerr=12 / np.sqrt(p['n']), fmt='s', mfc='white', color=COL[d],
                    ms=5, label=f'det{d} ry')
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xlabel('VMM the track hit', fontsize=8)
    ax.set_ylabel('median residual [mm]', fontsize=8)
    ax.set_title('3. per connector half (VMM)', fontsize=9, color=INK)
    ax.legend(fontsize=6, frameon=False, ncol=3)
    # widths
    for j, (lab, r, sv) in enumerate((('x', rx, sx), ('y', ry, sy))):
        ax = axs[1, 1 + j]
        style(ax)
        bins = np.arange(-a.tol, a.tol + 1, 3)
        ax.hist(r, bins=bins, histtype='step', color=INK, lw=2, label=f'data, core {core(r):.1f} mm')
        ax.hist(sv, bins=bins, histtype='step', color=MUTED, lw=1.5, ls='--',
                weights=np.full(len(sv), len(r) / len(sv)), label=f'pads alone, rms {np.std(sv):.1f} mm')
        ax.set_xlabel(f'r{lab} [mm]', fontsize=8)
        ax.set_title(f'width of r{lab}: data vs pad size', fontsize=9, color=INK)
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f'{a.title}: what widens the det2 residual (connector shift {a.conn_shift:+d})', fontsize=11)
    fig.savefig(a.out_stem + '_residual_diag.png', dpi=120)
    per.to_csv(a.out_stem + '_residual_per_vmm.csv', index=False)


if __name__ == '__main__':
    sys.exit(main())
