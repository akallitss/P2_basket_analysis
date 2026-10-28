#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Toy Monte Carlo of the cosmic-muon rates the three-chamber lab stand should see: singles, doubles and triples.

Geometry: the three P2 chambers lying flat, stacked in the SPS order and spacing (det3 P2_IN z = 320, det2 P2_MID
z = 630, det1 P2_OUT z = 940 mm), same footprint (no lateral offset unless --offset). The sensitive area of a chamber is
the union of its read-out pads: the P2 BASKET pad map (1280 pads, Gerber centroids pad_cx/pad_cy), a point belongs to
the nearest pad centre when it lies inside that pad's (rotated) rectangle. Read-out connectors c3-c5 (from run_112;
--connectors 4,5,6 for runs 1-111); masked or dead (vmm:ch) can be removed with --mask-from (an efficiency_tp.py json,
key dead_or_masked) so the acceptance matches the live pads of a given run.

Flux: I(theta) = I0 cos^n(theta) per m2 s sr; through a horizontal plane J = 2 pi I0 / (n + 2). Default J = 1 /cm2/min
(167 /m2/s, the usual sea-level figure, n = 2); --J to set your own (e.g. scaled from a DREAM/M3 bench rate). Muons
are thrown on the plane of the chamber the coincidence requires (uniform over the bounding box of its sensitive area,
cos(theta) = U^(1/(n+2)), phi uniform) and propagated as straight lines: rate = J x A_box x fraction accepted. That
is exact for any coincidence that contains that chamber (no muon can reach it from outside the box).

The geometric rates are for 100 % efficient chambers; the observed rate is rate x eps3 eps2 eps1. The script also
emulates the straight-track cut of straight_track_rate.py on the MC triples (pad centres, det2 within --tol mm of the
midpoint of the det3 and det1 pads) so the cut acceptance is known too.

Usage: cosmic_acceptance_mc.py OUT_DIR [--n 2e6] [--J 1.0] [--n-exp 2] [--connectors 3,4,5]
       [--mask-from eff124.json] [--offset P2_MID:dx,dy] [--observed 813:49] [--eff2 0.52:0.02]
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import vmm_qa_path                                          # noqa: F401
import vmm_stations as VS

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STATIONS = {'det3 P2_IN': ('P2_IN', 320.0), 'det2 P2_MID': ('P2_MID', 630.0), 'det1 P2_OUT': ('P2_OUT', 940.0)}
NAMES = list(STATIONS)
CM2_MIN = 1e4 / 60.0      # 1 /cm2/min in /m2/s


class Chamber:
    """Sensitive pads of one chamber, in its own frame (mm), with a lateral offset in the stand frame."""

    def __init__(self, pads, live, offset=(0.0, 0.0)):
        self.pads = pads
        self.live = np.asarray(live, bool)
        self.tree = cKDTree(pads[['pad_cx', 'pad_cy']].to_numpy())
        a = np.deg2rad(pads['pad_angle'].to_numpy())
        self.u = np.c_[np.cos(a), np.sin(a)]
        self.half = np.c_[pads['pad_w'].to_numpy(), pads['pad_h'].to_numpy()] / 2
        self.offset = np.asarray(offset, float)

    def pad_of(self, xy):
        """Pad index hit at stand-frame points xy (N, 2), -1 outside any live pad."""
        p = xy - self.offset
        _, i = self.tree.query(p, k=3)
        out = np.full(len(p), -1)
        for k in range(i.shape[1]):     # nearest centre first; the 2nd/3rd catch points near a rotated pad edge
            j = i[:, k]
            d = p - self.tree.data[j]
            lu = np.abs(d[:, 0] * self.u[j, 0] + d[:, 1] * self.u[j, 1])
            lv = np.abs(-d[:, 0] * self.u[j, 1] + d[:, 1] * self.u[j, 0])
            inside = (out < 0) & (np.minimum(np.maximum(lu / self.half[j, 0], lv / self.half[j, 1]),
                                             np.maximum(lu / self.half[j, 1], lv / self.half[j, 0])) <= 1.0)
            out[inside] = j[inside]
        hit = out >= 0
        out[hit & ~self.live[np.maximum(out, 0)]] = -1
        return out

    def box(self):
        xy = self.pads.loc[self.live, ['pad_cx', 'pad_cy']].to_numpy() + self.offset
        return xy.min(0) - 10, xy.max(0) + 10


def masked_pads(mask_json, pads_all, conn_shift, orders):
    """Pad rows (channel_id) of the dead/masked (vmm:ch) of an efficiency_tp.py json, per station."""
    import orientation_check as O
    m = O.Mapper(conn_shift)
    m.choice.update(O.parse_orders(orders))
    X, Y = m.xy_table()
    tree = cKDTree(pads_all[['pad_cx', 'pad_cy']].to_numpy())
    dead = json.load(open(mask_json))['dead_or_masked']
    out = {}
    for name, lst in dead.items():
        ids = set()
        for s in lst:
            v, c = (int(t) for t in s.split(':'))
            if np.isfinite(X[v, c]):
                ids.add(int(pads_all.index[tree.query([X[v, c], Y[v, c]])[1]]))
        out[name] = ids
    return out


def throw(rng, n, lo, hi, z0, n_exp):
    x = rng.uniform(lo[0], hi[0], n)
    y = rng.uniform(lo[1], hi[1], n)
    cth = rng.random(n) ** (1.0 / (n_exp + 2))
    tth = np.sqrt(1 - cth ** 2) / cth
    ph = rng.uniform(0, 2 * np.pi, n)
    return x, y, z0, tth * np.cos(ph), tth * np.sin(ph), np.degrees(np.arccos(cth))


def hits_at(ch, x, y, z0, tx, ty, z):
    return ch.pad_of(np.c_[x + (z - z0) * tx, y + (z - z0) * ty])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out_dir')
    ap.add_argument('--n', type=float, default=2e6, help='muons thrown per reference chamber')
    ap.add_argument('--J', type=float, default=1.0, help='flux through a horizontal plane, /cm2/min')
    ap.add_argument('--n-exp', type=float, default=2.0, help='cos^n zenith exponent')
    ap.add_argument('--connectors', default='3,4,5', help='read-out connectors (c3-c5 from run_112, 4,5,6 before)')
    ap.add_argument('--mask-from', help='efficiency_tp.py json: remove its dead_or_masked pads')
    ap.add_argument('--order', action='append', default=['P2_IN:3=bot:rev|top:lin', 'P2_IN:5=bot:rev|top:lin'],
                    help='channel order overrides used to place the masked channels (lab map from run_112)')
    ap.add_argument('--z', default='320,630,940', help='z of det3, det2, det1 in mm')
    ap.add_argument('--offset', action='append', default=[], help="lateral offset 'P2_MID:dx,dy' (mm)")
    ap.add_argument('--tol', type=float, default=30.0, help='straight-track cut, as straight_track_rate.py')
    ap.add_argument('--observed', default='813:49', help='observed straight tracks per hour, value:err (run_124)')
    ap.add_argument('--eff2', default='0.52:0.017', help='measured det2 efficiency value:err (eff124.json)')
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    n = int(a.n)
    J = a.J * CM2_MIN
    zs = dict(zip(NAMES, (float(v) for v in a.z.split(','))))
    conns = [int(c) for c in a.connectors.split(',')]
    offs = {}
    for s in a.offset:
        st, v = s.split(':')
        offs[st] = tuple(float(t) for t in v.split(','))

    pads_all = pd.read_csv(VS.MAP_CSV).set_index('channel_id').sort_index()
    read = pads_all[(pads_all.index // 128 + 1).isin(conns)]
    dead = masked_pads(a.mask_from, pads_all, conns[0] - 4, a.order) if a.mask_from else {}
    ch = {}
    for name in NAMES:
        st = STATIONS[name][0]
        live = ~read.index.isin(sorted(dead.get(name, ())))
        ch[name] = Chamber(read, live, offs.get(st, (0.0, 0.0)))
    full = Chamber(pads_all, np.ones(len(pads_all), bool))

    res = {'params': vars(a), 'J_per_m2_s': J, 'I0_per_m2_s_sr': J * (a.n_exp + 2) / (2 * np.pi)}

    # sensitive areas: MC over the bounding box, checked against the sum of the pad areas
    lo, hi = full.box()
    u = rng.uniform(lo, hi, (n, 2))
    box_area = np.prod(hi - lo) * 1e-6
    res['area_full_m2'] = {'mc': box_area * np.mean(full.pad_of(u) >= 0), 'pads': pads_all['pad_area'].sum() * 1e-6}
    res['area_readout_m2'] = {}
    for name in NAMES:
        lo, hi = ch[name].box()
        u = rng.uniform(lo, hi, (n, 2))
        res['area_readout_m2'][name] = {
            'mc': float(np.prod(hi - lo) * 1e-6 * np.mean(ch[name].pad_of(u) >= 0)),
            'pads': float(read.loc[ch[name].live, 'pad_area'].sum() * 1e-6),
            'live_pads': int(ch[name].live.sum())}

    # coincidences, each thrown on the plane of a chamber it contains
    combos = {'det3': ['det3 P2_IN'], 'det2': ['det2 P2_MID'], 'det1': ['det1 P2_OUT'],
              'det3&det2': ['det3 P2_IN', 'det2 P2_MID'], 'det2&det1': ['det2 P2_MID', 'det1 P2_OUT'],
              'det3&det1': ['det3 P2_IN', 'det1 P2_OUT'],
              'det3&det2&det1': ['det2 P2_MID', 'det3 P2_IN', 'det1 P2_OUT']}
    rates = {}
    keep = None
    for key, need in combos.items():
        ref = need[0]
        lo, hi = ch[ref].box()
        x, y, z0, tx, ty, th = throw(rng, n, lo, hi, zs[ref], a.n_exp)
        h = {nm: hits_at(ch[nm], x, y, z0, tx, ty, zs[nm]) for nm in need}
        ok = np.all([h[nm] >= 0 for nm in need], axis=0)
        A = np.prod(hi - lo) * 1e-6
        r = J * A * ok.mean()
        rates[key] = {'Hz': r, 'per_h': r * 3600, 'err_per_h': J * A * np.sqrt(ok.sum()) / n * 3600}
        if key == 'det3&det2&det1':
            keep = dict(x=x[ok], y=y[ok], tx=tx[ok], ty=ty[ok], th=th[ok], **{nm: h[nm][ok] for nm in need})
    res['geometric_rates'] = rates

    # straight-track cut on the pad centres (what straight_track_rate.py sees for a perfect chamber)
    P = read[['pad_cx', 'pad_cy']].to_numpy()
    p3, p2, p1 = (P[keep[nm]] + ch[nm].offset for nm in NAMES)
    mid = 0.5 * (p3 + p1)
    rr = np.hypot(*(p2 - mid).T)
    rr0 = np.hypot(*(p2 - mid - np.median(p2 - mid, 0)).T)
    f_cut = float(np.mean(rr0 < a.tol))
    res['straight_cut'] = {'tol_mm': a.tol, 'acceptance': f_cut,
                           'residual_core_mm': float(0.5 * (np.percentile(rr, 68.27))),
                           'theta_triples_deg': {'median': float(np.median(keep['th'])),
                                                 'p90': float(np.percentile(keep['th'], 90)),
                                                 'max': float(keep['th'].max())}}
    R3 = rates['det3&det2&det1']['per_h'] * f_cut

    # what the observed rate says about the efficiencies
    obs, obs_e = (float(v) for v in a.observed.split(':'))
    e2, e2_e = (float(v) for v in a.eff2.split(':'))
    prod = obs / R3
    prod_e = prod * np.hypot(obs_e / obs, rates['det3&det2&det1']['err_per_h'] / rates['det3&det2&det1']['per_h'])
    e31 = prod / e2
    e31_e = e31 * np.hypot(prod_e / prod, e2_e / e2)
    res['vs_observed'] = {'observed_per_h': obs, 'expected_straight_per_h_eff1': R3,
                          'eps_product': prod, 'eps_product_err': prod_e,
                          'eps3_eps1_given_eps2': e31, 'err': e31_e,
                          'eps_outer_if_equal': float(np.sqrt(e31)),
                          'eps_common_if_all_equal': float(prod ** (1 / 3)),
                          'table_per_h': {f'{e:.2f}': R3 * e ** 3 for e in (1.0, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5)}}
    json.dump(res, open(os.path.join(a.out_dir, 'acceptance_mc.json'), 'w'), indent=1, default=float)

    # ---- printout
    print(f"flux J = {a.J:g} /cm2/min = {J:.1f} /m2/s  (I0 = {res['I0_per_m2_s_sr']:.1f} /m2/s/sr, cos^{a.n_exp:g})")
    print(f"full P2 pad area: MC {res['area_full_m2']['mc']:.4f} m2, sum of pads {res['area_full_m2']['pads']:.4f} m2")
    for nm, v in res['area_readout_m2'].items():
        print(f"  {nm:12s} read-out c{a.connectors}: {v['live_pads']} live pads, MC {v['mc']:.4f} m2 "
              f"(pads {v['pads']:.4f})")
    print('geometric rates (100 % efficient chambers):')
    for k, v in rates.items():
        print(f"  {k:15s} {v['Hz']:8.3f} Hz  {v['per_h']:9.0f} +- {v['err_per_h']:.0f} /h  {v['per_h'] / 60:7.1f} /min")
    print(f"straight cut (det2 within {a.tol:g} mm of the det3-det1 midpoint, pad centres): acceptance {f_cut:.3f}; "
          f"triples theta median {res['straight_cut']['theta_triples_deg']['median']:.1f} deg, "
          f"p90 {res['straight_cut']['theta_triples_deg']['p90']:.1f} deg")
    print(f"expected straight triples at eff = 1: {R3:.0f} /h")
    for e, v in res['vs_observed']['table_per_h'].items():
        print(f"   eff {e} per chamber -> {v:7.0f} /h")
    print(f"observed {obs:.0f} +- {obs_e:.0f} /h -> eps3*eps2*eps1 = {prod:.3f} +- {prod_e:.3f} "
          f"(common eps {prod ** (1 / 3):.2f}); with eps2 = {e2:.2f}: eps3*eps1 = {e31:.2f} +- {e31_e:.2f} "
          f"(each ~{np.sqrt(e31):.2f})")

    # ---- figures
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    allp = pads_all[['pad_cx', 'pad_cy']].to_numpy()
    ax[0].scatter(*allp.T, s=2, c='0.85', label='P2 pads (all connectors)')
    for nm, c in zip(NAMES, ('C0', 'C1', 'C2')):
        q = P[ch[nm].live] + ch[nm].offset
        ax[0].scatter(*q.T, s=2, c=c, alpha=0.5, label=f'{nm} live read-out')
    sel = rng.choice(len(keep['x']), min(150, len(keep['x'])), replace=False)
    for i in sel:
        ax[0].plot([p3[i, 0], p1[i, 0]], [p3[i, 1], p1[i, 1]], '-', c='k', lw=0.3, alpha=0.5)
    ax[0].set_aspect('equal')
    ax[0].set_xlabel('x [mm]')
    ax[0].set_ylabel('y [mm]')
    ax[0].set_title(f'read-out region c{a.connectors} and 150 MC triples (det3 to det1 pad)')
    ax[0].legend(fontsize=7, markerscale=4, loc='upper left')
    ax[1].hist(keep['th'], bins=60, range=(0, 60), histtype='step', lw=1.5)
    ax[1].set_xlabel('zenith angle of triple-coincidence muons [deg]')
    ax[1].set_ylabel('MC tracks')
    ax[1].set_title(f'cos^{a.n_exp:g} flux, median {np.median(keep["th"]):.1f} deg')
    e = np.linspace(0.3, 1.0, 100)
    ax[2].plot(e, R3 * e ** 3, '-', label='all three chambers at eff')
    ax[2].plot(e, R3 * e2 * e ** 2, '--', label=f'det2 at {e2:.2f}, det3/det1 at eff')
    ax[2].axhspan(obs - obs_e, obs + obs_e, color='C3', alpha=0.25, label=f'observed run_124 {obs:.0f} +- {obs_e:.0f} /h')
    ax[2].set_xlabel('chamber efficiency')
    ax[2].set_ylabel('straight triple tracks per hour')
    ax[2].set_title(f'J = {a.J:g} /cm2/min; at eff 1: {R3:.0f} /h')
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out_dir, 'acceptance_mc.png'), dpi=130)
    print('wrote', a.out_dir)


if __name__ == '__main__':
    main()
