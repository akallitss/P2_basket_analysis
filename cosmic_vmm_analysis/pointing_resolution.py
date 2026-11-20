#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pointing resolution of the tag-and-probe prediction in each probe chamber, and the probe search radius / fiducial
edge margin that scale with it (20 Nov 2026).

Why: the prediction in the probe chamber is the straight line through the two tag cluster positions. For P2_MID it is an
INTERPOLATION (tags IN + OUT, f = 0.5), for P2_IN and P2_OUT an EXTRAPOLATION (f = -1 and 2: tags MID + OUT, IN + MID).
With z = 320 / 630 / 940 mm and the prediction p = (1 - f) a + f b, per axis

    sigma_point^2 = (1 - f)^2 sigma_A^2 + f^2 sigma_B^2                (pointing resolution: the line only)
    sigma_res^2   = sigma_point^2 + sigma_P^2                          (probe hit minus prediction)

    P2_MID: sigma_point = 0.71 s,  sigma_res = 1.22 s      (equal chambers, s = one chamber's resolution)
    P2_IN : sigma_point = 2.24 s,  sigma_res = 2.45 s
    P2_OUT: sigma_point = 2.24 s,  sigma_res = 2.45 s

so a fixed search radius (15 mm) is ~2.5 sigma for MID but only ~1.2 sigma for IN/OUT with s = 5 mm (the tracker
resolution measured at the SPS): the IN/OUT "far" class is mostly the pointing tail. Because the stand is equally
spaced, sigma_res(IN/OUT) = 2 sigma_res(MID) for ANY s_IN, s_MID, s_OUT (the residuals alone cannot separate the three
chamber resolutions); multiple scattering of cosmic muons (~1 mrad x 310 mm) is negligible next to 5 mm.

What it does:
  * prints the model (factors per probe, sigma for --sigma-tag);
  * fits the measured radial distance of the nearest probe hit (muon_funnel.py d_found + d_far, 5 mm bins, the
    140-310 ns window) per probe and group with a 2-D Gaussian (Rayleigh in r) + chance hits (a r + b), binned Poisson
    likelihood for r < --fit-max mm -> sigma_res measured, the ratio IN/MID, OUT/MID (model: 2) and the implied
    one-chamber resolution;
  * writes OUT.json with, per probe, the radius used to accept a probe hit: probe_r = max(--min-r, k sigma_res), with
    k = sqrt(-2 ln(1 - --contain)) (2-D containment; 0.99 -> 3.03), and the fiducial edge margin = --margin-nsig
    sigma_point (the predicted point must be that far inside the chamber, so that the true crossing point is inside).
    --source model uses --sigma-tag; --source data the fitted sigma_res of --group (sigma_point from it with the
    equal-chamber ratio).

Pass OUT.json to muon_funnel.py --pointing.

Usage: pointing_resolution.py --summary SUMMARY.json [--group REF] [--source model|data] [--sigma-tag 5]
                              [--contain 0.99] [--margin-nsig 2] --out pointing.json
"""
import argparse
import json
import sys

import numpy as np

Z = {'P2_IN': 320.0, 'P2_MID': 630.0, 'P2_OUT': 940.0}
TAGS = {'P2_MID': ('P2_IN', 'P2_OUT'), 'P2_OUT': ('P2_IN', 'P2_MID'), 'P2_IN': ('P2_MID', 'P2_OUT')}
PROBES = ['P2_IN', 'P2_MID', 'P2_OUT']


def factors(probe):
    """(w_A, w_B): prediction = w_A a + w_B b."""
    a, b = TAGS[probe]
    f = (Z[probe] - Z[a]) / (Z[b] - Z[a])
    return 1 - f, f


def sigma_model(probe, s):
    """s: dict chamber -> one-chamber resolution per axis (mm). -> (sigma_point, sigma_res)."""
    a, b = TAGS[probe]
    wa, wb = factors(probe)
    sp = np.hypot(wa * s[a], wb * s[b])
    return float(sp), float(np.hypot(sp, s[probe]))


def k_contain(c):
    """radius / sigma that contains the fraction c of a 2-D Gaussian (Rayleigh)."""
    return float(np.sqrt(-2 * np.log(1 - c)))


def hist(P):
    """d_found + d_far of one probe -> (lower bin edges, counts) for finite distances (5 mm bins, 300 = overflow)."""
    h = {}
    for k in ('d_found', 'd_far'):
        for b, n in P.get(k, {}).items():
            if b != 'inf':
                h[int(b)] = h.get(int(b), 0) + n
    e = np.arange(0, 300, 5)
    return e, np.array([h.get(int(x), 0) for x in e], float)


def fit(e, n, rmax, width=5.0):
    """binned Poisson fit of N_s Rayleigh(sigma) + chance (a r + b) on bins [e, e + width) with e + width <= rmax.
    -> dict sigma, err, n_sig, n_bkg_in_3sigma."""
    from scipy.optimize import minimize
    m = e + width <= rmax
    lo, hi, k = e[m], e[m] + width, n[m]

    def mu(q):
        ns, sig, a, b = np.exp(q[0]), np.exp(q[1]), np.exp(q[2]), np.exp(q[3])
        sgn = ns * (np.exp(-lo ** 2 / (2 * sig ** 2)) - np.exp(-hi ** 2 / (2 * sig ** 2)))
        bkg = a * (hi ** 2 - lo ** 2) / 2 + b * (hi - lo)
        return sgn + bkg

    def nll(q):
        u = np.maximum(mu(q), 1e-12)
        return float(np.sum(u - k * np.log(u)))

    tot = max(k.sum(), 1.0)
    tail = max(k[lo >= 40].mean() if np.any(lo >= 40) else 1.0, 0.5)
    best = None
    for s0 in (4.0, 7.0, 10.0, 14.0, 20.0):
        q0 = np.log([tot, s0, tail / (width * 60.0), tail / (2 * width)])
        r = minimize(nll, q0, method='Nelder-Mead', options={'maxiter': 20000, 'xatol': 1e-6, 'fatol': 1e-6})
        if best is None or r.fun < best.fun:
            best = r
    q = best.x
    # error on sigma from the curvature of the nll in log sigma (profile over the others by refitting)
    def prof(ls):
        r = minimize(lambda z: nll([z[0], ls, z[1], z[2]]), [q[0], q[2], q[3]], method='Nelder-Mead',
                     options={'maxiter': 20000, 'xatol': 1e-7, 'fatol': 1e-7})
        return r.fun
    h = 0.02
    f0, fp, fm = best.fun, prof(q[1] + h), prof(q[1] - h)
    curv = (fp + fm - 2 * f0) / h ** 2
    sig = float(np.exp(q[1]))
    err = float(sig / np.sqrt(curv)) if curv > 0 else float('nan')
    ns = float(np.exp(q[0]))
    r3 = 3 * sig
    bk3 = float(np.exp(q[2]) * r3 ** 2 / 2 + np.exp(q[3]) * r3)
    return {'sigma_res': round(sig, 2), 'err': round(err, 2), 'n_sig': round(ns, 1), 'n_bkg_r<3s': round(bk3, 1),
            'n_fit': int(tot), 'fit_max_mm': rmax}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--summary', help='muon_funnel_summary.py output (groups) or one muon_funnel.py run json')
    ap.add_argument('--group', default='REF', help='group used for --source data (default REF)')
    ap.add_argument('--groups', default='', help='comma list of groups to fit and print (default: all)')
    ap.add_argument('--source', choices=('model', 'data'), default='model')
    ap.add_argument('--sigma-tag', type=float, default=5.0, help='one-chamber resolution per axis, mm (SPS: ~5)')
    ap.add_argument('--contain', type=float, default=0.99, help='2-D containment of the probe search radius')
    ap.add_argument('--margin-nsig', type=float, default=2.0, help='fiducial edge margin in units of sigma_point')
    ap.add_argument('--min-r', type=float, default=0.0, help='floor of the probe search radius, mm')
    ap.add_argument('--fit-max', type=float, default=100.0, help='fit range of the radial distance, mm')
    ap.add_argument('--out', default='')
    a = ap.parse_args()

    k = k_contain(a.contain)
    s = {c: a.sigma_tag for c in PROBES}
    print(f'model, one-chamber resolution {a.sigma_tag} mm per axis (z = 320 / 630 / 940 mm):')
    mod = {}
    for p in PROBES:
        wa, wb = factors(p)
        sp, sr = sigma_model(p, s)
        mod[p] = {'w': [wa, wb], 'sigma_point': round(sp, 2), 'sigma_res': round(sr, 2)}
        print(f'  {p:7s} tags {TAGS[p][0]}+{TAGS[p][1]}  p = {wa:+.1f} a {wb:+.1f} b   sigma_point {sp:5.2f}  '
              f'sigma_res {sr:5.2f} mm   15 mm = {15 / sr:.2f} sigma_res   contains '
              f'{1 - np.exp(-15 ** 2 / (2 * sr ** 2)):.3f}')
    # equal-chamber ratios sigma_point / sigma_res
    rp = {p: mod[p]['sigma_point'] / mod[p]['sigma_res'] for p in PROBES}

    meas = {}
    if a.summary:
        S = json.load(open(a.summary))
        G = S['groups'] if 'groups' in S else {'run': S}
        want = [g for g in a.groups.split(',') if g] or list(G)
        print(f'\nmeasured radial distance of the nearest probe hit (Rayleigh + chance, r < {a.fit_max:.0f} mm):')
        for gname in want:
            if gname not in G:
                continue
            meas[gname] = {}
            line = []
            for p in PROBES:
                if p not in G[gname]['probes']:
                    continue
                e, n = hist(G[gname]['probes'][p])
                if n.sum() < 50:
                    continue
                r = fit(e, n, a.fit_max)
                r['contain_15mm'] = round(float(1 - np.exp(-15 ** 2 / (2 * r['sigma_res'] ** 2))), 3)
                meas[gname][p] = r
                line.append(f'{p[3:]} {r["sigma_res"]:5.2f}+-{r["err"]:.2f}')
            m_ = meas[gname]
            rat = ''
            if 'P2_MID' in m_:
                rat = '  ratio/MID ' + ' '.join(f'{p[3:]} {m_[p]["sigma_res"] / m_["P2_MID"]["sigma_res"]:.2f}'
                                               for p in ('P2_IN', 'P2_OUT') if p in m_)
                rat += f'  -> one chamber s = {m_["P2_MID"]["sigma_res"] / np.sqrt(1.5):.2f} mm (from MID)'
            print(f'  {gname:28s} ' + '  '.join(line) + rat)

    chosen = {}
    for p in PROBES:
        if a.source == 'data':
            if a.group not in meas or p not in meas[a.group]:
                sys.exit(f'no fit for {a.group} {p}')
            sr = meas[a.group][p]['sigma_res']
            sp = sr * rp[p]
        else:
            sr, sp = mod[p]['sigma_res'], mod[p]['sigma_point']
        chosen[p] = {'sigma_res': round(sr, 2), 'sigma_point': round(sp, 2),
                     'probe_r': round(max(a.min_r, k * sr), 1), 'margin_mm': round(a.margin_nsig * sp, 1)}
    print(f'\nchosen ({a.source}{" " + a.group if a.source == "data" else ""}; probe_r = max({a.min_r}, '
          f'{k:.2f} sigma_res) for {a.contain:.0%} containment, margin = {a.margin_nsig} sigma_point):')
    for p, c in chosen.items():
        print(f'  {p:7s} sigma_res {c["sigma_res"]:5.2f}  sigma_point {c["sigma_point"]:5.2f}  probe_r '
              f'{c["probe_r"]:5.1f} mm  margin {c["margin_mm"]:5.1f} mm')
    if a.out:
        json.dump({'params': vars(a), 'k_contain': k, 'model': mod, 'measured': meas, 'probes': chosen},
                  open(a.out, 'w'), indent=1)
        print('->', a.out)


if __name__ == '__main__':
    sys.exit(main())
