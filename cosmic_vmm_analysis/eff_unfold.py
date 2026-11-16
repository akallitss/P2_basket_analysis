#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ABSOLUTE single-channel efficiency with the tracking resolution fitted out (wk2 global-threshold scan).

Problem: the own-pad efficiency of efficiency_tp.py = eps x P(the prediction picked the pad the muon really crossed);
the second factor (resolution vs the ~12 mm pads) is unknown and dominates. It is NOT a degeneracy, because the pad
size is known: look at every live pad around the prediction, not only the nearest one, as a function of the offset
(u tangential, v radial, pad frame) of the predicted point from that pad's centre. The probability that the pad fires,

    P(u, v) = b + eps * A(u, v),
    A(u, v) = [Phi((W/2 - u)/sig) - Phi((-W/2 - u)/sig)] * [Phi((H/2 - v)/sig) - Phi((-H/2 - v)/sig)],
    W = w + 2 s, H = h + 2 s,

is a box of the known pad size (w, h; s = charge sharing / inclined tracks widening the response) smeared by the
prediction error sig, plus random hits b. Far from the pad P -> b, in the middle P -> b + eps A(0,0); the width of the
profile measures sig because w and h are known, its height then gives eps. Fitted per run and probe chamber on all
normal pads (blind ones excluded), off-time tags subtracted (on - r_scale * off). Per channel j:
eps_j = (K_j - b N_j) / sum_pairs A(u, v) with the fitted (sig, s) - the expected number of times it was crossed.

Inputs: tags_runN.npz (efficiency_tp.py --dump-tags), eff_runN.json (same run: live pads, title), the map/sdt files
of eff_vs_threshold_distance.py for the channel distance. Output OUT_DIR: fits.csv, channels.csv, profile figures,
eff_abs_vs_distance.png.

Usage: eff_unfold.py TAGS_DIR MAP_DIR OUT_DIR [--rmax 30]
"""
import argparse
import glob
import json
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from scipy.special import ndtr  # noqa: E402

import orientation_check as O  # noqa: E402
from event_display import style, INK, MUTED, COL, MRK  # noqa: E402
import vmm_stations as VS  # noqa: E402  (on the path through orientation_check)

ORDERS = ['P2_IN:3=bot:rev|top:lin', 'P2_IN:5=bot:rev|top:lin']
DETV = {3: range(2, 8), 2: range(8, 14), 1: range(14, 20)}
NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}


def pad_geometry():
    """(32, 64) arrays: x, y, phi, w, h of each read-out channel (same map as the efficiency runs)."""
    X, Y = O.Mapper(-1).xy_table(O.parse_orders(ORDERS))
    p = pd.read_csv(VS.MAP_CSV)
    tree = cKDTree(np.c_[p.pad_cx, p.pad_cy])
    G = {k: np.full((32, 64), np.nan) for k in ('phi', 'w', 'h')}
    ok = np.isfinite(X)
    d, i = tree.query(np.c_[X[ok], Y[ok]])
    assert d.max() < 0.01, d.max()
    for k, col in (('phi', 'phi'), ('w', 'pad_w'), ('h', 'pad_h')):
        G[k][ok] = p[col].to_numpy()[i]
    return X, Y, G


def box(u, half, sig):
    return ndtr((half - u) / sig) - ndtr((-half - u) / sig)


def acc(u, v, w, h, sig, s):
    return box(u, w / 2 + s, sig) * box(v, h / 2 + s, sig)


def pairs_for(tags, pads, rmax):
    """all (tag, pad) pairs within rmax: offsets in the pad frame, weight, fired flag."""
    px, py = tags['px'], tags['py']
    tree = cKDTree(np.c_[pads.x, pads.y])
    lst = tree.query_ball_point(np.c_[px, py], rmax)
    ti = np.repeat(np.arange(len(px)), [len(x) for x in lst])
    pj = np.concatenate([np.asarray(x, int) for x in lst]) if len(ti) else np.zeros(0, int)
    dx, dy = px[ti] - pads.x.to_numpy()[pj], py[ti] - pads.y.to_numpy()[pj]
    phi = pads.phi.to_numpy()[pj]
    v = dx * np.cos(phi) + dy * np.sin(phi)       # radial
    u = -dx * np.sin(phi) + dy * np.cos(phi)      # tangential
    # fired: the tag's fired pad centres, matched to pads by position
    key = lambda x, y: np.round(x * 10).astype(np.int64) * 100000 + np.round(y * 10).astype(np.int64)  # noqa: E731
    fired = set(zip(tags['fired_i'].tolist(), key(tags['fired_x'], tags['fired_y']).tolist()))
    pk = key(pads.x.to_numpy(), pads.y.to_numpy())[pj]
    f = np.fromiter(((a, b) in fired for a, b in zip(ti.tolist(), pk.tolist())), bool, len(ti))
    wt = np.where(tags['on'][ti], 1.0, -float(tags['r_scale']))
    return pd.DataFrame({'t': ti, 'j': pj, 'u': u, 'v': v, 'w': wt, 'f': f,
                         'pw': pads.w.to_numpy()[pj], 'ph': pads.h.to_numpy()[pj]})


def select_tags(tags, m):
    """keep the tags in mask m (fired-pad lists re-indexed)."""
    new = np.full(len(m), -1)
    new[m] = np.arange(m.sum())
    keep = m[tags['fired_i']]
    out = {k: (v[m] if np.ndim(v) and len(v) == len(m) else v) for k, v in tags.items()
           if not k.startswith('fired_')}
    out['fired_i'] = new[tags['fired_i'][keep]]
    out['fired_x'], out['fired_y'] = tags['fired_x'][keep], tags['fired_y'][keep]
    return out


def acc2(u, v, w, h, sig, s, sig2, ft):
    """core + tail resolution: (1 - ft) * box*G(sig) + ft * box*G(sig2)."""
    return (1 - ft) * acc(u, v, w, h, sig, s) + ft * acc(u, v, w, h, sig2, s)


def fit(pr, bw=1.5, rmax=30, s_fixed=None, tail=False):
    """binned least squares of P(u, v) = b + eps A on the off-time-subtracted pair counts (s fixed if given; with
    tail, the resolution is a core Gaussian sig plus a wider one sig2 with fraction ft)."""
    e = np.arange(-rmax, rmax + bw / 2, bw)
    iu, iv = np.digitize(pr.u, e) - 1, np.digitize(pr.v, e) - 1
    nb = len(e) - 1
    ok = (iu >= 0) & (iu < nb) & (iv >= 0) & (iv < nb)
    idx = iu[ok] * nb + iv[ok]
    w, f = pr.w.to_numpy()[ok], pr.f.to_numpy()[ok]
    N = np.bincount(idx, w, nb * nb)
    K = np.bincount(idx, w * f, nb * nb)
    N2 = np.bincount(idx, w ** 2, nb * nb)   # on + r^2 off: the variance weight of a bin
    c = (e[:-1] + e[1:]) / 2
    U, V = np.meshgrid(c, c, indexing='ij')
    U, V = U.ravel(), V.ravel()
    W, H = pr.pw.median(), pr.ph.median()
    m = N > 2
    names = ['eps', 'sig'] + (['s'] if s_fixed is None else []) + (['sig2', 'ft'] if tail else []) + ['b']
    p0 = {'eps': 0.6, 'sig': 3.0, 's': 1.0, 'sig2': 12.0, 'ft': 0.2, 'b': 0.01}
    lo = {'eps': 0, 'sig': 0.5, 's': -3, 'sig2': 3, 'ft': 0, 'b': 0}
    hi = {'eps': 1.5, 'sig': 30, 's': 15, 'sig2': 60, 'ft': 0.95, 'b': 1}

    def par(q):
        d = dict(zip(names, q))
        d.setdefault('s', s_fixed)
        d.setdefault('sig2', 1.0)
        d.setdefault('ft', 0.0)
        return d

    def res(q):
        # Pearson: binomial variance from the MODEL (the observed counts as variance bias the amplitude low)
        d = par(q)
        P = d['b'] + d['eps'] * acc2(U[m], V[m], W, H, d['sig'], d['s'], d['sig2'], d['ft'])
        return (K[m] - N[m] * P) / np.sqrt(N2[m] * np.clip(P * (1 - P), 1e-4, None))
    r = least_squares(res, [p0[n] for n in names], bounds=([lo[n] for n in names], [hi[n] for n in names]))
    npar = len(r.x)
    cov = np.linalg.pinv(r.jac.T @ r.jac) * max(1.0, (r.fun ** 2).sum() / max(m.sum() - npar, 1))
    d = par(r.x)
    err = {f'e_{n}': float(np.sqrt(cov[i, i])) for i, n in enumerate(names)}
    err.setdefault('e_s', 0.0)
    return dict(**d, **err, chi2ndf=(r.fun ** 2).sum() / max(m.sum() - npar, 1), W=W, H=H, prof=(U, V, N, K, m))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('tags_dir')
    ap.add_argument('map_dir')
    ap.add_argument('out_dir')
    ap.add_argument('--rmax', type=float, default=30)
    ap.add_argument('--scint', default='', help="LO:HI ns: keep only tags with a scintillator coincidence (H0 ch40) "
                    "LO..HI ns before the tag (muons: ~150:300)")
    ap.add_argument('--tail', action='store_true', help='core + tail (two-Gaussian) resolution')
    ap.add_argument('--s', default='mid', help="pad-response widening s (mm): 'mid' = weighted mean of the free "
                    "P2_MID fits (the interpolated chamber, where s and sig separate), 'free', or a number")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    X, Y, G = pad_geometry()
    mp = json.load(open(f'{a.map_dir}/map_Km18_i5_unmasked.json'))
    cal = json.load(open(f'{a.map_dir}/sdtcal.json'))['fit']
    rsdt = json.load(open(f'{a.map_dir}/run_sdt.json'))

    def dist(v, c, sdt):
        return (mp['threshold'][str(v)][str(c)] + (sdt - mp['sdt'][str(v)]) * cal[str(v)]['mV_per_sdt']
                - mp['baseline'][str(v)][str(c)])

    fits, chans, profs, store = [], [], {}, []
    for fz in sorted(glob.glob(f'{a.tags_dir}/tags_run*.npz')):
        run = re.search(r'run(\d+)', fz)[1]
        E = json.load(open(f'{a.tags_dir}/eff_run{run}.json'))
        tag = E['title'].split(' ', 1)[1]
        K = int(tag.split('_adwk2_K')[1].split('_')[0].replace('p', '+').replace('m', '-'))
        sdt = rsdt[f'p2b-config-cern-ext_gain4.5_peaktime200_opt_cosmic_{tag}.txt']['sdt']
        Z = np.load(fz)
        dead = {tuple(map(int, s.split(':'))) for lst in E['dead_or_masked'].values() for s in lst}
        for p in (3, 2, 1):
            tags = {k.split('_', 1)[1]: Z[k] for k in Z.files if k.startswith(f'p{p}_')}
            if a.scint:
                lo, hi = map(float, a.scint.split(':'))
                tags = select_tags(tags, (tags['dt_scint'] >= lo) & (tags['dt_scint'] <= hi))
            rows = [(v, c, X[v, c], Y[v, c], G['phi'][v, c], G['w'][v, c], G['h'][v, c], dist(v, c, sdt[str(v)]))
                    for v in DETV[p] for c in range(64) if np.isfinite(X[v, c]) and (v, c) not in dead]
            pads = pd.DataFrame(rows, columns=['vmm', 'ch', 'x', 'y', 'phi', 'w', 'h', 'dist'])
            pr = pairs_for(tags, pads, a.rmax)
            store.append((run, K, p, tags, pads, pr))
    # pass 1: s from the interpolated chamber
    if a.s in ('mid', 'free'):
        ss = []
        for run, K, p, tags, pads, pr in store:
            if p == 2:
                fr = fit(pr[(pads.dist >= 0).to_numpy()[pr.j]], rmax=a.rmax)
                ss.append((fr['s'], fr['e_s']))
                print(f'pass 1 run_{run} P2_MID free: s {fr["s"]:.2f}+-{fr["e_s"]:.2f} mm, sig {fr["sig"]:.2f} mm, '
                      f'eps {fr["eps"]:.3f}')
        ss = np.array(ss)
        s_mid = float(np.average(ss[:, 0], weights=1 / np.maximum(ss[:, 1], 0.05) ** 2))
        print(f'P2_MID s (weighted mean over {len(ss)} runs): {s_mid:.2f} mm')
    s_use = None if a.s == 'free' else (s_mid if a.s == 'mid' else float(a.s))
    for run, K, p, tags, pads, pr in store:
        if True:
            normal = (pads.dist >= 0).to_numpy()
            fr = fit(pr[normal[pr.j]], rmax=a.rmax, s_fixed=s_use, tail=a.tail)
            profs[(int(run), p)] = fr.pop('prof')
            fits.append(dict(run=int(run), K=K, det=NAME[p], tags_on=int(tags['on'].sum()),
                             tags_off=int((~tags['on']).sum()), **fr))
            print(f'run_{run} K={K:+d} {NAME[p]:6s}: eps {fr["eps"]:.3f}+-{fr["e_eps"]:.3f}  sig {fr["sig"]:.1f}+-'
                  f'{fr["e_sig"]:.1f} mm  s {fr["s"]:.1f} mm  tail {fr["ft"]:.2f} x {fr["sig2"]:.1f} mm  b {fr["b"]:.4f}  '
                  f'chi2/ndf {fr["chi2ndf"]:.2f}', flush=True)
            A = acc2(pr.u.to_numpy(), pr.v.to_numpy(), pr.pw.to_numpy(), pr.ph.to_numpy(), fr['sig'], fr['s'],
                     fr['sig2'], fr['ft'])
            g = pd.DataFrame({'j': pr.j, 'N': pr.w, 'K': pr.w * pr.f, 'A': pr.w * A,
                              'A2': pr.w ** 2 * A, 'V': pr.w ** 2 * pr.f}).groupby('j').sum()
            g = g.join(pads, how='left')
            g['eps'] = (g.K - fr['b'] * g.N) / g.A
            g['err'] = np.sqrt(g.V.clip(1)) / g.A.clip(1e-9)
            g['run'], g['K_scan'], g['det'] = int(run), K, NAME[p]
            gn = g[g.dist >= 0]
            fits[-1]['eps_count'] = float((gn.K - fr['b'] * gn.N).sum() / gn.A.sum())
            chans.append(g.reset_index(drop=True))
    F = pd.DataFrame(fits)
    F.to_csv(f'{a.out_dir}/fits.csv', index=False)
    C = pd.concat(chans)
    C.to_csv(f'{a.out_dir}/channels.csv', index=False)

    # profiles of one run per chamber (K=0, run 218): data vs model along u in the |v| < 3 mm band, and along v
    prun = 218 if (218, 3) in profs else min(profs)[0]
    fig, ax = plt.subplots(2, 3, figsize=(11, 5.5))
    for col, p in enumerate((3, 2, 1)):
        U, V, N, Kc, m = profs[(prun, p)]
        fr = F[(F.run == prun) & (F.det == NAME[p])].iloc[0]
        for row, (x, y) in enumerate(((U, V), (V, U))):
            band = np.abs(y) < 3
            xs = np.unique(x[band])
            dn = np.array([N[band & (x == q)].sum() for q in xs])
            dk = np.array([Kc[band & (x == q)].sum() for q in xs])
            mod = np.array([(N[band & (x == q)] * (fr.b + fr.eps * acc2(U[band & (x == q)], V[band & (x == q)],
                                                                          fr.W, fr.H, fr.sig, fr.s, fr.sig2, fr.ft))).sum() for q in xs])
            ok = dn > 20
            ax[row, col].errorbar(xs[ok], dk[ok] / dn[ok], np.sqrt(np.abs(dk[ok])) / dn[ok], fmt=MRK[p], color=COL[p],
                                  ms=3, lw=1, label='data')
            ax[row, col].plot(xs[ok], mod[ok] / dn[ok], color=INK, lw=1.2, label='fit')
            ax[row, col].set_xlabel(('u (tangential)' if row == 0 else 'v (radial)') + ' offset from pad centre (mm)',
                                    fontsize=8, color=INK)
            style(ax[row, col])
        ax[0, col].set_title(f'{NAME[p]} run_{prun}: eps {fr.eps:.2f}, sig {fr.sig:.1f} mm, s {fr.s:.1f} mm',
                             fontsize=8, color=INK, loc='left')
    ax[0, 0].set_ylabel('P(pad fires)', fontsize=8, color=INK)
    ax[1, 0].set_ylabel('P(pad fires)', fontsize=8, color=INK)
    ax[0, 0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(f'{a.out_dir}/profiles_run{prun}.png', dpi=150)

    # absolute efficiency vs threshold - baseline, per chamber (channels pooled: sum of (K - b N) / sum of A)
    C = C.merge(F[['run', 'det', 'b']], on=['run', 'det'])
    C['Kc'] = C.K - C.b * C.N
    bins = np.r_[-500, 0, 14, 18, 22, 26, 30, 34, 40, 46, 52, 58, 64, 70, 100]
    C['dbin'] = pd.cut(C.dist, bins)
    T = C.groupby(['det', 'dbin'], observed=True).agg(Kc=('Kc', 'sum'), A=('A', 'sum'), V=('V', 'sum'),
                                                      dist=('dist', 'mean'), ch=('ch', 'size')).reset_index()
    T['eps'] = T.Kc / T.A
    T['err'] = np.sqrt(T.V) / T.A
    T.astype({'dbin': str}).to_csv(f'{a.out_dir}/eff_abs_vs_distance.csv', index=False)
    pd.set_option('display.width', 200)
    print(T.pivot(index='dbin', columns='det', values='eps').round(3))
    print(T.pivot(index='dbin', columns='det', values='err').round(3))
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for p in (3, 2, 1):
        s = T[(T.det == NAME[p]) & (T.dist >= 0)]
        ax.errorbar(s.dist, s.eps, s.err, color=COL[p], marker=MRK[p], ms=4, lw=1.5, label=NAME[p])
    ax.set_xlabel('threshold - baseline of the channel (mV)', fontsize=8, color=INK)
    ax.set_ylabel('single-channel efficiency (resolution unfolded)', fontsize=8, color=INK)
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=7)
    style(ax)
    fig.tight_layout()
    fig.savefig(f'{a.out_dir}/eff_abs_vs_distance.png', dpi=150)


if __name__ == '__main__':
    main()
