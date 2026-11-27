#!/usr/bin/env python3
"""45_deform_local_z.py -- local effective charge-collection height z_eff of
the chamber, region by region, from the M3 reference tracks (det2 deformation
audit, purely geometric piece).

PRINCIPLE. The reference track is a straight line. If the chamber "sees" the
track at z_eff instead of the reference plane z0, the measured position minus
the track position at z0 is (z_eff - z0) * slope, in x and in y separately.
dz = z_eff - z0 is therefore the slope dependence of the position residual.
z_eff is the mean height of the ionisation a pad collects (~ mesh + gap/2), so
dz moves 1:1 with a bowed board and by half a local change of the drift gap.

SIGN CONVENTION. dz > 0 means z_eff is at LARGER M3 z than z0, i.e. towards
the upper M3 doublet (z = 1185/1302 mm; the lower doublet is at 24/144 mm and
the chamber at ~246 mm). The run_config gives det_orientation (0,0,0) and no
window side; it is measured here from the pad timing of crossing tracks (the
pad holding the larger-z ionisation fires ~50 ns later): the window is ABOVE
(+z), so dz > 0 = towards the window (thicker gap / board bowing up).

z0. The trimmed stage-34 alignment (alignment.npz, z = 245.98 mm) is used; it
already fitted ONE global plane z with the same residual-vs-slope estimator as
method 1, so the global centroid dz is ~0 by construction. Method 2 is an
independent measurement of the global z_eff.

TWO ESTIMATORS
 1. CENTROID: per region and per axis, OLS of (CC20 cluster centroid - track at
    z0) on the slope, trimmed |r| < 15 mm (as the alignment), bootstrap errors.
    Selection and region by the MEASURED max pad (a cut on the track position
    near dead pads biases dz by -12 %). The pad quantisation x the telescope
    slope-position correlation (corr -0.5) biases the raw dz by +0.4 mm
    (-0.8..+0.9 per reg24 region): corrected track by track with a paired
    geometric toy (same track, dz = 0, live tile containing the point).
 2. BOUNDARY: for every pair of touching live tiles the shared edge (midline of
    the two facing tile edges, overlap only, 1.5 mm kept from the ends) is built
    in the pad frame; the signed distance u of the track (at z0) to it (+ into
    pad B) and s_perp = slope . n are computed in M3 mm. The probability that B
    is the larger-amplitude pad (among pads within 15 mm of the track, the max
    must be A or B) is fitted by maximum likelihood with
        P = eps + (1 - 2 eps) Phi((u + dz s_perp - o . n) / w),
        w = sqrt(w0^2 + (k s_perp)^2),
    o = (ox, oy) a local in-plane translation of the region (sign-consistent
    with any pair ordering). (w0, k, eps) come from the whole-chamber fit and
    are fixed in the regions; (dz, ox, oy) free per region. Continuous variant
    (FRAC): when both pads fire, least squares of f = a_B/(a_A+a_B) to the same
    Phi.

VALIDATION: slope shuffle (null), re-projection to z0 + 3 mm (expect dz - 3),
toy injection of dz_true = +3 mm in both methods, splits in time and in track
angle, x vs y (centroid) and edge orientation (boundary), and the
slope-position (acceptance) correlation with a position-term variant.

GLOBAL MODELS on the per-reg48 dz: constant, plane, plane + isotropic bow,
full quadratic; chi2 and the significance of each added term.

Usage: python3 45_deform_local_z.py [run_key]
Out:   <run>/45_deform_local_z/{summary.txt, dz_regions.csv, dz_reg24.csv,
       global_models.csv, validation.csv, edges.csv, map_reg24.png,
       residual_vs_slope.png, boundary_step.png, methods_compare.png}
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import ndtr, log_ndtr
from scipy.spatial import cKDTree
from scipy.stats import chi2 as chi2d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_study_cache as sc     # noqa: E402
import p2_study_common as cm    # noqa: E402
import p2_align as pa           # noqa: E402

GRID = '#e6e6e6'
C1, C2, C3, C4 = '#9c9c9c', '#257de0', '#ec4647', '#1baf7a'
RNG = np.random.default_rng(20261127)
TRIM = 15.0          # mm, residual trim (as the alignment)
SD_MIN = 3.0         # mm inside the live area for the centroid method
U_MAX = 4.0          # mm, boundary window
END_MARGIN = 1.5     # mm kept away from the ends of a shared edge
NBOOT = 200
REG_TYPES = ['all', 'reg8', 'reg24', 'reg48', 'connector_N']


def inverse(T, xy):
    """M3 -> pad frame (inverse of p2_align.Transform.apply)."""
    return ((np.asarray(xy) - T.m3_mean) / T.s) @ T.R + T.pad_mean


# ======================================================================= centroid
def ols_dz(r, t, pos=None):
    """dz (and intercept) from r = a + dz t [+ b.pos], trimmed |r| < TRIM."""
    k = np.abs(r) < TRIM
    cols = [np.ones(k.sum()), t[k]] + ([] if pos is None else [p[k] for p in pos])
    A = np.column_stack(cols)
    p, *_ = np.linalg.lstsq(A, r[k], rcond=None)
    return p[1], p[0]


def boot_dz(r, t, nb=NBOOT):
    k = np.abs(r) < TRIM
    r, t = r[k], t[k]
    n = len(r)
    if n < 30:
        return np.nan, np.nan, n
    dz0 = np.cov(t, r)[0, 1] / np.var(t, ddof=1)
    idx = RNG.integers(0, n, (nb, n))
    tb, rb = t[idx], r[idx]
    tm, rm = tb.mean(1, keepdims=True), rb.mean(1, keepdims=True)
    dzb = ((tb - tm) * (rb - rm)).sum(1) / ((tb - tm) ** 2).sum(1)
    return dz0, dzb.std(ddof=1), n


def boot_dz_pos(r, t, X, Y, nb=NBOOT):
    """dz from r = a + bx (X-<X>) + by (Y-<Y>) + dz t (local alignment scale /
    rotation terms), trimmed, bootstrap error."""
    k = np.abs(r) < TRIM
    r, t, X, Y = r[k], t[k], X[k] - X[k].mean(), Y[k] - Y[k].mean()
    n = len(r)
    if n < 30:
        return np.nan, np.nan
    A = np.column_stack([np.ones(n), t, X, Y])
    p0 = np.linalg.lstsq(A, r, rcond=None)[0][1]
    pb = [np.linalg.lstsq(A[i], r[i], rcond=None)[0][1] for i in RNG.integers(0, n, (nb // 2, n))]
    return p0, np.std(pb, ddof=1)


def centroid_table(C, reg, pos=False):
    """C: dict with rx, ry, tx, ty, region columns. Returns per-region rows."""
    rows = []
    for rt in REG_TYPES:
        lab = np.zeros(len(C['rx']), int) if rt == 'all' else C[rt]
        for v in np.unique(lab):
            m = lab == v
            dx, ex, nx = boot_dz(C['rx'][m], C['tx'][m])
            dy, ey, ny = boot_dz(C['ry'][m], C['ty'][m])
            if not (np.isfinite(ex) and np.isfinite(ey)):
                continue
            w = np.array([1 / ex ** 2, 1 / ey ** 2])
            d = (w[0] * dx + w[1] * dy) / w.sum()
            row = dict(rtype=rt, region=int(v), n_c=int(m.sum()), dz_cx=dx, err_cx=ex,
                       dz_cy=dy, err_cy=ey, dz_c=d, err_c=1 / np.sqrt(w.sum()),
                       pull_xy=(dx - dy) / np.hypot(ex, ey))
            if pos:
                px, qx = boot_dz_pos(C['rx'][m], C['tx'][m], C['X'][m], C['Y'][m])
                py, qy = boot_dz_pos(C['ry'][m], C['ty'][m], C['X'][m], C['Y'][m])
                w2 = np.array([1 / qx ** 2, 1 / qy ** 2])
                row.update(dz_c_pos=(w2[0] * px + w2[1] * py) / w2.sum(), err_c_pos=1 / np.sqrt(w2.sum()))
            rows.append(row)
    return pd.DataFrame(rows)


# ======================================================================= boundary
def build_edges(tv, adj):
    """Shared edge of each touching tile pair (pad frame): centre c, direction d,
    normal n (A -> B), overlap interval [lo, hi] along d, gap, angle mismatch."""
    rows = []
    n_drop = 0
    for a, b in adj:
        qa, qb = tv[a], tv[b]
        best = None
        for i in range(4):          # facing edges: parallel, close, overlapping
            a0, a1 = qa[i], qa[(i + 1) % 4]
            da = (a1 - a0) / np.hypot(*(a1 - a0)); na = np.array([-da[1], da[0]])
            for j in range(4):
                b0, b1 = qb[j], qb[(j + 1) % 4]
                db = (b1 - b0) / np.hypot(*(b1 - b0))
                if abs(da @ db) < 0.995:
                    continue
                gp = abs((0.5 * (b0 + b1) - a0) @ na)
                ta = sorted([0.0, (a1 - a0) @ da]); tb = sorted([(b0 - a0) @ da, (b1 - a0) @ da])
                ov = min(ta[1], tb[1]) - max(ta[0], tb[0])
                if gp < 1.0 and ov > 2 * END_MARGIN + 1.0 and (best is None or gp < best[0]):
                    best = (gp, i, j)
        if best is None:
            n_drop += 1
            continue
        _, i, j = best
        ma = 0.5 * (qa + np.roll(qa, -1, 0)); mb = 0.5 * (qb + np.roll(qb, -1, 0))
        a0, a1 = qa[i], qa[(i + 1) % 4]; b0, b1 = qb[j], qb[(j + 1) % 4]
        da = (a1 - a0) / np.hypot(*(a1 - a0)); db = (b1 - b0) / np.hypot(*(b1 - b0))
        if da @ db < 0:
            db = -db
        mis = np.degrees(np.arccos(np.clip(da @ db, -1, 1)))
        d = (da + db) / np.hypot(*(da + db))
        c = 0.5 * (ma[i] + mb[j])
        n = np.array([-d[1], d[0]])
        if n @ (qb.mean(0) - qa.mean(0)) < 0:
            n = -n
        ta = sorted([(a0 - c) @ d, (a1 - c) @ d]); tb = sorted([(b0 - c) @ d, (b1 - c) @ d])
        lo, hi = max(ta[0], tb[0]), min(ta[1], tb[1])
        rows.append(dict(A=a, B=b, cx=c[0], cy=c[1], dx=d[0], dy=d[1], nx=n[0], ny=n[1],
                         lo=lo, hi=hi, gap=(mb[j] - ma[i]) @ n, mis_deg=mis))
    E = pd.DataFrame(rows)
    E.attrs['n_drop'] = n_drop
    return E


def edge_candidates(P, E):
    """(ray index, edge index, u_pad, along) for rays near a shared edge."""
    tree = cKDTree(E[['cx', 'cy']].to_numpy())
    rmax = float(np.max(np.maximum(-E.lo, E.hi))) + U_MAX + 1
    lists = tree.query_ball_point(P, rmax)
    ii = np.repeat(np.arange(len(P)), [len(l) for l in lists])
    ee = np.concatenate([np.asarray(l, int) for l in lists]) if len(ii) else np.empty(0, int)
    c = E[['cx', 'cy']].to_numpy()[ee]; d = E[['dx', 'dy']].to_numpy()[ee]
    n = E[['nx', 'ny']].to_numpy()[ee]
    q = P[ii] - c
    al = np.einsum('ij,ij->i', q, d); u = np.einsum('ij,ij->i', q, n)
    ok = (np.abs(u) < U_MAX + 2) & (al > E.lo.to_numpy()[ee] + END_MARGIN) & \
         (al < E.hi.to_numpy()[ee] - END_MARGIN)
    return ii[ok], ee[ok], u[ok]


def offset_design(B, affine=False):
    """In-plane offset of the edge crossing: n . (o + M (X - <X>)); o = local
    translation, M (affine) = local scale / rotation / shear of the alignment."""
    cols = [B['nx'], B['ny']]
    if affine:
        dX = (B['X'] - B['X'].mean()) / 100.0; dY = (B['Y'] - B['Y'].mean()) / 100.0
        cols += [B['nx'] * dX, B['nx'] * dY, B['ny'] * dX, B['ny'] * dY]
    return np.column_stack(cols)


def nll_step(p, u, s, D, y, fixed=None):
    m = D.shape[1]
    if fixed is None:
        dz, th, (w0, k, eps) = p[0], p[1:1 + m], p[1 + m:]
    else:
        dz, th = p[0], p[1:1 + m]
        w0, k, eps = fixed
    w = np.sqrt(w0 ** 2 + (k * s) ** 2)
    z = (u + dz * s - D @ th) / w
    eps = 0.5 / (1 + np.exp(-eps))          # logit -> (0, 0.5)
    P = eps + (1 - 2 * eps) * ndtr(z)
    P = np.clip(P, 1e-12, 1 - 1e-12)
    return -(y * np.log(P) + (1 - y) * np.log(1 - P)).sum()


def num_hess(f, p, h=None):
    p = np.asarray(p, float)
    h = np.full(len(p), 1e-3) if h is None else h
    H = np.zeros((len(p), len(p)))
    f0 = f(p)
    for i in range(len(p)):
        for j in range(i, len(p)):
            ei = np.zeros(len(p)); ej = np.zeros(len(p)); ei[i] = h[i]; ej[j] = h[j]
            if i == j:
                H[i, i] = (f(p + ei) - 2 * f0 + f(p - ei)) / h[i] ** 2
            else:
                H[i, j] = H[j, i] = (f(p + ei + ej) - f(p + ei - ej) - f(p - ei + ej)
                                     + f(p - ei - ej)) / (4 * h[i] * h[j])
    return H


def fit_step(B, fixed=None, p0=None, affine=False):
    u, s, y = B['u'], B['s'], B['y']
    D = offset_design(B, affine)
    m = D.shape[1]
    if p0 is None:
        p0 = [0.0] * (1 + m) + ([0.7, 2.0, -3.0] if fixed is None else [])
    elif len(p0) != 1 + m + (3 if fixed is None else 0):
        p0 = list(p0[:1]) + [0.0] * m + (list(p0[-3:]) if fixed is None else [])
    f = lambda p: nll_step(np.asarray(p), u, s, D, y, fixed)
    if fixed is None:
        r = minimize(f, p0, method='Nelder-Mead', options=dict(maxiter=20000, xatol=1e-5, fatol=1e-4))
        r = minimize(f, r.x, method='BFGS')
    else:
        r = minimize(f, p0, method='BFGS')
    H = num_hess(f, r.x)
    try:
        err = np.sqrt(np.clip(np.diag(np.linalg.inv(H)), 0, None))
    except np.linalg.LinAlgError:
        err = np.full(len(r.x), np.nan)
    return r.x, err


def lsq_frac(p, u, s, nx, ny, f, fixed):
    dz, ox, oy = p
    w0, k = fixed
    w = np.sqrt(w0 ** 2 + (k * s) ** 2)
    return ndtr((u + dz * s - ox * nx - oy * ny) / w) - f


def fit_frac(B, fixed=None):
    from scipy.optimize import least_squares
    u, s, nx, ny, f = (B[k] for k in ('u', 's', 'nx', 'ny', 'f'))
    if fixed is None:
        res = lambda p: lsq_frac(p[:3], u, s, nx, ny, f, p[3:])
        r = least_squares(res, [0, 0, 0, 0.7, 2.0])
    else:
        res = lambda p: lsq_frac(p, u, s, nx, ny, f, fixed)
        r = least_squares(res, [0, 0, 0])
    J = r.jac
    s2 = (r.fun ** 2).sum() / max(len(f) - len(r.x), 1)
    try:
        err = np.sqrt(np.diag(np.linalg.inv(J.T @ J) * s2))
    except np.linalg.LinAlgError:
        err = np.full(len(r.x), np.nan)
    return r.x, err


# ======================================================================= gap study helpers
GAP_EDGES = np.array([-np.inf, 3.5, 3.9, 4.3, 5.0, 6.0, np.inf])   # classes 0..5 of stage 47


def gap_class(g):
    g = np.asarray(g, float)
    return np.where(np.isfinite(g), np.digitize(g, GAP_EDGES) - 1, -1)


def nll_lin(p, u, s, Z, D, y, fixed):
    """Step likelihood with a position/gap-dependent dz = Z @ beta."""
    q = Z.shape[1]
    w0, k, eps = fixed
    w = np.sqrt(w0 ** 2 + (k * s) ** 2)
    z = (u + s * (Z @ p[:q]) - D @ p[q:]) / w
    e = 0.5 / (1 + np.exp(-eps))
    P = np.clip(e + (1 - 2 * e) * ndtr(z), 1e-12, 1 - 1e-12)
    return -(y * np.log(P) + (1 - y) * np.log(1 - P)).sum()


def fit_lin_b(B, Z, fixed, affine=True):
    """Boundary step with dz = Z @ beta (beta: a, b_gap, [gx, gy]); local affine
    alignment offsets as nuisance. Returns beta, cov(beta)."""
    D = offset_design(B, affine)
    q = Z.shape[1]
    f = lambda p: nll_lin(np.asarray(p), B['u'], B['s'], Z, D, B['y'], fixed)
    r = minimize(f, np.zeros(q + D.shape[1]), method='BFGS')
    H = num_hess(f, r.x)
    cov = np.linalg.inv(H)
    return r.x[:q], cov[:q, :q]


def fit_lin_f(F, Z, ffixed, affine=True):
    """Charge-fraction LSQ with dz = Z @ beta (local affine offsets free)."""
    from scipy.optimize import least_squares
    D = offset_design(F, affine)
    q = Z.shape[1]
    w0, k = ffixed
    w = np.sqrt(w0 ** 2 + (k * F['s']) ** 2)
    res = lambda p: ndtr((F['u'] + F['s'] * (Z @ p[:q]) - D @ p[q:]) / w) - F['f']
    r = least_squares(res, np.zeros(q + D.shape[1]))
    s2 = (r.fun ** 2).sum() / max(len(r.fun) - len(r.x), 1)
    cov = np.linalg.inv(r.jac.T @ r.jac) * s2
    return r.x[:q], cov[:q, :q]


def fit_lin_c(C, Zfun, axes=('x', 'y'), nb=100, rng=None):
    """Paired-centroid residuals: r_a = intercept(reg24, axis) + t_a * (Z @ beta),
    trimmed |r| < TRIM, stacked over the requested axes; bootstrap over tracks."""
    rng = RNG if rng is None else rng
    Z = Zfun(C)
    q = Z.shape[1]
    regs = np.unique(C['reg24'])
    ridx = np.searchsorted(regs, C['reg24'])
    blocks, rr, trk = [], [], []
    for ia, a in enumerate(axes):
        r = C['r' + a]; t = C['t' + a]
        k = np.abs(r) < TRIM
        A = np.zeros((k.sum(), len(regs) * len(axes) + q))
        A[np.arange(k.sum()), ia * len(regs) + ridx[k]] = 1
        A[:, len(regs) * len(axes):] = t[k, None] * Z[k]
        blocks.append(A); rr.append(r[k]); trk.append(np.nonzero(k)[0])
    A = np.vstack(blocks); r = np.concatenate(rr); trk = np.concatenate(trk)
    p = np.linalg.lstsq(A, r, rcond=None)[0][-q:]
    n = len(C['rx'])
    pb = []
    for _ in range(nb):
        w = np.bincount(rng.integers(0, n, n), minlength=n)[trk].astype(float)
        Aw = A * w[:, None]
        try:
            pb.append(np.linalg.solve(Aw.T @ A, Aw.T @ r)[-q:])
        except np.linalg.LinAlgError:
            pass
    return p, np.cov(np.array(pb).T)


def toy_clusters(n, G, cos_t, rng, ncl=2.5, nmax=60):
    """Primary ionisation clusters along a track in a gap of thickness G: number
    Poisson(ncl * G / cos theta), height fraction uniform, electrons with a
    1/k^2 tail (k = floor(1/U), <= 300). Returns (h_frac, electrons, valid)."""
    N = np.minimum(rng.poisson(ncl * G / cos_t), nmax)
    valid = np.arange(nmax)[None, :] < N[:, None]
    hf = rng.random((n, nmax))
    k = np.minimum(np.floor(1.0 / np.maximum(rng.random((n, nmax)), 1e-9)), 300.0)
    return hf, k * valid, valid


def amp_weight(d, Lg):
    """Amplitude weight of a cluster at drift distance d [mm]: 1 (pure charge)
    if Lg is None, else 1 for d < Lg and exp(-(d - Lg)/0.3 mm) beyond (a soft
    'only the early electrons count' gate, a stand-in for the observed fact that
    the two-pad amplitude ratio follows the mesh-side ionisation)."""
    if Lg is None:
        return np.ones_like(d)
    return np.where(d < Lg, 1.0, np.exp(-(d - Lg) / 0.3))


def toy_boundary(Bsrc, G, zm, Lg, thr, rng, n_rep=3, sig_p=0.55, sig_d=0.45, fire='amp'):
    """Boundary samples simulated on the real (u, s_perp, theta) of the data: the
    clusters sit at u + s*h (h = height relative to z0, mesh at zm, cathode at
    zm + G); pad B collects u > 0. Fired: amplitude > thr. Label: larger
    amplitude among fired pads (events with no fired pad dropped)."""
    u = np.tile(Bsrc['u'], n_rep); s = np.tile(Bsrc['s'], n_rep)
    ct = np.cos(np.radians(np.tile(Bsrc['theta'], n_rep)))
    n = len(u)
    hf, q, valid = toy_clusters(n, G, ct, rng)
    h = zm + G * hf
    pos = u[:, None] + s[:, None] * h + rng.normal(0, sig_p, n)[:, None] + rng.normal(0, sig_d, (n, hf.shape[1]))
    a = q * amp_weight(G * hf, Lg)
    aB = (a * (pos > 0)).sum(1); aA = (a * (pos <= 0)).sum(1)
    if fire == 'charge':          # a pad fires on its total charge, the amplitude is gated
        fA, fB = (q * (pos <= 0)).sum(1) > thr, (q * (pos > 0)).sum(1) > thr
    else:
        fA, fB = aA > thr, aB > thr
    keep = fA | fB
    y = np.where(fA & fB, aB > aA, fB).astype(float)
    both = fA & fB
    B = dict(u=u[keep], s=s[keep], y=y[keep], both=both[keep], nx=np.tile(Bsrc['nx'], n_rep)[keep],
             ny=np.tile(Bsrc['ny'], n_rep)[keep], X=np.tile(Bsrc['X'], n_rep)[keep],
             Y=np.tile(Bsrc['Y'], n_rep)[keep])
    F = {k_: v[B['both']] for k_, v in B.items()}
    with np.errstate(invalid='ignore', divide='ignore'):
        F['f'] = (aB / (aA + aB))[keep][B['both']]
    return B, F


def sub(B, m):
    return {k: v[m] for k, v in B.items()}


def boundary_table(B, fixed, Fq=None, ffixed=None, min_n=200):
    rows = []
    for rt in REG_TYPES:
        lab = np.zeros(len(B['u']), int) if rt == 'all' else B[rt]
        for v in np.unique(lab):
            m = lab == v
            if m.sum() < min_n:
                continue
            p, e = fit_step(sub(B, m), fixed)
            row = dict(rtype=rt, region=int(v), n_b=int(m.sum()), dz_b=p[0], err_b=e[0],
                       ox_b=p[1], oy_b=p[2])
            pa_, ea_ = fit_step(sub(B, m), fixed, affine=True)
            row.update(dz_b_aff=pa_[0], err_b_aff=ea_[0])
            m1 = m & ~B['both']
            if m1.sum() >= min_n:
                p1, e1 = fit_step(sub(B, m1), fixed)
                row.update(n_b1=int(m1.sum()), dz_b1=p1[0], err_b1=e1[0])
            if Fq is not None:
                mf = (Fq['all'] == 0) if rt == 'all' else (Fq[rt] == v)
                if mf.sum() >= 60:
                    pf, ef = fit_frac(sub(Fq, mf), ffixed)
                    row.update(n_f=int(mf.sum()), dz_f=pf[0], err_f=ef[0])
            rows.append(row)
    return pd.DataFrame(rows)


# ======================================================================= models
def model_fits(df, dzc, erc, xc, yc):
    d = df[np.isfinite(df[dzc]) & np.isfinite(df[erc]) & (df[erc] > 0)]
    X = (d.xm.to_numpy() - xc) / 100.0; Y = (d.ym.to_numpy() - yc) / 100.0   # per 100 mm
    z = d[dzc].to_numpy(); w = 1 / d[erc].to_numpy() ** 2
    one = np.ones(len(z))
    designs = {'constant': [one], 'plane': [one, X, Y], 'plane+bow': [one, X, Y, X * X + Y * Y],
               'quadratic': [one, X, Y, X * X, Y * Y, X * Y]}
    names = {'constant': ['c0'], 'plane': ['c0', 'gx', 'gy'], 'plane+bow': ['c0', 'gx', 'gy', 'bow'],
             'quadratic': ['c0', 'gx', 'gy', 'qxx', 'qyy', 'qxy']}
    out, prev = [], None
    for nm, cols in designs.items():
        A = np.column_stack(cols)
        Cv = np.linalg.inv(A.T @ (A * w[:, None]))
        p = Cv @ (A.T @ (w * z))
        chi = float((w * (z - A @ p) ** 2).sum()); ndf = len(z) - len(p)
        row = dict(model=nm, chi2=chi, ndf=ndf, p_value=float(chi2d.sf(chi, ndf)))
        if prev is not None:
            dchi = prev[0] - chi; dn = prev[1] - ndf
            row.update(dchi2=dchi, dndf=dn, p_extra=float(chi2d.sf(dchi, dn)))
        for k_, (pn, pv, pe) in enumerate(zip(names[nm], p, np.sqrt(np.diag(Cv)))):
            row[pn] = pv; row[pn + '_err'] = pe
        out.append(row); prev = (chi, ndf)
    return out


# ======================================================================= main
def main():
    run_key = sys.argv[1] if len(sys.argv) > 1 else 'det2_long_repos'
    hits, rays, g, cfg = sc.load(run_key)
    out = cfg.out_dir('45_deform_local_z')
    d34 = cfg.out_dir('34_clustering_matching')
    Al = np.load(f'{d34}/alignment.npz')
    T = cm.make_T(g, R=Al['R'], s=Al['s'], pad_mean=Al['pad_mean'], m3_mean=Al['m3_mean'])
    z0 = float(Al['z'])
    ra = pd.read_pickle(f'{d34}/rays_aligned.pkl')
    assert (ra.eventId.to_numpy() == rays.eventId.to_numpy()).all()
    ev = ra.eventId.to_numpy()
    tx, ty, theta = cm.ray_dirs(rays)
    x0, y0 = pa.project_to_z(rays, z0)
    xy0 = np.column_stack([x0, y0])
    assert np.abs(xy0 - ra[['x', 'y']].to_numpy()).max() < 1e-6
    sd = ra.sd.to_numpy()
    clu = pd.read_pickle(cfg.out_dir('34_study_cache') + '/clusters_cc20.pkl')
    reg = pd.read_csv(cfg.out_dir('43_deformation') + '/regions.csv').set_index('channel_id')
    pad_id = g['pad_id'].astype(int); tv = g['tile_verts']
    regs_of_row = {rt: reg.loc[pad_id, rt].to_numpy() for rt in REG_TYPES[1:]}
    L = [f'run {run_key}: {len(ev):,} M3 tracks; z0 = {z0:.3f} mm (stage-34 trimmed alignment; '
         f'nominal 246); rot {T.rotation_deg:.3f} deg, scale {float(T.s):.5f}',
         'SIGN: dz = z_eff - z0 > 0 means towards the upper M3 doublet (larger M3 z). The run_config '
         'records no window side (det_orientation 0,0,0); the window side is read from the pad timing below.']

    # time per event (for the time split)
    try:
        et = pa.load_event_times(cfg.m3_tracking_dir).set_index('eventId')['t_sec']
        tsec = et.reindex(ev).to_numpy()
        if np.isnan(tsec).mean() > 0.05:
            raise ValueError
        tsplit = tsec > np.nanmedian(tsec)
        tlab = 'M3 event time'
    except Exception:
        tsplit = ev > np.median(ev); tlab = 'eventId (time proxy)'

    # region of each ray: nearest live pad centre (pad frame) of the track at z0
    def region_of(xy):
        P = inverse(T, xy)
        _, row = cKDTree(g['pad_xy']).query(P)
        return {rt: regs_of_row[rt][row] for rt in REG_TYPES[1:]}, P, row
    rreg, P0, prow = region_of(xy0)
    # region centres in M3 (mean live-pad centre)
    padm3 = cm.apply_T(T, g['pad_xy'])
    cen = {rt: pd.DataFrame({'r': regs_of_row[rt], 'x': padm3[:, 0], 'y': padm3[:, 1]})
           .groupby('r').mean() for rt in REG_TYPES[1:]}
    xc, yc = padm3.mean(0)

    # hits lookup (pad rows), used by both methods
    h = hits[hits.channel_id.isin(pad_id)][['eventId', 'channel_id', 'amp', 'tot', 'time']].copy()
    h['channel_id'] = h.channel_id.astype(int)
    rowof = pd.Series(np.arange(len(pad_id)), index=pad_id)
    h['row'] = rowof.loc[h.channel_id.to_numpy()].to_numpy()

    # ================================================================ 1. centroid
    cxy = cm.apply_T(T, clu[['x_pad', 'y_pad']].to_numpy())
    dmin, j = cm.nearest(ev, xy0, cxy, clu['eventId'].to_numpy())
    # max-amplitude pad of each cluster (the selection and the region use the
    # MEASURED pad, not the track at z0: a cut on the track position next to a
    # dead pad / the edge removes a slope-dependent band of the pad and biases dz)
    cl = pd.DataFrame({'j': np.arange(len(clu)), 'eventId': clu.eventId.to_numpy(),
                       'amp': clu.qmax.to_numpy(), 'cx': clu.x_pad.to_numpy(), 'cy': clu.y_pad.to_numpy()})
    cl = cl.merge(h[['eventId', 'amp', 'row']], on=['eventId', 'amp'], how='left')
    cl['d'] = np.hypot(g['pad_xy'][cl.row.fillna(0).astype(int), 0] - cl.cx,
                       g['pad_xy'][cl.row.fillna(0).astype(int), 1] - cl.cy)
    cl = cl.sort_values('d').drop_duplicates('j').sort_values('j')
    crow = np.full(len(clu), -1); ok_ = cl.row.notna().to_numpy()
    crow[cl.j.to_numpy()[ok_]] = cl.row.to_numpy()[ok_].astype(int)
    selc = (dmin <= 20.0) & (j >= 0)
    selc[selc] &= crow[j[selc]] >= 0
    jj = j[selc]
    npad = clu.n_pad.to_numpy()[jj]
    creg = {rt: np.full(len(ev), -1) for rt in REG_TYPES[1:]}
    for rt in REG_TYPES[1:]:
        creg[rt][selc] = regs_of_row[rt][crow[jj]]

    def cent_inputs(xyref, slope_x=tx, slope_y=ty, m=None, rr=None):
        m = selc if m is None else m
        rr = creg if rr is None else rr
        r = cxy[j[m]] - xyref[m]
        C = dict(rx=r[:, 0], ry=r[:, 1], tx=slope_x[m], ty=slope_y[m],
                 X=xyref[m, 0], Y=xyref[m, 1])
        for rt in REG_TYPES[1:]:
            C[rt] = rr[rt][m]
        return C
    C0 = cent_inputs(xy0)
    L.append(f'\n== CENTROID method: {selc.sum():,} tracks matched (CC20 < 20 mm), region = max-amplitude '
             f'pad of the cluster, no track-position cut; {100*(npad == 1).mean():.1f} % single-pad '
             f'clusters; trim |r| < {TRIM:g} mm')
    TC = centroid_table(C0, None)

    # geometric toy: true point at z0 + dz_true*slope (+ 0.72 mm pointing) -> the LIVE
    # tile containing it (none: no hit) -> residual = pad centre - track at z0, |r| < 20
    pqd = __import__('importlib').import_module('34_clustering_matching').pt_quad_dist
    ptree = cKDTree(g['pad_xy'])

    def toy_centroid(dz_true, sig=0.72, track_cut=False):
        pt = xy0 + dz_true * np.column_stack([tx, ty]) + RNG.normal(0, sig, (len(ev), 2))
        Pp = inverse(T, pt)
        _, cand = ptree.query(Pp, k=4)
        rw = np.full(len(ev), -1)
        for kk in range(4):
            inside = (pqd(Pp, tv[cand[:, kk]]) == 0) & (rw < 0)
            rw[inside] = cand[inside, kk]
        m = rw >= 0
        if track_cut:
            m &= sd >= SD_MIN
        r = padm3[rw[m]] - xy0[m]
        k = np.hypot(r[:, 0], r[:, 1]) <= 20
        C = dict(rx=r[k, 0], ry=r[k, 1], tx=tx[m][k], ty=ty[m][k])
        for rt in REG_TYPES[1:]:
            C[rt] = regs_of_row[rt][rw[m][k]]
        return C
    V = []
    Ct0 = toy_centroid(0.0)
    TOY0 = centroid_table(Ct0, None)          # quantisation + acceptance bias per region
    toy_gain = []
    for dzt in (0.0, 3.0, -3.0):
        Ct = Ct0 if dzt == 0 else toy_centroid(dzt)
        dx, ex, _ = boot_dz(Ct['rx'], Ct['tx']); dy, ey, _ = boot_dz(Ct['ry'], Ct['ty'])
        toy_gain.append((dzt, 0.5 * (dx + dy)))
        V.append(dict(test=f'centroid toy (pad quantisation), dz_true {dzt:+g}', method='centroid',
                      dz=0.5 * (dx + dy), err=0.5 * np.hypot(ex, ey), expect=dzt,
                      note=f'x {dx:+.3f}+-{ex:.3f}, y {dy:+.3f}+-{ey:.3f}'))
    # old-style selection (track sd >= 3 mm, region of the track at z0) in the toy, for reference
    Cts = toy_centroid(3.0, track_cut=True)
    dx, ex, _ = boot_dz(Cts['rx'], Cts['tx']); dy, ey, _ = boot_dz(Cts['ry'], Cts['ty'])
    V.append(dict(test='centroid toy dz_true +3 with track sd >= 3 cut (old selection)', method='centroid',
                  dz=0.5 * (dx + dy), err=0.5 * np.hypot(ex, ey), expect=3.0,
                  note='shows why the selection is made on the measured pad'))
    # single-pad-only and multi-pad-only data subsets (whole chamber)
    for nm, mm in [('single-pad clusters', npad == 1), ('multi-pad clusters', npad > 1)]:
        dx, ex, _ = boot_dz(C0['rx'][mm], C0['tx'][mm]); dy, ey, _ = boot_dz(C0['ry'][mm], C0['ty'][mm])
        V.append(dict(test=f'centroid split: {nm}', method='centroid', dz=(dx / ex**2 + dy / ey**2) /
                      (1 / ex**2 + 1 / ey**2), err=1 / np.sqrt(1 / ex**2 + 1 / ey**2), expect=np.nan,
                      note=f'n {mm.sum()}, x {dx:+.3f}+-{ex:.3f}, y {dy:+.3f}+-{ey:.3f}'))
    # selection by the track position (sd >= 3 mm, as the efficiency studies) instead of the pad
    mm = sd[selc] >= SD_MIN
    dx, ex, _ = boot_dz(C0['rx'][mm], C0['tx'][mm]); dy, ey, _ = boot_dz(C0['ry'][mm], C0['ty'][mm])
    V.append(dict(test='centroid: track sd >= 3 mm selection', method='centroid',
                  dz=(dx / ex**2 + dy / ey**2) / (1 / ex**2 + 1 / ey**2),
                  err=1 / np.sqrt(1 / ex**2 + 1 / ey**2), expect=np.nan, note=f'n {mm.sum()}'))
    # PAIRED toy correction: for every selected data track the toy residual of the
    # SAME track (dz = 0, mean of K pointing replicas) is subtracted; quantisation
    # noise and the acceptance bias (slope-position correlation of the telescope
    # inside a pad) cancel track by track, the dz signal (data pad != toy pad near
    # pad edges) remains. Gain from the same pairing at dz_true = +-3 in the toy.
    def toy_track_resid(dz_true, K=8, sig=0.72):
        acc = np.zeros((len(ev), 2)); nn = np.zeros(len(ev))
        for _ in range(K):
            pt = xy0 + dz_true * np.column_stack([tx, ty]) + RNG.normal(0, sig, (len(ev), 2))
            Pp = inverse(T, pt)
            _, cand = ptree.query(Pp, k=4)
            rw = np.full(len(ev), -1)
            for kk in range(4):
                inside = (pqd(Pp, tv[cand[:, kk]]) == 0) & (rw < 0)
                rw[inside] = cand[inside, kk]
            m = rw >= 0
            acc[m] += padm3[rw[m]] - xy0[m]; nn[m] += 1
        with np.errstate(invalid='ignore'):
            return acc / nn[:, None]
    rt0 = toy_track_resid(0.0)[selc]
    okp = np.isfinite(rt0).all(1)
    Cp = {k: v[okp] for k, v in C0.items()}
    Cp['rx'] = C0['rx'][okp] - rt0[okp, 0]; Cp['ry'] = C0['ry'][okp] - rt0[okp, 1]
    gp = []
    for dzt in (1.0, -1.0):
        rtd = toy_track_resid(dzt)[selc][okp] - rt0[okp]
        dx, _, _ = boot_dz(rtd[:, 0], Cp['tx'], 20); dy, _, _ = boot_dz(rtd[:, 1], Cp['ty'], 20)
        gp.append(0.5 * (dx + dy))
    gain = (gp[0] - gp[1]) / 2.0
    TCP = centroid_table(Cp, None, pos=True)
    V.append(dict(test='centroid paired toy: gain (toy dz +-1, same tracks)', method='centroid', dz=gain,
                  err=np.nan, expect=1.0, note=f'+1 -> {gp[0]:+.3f}, -1 -> {gp[1]:+.3f}; '
                  f'{okp.mean()*100:.1f} % of data tracks paired'))
    # single-pad-only per region table (variant)
    TC1 = centroid_table({k: v[npad == 1] for k, v in C0.items()}, None)

    # acceptance: slope-position correlation inside regions (the paired toy carries it)
    acc_rows = []
    for v in np.unique(C0['reg24']):
        m = C0['reg24'] == v
        X = C0['X'][m] - C0['X'][m].mean(); Y = C0['Y'][m] - C0['Y'][m].mean()
        acc_rows.append(dict(region=int(v), corr_tx_x=np.corrcoef(C0['tx'][m], X)[0, 1],
                             corr_ty_y=np.corrcoef(C0['ty'][m], Y)[0, 1]))
    ACC = pd.DataFrame(acc_rows)

    # ================================================================ 2. boundary
    adj = np.array(__import__('importlib').import_module('34_clustering_matching').tile_adjacency(tv))
    E = build_edges(tv, adj)
    E.to_csv(f'{out}/edges.csv', index=False)
    L.append(f'\n== BOUNDARY method: {len(E)} shared edges (touching live tiles); gap median '
             f'{E.gap.median():.2f} mm (range {E.gap.min():.2f}..{E.gap.max():.2f}), facing-edge '
             f'angle mismatch median {E.mis_deg.median():.2f} deg (max {E.mis_deg.max():.2f}), overlap '
             f'median {(E.hi - E.lo).median():.1f} mm')

    # max-amp pad within 15 mm of each track; amplitude per (event, pad row)
    h['key'] = h.amp + 1e-4 * h.tot            # tot breaks ties of saturated pulses
    ray_df = pd.DataFrame({'i': np.arange(len(ev)), 'eventId': ev, 'px': P0[:, 0], 'py': P0[:, 1]})
    mh = ray_df.merge(h, on='eventId')
    pxy = g['pad_xy'][mh.row.to_numpy()]
    mh = mh[np.hypot(pxy[:, 0] - mh.px, pxy[:, 1] - mh.py) * float(T.s) < 15.0]
    best = mh.loc[mh.groupby('i')['key'].idxmax()]
    maxrow = np.full(len(ev), -1); maxrow[best.i.to_numpy()] = best.row.to_numpy()
    hk = pd.Series(h.key.to_numpy(), index=pd.MultiIndex.from_arrays([h.eventId.to_numpy(),
                                                                       h.row.to_numpy()]))
    hk = hk[~hk.index.duplicated()]
    ha = pd.Series(h.amp.to_numpy(), index=hk.index if False else pd.MultiIndex.from_arrays(
        [h.eventId.to_numpy(), h.row.to_numpy()]))
    ha = ha[~ha.index.duplicated()]
    ht = pd.Series(h.time.to_numpy(), index=pd.MultiIndex.from_arrays([h.eventId.to_numpy(),
                                                                       h.row.to_numpy()]))
    ht = ht[~ht.index.duplicated()]

    Rm = T.R

    def boundary_inputs(xyref, slx=tx, sly=ty, rr=None, Pref=None):
        rr = rreg if rr is None else rr
        Pp = inverse(T, xyref) if Pref is None else Pref
        ii, ee, upad = edge_candidates(Pp, E)
        A_ = E.A.to_numpy()[ee]; B_ = E.B.to_numpy()[ee]
        kA = hk.reindex(pd.MultiIndex.from_arrays([ev[ii], A_])).fillna(0).to_numpy()
        kB = hk.reindex(pd.MultiIndex.from_arrays([ev[ii], B_])).fillna(0).to_numpy()
        ok = ((kA > 0) | (kB > 0)) & ((maxrow[ii] == A_) | (maxrow[ii] == B_))
        ok &= np.abs(upad * float(T.s)) < U_MAX
        ii, ee, upad, kA, kB = ii[ok], ee[ok], upad[ok], kA[ok], kB[ok]
        n_pad = E[['nx', 'ny']].to_numpy()[ee]
        tm3 = np.column_stack([slx[ii], sly[ii]])
        tpad = tm3 @ Rm                                     # displacement M3 -> pad direction
        s_perp = np.einsum('ij,ij->i', tpad, n_pad)
        nm3 = n_pad @ Rm.T                                  # edge normal in the M3 frame
        Bd = dict(u=upad * float(T.s), s=s_perp, nx=nm3[:, 0], ny=nm3[:, 1], X=xyref[ii, 0], Y=xyref[ii, 1],
                  y=(kB > kA).astype(float), i=ii, e=ee, theta=theta[ii], tsplit=tsplit[ii],
                  all=np.zeros(len(ii), int))
        for rt in REG_TYPES[1:]:
            Bd[rt] = rr[rt][ii]
        both = (kA > 0) & (kB > 0)
        # charge fraction (amp only, tot tiebreak dropped)
        mA = pd.MultiIndex.from_arrays([ev[ii], A_[ok]]); mB = pd.MultiIndex.from_arrays([ev[ii], B_[ok]])
        aA = ha.reindex(mA).fillna(0).to_numpy(); aB = ha.reindex(mB).fillna(0).to_numpy()
        Bd['both'] = both
        Bd['dt'] = ht.reindex(mB).to_numpy() - ht.reindex(mA).to_numpy()
        Fq = {k: v[both] for k, v in Bd.items()}
        Fq['f'] = (aB / np.where(aA + aB > 0, aA + aB, 1))[both]
        return Bd, Fq
    B0, F0 = boundary_inputs(xy0)
    L.append(f'boundary samples: {len(B0["u"]):,} (track, edge) pairs with |u| < {U_MAX:g} mm, '
             f'{len(np.unique(B0["i"])):,} tracks; both pads fired in {len(F0["u"]):,}; '
             f's_perp rms {B0["s"].std():.3f}')
    # drift direction from the pulse times of the two pads of a crossing track:
    # for s_perp > 0 the track is on the B side at LARGER z; if B is later, the
    # larger-z ionisation drifts longer, i.e. the cathode (window) is at larger z.
    bb_ = B0['both'] & np.isfinite(B0['dt'])
    ori = []
    for lo, hi in [(-1, -0.2), (-0.2, -0.1), (-0.1, -0.03), (-0.03, 0.03), (0.03, 0.1), (0.1, 0.2), (0.2, 1)]:
        m = bb_ & (B0['s'] >= lo) & (B0['s'] < hi)
        ori.append((lo, hi, int(m.sum()), float(np.median(B0['s'][m])), float(np.median(B0['dt'][m]))))
    ORI = pd.DataFrame(ori, columns=['s_lo', 's_hi', 'n', 's_med', 'dt_med_ns'])
    ORI.to_csv(f'{out}/orientation_timing.csv', index=False, float_format='%.3f')
    m = bb_ & (np.abs(B0['s']) < 0.4) & (np.abs(B0['dt']) < 300)
    pfit = np.polyfit(B0['s'][m], B0['dt'][m], 1)
    sgn_dt = np.sign(ORI[ORI.s_med > 0.1].dt_med_ns.mean() - ORI[ORI.s_med < -0.1].dt_med_ns.mean())
    L.append('drift direction from pad timing (both pads fired; dt = t_B - t_A vs s_perp): ' + ', '.join(
        f's~{r.s_med:+.2f}: {r.dt_med_ns:+.0f} ns (n {r.n})' for r in ORI.itertuples()) +
        f'; linear slope {pfit[0]:+.0f} ns per unit slope')
    win_up = sgn_dt > 0
    L.append('  => ' + ('later signal on the larger-z side: electrons drift towards -z, the mesh is BELOW and '
                        'the drift cathode (window) is ABOVE (+z). So dz > 0 = towards the window.'
                        if win_up else 'later signal on the smaller-z side: the window is BELOW (-z). '
                        'So dz > 0 = towards the mesh/board side.'))
    pg, eg = fit_step(B0)
    fixed = tuple(pg[3:])
    eps = 0.5 / (1 + np.exp(-pg[5]))
    L.append(f'whole-chamber step fit: dz {pg[0]:+.3f} +- {eg[0]:.3f} mm, ox {pg[1]:+.3f}, oy {pg[2]:+.3f} mm, '
             f'w0 {pg[3]:.3f} mm, k {pg[4]:.2f} mm per unit slope, eps {100*eps:.2f} %')
    pf, ef = fit_frac(F0)
    ffixed = tuple(pf[3:])
    L.append(f'whole-chamber charge-fraction fit: dz {pf[0]:+.3f} +- {ef[0]:.3f} mm, w0 {pf[3]:.3f}, '
             f'k {pf[4]:.2f}')
    TB = boundary_table(B0, fixed, F0, ffixed)
    pa_g, ea_g = fit_step(B0, fixed, affine=True)
    L.append(f'whole-chamber step fit with a local affine offset (alignment scale/rotation/shear free): '
             f'dz {pa_g[0]:+.3f} +- {ea_g[0]:.3f} mm; affine terms [mm per 100 mm] ' +
             ' '.join(f'{v:+.3f}+-{e:.3f}' for v, e in zip(pa_g[3:], ea_g[3:])))

    # bootstrap check of the Hessian error (whole chamber and reg8, step fit)
    bt = []
    for nb in range(40):
        idx = RNG.integers(0, len(B0['u']), len(B0['u']))
        bt.append(fit_step(sub(B0, idx), fixed, p0=pg[:3])[0][0])
    L.append(f'bootstrap (40) of the whole-chamber step dz: rms {np.std(bt, ddof=1):.3f} mm vs Hessian '
             f'{TB[TB.rtype == "all"].err_b.iloc[0]:.3f} mm')

    # ================================================================ merge tables
    TBL = TC.merge(TB, on=['rtype', 'region'], how='outer')
    TBL = TBL.merge(TC1[['rtype', 'region', 'dz_c', 'err_c']].rename(columns={'dz_c': 'dz_c1', 'err_c': 'err_c1'}),
                    on=['rtype', 'region'], how='left')
    P_ = TCP[['rtype', 'region', 'dz_c', 'err_c', 'dz_cx', 'err_cx', 'dz_cy', 'err_cy', 'pull_xy',
              'dz_c_pos', 'err_c_pos']].rename(
        columns={'dz_c': 'dz_c_corr', 'err_c': 'err_c_corr', 'dz_cx': 'dz_cx_corr', 'err_cx': 'err_cx_corr',
                 'dz_cy': 'dz_cy_corr', 'err_cy': 'err_cy_corr', 'pull_xy': 'pull_xy_corr'})
    for c_ in ['dz_c_corr', 'err_c_corr', 'dz_cx_corr', 'err_cx_corr', 'dz_cy_corr', 'err_cy_corr',
               'dz_c_pos', 'err_c_pos']:
        P_[c_] = P_[c_] / gain
    TBL = TBL.merge(P_, on=['rtype', 'region'], how='left')
    TBL['toy0_c'] = TBL.dz_c - TBL.dz_c_corr * gain        # paired bias (quantisation + acceptance)
    xm, ym = [], []
    for rt, v in zip(TBL.rtype, TBL.region):
        if rt == 'all':
            xm.append(xc); ym.append(yc)
        else:
            xm.append(cen[rt].loc[v, 'x']); ym.append(cen[rt].loc[v, 'y'])
    TBL['xm'], TBL['ym'] = xm, ym
    TBL['z_eff_c'] = z0 + TBL.dz_c; TBL['z_eff_b'] = z0 + TBL.dz_b
    TBL['pull_cb'] = (TBL.dz_c_corr - TBL.dz_b) / np.hypot(TBL.err_c_corr, TBL.err_b)
    TBL.to_csv(f'{out}/dz_regions.csv', index=False, float_format='%.4f')
    r24 = TBL[TBL.rtype == 'reg24'].drop(columns=['rtype'])
    r24.to_csv(f'{out}/dz_reg24.csv', index=False, float_format='%.4f')

    def g_(r, k):
        v = r.get(k, np.nan)
        return np.nan if v is None else v

    def fmt_tab(rt):
        d = TBL[TBL.rtype == rt]
        s = [f'{"reg":>4s} {"n_c":>6s} {"c_x corr":>13s} {"c_y corr":>13s} {"pxy":>5s} {"c (x&y)":>13s} {"toy0":>6s} '
             f'{"c_corr":>13s} {"c 1-pad":>13s} {"n_b":>6s} {"b step":>13s} {"b 1-fired":>13s} '
             f'{"f frac":>13s} {"pull c-b":>8s} {"c +pos":>13s} {"b affine":>13s}']
        for _, r in d.iterrows():
            s.append(f'{int(r.region):4d} {g_(r, "n_c"):6.0f} {g_(r, "dz_cx_corr"):+5.2f} +- {g_(r, "err_cx_corr"):4.2f} '
                     f'{g_(r, "dz_cy_corr"):+5.2f} +- {g_(r, "err_cy_corr"):4.2f} {g_(r, "pull_xy_corr"):+5.1f} '
                     f'{g_(r, "dz_c"):+5.2f} +- {g_(r, "err_c"):4.2f} {g_(r, "toy0_c"):+6.2f} '
                     f'{g_(r, "dz_c_corr"):+5.2f} +- {g_(r, "err_c_corr"):4.2f} '
                     f'{g_(r, "dz_c1"):+5.2f} +- {g_(r, "err_c1"):4.2f} {g_(r, "n_b"):6.0f} '
                     f'{g_(r, "dz_b"):+5.2f} +- {g_(r, "err_b"):4.2f} {g_(r, "dz_b1"):+5.2f} +- {g_(r, "err_b1"):4.2f} '
                     f'{g_(r, "dz_f"):+5.2f} +- {g_(r, "err_f"):4.2f} {g_(r, "pull_cb"):+8.1f} '
                     f'{g_(r, "dz_c_pos"):+5.2f} +- {g_(r, "err_c_pos"):4.2f} {g_(r, "dz_b_aff"):+5.2f} +- {g_(r, "err_b_aff"):4.2f}')
        return s

    def chi_const(d, col, er):
        k = np.isfinite(d[col]) & np.isfinite(d[er]) & (d[er] > 0)
        w_ = 1 / d[er][k] ** 2; mu = np.average(d[col][k], weights=w_)
        return mu, float((w_ * (d[col][k] - mu) ** 2).sum()), int(k.sum()), float(d[er][k].median())
    for rt in REG_TYPES:
        L.append(f'\n-- per {rt} (dz [mm] w.r.t. z0 = {z0:.2f}; c = centroid (x, y, combined; toy0 = '
                 f'paired toy bias; c_corr = paired-toy-corrected, gain {gain:.3f}; c 1-pad = single-pad clusters '
                 f'only), b = boundary step (all; only one pad fired), f = charge fraction) --')
        L.extend(fmt_tab(rt))
        d = TBL[TBL.rtype == rt]
        if len(d) > 1:
            k1 = np.isfinite(d.pull_xy_corr); k2 = np.isfinite(d.pull_cb)
            cxy_ = float((d.pull_xy_corr[k1] ** 2).sum()); ccb = float((d.pull_cb[k2] ** 2).sum())
            L.append(f'   chi2 x-vs-y (centroid corr) {cxy_:.1f}/{k1.sum()} (p {chi2d.sf(cxy_, k1.sum()):.3f});  '
                     f'chi2 centroid(corr)-vs-boundary {ccb:.1f}/{k2.sum()} (p {chi2d.sf(ccb, k2.sum()):.3f})')
            for nm, col, er in [('centroid raw', 'dz_c', 'err_c'), ('centroid corr', 'dz_c_corr', 'err_c_corr'),
                                ('centroid +pos', 'dz_c_pos', 'err_c_pos'),
                                ('boundary', 'dz_b', 'err_b'), ('boundary affine', 'dz_b_aff', 'err_b_aff'),
                                ('boundary 1-fired', 'dz_b1', 'err_b1'),
                                ('fraction', 'dz_f', 'err_f')]:
                mu, ch, nn_, me = chi_const(d, col, er)
                L.append(f'   {nm:16s}: weighted mean {mu:+.3f}, chi2 vs constant {ch:.1f}/{nn_-1} '
                         f'(p {chi2d.sf(ch, nn_-1):.3g}); median error {me:.2f} mm')

    # ================================================================ global models
    GM = []
    for rt in ('reg48', 'reg24'):
        d = TBL[TBL.rtype == rt]
        for meth, col, er in [('centroid', 'dz_c_corr', 'err_c_corr'), ('cent+pos', 'dz_c_pos', 'err_c_pos'),
                              ('boundary', 'dz_b', 'err_b'), ('b_affine', 'dz_b_aff', 'err_b_aff'),
                              ('bound1', 'dz_b1', 'err_b1'), ('fraction', 'dz_f', 'err_f')]:
            for row in model_fits(d, col, er, xc, yc):
                row.update(rtype=rt, method=meth); GM.append(row)
    GM = pd.DataFrame(GM)
    GM.to_csv(f'{out}/global_models.csv', index=False, float_format='%.5g')
    L.append(f'\n== GLOBAL MODELS on per-region dz (coordinates per 100 mm from the live-pad centroid '
             f'({xc:.1f}, {yc:.1f}) mm M3; gx, gy in mm per 100 mm; bow in mm per (100 mm)^2) ==')
    for _, r in GM.iterrows():
        pars = ' '.join(f'{p}={r[p]:+.3f}+-{r[p+"_err"]:.3f}' for p in
                        ['c0', 'gx', 'gy', 'bow', 'qxx', 'qyy', 'qxy'] if p in r and np.isfinite(r[p]))
        ext = (f'  dchi2 {r.dchi2:.1f}/{r.dndf:.0f} (p {r.p_extra:.3g})' if np.isfinite(r.get('dchi2', np.nan))
               else '')
        L.append(f'{r.rtype:5s} {r.method:8s} {r.model:10s} chi2 {r.chi2:6.1f}/{r.ndf:<3.0f} '
                 f'(p {r.p_value:.3g}){ext}  {pars}')

    # ================================================================ validation
    def summarize_c(C):
        dx, ex, _ = boot_dz(C['rx'], C['tx'], 100); dy, ey, _ = boot_dz(C['ry'], C['ty'], 100)
        w = np.array([1 / ex ** 2, 1 / ey ** 2])
        return (w[0] * dx + w[1] * dy) / w.sum(), 1 / np.sqrt(w.sum()), f'x {dx:+.3f}+-{ex:.3f} y {dy:+.3f}+-{ey:.3f}'

    def summarize_b(Bd):
        p, e = fit_step(Bd, fixed)
        return p[0], e[0]

    def reg24_chi(C=None, Bd=None):
        """chi2 of per-reg24 dz against 0 (null tests)."""
        ch, n = 0.0, 0
        if C is not None:
            for v in np.unique(C['reg24']):
                m = C['reg24'] == v
                for a_ in ('x', 'y'):
                    d_, e_, _ = boot_dz(C['r' + a_][m], C['t' + a_][m], 100)
                    ch += (d_ / e_) ** 2; n += 1
        if Bd is not None:
            for v in np.unique(Bd['reg24']):
                m = Bd['reg24'] == v
                p, e = fit_step(sub(Bd, m), fixed)
                ch += (p[0] / e[0]) ** 2; n += 1
        return ch, n

    # (a) null: slopes shuffled between tracks (same permutation for tx, ty), N_PERM times
    nullc, nullb = [], []
    for ip in range(10):
        perm = RNG.permutation(len(ev))
        txs, tys = tx[perm], ty[perm]
        Cn = cent_inputs(xy0, txs, tys)
        Cnp = {k: v[okp] for k, v in Cn.items()}
        Cnp['rx'] = Cn['rx'][okp] - rt0[okp, 0]; Cnp['ry'] = Cn['ry'][okp] - rt0[okp, 1]
        nullc.append(summarize_c(Cnp)[:2])
        Bn, _ = boundary_inputs(xy0, txs, tys, Pref=P0)
        nullb.append(summarize_b(Bn))
        if ip == 0:
            chc, nc_ = reg24_chi(C=Cnp)
            chb, nb_ = reg24_chi(Bd=Bn)
    for nm, arr, ch, n_ in [('centroid', np.array(nullc), chc, nc_), ('boundary', np.array(nullb), chb, nb_)]:
        V.append(dict(test='null: slopes shuffled, 10 permutations', method=nm, dz=arr[:, 0].mean(),
                      err=arr[:, 1].mean() / np.sqrt(len(arr)), expect=0.0,
                      note=f'per-permutation rms {arr[:, 0].std(ddof=1):.3f} vs fit error {arr[:, 1].mean():.3f}; '
                           f'1st permutation reg24 chi2 vs 0 = {ch:.1f}/{n_} (p {chi2d.sf(ch, n_):.2f})'))
    # (b1) re-projection of the tracks to z0 + 3 mm: dz must move by -3
    x3, y3 = pa.project_to_z(rays, z0 + 3.0); xy3 = np.column_stack([x3, y3])
    rr3, _, _ = region_of(xy3)
    d_, e_, nt = summarize_c(cent_inputs(xy3))
    dref = TBL[(TBL.rtype == 'all')].dz_c.iloc[0]
    V.append(dict(test='reproject tracks to z0+3 (expect dz_ref-3)', method='centroid', dz=d_ - dref,
                  err=e_, expect=-3.0, note=nt))
    B3, _ = boundary_inputs(xy3, rr=rr3)
    d_, e_ = summarize_b(B3)
    dref_b = TBL[(TBL.rtype == 'all')].dz_b.iloc[0]
    V.append(dict(test='reproject tracks to z0+3 (expect dz_ref-3)', method='boundary', dz=d_ - dref_b,
                  err=e_, expect=-3.0, note=f'n {len(B3["u"])}'))
    # (b2) toy injection into the boundary labels: y ~ Bernoulli(P(u + 3 s)) with the fitted shape
    for dzt in (0.0, 3.0):
        w = np.sqrt(fixed[0] ** 2 + (fixed[1] * B0['s']) ** 2)
        Pt = eps + (1 - 2 * eps) * ndtr((B0['u'] + dzt * B0['s']) / w)
        Bt = dict(B0); Bt['y'] = (RNG.random(len(Pt)) < Pt).astype(float)
        d_, e_ = summarize_b(Bt)
        V.append(dict(test=f'boundary toy labels, dz_true {dzt:+g}', method='boundary', dz=d_, err=e_,
                      expect=dzt, note='real (u, s_perp), fitted w0/k/eps'))
    # (c) splits: time, angle, orientation
    for nm, mc, mb in [(f'time 1st half ({tlab})', ~tsplit[selc], ~B0['tsplit']),
                       (f'time 2nd half ({tlab})', tsplit[selc], B0['tsplit'])]:
        d_, e_, nt = summarize_c(sub(Cp, mc[okp]))
        V.append(dict(test=nm, method='centroid', dz=d_, err=e_, expect=np.nan, note=nt))
        d_, e_ = summarize_b(sub(B0, mb))
        V.append(dict(test=nm, method='boundary', dz=d_, err=e_, expect=np.nan,
                      note=f'both-fired fraction {B0["both"][mb].mean():.3f}'))
        d_, e_ = summarize_b(sub(B0, mb & ~B0['both']))
        V.append(dict(test=nm, method='boundary 1-fired', dz=d_, err=e_, expect=np.nan, note=''))
    th_c = theta[selc]
    for lo, hi in [(0, 10), (10, 17), (17, 90)]:
        mc = (th_c >= lo) & (th_c < hi); mb = (B0['theta'] >= lo) & (B0['theta'] < hi)
        d_, e_, nt = summarize_c(sub(Cp, mc[okp]))
        V.append(dict(test=f'theta {lo}-{hi} deg', method='centroid', dz=d_, err=e_, expect=np.nan, note=nt))
        d_, e_ = summarize_b(sub(B0, mb))
        V.append(dict(test=f'theta {lo}-{hi} deg', method='boundary', dz=d_, err=e_, expect=np.nan,
                      note=f'n {mb.sum()}, both-fired fraction {B0["both"][mb].mean():.3f}'))
        d_, e_ = summarize_b(sub(B0, mb & ~B0['both']))
        V.append(dict(test=f'theta {lo}-{hi} deg', method='boundary 1-fired', dz=d_, err=e_, expect=np.nan,
                      note=''))
    for nm, mb in [('edges with |n_x| > |n_y| (M3)', np.abs(B0['nx']) > np.abs(B0['ny'])),
                   ('edges with |n_y| > |n_x| (M3)', np.abs(B0['nx']) <= np.abs(B0['ny']))]:
        d_, e_ = summarize_b(sub(B0, mb))
        V.append(dict(test=nm, method='boundary', dz=d_, err=e_, expect=np.nan, note=f'n {mb.sum()}'))
    # free shape per region type reg8 (check that fixing w0/k/eps does not bias)
    fr = []
    for v in np.unique(B0['reg8']):
        p, e = fit_step(sub(B0, B0['reg8'] == v), None, p0=list(pg))
        fr.append((v, p[0], e[0], p[3], p[4]))
    # (d) is the region-to-region structure real? Independent halves (event parity):
    # per region dz from even and from odd events, same method and across methods.
    evc = ev[selc][okp] % 2 == 0; evb = ev[B0['i']] % 2 == 0
    HALF = []
    for rt in ('reg24', 'reg48'):
        res = {}
        for half, (mc, mb) in {'even': (evc, evb), 'odd': (~evc, ~evb)}.items():
            Tc = centroid_table(sub(Cp, mc), None)
            Tc = Tc[Tc.rtype == rt].set_index('region')
            res[('c', half)] = (Tc.dz_c / gain, Tc.err_c / gain)
            rows_b = []
            for v in np.unique(B0[rt]):
                m = mb & (B0[rt] == v)
                if m.sum() < 150:
                    continue
                p_, e_ = fit_step(sub(B0, m), fixed)
                m1 = m & ~B0['both']
                p1, e1 = fit_step(sub(B0, m1), fixed)
                rows_b.append((v, p_[0], e_[0], p1[0], e1[0]))
            tb = pd.DataFrame(rows_b, columns=['region', 'b', 'eb', 'b1', 'eb1']).set_index('region')
            res[('b', half)] = (tb.b, tb.eb); res[('b1', half)] = (tb.b1, tb.eb1)
        L_h = []
        for a_, b_ in [(('c', 'even'), ('c', 'odd')), (('b', 'even'), ('b', 'odd')), (('b1', 'even'), ('b1', 'odd')),
                       (('c', 'even'), ('b', 'odd')), (('b', 'even'), ('c', 'odd')),
                       (('c', 'even'), ('b1', 'odd')), (('b1', 'even'), ('c', 'odd'))]:
            x1, e1 = res[a_]; x2, e2 = res[b_]
            idx = x1.index.intersection(x2.index)
            x1, e1, x2, e2 = x1[idx], e1[idx], x2[idx], e2[idx]
            k = np.isfinite(x1) & np.isfinite(x2) & np.isfinite(e1) & np.isfinite(e2)
            x1, e1, x2, e2 = x1[k], e1[k], x2[k], e2[k]
            # weighted correlation of the deviations from each set's weighted mean
            d1 = x1 - np.average(x1, weights=1 / e1 ** 2); d2 = x2 - np.average(x2, weights=1 / e2 ** 2)
            r_ = float(np.corrcoef(d1, d2)[0, 1])
            same = a_[0] == b_[0]
            chi = float((((d1 - d2) ** 2) / (e1 ** 2 + e2 ** 2)).sum()) if same else np.nan
            HALF.append(dict(rtype=rt, a=f'{a_[0]}:{a_[1]}', b=f'{b_[0]}:{b_[1]}', n=int(k.sum()), corr=r_,
                             chi2_diff=chi))
            L_h.append(f'{a_[0]}:{a_[1]} vs {b_[0]}:{b_[1]} corr {r_:+.2f}' +
                       (f' (chi2 of the shape difference {chi:.1f}/{int(k.sum())-1})' if same else ''))
        L.append(f'independent halves, {rt} (corr of the region-to-region deviations; 0 if pure noise, '
                 f'+-{1/np.sqrt(len(x1)):.2f} 1-sigma): ' + '; '.join(L_h))
    pd.DataFrame(HALF).to_csv(f'{out}/halves_corr.csv', index=False, float_format='%.4f')
    VAL = pd.DataFrame(V)
    VAL['pull'] = (VAL.dz - VAL.expect) / VAL.err
    VAL.to_csv(f'{out}/validation.csv', index=False, float_format='%.4f')
    L.append('\n== VALIDATION ==')
    for _, r in VAL.iterrows():
        ex = '' if not np.isfinite(r.expect) else f'  expect {r.expect:+.1f}  pull {r.pull:+.1f}'
        L.append(f'{r.method:9s} {r.test:48s} dz {r.dz:+.3f} +- {r.err:.3f}{ex}   {r.note}')
    L.append('reg8 with the step shape free (dz, w0, k): ' + '; '.join(
        f'{int(v)}: {a:+.2f}+-{b:.2f} (w0 {c:.2f}, k {d:.2f})' for v, a, b, c, d in fr))
    t24 = TBL[TBL.rtype == 'reg24']
    L.append(f'\n== ACCEPTANCE (reg24): |corr(tx, x)| max {ACC.corr_tx_x.abs().max():.3f}, '
             f'|corr(ty, y)| max {ACC.corr_ty_y.abs().max():.3f} inside a region; paired-toy bias of the raw '
             f'centroid dz (quantisation x acceptance): whole chamber '
             f'{TBL[TBL.rtype == "all"].toy0_c.iloc[0]:+.3f} mm, reg24 range {t24.toy0_c.min():+.2f} .. '
             f'{t24.toy0_c.max():+.2f} (rms {t24.toy0_c.std():.2f}) mm. The boundary method fits P(B|u, s) '
             f'locally and is not affected by the track density.')
    gcorr = np.corrcoef(tx, x0)[0, 1], np.corrcoef(ty, y0)[0, 1]
    L.append(f'   whole sample corr(tx, x) {gcorr[0]:+.3f}, corr(ty, y) {gcorr[1]:+.3f}; mean slope '
             f'tx {tx.mean():+.4f} ty {ty.mean():+.4f}')
    ACC.merge(t24[['region', 'toy0_c']], on='region').to_csv(f'{out}/acceptance_reg24.csv', index=False,
                                                             float_format='%.4f')

    # ================================================================ 3. GAP STUDY
    # which surface moves? z_eff = mesh + gap/2 -> with an independent per-pad gap
    # (stage 47 drift time), window fixed / board moves: dz = -0.5 per mm of extra
    # gap; board fixed / window moves: +0.5; symmetric: 0.
    gpath = cfg.out_dir('47_deform_pulse_maps') + '/pad_gap_prelim.csv'
    GG = pd.read_csv(gpath).set_index('channel_id')
    gap_row = GG.gap_prelim_mm.reindex(pad_id).to_numpy()
    gcl_row = gap_class(gap_row)
    chk = GG.gap_class.reindex(pad_id).to_numpy()
    assert np.all((gcl_row == chk) | ~np.isfinite(chk))
    eA = E.A.to_numpy()[B0['e']]; eB = E.B.to_numpy()[B0['e']]
    B0['gap'] = 0.5 * (gap_row[eA] + gap_row[eB])
    B0['gcls'] = gap_class(B0['gap'])
    B0['gsame'] = (gcl_row[eA] == gcl_row[eB]) & (gcl_row[eA] >= 0)
    B0['gcls_same'] = np.where(B0['gsame'], gcl_row[eA], -1)
    fA_ = E.A.to_numpy()[F0['e']]; fB_ = E.B.to_numpy()[F0['e']]
    F0['gap'] = 0.5 * (gap_row[fA_] + gap_row[fB_]); F0['gcls'] = gap_class(F0['gap'])
    rowc = crow[jj][okp]
    Cp['gap'] = gap_row[rowc]; Cp['gcls'] = gcl_row[rowc]; Cp['npad'] = npad[okp]
    LG = [f'\n== GAP STUDY (per-pad drift gap from {gpath}; classes <3.5, 3.5-3.9, 3.9-4.3, 4.3-5, 5-6, '
          f'>6 mm; an edge gets the mean gap of its two pads) ==']

    def c_comb(C):
        dx, ex, _ = boot_dz(C['rx'], C['tx']); dy, ey, _ = boot_dz(C['ry'], C['ty'])
        w = np.array([1 / ex ** 2, 1 / ey ** 2])
        return (dx / gain, ex / gain, dy / gain, ey / gain,
                (w[0] * dx + w[1] * dy) / w.sum() / gain, 1 / np.sqrt(w.sum()) / gain)

    rows = []
    for c in range(6):
        mb = B0['gcls'] == c; ms = B0['gcls_same'] == c; mc = Cp['gcls'] == c
        row = dict(gcls=c, n_pads=int((gcl_row == c).sum()), gap_mean_b=float(B0['gap'][mb].mean()),
                   gap_mean_c=float(Cp['gap'][mc].mean()), n_b=int(mb.sum()), n_c=int(mc.sum()),
                   both_frac=float(B0['both'][mb].mean()))
        p_, e_ = fit_step(sub(B0, mb), fixed); row.update(b_all=p_[0], e_all=e_[0])
        p_, e_ = fit_step(sub(B0, mb & ~B0['both']), fixed); row.update(b_1f=p_[0], e_1f=e_[0])
        p_, e_ = fit_step(sub(B0, mb), fixed, affine=True); row.update(b_aff=p_[0], e_aff=e_[0])
        if ms.sum() > 150:
            p_, e_ = fit_step(sub(B0, ms), fixed); row.update(n_bs=int(ms.sum()), b_same=p_[0], e_same=e_[0])
            p_, e_ = fit_step(sub(B0, ms & ~B0['both']), fixed); row.update(b1_same=p_[0], e1_same=e_[0])
        mf = F0['gcls'] == c
        p_, e_ = fit_frac(sub(F0, mf), ffixed); row.update(n_f=int(mf.sum()), f_frac=p_[0], e_f=e_[0])
        dx, ex, dy, ey, dc, ec = c_comb(sub(Cp, mc))
        row.update(c_x=dx, e_cx=ex, c_y=dy, e_cy=ey, c_xy=dc, e_cxy=ec)
        m1 = mc & (Cp['npad'] == 1)
        r1 = c_comb(sub(Cp, m1)); row.update(c_1p=r1[4], e_c1p=r1[5])
        rows.append(row)
    GC = pd.DataFrame(rows)
    GC.to_csv(f'{out}/gap_classes.csv', index=False, float_format='%.4f')
    LG.append('per gap class (dz in mm w.r.t. z0; b = boundary step: all / single-fired / local-affine / both '
              'pads in the class; c = paired-toy-corrected centroid x / y / combined / single-pad clusters):')
    LG.append(f'{"cls":>3s} {"pads":>4s} {"<gap>":>5s} {"n_b":>5s} {"both":>5s} {"b all":>12s} {"b 1-fired":>12s} '
              f'{"b affine":>12s} {"b same":>12s} {"b1 same":>12s} {"n_c":>5s} {"c x":>12s} {"c y":>12s} '
              f'{"c x&y":>12s} {"c 1-pad":>12s} {"frac":>12s}')
    for r in GC.itertuples():
        f2 = lambda v, e: f'{v:+5.2f}+-{e:4.2f}' if np.isfinite(v) else f'{"n/a":>12s}'
        LG.append(f'{r.gcls:3d} {r.n_pads:4d} {r.gap_mean_b:5.2f} {r.n_b:5d} {r.both_frac:5.3f} '
                  f'{f2(r.b_all, r.e_all)} {f2(r.b_1f, r.e_1f)} {f2(r.b_aff, r.e_aff)} '
                  f'{f2(getattr(r, "b_same", np.nan), getattr(r, "e_same", np.nan))} '
                  f'{f2(getattr(r, "b1_same", np.nan), getattr(r, "e1_same", np.nan))} {r.n_c:5d} '
                  f'{f2(r.c_x, r.e_cx)} {f2(r.c_y, r.e_cy)} {f2(r.c_xy, r.e_cxy)} {f2(r.c_1p, r.e_c1p)} {f2(r.f_frac, r.e_f)}')
    # class-level chi2 for b fixed (a free)
    LG.append('class-level fit dz_class = a + b (<gap> - 4), chi2 (ndf 5) for b = -0.5 / 0 / +0.5 and free b:')
    for col, er in [('b_all', 'e_all'), ('b_1f', 'e_1f'), ('b_aff', 'e_aff'), ('b_same', 'e_same'),
                    ('c_xy', 'e_cxy'), ('c_x', 'e_cx'), ('c_y', 'e_cy'), ('c_1p', 'e_c1p'), ('f_frac', 'e_f')]:
        d = GC[np.isfinite(GC[col])]
        gx_ = (d.gap_mean_c if col.startswith('c') else d.gap_mean_b).to_numpy() - 4
        z_ = d[col].to_numpy(); w_ = 1 / d[er].to_numpy() ** 2
        out_ = []
        for b0 in (-0.5, 0.0, 0.5):
            a_ = np.average(z_ - b0 * gx_, weights=w_)
            out_.append(float((w_ * (z_ - a_ - b0 * gx_) ** 2).sum()))
        A_ = np.column_stack([np.ones(len(z_)), gx_])
        cv = np.linalg.inv(A_.T @ (A_ * w_[:, None])); pp = cv @ (A_.T @ (w_ * z_))
        chf = float((w_ * (z_ - A_ @ pp) ** 2).sum())
        LG.append(f'  {col:7s} chi2 {out_[0]:6.1f} / {out_[1]:6.1f} / {out_[2]:6.1f};  free b {pp[1]:+.3f} +- '
                  f'{np.sqrt(cv[1, 1]):.3f} (chi2 {chf:.1f}/{len(z_) - 2})')

    # pad-level fits: dz = a + b (gap - 4) [+ gx X + gy Y]
    def Zb(B, tilt):
        k = [np.ones(len(B['u'])), B['gap'] - 4.0]
        if tilt:
            k += [(B['X'] - xc) / 100.0, (B['Y'] - yc) / 100.0]
        return np.column_stack(k)

    def Zc_f(tilt):
        def f(C):
            k = [np.ones(len(C['rx'])), C['gap'] - 4.0]
            if tilt:
                k += [(C['X'] - xc) / 100.0, (C['Y'] - yc) / 100.0]
            return np.column_stack(k)
        return f
    okb = np.isfinite(B0['gap']); okc = np.isfinite(Cp['gap']); okf = np.isfinite(F0['gap'])
    apex_f = np.isin(F0['reg24'], [0, 1])
    apex_b = np.isin(B0['reg24'], [0, 1]); apex_c = np.isin(Cp['reg24'], [0, 1])
    FITS = []
    for tilt, excl in [(False, False), (True, False), (True, True)]:
        tag = ('tilt' if tilt else 'no tilt') + (', excl. reg24 0,1' if excl else '')
        mb = okb & (~apex_b if excl else True)
        mc = okc & (~apex_c if excl else True)
        ests = [('boundary all', 'b', mb), ('boundary 1-fired', 'b', mb & ~B0['both']),
                ('boundary all, same-class edges', 'b', mb & B0['gsame']),
                ('centroid x&y', 'c', ('x', 'y'), mc), ('centroid x', 'c', ('x',), mc),
                ('centroid y', 'c', ('y',), mc), ('centroid single-pad', 'c', ('x', 'y'), mc & (Cp['npad'] == 1)),
                ('charge fraction', 'f', okf & (~apex_f if excl else True))]
        for est in ests:
            if est[1] == 'b':
                Bs = sub(B0, est[2])
                beta, cov = fit_lin_b(Bs, Zb(Bs, tilt), fixed)
            elif est[1] == 'f':
                Fs = sub(F0, est[2])
                beta, cov = fit_lin_f(Fs, Zb(Fs, tilt), ffixed)
            else:
                beta, cov = fit_lin_c(sub(Cp, est[3]), Zc_f(tilt), axes=est[2])
                beta = beta / gain; cov = cov / gain ** 2
            eb = np.sqrt(cov[1, 1])
            row = dict(fit=tag, est=est[0], a=beta[0], ea=np.sqrt(cov[0, 0]), b=beta[1], eb=eb,
                       dchi2_m05=((beta[1] + 0.5) / eb) ** 2, dchi2_0=(beta[1] / eb) ** 2,
                       dchi2_p05=((beta[1] - 0.5) / eb) ** 2)
            if tilt:
                row.update(gx=beta[2], egx=np.sqrt(cov[2, 2]), gy=beta[3], egy=np.sqrt(cov[3, 3]),
                           corr_b_gx=cov[1, 2] / np.sqrt(cov[1, 1] * cov[2, 2]),
                           corr_b_gy=cov[1, 3] / np.sqrt(cov[1, 1] * cov[3, 3]))
            FITS.append(row)
    FT = pd.DataFrame(FITS)
    FT.to_csv(f'{out}/gap_fits.csv', index=False, float_format='%.4f')
    LG.append('\npad-level fit dz = a + b (gap - 4) [+ gx (X-xc)/100 + gy (Y-yc)/100]; dchi2 of b fixed at '
              '-0.5 / 0 / +0.5 relative to the best fit:')
    for r in FT.itertuples():
        t_ = (f'  gx {r.gx:+.3f}+-{r.egx:.3f} gy {r.gy:+.3f}+-{r.egy:.3f} corr(b,gx) {r.corr_b_gx:+.2f} '
              f'corr(b,gy) {r.corr_b_gy:+.2f}' if np.isfinite(getattr(r, 'gx', np.nan)) else '')
        LG.append(f'  [{r.fit:24s}] {r.est:32s} a {r.a:+.3f}+-{r.ea:.3f}  b {r.b:+.3f} +- {r.eb:.3f}  '
                  f'dchi2 {r.dchi2_m05:6.1f} / {r.dchi2_0:5.1f} / {r.dchi2_p05:6.1f}{t_}')
    # spatial layout of the gap (degeneracy with a tilt)
    LG.append(f'gap vs position (pads): corr(gap, X_M3) {np.corrcoef(gap_row[np.isfinite(gap_row)], padm3[np.isfinite(gap_row), 0])[0, 1]:+.2f}, '
              f'corr(gap, Y_M3) {np.corrcoef(gap_row[np.isfinite(gap_row)], padm3[np.isfinite(gap_row), 1])[0, 1]:+.2f}')

    # ---------------- toy: what each hypothesis predicts for each estimator ----------------
    rng_t = np.random.default_rng(4711)
    both_data = float(B0['both'].mean())
    b_src = {k_: B0[k_] for k_ in ('u', 's', 'theta', 'nx', 'ny', 'X', 'Y')}

    def calib_thr(Lg, fire):
        lo, hi = 0.2, 2000.0
        for _ in range(16):
            mid = np.sqrt(lo * hi)
            Bt, _ = toy_boundary(b_src, 4.0, -2.0, Lg, mid, np.random.default_rng(99), n_rep=1, fire=fire)
            if Bt['both'].mean() > both_data:
                lo = mid
            else:
                hi = mid
        return np.sqrt(lo * hi)

    def toy_cent(G, zm, Lg, thr, rng, fire, sig_p=0.7, sig_d=0.45, chunk=15000):
        """2D toy on the real selected tracks and tiles: clusters along the track,
        amplitudes per pad, fired pads (> thr), amplitude-weighted centroid,
        paired with the same dz = 0 single-tile residual as the data."""
        sel_i = np.nonzero(selc)[0][okp]
        res = np.full((len(sel_i), 2), np.nan); npd = np.zeros(len(sel_i), int)
        for c0 in range(0, len(sel_i), chunk):
            ii = sel_i[c0:c0 + chunk]; n = len(ii)
            ct = np.cos(np.radians(theta[ii]))
            hf, q, valid = toy_clusters(n, G, ct, rng)
            h = zm + G * hf
            pt = (xy0[ii, None, :] + np.stack([tx[ii], ty[ii]], -1)[:, None, :] * h[..., None]
                  + rng.normal(0, sig_p, (n, 1, 2)) + rng.normal(0, sig_d, (n, hf.shape[1], 2)))
            a = (q * amp_weight(G * hf, Lg))[valid]; qv = q[valid]
            ti = np.nonzero(valid)[0]
            Pp = inverse(T, pt[valid])
            _, cand = ptree.query(Pp, k=4)
            rw = np.full(len(Pp), -1)
            for kk in range(4):
                ins = (pqd(Pp, tv[cand[:, kk]]) == 0) & (rw < 0)
                rw[ins] = cand[ins, kk]
            m = rw >= 0
            key = ti[m].astype(np.int64) * 2000 + rw[m]
            uk, inv_ = np.unique(key, return_inverse=True)
            amp = np.bincount(inv_, weights=a[m])
            f = (np.bincount(inv_, weights=qv[m]) if fire == 'charge' else amp) > thr
            f &= amp > 0
            tr = (uk // 2000)[f]; pr = uk[f] % 2000; am = amp[f]
            W = np.bincount(tr, weights=am, minlength=n)
            cx_ = np.bincount(tr, weights=am * padm3[pr, 0], minlength=n)
            cy_ = np.bincount(tr, weights=am * padm3[pr, 1], minlength=n)
            ok_ = W > 0
            res[c0:c0 + n][ok_] = np.column_stack([cx_[ok_] / W[ok_], cy_[ok_] / W[ok_]]) - xy0[ii][ok_]
            npd[c0:c0 + n] = np.bincount(tr, minlength=n)
        r = res - rt0[okp]
        k = np.isfinite(r).all(1)
        out_ = {}
        for nm, mm in [('c_xy', k), ('c_1p', k & (npd == 1))]:
            dx, ex, _ = boot_dz(r[mm, 0], tx[sel_i][mm], 40); dy, ey, _ = boot_dz(r[mm, 1], ty[sel_i][mm], 40)
            w = np.array([1 / ex ** 2, 1 / ey ** 2])
            out_[nm] = ((w[0] * dx + w[1] * dy) / w.sum() / gain, 1 / np.sqrt(w.sum()) / gain)
        out_['multi_frac'] = float((npd[k] > 1).mean())
        return out_

    TOY = []
    HYP = [('base: gap 4, mesh -2', 4.0, -2.0), ('H_window: gap 7, mesh -2 (window moved up)', 7.0, -2.0),
           ('H_board: gap 7, mesh -5 (board moved down)', 7.0, -5.0), ('symmetric: gap 7, mesh -3.5', 7.0, -3.5)]
    MODELS = [('A charge', None, 'amp', True), ('B gate 1.5 (fire+amp)', 1.5, 'amp', True),
              ('C hybrid 2.0', 2.0, 'charge', False), ('C hybrid 1.5', 1.5, 'charge', True),
              ('C hybrid 1.0', 1.0, 'charge', False), ('C hybrid 0.5', 0.5, 'charge', False)]
    for mname, Lg, fire, do_c in MODELS:
        thr = calib_thr(Lg, fire)
        Bb, Fb = toy_boundary(b_src, 4.0, -2.0, Lg, thr, rng_t, fire=fire)
        pfree, _ = fit_step(Bb)
        fx = tuple(pfree[3:])
        pff, _ = fit_frac(Fb)
        ffx = tuple(pff[3:])
        for name, G, zm in HYP:
            Bt, Ft = toy_boundary(b_src, G, zm, Lg, thr, rng_t, fire=fire)
            r_ = dict(Lg=mname, thr=thr, hyp=name, G=G, zm=zm, geo_mid=zm + G / 2,
                      both_frac=float(Bt['both'].mean()))
            p_, e_ = fit_step(Bt, fx); r_.update(b_all=p_[0], e_all=e_[0])
            p_, e_ = fit_step(sub(Bt, ~Bt['both']), fx); r_.update(b_1f=p_[0], e_1f=e_[0])
            p_, e_ = fit_frac(Ft, ffx); r_.update(f_frac=p_[0], e_f=e_[0])
            if do_c:
                cc = toy_cent(G, zm, Lg, thr, rng_t, fire)
                r_.update(c_xy=cc['c_xy'][0], e_c=cc['c_xy'][1], c_1p=cc['c_1p'][0], e_c1p=cc['c_1p'][1],
                          multi_frac=cc['multi_frac'])
            TOY.append(r_)
    TY = pd.DataFrame(TOY)
    TY.to_csv(f'{out}/gap_toy.csv', index=False, float_format='%.4f')
    data4 = dict(b_all=float(TB[TB.rtype == 'all'].dz_b.iloc[0]), b_1f=float(TB[TB.rtype == 'all'].dz_b1.iloc[0]),
                 f_frac=float(pf[0]), both=both_data)
    LG.append('models: A = amplitude = charge; B = amplitude AND firing use only clusters within Lg of the mesh; '
              'C = a pad fires on its total charge, but the amplitude (max-pad choice, fraction, centroid weight) '
              'counts only clusters within Lg of the mesh.')
    LG.append(f'\ntoy (real boundary samples x3 / real centroid tracks; Ar-like clusters 2.5/mm, 1/k^2 sizes; threshold '
              f'tuned to the data both-fired fraction {both_data:.3f} at gap 4; amplitude gate Lg = drift distance '
              f'beyond which clusters stop counting in the amplitude). Data at the reference: b_all {data4["b_all"]:+.2f}, '
              f'b_1f {data4["b_1f"]:+.2f}, frac {data4["f_frac"]:+.2f}; frac - b_1f = {data4["f_frac"] - data4["b_1f"]:+.2f}.')
    for Lg_, d in TY.groupby('Lg', sort=False):
        b0 = d.iloc[0]
        LG.append(f' {Lg_}: base (geometric mid 0, both-fired {b0.both_frac:.3f}) b_all {b0.b_all:+.2f} b_1f '
                  f'{b0.b_1f:+.2f} frac {b0.f_frac:+.2f} (b_all - b_1f {b0.b_all - b0.b_1f:+.2f}, frac - b_1f '
                  f'{b0.f_frac - b0.b_1f:+.2f})' +
                  (f' c_xy {b0.c_xy:+.2f} c_1p {b0.c_1p:+.2f}' if np.isfinite(b0.get('c_xy', np.nan)) else ''))
        for r in d.iloc[1:].itertuples():
            pred = lambda col: (getattr(r, col) - getattr(b0, col)) / 3.0
            cpart = (f' c_xy {pred("c_xy"):+.2f} c_1p {pred("c_1p"):+.2f}'
                     if np.isfinite(getattr(r, 'c_xy', np.nan)) else '')
            LG.append(f'    {r.hyp:44s} predicted b = d(dz)/d(gap): geometric {(r.geo_mid - b0.geo_mid)/3:+.2f} | '
                      f'b_all {pred("b_all"):+.2f} b_1f {pred("b_1f"):+.2f} frac {pred("f_frac"):+.2f}{cpart} '
                      f'(both-fired {b0.both_frac:.3f} -> {r.both_frac:.3f}; toy stat ~0.03)')
    LG.append(f'data both-fired fraction per gap class: ' +
              ', '.join(f'{r.gap_mean_b:.2f} mm: {r.both_frac:.3f}' for r in GC.itertuples()))
    L.extend(LG)

    # ================================================================ figures
    # (1) residual vs slope, whole chamber
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    sb = np.linspace(-0.45, 0.45, 19)
    for ax, a_ in zip(axs, ('x', 'y')):
        r, t = C0['r' + a_], C0['t' + a_]
        k = np.abs(r) < TRIM
        c_ = np.digitize(t[k], sb)
        mu = [r[k][c_ == i].mean() if (c_ == i).sum() > 50 else np.nan for i in range(1, len(sb))]
        er = [r[k][c_ == i].std() / np.sqrt((c_ == i).sum()) if (c_ == i).sum() > 50 else np.nan
              for i in range(1, len(sb))]
        xs = 0.5 * (sb[1:] + sb[:-1])
        ax.errorbar(xs, mu, er, fmt='o', color=C2, label='data (mean, |r| < 15 mm)')
        dz_, a0 = ols_dz(r, t)
        ax.plot(xs, a0 + dz_ * xs, color=C3, lw=2, label=f'fit: dz = {dz_:+.2f} mm')
        ax.set(xlabel=f'track slope d{a_}/dz', ylabel=f'cluster - track at z0 in {a_} [mm]',
               title=f'Centroid residual vs slope ({a_}), whole chamber')
        ax.grid(color=GRID); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(f'{out}/residual_vs_slope.png', dpi=140, facecolor='white')
    plt.close(fig)

    # (2) boundary step in slope classes
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ub = np.linspace(-U_MAX, U_MAX, 33); xs = 0.5 * (ub[1:] + ub[:-1])
    for (lo, hi), col in zip([(-1, -0.15), (-0.05, 0.05), (0.15, 1)], (C2, C1, C3)):
        m = (B0['s'] >= lo) & (B0['s'] < hi)
        uu = B0['u'][m] - pg[1] * B0['nx'][m] - pg[2] * B0['ny'][m]
        c_ = np.digitize(uu, ub)
        pr = [B0['y'][m][c_ == i].mean() if (c_ == i).sum() > 20 else np.nan for i in range(1, len(ub))]
        ax.plot(xs, pr, 'o', color=col, ms=4, label=f's_perp in [{lo:+.2f}, {hi:+.2f}) (n {m.sum()})')
        sm = np.median(B0['s'][m]); w = np.sqrt(fixed[0] ** 2 + (fixed[1] * sm) ** 2)
        ax.plot(xs, eps + (1 - 2 * eps) * ndtr((xs + pg[0] * sm) / w), color=col, lw=1.5)
    ax.set(xlabel='signed distance u of the track at z0 to the shared edge [mm] (+ into pad B)',
           ylabel='P(pad B has the larger amplitude)',
           title=f'Boundary step by slope class (dz = {pg[0]:+.2f} mm)')
    ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f'{out}/boundary_step.png', dpi=140, facecolor='white')
    plt.close(fig)

    # (3) reg24 maps: centroid, boundary, difference (pull)
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.4))
    tm3 = cm.apply_T(T, tv)
    lab24 = regs_of_row['reg24']
    r24i = r24.set_index('region')
    vmax = np.nanmax(np.abs(np.r_[r24.dz_c_corr, r24.dz_b])); vmax = max(vmax, 0.5)
    for ax, col, ttl, cmap, vm in [(axs[0], 'dz_c_corr', 'centroid dz (toy-corrected) [mm]', 'RdBu_r', vmax),
                                   (axs[1], 'dz_b', 'boundary dz [mm]', 'RdBu_r', vmax),
                                   (axs[2], 'pull_cb', '(centroid - boundary)/sigma', 'PuOr', 3)]:
        vals = r24i[col].reindex(lab24).to_numpy()
        pc = PolyCollection(tm3, array=vals, cmap=cmap, edgecolors='none')
        pc.set_clim(-vm, vm); ax.add_collection(pc)
        for v, rr_ in r24i.iterrows():
            e_ = {'dz_c_corr': rr_.err_c_corr, 'dz_b': rr_.err_b}.get(col, None)
            txt = f'{rr_[col]:+.2f}' + (f'\n+-{e_:.2f}' if e_ is not None else '')
            ax.text(cen['reg24'].loc[v, 'x'], cen['reg24'].loc[v, 'y'], txt, ha='center', va='center',
                    fontsize=6.5)
        plt.colorbar(pc, ax=ax, label=ttl)
        ax.autoscale(); ax.set_aspect('equal')
        ax.set(xlabel='x M3 [mm]', ylabel='y M3 [mm]', title=f'reg24: {ttl} (z0 = {z0:.2f} mm)')
    fig.tight_layout(); fig.savefig(f'{out}/map_reg24.png', dpi=140, facecolor='white'); plt.close(fig)

    # (4) methods compare
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    for rt, col in [('reg24', C2), ('reg48', C1), ('reg8', C3)]:
        d = TBL[TBL.rtype == rt]
        axs[0].errorbar(d.dz_b, d.dz_c_corr, xerr=d.err_b, yerr=d.err_c_corr, fmt='o', color=col, ms=4, label=rt)
        axs[1].errorbar(d.dz_cx, d.dz_cy, xerr=d.err_cx, yerr=d.err_cy, fmt='o', color=col, ms=4, label=rt)
    for ax, (xl, yl, tt) in zip(axs, [('boundary dz [mm]', 'centroid dz, toy-corrected [mm]', 'centroid vs boundary'),
                                      ('centroid dz from x [mm]', 'centroid dz from y [mm]',
                                       'centroid: x vs y')]):
        lim = ax.get_xlim()
        ax.plot([-10, 10], [-10, 10], color=GRID, lw=1); ax.set_xlim(lim)
        ax.set(xlabel=xl, ylabel=yl, title=tt); ax.grid(color=GRID); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(f'{out}/methods_compare.png', dpi=140, facecolor='white')
    plt.close(fig)

    L.append(f'\nassumptions: trimmed stage-34 alignment (R, s, pad_mean, m3_mean, z0 = {z0:.3f}); '
             f'pad tiles from the 34 cache geometry; regions from 43_deformation/regions.csv; '
             f'centroid: CC20 cluster < 20 mm, region/selection by the cluster max pad (no track cut), paired-toy '
             f'correction (gain {gain:.3f}), trim {TRIM:g} mm, bootstrap {NBOOT}; '
             f'boundary: |u| < {U_MAX:g} mm, {END_MARGIN:g} mm from edge ends, max pad within 15 mm '
             f'must be A or B, (w0, k, eps) fixed from the chamber fit.')
    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
