#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pillar_stats.py -- offline analysis of the p2_pillars.py products.

Two independent things, sharing one pass over the data:

**Dead-area masking.**  Build a mask per station out of things that are known
independently of the efficiency being measured -- the gerber's medium and big
pillars, channels that never responded, and connector blocks that are weak as a
*block* -- then apply all three stations' masks to the stack at once and remeasure.
The masked numbers say how the stack performs where it is instrumented and
alive; the unmasked ones say what this particular stack delivered.  Both are
real, they answer different questions, and the gap between them is the part a
production detector would not inherit.

**The bulk pillar lattice.**  The 0.5 mm pillars sit on an exact 2.000 mm square
lattice at even-integer millimetres of the board frame, which is the frame the
fine map is binned in.  So the expected efficiency modulation is at
k = (pi, 0) and (0, pi) rad/mm with **phase zero and a negative sign** -- there
is nothing to fit, and a signal at the right |k| with the wrong phase is not a
detection.  The estimator is a periodogram of the map's high-pass residual; its
null is calibrated empirically from k well off the lattice, so no assumption
about per-track variance enters the significance.

The ratio of the harmonics measures the pointing resolution without knowing how
dead a pillar is: at fixed dead fraction the m-th harmonic carries
exp(-m^2 k0^2 sigma^2 / 2), so A2/A1 gives sigma and A1 then gives the dead
fraction.  That is the quantitative answer to "can the reference resolve them".
"""
import json
import os

import numpy as np

import pillar_geom as G

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
STATIONS = ("P2_IN", "P2_MID", "P2_OUT")
Z = {"P2_IN": 320.0, "P2_MID": 630.0, "P2_OUT": 940.0}
RUN = "eff_nominal_1"
K0 = 2 * np.pi / G.PITCH             # 3.14159 rad/mm, the lattice fundamental


def load(run=RUN):
    z = np.load(os.path.join(DATA, f"p2_pillars_{run}.npz"), allow_pickle=True)
    j = json.load(open(os.path.join(DATA, f"p2_pillars_{run}.json")))
    return z, j


def sample(z):
    cols = [str(c) for c in z["sample_cols"]]
    S = z["sample"].astype(np.float64)
    return {k: S[:, i] for i, k in enumerate(cols)}


def fine(z, s):
    """The 0.20 mm map of one station, with its bin centres."""
    b, half = float(z["fine_bin"][0]), float(z["fine_half"][0])
    cx, cy = z[f"centre_{s}"]
    n = z[f"fine_n_{s}"].astype(float)
    xc = cx - half + b * (0.5 + np.arange(n.shape[0]))
    yc = cy - half + b * (0.5 + np.arange(n.shape[1]))
    return dict(b=b, xc=xc, yc=yc, n=n,
                k=z[f"fine_k_{s}"].astype(float),
                asum=z[f"fine_asum_{s}"].astype(float),
                an=z[f"fine_an_{s}"].astype(float),
                csum=z[f"fine_csum_{s}"].astype(float))


# --------------------------------------------------------------------------- #
# the periodogram
# --------------------------------------------------------------------------- #
def _gauss_blur(a, sigma_bins):
    """Separable Gaussian blur, done in the Fourier domain so it is exact and
    wraps consistently for both the numerator and the denominator."""
    n0, n1 = a.shape
    f0 = np.exp(-0.5 * (2 * np.pi * np.fft.fftfreq(n0) * sigma_bins) ** 2)
    f1 = np.exp(-0.5 * (2 * np.pi * np.fft.fftfreq(n1) * sigma_bins) ** 2)
    return np.real(np.fft.ifft2(np.fft.fft2(a) * f0[:, None] * f1[None, :]))


def residual_map(f, what="eff", smooth_mm=4.0):
    """High-pass residual of a ratio map, and the weight that normalises it.

    The efficiency varies by several per cent across the beam spot (gain
    gradient) and that is a k -> 0 structure; subtracting a Gaussian-smoothed
    local mean removes it.  sigma = 4 mm suppresses k = pi by exp(-79), i.e. it
    takes nothing at all out of the signal.
    """
    if what == "eff":
        num, den = f["k"], f["n"]
    elif what == "amp":
        num, den = f["asum"], f["an"]
    else:
        num, den = f["csum"], f["an"]
    sb = smooth_mm / f["b"]
    sn, sd = _gauss_blur(num, sb), _gauss_blur(den, sb)
    with np.errstate(invalid="ignore", divide="ignore"):
        loc = np.where(sd > 0, sn / np.maximum(sd, 1e-9), 0.0)
    w = num - loc * den
    live = den > 0
    return w * live, loc * den * live, live


def dft(w, xc, yc, kx, ky, bin_mm):
    """S(k) = sum_b w_b exp(i k.x_b), separable, with the bin-width correction.

    A histogram is the true density convolved with the bin, so every amplitude
    read off one is low by sinc(k b / 2) per axis.  The correction is exact.
    """
    kx = np.atleast_1d(np.asarray(kx, float))
    ky = np.atleast_1d(np.asarray(ky, float))
    Ex = np.exp(1j * kx[:, None] * xc[None, :])        # (Nkx, Nx)
    Ey = np.exp(1j * ky[:, None] * yc[None, :])        # (Nky, Ny)
    S = (Ex @ w) @ Ey.T                                # (Nkx, Nky)
    h = bin_mm / 2
    cx = np.sinc(kx * h / np.pi)[:, None]
    cy = np.sinc(ky * h / np.pi)[None, :]
    return S / (cx * cy)


def modulation(f, kx, ky, what="eff", smooth_mm=4.0):
    """Complex modulation amplitude A e^{i phi} at each k.

    ratio(x) = local_mean(x) * (1 + Re[A e^{i(k.x + phi)}])  ->  A = 2|S| / K,
    with K the summed numerator.  Sign convention: a dip *on* the lattice node
    gives phase pi at k = (pi, 0), because the nodes are at even integers.
    """
    w, K, _ = residual_map(f, what, smooth_mm)
    S = dft(w, f["xc"], f["yc"], kx, ky, f["b"])
    return 2 * S / K.sum(), K.sum()


def null_scale(f, what="eff", smooth_mm=4.0, nk=24, seed=0):
    """Empirical null: |A|^2 averaged over k well away from any lattice line.

    Calibrating on the data itself rather than on a variance model means the
    significance does not depend on knowing the per-track amplitude spread, and
    it absorbs whatever the smoothing and the finite window do.
    """
    rng = np.random.default_rng(seed)
    ka = []
    while len(ka) < nk * nk:
        v = rng.uniform(-1.6 * K0, 1.6 * K0, 2)
        # keep away from both lattice lines and from k ~ 0
        if np.hypot(*v) < 0.6:
            continue
        if min(abs(abs(v[0]) - K0), abs(abs(v[1]) - K0)) < 0.35:
            continue
        if min(abs(v[0]), abs(v[1])) < 0.35 and np.hypot(*v) > 0.6:
            pass
        ka.append(v)
    ka = np.array(ka)
    A = np.array([modulation(f, [a[0]], [a[1]], what, smooth_mm)[0][0, 0]
                  for a in ka])
    # |A|^2 ~ Exp(mean), so the mean is the scale; sigma per component is
    # sqrt(mean/2) and A/sigma_A with sigma_A = sqrt(mean/2) is the 1-D width.
    return float(np.mean(np.abs(A) ** 2)), ka, A


def refine(f, k0, what="eff", span=0.06, n=61, smooth_mm=4.0):
    """Peak of |A| in a small box around k0, and the complex value there."""
    kx = k0[0] + np.linspace(-span, span, n)
    ky = k0[1] + np.linspace(-span, span, n)
    A, _ = modulation(f, kx, ky, what, smooth_mm)
    i = np.unravel_index(np.argmax(np.abs(A)), A.shape)
    return np.array([kx[i[0]], ky[i[1]]]), complex(A[i])


# the three first-order reciprocal vectors of a triangular lattice with rows
# along x, as measured: 30, 90 and -30 degrees, |G| = 4 pi / (sqrt3 d)
def hex_seeds(d=3.0):
    g = 4 * np.pi / (np.sqrt(3) * d)
    return [np.array([g * np.cos(np.radians(a)), g * np.sin(np.radians(a))])
            for a in (90.0, 30.0, -30.0)]


def dft_pairs(w, xc, yc, K, bin_mm, chunk=400):
    """S(k) for an arbitrary LIST of k vectors (not an outer-product grid)."""
    K = np.atleast_2d(np.asarray(K, float))
    out = np.empty(len(K), complex)
    for i in range(0, len(K), chunk):
        kk = K[i:i + chunk]
        Ex = np.exp(1j * kk[:, :1] * xc[None, :])
        Ey = np.exp(1j * kk[:, 1:2] * yc[None, :])
        out[i:i + chunk] = np.einsum("ia,ab,ib->i", Ex, w, Ey)
    h = bin_mm / 2
    cx = np.sinc(K[:, 0] * h / np.pi)
    cy = np.sinc(K[:, 1] * h / np.pi)
    return out / (cx * cy)


def _trim(w, xc, yc):
    """Drop the empty border of the map -- the beam lights maybe a third of the
    window and the rest is exactly zero, so carrying it costs time and nothing
    else."""
    ix = np.where(np.abs(w).sum(1) > 0)[0]
    iy = np.where(np.abs(w).sum(0) > 0)[0]
    if not len(ix) or not len(iy):
        return w, xc, yc
    sx = slice(ix[0], ix[-1] + 1)
    sy = slice(iy[0], iy[-1] + 1)
    return w[sx, sy], xc[sx], yc[sy]


def _hex_power(w, xc, yc, b, dv, tv, order):
    D, T = np.meshgrid(dv, np.radians(tv), indexing="ij")
    g = order * 4 * np.pi / (np.sqrt(3) * D)
    ks = np.stack([np.stack([g * np.cos(np.radians(a) + T),
                             g * np.sin(np.radians(a) + T)], -1)
                   for a in (90.0, 30.0, -30.0)])
    S = dft_pairs(w, xc, yc, ks.reshape(-1, 2), b)
    return ks, S.reshape(ks.shape[:-1])


def hex_scan(f, what="eff", smooth_mm=4.0, order=1,
             d0=(2.85, 3.15), th0=(-2.0, 2.0), n0=25):
    """Scan the lattice itself -- spacing and orientation -- not three
    independent wavevectors.

    Fitting the three first-order vectors separately lets them disagree by a
    few parts per thousand in |G|, and over a 64 mm window that is a radian of
    phase, which then puts the folded origin on the wrong site.  Scanning
    (d, theta) forces one lattice and the summed power over the three
    directions is the statistic.  Coarse pass, then a fine one around the peak.
    """
    w, K, _ = residual_map(f, what, smooth_mm)
    w, xc, yc = _trim(w, f["xc"], f["yc"])
    Ktot = K.sum()
    dv = np.linspace(*d0, n0)
    tv = np.linspace(*th0, n0)
    ks, S = _hex_power(w, xc, yc, f["b"], dv, tv, order)
    P = (np.abs(S) ** 2).sum(0)
    i = np.unravel_index(np.argmax(P), P.shape)
    hd = (dv[1] - dv[0]) * 1.5
    ht = (tv[1] - tv[0]) * 1.5
    dv2 = np.linspace(dv[i[0]] - hd, dv[i[0]] + hd, 31)
    tv2 = np.linspace(tv[i[1]] - ht, tv[i[1]] + ht, 31)
    ks2, S2 = _hex_power(w, xc, yc, f["b"], dv2, tv2, order)
    P2 = (np.abs(S2) ** 2).sum(0)
    j = np.unravel_index(np.argmax(P2), P2.shape)
    return dict(d=float(dv2[j[0]]), theta_deg=float(tv2[j[1]]),
                G=ks2[:, j[0], j[1], :], A=2 * S2[:, j[0], j[1]] / Ktot,
                power=(np.abs(2 * S2 / Ktot) ** 2).sum(0),
                d_grid=dv2, theta_grid=tv2)


def hex_at(f, d, theta_deg, what="eff", smooth_mm=4.0, order=1):
    """The three amplitudes of a GIVEN lattice -- used for the harmonics and
    for evaluating one station's lattice on another station's map."""
    w, K, _ = residual_map(f, what, smooth_mm)
    w, xc, yc = _trim(w, f["xc"], f["yc"])
    ks, S = _hex_power(w, xc, yc, f["b"], np.array([d]),
                       np.array([theta_deg]), order)
    return ks[:, 0, 0, :], 2 * S[:, 0, 0] / K.sum()


def lattice_basis(d, theta_deg):
    t = np.radians(theta_deg)
    a1 = d * np.array([np.cos(t), np.sin(t)])
    a2 = d * np.array([np.cos(t + np.pi / 3), np.sin(t + np.pi / 3)])
    return np.column_stack([a1, a2])


def _fold_centre(f, d, theta_deg, x0, what="eff", r=0.35):
    """The observable averaged within `r` of the lattice nodes -- the number
    that says whether a candidate origin is on a pillar."""
    if what == "eff":
        num, den = f["k"], f["n"]
    elif what == "amp":
        num, den = f["asum"], f["an"]
    else:
        num, den = f["csum"], f["an"]
    M = lattice_basis(d, theta_deg)
    X, Y = np.meshgrid(f["xc"], f["yc"], indexing="ij")
    p = np.stack([X.ravel() - x0[0], Y.ravel() - x0[1]])
    fr = np.linalg.solve(M, p)
    n0 = np.round(fr)
    best = None
    for du in (-1, 0, 1):
        for dv in (-1, 0, 1):
            off = p - M @ (n0 + np.array([[du], [dv]]))
            r2 = (off ** 2).sum(0)
            best = r2 if best is None else np.minimum(best, r2)
    m = best < r * r
    dn = den.ravel()[m].sum()
    return float(num.ravel()[m].sum() / dn) if dn > 0 else np.inf


def origin_from_fold(f, d, theta_deg, what="eff", smooth_mm=4.0, nb=30):
    """Where the pillar sits in the cell, taken from the data.

    Fold with an arbitrary origin, then the pillar is wherever the folded cell
    is lowest.  Doing it this way rather than from the harmonic phases avoids a
    branch and a sign convention, either of which lands the origin on one of
    the two interstitial sites -- which for a triangular lattice look enough
    like a real answer to be missed.
    """
    if what == "eff":
        num, den = f["k"], f["n"]
    elif what == "amp":
        num, den = f["asum"], f["an"]
    else:
        num, den = f["csum"], f["an"]
    M = lattice_basis(d, theta_deg)
    X, Y = np.meshgrid(f["xc"], f["yc"], indexing="ij")
    fr = np.linalg.solve(M, np.stack([X.ravel(), Y.ravel()])) % 1.0
    iu = np.minimum((fr[0] * nb).astype(int), nb - 1)
    iv = np.minimum((fr[1] * nb).astype(int), nb - 1)
    idx = iu * nb + iv
    N = np.bincount(idx, den.ravel(), nb * nb).reshape(nb, nb)
    K = np.bincount(idx, num.ravel(), nb * nb).reshape(nb, nb)
    # periodic 3x3 smoothing: the pillar is ~1/3 of a cell across, and a cell
    # split into 30 bins puts only ~1 % of the exposure in each
    sN = sum(np.roll(np.roll(N, i, 0), j, 1)
             for i in (-1, 0, 1) for j in (-1, 0, 1))
    sK = sum(np.roll(np.roll(K, i, 0), j, 1)
             for i in (-1, 0, 1) for j in (-1, 0, 1))
    with np.errstate(invalid="ignore", divide="ignore"):
        v = np.where(sN > 0, sK / np.maximum(sN, 1), np.inf)
    i = np.unravel_index(np.argmin(v), v.shape)
    # sub-bin: the deficit-weighted centroid over the 5x5 around the minimum,
    # so the answer is not quantised to a thirtieth of a cell
    plateau = np.nanmax(v[np.isfinite(v)])
    du = np.arange(-2, 3)
    iu = (i[0] + du[:, None]) % nb
    iv = (i[1] + du[None, :]) % nb
    w = np.clip(plateau - v[iu, iv], 0, None)
    if w.sum() > 0:
        cu = i[0] + float((w.sum(1) * du).sum() / w.sum())
        cv = i[1] + float((w.sum(0) * du).sum() / w.sum())
    else:
        cu, cv = i[0], i[1]
    return M @ np.array([(cu + 0.5) / nb, (cv + 0.5) / nb])


def d_profile(f, theta_deg, what="eff", smooth_mm=4.0,
              d=(2.90, 3.10, 101), order=1):
    """Modulation against lattice spacing at a fixed orientation.

    Only for the picture -- the fit is `hex_scan`.  A common d grid for the
    three stations, so the three curves can be drawn on one axis.
    """
    w, K, _ = residual_map(f, what, smooth_mm)
    w, xc, yc = _trim(w, f["xc"], f["yc"])
    dv = np.linspace(*d[:2], int(d[2]))
    _, S = _hex_power(w, xc, yc, f["b"], dv, np.array([theta_deg]), order)
    A = 2 * S[:, :, 0] / K.sum()
    return dv, np.sqrt((np.abs(A) ** 2).mean(0))


def hex_fit(f, what="eff", d_seed=3.0, smooth_mm=4.0):
    """Measure the pillar lattice: spacing, orientation, origin, amplitudes.

    Nothing here is taken from a drawing.  The three first-order vectors are
    refined independently, so their agreement on |G| and on 60 degrees between
    them is a check that the pattern really is triangular; the origin comes
    from the *phases*, over-determined by the third vector.
    """
    sc = hex_scan(f, what, smooth_mm=smooth_mm)
    Gs, As = sc["G"], list(sc["A"])
    mag = np.hypot(Gs[:, 0], Gs[:, 1])
    d = 4 * np.pi / (np.sqrt(3) * mag)
    # Origin, in two steps.  Rebuilding the modulation field from the three
    # measured harmonics and taking its minimum is the PRECISE answer -- it
    # uses the phases, so it is not quantised by any binning.  What it cannot
    # do on its own is prove it has landed on a pillar rather than on one of
    # the two trigonal interstitials, which for this lattice look the same to
    # a sign error.  So the field gives the position and a direct fold of the
    # data picks between the three candidates.
    x0_fold = origin_from_fold(f, d.mean(), sc["theta_deg"], what, smooth_mm)

    # the modulation field: M(x) = sum_i |A_i| cos(G_i.x - arg A_i)  Solving the phases algebraically has a
    # branch and a sign to get right and lands on an interstitial as readily as
    # on a pillar; the field has neither ambiguity.  Because a triangular
    # lattice has ONE node and TWO interstitials per cell, the depth of the
    # minimum relative to the maxima is itself the check that the pattern is a
    # lattice of dips rather than of bumps.
    ph = np.array([np.angle(a) for a in As])
    a1 = np.array([d.mean(), 0.0])
    a2 = np.array([d.mean() / 2, d.mean() * np.sqrt(3) / 2])
    t = np.linspace(0, 1, 361, endpoint=False)
    U, V = np.meshgrid(t, t, indexing="ij")
    P0 = U[..., None] * a1 + V[..., None] * a2
    M = np.zeros(U.shape)
    for g, a in zip(Gs, As):
        M += np.abs(a) * np.cos(P0 @ g - np.angle(a))
    i = np.unravel_index(np.argmin(M), M.shape)
    x0_field = P0[i]
    resid = float(M[i] / np.abs(As).sum())
    # pick between the node and the two interstitials on the data itself
    tri = (a1 + a2) / 3.0
    cands = [x0_field, x0_field + tri, x0_field - tri]
    depth = [_fold_centre(f, d.mean(), sc["theta_deg"], c, what) for c in cands]
    x0 = cands[int(np.argmin(depth))]
    da = x0 - x0_fold
    fr = np.linalg.solve(np.column_stack([a1, a2]), da)
    fr -= np.round(fr)
    origin_gap = float(np.hypot(*(np.column_stack([a1, a2]) @ fr)))      # -1 iff the three
    # phases agree on one site; a triangular lattice has G3 = G2 - G1, so this
    # is one real constraint, not three.
    out = dict(G=Gs, A=np.array(As), abs_A=np.abs(As), phase=ph,
               d_scan=sc["d"], theta_deg=sc["theta_deg"],
               d=d, d_mean=float(d.mean()), d_std=float(d.std()),
               angle_deg=np.degrees(np.arctan2(Gs[:, 1], Gs[:, 0])),
               x0=x0, x0_field=x0_field, x0_fold=x0_fold,
               origin_gap=origin_gap,
               field_min_ratio=resid, field=M,
               cell_area=float(np.sqrt(3) * d.mean() ** 2 / 2))
    out["scan"] = sc
    # second order: the SAME lattice, so d and theta are not refitted -- if
    # they were, the second harmonic could wander onto a noise peak and the
    # ratio that measures the resolution would be biased up.
    G2, A2 = hex_at(f, sc["d"], sc["theta_deg"], what, smooth_mm, order=2)
    out["G2"] = G2
    out["A2"] = A2
    out["abs_A2"] = np.abs(A2)
    return out


def fold_hex(f, lat, what="eff", nb=40, reach=2.0):
    """Every fine bin folded onto the Wigner-Seitz cell of the measured lattice.

    This is the average pillar: with 2.2 M tracks over a 128 mm window no single
    pillar has more than a handful of tracks on it, but every one of the ~20 000
    inside the beam spot contributes to this one picture.
    """
    if what == "eff":
        num, den = f["k"], f["n"]
    elif what == "amp":
        num, den = f["asum"], f["an"]
    else:
        num, den = f["csum"], f["an"]
    # the basis must carry the fitted orientation: theta is only a few tenths
    # of a degree, but over a 64 mm window that is a fifth of a millimetre of
    # drift, which is the width of the thing being folded
    M = lattice_basis(lat.get("d_scan", lat["d_mean"]),
                      lat.get("theta_deg", 0.0))
    X, Y = np.meshgrid(f["xc"], f["yc"], indexing="ij")
    p = np.stack([X.ravel() - lat["x0"][0], Y.ravel() - lat["x0"][1]])
    frac = np.linalg.solve(M, p)
    # nearest lattice node: round in fractional coords, then fix up over the
    # six neighbours, because rounding in a skew basis is not the nearest point
    best = None
    n0 = np.round(frac)
    for du in (-1, 0, 1):
        for dv in (-1, 0, 1):
            off = p - M @ (n0 + np.array([[du], [dv]]))
            r2 = (off ** 2).sum(0)
            if best is None:
                best, bestr = off, r2
            else:
                m = r2 < bestr
                best = np.where(m, off, best)
                bestr = np.where(m, r2, bestr)
    e = np.linspace(-reach, reach, nb + 1)
    iu = np.digitize(best[0], e) - 1
    iv = np.digitize(best[1], e) - 1
    ok = (iu >= 0) & (iu < nb) & (iv >= 0) & (iv < nb)
    idx = iu[ok] * nb + iv[ok]
    N = np.bincount(idx, den.ravel()[ok], nb * nb).reshape(nb, nb)
    K = np.bincount(idx, num.ravel()[ok], nb * nb).reshape(nb, nb)
    with np.errstate(invalid="ignore", divide="ignore"):
        val = np.where(N > 0, K / np.maximum(N, 1), np.nan)
    return e, val, N, np.sqrt(bestr).reshape(f["n"].shape)


def fold_profile(f, lat, what="eff", rmax=1.6, nb=32):
    """The same fold, as a radial profile about the pillar centre."""
    if what == "eff":
        num, den = f["k"], f["n"]
    elif what == "amp":
        num, den = f["asum"], f["an"]
    else:
        num, den = f["csum"], f["an"]
    _, _, _, r = fold_hex(f, lat, what)
    e = np.linspace(0.0, rmax, nb + 1)
    i = np.digitize(r.ravel(), e) - 1
    ok = (i >= 0) & (i < nb)
    N = np.bincount(i[ok], den.ravel()[ok], nb)
    K = np.bincount(i[ok], num.ravel()[ok], nb)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = np.where(N > 0, K / np.maximum(N, 1), np.nan)
    return 0.5 * (e[1:] + e[:-1]), v, N


# --------------------------------------------------------------------------- #
# what a lattice of dead discs looks like, and reading it backwards
# --------------------------------------------------------------------------- #
def _j1(x):
    """J1 by its series -- G a <= 2 here, where the series is exact to 1e-12."""
    x = np.asarray(x, float)
    s, term = np.zeros_like(x), x / 2
    for m in range(14):
        s = s + term
        term = term * -(x * x / 4) / ((m + 1) * (m + 2))
    return s


def disc_amplitude(G, a, sigma, cell_area, f_dead=1.0):
    """|A| of one reciprocal-lattice vector for a lattice of dead discs.

    A disc of radius a, one per cell of area A_c, gives a Fourier coefficient
    2 (pi a^2 / A_c) * 2 J1(G a)/(G a); a Gaussian pointing error sigma
    multiplies it by exp(-G^2 sigma^2 / 2).  Nothing else in the chain has a
    scale between 0.2 and 3 mm, which is what makes this invertible.
    """
    G = np.asarray(G, float)
    Ga = np.maximum(G * a, 1e-9)
    return (2 * np.pi * a ** 2 / cell_area * (2 * _j1(Ga) / Ga)
            * np.exp(-0.5 * G ** 2 * sigma ** 2) * f_dead)


def solve_lattice(A1, G1, A2, G2, cell_area, f_dead=1.0,
                  a_grid=None, s_grid=None):
    """(effective dead radius, pointing sigma) from the first two harmonics.

    Two measurements, two unknowns.  The ratio A2/A1 is where sigma comes from
    -- it does not depend on how big or how dead the pillar is, only on how
    much the smearing costs at twice the wavenumber -- and A1 then sets the
    radius.  `f_dead` is held at 1, so `a` is the radius of the fully dead disc
    that would produce the same modulation, not a claim about the artwork.
    """
    a_grid = np.arange(0.05, 1.201, 0.002) if a_grid is None else a_grid
    s_grid = np.arange(0.02, 1.501, 0.002) if s_grid is None else s_grid
    A, S = np.meshgrid(a_grid, s_grid, indexing="ij")
    m1 = disc_amplitude(G1, A, S, cell_area, f_dead)
    m2 = disc_amplitude(G2, A, S, cell_area, f_dead)
    cost = ((m1 - A1) / A1) ** 2 + ((m2 - A2) / A2) ** 2
    i = np.unravel_index(np.argmin(cost), cost.shape)
    a, sg = float(A[i]), float(S[i])
    return dict(a=a, sigma=sg, coverage=float(np.pi * a * a / cell_area),
                resid=float(np.sqrt(cost[i])),
                A1_model=float(disc_amplitude(G1, a, sg, cell_area, f_dead)),
                A2_model=float(disc_amplitude(G2, a, sg, cell_area, f_dead)))


# --------------------------------------------------------------------------- #
# dead-area masks
# --------------------------------------------------------------------------- #
def pad_frame(z, s, npad=None):
    """Per-pad ledger as plain arrays: n, k, mean amplitude, connector block."""
    n = z[f"pad_n_{s}"].astype(float)
    k = z[f"pad_k_{s}"].astype(float)
    an = z[f"pad_an_{s}"].astype(float)
    asum = z[f"pad_asum_{s}"].astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        eff = np.where(n > 0, k / np.maximum(n, 1), np.nan)
        amp = np.where(an > 0, asum / np.maximum(an, 1), np.nan)
    return dict(n=n, k=k, eff=eff, amp=amp,
                block=np.arange(len(n)) // 64)


def build_mask(z, s, min_n=100, block_frac=0.5, dead_max_k=0):
    """Which channel_ids to drop, and why.

    Three layers, none of which is a per-pad efficiency cut:

      dead    -- the channel never responded at all (k == 0 with real exposure).
                 That is a hardware statement, not a performance one.
      block   -- a whole 64-channel FEU input whose efficiency is below
                 `block_frac` x the station's median block.  The decision is one
                 number per cable, taken over 64 pads, and the whole cable goes;
                 it is the granularity at which the P2_IN connector-11 problem
                 actually lives.

    The gerber pillars are not in here: they are a position cut, not a channel
    cut, and are applied to the track in `stack_masks`.
    """
    p = pad_frame(z, s)
    ok = p["n"] >= min_n
    dead = ok & (p["k"] <= dead_max_k)
    beff = np.full(p["block"].max() + 1, np.nan)
    for b in range(len(beff)):
        m = ok & (p["block"] == b) & ~dead
        if m.sum() >= 3:
            beff[b] = p["k"][m].sum() / p["n"][m].sum()
    med = np.nanmedian(beff)
    weak_blocks = np.where(beff < block_frac * med)[0]
    weak = ok & np.isin(p["block"], weak_blocks)
    return dict(dead=np.where(dead)[0], weak=np.where(weak)[0],
                weak_blocks=weak_blocks.tolist(), block_eff=beff,
                median_block_eff=float(med),
                drop=np.union1d(np.where(dead)[0], np.where(weak)[0]))


def stack_masks(c, z, geom, masks, pillar_margin=1.0):
    """Per-station keep-flags, and the AND of them.

    A three-station track needs all three stations clean, so the stack mask is
    the intersection -- exactly as the vetoes are ANDed in the extraction.
    """
    keep = {}
    for s in STATIONS:
        x, y = c[f"ex_{s}"], c[f"ey_{s}"]
        k = G.big_pillar_mask(x, y, geom, margin=pillar_margin)
        pid = np.nan_to_num(c[f"pid_{s}"], nan=-1).astype(int)
        k &= ~np.isin(pid, masks[s]["drop"])
        keep[s] = k
    keep["all"] = np.logical_and.reduce([keep[s] for s in STATIONS])
    return keep


# --------------------------------------------------------------------------- #
# the stack, with and without the dead areas
# --------------------------------------------------------------------------- #
def station_eff(z, s, drop=()):
    """Per-station efficiency from the pad ledger, over ALL fiducial tracks.

    The ledger is indexed by channel_id, so a mask is an index list; it cannot
    express the gerber pillars (a position cut), which is why the stack numbers
    below are computed from the per-track sample instead.
    """
    p = pad_frame(z, s)
    keep = np.ones(len(p["n"]), bool)
    keep[np.asarray(drop, int)] = False
    n, k = p["n"][keep].sum(), p["k"][keep].sum()
    return dict(n=float(n), k=float(k), eff=float(k / max(n, 1)),
                n_pad_dropped=int((~keep & (p["n"] > 0)).sum()),
                n_track_dropped=float(p["n"][~keep].sum()))


def line_through(zs, u):
    zs = np.asarray(zs, float)
    dz = zs - zs.mean()
    sl = (u * dz).sum(1) / float((dz ** 2).sum())
    return sl, u.mean(1) - sl * zs.mean()


def iqr_sigma(v):
    q = np.percentile(v, [25, 75])
    return float((q[1] - q[0]) / 1.349)


def stack_block(c, keep, zs=(320.0, 630.0, 940.0), z_exit=940.0):
    """The three-station track against the reference, for a subset of tracks."""
    zs = np.asarray(zs, float)
    found = {s: (c[f"found_{s}"] > 0) for s in STATIONS}
    n3 = np.sum([found[s] for s in STATIONS], axis=0)
    out = {"n_fid_all3": int(keep.sum()),
           "n_station_frac": [float(((n3 == kk) & keep).mean()
                                    / max(keep.mean(), 1e-12))
                              for kk in range(4)]}
    for s in STATIONS:
        out[f"eff_{s}"] = float(found[s][keep].mean())
    m = keep & (n3 == 3)
    out["n_3of3"] = int(m.sum())
    out["frac_3of3"] = out["n_station_frac"][3]
    for ax in ("x", "y"):
        u = np.column_stack([c[f"u{ax}_{s}"][m] for s in STATIONS])
        sl, ic = line_through(zs, u)
        sl_ref, ic_ref = c[f"slope_{ax}"][m], c[f"{ax}f"][m]
        d = (ic + sl * z_exit) - (ic_ref + sl_ref * z_exit)
        a = (sl - sl_ref) * 1e3
        out[f"dpos_{ax}"] = dict(sigma_iqr=iqr_sigma(d), rms=float(d.std()),
                                 p95=float(np.percentile(np.abs(d), 95)),
                                 median=float(np.median(d)))
        out[f"dang_{ax}"] = dict(sigma_iqr=iqr_sigma(a), rms=float(a.std()),
                                 p95=float(np.percentile(np.abs(a), 95)),
                                 median=float(np.median(a)))
        # per-station residual to the reference, same subset
        for s in STATIONS:
            f = keep & found[s]
            r = c[f"u{ax}_{s}"][f] - (c[f"{ax}f"][f]
                                      + c[f"slope_{ax}"][f] * Z[s])
            r = r[np.abs(r) < 25.0]
            out[f"res_{ax}_{s}"] = dict(rms=float(r.std()),
                                        sigma_iqr=iqr_sigma(r), n=int(len(r)))
    return out


# --------------------------------------------------------------------------- #
# what the lattice costs, without a model
# --------------------------------------------------------------------------- #
def lattice_deficit(f, lat, what="eff", r_plateau=1.35):
    """Exposure-weighted deficit per lattice cell, straight from the fold.

    No disc, no Gaussian: take the value averaged over the part of the cell
    further than `r_plateau` from a pillar as the undisturbed level, and
    compare the cell average with it.  The Wigner-Seitz cell of a 3 mm
    triangular lattice reaches 1.73 mm, so r > 1.35 mm is a fifth of the cell
    and four smearing widths clear of the widest pillar footprint measured.
    """
    if what == "eff":
        num, den = f["k"], f["n"]
    elif what == "amp":
        num, den = f["asum"], f["an"]
    else:
        num, den = f["csum"], f["an"]
    _, _, _, r = fold_hex(f, lat, what)
    far = (r > r_plateau) & (den > 0)
    live = den > 0
    plateau = num[far].sum() / den[far].sum()
    mean = num[live].sum() / den[live].sum()
    return dict(plateau=float(plateau), mean=float(mean),
                deficit=float(plateau - mean),
                rel_deficit=float((plateau - mean) / plateau),
                far_frac=float(den[far].sum() / den[live].sum()))


try:                                     # scipy is present in the venv, but
    from scipy.special import erf as _erf  # this module should not need it
except ImportError:                       # pragma: no cover
    def _erf(x):
        """Abramowitz & Stegun 7.1.26 -- 1.5e-7 absolute, and it takes an
        array, which `math.erf` under `np.vectorize` does not."""
        x = np.asarray(x, float)
        t = 1.0 / (1.0 + 0.3275911 * np.abs(x))
        y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t
                    - 0.284496736) * t + 0.254829592) * t * np.exp(-x * x)
        return np.sign(x) * y


def edge_fit(z, s, ax="rw", tag="single", w_grid=None, s_grid=None):
    """Pointing resolution from the edge of the pad-residual box.

    A one-pad cluster reports its pad centre, so the residual to the reference
    is a box of the pad's own width convolved with the reference's pointing
    error and nothing else.  Fitting box width and edge width together means
    the answer does not depend on knowing the pitch -- which is just as well,
    because the pads are a fan and the pitch is not one number.
    """
    e = z["edge_edges"]
    c = 0.5 * (e[1:] + e[:-1])
    h = z[f"edge_{ax}_{s}_{tag}"].astype(float)
    if h.sum() < 100:
        return None
    w_grid = np.arange(9.5, 13.5, 0.02) if w_grid is None else w_grid
    s_grid = np.arange(0.03, 1.20, 0.005) if s_grid is None else s_grid
    W = w_grid[:, None, None]
    S = s_grid[None, :, None]
    X = c[None, None, :]
    mod = 0.5 * (_erf((W / 2 - X) / (S * np.sqrt(2)))
                 + _erf((W / 2 + X) / (S * np.sqrt(2))))
    amp = (mod * h).sum(-1) / np.maximum((mod * mod).sum(-1), 1e-9)
    chi = ((amp[..., None] * mod - h) ** 2 / np.maximum(h, 1)).sum(-1)
    i = np.unravel_index(np.argmin(chi), chi.shape)
    return dict(width=float(w_grid[i[0]]), sigma=float(s_grid[i[1]]),
                chi2_ndf=float(chi[i] / (len(c) - 3)), n=float(h.sum()),
                curve=amp[i] * mod[i[0], i[1]], centres=c, hist=h)
