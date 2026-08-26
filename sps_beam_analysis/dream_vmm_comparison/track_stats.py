#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""track_stats.py -- load and aggregate what p2_selftrack.py measured.

The NPZ carries counts, not ratios, so sub-runs add.  Every width quoted
downstream is therefore derived from the SUMMED histogram rather than averaged
over per-sub-run widths -- averaging widths would weight a short sub-run the
same as a long one and hide any drift instead of showing it.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

STATIONS = ("P2_IN", "P2_MID", "P2_OUT")
Z = {"P2_IN": 320.0, "P2_MID": 630.0, "P2_OUT": 940.0}
PITCH_X, PITCH_Y = 12.05, 11.60          # measured pad pitch [mm]
RUN = "eff_nominal_1"


def load(run=RUN):
    z = np.load(os.path.join(DATA, f"p2_selftrack_{run}.npz"), allow_pickle=True)
    j = json.load(open(os.path.join(DATA, f"p2_selftrack_{run}.json")))
    return z, j


# --------------------------------------------------------------------------- #
def centres(edges):
    return 0.5 * (edges[1:] + edges[:-1])


def hstats(counts, edges):
    """rms, IQR-sigma and the 95th percentile of |d|, from a histogram.

    Bin centres are exact enough here: the residual bins are 0.25 mm against a
    12 mm pad, and the percentile interpolates within the bin it lands in.
    """
    c = np.asarray(counts, float)
    n = c.sum()
    if n < 10:
        return dict(n=0)
    x = centres(edges)
    mean = float((c * x).sum() / n)
    rms = float(np.sqrt((c * (x - mean) ** 2).sum() / n))
    cum = np.concatenate([[0.0], np.cumsum(c)]) / n

    def q(p):
        return float(np.interp(p, cum, edges))

    p25, p75 = q(0.25), q(0.75)
    # |d| percentile: fold the histogram about zero first
    ax = np.abs(x)
    o = np.argsort(ax)
    ca = np.cumsum(c[o]) / n
    p95 = float(np.interp(0.95, ca, ax[o]))
    return dict(n=int(n), mean=mean, median=q(0.5), rms=rms,
                sigma_iqr=(p75 - p25) / 1.349, p95_abs=p95)


def sample(z):
    """The joined per-track table, as a dict of columns."""
    cols = [str(c) for c in z["sample_cols"]]
    S = z["sample"].astype(np.float64)
    return {k: S[:, i] for i, k in enumerate(cols)}, len(S)


def all3_mask(c):
    m = np.ones(len(next(iter(c.values()))), bool)
    for s in STATIONS:
        m &= (c[f"found_{s}"] > 0) & (c[f"fid_{s}"] > 0)
    return m


# --------------------------------------------------------------------------- #
def station_block(z, j):
    """Per-station efficiency and pointing, summed over sub-runs."""
    out = {}
    for s in STATIONS:
        k = sum(b["stations"][s]["n_found"] for b in j)
        n = sum(b["stations"][s]["n_fid"] for b in j)
        blk = dict(z=Z[s], k=k, n=n, eff=k / max(n, 1))
        for ax in ("x", "y"):
            hs = z[f"res_{ax}_{s}_single"]
            hm = z[f"res_{ax}_{s}_multi"]
            e = z["res_edges"]
            blk[f"res_{ax}"] = hstats(hs + hm, e)
            blk[f"res_{ax}_single"] = hstats(hs, e)
            blk[f"res_{ax}_multi"] = hstats(hm, e)
        blk["frac_multi"] = float(
            z[f"res_x_{s}_multi"].sum()
            / max(z[f"res_x_{s}_single"].sum() + z[f"res_x_{s}_multi"].sum(), 1))
        out[s] = blk
    return out


def selftrack_block(z, j, zs=(320.0, 630.0, 940.0)):
    """The P2-only track against the reference, from the sampled table.

    The sample is every k-th row of the fiducial table, so it is unbiased; the
    z-curve is computed from it rather than from the per-sub-run JSON because
    a percentile does not average.
    """
    c, _ = sample(z)
    m = all3_mask(c)
    d = {k: v[m] for k, v in c.items()}
    zs = np.asarray(zs, float)
    zbar, dzc = zs.mean(), zs - zs.mean()
    denom = float((dzc ** 2).sum())

    out = {"n": int(m.sum()), "n_fid_all3": int(len(m))}
    for ax in ("x", "y"):
        u = np.column_stack([d[f"u{ax}_{s}"] for s in STATIONS])
        sl = (u * dzc).sum(1) / denom
        ic = u.mean(1) - sl * zbar
        sl_ref, ic_ref = d[f"slope_{ax}"], d[f"{ax}f"]
        out[f"dang_{ax}"] = (sl - sl_ref) * 1e3          # mrad
        zgrid = np.arange(0.0, 3001.0, 50.0)
        out[f"zgrid"] = zgrid
        out[f"curve_{ax}"] = [(ic + sl * zz) - (ic_ref + sl_ref * zz)
                              for zz in zgrid]
        # the exit plane of the basket, quoted on its own: this is the number
        # a downstream detector inherits
        out[f"dpos_exit_{ax}"] = ((ic + sl * zs[2])
                                  - (ic_ref + sl_ref * zs[2]))
        # the two rulers on the middle station
        u_io = np.column_stack([d[f"u{ax}_P2_IN"], d[f"u{ax}_P2_OUT"]])
        out[f"mid_self_{ax}"] = d[f"u{ax}_P2_MID"] - u_io.mean(1)
        out[f"mid_ref_{ax}"] = d[f"u{ax}_P2_MID"] - (ic_ref + sl_ref * zs[1])

    # how often do all three report the SAME pad?  That is the whole reason
    # the self-consistency number is not a resolution.
    px = np.column_stack([d[f"ex_{s}"] + d[f"dx_{s}"] for s in STATIONS])
    py = np.column_stack([d[f"ey_{s}"] + d[f"dy_{s}"] for s in STATIONS])
    same = ((px.max(1) - px.min(1) < 1.0) & (py.max(1) - py.min(1) < 1.0))
    out["frac_same_pad"] = float(same.mean())
    out["same_pad"] = same

    nst = np.zeros(len(m), int)
    for s in STATIONS:
        nst += ((c[f"found_{s}"] > 0) & (c[f"fid_{s}"] > 0)).astype(int)
    out["n_station_frac"] = [float((nst == kk).mean()) for kk in range(4)]
    return out


def spread(v, p=(25, 75)):
    q = np.percentile(v, p)
    return float((q[1] - q[0]) / 1.349)


def padmap(z, s, what="eff", min_n=200):
    """Folded onto the face of a single pad: efficiency, mean leading-pad
    amplitude, or mean cluster size, against where the track landed relative
    to the centre of the pad it pointed at."""
    e = z["pad_edges"]
    n, k = z[f"padmap_n_{s}"], z[f"padmap_k_{s}"]
    if what == "eff":
        num, den = k, n
    elif what == "amp":
        num, den = z[f"padmap_asum_{s}"], z[f"padmap_an_{s}"]
    else:
        num, den = z[f"padmap_csum_{s}"], z[f"padmap_an_{s}"]
    with np.errstate(invalid="ignore", divide="ignore"):
        v = np.where(den >= min_n, num / np.maximum(den, 1), np.nan)
    return e, v


def padprofile(z, s, what="eff", nb=14, rmax=6.4):
    """Radial profile across the face of one pad.

    `eff` comes from the 1D counters; `amp` and `nclus` are re-binned from the
    2D maps, which is exact -- they are counts and sums, so summing them into
    an annulus is the same as having histogrammed in r to begin with.
    """
    e = z["pad_edges"]
    edges = np.linspace(0.0, rmax, nb + 1)
    c = 0.5 * (edges[1:] + edges[:-1])
    if what == "eff":
        n2, k2 = z[f"padmap_n_{s}"], z[f"padmap_k_{s}"]
    elif what == "amp":
        n2, k2 = z[f"padmap_an_{s}"], z[f"padmap_asum_{s}"]
    else:
        n2, k2 = z[f"padmap_an_{s}"], z[f"padmap_csum_{s}"]
    ctr = 0.5 * (e[1:] + e[:-1])
    X, Y = np.meshgrid(ctr, ctr, indexing="ij")
    R = np.hypot(X, Y).ravel()
    ib = np.digitize(R, edges) - 1
    ok = (ib >= 0) & (ib < nb)
    num = np.bincount(ib[ok], k2.ravel()[ok], nb)
    den = np.bincount(ib[ok], n2.ravel()[ok], nb)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = np.where(den > 0, num / np.maximum(den, 1), np.nan)
    return c, v, den


def effmap(z, s):
    """(x edges, y edges, efficiency, tracks) in the station's own pad frame."""
    xe, ye = z[f"map_edges_x_{s}"], z[f"map_edges_y_{s}"]
    n, k = z[f"map_n_{s}"], z[f"map_k_{s}"]
    with np.errstate(invalid="ignore", divide="ignore"):
        e = np.where(n > 0, k / np.maximum(n, 1), np.nan)
    return xe, ye, e, n


# --------------------------------------------------------------------------- #
def trim(run=RUN, rows=60000):
    """Shrink the committed NPZ by thinning its per-track sample.

    The sample is 99 % of the file and feeds only the self-track block; every
    exact number in the report comes from the histograms and maps, which are
    kept whole. The full sample stays on EOS under
    `analysis/selftrack/<run>/`. Thinning takes every k-th row, so it is
    deterministic and unbiased.
    """
    p = os.path.join(DATA, f"p2_selftrack_{run}.npz")
    z = np.load(p, allow_pickle=True)
    out = {k: z[k] for k in z.files}
    S = out["sample"]
    if len(S) > rows:
        out["sample"] = S[::int(np.ceil(len(S) / rows))]
    np.savez_compressed(p, **out)
    print(f"{p}: sample {len(S)} -> {len(out['sample'])} rows, "
          f"{os.path.getsize(p) / 1e6:.1f} MB")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--trim", type=int, metavar="ROWS",
                    help="thin the sample in data/p2_selftrack_<run>.npz")
    ap.add_argument("--run", default=RUN)
    a = ap.parse_args()
    if a.trim:
        trim(a.run, a.trim)
