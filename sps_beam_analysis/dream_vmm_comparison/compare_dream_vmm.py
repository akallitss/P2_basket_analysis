#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""compare_dream_vmm.py -- P2_OUT pad-by-pad pulse height, VMM vs DREAM.

Two readouts, the same 512-pad chamber, the same uRWELL reference, the same
track selection, the same operating point (mesh 450 V / drift 750 V):

    VMM     run_46 / cfg_gain4.5_peaktime200   1 Aug 04:53   sdt 224 on all six
    DREAM   eff_nominal_1 / eff_nominal_00..09  27 Jul 10:40

For every pad the two chains give the same pair of numbers -- how many uRWELL
tracks pointed at it, and the pulse-height distribution of the ones it
recorded -- so the pad-to-pad gain spread can be compared directly.

The comparison is in RATIO to each readout's own fleet median, because the two
ADCs share no scale: VMM is a 10-bit peak ADC whose spectrum is cut off from
below by the discriminator, DREAM a 12-bit one that saturates ~3700 at the top.

Variance decomposition
----------------------
If the two electronics chains fail independently, the covariance of the two
per-pad relative pulse heights is the DETECTOR's own pad-to-pad structure, and
what each readout shows on top of that is its own:

    var(V) = var_det + var_V,   var(D) = var_det + var_D,   cov(V,D) = var_det

Reported with a bootstrap interval, and repeated after regressing out pad area
(the fan pads differ in size by 2x, which is geometry, not gain).
"""
import json
import numpy as np
import pandas as pd

VMM_JSON = "data/report_run_46_rate.json"
DREAM_CSV = "data/dream_padadc_eff_nominal_1_P2_OUT.csv"
DREAM_NPZ = "data/dream_padadc_eff_nominal_1_P2_OUT_hist.npz"
GEOM_CSV = "data/pad_adc_run_46.csv"          # carries pad_area / radius
STATION = "P2_OUT"

MIN_TRACK = 500      # a pad the beam actually illuminated
MIN_AMP_N = 200      # enough recorded pulses for a per-pad median


def vmm_side():
    d = json.load(open(VMM_JSON))
    st = d["stations"][STATION]
    pp = st["per_pad"]
    hb = pp["adc_hist_bin"]
    nb = 1024 // hb
    n = len(pp["channel_id"])
    hist = np.zeros((n, nb), np.int64)
    for k, v in pp["adc_hist"].items():
        hist[int(k)] = v
    t = pd.DataFrame({
        "pad_id": pp["channel_id"], "x": pp["x"], "y": pp["y"],
        "vmm": pp["vmm"], "ch": pp["ch"],
        "n_track": pp["n"], "k_hit": pp["k"],
        "amp_n": pp["adc_n"], "amp_med": pp["adc_med_dnl"],
        "amp_p25": pp["adc_p25"], "amp_p75": pp["adc_p75"]})
    t["eff"] = np.where(t["n_track"] > 0, t["k_hit"] / t["n_track"].clip(lower=1),
                        np.nan)
    return t, hist, float(hb), st


def dream_side():
    """`_own` = the pad the track pointed at also led the event, which is the
    VMM's rule exactly (it takes the hit nearest the prediction, not the
    largest).  It is 87.9 % of the recorded pulses and moves the per-pad
    relative median by nothing (rms 0.391 -> 0.394, r = 0.998), so the choice
    is not load-bearing -- the strict one is used because it is the identical
    quantity, and that is the point of the comparison."""
    t = pd.read_csv(DREAM_CSV).rename(columns={"pad_cx": "x", "pad_cy": "y",
                                               "amp_med": "amp_med_lead",
                                               "amp_med_own": "amp_med",
                                               "amp_n": "amp_n_lead",
                                               "amp_n_own": "amp_n"})
    z = np.load(DREAM_NPZ)
    return t, z["hist_own"], float(z["amp_bin"])


def hist_quantile(h, q, binw):
    tot = h.sum()
    if tot == 0:
        return np.nan
    cum = np.cumsum(h) / tot
    i = int(np.searchsorted(cum, q, side="left"))
    below = cum[i - 1] if i > 0 else 0.0
    return float((i + (q - below) / max(cum[i] - below, 1e-12)) * binw)


def mpv(h, binw, smooth=3, hi=None):
    """Landau peak: argmax of a smoothed histogram, parabolically interpolated.

    `hi` blanks the top of the range first.  Both ADCs pile up in their last
    bin -- 1.0 % of the VMM entries at 1016, 1.7 % of the DREAM ones at 4092 --
    and on a pad whose Landau peak is low that overflow spike is the tallest
    bin in the histogram, so an unrestricted argmax returns the saturation
    edge instead of the MPV.
    """
    if h.sum() == 0:
        return np.nan
    h = h.astype(float).copy()
    if hi is not None:
        h[hi:] = 0.0
    k = np.ones(smooth) / smooth
    s = np.convolve(h, k, mode="same")
    i = int(s.argmax())
    if 0 < i < len(s) - 1:
        a, b, c = s[i - 1], s[i], s[i + 1]
        den = a - 2 * b + c
        d = 0.5 * (a - c) / den if den != 0 else 0.0
        d = float(np.clip(d, -0.5, 0.5))
    else:
        d = 0.0
    return (i + 0.5 + d) * binw


def relspread(r):
    """Robust and plain spread of a per-pad ratio, both as a fraction."""
    q = np.percentile(r, [25, 75])
    return dict(sigma_iqr=float((q[1] - q[0]) / 1.349), rms=float(np.std(r, ddof=1)),
                p05=float(np.percentile(r, 5)), p95=float(np.percentile(r, 95)),
                lo=float(r.min()), hi=float(r.max()))


def decompose(v, d, nboot=2000, seed=0):
    """var_det = cov(v, d); the rest belongs to each chain.  Bootstrap over pads."""
    def one(iv):
        a, b = v[iv], d[iv]
        cov = float(np.cov(a, b, ddof=1)[0, 1])
        return np.var(a, ddof=1), np.var(b, ddof=1), cov
    n = len(v)
    va, vb, cov = one(np.arange(n))
    rng = np.random.default_rng(seed)
    bs = np.array([one(rng.integers(0, n, n)) for _ in range(nboot)])
    def ci(x):
        return [float(np.percentile(x, 16)), float(np.percentile(x, 84))]
    return dict(
        var_vmm=float(va), var_dream=float(vb), cov=float(cov),
        r=float(cov / np.sqrt(va * vb)),
        sigma_vmm=float(np.sqrt(va)), sigma_dream=float(np.sqrt(vb)),
        sigma_shared=float(np.sqrt(max(cov, 0.0))),
        sigma_vmm_only=float(np.sqrt(max(va - cov, 0.0))),
        sigma_dream_only=float(np.sqrt(max(vb - cov, 0.0))),
        ci_sigma_shared=[float(np.sqrt(max(x, 0))) for x in ci(bs[:, 2])],
        ci_r=ci(bs[:, 2] / np.sqrt(bs[:, 0] * bs[:, 1])),
        n_pads=int(n))


def build():
    v, vh, vbin, vst = vmm_side()
    d, dh, dbin = dream_side()
    geom = pd.read_csv(GEOM_CSV)
    geom = geom[geom["station"] == STATION][
        ["channel_id", "pad_area", "radius", "connector_N", "half", "masked"]]
    geom = geom.rename(columns={"channel_id": "pad_id"})

    # index of each pad_id in its own histogram array, so spectra can be pulled
    v = v.reset_index().rename(columns={"index": "vi"})
    d = d.reset_index().rename(columns={"index": "di"})
    m = v.merge(d, on="pad_id", suffixes=("_v", "_d")).merge(geom, on="pad_id",
                                                             how="left")
    # the two chains must agree on where the pad is, or the ids mean different
    # pads and every number below is a coincidence
    dx = np.hypot(m["x_v"] - m["x_d"], m["y_v"] - m["y_d"])
    assert float(np.nanmax(dx)) < 0.5, f"pad positions disagree by {dx.max():.2f} mm"

    # blank the saturation shoulder: VMM above 1008 ADC, DREAM above 3660
    m["mpv_v"] = [mpv(vh[i], vbin, smooth=3, hi=int(1008 / vbin)) for i in m["vi"]]
    m["mpv_d"] = [mpv(dh[i], dbin, smooth=9, hi=int(3660 / dbin)) for i in m["di"]]
    m["use"] = ((m["n_track_v"] >= MIN_TRACK) & (m["n_track_d"] >= MIN_TRACK)
                & (m["amp_n_v"] >= MIN_AMP_N) & (m["amp_n_d"] >= MIN_AMP_N)
                & ~m["masked"].fillna(False))
    return m, vh, vbin, dh, dbin, vst


def main():
    m, vh, vbin, dh, dbin, vst = build()
    g = m[m["use"]].copy()
    print(f"{len(m)} pads in common, {len(g)} used "
          f"(>= {MIN_TRACK} tracks and >= {MIN_AMP_N} pulses in BOTH)")

    for tag, cv, cd in (("median", "amp_med_v", "amp_med_d"),
                        ("MPV", "mpv_v", "mpv_d")):
        rv = g[cv] / g[cv].median()
        rd = g[cd] / g[cd].median()
        sv, sd = relspread(rv.to_numpy()), relspread(rd.to_numpy())
        dec = decompose(rv.to_numpy(), rd.to_numpy())
        print(f"\n--- per-pad {tag}, relative to each readout's fleet {tag} ---")
        print(f"  VMM    sigma_iqr {sv['sigma_iqr']:.3f}  rms {sv['rms']:.3f}  "
              f"5-95% {sv['p05']:.2f}-{sv['p95']:.2f}  full {sv['lo']:.2f}-{sv['hi']:.2f}")
        print(f"  DREAM  sigma_iqr {sd['sigma_iqr']:.3f}  rms {sd['rms']:.3f}  "
              f"5-95% {sd['p05']:.2f}-{sd['p95']:.2f}  full {sd['lo']:.2f}-{sd['hi']:.2f}")
        print(f"  ratio of spreads (rms) VMM/DREAM = {sv['rms'] / sd['rms']:.2f}")
        print(f"  r(VMM, DREAM) = {dec['r']:+.3f}  "
              f"[{dec['ci_r'][0]:+.3f}, {dec['ci_r'][1]:+.3f}]")
        print(f"  shared (detector)  sigma {dec['sigma_shared']:.3f} "
              f"[{dec['ci_sigma_shared'][0]:.3f}, {dec['ci_sigma_shared'][1]:.3f}]")
        print(f"  VMM-only   sigma {dec['sigma_vmm_only']:.3f}")
        print(f"  DREAM-only sigma {dec['sigma_dream_only']:.3f}")

    # pad area is 2:1 across the fan and is geometry, not gain: report the same
    # numbers with it regressed out of both
    A = np.column_stack([np.ones(len(g)), g["pad_area"].to_numpy()])
    print("\n--- with pad area regressed out (both readouts) ---")
    for tag, cv, cd in (("median", "amp_med_v", "amp_med_d"),
                        ("MPV", "mpv_v", "mpv_d")):
        out = {}
        for lab, c in (("V", cv), ("D", cd)):
            y = (g[c] / g[c].median()).to_numpy()
            b, *_ = np.linalg.lstsq(A, y, rcond=None)
            out[lab] = y - A @ b + 1.0
            print(f"  r(pad_area, {tag} {lab}) = "
                  f"{np.corrcoef(g['pad_area'], y)[0, 1]:+.3f}")
        dec = decompose(out["V"], out["D"])
        print(f"  {tag}: VMM rms {np.std(out['V'], ddof=1):.3f}  "
              f"DREAM rms {np.std(out['D'], ddof=1):.3f}  "
              f"r {dec['r']:+.3f}  shared sigma {dec['sigma_shared']:.3f}")

    # ---------------------------------------------------------------- #
    # The one control that decides the headline.  The VMM spectrum is cut off
    # from below at ~64 ADC (its pooled MPV is 112), and a fixed cut RAISES the
    # median of a low-gain pad more than that of a high-gain one -- so it
    # compresses the pad-to-pad spread.  Put the same relative cut on DREAM and
    # see whether its spread collapses onto the VMM's.
    v_thr, v_mpv, d_mpv = 64.0, 111.8, 237.5
    scale = v_mpv / d_mpv                       # VMM ADC per DREAM ADC
    matched = v_thr / scale
    print(f"\n--- DREAM with a low-amplitude cut (VMM MPV {v_mpv:.0f} / DREAM "
          f"MPV {d_mpv:.0f} -> {scale:.3f} VMM ADC per DREAM ADC;")
    print(f"    the VMM's ~{v_thr:.0f} ADC turn-on is {matched:.0f} DREAM ADC) ---")
    rv = (g["amp_med_v"] / g["amp_med_v"].median()).to_numpy()
    di = g["di"].to_numpy()
    print(f"  {'cut':>6}  {'rms':>6}  {'sig_iqr':>7}  {'r(V,D)':>7}  "
          f"{'slope V~D':>9}  eff_kept")
    for cut in (0, 40, 80, matched, 160, 200, 240):
        b0 = int(np.ceil(cut / dbin))
        hh = dh[di][:, b0:]
        med = np.array([hist_quantile(h, 0.5, dbin) for h in hh]) + b0 * dbin
        keep = hh.sum() / dh[di].sum()
        r = med / np.median(med)
        sl = np.polyfit(r, rv, 1)[0]
        tag = "  <-- matched" if abs(cut - matched) < 1e-9 else ""
        print(f"  {cut:6.0f}  {np.std(r, ddof=1):6.3f}  "
              f"{(np.percentile(r, 75) - np.percentile(r, 25)) / 1.349:7.3f}  "
              f"{np.corrcoef(rv, r)[0, 1]:+7.3f}  {sl:9.3f}  {keep:7.3f}{tag}")
    print(f"  VMM for comparison: rms {np.std(rv, ddof=1):.3f}")

    print("\n--- efficiency ---")
    print(f"  VMM   {g['k_hit_v'].sum() / g['n_track_v'].sum():.4f} on "
          f"{int(g['n_track_v'].sum())} tracks")
    print(f"  DREAM {g['k_hit_d'].sum() / g['n_track_d'].sum():.4f} on "
          f"{int(g['n_track_d'].sum())} tracks")

    m.to_csv("data/compare_dream_vmm_P2_OUT.csv", index=False)
    print("\nwrote data/compare_dream_vmm_P2_OUT.csv")


if __name__ == "__main__":
    main()
