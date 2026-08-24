#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
decide.py -- beam-test gas-change decision table for the P2 Micromegas.

analyze.py answers "which mesh HV reproduces my gain in the new gas".  This
script answers the wider question asked at the beam: for each candidate gas,
what do I get in GAIN / HV / DRIFT SPEED / STABILITY / EFFICIENCY / TIMING, all
in one table, with the timing column referred to the < 20 ns goal.

It is deliberately tolerant of PARTIAL scans: the drift rows are written first
by p2_gas_scan, so the drift/timing/efficiency columns become available a couple
of minutes into a run, while the (slow) Townsend rows are still being computed.
Missing quantities print as "--" rather than aborting.

Timing model (gas-limited floor, 3 mm conversion gap)
-----------------------------------------------------
A track deposits primary clusters at an exponentially distributed depth with
mean free path 1/n_cl.  If the front-end fires on the FIRST arriving cluster,
the arrival-time jitter is the drift time of that exponential depth:

    sigma_first = 1 / (n_cl * v_d)

If it needs k clusters to cross threshold, sigma grows as sqrt(k)/(n_cl*v_d).
Longitudinal diffusion adds, over a full-gap drift L:

    sigma_diff  = D_L * sqrt(L) / v_d          (Magboltz D_L in cm^1/2)

and the two add in quadrature.  Both terms scale as 1/v_d, so the drift
velocity of the gas is the single lever the mixture gives you on timing.  This
is a FLOOR: it contains no electronics jitter, no time-walk, no signal-shaping
term, which is why the measured resolution can sit well above it.

Usage
-----
  python3 decide.py                              # all gases with a CSV
  python3 decide.py --ref ar_co2_iso_93_5_2 --ref-vmesh 410 --drift-dv 200
  python3 decide.py --n-thresh 3                 # threshold at the 3rd cluster
"""
import os
import math
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import matplotlib as _mpl
# Slide typography, matched to the rest of the MPGD26 figure set.
_mpl.rcParams.update({
    'font.size': 13, 'axes.titlesize': 14, 'axes.labelsize': 17,
    'xtick.labelsize': 16, 'ytick.labelsize': 16, 'legend.fontsize': 12,
})


import gases as G
from analyze import load, ln_gain_at, invert_gain, RESULTS, FIGS, COLORS

# Goal set by the P2 timing study: sub-20 ns on the track time stamp.
TIMING_GOAL_NS = 20.0


def _interp(x, xp, fp):
    """np.interp with NaN outside the scanned range (no silent extrapolation)."""
    out = np.interp(x, xp, fp, left=np.nan, right=np.nan)
    return out


def drift_metrics(d, dv_op, n_cl, n_thresh):
    """Timing / transport numbers for gas `d` at operating HV difference dv_op."""
    dV, vd = d["dV"], d["vd"]                     # [V], [cm/us]
    m = {}
    m["vd_op"] = float(_interp(dv_op, dV, vd))    # cm/us at the operating point
    i_best = int(np.nanargmax(vd))
    m["vd_max"] = float(vd[i_best])
    m["dv_best"] = float(dV[i_best])
    # How flat is the plateau? spread of v_d over +-50 V around the operating
    # point, relative to v_d itself -> sensitivity to HV / pressure drift.
    win = (dV >= dv_op - 50) & (dV <= dv_op + 50)
    m["vd_spread_pct"] = (float(np.ptp(vd[win])) / m["vd_op"] * 100.0
                          if win.sum() >= 2 and np.isfinite(m["vd_op"]) else np.nan)

    L = G.D_DRIFT_CM
    vd_cm_ns = m["vd_op"] * 1e-3                  # cm/us -> cm/ns
    m["t_drift"] = L / vd_cm_ns if vd_cm_ns > 0 else np.nan   # full-gap time [ns]

    # longitudinal diffusion at the operating field, in ns
    dl = _interp(dv_op, dV, d["dl"])              # cm^1/2
    m["sig_diff"] = float(dl * np.sqrt(L) / vd_cm_ns) if vd_cm_ns > 0 else np.nan
    # first / k-th cluster arrival jitter
    m["sig_first"] = float(np.sqrt(n_thresh) / (n_cl * vd_cm_ns)) if vd_cm_ns > 0 else np.nan
    m["sig_gas"] = float(np.hypot(m["sig_first"], m["sig_diff"]))

    # attachment in the drift gap -> electrons lost before reaching the mesh
    eta = _interp(dv_op, dV, d["eta_drift"])
    m["eta_op"] = float(eta)
    m["surv"] = float(np.exp(-eta * L)) if np.isfinite(eta) else np.nan
    return m


def sigma_curve(d, n_cl, n_thresh):
    """sigma_t(gas floor) as a function of the drift HV difference, for plotting."""
    vd_cm_ns = d["vd"] * 1e-3
    ok = vd_cm_ns > 0
    sig = np.full(d["dV"].shape, np.nan)
    s_first = np.sqrt(n_thresh) / (n_cl * vd_cm_ns[ok])
    s_diff = d["dl"][ok] * np.sqrt(G.D_DRIFT_CM) / vd_cm_ns[ok]
    sig[ok] = np.hypot(s_first, s_diff)
    return sig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="ar_co2_iso_93_5_2",
                    help="gas currently in the detector (baseline for the table)")
    ap.add_argument("--ref-vmesh", type=float, default=410.0,
                    help="mesh working point in the reference gas [V]")
    ap.add_argument("--drift-dv", type=float, default=200.0,
                    help="mesh->drift HV difference across the 3 mm gap [V]")
    ap.add_argument("--n-thresh", type=int, default=1,
                    help="clusters needed to cross the front-end threshold")
    ap.add_argument("--measured-ns", type=float, default=20.0,
                    help="measured time resolution in the present gas [ns], for context. "
                         "Beam value 2026-08-01 is ~20 ns; the 32 ns of the earlier "
                         "timing study was a different (slower) working point.")
    args = ap.parse_args()

    os.makedirs(FIGS, exist_ok=True)

    keys = [args.ref] + [k for k in G.GASES if k != args.ref]
    data = [load(k) for k in keys]
    data = [d for d in data if d is not None and len(d["dV"]) > 0]
    if not data:
        raise SystemExit("no usable CSVs in results/ -- run ./run_lxplus.sh first.")

    ref = data[0]
    have_ref_gain = len(ref["v_mesh"]) >= 2
    lg_ref = ln_gain_at(ref, args.ref_vmesh) if have_ref_gain else np.nan

    rows = []
    for i, d in enumerate(data):
        n_cl = G.cluster_density(d["key"])
        m = drift_metrics(d, args.drift_dv, n_cl, args.n_thresh)
        m["key"], m["label"], m["n_cl"] = d["key"], d["label"], n_cl
        # efficiency floor: at least n_thresh clusters in the 3 mm gap (Poisson)
        mu = n_cl * G.D_DRIFT_CM
        p = 1.0 - np.exp(-mu) if args.n_thresh <= 1 else 1.0 - sum(
            np.exp(-mu) * mu ** j / math.factorial(j) for j in range(args.n_thresh))
        m["eff_geom"] = p * (m["surv"] if np.isfinite(m["surv"]) else 1.0)
        # gain / HV: equal-gain mesh voltage, and the gain slope at that point
        m["v_equal"] = np.nan
        m["slope_pct_per_V"] = np.nan
        if have_ref_gain and np.isfinite(lg_ref) and len(d["v_mesh"]) >= 2:
            v_eq = float(np.atleast_1d(invert_gain(d, lg_ref))[0])
            m["v_equal"] = v_eq
            if np.isfinite(v_eq):
                dv = 5.0
                lo, hi = ln_gain_at(d, v_eq - dv), ln_gain_at(d, v_eq + dv)
                if np.isfinite(lo) and np.isfinite(hi):
                    m["slope_pct_per_V"] = (hi - lo) / (2 * dv) * 100.0
        m["color"] = COLORS[i % len(COLORS)]
        rows.append(m)

    # ---- Fig: gas-limited timing floor vs drift HV ------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.6, 8.2), sharex=True)
    for d, m in zip(data, rows):
        ax1.plot(d["dV"], d["vd"], "-o", ms=3, color=m["color"], label=m["label"])
        ax2.plot(d["dV"], sigma_curve(d, m["n_cl"], args.n_thresh), "-o", ms=3,
                 color=m["color"], label=m["label"])
    for ax in (ax1, ax2):
        ax.axvline(args.drift_dv, color="0.5", ls=":", lw=1)
        ax.grid(alpha=0.3)
    ax2.axhline(TIMING_GOAL_NS, color="k", ls="--", lw=1.2)
    ax2.text(ax2.get_xlim()[1], TIMING_GOAL_NS, " %.0f ns goal" % TIMING_GOAL_NS,
             va="bottom", ha="right", fontsize=9)
    ax1.set_ylabel(r"drift velocity  $v_d$  [cm/$\mu$s]")
    ax1.set_title("P2 gas comparison: drift speed and gas-limited timing floor\n"
                  "(3 mm gap, threshold at cluster #%d)" % args.n_thresh)
    ax1.legend(fontsize=8)
    ax2.set_yscale("log")
    ax2.set_ylabel(r"$\sigma_t$ gas floor  [ns]")
    ax2.set_xlabel(r"mesh$\to$drift HV difference across 3 mm  [V]")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "timing_floor_vs_driftHV.png"), dpi=140)
    plt.close(fig)

    # ---- decision table ---------------------------------------------------
    L = []
    L.append("=" * 108)
    L.append("P2 Micromegas gas-change decision table  (150 um amp gap, 3 mm drift gap, Magboltz 20 C / 760 Torr)")
    L.append("baseline gas: %s   mesh %.0f V   drift dV %.0f V (E = %.0f V/cm)"
             % (ref["label"], args.ref_vmesh, args.drift_dv, args.drift_dv / G.D_DRIFT_CM))
    L.append("=" * 108)
    L.append("")
    hdr = ("%-24s %7s %7s %7s %8s %8s %8s %7s %8s %7s"
           % ("gas", "n_cl", "v_d", "v_d,max", "t_drift", "sig_gas", "V_equal",
              "dG/dV", "e-surv", "eff"))
    L.append(hdr)
    L.append("%-24s %7s %7s %7s %8s %8s %8s %7s %8s %7s"
             % ("", "[1/cm]", "[cm/us]", "@dV[V]", "[ns]", "[ns]", "[V]",
                "[%/V]", "3mm", "[%]"))
    L.append("-" * 108)

    def f(x, fmt):
        return ("%" + fmt) % x if np.isfinite(x) else "--".rjust(int(fmt.split(".")[0] or 6))

    for m in rows:
        L.append("%-24s %7.1f %7.2f %4.1f@%3.0f %8.1f %8.1f %8s %7s %8.3f %6.1f%s"
                 % (m["label"], m["n_cl"], m["vd_op"], m["vd_max"], m["dv_best"],
                    m["t_drift"], m["sig_gas"],
                    f(m["v_equal"], "8.0f"), f(m["slope_pct_per_V"], "7.2f"),
                    m["surv"], 100 * m["eff_geom"],
                    "   <- baseline" if m["key"] == ref["key"] else ""))
    L.append("-" * 108)
    L.append("")
    L.append("columns")
    L.append("  n_cl     primary ionisation clusters per cm (partial-pressure weighted, PDG/Sauli)")
    L.append("  v_d      drift velocity at the operating drift HV (dV = %.0f V over 3 mm)" % args.drift_dv)
    L.append("  v_d,max  peak drift velocity in the scan, and the drift dV where it peaks")
    L.append("  t_drift  full-gap drift time = 3 mm / v_d  (the time-walk budget)")
    L.append("  sig_gas  GAS-LIMITED time resolution floor = sqrt( (sqrt(k)/(n_cl v_d))^2 + (D_L sqrt(L)/v_d)^2 ), k=%d" % args.n_thresh)
    L.append("  V_equal  mesh voltage giving the SAME gain as the baseline at %.0f V" % args.ref_vmesh)
    L.append("  dG/dV    gain slope at that voltage -- higher = twitchier HV, harder to hold gain stable")
    L.append("  e-surv   fraction of drift electrons surviving attachment over 3 mm")
    L.append("  eff      P(>= %d cluster(s) in 3 mm) x e-surv -- an intrinsic efficiency ceiling, not the measured one" % args.n_thresh)
    L.append("")

    base = next(m for m in rows if m["key"] == ref["key"])
    L.append("timing verdict  (goal %.0f ns; measured %.0f ns in the present gas)"
             % (TIMING_GOAL_NS, args.measured_ns))
    L.append("-" * 108)
    for m in rows:
        if not np.isfinite(m["sig_gas"]):
            continue
        gain_factor = base["sig_gas"] / m["sig_gas"] if m["sig_gas"] > 0 else np.nan
        # Scale the MEASURED resolution by the gas-floor improvement, assuming the
        # non-gas (electronics/shaping) term is unchanged and adds in quadrature.
        other2 = args.measured_ns ** 2 - base["sig_gas"] ** 2
        proj = np.sqrt(max(other2, 0.0) + m["sig_gas"] ** 2)
        L.append("  %-24s floor %5.1f ns  (x%.2f faster)   projected measured -> %5.1f ns  %s"
                 % (m["label"], m["sig_gas"], gain_factor, proj,
                    "OK" if proj < TIMING_GOAL_NS else "still above goal"))
    L.append("")
    L.append("  NB the projection assumes the non-gas contribution (%.1f ns, from electronics /"
             % np.sqrt(max(args.measured_ns ** 2 - base["sig_gas"] ** 2, 0.0)))
    L.append("  shaping / threshold walk) is unchanged by the gas.  If that term dominates, no")
    L.append("  mixture reaches the goal on its own -- see results/gain_summary.txt and the")
    L.append("  P2 timing study for the electronics side.")
    L.append("")

    txt = "\n".join(L)
    with open(os.path.join(RESULTS, "gas_decision_table.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)
    print("figure -> figs/timing_floor_vs_driftHV.png")
    print("table  -> results/gas_decision_table.txt")


if __name__ == "__main__":
    main()
