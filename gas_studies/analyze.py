#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze.py -- turn the Magboltz scans (results/<gas>.csv, produced on lxplus by
run_lxplus.sh) into the two deliverables:

  1. GAIN EQUIVALENCE between gases:  at what mesh HV does a candidate gas give
     the same gas gain as the reference gas?  We use the parallel-plate gain
        G(V_mesh) = exp[ (alpha - eta)(E) * d_amp ],   E = V_mesh / d_amp
     and, for every reference mesh voltage, invert G_candidate to find the mesh
     voltage that reproduces the same gain.

  2. DRIFT VELOCITY vs the mesh<->drift HV difference across the 3 mm gap:
        dV = E * d_drift,   v_d(E)   for each gas overlaid.

Outputs (figs/ and results/):
  figs/gain_vs_vmesh.png            G(V_mesh) both gases, working point marked
  figs/hv_equivalence.png           candidate mesh HV vs reference mesh HV (+dV)
  figs/drift_velocity_vs_dV.png     v_d vs mesh-drift HV difference
  figs/townsend_vs_E.png            (alpha-eta)(E) diagnostic
  results/hv_equivalence.csv        the mapping table (ref V -> candidate V)
  results/gain_summary.txt          headline numbers at the working point(s)

Usage:
  python3 analyze.py [--ref-vmesh 415] [--ref ar_iso_95_5]
                     [--working-points 375 400 415 450]
"""
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import gases as G

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
FIGS = os.path.join(HERE, "figs")

# consistent per-gas colours
COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf"]


def load(key):
    """Load a scan CSV -> dict with amp/drift arrays. Returns None if missing."""
    path = os.path.join(RESULTS, key + ".csv")
    if not os.path.exists(path):
        return None
    rows = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    reg = np.array([str(r) for r in rows["region"]])
    amp = rows[reg == "amp"]
    drift = rows[reg == "drift"]
    # amplification: mesh voltage and effective Townsend
    v_mesh = amp["E_Vcm"] * G.D_AMP_CM
    alpha_eff = amp["alpha_cm"] - amp["eta_cm"]          # 1/cm
    ln_gain = alpha_eff * G.D_AMP_CM                     # ln(G) = (alpha-eta) d
    order = np.argsort(v_mesh)
    # drift: HV difference across the 3 mm gap and drift velocity
    dV = drift["E_Vcm"] * G.D_DRIFT_CM
    vd = drift["vz_cm_ns"] * 1.0e3                       # cm/ns -> cm/us
    dorder = np.argsort(dV)
    return {
        "key": key,
        "label": G.GASES[key]["label"] if key in G.GASES else key,
        "role": G.GASES[key].get("role", "") if key in G.GASES else "",
        "v_mesh": v_mesh[order],
        "ln_gain": ln_gain[order],
        "alpha_eff": alpha_eff[order],
        "E_amp": amp["E_Vcm"][order],
        "dV": dV[dorder],
        "vd": vd[dorder],
        "E_drift": drift["E_Vcm"][dorder],
    }


def ln_gain_at(gas, v_mesh):
    """ln(gain) at an arbitrary mesh voltage (interp in V; extrapolate flagged)."""
    return np.interp(v_mesh, gas["v_mesh"], gas["ln_gain"],
                     left=np.nan, right=np.nan)


def invert_gain(gas, target_ln_gain):
    """Mesh voltage in `gas` that reproduces target ln(gain). NaN if out of range.

    ln(gain) rises monotonically with V over the amplification scan; we sort by
    ln(gain) before interpolating so np.interp's increasing-xp rule holds even if
    Magboltz statistics make two neighbouring points wobble.
    """
    x, y = gas["v_mesh"], gas["ln_gain"]
    order = np.argsort(y)
    ys, xs = y[order], x[order]
    tg = np.atleast_1d(target_ln_gain).astype(float)
    out = np.full(tg.shape, np.nan)
    inside = (tg >= ys.min()) & (tg <= ys.max())
    out[inside] = np.interp(tg[inside], ys, xs)
    return out if out.size > 1 else float(out[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="ar_iso_95_5", help="reference gas key")
    ap.add_argument("--ref-vmesh", type=float, default=415.0,
                    help="reference working-point mesh voltage [V]")
    ap.add_argument("--working-points", type=float, nargs="*",
                    default=[375, 400, 415, 450, 475, 500],
                    help="reference mesh voltages to tabulate")
    args = ap.parse_args()

    os.makedirs(FIGS, exist_ok=True)

    # load every gas that has a CSV (reference first)
    keys = [args.ref] + [k for k in G.GASES if k != args.ref]
    data = [load(k) for k in keys]
    data = [d for d in data if d is not None]
    if not data:
        raise SystemExit("no result CSVs found in results/. Run ./run_lxplus.sh first.")
    ref = data[0]
    if ref["key"] != args.ref:
        print("WARNING: reference gas %s has no CSV; using %s" % (args.ref, ref["key"]))
    cands = data[1:]

    # ---- Fig 1: gain vs mesh voltage --------------------------------------
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for i, d in enumerate(data):
        ax.plot(d["v_mesh"], d["ln_gain"] / np.log(10.0), "-o", ms=3,
                color=COLORS[i % len(COLORS)],
                label="%s (%s)" % (d["label"], d["role"]))
    ax.axvline(args.ref_vmesh, color="0.5", ls="--", lw=1)
    ax.set_xlabel("mesh voltage across 150 um gap  [V]")
    ax.set_ylabel(r"$\log_{10}$ gas gain   $= (\alpha-\eta)\,d_{amp}/\ln 10$")
    ax.set_title("P2 Micromegas gas gain vs mesh HV (Magboltz, 150 um gap)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "gain_vs_vmesh.png"), dpi=140)
    plt.close(fig)

    # ---- Fig 2: HV equivalence mapping ------------------------------------
    # over the reference voltage grid, find candidate V giving the same gain
    v_ref_grid = np.linspace(max(275, ref["v_mesh"].min()),
                             min(650, ref["v_mesh"].max()), 61)
    lg_ref = ln_gain_at(ref, v_ref_grid)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 7.6), sharex=True,
                                   gridspec_kw={"height_ratios": [2, 1]})
    mapping = {"V_ref": v_ref_grid}
    for i, d in enumerate(cands):
        v_equiv = invert_gain(d, lg_ref)
        mapping["V_" + d["key"]] = v_equiv
        c = COLORS[(i + 1) % len(COLORS)]
        ax1.plot(v_ref_grid, v_equiv, "-", color=c, lw=2, label=d["label"])
        ax2.plot(v_ref_grid, v_equiv - v_ref_grid, "-", color=c, lw=2)
    ax1.plot(v_ref_grid, v_ref_grid, "k--", lw=1, label="unity (%s)" % ref["label"])
    ax1.axvline(args.ref_vmesh, color="0.5", ls=":", lw=1)
    ax2.axvline(args.ref_vmesh, color="0.5", ls=":", lw=1)
    ax2.axhline(0, color="k", lw=0.8)
    ax1.set_ylabel("candidate mesh HV for equal gain  [V]")
    ax1.set_title("Equal-gain mesh HV:  candidate gas vs %s" % ref["label"])
    ax1.grid(alpha=0.3); ax1.legend()
    ax2.set_ylabel(r"$\Delta V_{mesh}$  [V]")
    ax2.set_xlabel("reference (%s) mesh HV  [V]" % ref["label"])
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "hv_equivalence.png"), dpi=140)
    plt.close(fig)

    # ---- Fig 3: drift velocity vs mesh-drift HV difference ----------------
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for i, d in enumerate(data):
        ax.plot(d["dV"], d["vd"], "-o", ms=3, color=COLORS[i % len(COLORS)],
                label=d["label"])
    ax.set_xlabel("mesh$\\to$drift HV difference across 3 mm gap  [V]")
    ax.set_ylabel(r"drift velocity  $v_d$  [cm/$\mu$s]")
    ax.set_title("P2 drift velocity vs drift HV (Magboltz, 3 mm gap)")
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "drift_velocity_vs_dV.png"), dpi=140)
    plt.close(fig)

    # ---- Fig 4: Townsend diagnostic ---------------------------------------
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for i, d in enumerate(data):
        ax.plot(d["E_amp"] / 1e3, d["alpha_eff"], "-o", ms=3,
                color=COLORS[i % len(COLORS)], label=d["label"])
    ax.set_xlabel("amplification field  E  [kV/cm]")
    ax.set_ylabel(r"effective Townsend  $\alpha-\eta$  [cm$^{-1}$]")
    ax.set_title("Effective Townsend coefficient vs field")
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "townsend_vs_E.png"), dpi=140)
    plt.close(fig)

    # ---- write mapping CSV ------------------------------------------------
    cols = list(mapping.keys())
    arr = np.column_stack([mapping[c] for c in cols])
    hdr = ",".join(cols)
    np.savetxt(os.path.join(RESULTS, "hv_equivalence.csv"), arr, delimiter=",",
               header=hdr, comments="", fmt="%.3f")

    # ---- headline summary -------------------------------------------------
    lines = []
    lines.append("P2 gas gain-equivalence summary (Magboltz, 150 um amp gap)")
    lines.append("reference gas: %s  (working point mesh = %.0f V)"
                 % (ref["label"], args.ref_vmesh))
    lines.append("gain model: G = exp[(alpha-eta) * d_amp], d_amp = %.4f cm" % G.D_AMP_CM)
    lines.append("")
    for wp in sorted(args.working_points):
        lg = ln_gain_at(ref, wp)
        if not np.isfinite(lg):
            continue
        lines.append("--- reference mesh %.0f V : gain = %.3g (log10 G = %.2f) ---"
                     % (wp, np.exp(lg), lg / np.log(10)))
        for d in cands:
            v_eq = invert_gain(d, lg)
            if np.isfinite(v_eq):
                lines.append("    %-24s equal gain at mesh = %6.1f V   (%+.1f V)"
                             % (d["label"], v_eq, v_eq - wp))
            else:
                lines.append("    %-24s equal gain OUTSIDE scanned range" % d["label"])
        lines.append("")

    txt = "\n".join(lines)
    with open(os.path.join(RESULTS, "gain_summary.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)
    print("figures -> figs/ , tables -> results/")


if __name__ == "__main__":
    main()
