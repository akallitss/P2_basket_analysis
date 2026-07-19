#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_pdf.py

Bundle the gas-study products into a single styled PDF:
  page 1  headline + equal-gain HV table + the HV-equivalence plot
  page 2  gas gain vs mesh HV, and effective Townsend vs field
  page 3  drift velocity vs the mesh->drift HV difference (3 mm gap)

Replots from results/*.csv (vector, no blur), reusing analyze.py's loaders, so
the PDF always matches the current Magboltz scans. Run analyze.py first (or just
run this -- it reads the same CSVs).

Usage:
  python3 build_pdf.py [--ref ar_iso_95_5] [--ref-vmesh 415]
          [--working-points 375 400 415 450 475 500] [--out PATH]
"""
import os
import argparse
import datetime

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import gases as G
import analyze as A   # reuse load / ln_gain_at / invert_gain (no side effects on import)

INK = "#1f3a5f"
HERE = os.path.dirname(os.path.abspath(__file__))
COLORS = A.COLORS


def header(fig, title, sub):
    fig.text(0.06, 0.955, title, fontsize=16, weight="bold", color=INK)
    fig.text(0.06, 0.930, sub, fontsize=9.5, color="0.35")
    fig.text(0.94, 0.955, datetime.date.today().isoformat(), fontsize=8,
             color="0.5", ha="right")
    fig.add_artist(plt.Line2D([0.06, 0.94], [0.918, 0.918], color=INK,
                              lw=1.2, transform=fig.transFigure))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="ar_iso_95_5")
    ap.add_argument("--ref-vmesh", type=float, default=415.0)
    ap.add_argument("--working-points", type=float, nargs="*",
                    default=[375, 400, 415, 450, 475, 500])
    ap.add_argument("--out", default=os.path.join(HERE, "gas_hv_equivalence.pdf"))
    args = ap.parse_args()

    keys = [args.ref] + [k for k in G.GASES if k != args.ref]
    data = [d for d in (A.load(k) for k in keys) if d is not None]
    if not data:
        raise SystemExit("no result CSVs in results/. Run ./run_lxplus.sh then analyze.py.")
    ref, cands = data[0], data[1:]

    with PdfPages(args.out) as pdf:
        # ---- page 1: headline + table + HV-equivalence -------------------
        fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
        header(fig, "P2 Micromegas -- gas gain / HV equivalence",
               "Magboltz (Garfield++), 150 um amplification gap, 20 C / 760 Torr")

        # headline
        lg_wp = A.ln_gain_at(ref, args.ref_vmesh)
        msg = "Reference: %s at mesh %.0f V (gain %.0f)." % (
            ref["label"], args.ref_vmesh, np.exp(lg_wp))
        for d in cands:
            v = A.invert_gain(d, lg_wp)
            if np.isfinite(v):
                msg += "\n%s needs mesh %.0f V (%+.1f V) for equal gain." % (
                    d["label"], v, v - args.ref_vmesh)
        fig.text(0.06, 0.885, msg, fontsize=11, color=INK, va="top")

        # equal-gain table
        rows = []
        for wp in sorted(args.working_points):
            lg = A.ln_gain_at(ref, wp)
            if not np.isfinite(lg):
                continue
            row = [f"{wp:.0f}", f"{np.exp(lg):.0f}"]
            for d in cands:
                v = A.invert_gain(d, lg)
                row.append(f"{v:.0f} ({v-wp:+.1f})" if np.isfinite(v) else "--")
            rows.append(row)
        col_labels = [f"{ref['label']}\nmesh [V]", "gain"] + \
                     [f"{d['label']}\nmesh [V] (dV)" for d in cands]
        ax_t = fig.add_axes([0.06, 0.60, 0.88, 0.20]); ax_t.axis("off")
        tbl = ax_t.table(cellText=rows, colLabels=col_labels, loc="center",
                         cellLoc="center")
        tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1, 1.6)
        for (r, c), cell in tbl.get_celld().items():
            cell.set_edgecolor("0.8")
            if r == 0:
                cell.set_facecolor(INK); cell.set_text_props(color="w", weight="bold")
        ax_t.set_title("Mesh HV for equal gas gain", color=INK, fontsize=11,
                       weight="bold", loc="left", pad=8)

        # HV-equivalence curve (candidate HV vs reference HV)
        ax = fig.add_axes([0.12, 0.10, 0.78, 0.42])
        vg = np.linspace(max(ref["v_mesh"].min(), 300),
                         min(ref["v_mesh"].max(), 640), 61)
        lgr = A.ln_gain_at(ref, vg)
        for i, d in enumerate(cands):
            ax.plot(vg, A.invert_gain(d, lgr), "-", lw=2,
                    color=COLORS[(i + 1) % len(COLORS)], label=d["label"])
        ax.plot(vg, vg, "k--", lw=1, label="unity (%s)" % ref["label"])
        ax.axvline(args.ref_vmesh, color="0.5", ls=":", lw=1)
        ax.set_xlabel("reference mesh HV  [V]"); ax.set_ylabel("equal-gain mesh HV  [V]")
        ax.set_title("Equal-gain mesh HV", color=INK, loc="left")
        ax.grid(alpha=0.3); ax.legend(fontsize=9)
        pdf.savefig(fig); plt.close(fig)

        # ---- page 2: gain vs V, Townsend vs E ----------------------------
        fig = plt.figure(figsize=(8.27, 11.69))
        header(fig, "Gas gain and Townsend coefficient",
               "G = exp[(alpha - eta) * d_amp],  d_amp = 150 um,  E = V_mesh / d_amp")
        ax1 = fig.add_axes([0.12, 0.55, 0.78, 0.33])
        for i, d in enumerate(data):
            ax1.plot(d["v_mesh"], d["ln_gain"] / np.log(10), "-o", ms=3,
                     color=COLORS[i % len(COLORS)], label=d["label"])
        ax1.axvline(args.ref_vmesh, color="0.5", ls="--", lw=1)
        ax1.set_xlabel("mesh voltage  [V]"); ax1.set_ylabel(r"$\log_{10}$ gas gain")
        ax1.set_title("Gas gain vs mesh HV", color=INK, loc="left")
        ax1.grid(alpha=0.3); ax1.legend(fontsize=9)
        ax2 = fig.add_axes([0.12, 0.10, 0.78, 0.33])
        for i, d in enumerate(data):
            ax2.plot(d["E_amp"] / 1e3, d["alpha_eff"], "-o", ms=3,
                     color=COLORS[i % len(COLORS)], label=d["label"])
        ax2.set_xlabel("amplification field  [kV/cm]")
        ax2.set_ylabel(r"$\alpha-\eta$  [cm$^{-1}$]")
        ax2.set_title("Effective Townsend coefficient", color=INK, loc="left")
        ax2.grid(alpha=0.3); ax2.legend(fontsize=9)
        pdf.savefig(fig); plt.close(fig)

        # ---- page 3: drift velocity --------------------------------------
        fig = plt.figure(figsize=(8.27, 11.69))
        header(fig, "Drift velocity vs mesh->drift HV difference",
               "3 mm conversion gap,  HV difference dV = E * d_drift")
        ax = fig.add_axes([0.12, 0.30, 0.78, 0.50])
        for i, d in enumerate(data):
            ax.plot(d["dV"], d["vd"], "-o", ms=3, color=COLORS[i % len(COLORS)],
                    label=d["label"])
        ax.set_xlabel(r"mesh$\to$drift HV difference across 3 mm gap  [V]")
        ax.set_ylabel(r"drift velocity  [cm/$\mu$s]")
        ax.set_title("Drift velocity", color=INK, loc="left")
        ax.grid(alpha=0.3); ax.legend(fontsize=10)
        fig.text(0.12, 0.20,
                 "Both gases peak near ~4 cm/us at a few hundred V/cm; the NSW gas\n"
                 "stays faster and flatter at higher drift field (better timing uniformity).",
                 fontsize=9.5, color="0.3", va="top")
        pdf.savefig(fig); plt.close(fig)

    print("wrote", args.out)


if __name__ == "__main__":
    main()
