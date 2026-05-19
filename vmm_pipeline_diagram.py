#!/usr/bin/env python3
"""
Generate a pipeline diagram for the VMM trigger rate analysis.
Run standalone: python vmm_pipeline_diagram.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ── Colours ────────────────────────────────────────────────────────────────
C_INPUT   = "#D0E8F2"   # light blue   – data source
C_PROC    = "#D4EAD0"   # light green  – computation / algorithm
C_QC      = "#FFF3CC"   # light yellow – quality selection
C_OUT     = "#FDDCB5"   # light orange – plot / file output
C_MASK    = "#E8D5F0"   # light purple – derived object (good channel mask)
C_ARROW   = "#444444"
C_PARAM   = "#555555"
C_BORDER  = "#888888"


def box(ax, x, y, w, h, label, sublabel=None, color=C_PROC, fontsize=9.5):
    """Draw a rounded rectangle with a bold label and optional sublabel."""
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.07",
        facecolor=color, edgecolor=C_BORDER, linewidth=1.2, zorder=3
    )
    ax.add_patch(patch)
    if sublabel:
        ax.text(x, y + 0.09, label, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", zorder=4,
                linespacing=1.3)
        ax.text(x, y - 0.20, sublabel, ha="center", va="center",
                fontsize=fontsize - 1.2, color=C_PARAM, zorder=4,
                linespacing=1.3)
    else:
        ax.text(x, y, label, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", zorder=4,
                linespacing=1.3)


def arrow(ax, x1, y1, x2, y2, label=None, label_x=None, label_y=None,
          color=C_ARROW, lw=1.4, style="->"):
    """Draw an annotated arrow."""
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color,
                                lw=lw, connectionstyle="arc3,rad=0.0"),
                zorder=2)
    if label:
        lx = label_x if label_x is not None else (x1 + x2) / 2 + 0.12
        ly = label_y if label_y is not None else (y1 + y2) / 2
        ax.text(lx, ly, label, ha="left", va="center",
                fontsize=8, color=C_PARAM, style="italic", zorder=5)


def out_box(ax, x, y, w, h, text, color=C_OUT):
    """Small output file / plot box."""
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.05",
        facecolor=color, edgecolor=C_BORDER, linewidth=0.9, zorder=3,
        linestyle="--"
    )
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center",
            fontsize=7.8, zorder=4, linespacing=1.4)


def param_note(ax, x, y, text):
    ax.text(x, y, text, ha="left", va="center",
            fontsize=7.6, color=C_PARAM, zorder=4,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.0))


# ── Canvas ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 13))
ax.set_xlim(0, 10)
ax.set_ylim(0, 13)
ax.axis("off")

fig.patch.set_facecolor("white")

# ── Layout constants ───────────────────────────────────────────────────────
XL  = 5.0    # centre x of pipeline boxes
WL  = 8.8    # width of main boxes
HL  = 0.72   # height of main boxes

# y positions (top → bottom)
Y = {
    "input" : 12.1,
    "s1"    : 10.8,
    "s2"    : 9.35,
    "s3a"   : 7.90,
    "s3b"   : 6.45,
    "s3c"   : 5.00,
    "s4"    : 3.55,
    "s5"    : 2.10,
}

# ── Main pipeline boxes ────────────────────────────────────────────────────
box(ax, XL, Y["input"], WL, HL,
    "ROOT files  (sng=0 signal runs)",
    sublabel="VMM hit data: time, vmm_id, channel  |  one file per run loaded",
    color=C_INPUT)

box(ax, XL, Y["s1"], WL, HL,
    "Step 1 — Trigger rate time series",
    sublabel="Bin VMM 0 ch 40 timestamps into 1 ms windows  →  rate vs time",
    color=C_PROC)

box(ax, XL, Y["s2"], WL, HL,
    "Step 2 — Spill mask construction",
    sublabel="threshold = 1 kHz  |  hole-fill ≤ 2 s  |  min spill = 1 s",
    color=C_PROC)

box(ax, XL, Y["s3a"], WL, HL,
    "Step 3a — Channel hit counts  (dead / noisy identification)",
    sublabel="dead < 10 % × median  |  noisy (diagnostic) > 5 × median  |  3 files, skip file 0",
    color=C_QC)

box(ax, XL, Y["s3b"], WL, HL,
    "Step 3b — Per-channel spill-on / spill-off rates",
    sublabel="noisy if off-spill rate > 3 × median  |  max across all diagnostic runs  |  5 files",
    color=C_QC)

box(ax, XL, Y["s3c"], WL, HL,
    "Step 3c — Good channel mask",
    sublabel="good = connected  −  dead  −  noisy off-spill  (applied to all configs)",
    color=C_MASK)

box(ax, XL, Y["s4"], WL, HL,
    "Step 4 — Spill-on / spill-off rate per configuration",
    sublabel="rate_on, rate_off, ±1σ across channels  |  normalised by n_good_channels",
    color=C_PROC)

box(ax, XL, Y["s5"], WL, HL,
    "Step 5 — Summary plots & CSV output",
    sublabel="avg rate per channel vs (sg, snt)  |  ±1σ band  |  trigger reference overlay",
    color=C_PROC)

# ── Vertical arrows (main flow) ────────────────────────────────────────────
for ya, yb in [
    (Y["input"] - HL/2, Y["s1"]  + HL/2),
    (Y["s1"]   - HL/2, Y["s2"]  + HL/2),
    (Y["s2"]   - HL/2, Y["s3a"] + HL/2),
    (Y["s3a"]  - HL/2, Y["s3b"] + HL/2),
    (Y["s3b"]  - HL/2, Y["s3c"] + HL/2),
    (Y["s3c"]  - HL/2, Y["s4"]  + HL/2),
    (Y["s4"]   - HL/2, Y["s5"]  + HL/2),
]:
    arrow(ax, XL, ya, XL, yb)

# ── Feedback arrow: mask applied in Step 4 ────────────────────────────────
# curved arrow from s3c left side down to s4 left side
ax.annotate("", xy=(XL - WL/2, Y["s4"]), xytext=(XL - WL/2, Y["s3c"] - HL/2),
            arrowprops=dict(arrowstyle="->", color="#9966CC", lw=1.6,
                            connectionstyle="arc3,rad=-0.5"),
            zorder=2)
ax.text(XL - WL/2 - 0.22, (Y["s3c"] + Y["s4"]) / 2,
        "good channel\nmask applied", ha="right", va="center",
        fontsize=7.8, color="#9966CC", style="italic")

# ── Legend ─────────────────────────────────────────────────────────────────
legend_items = [
    mpatches.Patch(facecolor=C_INPUT,  edgecolor=C_BORDER, label="Data source"),
    mpatches.Patch(facecolor=C_PROC,   edgecolor=C_BORDER, label="Computation / algorithm"),
    mpatches.Patch(facecolor=C_QC,     edgecolor=C_BORDER, label="Quality selection"),
    mpatches.Patch(facecolor=C_MASK,   edgecolor=C_BORDER, label="Derived object (mask)"),
    mpatches.Patch(facecolor=C_OUT,    edgecolor=C_BORDER,
                   label="Output (plot / CSV)", linestyle="--"),
]
ax.legend(handles=legend_items, loc="lower left",
          bbox_to_anchor=(0.01, 0.01), fontsize=8.5, framealpha=0.9,
          edgecolor=C_BORDER)

# ── Title ──────────────────────────────────────────────────────────────────
ax.text(5.0, 12.82,
        "VMM Trigger Rate Analysis — Pipeline Overview",
        ha="center", va="center", fontsize=13, fontweight="bold")

ax.text(5.0, 12.55,
        "P2 Basket Detector  |  SPS Beam Test",
        ha="center", va="center", fontsize=9.5, color="#555555")

plt.tight_layout(pad=0.3)

out_dir = "/drf/projets/clas12/P2/akallits/plots_trigger"
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(out_dir, f"analysis_pipeline.{ext}"),
                bbox_inches="tight", dpi=180)
print("Saved analysis_pipeline.pdf / .png")
plt.show()
plt.close(fig)