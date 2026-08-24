#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_bench_schematic.py -- the Saclay cosmic bench tower, slide 10.

Replaces the inline SVG in the draft deck, which had two problems:

  * its viewBox was 420 units wide while the right-hand plane labels started at
    x=368 and ran to ~443, so every "M3 z = ..." caption was clipped at the
    frame edge;
  * it drew the four M3 planes and the two devices under test but not the
    trigger, even though the trigger is what defines the sample -- two
    scintillator counters, one above the top M3 station and one below the
    bottom one, in coincidence.

Both are fixed here by drawing the whole thing from the geometry instead of
from hand-placed pixels, and letting a tight bounding box size the canvas.

Geometry is the measured bench, not a sketch:
  M3 telescope planes   z = 24, 144, 1185, 1302 mm
  devices under test    z = 232 (P1 position) and 712 (P2 position)
  trigger               2-PMT scintillator coincidence, 5 ns each

The scintillators carry NO z label on purpose: the cosmic-bench report records
them as "above the top and below the bottom M3 stations" and never gives a
coordinate, so putting a number on them would be inventing one.

Sources: reports/cosmic_bench_2026-07/p2_cosmic_bench_report.tex (trigger
section), cosmic_bench_analysis/13_timing_waveforms.py (2-PMT coincidence).

    python3 make_bench_schematic.py            # -> figures/3_act2_bench/
    python3 make_bench_schematic.py -o DIR

@author: ak271430 Alexandra Kallitsopoulou
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

# Deck palette, kept so the slide still reads as part of the same deck.
C_M3_FILL, C_M3_EDGE = "#c9d6e8", "#5a7ca8"
C_DUT_FILL, C_DUT_EDGE, C_DUT_TEXT = "#f6d7c4", "#eb6834", "#c1521f"
C_SC_FILL, C_SC_EDGE, C_SC_TEXT = "#e7dff2", "#7a5ea8", "#5b4483"
C_MU = "#b8860b"
C_INK, C_GREY = "#0b0b0b", "#52514e"

M3_Z = [24.0, 144.0, 1185.0, 1302.0]
# Deck's own names for the two tower positions holding a device under test.
DUT = [(232.0, "P1 plane"), (712.0, "P2 plane")]

# Where the scintillators sit on the drawing. They are outside the M3 span by
# construction; the offset is a drawing choice, which is exactly why they are
# labelled without a z.
SC_GAP = 118.0          # mm above/below the outer M3 planes (drawing choice)
Z_SCALE = 0.36          # drawing units per mm, as in the deck original

BASE = os.path.dirname(os.path.abspath(__file__))
SACLAY = os.path.abspath(os.path.join(BASE, "..", ".."))

# Write to the canonical bench-figure tree, beside the other conference figures
# from 22_conference_figures.py, and let conference/gather_figures.py copy it
# into figures/ like everything else. conference/figures/ is regenerated output
# that gets rebuilt from the manifest -- an unregistered figure written straight
# into it does not survive that rebuild, which is exactly what happened to the
# first version of this one.
OUT_DEFAULT = os.path.join(SACLAY, "data", "Cosmic_Bench", "Analysis",
                           "conference")
# No slide prefix: the manifest owns the slide number and adds `s10_`.
STEM = "cosmic_bench_schematic"


def draw(ax):
    """Everything is placed from the geometry in one coordinate system.

    z is compressed by Z_SCALE for the drawing, exactly as the deck's original
    SVG did: 1.3 m of tower against a 0.3 m plane is a 4:1 sliver at true
    scale, unreadable on a slide. Every plane still carries its real z, so the
    compression costs no information.
    """
    zs = lambda z: z * Z_SCALE

    x0, x1 = 0.0, 300.0            # M3 planes
    dx0, dx1 = 22.0, 278.0         # devices under test, visibly narrower
    sx0, sx1 = -16.0, 316.0        # scintillators, widest of all
    x_lab = 344.0                  # one left-aligned label column
    x_and = -104.0                 # trigger logic column
    x_axis = -178.0                # z axis, outboard of the logic

    y_top, y_bot = zs(M3_Z[-1] + SC_GAP), zs(M3_Z[0] - SC_GAP)

    # ---- cosmic muon, behind everything ---------------------------------
    ax.annotate("", xy=(212, y_bot - 30), xytext=(150, y_top + 30),
                arrowprops=dict(arrowstyle="-|>", color=C_MU, lw=2.0,
                                linestyle=(0, (6, 4)), shrinkA=0, shrinkB=0),
                zorder=1)
    ax.text(146, y_top + 34, "cosmic µ", color=C_MU, fontsize=12.5,
            fontweight="bold", ha="right", va="bottom", zorder=5)

    # ---- trigger scintillators ------------------------------------------
    for z, tag in ((M3_Z[-1] + SC_GAP, "top"), (M3_Z[0] - SC_GAP, "bottom")):
        y, h, pw = zs(z), 11.0, 26.0
        ax.add_patch(Rectangle((sx0, y - h), sx1 - sx0, 2 * h,
                               facecolor=C_SC_FILL, edgecolor=C_SC_EDGE,
                               lw=1.8, zorder=3))
        # PMTs sit OUTSIDE the bar, so they never sit under its label
        for xp in (sx0 - pw, sx1):
            ax.add_patch(FancyBboxPatch((xp, y - 8), pw, 16,
                                        boxstyle="round,pad=0,rounding_size=4",
                                        facecolor=C_SC_EDGE, edgecolor="none",
                                        alpha=0.9, zorder=3))
        ax.text((sx0 + sx1) / 2, y, "scintillator", ha="center", va="center",
                fontsize=11.5, color=C_SC_TEXT, fontweight="bold", zorder=4)
        ax.text(x_lab, y, f"{tag} trigger counter\nplastic + 2 PMT",
                ha="left", va="center", fontsize=11.5, color=C_SC_TEXT,
                linespacing=1.35)

    # ---- M3 telescope planes --------------------------------------------
    for z in M3_Z:
        ax.add_patch(Rectangle((x0, zs(z) - 5), x1 - x0, 10,
                               facecolor=C_M3_FILL, edgecolor=C_M3_EDGE,
                               lw=1.4, zorder=3))
        ax.text(x_lab, zs(z), f"M3    z = {z:.0f}", ha="left", va="center",
                fontsize=12, color=C_GREY)

    # ---- devices under test ---------------------------------------------
    for z, name in DUT:
        ax.add_patch(Rectangle((dx0, zs(z) - 9), dx1 - dx0, 18,
                               facecolor=C_DUT_FILL, edgecolor=C_DUT_EDGE,
                               lw=2.2, zorder=3))
        # One line, not two: z = 232 and the M3 at z = 144 are only 88 mm
        # apart, and a two-line caption on the first overlaps the second.
        ax.text(x_lab, zs(z), f"{name}  —  P2 under test,  z = {z:.0f}",
                ha="left", va="center", fontsize=12, color=C_DUT_TEXT,
                fontweight="bold")

    # ---- the coincidence that defines the trigger ------------------------
    y_and = (y_top + y_bot) / 2
    for y in (y_top, y_bot):
        ax.plot([sx0 - 26, x_and], [y, y], color=C_SC_EDGE, lw=1.5, zorder=2)
    ax.plot([x_and, x_and], [y_bot, y_and - 22], color=C_SC_EDGE, lw=1.5, zorder=2)
    ax.plot([x_and, x_and], [y_and + 22, y_top], color=C_SC_EDGE, lw=1.5, zorder=2)
    ax.add_patch(FancyBboxPatch((x_and - 40, y_and - 22), 80, 44,
                                boxstyle="round,pad=0,rounding_size=8",
                                facecolor="white", edgecolor=C_SC_EDGE,
                                lw=1.8, zorder=4))
    ax.text(x_and, y_and, "AND", ha="center", va="center", fontsize=13.5,
            fontweight="bold", color=C_SC_TEXT, zorder=5)
    ax.annotate("", xy=(x_and, y_bot - 46), xytext=(x_and, y_and - 22),
                arrowprops=dict(arrowstyle="-|>", color=C_SC_EDGE, lw=1.6))
    ax.text(x_and, y_bot - 54, "trigger to all FEUs\n2-PMT coincidence, 5 ns each",
            ha="center", va="top", fontsize=11, color=C_SC_TEXT,
            linespacing=1.35)

    # ---- z axis ----------------------------------------------------------
    ax.annotate("", xy=(x_axis, zs(M3_Z[-1]) + 16), xytext=(x_axis, zs(M3_Z[0]) - 16),
                arrowprops=dict(arrowstyle="-|>", color=C_GREY, lw=1.4))
    ax.text(x_axis - 17, zs((M3_Z[0] + M3_Z[-1]) / 2), "z  [mm]", rotation=90,
            ha="center", va="center", fontsize=12, color=C_GREY)

    ax.set_title("The M3 telescope brackets the devices under test;\n"
                 "the trigger is a top–bottom scintillator coincidence",
                 fontsize=13.5, color=C_INK, fontweight="bold", pad=14)

    ax.set_xlim(x_axis - 34, x_lab + 210)
    ax.set_ylim(y_bot - 96, y_top + 62)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out-dir", default=OUT_DEFAULT)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(11.0, 8.6))
    draw(ax)
    for ext in ("svg", "pdf", "png"):
        p = os.path.join(args.out_dir, f"{STEM}.{ext}")
        # tight bbox with a real pad is what guarantees nothing is clipped --
        # the bug this figure exists to fix.
        fig.savefig(p, dpi=200, bbox_inches="tight", pad_inches=0.28,
                    facecolor="white")
        print(f"  {p}")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
