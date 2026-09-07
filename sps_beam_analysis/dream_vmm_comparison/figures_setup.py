#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figures_setup.py -- the two backing slides the three-slide deck assumes.

The deck says "same chamber, same beam, same working point" and "one
discriminator level, all six chips" and then moves on.  Both are claims about
the RUN, and neither is visible on any of the three slides.  These two answer
them:

  setup_conditions.png   which runs, at what voltages, with which chip
                         settings -- and DREAM's OWN per-pad threshold on the
                         pad plane, the number the VMM's is being compared to
  setup_groups.png       what a ridgeline row actually is: which pads are in
                         it and what they were sorted on

Both are projector-sized and use figures_deck's type scale, so they can be
dropped into the same talk.

WHY DREAM'S THRESHOLD BELONGS ON THE CONDITIONS SLIDE.  The whole argument is
that the VMM's discriminator sits inside the Landau where DREAM's does not, so
the reader is entitled to ask where DREAM's actually is.  It is 27 ADC -- and
it is 27 ADC on every pad, +/- 1, because it is set per channel as 5 sigma of
that channel's own pedestal noise.  That is the control: a per-channel
threshold, calibrated per channel, comes out uniform; the VMM's single global
DAC setting comes out with a 26 % spread once expressed in charge.  The map is
worth drawing precisely because it is boring.
"""

import os
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import figures as F
import figures_slide as S
import figures_deck as D
import threshold_model as M

HERE = os.path.dirname(os.path.abspath(__file__))
FIG, DATA = os.path.join(HERE, "figures"), os.path.join(HERE, "data")

NBAND = D.NGROUP        # the deck's six bands, so the map matches slide 3


def load():
    g, H, bw = M.load()
    n = json.load(open(os.path.join(DATA, "threshold_model_P2_OUT.json")))
    c = json.load(open(os.path.join(DATA, "run_conditions_P2_OUT.json")))
    thr = pd.read_csv(os.path.join(
        DATA, "dream_zs_threshold_eff_nominal_1_P2_OUT.csv"))
    # pad_id IS the DREAM channel_id -- checked against the pad centres, which
    # agree to the last digit printed in both tables
    g = g.merge(thr[["channel_id", "thr_adc", "ped_mean"]],
                left_on="pad_id", right_on="channel_id", how="left")
    assert g["thr_adc"].notna().all(), "a used pad has no DREAM threshold"
    return g, H, bw, n, c


def padnum(ax, x, y, val, c, fmt="{:.0f}", fs=1.0, txt=F.SURFACE):
    """A pad map with the value written inside the pad.  The circles are drawn
    big enough for two digits and no bigger -- at this size the map is a
    lookup table, not a picture of a gradient, so the numbers have to win."""
    ax.scatter(x, y, c=c, s=560, edgecolor=F.SURFACE, linewidth=1.2, zorder=3)
    for xi, yi, v in zip(x, y, val):
        ax.text(xi, yi, fmt.format(v), fontsize=7.4 * fs, color=txt,
                fontweight="bold", ha="center", va="center", zorder=4)
    ax.set_aspect("equal")
    ax.set_xlim(*D.XLIM)
    ax.set_ylim(*D.YLIM)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_color(F.GRID)


# --------------------------------------------------------------------------- #
# The card's two columns.  Fixed x, so the labels and the values each line up
# down the whole slide; a matplotlib table would centre them and lose that.
# Measured, not guessed: the widest label is 0.125 of the figure and the widest
# value 0.355, which is what puts the value column at 0.190 and the map at
# 0.590.
XLAB, XVAL, XMAP = 0.045, 0.190, 0.590


def conditions(g, n, c):
    """Left: the run card.  Right: DREAM's own per-pad threshold, on the pads."""
    d, v, hv = c["dream"], c["vmm"], c["hv"]
    fig = plt.figure(figsize=D.SLIDE)

    rows = [
        ("§", "Both readouts, same chamber", ""),
        ("", "detector", f"{c['station']} · {n['npad']} pads in the fit"),
        ("", "gas", c["gas"]),
        ("", "drift / mesh", f"{hv['p2_drift_v']} V / {hv['p2_mesh_v']} V"
                             "   (identical in both runs)"),
        ("", "reference tracker", f"EIC uRWELL, drift {hv['urwell_drift_v']} V"
                                  f" / resistive {hv['urwell_resistive_v']} V"),
        ("§", "DREAM", ""),
        ("", "run", f"{d['run']} · {d['sub_runs']}"
                    f"  ({d['n_sub_runs']}×{d['sub_run_minutes']} min)"),
        ("", "taken", d["start"]),
        ("", "trigger", d["trigger"]),
        ("", "zero suppression", f"{d['zs_sigma']:.0f}σ of each channel's own "
                                 f"pedestal → {d['thr_med_adc']} ADC"),
        ("", "tracks / efficiency", f"{d['n_probe_tracks'] / 1e6:.2f} M "
                                    f"· {d['eff'] * 100:.1f} %"),
        ("§", "VMM3a", ""),
        ("", "run", f"{v['run']} · {v['sub_run']}"),
        ("", "taken", v["start"]),
        ("", "trigger", v["trigger"]),
        ("", "gain / peaking", f"{v['gain_mv_per_fc']} mV/fC · "
                               f"{v['peaking_ns']} ns · {v['polarity']} "
                               "polarity"),
        ("", "threshold", f"sdt = {v['sdt']} on all six chips, {v['chips']}"),
        ("", "tracks / efficiency", f"{v['n_good_tracks'] / 1e6:.2f} M "
                                    f"· {v['eff'] * 100:.1f} %"),
    ]
    fs = D.NOTE * 0.80
    y, dy = 0.762, 0.0365
    for kind, label, value in rows:
        if kind == "§":
            y -= 0.014
            fig.text(XLAB - 0.012, y, label, fontsize=fs * 1.06, color=F.INK,
                     fontweight="bold", ha="left", va="center")
        else:
            fig.text(XLAB, y, label, fontsize=fs, color=F.INK2,
                     ha="left", va="center")
            fig.text(XVAL, y, value, fontsize=fs, color=F.INK,
                     ha="left", va="center")
        y -= dy

    # -- the map ------------------------------------------------------------- #
    # width chosen so the equal-aspect box exactly fills it: the map is
    # 122 x 119 mm, so 0.520 of the height needs 0.300 of the width.  Any
    # wider and matplotlib centres the box, and the title and the note below
    # no longer line up with the map's own left edge.
    ax = fig.add_axes([XMAP, 0.255, 0.300, 0.520])
    t = g["thr_adc"].to_numpy(float)
    padnum(ax, g["x_d"].to_numpy(), g["y_d"].to_numpy(), t,
           c=[D.AMBER(0.55)] * len(t), fs=D.SC)
    D.scalebar(ax, x0=434.0, y0=166.0)
    ax.set_title("DREAM's own threshold, per pad  [ADC]", loc="left",
                 color=F.INK, fontsize=D.PTITLE * 0.70, fontweight="bold",
                 pad=6 * D.SC)

    fig.text(XMAP, 0.215,
             f"{d['zs_sigma']:.0f}σ of each channel's own pedestal noise, so "
             "it is\ncalibrated per channel — and it comes out flat: "
             f"{t.min():.0f}–{t.max():.0f} ADC\nover the {len(t)} pads.  The "
             "VMM's level, in these same units,\nis "
             f"{n['T']:.0f} ADC — {n['T'] / np.median(t):.0f}× higher, and "
             f"scattering {n['perpad']['T_relrms'] * 100:.0f} % pad to pad.",
             fontsize=D.NOTE * 0.72, color=F.INK2, ha="left", va="top",
             linespacing=1.5)

    D.chrome(fig, "The two runs this comparison is made of",
             "Same chamber, same gas, same voltages; the two readouts were "
             "taken five days apart.",
             "P2 SPS July 2026 · P2_OUT")
    fig.savefig(f"{FIG}/setup_conditions.png", dpi=D.DPI)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def groups(g, H, bw, n):
    """What a ridgeline row is.  Left: the pads, coloured by which band they
    landed in.  Right: the sort itself, with the cuts drawn on it.

    The point of the right panel is that the bands are cuts on a CONTINUUM,
    not clusters: nothing separates band 3 from band 4 except the requirement
    that the bands hold equal numbers of pads."""
    gr = S.groups(g, H, bw, ngroup=NBAND)
    band = np.empty(len(g), int)
    for k, G in enumerate(gr):
        band[G["idx"]] = k
    amp = g["amp_med_d"].to_numpy(float)

    fig = plt.figure(figsize=D.SLIDE)
    # as on the conditions slide, the map axes is sized so the equal-aspect
    # box exactly fills it and nothing gets centred out from under its title
    ax = fig.add_axes([0.035, 0.265, 0.274, 0.475])
    axr = fig.add_axes([0.415, 0.265, 0.555, 0.475])

    # bands run weak -> strong, and so does the ramp: same encoding as slide 2
    col = [D.AMBER(0.10 + 0.86 * b / (NBAND - 1)) for b in band]
    padnum(ax, g["x_d"].to_numpy(), g["y_d"].to_numpy(), band + 1,
           c=col, fs=D.SC)
    D.scalebar(ax, x0=434.0, y0=167.0)
    ax.set_title("which band each pad is in", loc="left", color=F.INK,
                 fontsize=D.PTITLE * 0.70, fontweight="bold", pad=6 * D.SC)

    # -- the sort ------------------------------------------------------------ #
    o = np.argsort(amp)
    axr.scatter(np.arange(len(o)) + 1, amp[o], s=110,
                c=[D.AMBER(0.10 + 0.86 * band[i] / (NBAND - 1)) for i in o],
                edgecolor=F.SURFACE, linewidth=0.9, zorder=3)
    edge = np.cumsum([len(G["idx"]) for G in gr])[:-1]
    for e in edge:
        axr.axvline(e + 0.5, color=F.MUTED, lw=1.0 * D.SC, ls=(0, (4, 3)),
                    zorder=2)
    for k, G in enumerate(gr):
        mid = np.mean([np.where(o == i)[0][0] for i in G["idx"]]) + 1
        axr.text(mid, amp.max() * 1.045, f"{k + 1}", fontsize=D.NOTE,
                 color=F.INK2, fontweight="bold", ha="center", va="bottom")
    axr.set_xlim(0, len(o) + 1)
    axr.set_ylim(0, amp.max() * 1.16)
    axr.set_xlabel("pads, sorted by DREAM pulse height")
    # the axes is only 3.6 in tall; the label with "DREAM" in it is
    # longer than that and overruns both ends.  The x label already says it.
    axr.set_ylabel("pad median pulse height  [ADC]")
    axr.grid(True, lw=0.6, alpha=0.7)
    axr.set_axisbelow(True)
    axr.set_title("the sort, and where the six cuts fall", loc="left",
                  color=F.INK, fontsize=D.PTITLE * 0.70, fontweight="bold",
                  pad=6 * D.SC)

    fig.text(0.035, 0.150,
             "SORTED ON  the pad's median DREAM pulse height on tracked "
             "events — the same quantity the VMM reports, from the readout "
             "not on trial.\nNot efficiency, not position, and not the VMM's "
             "own median: the threshold under test has already biased that "
             f"one.\nCUT INTO  {NBAND} bands of equal pad count "
             f"({', '.join(str(len(G['idx'])) for G in gr)}).  They come out "
             "spatially banded only because the gain is one smooth gradient.",
             fontsize=D.NOTE * 0.75, color=F.INK2, ha="left", va="top",
             linespacing=1.55)

    D.chrome(fig, "A ridgeline row is a gain band, not a pad",
             f"The {n['npad']} pads sorted by pulse height and cut into "
             f"{NBAND} equal groups; each row of slide 3 is one group's "
             "average.",
             "P2 SPS July 2026 · P2_OUT")
    fig.savefig(f"{FIG}/setup_groups.png", dpi=D.DPI)
    plt.close(fig)


def main():
    g, H, bw, n, c = load()
    with plt.rc_context(D.RC):
        conditions(g, n, c)
        groups(g, H, bw, n)
    print("wrote figures/setup_{conditions,groups}.png")


if __name__ == "__main__":
    main()
