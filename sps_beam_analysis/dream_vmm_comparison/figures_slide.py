#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figures_slide.py -- the MPGD2026 slide figures for P2_OUT: one conclusion,
readable in a minute.

The conclusion: the pad-to-pad gain variation belongs to the DETECTOR -- both
readouts measure the same map -- and the whole difference between them is that
the VMM's discriminator cuts into the Landau, so on the low-gain pads it loses
the small pulses DREAM still records.

Four figures, three of them stand-alone slides and one composite:
  slide_ridge.png  the mechanism, in one panel: the Landau slides across a
                   fixed threshold line as the pad gain falls
  slide_proof.png  that this is quantitatively the whole story
  slide_fix.png    what it costs to fix
  slide_full.png   all three on one 16:9 slide

Convention carried from the rest of the deck: DREAM = orange, VMM = blue.  In
the ridge panel saturation carries the second meaning -- faint orange is what
DREAM records, solid orange is the part of it below the VMM's threshold.
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import figures as F
import threshold_model as M

HERE = os.path.dirname(os.path.abspath(__file__))
FIG, DATA = os.path.join(HERE, "figures"), os.path.join(HERE, "data")

VMM_C, DREAM_C = F.C1, F.C2
NGROUP = 8


def _style(ax):
    ax.grid(True, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)


def load():
    g, H, bw = M.load()
    S = M.Spectra(H, bw)
    n = json.load(open(os.path.join(DATA, "threshold_model_P2_OUT.json")))
    return g, H, bw, S, n


def groups(g, H, bw, ngroup=None):
    """Pads sorted by gain, cut into `ngroup` bands (NGROUP by default).  Each band's spectrum is the
    track-weighted mean of its pads' UNIT-AREA spectra, not the pooled raw
    counts -- that way the band's sub-threshold fraction is exactly the mean of
    its pads', which is the number the efficiency model uses."""
    order = np.argsort(g["amp_med_d"].to_numpy())
    out = []
    for idx in np.array_split(order, ngroup or NGROUP):
        w = g["n_track_v"].to_numpy()[idx].astype(float)
        p = H[idx] / H[idx].sum(1, keepdims=True)
        out.append(dict(
            idx=idx, h=np.average(p, axis=0, weights=w),
            gain=float(np.average(g["amp_med_d"].to_numpy()[idx], weights=w)
                       / np.median(g["amp_med_d"])),
            eff_v=float(np.average(g["eff_v"].to_numpy()[idx], weights=w)),
            eff_d=float(np.average(g["eff_d"].to_numpy()[idx],
                                   weights=g["n_track_d"].to_numpy()[idx]))))
    return out


# --------------------------------------------------------------------------- #
def ridge(ax, gr, bw, T, xlim=(28.0, 3000.0), h=1.02, fs=1.0):
    """The hero panel.  One filled curve per gain band, stacked bottom (weak)
    to top (strong), with the single discriminator level drawn straight through
    all of them.

    The x axis is LOGARITHMIC on purpose: a gain difference is a multiplicative
    factor, so on a log axis the eight bands are the same curve TRANSLATED, and
    what changes from row to row is only how much of it falls left of the fixed
    threshold line.  (Plotted as A x dN/dA, the density per decade, so the areas
    on screen are the real fractions.)"""
    c = (np.arange(gr[0]["h"].size) + 0.5) * bw
    # 3600 up: the DREAM amplitude saturates and the last bins are an overflow
    # pile-up.  Left in, that spike sets the normalisation and flattens the
    # high-gain rows into lines.
    k = (c >= xlim[0]) & (c <= min(xlim[1] * 1.35, 3600.0))
    for j, G in enumerate(gr):
        y = np.convolve(G["h"], np.ones(3) / 3, "same")[k] * c[k]
        y = y / y.max() * h
        base = float(j)
        ax.fill_between(c[k], base, base + y, color=DREAM_C, alpha=0.15,
                        lw=0, zorder=3 + j)
        lost = c[k] <= T
        ax.fill_between(c[k][lost], base, (base + y)[lost], color=DREAM_C,
                        alpha=0.80, lw=0, zorder=3 + j)
        ax.plot(c[k], base + y, lw=1.6, color=DREAM_C, zorder=3.5 + j)
        ax.plot(xlim, [base, base], lw=0.7, color=F.GRID, zorder=1)
        ax.text(xlim[0] * 0.90, base + 0.30, f"{G['gain']:.2f}x",
                fontsize=9 * fs, color=F.INK2, ha="right", va="center")

    ax.axvline(T, lw=2.2, color=VMM_C, zorder=40)
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.42, len(gr) + 0.35)
    ax.set_yticks([])
    ax.set_xticks([50, 100, 200, 500, 1000, 2000])
    ax.set_xticklabels(["50", "100", "200", "500", "1000", "2000"])
    ax.set_xlabel("pulse height on the pad, as DREAM records it  [ADC]")
    ax.spines["left"].set_visible(False)
    ax.grid(True, axis="x", lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    return 1.0


def effbars(ax, gr, n, fs=1.0, bh=0.26):
    """The consequence, on the same rows: what each band's pads actually
    recorded."""
    for j, G in enumerate(gr):
        # pair centred at j+0.45 with a 0.04 surface gap, so bh=0.26 lands on
        # the 0.30 / 0.60 rows this figure has always used
        d = (bh + 0.04) / 2
        yd, yv = j + 0.45 - d, j + 0.45 + d
        ax.barh(yd, G["eff_d"], height=bh, color=DREAM_C, alpha=0.85, zorder=3)
        ax.barh(yv, G["eff_v"], height=bh, color=VMM_C, zorder=3)
        ax.text(G["eff_v"] - 0.010, yv, f"{G['eff_v'] * 100:.0f}",
                ha="right", va="center", fontsize=8.5 * fs, color=F.SURFACE,
                fontweight="bold", zorder=4)
        ax.text(G["eff_d"] - 0.010, yd, f"{G['eff_d'] * 100:.0f}",
                ha="right", va="center", fontsize=8.5 * fs, color=F.SURFACE,
                fontweight="bold", zorder=4)
    ax.set_xlim(0.40, 1.0)
    ax.set_ylim(-0.42, len(gr) + 0.35)
    ax.set_yticks([])
    ax.set_xticks([0.5, 0.75, 1.0])
    ax.set_xticklabels(["50 %", "75 %", "100 %"])
    ax.spines["left"].set_visible(False)
    ax.grid(True, axis="x", lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)


# --------------------------------------------------------------------------- #
def fig_ridge(g, H, bw, S, n):
    gr = groups(g, H, bw)
    T = n["T"]
    fig, (ax, axb) = plt.subplots(
        1, 2, figsize=(12.2, 6.4), sharey=True,
        gridspec_kw=dict(width_ratios=[2.7, 1], wspace=0.05))
    ridge(ax, gr, bw, T)
    effbars(axb, gr, n)

    ax.annotate("the VMM's discriminator — one level, all six chips",
                xy=(T, NGROUP + 0.02), xytext=(T * 1.28, NGROUP + 0.26),
                fontsize=9.5, color=VMM_C, fontweight="bold", va="center",
                arrowprops=dict(arrowstyle="-", color=VMM_C, lw=1.0))
    ax.text(T * 0.52, -0.30, "lost", fontsize=10.5, color=DREAM_C,
            fontweight="bold", ha="center", va="center")
    ax.text(T * 2.1, -0.30, "recorded", fontsize=10.5, color=F.INK2,
            ha="center", va="center")
    ax.annotate("", xy=(T * 0.99, -0.30), xytext=(T * 0.72, -0.30),
                arrowprops=dict(arrowstyle="<-", color=DREAM_C, lw=1.1))
    ax.text(28 * 0.90, NGROUP + 0.10, "pad gain", fontsize=9, color=F.INK2,
            ha="right", va="bottom", style="italic")
    axb.text(0.405, NGROUP + 0.06, "of the tracks that point at the pad,\n"
             "the % each readout records",
             fontsize=8.5, color=F.INK2, ha="left", va="bottom")
    axb.text(0.50, -0.30, "VMM", fontsize=9.5, color=VMM_C,
             fontweight="bold", ha="center", va="center")
    axb.text(0.72, -0.30, "DREAM", fontsize=9.5, color=DREAM_C,
             fontweight="bold", ha="center", va="center")

    fig.suptitle("On the weak pads the Landau slides down onto the VMM's "
                 "threshold — and that is the whole difference",
                 x=0.006, ha="left", fontsize=12.5, fontweight="bold",
                 color=F.INK)
    ax.set_title("P2_OUT · 53 pads under the beam, sorted by gain into 8 bands "
                 "· log axis, so a gain factor is a sideways shift",
                 loc="left", color=F.INK2, fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.945])
    fig.savefig(f"{FIG}/slide_ridge.png", dpi=170)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def fig_proof(g, H, bw, S, n):
    T = n["T"]
    pred = g["eff_d"].to_numpy() * S.frac_above(T)
    obs = g["eff_v"].to_numpy()

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.6, 5.4),
                                 gridspec_kw=dict(width_ratios=[1, 1.12]))
    _style(a1)
    lim = (0.25, 1.0)
    a1.plot(lim, lim, ls=(0, (5, 4)), lw=1.3, color=F.MUTED, zorder=2)
    a1.scatter(pred, obs, s=np.clip(g["n_track_v"] / 380, 16, 120),
               facecolor=VMM_C, alpha=0.55, edgecolor=F.SURFACE, linewidth=0.7,
               zorder=4)
    a1.set_xlim(*lim)
    a1.set_ylim(*lim)
    a1.set_aspect("equal")
    a1.set_xlabel("predicted: DREAM's spectrum for that pad, cut at "
                  f"{T:.0f} ADC")
    a1.set_ylabel("measured VMM efficiency of that pad")
    a1.text(0.285, 0.955, f"53 pads, ONE fitted number\n"
                          f"r = {n['r_pred']:+.2f}   ·   residual "
                          f"{n['resid_rms'] * 100:.1f} points",
            fontsize=9.5, color=F.INK, va="top")
    dead = (g["pad_id"] == 635).to_numpy()
    if dead.any():
        a1.annotate("pad 635 — dead in BOTH readouts,\nnot a threshold effect",
                    xy=(pred[dead][0], obs[dead][0]),
                    xytext=(0.30, 0.42), fontsize=8.5, color=F.INK2,
                    va="center",
                    arrowprops=dict(arrowstyle="-", color=F.INK2, lw=0.8))
    a1.set_title("One cut reproduces the VMM's efficiency, pad by pad",
                 loc="left", color=F.INK, fontsize=10, fontweight="bold")

    # -- the pulse-height half, as a waterfall ------------------------------- #
    _style(a2)
    vals = [n["rms_dream_raw"], n["rms_dream_cut"], n["rms_vmm"]]
    labs = ["what the detector does\n(DREAM, everything kept)",
            "…after the VMM's\nthreshold cuts the Landau",
            "…and after the VMM ADC's\n+43-count offset"]
    cols = [DREAM_C, "#b0783f", VMM_C]
    xs = np.arange(3)
    a2.bar(xs, vals, width=0.56, color=cols, zorder=3)
    for x, v in zip(xs, vals):
        a2.text(x, v + 0.012, f"{v:.2f}", ha="center", fontsize=11,
                fontweight="bold", color=F.INK)
    for i in range(2):
        a2.annotate("", xy=(xs[i + 1] - 0.28, vals[i + 1]),
                    xytext=(xs[i] + 0.28, vals[i]),
                    arrowprops=dict(arrowstyle="->", color=F.INK2, lw=1.1,
                                    connectionstyle="arc3,rad=-0.25"))
    a2.text(0.5, 0.455, "the threshold\ncuts the low tail", fontsize=8.5,
            color=F.INK2, ha="center")
    a2.text(1.5, 0.395, "and the ADC counts\nfrom a pedestal", fontsize=8.5,
            color=F.INK2, ha="center")
    a2.set_xticks(xs)
    a2.set_xticklabels(labs, fontsize=8.5, color=F.INK2)
    a2.set_ylim(0, 0.54)
    a2.set_ylabel("pad-to-pad spread of the recorded\npulse height  (rms, "
                  "relative)")
    a2.set_title("The VMM only LOOKS more uniform",
                 loc="left", color=F.INK, fontsize=10, fontweight="bold")

    fig.suptitle("One number — the discriminator level — accounts for the "
                 "whole difference between the two readouts",
                 x=0.006, ha="left", fontsize=12, fontweight="bold",
                 color=F.INK)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"{FIG}/slide_proof.png", dpi=170)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def fig_fix(g, H, bw, S, n):
    T = n["T"]
    fac = np.geomspace(1.0, 4.0, 40)
    w = g["n_track_v"].to_numpy().astype(float)
    ed = g["eff_d"].to_numpy()
    allc, minc = [], []
    for f in fac:
        p = ed * S.frac_above(T / f)
        allc.append(np.average(p, weights=w))
        minc.append(np.sort(p)[1])
    allc, minc = np.array(allc), np.array(minc)

    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    _style(ax)
    ax.axhline(n["eff_dream_all"], lw=1.6, ls=(0, (5, 4)), color=DREAM_C,
               zorder=2)
    ax.text(3.98, n["eff_dream_all"] + 0.006, "DREAM, as measured", fontsize=9,
            color=DREAM_C, fontweight="bold", ha="right")
    ax.plot(fac, allc, lw=2.4, color=VMM_C, zorder=4)
    ax.plot(fac, minc, lw=2.0, color=VMM_C, ls=(0, (4, 3)), zorder=4)
    ax.scatter([1.0], [n["eff_obs_all"]], s=70, color=VMM_C, zorder=6,
               edgecolor=F.SURFACE, linewidth=1.0)
    ax.annotate(f"today  {n['eff_obs_all'] * 100:.0f} %", xy=(1.0, n["eff_obs_all"]),
                xytext=(1.12, 0.80), fontsize=9.5, color=VMM_C,
                fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=VMM_C, lw=0.9))
    ax.text(2.55, allc[np.searchsorted(fac, 2.55)] - 0.028, "all pads",
            fontsize=9.5, color=VMM_C, fontweight="bold", ha="center")
    ax.text(2.55, minc[np.searchsorted(fac, 2.55)] - 0.030,
            "the weakest pad", fontsize=9.5, color=VMM_C, fontweight="bold")
    for f in (1.5, 2.0):
        i = np.searchsorted(fac, f)
        ax.plot([f, f], [0.55, allc[i]], lw=1.0, color=F.MUTED, ls=(0, (2, 3)),
                zorder=2)
        ax.text(f * 1.035, 0.565, f"x{f:g}\n{allc[i] * 100:.0f} %", fontsize=9,
                color=F.INK2, ha="left")
    ax.set_xscale("log")
    ax.set_xticks([1, 1.5, 2, 3, 4])
    ax.set_xticklabels(["1", "1.5", "2", "3", "4"])
    ax.set_xlim(0.97, 4.15)
    ax.set_ylim(0.55, 1.0)
    ax.set_xlabel("signal-over-threshold, relative to today\n"
                  "(raise the mesh gain, lower sdt, or raise the VMM's "
                  "front-end gain — the model cannot tell them apart)")
    ax.set_ylabel("efficiency")
    ax.set_title("Halving the threshold — or doubling the gain — buys back "
                 "8 points overall,\nand more than half the pad-to-pad "
                 "efficiency spread",
                 loc="left", color=F.INK, fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0.085, 1, 1])
    fig.text(0.012, 0.018, "How far sdt can actually drop is set by the noise, "
             "which this measurement does not constrain.\nThe curve says what "
             "the move is worth — not that it is free.",
             ha="left", va="bottom", fontsize=8.5, color=F.INK2)
    fig.savefig(f"{FIG}/slide_fix.png", dpi=170)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def fig_full(g, H, bw, S, n):
    """One 16:9 slide: the mechanism and its consequence on the left, the proof
    and the fix stacked on the right."""
    gr = groups(g, H, bw)
    T = n["T"]
    fig = plt.figure(figsize=(13.33, 7.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[2.05, 0.72, 1.15],
                          height_ratios=[1, 1],
                          left=0.068, right=0.975, top=0.845, bottom=0.085,
                          wspace=0.30, hspace=0.50)
    ax = fig.add_subplot(gs[:, 0])
    axb = fig.add_subplot(gs[:, 1], sharey=ax)
    ridge(ax, gr, bw, T)
    effbars(axb, gr, n)
    ax.annotate("VMM threshold", xy=(T, NGROUP - 0.30),
                xytext=(T * 1.30, NGROUP - 0.05), fontsize=9.5, color=VMM_C,
                fontweight="bold", va="center",
                arrowprops=dict(arrowstyle="-", color=VMM_C, lw=1.0))
    ax.text(T * 0.52, -0.30, "lost", fontsize=10, color=DREAM_C,
            fontweight="bold", ha="center", va="center")
    ax.text(28 * 0.97, NGROUP + 0.08, "pad gain", fontsize=8.5, color=F.INK2,
            ha="right", va="bottom", style="italic")
    ax.set_title("On a weak pad the Landau slides onto the threshold",
                 loc="left", color=F.INK, fontsize=10.5, fontweight="bold")
    axb.set_title("what each\nreadout records", loc="left", color=F.INK,
                  fontsize=9.5)
    axb.text(0.50, -0.30, "VMM", fontsize=9, color=VMM_C, fontweight="bold",
             ha="center", va="center")
    axb.text(0.78, -0.30, "DREAM", fontsize=9, color=DREAM_C,
             fontweight="bold", ha="center", va="center")

    a1 = fig.add_subplot(gs[0, 2])
    _style(a1)
    pred = g["eff_d"].to_numpy() * S.frac_above(T)
    a1.plot([0.25, 1], [0.25, 1], ls=(0, (5, 4)), lw=1.2, color=F.MUTED)
    a1.scatter(pred, g["eff_v"], s=np.clip(g["n_track_v"] / 900, 10, 55),
               facecolor=VMM_C, alpha=0.55, edgecolor=F.SURFACE, linewidth=0.6,
               zorder=4)
    a1.set_xlim(0.25, 1)
    a1.set_ylim(0.25, 1)
    a1.set_aspect("equal")
    a1.set_xticks([0.4, 0.6, 0.8, 1.0])
    a1.set_yticks([0.4, 0.6, 0.8, 1.0])
    a1.set_xlabel("predicted from DREAM + one threshold", fontsize=8.5)
    a1.set_ylabel("VMM, measured", fontsize=8.5)
    a1.set_title(f"That one number predicts\nevery pad   (r = {n['r_pred']:+.2f})",
                 loc="left", color=F.INK, fontsize=10, fontweight="bold")

    a2 = fig.add_subplot(gs[1, 2])
    _style(a2)
    fac = np.geomspace(1.0, 4.0, 40)
    wv = g["n_track_v"].to_numpy().astype(float)
    allc = np.array([np.average(g["eff_d"].to_numpy() * S.frac_above(T / f),
                                weights=wv) for f in fac])
    a2.axhline(n["eff_dream_all"], lw=1.4, ls=(0, (5, 4)), color=DREAM_C)
    a2.text(3.95, n["eff_dream_all"] + 0.008, "DREAM", fontsize=9,
            color=DREAM_C, fontweight="bold", ha="right")
    a2.plot(fac, allc, lw=2.4, color=VMM_C, zorder=4)
    a2.scatter([1.0], [n["eff_obs_all"]], s=55, color=VMM_C, zorder=6,
               edgecolor=F.SURFACE, linewidth=1.0)
    a2.text(1.06, n["eff_obs_all"] - 0.052,
            f"today {n['eff_obs_all'] * 100:.0f} %", fontsize=9, color=VMM_C,
            fontweight="bold")
    i2 = np.searchsorted(fac, 2.0)
    a2.text(2.05, allc[i2] - 0.052, f"x2 -> {allc[i2] * 100:.0f} %", fontsize=9,
            color=VMM_C, fontweight="bold")
    a2.set_xscale("log")
    a2.set_xticks([1, 1.5, 2, 3, 4])
    a2.set_xticklabels(["1", "1.5", "2", "3", "4"])
    a2.set_xlim(0.97, 4.15)
    a2.set_ylim(0.62, 1.0)
    a2.set_xlabel("signal over threshold, relative to today", fontsize=8.5)
    a2.set_ylabel("efficiency", fontsize=8.5)
    a2.set_title("The fix is one factor of two", loc="left", color=F.INK,
                 fontsize=10, fontweight="bold")

    fig.text(0.008, 0.955, "P2_OUT: the gain spread belongs to the CHAMBER — "
             "what the VMM adds is a threshold inside the Landau",
             ha="left", va="bottom", fontsize=14, fontweight="bold",
             color=F.INK)
    fig.text(0.008, 0.912,
             "53 pads, uRWELL tracks, matched runs (VMM run_46 / DREAM "
             "eff_nominal_1, mesh 450 V, drift 750 V).  The two readouts "
             "measure the same per-pad gain map, r = +0.94.",
             fontsize=9.5, color=F.INK2, ha="left", va="bottom")
    fig.savefig(f"{FIG}/slide_full.png", dpi=170)
    plt.close(fig)


def main():
    g, H, bw, S, n = load()
    fig_ridge(g, H, bw, S, n)
    fig_proof(g, H, bw, S, n)
    fig_fix(g, H, bw, S, n)
    fig_full(g, H, bw, S, n)
    print("wrote figures/slide_{ridge,proof,fix,full}.png")


if __name__ == "__main__":
    main()
