#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figures_slide.py -- the MPGD2026 slide figures for P2_OUT: one conclusion,
readable in a minute.

The conclusion: the pad-to-pad gain variation belongs to the DETECTOR -- both
readouts measure the same map -- and the whole difference between them is that
the VMM's discriminator cuts into the Landau, so on the low-gain pads it loses
the small pulses DREAM still records.

Five figures, four of them stand-alone slides and one composite:
  slide_ridge.png  the mechanism, in one panel: the Landau slides across a
                   fixed threshold line as the pad gain falls
  slide_ridge_perpad.png
                   the same panel with each pad's own fitted threshold ticked
                   on it in green, and both colours named on the figure
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
TICK_C = F.C3


def ridge(ax, gr, bw, T, xlim=(28.0, 3000.0), h=1.02, fs=1.0, Tpad_rows=None,
          tick=(0.03, 0.30, 1.2, 0.75)):
    """The hero panel.  One filled curve per gain band, stacked bottom (weak)
    to top (strong), with the single discriminator level drawn straight through
    all of them.

    The x axis is LOGARITHMIC on purpose: a gain difference is a multiplicative
    factor, so on a log axis the eight bands are the same curve TRANSLATED, and
    what changes from row to row is only how much of it falls left of the fixed
    threshold line.  (Plotted as A x dN/dA, the density per decade, so the areas
    on screen are the real fractions.)

    Tpad_rows, if given, is one array per band of the INDEPENDENT per-pad
    threshold fits (threshold_model.fit_threshold_per_pad) for the pads in
    that band, drawn as short ticks against the single global line -- the
    test of whether the global-threshold claim actually holds pad by pad.
    `tick` is (bottom, top, linewidth, alpha) for those ticks, in row units:
    the default is the unobtrusive version for the two-panel check figure,
    slide_ridge_perpad.png asks for taller and more solid ones because there
    they carry half the message."""
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
        if Tpad_rows is not None:
            tp = np.asarray(Tpad_rows[j])
            tp = tp[(tp >= xlim[0]) & (tp <= xlim[1])]
            ax.vlines(tp, base + tick[0], base + tick[1], color=TICK_C,
                      lw=tick[2], alpha=tick[3], zorder=44)

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
def fig_ridge_perpad(g, H, bw, S, n):
    """slide_ridge.png with the per-pad fits drawn on it -- the same eight
    bands and the same efficiency bars, so the slide reads exactly as before,
    plus one green tick per pad at the threshold fitted from that pad alone.

    slide_ridge_check.png answers the same question but spends half the figure
    on the T_pad-vs-gain scatter; this one keeps the original composition and
    only says, in the picture itself, what the two colours are: blue is one
    number fitted to all 53 pads at once, green is that fit repeated 53 times
    independently."""
    gr = groups(g, H, bw)
    T = n["T"]
    Tpad, _, lo_clip, hi_clip = M.fit_threshold_per_pad(g, S)
    Tpad_rows = [Tpad[G["idx"]] for G in gr]
    p = n["perpad"]

    fig, (ax, axb) = plt.subplots(
        1, 2, figsize=(12.2, 6.8), sharey=True,
        gridspec_kw=dict(width_ratios=[2.7, 1], wspace=0.05))
    ridge(ax, gr, bw, T, Tpad_rows=Tpad_rows, tick=(0.02, 0.50, 2.0, 0.95))
    effbars(axb, gr, n)
    # room above for the blue key, below for the green one; sharey, so this
    # sets both panels
    ax.set_ylim(-0.80, NGROUP + 0.72)

    # -- the key: what each of the two colours is ---------------------------- #
    ax.annotate("the VMM's discriminator — one level, all six chips\n"
                "blue = the single global fit, all 53 pads at once",
                xy=(T, NGROUP + 0.02), xytext=(T * 1.28, NGROUP + 0.34),
                fontsize=9.5, color=VMM_C, fontweight="bold", va="center",
                arrowprops=dict(arrowstyle="-", color=VMM_C, lw=1.0))
    ax.vlines([29.6], -0.75, -0.57, color=TICK_C, lw=2.0, zorder=44)
    # zorder above the global line, which spans the whole axes and would
    # otherwise be drawn straight through the sentence
    ax.text(31.5, -0.66, "green = the same fit run on each pad ALONE, one "
            "tick per pad (6–7 per band)", fontsize=9.5, color=TICK_C,
            fontweight="bold", ha="left", va="center", zorder=45)

    ax.text(28 * 0.90, NGROUP + 0.10, "pad gain", fontsize=9, color=F.INK2,
            ha="right", va="bottom", style="italic")
    ax.text(T * 0.52, -0.30, "lost", fontsize=10.5, color=DREAM_C,
            fontweight="bold", ha="center", va="center")
    ax.text(T * 2.1, -0.30, "recorded", fontsize=10.5, color=F.INK2,
            ha="center", va="center")
    ax.annotate("", xy=(T * 0.99, -0.30), xytext=(T * 0.72, -0.30),
                arrowprops=dict(arrowstyle="<-", color=DREAM_C, lw=1.1))
    axb.text(0.405, NGROUP + 0.06, "of the tracks that point at the pad,\n"
             "the % each readout records",
             fontsize=8.5, color=F.INK2, ha="left", va="bottom")
    axb.text(0.50, -0.30, "VMM", fontsize=9.5, color=VMM_C,
             fontweight="bold", ha="center", va="center")
    axb.text(0.72, -0.30, "DREAM", fontsize=9.5, color=DREAM_C,
             fontweight="bold", ha="center", va="center")

    fig.suptitle("On the weak pads the Landau slides down onto the VMM's "
                 "threshold — and every pad, fitted on its own, agrees",
                 x=0.006, ha="left", fontsize=12.5, fontweight="bold",
                 color=F.INK)
    ax.set_title("P2_OUT · 53 pads under the beam, sorted by gain into 8 bands "
                 "· log axis, so a gain factor is a sideways shift",
                 loc="left", color=F.INK2, fontsize=9)
    n_clip = int((lo_clip | hi_clip).sum())
    fig.text(0.006, 0.012,
             f"The {p['n_ok']} independent fits centre on the global level — "
             f"median {p['T_med']:.0f} against {T:.0f} DREAM ADC — and "
             f"scatter {p['T_rms']:.0f} ADC rms ({p['T_relrms'] * 100:.0f} %).\n"
             f"At {p['adc_per_point']:.0f} ADC per efficiency point that is "
             f"the {n['resid_rms'] * 100:.1f}-point residual the global fit "
             "already carries, not a second effect."
             + (f"  {n_clip} pad(s) off the scanned range, not drawn."
                if n_clip else ""),
             fontsize=8.5, color=F.INK2, ha="left", va="bottom")
    fig.tight_layout(rect=[0, 0.052, 1, 0.945])
    fig.savefig(f"{FIG}/slide_ridge_perpad.png", dpi=170)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def fig_ridge_check(g, H, bw, S, n):
    """The test the global fit doesn't otherwise get: fit T pad by pad,
    independently, and see whether the 53 answers cluster on the global line
    or scatter away from it (with gain, or just noisily)."""
    gr = groups(g, H, bw)
    T = n["T"]
    Tpad, sigT, lo_clip, hi_clip = M.fit_threshold_per_pad(g, S)
    Tpad_rows = [Tpad[G["idx"]] for G in gr]
    ok = ~(lo_clip | hi_clip) & np.isfinite(sigT)

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(13.6, 6.4),
        gridspec_kw=dict(width_ratios=[2.0, 1.35], wspace=0.28))

    ridge(ax, gr, bw, T, Tpad_rows=Tpad_rows)
    ax.annotate("global fit — all 53 pads together",
                xy=(T, NGROUP - 0.05), xytext=(T * 1.35, NGROUP + 0.20),
                fontsize=9.5, color=VMM_C, fontweight="bold", va="center",
                arrowprops=dict(arrowstyle="-", color=VMM_C, lw=1.0))
    ax.text(28 * 0.90, -0.34, "ticks: each pad's OWN threshold, fit from that "
            "pad alone", fontsize=8.5, color=TICK_C, fontweight="bold",
            ha="left", va="center")
    ax.set_title("P2_OUT · same 8 gain bands · one tick per pad, from an "
                 "independent per-pad fit", loc="left", color=F.INK2,
                 fontsize=9)

    # -- the actual test: does T_pad cluster on T, or trend/scatter? -------- #
    _style(ax2)
    gain = (g["amp_med_d"].to_numpy() / np.median(g["amp_med_d"].to_numpy()))
    ax2.axhspan(n["T_lo"], n["T_hi"], color=VMM_C, alpha=0.12, zorder=1)
    ax2.axhline(T, lw=2.0, color=VMM_C, zorder=3)
    ax2.errorbar(gain[ok], Tpad[ok], yerr=sigT[ok], fmt="o", ms=5.5,
                color=TICK_C, ecolor=TICK_C, elinewidth=1.0, capsize=2,
                alpha=0.85, zorder=5, label="in-range pads")
    bad = lo_clip | hi_clip
    if bad.any():
        ax2.scatter(gain[bad], Tpad[bad], marker="x", s=42, color=F.MUTED,
                    zorder=4, label="off scan range (ratio ill-defined)")
    ax2.set_xscale("log")
    ax2.set_xlabel("pad gain (DREAM median / fleet median)")
    ax2.set_ylabel("independent per-pad threshold fit  [DREAM ADC]")
    ax2.legend(loc="upper right", fontsize=8, frameon=False)

    n_ok = int(ok.sum())
    med, iqr_lo, iqr_hi = (np.median(Tpad[ok]), *np.percentile(Tpad[ok], [25, 75])) \
        if n_ok else (np.nan, np.nan, np.nan)
    chi2 = float(np.sum(((Tpad[ok] - T) / sigT[ok]) ** 2)) if n_ok else np.nan
    dof = max(n_ok - 1, 1)
    ax2.text(0.03, 0.03, f"{n_ok}/{len(g)} pads in range\n"
             f"median T_pad = {med:.0f} ADC  (IQR {iqr_lo:.0f}–{iqr_hi:.0f})\n"
             f"global T = {T:.0f} ADC\n"
             f"$\\chi^2$/dof vs global = {chi2 / dof:.2f}  ({n_ok} pads, "
             "counting-noise error bars)",
             transform=ax2.transAxes, fontsize=8.5, color=F.INK,
             va="bottom", ha="left",
             bbox=dict(boxstyle="round,pad=0.35", fc=F.SURFACE, ec=F.GRID))
    ax2.set_title("Independent per-pad fits vs. the global line",
                 loc="left", color=F.INK, fontsize=10, fontweight="bold")

    fig.suptitle("Testing the global-threshold claim: fit each pad on its own, "
                 "not jointly — do the 53 answers agree?",
                 x=0.006, ha="left", fontsize=12.5, fontweight="bold",
                 color=F.INK)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"{FIG}/slide_ridge_check.png", dpi=170)
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
    fig_ridge_perpad(g, H, bw, S, n)
    fig_ridge_check(g, H, bw, S, n)
    fig_proof(g, H, bw, S, n)
    fig_fix(g, H, bw, S, n)
    fig_full(g, H, bw, S, n)
    print("wrote figures/slide_{ridge,ridge_perpad,ridge_check,proof,fix,"
          "full}.png")


if __name__ == "__main__":
    main()
