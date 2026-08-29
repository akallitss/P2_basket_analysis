#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figures_perpad.py -- the hero ridge with the averaging taken OFF.

slide_ridge.png draws eight curves, but P2_OUT has 53 pads: each curve is a
gain band, the track-weighted mean of six or seven pads' unit-area spectra.
That averaging is honest for the model -- a band's sub-threshold fraction is
exactly the mean of its pads' -- but it is not neutral for the EYE.  Averaging
six spectra that sit at slightly different places produces one smooth curve
that sits where none of them is, and the question "do the cut-offs line up?"
cannot be asked of a picture that has already answered it by construction.

So: same plot, 53 rows, no averaging, each pad carrying its own independently
fitted threshold next to the one global line.  Then three panels that say what
the spread in those 53 answers is worth --

  perpad_ridge.png    every pad on its own row, its own fit, the global line
  perpad_check.png    (a) where the spread lives -- within a gain band, or
                          between bands?
                      (b) is it wider than counting noise alone would make it?
                      (c) what it is in the units that matter: efficiency

The conclusion the three panels reach is that the per-pad threshold spread and
the 5-point efficiency residual already quoted on slide_proof.png are the same
fact in two different units, not two separate problems.
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import figures as F
import figures_slide as SL
import threshold_model as M

HERE = os.path.dirname(os.path.abspath(__file__))
FIG, DATA = os.path.join(HERE, "figures"), os.path.join(HERE, "data")

VMM_C, DREAM_C, TICK_C = F.C1, F.C2, F.C3
NBAND = SL.NGROUP


def load():
    g, H, bw = M.load()
    S = M.Spectra(H, bw)
    n = json.load(open(os.path.join(DATA, "threshold_model_P2_OUT.json")))
    Tpad, sigT, lo, hi = M.fit_threshold_per_pad(g, S)
    return g, H, bw, S, n, Tpad, sigT, lo | hi


def bands(g, nband=NBAND):
    """Which gain band each pad falls in -- the same split slide_ridge.png
    averages over, kept here only so the spread can be decomposed against it."""
    order = np.argsort(g["amp_med_d"].to_numpy())
    b = np.empty(len(g), int)
    for k, idx in enumerate(np.array_split(order, nband)):
        b[idx] = k
    return b, order


# --------------------------------------------------------------------------- #
def fig_perpad_ridge(g, H, bw, S, n, Tpad, clip):
    """53 rows, one pad each, weakest at the bottom -- slide_ridge.png with the
    band averaging removed."""
    T = n["T"]
    xlim = (28.0, 3000.0)
    order = np.argsort(g["amp_med_d"].to_numpy())
    gain = g["amp_med_d"].to_numpy() / np.median(g["amp_med_d"].to_numpy())
    pred = g["eff_d"].to_numpy() * S.frac_above(T)

    c = (np.arange(H.shape[1]) + 0.5) * bw
    k = (c >= xlim[0]) & (c <= 3600.0)          # above 3600 DREAM piles up

    fig, (ax, axb) = plt.subplots(
        1, 2, figsize=(12.6, 15.0), sharey=True,
        gridspec_kw=dict(width_ratios=[2.6, 1], wspace=0.04))
    # tight_layout cannot handle 53 stacked rows with a shared y; place the
    # axes by hand so the rows get the whole page
    fig.subplots_adjust(left=0.068, right=0.985, top=0.955, bottom=0.042)

    for row, j in enumerate(order):
        h = H[j] / max(H[j].sum(), 1.0)
        # 5 bins, not the 3 the band plot uses: a single pad has ~1e3-7e4
        # pulses spread over 450 bins and is visibly grainy at 3
        y = np.convolve(h, np.ones(5) / 5, "same")[k] * c[k]
        y = y / max(y.max(), 1e-12) * 1.05
        base = float(row)
        ax.fill_between(c[k], base, base + y, color=DREAM_C, alpha=0.13, lw=0,
                        zorder=3 + row * 0.01)
        lost = c[k] <= T
        ax.fill_between(c[k][lost], base, (base + y)[lost], color=DREAM_C,
                        alpha=0.70, lw=0, zorder=3 + row * 0.01)
        ax.plot(c[k], base + y, lw=0.9, color=DREAM_C, zorder=3.5 + row * 0.01)
        ax.plot(xlim, [base, base], lw=0.5, color=F.GRID, zorder=1)
        ax.text(xlim[0] * 0.93, base + 0.30, f"{gain[j]:.2f}x", fontsize=6.4,
                color=F.INK2, ha="right", va="center")
        # the pad's own answer, fitted from that pad alone
        if clip[j]:
            ax.plot([np.clip(Tpad[j], *xlim)], [base + 0.20], marker="x",
                    ms=4.0, color=F.MUTED, zorder=44)
        else:
            ax.vlines(Tpad[j], base + 0.02, base + 0.60, color=TICK_C, lw=2.0,
                      zorder=44)

    ax.axvline(T, lw=2.0, color=VMM_C, zorder=40)
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.9, len(g) + 0.6)
    ax.set_yticks([])
    ax.set_xticks([50, 100, 200, 500, 1000, 2000])
    ax.set_xticklabels(["50", "100", "200", "500", "1000", "2000"])
    ax.set_xlabel("pulse height on the pad, as DREAM records it  [ADC]")
    ax.spines["left"].set_visible(False)
    ax.grid(True, axis="x", lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.annotate("global fit,\nall 53 pads at once", xy=(T, len(g) - 0.2),
                xytext=(T * 1.45, len(g) + 0.05), fontsize=9, color=VMM_C,
                fontweight="bold", va="center",
                arrowprops=dict(arrowstyle="-", color=VMM_C, lw=1.0))
    # the x marker only exists if some pad's ratio put it off the scan range;
    # with the July data none do, so do not advertise a key that has no entry
    ax.text(xlim[0] * 0.93, -0.62, "green tick = that pad's OWN threshold, "
            "fitted from that pad alone"
            + ("   ·   x = solution outside the scanned range"
               if clip.any() else ""),
            fontsize=8.2, color=TICK_C, fontweight="bold", ha="left",
            va="center")

    # -- right: what the pad actually recorded, and what the one cut predicts - #
    for row, j in enumerate(order):
        yd, yv = row + 0.28, row + 0.62
        axb.barh(yd, g["eff_d"].to_numpy()[j], height=0.30, color=DREAM_C,
                 alpha=0.80, zorder=3)
        axb.barh(yv, g["eff_v"].to_numpy()[j], height=0.30, color=VMM_C,
                 zorder=3)
        axb.plot([pred[j]], [yv], marker="|", ms=9, mew=1.8, color=F.INK,
                 zorder=6)
    axb.set_xlim(0.0, 1.0)
    axb.set_yticks([])
    axb.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    axb.set_xticklabels(["0", "25 %", "50 %", "75 %", "100 %"])
    axb.spines["left"].set_visible(False)
    axb.grid(True, axis="x", lw=0.6, alpha=0.7)
    axb.set_axisbelow(True)
    axb.text(0.02, -0.62, "bars: measured efficiency   ·   | : predicted "
             "from the global cut", fontsize=8.2, color=F.INK2, ha="left",
             va="center")
    axb.text(0.60, len(g) + 0.05, "VMM", fontsize=9, color=VMM_C,
             fontweight="bold", ha="center")
    axb.text(0.86, len(g) + 0.05, "DREAM", fontsize=9, color=DREAM_C,
             fontweight="bold", ha="center")

    fig.suptitle("The same figure with the band averaging taken off: 53 pads, "
                 "53 independent thresholds",
                 x=0.006, y=0.988, ha="left", fontsize=13, fontweight="bold",
                 color=F.INK)
    ax.set_title("P2_OUT · every pad on its own row, sorted by gain · log "
                 "axis, so a gain factor is a sideways shift",
                 loc="left", color=F.INK2, fontsize=9)
    fig.savefig(f"{FIG}/perpad_ridge.png", dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def fig_perpad_check(g, H, bw, S, n, Tpad, sigT, clip):
    """What the spread in those 53 answers is worth: where it lives, whether
    it is bigger than counting noise, and what it costs in efficiency."""
    T = n["T"]
    ok = ~clip & np.isfinite(sigT)
    b, _ = bands(g)
    gain = g["amp_med_d"].to_numpy() / np.median(g["amp_med_d"].to_numpy())
    effv, effd = g["eff_v"].to_numpy(), g["eff_d"].to_numpy()
    pred = effd * S.frac_above(T)

    # variance of T_pad split into a part that separates the gain bands and a
    # part that lives inside them -- i.e. how much of the spread the eight-row
    # figure averages away
    t, bb = Tpad[ok], b[ok]
    gm = t.mean()
    betw = sum((bb == q).sum() * (t[bb == q].mean() - gm) ** 2
               for q in range(NBAND) if (bb == q).any())
    with_ = sum(((t[bb == q] - t[bb == q].mean()) ** 2).sum()
                for q in range(NBAND) if (bb == q).any())
    f_within = 100 * with_ / (with_ + betw)

    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15.4, 5.5))

    # -- (a) within a band vs between bands --------------------------------- #
    SL._style(a1)
    a1.axhspan(n["T_lo"], n["T_hi"], color=VMM_C, alpha=0.12, zorder=1)
    a1.axhline(T, lw=2.0, color=VMM_C, zorder=3)
    rng = np.random.default_rng(7)
    for q in range(NBAND):
        m = ok & (b == q)
        if not m.any():
            continue
        x = q + rng.uniform(-0.20, 0.20, m.sum())
        a1.errorbar(x, Tpad[m], yerr=sigT[m], fmt="o", ms=4.5, color=TICK_C,
                    ecolor=TICK_C, elinewidth=0.8, capsize=0, alpha=0.75,
                    zorder=5)
        a1.plot([q - 0.36, q + 0.36], [np.median(Tpad[m])] * 2, lw=2.4,
                color=F.INK, zorder=6)
    a1.set_xticks(range(NBAND))
    a1.set_xticklabels([f"{np.average(gain[b == q]):.2f}x" for q in range(NBAND)],
                       fontsize=8)
    a1.set_xlabel("gain band  (the eight rows of slide_ridge.png)")
    a1.set_ylabel("independent per-pad threshold  [DREAM ADC]")
    a1.text(0.03, 0.03, f"{f_within:.0f} % of the spread is WITHIN a band\n"
            f"{100 - f_within:.0f} % separates the bands\n"
            "(black bar = band median)",
            transform=a1.transAxes, fontsize=8.5, color=F.INK, va="bottom",
            bbox=dict(boxstyle="round,pad=0.35", fc=F.SURFACE, ec=F.GRID))
    a1.set_title("The averaging hides the spread, it does not remove it",
                 loc="left", color=F.INK, fontsize=10, fontweight="bold")

    # -- (b) is it wider than counting noise? -------------------------------- #
    SL._style(a2)
    edges = np.arange(60, 320, 15.0)
    cnt, _, _ = a2.hist(Tpad[ok], bins=edges, color=TICK_C, alpha=0.55,
                        zorder=3)
    a2.axvline(T, lw=2.0, color=VMM_C, zorder=5)
    # what the 53 answers would have looked like if one common threshold were
    # the whole truth and only track/pulse counting noise moved them.  That
    # curve is ~10x narrower than the histogram, so it is drawn as the band it
    # is rather than as a spike that would own the y axis on its own.
    s = float(np.median(sigT[ok]))
    a2.axvspan(T - s, T + s, color=VMM_C, alpha=0.30, zorder=4)
    a2.set_ylim(0, cnt.max() * 1.35)
    a2.annotate(f"if one common threshold were the whole\ntruth, all 53 would "
                f"sit in here (±{s:.0f} ADC)",
                xy=(T + s, cnt.max() * 0.62),
                xytext=(T + 42, cnt.max() * 0.78), fontsize=8.5, color=VMM_C,
                fontweight="bold", va="center",
                arrowprops=dict(arrowstyle="->", color=VMM_C, lw=1.1))
    a2.set_xlabel("independent per-pad threshold  [DREAM ADC]")
    a2.set_ylabel("pads")
    a2.text(0.025, 0.985, f"observed: median {np.median(Tpad[ok]):.0f}, "
            f"rms {np.std(Tpad[ok], ddof=1):.0f} ADC\n"
            f"counting noise alone: {s:.0f} ADC\n"
            f"$\\chi^2$/dof vs one common T = "
            f"{np.sum(((Tpad[ok] - T) / sigT[ok]) ** 2) / max(ok.sum() - 1, 1):.0f}",
            transform=a2.transAxes, fontsize=8.5, color=F.INK, va="top",
            ha="left",
            bbox=dict(boxstyle="round,pad=0.35", fc=F.SURFACE, ec=F.GRID))
    a2.set_title("Wider than statistics — a common threshold would\n"
                 "have put all 53 inside the blue band", loc="left",
                 color=F.INK, fontsize=10, fontweight="bold")

    # -- (c) the same spread, in efficiency ---------------------------------- #
    SL._style(a3)
    dT = Tpad - T
    resid = (effv - pred) * 100.0          # points; = effd*(F(T_pad) - F(T))
    a3.axhline(0, lw=1.0, color=F.MUTED, zorder=2)
    a3.axvline(0, lw=1.0, color=F.MUTED, zorder=2)
    a3.scatter(dT[ok], resid[ok], s=np.clip(g["n_track_v"].to_numpy()[ok] / 400,
                                            14, 90),
               facecolor=TICK_C, alpha=0.60, edgecolor=F.SURFACE, linewidth=0.7,
               zorder=4)
    if (~ok).any():
        a3.scatter(np.clip(dT[~ok], -140, 160), resid[~ok], marker="x", s=34,
                   color=F.MUTED, zorder=4)
    rms = n["resid_rms"] * 100.0
    a3.axhspan(-rms, rms, color=VMM_C, alpha=0.10, zorder=1)
    a3.set_xlabel("pad's own threshold minus the global one  [DREAM ADC]")
    a3.set_ylabel("efficiency the global cut gets wrong  [points]")
    a3.text(0.03, 0.03, "the two axes are the same number twice:\n"
            f"~10 ADC of threshold = 1 point of efficiency\n"
            f"shaded: the {rms:.1f}-point residual already\n"
            "quoted for the global fit",
            transform=a3.transAxes, fontsize=8.5, color=F.INK, va="bottom",
            bbox=dict(boxstyle="round,pad=0.35", fc=F.SURFACE, ec=F.GRID))
    a3.set_title("…and in efficiency it is the residual\nwe already had",
                 loc="left", color=F.INK, fontsize=10, fontweight="bold")

    fig.suptitle("What the per-pad threshold spread is, where it lives, and "
                 "what it costs",
                 x=0.006, ha="left", fontsize=12.5, fontweight="bold",
                 color=F.INK)
    fig.tight_layout(rect=[0, 0, 1, 0.925])
    fig.savefig(f"{FIG}/perpad_check.png", dpi=170)
    plt.close(fig)


def main():
    g, H, bw, S, n, Tpad, sigT, clip = load()
    fig_perpad_ridge(g, H, bw, S, n, Tpad, clip)
    fig_perpad_check(g, H, bw, S, n, Tpad, sigT, clip)
    print("wrote figures/perpad_{ridge,check}.png")


if __name__ == "__main__":
    main()
