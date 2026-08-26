#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figures_dv.py -- P2_OUT, VMM vs DREAM per-pad pulse height.

Reads data/compare_dream_vmm_P2_OUT.csv plus the two per-pad histogram sets,
writes figures/dv_*.png.  Two series throughout: VMM = categorical slot 1,
DREAM = slot 2 (validated as a pair, light surface).
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap

import figures as F                       # design tokens + rcParams

HERE = os.path.dirname(os.path.abspath(__file__))
FIG, DATA = os.path.join(HERE, "figures"), os.path.join(HERE, "data")

VMM_C, DREAM_C = F.C1, F.C2
VMM_LBL = "VMM3a  ·  run_46"
DREAM_LBL = "DREAM  ·  eff_nominal_1"


# Blues, but starting at step ~250 rather than at the surface: on an ORDINAL
# map every pad must stay visible, and the bottom of the full ramp is under
# 2:1 against #fcfcfb.
_B = plt.get_cmap("Blues")
BLUES = LinearSegmentedColormap.from_list(
    "blues_vis", _B(np.linspace(0.24, 1.0, 256)))


def _style(ax):
    ax.grid(True, axis="both", linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)


def load():
    m = pd.read_csv(os.path.join(DATA, "compare_dream_vmm_P2_OUT.csv"))
    g = m[m["use"]].reset_index(drop=True)
    d = json.load(open(os.path.join(DATA, "report_run_46_rate.json")))
    pp = d["stations"]["P2_OUT"]["per_pad"]
    hb = pp["adc_hist_bin"]
    VH = np.zeros((len(pp["channel_id"]), 1024 // hb), np.int64)
    for k, v in pp["adc_hist"].items():
        VH[int(k)] = v
    z = np.load(os.path.join(DATA,
                             "dream_padadc_eff_nominal_1_P2_OUT_hist.npz"))
    return g, VH, float(hb), z["hist_own"], float(z["amp_bin"])


def rel(g, col):
    return (g[col] / g[col].median()).to_numpy()


# --------------------------------------------------------------------------- #
def fig_headline(g, DH, dbin):
    rv, rd = rel(g, "amp_med_v"), rel(g, "amp_med_d")
    b = np.polyfit(rd, rv, 1)
    r = np.corrcoef(rd, rv)[0, 1]

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(11.4, 4.9), gridspec_kw=dict(width_ratios=[1.15, 1]))

    _style(ax)
    lim = (0.4, 2.1)
    ax.plot(lim, lim, ls=(0, (5, 4)), lw=1.4, color=F.MUTED, zorder=2)
    ax.annotate("if the two readouts\nmeasured the same spread",
                xy=(1.72, 1.72), xytext=(1.30, 1.93), fontsize=8,
                color=F.MUTED, ha="center",
                arrowprops=dict(arrowstyle="-", color=F.MUTED, lw=0.8))
    xs = np.linspace(*lim, 50)
    ax.plot(xs, b[0] * xs + b[1], lw=2.0, color=VMM_C, zorder=3)
    ax.scatter(rd, rv, s=np.clip(g["n_track_v"] / 380, 14, 110),
               facecolor=VMM_C, alpha=0.55, edgecolor=F.SURFACE, linewidth=0.7,
               zorder=4)
    ax.annotate(f"slope {b[0]:.2f}", xy=(1.60, b[0] * 1.60 + b[1]),
                xytext=(1.72, 0.78), fontsize=9.5,
                color=VMM_C, fontweight="bold", ha="center",
                arrowprops=dict(arrowstyle="-", color=VMM_C, lw=0.9))
    dead = g["pad_id"] == 635
    if dead.any():
        ax.scatter(rd[dead], rv[dead], s=95, facecolor="none", edgecolor=F.INK,
                   linewidth=1.3, zorder=5)
        ax.annotate("pad 635 (the dead pad\nboth readouts find)",
                    xy=(rd[dead][0], rv[dead][0]),
                    xytext=(rd[dead][0] + 0.30, rv[dead][0] - 0.22),
                    fontsize=8, color=F.INK2,
                    arrowprops=dict(arrowstyle="-", color=F.INK2, lw=0.8))
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal")
    ax.set_xlabel("DREAM pulse height  /  its own fleet median")
    ax.set_ylabel("VMM pulse height  /  its own fleet median")
    ax.set_title(f"53 pads, same chamber, same tracks   ·   r = {r:+.3f}",
                 loc="left", color=F.INK, fontsize=9.5)

    # -- the control: cut the DREAM spectra from below and watch the spread -- #
    _style(ax2)
    di = g["di"].to_numpy()
    cuts = np.arange(0, 340, 20.0)
    rms = []
    keep = []
    for c in cuts:
        b0 = int(np.ceil(c / dbin))
        hh = DH[di][:, b0:]
        med = np.array([_hq(h, dbin) for h in hh]) + b0 * dbin
        rms.append(np.std(med / np.median(med), ddof=1))
        keep.append(hh.sum() / DH[di].sum())
    rms = np.array(rms)
    ax2.plot(cuts, rms, lw=2.0, color=DREAM_C, zorder=4)
    ax2.axhline(np.std(rv, ddof=1), lw=2.0, color=VMM_C, zorder=3)
    ax2.annotate(f"VMM measures {np.std(rv, ddof=1):.2f}",
                 xy=(30, np.std(rv, ddof=1)), xytext=(18, 0.268),
                 fontsize=9.5, color=VMM_C, fontweight="bold", ha="left")
    ax2.annotate("DREAM, uncut", xy=(4, rms[0]), xytext=(28, rms[0] + 0.012),
                 fontsize=9, color=DREAM_C, fontweight="bold")
    ax2.axvline(136, lw=1.2, ls=(0, (4, 3)), color=F.MUTED, zorder=2)
    ax2.annotate("where the VMM's\nturn-on sits (136)\n— 9 % of pulses",
                 xy=(136, 0.44), xytext=(150, 0.435), fontsize=8, color=F.INK2)
    ax2.set_xlabel("low-amplitude cut applied to the DREAM spectra  [DREAM ADC]")
    ax2.set_ylabel("spread of the per-pad relative pulse height  (rms)")
    ax2.set_ylim(0.24, 0.46)
    ax2.set_xlim(-6, 330)
    ax2.set_title("A threshold shrinks the measured spread — but not by enough",
                  loc="left", color=F.INK, fontsize=9.5)

    fig.suptitle("The VMM sees the SAME pad-to-pad structure as DREAM, "
                 "compressed to 0.6x by its threshold",
                 x=0.006, ha="left", fontsize=11.5, fontweight="bold",
                 color=F.INK)
    fig.tight_layout(rect=[0, 0, 1, 0.935])
    fig.savefig(f"{FIG}/dv_headline.png", dpi=170)
    plt.close(fig)
    return dict(slope=float(b[0]), intercept=float(b[1]), r=float(r),
                rms_vmm=float(np.std(rv, ddof=1)),
                rms_dream=float(rms[0]),
                rms_dream_at_matched=float(np.interp(136, cuts, rms)),
                keep_at_matched=float(np.interp(136, cuts, keep)))


def _hq(h, binw, q=0.5):
    t = h.sum()
    if t == 0:
        return np.nan
    cum = np.cumsum(h) / t
    i = int(np.searchsorted(cum, q, side="left"))
    below = cum[i - 1] if i > 0 else 0.0
    return float((i + (q - below) / max(cum[i] - below, 1e-12)) * binw)


# --------------------------------------------------------------------------- #
def fig_spectra(g, VH, vbin, DH, dbin):
    """Pooled, then one low-gain and one high-gain pad, on a common axis."""
    vi, di = g["vi"].to_numpy(), g["di"].to_numpy()
    vmpv, dmpv = 111.8, 237.5          # pooled Landau peaks, both readouts
    order = np.argsort(rel(g, "amp_med_d"))
    picks = [("pooled  ·  all 53 pads", None),
             (f"a low-gain pad  ·  {int(g['pad_id'][order[2]])}", order[2]),
             (f"a high-gain pad  ·  {int(g['pad_id'][order[-1]])}", order[-1])]

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.3), sharey=True)
    for ax, (title, k) in zip(axes, picks):
        _style(ax)
        hv = VH[vi].sum(0) if k is None else VH[vi[k]]
        hd = DH[di].sum(0) if k is None else DH[di[k]]
        cv = (np.arange(len(hv)) + 0.5) * vbin / vmpv
        cd = (np.arange(len(hd)) + 0.5) * dbin / dmpv
        ax.step(cv, hv / hv.sum() / (vbin / vmpv), where="mid", lw=2.0,
                color=VMM_C, zorder=4)
        ax.step(cd, hd / hd.sum() / (dbin / dmpv), where="mid", lw=2.0,
                color=DREAM_C, zorder=3)
        ax.set_xlim(0, 4.2)
        ax.set_yscale("log")
        ax.set_ylim(2e-3, 3)
        ax.set_xlabel("pulse height  /  fleet MPV of that readout")
        ax.set_title(title, loc="left", color=F.INK, fontsize=9.5)
        # the medians, which are the numbers every other figure is built on:
        # on a per-pad panel the gap between the two ticks IS the compression
        if k is not None:
            for h, c, mp, bw, yy in ((hv, VMM_C, vmpv, vbin, 1.7),
                                     (hd, DREAM_C, dmpv, dbin, 2.4)):
                mm = _hq(h, bw) / mp
                ax.plot([mm, mm], [2e-3, yy], lw=1.2, ls=(0, (2, 2)), color=c,
                        zorder=5)
                far = mm > 3.4
                ax.annotate(f"{mm:.2f}", xy=(mm, yy),
                            xytext=(mm - 0.08 if far else mm + 0.08, yy),
                            ha="right" if far else "left",
                            fontsize=8.5, color=c, fontweight="bold")
    axes[0].set_ylabel("probability density")
    axes[0].axvspan(0, 64 / vmpv, color=VMM_C, alpha=0.08, zorder=1)
    axes[0].annotate("nothing below the\nVMM discriminator", xy=(0.28, 1.1),
                     xytext=(0.62, 1.35), fontsize=8, color=VMM_C,
                     arrowprops=dict(arrowstyle="->", color=VMM_C, lw=0.9))
    axes[0].annotate(VMM_LBL, xy=(1.0, 0.55), xytext=(1.9, 0.9), fontsize=9,
                     color=VMM_C, fontweight="bold")
    axes[0].annotate(DREAM_LBL, xy=(1.0, 0.30), xytext=(1.9, 0.30), fontsize=9,
                     color=DREAM_C, fontweight="bold")
    fig.suptitle("Same pad, two readouts: DREAM records the whole Landau, "
                 "the VMM only what clears its threshold",
                 x=0.006, ha="left", fontsize=11.5, fontweight="bold",
                 color=F.INK)
    fig.tight_layout(rect=[0, 0, 1, 0.925])
    fig.savefig(f"{FIG}/dv_spectra.png", dpi=170)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def fig_maps(g):
    """The gradient, as two maps and as one profile along its own axis."""
    rv, rd = rel(g, "amp_med_v"), rel(g, "amp_med_d")
    x, y = g["x_v"].to_numpy(), g["y_v"].to_numpy()
    A = np.column_stack([np.ones(len(g)), x, y])
    bd = np.linalg.lstsq(A, rd, rcond=None)[0]
    ang = np.arctan2(bd[2], bd[1])
    s = (x - x.mean()) * np.cos(ang) + (y - y.mean()) * np.sin(ang)

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.4),
                             gridspec_kw=dict(width_ratios=[1, 1, 1.15]))
    for ax, r, lab in ((axes[0], rv, VMM_LBL), (axes[1], rd, DREAM_LBL)):
        _style(ax)
        sc = ax.scatter(x, y, c=r, s=95, cmap=BLUES,
                        norm=Normalize(0.5, 2.0), edgecolor=F.SURFACE,
                        linewidth=0.8, zorder=3)
        ax.set_aspect("equal")
        ax.set_xlabel("x  [mm]")
        ax.set_title(lab, loc="left", color=F.INK, fontsize=9.5)
        cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
        cb.outline.set_edgecolor(F.MUTED)
        cb.ax.tick_params(labelsize=8, color=F.MUTED)
    axes[0].set_ylabel("y  [mm]")
    for ax in axes[:2]:
        ax.annotate("", xy=(x.mean() + 34 * np.cos(ang),
                            y.mean() + 34 * np.sin(ang)),
                    xytext=(x.mean() - 34 * np.cos(ang),
                            y.mean() - 34 * np.sin(ang)),
                    arrowprops=dict(arrowstyle="-|>", color=F.INK, lw=1.6),
                    zorder=6)
    axes[0].text(x.mean() - 40 * np.cos(ang), y.mean() - 40 * np.sin(ang) - 9,
                 "gain gradient", fontsize=8.5, color=F.INK, ha="left")

    ax = axes[2]
    _style(ax)
    for r, c, lab in ((rd, DREAM_C, DREAM_LBL), (rv, VMM_C, VMM_LBL)):
        ax.scatter(s, r, s=26, color=c, alpha=0.6, edgecolor=F.SURFACE,
                   linewidth=0.5, zorder=3)
        bb = np.polyfit(s, r, 1)
        ss = np.linspace(s.min(), s.max(), 20)
        ax.plot(ss, bb[0] * ss + bb[1], lw=2.0, color=c, zorder=4)
        ax.annotate(f"{lab}\n{bb[0] * 10 * 100:+.0f} % per 10 mm",
                    xy=(ss[-1], bb[0] * ss[-1] + bb[1]),
                    xytext=(s.min() + 4, 2.05 if c == DREAM_C else 1.62),
                    fontsize=8.5, color=c, fontweight="bold")
    ax.set_xlabel("distance along the gradient  [mm]")
    ax.set_ylabel("pulse height  /  fleet median")
    ax.set_ylim(0.35, 2.3)
    ax.set_title("Same gradient, VMM at 0.6x the contrast",
                 loc="left", color=F.INK, fontsize=9.5)

    fig.suptitle("P2_OUT's pad-to-pad spread is one smooth gain gradient "
                 "across the beam spot",
                 x=0.006, ha="left", fontsize=11.5, fontweight="bold",
                 color=F.INK)
    fig.tight_layout(rect=[0, 0, 1, 0.925])
    fig.savefig(f"{FIG}/dv_maps.png", dpi=170)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def fig_efficiency(g):
    rv, rd = rel(g, "amp_med_v"), rel(g, "amp_med_d")
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    _style(ax)
    for r, e, c, lab in ((rd, g["eff_d"], DREAM_C, DREAM_LBL),
                         (rv, g["eff_v"], VMM_C, VMM_LBL)):
        ax.scatter(r, e, s=32, color=c, alpha=0.6, edgecolor=F.SURFACE,
                   linewidth=0.5, zorder=3)
    ax.text(1.30, 0.72, f"{DREAM_LBL}\n0.956 overall   ·   0.93–0.99 per pad",
            fontsize=9, color=DREAM_C, fontweight="bold")
    ax.text(1.30, 0.56, f"{VMM_LBL}\n0.853 overall   ·   0.50–0.96 per pad",
            fontsize=9, color=VMM_C, fontweight="bold")
    d = g["pad_id"] == 635
    ax.scatter(rd[d], g["eff_d"][d], s=95, facecolor="none", edgecolor=F.INK,
               linewidth=1.3, zorder=5)
    ax.scatter(rv[d], g["eff_v"][d], s=95, facecolor="none", edgecolor=F.INK,
               linewidth=1.3, zorder=5)
    ax.annotate("pad 635 — the one pad that is\ndead in both readouts",
                xy=(rv[d][0], g["eff_v"][d].iloc[0]),
                xytext=(0.86, 0.30), fontsize=8, color=F.INK2, va="center",
                arrowprops=dict(arrowstyle="-", color=F.INK2, lw=0.8))
    ax.set_xlabel("pulse height  /  fleet median of that readout")
    ax.set_ylabel("per-pad efficiency")
    ax.set_ylim(0.25, 1.02)
    ax.set_title("One gain map, two readouts: the VMM's per-pad efficiency "
                 "spreads eight times wider",
                 loc="left", color=F.INK, fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    fig.savefig(f"{FIG}/dv_efficiency.png", dpi=170)
    plt.close(fig)


def main():
    g, VH, vbin, DH, dbin = load()
    s = fig_headline(g, DH, dbin)
    fig_spectra(g, VH, vbin, DH, dbin)
    fig_maps(g)
    fig_efficiency(g)
    print(json.dumps(s, indent=1))
    print("wrote figures/dv_*.png")


if __name__ == "__main__":
    main()
