#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figures_p2out.py -- P2_OUT only: the per-pad pulse height / efficiency result."""

import json, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy import stats

import figures as F                      # tokens + loader

FIG, DATA = F.FIG, F.DATA
ST = "P2_OUT"
C1, C2, INK, INK2, MUTED, SURFACE = F.C1, F.C2, F.INK, F.INK2, F.MUTED, F.SURFACE


def get():
    df = F.load()
    return df[df.station == ST].copy()


def fig_main(g):
    u = g[g.use]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3),
                             gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]; F._style(ax)
    ax.errorbar(u.adc_med, u.eff,
                yerr=np.sqrt(np.clip(u.eff*(1-u.eff), 0, None)/u.n_track),
                fmt="none", ecolor=MUTED, elinewidth=0.8, zorder=2)
    ax.scatter(u.adc_med, u.eff, s=10 + 46*np.sqrt(u.n_track/u.n_track.max()),
               c=C1, edgecolor=SURFACE, linewidth=1.0, zorder=3)
    b = np.polyfit(u.adc_med, u.eff, 1)
    xs = np.array([u.adc_med.min(), u.adc_med.max()])
    ax.plot(xs, np.clip(np.polyval(b, xs), None, 1.0), color=C1, lw=2.2, zorder=4)
    r, p = stats.pearsonr(u.adc_med, u.eff)
    rho, _ = stats.spearmanr(u.adc_med, u.eff)
    ax.axhline(0.9604, color=C2, lw=1.4, ls=(0, (5, 3)), zorder=5)
    ax.text(u.adc_med.min(), 0.968, "DREAM, same chamber  0.960", color=C2,
            fontsize=8.5, ha="left", va="bottom")
    ax.text(0.035, 0.05,
            f"Pearson r = {r:+.2f}   (p = {p:.0e})\nSpearman ρ = {rho:+.2f}\n"
            f"slope {b[0]*1000:+.1f} pts / 10 ADC",
            transform=ax.transAxes, fontsize=9, color=INK2, va="bottom")
    ax.set_xlabel("pad median pulse height  [VMM ADC, DNL-corrected]")
    ax.set_ylabel("pad efficiency")
    ax.set_ylim(0.55, 1.02)
    ax.set_title(f"{ST}: {len(u)} illuminated pads, all chips at sdt 224",
                 loc="left", color=INK, fontweight="bold", fontsize=10)

    # track-weighted quintiles
    ax = axes[1]; F._style(ax)
    s = u.sort_values("adc_med")
    cw = np.cumsum(s.n_track.to_numpy())
    cut = np.searchsorted(cw, np.linspace(0, cw[-1], 6)[1:-1])
    xs, ys, es, lab = [], [], [], []
    for q in np.split(np.arange(len(s)), cut):
        z = s.iloc[q]
        e = z.k_track.sum()/z.n_track.sum()
        xs.append(np.average(z.adc_med, weights=z.n_track)); ys.append(e)
        es.append(np.sqrt(e*(1-e)/z.n_track.sum()))
        lab.append(f"{z.adc_med.min():.0f}–{z.adc_med.max():.0f}")
    ax.errorbar(xs, ys, yerr=es, color=C1, lw=2.2, marker="o", ms=8,
                mfc=C1, mec=SURFACE, mew=1.4, capsize=3, zorder=3)
    for x, y, t in zip(xs, ys, ys):
        ax.annotate(f"{y:.3f}", (x, y), textcoords="offset points",
                    xytext=(0, 11), ha="center", fontsize=9, color=INK2)
    ax.set_xlabel("track-weighted mean pad median ADC")
    ax.set_ylabel("efficiency of that quintile")
    ax.set_ylim(0.70, 0.99)
    ax.set_title("Quintiles of pads by pulse height", loc="left", color=INK,
                 fontweight="bold", fontsize=10)
    fig.suptitle("P2_OUT: pad efficiency is set by the pad's pulse height — "
                 "the best pads already reach the DREAM number",
                 x=0.008, ha="left", fontsize=11.5, fontweight="bold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(f"{FIG}/p2out_headline.png", dpi=170); plt.close(fig)


def fig_spectra_map(g):
    z = np.load(os.path.join(DATA, "adc_vs_ch_run_46.npz"))
    H = z["adc_vs_ch"].astype(float)
    ctr = np.arange(128)*8 + 4
    u = g[g.use]

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0),
                             gridspec_kw={"width_ratios": [1.15, 1, 1]})
    ax = axes[0]; F._style(ax)
    s = u.sort_values("eff"); n3 = len(s)//3
    for q, col, lab in ((s.iloc[:n3], C2, "lowest third"),
                        (s.iloc[-n3:], C1, "highest third")):
        h = H[q.vmm.to_numpy(), q.ch.to_numpy()].sum(0); h = h/h.sum()
        ax.step(ctr, h, where="mid", color=col, lw=1.9, zorder=3,
                label=f"{lab}   eff {q.k_track.sum()/q.n_track.sum():.3f}")
    ax.set_xlim(0, 400); ax.set_xlabel("VMM ADC")
    ax.set_ylabel("fraction of hits")
    ax.legend(loc="upper right", fontsize=8.5, labelcolor=INK2)
    ax.annotate("same threshold onset", xy=(56, 0.004), xytext=(120, 0.028),
                color=INK2, fontsize=9,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.0))
    ax.set_title("Gain differs, threshold does not", loc="left", color=INK,
                 fontweight="bold", fontsize=10)

    for ax, col, cmap, lab, vlim in ((axes[1], "adc_med", "Blues",
                                      "pad median ADC", (90, 270)),
                                     (axes[2], "eff", "Blues",
                                      "pad efficiency", (0.7, 1.0))):
        F._style(ax)
        sc = ax.scatter(u.x, u.y, c=u[col], s=64, cmap=cmap,
                        norm=Normalize(*vlim), edgecolor=SURFACE,
                        linewidth=0.8, zorder=3)
        ax.set_aspect("equal"); ax.set_xlabel("x  [mm]")
        ax.set_title(lab, loc="left", color=INK, fontweight="bold", fontsize=10)
        cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
        cb.outline.set_edgecolor(MUTED); cb.ax.tick_params(labelsize=8)
    axes[1].set_ylabel("y  [mm]")
    fig.suptitle("The pulse-height map and the efficiency map are the same map",
                 x=0.008, ha="left", fontsize=11.5, fontweight="bold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(f"{FIG}/p2out_spectra_map.png", dpi=170); plt.close(fig)


def main():
    g = get()
    fig_main(g); fig_spectra_map(g)
    print("wrote p2out_headline.png, p2out_spectra_map.png")


if __name__ == "__main__":
    main()
