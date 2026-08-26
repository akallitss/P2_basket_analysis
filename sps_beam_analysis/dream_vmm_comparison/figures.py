#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figures.py -- figures for the per-pad ADC / efficiency study (run_46).

Reads data/pad_adc_run_46.csv (written by pad_adc.py) and writes figures/*.png.
"""

import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
DATA = os.path.join(HERE, "data")
os.makedirs(FIG, exist_ok=True)

# --- design tokens (validated: see notes in report) ------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a887f"
GRID = "#e6e5e0"
# sequential ramp for the discriminator threshold (ordinal): light -> dark
SDT_RAMP = {224: "#84b4e8", 256: "#4a8ad6", 288: "#2765b4", 300: "#11406c"}
# categorical slots 1..3
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"

STATIONS = ["P2_IN", "P2_MID", "P2_OUT"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.8, "grid.color": GRID, "grid.linewidth": 0.6,
    "font.size": 9, "axes.titlesize": 10, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
})


def load():
    """Per-pad table: TRACKED pulse height (the faithful quantity) joined to the
    untracked proxy, the chip threshold and the pad geometry."""
    u = pd.read_csv(os.path.join(DATA, "pad_adc_run_46.csv"))
    t = pd.read_csv(os.path.join(DATA, "pad_adc_tracked_run_46.csv"))
    t = t[["station", "vmm", "ch", "n", "k", "eff", "adc_med", "adc_med_dnl",
           "adc_p25", "adc_p75", "adc_n", "adc_med_untracked"]]
    geo = u[["station", "vmm", "ch", "x", "y", "radius", "pad_area",
             "n_hits", "masked"]]
    df = t.merge(geo, on=["station", "vmm", "ch"], how="left")

    thr = json.load(open(os.path.join(DATA, "thresholds_run_46.json")))
    sdt = {int(v): tt for st in thr for v, tt in thr[st].items()}
    df["sdt"] = df.vmm.map(sdt)
    df["dead"] = df.n_hits == 0
    df["n_track"] = df["n"]
    df["k_track"] = df["k"]
    # `adc_med` below always means the DNL-clean tracked median
    df["adc_med"] = df["adc_med_dnl"]
    df["use"] = (df.n_track >= 500) & (df.adc_n >= 50) & df.adc_med.notna() \
        & ~df.dead & ~df.masked
    df["illuminated"] = df.use | (df.dead & (df.n_track >= 500))
    return df


def _style(ax):
    ax.grid(True, axis="both", alpha=0.7, zorder=0)
    ax.set_axisbelow(True)


def fig_eff_vs_adc(df):
    """The headline: per-pad efficiency against per-pad median pulse height."""
    XLO, XHI = 60, 285
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.1), sharey=True, sharex=True)
    for ax, st in zip(axes, STATIONS):
        g = df[(df.station == st) & df.use]
        _style(ax)
        multi = g.sdt.nunique() > 1
        for sdt_v, q in g.groupby("sdt"):
            ax.scatter(q.adc_med.clip(XLO, XHI), q.eff,
                       s=8 + 40 * np.sqrt(q.n_track / q.n_track.max()),
                       c=SDT_RAMP[int(sdt_v)], edgecolor=SURFACE, linewidth=1.0,
                       alpha=0.95, zorder=3,
                       label=f"sdt {int(sdt_v)}" if multi else None)
        # per-chip trend lines -- threshold is fixed within a chip
        for _, q in g.groupby("vmm"):
            if len(q) < 8:
                continue
            b = np.polyfit(q.adc_med, q.eff, 1)
            xs = np.linspace(max(q.adc_med.min(), XLO), min(q.adc_med.max(), XHI), 2)
            ys = np.clip(np.polyval(b, xs), 0, 1)
            ax.plot(xs, ys, color=SDT_RAMP[int(q.sdt.iloc[0])], lw=1.7,
                    alpha=0.8, zorder=2)
        off = int((g.adc_med > XHI).sum())
        r, _ = stats.pearsonr(g.adc_med, g.eff)
        gg = g.assign(a=g.adc_med - g.groupby("vmm").adc_med.transform("mean"),
                      e=g.eff - g.groupby("vmm").eff.transform("mean"))
        rw, _ = stats.pearsonr(gg.a, gg.e)
        ax.set_title(st, color=INK, loc="left", fontweight="bold")
        note = f"pooled  r = {r:+.2f}\nwithin chip  r = {rw:+.2f}"
        if off:
            note += f"\n{off} pad(s) off scale (ADC > {XHI})"
        ax.text(0.97, 0.03, note, transform=ax.transAxes, fontsize=8.5,
                color=INK2, va="bottom", ha="right", zorder=6,
                bbox=dict(fc=SURFACE, ec="none", alpha=0.85, pad=2.0))
        ax.set_xlabel("pad median pulse height  [VMM ADC]")
        if multi:
            ax.legend(loc="upper left", fontsize=8, handletextpad=0.3,
                      labelcolor=INK2, borderpad=0.2)
        else:
            ax.text(0.03, 0.955, f"all chips at sdt {int(g.sdt.iloc[0])}",
                    transform=ax.transAxes, fontsize=8.5, color=INK2, va="top")
    axes[0].set_ylabel("pad efficiency")
    axes[0].set_ylim(-0.04, 1.06)
    axes[0].set_xlim(XLO - 6, XHI + 6)
    fig.suptitle("Per-pad efficiency rises with pulse height WITHIN a chip; "
                 "across chips the threshold sets both",
                 x=0.008, ha="left", fontsize=11, fontweight="bold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"{FIG}/eff_vs_adc.png", dpi=170)
    plt.close(fig)


def fig_centered(df):
    """Chip mean removed: the gain-side trend, with the threshold divided out."""
    XLO, XHI = -80, 110
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    _style(ax)
    handles = []
    for st, col in zip(STATIONS, (C3, C2, C1)):
        g = df[(df.station == st) & df.use].copy()
        g["a"] = g.adc_med - g.groupby("vmm").adc_med.transform("mean")
        g["e"] = g.eff - g.groupby("vmm").eff.transform("mean")
        ax.scatter(g.a.clip(XLO, XHI), g.e, s=14, c=col, edgecolor=SURFACE,
                   linewidth=0.9, alpha=0.9, zorder=3)
        b = np.polyfit(g.a, g.e, 1)
        xs = np.linspace(max(g.a.min(), XLO), min(g.a.max(), XHI), 2)
        ax.plot(xs, np.polyval(b, xs), color=col, lw=2.2, zorder=4)
        r, _ = stats.pearsonr(g.a, g.e)
        # b[0] is d(efficiency)/d(ADC); x1000 -> percentage points per 10 ADC
        handles.append(Line2D([], [], color=col, lw=2, marker="o", ms=5,
                              label=f"{st}   r = {r:+.2f}   "
                                    f"{b[0]*1000:+.1f} pts / 10 ADC"))
    ax.axhline(0, color=MUTED, lw=0.8, zorder=1)
    ax.axvline(0, color=MUTED, lw=0.8, zorder=1)
    ax.set_xlim(XLO - 6, XHI + 6)
    ax.set_xlabel("pad median pulse height  −  its chip's mean   [VMM ADC]")
    ax.set_ylabel("pad efficiency  −  its chip's mean")
    ax.legend(handles=handles, loc="lower right", fontsize=8.5, labelcolor=INK2)
    ax.set_title("With the per-chip threshold removed, the two working stations\n"
                 "share one positive slope — P2_IN still runs backwards",
                 loc="left", color=INK, fontweight="bold", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(f"{FIG}/eff_vs_adc_centered.png", dpi=170)
    plt.close(fig)


def fig_spectra(df):
    """The ADC distributions themselves, low- vs high-efficiency pads."""
    z = np.load(os.path.join(DATA, "adc_vs_ch_run_46.npz"))
    H = z["adc_vs_ch"].astype(float)
    edges = np.arange(129) * 8
    ctr = edges[:-1] + 4

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9))
    for ax, st in zip(axes, ("P2_OUT", "P2_MID")):
        _style(ax)
        g = df[(df.station == st) & df.use].sort_values("eff")
        n3 = max(len(g) // 3, 1)
        lo, hi = g.iloc[:n3], g.iloc[-n3:]
        for q, col, lab in ((lo, C2, "lowest third"), (hi, C1, "highest third")):
            h = H[q.vmm.to_numpy(), q.ch.to_numpy()].sum(0)
            h = h / h.sum()
            ax.step(ctr, h, where="mid", color=col, lw=1.8, zorder=3,
                    label=f"{lab}   eff {q.k_track.sum()/q.n_track.sum():.3f}")
            med = np.interp(0.5, np.cumsum(h), ctr)
            ax.axvline(med, color=col, lw=1.0, ls=(0, (4, 3)), zorder=2)
            ax.annotate(f"median {med:.0f}", xy=(med, h.max() * 0.62),
                        xytext=(med + 42, h.max() * (0.78 if col == C2 else 0.62)),
                        color=col, fontsize=8.5,
                        arrowprops=dict(arrowstyle="-", color=col, lw=0.9))
        onset = np.flatnonzero(H[q.vmm.to_numpy(), q.ch.to_numpy()].sum(0) > 0)[0] * 8
        ax.annotate("same threshold onset", xy=(onset + 6, h.max() * 0.10),
                    xytext=(onset + 74, h.max() * 0.30), color=INK2, fontsize=8.5,
                    arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.9))
        ax.set_xlim(0, 400)
        ax.set_xlabel("VMM ADC  (8-code bins, all recorded hits)")
        ax.set_title(st, loc="left", color=INK, fontweight="bold")
        ax.legend(loc="upper right", fontsize=8.5, labelcolor=INK2)
    axes[0].set_ylabel("fraction of hits")
    fig.suptitle("Same threshold onset, smaller pulse: the low-efficiency pads "
                 "differ in GAIN, not in threshold",
                 x=0.008, ha="left", fontsize=11, fontweight="bold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(f"{FIG}/spectra_low_high.png", dpi=170)
    plt.close(fig)


def fig_maps(df):
    """Where on the chamber the pulse height and the efficiency sit."""
    from matplotlib.colors import Normalize
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 6.6))
    for j, st in enumerate(STATIONS):
        g = df[(df.station == st) & df.illuminated]
        for i, (col, cmap, lab, vmin, vmax) in enumerate((
                ("adc_med", "Blues", "pad median ADC", 60, 220),
                ("eff", "Blues", "pad efficiency", 0, 1))):
            ax = axes[i, j]
            _style(ax)
            live = g[~g.dead]
            sc = ax.scatter(live.x, live.y, c=live[col], s=52, cmap=cmap,
                            norm=Normalize(vmin, vmax), edgecolor=SURFACE,
                            linewidth=0.8, zorder=3)
            dead = g[g.dead]
            if len(dead):
                ax.scatter(dead.x, dead.y, s=52, facecolor="none",
                           edgecolor=C2, linewidth=1.4, zorder=4,
                           label=f"{len(dead)} dead")
                ax.legend(loc="upper left", fontsize=8, labelcolor=INK2)
            ax.set_aspect("equal")
            ax.set_title(f"{st} — {lab}", loc="left", color=INK, fontsize=9.5)
            cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
            cb.outline.set_edgecolor(MUTED)
            cb.ax.tick_params(labelsize=8, color=MUTED)
            if i == 1:
                ax.set_xlabel("x  [mm]")
            if j == 0:
                ax.set_ylabel("y  [mm]")
    fig.suptitle("Pulse height varies smoothly across the chamber; the dead "
                 "channels are a separate, chip-local defect",
                 x=0.008, ha="left", fontsize=11, fontweight="bold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(f"{FIG}/chamber_maps.png", dpi=170)
    plt.close(fig)


def fig_dnl(df):
    """The period-16 differential non-linearity of the VMM3a ADC."""
    z = np.load(os.path.join(DATA, "adc_vs_ch_run_46.npz"))
    A = z["adc_per_vmm"].astype(float)
    H = z["adc_vs_ch"].astype(float)
    sys.path.insert(0, "/home/dylan/PycharmProjects/DAQ_Control_VMM_Beam/vmm_qa")
    import vmm_stations as VS

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9))
    ax = axes[0]
    _style(ax)
    lo, hi = 48, 512
    for st, col in zip(("P2_MID", "P2_OUT"), (C1, C2)):
        bad = set(df[(df.station == st) & df.masked].vmm.unique())
        clean = [v for v in VS.STATION_VMMS[st] if v not in bad]
        s = A[clean].sum(0)[lo:hi]
        ph = np.arange(lo, hi) % 16
        w = np.array([s[ph == p].sum() for p in range(16)])
        w = w / w.mean()
        ax.step(np.arange(16), w, where="mid", color=col, lw=1.8, zorder=3,
                label=f"{st}   {100*s[ph==0].sum()/s.sum():.0f}% on phase 0")
    ax.axhline(1.0, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.text(15.6, 1.03, "flat", color=MUTED, fontsize=8, ha="right")
    ax.set_xlabel("ADC code mod 16")
    ax.set_ylabel("relative code width")
    ax.set_title("Phase 0 is always too wide (P2_OUT carries extra structure)",
                 loc="left", color=INK, fontweight="bold", fontsize=9.5)
    ax.legend(loc="upper center", fontsize=8.5, labelcolor=INK2)

    ax = axes[1]
    _style(ax)
    g = df[(df.station == "P2_OUT") & df.use]
    h = H[g.vmm.to_numpy(), g.ch.to_numpy()].sum(0)
    ctr8 = np.arange(128) * 8 + 4
    h16 = h.reshape(64, 2).sum(1)
    ctr16 = np.arange(64) * 16 + 8
    ax.step(ctr8, h / h.sum(), where="mid", color=MUTED, lw=1.2, zorder=2,
            label="8-code bins")
    ax.step(ctr16, h16 / h16.sum() / 2, where="mid", color=C1, lw=1.8, zorder=3,
            label="16-code bins (one full DNL period)")
    mpv = ctr16[np.argmax(h16[2:]) + 2]
    ax.axvline(mpv, color=C1, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.annotate(f"MPV {mpv:.0f}", xy=(mpv, ax.get_ylim()[1] * 0.7),
                xytext=(mpv + 55, ax.get_ylim()[1] * 0.8), color=C1,
                fontsize=9, arrowprops=dict(arrowstyle="-", color=C1, lw=1))
    ax.set_xlim(0, 420)
    ax.set_xlabel("VMM ADC")
    ax.set_ylabel("fraction of hits per code")
    ax.set_title("P2_OUT illuminated pads, DNL-clean binning", loc="left",
                 color=INK, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8.5, labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(f"{FIG}/dnl.png", dpi=170)
    plt.close(fig)


def fig_validation(df):
    """Does the untracked proxy measure the same thing as the tracked ADC?"""
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9))
    for ax, st in zip(axes, STATIONS):
        g = df[(df.station == st) & df.use].dropna(subset=["adc_med_untracked"])
        _style(ax)
        lim = (60, 300)
        ax.plot(lim, lim, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=2)
        ax.scatter(g.adc_med, g.adc_med_untracked.clip(*lim), s=16, c=C1,
                   edgecolor=SURFACE, linewidth=0.9, zorder=3)
        r, _ = stats.pearsonr(g.adc_med, g.adc_med_untracked)
        d = g.adc_med_untracked - g.adc_med
        ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_aspect("equal")
        ax.set_title(st, loc="left", color=INK, fontweight="bold")
        ax.text(0.96, 0.06, f"r = {r:+.3f}\nspread {d.std():.0f} ADC",
                transform=ax.transAxes, fontsize=8.5, color=INK2,
                ha="right", va="bottom")
        ax.set_xlabel("tracked median  [ADC]")
    axes[0].set_ylabel("untracked median  [ADC]")
    fig.suptitle("The untracked per-channel histogram measures the same pulse "
                 "height as the tracked hits — except on P2_IN",
                 x=0.008, ha="left", fontsize=11, fontweight="bold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(f"{FIG}/tracked_vs_untracked.png", dpi=170)
    plt.close(fig)


def main():
    df = load()
    fig_validation(df)
    fig_eff_vs_adc(df)
    fig_centered(df)
    fig_spectra(df)
    fig_maps(df)
    fig_dnl(df)
    print("wrote:", ", ".join(sorted(os.listdir(FIG))))


if __name__ == "__main__":
    main()
