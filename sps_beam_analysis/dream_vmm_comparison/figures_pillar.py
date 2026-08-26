#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figures_pillar.py -- the dead areas of the P2 stack, from a 6 mm support
pillar down to the 0.75 mm bulk pillars.

  pill_1_bigpillar.png  the direct image: one Ohm-scale support pillar, mapped
                        at 0.2 mm from reference tracks alone
  pill_2_kplane.png     the whole k-plane of the efficiency map, and the six
                        first-order spots of a triangular lattice in it
  pill_3_avgpillar.png  the average bulk pillar: every one in the beam spot,
                        folded onto one cell
  pill_4_resolution.png two independent measurements of the reference's
                        pointing resolution, and the MTF they imply
  pill_5_cost.png       what the pillars cost, and why it is not the same on
                        the three stations
  pill_6_mask.png       the dead areas that CAN be masked, and what masking
                        them does to the stack

Colour as everywhere else in this package: the three stations are the
categorical triple, mapped efficiency is the deck's purple, mapped amplitude
its amber.
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

import figures as F
import pillar_stats as P
import pillar_geom as G

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
DATA = os.path.join(HERE, "data")

SC = {"P2_IN": F.C1, "P2_MID": F.C2, "P2_OUT": F.C3}
SHORT = {"P2_IN": "IN", "P2_MID": "MID", "P2_OUT": "OUT"}
STATIONS = P.STATIONS

PURPLE = LinearSegmentedColormap.from_list(
    "eff", ["#b9a3da", "#a288cc", "#8a6bbc", "#7250a4", "#583586", "#3a1c60"])
AMBER = LinearSegmentedColormap.from_list(
    "amp", ["#f3d9a4", "#eec27a", "#e3a54e", "#cf8730", "#ad6a1e", "#7d4a12"])

RC = {"font.size": 9.2, "axes.labelsize": 9.8, "axes.titlesize": 10.5,
      "xtick.labelsize": 8.8, "ytick.labelsize": 8.8,
      "figure.facecolor": F.SURFACE, "axes.facecolor": F.SURFACE,
      "savefig.facecolor": F.SURFACE, "axes.edgecolor": F.INK2,
      "text.color": F.INK, "axes.labelcolor": F.INK,
      "xtick.color": F.INK2, "ytick.color": F.INK2}


def _style(ax, grid=True):
    if grid:
        ax.grid(True, color=F.GRID, linewidth=0.7)
        ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def load(run=P.RUN):
    n = json.load(open(os.path.join(DATA, f"pillar_numbers_{run}.json")))
    a = np.load(os.path.join(DATA, f"pillar_arrays_{run}.npz"))
    return n, a


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"  {name}  ({os.path.getsize(p) / 1024:.0f} kB)")


# --------------------------------------------------------------------------- #
def _rebin(a, g):
    nx, ny = a.shape[0] // g * g, a.shape[1] // g * g
    return a[:nx, :ny].reshape(nx // g, g, ny // g, g).sum((1, 3))


def fig_bigpillar(num, arr):
    """One 6.15 mm bulk support pillar, imaged.  Nothing is folded or fitted
    here: this is where the reference tracks landed and whether the chamber
    responded."""
    fig, axes = plt.subplots(1, 4, figsize=(14.6, 3.9),
                             gridspec_kw={"width_ratios": [1, 1, 1, 1.35],
                                          "wspace": 0.52})
    R = 0.5 * num["big_pillar"][STATIONS[0]]["dia"]
    for ax, s in zip(axes[:3], STATIONS):
        g = 6                                    # 1.2 mm display bins
        n = _rebin(arr[f"bigmap_n_{s}"].astype(float), g)
        k = _rebin(arr[f"bigmap_k_{s}"].astype(float), g)
        with np.errstate(invalid="ignore"):
            e = np.where(n >= 40, k / np.maximum(n, 1), np.nan)
        x = arr[f"bigmap_x_{s}"][:n.shape[0] * g].reshape(-1, g).mean(1)
        y = arr[f"bigmap_y_{s}"][:n.shape[1] * g].reshape(-1, g).mean(1)
        im = ax.pcolormesh(x, y, e.T, cmap=PURPLE, vmin=0.0, vmax=1.0,
                           shading="nearest")
        ax.add_patch(plt.Circle((0, 0), R, fill=False, ec=F.INK, lw=1.7,
                                ls="--"))
        ax.set_aspect("equal")
        ax.set_xlim(-9.5, 9.5)
        ax.set_ylim(-9.5, 9.5)
        ax.set_title(f"{SHORT[s]}   (pale = dead)", loc="left", color=SC[s],
                     fontweight="bold", fontsize=9.4)
        ax.set_xlabel("x $-$ pillar [mm]")
        if s == STATIONS[0]:
            ax.set_ylabel("y $-$ pillar [mm]")
        bp = num["big_pillar"][s]
        ax.text(0.035, 0.035, f"core {bp['eff_core']:.3f}",
                transform=ax.transAxes, color=F.SURFACE, fontsize=8.6,
                fontweight="bold",
                bbox=dict(fc=F.INK, ec="none", pad=2.2, alpha=0.8))
    cax = axes[2].inset_axes([1.07, 0.0, 0.055, 1.0])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.set_title("eff.", fontsize=8.2, color=F.INK2, pad=5)
    cb.ax.tick_params(labelsize=8)

    ax = axes[3]
    for s in STATIONS:
        r = arr[f"big_r_{s}"]
        n, k = arr[f"big_n_{s}"].astype(float), arr[f"big_k_{s}"].astype(float)
        m = n >= 30
        with np.errstate(invalid="ignore"):
            e = k / np.maximum(n, 1)
        ax.plot(r[m], e[m], color=SC[s], lw=2.0)
        ax.text(0.68, 0.36 - 0.088 * list(STATIONS).index(s), SHORT[s],
                transform=ax.transAxes, color=SC[s], fontweight="bold",
                fontsize=9.6)
    ax.axvspan(0, R, color=F.GRID, alpha=0.85, zorder=0)
    ax.text(R * 0.5, 0.60, "the pillar,\nfrom the\ngerber", ha="center",
            va="center", color=F.INK2, fontsize=8.4)
    ax.set_xlim(0, 9.5)
    ax.set_ylim(-0.02, 1.16)
    ax.set_xlabel("distance from the pillar centre [mm]")
    ax.set_ylabel("efficiency")
    ax.set_title(f"$\\varnothing${2 * R:.2f} mm, and the chamber is blind "
                 f"across it", loc="left", color=F.INK2, fontsize=9.4)
    _style(ax)
    fig.suptitle("A support pillar, imaged by the reference telescope alone",
                 x=0.005, ha="left", fontsize=12.5, fontweight="bold", y=1.05)
    save(fig, "pill_1_bigpillar.png")


# --------------------------------------------------------------------------- #
def fig_kplane(num, arr):
    """The whole k-plane, then the lattice scan.  Nothing about a pitch is
    assumed anywhere in this figure."""
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.0),
                             gridspec_kw={"width_ratios": [1.1, 1, 1.05],
                                          "wspace": 0.40})
    s0 = "P2_OUT"
    ax = axes[0]
    kp = arr[f"kplane_{s0}"]
    ka = arr["kplane_axis"]
    im = ax.pcolormesh(ka, ka, (100 * kp).T, cmap=PURPLE, shading="nearest",
                       vmin=0, vmax=np.percentile(100 * kp, 99.97))
    d = num["stations"][s0]["lattice"]["d"]
    g = 4 * np.pi / (np.sqrt(3) * d)
    for a in range(0, 360, 60):
        ax.plot(g * np.cos(np.radians(a + 90)), g * np.sin(np.radians(a + 90)),
                "o", mfc="none", mec=F.C2, mew=1.7, ms=15)
    ax.set_aspect("equal")
    ax.set_xlabel("$k_x$ [rad/mm]")
    ax.set_ylabel("$k_y$ [rad/mm]")
    ax.set_title(f"{SHORT[s0]}: |modulation| over the k-plane", loc="left",
                 color=F.INK2, fontsize=9.4)
    cax = ax.inset_axes([1.05, 0.0, 0.05, 1.0])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("amplitude [%]", fontsize=8.4)
    cb.ax.tick_params(labelsize=8)
    ax.text(0.03, 0.965, "six spots, 60$^\\circ$ apart",
            transform=ax.transAxes, color=F.C2, fontsize=8.8,
            fontweight="bold", va="top")
    ax.text(0.03, 0.055, "the faint diagonals are the reference's own\n"
                         "1 mm strip pitch \u2014 a different frame, 60$^\\circ$ away",
            transform=ax.transAxes, color=F.INK2, fontsize=7.8, va="bottom")

    ax = axes[1]
    dg = arr["dscan_d"]
    for i, s in enumerate(STATIONS):
        ax.plot(dg, 100 * arr[f"dscan_amp_{s}"], color=SC[s], lw=2.0)
        b = num["stations"][s]["lattice"]
        ax.text(0.035, 0.94 - 0.082 * i, f"{SHORT[s]}   d = {b['d']:.4f} mm",
                transform=ax.transAxes, color=SC[s], fontweight="bold",
                fontsize=9.0)
    ax.axvline(3.0, color=F.INK2, lw=1.0, ls=":")
    ax.set_xlabel("triangular lattice spacing $d$ [mm]")
    ax.set_ylabel("amplitude modulation per direction [%]")
    ax.set_xlim(dg[0], dg[-1])
    ax.set_ylim(0, max(100 * arr[f"dscan_amp_{s}"].max() for s in STATIONS)
                * 1.55)
    ax.text(3.0, ax.get_ylim()[1] * 0.045, " 3.000 mm", color=F.INK2,
            fontsize=8.6)
    ax.set_title("scanning the lattice, not a wavevector", loc="left",
                 color=F.INK2, fontsize=9.4)
    _style(ax)

    ax = axes[2]
    lbl, val, err, col = [], [], [], []
    for s in STATIONS:
        for what, tag in (("eff", "efficiency"), ("amp", "amplitude")):
            b = num["stations"][s][what]
            lbl.append(f"{SHORT[s]}  {tag}")
            val.append(100 * b["A1"])
            err.append(100 * b["sigma_A"])
            col.append(SC[s])
    y = np.arange(len(lbl))
    ax.barh(y, val, xerr=err, color=col, height=0.62,
            error_kw=dict(ecolor=F.INK2, lw=1.2, capsize=2.5))
    for i, (v, e) in enumerate(zip(val, err)):
        ax.text(v + e + 0.4, i, f"{v / e:.0f}$\\sigma$", va="center",
                color=F.INK2, fontsize=8.8, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(lbl, fontsize=8.4)
    ax.invert_yaxis()
    ax.set_xlabel("first-order modulation amplitude [%]")
    ax.set_xlim(0, max(val) * 1.32)
    ax.set_title("against a null taken off the lattice", loc="left",
                 color=F.INK2, fontsize=9.4)
    _style(ax)
    ax.spines["left"].set_visible(False)
    sg = [num["stations"][s]["eff"]["sig1"] for s in STATIONS]
    fig.suptitle(f"The bulk pillars are a lattice, and the lattice is in the "
                 f"data at {min(sg):.0f}$-${max(sg):.0f}$\\sigma$", x=0.005,
                 ha="left", fontsize=12.5, fontweight="bold", y=1.045)
    save(fig, "pill_2_kplane.png")


# --------------------------------------------------------------------------- #
def fig_avgpillar(num, arr):
    """Every pillar in the beam spot, folded onto one cell."""
    fig, axes = plt.subplots(1, 5, figsize=(15.0, 3.9),
                             gridspec_kw={"width_ratios":
                                          [1, 1, 1, 1.3, 1.3],
                                          "wspace": 0.58})
    for ax, s in zip(axes[:3], STATIONS):
        e = arr[f"fold_edges_eff_{s}"]
        v = arr[f"fold_eff_{s}"]
        n = arr[f"foldn_eff_{s}"]
        plateau = num["stations"][s]["eff"]["deficit"]["plateau"]
        v = np.where(n >= 400, v / plateau, np.nan)
        im = ax.pcolormesh(e, e, v.T, cmap=PURPLE, shading="flat",
                           vmin=0.40, vmax=1.0)
        b = num["stations"][s]["eff"]
        ax.add_patch(plt.Circle((0, 0), b["solve"]["a"], fill=False,
                                ec=F.SURFACE, lw=1.6, ls="--"))
        ax.set_aspect("equal")
        ax.set_xlabel("x $-$ pillar [mm]")
        if s == STATIONS[0]:
            ax.set_ylabel("y $-$ pillar [mm]")
        ax.set_title(SHORT[s], loc="left", color=SC[s], fontweight="bold")
    cax = axes[2].inset_axes([1.07, 0.0, 0.055, 1.0])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.set_title("eff. /\nplateau", fontsize=8.0, color=F.INK2, pad=5)
    cb.ax.tick_params(labelsize=8)

    for ax, what, lab, ttl, ylo in (
            (axes[3], "eff", "efficiency / plateau", "the average pillar", 0.10),
            (axes[4], "amp", "charge / plateau",
             "and what it takes out of the charge", 0.50)):
        for s in STATIONS:
            r = arr[f"prof_r_{what}_{s}"]
            v = arr[f"prof_v_{what}_{s}"] / np.nanmean(
                arr[f"prof_v_{what}_{s}"][arr[f"prof_r_{what}_{s}"] > 1.35])
            # the charge is only defined for tracks the pad actually recorded,
            # so at the very centre it is the surviving high-charge tail and
            # reads too high; that bin is not plotted and the dip below is a
            # lower bound on the real charge loss
            m = r > (0.14 if what == "amp" else 0.0)
            ax.plot(r[m], v[m], color=SC[s], lw=2.0)
        ax.set_xlabel("distance from the pillar centre [mm]")
        ax.set_ylabel(lab)
        ax.set_xlim(0, 1.73)
        ax.set_ylim(ylo, 1.09)
        _style(ax)
        ax.set_title(ttl, loc="left", color=F.INK2, fontsize=9.4)
    axes[4].text(0.03, 0.06, "recorded charge only, so this is a\n"
                             "LOWER bound on what a pillar takes",
                 transform=axes[4].transAxes, color=F.INK2, fontsize=7.8,
                 va="bottom")
    for i, s in enumerate(STATIONS):
        axes[3].text(0.66, 0.30 - 0.088 * i, SHORT[s],
                     transform=axes[3].transAxes, color=SC[s],
                     fontweight="bold", fontsize=9.6)
    fig.suptitle("No single bulk pillar has enough tracks on it \u2014 all "
                 "\u224820 000 of them do", x=0.005, ha="left", fontsize=12.5,
                 fontweight="bold", y=1.05)
    save(fig, "pill_3_avgpillar.png")


# --------------------------------------------------------------------------- #
def fig_resolution(num, arr):
    """Two independent measurements of the same thing."""
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 3.9),
                             gridspec_kw={"wspace": 0.34})

    ax = axes[0]
    s0 = "P2_OUT"
    c = arr[f"edge_c_rw_{s0}"]
    h = arr[f"edge_h_rw_{s0}"].astype(float)
    fcur = arr[f"edge_f_rw_{s0}"]
    b = num["stations"][s0]["edge"]["rw"]
    lo, hi = b["width"] / 2 - 1.6, b["width"] / 2 + 1.6
    m = (c > lo) & (c < hi)
    ax.step(c[m], h[m], where="mid", color=F.INK2, lw=1.0, alpha=0.85)
    ax.plot(c[m], fcur[m], color=SC[s0], lw=2.4)
    ax.set_xlim(lo, hi)
    ax.set_ylim(0, max(np.nanmax(fcur[m]), np.nanmax(h[m])) * 1.30)
    ax.set_xlabel("track $-$ pad centre, along the pad [mm]")
    ax.set_ylabel("one-pad clusters / 0.02 mm")
    ax.set_title(f"{SHORT[s0]}: the edge of the pad box", loc="left",
                 color=F.INK2, fontsize=9.4)
    ax.text(0.035, 0.95, f"box {b['width']:.2f} mm\n"
                         f"edge $\\sigma$ = {b['sigma']:.3f} mm",
            transform=ax.transAxes, color=F.INK, fontsize=9.6,
            fontweight="bold", va="top")
    ax.text(0.035, 0.055, "the fit is to the EDGE; the interior carries\n"
                          "the pad fan and the beam profile",
            transform=ax.transAxes, color=F.INK2, fontsize=8.4, va="bottom")
    _style(ax)

    ax = axes[1]
    lbl, val, col, hat = [], [], [], []
    for s in STATIONS:
        b = num["stations"][s]
        for tag, v in (("box\nedge", float(np.mean(
                            [b["edge"][a]["sigma"] for a in b["edge"]]))),
                       ("lattice\nharmonics", b["eff"]["solve"]["sigma"])):
            lbl.append(f"{SHORT[s]}\n{tag}")
            val.append(v)
            col.append(SC[s])
            hat.append("" if "box" in tag else "///")
    x = np.arange(len(val))
    for i in range(len(val)):
        ax.bar(x[i], val[i], color=col[i], width=0.68, hatch=hat[i],
               edgecolor=F.SURFACE, linewidth=1.4)
        ax.text(x[i], val[i] + 0.008, f"{val[i]:.2f}", ha="center",
                color=F.INK2, fontsize=8.8, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(lbl, fontsize=7.4)
    ax.set_ylabel("reference pointing $\\sigma$ [mm]")
    ax.set_ylim(0, max(val) * 1.45)
    ax.set_title("two methods, nothing in common but the tracks", loc="left",
                 color=F.INK2, fontsize=9.4)
    _style(ax)

    ax = axes[2]
    per = np.linspace(0.6, 14.0, 400)
    k = 2 * np.pi / per
    for s in STATIONS:
        sg = num["stations"][s]["eff"]["solve"]["sigma"]
        ax.plot(per, np.exp(-0.5 * k ** 2 * sg ** 2), color=SC[s], lw=2.0)
    for v, lab in ((3.0, "bulk pillars\n3.00 mm"), (12.0, "pad pitch\n12 mm")):
        ax.axvline(v, color=F.INK2, lw=1.0, ls=":")
        ax.text(v * 1.06, 0.07, lab, color=F.INK2, fontsize=8.4)
    ax.set_xlabel("period [mm]")
    ax.set_ylabel("contrast retained")
    ax.set_xlim(0.6, 14)
    ax.set_ylim(0, 1.08)
    ax.set_title("what a 0.2 mm pointing error costs at each scale",
                 loc="left", color=F.INK2, fontsize=9.4)
    _style(ax)
    sg = [num["stations"][s]["eff"]["solve"]["sigma"] for s in STATIONS]
    fig.suptitle(f"The reference points to {min(sg):.2f}$-${max(sg):.2f} mm at "
                 f"the P2 planes \u2014 measured two ways", x=0.005, ha="left",
                 fontsize=12.5, fontweight="bold", y=1.045)
    save(fig, "pill_4_resolution.png")


# --------------------------------------------------------------------------- #
def fig_cost(num, arr, others=()):
    """What the lattice costs, and why the three stations differ."""
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 3.9),
                             gridspec_kw={"wspace": 0.34})

    ax = axes[0]
    plat = [num["stations"][s]["eff"]["deficit"]["plateau"] for s in STATIONS]
    mean = [num["stations"][s]["eff"]["deficit"]["mean"] for s in STATIONS]
    for i, s in enumerate(STATIONS):
        ax.bar(i, mean[i], color=SC[s], width=0.62)
        ax.bar(i, plat[i] - mean[i], bottom=mean[i], color=SC[s], alpha=0.30,
               width=0.62, edgecolor=F.SURFACE, linewidth=1.6)
        ax.text(i, plat[i] + 0.005, f"$-${100 * (plat[i] - mean[i]):.1f} pp",
                ha="center", color=F.INK2, fontsize=9.2, fontweight="bold")
        ax.text(i, mean[i] - 0.032, f"{mean[i]:.3f}", ha="center",
                color=F.SURFACE, fontsize=9.2, fontweight="bold")
    ax.set_xticks(range(3))
    ax.set_xticklabels([SHORT[s] for s in STATIONS])
    ax.set_ylim(0.80, 1.028)
    ax.set_ylabel("efficiency")
    ax.set_title("solid: measured.  faint: what the pillars take",
                 loc="left", color=F.INK2, fontsize=9.4)
    _style(ax)

    ax = axes[1]
    for i, s in enumerate(STATIONS):
        b = num["stations"][s]["eff"]
        r = arr[f"prof_r_eff_{s}"]
        v = arr[f"prof_v_eff_{s}"]
        ax.plot(r, 1 - v / b["deficit"]["plateau"], color=SC[s], lw=2.0)
        ax.text(0.56, 0.94 - 0.088 * i,
                f"{SHORT[s]}  $\\varnothing_\\mathrm{{eff}}$ "
                f"{2 * b['solve']['a']:.2f} mm", transform=ax.transAxes,
                color=SC[s], fontweight="bold", fontsize=9.0)
    ax.set_xlabel("distance from the pillar centre [mm]")
    ax.set_ylabel("efficiency lost")
    ax.set_xlim(0, 1.73)
    ax.set_ylim(0, 0.95)
    ax.set_title("the footprint is not the pillar", loc="left", color=F.INK2,
                 fontsize=9.4)
    _style(ax)

    ax = axes[2]
    sets = [(num, arr, "o", "s", 1.0)]
    for n2, a2 in others:
        sets.append((n2, a2, "^", "D", 0.6))
    lo, hi = [], []
    for n2, a2, m1, m2, al in sets:
        for s in STATIONS:
            b = n2["stations"][s]
            eff = b["eff"]["deficit"]["mean"]
            ax.plot(eff, 2 * b["eff"]["solve"]["a"], m1, ms=11, color=SC[s],
                    alpha=al)
            ax.plot(eff, 2 * b["amp"]["solve"]["a"], m2, ms=9, mfc="none",
                    mec=SC[s], mew=1.8, alpha=al)
            lo.append(2 * b["amp"]["solve"]["a"])
            hi.append(2 * b["eff"]["solve"]["a"])
    ax.axhspan(min(lo), max(lo), color=F.GRID, alpha=0.75, zorder=0)
    ax.text(0.03, (max(lo) - 0.35) / (1.22 - 0.35) + 0.03,
            f"charge: {min(lo):.2f}\u2013{max(lo):.2f} mm, everywhere",
            transform=ax.transAxes, ha="left", color=F.INK2, fontsize=8.4)
    ax.set_xlabel("station efficiency")
    ax.set_ylabel("effective dead diameter [mm]")
    ax.set_ylim(0.35, 1.22)
    ax.set_title("filled: efficiency footprint.  open: charge footprint",
                 loc="left", color=F.INK2, fontsize=9.0)
    ax.text(0.03, 0.06, "five of the six working points sit at "
                        f"{min(hi):.2f}\u2013{max(sorted(hi)[:-1]):.2f} mm;\n"
                        "the sixth is the one station that is only 89 % "
                        "efficient\n(circles/squares: this run, "
                        "triangles/diamonds: highstat_eff_1)",
            transform=ax.transAxes, color=F.INK2, fontsize=7.8, va="bottom")
    _style(ax)
    cov = 100 * float(np.mean([num["stations"][s]["amp"]["solve"]["coverage"]
                               for s in STATIONS]))
    fig.suptitle(f"A {cov:.1f} % dead area, and what it actually costs each "
                 f"station", x=0.005, ha="left", fontsize=12.5,
                 fontweight="bold", y=1.045)
    save(fig, "pill_5_cost.png")


# --------------------------------------------------------------------------- #
def fig_mask(num, arr, run):
    """The dead areas that can be masked -- and what masking them buys."""
    z, _ = P.load(run)
    geom = G.load()
    fig, axes = plt.subplots(1, 4, figsize=(15.2, 3.9),
                             gridspec_kw={"width_ratios": [1, 1, 1, 1.6],
                                          "wspace": 0.62})
    for ax, s in zip(axes[:3], STATIONS):
        f = P.fine(z, s)
        g = 10                                   # 2 mm display bins
        n = _rebin(f["n"], g)
        k = _rebin(f["k"], g)
        with np.errstate(invalid="ignore"):
            e = np.where(n >= 60, k / np.maximum(n, 1), np.nan)
        x = f["xc"][:n.shape[0] * g].reshape(-1, g).mean(1)
        y = f["yc"][:n.shape[1] * g].reshape(-1, g).mean(1)
        ok = np.where(n >= 60)
        im = ax.pcolormesh(x, y, e.T, cmap=PURPLE, shading="nearest",
                           vmin=0.75, vmax=1.0)
        for tag, col in (("big", F.C2), ("medium", F.INK2)):
            p, dd = geom[f"{tag}_xy"], float(geom[f"{tag}_dia"][0])
            for px, py in p:
                if x[ok[0].min()] < px < x[ok[0].max()] \
                        and y[ok[1].min()] < py < y[ok[1].max()]:
                    ax.add_patch(plt.Circle((px, py), dd / 2, fill=False,
                                            ec=col, lw=2.0))
        ax.set_xlim(x[ok[0].min()] - 4, x[ok[0].max()] + 4)
        ax.set_ylim(y[ok[1].min()] - 4, y[ok[1].max()] + 4)
        ax.set_aspect("equal")
        ax.set_xlabel("x [mm]")
        if s == STATIONS[0]:
            ax.set_ylabel("y [mm]")
        m = num["stations"][s]["mask"]
        ax.set_title(f"{SHORT[s]}   {m['n_dead']} dead, {m['n_weak']} weak",
                     loc="left", color=SC[s], fontweight="bold", fontsize=9.4)
    cax = axes[2].inset_axes([1.07, 0.0, 0.055, 1.0])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.set_title("eff.", fontsize=8.2, color=F.INK2, pad=5)
    cb.ax.tick_params(labelsize=8)

    ax = axes[3]
    a, b = num["stack"]["all"], num["stack"]["masked"]
    rows = [("efficiency IN", a["eff_P2_IN"], b["eff_P2_IN"], "%"),
            ("efficiency MID", a["eff_P2_MID"], b["eff_P2_MID"], "%"),
            ("efficiency OUT", a["eff_P2_OUT"], b["eff_P2_OUT"], "%"),
            ("3-of-3 tracks", a["frac_3of3"], b["frac_3of3"], "%"),
            ("exit core x [mm]", a["dpos_x"]["sigma_iqr"],
             b["dpos_x"]["sigma_iqr"], ""),
            ("exit core y [mm]", a["dpos_y"]["sigma_iqr"],
             b["dpos_y"]["sigma_iqr"], ""),
            ("angle core x [mrad]", a["dang_x"]["sigma_iqr"],
             b["dang_x"]["sigma_iqr"], "")]
    for i, (lab, v0, v1, unit) in enumerate(rows):
        rel = (v1 - v0) / v0 * 100
        ax.barh(i, rel, color=F.C3 if rel > 0 else F.C2, height=0.58)
        f0 = f"{100 * v0:.2f}" if unit == "%" else f"{v0:.3f}"
        f1 = f"{100 * v1:.2f}" if unit == "%" else f"{v1:.3f}"
        ax.text(rel + (0.05 if rel >= 0 else -0.05), i,
                f"{f0} \u2192 {f1}", va="center",
                ha="left" if rel >= 0 else "right", color=F.INK2, fontsize=8.6)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=8.6)
    ax.yaxis.tick_right()
    ax.invert_yaxis()
    ax.axvline(0, color=F.INK2, lw=1.0)
    ax.set_xlabel("change when the maskable dead areas are cut [%]")
    ax.set_xlim(-1.0, 2.4)
    ax.set_title(f"{100 * (1 - num['stack']['keep_frac']):.2f} % of tracks "
                 f"removed", loc="left", color=F.INK2, fontsize=9.4)
    _style(ax)
    ax.spines["left"].set_visible(False)
    fig.suptitle("The maskable dead area is one support pillar; the rest is "
                 "everywhere", x=0.005, ha="left", fontsize=12.5,
                 fontweight="bold", y=1.045)
    save(fig, "pill_6_mask.png")


# --------------------------------------------------------------------------- #
def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=P.RUN)
    ap.add_argument("--also", action="append", default=["highstat_eff_1"],
                    help="extra runs to overlay on the gain-margin panel")
    a = ap.parse_args()
    num, arr = load(a.run)
    others = []
    for r in a.also:
        if r == a.run:
            continue
        try:
            others.append(load(r))
        except FileNotFoundError:
            print(f"  (no products for {r}, skipped)")
    with plt.rc_context(RC):
        fig_bigpillar(num, arr)
        fig_kplane(num, arr)
        fig_avgpillar(num, arr)
        fig_resolution(num, arr)
        fig_cost(num, arr, others)
        fig_mask(num, arr, a.run)


if __name__ == "__main__":
    main()
