#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figures_track.py -- how well the three P2 stations track, against truth.

Five figures, in the order the argument runs:

  track_1_pointing.png   what ONE station knows: the residual to the reference
                         track, and why it is a box and not a Gaussian
  track_2_selftrack.png  what the THREE together know: the P2-only track
                         against the reference, and how the error grows with
                         the lever arm
  track_3_illusion.png   why P2 cannot check its own tracking: 69 % of the
                         time all three stations report the same pad, so the
                         self-consistency residual is 25x too good
  track_4_maps.png       where on each chamber the tracks are found
  track_5_purity.png     how much of the efficiency is an accidental match

Colour follows the convention already in this package: the three stations are
a categorical triple (validated: worst adjacent CVD dE 9.2 deutan, all three
inside the lightness band; the green sits at 2.7:1 on the surface, so every
series is directly labelled rather than left to the legend), and anything
MAPPED uses the sequential purple that means efficiency in the deck.
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

import figures as F
import track_stats as T

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")

# station identity -- categorical, fixed order, never cycled
SC = {"P2_IN": F.C1, "P2_MID": F.C2, "P2_OUT": F.C3}
SHORT = {"P2_IN": "IN", "P2_MID": "MID", "P2_OUT": "OUT"}

# efficiency, mapped -- the deck's purple, same ramp, same meaning
PURPLE = LinearSegmentedColormap.from_list(
    "eff", ["#b9a3da", "#a288cc", "#8a6bbc", "#7250a4", "#583586", "#3a1c60"])

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


def step(ax, counts, edges, colour, label=None, lw=2.0, norm=True, **kw):
    c = np.asarray(counts, float)
    if norm and c.sum():
        c = c / c.sum() / np.diff(edges)
    ax.step(edges[:-1], c, where="post", color=colour, lw=lw, label=label, **kw)
    return c


def curve_from_samples(v, lo, hi, nb=200):
    h, e = np.histogram(v, np.linspace(lo, hi, nb + 1))
    return h, e


# --------------------------------------------------------------------------- #
def fig_pointing(z, j, st):
    """One station against the reference.  The residual is the pad, not the
    detector: a 12 mm box, the same on all three, and a cluster spread over two
    pads does not improve it."""
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9))
    e = z["res_edges"]

    for ax, coord, pitch in ((axes[0], "x", T.PITCH_X), (axes[1], "y", T.PITCH_Y)):
        for s in T.STATIONS:
            h = z[f"res_{coord}_{s}_single"] + z[f"res_{coord}_{s}_multi"]
            c = step(ax, h, e, SC[s], lw=2.0)
            # direct label at the top of each curve, staggered
            k = list(T.STATIONS).index(s)
            ax.text(0.04, 0.93 - 0.085 * k, SHORT[s], transform=ax.transAxes,
                    color=SC[s], fontweight="bold", fontsize=10)
            ax.text(0.155, 0.93 - 0.085 * k,
                    f"rms {st[s][f'res_{coord}']['rms']:.2f} mm",
                    transform=ax.transAxes, color=F.INK2, fontsize=9)
        ax.text(0.04, 0.655, "the three are superimposed:\none PCB design, one pad",
                transform=ax.transAxes, color=F.INK2, fontsize=8.6, va="top")
        # the box a 12 mm pad gives with no sharing at all
        half = pitch / 2
        ax.axvspan(-half, half, color=F.GRID, alpha=0.65, zorder=0)
        ax.text(0.0, ax.get_ylim()[1] * 0.06, f"one pad, {pitch:.1f} mm",
                ha="center", color=F.INK2, fontsize=8.4)
        ax.set_xlim(-18, 18)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.30)
        ax.set_xlabel(f"track $-$ cluster, {coord} [mm]")
        ax.set_title(f"{coord}: pitch/$\\sqrt{{12}}$ = {pitch / np.sqrt(12):.2f} mm",
                     loc="left", color=F.INK2, fontsize=9.6)
        _style(ax)
    axes[0].set_ylabel("tracks / mm (normalised)")

    # -- does a two-pad cluster do better?  no. ------------------------------ #
    ax = axes[2]
    hs = sum(z[f"res_x_{s}_single"] for s in T.STATIONS)
    hm = sum(z[f"res_x_{s}_multi"] for s in T.STATIONS)
    step(ax, hs, e, F.INK, lw=2.0)
    step(ax, hm, e, F.C2, lw=2.0)
    ss, sm = T.hstats(hs, e), T.hstats(hm, e)
    fm = hm.sum() / (hs.sum() + hm.sum())
    ax.text(0.04, 0.93, f"1 pad   rms {ss['rms']:.2f} mm",
            transform=ax.transAxes, color=F.INK, fontweight="bold", fontsize=9.6)
    ax.text(0.04, 0.845, f"$\\geq$2 pads   rms {sm['rms']:.2f} mm",
            transform=ax.transAxes, color=F.C2, fontweight="bold", fontsize=9.6)
    ax.text(0.04, 0.775, f"({fm:.0%} of clusters)",
            transform=ax.transAxes, color=F.C2, fontsize=8.8)
    ax.text(0.04, 0.655, "a centroid over two pads is not\nbetter — it is slightly worse",
            transform=ax.transAxes, color=F.INK2, fontsize=8.8, va="top")
    ax.set_xlim(-18, 18)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.48)
    ax.set_xlabel("track $-$ cluster, x [mm]")
    ax.set_title("split by cluster size, all stations",
                 loc="left", color=F.INK2, fontsize=9.6)
    _style(ax)

    fig.suptitle("What one station knows about a track", x=0.008, ha="left",
                 fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = os.path.join(FIG, "track_1_pointing.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
def fig_selftrack(sb):
    """The three stations as a tracker, judged against the reference."""
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9))

    # -- (a) offset at the exit plane ---------------------------------------- #
    ax = axes[0]
    for coord, colour in (("x", F.C1), ("y", F.C2)):
        d = sb[f"dpos_exit_{coord}"]
        h, e = curve_from_samples(d, -30, 30, 240)
        step(ax, h, e, colour, lw=1.9)
        s = T.spread(d)
        ax.text(0.03, 0.93 if coord == "x" else 0.84,
                f"$\\Delta${coord}   {s:.2f} mm core,  {d.std():.1f} mm rms",
                transform=ax.transAxes, color=colour, fontweight="bold",
                fontsize=9.6)
    ax.set_yscale("log")
    ax.set_xlim(-30, 30)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * 60)
    ax.set_xlabel("P2 track $-$ reference, at the exit plane [mm]")
    ax.set_ylabel("tracks / mm (normalised)")
    ax.set_title("where the P2-only track says the particle was",
                 loc="left", color=F.INK2, fontsize=9.6)
    _style(ax)

    # -- (b) angle ------------------------------------------------------------ #
    ax = axes[1]
    for coord, colour in (("x", F.C1), ("y", F.C2)):
        d = sb[f"dang_{coord}"]
        h, e = curve_from_samples(d, -30, 30, 240)
        step(ax, h, e, colour, lw=1.9)
        ax.text(0.03, 0.93 if coord == "x" else 0.84,
                f"$\\Delta\\theta_{coord}$   {T.spread(d):.2f} mrad core,  "
                f"{np.percentile(np.abs(d), 95):.0f} mrad at 95 %",
                transform=ax.transAxes, color=colour, fontweight="bold",
                fontsize=9.6)
    ax.set_yscale("log")
    ax.set_xlim(-30, 30)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * 60)
    # one pad of disagreement between the first and last station, over the
    # 620 mm between them, IS an angle -- and it is where the satellites sit
    q = T.PITCH_X / (T.Z["P2_OUT"] - T.Z["P2_IN"]) * 1e3
    for sgn in (-1, 1):
        ax.axvline(sgn * q, color=F.INK, lw=1.0, ls=(0, (3, 3)), alpha=0.7)
    ax.text(q + 1.2, 0.55, f"$\\pm$1 pad over the 620 mm\nlever arm = {q:.0f} mrad",
            transform=ax.get_xaxis_transform(), color=F.INK, fontsize=8.4,
            ha="left", va="center")
    ax.set_xlabel("P2 track $-$ reference, angle [mrad]")
    ax.set_title("and which way it was going", loc="left", color=F.INK2,
                 fontsize=9.6)
    ax.text(0.03, 0.70, "the satellites are the tail:\none station picking the\n"
                        "neighbouring pad",
            transform=ax.transAxes, color=F.INK2, fontsize=8.8, va="top")
    _style(ax)

    # -- (c) error against the lever arm -------------------------------------- #
    ax = axes[2]
    zz = sb["zgrid"]
    ymax = 0.0
    for coord, colour in (("x", F.C1), ("y", F.C2)):
        core = np.array([T.spread(b) for b in sb[f"curve_{coord}"]])
        p95 = np.array([np.percentile(np.abs(b), 95)
                        for b in sb[f"curve_{coord}"]])
        ax.fill_between(zz, core, p95, color=colour, alpha=0.10, lw=0)
        ax.plot(zz, core, color=colour, lw=2.2)
        ax.plot(zz, p95, color=colour, lw=1.4, ls=(0, (4, 3)))
        ymax = max(ymax, p95.max())
    ax.set_xlim(0, 3000)
    ax.set_ylim(0, ymax * 1.16)
    cx = np.array([T.spread(b) for b in sb["curve_x"]])
    px = np.array([np.percentile(np.abs(b), 95) for b in sb["curve_x"]])
    ax.text(2940, cx[-1] + ymax * 0.055, "core (x and y, superimposed)",
            color=F.INK, ha="right", fontweight="bold", fontsize=9.2)
    ax.text(2940, px[-1] - ymax * 0.03, "95th percentile", color=F.INK,
            ha="right", va="top", fontweight="bold", fontsize=9.2)
    ax.text(0.03, 0.20, "$\\Delta$x", transform=ax.transAxes, color=F.C1,
            fontweight="bold", fontsize=10)
    ax.text(0.105, 0.20, "$\\Delta$y", transform=ax.transAxes, color=F.C2,
            fontweight="bold", fontsize=10)
    ax.axvspan(940, 3000, color=F.GRID, alpha=0.55, zorder=0)
    for k, s in enumerate(T.STATIONS):
        ax.axvline(T.Z[s], color=SC[s], lw=1.2, alpha=0.8, zorder=1)
        ax.text(T.Z[s], ymax * (1.145 - 0.055 * (k % 2)), SHORT[s],
                color=SC[s], fontsize=8.8, ha="center", va="top",
                fontweight="bold")
    ax.text(1970, ymax * 0.30, "extrapolated beyond\nthe basket", ha="center",
            color=F.INK2, fontsize=8.8)
    ax.set_xlabel("z along the beam [mm]   (reference front = 0)")
    ax.set_ylabel("pointing error [mm]")
    ax.set_title("core (solid) and 95th percentile (dashed)",
                 loc="left", color=F.INK2, fontsize=9.6)
    _style(ax)

    fig.suptitle(f"What the three of them know together — "
                 f"{sb['n']:,} three-station tracks",
                 x=0.008, ha="left", fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = os.path.join(FIG, "track_2_selftrack.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
def fig_illusion(sb):
    """The headline: P2 checking itself gets an answer 25x too good."""
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9))

    # -- (a) two rulers, same events, same clusters --------------------------- #
    ax = axes[0]
    self_d, ref_d = sb["mid_self_x"], sb["mid_ref_x"]
    for d, colour, lab in ((self_d, F.C3, "checked against the other two P2"),
                           (ref_d, F.C2, "checked against the reference")):
        h, e = curve_from_samples(d - np.median(d), -20, 20, 320)
        step(ax, h, e, colour, lw=2.0)
    ax.set_yscale("log")
    ax.set_xlim(-20, 20)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * 120)
    ax.text(0.03, 0.95, f"self-checked   {T.spread(self_d):.2f} mm",
            transform=ax.transAxes, color=F.C3, fontweight="bold", fontsize=10,
            va="top")
    ax.text(0.03, 0.885, f"truth-checked   {T.spread(ref_d):.2f} mm",
            transform=ax.transAxes, color=F.C2, fontweight="bold", fontsize=10,
            va="top")
    ax.text(0.03, 0.80,
            f"the same events and the same clusters.\n"
            f"only the ruler changes — and it flatters\n"
            f"P2_MID by {T.spread(ref_d) / T.spread(self_d):.0f}$\\times$.",
            transform=ax.transAxes, color=F.INK2, fontsize=8.8, va="top")
    ax.set_xlabel("P2_MID residual [mm]")
    ax.set_ylabel("tracks / mm (normalised)")
    ax.set_title("one station, two rulers", loc="left", color=F.INK2,
                 fontsize=9.6)
    _style(ax)

    # -- (b) why: the same pad fires in all three ----------------------------- #
    ax = axes[1]
    f = sb["frac_same_pad"]
    ax.barh([0.78], [f], height=0.10, color=F.C3, edgecolor=F.SURFACE, lw=2)
    ax.barh([0.78], [1 - f], left=[f], height=0.10, color=F.GRID,
            edgecolor=F.SURFACE, lw=2)
    ax.text(f / 2, 0.78, f"{f:.0%}", ha="center", va="center", color="white",
            fontweight="bold", fontsize=13)
    ax.text(f + (1 - f) / 2, 0.78, f"{1 - f:.0%}", ha="center", va="center",
            color=F.INK2, fontweight="bold", fontsize=11)
    ax.text(0.0, 0.90, "all three stations report the SAME pad",
            color=F.C3, fontweight="bold", fontsize=10, va="center")
    ax.text(0.0, 0.60,
            "The stations are the same PCB, mounted parallel,\n"
            "and the beam is near-parallel to z.  So the pad the\n"
            "particle lands on is usually the same pad in all\n"
            "three — and a straight line through three copies of\n"
            "the same number has no residual to show.\n\n"
            "That is not resolution.  It is the same rounding\n"
            "error three times.",
            color=F.INK2, fontsize=8.8, va="top", linespacing=1.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("fraction of three-station tracks")
    ax.set_title("why the self-check cannot see anything", loc="left",
                 color=F.INK2, fontsize=9.6)
    _style(ax, grid=False)
    ax.spines["left"].set_visible(False)

    # -- (c) what a self-referenced efficiency does to the same detectors ----- #
    ax = axes[2]
    n = sb["n_station_frac"]
    ax.bar(range(4), n, color=[F.GRID, F.GRID, "#c9c7bf", F.C1],
           edgecolor=F.SURFACE, lw=2, width=0.68)
    for k, v in enumerate(n):
        ax.text(k, v + 0.015, f"{v:.1%}", ha="center", color=F.INK2,
                fontsize=9)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["none", "1 of 3", "2 of 3", "all 3"])
    ax.set_ylim(0, max(n) * 1.22)
    ax.set_xlabel("stations that found the particle")
    ax.set_ylabel("fraction of reference tracks")
    ax.set_title("of tracks inside all three chambers", loc="left",
                 color=F.INK2, fontsize=9.6)
    ax.text(0.04, 0.72,
            f"A three-station track exists for\n{n[3]:.0%} of the particles the\n"
            f"reference says went through.\n\n"
            f"{n[2]:.0%} more give a two-station\ntrack: enough for a position,\n"
            f"not for an angle.",
            transform=ax.transAxes, color=F.INK2, fontsize=8.8, va="top",
            linespacing=1.5)
    _style(ax)

    fig.suptitle("Why P2 cannot measure its own tracking", x=0.008,
                 ha="left", fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = os.path.join(FIG, "track_3_illusion.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
def fig_maps(z, st):
    """Where on each chamber the reference tracks are found."""
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.3))
    norm = Normalize(0.80, 1.0)
    sc = None
    for ax, s in zip(axes, T.STATIONS):
        xe, ye, e, n = T.effmap(z, s)
        e = np.where(n >= 60, e, np.nan)
        sc = ax.pcolormesh(xe, ye, e.T, cmap=PURPLE, norm=norm,
                           shading="flat")
        # the beam lights up a fraction of the chamber; showing the whole pad
        # plane would shrink the part that carries the answer to nothing
        ix, iy = np.where(n >= 60)
        pad = 3
        ax.set_xlim(xe[max(ix.min() - pad, 0)],
                    xe[min(ix.max() + 1 + pad, len(xe) - 1)])
        ax.set_ylim(ye[max(iy.min() - pad, 0)],
                    ye[min(iy.max() + 1 + pad, len(ye) - 1)])
        ax.set_aspect("equal")
        ax.set_title(f"{s}    {st[s]['eff']:.3f}", loc="left", color=SC[s],
                     fontweight="bold", fontsize=11)
        ax.set_xlabel("pad-frame x [mm]")
        ax.tick_params(labelsize=8)
        for side in ax.spines.values():
            side.set_color(F.GRID)
    axes[0].set_ylabel("pad-frame y [mm]")
    cb = fig.colorbar(sc, ax=axes, orientation="horizontal", fraction=0.055,
                      pad=0.14, aspect=48)
    cb.set_label("efficiency against the reference track  "
                 "(bins with $\\geq$60 tracks)", fontsize=9.4)
    cb.outline.set_visible(False)
    fig.suptitle("Where each chamber finds the particle", x=0.008, ha="left",
                 fontsize=12.5, fontweight="bold")
    out = os.path.join(FIG, "track_4_maps.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
def fig_purity(z, st, probe_r=15.0):
    """How much of the efficiency is a real hit and how much is an accident."""
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.9))
    e = z["dmin_edges"]
    ctr = T.centres(e)

    ax = axes[0]
    for s in T.STATIONS:
        h = np.asarray(z[f"dmin_{s}"], float)
        # normalised by every track in the fiducial, not by the ones that found
        # a cluster: a track with no cluster at all has dmin = infinity and
        # falls outside the histogram, and dividing it out would hide the
        # inefficiency this figure is about
        nfid = st[s]["n"]
        # per unit AREA: an accidental match is uniform in the plane, so a flat
        # tail here -- not a flat tail in dmin -- is what "accidental" looks like
        area = np.pi * (e[1:] ** 2 - e[:-1] ** 2)
        ax.step(e[:-1], h / nfid / area, where="post", color=SC[s], lw=2.0)
        k = list(T.STATIONS).index(s)
        ax.text(0.62, 0.93 - 0.09 * k, SHORT[s], transform=ax.transAxes,
                color=SC[s], fontweight="bold", fontsize=10)
    ax.axvline(probe_r, color=F.INK, lw=1.4, ls=(0, (4, 3)))
    ax.text(probe_r + 0.7, 0.5, f"probe radius {probe_r:g} mm", rotation=90,
            transform=ax.get_xaxis_transform(), color=F.INK, fontsize=8.6,
            va="center")
    ax.set_yscale("log")
    ax.set_xlim(0, 40)
    ax.set_xlabel("distance, track to nearest cluster [mm]")
    ax.set_ylabel("fraction of fiducial tracks / mm$^2$")
    ax.set_title("a real hit, or an accident?", loc="left", color=F.INK2,
                 fontsize=9.6)
    _style(ax)

    ax = axes[1]
    for s in T.STATIONS:
        h = np.asarray(z[f"dmin_{s}"], float)
        cum = np.cumsum(h) / st[s]["n"]
        ax.plot(e[1:], cum, color=SC[s], lw=2.2)
        k = list(T.STATIONS).index(s)
        ax.text(0.62, 0.30 - 0.08 * k, SHORT[s], transform=ax.transAxes,
                color=SC[s], fontweight="bold", fontsize=10)
        # the accidental slope: the plateau's rise from the probe radius out
        i15 = np.searchsorted(e[1:], probe_r)
        slope = (cum[-1] - cum[i15]) / (e[-1] - probe_r) * 10
        ax.text(0.62 + 0.10, 0.30 - 0.08 * k,
                f"+{slope * 100:.2f} % per 10 mm", transform=ax.transAxes,
                color=F.INK2, fontsize=8.8)
    ax.axvline(probe_r, color=F.INK, lw=1.4, ls=(0, (4, 3)))
    ax.set_xlim(0, 40)
    ax.set_ylim(0.5, 1.0)
    ax.axhline(1.0, color=F.GRID, lw=1.0)
    ax.set_xlabel("probe radius [mm]")
    ax.set_ylabel("efficiency")
    ax.set_title("efficiency against the one knob that moves it",
                 loc="left", color=F.INK2, fontsize=9.6)
    _style(ax)

    fig.suptitle("How much of the efficiency is real", x=0.008, ha="left",
                 fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(FIG, "track_5_purity.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
AMBER = LinearSegmentedColormap.from_list(
    "amp", ["#e0a94f", "#d4903a", "#c47829", "#ac5f1a", "#8e470f", "#6c3306"])


def fig_inpad(z, st):
    """Every track folded onto the face of the pad it pointed at.

    The question this was built to answer was whether a pad loses hits at its
    own edge, where the avalanche is split with its neighbour.  It does not:
    the efficiency is flat across the face and if anything a shade HIGHER at
    the rim, because there two pads each get a chance at the same charge.
    What the rim does show is the sharing itself -- the leading pad keeps ~5 %
    less charge and the cluster grows -- and that it is only ~5 % is exactly
    why the residual is the full pitch/sqrt(12) box.
    """
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2))
    e = z["pad_edges"]
    ref = "P2_OUT"        # one station, not an average: the three sit at
    # different efficiencies, and averaging maps that are NaN in different
    # bins manufactures structure that is not in any of them

    for ax, what, cmap, lab in (
            (axes[0], "eff", PURPLE, "efficiency"),
            (axes[1], "amp", AMBER, "mean leading-pad amplitude [ADC]")):
        v = T.padmap(z, ref, what)[1]
        lo, hi = np.nanpercentile(v, [3, 97])
        m = ax.pcolormesh(e, e, v.T, cmap=cmap, norm=Normalize(lo, hi))
        cb = fig.colorbar(m, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(lab, fontsize=9)
        cb.outline.set_visible(False)
        ax.set_aspect("equal")
        ax.set_xlim(-7, 7)
        ax.set_ylim(-7, 7)
        ax.set_xlabel("track $-$ pad centre, x [mm]")
        for side in ax.spines.values():
            side.set_color(F.GRID)
        ax.add_patch(plt.Circle((0, 0), T.PITCH_Y / 2, fill=False,
                                edgecolor=F.INK, lw=1.2, ls=(0, (4, 3))))
    axes[0].set_ylabel("track $-$ pad centre, y [mm]")
    axes[0].set_title(f"efficiency across one pad — {ref}", loc="left",
                      color=F.INK2, fontsize=9.6)
    axes[1].set_title("the charge the leading pad keeps", loc="left",
                      color=F.INK2, fontsize=9.6)
    axes[1].set_xlabel("track $-$ pad centre, x [mm]     "
                       "(dashed: half a pad)")

    # -- (c) everything relative to the middle of the pad -------------------- #
    # Three different quantities on one axis is only legitimate because each is
    # divided by its own value at the pad centre: the axis is "relative to the
    # centre", not three units pretending to be one.
    ax = axes[2]
    style = {"eff": (F.C1, "-", "efficiency"),
             "amp": (F.C2, "-", "leading-pad charge"),
             "nclus": (F.C3, "-", "pads in the cluster")}
    for what, (colour, ls, lab) in style.items():
        prof = []
        for s in T.STATIONS:
            c, v, den = T.padprofile(z, s, what)
            # normalise on the interior plateau, not on the innermost bin: at
            # small r only a handful of 0.4 mm cells contribute, so that bin is
            # an average over a few particular pads rather than over all of them
            flat = (c > 1.0) & (c < 4.5)
            prof.append(v / np.mean(v[flat]))
        y = np.mean(prof, axis=0)
        keep = c >= 1.0
        ax.plot(c[keep], y[keep], color=colour, lw=2.3, ls=ls)
        ax.text(c[-1] + 0.18, y[-1], lab, color=colour,
                fontweight="bold", fontsize=9.2, va="center")
        ax.plot([c[-1]], [y[-1]], "o", color=colour, ms=6,
                markeredgecolor=F.SURFACE, markeredgewidth=1.5)
    ax.axhline(1.0, color=F.GRID, lw=1.2)
    ax.axvline(T.PITCH_Y / 2, color=F.INK, lw=1.1, ls=(0, (4, 3)))
    ax.text(T.PITCH_Y / 2 - 0.15, 0.04, "half a pad", rotation=90,
            transform=ax.get_xaxis_transform(), color=F.INK, fontsize=8.4,
            ha="right", va="bottom")
    ax.set_xlim(0.6, 10.2)
    ax.set_xticks([1, 2, 3, 4, 5, 6])
    ax.set_xlabel("distance from the pad centre [mm]")
    ax.set_ylabel("relative to the pad centre")
    ax.set_title("mean of the three stations", loc="left", color=F.INK2,
                 fontsize=9.6)
    ax.text(0.04, 0.60, "flat inside 4 mm; the sharing\n"
                        "turns on in the outer 10 %\nof the pad face",
            transform=ax.transAxes, color=F.INK2, fontsize=8.8, va="top")
    _style(ax)

    fig.suptitle("Inside one pad: the sharing is real, and it is small",
                 x=0.008, ha="left", fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(FIG, "track_6_inpad.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return out


def main():
    os.makedirs(FIG, exist_ok=True)
    z, j = T.load()
    st = T.station_block(z, j)
    sb = T.selftrack_block(z, j)
    with plt.rc_context(RC):
        for f in (fig_pointing(z, j, st), fig_selftrack(sb),
                  fig_illusion(sb), fig_maps(z, st), fig_purity(z, st),
                  fig_inpad(z, st)):
            print("wrote", f)


if __name__ == "__main__":
    main()
