#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figures_deck.py -- the three-slide MPGD2026 sequence for the VMM efficiency
deficit on P2_OUT.

The talk arrives here having already shown DREAM working (efficiency, timing,
space) and the VMM matching it on timing and space.  These three slides are the
turn: the VMM is *less efficient*, that deficit is a place on the chamber, and
it is one discriminator level sitting inside the Landau.

  deck_1_deficit.png   the observation -- same pads, same beam, same working
                       point, 96 % vs 85 %, and the loss is a corner
  deck_2_gainmap.png   the cause is upstream of both readouts -- a factor 3.9
                       gas-gain gradient that BOTH measure, pad for pad
  deck_3_threshold.png the mechanism -- the Landau slides across a fixed
                       discriminator line; the closer.  Carries the audit on
                       it: the one global level in blue, and each pad's own
                       independently fitted threshold as a green tick, for
                       the audience that asks whether it really is ONE level.
                       The only one of the three with no chrome -- it goes in
                       the talk's own template, which supplies the title

BUILT FOR A PROJECTOR, not for reading at a desk.  Everything is sized off
`SC`, roughly twice what a printed figure would use, and that costs about half
the information density: the headline is one short sentence, the sub-line two,
and every in-panel note is at most two short lines.  The pad maps carry no tick
labels at all -- a scale bar says what 20 mm is, which is the only thing the
absolute coordinates were ever used for.

COLOUR.  Two jobs, two encodings, kept apart:
  * WHICH READOUT is categorical -- DREAM orange, VMM blue -- and is carried by
    the series marks (dots, bars, lines) and the panel names.
  * WHAT IS MAPPED is sequential, and each slide owns a hue so the audience can
    tell at a glance which quantity is on screen: slide 1 EFFICIENCY is purple,
    slide 2 GAIN is amber.  Both maps on a slide share one ramp and one
    colorbar -- that is what makes them comparable, and it is also why nobody
    can read the ramp as a readout identity.
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Ellipse

import figures as F
import figures_slide as S
import threshold_model as M

HERE = os.path.dirname(os.path.abspath(__file__))
FIG, DATA = os.path.join(HERE, "figures"), os.path.join(HERE, "data")

VMM_C, DREAM_C = F.C1, F.C2
VMM_LBL, DREAM_LBL = "VMM3a", "DREAM"

# --- type scale ------------------------------------------------------------ #
# One knob.  SC = 1.0 is the desk-reading size these figures started at; the
# deck runs at ~2x so the back of a conference hall can read it.
SC = 1.70
HEAD, SUB, STAMP = 15.5 * SC, 9.8 * SC, 7.5 * SC
PTITLE, TITLE, NOTE = 13.5 * SC, 12.0 * SC, 10.2 * SC
# Axis numbers and axis labels are the last thing to survive a projector, and
# at 9.2/9.8 they were the smallest type on the slide -- below the in-panel
# notes they sit next to.  Raised to sit just under NOTE.
AXL, TICK, CBL = 11.6 * SC, 11.0 * SC, 11.2 * SC

RC = {"font.size": TICK, "axes.labelsize": AXL, "axes.titlesize": TITLE,
      "xtick.labelsize": TICK, "ytick.labelsize": TICK,
      "xtick.major.size": 3.5 * SC, "ytick.major.size": 3.5 * SC,
      "xtick.major.width": 0.8 * SC, "ytick.major.width": 0.8 * SC,
      "axes.linewidth": 0.8 * SC}

# --- the two sequential ramps ---------------------------------------------- #
# Both pass the ordinal checks against this surface: monotone in OKLab
# lightness, every adjacent gap >= 0.06, one hue, and -- the one that matters on
# a projector -- the LIGHT END clears 2:1 against the slide background, so the
# pale steps do not wash out to paper.  That floor is why neither ramp starts at
# a tint: #b9a3da is 2.20:1 and #e0a94f is 2.05:1, and anything lighter fails.
PURPLE = LinearSegmentedColormap.from_list(   # slide 1: efficiency
    "eff", ["#b9a3da", "#a288cc", "#8a6bbc", "#7250a4", "#583586", "#3a1c60"])
AMBER = LinearSegmentedColormap.from_list(    # slide 2: gain
    "gain", ["#e0a94f", "#d4903a", "#c47829", "#ac5f1a", "#8e470f", "#6c3306"])

SLIDE = (13.33, 7.5)
DPI = 160
XLIM, YLIM = (343, 465), (161, 280)
NGROUP = 6          # fewer, taller bands than the desk figure: 53 pads / 6


def chrome(fig, head, sub, stamp):
    """Headline, sub-line and source stamp -- identical on all three slides."""
    fig.text(0.010, 0.930, head, ha="left", va="bottom", fontsize=HEAD,
             fontweight="bold", color=F.INK)
    fig.text(0.010, 0.822, sub, ha="left", va="bottom", fontsize=SUB,
             color=F.INK2, linespacing=1.45)
    fig.text(0.992, 0.010, stamp, ha="right", va="bottom", fontsize=STAMP,
             color=F.MUTED)


def scalebar(ax, x0=434.0, y0=173.0, mm=20.0):
    ax.plot([x0, x0 + mm], [y0, y0], lw=2.6 * SC, color=F.INK2,
            solid_capstyle="butt", zorder=8)
    ax.text(x0 + mm / 2, y0 + 3.5, f"{mm:.0f} mm", fontsize=NOTE * 0.86,
            color=F.INK2, ha="center", va="bottom", zorder=8)


PAD_MM = 11.88          # measured pad-centre spacing on these 53 pads


def _cell_mask(x, y, GX, GY):
    """True where a grid point falls inside some pad's own square cell.

    A radius cut around each pad centre scallops the outline -- the discs of
    the outermost pads bulge past the footprint and bite into each other -- and
    that scalloping reads as structure the detector does not have.  The pads
    are a regular lattice (11.88 mm, one spacing to 0.02 mm), so give each pad
    the square cell it actually occupies: rotate into the lattice frame, and
    keep grid points within half a pitch in BOTH lattice directions.  The union
    is then the true stair-stepped pad footprint, edge pads included.
    """
    from scipy.spatial import cKDTree
    p = np.column_stack([x, y])
    # lattice orientation from the shortest pad-to-pad vector
    d, j = cKDTree(p).query(p, k=2)
    v = p[j[:, 1]] - p
    v = v[np.argmin(d[:, 1])]
    a = np.arctan2(v[1], v[0])
    R = np.array([[np.cos(-a), -np.sin(-a)], [np.sin(-a), np.cos(-a)]])
    pr = p @ R.T
    gr = np.column_stack([GX.ravel(), GY.ravel()]) @ R.T
    # Chebyshev distance in the rotated frame == "inside a square cell"
    dd, _ = cKDTree(pr).query(gr, p=np.inf)
    return (dd <= 0.5 * PAD_MM).reshape(GX.shape)


def padsurface(ax, x, y, c, norm, cmap, label, colour):
    """The same map as a continuous surface instead of 53 discrete dots.

    What the slide is claiming is a GRADIENT -- one smooth roll-off across the
    chamber -- and a field of separate dots makes the eye read 53 independent
    measurements instead of one shape.  Interpolated linearly between pad
    centres and masked back to within 0.75 of a pad pitch of a real pad, so
    nothing is painted where no pad was read and the footprint keeps its own
    ragged edge rather than being squared off by the grid.

    The cost is honest and worth stating: between two pad centres the surface
    shows a value nobody measured, and a single anomalous pad is spread over
    its neighbourhood instead of standing alone.  On slide 1 that is why the
    one dead pad is called out in the third panel, where it is still a point.
    """
    from scipy.interpolate import griddata
    gx = np.linspace(x.min() - PAD_MM, x.max() + PAD_MM, 340)
    gy = np.linspace(y.min() - PAD_MM, y.max() + PAD_MM, 340)
    GX, GY = np.meshgrid(gx, gy)
    Z = griddata((x, y), c, (GX, GY), method="linear")
    # nearest-neighbour fill first, so the MASK decides the outline; linear
    # alone stops at the triangulation hull and clips the outer half-pads
    Zn = griddata((x, y), c, (GX, GY), method="nearest")
    Z = np.where(np.isfinite(Z), Z, Zn)
    Z = np.where(_cell_mask(x, y, GX, GY), Z, np.nan)
    sc = ax.imshow(Z, origin="lower", cmap=cmap, norm=norm, aspect="equal",
                   extent=[gx[0], gx[-1], gy[0], gy[-1]],
                   interpolation="bilinear", zorder=3)
    _map_frame(ax, label, colour)
    return sc


def padmap(ax, x, y, c, norm, cmap, label, colour):
    """The discrete version: one mark per pad, kept for the desk figures."""
    sc = ax.scatter(x, y, c=c, s=420, cmap=cmap, norm=norm,
                    edgecolor=F.SURFACE, linewidth=1.4, zorder=3)
    _map_frame(ax, label, colour)
    return sc


def _map_frame(ax, label, colour):
    """Frame, limits, title and scale bar -- shared by both renderings.  No
    tick labels: the absolute mm never mattered, only the pattern and the
    scale, and the scale is a bar."""
    ax.set_aspect("equal")
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_color(F.GRID)
    ax.set_title(label, loc="left", color=colour, fontsize=PTITLE,
                 fontweight="bold", pad=8 * SC)
    scalebar(ax)


def hbar(fig, sc, label, ticks, rect=(0.050, 0.093, 0.38, 0.028)):
    """One horizontal colorbar under the two maps -- they share a scale, so
    they share a legend."""
    cax = fig.add_axes(rect)
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.outline.set_edgecolor(F.MUTED)
    cb.ax.tick_params(labelsize=TICK, color=F.MUTED)
    cb.set_ticks(ticks)
    cb.set_label(label, fontsize=CBL, color=F.INK2, labelpad=4)
    return cb


CORNER = dict(xy=(370.0, 214.0), width=44.0, height=66.0)


def corner(ax):
    """The low-gain / low-efficiency corner, marked identically on slide 1 and
    slide 2 so the eye carries it from one to the other."""
    ax.add_patch(Ellipse(**CORNER, facecolor="none", edgecolor=F.INK,
                         lw=1.4 * SC, ls=(0, (5, 4)), zorder=7))


def gradient(g):
    """The plane fit both readouts share: direction, per-readout slope along it,
    and how much of the pad-to-pad variance it takes."""
    x, y = g["x_v"].to_numpy(), g["y_v"].to_numpy()
    A = np.column_stack([np.ones(len(g)), x, y])
    out = {}
    for key, col in (("d", "amp_med_d"), ("v", "amp_med_v")):
        r = (g[col] / g[col].median()).to_numpy()
        b = np.linalg.lstsq(A, r, rcond=None)[0]
        out[key] = dict(b=b, ang=np.arctan2(b[2], b[1]),
                        r2=1 - np.var(r - A @ b) / np.var(r), rel=r)
    ang = out["d"]["ang"]
    s = (x - x.mean()) * np.cos(ang) + (y - y.mean()) * np.sin(ang)
    for key in ("d", "v"):
        out[key]["slope10"] = np.polyfit(s, out[key]["rel"], 1)[0] * 10
    out["s"], out["ang"] = s, ang
    out["dang"] = np.degrees(abs(np.arctan2(
        np.sin(out["v"]["ang"] - ang), np.cos(out["v"]["ang"] - ang))))
    return out


def _maps_grid(fig):
    # The maps are surfaces now, not fields of dots, so they carry their
    # pattern at a smaller size; the width that frees goes to the third panel,
    # which needs it for a y axis at the raised type size.
    return fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.45],
                            left=0.015, right=0.992, top=0.735, bottom=0.205,
                            wspace=0.13)


def _shift_right(ax, dx=0.036):
    """Push the third panel clear of the middle map.

    The two maps carry no tick labels and sit happily close; the third panel
    carries a y axis -- numbers AND a label -- which at this type size reaches
    ~5 % of the figure width to its left and lands on the map.  `wspace` is one
    number for every gap in the row, so widening it would also pull the two
    maps apart and run their captions into each other.  Moving this one axes
    after the fact is the only change that touches nothing else.
    """
    b = ax.get_position()
    ax.set_position([b.x0 + dx, b.y0, b.width - dx, b.height])


def _panel(ax, title):
    ax.grid(True, lw=0.6 * SC, alpha=0.7)
    ax.set_axisbelow(True)
    ax.set_title(title, loc="left", color=F.INK, fontsize=TITLE,
                 fontweight="bold", pad=8 * SC)


# --------------------------------------------------------------------------- #
def slide_1(g, n):
    """The observation.  Two efficiency maps on one purple ramp, then the
    pad-by-pad pairing that says the deficit is not a uniform scale factor."""
    x, y = g["x_v"].to_numpy(), g["y_v"].to_numpy()
    ev, ed = g["eff_v"].to_numpy(), g["eff_d"].to_numpy()

    fig = plt.figure(figsize=SLIDE)
    gs = _maps_grid(fig)
    norm = Normalize(0.45, 1.0)

    axd = fig.add_subplot(gs[0, 0])
    sc = padsurface(axd, x, y, ed, norm, PURPLE,
                    f"{DREAM_LBL} · {n['eff_dream_all'] * 100:.1f} %", DREAM_C)
    axv = fig.add_subplot(gs[0, 1])
    padsurface(axv, x, y, ev, norm, PURPLE,
               f"{VMM_LBL} · {n['eff_obs_all'] * 100:.1f} %", VMM_C)
    hbar(fig, sc, "efficiency, per pad", [0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    corner(axv)
    # the callout goes in the empty xlabel slot, not over the pads: at this
    # type size there is no clear space left inside a map.
    axd.set_xlabel("uniform across the chamber", fontsize=NOTE,
                   color=F.INK, fontweight="bold", labelpad=6 * SC)
    axv.set_xlabel("all its losses are in the ring", fontsize=NOTE,
                   color=F.INK, fontweight="bold", labelpad=6 * SC)

    # -- pad-by-pad pairing -------------------------------------------------- #
    ax = fig.add_subplot(gs[0, 2])
    _shift_right(ax)
    order = np.argsort(ev)
    r = np.arange(len(order))
    ax.vlines(r, ev[order], ed[order], color=F.MUTED, lw=1.1 * SC, alpha=0.65,
              zorder=2)
    ax.scatter(r, ed[order], s=70, color=DREAM_C, zorder=4,
               edgecolor=F.SURFACE, linewidth=1.1)
    ax.scatter(r, ev[order], s=70, color=VMM_C, zorder=4,
               edgecolor=F.SURFACE, linewidth=1.1)
    ax.grid(True, axis="y", lw=0.6 * SC, alpha=0.7)
    ax.set_axisbelow(True)
    ax.set_xticks([])
    ax.set_xlim(-3.0, len(order) + 2.0)
    ax.set_ylim(0.24, 1.13)
    ax.set_yticks([0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("the 53 pads under the beam")
    ax.set_ylabel("efficiency")
    _panel(ax, "A place, not a scale factor")
    ax.text(len(order) + 1, 1.005, DREAM_LBL, fontsize=PTITLE * 0.82,
            color=DREAM_C, fontweight="bold", ha="right", va="bottom")
    ax.text(34, 0.855, VMM_LBL, fontsize=PTITLE * 0.82, color=VMM_C,
            fontweight="bold", ha="center", va="top")
    ax.annotate("half the tracks lost\non the weakest pad",
                xy=(1.4, ev[order][1]), xytext=(8, 0.395), fontsize=NOTE,
                color=F.INK2, va="center", ha="left",
                arrowprops=dict(arrowstyle="-", color=F.INK2, lw=0.9 * SC))
    ax.scatter([0], [ev[order][0]], s=260, facecolor="none",
               edgecolor=F.INK, linewidth=1.2 * SC, zorder=5)
    ax.text(2.8, ev[order][0], "dead in both", fontsize=NOTE * 0.86,
            color=F.INK2, va="center")

    chrome(fig,
           "DREAM 95.6 %, VMM 85.3 % — same pads, same beam",
           "P2_OUT.  Efficiency of the uRWELL tracks pointing at each pad, "
           "measured the same way on both.\nMatched runs at mesh 450 V / "
           "drift 750 V: VMM run_46, DREAM eff_nominal_1.",
           "P2 SPS July 2026 · 53 pads · 0.84 M (VMM) / 1.87 M (DREAM) tracks")
    fig.savefig(f"{FIG}/deck_1_deficit.png", dpi=DPI)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def slide_2(g, n):
    """The cause is upstream of the electronics: one gas-gain gradient, in
    amber, and the two readouts return the same map pad for pad."""
    x, y = g["x_v"].to_numpy(), g["y_v"].to_numpy()
    G = gradient(g)
    rd, rv = G["d"]["rel"], G["v"]["rel"]

    fig = plt.figure(figsize=SLIDE)
    gs = _maps_grid(fig)
    norm = Normalize(0.5, 2.0)

    axd = fig.add_subplot(gs[0, 0])
    sc = padsurface(axd, x, y, rd, norm, AMBER, DREAM_LBL, DREAM_C)
    axv = fig.add_subplot(gs[0, 1])
    padsurface(axv, x, y, rv, norm, AMBER, VMM_LBL, VMM_C)
    hbar(fig, sc, "pad gain  /  that readout's own median",
         [0.5, 1.0, 1.5, 2.0])

    ang = G["ang"]
    for ax, key in ((axd, "d"), (axv, "v")):
        corner(ax)
        ax.annotate("", xy=(x.mean() + 44 * np.cos(ang),
                            y.mean() + 44 * np.sin(ang)),
                    xytext=(x.mean() - 44 * np.cos(ang),
                            y.mean() - 44 * np.sin(ang)),
                    arrowprops=dict(arrowstyle="-|>", color=F.INK,
                                    lw=1.7 * SC), zorder=6)
        # in the xlabel slot, as on slide 1 -- inside the map it lands on the
        # bottom-left pads
        ax.set_xlabel(f"{G[key]['slope10'] * 100:+.0f} % per 10 mm",
                      fontsize=NOTE, color=F.INK, fontweight="bold",
                      labelpad=6 * SC)
    axv.text(346, 278, "the same corner", fontsize=NOTE, color=F.INK,
             fontweight="bold", va="top", ha="left")

    # -- pad for pad --------------------------------------------------------- #
    ax = fig.add_subplot(gs[0, 2])
    _shift_right(ax)
    lim = (0.40, 2.18)
    ax.plot(lim, lim, ls=(0, (5, 4)), lw=1.2 * SC, color=F.MUTED, zorder=2)
    ax.scatter(rd, rv, s=np.clip(g["n_track_v"] / 210, 40, 280),
               facecolor=VMM_C, alpha=0.52, edgecolor=F.SURFACE,
               linewidth=1.0, zorder=4)
    b = np.polyfit(rd, rv, 1)
    xs = np.linspace(*lim, 20)
    ax.plot(xs, b[0] * xs + b[1], lw=2.2 * SC, color=VMM_C, zorder=5)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal")
    ax.set_xticks([0.5, 1.0, 1.5, 2.0])
    ax.set_yticks([0.5, 1.0, 1.5, 2.0])
    ax.set_xlabel("relative pad gain, DREAM")
    ax.set_ylabel("relative pad gain, VMM")
    _panel(ax, f"Pad for pad:  r = {np.corrcoef(rd, rv)[0, 1]:+.2f}")
    ax.text(0.47, 2.10, "one map,\nmeasured twice", fontsize=NOTE,
            color=F.INK, va="top")
    ax.text(1.02, 0.94, f"the VMM's copy is\nFLATTER: slope {b[0]:.2f}",
            fontsize=NOTE, color=VMM_C, fontweight="bold", va="top")
    ax.text(2.13, 2.04, "1:1", fontsize=NOTE * 0.86, color=F.MUTED,
            ha="right", va="top")

    chrome(fig,
           "The chamber's gain rolls off ×3.9 — and both readouts see it",
           f"The same corner, now in pulse height.  One plane in (x, y) "
           f"takes {G['d']['r2'] * 100:.0f} % (DREAM) / "
           f"{G['v']['r2'] * 100:.0f} % (VMM) of the variance,\nin the same "
           f"direction to {G['dang']:.0f}°.  The variation is the chamber's, "
           f"not the electronics'.",
           "P2 SPS July 2026 · P2_OUT · leading-pad pulse height on tracked "
           "events")
    fig.savefig(f"{FIG}/deck_2_gainmap.png", dpi=DPI)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def slide_3(g, H, bw, Sp, n):
    """The closer: the spectra themselves, banded by gain, against the
    discriminator line -- and what each band actually records.

    NO CHROME.  This is the one slide that goes into the talk's own template,
    so the headline, the sub-line and the source stamp are the deck's job and
    not the figure's; what is left on the picture is only what it cannot be
    read without.  Everything that would have been said in prose is said as a
    label instead: BLUE is the single level fitted to all 53 pads at once,
    GREEN is that same model inverted for each pad ALONE, one tick per pad.
    Only the blue one is keyed on the figure -- the green key goes on in
    PowerPoint, so nothing here may occupy the strip under the rows.

    With the chrome gone the panels take the whole slide, and the space that
    freed up at the top is what the blue key sits in -- inside the axes, not
    in a legend box, because it has to sit against the line it names."""
    gr = S.groups(g, H, bw, ngroup=NGROUP)
    T = n["T"]
    Tpad, _, lo_clip, hi_clip = M.fit_threshold_per_pad(g, Sp)
    rows = [Tpad[G["idx"]] for G in gr]

    fig = plt.figure(figsize=SLIDE)
    gs = fig.add_gridspec(1, 2, width_ratios=[2.45, 1],
                          left=0.068, right=0.962, top=0.945, bottom=0.145,
                          wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1], sharey=ax)
    S.ridge(ax, gr, bw, T, fs=SC, Tpad_rows=rows,
            tick=(0.03, 0.58, 1.7 * SC, 0.95))
    S.effbars(axb, gr, n, fs=SC, bh=0.38)

    # The lost/recorded and VMM3a/DREAM keys live BELOW the lowest band, in the
    # margin the y limit leaves under it.  The desk figure's -0.42 leaves ~0.08
    # of clear space, which at this type size put them on the axis line and
    # into the tick numbers underneath.  Deepen that margin and hang the keys
    # in the middle of it: no label moves relative to the data, they simply
    # stop sharing a row with the axis.  Shared y, so one call does both.
    ax.set_ylim(-0.92, NGROUP + 0.35)
    KEY = -0.46           # the key row, clear of both spine and lowest band

    ax.annotate("the VMM's discriminator", xy=(T, NGROUP + 0.01),
                xytext=(T * 1.32, NGROUP + 0.26), fontsize=NOTE,
                color=VMM_C, fontweight="bold", va="center",
                arrowprops=dict(arrowstyle="-", color=VMM_C, lw=1.0 * SC))
    ax.text(T * 0.46, KEY, "lost", fontsize=NOTE * 1.15, color=DREAM_C,
            fontweight="bold", ha="center", va="center")
    ax.text(T * 2.4, KEY, "recorded", fontsize=NOTE * 1.15, color=F.INK2,
            ha="center", va="center")
    ax.annotate("", xy=(T * 0.98, KEY), xytext=(T * 0.66, KEY),
                arrowprops=dict(arrowstyle="<-", color=DREAM_C, lw=1.1 * SC))
    ax.text(30.5, NGROUP + 0.08, "pad gain", fontsize=NOTE, color=F.INK2,
            ha="left", va="bottom", style="italic")
    ax.set_xlabel("pulse height on the pad  [DREAM ADC]", labelpad=7 * SC)
    ax.set_title("log axis — so a gain factor is a sideways shift",
                 loc="left", color=F.INK2, fontsize=TITLE * 0.86, pad=8 * SC)

    axb.text(0.405, NGROUP + 0.08, "% of those tracks recorded",
             fontsize=NOTE, color=F.INK2, ha="left", va="bottom")
    axb.text(0.50, KEY, VMM_LBL, fontsize=PTITLE * 0.82, color=VMM_C,
             fontweight="bold", ha="center", va="center")
    axb.text(0.74, KEY, DREAM_LBL, fontsize=PTITLE * 0.82, color=DREAM_C,
             fontweight="bold", ha="center", va="center")

    fig.savefig(f"{FIG}/deck_3_threshold.png", dpi=DPI)
    plt.close(fig)


def main():
    g, H, bw, Sp, n = S.load()
    with plt.rc_context(RC):
        slide_1(g, n)
        slide_2(g, n)
        slide_3(g, H, bw, Sp, n)
    print("wrote figures/deck_{1_deficit,2_gainmap,3_threshold}.png")


if __name__ == "__main__":
    main()
