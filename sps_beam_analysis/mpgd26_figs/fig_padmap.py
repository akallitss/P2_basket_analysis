#!/usr/bin/env python3
"""Per-pad efficiency maps at the working point + cold-connector highlight.

Two rows, because they answer two different questions:

  top     the pads themselves, on the fixed 0.5-1.0 talk scale, with the cold
          pads ringed -- "are any pads dead?"
  bottom  the same numbers as a smooth surface on a stretched scale around the
          bulk -- "is what is left uniform?"  At 0.5-1.0 every healthy pad is
          the same dark blue and a 3-point gradient across the chamber is
          invisible; the auto-stretched scale is the only way that structure
          reaches the back of the room.  The colourbars are separate and
          labelled, so the two rows are never read as one scale.
"""
import os, sys, glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2style as st
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy.interpolate import griddata

from paths import S, A, URW, RD, OUT  # noqa: E402
DETS = ['P2_IN', 'P2_MID', 'P2_OUT']
PAD_MM = 12.0

# the five big mesh-support pillars, from the insulation-mask Gerber -- same
# board frame as pad_cx/pad_cy, so they drop straight onto these axes.  One of
# them sits at (411, 238), inside the beam spot on all three stations: it is
# the white hole that shows up in the middle of every surface panel.
_CBA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), 'cosmic_bench_analysis')
sys.path.insert(0, _CBA)
try:
    import p2_qa_config as _qa
    import p2_mapping as _pm
    _pil = _pm.load_pillars(_qa.MASK_GBR_PATH)
    BIG = _pil[_pil.big] if len(_pil) else None
except Exception as _e:                      # figure still stands without them
    print('  ! pillars unavailable:', _e)
    BIG = None


def draw_big_pillars(ax, label=None):
    from matplotlib.patches import Circle
    if BIG is None:
        return
    first = True
    for xi, yi, ri in zip(BIG.x, BIG.y, BIG.r):
        ax.add_patch(Circle((xi, yi), ri, facecolor=st.C_RED,
                            edgecolor=st.C_RED, zorder=6))
        ax.add_patch(Circle((xi, yi), 9.0, facecolor='none', edgecolor=st.C_RED,
                            ls='--', lw=1.5, zorder=6,
                            label=(label if first else None)))
        first = False

fig, axes = plt.subplots(2, 3, figsize=(16.2, 10.4))
last, last_s = None, None
frames = {}

for col, det in enumerate(DETS):
    # working-point map: highstat_eff_1 aggregate (scan dir holds per-sub_run
    # maps; use beam_commissioning_00)
    cands = sorted(glob.glob(
        f'{A}/{det}/highstat_eff_1/*/22_tag_probe_efficiency/'
        f'eff_map_{det}_beam_commissioning_00*.csv'))
    if not cands:
        cands = sorted(glob.glob(
            f'{A}/{det}/eff_nominal_1/scan/22_tag_probe_efficiency/'
            f'eff_map_{det}_*.csv'))
    df = pd.read_csv(cands[0])
    df = df[df.n_tag >= 2000].copy()
    df['eff'] = pd.to_numeric(df['eff'], errors='coerce')
    df = df.dropna(subset=['eff'])
    frames[det] = df

# one stretched scale for the whole bottom row, so the three stations stay
# comparable to each other even though they are not comparable to the top row
allv = np.concatenate([f.eff.to_numpy() for f in frames.values()])
lo = max(0.80, float(np.percentile(allv, 2)))
hi = min(1.0, float(np.percentile(allv, 99.5)))

for col, det in enumerate(DETS):
    df = frames[det]
    ax, axs = axes[0, col], axes[1, col]

    # ---- top: the pads, fixed talk scale ---------------------------------- #
    sc = ax.scatter(df.pad_cx, df.pad_cy, c=df.eff, cmap=st.SEQ_CMAP,
                    vmin=0.5, vmax=1.0, s=62, marker='s', edgecolors='none')
    last = sc
    med = df.eff.median()
    cold = df[df.eff < 0.75 * med]
    ax.scatter(cold.pad_cx, cold.pad_cy, facecolors='none',
               edgecolors=st.C_RED, s=110, marker='s', linewidths=1.6,
               label=f'{len(cold)} cold pad' + ('s' if len(cold) != 1 else ''))
    draw_big_pillars(ax)
    ax.set_title(f'{det}  —  median eff {med:.3f}', color=st.DET_COLOR[det],
                 fontsize=15)
    if len(cold):
        ax.legend(loc='lower left', fontsize=11)

    # ---- bottom: the same numbers as a surface ---------------------------- #
    # linear interpolation between pad centres, masked back to the illuminated
    # footprint so nothing is drawn where no pad was read.
    x, y, v = df.pad_cx.to_numpy(), df.pad_cy.to_numpy(), df.eff.to_numpy()
    gx = np.linspace(x.min() - PAD_MM, x.max() + PAD_MM, 260)
    gy = np.linspace(y.min() - PAD_MM, y.max() + PAD_MM, 260)
    GX, GY = np.meshgrid(gx, gy)
    Z = griddata((x, y), v, (GX, GY), method='linear')
    # keep only grid points within one pad of a real pad centre
    from scipy.spatial import cKDTree
    dist, _ = cKDTree(np.column_stack([x, y])).query(
        np.column_stack([GX.ravel(), GY.ravel()]))
    Z = np.where(dist.reshape(Z.shape) <= 0.85 * PAD_MM, Z, np.nan)
    last_s = axs.imshow(Z, origin='lower', cmap=st.SEQ_CMAP, vmin=lo, vmax=hi,
                        extent=[gx[0], gx[-1], gy[0], gy[-1]], aspect='equal',
                        interpolation='bilinear')
    axs.scatter(cold.pad_cx, cold.pad_cy, facecolors='none',
                edgecolors=st.C_RED, s=110, marker='s', linewidths=1.6)
    draw_big_pillars(axs, label='big mesh pillar' if col == 0 else None)
    p16, p84 = np.percentile(v, [16, 84])
    axs.set_title(f'spread p16–p84 = {100 * (p84 - p16):.1f} points',
                  color=st.DET_COLOR[det], fontsize=15)
    if col == 0:
        axs.legend(loc='lower left', fontsize=11)

    for a in (ax, axs):
        a.set_xlabel('pad x [mm]')
        a.set_aspect('equal')
        a.grid(False)
        # one frame for all six panels: the stations only compare if the axes do
        a.set_xlim(292, 528)
        a.set_ylim(132, 344)

axes[0, 0].set_ylabel('pad y [mm]')
axes[1, 0].set_ylabel('pad y [mm]')

cb1 = fig.colorbar(last, ax=axes[0, :], fraction=0.022, pad=0.015)
st.bold_cbar(cb1, 'per-pad tag-probe efficiency\n(fixed 0.5–1.0 talk scale)',
             size=13)
cb2 = fig.colorbar(last_s, ax=axes[1, :], fraction=0.022, pad=0.015)
st.bold_cbar(cb2, f'same efficiency, surface\n(stretched to {lo:.2f}–{hi:.2f})',
             size=13)

fig.suptitle('Per-pad efficiency at the working point (mesh 450 / drift 700, '
             'beam-illuminated pads, n_tag $\\geq$ 2000)\n'
             'top: every pad on the talk scale — bottom: the same numbers as a '
             'surface, stretched, so the uniformity is visible\n'
             'the white hole at (411, 238) in all three surfaces is a big mesh '
             'pillar, from the Gerber — not a fit, not a dead channel',
             fontweight='bold', fontsize=16)
st.finish(fig, f'{OUT}/padmap_working_point.png', tight=False)
print('wrote', f'{OUT}/padmap_working_point.png')
