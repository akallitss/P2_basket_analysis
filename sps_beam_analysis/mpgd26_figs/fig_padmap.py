#!/usr/bin/env python3
"""Per-pad efficiency maps at the working point + cold-connector highlight."""
import os, sys, glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2style as st
import matplotlib.pyplot as plt

from paths import S, A, URW, RD, OUT  # noqa: E402
DETS = ['P2_IN', 'P2_MID', 'P2_OUT']

fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))
last = None
for ax, det in zip(axes, DETS):
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
    sc = ax.scatter(df.pad_cx, df.pad_cy, c=df.eff, cmap=st.SEQ_CMAP,
                    vmin=0.5, vmax=1.0, s=52, marker='s', edgecolors='none')
    last = sc
    # ring the cold pads (eff < 0.5 * healthy median)
    med = df.eff.median()
    cold = df[df.eff < 0.75 * med]
    ax.scatter(cold.pad_cx, cold.pad_cy, facecolors='none',
               edgecolors=st.C_RED, s=90, marker='s', linewidths=1.4,
               label=f'{len(cold)} cold pads')
    ax.set_title(f'{det}   (median eff {med:.3f})', color=st.DET_COLOR[det])
    ax.set_xlabel('pad x [mm]')
    ax.set_aspect('equal')
    ax.grid(False)
    if len(cold):
        ax.legend(loc='upper right', fontsize=8.5)
axes[0].set_ylabel('pad y [mm]')
fig.colorbar(last, ax=axes, label='per-pad tag-probe efficiency',
             fraction=0.02)
fig.suptitle('Per-pad efficiency at the working point (mesh 450 / drift 700, '
             'beam-illuminated pads, n_tag $\\geq$ 2000)', fontweight='bold')
st.finish(fig, f'{OUT}/padmap_working_point.png', tight=False)
print('wrote', f'{OUT}/padmap_working_point.png')
