#!/usr/bin/env python3
"""Efficiency figures: HV scan curves, method comparison, 2D maps, P2_IN swap."""
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2style as st
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import chamber_history as ch

from paths import S, A, URW, RD, OUT  # noqa: E402

DETS = ['P2_IN', 'P2_MID', 'P2_OUT']
DY = {'P2_IN': 0, 'P2_MID': 9, 'P2_OUT': -9}
u = pd.read_csv(f'{URW}/drift_mesh_scan_1/urw_p2_efficiency_drift_mesh_scan_1.csv')

# ---------------------------------------------------------------- F1: eps(mesh)
# The mesh turn-on is spread over three runs, deliberately: drift_mesh_scan_1
# (25 Jul) covers 390-450 V, low_mesh_scan_1 (26 Jul) walks MID/OUT down to
# 330 V, and the p2in_hvrange_* pair (28 Jul) brings the P2_IN station up from
# 200 V.  Joined per station ONLY where it is the same physical chamber: P2_MID
# and P2_OUT never changed (det1, det3), but P2_IN held det4 on 25 Jul and the
# CERN-built chamber from 28 Jul, so those are two separate series.
# See ../chamber_history.py.
def _urw(run):
    p = f'{URW}/{run}/urw_p2_efficiency_{run}.csv'
    return pd.read_csv(p) if os.path.exists(p) else pd.DataFrame()


m = u[u.sub_run.str.startswith('meshscan')].copy()
low = _urw('low_mesh_scan_1')
inA, inB = _urw('p2in_hvrange_1'), _urw('p2in_hvrange_2')


def _series(det, frames):
    """Concatenate scan points of one station, keeping the drift gap."""
    parts = [f[f.station == det] for f in frames if len(f)]
    if not parts:
        return pd.DataFrame()
    g = pd.concat(parts, ignore_index=True)
    # the p2in runs park the other two stations at their working point; those
    # single points are not scan points and would draw a spurious spur
    g = g[g.groupby('mesh_hv').mesh_hv.transform('size') > 0]
    return g.sort_values('mesh_hv')


fig, ax = plt.subplots(figsize=(9.4, 5.8))
for det in DETS:
    frames = [m, low] if det != 'P2_IN' else [m]
    g = _series(det, frames)
    if not len(g):
        continue
    gap = (g.drift_hv - g.mesh_hv)
    gaps = sorted(gap.unique())
    gl = (f'{gaps[0]:.0f} V' if len(gaps) == 1
          else f'{min(gaps):.0f}–{max(gaps):.0f} V')
    lab = ch.label(det, 'drift_mesh_scan_1', window=(det == 'P2_IN'))
    ax.errorbar(g.mesh_hv, g.eff, yerr=[g.eff - g.lo, g.hi - g.eff],
                color=st.DET_COLOR[det], marker='o', capsize=2.5, lw=2,
                label=f'{lab} — drift = mesh + {gl}')
    st.direct_label(ax, g.mesh_hv.iloc[-1], g.eff.iloc[-1], det,
                    st.DET_COLOR[det], dy=DY[det])

# the replacement P2_IN: a different chamber, so a different series
rep = _series('P2_IN', [inA, inB])
if len(rep):
    ax.errorbar(rep.mesh_hv, rep.eff, yerr=[rep.eff - rep.lo, rep.hi - rep.eff],
                color=st.C_VIO, marker='^', ls='--', capsize=2.5, lw=2,
                label=ch.label('P2_IN', 'p2in_hvrange_1', window=True)
                      + ' — drift = mesh + 300 V')
    st.direct_label(ax, rep.mesh_hv.iloc[-1], rep.eff.iloc[-1],
                    'P2_IN (det5)', st.C_VIO, dx=4, dy=-20)

ax.set_xlabel('mesh voltage [V]')
ax.set_ylabel('absolute efficiency')
ax.set_ylim(0, 1.02)
ax.axhline(0.95, color=st.GRID, lw=1, ls=':')
ax.legend(loc='upper left', fontsize=8.5)
ax.set_title('P2 efficiency vs mesh HV — SPS muons, Ar/CO$_2$/iso 93/5/2\n'
             'three runs stitched into one turn-on: drift_mesh_scan_1 + '
             'low_mesh_scan_1 for MID/OUT (det1, det3),\nand p2in_hvrange_1/2 '
             'for det5 = the CERN-built chamber that took over P2_IN on 28 Jul',
             fontsize=10.5)
st.finish(fig, f'{OUT}/eff_vs_mesh_urw.png')

# --------------------------------------------------------------- F2: eps(drift)
d = u[u.sub_run.str.startswith('drift_')].copy()
# P2_IN is deliberately absent: during the drift half of this scan it sat
# parked at a fixed 430/630 V, so its points are not scan points at all --
# plotting them would draw a flat line that looks like a measurement.
d = d[~((d.station == 'P2_IN'))]
fig, ax = plt.subplots(figsize=(8.8, 5.6))
for det in DETS:
    g = d[d.station == det].sort_values('drift_hv')
    if not len(g):
        continue
    ax.errorbar(g.drift_hv, g.eff, yerr=[g.eff - g.lo, g.hi - g.eff],
                color=st.DET_COLOR[det], marker='o', capsize=2.5, lw=2,
                label=ch.label(det, 'drift_mesh_scan_1'))
    st.direct_label(ax, g.drift_hv.iloc[-1], g.eff.iloc[-1], det,
                    st.DET_COLOR[det], dy=DY[det])
ax.plot([], [], ' ',
        label=ch.label('P2_IN', 'drift_mesh_scan_1') + ' — parked at '
              '430/630 V, not scanned')
ax.set_xlabel('drift voltage [V] (mesh fixed at station working point)')
ax.set_ylabel('absolute efficiency')
ax.set_ylim(0, 1.02)
ax.legend(loc='lower right')
ax.set_title('P2 efficiency vs drift HV — drift_450 point = zero drift field\n'
             'only P2_MID (det1) and P2_OUT (det3) were scanned in drift; '
             'P2_IN (det4) was held at its working point', fontsize=11)
st.finish(fig, f'{OUT}/eff_vs_drift_urw.png')

# ------------------------------------------- F3: method comparison on one scan
tp = pd.read_csv(f'{RD}/dream_tag_probe.csv')
rs = pd.read_csv(f'{RD}/dream_raw_stream.csv')
fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.0), sharey=True)
for ax, det in zip(axes, DETS):
    g = m[m.station == det].sort_values('mesh_hv')
    ax.plot(g.mesh_hv, g.eff, color=st.DET_COLOR[det], marker='o', lw=2,
            label='uRWELL tracks (absolute)')
    t = tp[(tp.run == 'drift_mesh_scan_1') & (tp.probe == det)
           & (tp.prod_sub == 'scan')
           & (tp.sub_run.str.startswith('meshscan'))].sort_values('hv')
    if len(t):
        ax.plot(t.hv, t.eff_corr, color=st.TEXT2, marker='s', lw=1.6, ls='--',
                label='tag-and-probe (2-of-3)')
    # The raw zero-suppressed stream at a fixed 100 ADC threshold used to be a
    # third series here.  It was dropped (2026-08-19): it is not an efficiency
    # measurement of the chamber but of that threshold -- it sits 20-50 points
    # low over the whole scan and moves with the gain, so it only invited the
    # question "which one is right?" when the answer is that the first two are.
    ax.set_title(ch.label(det, 'drift_mesh_scan_1'),
                 color=st.DET_COLOR[det])
    ax.set_xlabel('mesh voltage [V]')
    ax.set_ylim(0, 1.02)
axes[0].set_ylabel('efficiency')
axes[0].legend(loc='upper left', fontsize=9)
fig.suptitle('Two independent efficiency methods on the same mesh scan '
             '(drift_mesh_scan_1, 25 Jul)\n'
             'uRWELL-referenced needs an external tracker; tag-and-probe needs '
             'only the other two P2 planes — they agree to 1–2 points',
             fontweight='bold')
st.finish(fig, f'{OUT}/eff_methods_comparison.png')

# ----------------------------------------------- F4: 2D eff(mesh,drift) maps
hv = pd.read_csv(f'{S}/eos_inventory/hv_setpoints.csv').rename(
    columns={'mesh_or_resist': 'mesh_v', 'drift': 'drift_v'})
for run2d in ['drift_mesh_2d_2']:
    t = tp[(tp.run == run2d) & (tp.prod_sub == 'scan')]
    if not len(t):
        continue
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
    for ax, det in zip(axes, DETS):
        g = t[t.probe == det].copy()
        j = g.merge(hv[hv.det == det][['run', 'sub_run', 'mesh_v', 'drift_v']],
                    on=['run', 'sub_run'], how='left')
        j = j.dropna(subset=['mesh_v', 'drift_v'])
        j['mesh_v'] = j.mesh_v.astype(float)
        j['drift_v'] = j.drift_v.astype(float)
        j['gap'] = j.drift_v - j.mesh_v
        piv = j.pivot_table(index='gap', columns='mesh_v', values='eff_corr',
                            aggfunc='mean')
        im = ax.pcolormesh(piv.columns, piv.index, piv.values,
                           cmap=st.SEQ_CMAP, vmin=0, vmax=1,
                           edgecolors=st.SURFACE, linewidth=1.5)
        for yy in piv.index:
            for xx in piv.columns:
                v = piv.loc[yy, xx]
                if np.isfinite(v):
                    ax.text(xx, yy, f'{v:.2f}', ha='center', va='center',
                            fontsize=7.2,
                            color='white' if v > 0.55 else st.TEXT)
        ax.set_title(ch.label(det, run2d), color=st.DET_COLOR[det])
        ax.set_xlabel('mesh voltage [V]')
        ax.grid(False)
    axes[0].set_ylabel('drift $-$ mesh voltage [V]')
    fig.colorbar(im, ax=axes, label='tag-probe efficiency', fraction=0.02)
    fig.suptitle(f'Efficiency over the (mesh, drift) plane — {run2d} '
                 f'(7$\\times$7 grid)', fontweight='bold')
    st.finish(fig, f'{OUT}/eff_2d_{run2d}.png', tight=False)
    print('wrote', f'{OUT}/eff_2d_{run2d}.png')

# ------------------------------- F4b: the same surface as curves, not boxes
# A 7x7 heatmap answers "what is the number here"; the talk needs "what does
# the knob do".  One curve per drift gap, so the gain turn-on and the gap
# ordering are both readable at a glance.
for run2d in ['drift_mesh_2d_2']:
    t = tp[(tp.run == run2d) & (tp.prod_sub == 'scan')]
    if not len(t):
        continue
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), sharey=True)
    gaps = None
    for ax, det in zip(axes, DETS):
        g = t[t.probe == det].copy()
        j = g.merge(hv[hv.det == det][['run', 'sub_run', 'mesh_v', 'drift_v']],
                    on=['run', 'sub_run'], how='left')
        j = j.dropna(subset=['mesh_v', 'drift_v'])
        j['mesh_v'] = j.mesh_v.astype(float)
        j['gap'] = j.drift_v.astype(float) - j.mesh_v
        gaps = sorted(j.gap.unique())
        cmap = plt.get_cmap('viridis')
        for i, gap in enumerate(gaps):
            s = j[j.gap == gap].sort_values('mesh_v')
            col = cmap(i / max(len(gaps) - 1, 1))
            ax.plot(s.mesh_v, s.eff_corr, marker='o', ms=4.5, lw=1.8,
                    color=col, label=f'{gap:.0f} V')
            if len(s):
                st.direct_label(ax, s.mesh_v.iloc[-1], s.eff_corr.iloc[-1],
                                f'{gap:.0f}', col, dx=4, fontsize=8)
        ax.set_title(ch.label(det, run2d), color=st.DET_COLOR[det])
        ax.set_xlabel('mesh voltage [V]')
        ax.set_ylim(0, 1.02)
    axes[0].set_ylabel('tag-probe efficiency')
    axes[-1].legend(title='drift $-$ mesh [V]', fontsize=8, title_fontsize=8,
                    loc='lower right', ncol=2)
    fig.suptitle('Efficiency over the (mesh, drift) plane, as curves — '
                 f'{run2d}\ngain sets the turn-on; the drift gap only moves '
                 'the curve sideways once it is above ~250 V',
                 fontweight='bold')
    st.finish(fig, f'{OUT}/eff_2d_curves_{run2d}.png')

# ------------------------------------ F5: every chamber that sat at P2_IN
# Rewritten 2026-08-19.  The old version drew TWO CERN series, one of them
# `p2_mesh_drift_eff_1` -- a 2D (mesh, drift) scan plotted against mesh alone,
# so every mesh voltage appeared once per drift gap and the curve grew a
# vertical stack of duplicate points at 440 V.  It also named det2 in the
# title without ever plotting it.
#
# Now: one clean series per chamber, each at a stated drift gap, so the three
# chambers that occupied this one station are actually comparable.
hvset = pd.read_csv(f'{S}/eos_inventory/hv_setpoints.csv')


def _gap_slice(run, gap):
    """sub_runs of `run` where P2_IN sat at drift = mesh + gap."""
    h = hvset[(hvset.run == run) & (hvset.det == 'P2_IN')].copy()
    h['gap'] = h.drift - h.mesh_or_resist
    return set(h.loc[h.gap == gap, 'sub_run'])


fig, ax = plt.subplots(figsize=(9.0, 5.8))

# det4, 24-26 Jul -- uRWELL-referenced, drift = mesh + 200 V
g = m[m.station == 'P2_IN'].sort_values('mesh_hv')
ax.errorbar(g.mesh_hv, g.eff, yerr=[g.eff - g.lo, g.hi - g.eff],
            color=st.C_MID, marker='o', lw=2.2, capsize=2.5,
            label=ch.label('P2_IN', 'drift_mesh_scan_1', window=True)
                  + ' — uRWELL-referenced, drift = mesh + 200 V')

# det2, 26-27 Jul -- tag-probe, same 200 V gap out of the 7x7 scan
keep = _gap_slice('drift_mesh_2d_2', 200)
t2 = tp[(tp.run == 'drift_mesh_2d_2') & (tp.probe == 'P2_IN')
        & tp.sub_run.isin(keep)].sort_values('hv')
if len(t2):
    ax.plot(t2.hv, t2.eff_corr, color=st.C_OUT, marker='s', lw=2.2, ls='--',
            label=ch.label('P2_IN', 'drift_mesh_2d_2', window=True)
                  + ' — tag-and-probe, drift = mesh + 200 V')

# det5 = the CERN-built chamber, from 28 Jul -- uRWELL-referenced
rep5 = _series('P2_IN', [inA, inB])
if len(rep5):
    ax.errorbar(rep5.mesh_hv, rep5.eff,
                yerr=[rep5.eff - rep5.lo, rep5.hi - rep5.eff],
                color=st.C_VIO, marker='^', lw=2.2, ls='-.', capsize=2.5,
                label=ch.label('P2_IN', 'p2in_hvrange_1', window=True)
                      + ' — uRWELL-referenced, drift = mesh + 300 V')

ax.set_xlabel('mesh voltage [V]')
ax.set_ylabel('efficiency')
ax.set_ylim(0, 1.02)
ax.axhline(0.95, color=st.GRID, lw=1, ls=':')
ax.legend(loc='upper left', fontsize=10.5)
ax.set_title('One station, three chambers — everything that sat at P2_IN\n'
             'det2 (22–23 and 26–27 Jul; dead on the first stint, repaired '
             'for the second), det4 (24–26 Jul,\nleaky drift frame) and '
             'det5 = the CERN-built chamber (from 28 Jul)', fontsize=11.5)
st.finish(fig, f'{OUT}/eff_p2in_swap.png')
print('done')
