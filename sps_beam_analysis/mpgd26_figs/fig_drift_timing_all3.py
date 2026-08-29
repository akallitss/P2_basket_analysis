#!/usr/bin/env python3
"""Timing resolution vs drift field with ALL THREE stations, on both arms.

The two published drift-timing figures each show only two stations, for two
different reasons, and neither is a data limitation:

  DREAM  `timing_vs_drift_magboltz.png` uses drift_mesh_scan_1, where the drift
         half PARKS P2_IN at a fixed 430/630 -- it has no drift axis in that
         run at all (sigma sits at 26.8 ns for all eleven points).
  VMM    `vmm_timing.png` panel (a) omits P2_IN as "raised thresholds, few
         captures, fits unstable".

Two DREAM runs did sweep the drift on all three stations at once, and the VMM
by-HV table does carry P2_IN, so both arms can be shown complete:

  D1  mesh_drift_scan_up_1  Sun 26 Jul, mesh 450 V on every station, drift
      460 -> 875 V in 25 steps: the full turn-on, from the field being too low
      to pull the primaries out of the drift gap to the plateau.
  D2  p2_mesh_drift_eff_1   Tue 28 Jul, the same sweep after the P2_IN swap to
      det5 (the CERN-built chamber), coarser.  P2_IN ran at mesh 440 with its
      drift 10 V lower, so the drift-mesh GAP is matched across the three
      stations point for point (the sub_run names g150...g450 are that gap).
      The eleven eff_* sub_runs are one fixed working point held for the
      stability run, so they are collapsed to their median.
  V   vmm_trigger_timing_by_hv.csv at mesh 450, both gases.  P2_IN is kept.
      The two-capture gas-A point at 600 V (sigma 6.5 ns with contrast 5 and
      zero efficiency -- a fit to noise) is dropped by the same n >= 5 cut the
      published figure uses.

Products (into the report's figs/, gathered to the deck from there):
    timing_vs_drift_all3.png       DREAM, two runs, three stations each
    vmm_timing_vs_drift_all3.png   VMM, two gases, three stations each

Usage:  python3 fig_drift_timing_all3.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2style as st
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import chamber_history as ch

from paths import RD, OUT, REPO  # noqa: E402

DETS = ['P2_IN', 'P2_MID', 'P2_OUT']
GOAL = 20.0                       # ns, the P2 timing goal
MIN_CAPT = 5                      # VMM: below this the Gaussian fit is noise
VMM_BY_HV = os.path.join(REPO, 'sps_beam_analysis', 'vmm_dream_matching',
                         'vmm_trigger_timing_by_hv.csv')
# Every mesh-450 VMM capture that carries a P2_IN drift point was taken
# 30 Jul - 3 Aug, so P2_IN is det5 for the whole VMM panel; any VMM run number
# resolves to that, and run_36 is the nominal working point.
VMM_RUN = 'run_36'
# below ~555 V the DREAM sigma saturates at the readout-window scale -- the
# drift field is too low for the arrival spread to fit in the window at all.
# The points are kept (they are the turn-on), but no longer called out.

# Marker/colour contract shared by all three figures in this file, so the set
# reads as one set: the chamber sets the colour (det5, the CERN-built one, gets
# its own hue instead of inheriting P2_IN's), and the marker says which readout
# and which campaign a point came from.
#   circles + solid   DREAM, 26 Jul      triangles + dashed  DREAM, 28 Jul
#   squares + solid   VMM3a
DET5_COLOR = st.C_VIO
DREAM26 = dict(ls='-', marker='o', ms=6.5)
DREAM28 = dict(ls='--', marker='^', ms=8.0)
VMM = dict(ls='-', marker='s', ms=7.5)


def det_color(det, run):
    """Station hue, except at P2_IN, where the chamber changed mid-campaign."""
    if det == 'P2_IN' and ch.chamber('P2_IN', run) == 'det5':
        return DET5_COLOR
    return st.DET_COLOR[det]


def log_sigma_axis(ax, ticks):
    ax.set_yscale('log')
    ax.set_yticks(ticks)
    ax.set_yticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.axhline(GOAL, color=st.C_RED, lw=1.4, ls=':')


def goal_note(ax, x=0.02):
    ax.annotate(f'P2 goal: {GOAL:.0f} ns', xy=(x, GOAL),
                xycoords=('axes fraction', 'data'), xytext=(0, 6),
                textcoords='offset points', fontsize=11, color=st.C_RED,
                fontweight='bold')


# =============================================================== DREAM ======
def dream_curve(g, det):
    """(drift V, sigma) for one station, sorted, with repeats collapsed.

    p2_mesh_drift_eff_1 holds one setpoint for eleven stability sub_runs; a
    median over them is the honest single point (they scatter by 0.3 ns).
    """
    d = pd.DataFrame({'dv': g[f'drift_v_{det}'].astype(float),
                      'mv': g[f'mesh_v_{det}'].astype(float),
                      'sg': g[f'{det}_sigma'].astype(float)}).dropna()
    d = d[d.dv > d.mv]                       # drop any zero/reverse-field point
    if not len(d):
        return np.array([]), np.array([]), np.nan
    m = d.groupby('dv').sg.median().sort_index()
    return m.index.to_numpy(), m.to_numpy(), float(d.mv.iloc[0])


def fig_dream(tm):
    """Both drift runs on one axis.

    Colour stays station identity for the two chambers that never moved, so the
    28 Jul triangles land on top of the 26 Jul circles and the retrace is the
    thing you see.  P2_IN is the exception: it is a DIFFERENT CHAMBER on the two
    dates (det4 -> det5, the CERN-built one), so it gets its own colour rather
    than pretending to be a repeat of itself.
    """
    RUNS = [('mesh_drift_scan_up_1', '26 Jul', DREAM26, 1.0),
            ('p2_mesh_drift_eff_1', '28 Jul', DREAM28, 0.95)]
    END_LABEL = {('P2_IN', 'mesh_drift_scan_up_1'): ('det4', 9),
                 ('P2_IN', 'p2_mesh_drift_eff_1'): ('det5 = CERN', 0),
                 ('P2_MID', 'mesh_drift_scan_up_1'): ('P2_MID', -10),
                 # hung off the 28 Jul curve, whose last point is 25 V
                 # further right: on the 26 Jul one it collides with the green
                 # dashes, and lifting it far enough hits the 20 ns line
                 ('P2_OUT', 'p2_mesh_drift_eff_1'): ('P2_OUT', 4)}

    fig, ax = plt.subplots(figsize=(10.8, 6.6))
    for run, when, style, alpha in RUNS:
        g = tm[(tm.run == run) & (tm.axis == 'drift')]
        for det in DETS:
            x, y, mv = dream_curve(g, det)
            if len(x) < 2:
                continue
            c = det_color(det, run)
            ax.plot(x, y, color=c, lw=2.2, alpha=alpha, **style,
                    label=f'{ch.label(det, run)} — {when}, mesh {mv:.0f} V')
            lab = END_LABEL.get((det, run))
            if lab:
                st.direct_label(ax, x[-1], y[-1], lab[0], c, dy=lab[1])

    log_sigma_axis(ax, [15, 20, 30, 50, 100, 200])
    goal_note(ax)
    ax.set_xlim(440, 960)
    ax.set_xlabel('drift voltage [V]')
    ax.set_ylabel('single-station time resolution $\\sigma$ [ns]')
    ax.legend(loc='upper right', fontsize=10.5, ncol=2, columnspacing=1.0)
    ax.set_title(
        'Timing vs drift HV, all three stations, both drift campaigns\n'
        'circles 26 Jul, triangles 28 Jul: P2_MID (det1) and P2_OUT (det3) are '
        'the same chambers on both dates and retrace each other;\n'
        'only P2_IN changed hardware, det4 → det5 = the CERN-built chamber, '
        'which ran at mesh 440 V so its drift$-$mesh gap still matches',
        fontsize=11.5)
    st.finish(fig, f'{OUT}/timing_vs_drift_all3.png')


# ================================================================= VMM ======
def fig_vmm():
    d = pd.read_csv(VMM_BY_HV)
    d = d[(d.mesh == 450) & (d.n_captures >= MIN_CAPT)]
    panels = [(st.GAS_A, '(a) gas A — Ar/CO$_2$/iC$_4$H$_{10}$ 93/5/2'),
              (st.GAS_B, '(b) gas B — Ar/CF$_4$/iC$_4$H$_{10}$ 88/10/2')]
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.6), sharey=True)
    for ax, (gas, title) in zip(axes, panels):
        for det in DETS:
            s = d[(d.station == det) & (d.gas == gas)].sort_values('drift')
            if len(s) < 2:
                continue
            c = det_color(det, VMM_RUN)
            ax.plot(s.drift, s.sigma_ns, color=c, lw=2.2, ls=VMM['ls'],
                    zorder=2, label=ch.label(det, VMM_RUN))
            ax.plot(s.drift, s.sigma_ns, ls='none', marker=VMM['marker'],
                    ms=VMM['ms'], color=c, zorder=3)
            st.direct_label(ax, s.drift.iloc[-1], s.sigma_ns.iloc[-1], det, c,
                            dy={'P2_IN': 0, 'P2_MID': -9, 'P2_OUT': 9}[det])
        ax.axhline(GOAL, color=st.C_RED, lw=1.4, ls=':')
        goal_note(ax)
        ax.set_xlim(575, 900)
        ax.set_ylim(0, 105)
        ax.set_xlabel('drift voltage [V] (mesh 450 V)')
        ax.set_title(title, fontsize=12.5)
        ax.legend(loc='upper right', fontsize=11)
    axes[0].set_ylabel('coincidence width $\\sigma$ [ns]')
    axes[1].text(0.97, 0.06, 'gas B was never scanned above 750 V',
                 transform=axes[1].transAxes, fontsize=10.5, color=st.TEXT2,
                 ha='right')
    fig.suptitle('VMM3a trigger-referenced timing vs drift field, all three '
                 'stations\nP2_IN (det5, the CERN chamber) ran at raised '
                 'thresholds — it is the widest, but it does follow the field',
                 fontsize=12.5, y=0.995)
    fig.tight_layout(rect=(0.01, 0, 1, 0.945))
    st.finish(fig, f'{OUT}/vmm_timing_vs_drift_all3.png', tight=False)


# ==================================================== DREAM vs VMM ==========
def fig_compare(tm):
    """The two readouts side by side on one axis, same chambers, same gas.

    DREAM is the 28 Jul drift sweep (p2_mesh_drift_eff_1), the same run as the
    triangles in timing_vs_drift_all3; VMM is gas A at mesh 450, which is the
    same Ar/CO2/iso 93/5/2 the July DREAM data was taken on (the CF4 changeover
    is run_60, 2 Aug).  Both arms therefore see det1 / det3 / det5 in the same
    three slots at the same working point.

    The two sigmas are close cousins, not the same number, and the legend says
    which is which: DREAM's is pair-derived, so the trigger jitter cancels out
    of it, while the VMM one is referenced to the trigger and still carries the
    trigger channel.  That is the honest reason the VMM curves sit a little
    higher, and it is left as a legend label rather than an argument on the
    slide.
    """
    dr = tm[(tm.run == 'p2_mesh_drift_eff_1') & (tm.axis == 'drift')]
    vm = pd.read_csv(VMM_BY_HV)
    vm = vm[(vm.mesh == 450) & (vm.gas == st.GAS_A)
            & (vm.n_captures >= MIN_CAPT)]

    LW = 2.4

    fig, ax = plt.subplots(figsize=(10.8, 6.6))
    for det in DETS:
        c = det_color(det, 'p2_mesh_drift_eff_1')
        x, y, _ = dream_curve(dr, det)
        if len(x):
            ax.plot(x, y, color=c, lw=LW, zorder=3, **DREAM28,
                    label=f'{ch.label(det, "p2_mesh_drift_eff_1")} — DREAM')
        s = vm[vm.station == det].sort_values('drift')
        if len(s):
            ax.plot(s.drift, s.sigma_ns, color=c, ls=VMM['ls'], lw=LW,
                    zorder=2, label=f'{ch.label(det, VMM_RUN)} — VMM3a')
            ax.plot(s.drift, s.sigma_ns, ls='none', marker=VMM['marker'],
                    ms=VMM['ms'], color=c, zorder=3)

    log_sigma_axis(ax, [15, 20, 30, 50, 100])
    goal_note(ax)
    ax.set_xlim(575, 925)
    ax.set_xlabel('drift voltage [V] (mesh 450 V, Ar/CO$_2$/iC$_4$H$_{10}$ '
                  '93/5/2)')
    ax.set_ylabel('time resolution $\\sigma$ [ns]')
    ax.legend(loc='upper right', fontsize=10.5, ncol=2, columnspacing=1.0,
              title='dashed + triangles: DREAM   |   solid + squares: VMM3a',
              title_fontsize=10.5)
    ax.set_title('DREAM and VMM3a on the same three chambers, same gas, '
                 'same drift sweep', fontsize=13)
    # under the axes, not inside them: at 15-20 ns the plot floor is where the
    # det1 curves live and the line ran straight through them
    fig.text(0.5, -0.02, 'DREAM $\\sigma$: pair-derived, trigger jitter '
             'cancels.    VMM3a $\\sigma$: trigger-referenced, trigger '
             'channel included.',
             ha='center', va='top', fontsize=10.5, color=st.TEXT2)
    st.finish(fig, f'{OUT}/timing_vs_drift_dream_vs_vmm.png')


if __name__ == '__main__':
    tm = pd.read_csv(f'{RD}/dream_timing_scans.csv')
    fig_dream(tm)
    fig_vmm()
    fig_compare(tm)
