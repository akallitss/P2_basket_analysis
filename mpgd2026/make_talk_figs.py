#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_talk_figs.py -- the figures the MPGD26 deck needs and that no pipeline
stage produces yet.  Everything else in the deck is copied from an existing
stage/report output; only these five are new, and each closes a gap the
conference roadmap lists as open:

  bench_beam_mesh        WP-C: the cosmic bench and the beam on one eps-vs-mesh axis
  dream_vs_vmm           WP-B: DREAM and VMM efficiency, same uRWELL reference
  vmm_threshold          WP-B: why the VMM number is lower -- the discriminator
  snr_matrix             WP-E: the Nov-2025 SNR(gain, peaking) matrix, frozen
  timing_campaigns       WP-G: the timing ladder across all three campaigns

Titles carry the conditions, because a slide has no caption.

Usage:  python3 make_talk_figs.py [-o OUTDIR]
"""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- paths ---- #
# the bench pipeline owns the pad map and the insulation-mask pillars
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'cosmic_bench_analysis'))

BENCH = ('/local/home/ak271430/Documents/PostDocSaclay/data/'
         'Cosmic_Bench/Analysis')
BEAM = ('/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/'
        'SPS_Beam_Test/TB_July2026_H4/analysis/urw_referenced_efficiency')
VMM = ('/local/home/ak271430/Documents/PostDocSaclay/P2_basket_analysis/'
       'sps_beam_analysis/vmm_dream_matching')
NOV = ('/local/home/ak271430/Documents/PostDocSaclay/data/'
       'SPS_Beam_Test/VMM-alinx-data')
WS = ('/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/'
      'mpgd26_workspace')
RD = f'{WS}/report_data'
HVSET = f'{WS}/eos_inventory/hv_setpoints.csv'

# Station / detector colours, kept identical to the pipeline figures so the
# deck reads as one system.
C = {'P2_IN': '#1f77b4', 'P2_MID': '#ff7f0e', 'P2_OUT': '#2ca02c',
     'det1': '#1f77b4', 'det2': '#2ca02c', 'det3': '#d62728',
     'det4': '#ff7f0e'}
SIZE = (9.2, 5.6)

# Chamber identity across the two campaigns (hardware logbook, 26 + 28 Jul;
# full table in sps_beam_analysis/chamber_history.py).  ALL FOUR Saclay
# chambers travelled to H4.  P2_MID = det1 and P2_OUT = det3 for the whole
# test, never touched; every swap happened at P2_IN -- det2 (22-23 Jul, ~20x
# low signal from a drift HV contact), det4 (24-25 Jul, leaky drift frame),
# det2 again (26-27 Jul, contact believed fixed), then the CERN-built chamber
# from 28 Jul.
BENCH_TO_STATION = {'det1': 'P2_MID', 'det3': 'P2_OUT'}
BENCH_LABEL = {
    'det1': 'bench det1 (= beam P2_MID) — cosmics, M3 tracks, '
            'Ar/iC$_4$H$_{10}$ 95/5',
    'det2': 'bench det2 (= beam P2_IN, 22–23 and 26–27 Jul) — cosmics, '
            'M3 tracks',
    'det3': 'bench det3 (= beam P2_OUT) — cosmics, M3 tracks, '
            'Ar/iC$_4$H$_{10}$ 95/5',
    'det4': 'bench det4 (= beam P2_IN, 24–25 Jul — leaky drift frame) — '
            'cosmics, M3 tracks',
}

# Slide typography (2026-08-18): projected, not read at arm's length -- axis
# labels large AND bold everywhere, matching sps_beam_analysis/mpgd26_figs/
# p2style.py so the whole deck reads as one system.
plt.rcParams.update({
    'font.size': 13, 'axes.titlesize': 14, 'axes.labelsize': 16,
    'axes.labelweight': 'bold', 'axes.titleweight': 'bold',
    'xtick.labelsize': 13.5, 'ytick.labelsize': 13.5,
    'legend.fontsize': 12, 'figure.facecolor': 'white',
    'figure.titlesize': 15, 'figure.titleweight': 'bold',
    'axes.grid': True, 'grid.alpha': 0.3, 'axes.axisbelow': True,
})


def save(fig, out, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(out, f'{name}.{ext}'), dpi=170,
                    bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  [{name}]')


# ==================================================================== 1 ==== #
def fig_bench_beam_mesh(out):
    """eps vs mesh on the bench (cosmics/DREAM/M3) and on the beam
    (muons/DREAM/uRWELL).

    Chamber identity, from the 2026-07-28 logbook: the beam telescope ran
    P2_MID = det1 and P2_OUT = det3 for the WHOLE test, and both are the same
    physical chambers the cosmic bench characterised.  So bench det1 vs beam
    P2_MID is a like-for-like pair; only the gas and the probe differ.  det2
    (drawn here) went to the beam too, at P2_IN."""
    bench = {
        'det1': (f'{BENCH}/det1/p2_det1_long_run_mesh_scan_7-19-26/mesh_scan/'
                 '11_hv_scan_efficiency/'
                 'efficiency_vs_hv_without_connectors_1_2_10_spark_vetoed.csv'),
        'det2': (f'{BENCH}/det2/p2_det1_det2_long_run_mesh_scan_7-9-26/hv_scan/'
                 '11_hv_scan_efficiency/'
                 'efficiency_vs_hv_without_connectors_1_8_9_10_spark_vetoed.csv'),
    }
    fig, ax = plt.subplots(figsize=SIZE)

    for det, path in bench.items():
        if not os.path.exists(path):
            print(f'  ! missing {path}')
            continue
        d = pd.read_csv(path).sort_values('hv')
        ax.errorbar(d['hv'], d['eff_reco'], yerr=d['eff_reco_err'],
                    marker='s', ms=5, lw=1.6, ls='--', color=C[det],
                    label=BENCH_LABEL[det])

    b = pd.read_csv(f'{BEAM}/drift_mesh_scan_1/'
                    'urw_p2_efficiency_drift_mesh_scan_1.csv')
    # The CSV tags every row scan_axis='mesh'; the mesh scan proper is the
    # meshscan_* sub_runs, the drift_* ones hold mesh fixed at 450 V.
    b = b[b['sub_run'].str.startswith(('meshscan', 'nominal'))]
    for st in ('P2_IN', 'P2_MID', 'P2_OUT'):
        s = b[b['station'] == st].sort_values('mesh_hv')
        if not len(s):
            continue
        ax.plot(s['mesh_hv'], s['eff'], marker='o', ms=5.5, lw=2.2,
                color=C[st],
                label=f'beam {st} — SPS muons, uRWELL tracks')

    ax.axhline(0.95, color='0.5', ls=':', lw=1.2)
    ax.text(0.012, 0.952, '95 %', transform=ax.get_yaxis_transform(),
            fontsize=9.5, color='0.35', va='bottom')
    ax.annotate('working point\n(top of the scanned range —\nno plateau reached)',
                xy=(449, 0.955), xytext=(432, 0.60), fontsize=9.5,
                ha='center', color='#333',
                arrowprops=dict(arrowstyle='->', color='#777', lw=1.2),
                bbox=dict(fc='white', ec='#ccc', boxstyle='round,pad=0.35'))
    ax.set_xlabel('mesh voltage [V]')
    ax.set_ylabel('efficiency')
    ax.set_ylim(0.25, 1.02)
    ax.set_title('The bench predicted the beam — efficiency vs mesh voltage\n'
                 'bench det1 IS the beam\'s P2_MID chamber: same hardware, '
                 'cosmics + M3 vs SPS muons + uRWELL,\n'
                 'Ar/iC$_4$H$_{10}$ 95/5 vs Ar/CO$_2$/iC$_4$H$_{10}$ 93/5/2 — '
                 'the gas shifts the curve, the shape and the working point do '
                 'not move', fontsize=11.5)
    ax.legend(loc='lower right', framealpha=0.95, fontsize=9.5)
    save(fig, out, 'bench_beam_mesh')


# =================================================================== 1b ==== #
def fig_bench_beam_drift(out):
    """The transport half of WP-C, and the strongest form of it: bench det1 and
    det3 are the very chambers that ran the beam as P2_MID and P2_OUT (logbook
    26 + 28 Jul; neither station was touched during the test).  So all four
    curves here are two chambers measured twice, in two labs, with two probes
    and two gases.  The x axis has to be the drift FIELD rather than the
    voltage, because the campaigns sit at different mesh settings (415/420 V on
    the bench, 450 V on the beam) and the drift electrode reads mesh + gap in
    both."""
    GAP_CM = 0.3          # P2 drift gap, 3 mm (bench CSVs: dV/E = 0.3 cm)
    bench = {
        'det1': (f'{BENCH}/det1/p2_det1_drift_scan_7-19-26/drift_scan/'
                 '16_drift_scan_efficiency/'
                 'efficiency_vs_drift_without_connectors_1_2_10_spark_vetoed.csv'),
        'det3': (f'{BENCH}/det3/p2_det3_det4_drift_scan_7-16-26/drift_scan/'
                 '16_drift_scan_efficiency/'
                 'efficiency_vs_drift_without_connectors_1_8_9_10_spark_vetoed.csv'),
    }
    fig, ax = plt.subplots(figsize=SIZE)

    for det, path in bench.items():
        if not os.path.exists(path):
            print(f'  ! missing {path}')
            continue
        d = pd.read_csv(path).sort_values('drift')
        e = (d['drift'] - d['mesh']) / GAP_CM
        ax.errorbar(e, d['eff_reco'], yerr=d['eff_reco_err'],
                    marker='s', ms=5, lw=1.6, ls='--', color=C[det],
                    label=f'bench {det} (= beam {BENCH_TO_STATION[det]}) — '
                          f'cosmics, M3 tracks, mesh '
                          f'{int(d["mesh"].iloc[0])} V')

    b = pd.read_csv(f'{BEAM}/drift_mesh_scan_1/'
                    'urw_p2_efficiency_drift_mesh_scan_1.csv')
    # the drift_* sub_runs hold mesh at 450 V; drift_450 is zero drift field
    b = b[b['sub_run'].str.startswith('drift_')]
    for st in ('P2_MID', 'P2_OUT'):          # P2_IN was parked at 430/630
        s = b[b['station'] == st].sort_values('drift_hv')
        if not len(s):
            continue
        e = (s['drift_hv'] - s['mesh_hv']) / GAP_CM
        ax.errorbar(e, s['eff'], yerr=[s['eff'] - s['lo'], s['hi'] - s['eff']],
                    marker='o', ms=5.5, lw=2.2, color=C[st],
                    label=f'beam {st} — SPS muons, uRWELL tracks, mesh 450 V')

    ax.axvspan(1000, 1167, color='#999', alpha=0.13, zorder=0)
    ax.text(1083, 0.32, 'beam working point\n700 V drift / 450 V mesh\n'
            '(833 V/cm) sits just below\nthe 750–800 V plateau',
            fontsize=8.8, ha='center', va='center', color='#444')
    ax.axhline(0.95, color='0.5', ls=':', lw=1.2)
    ax.set_xlabel('drift field [V/cm]   (3 mm gap; electrode = mesh + gap V)')
    ax.set_ylabel('efficiency')
    ax.set_ylim(0.0, 1.02)
    ax.set_title('The bench predicted the beam — the transport half\n'
                 'TWO chambers measured twice: det1 = P2_MID and det3 = '
                 'P2_OUT, on the Saclay bench and at H4\n'
                 'sharp rise to ~0.9 by 170 V/cm, saturation by ~1000 V/cm, '
                 'and the two campaigns agree to 2–4 points\nat the plateau '
                 'despite different gas and different tracking',
                 fontsize=11.5)
    ax.legend(loc='lower right', framealpha=0.95, fontsize=9.5)
    save(fig, out, 'bench_beam_drift')


# =================================================================== 1c ==== #
def fig_fe55_bench_beam(out):
    """PARKED 2026-08-17 -- the chamber identity underneath it is contradictory.

    The premise was: the 7-18-26 Fe55 scan measured the chamber that then ran as
    P2_OUT.  That came from `p2_qa_config.py:624` ("det2 = P2_OUT ... det3 =
    P2_MID"), which the 2026-07-28 logbook contradicts: the beam ran P2_MID =
    det1 and P2_OUT = det3, and det2 never left Saclay (leaky drift frame).
    Under the logbook the healthy Fe55 curve belongs to a chamber that never saw
    the beam, and the dead one (det3) is the beam's P2_OUT -- which would need a
    repair between 18 and 23 Jul to explain 96 % in the beam.  Not drawn until
    the hardware record settles which chamber sat on the Fe55 bench that day.
    Kept because the analysis is right once the label is; call it explicitly to
    regenerate."""
    fe = pd.read_csv(f'{BENCH}/det2/p2_fe55_det2_det3_mesh_scan_7-18-26/'
                     'fe55_scan/18_fe55_spectra/fe55_gain_vs_hv_spark_vetoed.csv')
    fe = fe[np.isfinite(fe['peak_adc']) & (~fe['near_floor'])].sort_values('hv')

    fig, (ax, axr) = plt.subplots(1, 2, figsize=(15.0, 5.6),
                                  gridspec_kw=dict(width_ratios=[1.35, 1],
                                                   wspace=0.42))
    ax.errorbar(fe['hv'], fe['peak_adc'], yerr=fe['peak_err'], marker='s',
                ms=6, lw=1.8, ls='--', color='#7f4fc4',
                label='bench: Fe$^{55}$ photopeak, DREAM self-trigger '
                      '(18 Jul, Saclay)')
    ax.set_yscale('log')
    ax.set_ylabel('Fe$^{55}$ 5.9 keV photopeak [ADC]  $\\propto$ gas gain',
                  color='#7f4fc4')
    ax.tick_params(axis='y', labelcolor='#7f4fc4')

    # gain e-folding: peak = A exp(V / V0)
    sl, ic = np.polyfit(fe['hv'], np.log(fe['peak_adc']), 1)
    v0 = 1.0 / sl
    xs = np.linspace(fe['hv'].min() - 5, fe['hv'].max() + 5, 50)
    ax.plot(xs, np.exp(ic + sl * xs), color='#7f4fc4', lw=1.0, alpha=0.5)
    ax.text(0.03, 0.93, f'gain e-folds every {v0:.1f} V\n'
                        f'(doubles every {v0 * np.log(2):.1f} V)',
            transform=ax.transAxes, fontsize=10.5, color='#7f4fc4',
            va='top', bbox=dict(fc='white', ec='#ddd',
                                boxstyle='round,pad=0.3'))

    ax2 = ax.twinx()
    b = pd.read_csv(f'{BEAM}/drift_mesh_scan_1/'
                    'urw_p2_efficiency_drift_mesh_scan_1.csv')
    b = b[b['sub_run'].str.startswith(('meshscan', 'nominal'))]
    s = b[b['station'] == 'P2_OUT'].sort_values('mesh_hv')
    ax2.errorbar(s['mesh_hv'], s['eff'],
                 yerr=[s['eff'] - s['lo'], s['hi'] - s['eff']],
                 marker='o', ms=6, lw=2.4, color=C['P2_OUT'],
                 label='beam: same chamber as P2_OUT, uRWELL-referenced '
                       'efficiency (26 Jul, H4)')
    low = f'{BEAM}/low_mesh_scan_1/urw_p2_efficiency_low_mesh_scan_1.csv'
    if os.path.exists(low):
        lo = pd.read_csv(low)
        lo = lo[lo['station'] == 'P2_OUT'].sort_values('mesh_hv')
        ax2.errorbar(lo['mesh_hv'], lo['eff'],
                     yerr=[lo['eff'] - lo['lo'], lo['hi'] - lo['eff']],
                     marker='o', ms=6, lw=2.4, ls=':', color=C['P2_OUT'],
                     label='beam: low-mesh extension (low_mesh_scan_1)')
    ax2.set_ylabel('absolute efficiency (uRWELL-referenced)',
                   color=C['P2_OUT'])
    ax2.tick_params(axis='y', labelcolor=C['P2_OUT'])
    ax2.set_ylim(0, 1.03)
    ax2.grid(False)

    ax.set_xlabel('mesh voltage [V]')
    ax.set_xlim(360, 455)
    ax.set_title('(a) the two campaigns on one voltage axis', fontsize=11.5)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc='lower right', fontsize=8.6,
              framealpha=0.95)

    # (b) eliminate the voltage: beam efficiency against the gain the bench
    # measured at the same mesh setting.  Nothing is fitted across campaigns --
    # the bench gain is interpolated in log space at each beam scan point.
    beam = pd.concat([s, lo] if os.path.exists(low) else [s], ignore_index=True)
    beam = beam.sort_values('mesh_hv')
    ok = ((beam['mesh_hv'] >= fe['hv'].min()) &
          (beam['mesh_hv'] <= fe['hv'].max()))
    bb = beam[ok]
    gain = np.exp(np.interp(bb['mesh_hv'], fe['hv'], np.log(fe['peak_adc'])))
    axr.errorbar(gain, bb['eff'],
                 yerr=[bb['eff'] - bb['lo'], bb['hi'] - bb['eff']],
                 marker='o', ms=7, lw=2.2, color=C['P2_OUT'])
    for g_, e_, v_ in zip(gain, bb['eff'], bb['mesh_hv']):
        axr.annotate(f'{v_:.0f} V', xy=(g_, e_), xytext=(6, -10),
                     textcoords='offset points', fontsize=8.5, color='#555')
    axr.set_xscale('log')
    axr.set_xlabel('gas gain the bench measured at that mesh setting\n'
                   '[Fe$^{55}$ photopeak, ADC]')
    axr.set_ylabel('beam efficiency (uRWELL-referenced)')
    axr.set_ylim(0, 1.02)
    axr.set_title('(b) the voltage divided out: efficiency vs gain',
                  fontsize=11.5)
    axr.text(0.04, 0.30,
             f'over the range the Fe$^{{55}}$ fit reaches\n'
             f'({bb.mesh_hv.min():.0f}–{bb.mesh_hv.max():.0f} V), gain '
             f'$\\times${gain.max() / gain.min():.1f} takes the beam\n'
             f'efficiency {bb.eff.iloc[0]:.2f} $\\rightarrow$ '
             f'{bb.eff.iloc[-1]:.2f}. Below 380 V the photopeak\n'
             'sits in the noise bulge and is not fitted.',
             transform=axr.transAxes, fontsize=9.5, va='top', color='#444',
             bbox=dict(fc='white', ec='#ddd', boxstyle='round,pad=0.3'))

    fig.suptitle('The same chamber, on the bench and in the beam — '
                 'P2_OUT (bulked 25 Jun)\n'
                 'Fe$^{55}$ gain at Saclay on 18 Jul, SPS muons at H4 from '
                 '23 Jul: the bench gain curve is the beam turn-on',
                 fontsize=12.5, y=1.02)
    save(fig, out, 'fe55_bench_beam_P2_OUT')


# =================================================================== 1d ==== #
# ------------------------------------------------------- geometry helpers -- #
def _board_to_bench():
    """Similarity transform from the Gerber board frame (pad_cx/pad_cy, the
    frame the beam per-pad maps live in) into the bench M3 reference frame
    (the frame the bench efficiency map lives in).

    The two pillar tables -- p2_mapping.load_pillars() in the board frame and
    the pipeline's pillars_m3 CSV -- are the SAME 11,683 circles row for row,
    so the map is over-determined by four orders of magnitude: fitting a full
    affine to it closes to ~1e-13 mm and comes out a pure rotation (89.25 deg)
    plus the alignment's 0.9743 scale, no reflection.  That is what lets the
    beam map be drawn on top of the bench map instead of beside it.
    """
    import p2_qa_config as _qa
    import p2_mapping as _pm
    brd = _pm.load_pillars(_qa.MASK_GBR_PATH)
    m3 = pd.read_csv(f'{BENCH}/det1/p2_det1_long_run_efficiency_7-19-26/'
                     'long_run_det1_415_615/06_efficiency/'
                     'pillars_m3_without_connectors_1_2_10_spark_vetoed.csv')
    if not len(brd) or len(brd) != len(m3):
        return None, None, None
    A = np.column_stack([brd['x'], brd['y'], np.ones(len(brd))])
    M, *_ = np.linalg.lstsq(A, np.column_stack([m3['x'], m3['y']]), rcond=None)
    res = float(np.hypot(*(A @ M - np.column_stack([m3['x'], m3['y']])).T).max())

    def xf(x, y):
        v = np.column_stack([np.asarray(x, float), np.asarray(y, float),
                             np.ones(np.size(x))]) @ M
        return v[:, 0], v[:, 1]

    scale = float(np.sqrt(abs(np.linalg.det(M[:2].T))))
    return xf, scale, res


def _pad_polys(xf=None):
    """Every mapped pad as a real rectangle (Gerber tile geometry), optionally
    carried into the bench frame.  Returns (channel_id array, (n,4,2) verts)."""
    import p2_qa_config as _qa
    import p2_mapping as _pm
    m = _pm.load_pad_map(_qa.MAP_CSV_PATH)
    a = np.radians(m['pad_angle'].to_numpy())
    w, h = m['pad_w'].to_numpy(), m['pad_h'].to_numpy()
    cx, cy = m['pad_cx'].to_numpy(), m['pad_cy'].to_numpy()
    ca, sa = np.cos(a), np.sin(a)
    verts = np.stack([np.stack([cx + dx * w * ca - dy * h * sa,
                                cy + dx * w * sa + dy * h * ca], axis=1)
                      for dx, dy in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))],
                     axis=1)
    if xf is not None:
        fx, fy = xf(verts[..., 0].ravel(), verts[..., 1].ravel())
        verts = np.stack([fx, fy], axis=1).reshape(verts.shape)
    return m.index.to_numpy(), verts


def _big_pillars(ax, x, y, r, label=None, color='#d62728'):
    from matplotlib.patches import Circle
    first = True
    for xi, yi, ri in zip(np.atleast_1d(x), np.atleast_1d(y), np.atleast_1d(r)):
        ax.add_patch(Circle((xi, yi), ri, facecolor=color, edgecolor=color,
                            zorder=6))
        ax.add_patch(Circle((xi, yi), 9.0, facecolor='none', edgecolor=color,
                            ls='--', lw=1.6, zorder=6,
                            label=(label if first else None)))
        first = False


def fig_bench_beam_maps(out):
    """Wish-list 5a.4: the efficiency map on the bench beside the map on the
    beam -- same colour scale, and (since 2026-08-18) the same FRAME.

    Two different probes on the same physical chamber (det1 = P2_MID): 18 h of
    cosmics through an M3 telescope over the whole active area, vs 1.15 M SPS
    muons through a beam spot that covers a fifth of it.  Drawing them in one
    frame is the point of the figure -- it says which part of the left panel
    the right panel is actually measuring, which the old side-by-side version,
    with two unrelated axis ranges, could not.

    Both panels are surfaces: the bench map is the 5 mm sliding window, the
    beam map is every illuminated pad filled at its true Gerber outline.  The
    five big mesh-support pillars are marked in both, in the same places.
    """
    bench_npz = (f'{BENCH}/det1/p2_det1_long_run_efficiency_7-19-26/'
                 'long_run_det1_415_615/06_efficiency/'
                 'efficiency_map_sliding_without_connectors_1_2_10_'
                 'spark_vetoed.npz')
    beam_glob = (f'{BEAM}/../P2_MID/highstat_eff_1/*/22_tag_probe_efficiency/'
                 'eff_map_P2_MID_beam_commissioning_00*.csv')
    import glob as _glob
    from matplotlib.collections import PolyCollection
    beam_csv = sorted(_glob.glob(beam_glob))
    if not os.path.exists(bench_npz) or not beam_csv:
        print('  ! bench npz or beam map missing')
        return

    z = np.load(bench_npz, allow_pickle=True)
    grid = np.where(z['counts'] >= int(z['min_rays']), z['eff_within'], np.nan)
    extent = z['extent']

    d = pd.read_csv(beam_csv[0])
    lit = d[(d['n_tag'] >= 2000) & np.isfinite(d['eff'])]

    xf, scale, res = _board_to_bench()
    if xf is None:
        print('  ! pillar correspondence unavailable — cannot co-register')
        return
    print(f'  board->bench frame: scale {scale:.4f}, closure {res:.2e} mm')

    ids, verts = _pad_polys(xf)                      # all 1280 pads, bench frame
    idx = {c: i for i, c in enumerate(ids)}
    keep = [idx[c] for c in lit['channel_id'] if c in idx]
    lit_verts = verts[keep]
    lit_eff = lit['eff'].to_numpy()[[k in idx for k in lit['channel_id']]]

    # the beam spot as one outline, for the bench panel
    from scipy.spatial import ConvexHull
    pts = lit_verts.reshape(-1, 2)
    hull = pts[np.append(ConvexHull(pts).vertices, ConvexHull(pts).vertices[0])]

    pil = pd.read_csv(f'{BENCH}/det1/p2_det1_long_run_efficiency_7-19-26/'
                      'long_run_det1_415_615/06_efficiency/'
                      'pillars_m3_without_connectors_1_2_10_spark_vetoed.csv')
    big = pil[pil['big']]

    fig, axes = plt.subplots(1, 2, figsize=(15.4, 7.4))
    vmin, vmax = 0.5, 1.0
    lim_x = (float(extent[0]) - 6, float(extent[1]) + 6)
    lim_y = (float(extent[2]) - 6, float(extent[3]) + 6)

    # ---- left: the bench, whole active area ------------------------------ #
    im = axes[0].imshow(grid.T, origin='lower', extent=extent, cmap='viridis',
                        vmin=vmin, vmax=vmax, aspect='equal', zorder=1)
    axes[0].plot(hull[:, 0], hull[:, 1], color='#ffffff', lw=3.4, zorder=5)
    axes[0].plot(hull[:, 0], hull[:, 1], color='#d62728', lw=2.0, ls='--',
                 zorder=5, label='SPS beam spot (right panel)')
    _big_pillars(axes[0], big['x'], big['y'], big['r'], label='big pillars (5)')
    axes[0].set_title('BENCH — det1, 18.3 h of cosmics, M3 tracks\n'
                      'sliding 5 mm window, mesh 415 / drift 615 V\n'
                      f'map median {np.nanmedian(grid):.3f} over the whole '
                      'active area', fontsize=12.5)

    # ---- right: the beam, same frame, same scale -------------------------- #
    # three tiers, because they are three different things: the pad plane, the
    # part of it the beam DAQ read out (connectors 3-6 = 512 pads), and the
    # part the beam actually lit.
    read = [idx[c] for c in d['channel_id'] if c in idx]
    axes[1].add_collection(PolyCollection(
        verts, facecolors='#f5f5f3', edgecolors='#e4e3df', linewidths=0.3,
        zorder=1))
    axes[1].add_collection(PolyCollection(
        verts[read], facecolors='#e2e1dc', edgecolors='#c9c8c3',
        linewidths=0.35, zorder=2,
        label=f'read out at the beam ({len(read)} pads, connectors 3–6)'))
    axes[1].add_collection(PolyCollection(
        lit_verts, array=lit_eff, cmap='viridis', clim=(vmin, vmax),
        edgecolors='face', linewidths=0.0, zorder=3))
    axes[1].plot(hull[:, 0], hull[:, 1], color='#d62728', lw=2.0, ls='--',
                 zorder=5, label=f'lit by the beam ({len(lit_verts)} pads, '
                                 f'{100 * len(lit_verts) / len(read):.0f} % '
                                 'of what was read out)')
    _big_pillars(axes[1], big['x'], big['y'], big['r'], label='big pillars (5)')
    axes[1].set_title('BEAM — P2_MID (the same chamber), 1.15 M SPS muons, '
                      'uRWELL tracks\nper pad, only pads with $\\geq$ 2000 tags, '
                      'on the full pad layout\n'
                      f'mesh 450 / drift 700 V, median {lit["eff"].median():.3f} '
                      'inside the spot', fontsize=12.5)

    for ax in axes:
        ax.set_aspect('equal')
        ax.set_xlim(*lim_x); ax.set_ylim(*lim_y)
        ax.set_xlabel('detector-plane x [mm]')
        ax.legend(loc='lower left', fontsize=11.5, framealpha=0.92)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel('detector-plane y [mm]')

    cb = fig.colorbar(im, ax=axes, fraction=0.021, pad=0.02)
    cb.set_label('efficiency', fontsize=16, fontweight='bold')
    cb.ax.tick_params(labelsize=13)
    fig.suptitle('The same chamber, two probes, one frame — wish-list 5a.4\n'
                 f'the beam only ever lights the red patch '
                 f'({100 * len(lit_verts) / len(verts):.0f} % of the pad '
                 'plane); the bench measures all of it\n'
                 'which is why the bench pays the pillars as a 3.2 % no-hit '
                 'floor and the narrow spot only 0.07 points',
                 fontsize=14.5, y=1.10)
    save(fig, out, 'bench_beam_maps')


# =================================================================== 1e ==== #
# The extensive bench<->beam comparison: every chamber, both campaigns.
#
# Five chambers went through the two campaigns and not one of them was measured
# both ways in every axis, so the honest figure shows the gaps as well as the
# overlaps:
#
#            bench mesh   bench drift   beam mesh   beam drift
#   det1        yes          yes           yes         yes
#   det2        yes           -            yes*         -      (* tag-probe)
#   det3        yes          yes           yes         yes
#   det4         -           yes           yes          -      (P2_IN parked)
#   det5         -            -            yes          -      (CERN-built)
#
# Beam efficiencies are uRWELL-referenced wherever that product exists; det2's
# beam curve comes from 2-of-3 tag-and-probe, which the method-comparison
# figure shows agrees with the uRWELL reference to 1-2 points.
# --------------------------------------------------------------------------- #
CHAMBER_C = {'det1': '#1f77b4', 'det2': '#2ca02c', 'det3': '#d62728',
             'det4': '#ff7f0e', 'det5': '#7d3ac1'}
CHAMBERS = ['det1', 'det2', 'det3', 'det4', 'det5']
BENCH_MESH_KEY = {'det1': 'det1_meshscan1', 'det2': 'det2_hvscan1',
                  'det3': 'det3_meshscan1'}
BENCH_DRIFT_KEY = {'det1': 'det1_driftscan2', 'det3': 'det3_driftscan1',
                   'det4': 'det4_driftscan1'}
BEAM_STATION = {'det1': 'P2_MID', 'det2': 'P2_IN', 'det3': 'P2_OUT',
                'det4': 'P2_IN', 'det5': 'P2_IN'}
# caveats that belong on the series, not in a caption nobody reads
BENCH_NOTE = {'det3': ' (scan stopped at 420 V, below the plateau)',
              'det2': ' (scan starts at 395 V)'}
BEAM_RUN = {'det1': 'drift_mesh_scan_1', 'det3': 'drift_mesh_scan_1',
            'det4': 'drift_mesh_scan_1', 'det5': 'p2in_hvrange_1',
            'det2': 'drift_mesh_2d_2'}


def _ch():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), 'sps_beam_analysis'))
    import chamber_history as ch
    return ch


def _bench_scan(key, stage, stem):
    import glob as _g
    import p2_qa_config as _qa
    try:
        c = _qa.get_config(key)
    except Exception:
        return None
    f = [q for q in sorted(_g.glob(os.path.join(c.OUT_BASE, stage,
                                                f'{stem}*spark_vetoed.csv')))
         if 'last' not in os.path.basename(q)]
    return pd.read_csv(f[-1]) if f else None


def bench_mesh(det):
    k = BENCH_MESH_KEY.get(det)
    d = _bench_scan(k, '11_hv_scan_efficiency', 'efficiency_vs_hv') if k else None
    return None if d is None else d.sort_values('hv')


def bench_drift(det):
    k = BENCH_DRIFT_KEY.get(det)
    d = (_bench_scan(k, '16_drift_scan_efficiency', 'efficiency_vs_drift')
         if k else None)
    return None if d is None else d.sort_values('drift')


def _urw(run):
    p = f'{BEAM}/{run}/urw_p2_efficiency_{run}.csv'
    return pd.read_csv(p) if os.path.exists(p) else pd.DataFrame()


def beam_mesh(det):
    """Beam eps(mesh) for one chamber -> (frame, method, drift gap)."""
    stn = BEAM_STATION[det]
    if det in ('det1', 'det3'):
        f = _urw('drift_mesh_scan_1')
        parts = [f[(f.station == stn) & f.sub_run.str.startswith('meshscan')]] \
            if len(f) else []
        low = _urw('low_mesh_scan_1')
        if len(low):
            parts.append(low[low.station == stn])
        g = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        return g.sort_values('mesh_hv'), 'uRWELL-referenced', '200 V'
    if det == 'det4':
        f = _urw('drift_mesh_scan_1')
        g = f[(f.station == stn) & f.sub_run.str.startswith('meshscan')] \
            if len(f) else pd.DataFrame()
        return (g.sort_values('mesh_hv') if len(g) else g,
                'uRWELL-referenced', '200 V')
    if det == 'det5':
        parts = [f[f.station == stn] for f in
                 (_urw('p2in_hvrange_1'), _urw('p2in_hvrange_2')) if len(f)]
        g = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        return (g.sort_values('mesh_hv') if len(g) else g,
                'uRWELL-referenced', '300 V')
    if det == 'det2':
        # tag-and-probe, drift_mesh_2d_2 sliced at the same 200 V drift gap the
        # uRWELL mesh scan used, so the two curves are directly comparable
        tp = pd.read_csv(f'{RD}/dream_tag_probe.csv')
        hv = pd.read_csv(HVSET)
        h = hv[(hv.run == 'drift_mesh_2d_2') & (hv.det == stn)].copy()
        h['gap'] = h.drift - h.mesh_or_resist
        keep = set(h.loc[h.gap == 200, 'sub_run'])
        g = tp[(tp.run == 'drift_mesh_2d_2') & (tp.probe == stn) & tp.vetoed
               & tp.sub_run.isin(keep)].copy()
        g = g.rename(columns={'hv': 'mesh_hv', 'eff_lo': 'lo', 'eff_hi': 'hi'})
        g = g[np.isfinite(g.mesh_hv) & np.isfinite(g.eff)]
        return g.sort_values('mesh_hv'), 'tag-and-probe (2-of-3)', '200 V'
    return pd.DataFrame(), '', ''


def beam_drift(det):
    """Beam eps(drift). Only MID/OUT were scanned -- P2_IN sat parked."""
    stn = BEAM_STATION[det]
    if det not in ('det1', 'det3'):
        return pd.DataFrame()
    f = _urw('drift_mesh_scan_1')
    if not len(f):
        return pd.DataFrame()
    g = f[(f.station == stn) & f.sub_run.str.startswith('drift_')]
    return g.sort_values('drift_hv')


def fig_bench_beam_all(out):
    """One overlay panel with everything, plus a per-chamber grid as backup."""
    ch = _ch()

    for axis, xlab, fname in (
            ('mesh', 'mesh voltage [V]', 'bench_beam_mesh_all'),
            ('drift', 'drift voltage [V]', 'bench_beam_drift_all')):
        fig, ax = plt.subplots(figsize=(13.8, 6.8))
        n_b = n_m = 0
        for det in CHAMBERS:
            col = CHAMBER_C[det]
            b = bench_mesh(det) if axis == 'mesh' else bench_drift(det)
            if b is not None and len(b):
                x = b['hv'] if axis == 'mesh' else b['drift']
                ax.errorbar(x, b['eff_reco'], yerr=b.get('eff_reco_err'),
                            marker='s', ms=6, mfc='white', mew=1.8, lw=1.8,
                            ls='--', color=col, capsize=2,
                            label=f'{det} — cosmic bench{BENCH_NOTE.get(det, "")}')
                n_b += 1
            if axis == 'mesh':
                g, meth, gap = beam_mesh(det)
            else:
                g, meth, gap = beam_drift(det), 'uRWELL-referenced', ''
            xk = 'mesh_hv' if axis == 'mesh' else 'drift_hv'
            if len(g):
                lab = ch.label(BEAM_STATION[det], BEAM_RUN[det],
                               window=(BEAM_STATION[det] == 'P2_IN'))
                extra = '' if meth.startswith('uRWELL') else f', {meth}'
                ax.errorbar(g[xk], g['eff'],
                            yerr=[g['eff'] - g['lo'], g['hi'] - g['eff']],
                            marker='o', ms=6.5, lw=2.4, color=col, capsize=2,
                            label=f'{det} — SPS beam, {lab}{extra}')
                n_m += 1
        ax.set_xlabel(xlab)
        ax.set_ylabel('efficiency')
        ax.set_ylim(0, 1.03)
        ax.axhline(0.95, color='0.6', ls=':', lw=1.2)
        ax.grid(alpha=.3)
        ax.legend(loc='center left', bbox_to_anchor=(1.01, 0.5),
                  fontsize=11.5, framealpha=.95)
        ax.set_title(
            f'Every chamber, both campaigns — efficiency vs {axis} voltage\n'
            'open squares + dashed = cosmic bench (M3 tracks, '
            'Ar/iC$_4$H$_{10}$ 95/5)\n'
            'filled circles + solid = SPS beam '
            '(Ar/CO$_2$/iC$_4$H$_{10}$ 93/5/2) — the gas shifts the curve '
            'along the voltage axis, the shape does not move', fontsize=12)
        save(fig, out, fname)
        print(f'    {axis}: {n_b} bench series, {n_m} beam series')

    fig, axes = plt.subplots(2, 5, figsize=(24.0, 10.0), sharey=True)
    for j, det in enumerate(CHAMBERS):
        col = CHAMBER_C[det]
        for i, axis in enumerate(('mesh', 'drift')):
            ax = axes[i, j]
            b = bench_mesh(det) if axis == 'mesh' else bench_drift(det)
            if axis == 'mesh':
                g, meth, gap = beam_mesh(det)
            else:
                g, meth = beam_drift(det), 'uRWELL-referenced'
            xk = 'mesh_hv' if axis == 'mesh' else 'drift_hv'
            have = []
            if b is not None and len(b):
                x = b['hv'] if axis == 'mesh' else b['drift']
                ax.errorbar(x, b['eff_reco'], yerr=b.get('eff_reco_err'),
                            marker='s', ms=5, mfc='white', mew=1.6, lw=1.6,
                            ls='--', color=col, capsize=2, label='cosmic bench')
                have.append('bench')
            if len(g):
                ax.errorbar(g[xk], g['eff'],
                            yerr=[g['eff'] - g['lo'], g['hi'] - g['eff']],
                            marker='o', ms=5.5, lw=2.0, color=col, capsize=2,
                            label='SPS beam')
                have.append('beam')
            if not have:
                ax.text(.5, .5, 'not measured\non either setup',
                        transform=ax.transAxes, ha='center', va='center',
                        fontsize=14, color='0.45', style='italic')
            else:
                miss = {'bench', 'beam'} - set(have)
                if miss:
                    # as a legend entry, not free text -- free text lands on
                    # top of whichever curve happens to be there
                    why = ('no bench scan for this chamber' if 'bench' in miss
                           else 'P2_IN parked, never scanned' if axis == 'drift'
                           else 'no beam scan')
                    ax.plot([], [], ' ', label=why)
                ax.legend(loc='upper left', fontsize=10.5, handlelength=1.6,
                          labelspacing=.3)
            ax.set_ylim(0, 1.03)
            ax.axhline(0.95, color='0.6', ls=':', lw=1.1)
            ax.grid(alpha=.3)
            ax.set_xlabel(f'{axis} voltage [V]')
            if i == 0:
                stn = BEAM_STATION[det]
                extra = ('' if det in ('det1', 'det3')
                         else f'\n{ch.P2IN_WINDOW.get(det, "")}')
                ax.set_title(f'{det}   (beam station {stn}){extra}', color=col,
                             fontsize=14.5)
    axes[0, 0].set_ylabel('efficiency')
    axes[1, 0].set_ylabel('efficiency')
    fig.suptitle('Bench and beam, chamber by chamber — top row vs mesh, '
                 'bottom row vs drift\n'
                 'P2_MID = det1 and P2_OUT = det3 all campaign; P2_IN held '
                 'det2 (22–23 and 26–27 Jul), det4 (24–26 Jul) and '
                 'det5 = the CERN-built chamber (from 28 Jul)',
                 fontsize=15, y=1.01)
    save(fig, out, 'bench_beam_grid')


# ==================================================================== 2 ==== #
def _eff_table():
    d = pd.read_csv(f'{VMM}/efficiency_table.csv')
    d['cfg'] = d['run'].astype(str) + ' / ' + d['sub'].astype(str)
    return d


def fig_dream_vs_vmm(out):
    """The like-for-like plot: both readouts on the same three detectors,
    the same uRWELL tracks, the same cuts.  This is only legitimate because
    the VMM efficiency was re-derived against that reference (Aug 2026)."""
    dream = {'P2_IN': 0.9649, 'P2_MID': 0.9706, 'P2_OUT': 0.9604}
    d = _eff_table()

    def row(run, sub):
        m = d[(d['run'] == run) & (d['sub'] == sub)]
        return {st: float(m[f'{st}_eff'].iloc[0]) if len(m) else np.nan
                for st in dream}

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(14.6, 5.6), gridspec_kw={'width_ratios': [1.05, 1]})

    # ---- left: the like-for-like comparison, station by station -------- #
    series = [
        ('DREAM (25 Jul, highstat_eff_1)', dream, '#111111', 'o'),
        ('VMM, best config of the campaign\n(run_46: gain 4.5 mV/fC, 200 ns)',
         row('run_46', 'cfg_gain4.5_peaktime200'), '#2ca02c', 's'),
        ('VMM, same config file with the per-chip\nthreshold lines removed '
         '(run_47)',
         row('run_47', 'cfg_gain4.5_peaktime200_deflt'), '#d62728', 'v'),
    ]
    sts = ['P2_IN', 'P2_MID', 'P2_OUT']
    x = np.arange(3)
    for lbl, vals, col, mk in series:
        y = [vals[s] for s in sts]
        ax.plot(x, y, marker=mk, ms=11, lw=2.2, color=col, label=lbl)
        for xi, yi in zip(x, y):
            if np.isfinite(yi):
                ax.annotate(f'{yi:.3f}', (xi, yi), textcoords='offset points',
                            xytext=(0, 13), ha='center', fontsize=10.5,
                            color=col, fontweight='bold')

    ax.axhspan(0.95, 0.975, color='#111111', alpha=0.07, zorder=0)
    ax.text(-0.42, 0.962, 'DREAM band', fontsize=9, color='#555',
            va='center')
    ax.annotate('the best individual VMM pads of P2_OUT\n'
                'already read 0.962 / 0.957 / 0.954',
                xy=(2.0, 0.854), xytext=(0.75, 0.72), fontsize=9.8,
                color='#1a1a1a', ha='center',
                arrowprops=dict(arrowstyle='->', color='#555', lw=1.3),
                bbox=dict(fc='#f4faf4', ec='#8bbf8b',
                          boxstyle='round,pad=0.4'))
    ax.text(0.0, 0.235, 'P2_IN: lower bound only\n(raised thresholds +\n'
            'suspect pad map)', fontsize=8.4, color='#7a3030', ha='center',
            va='bottom')
    ax.set_xticks(x)
    ax.set_xticklabels([f'{s}\nz = {z} mm' for s, z in
                        zip(sts, (320, 630, 940))])
    ax.set_xlim(-0.45, 2.45)
    ax.set_ylim(0.0, 1.10)
    ax.set_ylabel('absolute efficiency (uRWELL-referenced)')
    ax.set_title('(a) same detectors, same tracks, same cuts', fontsize=11.5)
    ax.legend(loc='lower right', framealpha=0.96, fontsize=8.8)

    # ---- right: the VMM configuration ladder --------------------------- #
    ladder = [(25, 'run_39', 'cfg_gain3.0_peaktime25', 3.0),
              (50, 'run_40', 'cfg_gain3.0_peaktime50', 3.0),
              (100, 'run_38', 'cfg_gain3.0_peaktime100', 3.0),
              (200, 'run_41', 'cfg_gain3.0_peaktime200', 3.0),
              (25, 'run_44', 'cfg_gain4.5_peaktime25', 4.5),
              (50, 'run_54', 'cfg_gain4.5_peaktime50', 4.5),
              (100, 'run_45', 'cfg_gain4.5_peaktime100', 4.5),
              (200, 'run_46', 'cfg_gain4.5_peaktime200', 4.5)]
    for gain, ls, mk in ((3.0, '--', 'o'), (4.5, '-', 's')):
        for st in ('P2_MID', 'P2_OUT'):
            xs, ys = [], []
            for pt, run, sub, g in ladder:
                if g != gain:
                    continue
                m = d[(d['run'] == run) & (d['sub'] == sub)]
                if len(m):
                    xs.append(pt)
                    ys.append(float(m[f'{st}_eff'].iloc[0]))
            ax2.plot(xs, ys, marker=mk, ms=8, lw=2.0, ls=ls, color=C[st],
                     label=f'{st}, gain {gain} mV/fC')
    for st in ('P2_MID', 'P2_OUT'):
        ax2.axhline(dream[st], color=C[st], ls=':', lw=1.6, alpha=0.8)
        dy = 0.022 if st == 'P2_MID' else -0.030
        ax2.text(24, dream[st] + dy, f'DREAM {st} = {dream[st]:.3f}',
                 fontsize=9, color=C[st],
                 va='bottom' if dy > 0 else 'top', ha='left')
    ax2.set_xscale('log')
    ax2.set_xticks([25, 50, 100, 200])
    ax2.set_xticklabels(['25', '50', '100', '200'])
    ax2.set_xlim(22, 235)
    ax2.set_ylim(0, 1.05)
    ax2.set_xlabel('VMM3a peaking time [ns]')
    ax2.set_ylabel('absolute efficiency (uRWELL-referenced)')
    ax2.set_title('(b) every electronics knob still moves the efficiency\n'
                  '— the readout is gain-starved, the chamber is not',
                  fontsize=11.5)
    ax2.legend(loc='lower right', framealpha=0.96, fontsize=9)

    fig.suptitle('DREAM vs VMM3a on the same three chambers — the deficit is '
                 'the discriminator threshold, not the detector', fontsize=13.5,
                 y=1.02)
    fig.tight_layout()
    save(fig, out, 'dream_vs_vmm')


# ==================================================================== 3 ==== #
def fig_vmm_threshold(out):
    """Three independent handles that all say 'threshold', on one page."""
    d = _eff_table()
    fig, axes = plt.subplots(1, 3, figsize=(15.4, 5.0))

    # (a) the mesh scan and the Landau-survival model -------------------- #
    ax = axes[0]
    with open(f'{VMM}/model_P2_OUT.json') as fh:
        m = json.load(fh)
    pts = pd.DataFrame(m['points'])
    ax.plot(pts['dv'], pts['eff'], 'o', ms=7, color=C['P2_OUT'],
            label='P2_OUT, VMM mesh scan (run_32/33)')
    ax.plot(pts['dv'], pts['model'], '-', lw=2, color='#333',
            label=(r'Landau swept past a fixed threshold,'
                   '\n' r'fitted $V_0$ = %.1f V, thr = %.2f$\times$MPV'
                   % (m['V0_volts'], m['r0_over_mpv'])))
    g = m['amplifier_gain_step']
    ax.plot([g['equivalent_mesh_volts']], [g['observed_eff']], '*', ms=20,
            color='#d62728', zorder=5,
            label=(r'$\times$1.5 amplifier gain (run_41$\to$48):'
                   '\nmeasured %.3f, model %.3f'
                   % (g['observed_eff'], g['predicted_eff'])))
    ax.set_xlabel('mesh voltage relative to working point [V]')
    ax.set_ylabel('VMM efficiency, P2_OUT')
    ax.set_yscale('log')
    ax.set_title('(a) one curve fits the HV scan AND the\n'
                 r'amplifier-gain step — $V_0$ = 22.4 V is the'
                 '\nbulk-Micromegas value, not a fit input', fontsize=11)
    ax.legend(fontsize=8.6, loc='upper left')

    # (b) the threshold DAC alone --------------------------------------- #
    ax = axes[1]
    pairs = [('run_47', 'cfg_gain4.5_peaktime200_deflt',
              'run_48', 'cfg_gain4.5_peaktime200_opt', '4.5 / 200, 1 Aug'),
             ('run_42', 'cfg_gain3.0_peaktime200_deflt',
              'run_43', 'cfg_gain3.0_peaktime200_opt', '3.0 / 200, 1 Aug'),
             ('run_51', 'cfg_gain3.0_peaktime200_deflt',
              'run_50', 'cfg_gain3.0_peaktime200_opt', '3.0 / 200, later'),
             ('run_67', 'cfg_gain4.5_peaktime200_deflt',
              'run_66', 'cfg_gain4.5_peaktime200_opt', '4.5 / 200, gas B')]
    w, off = 0.19, -0.30
    for i, (rd, sd, ro, so, lbl) in enumerate(pairs):
        a = d[(d['run'] == rd) & (d['sub'] == sd)]
        b = d[(d['run'] == ro) & (d['sub'] == so)]
        if not len(a) or not len(b):
            continue
        for j, st in enumerate(('P2_MID', 'P2_OUT')):
            lo, hi = float(a[f'{st}_eff'].iloc[0]), float(b[f'{st}_eff'].iloc[0])
            xp = j + off + i * w
            ax.vlines(xp, lo, hi, color=C[st], lw=2.4, alpha=0.8)
            ax.plot([xp], [lo], 'v', ms=8, color='#d62728', zorder=4,
                    label='base thresholds' if (i == 0 and j == 0) else None)
            ax.plot([xp], [hi], '^', ms=8, color='#2ca02c', zorder=4,
                    label='tuned per chip' if (i == 0 and j == 0) else None)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['P2_MID', 'P2_OUT'])
    ax.set_ylabel('VMM efficiency')
    ax.set_ylim(0, 1.0)
    ax.set_title('(b) the SAME config file, per-chip threshold\n'
                 'lines set or commented out.\n'
                 'Four pairs, 36 min apart: nothing else changed',
                 fontsize=11)
    ax.legend(loc='lower left', fontsize=9, framealpha=0.95)

    # (c) per-pad pulse height vs efficiency ---------------------------- #
    ax = axes[2]
    q = pd.DataFrame({
        'label': ['90–127', '127–139', '139–168', '168–202', '202–271'],
        'eff': [0.752, 0.816, 0.892, 0.920, 0.921],
        'n': [151914, 241640, 163012, 174252, 106159]})
    ax.bar(q['label'], q['eff'], color=C['P2_OUT'], alpha=0.85, width=0.68)
    ax.axhline(0.9604, color='#111', ls='--', lw=1.8)
    ax.text(0.05, 0.972, 'DREAM, same chamber: 0.960', va='bottom',
            ha='left', fontsize=9.5, color='#111')
    for i, (e, n) in enumerate(zip(q['eff'], q['n'])):
        ax.text(i, e + 0.012, f'{e:.3f}', ha='center', fontsize=10,
                fontweight='bold')
    ax.set_xlabel('pad median pulse height [ADC], quintiles')
    ax.set_ylabel('VMM efficiency of those pads')
    ax.set_ylim(0, 1.06)
    ax.set_title('(c) pad by pad, efficiency tracks PULSE HEIGHT\n'
                 r'($r$ = +0.55 over 69 illuminated pads).'
                 '\nA selection effect would run the other way', fontsize=11)

    fig.suptitle('Why the VMM reads 85 % where DREAM reads 96 % — three '
                 'independent handles, all pointing at the discriminator',
                 fontsize=13.5, y=1.005)
    fig.tight_layout()
    save(fig, out, 'vmm_threshold')


# ==================================================================== 4 ==== #
def fig_snr_matrix(out):
    """Nov-2025 VMM3a shaping scan, frozen for the talk.  The scan is NOT a
    full grid -- only 6 of the 9 (gain, peaking) cells were taken -- so it is
    drawn as a grid with the untaken cells left blank, never as lines."""
    d = pd.read_csv(f'{NOV}/vmm_snr_results.csv')
    geo = {8: 'P2 Small det 3', 9: 'P2 Small det 3', 10: 'P2 Small det 1',
           11: 'P2 Small det 1', 12: 'P2 Large det', 13: 'P2 Large det',
           14: 'P2 Large det', 15: 'P2 Large det'}
    d['geo'] = d['vmm_id'].map(geo)
    d = d[d['noise_quality'] == 'ok']

    gains = [3.0, 4.5, 6.0]
    pts = [25, 50, 100, 200]
    geos = ['P2 Large det', 'P2 Small det 1', 'P2 Small det 3']

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.5), sharey=True)
    vmin, vmax = 12, 34
    for ax, g in zip(axes, geos):
        grid = np.full((len(gains), len(pts)), np.nan)
        for i, sg in enumerate(gains):
            for j, pt in enumerate(pts):
                v = d[(d['geo'] == g) & (d['sg'] == sg) & (d['snt'] == pt)]
                if len(v):
                    grid[i, j] = v['snr'].mean()
        im = ax.imshow(grid, cmap='viridis', vmin=vmin, vmax=vmax,
                       aspect='auto', origin='lower')
        for i in range(len(gains)):
            for j in range(len(pts)):
                if np.isnan(grid[i, j]):
                    ax.text(j, i, 'not\ntaken', ha='center', va='center',
                            fontsize=8.5, color='#999')
                else:
                    best = grid[i, j] == np.nanmax(grid)
                    ax.text(j, i, f'{grid[i, j]:.1f}', ha='center',
                            va='center', fontsize=13,
                            fontweight='bold' if best else 'normal',
                            color='white' if grid[i, j] < 26 else 'black')
                    if best:
                        ax.add_patch(plt.Rectangle(
                            (j - 0.5, i - 0.5), 1, 1, fill=False,
                            ec='#d62728', lw=3.0))
        ax.set_xticks(range(len(pts)))
        ax.set_xticklabels(pts)
        ax.set_yticks(range(len(gains)))
        ax.set_yticklabels(gains)
        ax.set_xlabel('peaking time [ns]')
        ax.set_title(g, fontsize=12)
        ax.grid(False)
    axes[0].set_ylabel('amplifier gain [mV/fC]')
    fig.colorbar(im, ax=axes, fraction=0.022, pad=0.015,
                 label='SNR = signal MPV / noise σ (MAD)')
    fig.suptitle('VMM3a shaping scan — SPS Nov 2025, two pad geometries, '
                 '~5 kHz muons\n'
                 '100 ns peaking wins on all three geometries (red box); the '
                 'gain choice is within the VMM-to-VMM spread.\n'
                 'Plurality vote over all 8 VMMs, VMM and channel level '
                 'agreeing: gain 3.0 mV/fC / 100 ns',
                 fontsize=12.5, y=1.10)
    save(fig, out, 'snr_matrix')


# ==================================================================== 5 ==== #
def fig_timing_campaigns(out):
    """The timing narrative: three campaigns, one physics floor, one goal."""
    items = [
        ('Garfield++ / Magboltz\nphysics floor\n(Ar/iC$_4$H$_{10}$ 95/5)',
         5.0, 3.0, 7.0, '#7f7f7f'),
        ('Cosmic bench\nbest conditions\n(DREAM waveforms)',
         28.8, 28.8, 32.0, '#1f77b4'),
        ('SPS beam\nP2_MID', 15.5, None, None, '#ff7f0e'),
        ('SPS beam\nP2_OUT', 18.4, None, None, '#2ca02c'),
        ('SPS beam\nP2_IN', 22.4, None, None, '#d62728'),
        ('Next campaign,\nAr/CF$_4$/iC$_4$H$_{10}$ 88/10/2\n(model)',
         13.5, 13.0, 14.0, '#9467bd'),
    ]
    fig, ax = plt.subplots(figsize=(10.4, 5.6))
    y = np.arange(len(items))[::-1]
    for yi, (lbl, v, lo, hi, col) in zip(y, items):
        if lo is not None and hi is not None and hi > lo:
            ax.barh(yi, hi - lo, left=lo, height=0.46, color=col, alpha=0.35)
        ax.plot([v], [yi], 'o', ms=13, color=col, zorder=5)
        ax.text(37.5, yi, f'{v:.1f} ns', va='center', ha='right',
                fontsize=12, fontweight='bold', color=col)

    ax.axvline(20, color='#111', ls='--', lw=2)
    ax.text(19.4, len(items) - 0.4, 'P2 goal: 20 ns', rotation=90,
            va='top', ha='right', fontsize=11, fontweight='bold')
    ax.set_yticks(y)
    ax.set_yticklabels([i[0] for i in items], fontsize=10.5)
    ax.set_xlabel('time resolution σ [ns]')
    ax.set_xlim(0, 38)
    ax.set_title('Time resolution across the programme — the 20 ns goal is met '
                 'at two of three stations\n'
                 'the bench was drift-geometry limited; the beam is '
                 'gas + walk limited; the new gas removes the gas term',
                 fontsize=12)
    ax.grid(axis='y', alpha=0)
    save(fig, out, 'timing_campaigns')


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--out', default='figs')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    print(f'writing to {a.out}')
    for f in (fig_bench_beam_mesh, fig_bench_beam_drift,
              fig_bench_beam_all, fig_bench_beam_maps,        # fig_fe55_bench_beam: parked, see above
              fig_dream_vs_vmm, fig_vmm_threshold,
              fig_snr_matrix, fig_timing_campaigns):
        try:
            f(a.out)
        except Exception as exc:                      # keep going, report
            print(f'  !! {f.__name__}: {type(exc).__name__}: {exc}')


if __name__ == '__main__':
    main()
