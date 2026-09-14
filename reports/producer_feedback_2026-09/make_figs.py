#!/usr/bin/env python3
"""Figures for the two chamber-feedback decks (2026-09): one for the Saclay
bulk production (det1-det4), one for the CERN MPT workshop (det5).

The decks are deliberately separate -- no Saclay/CERN comparison anywhere --
and DREAM-only (no VMM). Minimal style: no titles, the slide text carries the
message.

Two efficiency definitions appear and are never mixed on one curve:
  * beam, uRWELL-referenced: tracks from the two EIC uRWELL planes
    (urw_reference/urw_p2_efficiency.py)
  * beam, P2 self-tracking ("tag-and-probe"): the other two P2 stations tag
    the track (stage 22). This is how P2 itself will run, without an external
    reference; it reads a few points LOW because the tag is built from 12 mm
    pads.
  * bench: M3-telescope cosmics, whole active area, reconstructed within
    R = 40 mm, spark-vetoed.

Bench gas Ar/iC4H10 95/5, beam gas Ar/CO2/iC4H10 93/5/2 -- the working
voltages are not comparable between the two.

Usage:
  python3 make_figs.py --saclay-out <dir> --cern-out <dir>
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HOME = os.path.expanduser('~')
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'sps_beam_analysis'))
import chamber_history as ch  # noqa: E402

WS = os.path.join(HOME, 'Documents', 'PostDocSaclay', 'data', 'SPS_Beam_Test',
                  'mpgd26_workspace')
URW = os.path.join(WS, 'products', 'urw_local', 'urw_referenced_efficiency')
BENCH = os.path.join(HOME, 'Documents', 'PostDocSaclay', 'data', 'Cosmic_Bench',
                     'Analysis')

COL = {'det1': '#1f6fb4', 'det2': '#2e9b48', 'det3': '#c8322f',
       'det4': '#e3802a', 'det5': '#4a3aa7'}

plt.rcParams.update({
    'font.size': 15, 'axes.labelsize': 16, 'legend.fontsize': 13,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': True, 'grid.alpha': 0.25, 'legend.frameon': False,
    'lines.linewidth': 2.2, 'lines.markersize': 7,
    'savefig.dpi': 170, 'savefig.bbox': 'tight',
    'figure.facecolor': 'white', 'axes.facecolor': 'white'})


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, f'{name}.png'))
    plt.close(fig)
    print('  wrote', os.path.join(out, f'{name}.png'))


def pct_axis(ax, lo=0, hi=100):
    ax.set_ylim(lo, hi)
    ax.set_ylabel('efficiency [%]')


# ----------------------------------------------------------------- beam data --
def urw(run):
    f = os.path.join(URW, run, f'urw_p2_efficiency_{run}.csv')
    return pd.read_csv(f) if os.path.isfile(f) else None


def urw_by_hv(d, station):
    """uRWELL efficiency per (mesh, drift), repeated settings combined by counts."""
    s = d[d.station == station].copy()
    s['k'] = s.eff * s.n
    g = s.groupby(['mesh_hv', 'drift_hv']).agg(k=('k', 'sum'), n=('n', 'sum'))
    g['eff'] = 100 * g.k / g.n
    return g.reset_index().rename(columns={'mesh_hv': 'mesh', 'drift_hv': 'drift'})


_TP = None


def tag_probe(run, station):
    """P2 self-tracking efficiency per (mesh, drift) for one station and run."""
    global _TP
    if _TP is None:
        t = pd.read_csv(os.path.join(WS, 'report_data', 'dream_tag_probe.csv'))
        t = t[t.vetoed == True]  # noqa: E712
        hv = pd.read_csv(os.path.join(WS, 'eos_inventory', 'hv_setpoints.csv'))
        hv = hv.rename(columns={'mesh_or_resist': 'mesh'})
        _TP = t.merge(hv, on=['run', 'sub_run', 'det'], how='left')
    s = _TP[(_TP.run == run) & (_TP.det == station)].copy()
    s['k'] = s.eff_corr * s.n_tag
    g = s.groupby(['mesh', 'drift']).agg(k=('k', 'sum'), n=('n_tag', 'sum'))
    g['eff'] = 100 * g.k / g.n
    return g.reset_index()


def gap(d, dv):
    return d[(d.drift - d.mesh).round() == dv].sort_values('mesh')


# ---------------------------------------------------------------- bench data --
def bench_csv(pattern):
    fs = sorted(glob.glob(os.path.join(BENCH, pattern)))
    if not fs:
        raise FileNotFoundError(pattern)
    return pd.read_csv(fs[0])


def eff_map(ax, npz, vmin=80):
    z = np.load(npz, allow_pickle=True)
    e = np.array(z['eff_any'], dtype=float) * 100
    e[np.array(z['counts']) < float(z['min_rays'])] = np.nan
    x0, x1, y0, y1 = z['extent']
    im = ax.imshow(e, origin='lower', extent=(x0, x1, y0, y1), cmap='viridis',
                   vmin=vmin, vmax=100, interpolation='nearest')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_aspect('equal')
    ax.grid(False)
    return im


# ================================================================== SACLAY ===
def saclay(out):
    print('Saclay deck ->', out)
    # --- beam: efficiency vs mesh --------------------------------------------
    dms = urw('drift_mesh_scan_1')
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    for det, st, dv in [('det1', 'P2_MID', 250), ('det3', 'P2_OUT', 250),
                        ('det4', 'P2_IN', 200)]:
        u = gap(urw_by_hv(dms, st), dv)
        ax.plot(u.mesh, u.eff, 'o-', color=COL[det], label=f'{det}')
        t = gap(tag_probe('drift_mesh_scan_1', st), dv)
        ax.plot(t.mesh, t.eff, 'o--', color=COL[det], mfc='white', alpha=.8)
    dm2 = urw('drift_mesh_2d_2')
    if dm2 is not None:
        u = gap(urw_by_hv(dm2, 'P2_IN'), 300)
        ax.plot(u.mesh, u.eff, 'o-', color=COL['det2'], label='det2 (repaired)')
    t = gap(tag_probe('drift_mesh_2d_2', 'P2_IN'), 300)
    ax.plot(t.mesh, t.eff, 'o--', color=COL['det2'], mfc='white', alpha=.8,
            label=None if dm2 is not None else 'det2 (repaired)')
    ax.plot([], [], 'k-', label='uRWELL tracks')
    ax.plot([], [], 'k--', label='P2 self-tracking')
    ax.set_xlabel('mesh voltage [V]')
    pct_axis(ax)
    ax.legend(loc='lower right', ncol=2)
    save(fig, out, 'beam_mesh')

    # --- beam: efficiency vs drift at mesh 450 V -----------------------------
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    for det, st in [('det1', 'P2_MID'), ('det3', 'P2_OUT')]:
        u = urw_by_hv(dms, st)
        u = u[(u.mesh == 450) & (u.drift > 450)].sort_values('drift')
        ax.plot(u.drift, u.eff, 'o-', color=COL[det], label=det)
        t = tag_probe('drift_mesh_scan_1', st)
        t = t[(t.mesh == 450) & (t.drift > 450)].sort_values('drift')
        ax.plot(t.drift, t.eff, 'o--', color=COL[det], mfc='white', alpha=.8)
    if dm2 is not None:
        u = urw_by_hv(dm2, 'P2_IN')
        u = u[u.mesh == 450].sort_values('drift')
        ax.plot(u.drift, u.eff, 'o-', color=COL['det2'], label='det2 (repaired)')
    t = tag_probe('drift_mesh_2d_2', 'P2_IN')
    t = t[t.mesh == 450].sort_values('drift')
    ax.plot(t.drift, t.eff, 'o--', color=COL['det2'], mfc='white', alpha=.8,
            label=None if dm2 is not None else 'det2 (repaired)')
    ax.plot([], [], 'k-', label='uRWELL tracks')
    ax.plot([], [], 'k--', label='P2 self-tracking')
    ax.set_xlabel('drift voltage [V]   (mesh 450 V)')
    pct_axis(ax, 70, 100)
    ax.legend(loc='lower right', ncol=2)
    save(fig, out, 'beam_drift')

    # --- bench: efficiency vs mesh --------------------------------------------
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    d1 = bench_csv('det1/p2_det1_long_run_mesh_scan_7-19-26/mesh_scan/'
                   '16_mesh_scan_efficiency/efficiency_vs_mesh_without_connectors_1_2_10_spark_vetoed.csv')
    ax.plot(d1.mesh, 100 * d1.eff_reco, 's-', color=COL['det1'],
            label='det1  (drift = mesh + 310 V)')
    d2 = bench_csv('det2/p2_det1_det2_long_run_mesh_scan_7-9-26/hv_scan/'
                   '11_hv_scan_efficiency/efficiency_vs_hv_without_connectors_1_8_9_10_spark_vetoed.csv')
    ax.plot(d2.hv, 100 * d2.eff_reco, 's-', color=COL['det2'],
            label='det2  (drift = mesh + 170 V)')
    ax.set_xlabel('mesh voltage [V]')
    pct_axis(ax, 60, 100)
    ax.legend(loc='lower right')
    save(fig, out, 'bench_mesh')

    # --- bench: efficiency vs drift -------------------------------------------
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    d1 = bench_csv('det1/p2_det1_drift_scan_7-19-26/drift_scan/16_drift_scan_efficiency/'
                   'efficiency_vs_drift_without_connectors_1_2_10_spark_vetoed.csv')
    ax.plot(d1.drift - d1.mesh, 100 * d1.eff_reco, 's-', color=COL['det1'],
            label='det1  (mesh 415 V)')
    d3 = bench_csv('det3/p2_det3_det4_drift_scan_7-16-26/drift_scan/16_drift_scan_efficiency/'
                   'efficiency_vs_drift_without_connectors_1_8_9_10_spark_vetoed.csv')
    ax.plot(d3.drift - d3.mesh, 100 * d3.eff_reco, 's-', color=COL['det3'],
            label='det3  (mesh 420 V)')
    ax.set_xlabel('drift − mesh voltage [V]')
    pct_axis(ax)
    ax.legend(loc='lower right')
    save(fig, out, 'bench_drift')

    # --- bench: efficiency vs time (det4 charges up, det3 loses its drift HV) -
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    d4 = bench_csv('det4/p2_det4_long_run_7-20-26/long_run_det4_410_610/20_charging_up/'
                   'charging_vs_time_without_connectors_1_2_7_10_spark_vetoed.csv')
    ax.errorbar(d4.h, d4.eff, d4.eff_err, fmt='s-', color=COL['det4'],
                label='det4  (mesh 410 / drift 610 V)')
    d3 = bench_csv('det3/p2_det3_mesh_scan_det4_initial_7-16-26/initial_run_det3_420_820_det4_430_830/'
                   '20_charging_up/charging_vs_time_without_connectors_1_8_9_10_spark_vetoed.csv')
    ax.errorbar(d3.h, d3.eff, d3.eff_err, fmt='s-', color=COL['det3'],
                label='det3  (mesh 420 / drift 820 V)')
    ax.set_xlabel('time since HV on [h]')
    pct_axis(ax, 60, 100)
    ax.legend(loc='lower right')
    save(fig, out, 'bench_time')

    # --- bench: det1 efficiency map -------------------------------------------
    f = glob.glob(os.path.join(BENCH, 'det1/p2_det1_long_run_efficiency_7-19-26/long_run_det1_415_615/'
                                      '06_efficiency/efficiency_map_sliding_*.npz'))[0]
    fig, ax = plt.subplots(figsize=(7.4, 6.6))
    im = eff_map(ax, f)
    fig.colorbar(im, ax=ax, shrink=.85, label='efficiency (any pad fired) [%]')
    save(fig, out, 'bench_map_det1')


# ==================================================================== CERN ===
def cern(out):
    print('CERN deck ->', out)
    det = 'det5'
    c = COL[det]
    # --- beam: efficiency vs mesh --------------------------------------------
    hv2 = urw('p2in_hvrange_2')
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    u = gap(urw_by_hv(hv2, 'P2_IN'), 300)
    ax.plot(u.mesh, u.eff, 'o-', color=c, label='uRWELL tracks')
    t = gap(tag_probe('p2in_hvrange_2', 'P2_IN'), 300)
    ax.plot(t.mesh, t.eff, 'o--', color=c, mfc='white', label='P2 self-tracking')
    ax.set_xlabel('mesh voltage [V]   (drift = mesh + 300 V)')
    pct_axis(ax)
    ax.legend(loc='lower right')
    save(fig, out, 'beam_mesh')

    # --- beam: efficiency vs drift (self-tracking only: no uRWELL for this run)
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    t = tag_probe('p2_mesh_drift_2d_1', 'P2_IN')
    for mesh, a in [(450, 1.0), (440, .75), (430, .5)]:
        s = t[t.mesh == mesh].sort_values('drift')
        ax.plot(s.drift, s.eff, 'o--', color=c, mfc='white', alpha=a,
                label=f'mesh {mesh} V')
    ax.set_xlabel('drift voltage [V]')
    pct_axis(ax, 70, 100)
    ax.legend(loc='lower right', title='P2 self-tracking')
    save(fig, out, 'beam_drift')

    base = os.path.join(BENCH, 'det5')
    sfx = '_without_connectors_9_10_spark_vetoed'

    def scan(kind):
        d = os.path.join(base, 'p2_det5_alignment_mesh_drift_scan_9-8-26',
                         f'{kind}_scan', f'16_{kind}_scan_efficiency')
        f40 = os.path.join(d, f'efficiency_vs_{kind}{sfx}_r40.csv')
        return pd.read_csv(f40 if os.path.isfile(f40)
                           else os.path.join(d, f'efficiency_vs_{kind}{sfx}.csv'))

    # --- bench: mesh ------------------------------------------------------------
    m = scan('mesh')
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    ax.plot(m.mesh, 100 * m.eff_reco, 's-', color=c, label='reconstructed')
    ax.plot(m.mesh, 100 * m.eff_anyhit, 's:', color=c, mfc='white',
            label='any pad fired')
    ax.set_xlabel('mesh voltage [V]   (drift = mesh + 250 V)')
    pct_axis(ax, 40, 100)
    ax.legend(loc='lower right')
    save(fig, out, 'bench_mesh')

    # --- bench: drift -----------------------------------------------------------
    d = scan('drift')
    d = d[d.drift > d.mesh]
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    ax.plot(d.drift - d.mesh, 100 * d.eff_reco, 's-', color=c, label='reconstructed')
    ax.plot(d.drift - d.mesh, 100 * d.eff_anyhit, 's:', color=c, mfc='white',
            label='any pad fired')
    ax.set_xlabel('drift − mesh voltage [V]   (mesh 420 V)')
    pct_axis(ax, 70, 100)
    ax.legend(loc='lower right')
    save(fig, out, 'bench_drift')

    # --- bench: 12 h stability ---------------------------------------------------
    L = os.path.join(base, 'p2_det5_long_run_9-9-26', 'long_run_det5_420_820')
    e = pd.read_csv(os.path.join(L, '06_efficiency',
                                 f'efficiency_vs_time{sfx}_r40.csv'))
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.plot(e.t_h, 100 * e.eff_within, 's-', color=c, label='reconstructed')
    ax.plot(e.t_h, 100 * e.eff_any, 's:', color=c, mfc='white', label='any pad fired')
    ax.set_xlabel('time [h]   (mesh 420 / drift 820 V)')
    pct_axis(ax, 70, 100)
    ax.legend(loc='lower right')
    save(fig, out, 'bench_time')

    # --- bench: efficiency map ----------------------------------------------------
    f = glob.glob(os.path.join(L, '06_efficiency', 'efficiency_map_sliding_*.npz'))[0]
    fig, ax = plt.subplots(figsize=(7.4, 6.6))
    im = eff_map(ax, f)
    fig.colorbar(im, ax=ax, shrink=.85, label='efficiency (any pad fired) [%]')
    save(fig, out, 'bench_map')

    # --- bench: drift-time surface (possible bow) ---------------------------------
    ts = pd.read_csv(os.path.join(L, '19_timing_surface', f'timing_surface{sfx}.csv'))
    ts = ts.dropna(subset=['median_ns'])
    ts = ts[ts.n >= 20]
    off = ts.median_ns - ts.median_ns.median()
    # per-connector constant + one quadratic surface, fitted together; plot the
    # measurement with only the connector constants removed (stage 23 method)
    x = (ts.pad_cx - ts.pad_cx.mean()) / 100
    y = (ts.pad_cy - ts.pad_cy.mean()) / 100
    conns = sorted(ts.connector_N.unique())
    C = np.column_stack([(ts.connector_N == k).astype(float) for k in conns])
    Q = np.column_stack([x, y, x * x, x * y, y * y])
    coef, *_ = np.linalg.lstsq(np.hstack([C, Q]), off.values, rcond=None)
    corr = off.values - C @ coef[:len(conns)]
    corr -= np.median(corr)
    b = 36.0
    ts = ts.assign(v=corr, bx=(ts.pad_cx // b) * b + b / 2, by=(ts.pad_cy // b) * b + b / 2)
    g = ts.groupby(['bx', 'by']).agg(v=('v', 'median'), n=('v', 'size')).reset_index()
    g = g[g.n >= 4]
    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    lim = np.nanpercentile(np.abs(g.v), 95)
    sc = ax.scatter(g.bx, g.by, c=g.v, cmap='coolwarm', vmin=-lim, vmax=lim,
                    marker='s', s=330, edgecolors='none')
    fig.colorbar(sc, ax=ax, shrink=.85, label='arrival time offset [ns]')
    ax.set_xlabel('board x [mm]')
    ax.set_ylabel('board y [mm]')
    ax.set_aspect('equal')
    ax.grid(False)
    save(fig, out, 'bench_timing_surface')

    v = pd.read_csv(os.path.join(L, '23_drift_gap_bow', f'bow_vd_scaling{sfx}.csv'))
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.plot(v.drift, v.dz_um / 1000, 's-', color=c)
    ax.axhline(np.mean(v.dz_um) / 1000, color='0.4', ls='--', lw=1.2)
    ax.set_xlabel('drift voltage [V]   (mesh 420 V)')
    ax.set_ylabel('implied gap variation [mm]')
    ax.set_ylim(0, 1.6)
    save(fig, out, 'bench_bow_scaling')


# ===================================================== shared / talk figures ===
class _NoTitles:
    """Silence set_title / suptitle, so the MPGD26 talk producers (whose titles
    carry the talk's narrative) draw bare panels for the decks."""
    def __enter__(self):
        from matplotlib.axes import Axes
        from matplotlib.figure import Figure
        self._a, self._f = Axes.set_title, Figure.suptitle
        Axes.set_title = lambda self, *k, **kw: None
        Figure.suptitle = lambda self, *k, **kw: None

    def __exit__(self, *exc):
        from matplotlib.axes import Axes
        from matplotlib.figure import Figure
        Axes.set_title, Figure.suptitle = self._a, self._f


def setup_figs(out):
    """Test-setup pictures: the H4 three-station render and the bench schematic."""
    import shutil
    os.makedirs(out, exist_ok=True)
    shutil.copy(os.path.join(REPO, 'conference', 'figures', '4_act3_beam',
                             's19_coincidence_hero.png'),
                os.path.join(out, 'setup_beam.png'))
    sys.path.insert(0, os.path.join(REPO, 'conference'))
    import make_bench_schematic as mbs
    with plt.rc_context({'axes.grid': False}), _NoTitles():
        fig, ax = plt.subplots(figsize=(11.0, 8.6))
        mbs.draw(ax)
        fig.savefig(os.path.join(out, 'setup_bench.png'), dpi=150,
                    bbox_inches='tight', pad_inches=0.2, facecolor='white')
        plt.close(fig)
    print('  wrote setup_beam.png, setup_bench.png')


def _timing_curve(tm, run, station):
    """(drift V, sigma ns, mesh V) for one station in one DREAM drift run;
    repeated setpoints collapsed to their median (fig_drift_timing_all3 rule)."""
    g = tm[(tm.run == run) & (tm.axis == 'drift')]
    d = pd.DataFrame({'dv': g[f'drift_v_{station}'], 'mv': g[f'mesh_v_{station}'],
                      'sg': g[f'{station}_sigma']}).astype(float).dropna()
    d = d[d.dv > d.mv]
    m = d.groupby('dv').sg.median().sort_index()
    return m.index.to_numpy(), m.to_numpy(), (d.mv.iloc[0] if len(d) else np.nan)


def timing_fig(out, series):
    """DREAM single-station time resolution vs drift voltage.
    series: [(run, station, det, legend label, style)]."""
    tm = pd.read_csv(os.path.join(WS, 'report_data', 'dream_timing_scans.csv'))
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    for run, station, det, lab, style in series:
        x, y, mv = _timing_curve(tm, run, station)
        if len(x) < 2:
            continue
        ax.plot(x, y, style, color=COL[det], mfc='white' if '--' in style else None,
                label=f'{lab}, mesh {mv:.0f} V')
        print(f'    timing {det} {run}: sigma {y[-1]:.1f} ns at drift {x[-1]:.0f} V '
              f'(min {y.min():.1f} ns)')
    ax.set_yscale('log')
    ticks = [15, 20, 30, 50, 100, 200]
    ax.set_yticks(ticks)
    ax.set_yticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.axhline(20, color='#c8322f', lw=1.2, ls=':')
    ax.text(0.02, 20, 'P2 goal 20 ns', transform=ax.get_yaxis_transform(),
            va='bottom', color='#c8322f', fontsize=12)
    ax.set_xlabel('drift voltage [V]')
    ax.set_ylabel('time resolution σ [ns]')
    ax.legend(loc='upper right')
    save(fig, out, 'beam_timing')


def bench_beam_field(out):
    """The same two chambers on the bench and in the beam vs drift FIELD
    (4 mm gap), which lines up the different mesh settings of the two setups."""
    GAP_CM = 0.4
    dms = urw('drift_mesh_scan_1')
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    for det, st, pat in [
            ('det1', 'P2_MID', 'det1/p2_det1_drift_scan_7-19-26/drift_scan/16_drift_scan_efficiency/'
                               'efficiency_vs_drift_without_connectors_1_2_10_spark_vetoed.csv'),
            ('det3', 'P2_OUT', 'det3/p2_det3_det4_drift_scan_7-16-26/drift_scan/16_drift_scan_efficiency/'
                               'efficiency_vs_drift_without_connectors_1_8_9_10_spark_vetoed.csv')]:
        b = bench_csv(pat).sort_values('drift')
        ax.plot((b.drift - b.mesh) / GAP_CM, 100 * b.eff_reco, 's--', color=COL[det],
                mfc='white', label=f'{det} bench (mesh {int(b.mesh.iloc[0])} V)')
        u = urw_by_hv(dms, st)
        u = u[(u.mesh == 450) & (u.drift >= 450)].sort_values('drift')
        ax.plot((u.drift - u.mesh) / GAP_CM, u.eff, 'o-', color=COL[det],
                label=f'{det} beam (mesh 450 V)')
    ax.set_xlabel('drift field [V/cm]   (4 mm drift gap)')
    pct_axis(ax)
    ax.legend(loc='lower right', ncol=2)
    save(fig, out, 'bench_beam_field')


def saclay_extra(out):
    print('Saclay deck, talk figures ->', out)
    setup_figs(out)
    bench_beam_field(out)
    timing_fig(out, [
        ('mesh_drift_scan_up_1', 'P2_MID', 'det1', 'det1, 26 Jul', 'o-'),
        ('p2_mesh_drift_eff_1', 'P2_MID', 'det1', 'det1, 28 Jul', '^--'),
        ('mesh_drift_scan_up_1', 'P2_OUT', 'det3', 'det3, 26 Jul', 'o-'),
        ('p2_mesh_drift_eff_1', 'P2_OUT', 'det3', 'det3, 28 Jul', '^--'),
        # chamber_history flags this run's P2_IN as det4 by the clock, UNCONFIRMED
        ('mesh_drift_scan_up_1', 'P2_IN', 'det4', 'det4*, 26 Jul', 'o-')])
    # the talk's all-chamber overlays and the bench/beam area maps, Saclay
    # chambers only (the CERN chamber is its own deck), talk titles removed
    with plt.rc_context():
        sys.path.insert(0, os.path.join(REPO, 'mpgd2026'))
        import make_talk_figs as mtf
        mtf.CHAMBERS = ['det1', 'det2', 'det3', 'det4']
        # the det3 note ("scan stopped at 420 V") describes its MESH scan, which
        # the overlay drops; on the drift overlay it would mislabel the drift scan
        mtf.BENCH_NOTE.pop('det3', None)
        with _NoTitles():
            mtf.fig_bench_beam_all(out)
            mtf.fig_bench_beam_maps(out)
    for f in glob.glob(os.path.join(out, '*.pdf')) + [os.path.join(out, 'bench_beam_grid.png')]:
        if os.path.exists(f):
            os.remove(f)


def cern_extra(out, saclay_out=None):
    import shutil
    print('CERN deck, talk figures ->', out)
    setup_figs(out)
    timing_fig(out, [('p2_mesh_drift_eff_1', 'P2_IN', 'det5', 'det5, 28 Jul', '^-')])
    # the measurement-area illustration is drawn once (det1 at P2_MID) and shared
    src = os.path.join(saclay_out, 'bench_beam_maps.png') if saclay_out else None
    if src and os.path.isfile(src):
        shutil.copy(src, os.path.join(out, 'bench_beam_maps.png'))
        print('  copied bench_beam_maps.png')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--saclay-out')
    ap.add_argument('--cern-out')
    a = ap.parse_args()
    if a.saclay_out:
        saclay(a.saclay_out)
        saclay_extra(a.saclay_out)
    if a.cern_out:
        cern(a.cern_out)
        cern_extra(a.cern_out, a.saclay_out)


if __name__ == '__main__':
    main()
