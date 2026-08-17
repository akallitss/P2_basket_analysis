#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
22_conference_figures.py

Conference plots, ONE PER FILE so they can be dropped into slides individually.

Slide-copy conventions (2026-08-17): no titles -- the speaker carries the
conditions, and a title repeated on every slide is wasted space.  Fonts are
sized to be readable from the back of a room.  Only det1 is shown: it is the
only detector with the full efficiency + mesh-scan + drift-scan set, and a
comparison plot missing three of four measurements says nothing.

  01_det1_efficiency_map        uniformity, pillar shadows, dead regions
  02_det1_mesh_scan             gain turn-on; operating point on the plateau
  03_det1_drift_scan            drift turn-on + the zero-field null test
  04_det1_timing_vs_drift       sigma_t vs field, plateau band annotated
  05_all_efficiency             charged-up efficiency
  06_coverage_table             which measurement exists
  07_all_mesh_scans             gain scan
  08_all_drift_scans            drift scan
  09_all_timing_vs_drift        timing vs drift, log axis
  10_det1_loss_budget           reco / fired-not-reco / no-hit
  11_pillar_accounting          measured no-hit vs Gerber pillar area
  12_charging_up                efficiency vs time

The 05/07/08/09/12 file names are kept even though they now carry a single
detector: they are already referenced from the slide deck.

Conventions applied throughout, because both would otherwise mislead:
  * efficiency is the CHARGED-UP (late-window) value where a run contains
    charging-up.
  * timing points from too few events are drawn hollow, never as a resolution.

Usage:
  python3 22_conference_figures.py [-o OUTDIR] [--fmt png,pdf] [--dets det1,...]
"""

import argparse
import glob
import os
import re
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, ScalarFormatter
import numpy as np
import pandas as pd

import p2_qa_config as qa

COL = {'det1': '#1f77b4', 'det2': '#2ca02c', 'det3': '#d62728', 'det4': '#ff7f0e'}
ALL_STD = {
    'det1': ('det1_long5', 'det1_meshscan1', 'det1_driftscan2'),
    'det2': ('det2_long1', 'det2_hvscan1', None),
    'det3': ('det3_initial1', None, 'det3_driftscan1'),
    'det4': ('det4_long2', None, 'det4_driftscan1'),
}
STD = {}                   # filled in main() from --dets
MIN_EV_TIMING = 500
SIZE = (8.6, 5.9)          # a single slide plot
DRIFT_PLATEAU = 500.0      # V/cm; where both eff and sigma_t stop moving

# Slide typography: everything one size up from the paper default, because the
# figure is projected, not read at arm's length.
plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'xtick.labelsize': 17,
    'ytick.labelsize': 17,
    'legend.fontsize': 16,
    'axes.linewidth': 1.3,
    'lines.linewidth': 2.4,
    'lines.markersize': 9,
    'xtick.major.width': 1.3,
    'ytick.major.width': 1.3,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
})


def find(cfg, stage, stem, ext):
    f = [p for p in sorted(glob.glob(os.path.join(cfg.OUT_BASE, stage,
                                                  f'{stem}*spark_vetoed.{ext}')))
         if 'last' not in os.path.basename(p)]
    return f[-1] if f else None


def charging(cfg):
    f = sorted(glob.glob(os.path.join(cfg.OUT_BASE, '20_charging_up',
                                      'charging_summary*.txt')))
    if not f:
        return np.nan, np.nan, ''
    t = open(f[-1]).read()
    def g(p):
        m = re.search(p, t)
        return float(m.group(1)) if m else np.nan
    ver = ('DEGRADING' if 'DEGRADING' in t else
           'RISING' if 'STILL RISING' in t else 'PLATEAUED')
    return (g(r'FULL RUN \(as quoted\)\s*: eff ([\d.]+)'),
            g(r'late\s+\([^)]*\)\s*: eff ([\d.]+)'), ver)


def save(fig, out, name, fmts):
    for f in fmts:
        fig.savefig(os.path.join(out, f'{name}.{f}'), dpi=200,
                    bbox_inches='tight')
    plt.close(fig)
    print(f'  {name}  ({", ".join(fmts)})')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--out', default=None)
    ap.add_argument('--fmt', default='png,pdf')
    ap.add_argument('--dets', default='det1',
                    help='comma-separated detectors to show (default: det1 only)')
    args = ap.parse_args()
    out = args.out or os.path.join(qa.ANALYSIS_ROOT, 'conference')
    os.makedirs(out, exist_ok=True)
    fmts = [f.strip() for f in args.fmt.split(',') if f.strip()]

    STD.clear()
    STD.update({d: ALL_STD[d] for d in
                (x.strip() for x in args.dets.split(',')) if d in ALL_STD})

    c1, m1, d1 = (qa.get_config(k) for k in STD['det1'])

    # -- 01 efficiency map. The stage-10 product is a 3-panel diagnostic with the
    #    run path in the suptitle; for a slide only the efficiency panel is the
    #    figure, so re-draw it alone from the saved grids.
    p = find(c1, '06_efficiency', 'efficiency_map_sliding', 'npz')
    if p:
        z = np.load(p)
        eff, ext = z['eff_within'], z['extent']
        cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('lightgrey')
        fig, ax = plt.subplots(figsize=(8.0, 6.6))
        im = ax.imshow(100 * eff.T, origin='lower', extent=ext, aspect='equal',
                       cmap=cmap, vmin=0, vmax=100)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label('efficiency [%]', fontsize=20)
        cb.ax.tick_params(labelsize=17)
        try:
            from scipy.spatial import ConvexHull
            h = ConvexHull(np.column_stack([z['fp_x'], z['fp_y']]))
            v = np.append(h.vertices, h.vertices[0])
            ax.plot(z['fp_x'][v], z['fp_y'][v], color='red', lw=2.0, alpha=.85)
        except Exception:
            pass
        pcsv = p.replace('efficiency_map_sliding', 'pillars_m3').replace('.npz', '.csv')
        if os.path.isfile(pcsv):
            import p2_mapping as pmap
            pmap.draw_pillars(ax, pd.read_csv(pcsv), small=False)
            ax.legend(loc='lower left', framealpha=.9, fontsize=15)
        ax.set_xlabel('reference X [mm]'); ax.set_ylabel('reference Y [mm]')
        save(fig, out, '01_det1_efficiency_map', fmts)
    else:
        p = find(c1, '06_efficiency', 'efficiency_map_sliding', 'png')
        if p:
            shutil.copyfile(p, os.path.join(out, '01_det1_efficiency_map.png'))
            print('  01_det1_efficiency_map  (copied — re-run stage 10 for the '
                  'single-panel version)')

    # -- 02 det1 mesh scan
    f = find(m1, '11_hv_scan_efficiency', 'efficiency_vs_hv', 'csv')
    if f:
        d = pd.read_csv(f).sort_values('hv')
        fig, ax = plt.subplots(figsize=SIZE)
        ax.axvspan(400, d['hv'].max(), color='seagreen', alpha=.12, zorder=0)
        ax.errorbar(d['hv'], 100 * d['eff_reco'],
                    yerr=100 * d.get('eff_reco_err', 0), fmt='o-',
                    color=COL['det1'], capsize=3, label='reconstructed within R')
        ax.plot(d['hv'], 100 * d['eff_anyhit'], 's--', ms=6, alpha=.75,
                color='grey', label='any pad fired')
        pl = 100 * d.loc[d['hv'] >= 400, 'eff_reco']
        ax.text(0.03, 0.97,
                f'plateau\n$\\varepsilon$ = {pl.mean():.1f} %',
                transform=ax.transAxes, ha='left', va='top',
                fontsize=19, weight='bold', color='seagreen',
                bbox=dict(boxstyle='round,pad=0.35', fc='white',
                          ec='seagreen', alpha=.85))
        ax.set_xlabel('mesh HV [V]'); ax.set_ylabel('efficiency [%]')
        ax.grid(alpha=.3); ax.legend(loc='lower right')
        save(fig, out, '02_det1_mesh_scan', fmts)

    # -- 03/04 det1 drift scan + timing
    f = find(d1, '16_drift_scan_efficiency', 'efficiency_vs_drift', 'csv')
    dd = pd.read_csv(f).sort_values('e_drift_Vcm') if f else None
    if dd is not None:
        plat = dd['e_drift_Vcm'] >= DRIFT_PLATEAU
        eff_plateau = 100 * dd.loc[plat, 'eff_reco'].mean()
        xhi = dd['e_drift_Vcm'].max() * 1.04

        fig, ax = plt.subplots(figsize=SIZE)
        ax.axvspan(DRIFT_PLATEAU, xhi, color='seagreen', alpha=.12, zorder=0)
        ax.errorbar(dd['e_drift_Vcm'], 100 * dd['eff_reco'],
                    yerr=100 * dd.get('eff_reco_err', 0), fmt='o-',
                    color=COL['det1'], capsize=3)
        ax.annotate('zero drift field:\nprimaries never reach the mesh\n'
                    '(null test — it must collapse)',
                    xy=(dd['e_drift_Vcm'].iloc[0], 100 * dd['eff_reco'].iloc[0]),
                    xytext=(230, 44), fontsize=15, color='crimson',
                    arrowprops=dict(arrowstyle='->', color='crimson', lw=2))
        ax.text(0.97, 0.42,
                f'plateau\n$\\varepsilon$ = {eff_plateau:.1f} %',
                transform=ax.transAxes, ha='right', va='top',
                fontsize=20, weight='bold', color='seagreen',
                bbox=dict(boxstyle='round,pad=0.35', fc='white',
                          ec='seagreen', alpha=.85))
        ax.set_xlim(right=xhi)
        ax.set_xlabel('drift field [V/cm]'); ax.set_ylabel('efficiency [%]')
        ax.grid(alpha=.3)
        save(fig, out, '03_det1_drift_scan', fmts)

        if 'time_res_ns' in dd.columns:
            ok = dd['n_active'] >= MIN_EV_TIMING if 'n_active' in dd else slice(None)
            e, t = dd.loc[ok, 'e_drift_Vcm'], dd.loc[ok, 'time_res_ns']
            fig, ax = plt.subplots(figsize=SIZE)
            ax.axvspan(DRIFT_PLATEAU, xhi, color='seagreen', alpha=.12, zorder=0)
            ax.plot(e, t, 'o-', color=COL['det1'])
            pl = t[e >= DRIFT_PLATEAU]
            if len(pl):
                v = float(np.median(pl))
                ax.axhline(v, color='crimson', ls='--', lw=2)
                # the plateau is the same region in both observables -- say so
                # here, so the timing slide does not need the efficiency slide.
                ax.text(0.97, 0.93,
                        f'plateau  ($E_d \\geq$ {DRIFT_PLATEAU:.0f} V/cm)\n'
                        f'$\\varepsilon$ = {eff_plateau:.1f} %\n'
                        f'$\\sigma_t$ = {v:.0f} ns',
                        transform=ax.transAxes, ha='right', va='top',
                        fontsize=19, weight='bold', color='crimson',
                        bbox=dict(boxstyle='round,pad=0.4', fc='white',
                                  ec='crimson', lw=1.6, alpha=.9))
            ax.set_xlim(right=xhi)
            ax.set_xlabel('drift field [V/cm]')
            ax.set_ylabel(r'time resolution $\sigma_t$ [ns]')
            ax.grid(alpha=.3)
            save(fig, out, '04_det1_timing_vs_drift', fmts)

    # -- 05 efficiency, all detectors
    rows = []
    for det, (lk, mk, dk) in STD.items():
        full, late, ver = charging(qa.get_config(lk))
        rows.append(dict(det=det, full=full, late=late, ver=ver,
                         mesh=mk is not None, drift=dk is not None))
    s = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(5.2, 5.9) if len(s) == 1 else SIZE)
    x = np.arange(len(s))
    for i, r in s.iterrows():
        deg = r['ver'] == 'DEGRADING'
        v = r['late'] if np.isfinite(r['late']) else r['full']
        ax.bar(x[i], v, .55, color=COL[r['det']], alpha=.35 if deg else .95,
               hatch='//' if deg else None)
        ax.text(x[i], v + 2.0, f'{v:.1f} %', ha='center', fontsize=22,
                weight='bold')
        if deg:
            ax.text(x[i], v / 2, 'DEGRADING\nnot a stable\nmeasurement',
                    ha='center', fontsize=14, color='crimson', weight='bold')
    ax.set_xticks(x); ax.set_xticklabels(s['det'])
    ax.set_xlim(-.75, len(s) - .25)
    ax.set_ylim(0, 112); ax.set_ylabel('efficiency [%]')
    ax.grid(axis='y', alpha=.3)
    save(fig, out, '05_all_efficiency', fmts)

    # -- 06 coverage table
    fig, ax = plt.subplots(figsize=(9.6, 0.9 + 0.7 * len(s))); ax.axis('off')
    tab = [[r['det'], '✓', '✓' if r['mesh'] else '—', '✓' if r['drift'] else '—',
            'degrading' if r['ver'] == 'DEGRADING' else f"{r['late']:.1f} %"]
           for _, r in s.iterrows()]
    t = ax.table(cellText=tab,
                 colLabels=['', 'long run', 'mesh scan', 'drift scan', 'efficiency'],
                 colWidths=[.14, .20, .22, .22, .22],
                 loc='center', cellLoc='center')
    t.auto_set_font_size(False); t.set_fontsize(18); t.scale(1, 2.6)
    for (row, _), cell in t.get_celld().items():
        cell.set_linewidth(1.2)
        if row == 0:
            cell.set_text_props(weight='bold')
            cell.set_facecolor('#eeeeee')
    save(fig, out, '06_coverage_table', fmts)

    # -- 07 mesh scans, all
    fig, ax = plt.subplots(figsize=SIZE)
    for det, (lk, mk, dk) in STD.items():
        if not mk:
            continue
        f = find(qa.get_config(mk), '11_hv_scan_efficiency', 'efficiency_vs_hv', 'csv')
        if f:
            d = pd.read_csv(f).sort_values('hv')
            ax.errorbar(d['hv'], 100 * d['eff_reco'], fmt='o-', color=COL[det],
                        label=det, capsize=3)
    ax.set_xlabel('mesh HV [V]'); ax.set_ylabel('efficiency [%]')
    ax.grid(alpha=.3); ax.legend()
    save(fig, out, '07_all_mesh_scans', fmts)

    # -- 08/09 drift scans + timing, all
    fig, ax = plt.subplots(figsize=SIZE)
    figt, axt = plt.subplots(figsize=SIZE)
    for det, (lk, mk, dk) in STD.items():
        if not dk:
            continue
        f = find(qa.get_config(dk), '16_drift_scan_efficiency',
                 'efficiency_vs_drift', 'csv')
        if not f:
            continue
        d = pd.read_csv(f).sort_values('e_drift_Vcm')
        ax.errorbar(d['e_drift_Vcm'], 100 * d['eff_reco'], fmt='o-',
                    color=COL[det], label=det, capsize=3)
        if 'time_res_ns' in d.columns:
            low = d['n_active'] < MIN_EV_TIMING if 'n_active' in d.columns \
                else pd.Series(False, index=d.index)
            axt.plot(d.loc[~low, 'e_drift_Vcm'], d.loc[~low, 'time_res_ns'],
                     'o-', color=COL[det], label=det)
            if low.any():
                axt.plot(d.loc[low, 'e_drift_Vcm'], d.loc[low, 'time_res_ns'],
                         'o:', color=COL[det], mfc='none', alpha=.55,
                         label=f'{det} (<{MIN_EV_TIMING} evt — not a resolution)')
    ax.set_xlabel('drift field [V/cm]'); ax.set_ylabel('efficiency [%]')
    ax.grid(alpha=.3); ax.legend()
    save(fig, out, '08_all_drift_scans', fmts)
    axt.set_yscale('log')
    # a decade axis over 20-200 ns labels exactly one tick by default
    axt.set_yticks([20, 30, 50, 100, 200])
    axt.yaxis.set_major_formatter(ScalarFormatter())
    axt.yaxis.set_minor_formatter(NullFormatter())
    axt.set_xlabel('drift field [V/cm]')
    axt.set_ylabel(r'time resolution $\sigma_t$ [ns]')
    axt.grid(alpha=.3, which='both'); axt.legend(fontsize=14)
    save(figt, out, '09_all_timing_vs_drift', fmts)

    # -- 10 loss budget
    fig, ax = plt.subplots(figsize=SIZE)
    vals = [93.4, 3.4, 3.2]
    ax.bar(['reconstructed', 'fired, not\nreconstructed', 'no hit'], vals,
           color=['#2ca02c', '#ff7f0e', '#d62728'])
    for i, v in enumerate(vals):
        ax.text(i, v + 2.0, f'{v:.1f} %', ha='center', fontsize=20, weight='bold')
    ax.set_ylabel('fraction of reference tracks [%]'); ax.set_ylim(0, 112)
    ax.grid(axis='y', alpha=.3)
    save(fig, out, '10_det1_loss_budget', fmts)

    # -- 11 pillar accounting
    fig, ax = plt.subplots(figsize=(6.8, 5.9))
    ax.bar(['measured\nno-hit', 'mesh pillar area\n(from the Gerber)'],
           [3.17, 3.43], width=.55, color=['#d62728', '#7f7f7f'])
    for i, v in enumerate([3.17, 3.43]):
        ax.text(i, v + .14, f'{v:.2f} %', ha='center', fontsize=20, weight='bold')
    ax.set_ylabel('fraction of active area [%]'); ax.set_ylim(0, 4.6)
    ax.grid(axis='y', alpha=.3)
    save(fig, out, '11_pillar_accounting', fmts)

    # -- 12 charging-up
    fig, ax = plt.subplots(figsize=SIZE)
    for det, (lk, mk, dk) in STD.items():
        f = sorted(glob.glob(os.path.join(qa.get_config(lk).OUT_BASE,
                                          '20_charging_up',
                                          'charging_vs_time*.csv')))
        if f:
            d = pd.read_csv(f[-1])
            ax.plot(d['h'], d['eff'], 'o-', ms=6, color=COL[det], label=det)
    ax.set_xlabel('time since run start [h]'); ax.set_ylabel('efficiency [%]')
    ax.grid(alpha=.3); ax.legend()
    save(fig, out, '12_charging_up', fmts)

    print(f'\nwrote to {out}')


if __name__ == '__main__':
    main()
