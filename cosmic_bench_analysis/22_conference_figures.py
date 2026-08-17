#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
22_conference_figures.py

Conference plots, ONE PER FILE so they can be dropped into slides individually.
Each is self-contained: the title carries the conditions, because a plot on a
slide has no caption and no neighbouring panel to lean on.

  01_det1_efficiency_map        uniformity, pillar shadows, dead regions
  02_det1_mesh_scan             gain turn-on; operating point on the plateau
  03_det1_drift_scan            drift turn-on + the zero-field null test
  04_det1_timing_vs_drift       sigma_t vs field, plateau annotated
  05_all_efficiency             every detector, charged-up
  06_coverage_table             which measurement exists for which detector
  07_all_mesh_scans             gain scans overlaid
  08_all_drift_scans            drift scans overlaid (det4 flat = no drift field)
  09_all_timing_vs_drift        timing overlaid; sparse points hollow
  10_det1_loss_budget           reco / fired-not-reco / no-hit
  11_pillar_accounting          measured no-hit vs Gerber pillar area
  12_charging_up                efficiency vs time, all detectors

Conventions applied throughout, because both would otherwise mislead:
  * efficiency is the CHARGED-UP (late-window) value where a run contains
    charging-up; det4 differs by 4 points between the two conventions.
  * det3 is shown but marked DEGRADING -- averaging a collapsing detector into
    a performance plot would be the one genuinely misleading thing here.
  * timing points from too few events are drawn hollow, never as a resolution.

Usage:
  python3 22_conference_figures.py [-o OUTDIR] [--fmt png,pdf]
"""

import argparse
import glob
import os
import re
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import p2_qa_config as qa

COL = {'det1': '#1f77b4', 'det2': '#2ca02c', 'det3': '#d62728', 'det4': '#ff7f0e'}
STD = {
    'det1': ('det1_long5', 'det1_meshscan1', 'det1_driftscan2'),
    'det2': ('det2_long1', 'det2_hvscan1', None),
    'det3': ('det3_initial1', None, 'det3_driftscan1'),
    'det4': ('det4_long2', None, 'det4_driftscan1'),
}
MIN_EV_TIMING = 500
SIZE = (8.0, 5.4)          # a single slide plot


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
    args = ap.parse_args()
    out = args.out or os.path.join(qa.ANALYSIS_ROOT, 'conference')
    os.makedirs(out, exist_ok=True)
    fmts = [f.strip() for f in args.fmt.split(',') if f.strip()]

    c1, m1, d1 = (qa.get_config(k) for k in STD['det1'])

    # -- 01 efficiency map: copy the stage-06 product rather than re-render it,
    #    which would only lose resolution.
    p = find(c1, '06_efficiency', 'efficiency_map_sliding', 'png')
    if p:
        shutil.copyfile(p, os.path.join(out, '01_det1_efficiency_map.png'))
        print('  01_det1_efficiency_map  (copied from stage 06)')

    # -- 02 det1 mesh scan
    f = find(m1, '11_hv_scan_efficiency', 'efficiency_vs_hv', 'csv')
    if f:
        d = pd.read_csv(f).sort_values('hv')
        fig, ax = plt.subplots(figsize=SIZE)
        ax.errorbar(d['hv'], 100 * d['eff_reco'],
                    yerr=100 * d.get('eff_reco_err', 0), fmt='o-',
                    color=COL['det1'], capsize=2, label='reconstructed within R')
        ax.plot(d['hv'], 100 * d['eff_anyhit'], 's--', ms=4, alpha=.75,
                color='grey', label='any pad fired')
        ax.axvspan(400, d['hv'].max(), color='seagreen', alpha=.10)
        ax.text(0.98, 0.05, 'plateau', transform=ax.transAxes, ha='right',
                fontsize=10, color='seagreen')
        ax.set_xlabel('mesh HV [V]'); ax.set_ylabel('efficiency [%]')
        ax.grid(alpha=.3); ax.legend(fontsize=9, loc='lower right')
        ax.set_title('det1 — efficiency vs mesh voltage\n'
                     'drift stepped in tandem; operating point 415 V is on the '
                     'plateau', fontsize=11)
        save(fig, out, '02_det1_mesh_scan', fmts)

    # -- 03/04 det1 drift scan + timing
    f = find(d1, '16_drift_scan_efficiency', 'efficiency_vs_drift', 'csv')
    dd = pd.read_csv(f).sort_values('e_drift_Vcm') if f else None
    if dd is not None:
        fig, ax = plt.subplots(figsize=SIZE)
        ax.errorbar(dd['e_drift_Vcm'], 100 * dd['eff_reco'],
                    yerr=100 * dd.get('eff_reco_err', 0), fmt='o-',
                    color=COL['det1'], capsize=2)
        ax.annotate('zero drift field:\nprimaries never reach the mesh\n'
                    '(null test — it must collapse)',
                    xy=(dd['e_drift_Vcm'].iloc[0], 100 * dd['eff_reco'].iloc[0]),
                    xytext=(300, 42), fontsize=9, color='crimson',
                    arrowprops=dict(arrowstyle='->', color='crimson'))
        ax.set_xlabel('drift field [V/cm]'); ax.set_ylabel('efficiency [%]')
        ax.grid(alpha=.3)
        ax.set_title('det1 — efficiency vs drift field\n'
                     'mesh fixed at 415 V; plateau from ~500 V/cm', fontsize=11)
        save(fig, out, '03_det1_drift_scan', fmts)

        if 'time_res_ns' in dd.columns:
            ok = dd['n_active'] >= MIN_EV_TIMING if 'n_active' in dd else slice(None)
            e, t = dd.loc[ok, 'e_drift_Vcm'], dd.loc[ok, 'time_res_ns']
            fig, ax = plt.subplots(figsize=SIZE)
            ax.plot(e, t, 'o-', color=COL['det1'])
            pl = t[e >= 500]
            if len(pl):
                v = float(np.median(pl))
                ax.axhline(v, color='crimson', ls='--')
                ax.text(0.98, 0.90, f'plateau  $\\sigma_t$ = {v:.0f} ns',
                        transform=ax.transAxes, ha='right', color='crimson',
                        fontsize=12, weight='bold')
            ax.set_xlabel('drift field [V/cm]')
            ax.set_ylabel(r'time resolution $\sigma_t$ [ns]')
            ax.grid(alpha=.3)
            ax.set_title('det1 — time resolution vs drift field\n'
                         r'leading pad, $\sigma=(p_{84.1}-p_{15.9})/2$; the '
                         'operating point is on a plateau, not a minimum',
                         fontsize=11)
            save(fig, out, '04_det1_timing_vs_drift', fmts)

    # -- 05 efficiency, all detectors
    rows = []
    for det, (lk, mk, dk) in STD.items():
        full, late, ver = charging(qa.get_config(lk))
        rows.append(dict(det=det, full=full, late=late, ver=ver,
                         mesh=mk is not None, drift=dk is not None))
    s = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=SIZE)
    x = np.arange(len(s))
    for i, r in s.iterrows():
        deg = r['ver'] == 'DEGRADING'
        v = r['late'] if np.isfinite(r['late']) else r['full']
        ax.bar(x[i], v, .6, color=COL[r['det']], alpha=.35 if deg else .95,
               hatch='//' if deg else None)
        ax.text(x[i], v + 1.5, f'{v:.1f}', ha='center', fontsize=12,
                weight='bold')
        if deg:
            ax.text(x[i], v / 2, 'DEGRADING\nnot a stable\nmeasurement',
                    ha='center', fontsize=8.5, color='crimson', weight='bold')
    ax.set_xticks(x); ax.set_xticklabels(s['det'], fontsize=12)
    ax.set_ylim(0, 108); ax.set_ylabel('efficiency [%]')
    ax.grid(axis='y', alpha=.3)
    ax.set_title('P2 cosmic bench — efficiency per detector\n'
                 'charged-up value (late window of each long run)', fontsize=11)
    save(fig, out, '05_all_efficiency', fmts)

    # -- 06 coverage table
    fig, ax = plt.subplots(figsize=(7.6, 2.9)); ax.axis('off')
    tab = [[r['det'], '✓', '✓' if r['mesh'] else '—', '✓' if r['drift'] else '—',
            'degrading' if r['ver'] == 'DEGRADING' else f"{r['late']:.1f} %"]
           for _, r in s.iterrows()]
    t = ax.table(cellText=tab,
                 colLabels=['', 'long run', 'mesh scan', 'drift scan', 'efficiency'],
                 loc='center', cellLoc='center')
    t.auto_set_font_size(False); t.set_fontsize(11); t.scale(1, 1.8)
    ax.set_title('Measurement coverage', fontsize=11)
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
                        label=det, capsize=2)
    ax.set_xlabel('mesh HV [V]'); ax.set_ylabel('efficiency [%]')
    ax.grid(alpha=.3); ax.legend(fontsize=10)
    ax.set_title('Gain scans — det1 plateaus by 400 V; det2 is still rising at '
                 'the top of its scan', fontsize=10.5)
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
                    color=COL[det], label=det, capsize=2)
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
    ax.grid(alpha=.3); ax.legend(fontsize=10)
    ax.set_title('Drift scans — det4 stays flat: no drift field was established\n'
                 '(HV contact to the drift foil)', fontsize=10.5)
    save(fig, out, '08_all_drift_scans', fmts)
    axt.set_yscale('log')
    axt.set_xlabel('drift field [V/cm]')
    axt.set_ylabel(r'time resolution $\sigma_t$ [ns]')
    axt.grid(alpha=.3); axt.legend(fontsize=9)
    axt.set_title('Time resolution vs drift field\n'
                  'hollow points carry too few events to be a resolution',
                  fontsize=10.5)
    save(figt, out, '09_all_timing_vs_drift', fmts)

    # -- 10 loss budget
    fig, ax = plt.subplots(figsize=SIZE)
    vals = [93.4, 3.4, 3.2]
    ax.bar(['reconstructed', 'fired, not\nreconstructed', 'no hit'], vals,
           color=['#2ca02c', '#ff7f0e', '#d62728'])
    for i, v in enumerate(vals):
        ax.text(i, v + 1.8, f'{v:.1f} %', ha='center', fontsize=13, weight='bold')
    ax.set_ylabel('fraction of reference tracks [%]'); ax.set_ylim(0, 108)
    ax.grid(axis='y', alpha=.3)
    ax.set_title('det1 — where the reference tracks go\n'
                 'the 3.4 % is centroid confusion, recoverable in software',
                 fontsize=11)
    save(fig, out, '10_det1_loss_budget', fmts)

    # -- 11 pillar accounting
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    ax.bar(['measured\nno-hit', 'mesh pillar area\n(from the Gerber)'],
           [3.17, 3.43], color=['#d62728', '#7f7f7f'])
    for i, v in enumerate([3.17, 3.43]):
        ax.text(i, v + .12, f'{v:.2f} %', ha='center', fontsize=13, weight='bold')
    ax.set_ylabel('fraction of active area [%]'); ax.set_ylim(0, 4.6)
    ax.grid(axis='y', alpha=.3)
    ax.set_title('det1 — the inefficiency is the mesh pillars\n'
                 '8,661 small + 5 big inside the active area; no fitted parameter',
                 fontsize=11)
    save(fig, out, '11_pillar_accounting', fmts)

    # -- 12 charging-up
    fig, ax = plt.subplots(figsize=SIZE)
    for det, (lk, mk, dk) in STD.items():
        f = sorted(glob.glob(os.path.join(qa.get_config(lk).OUT_BASE,
                                          '20_charging_up',
                                          'charging_vs_time*.csv')))
        if f:
            d = pd.read_csv(f[-1])
            ax.plot(d['h'], d['eff'], 'o-', ms=3.5, color=COL[det], label=det)
    ax.set_xlabel('time since run start [h]'); ax.set_ylabel('efficiency [%]')
    ax.grid(alpha=.3); ax.legend(fontsize=10)
    ax.set_title('Charging-up — det4 gains 17.6 points over 10 h\n'
                 'det3 falls: it was degrading, not charging', fontsize=11)
    save(fig, out, '12_charging_up', fmts)

    print(f'\nwrote to {out}')


if __name__ == '__main__':
    main()
