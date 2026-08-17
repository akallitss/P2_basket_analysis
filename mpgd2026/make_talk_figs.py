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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- paths ---- #
BENCH = ('/local/home/ak271430/Documents/PostDocSaclay/data/'
         'Cosmic_Bench/Analysis')
BEAM = ('/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/'
        'SPS_Beam_Test/TB_July2026_H4/analysis/urw_referenced_efficiency')
VMM = ('/local/home/ak271430/Documents/PostDocSaclay/P2_basket_analysis/'
       'sps_beam_analysis/vmm_dream_matching')
NOV = ('/local/home/ak271430/Documents/PostDocSaclay/data/'
       'SPS_Beam_Test/VMM-alinx-data')

# Station / detector colours, kept identical to the pipeline figures so the
# deck reads as one system.
C = {'P2_IN': '#1f77b4', 'P2_MID': '#ff7f0e', 'P2_OUT': '#2ca02c',
     'det1': '#1f77b4', 'det2': '#2ca02c', 'det3': '#d62728',
     'det4': '#ff7f0e'}
SIZE = (9.2, 5.6)

plt.rcParams.update({
    'font.size': 12, 'axes.titlesize': 13, 'axes.labelsize': 12,
    'legend.fontsize': 10.5, 'figure.facecolor': 'white',
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
    (muons/DREAM/uRWELL).  Different detectors and different gas -- the
    claim is the SHAPE and the working point, never the same number."""
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
                    label=f'bench {det} — cosmics, M3 tracks, Ar/iC$_4$H$_{{10}}$ 95/5')

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
                 'bench: cosmics + M3, Ar/iC$_4$H$_{10}$ 95/5   |   '
                 'beam: SPS muons + uRWELL, Ar/CO$_2$/iC$_4$H$_{10}$ 93/5/2\n'
                 'different chambers, different gas — the SHAPE and the '
                 'working point transfer, the value does not',
                 fontsize=11.5)
    ax.legend(loc='lower right', framealpha=0.95, fontsize=9.5)
    save(fig, out, 'bench_beam_mesh')


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
    for f in (fig_bench_beam_mesh, fig_dream_vs_vmm, fig_vmm_threshold,
              fig_snr_matrix, fig_timing_campaigns):
        try:
            f(a.out)
        except Exception as exc:                      # keep going, report
            print(f'  !! {f.__name__}: {type(exc).__name__}: {exc}')


if __name__ == '__main__':
    main()
