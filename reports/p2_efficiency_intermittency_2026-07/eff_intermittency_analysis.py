#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Efficiency-intermittency investigation: correlate the P2 efficiency-vs-time
behaviour with the mesh HV current (baseline + sparks) and with the pad
amplitude spectra, across all det1 long runs and the det2 long run.

Figures land in  <DATA_ROOT>/Analysis/reports/p2_efficiency_intermittency/figs
and a numbers file (for the LaTeX report) next to them.
"""
import os
import sys
import glob
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '/local/home/ak271430/Documents/PostDocSaclay/P2_basket_analysis/cosmic_bench_analysis')
import p2_qa_config as qa
import p2_mapping as pmap

C_WITHIN = '#2a9d8f'   # efficiency (within match_r)  - green
C_ANY    = '#e9a13b'   # efficiency (has_any)         - orange
C_IMON   = '#7a1717'   # mesh current                 - dark red
C_SPARK  = '#d62757'   # spark markers                - crimson
C_AMP    = '#39558c'   # amplitude                    - blue

REPORT = os.path.join(qa.DATA_ROOT, 'Analysis', 'reports',
                      'p2_efficiency_intermittency')
FIGS = os.path.join(REPORT, 'figs')
os.makedirs(FIGS, exist_ok=True)

RUNS = [
    # key, label, mesh V, efficiency dir stage tag
    ('det1_long',  'P2_1  6-30-26  (mesh 420 V)', 420),
    ('det1_long2', 'P2_1  7-4-26  (mesh 440 V)', 440),
    ('det1_long3', 'P2_1  7-7-26  (mesh 440 V)', 440),
    ('det1_long4', 'P2_1  7-9-26  (mesh 430 V)', 430),
    ('det2_long1', 'P2_2  7-9-26  (mesh 430 V)', 430),
]

lines_out = []


def log(s=''):
    print(s)
    lines_out.append(s)


def run_start(cfg):
    with open(cfg.run_config_path) as f:
        return pd.to_datetime(json.load(f)['start_time'])


def load_eff(cfg):
    d = os.path.join(cfg.OUT_BASE, '06_efficiency')
    fp = os.path.join(d, f'efficiency_vs_time{cfg.dead_suffix}.csv')
    return pd.read_csv(fp)


def load_mesh_hv(cfg):
    hv = pd.read_csv(cfg.hv_monitor_csv)
    ts = pd.to_datetime(hv['timestamp'])
    t0 = run_start(cfg)
    t_h = (ts - t0).dt.total_seconds() / 3600.0
    return pd.DataFrame({'t_h': t_h,
                         'vmon': hv[f'{cfg.SPARK_CHANNEL} vmon'],
                         'imon': hv[f'{cfg.SPARK_CHANNEL} imon']})


def baseline_curve(hv, bin_h=0.5):
    """Per-bin median (robust vs spark spikes) and p90 of imon."""
    h = hv[hv['t_h'] >= 0].copy()
    h['_b'] = (h['t_h'] / bin_h).astype(int)
    g = h.groupby('_b')['imon']
    return pd.DataFrame({'t_h': g.median().index * bin_h + bin_h / 2,
                         'med': g.median().to_numpy(),
                         'p90': g.quantile(0.90).to_numpy(),
                         'n': g.size().to_numpy()})


def spark_times(hv, thr=2.0):
    s = hv[hv['imon'] >= thr]
    return s['t_h'].to_numpy()


def overlay_figure(key, label, meshV):
    cfg = qa.get_config(key)
    eff = load_eff(cfg)
    hv = load_mesh_hv(cfg)
    bl = baseline_curve(hv)
    sp = spark_times(hv, cfg.SPARK_IMON_THR)
    t0 = run_start(cfg)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6.4), sharex=True,
                                   gridspec_kw={'height_ratios': [1.15, 1],
                                                'hspace': 0.08})
    ax1.plot(eff['t_h'], 100 * eff['eff_any'], 'o-', ms=4, lw=1.6,
             color=C_ANY, label='has any hit')
    ax1.plot(eff['t_h'], 100 * eff['eff_within'], 'o-', ms=4, lw=1.6,
             color=C_WITHIN, label=f'within {cfg.MATCH_R:.0f} mm')
    for t in sp:
        ax1.axvline(t, color=C_SPARK, alpha=0.25, lw=0.8, zorder=0)
    ax1.set_ylabel('efficiency [%]')
    ax1.set_ylim(0, 100)
    ax1.legend(loc='upper right', fontsize=9, framealpha=0.95)
    ax1.grid(alpha=0.25)

    ax2.plot(hv['t_h'], hv['imon'], lw=0.4, color=C_IMON, alpha=0.45,
             label='mesh imon (raw)')
    ax2.plot(bl['t_h'], bl['med'], 'o-', ms=3.5, lw=1.6, color=C_IMON,
             label='30-min median')
    ax2.axhline(cfg.SPARK_IMON_THR, color='k', ls='--', lw=0.9,
                label=f'spark threshold {cfg.SPARK_IMON_THR:.0f} µA')
    ax2.set_ylim(-0.05, 2.4)
    ax2.set_ylabel('mesh current [µA]')
    ax2.set_xlabel(f'time since run start [h]   (run start {t0})')
    ax2.legend(loc='upper right', fontsize=9, framealpha=0.95)
    ax2.grid(alpha=0.25)

    fig.suptitle(f'{label} — efficiency vs mesh current '
                 f'(sparks marked in pink)', y=0.97)
    out = os.path.join(FIGS, f'overlay_{key}.pdf')
    fig.savefig(out, bbox_inches='tight')
    fig.savefig(out.replace('.pdf', '.png'), dpi=140, bbox_inches='tight')
    plt.close(fig)

    # ---- numbers ---------------------------------------------------------
    e = eff.dropna(subset=['eff_any']).copy()
    med_any = e['eff_any'].median()
    win = e[e['eff_any'] > max(3 * med_any, 0.15)]
    log(f'--- {key}  ({label}) ---')
    log(f'  run start            : {t0}')
    log(f'  duration             : {e["t_h"].max() + 0.25:.1f} h')
    log(f'  median eff_any       : {100*med_any:.1f} %   '
        f'median eff_within: {100*e["eff_within"].median():.1f} %')
    log(f'  peak eff_any         : {100*e["eff_any"].max():.1f} %   '
        f'peak eff_within: {100*e["eff_within"].max():.1f} %')
    if len(win) and len(win) < 0.7 * len(e):
        lo, hi = win['t_h'].min() - 0.25, win['t_h'].max() + 0.25
        wall_lo = (t0 + pd.Timedelta(hours=lo)).strftime('%H:%M')
        wall_hi = (t0 + pd.Timedelta(hours=hi)).strftime('%H:%M')
        log(f'  active window        : t = {lo:.2f}-{hi:.2f} h '
        f'({wall_lo}-{wall_hi} wall clock), {len(win)} bins')
        in_w = (sp >= lo) & (sp <= hi)
        log(f'  sparks in window     : {int(in_w.sum())} / {len(sp)}')
        bw = bl[(bl['t_h'] >= lo) & (bl['t_h'] <= hi)]
        bo = bl[(bl['t_h'] < lo) | (bl['t_h'] > hi)]
        log(f'  baseline imon median : window {bw["med"].median():.3f} µA  '
            f'vs outside {bo["med"].median():.3f} µA')
    else:
        log(f'  active window        : none / whole run  '
            f'(bins above 3x median: {len(win)})')
        log(f'  sparks total         : {len(sp)}')
    # per-bin correlation eff_any vs baseline imon
    m = pd.merge_asof(e.sort_values('t_h'), bl.sort_values('t_h'),
                      on='t_h', direction='nearest', tolerance=0.3).dropna()
    if len(m) > 4:
        r = np.corrcoef(m['eff_any'], m['med'])[0, 1]
        log(f'  corr(eff_any, imon)  : r = {r:+.2f}  ({len(m)} bins)')
    log()
    return eff, bl, sp, t0


def amplitude_vs_time(key, tag, bin_h=0.5, veto_thr=None):
    """Median pad-hit amplitude + hit rate per time bin from combined_hits."""
    cfg = qa.get_config(key)
    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy='reverse',
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    feu_set = set(ct.attrs['feus'])
    files = sorted(glob.glob(os.path.join(cfg.combined_hits_dir, '*.root')))
    import uproot
    parts = []
    for fp in files:
        a = uproot.open(f'{fp}:hits').arrays(
            ['eventId', 'trigger_timestamp_ns', 'channel', 'amplitude', 'feu'],
            library='pd')
        a = a[a['feu'].isin(feu_set)]
        parts.append(a)
    df = pd.concat(parts, ignore_index=True)
    df = pmap.attach_pads_to_hits(df, ct)
    df = df[df['mapped'].fillna(False).astype(bool)].copy()
    # drop burst events (P2-internal discharges: >= burst_npads pads at once)
    npads = df.groupby('eventId')['channel'].transform('size')
    df = df[npads < cfg.BURST_NPADS]
    df['t_h'] = df['trigger_timestamp_ns'] / 1e9 / 3600.0
    df['_b'] = (df['t_h'] / bin_h).astype(int)
    g = df.groupby('_b')
    curve = pd.DataFrame({'t_h': g['t_h'].first().index * bin_h + bin_h / 2,
                          'amp_med': g['amplitude'].median().to_numpy(),
                          'amp_p90': g['amplitude'].quantile(.9).to_numpy(),
                          'rate': (g.size() / (bin_h * 3600.0)).to_numpy()})
    return cfg, df, curve


def amp_figure(key, label, window):
    """2-panel: median amplitude vs time + pad-hit rate vs time; and a spectra
    comparison inside/outside the efficiency window."""
    cfg, df, cu = amplitude_vs_time(key, label)
    t0 = run_start(cfg)
    lo, hi = window

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6.0), sharex=True,
                                   gridspec_kw={'hspace': 0.08})
    ax1.plot(cu['t_h'], cu['amp_med'], 'o-', ms=3.5, lw=1.5, color=C_AMP,
             label='median')
    ax1.plot(cu['t_h'], cu['amp_p90'], 'o--', ms=3, lw=1.2, color=C_AMP,
             alpha=0.55, label='90th percentile')
    ax1.set_ylabel('pad-hit amplitude [ADC]')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(alpha=0.25)
    ax2.plot(cu['t_h'], cu['rate'], 'o-', ms=3.5, lw=1.5, color='0.25')
    ax2.set_ylabel('pad-hit rate [Hz]')
    ax2.set_xlabel(f'time since run start [h]   (run start {t0})')
    ax2.grid(alpha=0.25)
    if lo is not None:
        for ax in (ax1, ax2):
            ax.axvspan(lo, hi, color=C_ANY, alpha=0.15, zorder=0)
        ax1.text(0.5 * (lo + hi), ax1.get_ylim()[1] * 0.93, 'efficiency\nwindow',
                 ha='center', va='top', fontsize=8, color='#8a5a00')
    fig.suptitle(f'{label} — pad amplitude and hit rate vs time', y=0.96)
    out = os.path.join(FIGS, f'amplitude_{key}.pdf')
    fig.savefig(out, bbox_inches='tight')
    fig.savefig(out.replace('.pdf', '.png'), dpi=140, bbox_inches='tight')
    plt.close(fig)

    if lo is not None:
        inw = df[(df['t_h'] >= lo) & (df['t_h'] <= hi)]['amplitude']
        out_w = df[(df['t_h'] < lo) | (df['t_h'] > hi)]['amplitude']
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        bins = np.linspace(0, 4200, 105)
        ax.hist(out_w, bins=bins, density=True, histtype='step', lw=1.8,
                color='0.35', label=f'outside window ({len(out_w):,} hits)')
        ax.hist(inw, bins=bins, density=True, histtype='step', lw=1.8,
                color=C_ANY, label=f'inside window ({len(inw):,} hits)')
        ax.set_yscale('log')
        ax.set_xlabel('pad-hit amplitude [ADC]')
        ax.set_ylabel('normalised entries')
        ax.set_title(f'{label} — amplitude spectrum in / out of the '
                     f'efficiency window')
        ax.legend(fontsize=9)
        ax.grid(alpha=0.25)
        out = os.path.join(FIGS, f'amp_spectra_{key}.pdf')
        fig.savefig(out, bbox_inches='tight')
        fig.savefig(out.replace('.pdf', '.png'), dpi=140, bbox_inches='tight')
        plt.close(fig)
        log(f'  {key}: median amp inside window {inw.median():.0f} ADC, '
            f'outside {out_w.median():.0f} ADC; '
            f'rate inside {len(inw):,}, outside {len(out_w):,}')


def main():
    if '--amp-only' not in sys.argv:
        for key, label, meshV in RUNS:
            overlay_figure(key, label, meshV)

    log('=== amplitude analyses ===')
    # det1 7-7 run: window found above (hard-coded from the numbers printed)
    amp_figure('det1_long3', 'P2_1  7-7-26  (mesh 440 V)', (4.0, 6.25))
    # det2 7-9 run: turn-on = first ~1.8 h
    amp_figure('det2_long1', 'P2_2  7-9-26  (mesh 430 V)', (0.0, 1.8))
    # det1 6-30 run: death at ~8.9 h -> compare before/after via window=(0,8.8)
    amp_figure('det1_long', 'P2_1  6-30-26  (mesh 420 V)', (0.0, 8.8))

    with open(os.path.join(REPORT, 'numbers.txt'), 'w') as f:
        f.write('\n'.join(lines_out) + '\n')
    print(f'\nfigures -> {FIGS}')


if __name__ == '__main__':
    main()
