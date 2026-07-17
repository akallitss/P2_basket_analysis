#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
15_timing_kinematics.py

Where does the P2 timing resolution go? Differential views of the stage-13
per-hit TOA table (toa_hits<sfx>.csv):

  sigma vs amplitude    : leading-edge algorithms should improve ~1/A until
                          the drift-time floor takes over
  sigma vs cosmic angle : inclined tracks illuminate a longer drift column ->
                          earliest-electron time tightens but charge arrives
                          over a wider window; M3 gives the zenith angle
  sigma vs run/sub-run  : stability across long-run chunks and across the
                          mesh-HV scan points (stage 13 must have been run per
                          sub-run first, e.g. via --sub-run scan_mesh_...)

The slewing correction is re-derived on the full long-run sample and applied
before every differential split (so the amplitude trend shows the residual,
not the raw walk).

Products (<Analysis>/<det>/<run>/<sub_run>/13_timing/):
  sigma_vs_amplitude<sfx>.png
  sigma_vs_angle<sfx>.png
  sigma_vs_run<sfx>.png
  timing_kinematics<sfx>.csv / timing_kinematics_summary<sfx>.txt

Usage: python3 15_timing_kinematics.py [run_key] [--algo thr5sig]
"""

import os
import glob
import argparse
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import p2_qa_config as qa
import p2_align as pa
import p2_waveforms as pw


def sigma_err(sig, n):
    """Approximate statistical error on a robust sigma."""
    return sig / np.sqrt(max(2 * (n - 1), 1))


def bin_sigma(x):
    """Core sigma for a (possibly contaminated) bin: the smaller of the
    clipped-Gaussian sigma and the 68% quantile half-width. The clipped fit
    can run away when a bin holds a secondary population (e.g. late noise
    peaks at low amplitude); the quantile width cannot."""
    s1, s2 = pw.robust_sigma(x), pw.q68_half_width(x)
    vals = [v for v in (s1, s2) if np.isfinite(v)]
    return min(vals) if vals else np.nan


def main():
    ap = argparse.ArgumentParser(description='P2 timing kinematics.')
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--algo', default='thr5sig',
                    help='TOA algorithm column to analyse (default thr5sig, '
                         'the stage-13 winner).')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True)
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    out_dir = cfg.out_dir('13_timing')
    sfx = cfg.product_suffix(args.veto_sparks)
    col = f't_{args.algo}'

    toa_csv = f'{out_dir}/toa_hits{sfx}.csv'
    d = pd.read_csv(toa_csv)
    if col not in d.columns:
        raise SystemExit(f'{col} not in {toa_csv} — available: '
                         f'{[c for c in d.columns if c.startswith("t_")]}')
    print(f'  {len(d):,} hits from {toa_csv}')

    # global slewing correction, applied everywhere below
    t_corr, _, _ = pw.slewing_correction(d[col].to_numpy(),
                                         d['amp_max'].to_numpy())
    d['t'] = t_corr
    d = d[np.isfinite(d['t'])]

    # noisy-timing channel mask: a real pad's TOA sits in the trigger-
    # synchronous core; a noisy pad fires at random times (flat t, flat ftst).
    # Drop channels with >=100 hits of which >30% sit far (>200 ns) from the
    # core — e.g. det2's channel 833 (FEU 7), which otherwise fakes a huge
    # sigma in its amplitude band.
    g = d.groupby('channel_id')['t']
    tail_frac = g.apply(lambda x: (np.abs(x) > 200).mean())
    n_ch = g.size()
    noisy = tail_frac.index[(n_ch >= 100) & (tail_frac > 0.3)]
    if len(noisy):
        n0 = len(d)
        d = d[~d['channel_id'].isin(noisy)]
        print(f'  noisy-timing channels dropped: {list(noisy)} '
              f'(-{n0 - len(d):,} hits)')

    # ---------------- M3 zenith angle per event ----------------
    ep = pa.load_m3_endpoints(cfg.m3_tracking_dir)
    dx = ep['X_Up'] - ep['X_Down']
    dy = ep['Y_Up'] - ep['Y_Down']
    dz = np.abs(ep['Z_Up'] - ep['Z_Down'])
    ep['zenith_deg'] = np.degrees(np.arctan2(np.hypot(dx, dy), dz))
    d = d.merge(ep[['eventId', 'zenith_deg']], on='eventId', how='left')
    n_ang = int(d['zenith_deg'].notna().sum())
    print(f'  hits with an M3 track angle: {n_ang:,} '
          f'({100*n_ang/len(d):.0f}%)')

    rows = []

    # ---------------- sigma vs amplitude ----------------
    la = np.log10(d['amp_max'])
    edges = np.percentile(la, np.linspace(0, 100, 13))
    edges = np.unique(edges)
    ctr, sig, err, cnt = [], [], [], []
    for i in range(len(edges) - 1):
        m = (la >= edges[i]) & (la < edges[i + 1])
        if m.sum() < 60:
            continue
        s = bin_sigma(d.loc[m, "t"])
        ctr.append(10 ** (0.5 * (edges[i] + edges[i + 1])))
        sig.append(s); err.append(sigma_err(s, int(m.sum()))); cnt.append(int(m.sum()))
        rows.append(dict(split='amplitude', bin=f'{ctr[-1]:.0f}',
                         sigma_ns=s, n=int(m.sum())))
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.errorbar(ctr, sig, yerr=err, fmt='o-', color='steelblue', lw=1.5)
    ax.set_xscale('log')
    ax.set_xlabel('pulse amplitude [ADC]')
    ax.set_ylabel(f'sigma_t [ns]  ({args.algo}, walk-corrected)')
    ax.set_title(f'{cfg.DET_NAME} timing resolution vs amplitude — '
                 f'{cfg.RUN}/{cfg.SUB_RUN}')
    ax.grid(True, alpha=0.3, which='both')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/sigma_vs_amplitude{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    # ---------------- sigma vs zenith angle ----------------
    da = d[np.isfinite(d['zenith_deg'])]
    aedges = [0, 5, 10, 15, 20, 25, 30, 40, 55]
    actr, asig, aerr = [], [], []
    for i in range(len(aedges) - 1):
        m = (da['zenith_deg'] >= aedges[i]) & (da['zenith_deg'] < aedges[i + 1])
        if m.sum() < 60:
            continue
        s = bin_sigma(da.loc[m, "t"])
        actr.append(0.5 * (aedges[i] + aedges[i + 1]))
        asig.append(s); aerr.append(sigma_err(s, int(m.sum())))
        rows.append(dict(split='zenith', bin=f'{actr[-1]:.0f}deg',
                         sigma_ns=s, n=int(m.sum())))
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.errorbar(actr, asig, yerr=aerr, fmt='s-', color='seagreen', lw=1.5)
    ax.set_xlabel('M3 zenith angle [deg]')
    ax.set_ylabel(f'sigma_t [ns]  ({args.algo}, walk-corrected)')
    ax.set_title(f'{cfg.DET_NAME} timing resolution vs cosmic angle — '
                 f'{cfg.RUN}/{cfg.SUB_RUN}')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/sigma_vs_angle{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    # ---------------- sigma vs run (chunks + HV-scan points) ----------------
    labels, sigs, errs, ns = [], [], [], []
    for ch, g in d.groupby('chunk'):
        s = bin_sigma(g["t"])
        labels.append(f'long run\nchunk {ch:03d}')
        sigs.append(s); errs.append(sigma_err(s, len(g))); ns.append(len(g))
        rows.append(dict(split='chunk', bin=f'{ch:03d}', sigma_ns=s, n=len(g)))

    # HV-scan points: read their stage-13 toa tables (run 13 --sub-run first)
    scan_dirs = sorted(glob.glob(os.path.join(
        cfg.ANALYSIS_ROOT, cfg.DET_TAG, cfg.RUN,
        'scan_mesh_*', '13_timing', f'toa_hits{sfx}.csv')))
    for fp in scan_dirs:
        sub = fp.split(os.sep)[-3]
        mesh = sub.split('_')[2]
        ds = pd.read_csv(fp)
        if col not in ds.columns or len(ds) < 100:
            continue
        if len(noisy):
            ds = ds[~ds['channel_id'].isin(noisy)]
        tsc, _, _ = pw.slewing_correction(ds[col].to_numpy(),
                                          ds['amp_max'].to_numpy())
        s = bin_sigma(tsc)
        labels.append(f'scan\n{mesh}')
        sigs.append(s); errs.append(sigma_err(s, len(ds))); ns.append(len(ds))
        rows.append(dict(split='hv_scan', bin=mesh, sigma_ns=s, n=len(ds)))

    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(labels)), 5))
    x = np.arange(len(labels))
    ax.bar(x, sigs, yerr=errs, color=['tab:blue'] * sum(l.startswith('long') for l in labels)
           + ['tab:orange'] * sum(l.startswith('scan') for l in labels))
    for i, (s, n) in enumerate(zip(sigs, ns)):
        ax.text(i, s + 0.5, f'{s:.1f}\n(n={n:,})', ha='center', fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(f'sigma_t [ns]  ({args.algo}, walk-corrected)')
    ax.set_title(f'{cfg.DET_NAME} timing resolution by run — {cfg.RUN}')
    ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/sigma_vs_run{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    # ------- per-run differentials: sigma vs amplitude / angle BY RUN -------
    # one curve per sample: the three long-run chunks + every HV-scan point
    # (scan points need their stage-13 toa tables and local m3_tracking).
    def diff_curves(dd, zen_col='zenith_deg'):
        out = {}
        laa = np.log10(dd['amp_max'])
        aed = np.percentile(laa, np.linspace(0, 100, 9))
        aed = np.unique(aed)
        pts = []
        for i in range(len(aed) - 1):
            m = (laa >= aed[i]) & (laa < aed[i + 1])
            if m.sum() >= 50:
                pts.append((10 ** (0.5 * (aed[i] + aed[i + 1])),
                            bin_sigma(dd.loc[m, 't']), int(m.sum())))
        out['amp'] = pts
        pts = []
        if zen_col in dd and dd[zen_col].notna().any():
            for lo, hi_ in [(0, 6), (6, 12), (12, 18), (18, 26), (26, 40)]:
                m = (dd[zen_col] >= lo) & (dd[zen_col] < hi_)
                if m.sum() >= 50:
                    pts.append((0.5 * (lo + hi_),
                                bin_sigma(dd.loc[m, 't']), int(m.sum())))
        out['zen'] = pts
        return out

    samples = {}
    for ch, gch in d.groupby('chunk'):
        samples[f'long chunk {ch:03d}'] = diff_curves(gch)
    for fp in scan_dirs:
        sub = fp.split(os.sep)[-3]
        mesh = sub.split('_')[2]
        ds = pd.read_csv(fp)
        if col not in ds.columns or len(ds) < 500:
            continue
        if len(noisy):
            ds = ds[~ds['channel_id'].isin(noisy)]
        tsc, _, _ = pw.slewing_correction(ds[col].to_numpy(),
                                          ds['amp_max'].to_numpy())
        ds = ds.assign(t=tsc)
        ds = ds[np.isfinite(ds['t'])]
        try:
            eps = pa.load_m3_endpoints(os.path.join(cfg.run_dir, sub,
                                                    'm3_tracking_root'))
            eps['zenith_deg'] = np.degrees(np.arctan2(
                np.hypot(eps['X_Up'] - eps['X_Down'],
                         eps['Y_Up'] - eps['Y_Down']),
                np.abs(eps['Z_Up'] - eps['Z_Down'])))
            ds = ds.merge(eps[['eventId', 'zenith_deg']], on='eventId',
                          how='left')
        except Exception as e:
            print(f'  [by-run] no M3 angles for {sub}: {e}')
        samples[f'scan {mesh}'] = diff_curves(ds)

    cmap = plt.get_cmap('viridis')
    for key, xlabel, fname, logx in [
            ('amp', 'pulse amplitude [ADC]', 'sigma_vs_amplitude_by_run', True),
            ('zen', 'M3 zenith angle [deg]', 'sigma_vs_angle_by_run', False)]:
        fig, ax = plt.subplots(figsize=(9.5, 6))
        n_s = max(len(samples) - 1, 1)
        for i, (lab, cur) in enumerate(samples.items()):
            pts = cur[key]
            if not pts:
                continue
            xx, yy, nn = zip(*pts)
            ls = '--' if lab.startswith('long') else '-'
            ax.plot(xx, yy, 'o' + ls, ms=4, lw=1.4, color=cmap(i / n_s),
                    label=lab)
            for x_, y_, n_ in pts:
                rows.append(dict(split=f'{fname}:{lab}', bin=f'{x_:.0f}',
                                 sigma_ns=y_, n=n_))
        if logx:
            ax.set_xscale('log')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f'sigma_t [ns]  ({args.algo}, walk-corrected)')
        ax.set_title(f'{cfg.DET_NAME} timing resolution vs '
                     f'{"amplitude" if key == "amp" else "cosmic angle"} '
                     f'by run — {cfg.RUN}')
        ax.grid(True, alpha=0.3, which='both')
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(f'{out_dir}/{fname}{sfx}.png', dpi=150,
                    bbox_inches='tight')
        plt.close(fig)

    # ---------------- outputs ----------------
    res = pd.DataFrame(rows)
    res.to_csv(f'{out_dir}/timing_kinematics{sfx}.csv', index=False)
    lines = [f'P2 timing kinematics — {cfg.DET_TAG} {cfg.RUN}/{cfg.SUB_RUN} '
             f'[{args.algo}, walk-corrected]',
             f'  overall sigma: {pw.robust_sigma(d["t"]):.1f} ns '
             f'({len(d):,} hits)']
    for split in ['amplitude', 'zenith', 'chunk', 'hv_scan']:
        sub = res[res['split'] == split]
        if not len(sub):
            continue
        lines.append(f'  {split}:')
        for _, r in sub.iterrows():
            lines.append(f'    {r["bin"]:>10s}: {r["sigma_ns"]:5.1f} ns '
                         f'(n={r["n"]:,})')
    txt = '\n'.join(lines)
    print('\n' + txt)
    with open(f'{out_dir}/timing_kinematics_summary{sfx}.txt', 'w') as f:
        f.write(txt + '\n')
    print(f'\nWritten to: {out_dir}')


if __name__ == '__main__':
    main()
