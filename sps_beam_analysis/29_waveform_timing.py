#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
29_waveform_timing.py

Waveform-level timing for the P2 telescope at the SPS -- the beam-arm
counterpart of the cosmic bench 13_timing_waveforms.py, working on the per-FEU
decoded_root files (nt tree). Where 28_timing_qa uses the one `time` number the
standard processing wrote per hit, this stage recomputes the time of arrival
from the raw 16-sample waveforms with several algorithms and applies the two
corrections that matter:

  ftst        the trigger phase within the coarse DREAM clock. The correct
              ns-per-ftst-unit depends on the clock configuration, so it is
              FITTED from the data (slope of mean TOA vs ftst) rather than
              assumed, and validated by the flatness after correction.
  time-walk   leading-edge slewing vs amplitude, corrected with the median
              TOA-vs-log(amplitude) profile.

Beam data is zero-suppressed with on-FEU pedestal subtraction, so waveforms
arrive baseline-subtracted and only fired channels are present -- no offline
pedestal/CNS step (unlike the bench, where the raw non-ZS stream needs both).

The headline number needs no external reference: for events with a hit in BOTH
P2_MID and P2_OUT, the pair time difference cancels the trigger jitter, and
sigma(dt)/sqrt(2) is the single-station resolution (equal-resolution
assumption).

Products:
  per station  (<Analysis>/<det_tag>/<run>/<sub_run>/29_waveform_timing/):
    toa_algorithms_<det>.png   TOA distribution per algorithm (raw / +ftst)
    ftst_validation_<det>.png  mean TOA vs ftst before/after the fitted corr.
    time_walk_<det>.png        TOA vs amplitude + median profile + corrected
  telescope    (<Analysis>/telescope/<run>/<sub_run>/29_waveform_timing/):
    pair_dt_<A>_<B>.png        pair delta-t per algorithm, walk-corrected
    waveform_timing_summary.json

Usage:
  SPS_DATA_ROOT=.../runs SPS_ANALYSIS_ROOT=.../analysis SPS_RUN=<run> \
    python3 29_waveform_timing.py live --sub-run nominal_00 \
    [--amp-min 200] [--amp-max 3500] [--max-events N]
"""

import os
import json
import glob
import argparse

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sps_config as sc
import p2_waveforms as pw

N_CH = 512


def extract_zs_hits(path, n_samples, amp_min, max_events=None, step=4000):
    """Hit waveforms from a ZS decoded file.

    Returns (df, waves): df(eventId, ftst, channel, amp) and waves(n, n_samples)
    baseline-subtracted (the FEU already subtracted pedestals; a residual
    per-waveform baseline = min of the first two samples is removed).
    Only channels PRESENT in the ZS stream are considered -- presence is
    tracked explicitly so genuine zero samples are not confused with absence.
    """
    import uproot
    import awkward as ak
    dfs, wlist = [], []
    seen = 0
    with uproot.open(f'{path}:nt') as t:
        for arrs in t.iterate(['eventId', 'ftst', 'sample', 'channel',
                               'amplitude'], step_size=step, library='ak'):
            ev = np.asarray(arrs['eventId'], dtype=np.int64)
            ft = np.asarray(arrs['ftst'], dtype=np.int64)
            n = len(ev)
            counts = np.asarray(ak.num(arrs['channel']))
            ch = np.asarray(ak.flatten(arrs['channel']), dtype=np.int64)
            sa = np.asarray(ak.flatten(arrs['sample']), dtype=np.int64)
            am = np.asarray(ak.flatten(arrs['amplitude']), dtype=np.float32)
            iev = np.repeat(np.arange(n, dtype=np.int64), counts)
            ok = (ch >= 0) & (ch < N_CH) & (sa >= 0) & (sa < n_samples)
            wf = np.zeros((n, N_CH, n_samples), dtype=np.float32)
            present = np.zeros((n, N_CH), dtype=bool)
            wf[iev[ok], ch[ok], sa[ok]] = am[ok]
            present[iev[ok], ch[ok]] = True
            peak = wf.max(axis=2)
            sel = present & (peak >= amp_min)
            ie, ic = np.nonzero(sel)
            if len(ie):
                # ZS waveforms are re-centred at 256 by the FEU; absent samples
                # (zeros) are below the ZS threshold -> map them to baseline 0.
                w = wf[ie, ic, :]
                w = np.where(w > 0, w - 256.0, 0.0).astype(np.float32)
                dfs.append(pd.DataFrame({'eventId': ev[ie], 'ftst': ft[ie],
                                         'channel': ic,
                                         'amp': w.max(axis=1)}))
                wlist.append(w)
            seen += n
            if max_events is not None and seen >= max_events:
                break
    if not dfs:
        return (pd.DataFrame(columns=['eventId', 'ftst', 'channel', 'amp']),
                np.empty((0, n_samples), np.float32))
    return (pd.concat(dfs, ignore_index=True),
            np.vstack(wlist).astype(np.float32))


ALGOS = {
    'frac30': lambda w: pw.t_frac(w, 0.3),
    'frac50': lambda w: pw.t_frac(w, 0.5),
    'parabola': pw.t_parabola,
    'dcfd': lambda w: pw.t_dcfd(w, frac=0.5, delay=2),
}


def fit_ftst_slope(t_ns, ftst):
    """Slope [ns / ftst unit] of mean TOA vs ftst (the fine-timestamp scale,
    fitted rather than assumed). Returns (slope, ftst_values, mean_t, sem_t).
    NaN TOAs (failed edge crossings) are dropped up front."""
    fin = np.isfinite(t_ns)
    t_ns, ftst = t_ns[fin], ftst[fin]
    vals = np.unique(ftst)
    mt, st, keep = [], [], []
    for v in vals:
        m = ftst == v
        if m.sum() < 50:
            continue
        keep.append(v)
        mt.append(np.mean(t_ns[m]))
        st.append(np.std(t_ns[m]) / np.sqrt(m.sum()))
    keep, mt, st = np.array(keep, float), np.array(mt), np.array(st)
    if len(keep) < 2 or not np.isfinite(mt).all():
        return 0.0, keep, mt, st
    slope = np.polyfit(keep, mt, 1, w=1 / np.maximum(st, 1e-3))[0]
    return float(slope), keep, mt, st


def walk_correct(t_ns, amp, nbins=30):
    """Median-profile time-walk correction in log-amplitude bins.
    Returns (t_corrected, prof_amp, prof_t)."""
    good = np.isfinite(t_ns) & (amp > 0)
    if good.sum() < 100:
        return t_ns, np.array([]), np.array([])
    edges = np.logspace(np.log10(amp[good].min()),
                        np.log10(amp[good].max()) + 1e-6, nbins + 1)
    idx = np.clip(np.digitize(amp, edges) - 1, 0, nbins - 1)
    med_global = np.nanmedian(t_ns[good])
    corr = np.zeros_like(t_ns)
    px, py = [], []
    for b in range(nbins):
        m = good & (idx == b)
        if m.sum() > 30:
            mb = np.nanmedian(t_ns[m])
            corr[idx == b] = mb - med_global
            px.append(np.sqrt(edges[b] * edges[b + 1]))
            py.append(mb)
    return t_ns - corr, np.array(px), np.array(py)


def analyse_station(cfg, det, sub_run, args, n_samples, tps):
    """Full per-station chain. Returns (per-hit DataFrame with corrected TOAs,
    summary dict) or (None, None)."""
    dec_dir = os.path.join(cfg.subrun_dir(sub_run), 'decoded_root')
    files = sorted(glob.glob(os.path.join(dec_dir, f'*_{det.feus[0]:02d}.root')))
    if not files:
        print(f'  {det.name}: no decoded file for FEU {det.feus[0]} in '
              f'{dec_dir} -- fetch a decoded chunk first')
        return None, None
    df, waves = extract_zs_hits(files[0], n_samples, args.amp_min,
                                max_events=args.max_events)
    print(f'  {det.name}: {len(df)} hit waveforms from '
          f'{os.path.basename(files[0])}')
    if len(df) < 500:
        return None, None
    # benchmark selection: clean MIP-like pulses, no saturation
    sel = (df['amp'] >= args.amp_min) & (df['amp'] <= args.amp_max)
    out = cfg.out_dir(det.det_tag, sub_run, '29_waveform_timing')
    summary = {'detector': det.name, 'n_hits': int(len(df)),
               'n_benchmark': int(sel.sum()), 'algorithms': {}}

    toas = {}
    for name, fn in ALGOS.items():
        t_smp = fn(waves)
        t_ns = t_smp * tps
        slope, fv, fm, fs = fit_ftst_slope(t_ns[sel.to_numpy()],
                                           df.loc[sel, 'ftst'].to_numpy())
        t_ftst = t_ns + slope * (np.median(df['ftst']) - df['ftst'])
        t_walk, px, py = walk_correct(t_ftst[sel.to_numpy()],
                                      df.loc[sel, 'amp'].to_numpy())
        s_raw = pw.robust_sigma(t_ns[sel.to_numpy()][np.isfinite(
            t_ns[sel.to_numpy()])])
        s_ftst = pw.robust_sigma(t_ftst[sel.to_numpy()][np.isfinite(
            t_ftst[sel.to_numpy()])])
        s_walk = pw.robust_sigma(t_walk[np.isfinite(t_walk)])
        toas[name] = dict(t_ns=t_ns, t_ftst=t_ftst, slope=slope,
                          fv=fv, fm=fm, fs=fs, t_walk=t_walk, px=px, py=py,
                          s_raw=s_raw, s_ftst=s_ftst, s_walk=s_walk)
        summary['algorithms'][name] = dict(
            ftst_slope_ns=round(slope, 3),
            sigma_raw_ns=round(float(s_raw), 2),
            sigma_ftst_ns=round(float(s_ftst), 2),
            sigma_walk_ns=round(float(s_walk), 2))
        print(f'    {name:9s} ftst {slope:+.2f} ns/unit   sigma '
              f'{s_raw:.1f} -> {s_ftst:.1f} (ftst) -> {s_walk:.1f} ns (walk)')

    # ---- plots ---------------------------------------------------------- #
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    for name, d in toas.items():
        for ax, key in ((axes[0], 't_ns'), (axes[1], 't_ftst')):
            v = d[key][sel.to_numpy()]
            v = v[np.isfinite(v)]
            med = np.median(v)
            ax.hist(v - med, bins=np.linspace(-150, 150, 121),
                    histtype='step', lw=1.5,
                    label=f'{name} (sigma {pw.robust_sigma(v):.0f} ns)')
    axes[0].set_title('raw TOA (median-centred)')
    axes[1].set_title('after fitted ftst correction')
    for ax in axes:
        ax.set_xlabel('TOA - median [ns]'); ax.set_ylabel('hits')
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.suptitle(f'{det.name} TOA algorithms — {sub_run} '
                 f'(amp {args.amp_min:g}-{args.amp_max:g} ADC)', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(out, f'toa_algorithms_{det.det_tag}.png'),
                dpi=150, bbox_inches='tight'); plt.close(fig)

    best = min(toas, key=lambda k: toas[k]['s_walk'])
    d = toas[best]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    axes[0].errorbar(d['fv'], d['fm'] - np.nanmean(d['fm']), yerr=d['fs'],
                     fmt='o-', label=f'raw (slope {d["slope"]:+.2f} ns/unit)')
    tf_sel = d['t_ftst'][sel.to_numpy()]
    ft_sel = df.loc[sel, 'ftst'].to_numpy()
    _, fv2, fm2, fs2 = fit_ftst_slope(tf_sel, ft_sel)
    axes[0].errorbar(fv2, fm2 - np.nanmean(fm2), yerr=fs2, fmt='s--',
                     label='after correction')
    axes[0].set_xlabel('ftst'); axes[0].set_ylabel('mean TOA - avg [ns]')
    axes[0].set_title(f'ftst validation ({best})')
    axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)

    a = df.loc[sel, 'amp'].to_numpy()
    axes[1].scatter(a, tf_sel, s=2, alpha=0.1, c='slategrey', linewidths=0)
    axes[1].plot(d['px'], d['py'], 'o-', color='crimson', ms=4,
                 label='median profile (subtracted)')
    axes[1].set_xscale('log'); axes[1].set_xlabel('amplitude [ADC]')
    axes[1].set_ylabel('TOA (ftst-corr.) [ns]')
    ylo, yhi = np.nanpercentile(tf_sel, [1, 99])
    axes[1].set_ylim(ylo - 20, yhi + 20)
    axes[1].set_title(f'time-walk ({best}): sigma {d["s_ftst"]:.1f} -> '
                      f'{d["s_walk"]:.1f} ns')
    axes[1].legend(fontsize=8); axes[1].grid(alpha=0.3)
    fig.suptitle(f'{det.name} — {sub_run}', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(out, f'time_walk_{det.det_tag}.png'),
                dpi=150, bbox_inches='tight'); plt.close(fig)

    summary['best_algorithm'] = best
    # per-hit table for the pair stage: best algo, ftst+walk corrected
    tbl = df.loc[sel, ['eventId', 'amp']].copy()
    tbl['t'] = d['t_walk']
    return tbl, summary


def pair_dt(tblA, tblB, nameA, nameB, sub_run, out_dir):
    """Event-matched pair time difference; max-amp hit per event per station."""
    a = tblA.sort_values('amp').drop_duplicates('eventId', keep='last')
    b = tblB.sort_values('amp').drop_duplicates('eventId', keep='last')
    m = a.merge(b, on='eventId', suffixes=('_a', '_b'))
    m = m[np.isfinite(m['t_a']) & np.isfinite(m['t_b'])]
    if len(m) < 200:
        print(f'  pair {nameA}-{nameB}: only {len(m)} matched events, skipping')
        return None
    dt = (m['t_a'] - m['t_b']).to_numpy()
    dt = dt - np.median(dt)
    sig = pw.robust_sigma(dt)
    single = sig / np.sqrt(2)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.hist(dt, bins=np.linspace(-120, 120, 121), histtype='stepfilled',
            alpha=0.7, color='steelblue')
    ax.set_xlabel(f't({nameA}) - t({nameB})  [ns]')
    ax.set_ylabel('events')
    ax.set_title(f'{nameA}-{nameB} pair dt — {len(m)} events\n'
                 f'sigma(dt) {sig:.1f} ns  ->  single-station '
                 f'{single:.1f} ns (trigger jitter cancelled)')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'pair_dt_{nameA}_{nameB}.png'),
                dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f'  pair {nameA}-{nameB}: {len(m)} events  sigma(dt) {sig:.1f} ns '
          f'-> single {single:.1f} ns')
    return dict(pair=f'{nameA}-{nameB}', n_events=int(len(m)),
                sigma_dt_ns=round(float(sig), 2),
                sigma_single_ns=round(float(single), 2))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_key', nargs='?', default=None)
    ap.add_argument('--sub-run', required=True)
    ap.add_argument('--amp-min', type=float, default=200.0)
    ap.add_argument('--amp-max', type=float, default=3500.0)
    ap.add_argument('--max-events', type=int, default=None)
    args = ap.parse_args()

    cfg = sc.get_config(args.run_key)
    print(cfg)
    daq = cfg.daq_info()
    n_samples = int(daq.get('n_samples_per_waveform') or 16)
    tps = float(daq.get('sample_period') or 60)
    print(f'  {n_samples} samples x {tps:g} ns')
    dets = cfg.mappable_detectors()

    tables, summaries = {}, []
    for det in dets:
        tbl, s = analyse_station(cfg, det, args.sub_run, args, n_samples, tps)
        if tbl is not None:
            tables[det.name] = tbl
            summaries.append(s)

    tele_out = cfg.out_dir(sc.TELESCOPE_TAG, args.sub_run,
                           '29_waveform_timing')
    pairs = []
    names = list(tables)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            r = pair_dt(tables[names[i]], tables[names[j]],
                        names[i], names[j], args.sub_run, tele_out)
            if r:
                pairs.append(r)
    with open(os.path.join(tele_out, 'waveform_timing_summary.json'),
              'w') as fh:
        json.dump({'run': cfg.RUN, 'sub_run': args.sub_run,
                   'amp_window': [args.amp_min, args.amp_max],
                   'stations': summaries, 'pairs': pairs}, fh, indent=2)
    print(f'  -> {tele_out}')


if __name__ == '__main__':
    main()
