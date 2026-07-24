#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
30_raw_stream_efficiency.py

Efficiency-per-trigger straight from the RAW decoded ZS stream, bypassing the
offline processing entirely. Motivation (2026-07-24): at nominal HV the raw
stream shows P2_OUT/P2_MID registering an in-time >=200 ADC waveform in 97.6%
of triggers, while the processed combined-hits chain retains only 74% / 52% --
i.e. the tag-probe "turn-on that never plateaus" is dominated by an
amplitude-dependent OFFLINE PROCESSING loss, not by detector gain. This stage
measures the turn-on the processing cannot distort.

Definition per station and sub_run:
    eff_raw(thr) = fraction of recorded triggers with >=1 ZS waveform whose
                   baseline-subtracted peak >= thr ADC AND whose peak sample
                   falls in the in-time window (default samples 4..11 -- the
                   signal sits at 6-7 for latency 32).
The in-time requirement kills accidental/noise ZS packets; ZS itself already
applies a ~5-sigma channel threshold on the FEU. Several thr values are
reported so the low-gain end of a scan is not biased by the cut itself.

Needs the per-FEU decoded chunk(s) locally -- one chunk (~2M triggers) per
sub_run is plenty. fetch e.g.:  rsync banco:...<sub_run>/decoded_root/*_000_0{4,5}.root

Products (<Analysis>/<det_tag>/<run>/scan/30_raw_stream_efficiency/):
  raw_eff_vs_hv_<det>.png    raw turn-on per threshold, processed curve overlaid
                             when the stage-22 CSV exists
  raw_stream_efficiency_<det>.csv

Usage:
  SPS_DATA_ROOT=.../runs SPS_ANALYSIS_ROOT=.../analysis SPS_RUN=<run> \
    python3 30_raw_stream_efficiency.py live [--thr 50 100 200]
    [--t-window 4 11] [--max-events 600000] [--sub-run NAME]
"""

import os
import csv
import glob
import argparse

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sps_config as sc

N_CH = 512
ZS_BASELINE = 256.0     # the FEU re-centres ZS pedestal-subtracted waveforms here


def raw_share(path, n_samples, thresholds, t_lo, t_hi, max_events):
    """Fraction of triggers with an in-time waveform above each threshold."""
    import uproot
    import awkward as ak
    n_seen = 0
    fired = {t: 0 for t in thresholds}
    with uproot.open(f'{path}:nt') as t:
        for arrs in t.iterate(['eventId', 'sample', 'channel', 'amplitude'],
                              step_size=4000, library='ak'):
            ev = np.asarray(arrs['eventId'])
            n = len(ev)
            n_seen += n
            counts = np.asarray(ak.num(arrs['channel']))
            ch = np.asarray(ak.flatten(arrs['channel']), np.int64)
            sa = np.asarray(ak.flatten(arrs['sample']), np.int64)
            am = np.asarray(ak.flatten(arrs['amplitude']), np.float32)
            iev = np.repeat(np.arange(n), counts)
            ok = (sa >= 0) & (sa < n_samples) & (ch >= 0) & (ch < N_CH)
            wf = np.zeros((n, N_CH, n_samples), np.float32)
            pres = np.zeros((n, N_CH), bool)
            wf[iev[ok], ch[ok], sa[ok]] = am[ok]
            pres[iev[ok], ch[ok]] = True
            # ZS + on-FEU pedestal subtraction: every present sample is
            # re-centred at 256, and absent samples are BELOW the ZS threshold
            # by construction. Signal = present amplitude - 256; absent = 0.
            w = np.where(wf > 0, wf - ZS_BASELINE, 0.0).astype(np.float32)
            peak = w.max(axis=2)
            pk_s = w.argmax(axis=2)
            intime = (pk_s >= t_lo) & (pk_s <= t_hi)
            for thr in thresholds:
                fired[thr] += int((pres & intime & (peak >= thr))
                                  .any(axis=1).sum())
            if max_events and n_seen >= max_events:
                break
    return n_seen, fired


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_key', nargs='?', default=None)
    ap.add_argument('--sub-run', default=None)
    ap.add_argument('--thr', type=float, nargs='+', default=[50, 100, 200])
    ap.add_argument('--t-window', type=int, nargs=2, default=[4, 11],
                    metavar=('LO', 'HI'),
                    help='in-time peak-sample window (default 4 11).')
    ap.add_argument('--max-events', type=int, default=600000)
    args = ap.parse_args()

    cfg = sc.get_config(args.run_key)
    print(cfg)
    n_samples = int(cfg.daq_info().get('n_samples_per_waveform') or 16)
    dets = cfg.mappable_detectors()
    subs = [args.sub_run] if args.sub_run else sorted(
        os.path.basename(d.rstrip('/')) for d in
        glob.glob(os.path.join(cfg.run_dir, '*/'))
        if glob.glob(os.path.join(d, 'decoded_root', '*.root')))
    if not subs:
        raise SystemExit('No sub_run with local decoded_root chunks.')
    print(f'  sub_runs with decoded chunks: {subs}')
    t_lo, t_hi = args.t_window

    for det in dets:
        rows = []
        print(f'\n== {det.name} (FEU {det.feus[0]})')
        for sub in subs:
            files = sorted(glob.glob(os.path.join(
                cfg.subrun_dir(sub), 'decoded_root',
                f'*_{det.feus[0]:02d}.root')))
            if not files:
                continue
            n, fired = raw_share(files[0], n_samples, args.thr, t_lo, t_hi,
                                 args.max_events)
            mesh_v = cfg.subrun_mesh_hv(sub, det)
            row = dict(sub_run=sub, mesh_v=mesh_v, n_triggers=n)
            for thr in args.thr:
                row[f'eff_thr{thr:g}'] = fired[thr] / n if n else None
            rows.append(row)
            effs = '  '.join(f'thr{thr:g}={fired[thr]/n:.3f}'
                             for thr in args.thr)
            print(f'  {sub}: mesh {mesh_v} V  {n} trig  {effs}')
        if not rows:
            continue

        out = cfg.out_dir(det.det_tag, 'scan', '30_raw_stream_efficiency')
        with open(os.path.join(out,
                               f'raw_stream_efficiency_{det.det_tag}.csv'),
                  'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(sorted(rows, key=lambda r: (r['mesh_v'] or 0)))

        r = [x for x in rows if x['mesh_v'] is not None]
        if len(r) >= 2:
            r.sort(key=lambda x: x['mesh_v'])
            v = np.array([x['mesh_v'] for x in r], float)
            fig, ax = plt.subplots(figsize=(8.5, 5.5))
            for thr in args.thr:
                y = [x[f'eff_thr{thr:g}'] for x in r]
                ax.plot(v, y, 'o-', lw=2, ms=7,
                        label=f'raw stream, in-time peak >= {thr:g} ADC')
            # processed tag-probe curve for comparison, if present
            tp = glob.glob(os.path.join(
                cfg.ANALYSIS_ROOT, det.det_tag, cfg.RUN, 'scan',
                '22_tag_probe_efficiency', 'tag_probe_efficiency*.csv'))
            if tp:
                t22 = pd.read_csv(tp[0])
                col = 'eff_corr' if t22.get('eff_corr') is not None and \
                    t22['eff_corr'].notna().any() else 'eff'
                t22 = t22.dropna(subset=['hv']).sort_values('hv')
                ax.plot(t22['hv'], t22[col], 's--', lw=1.5, ms=6, color='k',
                        alpha=0.6, label='processed tag-probe (stage 22)')
            ax.set_xlabel('mesh HV [V]')
            ax.set_ylabel('fraction of triggers with in-time signal')
            ax.set_ylim(0, 1.02)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
            ax.set_title(f'{det.name} -- RAW-stream efficiency per trigger '
                         f'vs mesh HV\n(peak sample in [{t_lo},{t_hi}]; '
                         f'gap between raw and processed = offline '
                         f'processing loss)', fontsize=10)
            fig.tight_layout()
            fig.savefig(os.path.join(out, f'raw_eff_vs_hv_{det.det_tag}.png'),
                        dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f'  -> {out}')


if __name__ == '__main__':
    main()
