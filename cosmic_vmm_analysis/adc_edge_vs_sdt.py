#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Per-chip conversion between the threshold register sdt and the hit ADC, from the lower edge of the ADC spectra.

Why: the threshold is set in sdt (DAC) steps, the trims in their own steps, and the hits are recorded in ADC counts of
the VMM's own 10-bit ADC. Every over-threshold hit has a peak above the channel's threshold, so the lower edge of a chip's
hit-ADC spectrum sits at that chip's threshold, in ADC units. Measuring the edge at several sdt values gives, per chip,
    edge_ADC = intercept + slope * sdt          (the July beam data gave one shared slope of ~0.9 ADC per sdt step)
With neighbour-read-out baselines (neighbor_baseline.py, same ADC units) it also gives the margin above the baseline,
edge - baseline, per chip and, with --per-channel, per channel.

The edge does not depend on the chamber gain (a lower gain moves the muon peak, not the threshold), so runs with a gain
problem (the weekend scan, runs 133-169) are still usable here. Only over-threshold hits are used (read_times keeps
flag = 1), so neighbour read-outs of sng 1 runs do not pull the edge down.

Edge estimators, per chip and sub-run (all in ADC counts):
    p1, p2, p5  low percentiles of the hit ADC (p1 is the default for the fit: p0.5 is pulled down by a few isolated hits)
    rise        lowest ADC bin whose count reaches 10 % of the mean bin count between p2 and p10 (insensitive to
                a sparse tail below the edge)
The chip edge is set by its lowest-threshold channels; the per-channel table (--per-channel) gives each channel's own p5.

Near and below the baseline the edge stops following sdt (the spectrum is cut by the baseline itself, and chips held at a
fixed sdt repeat a value), so with --baselines only points whose edge is at least --fit-margin ADC above the baseline
enter the fits.

The spread of margin_ref inside a chip is the spread of the channel thresholds above their baselines: what per-channel
trims would remove (in ADC counts; divide by the slope for sdt steps; ~1.4 sdt steps per trim step).

Inputs: RUNS_DIR as on the DAQ machine / EOS (run_N/p2b-config-*.txt, run_N/cfg_*/raw_daq_data/*.pcapng). sdt per chip
is read from the applied config with inventory_runs.chip_settings. A sample of each sub-run is decoded: --max-mb from
the middle of its capture files.

Outputs (OUT_DIR):
    adc_edge_points.csv     run, subrun, K label, vmm, sdt, hits, p1, p2, p5, rise, [baseline, margin_p1]
    adc_edge_fit.csv        per vmm: points, sdt range, slope and intercept (ADC vs sdt), rms of residuals,
                            nominal mV per ADC count (0.801 mV per sdt step / slope), [baseline median, sdt at baseline]
    adc_edge_channels.csv   (--per-channel) run, subrun, vmm, ch, sdt, hits, p5, [baseline, margin_p5]
    adc_edge_channel_fit.csv (--per-channel) per channel with >= --min-ch-sdt distinct sdt: slope, intercept, rms,
                            edge at the chip's reference sdt (its K = 0 value if scanned, else the median), and with
                            baselines: margin_ref = edge_ref - baseline (the channel's threshold above its own baseline,
                            in ADC counts; the ADC offset of the channel cancels) and sdt_at_baseline
    adc_edge_vs_sdt.png     (--plot) edge vs sdt per chip with the fit

Usage:
    adc_edge_vs_sdt.py RUNS_DIR OUT_DIR [--runs 133-169] [--max-mb 60] [--estimator p1] [--min-hits 1000]
        [--baselines analysis_x/neighbor_baseline.csv] [--fit-margin 10] [--per-channel] [--min-ch-sdt 4] [--plot]
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coincidence_check import read_times, CHAMBER          # noqa: E402
from inventory_runs import chip_settings, parse_range     # noqa: E402

VMMS = list(range(2, 20))
MV_PER_SDT_NOMINAL = 0.801          # vmmsc ThresholdDAC_to_mV slope; each chip has its own, to be measured (monitor ADC)
ESTIMATORS = ('p1', 'p2', 'p5', 'rise')


def k_label(subrun):
    m = re.search(r'_K(p|m)?(\d+)', subrun)
    if not m:
        return ''
    return ('-' if m[1] == 'm' else '+' if m[1] == 'p' else '') + m[2]


def edge_stats(adc):
    adc = adc[adc < 1023]                       # saturated hits carry no edge information
    if len(adc) == 0:
        return None
    p1, p2, p5, p10 = np.percentile(adc, [1, 2, 5, 10])
    h = np.bincount(adc, minlength=1024)
    lo, hi = int(p2), max(int(p10), int(p2) + 1)
    level = h[lo:hi + 1].mean()
    above = np.nonzero(h >= 0.1 * level)[0]
    rise = int(above[0]) if len(above) else None
    return {'p1': float(p1), 'p2': float(p2), 'p5': float(p5), 'rise': rise}


def sample_files(files, max_bytes):
    """capture files from the middle of the sub-run until max_bytes"""
    if not files:
        return []
    mid = len(files) // 2
    order = sorted(files, key=lambda f: abs(files.index(f) - mid))
    out, tot = [], 0
    for f in order:
        out.append(f)
        tot += os.path.getsize(f)
        if tot >= max_bytes:
            break
    return sorted(out)


def decode(files, max_bytes):
    V, C, A = [], [], []
    for f in files:
        state = {}                               # markers are per file stream
        for v, _t, ch, adc in read_times(f, state, with_adc=True, max_bytes=max_bytes):
            V.append(v); C.append(ch); A.append(adc)
    if not V:
        return None
    return np.concatenate(V), np.concatenate(C), np.concatenate(A)


def load_baselines(path):
    """neighbor_baseline.csv -> {(vmm, ch): adc_nb_median}"""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            try:
                out[(int(r['vmm']), int(r['ch']))] = float(r['adc_nb_median'])
            except (KeyError, ValueError):
                continue
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('runs_dir')
    ap.add_argument('out_dir')
    ap.add_argument('--runs', default=None, help='range, e.g. 133-169')
    ap.add_argument('--max-mb', type=float, default=60, help='decoded per sub-run, from the middle of its files')
    ap.add_argument('--estimator', choices=ESTIMATORS, default='p1', help='edge used in the fit')
    ap.add_argument('--min-hits', type=int, default=1000, help='hits of a chip needed to use a point')
    ap.add_argument('--min-ch-hits', type=int, default=200, help='hits of a channel needed for --per-channel')
    ap.add_argument('--baselines', default=None, help='neighbor_baseline.csv (same ADC units)')
    ap.add_argument('--fit-margin', type=float, default=10, help='with baselines: fit only edges >= baseline + this')
    ap.add_argument('--min-ch-sdt', type=int, default=4, help='distinct sdt values needed for a channel fit')
    ap.add_argument('--per-channel', action='store_true')
    ap.add_argument('--plot', action='store_true')
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    keep = parse_range(a.runs) if a.runs else None
    base = load_baselines(a.baselines) if a.baselines else {}
    chip_base = {v: float(np.median([b for (vv, _c), b in base.items() if vv == v]))
                 for v in VMMS if any(vv == v for vv, _c in base)}

    runs = sorted((int(d.split('_')[-1]), d) for d in glob.glob(os.path.join(a.runs_dir, 'run_*'))
                  if d.split('_')[-1].isdigit())
    points, chans = [], []
    for n, rdir in runs:
        if keep and n not in keep:
            continue
        for sdir in sorted(glob.glob(os.path.join(rdir, 'cfg_*'))):
            sub = os.path.basename(sdir)
            cfgs = glob.glob(os.path.join(rdir, 'p2b-config-*' + sub[len('cfg'):] + '.txt'))
            if not cfgs:
                print(f'run_{n} {sub}: no applied config file, skipped', flush=True)
                continue
            chips = chip_settings(cfgs[0])
            files = sorted(glob.glob(os.path.join(sdir, 'raw_daq_data', '*.pcapng')))
            d = decode(sample_files(files, a.max_mb * 1e6), int(a.max_mb * 1e6))
            if d is None:
                print(f'run_{n} {sub}: no hits decoded', flush=True)
                continue
            V, C, A = d
            k = k_label(sub)
            line = []
            for v in VMMS:
                sdt = chips.get(v, {}).get('sdt')
                m = V == v
                if sdt is None or m.sum() == 0:
                    continue
                e = edge_stats(A[m])
                if e is None:
                    continue
                row = {'run': n, 'subrun': sub[len('cfg_'):], 'K': k, 'vmm': v, 'chamber': f'det{CHAMBER[v]}',
                       'sdt': sdt, 'hits': int(m.sum()), **{q: e[q] for q in ESTIMATORS}}
                if v in chip_base:
                    row['baseline'] = round(chip_base[v], 1)
                    row['margin_p1'] = round(e['p1'] - chip_base[v], 1)
                points.append(row)
                line.append(f'{v}:{sdt}->{e[a.estimator]:.0f}')
                if a.per_channel:
                    for ch in np.unique(C[m]):
                        aa = A[m][C[m] == ch]
                        aa = aa[aa < 1023]
                        if len(aa) < a.min_ch_hits:
                            continue
                        r = {'run': n, 'subrun': sub[len('cfg_'):], 'vmm': v, 'ch': int(ch), 'sdt': sdt,
                             'hits': int(len(aa)), 'p5': float(np.percentile(aa, 5))}
                        b = base.get((v, int(ch)))
                        if b is not None:
                            r['baseline'] = b
                            r['margin_p5'] = round(r['p5'] - b, 1)
                        chans.append(r)
            print(f'run_{n} K{k or "?"} ' + ' '.join(line), flush=True)

    def write(name, rows):
        if not rows:
            return
        cols = []
        for r in rows:
            cols += [c for c in r if c not in cols]
        with open(os.path.join(a.out_dir, name), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)

    write('adc_edge_points.csv', points)
    write('adc_edge_channels.csv', chans)

    fits = []
    print(f'\nfit of the {a.estimator} edge against sdt (points with >= {a.min_hits} hits):')
    print(f'{"VMM":>4} {"pts":>4} {"sdt range":>11} {"ADC/sdt":>8} {"intercept":>10} {"rms":>6} {"mV/ADC*":>8}'
          + (f' {"baseline":>9} {"sdt@base":>9}' if chip_base else ''))
    for v in VMMS:
        pts = [p for p in points if p['vmm'] == v and p['hits'] >= a.min_hits and p[a.estimator] is not None
               and (v not in chip_base or p[a.estimator] >= chip_base[v] + a.fit_margin)]
        sdts = sorted({p['sdt'] for p in pts})
        if len(sdts) < 2:
            continue
        x = np.array([p['sdt'] for p in pts], float)
        y = np.array([p[a.estimator] for p in pts], float)
        slope, icpt = np.polyfit(x, y, 1)
        rms = float(np.sqrt(np.mean((y - (icpt + slope * x)) ** 2)))
        f = {'vmm': v, 'chamber': f'det{CHAMBER[v]}', 'points': len(pts), 'distinct_sdt': len(sdts),
             'sdt_min': sdts[0], 'sdt_max': sdts[-1], 'slope_adc_per_sdt': round(slope, 4),
             'intercept_adc': round(icpt, 2), 'rms_adc': round(rms, 2),
             'mV_per_adc_nominal': round(MV_PER_SDT_NOMINAL / slope, 3) if slope > 0 else None}
        if v in chip_base:
            f['baseline_median'] = round(chip_base[v], 1)
            f['sdt_at_baseline'] = round((chip_base[v] - icpt) / slope, 1) if slope > 0 else None
        fits.append(f)
        print(f'{v:>4} {len(pts):>4} {sdts[0]:>5}-{sdts[-1]:<5} {slope:>8.3f} {icpt:>10.1f} {rms:>6.1f} '
              f'{f["mV_per_adc_nominal"] if f["mV_per_adc_nominal"] else float("nan"):>8.3f}'
              + (f' {f.get("baseline_median", float("nan")):>9.1f} {f.get("sdt_at_baseline") or float("nan"):>9.1f}'
                 if chip_base else ''))
    write('adc_edge_fit.csv', fits)
    print('* nominal: 0.801 mV per sdt step / slope; the real mV per sdt step of each chip needs the monitor ADC.')

    if a.per_channel and chans:
        ref = {}
        for v in VMMS:
            k0 = [p['sdt'] for p in points if p['vmm'] == v and p['K'] in ('0', '+0', '-0')]
            allv = [p['sdt'] for p in points if p['vmm'] == v]
            if allv:
                ref[v] = int(np.median(k0)) if k0 else int(np.median(allv))
        groups = {}
        for r in chans:
            groups.setdefault((r['vmm'], r['ch']), []).append(r)
        cfits = []
        for (v, ch), rs in sorted(groups.items()):
            b = base.get((v, ch))
            rs = [r for r in rs if b is None or r['p5'] >= b + a.fit_margin]
            if len({r['sdt'] for r in rs}) < a.min_ch_sdt:
                continue
            x = np.array([r['sdt'] for r in rs], float)
            y = np.array([r['p5'] for r in rs], float)
            slope, icpt = np.polyfit(x, y, 1)
            f = {'vmm': v, 'ch': ch, 'n_sdt': len(set(x)), 'sdt_min': int(x.min()), 'sdt_max': int(x.max()),
                 'slope_adc_per_sdt': round(slope, 4), 'intercept_adc': round(icpt, 2),
                 'rms_adc': round(float(np.sqrt(np.mean((y - icpt - slope * x) ** 2))), 2),
                 'ref_sdt': ref.get(v), 'edge_ref': round(icpt + slope * ref[v], 1) if v in ref else None}
            if b is not None:
                f['baseline'] = b
                if f['edge_ref'] is not None:
                    f['margin_ref'] = round(f['edge_ref'] - b, 1)
                f['sdt_at_baseline'] = round((b - icpt) / slope, 1) if slope > 0 else None
            cfits.append(f)
        write('adc_edge_channel_fit.csv', cfits)
        print(f'\nper-channel fits (p5 edge, >= {a.min_ch_sdt} distinct sdt, >= {a.min_ch_hits} hits per point):')
        print(f'{"VMM":>4} {"ref sdt":>7} {"chans":>6} {"slope med":>9}' + (
            f' {"margin med":>10} {"margin std":>10} {"std in sdt":>10}' if base else f' {"edge std":>9}'))
        for v in VMMS:
            fs = [f for f in cfits if f['vmm'] == v]
            if not fs:
                continue
            sl = float(np.median([f['slope_adc_per_sdt'] for f in fs]))
            if base:
                mg = np.array([f['margin_ref'] for f in fs if f.get('margin_ref') is not None])
                if len(mg) < 3:
                    continue
                print(f'{v:>4} {ref.get(v, ""):>7} {len(mg):>6} {sl:>9.3f} {np.median(mg):>10.1f} {np.std(mg):>10.1f}'
                      f' {np.std(mg) / sl if sl > 0 else float("nan"):>10.1f}')
            else:
                ed = np.array([f['edge_ref'] for f in fs if f['edge_ref'] is not None])
                print(f'{v:>4} {ref.get(v, ""):>7} {len(ed):>6} {sl:>9.3f} {np.std(ed):>9.1f}')
        print('margin = channel threshold above its own baseline at the reference sdt (ADC counts).')

    if a.plot and fits:
        try:
            import matplotlib
        except ImportError:
            print('--plot: matplotlib not available here (on lxplus: source the LCG view), no plot written')
            return
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(3, 6, figsize=(18, 9), sharex=False, sharey=False)
        for ax, v in zip(axs.flat, VMMS):
            pts = [p for p in points if p['vmm'] == v and p['hits'] >= a.min_hits]
            ax.set_title(f'VMM {v} (det{CHAMBER[v]})', fontsize=9)
            if not pts:
                ax.axis('off')
                continue
            ax.plot([p['sdt'] for p in pts], [p[a.estimator] for p in pts], 'o', ms=4, label=a.estimator)
            ax.plot([p['sdt'] for p in pts], [p['rise'] for p in pts], 'x', ms=4, alpha=.6, label='rise')
            f = next((f for f in fits if f['vmm'] == v), None)
            if f:
                xx = np.array([f['sdt_min'], f['sdt_max']])
                ax.plot(xx, f['intercept_adc'] + f['slope_adc_per_sdt'] * xx, '-', lw=1,
                        label=f'{f["slope_adc_per_sdt"]:.2f} ADC/sdt')
            if v in chip_base:
                ax.axhline(chip_base[v], ls=':', c='gray', lw=1, label='nb baseline')
            ax.tick_params(labelsize=7)
            ax.legend(fontsize=6)
        fig.supxlabel('sdt (DAC)')
        fig.supylabel('lower edge of the hit ADC spectrum (ADC counts)')
        fig.tight_layout()
        fig.savefig(os.path.join(a.out_dir, 'adc_edge_vs_sdt.png'), dpi=110)
        print('plot:', os.path.join(a.out_dir, 'adc_edge_vs_sdt.png'))


if __name__ == '__main__':
    main()
