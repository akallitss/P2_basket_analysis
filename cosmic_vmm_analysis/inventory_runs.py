#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Inventory of the cosmic campaign: one row per sub-run (run_N/cfg_*), from what is on disk only.

Per sub-run: start (first capture file name), minutes (first to last file + one capture length), capture files and
GB, chip settings read from the applied ext config (gain sg, peaking time st, neighbour logic sng, threshold
sdt = spi2*64 + spi1 per VMM, channel mask lines split by value sm 1 / sm 0), the read-out connectors and the HV
declaration from run_config.json (the HV is set by hand and not read back: the declaration is the only record; each
distinct text gets a label E0, E1, ... in order of appearance). With --decode, the middle capture file is decoded
(at most --max-gb of it): over-threshold hit rate per VMM, number of channels with hits, fraction of hits at
ADC 1023.

Output: OUT.csv (one row per sub-run) and OUT.json (the same plus per-VMM rates and the HV texts).

Usage: inventory_runs.py RUNS_DIR OUT_STEM [--decode] [--max-gb 0.3] [--runs 1-121]
"""
import argparse
import csv
import glob
import json
import os
import re
from datetime import datetime

import numpy as np

from coincidence_check import read_times, CHAMBER

VMMS = list(range(2, 20))
SG = {0: 0.5, 1: 1, 2: 3, 3: 4.5, 4: 6, 5: 9, 6: 12, 7: 16}
ST = {0: 200, 1: 100, 2: 50, 3: 25}
LINE = re.compile(r'Hyb\s+(\d+)\s+Vmm\s+(\d+)\s+Reg\s+(\S+)\s+Field\s+(\w+)\s+(\d+)')


def chip_settings(path):
    spi = {}
    out = {v: {} for v in VMMS}
    masks = {v: {'sm1': 0, 'sm0': 0} for v in VMMS}
    for line in open(path):
        line = line.split('#')[0]
        m = LINE.search(line)
        if not m:
            continue
        hyb, vmm, reg, field, val = int(m[1]), int(m[2]), m[3], m[4], int(m[5])
        v = 2 * hyb + vmm
        if v not in out:
            continue
        if field == 'sdt':
            spi[(v, reg)] = val
        elif field == 'sm':
            masks[v]['sm1' if val else 'sm0'] += 1
        elif field in ('sg', 'st', 'sng'):
            out[v][field] = val
    for v in VMMS:
        if (v, 'bank1spi1') in spi and (v, 'bank1spi2') in spi:
            out[v]['sdt'] = spi[(v, 'bank1spi2')] * 64 + spi[(v, 'bank1spi1')]
        out[v].update(masks[v])
    return out


def file_time(name):
    m = re.search(r'_(\d{14})\.pcapng$', name)
    return datetime.strptime(m[1], '%Y%m%d%H%M%S') if m else None


def decode(path, max_bytes):
    V, T, C, A = [], [], [], []
    for v, t, ch, adc in read_times(path, {}, with_adc=True, max_bytes=max_bytes):
        V.append(v); T.append(t); C.append(ch); A.append(adc)
    if not V:
        return {}
    V, T, C, A = (np.concatenate(x) for x in (V, T, C, A))
    lo, hi = np.percentile(T, [0.1, 99.9])
    span = (hi - lo) / 0.998 * 1e-9
    res = {'span_s': span}
    for v in VMMS:
        m = V == v
        n = int(m.sum())
        res[v] = {'hz': n / span if span > 0 else None, 'n_ch': int(len(np.unique(C[m]))),
                  'sat': float(np.mean(A[m] == 1023)) if n else None}
    return res


def parse_range(s):
    a, b = s.split('-')
    return set(range(int(a), int(b) + 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('runs_dir')
    ap.add_argument('out_stem')
    ap.add_argument('--decode', action='store_true')
    ap.add_argument('--max-gb', type=float, default=0.3)
    ap.add_argument('--runs', default=None)
    a = ap.parse_args()
    keep = parse_range(a.runs) if a.runs else None

    runs = sorted((int(d.split('_')[-1]), d) for d in glob.glob(os.path.join(a.runs_dir, 'run_*'))
                  if d.split('_')[-1].isdigit())
    hv_labels, rows = {}, []
    for n, rdir in runs:
        if keep and n not in keep:
            continue
        try:
            rc = json.load(open(os.path.join(rdir, 'run_config.json')))
        except (OSError, ValueError):
            rc = {}
        hv = rc.get('lab_hv_state_declared', '(none declared)')
        hv_labels.setdefault(hv, f'E{len(hv_labels)}')
        conns = []
        for det in rc.get('detectors', []):
            for m in det.get('vmm_map', {}).values():
                conns.append(det['name'] + ':' + '/'.join(m.get('connectors', {})))
        for sdir in sorted(glob.glob(os.path.join(rdir, 'cfg_*'))):
            sub = os.path.basename(sdir)
            files = sorted(glob.glob(os.path.join(sdir, 'raw_daq_data', '*.pcapng')))
            sizes = [os.path.getsize(f) for f in files]
            times = [t for t in (file_time(f) for f in files) if t]
            minutes = ((times[-1] - times[0]).total_seconds()
                       + rc.get('vmm_daq_info', {}).get('capture_duration_s', 120)) / 60 if times else 0
            cfgs = glob.glob(os.path.join(rdir, 'p2b-config-*' + sub[len('cfg'):] + '.txt'))
            chips = chip_settings(cfgs[0]) if cfgs else {}
            row = {'run': n, 'subrun': sub[len('cfg_'):], 'start': times[0].isoformat(' ') if times else '',
                   'minutes': round(minutes, 1), 'files': len(files), 'GB': round(sum(sizes) / 1e9, 3),
                   'GB_per_h': round(sum(sizes) / 1e9 / (minutes / 60), 2) if minutes else None,
                   'hv': hv_labels[hv], 'connectors': ' '.join(conns)}
            if chips:
                for key, conv in (('sg', SG), ('st', ST), ('sng', None)):
                    vals = sorted({chips[v].get(key) for v in VMMS if key in chips[v]})
                    row[key] = '/'.join(str(conv.get(x, x) if conv else x) for x in vals)
                row['sm1'] = sum(chips[v]['sm1'] for v in VMMS)
                row['sm0'] = sum(chips[v]['sm0'] for v in VMMS)
                for v in VMMS:
                    row[f'sdt{v}'] = chips[v].get('sdt')
            if a.decode and files:
                f = files[len(files) // 2]
                d = decode(f, int(a.max_gb * 1e9))
                if d:
                    row['decoded'] = os.path.basename(f)
                    for c in (3, 2, 1):
                        row[f'hz_det{c}'] = round(sum(d[v]['hz'] for v in VMMS if CHAMBER[v] == c), 1)
                    row['vmm'] = {v: d[v] for v in VMMS}
            rows.append(row)
            print(n, row['subrun'], row['minutes'], 'min', row['GB'], 'GB', row['hv'],
                  *(f"det{c} {row.get(f'hz_det{c}')} Hz" for c in (3, 2, 1) if f'hz_det{c}' in row), flush=True)

    json.dump({'hv_declarations': {lab: txt for txt, lab in hv_labels.items()}, 'rows': rows},
              open(a.out_stem + '.json', 'w'), indent=1, default=str)
    cols = []
    for r in rows:
        cols += [k for k in r if k != 'vmm' and k not in cols]
    with open(a.out_stem + '.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    for lab, txt in ((lab, txt) for txt, lab in hv_labels.items()):
        print(lab, ':', txt)


if __name__ == '__main__':
    main()
