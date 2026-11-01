#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Chamber response to the scintillator trigger in time slices of one run: separates a gain drift from the threshold.

The capture files (one per minute, wall time in the file name) are read --slice files at a time. Per slice and chamber:
the fraction of scintillator triggers with a hit anywhere in the chamber within --half-win ns of the latency (random
hits from a window --sideband ns earlier subtracted), and the median ADC of those hits (the largest per trigger): the
pulse height of muons, a gain monitor that does not depend on the chamber thresholds as long as the pulses sit well
above them. Trigger: the hardware triple (--trig) or the offline (--top OR) AND --bottom coincidence.

Usage: scint_timeline.py OUT.json run_N:f1.pcapng,f2,... [run_M:...] [--slice 30] [--top 0:41,0:45 --bottom 0:44]
"""
import argparse
import json
import re
import sys
from datetime import datetime

import numpy as np

from coincidence_check import read_times

DET = np.zeros(32, int)
DET[2:8], DET[8:14], DET[14:20] = 3, 2, 1
NAME = {3: 'det3 P2_IN', 2: 'det2 P2_MID', 1: 'det1 P2_OUT'}


def ftime(p):
    m = re.search(r'_(\d{14})\.pcapng$', p)
    return datetime.strptime(m[1], '%Y%m%d%H%M%S')


def triggers(V, T, C, a):
    if a.top:
        def times(spec):
            return np.sort(np.concatenate([T[(V == int(s.split(':')[0])) & (C == int(s.split(':')[1]))] for s in spec.split(',')]))
        top, bot = times(a.top), times(a.bottom)
        s = np.searchsorted(bot, top - a.coinc_win)
        e = np.searchsorted(bot, top + a.coinc_win)
        tr = top[e > s]
        return tr[np.concatenate([[True], np.diff(tr) > 300])] if len(tr) else tr
    tv, tc = (int(x) for x in a.trig.split(':'))
    return T[(V == tv) & (C == tc)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out_json')
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--slice', type=int, default=30, help='capture files per slice (1 file = 1 minute)')
    ap.add_argument('--trig', default='0:40')
    ap.add_argument('--top', default='')
    ap.add_argument('--bottom', default='')
    ap.add_argument('--coinc-win', type=float, default=60)
    ap.add_argument('--latency', type=float, default=230)
    ap.add_argument('--half-win', type=float, default=150)
    ap.add_argument('--sideband', type=float, default=3000)
    a = ap.parse_args()
    out = {'params': vars(a), 'slices': []}
    for spec in a.runs:
        name, files = spec.split(':', 1)
        files = sorted(files.split(','))
        for i in range(0, len(files), a.slice):
            part = files[i:i + a.slice]
            state = {}
            V, T, C, A = [], [], [], []
            for p in part:
                for v, t, c, adc in read_times(p, state, with_adc=True):
                    V.append(v); T.append(t); C.append(c); A.append(adc)
            if not V:
                continue
            V, T, C, A = (np.concatenate(x) for x in (V, T, C, A))
            o = np.argsort(T, kind='stable')
            V, T, C, A = V[o], T[o], C[o], A[o]
            tr = triggers(V, T, C, a)
            row = {'run': name, 'start': ftime(part[0]).isoformat(' '), 'files': len(part), 'n_trig': int(len(tr))}
            for d in (3, 2, 1):
                m = DET[V] == d
                ht, ha = T[m], A[m]
                lo, hi = a.latency - a.half_win, a.latency + a.half_win
                s = np.searchsorted(ht, tr + lo)
                e = np.searchsorted(ht, tr + hi, side='right')
                s2 = np.searchsorted(ht, tr + lo - a.sideband)
                e2 = np.searchsorted(ht, tr + hi - a.sideband, side='right')
                on, acc = (np.mean(e > s), np.mean(e2 > s2)) if len(tr) else (np.nan, np.nan)
                eff = (on - acc) / (1 - acc) if len(tr) else np.nan
                amax = [int(ha[x:y].max()) for x, y in zip(s, e) if y > x]
                row[NAME[d]] = {'resp': float(eff), 'err': float(np.sqrt(max(eff * (1 - eff), 0) / len(tr))) if len(tr) else None,
                                'adc_median': float(np.median(amax)) if amax else None, 'n_hit': len(amax),
                                'rate_hz': float(m.sum() / (len(part) * 60))}
            out['slices'].append(row)
            print(f'{name} {row["start"][5:16]} trig {len(tr):5d}  ' + '  '.join(
                f'{NAME[d].split()[0]} {row[NAME[d]]["resp"]:.3f} adc {row[NAME[d]]["adc_median"] or 0:5.0f}' for d in (3, 2, 1)), flush=True)
    json.dump(out, open(a.out_json, 'w'), indent=1, default=float)


if __name__ == '__main__':
    sys.exit(main())
