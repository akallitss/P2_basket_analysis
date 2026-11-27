#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Probe hit timing before/after the DAQ PC change and the P2_MID <-> P2_OUT read-out swap (25 Nov 2026).

Old: night 23-24 Nov (runs 391-406, old box, default read-out); new: 24-25 Nov (runs 410-424, new box, MID <-> OUT
read-out swapped). REF runs only. From the tag dumps (tags_run*.npz of efficiency_tp --dump-diag):
  * <p>_latency    : per chunk, the probe latency efficiency_tp measured (median probe time - tag time, ns);
  * <p>_fired_nr_* : every live probe hit within --diag-r (15) mm of the prediction and +-2000 ns of tag + latency,
                     with its time from tag + latency (dt) and its channel key (vmm*64+ch).
Per probe and chip: the share of these near hits inside the efficiency window (|dt| <= win, 100 ns), the late share
(win < dt <= 2000), early (dt < -win), the core median and 16-84 % half width of |dt| <= win.
Usage: in_timing_check.py AREA_OLD RUNS_OLD AREA_NEW RUNS_NEW --out JSON   (RUNS as 391-406)
@author: Alexandra Kallitsopoulou
"""
import argparse
import json
import os

import numpy as np

B = '/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus'
PROBE = {'p3': 'P2_IN', 'p2': 'P2_MID', 'p1': 'P2_OUT'}


def load(area, runs):
    lo, hi = map(int, runs.split('-'))
    d = json.load(open(f'{B}/{area}/pseudo/chunks.json'))
    keys = [k for k, v in sorted(d.items()) if lo <= v['run'] <= hi and v['cfg'].endswith('trimmed_Km18_v1_hvm')]
    out = {p: dict(dt=[], key=[], lat=[], win=None) for p in PROBE}
    for k in keys:
        f = f'{B}/{area}/tags/tags_run{k}.npz'
        if not os.path.exists(f):
            continue
        z = np.load(f)
        for p in PROBE:
            if f'{p}_fired_nr_dt' not in z.files:
                continue
            out[p]['dt'].append(z[f'{p}_fired_nr_dt'])
            out[p]['key'].append(z[f'{p}_fired_nr_key'])
            out[p]['lat'].append(float(z[f'{p}_latency']))
            out[p]['win'] = float(z[f'{p}_win'])
    for p in PROBE:
        o = out[p]
        o['dt'] = np.concatenate(o['dt']) if o['dt'] else np.zeros(0)
        o['key'] = np.concatenate(o['key']) if o['key'] else np.zeros(0, int)
    return out, len(keys)


def stats(dt, win):
    if len(dt) < 30:
        return None
    core = dt[np.abs(dt) <= win]
    q16, q50, q84 = np.percentile(core, [16, 50, 84]) if len(core) else (np.nan,) * 3
    return dict(n=int(len(dt)), in_win=float(np.mean(np.abs(dt) <= win)), late=float(np.mean(dt > win)),
                early=float(np.mean(dt < -win)), core_median=float(q50), core_hw=float((q84 - q16) / 2))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('area_old'), ap.add_argument('runs_old'), ap.add_argument('area_new'), ap.add_argument('runs_new')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    O, no = load(a.area_old, a.runs_old)
    N, nn = load(a.area_new, a.runs_new)
    res = {'old': f'{a.area_old} {a.runs_old} ({no} chunks)', 'new': f'{a.area_new} {a.runs_new} ({nn} chunks)',
           'probes': {}}
    for p, name in PROBE.items():
        win = O[p]['win'] or N[p]['win'] or 100.0
        r = {'win': win, 'chips': {}}
        for lab, D in (('old', O), ('new', N)):
            lat = np.array(D[p]['lat'])
            r[f'latency_{lab}'] = dict(median=float(np.median(lat)), p16=float(np.percentile(lat, 16)),
                                       p84=float(np.percentile(lat, 84))) if len(lat) else None
            r[lab] = stats(D[p]['dt'], win)
            vm = D[p]['key'] // 64
            for c in sorted(set(vm.tolist())):
                r['chips'].setdefault(str(c), {})[lab] = stats(D[p]['dt'][vm == c], win)
        res['probes'][name] = r
        lo, ln = r['latency_old'], r['latency_new']
        print(f'== {name}: latency old {lo["median"]:.1f} [{lo["p16"]:.0f},{lo["p84"]:.0f}]  new {ln["median"]:.1f} '
              f'[{ln["p16"]:.0f},{ln["p84"]:.0f}] ns, window +-{win:.0f} ns')
        for lab in ('old', 'new'):
            s = r[lab]
            if s:
                print(f'   {lab}: near hits {s["n"]:6d}  in window {100 * s["in_win"]:5.1f} %  late {100 * s["late"]:5.1f}'
                      f'  early {100 * s["early"]:5.1f}  core median {s["core_median"]:6.1f}  hw {s["core_hw"]:5.1f} ns')
        for c, d in r['chips'].items():
            so, sn = d.get('old'), d.get('new')
            f = lambda s: (f'in {100 * s["in_win"]:5.1f} late {100 * s["late"]:5.1f} med {s["core_median"]:6.1f} '
                           f'(n {s["n"]})') if s else '-'
            print(f'   chip {c:>2}: {f(so)}  ->  {f(sn)}')
    json.dump(res, open(a.out, 'w'), indent=1)
    print('written', a.out)


if __name__ == '__main__':
    main()
