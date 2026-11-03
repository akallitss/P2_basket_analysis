#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Bad-time intervals of a lab run from the HV and LV monitor files (sparks, trips, LV glitches), and a check of
those intervals against the chamber data.

Inputs (written by the DAQ in every run since run_170, November 2026):
  <run>/cfg_*/hv_monitor.csv   every 5 s: per electrode power, v0, vmon [V], imon [uA]; `timestamp` = DAQ box clock
                               (Europe/Paris), `remote_timestamp` = banco clock (about 30 s ahead; not used)
  <run>/tdk_lv_monitor.csv     every ~11 s: ALINX 12V, Hybrids 3.3V, Hybrids 2.2V volt_v / curr_a

Flags (each becomes an interval [t - one sampling period - --pre-extra s, t + --post s]: the event happened between
two readings, and the chamber activity it causes lasts a few seconds):
  hv_spark    imon above the run median of that electrode by more than --imon-big uA      (discharge)
  hv_small    imon above the median by more than --imon-small uA, below --imon-big          (small: to be validated)
  hv_vdev     |vmon - v0| > --vdev V while powered                                          (dip, ramp, trip)
  hv_off      power flag 0 on any electrode
  lv_glitch   LV current off its running median (+-5 readings) by more than --lv-di A (2.2V, 3.3V) or 0.5 % (12V),
              or a voltage off by more than 0.5 %
--use selects which flags make the veto (default: hv_spark,hv_small,hv_vdev,hv_off).

Clock: the hit times of read_times() are SRS-clock ns. Each capture file is anchored to wall time with its first
packet timestamp (pcapng EPB), so a hit's wall time = file anchor + (t - first hit t of that file).

Usage:
  spark_veto.py RUN_DIR [RUN_DIR...] --out DIR                 intervals: DIR/veto_<run>.json, DIR/veto_<run>.png
  spark_veto.py RUN_DIR --out DIR --check                      + decode the pcapngs: chamber hit rate and saturated
                                                                 (ADC >= 1000) hits per second around every flag class
                                                                 against random times (lxplus, LCG view)
In other scripts:  from spark_veto import load_veto, WallClock   (see scint_efficiency.py --bad-times)
"""

import argparse
import glob
import json
import os
import struct
from datetime import datetime
from zoneinfo import ZoneInfo

import numpy as np

TZ = ZoneInfo('Europe/Paris')
ELECTRODES = ('P2_IN_Drift', 'P2_IN_Mesh', 'P2_MID_Drift', 'P2_MID_Mesh', 'P2_OUT_Drift', 'P2_OUT_Mesh')
DEFAULT_USE = 'hv_spark,hv_small,hv_vdev,hv_off'


def _epoch(s):
    return datetime.strptime(s.strip(), '%Y-%m-%d %H:%M:%S').replace(tzinfo=TZ).timestamp()


def _read_csv(path):
    import csv
    with open(path) as f:
        return list(csv.DictReader(f))


def _merge(iv):
    iv = sorted(iv)
    out = []
    for a, b, why in iv:
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
            if why not in out[-1][2]:
                out[-1][2].append(why)
        else:
            out.append([a, b, [why]])
    return out


def hv_flags(path, a):
    rows = _read_csv(path)
    t = np.array([_epoch(r['timestamp']) for r in rows])
    flags = []
    for el in ELECTRODES:
        if f'{el} imon' not in rows[0]:
            continue
        num = lambda k: np.array([float(r[f'{el} {k}']) if r[f'{el} {k}'] not in ('', 'nan', 'None') else np.nan for r in rows])
        pw, v0, vm, im = num('power'), num('v0'), num('vmon'), num('imon')
        on = pw == 1
        base = np.nanmedian(im[on]) if on.any() else 0.0
        d = im - base
        for k in np.nonzero(on & (d > a.imon_big))[0]:
            flags.append((t[k], 'hv_spark', el, float(im[k])))
        for k in np.nonzero(on & (d > a.imon_small) & (d <= a.imon_big))[0]:
            flags.append((t[k], 'hv_small', el, float(im[k])))
        for k in np.nonzero(on & (np.abs(vm - v0) > a.vdev))[0]:
            flags.append((t[k], 'hv_vdev', el, float(vm[k] - v0[k])))
        for k in np.nonzero(pw == 0)[0]:
            flags.append((t[k], 'hv_off', el, 0.0))
    period = float(np.median(np.diff(t))) if len(t) > 1 else 5.0
    return flags, period, (float(t[0]), float(t[-1])) if len(t) else None


def lv_flags(path, a):
    rows = _read_csv(path)
    flags, period = [], 11.0
    for name in sorted({r['name'] for r in rows}):
        rr = [r for r in rows if r['name'] == name]
        t = np.array([_epoch(r['timestamp']) for r in rr])
        i = np.array([float(r['curr_a']) for r in rr])
        v = np.array([float(r['volt_v']) for r in rr])
        if len(t) > 1:
            period = float(np.median(np.diff(t)))
        n = len(i)
        med_i = np.array([np.median(i[max(0, k - 5):k + 6]) for k in range(n)])
        med_v = np.array([np.median(v[max(0, k - 5):k + 6]) for k in range(n)])
        di = a.lv_di if '12V' not in name else 0.005 * np.median(i)
        bad = (np.abs(i - med_i) > di) | (np.abs(v - med_v) > 0.005 * np.abs(med_v))
        for k in np.nonzero(bad)[0]:
            flags.append((t[k], 'lv_glitch', name, float(i[k] - med_i[k])))
    return flags, period


def build(run_dir, a):
    hv = glob.glob(os.path.join(run_dir, 'cfg_*', 'hv_monitor.csv'))
    lv = os.path.join(run_dir, 'tdk_lv_monitor.csv')
    flags, iv, span = [], [], None
    use = set(a.use.split(','))
    if hv:
        f, per, span = hv_flags(hv[0], a)
        flags += f
        iv += [(t - per - a.pre_extra, t + a.post, w) for t, w, _, _ in f if w in use]
    if os.path.exists(lv):
        f, per = lv_flags(lv, a)
        flags += f
        iv += [(t - per - a.pre_extra, t + a.post, w) for t, w, _, _ in f if w in use]
    iv = _merge(iv)
    by = {}
    for t, w, el, val in flags:
        by.setdefault(w, {}).setdefault(el, 0)
        by[w][el] += 1
    res = {'run': os.path.basename(os.path.normpath(run_dir)), 'hv_file': hv[0] if hv else None,
           'lv_file': lv if os.path.exists(lv) else None, 'monitor_span': span, 'use': sorted(use),
           'params': {k: getattr(a, k) for k in ('imon_big', 'imon_small', 'vdev', 'lv_di', 'pre_extra', 'post')},
           'flag_counts': by, 'flags': [[t, w, el, val] for t, w, el, val in sorted(flags)],
           'intervals': [[x, y, w] for x, y, w in iv]}
    if span:
        tot = span[1] - span[0]
        clip = sum(max(0.0, min(y, span[1]) - max(x, span[0])) for x, y, _ in iv)
        res['vetoed_s'], res['vetoed_frac'] = clip, clip / tot if tot > 0 else 0.0
    return res


def load_veto(path):
    """intervals [[t0, t1], ...] (epoch s, sorted, merged) from a veto_<run>.json."""
    d = json.load(open(path))
    iv = np.array([[x, y] for x, y, _ in d['intervals']], float).reshape(-1, 2)
    return iv


def in_intervals(t_s, iv):
    """bool array: epoch seconds t_s inside any of the merged sorted intervals iv."""
    if len(iv) == 0:
        return np.zeros(len(t_s), bool)
    k = np.searchsorted(iv[:, 0], t_s, side='right') - 1
    return (k >= 0) & (t_s <= iv[np.maximum(k, 0), 1])


def first_packet_epoch(path):
    """epoch seconds of the first enhanced packet block of a pcapng (honours if_tsresol)."""
    resol = 1e-6
    with open(path, 'rb') as f:
        while True:
            hdr = f.read(8)
            if len(hdr) < 8:
                return None
            bt, bl = struct.unpack('<II', hdr)
            body = f.read(bl - 8)
            if bt == 1:                               # interface description block: look for if_tsresol (9)
                p = 8
                while p + 4 <= len(body) - 4:
                    code, ln = struct.unpack('<HH', body[p:p + 4])
                    if code == 0:
                        break
                    if code == 9 and ln >= 1:
                        r = body[p + 4]
                        resol = 2.0 ** -(r & 0x7F) if r & 0x80 else 10.0 ** -r
                    p += 4 + ((ln + 3) & ~3)
            elif bt == 6:
                hi, lo = struct.unpack('<II', body[4:12])
                return ((hi << 32) | lo) * resol


class WallClock:
    """maps SRS-clock hit times to epoch seconds, one anchor per capture file."""

    def __init__(self):
        self.anchors = []      # (first hit t_ns of the file, epoch s of its first packet)

    def add_file(self, path, first_t_ns):
        e = first_packet_epoch(path)
        if e is not None and first_t_ns is not None:
            self.anchors.append((float(first_t_ns), e))

    def __call__(self, t_ns):
        if not self.anchors:
            raise RuntimeError('WallClock: no anchors')
        A = np.array(sorted(self.anchors))
        k = np.clip(np.searchsorted(A[:, 0], t_ns, side='right') - 1, 0, len(A) - 1)
        return A[k, 1] + (np.asarray(t_ns, float) - A[k, 0]) * 1e-9


def read_run_hits(files, with_adc=True):
    """V, T, C, (A), wall clock. Files are read in order; each one anchored by its first packet."""
    from coincidence_check import read_times
    state, out, clock = {}, [], WallClock()
    for p in files:
        recs = read_times(p, state, with_adc=with_adc)
        ts = [r[1] for r in recs if len(r[1])]
        clock.add_file(p, min(x.min() for x in ts) if ts else None)
        out += recs
    cols = [np.concatenate([r[i] for r in out]) for i in range(4 if with_adc else 3)]
    o = np.argsort(cols[1], kind='stable')
    return [c[o] for c in cols], clock


def check(run_dir, res, out_png, out_json, win=60):
    """per-second chamber hit rate and saturated hits around each flag class vs random times."""
    files = sorted(glob.glob(os.path.join(run_dir, 'cfg_*', 'raw_daq_data', '*.pcapng')))
    (V, T, C, A), clock = read_run_hits(files)
    det = np.zeros(32, int)
    det[2:8], det[8:14], det[14:20] = 3, 2, 1
    ws = clock(T)
    t0, t1 = np.floor(ws.min()), np.ceil(ws.max())
    edges = np.arange(t0, t1 + 1)
    series = {}
    for d, nm in ((3, 'P2_IN'), (2, 'P2_MID'), (1, 'P2_OUT')):
        m = det[V] == d
        series[f'{nm} hits'] = np.histogram(ws[m], edges)[0]
        series[f'{nm} sat'] = np.histogram(ws[m & (A >= 1000)], edges)[0]
    classes = {}
    for t, w, el, _ in res['flags']:
        if t0 + win < t < t1 - win:
            classes.setdefault(f'{w} {el.split("_")[1] if el.startswith("P2_") else el}', []).append(t)
    rng = np.random.default_rng(1)
    classes['random'] = list(rng.uniform(t0 + win, t1 - win, 300))
    lag = np.arange(-win, win)
    out = {'lag_s': lag.tolist(), 'classes': {}}
    for cl, ts in classes.items():
        k = (np.array(ts) - t0).astype(int)
        out['classes'][cl] = {'n': len(ts)}
        for s, y in series.items():
            prof = np.array([y[kk + lag] for kk in k if kk + lag[0] >= 0 and kk + lag[-1] < len(y)]) if len(k) else np.zeros((0, len(lag)))
            out['classes'][cl][s] = prof.mean(0).tolist() if len(prof) else []
    # peak excess of saturated hits within +-15 s of the flag, relative to random
    print(f'{"class":28s} {"n":>4}  saturated hits/s in [-15,+15] s around the flag (random level) -- P2_IN / P2_MID / P2_OUT')
    near = (lag >= -15) & (lag <= 15)
    for cl, d in out['classes'].items():
        if not d.get('P2_IN sat'):
            continue
        vals = ' / '.join(f'{np.max(np.array(d[f"{n} sat"])[near]):7.1f} ({np.mean(out["classes"]["random"][f"{n} sat"]):5.1f})' for n in ('P2_IN', 'P2_MID', 'P2_OUT'))
        print(f'{cl:28s} {d["n"]:4d}  {vals}')
    json.dump(out, open(out_json, 'w'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    show = [c for c in out['classes'] if out['classes'][c].get('P2_IN sat')]
    fig, axes = plt.subplots(2, 3, figsize=(15, 7), sharex=True)
    for j, n in enumerate(('P2_IN', 'P2_MID', 'P2_OUT')):
        for i, s in enumerate(('hits', 'sat')):
            ax = axes[i, j]
            for cl in show:
                ax.plot(lag, out['classes'][cl][f'{n} {s}'], label=f'{cl} ({out["classes"][cl]["n"]})',
                        lw=2 if cl == 'random' else 1, color='k' if cl == 'random' else None)
            ax.set_title(f'{n}: {"all hits" if s == "hits" else "hits with ADC >= 1000"} per s')
            ax.set_yscale('log')
            ax.axvline(0, color='0.6', lw=0.8)
            ax.grid(alpha=0.3)
    axes[1, 1].set_xlabel('time from the monitor reading [s] (DAQ-box clock)')
    axes[0, 0].legend(fontsize=7)
    fig.suptitle(f'{res["run"]}: chamber activity around the HV / LV monitor flags, against random times')
    fig.tight_layout()
    fig.savefig(out_png, dpi=100)


def plot(res, out_png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = _read_csv(res['hv_file']) if res['hv_file'] else []
    lv = _read_csv(res['lv_file']) if res['lv_file'] else []
    fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True)
    if rows:
        t = np.array([_epoch(r['timestamp']) for r in rows])
        T0 = t[0]
        for el in ELECTRODES:
            axes[0].plot((t - T0) / 60, [float(r[f'{el} imon']) for r in rows], lw=0.8, label=el)
            axes[1].plot((t - T0) / 60, [float(r[f'{el} vmon']) - float(r[f'{el} v0']) for r in rows], lw=0.8, label=el)
    else:
        T0 = min(_epoch(r['timestamp']) for r in lv)
    for name in sorted({r['name'] for r in lv}):
        rr = [r for r in lv if r['name'] == name]
        i = np.array([float(r['curr_a']) for r in rr])
        axes[2].plot((np.array([_epoch(r['timestamp']) for r in rr]) - T0) / 60, i - np.median(i), lw=0.8, label=name)
    for ax in axes:
        for x, y, _ in res['intervals']:
            ax.axvspan((x - T0) / 60, (y - T0) / 60, color='#d62728', alpha=0.15, lw=0)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=6)
    axes[0].set_ylabel('imon [uA]')
    axes[0].set_yscale('symlog', linthresh=0.05)
    axes[1].set_ylabel('vmon - v0 [V]')
    axes[1].set_yscale('symlog', linthresh=1)
    axes[2].set_ylabel('LV I - median [A]')
    axes[2].set_xlabel('minutes from the first HV reading')
    fig.suptitle(f'{res["run"]}: vetoed {res.get("vetoed_frac", 0) * 100:.1f} % of the time ({", ".join(res["use"])})')
    fig.tight_layout()
    fig.savefig(out_png, dpi=100)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--out', required=True)
    ap.add_argument('--imon-big', type=float, default=1.0)
    ap.add_argument('--imon-small', type=float, default=0.15)
    ap.add_argument('--vdev', type=float, default=1.5)
    ap.add_argument('--lv-di', type=float, default=0.01)
    ap.add_argument('--post', type=float, default=5.0)
    ap.add_argument('--pre-extra', type=float, default=1.0)
    ap.add_argument('--use', default=DEFAULT_USE)
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--no-plot', action='store_true')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for r in a.runs:
        res = build(r, a)
        name = res['run']
        json.dump(res, open(os.path.join(a.out, f'veto_{name}.json'), 'w'), indent=1)
        print(f'{name}: {len(res["intervals"])} intervals, vetoed {res.get("vetoed_s", 0):.0f} s '
              f'({res.get("vetoed_frac", 0) * 100:.2f} %); flags {json.dumps(res["flag_counts"])}')
        if not a.no_plot:
            plot(res, os.path.join(a.out, f'veto_{name}.png'))
        if a.check:
            check(r, res, os.path.join(a.out, f'check_{name}.png'), os.path.join(a.out, f'check_{name}.json'))


if __name__ == '__main__':
    main()
