#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Dead time and sparks around the scintillator-tagged muons of ONE chunk (pid = run x 100 + chunk).

Reads the chunk's tag dump (TAGS/tags_run<pid>.npz, efficiency_tp.py --dump-diag) and re-decodes its raw captures with
efficiency_tp.load(), the same function and file order that made the dump, so the hit times share the tag clock.

Per muon tag of each probe chamber (scintillator window --scint, on and off tags, t0 = tag time + latency):
  spark_prev_<k> / spark_next_<k>   time (ns) from the end of the last / to the start of the next SPARK EVENT of the
                                    probe chamber, outside t0 +-2 us (the muon's own large signal is not counted).
                                    Spark event: saturated hits (ADC >= --sat-adc) of the chamber grouped with a 1 us
                                    gap; class k = at least k saturated hits in the group (1, 3, 10).
  spark_any_prev_3                  the same for spark events (k >= 3) in ANY of the three chambers
  dark_<scope>                      length (ns) of the interval without any hit of that scope that contains t0, the
                                    window t0 +-win removed first (silent muons have no live hit there by definition):
                                    t0 - win - last hit before t0 - win, plus first hit after t0 + win - (t0 + win).
                                    scope: hyb (the two chips of the predicted pad's hybrid), cham (all chips of the
                                    probe chamber, live and dead), all3 (the three chambers)
The same dark_* lengths at --n-random random times of the chunk give the reference distribution.

Usage: spark_gap_chunk.py PID TAGS_DIR PSEUDO_DIR OUT_DIR [--scint 140:310] [--sat-adc 1000]
(lxplus, LCG view, VMM_QA_DIR set: source pipeline/env.sh; PYTHONPATH = pipeline/cva)
"""
import argparse
import glob
import os
import sys
import time

import numpy as np

import efficiency_tp as E

DETV = {3: list(range(2, 8)), 2: list(range(8, 14)), 1: list(range(14, 20))}
KS = (1, 3, 10)
EXCL = 2000.0          # ns around t0 where spark events are not counted (the muon's own signal)


def groups(t, gap):
    """start / end / size of the groups of sorted times t with gaps <= gap."""
    if not len(t):
        return np.zeros(0), np.zeros(0), np.zeros(0, int)
    cut = np.nonzero(np.diff(t) > gap)[0]
    s = np.r_[0, cut + 1]
    e = np.r_[cut, len(t) - 1]
    return t[s], t[e], e - s + 1


def prev_end(ends, starts, t0):
    """t0 - end of the last event ending before t0 - EXCL (inf if none)."""
    k = np.searchsorted(ends, t0 - EXCL) - 1
    out = np.full(len(t0), np.inf)
    ok = k >= 0
    out[ok] = t0[ok] - ends[k[ok]]
    return out


def next_start(starts, t0):
    k = np.searchsorted(starts, t0 + EXCL, side='right')
    out = np.full(len(t0), np.inf)
    ok = k < len(starts)
    out[ok] = starts[k[ok]] - t0[ok]
    return out


def dark(tt, t0, win):
    """length of the hit-free interval of the sorted times tt containing [t0 - win, t0 + win], window excluded."""
    a = np.searchsorted(tt, t0 - win) - 1
    b = np.searchsorted(tt, t0 + win, side='right')
    before = np.where(a >= 0, t0 - win - tt[np.maximum(a, 0)], np.inf)
    after = np.where(b < len(tt), tt[np.minimum(b, len(tt) - 1)] - (t0 + win), np.inf)
    return before + after


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pid', type=int)
    ap.add_argument('tags'); ap.add_argument('pseudo'); ap.add_argument('out')
    ap.add_argument('--scint', default='140:310')
    ap.add_argument('--sat-adc', type=int, default=1000)
    ap.add_argument('--n-random', type=int, default=20000)
    a = ap.parse_args()
    t_start = time.time()
    out = f'{a.out}/sg_run{a.pid}.npz'
    if os.path.exists(out):
        print(f'{out} exists'); return
    Z = np.load(f'{a.tags}/tags_run{a.pid}.npz')
    d = sorted(glob.glob(f'{a.pseudo}/run_{a.pid}/cfg_*'))
    files = sorted(glob.glob(f'{d[0]}/raw_daq_data/*.pcapng'))
    V, T, C, A = E.load(','.join(files))
    DET = E.DET
    lo, hi = map(float, a.scint.split(':'))
    res = {}
    rng = np.random.default_rng(a.pid)
    t_rand = np.sort(rng.uniform(T[0] + 1e9, T[-1] - 1e9, a.n_random)) if T[-1] - T[0] > 3e9 else np.zeros(0)
    win = float(Z['p3_win'])
    in3 = DET[V] > 0
    T_all3 = T[in3]
    # spark events per chamber and class
    ev = {}
    for p in (3, 2, 1):
        ts = T[(DET[V] == p) & (A >= a.sat_adc)]
        s, e, n = groups(ts, 1000.0)
        for k in KS:
            m = n >= k
            ev[(p, k)] = (s[m], e[m])
        res[f'ev_p{p}_n'] = n.astype(np.int32)
        res[f'ev_p{p}_t'] = s
    sa = np.concatenate([ev[(p, 3)][0] for p in (3, 2, 1)]); ea = np.concatenate([ev[(p, 3)][1] for p in (3, 2, 1)])
    o = np.argsort(ea); sa, ea = sa[o], ea[o]
    T_cham = {p: T[DET[V] == p] for p in (3, 2, 1)}
    T_hyb = {}
    for h in range(1, 10):
        T_hyb[h] = T[(V == 2 * h) | (V == 2 * h + 1)]
    for p in (3, 2, 1):
        g = lambda k: Z[f'p{p}_{k}']
        dts = g('dt_scint')
        idx = np.nonzero((dts >= lo) & (dts <= hi))[0]
        t0 = g('t_tag')[idx] + float(g('latency'))
        vmm = g('vmm')[idx]
        r = dict(idx=idx, on=g('on')[idx], vmm=vmm, n_any=g('n_any')[idx], prev_chip=g('prev_chip')[idx],
                 wall=g('wall')[idx], r_scale=np.array(float(g('r_scale'))))
        for k in KS:
            s, e = ev[(p, k)]
            r[f'spark_prev_{k}'] = prev_end(e, s, t0)
            r[f'spark_next_{k}'] = next_start(s, t0)
        r['spark_any_prev_3'] = prev_end(ea, sa, t0)
        dh = np.full(len(t0), np.inf)
        for h in np.unique(vmm // 2):
            m = vmm // 2 == h
            dh[m] = dark(T_hyb[int(h)], t0[m], win)
        r['dark_hyb'] = dh
        r['dark_cham'] = dark(T_cham[p], t0, win)
        r['dark_all3'] = dark(T_all3, t0, win)
        # reference: random times, per hybrid of this chamber, the chamber and all three
        if len(t_rand):
            r['rnd_cham'] = dark(T_cham[p], t_rand, win)
            r['rnd_all3'] = dark(T_all3, t_rand, win)
            for v in DETV[p][::2]:
                r[f'rnd_hyb{v // 2}'] = dark(T_hyb[v // 2], t_rand, win)
        for kk, vv in r.items():
            res[f'p{p}_{kk}'] = vv
    res['span_s'] = np.array((T[-1] - T[0]) * 1e-9)
    res['n_hits'] = np.array(len(T))
    np.savez_compressed(out + '.tmp.npz', **res)
    os.replace(out + '.tmp.npz', out)
    print(f'pid {a.pid}: {len(T)} hits, {res["span_s"]:.0f} s, events k>=3 per chamber',
          {p: len(ev[(p, 3)][0]) for p in (3, 2, 1)}, f'{time.time() - t_start:.0f} s', flush=True)


if __name__ == '__main__':
    main()
