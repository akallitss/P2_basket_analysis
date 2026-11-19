#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Is a chamber blind after a spark? Hit rate of the chamber in log-spaced bins of time AFTER each spark event and,
as the reference, the same bins BEFORE it (mirror), summed over the spark events of a chunk.

Spark event: saturated hits (ADC >= 1000) of one chamber grouped with a 1 us gap, at least K of them (K = 3, 10).
Counted after the event end / before the event start, the event's own hits excluded:
  cham   all hits of that chamber (live and dead channels), the pulsed channels of a pulser run excluded
  pulse  (pulser runs, *_tp*) test pulses of the chamber's pulsed chips (channel 31 of chips 2 5 / 8 11 / 14 17)
  other  all hits of the OTHER two chambers (is the blindness local to the sparking chamber?)
Windows that run past the chunk edges are dropped per bin (exposure counted per bin).
Output: CHUNKS/sp_run<pid>.npz with counts[after/before][K][scope][det] per bin and the exposures.

Usage: spark_profile_chunk.py PID PSEUDO_DIR OUT_DIR
"""
import argparse
import glob
import time

import numpy as np

import efficiency_tp as E

EDGES = np.array([0, 300, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9])
PULSED = {3: (2, 5), 2: (8, 11), 1: (14, 17)}


def groups(t, gap):
    if not len(t):
        return np.zeros(0), np.zeros(0), np.zeros(0, int)
    cut = np.nonzero(np.diff(t) > gap)[0]
    s = np.r_[0, cut + 1]; e = np.r_[cut, len(t) - 1]
    return t[s], t[e], e - s + 1


def profile(tt, ref, sign, t_lo, t_hi):
    """counts of the sorted times tt in [ref + sign*EDGES[i], ref + sign*EDGES[i+1]) per bin, and the number of events
    whose bin lies fully inside [t_lo, t_hi]."""
    nb = len(EDGES) - 1
    cnt = np.zeros(nb); expo = np.zeros(nb)
    for i in range(nb):
        a, b = EDGES[i], EDGES[i + 1]
        if sign > 0:
            lo, hi = ref + a, ref + b
        else:
            lo, hi = ref - b, ref - a
        ok = (lo >= t_lo) & (hi <= t_hi)
        cnt[i] = (np.searchsorted(tt, hi[ok]) - np.searchsorted(tt, lo[ok])).sum()
        expo[i] = ok.sum() * (b - a)
    return cnt, expo


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pid', type=int); ap.add_argument('pseudo'); ap.add_argument('out')
    a = ap.parse_args()
    t_start = time.time()
    d = sorted(glob.glob(f'{a.pseudo}/run_{a.pid}/cfg_*'))[0]
    pulser = '_tp' in d
    V, T, C, A = E.load(','.join(sorted(glob.glob(f'{d}/raw_daq_data/*.pcapng'))))
    DET = E.DET[V]
    pulsed_hit = np.zeros(len(V), bool)
    if pulser:
        for p, chips in PULSED.items():
            for v in chips:
                pulsed_hit |= (V == v) & (C == 31)
    sat = A >= 1000
    t_lo, t_hi = T[0], T[-1]
    res = dict(edges=EDGES, pulser=np.array(pulser))
    for p in (3, 2, 1):
        s, e, n = groups(T[(DET == p) & sat], 1000.0)
        for K in (3, 10):
            m = n >= K
            S, Ee = s[m], e[m]
            res[f'n_ev_K{K}_p{p}'] = np.array(m.sum())
            scopes = {'cham': T[(DET == p) & ~pulsed_hit], 'other': T[(DET > 0) & (DET != p) & ~pulsed_hit]}
            if pulser:
                tp = np.sort(T[np.isin(V, PULSED[p]) & (C == 31)])
                scopes['pulse'] = tp
            for sc, tt in scopes.items():
                res[f'mean_rate_{sc}_p{p}'] = np.array(len(tt) / (t_hi - t_lo))
                ca, xa = profile(tt, Ee, +1, t_lo, t_hi)
                cb, xb = profile(tt, S, -1, t_lo, t_hi)
                res[f'after_K{K}_{sc}_p{p}'], res[f'expo_after_K{K}_{sc}_p{p}'] = ca, xa
                res[f'before_K{K}_{sc}_p{p}'], res[f'expo_before_K{K}_{sc}_p{p}'] = cb, xb
    np.savez_compressed(f'{a.out}/sp_run{a.pid}.npz', **res)
    print(f'pid {a.pid} (pulser {pulser}): sparks K>=3', {p: int(res[f"n_ev_K3_p{p}"]) for p in (3, 2, 1)},
          f'{time.time() - t_start:.0f} s', flush=True)


if __name__ == '__main__':
    main()
