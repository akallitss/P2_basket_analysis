#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Keep only the hits of small multi-chamber time groups, for the orientation / mapping studies.

Hits (over-threshold) are grouped in time (gap > --gap ns starts a new group); a group is kept when it has at
most --burst hits and touches at least --min-chambers chambers. Output .npz: vmm, ch, adc, t_ns (float64), group id
(unique across the files of one call), chamber (3 = P2_IN, 2 = P2_MID, 1 = P2_OUT), and the span in s.

Usage: dump_candidates.py f1.pcapng[,f2,...] OUT.npz [--gap 100] [--burst 36] [--min-chambers 2]
"""
import argparse
import sys

import numpy as np

from coincidence_check import read_times, CHAMBER


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('files')
    ap.add_argument('out_npz')
    ap.add_argument('--gap', type=float, default=100)
    ap.add_argument('--burst', type=int, default=36)
    ap.add_argument('--min-chambers', type=int, default=2)
    a = ap.parse_args()
    state = {}
    cols = [[], [], [], []]
    for p in a.files.split(','):
        for rec in read_times(p, state, with_adc=True):
            for k in range(4):
                cols[k].append(rec[k])
    V, T, C, A = (np.concatenate(c) for c in cols)
    o = np.argsort(T, kind='stable')
    V, T, C, A = V[o], T[o], C[o], A[o]
    d = np.array([CHAMBER.get(int(v), 0) for v in V], dtype=np.int8)
    ev = np.concatenate([[0], np.cumsum(np.diff(T) > a.gap)])
    ne = int(ev[-1]) + 1
    n = np.bincount(ev, minlength=ne)
    nch = sum((np.bincount(ev[d == c], minlength=ne) > 0).astype(int) for c in (1, 2, 3))
    keep = ((n <= a.burst) & (nch >= a.min_chambers))[ev]
    np.savez_compressed(a.out_npz, vmm=V[keep].astype(np.int8), ch=C[keep].astype(np.int8), adc=A[keep].astype(np.int16),
                        t=T[keep], ev=ev[keep], det=d[keep], span_s=(T[-1] - T[0]) / 1e9)
    print(f'{len(T)} hits, kept {int(keep.sum())} in {len(np.unique(ev[keep]))} groups')


if __name__ == '__main__':
    sys.exit(main())
