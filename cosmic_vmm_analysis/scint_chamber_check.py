#!/usr/bin/env python3
"""Quick chamber-level check of scintillator-triggered cosmics (no pad map, no tracking).

For every hit of the scintillator coincidence channel (H0 VMM 0 ch 40 by default) it asks whether each chamber has
ANY hit 140-310 ns later (on-time window) and 4140-4310 ns later (off-time window = random hits). Prints per run:
  - P(hit | scint): depends on trigger purity and on how many triggered muons cross the read-out area;
  - P(probe | scint & the other two chambers): closer to a chamber efficiency, but still contains the geometric
    acceptance of the probe's read-out region. For the outer chambers (P2_IN, P2_OUT) a muon through the other two
    can leave the read-out footprint (extrapolation), for P2_MID it is bracketed (interpolation), so the outer
    numbers are efficiency x acceptance;
  - 3-chamber triples per minute.
This is NOT the single-channel efficiency (channel_efficiency.py, CHANNEL_EFFICIENCY.md); it saturates with
threshold (any pad of the chamber counts) and cannot resolve a few mV of threshold. Use it to check that a run has
scintillator tags and to compare trigger conditions.

    python3 scint_chamber_check.py --runs-dir /local/p2/p2data/cosmics_3det_202610/runs --files 15 231 232 234 235
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coincidence_check import read_times  # noqa: E402

DET = {'IN': range(2, 8), 'MID': range(8, 14), 'OUT': range(14, 20)}


def check_run(files, scint, windows):
    dur, nsc, acc = 0.0, 0, {w: [] for w in windows}
    for f in files:
        V, T, C = [], [], []
        for v, t, c in read_times(f, {}):
            V.append(v); T.append(t); C.append(c)
        if not V:
            continue
        V, T, C = map(np.concatenate, (V, T, C))
        dur += (T.max() - T.min()) * 1e-9
        sc = np.sort(T[(V == scint[0]) & (C == scint[1])])
        nsc += len(sc)
        for w, (lo, hi) in windows.items():
            fl = []
            for r in DET.values():
                ht = np.sort(T[np.isin(V, list(r))])
                fl.append(np.searchsorted(ht, sc + hi) > np.searchsorted(ht, sc + lo))
            acc[w].append(np.array(fl).reshape(3, -1))
    return dur, nsc, {w: np.concatenate(a, axis=1) if a else np.zeros((3, 0), bool) for w, a in acc.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--runs-dir', required=True, help='directory holding run_N/cfg_*/raw_daq_data/*.pcapng')
    ap.add_argument('--files', type=int, default=15, help='capture files per run, spread over the run (0 = all)')
    ap.add_argument('--scint-ch', default='0:40')
    ap.add_argument('--window', default='140:310', help='on-time window (ns after the scintillator hit)')
    ap.add_argument('--offset', type=float, default=4000, help='off-time window shift (ns)')
    a = ap.parse_args()
    scint = tuple(map(int, a.scint_ch.split(':')))
    lo, hi = map(float, a.window.split(':'))
    windows = {'on': (lo, hi), 'off': (lo + a.offset, hi + a.offset)}
    for run in a.runs:
        files = sorted(glob.glob(f'{a.runs_dir}/run_{run}/cfg_*/raw_daq_data/*.pcapng'))
        if a.files and len(files) > a.files:
            step = len(files) // a.files
            files = files[1::step][:a.files]      # skip the first file (run start)
        dur, nsc, acc = check_run(files, scint, windows)
        if dur == 0:
            print(f'run_{run}: no data'); continue
        print(f'run_{run}: {len(files)} files, {dur / 60:.1f} min, scint {nsc / dur:.2f} Hz')
        for w, F in acc.items():
            single = ' '.join(f'{d} {F[i].mean():.3f}' for i, d in enumerate(DET)) if F.shape[1] else '-'
            tp = []
            for i, d in enumerate(DET):
                o = [j for j in range(3) if j != i]
                m = F[o[0]] & F[o[1]]
                n, k = int(m.sum()), int(F[i][m].sum())
                p = k / n if n else float('nan')
                err = np.sqrt(p * (1 - p) / n) if n else float('nan')
                tp.append(f'{d} {k}/{n} = {p:.3f} +- {err:.3f}')
            print(f'  {w:3s} P(hit|scint) {single} | triples {F.all(0).sum() / dur * 60:.1f}/min'
                  f' | P(probe|scint & other two) {"  ".join(tp)}')


if __name__ == '__main__':
    main()
