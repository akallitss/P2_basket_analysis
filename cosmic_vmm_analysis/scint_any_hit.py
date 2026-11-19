#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reconstruction-free chamber efficiency with the scintillator coincidence as the only reference.

For every hit of the scintillator coincidence channel (H0 VMM 0 ch 40) and every chamber: is there ANY hit of that
chamber (any channel, live or masked-but-readout, any ADC) in a window after the scintillator? No pad map, no
clusters, no tracks. Random hits: the same window shifted by --offset ns (no muon there).

  single   P(chamber hit | scint), random-corrected: (p_on - p_off) / (1 - p_off). Contains the geometric acceptance
           (muons through the paddles that miss the chamber's read-out area).
  probe    P(probe hit | scint and the two OTHER chambers hit in their on-time windows), random-corrected with the
           probe's off-time rate on the same triggers. Removes most of the acceptance for P2_MID (bracketed); for the
           outer chambers some acceptance remains (extrapolation).
Windows (ns after the scintillator hit): several widths around the muon time, the widest first-to-last.

Usage: scint_any_hit.py RUNS_DIR RUN [RUN ...] [--files 15] [--offset 4000]
"""
import argparse
import glob

import numpy as np

from coincidence_check import read_times

DET = {'P2_IN': range(2, 8), 'P2_MID': range(8, 14), 'P2_OUT': range(14, 20)}
WINDOWS = {'140:310 (now)': (140, 310), '100:400': (100, 400), '0:600': (0, 600), '-200:1000': (-200, 1000),
           '-500:2000': (-500, 2000)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs_dir'); ap.add_argument('runs', nargs='+')
    ap.add_argument('--files', type=int, default=15)
    ap.add_argument('--offset', type=float, default=4000)
    a = ap.parse_args()
    for run in a.runs:
        files = sorted(glob.glob(f'{a.runs_dir}/run_{run}/cfg_*/raw_daq_data/*.pcapng'))
        if a.files and len(files) > a.files:
            step = len(files) // a.files
            files = files[1::step][:a.files]
        cfg = glob.glob(f'{a.runs_dir}/run_{run}/cfg_*')[0].split('cfg_')[-1]
        on = {w: [] for w in WINDOWS}; off = {w: [] for w in WINDOWS}
        dur = 0.0; nsc = 0
        dts = {d: [] for d in DET}
        for f in files:
            V, T, C = [], [], []
            for v, t, c in read_times(f, {}):
                V.append(v); T.append(t); C.append(c)
            V, T, C = np.concatenate(V), np.concatenate(T), np.concatenate(C)
            dur += (T.max() - T.min()) * 1e-9
            sc = np.sort(T[(V == 0) & (C == 40)])
            sc = sc[(sc > T.min() + 1e5) & (sc < T.max() - 1e5)]
            nsc += len(sc)
            ht = {d: np.sort(T[np.isin(V, list(r))]) for d, r in DET.items()}
            for w, (lo, hi) in WINDOWS.items():
                on[w].append(np.array([np.searchsorted(ht[d], sc + hi) > np.searchsorted(ht[d], sc + lo) for d in DET]))
                off[w].append(np.array([np.searchsorted(ht[d], sc + hi + a.offset) >
                                        np.searchsorted(ht[d], sc + lo + a.offset) for d in DET]))
            for d in DET:      # first chamber hit after the scintillator, for the timing spectrum
                k = np.searchsorted(ht[d], sc - 500)
                ok = k < len(ht[d])
                dts[d].append(ht[d][k[ok]] - sc[ok])
        print(f'\nrun_{run} {cfg}: {len(files)} files, {dur / 60:.1f} min, scintillator {nsc} ({nsc / dur:.2f} Hz)')
        for d in DET:
            x = np.concatenate(dts[d]); x = x[(x > -500) & (x < 2000)]
            print(f'  {d}: first hit after scint, quantiles 5/25/50/75/95 % (ns):',
                  np.round(np.percentile(x, [5, 25, 50, 75, 95])).astype(int).tolist())
        print(f'  {"window (ns)":16s} | ' + ' | '.join(f'{d:^25s}' for d in DET))
        print(f'  {"":16s} | ' + ' | '.join(f'{"single":>8s} {"probe":>8s} {"rand":>6s}' for _ in DET))
        for w in WINDOWS:
            F, G = np.concatenate(on[w], axis=1), np.concatenate(off[w], axis=1)
            cells = []
            for i, d in enumerate(DET):
                p_on, p_off = F[i].mean(), G[i].mean()
                single = (p_on - p_off) / (1 - p_off)
                o = [j for j in range(3) if j != i]
                m = F[o[0]] & F[o[1]]
                q_on, q_off = F[i][m].mean(), G[i][m].mean()
                probe = (q_on - q_off) / (1 - q_off)
                cells.append(f'{single:8.3f} {probe:8.3f} {p_off:6.3f}')
            print(f'  {w:16s} | ' + ' | '.join(cells) + f'   (n probe: {int(m.sum())})')


if __name__ == '__main__':
    main()
