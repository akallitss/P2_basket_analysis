#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read-out efficiency from the internal test pulses of the weekend pulser runs (19-23 Nov 2026).

The *_tpA / *_tpB configurations pulse ONE channel on six chips (maps/pulsed_channels.json: {config: [[chip, ch],
...]}) once per 4096-BC frame (92.16 us, ~10.85 kHz), during cosmics. A pulse that does not arrive is a read-out loss
(the pulse amplitude is far above threshold: 100 % efficiency on the bench).

Per chunk (pseudo-run) of a pulser configuration:
  pulses    per pulsed channel: the pulse lattice (period fitted, phase from the hits), pulses received vs expected
            frames between the chunk's first and last hit, missing frames, the gap distribution (1 = no loss) and the
            run lengths of consecutive missing frames (random vs bursts)    -> pulser_channels.csv, pulser_gaps.csv
  together  frames missing on 1..6 pulsed chips at once, against the independent expectation (common cause: FEC /
            link / DAQ loss vs per-chip)                                     -> pulser_coincidence.csv
  muons     for every scintillator-confirmed muon (ON tags, --scint window) whose probe chamber holds a pulsed chip:
            were the two lattice pulses around the muon time received on that chamber's pulsed chips? Compared for
            silent muons (no live hit in the chamber, n_any == 0) and detected ones  -> pulser_muons.csv
Hits of the pulsed channel off the lattice (noise) are counted separately.

Usage: pulser_readout.py PSEUDO_DIR TAGS_DIR PULSED_JSON OUT_DIR [--pids P1,P2] [--scint 140:310] [--tol-ns 1000]
(lxplus / HTCondor: decodes the capture files of each chunk; ~2-3 GB per chunk.)
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
FRAME_NS = 4096 * 22.5          # 92160 ns
DETP = {3: range(2, 8), 2: range(8, 14), 1: range(14, 20)}
NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}


def lattice(t, P0, tol):
    """(t0, P, k of the on-lattice hits, mask of on-lattice hits) for one pulsed channel."""
    if len(t) < 20:
        return None
    # start from the NOMINAL period: the pulse comes once per 4096-BC frame, so P0 is exact. A period from the median
    # interval is biased by up to ~5 ns (noise pairs of a noisy pulsed channel, chip 8 in run_270) or ~0.2 ns (ns
    # rounding of the hit times, chip 11), and over 650 000 frames the lattice walks off the pulses. The least-squares
    # step below still refines P on the selected pulses.
    P = float(P0)
    t0 = float(t[0])
    for _ in range(3):
        k = np.round((t - t0) / P)
        r = t - t0 - k * P
        # phase = peak of the residuals (200 ns bins), not their median: noise hits spread over the whole frame pull
        # the median away from the pulses
        h, e = np.histogram(r, bins=np.arange(-P / 2, P / 2 + 200, 200))
        r0 = 0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])
        # re-centre the lattice on the pulses (t[0] may be a noise hit half a frame off: residuals would wrap at +-P/2)
        t0 += r0
        k = np.round((t - t0) / P)
        r = t - t0 - k * P
        on = np.abs(r) < tol
        if on.sum() < 10:
            return None
        A = np.vstack([np.ones(on.sum()), k[on]]).T
        c = np.linalg.lstsq(A, t[on], rcond=None)[0]
        t0, P = float(c[0]), float(c[1])
    k = np.round((t - t0) / P).astype(np.int64)
    on = np.abs(t - t0 - k * P) < tol
    return t0, P, k, on


def runs_of(missing):
    """lengths of runs of consecutive integers in a sorted array."""
    if not len(missing):
        return np.zeros(0, int)
    br = np.nonzero(np.diff(missing) != 1)[0]
    starts = np.r_[0, br + 1]
    ends = np.r_[br, len(missing) - 1]
    return ends - starts + 1


def one_chunk(pid, info, pulsed, a, tagsdir):
    from efficiency_tp import load
    d = sorted(glob.glob(f'{a.pseudo}/run_{pid}/cfg_*'))
    files = sorted(glob.glob(f'{d[0]}/raw_daq_data/*.pcapng'))
    V, T, C, A = load(','.join(files))
    tmin, tmax = float(T.min()), float(T.max())
    ch_rows, gap_rows, miss = [], [], {}
    lat = {}
    for v, c in pulsed:
        t = T[(V == v) & (C == c)].astype(float)
        L = lattice(t, FRAME_NS, a.tol_ns)
        if L is None:
            ch_rows.append(dict(pid=pid, run=info['run'], cfg=info['cfg'], chip=v, ch=c, hits=len(t), pulses=0,
                                expected=np.nan, eff=np.nan))
            continue
        t0, P, k, on = L
        kk = np.unique(k[on])
        k0, k1 = int(np.ceil((tmin - t0) / P)), int(np.floor((tmax - t0) / P))
        kk = kk[(kk >= k0) & (kk <= k1)]
        exp = k1 - k0 + 1
        allk = np.arange(k0, k1 + 1)
        m = np.setdiff1d(allk, kk, assume_unique=True)
        miss[v] = set(m.tolist())
        lat[v] = (t0, P, set(kk.tolist()))
        rl = runs_of(m)
        gaps = np.diff(kk)
        ch_rows.append(dict(pid=pid, run=info['run'], cfg=info['cfg'], chip=v, ch=c, hits=len(t), pulses=len(kk),
                            expected=exp, eff=round(len(kk) / exp, 5), missing=len(m), off_lattice_hits=int((~on).sum()),
                            dup_hits=int(on.sum() - len(np.unique(k[on]))), period_ns=round(P, 3),
                            miss_runs=len(rl), miss_run_max=int(rl.max()) if len(rl) else 0,
                            miss_in_runs_ge5=int(rl[rl >= 5].sum()) if len(rl) else 0,
                            live_s=round((tmax - tmin) * 1e-9, 1)))
        for g, n in zip(*np.unique(gaps, return_counts=True)):
            gap_rows.append(dict(pid=pid, chip=v, gap_frames=int(g), n=int(n)))
    # frames missing on several chips at once (common lattice: the first channel's; chips share the FEC pulser)
    co_rows = []
    if lat:
        ref = next(iter(lat))
        t0r, Pr, _ = lat[ref]
        k0, k1 = int(np.ceil((tmin - t0r) / Pr)), int(np.floor((tmax - t0r) / Pr))
        nfr = k1 - k0 + 1
        cnt = np.zeros(nfr, int)
        ps = []
        for v, (t0, P, got) in lat.items():
            off = int(round((t0 - t0r) / Pr))      # same physical frame, other chip
            mk = np.array(sorted(miss[v]), int) + off - k0
            mk = mk[(mk >= 0) & (mk < nfr)]
            cnt[mk] += 1
            ps.append(len(mk) / nfr)
        obs = np.bincount(cnt, minlength=len(lat) + 1)
        # independent expectation: Poisson-binomial over the chips' loss probabilities
        pb = np.zeros(len(ps) + 1); pb[0] = 1
        for p in ps:
            pb = np.r_[pb * (1 - p), 0] + np.r_[0, pb * p]
        for n in range(len(obs)):
            co_rows.append(dict(pid=pid, run=info['run'], cfg=info['cfg'], n_chips_missing=n, frames=int(obs[n]),
                                expected_indep=round(float(pb[n] * nfr), 2)))
    # muons: pulses around the muon time on the pulsed chips of the probe chamber
    mu_rows = []
    fz = f'{tagsdir}/tags_run{pid}.npz'
    if lat and os.path.exists(fz):
        Z = np.load(fz)
        lo, hi = map(float, a.scint.split(':'))
        for p, chips in DETP.items():
            pc = [v for v in lat if v in chips]
            if not pc or f'p{p}_t_tag' not in Z.files:
                continue
            dts = Z[f'p{p}_dt_scint']
            sel = (dts >= lo) & (dts <= hi) & Z[f'p{p}_on'].astype(bool)
            tt = Z[f'p{p}_t_tag'][sel].astype(float) + float(Z[f'p{p}_latency'])
            silent = Z[f'p{p}_n_any'][sel] == 0
            for v in pc:
                t0, P, got = lat[v]
                kb = np.floor((tt - t0) / P).astype(np.int64)
                rec = np.array([(k in got) and (k + 1 in got) for k in kb])
                for s in (True, False):
                    q = silent == s
                    mu_rows.append(dict(pid=pid, run=info['run'], cfg=info['cfg'], det=NAME[p], chip=v,
                                        muons='silent' if s else 'detected', n=int(q.sum()),
                                        both_pulses_received=int(rec[q].sum())))
    return ch_rows, gap_rows, co_rows, mu_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('pseudo'); ap.add_argument('tags'); ap.add_argument('pulsed'); ap.add_argument('out')
    ap.add_argument('--pids', default='')
    ap.add_argument('--scint', default='140:310')
    ap.add_argument('--tol-ns', type=float, default=1000.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    chunks = {int(k): v for k, v in json.load(open(f'{a.pseudo}/chunks.json')).items()}
    pulsed = {k: [tuple(x) for x in v] for k, v in json.load(open(a.pulsed)).items() if not k.startswith('_')}
    want = {int(x) for x in a.pids.split(',') if x}
    todo = []
    for pid, info in sorted(chunks.items()):
        name = re.sub(r'^cfg_', '', info['cfg'])
        if name in pulsed and (not want or pid in want):
            todo.append((pid, info, pulsed[name]))
    print(f'{len(todo)} pulser chunks', flush=True)
    R = [[], [], [], []]
    for pid, info, pl in todo:
        out = one_chunk(pid, info, pl, a, a.tags)
        for i, x in enumerate(out):
            R[i] += x
        c = pd.DataFrame(out[0])
        print(f'{pid} {info["cfg"]}: ' + ', '.join(f'{r.chip}:{r.ch} {r.eff:.4f}' for r in c.itertuples()), flush=True)
        # write as we go (a long job keeps partial results)
        for name, rows in zip(('channels', 'gaps', 'coincidence', 'muons'), R):
            pd.DataFrame(rows).to_csv(f'{a.out}/pulser_{name}.csv', index=False)
    if not todo:
        return
    ch = pd.DataFrame(R[0])
    g = ch.groupby(['cfg', 'chip', 'ch'])[['pulses', 'expected', 'missing', 'miss_in_runs_ge5']].sum()
    g['eff'] = (g.pulses / g.expected).round(5)
    print('\n== per configuration and pulsed channel (all chunks)')
    print(g.to_string())
    if R[3]:
        mu = pd.DataFrame(R[3]).groupby(['cfg', 'det', 'chip', 'muons'])[['n', 'both_pulses_received']].sum()
        mu['frac_received'] = (mu.both_pulses_received / mu.n).round(4)
        print('\n== muons on chambers with a pulsed chip: both lattice pulses around the muon received?')
        print(mu.to_string())


if __name__ == '__main__':
    main()
