#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""How many scintillator-confirmed muons does the tag-and-probe throw away because a tag chamber's cluster is not
"clean" (handoff 20261119d, finding 11 / next step 4)?

efficiency_tp.chamber_clusters groups the hits of a chamber in TIME (gap --gap ns) and keeps a group as a tag only if
it is clean: 1..--max-tag hits, all on live pads, within --cluster-r mm of their mean, not inside a stack burst
(> --burst hits in the whole-stack group). A muon whose tag group fails is not counted as a miss, it just leaves the
sample. Here every scintillator coincidence (--scint-ch) is taken as a muon and, per chamber, the time groups with a hit
in [t_scint + lo, t_scint + hi] are classified:

  clean      a clean group overlaps the window (it would be a tag)
  rejected   only non-clean groups: reasons (not exclusive) burst / dead (a hit on a dead/masked/unmapped pad) /
             npads (> max_tag hits) / wide (spread > cluster_r); and the number of SPATIAL sub-clusters of the group
             (single linkage within --cluster-r mm on live mapped pads): 1 = one big blob, >= 2 = muon + something else
  none       no hit of the chamber in the window
The same at t_scint + --offset (no muon there) gives the random rates; "corr" columns are on - off.

Per probe chamber (tags = the other two): muons with both tag chambers active, split into USED (both clean) and LOST
(>= 1 rejected), and P(probe any hit in its window) for both samples - if the lost muons had a different probe
efficiency, dropping them would bias the result.

Usage: tag_rejects.py EFF_JSON [EFF_JSON ...] --out OUT.json   (EFF_JSON = efficiency_tp.py json of a (pseudo-)run:
       files and all mapping / cut parameters are taken from its 'params', so the classification uses the same cuts)
"""
import argparse
import json
import sys
from collections import Counter

import numpy as np

import orientation_check as O
from efficiency_tp import DET, load, split_groups
from straight_track_rate import parse_veto

NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
TAGS = {2: (3, 1), 1: (3, 2), 3: (2, 1)}


def n_subclusters(x, y, r):
    """single-linkage clusters of points within r mm."""
    n = len(x)
    lab = -np.ones(n, int)
    k = 0
    for i in range(n):
        if lab[i] >= 0:
            continue
        stack, lab[i] = [i], k
        while stack:
            j = stack.pop()
            nb = np.nonzero((lab < 0) & (np.hypot(x - x[j], y - y[j]) <= r))[0]
            lab[nb] = k
            stack += nb.tolist()
        k += 1
    return k


def groups_of(d, V, T, C, X, Y, live, bad_hit, burst_hit, p):
    """all time groups of chamber d: start, end, clean flag, reason flags, spatial sub-clusters (lazy)."""
    m = DET[V] == d
    idx = np.nonzero(m)[0]
    t, v, c = T[m], V[m], C[m]
    x, y = X[v, c], Y[v, c]
    ok = live[v, c]
    g = split_groups(t, p['gap'])
    ng = int(g[-1]) + 1 if len(g) else 0
    n_all = np.bincount(g, minlength=ng)
    n_live = np.bincount(g, weights=ok, minlength=ng)
    sx = np.bincount(g, weights=np.where(ok, x, 0), minlength=ng)
    sy = np.bincount(g, weights=np.where(ok, y, 0), minlength=ng)
    with np.errstate(invalid='ignore', divide='ignore'):
        mx, my = sx / n_live, sy / n_live
    r = np.hypot(x - mx[g], y - my[g])
    spread = np.zeros(ng)
    np.maximum.at(spread, g, np.where(ok, r, 0))
    burst = np.bincount(g, weights=burst_hit[m], minlength=ng) > 0
    veto = np.bincount(g, weights=bad_hit[m], minlength=ng) > 0
    t0 = np.full(ng, np.inf); np.minimum.at(t0, g, t)
    t1 = np.full(ng, -np.inf); np.maximum.at(t1, g, t)
    R = dict(burst=burst | veto, dead=n_live < n_all, npads=n_all > p['max_tag'], wide=spread > p['cluster_r'],
             nolive=n_live < 1)
    clean = (n_live >= 1) & (n_all <= p['max_tag']) & (n_live == n_all) & (spread <= p['cluster_r']) & ~burst & ~veto
    starts = np.searchsorted(g, np.arange(ng))
    return dict(t0=t0, t1=t1, clean=clean, R=R, n_all=n_all, g=g, x=x, y=y, ok=ok, starts=starts,
                ends=np.r_[starts[1:], len(g)])


def classify(G, ts, lo, hi, r_sub):
    """state per scintillator time: 0 none, 1 clean, 2 rejected; reasons and sub-clusters of the rejected."""
    n = len(ts)
    state = np.zeros(n, int)
    reasons, subs, npads = [None] * n, np.zeros(n, int), np.zeros(n, int)
    # groups overlapping [ts+lo, ts+hi]: t0 <= ts+hi and t1 >= ts+lo; groups are time-ordered
    a = np.searchsorted(G['t1'], ts + lo)            # first group ending after the window start
    b = np.searchsorted(G['t0'], ts + hi, side='right')  # groups starting before the window end
    for k in np.nonzero(b > a)[0]:
        gs = np.arange(a[k], b[k])
        if G['clean'][gs].any():
            state[k] = 1
            continue
        state[k] = 2
        q = gs[np.argmax(G['n_all'][gs])]           # the biggest overlapping group
        reasons[k] = tuple(r for r in ('burst', 'dead', 'npads', 'wide', 'nolive') if G['R'][r][q])
        s0, s1 = G['starts'][q], G['ends'][q]
        ok = G['ok'][s0:s1]
        subs[k] = n_subclusters(G['x'][s0:s1][ok], G['y'][s0:s1][ok], r_sub) if ok.any() else 0
        npads[k] = G['n_all'][q]
    return state, reasons, subs, npads


def summarise(st, rs, sb, npd, n):
    out = {'n_scint': int(n), 'clean': int((st == 1).sum()), 'rejected': int((st == 2).sum()),
           'none': int((st == 0).sum())}
    rej = st == 2
    out['reasons'] = dict(Counter('+'.join(r) for r, x in zip(rs, rej) if x and r is not None))
    out['reason_any'] = {k: int(sum(1 for r, x in zip(rs, rej) if x and r and k in r))
                         for k in ('burst', 'dead', 'npads', 'wide', 'nolive')}
    out['subclusters'] = {str(k): int(v) for k, v in Counter(sb[rej].tolist()).items()}
    out['npads_rejected'] = {str(k): int(v) for k, v in Counter(np.minimum(npd[rej], 40).tolist()).items()}
    return out


def one_run(path, a):
    J = json.load(open(path))
    p = J['params']
    X, Y = O.Mapper(p['conn_shift']).xy_table(O.parse_orders(p['order']))
    V, T, C, A = load(p['files'])
    counts = np.zeros((32, 64), int)
    np.add.at(counts, (V, C), 1)
    live = np.isfinite(X) & (counts >= p['min_hits'])
    for v in {int(v) for v in str(p['dead_vmm']).split(',') if v}:
        live[v] = False
    veto = parse_veto(p['veto'])
    vmask = np.zeros((32, 64), bool)
    for v, c in veto:
        live[v, c] = False
        vmask[v, c] = True
    g = split_groups(T, p['gap'])
    burst_hit = np.bincount(g)[g] > p['burst']
    bad_hit = vmask[V, C]
    sv, sc = map(int, p['scint_ch'].split(':'))
    ts = np.sort(T[(V == sv) & (C == sc)])
    ts = ts[(ts > T[0] + 1e4) & (ts < T[-1] - 1e4)]
    res = {'title': p['title'], 'span_s': float(T[-1] - T[0]) * 1e-9, 'chambers': {}, 'probes': {}}
    S = {}
    for d in (3, 2, 1):
        G = groups_of(d, V, T, C, X, Y, live, bad_hit, burst_hit, p)
        S[d] = {w: classify(G, ts + off, a.lo, a.hi, p['cluster_r']) for w, off in (('on', 0.0), ('off', a.offset))}
        res['chambers'][NAME[d]] = {w: summarise(*S[d][w], len(ts)) for w in ('on', 'off')}
    for pr, (ta, tb) in TAGS.items():
        out = {}
        for w in ('on', 'off'):
            sa, sb_, sp = S[ta][w][0], S[tb][w][0], S[pr][w][0]
            act = (sa > 0) & (sb_ > 0)
            used = act & (sa == 1) & (sb_ == 1)
            lost = act & ~used
            out[w] = {'both_active': int(act.sum()), 'used': int(used.sum()), 'lost': int(lost.sum()),
                      'lost_one_clean': int((lost & ((sa == 1) | (sb_ == 1))).sum()),
                      'probe_hit_used': int((used & (sp > 0)).sum()), 'probe_hit_lost': int((lost & (sp > 0)).sum())}
        res['probes'][NAME[pr]] = out
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('jsons', nargs='+')
    ap.add_argument('--out', required=True)
    ap.add_argument('--lo', type=float, default=140)
    ap.add_argument('--hi', type=float, default=310)
    ap.add_argument('--offset', type=float, default=4000)
    a = ap.parse_args()
    allres = {}
    for j in a.jsons:
        r = one_run(j, a)
        allres[j] = r
        print(r['title'], {k: (v['on']['clean'], v['on']['rejected'], v['on']['none']) for k, v in r['chambers'].items()},
              flush=True)
    json.dump(allres, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    sys.exit(main())
