#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""det5, pad by pad: DREAM cosmic bench (long run 9 Sep, 420/820 V) vs VMM lab stand at P2_MID (25 Nov 2026).

Are the dead / weak / silent pads the same in the two set-ups? Connectors c3-c5 only (the stand's read-out area).
  DREAM: stage-02 hit counts per pad (dead = 0 hits, weak = < 30 % of its half-connector median) and the per-pad
         efficiency of dream_pad_export.py (within MATCH_R, any hit). The even DREAM channels of FEU 6 / DREAM
         connector 7 are flagged BENCH: that pin row is dead on the bench for every chamber plugged there since
         Sep (det5 9 Sep on its c4, det2 24 Nov on its c5; alive in July for det1-det4).
  VMM:   silent_tests.py per-pad rows [vmm, ch, x, y, n, silent] of the REF group, P2_MID, before the read-out swap
         (night 23-24 Nov, det5 read by VMM 8-13) and after (24-25 Nov, VMM 14-19); stand-dead channels (no hit in
         a whole run: dead_or_masked lists) placed with orientation_check.Mapper(-1) positions (VMM_XY json).
Pads are matched by position (board frame, < 2 mm).
Usage: det5_dream_vmm_pads.py --stage02 CSV --dream CSV --night JSON --swap JSON --xy JSON [--out JSON]
@author: Alexandra Kallitsopoulou
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

STAND_DEAD = {'night': [(9, 0), (9, 35), (11, 27), (11, 52)],          # dead_or_masked det5-only (not on det2)
              'swap': []}
MASKED_ALWAYS = [(8, 26), (10, 22), (10, 45), (13, 25)]                # in the det2 weekend list too (channels)


def vmm_pads(st, vmin):
    rows = st['REF']['P2_MID']['pads']
    d = pd.DataFrame(rows, columns=['vmm', 'ch', 'x', 'y', 'n', 'silent'][:len(rows[0])])
    return d[(d.vmm >= vmin) & (d.vmm < vmin + 6)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for k in ('stage02', 'dream', 'night', 'swap', 'xy'):
        ap.add_argument('--' + k, required=True)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    m = pd.read_csv(a.stage02)
    m = m[m.mapped].drop_duplicates('channel_id')
    m['rel'] = m.n_hits / m.groupby(['connector_N', 'half']).n_hits.transform('median')
    m['bench'] = (m.feu == 6) & (m.dream_conn == 7) & (m.channel % 2 == 0)
    e = pd.read_csv(a.dream).rename(columns={'within': 'eff_within', 'any': 'eff_any', 'n': 'rays'})
    D = m.merge(e[['channel_id', 'rays', 'eff_within', 'eff_any', 'pad_cx', 'pad_cy']], on='channel_id',
                suffixes=('_02', ''))
    D = D[D.connector_N.between(3, 5)].reset_index(drop=True)
    D['dream'] = np.where(D.bench, 'BENCH', np.where(D.n_hits == 0, 'dead', np.where(D.rel < 0.3, 'weak', 'ok')))
    tree = cKDTree(D[['pad_cx', 'pad_cy']].to_numpy())

    def attach(v):
        dist, i = tree.query(v[['x', 'y']].to_numpy())
        v = v.assign(dist=dist, idx=i)
        return v[v.dist < 2.0]

    N = attach(vmm_pads(json.load(open(a.night)), 8))
    S = attach(vmm_pads(json.load(open(a.swap)), 14))
    D['vmm_night'] = ''
    D['vmm_swap'] = ''
    for lab, V in (('night', N), ('swap', S)):
        D.loc[V.idx, f'n_{lab}'] = V.n.to_numpy()
        D.loc[V.idx, f'silent_{lab}'] = V.silent.to_numpy()
        D.loc[V.idx, f'vmm_{lab}'] = [f'{int(r.vmm)}:{int(r.ch)}' for r in V.itertuples()]
    xy = json.load(open(a.xy))
    print(f'det5 c3-c5: {len(D)} pads; matched VMM pads night {len(N)}, swap {len(S)}')
    print(f'DREAM status: {D.dream.value_counts().to_dict()}')

    print('\n1. Stand-dead det5 channels (no hit in a whole night), and what DREAM saw on the same pad:')
    for v, c in STAND_DEAD['night'] + MASKED_ALWAYS:
        x, y = xy['X'][v][c], xy['Y'][v][c]
        if x < -900:
            print(f'   VMM {v} ch {c}: not mapped')
            continue
        dist, i = tree.query([x, y])
        r = D.iloc[i]
        tag = 'det5 only' if (v, c) in STAND_DEAD['night'] else 'masked/dead on every chamber'
        print(f'   VMM {v:2d} ch {c:2d} ({tag}) at ({x:6.1f},{y:6.1f}) -> DREAM ch_id {r.channel_id} c{r.connector_N} '
              f'{r.half} strip {r.strip} [{r.dream}] hits {r.n_hits} (rel {r.rel:.2f}), bench eff within '
              f'{r.eff_within:.2f} any {r.eff_any:.2f} (rays {r.rays:.0f}); after swap VMM {r.vmm_swap or "-"} '
              f'silent {r.get("silent_swap", np.nan):.2f}')

    print('\n2. DREAM dead / weak pads (not the bench pin row), and the stand on the same pad:')
    for r in D[D.dream.isin(['dead', 'weak'])].itertuples():
        print(f'   ch_id {r.channel_id} c{r.connector_N} {r.half} strip {r.strip} [{r.dream}] hits {r.n_hits} '
              f'-> stand night {r.vmm_night or "-"} silent {r.silent_night:.2f} (n {r.n_night:.0f}), '
              f'swap {r.vmm_swap or "-"} silent {r.silent_swap:.2f} (n {r.n_swap:.0f})')

    b = D[D.dream == 'BENCH']
    print(f'\n3. Bench pin row (FEU 6 / conn 7 even, {len(b)} pads of det5 c{sorted(b.connector_N.unique())} '
          f'{sorted(b.half.unique())}): on the stand these pads are silent '
          f'{np.nansum(b.silent_night * b.n_night) / np.nansum(b.n_night):.3f} (night) / '
          f'{np.nansum(b.silent_swap * b.n_swap) / np.nansum(b.n_swap):.3f} (swap); the other pads of that half '
          f'{np.nansum(D[(D.connector_N == 4) & (D.half == "bot") & ~D.bench].silent_night * D[(D.connector_N == 4) & (D.half == "bot") & ~D.bench].n_night) / np.nansum(D[(D.connector_N == 4) & (D.half == "bot") & ~D.bench].n_night):.3f} (night)')

    G = D[(D.dream == 'ok')]
    print('\n4. Per half-connector, pads DREAM sees as ok (bench any / within vs stand silent before / after swap):')
    rows = []
    for (c, h), g in G.groupby(['connector_N', 'half']):
        w = lambda col, n: np.nansum(g[col] * g[n]) / np.nansum(g[n])
        r = dict(conn=int(c), half=h, pads=len(g), dream_any=float(np.average(g.eff_any, weights=g.rays)),
                 dream_within=float(np.average(g.eff_within, weights=g.rays)), silent_night=float(w('silent_night', 'n_night')),
                 silent_swap=float(w('silent_swap', 'n_swap')))
        rows.append(r)
        print(f'   c{c} {h}: {len(g):3d} pads  DREAM any {r["dream_any"]:.3f} within {r["dream_within"]:.3f}  |  '
              f'stand silent {r["silent_night"]:.3f} -> {r["silent_swap"]:.3f}')

    print('\n5. Pad-to-pad correlation on DREAM-ok pads (>= 30 VMM muons in both periods, >= 30 bench rays):')
    H = G[(G.n_night >= 30) & (G.n_swap >= 30) & (G.rays >= 30)]
    c_ns = np.corrcoef(H.silent_night, H.silent_swap)[0, 1]
    c_nd = np.corrcoef(H.silent_night, 1 - H.eff_any)[0, 1]
    c_sd = np.corrcoef(H.silent_swap, 1 - H.eff_any)[0, 1]
    print(f'   {len(H)} pads: corr(stand night, stand swap) = {c_ns:.2f}; corr(stand night, DREAM 1-any) = {c_nd:.2f}; '
          f'corr(stand swap, DREAM 1-any) = {c_sd:.2f}')
    # expected scatter of a pad's silent share from counting alone, to judge the correlations
    print(f'   median pad: n night {H.n_night.median():.0f}, n swap {H.n_swap.median():.0f}; spread of silent night '
          f'{H.silent_night.std():.3f}, binomial expectation {np.sqrt(np.mean(H.silent_night * (1 - H.silent_night) / H.n_night)):.3f}')
    top = H.assign(s=(H.silent_night * H.n_night + H.silent_swap * H.n_swap) / (H.n_night + H.n_swap)).nlargest(12, 's')
    print('\n6. The 12 most silent det5 pads on the stand (both periods together), with the bench:')
    for r in top.itertuples():
        print(f'   ch_id {r.channel_id} c{r.connector_N} {r.half} strip {r.strip:3d} ({r.pad_cx:6.1f},{r.pad_cy:6.1f})  '
              f'stand {r.s:.2f} (night {r.silent_night:.2f}/{r.n_night:.0f}, swap {r.silent_swap:.2f}/{r.n_swap:.0f})  '
              f'DREAM any {r.eff_any:.2f} within {r.eff_within:.2f} (rays {r.rays:.0f}, hits {r.n_hits})')
    if a.out:
        json.dump(dict(halves=rows, corr=dict(n=len(H), night_swap=c_ns, night_dream=c_nd, swap_dream=c_sd)),
                  open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
