#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Why are the silent muons silent? Three tests on the per-tag dumps of the weekend runs (20 Nov 2026).

Tags: scintillator-confirmed (dt_scint 140-310 ns), weight on = 1, off = -r_scale (random tag pairs subtracted), the
dump's fiducial (nearest pad live within 9 mm) plus the pointing edge margin of pointing_resolution.py (--pointing).
Silent = no live probe hit within +-win of tag time + latency (n_any == 0); 'dead1us' = in addition no hit at all on
any chip of the probe chamber within +-1 us (chip_n1us). Found = a live hit within the probe's search radius in the
window (d_win <= probe_r).

  1. ANGLE: silent fraction vs tan(theta) of the tag pair. A hole in the amplification (pillar, dead area) kills
     mostly near-vertical tracks: with a 4 mm drift gap a track spreads its charge over 4 tan(theta) mm on the mesh,
     so a 0.8 mm pillar can only swallow it whole when tan(theta) << 0.2; a signal-size effect (Landau tail under
     threshold) falls only gently (path length 1 / cos(theta): +12 % at tan(theta) = 0.5).
  2. WHERE: silent fraction per predicted pad (n >= --min-pad), its spread against the binomial expectation
     (chi2 / ndf: 1 = the same loss everywhere), the worst pads, the fraction within --big-r mm of the five big
     6.15 mm support pillars of mask V1 (Saclay bulk, all three stand chambers), and vs the distance to the nearest
     dead / masked pad.
  3. TAG PULSE HEIGHT: silent fraction vs the smaller of the two tags' largest ADC (tag_adc). If low-ADC tags are
     more often followed by a silent probe, some tags are not muons that crossed the probe (noise in time with the
     scintillator, or tracks that did not go straight through).

Usage: silent_tests.py --tags TAGS_DIR --chunks pseudo/chunks.json --pointing pointing.json --out OUT.json
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np

from efficiency_tp import DET

NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
# big support pillars of mask V1 (P2_BASKET-Mask_M2_V1.gbr, pad frame, mm), diameter 6.15 mm
BIG = np.array([[206.40, 118.98], [312.39, 180.36], [411.36, 237.50], [339.74, 331.97], [457.36, 128.25]])
SLOPE = [0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 10]
TADC = [0, 100, 150, 200, 300, 400, 600, 900, 2000]
DDEAD = [0, 12, 18, 25, 40, 80, 1e9]
_ang = np.linspace(0, 2 * np.pi, 16, endpoint=False)
RING = np.r_[np.c_[np.cos(_ang), np.sin(_ang)], 0.5 * np.c_[np.cos(_ang), np.sin(_ang)]]


def wfrac(w, m):
    """weighted fraction of m, its error (binomial on the on-time count), sum of weights."""
    s = w.sum()
    if s <= 0:
        return None, None, 0.0
    f = float((w * m).sum() / s)
    n_on = max((w > 0).sum(), 1)
    return round(f, 4), round(float(np.sqrt(max(f * (1 - f), 0) / n_on)), 4), round(float(s), 1)


def binned(x, edges, w, m):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (x >= lo) & (x < hi)
        f, e, n = wfrac(w[k], m[k])
        out.append({'lo': lo, 'hi': hi, 'f': f, 'err': e, 'n': n})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--tags', required=True)
    ap.add_argument('--chunks', required=True)
    ap.add_argument('--pointing', required=True)
    ap.add_argument('--fid-r', type=float, default=9.0)
    ap.add_argument('--min-pad', type=float, default=40)
    ap.add_argument('--big-r', type=float, default=8.0)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    import orientation_check as O
    chunks = {int(k): v for k, v in json.load(open(a.chunks)).items()}
    PT = json.load(open(a.pointing))['probes']
    acc = {}           # group -> probe -> list of column arrays
    pads_of = {}
    for fz in sorted(glob.glob(f'{a.tags}/tags_run*.npz')):
        pid = int(re.search(r'run(\d+)\.npz', fz)[1])
        if pid not in chunks or '_tp' in chunks[pid]['cfg']:
            continue
        ej = fz.replace('tags_run', 'eff_run').replace('.npz', '.json')
        if not os.path.exists(ej):
            continue
        J = json.load(open(ej))
        p_ = J['params']
        if not pads_of:
            X, Y = O.Mapper(p_['conn_shift']).xy_table(O.parse_orders(p_['order']))
            vv, cc = np.nonzero(np.isfinite(X))
            for d in (3, 2, 1):
                k = DET[vv] == d
                pads_of[d] = (X[vv[k], cc[k]], Y[vv[k], cc[k]])
            XY = (X, Y)
        dead = {}
        for nm, lst in J.get('dead_or_masked', {}).items():
            d = {'P2_IN': 3, 'P2_MID': 2, 'P2_OUT': 1}[nm.split()[-1]]
            dead[d] = np.array([[XY[0][int(s.split(':')[0]), int(s.split(':')[1])],
                                 XY[1][int(s.split(':')[0]), int(s.split(':')[1])]] for s in lst]).reshape(-1, 2)
        Z = np.load(fz)
        if 'p2_chip_n1us' not in Z.files:
            continue
        cfg = chunks[pid]['cfg']
        groups = ['ALL', 'cfg:' + cfg] + (['REF'] if 'trimmed_Km18' in cfg else [])
        for p in (3, 2, 1):
            g = lambda k: Z[f'p{p}_{k}']
            dts = g('dt_scint')
            sel = np.nonzero((dts >= 140) & (dts <= 310))[0]
            if not len(sel):
                continue
            px, py = g('px')[sel], g('py')[sel]
            cx, cy = pads_of[p]
            mrg = PT[NAME[p]]['margin_mm']
            ok = np.ones(len(sel), bool)
            for rx, ry in RING * mrg:
                d2 = (px[:, None] + rx - cx[None]) ** 2 + (py[:, None] + ry - cy[None]) ** 2
                ok &= d2.min(1) < a.fid_r ** 2
            sel, px, py = sel[ok], px[ok], py[ok]
            w = np.where(g('on')[sel], 1.0, -float(g('r_scale')))
            n_any = g('n_any')[sel]
            c1 = g('chip_n1us')[sel].sum(1)
            dwin = g('d_win')[sel]
            dd = dead.get(p, np.zeros((0, 2)))
            ddead = (np.sqrt(((px[:, None] - dd[None, :, 0]) ** 2 + (py[:, None] - dd[None, :, 1]) ** 2).min(1))
                     if len(dd) else np.full(len(sel), 1e9))
            dbig = np.sqrt(((px[:, None] - BIG[None, :, 0]) ** 2 + (py[:, None] - BIG[None, :, 1]) ** 2).min(1))
            cols = dict(w=w, sil=(n_any == 0).astype(float), dead1us=((n_any == 0) & (c1 == 0)).astype(float),
                        found=(dwin <= PT[NAME[p]]['probe_r']).astype(float), slope=g('slope')[sel],
                        tadc=g('tag_adc')[sel].astype(float), key=(g('vmm')[sel] * 64 + g('ch')[sel]).astype(float),
                        ddead=ddead, dbig=dbig, cx=g('cx')[sel], cy=g('cy')[sel])
            for gr in groups:
                L = acc.setdefault(gr, {}).setdefault(NAME[p], {k: [] for k in cols})
                for k, v in cols.items():
                    L[k].append(v)
    out = {}
    for gr, byp in acc.items():
        out[gr] = {}
        for pn, L in byp.items():
            C = {k: np.concatenate(v) for k, v in L.items()}
            w, sil = C['w'], C['sil']
            r = {'n_tags': wfrac(w, sil)[2], 'silent': wfrac(w, sil), 'dead1us': wfrac(w, C['dead1us']),
                 'found': wfrac(w, C['found'])}
            r['vs_slope'] = binned(C['slope'], SLOPE, w, sil)
            r['found_vs_slope'] = binned(C['slope'], SLOPE, w, C['found'])
            r['vs_tag_adc'] = binned(C['tadc'], TADC, w, sil)
            r['vs_dead_dist'] = binned(C['ddead'], DDEAD, w, sil)
            r['vs_big_pillar'] = binned(C['dbig'], [0, a.big_r, 20, 1e9], w, sil)
            pads = []
            f0 = r['silent'][0]
            chi2, ndf = 0.0, 0
            for k in np.unique(C['key']):
                m = C['key'] == k
                f, e, n = wfrac(w[m], sil[m])
                if n >= a.min_pad and f is not None:
                    non = max((w[m] > 0).sum(), 1)
                    chi2 += (f - f0) ** 2 / (f0 * (1 - f0) / non)
                    ndf += 1
                    pads.append([int(k) // 64, int(k) % 64, round(float(C['cx'][m][0]), 1),
                                 round(float(C['cy'][m][0]), 1), n, f])
            r['pads'] = pads
            r['pad_chi2_ndf'] = [round(chi2, 1), ndf]
            out[gr][pn] = r
            s = r['silent']
            if gr.startswith('cfg:'):
                print(f'{gr} {pn}: silent {s[0]} pad chi2/ndf {chi2:.0f}/{ndf}')
                continue
            print(f'{gr} {pn}: tags {r["n_tags"]:.0f}  silent {s[0]:.3f}+-{s[1]:.3f}  no chip hit +-1us '
                  f'{r["dead1us"][0]:.3f}  found {r["found"][0]:.3f}  pad chi2/ndf {chi2:.0f}/{ndf}')
            print('   vs tan(theta):', [(b['lo'], b['f'], int(b['n'])) for b in r['vs_slope']])
            print('   vs tag ADC   :', [(b['lo'], b['f'], int(b['n'])) for b in r['vs_tag_adc']])
            print('   vs dead pad  :', [(b['lo'], b['f'], int(b['n'])) for b in r['vs_dead_dist']])
            print('   vs big pillar:', [(b['lo'], b['f'], int(b['n'])) for b in r['vs_big_pillar']])
    json.dump(out, open(a.out, 'w'))


if __name__ == '__main__':
    sys.exit(main())
