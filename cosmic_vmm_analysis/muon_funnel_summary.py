#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Merge the per-run muon_funnel.py outputs into configuration groups for the explanation page.

Groups: REF (g4.5 trimmed_Km18), ALL, and every configuration (label 'g<gain> <cfg>'). Per group: number of
scintillator coincidences, per probe the class counts on / off (random, probe window shifted) and the per-chip
classes, nearest-hit distance histograms; per chamber the tag-group states, rejection reasons, number of spatial
pieces and histograms of the extra pieces (separation, same chip, channel difference, time difference, size, ADC).
Example events (all hits of the three chambers) are kept from the REF runs, the pad tables once.

Usage (from $W): muon_funnel_summary.py analysis/muon_funnel OUT.json
"""
import argparse
import glob
import json
import os
import re
from collections import Counter, defaultdict

import numpy as np


def label(cfg):
    g = re.search(r'gain([0-9.]+)', cfg)[1]
    return f"g{g} " + re.sub(r'^cfg_gain[0-9.]+_peaktime200_opt_cosmic_', '', cfg)


def hist(vals, edges):
    return np.histogram(np.asarray(vals, float), edges)[0].tolist() if len(vals) else [0] * (len(edges) - 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('indir'); ap.add_argument('out')
    ap.add_argument('--examples', type=int, default=6)
    a = ap.parse_args()
    chunks = json.load(open('pseudo/chunks.json'))
    cfg_of_run = {v['run']: label(v['cfg']) for v in chunks.values()}
    G = {}
    examples, pads = defaultdict(lambda: defaultdict(list)), None

    def grp(name):
        return G.setdefault(name, {'n_scint': 0, 'runs': [], 'probes': {}, 'rejects': {}, 'xtalk': {},
                                   '_extra': defaultdict(list)})

    for f in sorted(glob.glob(os.path.join(a.indir, 'run*.json'))):
        run = int(re.search(r'run(\d+)\.json', f)[1])
        R = json.load(open(f))
        cfg = cfg_of_run[run]
        if '_tp' in cfg:          # pulser runs: test pulses on some channels, not used
            continue
        is_ref = 'trimmed_Km18' in cfg and 'tp' not in cfg
        pads = pads or R.get('pads')
        for name in {cfg, 'ALL', *(['REF'] if is_ref else [])}:
            g = grp(name)
            g['n_scint'] += R['n_scint']
            g['runs'].append(run)
            for p, v in R['probes'].items():
                q = g['probes'].setdefault(p, {'on': Counter(), 'off': Counter(), 'chip': defaultdict(Counter),
                                               'd_far': Counter(), 'd_found': Counter(),
                                               'paddle': defaultdict(Counter)})
                for k in ('on', 'off', 'd_far', 'd_found', 'margin_cut', 'margin_cut_off'):
                    q.setdefault(k, Counter()).update(v.get(k, {}))
                for k in ('probe_r', 'margin_mm'):          # muon_funnel.py --pointing: per-probe radius, margin
                    if k in v:
                        q.setdefault(k + '_set', set()).add(v[k])
                for c, cc in v['chip'].items():
                    q['chip'][c].update(cc)
                for c, cc in v.get('paddle', {}).items():
                    q['paddle'][c].update(cc)
            for d, v in R.get('xtalk', {}).items():
                q = g['xtalk'].setdefault(d, {'n': [0] * 11, 'on': [[0] * 11 for _ in range(3)],
                                              'off': [[0] * 11 for _ in range(3)]})
                q['n'] = [x + y for x, y in zip(q['n'], v['n'])]
                for w in ('on', 'off'):
                    q[w] = [[x + y for x, y in zip(a_, b_)] for a_, b_ in zip(q[w], v[w])]
            for d, v in R['rejects'].items():
                q = g['rejects'].setdefault(d, {'states': Counter(), 'states_off': Counter(), 'reasons': Counter(),
                                                'pieces': Counter()})
                for k in ('states', 'states_off', 'reasons', 'pieces'):
                    q[k].update(v[k])
                g['_extra'][d] += v['extra']
        if is_ref:
            for p, byc in R['examples'].items():
                for c, lst in byc.items():
                    examples[p][c] += lst
    for name, g in G.items():
        for d, ex in g.pop('_extra').items():
            E = np.array(ex, float).reshape(-1, 6) if ex else np.zeros((0, 6))
            same = E[:, 3] > 0
            g['rejects'][d]['extra'] = {
                'n': int(len(E)), 'same_chip': int(same.sum()),
                'sep_edges': list(range(0, 520, 20)), 'sep': hist(E[:, 0], range(0, 520, 20)),
                'sep_same': hist(E[same, 0], range(0, 520, 20)),
                'dch_edges': list(range(0, 66, 2)), 'dch': hist(E[same, 4], range(0, 66, 2)),
                'dt_edges': list(range(-200, 220, 20)), 'dt': hist(E[:, 5], range(-200, 220, 20)),
                'size': dict(Counter(np.minimum(E[:, 1], 5).astype(int).astype(str).tolist())),
                'adc_edges': list(range(0, 1100, 50)), 'adc': hist(E[:, 2], range(0, 1100, 50))}
        for p in g['probes'].values():
            p['chip'] = {c: dict(v) for c, v in p['chip'].items()}
            p['paddle'] = {c: dict(v) for c, v in p['paddle'].items()}
            for k in ('probe_r', 'margin_mm'):
                if k + '_set' in p:
                    p[k] = sorted(p.pop(k + '_set'))
    out = {'groups': G, 'pads': pads,
           'examples': {p: {c: lst[:a.examples] for c, lst in byc.items()} for p, byc in examples.items()}}
    json.dump(out, open(a.out, 'w'))
    for name in ('REF', 'ALL'):
        g = G[name]
        print(name, g['n_scint'], sorted(set(g['runs'])))
        for p, v in g['probes'].items():
            print('  ', p, dict(v['on']), 'off', dict(v['off']), 'probe_r', v.get('probe_r'), 'margin',
                  v.get('margin_mm'), 'margin_cut', dict(v.get('margin_cut', {})))
        for d, v in g['rejects'].items():
            print('  ', d, dict(v['states']), dict(v['reasons']), dict(v['pieces']), 'same chip',
                  v['extra']['same_chip'], '/', v['extra']['n'])
        for p, v in g['probes'].items():
            print('  paddle', p, {b: dict(c) for b, c in sorted(v['paddle'].items())})
        for d, v in g['xtalk'].items():
            print('  xtalk', d, 'n', v['n'], 'on k=1', v['on'][0], 'off k=1', v['off'][0])


if __name__ == '__main__':
    main()
