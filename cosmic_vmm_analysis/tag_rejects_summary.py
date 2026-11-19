#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summary of tag_rejects.py over the weekend runs, per configuration and chamber.

  chamber table  per scintillator coincidence: P(clean tag group), P(rejected group), P(none), random-corrected
                 (on - off), the rejected share of the active muons, the reasons and the spatial sub-clusters
  probe table    muons with both tag chambers active: share LOST (>= 1 tag group rejected), of them with one clean
                 tag; P(probe any hit) for USED and for LOST muons

Usage (from $W, env.sh): tag_rejects_summary.py analysis/tag_rejects OUT_PREFIX
"""
import argparse
import glob
import json
import os
import re
from collections import Counter, defaultdict

import pandas as pd

ORDER = ['P2_IN', 'P2_MID', 'P2_OUT']


def label(cfg):
    g = re.search(r'gain([0-9.]+)', cfg)[1]
    name = re.sub(r'^cfg_gain[0-9.]+_peaktime200_opt_cosmic_', '', cfg)
    return f'g{g} {name}'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('indir'); ap.add_argument('out')
    a = ap.parse_args()
    chunks = json.load(open('pseudo/chunks.json'))
    C = defaultdict(Counter)            # (cfg, det, window) -> counts
    Rs = defaultdict(Counter)           # (cfg, det) -> reason combos (on)
    Sub = defaultdict(Counter)          # (cfg, det) -> sub-cluster counts (on)
    P = defaultdict(Counter)            # (cfg, probe, window) -> probe counts
    for f in sorted(glob.glob(os.path.join(a.indir, 'run*.json'))):
        for path, r in json.load(open(f)).items():
            pid = re.search(r'eff_run(\d+)\.json', path)[1]
            cfg = label(chunks[pid]['cfg'])
            for key in (cfg, 'ALL'):
                for d, v in r['chambers'].items():
                    for w in ('on', 'off'):
                        C[(key, d, w)].update({k: v[w][k] for k in ('n_scint', 'clean', 'rejected', 'none')})
                    Rs[(key, d)].update(v['on']['reasons'])
                    Sub[(key, d)].update(v['on']['subclusters'])
                for d, v in r['probes'].items():
                    for w in ('on', 'off'):
                        P[(key, d, w)].update(v[w])
    rows = []
    for (cfg, d, w) in sorted(C):
        if w != 'on':
            continue
        on, off = C[(cfg, d, 'on')], C[(cfg, d, 'off')]
        n = on['n_scint']
        cl = (on['clean'] - off['clean']) / n
        rj = (on['rejected'] - off['rejected']) / n
        sub = Sub[(cfg, d)]
        nrej = sum(sub.values())
        rows.append(dict(config=cfg, det=d, scint=n, clean=cl, rejected=rj, none=1 - (on['clean'] + on['rejected']) / n,
                         rejected_of_active=rj / (cl + rj) if cl + rj > 0 else float('nan'),
                         rand_rejected=off['rejected'] / n,
                         sub1=sub.get('1', 0) / nrej if nrej else float('nan'),
                         sub2=sub.get('2', 0) / nrej if nrej else float('nan'),
                         sub3plus=sum(v for k, v in sub.items() if int(k) >= 3) / nrej if nrej else float('nan'),
                         reasons='; '.join(f'{k} {v / nrej:.2f}' for k, v in Rs[(cfg, d)].most_common(4)) if nrej else ''))
    T = pd.DataFrame(rows)
    T['det'] = pd.Categorical(T.det, ORDER)
    T = T.sort_values(['config', 'det'])
    T.to_csv(a.out + '_chambers.csv', index=False, float_format='%.4f')
    rows = []
    for (cfg, d, w) in sorted(P):
        if w != 'on':
            continue
        on, off = P[(cfg, d, 'on')], P[(cfg, d, 'off')]
        act = on['both_active'] - off['both_active']
        used, lost = on['used'] - off['used'], on['lost'] - off['lost']
        rows.append(dict(config=cfg, probe=d, both_active=act, lost_share=lost / act if act else float('nan'),
                         lost_with_one_clean=(on['lost_one_clean'] - off['lost_one_clean']) / lost if lost else float('nan'),
                         probe_any_used=(on['probe_hit_used'] - off['probe_hit_used']) / used if used else float('nan'),
                         probe_any_lost=(on['probe_hit_lost'] - off['probe_hit_lost']) / lost if lost else float('nan'),
                         n_lost=lost))
    Q = pd.DataFrame(rows)
    Q['probe'] = pd.Categorical(Q.probe, ORDER)
    Q = Q.sort_values(['config', 'probe'])
    Q.to_csv(a.out + '_probes.csv', index=False, float_format='%.4f')
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
    print('== per chamber, per scintillator coincidence (random-corrected)')
    print(T.round(3).to_string(index=False))
    print('\n== per probe: muons with both tag chambers active')
    print(Q.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
