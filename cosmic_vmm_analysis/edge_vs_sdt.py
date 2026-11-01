#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Threshold edge of every chip against its threshold register, run by run (lab runs on EOS).

For each run: the first capture file is read (at most --max-mb MB), the over-threshold hits only; per chip the
ADC spectrum, its edge (the ADC below which --q of the hits lie, as in the SPS analysis of 21 October:
sps_beam_analysis/vmm_threshold_tradeoff/note_figs4.py) and the number of hits; the chip's sdt and masks are read
from the configuration file of the run (inventory_runs.chip_settings). Chips with fewer than --min-hits hits are
kept with edge None.

Output OUT.json: {run: {label, sdt: {vmm: sdt}, edge: {vmm: adc}, n: {vmm: hits}, sm1: {vmm: masked}}}.

Usage: edge_vs_sdt.py OUT.json RUN [RUN ...] [--max-mb 200] [--q 0.01] [--min-hits 500]
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

from coincidence_check import read_times
from inventory_runs import chip_settings

RUNS_DIR = '/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/runs'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out_json')
    ap.add_argument('runs', nargs='+', type=int)
    ap.add_argument('--max-mb', type=float, default=200)
    ap.add_argument('--q', type=float, default=0.01)
    ap.add_argument('--min-hits', type=int, default=500)
    ap.add_argument('--runs-dir', default=RUNS_DIR)
    a = ap.parse_args()
    out = {}
    for r in a.runs:
        rdir = os.path.join(a.runs_dir, f'run_{r}')
        subs = sorted(glob.glob(os.path.join(rdir, 'cfg_*')))
        if not subs:
            continue
        sub = os.path.basename(subs[0])
        files = sorted(glob.glob(os.path.join(subs[0], 'raw_daq_data', '*.pcapng')))
        cfgs = glob.glob(os.path.join(rdir, 'p2b-config-*' + sub[len('cfg'):] + '.txt'))
        if not files or not cfgs:
            continue
        chips = chip_settings(cfgs[0])
        h = np.zeros((32, 1024), np.int64)
        for v, t, ch, adc in read_times(files[0], {}, with_adc=True, max_bytes=int(a.max_mb * 1e6)):
            np.add.at(h, (v, adc), 1)
        row = {'label': sub.replace('cfg_gain4.5_peaktime200_opt_cosmic_', ''), 'sdt': {}, 'edge': {}, 'n': {}, 'sm1': {}}
        for v in range(2, 20):
            n = int(h[v].sum())
            row['sdt'][v] = chips.get(v, {}).get('sdt')
            row['sm1'][v] = chips.get(v, {}).get('sm1')
            row['n'][v] = n
            row['edge'][v] = int(np.searchsorted(np.cumsum(h[v]) / n, a.q)) if n >= a.min_hits else None
        out[r] = row
        print(r, row['label'], ' '.join(f'{v}:{row["sdt"][v]}/{row["edge"][v]}' for v in range(2, 20)), flush=True)
    json.dump(out, open(a.out_json, 'w'), indent=1)


if __name__ == '__main__':
    sys.exit(main())
