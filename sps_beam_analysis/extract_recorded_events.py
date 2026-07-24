#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_recorded_events.py

Per-FEU RECORDED-event extraction for the DAQ-overlap-corrected efficiency.

With zero suppression, a trigger in which a FEU had no hit leaves NO row in the
combined hits, so "the FEU recorded this trigger but saw nothing" and "the FEU
never received this trigger" are indistinguishable there. The decoded per-FEU
files resolve it: their nt tree has one entry per RECORDED trigger, including
events with zero fired channels. This script reduces those (large) files to a
tiny per-sub_run summary that travels with the run:

    <sub_run>/recorded_events.npz
        feu<N>_range   = [min_eventId, max_eventId]
        feu<N>_missing = sorted eventIds absent within that range (uint64)
        feu<N>_n       = number of recorded events

The recorded set is reconstructed locally as range minus missing -- the lists
are nearly contiguous, so `missing` is tiny even for multi-million-event runs.
Run this ON THE DAQ HOST (where decoded_root lives); the output is then picked
up by fetch_run_partial.sh and consumed by 22_tag_probe_efficiency.

Usage:  python3 extract_recorded_events.py <run_dir> [<run_dir> ...] [--force]
        (a run_dir = .../runs/<run_name>; every sub_run with decoded_root done)
"""

import os
import re
import sys
import glob
import argparse

import numpy as np
import uproot

_FEU_RE = re.compile(r'_(\d{3})_(\d{2})\.root$')


def subrun_recorded(decoded_dir):
    """{feu: sorted unique eventId array} over all chunks of one sub_run."""
    per_feu = {}
    for fp in sorted(glob.glob(os.path.join(decoded_dir, '*.root'))):
        m = _FEU_RE.search(os.path.basename(fp))
        if not m:
            continue
        feu = int(m.group(2))
        try:
            with uproot.open(f'{fp}:nt') as t:
                ev = t['eventId'].array(library='np').astype(np.uint64)
        except Exception as e:  # truncated/foreign file: report, keep going
            print(f'    [warn] {os.path.basename(fp)}: {e}')
            continue
        per_feu.setdefault(feu, []).append(ev)
    return {f: np.unique(np.concatenate(a)) for f, a in per_feu.items()}


def process_subrun(sub_dir, force=False):
    out = os.path.join(sub_dir, 'recorded_events.npz')
    dec = os.path.join(sub_dir, 'decoded_root')
    if not os.path.isdir(dec):
        return None
    if os.path.isfile(out) and not force:
        print(f'  {os.path.basename(sub_dir)}: exists, skipping (--force to redo)')
        return out
    rec = subrun_recorded(dec)
    if not rec:
        return None
    payload = {}
    line = []
    for feu, ev in sorted(rec.items()):
        lo, hi = int(ev[0]), int(ev[-1])
        full = np.arange(lo, hi + 1, dtype=np.uint64)
        missing = np.setdiff1d(full, ev, assume_unique=True)
        payload[f'feu{feu}_range'] = np.array([lo, hi], dtype=np.uint64)
        payload[f'feu{feu}_missing'] = missing
        payload[f'feu{feu}_n'] = np.uint64(len(ev))
        line.append(f'FEU{feu} {len(ev)} [{lo}..{hi}] miss {len(missing)}')
    np.savez_compressed(out, **payload)
    print(f'  {os.path.basename(sub_dir)}: ' + '; '.join(line))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('run_dirs', nargs='+')
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args()
    for rd in args.run_dirs:
        print(f'== {rd}')
        subs = sorted(d for d in glob.glob(os.path.join(rd, '*'))
                      if os.path.isdir(os.path.join(d, 'decoded_root')))
        if not subs:
            print('  (no sub_runs with decoded_root)')
        for s in subs:
            process_subrun(s, force=args.force)


if __name__ == '__main__':
    main()
