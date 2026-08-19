#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
decode_to_store.py -- rebuild the reduced column store for a sub_run that was
recorded with --drop-columns (or never reduced at all), so the timing scripts
can read it unchanged.

Most of the campaign kept only counts/meta/scalars, and the whole gas-B period
has no hits_store at all, so the raw pcapng is the only way to the hit-level
columns.  Only the five columns the coincidence needs are kept, and only for
the detector and trigger VMMs, which is what makes this fit in memory: one
chunk at a time, filtered before it is stored.

Usage:
  python3 decode_to_store.py <raw_daq_data dir> <out dir> [--captures N]
"""

import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vmm_decode                                            # noqa: E402

KEEP = ('vmm', 'ch', 'bcid', 'tdc', 'srs_timestamp')
DTYPE = dict(vmm=np.int8, ch=np.int8, bcid=np.uint16, tdc=np.uint8,
             srs_timestamp=np.uint64)
WANTED_VMMS = np.arange(0, 20)          # trigger hybrid (0,1) + the three stations


def decode_one(pcap, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    acc = {k: [] for k in KEEP}
    n_raw = 0
    for chunk in vmm_decode.iter_chunks(pcap):
        v = np.asarray(chunk['vmm'])
        n_raw += v.size
        m = np.isin(v, WANTED_VMMS)
        if not m.any():
            continue
        for k in KEEP:
            acc[k].append(np.asarray(chunk[k])[m].astype(DTYPE[k], copy=False))
    if not acc['vmm']:
        return 0, n_raw
    n = 0
    for k in KEEP:
        a = np.concatenate(acc[k])
        np.save(os.path.join(out_dir, f'{k}.npy'), a)
        n = a.size
        acc[k] = None
    return n, n_raw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('raw')
    ap.add_argument('out')
    ap.add_argument('--captures', type=int, default=10)
    a = ap.parse_args()

    caps = sorted(glob.glob(os.path.join(a.raw, '*.pcapng')))[:a.captures]
    print(f'{len(caps)} captures -> {a.out}')
    for i, p in enumerate(caps, 1):
        tag = os.path.splitext(os.path.basename(p))[0]
        d = os.path.join(a.out, tag)
        if os.path.exists(os.path.join(d, 'tdc.npy')):
            print(f'  [{i}/{len(caps)}] {tag}: already done')
            continue
        n, n_raw = decode_one(p, d)
        print(f'  [{i}/{len(caps)}] {tag}: {n_raw} hits -> {n} kept',
              flush=True)


if __name__ == '__main__':
    main()
