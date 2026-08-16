#!/usr/bin/env python3
"""Extract the VMM trigger-channel times of one sub_run into a small npz.

The DREAM trigger is fanned into the VMM/SRS DAQ on VMM 0, channel 44, so the
only thing stream matching needs from the VMM side is the arrival time of that
one channel -- a few 10^5 numbers out of 10^8-10^9 hits. Two sources:

  hits_store/<capture>/{vmm,ch,offset,bcid,tdc,srs_timestamp}.npy
      the reduced column store, when the online pass kept the hit columns;

  raw_daq_data/<capture>.pcapng
      streamed through vmm_decode.iter_chunks otherwise. Most of the campaign
      was reduced with --drop-columns, so for those sub_runs this is the only
      way -- and it is cheap, because every chunk is filtered down to channel
      44 and thrown away immediately.

Time base (see match_streams.py):
    t = srs_timestamp*22.5 + (offset*4096 + bcid)*22.5 + (22.5 - tdc*60/255)
The SRS timestamp tick is 22.5 ns, NOT the 25 ns vmm_decode.derive assumes.

    ./extract_vmm_triggers.py run_57 driftscan_gap150V --out trig.npz
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

BASE = "/eos/experiment/ntof/data/x17/p2_sps_july"
CLOCK_PERIOD_NS = 22.5
SRS_TICK_NS = 22.5
TAC_SLOPE_NS = 60.0
TDC_RANGE = 255
TRIG_VMM, TRIG_CH = 0, 44
DEDUP_NS = 1500.0


def trigger_times(vmm, ch, off, bcid, tdc, srs):
    """t_ns of the trigger-channel hits in one chunk (empty array if none)."""
    m = (vmm == TRIG_VMM) & (ch == TRIG_CH)
    if not m.any():
        return np.empty(0)
    return (srs[m].astype(np.float64) * SRS_TICK_NS
            + (off[m].astype(np.float64) * 4096
               + bcid[m].astype(np.float64)) * CLOCK_PERIOD_NS
            + (CLOCK_PERIOD_NS - tdc[m].astype(np.float64)
               * TAC_SLOPE_NS / TDC_RANGE))


def from_store(sub_dir):
    """Trigger times from the reduced column store, or None if it has none."""
    caps = sorted(glob.glob(f"{sub_dir}/hits_store/*/"))
    if not caps or not os.path.exists(caps[0] + "vmm.npy"):
        return None, []
    ts, used = [], []
    for cap in caps:
        try:
            cols = {n: np.load(cap + n + ".npy") for n in
                    ("vmm", "ch", "offset", "bcid", "tdc", "srs_timestamp")}
        except Exception as e:                      # partial/corrupt capture
            print(f"  skip {os.path.basename(cap.rstrip('/'))}: {e}")
            continue
        t = trigger_times(cols["vmm"], cols["ch"], cols["offset"],
                          cols["bcid"], cols["tdc"], cols["srs_timestamp"])
        if len(t):
            ts.append(t)
            used.append(os.path.basename(cap.rstrip("/")))
    if not ts:
        return None, []
    return np.concatenate(ts), used


def from_pcapng(sub_dir, max_captures=None):
    """Trigger times decoded straight from the raw captures."""
    try:
        import vmm_decode
    except ImportError:
        sys.exit("vmm_decode.py not importable -- put it on PYTHONPATH "
                 "(it lives in P2_basket_online_analysis)")
    files = sorted(glob.glob(f"{sub_dir}/raw_daq_data/*.pcapng"))
    if not files:
        return None, []
    if max_captures:
        files = files[:max_captures]
    ts, used = [], []
    for fn in files:
        got = []
        try:
            for c in vmm_decode.iter_chunks(fn):
                t = trigger_times(c["vmm"], c["ch"], c["offset"], c["bcid"],
                                  c["tdc"], c["srs_timestamp"])
                if len(t):
                    got.append(t)
                del c
        except Exception as e:                      # truncated final capture
            print(f"  {os.path.basename(fn)}: decode stopped ({e})")
        if got:
            ts.append(np.concatenate(got))
            used.append(os.path.basename(fn).replace(".pcapng", ""))
        print(f"  {os.path.basename(fn)}: "
              f"{sum(len(g) for g in got)} trigger hits")
    if not ts:
        return None, []
    return np.concatenate(ts), used


def dedup(t):
    """One trigger pulse can fire several VMM samples: keep the first."""
    t = np.sort(t)
    keep = np.ones(len(t), bool)
    keep[1:] = np.diff(t) > DEDUP_NS
    return t[keep]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("sub")
    ap.add_argument("--base", default=BASE)
    ap.add_argument("--sub-dir", default=None,
                    help="override <base>/vmm/runs/<run>/<sub>")
    ap.add_argument("--out", default=None)
    ap.add_argument("--source", choices=["auto", "store", "pcapng"],
                    default="auto")
    ap.add_argument("--max-captures", type=int, default=None)
    args = ap.parse_args()

    sub_dir = args.sub_dir or f"{args.base}/vmm/runs/{args.run}/{args.sub}"
    out = args.out or f"vmm_triggers_{args.run}_{args.sub}.npz"

    t, used, src = None, [], None
    if args.source in ("auto", "store"):
        t, used = from_store(sub_dir)
        src = "hits_store" if t is not None else None
    if t is None and args.source in ("auto", "pcapng"):
        t, used = from_pcapng(sub_dir, args.max_captures)
        src = "pcapng" if t is not None else None
    if t is None:
        sys.exit(f"no trigger hits found for {args.run}/{args.sub} "
                 f"(source={args.source}) under {sub_dir}")

    n_raw = len(t)
    t = dedup(t)
    meta = dict(run=args.run, sub=args.sub, source=src, n_captures=len(used),
                captures=used, n_trigger_hits=int(n_raw),
                n_triggers=int(len(t)),
                span_s=float((t[-1] - t[0]) / 1e9) if len(t) > 1 else 0.0)
    np.savez_compressed(out, t_ns=t, meta=json.dumps(meta))
    print(json.dumps({k: v for k, v in meta.items() if k != "captures"},
                     indent=1))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
