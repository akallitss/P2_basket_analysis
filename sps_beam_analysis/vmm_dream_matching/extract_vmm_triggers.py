#!/usr/bin/env python3
"""Extract the VMM trigger-channel times of one sub_run into a small npz.

The DREAM trigger is fanned into the VMM/SRS DAQ on VMM 0, channel 44 (for
most of the campaign -- see find_trigger_channel.py), so the only thing stream
matching needs from the VMM side is the arrival time of that one channel: a
few 10^5 numbers out of 10^8-10^9 hits. Two sources:

  hits_store/<capture>/{vmm,ch,offset,bcid,tdc,srs_timestamp}.npy
      the reduced column store, when the online pass kept the hit columns;

  raw_daq_data/<capture>.pcapng
      streamed through vmm_decode.iter_chunks otherwise. Most of the campaign
      was reduced with --drop-columns, so for those sub_runs this is the only
      way -- and it is cheap, because every chunk is filtered down to the
      wanted channel and thrown away immediately.

Time base (see match_streams.py):
    t = srs_timestamp*22.5 + (offset*4096 + bcid)*22.5 + (22.5 - tdc*60/255)
The SRS timestamp tick is 22.5 ns, NOT the 25 ns vmm_decode.derive assumes.

    ./extract_vmm_triggers.py run_57 driftscan_gap150V --out trig.npz
    ./extract_vmm_triggers.py run_25 meshscan_m00V --scan     # busy channels
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
N_VMM, N_CH = 8, 64


def sub_dir_of(run, sub, base=BASE):
    return f"{base}/vmm/runs/{run}/{sub}"


def iter_columns(sub_dir, max_captures=None, source="auto", quiet=False):
    """Yield (capture_name, columns) for every capture of the sub_run.

    Prefers the reduced column store and falls back to decoding the raw
    captures; the two give identical hit columns (verified bit for bit on
    run_32/meshscan_m00V).
    """
    caps = sorted(glob.glob(f"{sub_dir}/hits_store/*/"))
    have_store = bool(caps) and os.path.exists(caps[0] + "vmm.npy")
    if source == "store" and not have_store:
        return
    if have_store and source != "pcapng":
        for cap in caps[:max_captures]:
            name = os.path.basename(cap.rstrip("/"))
            try:
                yield name, {n: np.load(cap + n + ".npy") for n in
                             ("vmm", "ch", "offset", "bcid", "tdc",
                              "srs_timestamp")}
            except Exception as e:               # partial/corrupt capture
                if not quiet:
                    print(f"  skip {name}: {e}")
        return
    try:
        import vmm_decode
    except ImportError:
        sys.exit("vmm_decode.py not importable -- put it on PYTHONPATH "
                 "(it lives in P2_basket_online_analysis)")
    for fn in sorted(glob.glob(f"{sub_dir}/raw_daq_data/*.pcapng"))[:max_captures]:
        name = os.path.basename(fn).replace(".pcapng", "")
        try:
            for c in vmm_decode.iter_chunks(fn):
                yield name, c
        except Exception as e:                   # truncated final capture
            if not quiet:
                print(f"  {name}: decode stopped ({e})")


def hit_times(cols, m):
    return (cols["srs_timestamp"][m].astype(np.float64) * SRS_TICK_NS
            + (cols["offset"][m].astype(np.float64) * 4096
               + cols["bcid"][m].astype(np.float64)) * CLOCK_PERIOD_NS
            + (CLOCK_PERIOD_NS - cols["tdc"][m].astype(np.float64)
               * TAC_SLOPE_NS / TDC_RANGE))


def collect(sub_dir, channels, max_captures=None, source="auto", quiet=False,
            counts=None):
    """Hit times of each (vmm, ch) in `channels`, in one pass over the data.

    `counts`, if given, is an (8, 64) array filled with the per-channel hit
    counts of the same pass -- free, and what the trigger-channel search
    needs to pick its candidates.
    """
    acc = {c: [] for c in channels}
    used = []
    for name, cols in iter_columns(sub_dir, max_captures, source, quiet):
        if counts is not None:
            np.add.at(counts, (np.clip(cols["vmm"].astype(np.int64), 0,
                                       N_VMM - 1),
                               np.clip(cols["ch"].astype(np.int64), 0,
                                       N_CH - 1)), 1)
        for (v, c) in channels:
            m = (cols["vmm"] == v) & (cols["ch"] == c)
            if m.any():
                acc[(v, c)].append(hit_times(cols, m))
        if name not in used:
            used.append(name)
    out = {}
    for c, parts in acc.items():
        out[c] = np.concatenate(parts) if parts else np.empty(0)
    return out, used


def dedup(t):
    """One trigger pulse can fire several VMM samples: keep the first."""
    t = np.sort(t)
    keep = np.ones(len(t), bool)
    keep[1:] = np.diff(t) > DEDUP_NS
    return t[keep]


def source_of(sub_dir):
    caps = sorted(glob.glob(f"{sub_dir}/hits_store/*/"))
    if caps and os.path.exists(caps[0] + "vmm.npy"):
        return "hits_store"
    if glob.glob(f"{sub_dir}/raw_daq_data/*.pcapng"):
        return "pcapng"
    return "none"


def write_npz(out, t, run, sub, src, used, vmm, ch, n_raw):
    meta = dict(run=run, sub=sub, source=src, trigger_vmm=vmm, trigger_ch=ch,
                n_captures=len(used), captures=used,
                n_trigger_hits=int(n_raw), n_triggers=int(len(t)),
                span_s=float((t[-1] - t[0]) / 1e9) if len(t) > 1 else 0.0)
    np.savez_compressed(out, t_ns=t, meta=json.dumps(meta))
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("sub")
    ap.add_argument("--base", default=BASE)
    ap.add_argument("--sub-dir", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--source", choices=["auto", "store", "pcapng"],
                    default="auto")
    ap.add_argument("--vmm", type=int, default=TRIG_VMM)
    ap.add_argument("--ch", type=int, default=TRIG_CH)
    ap.add_argument("--max-captures", type=int, default=None)
    ap.add_argument("--scan", action="store_true",
                    help="just report the busiest channels")
    args = ap.parse_args()

    sub_dir = args.sub_dir or sub_dir_of(args.run, args.sub, args.base)
    src = source_of(sub_dir)
    if src == "none":
        sys.exit(f"no VMM data for {args.run}/{args.sub} under {sub_dir}")

    if args.scan:
        counts = np.zeros((N_VMM, N_CH), np.int64)
        collect(sub_dir, [], args.max_captures or 2, args.source,
                counts=counts)
        tot = counts.sum()
        print(f"{args.run}/{args.sub} [{src}] {tot} hits")
        for k in np.argsort(counts.ravel())[::-1][:10]:
            v, c = divmod(int(k), N_CH)
            print(f"  vmm {v} ch {c:2d}: {counts[v, c]:10d} "
                  f"({100*counts[v, c]/max(tot, 1):5.2f}%)")
        return

    out = args.out or f"vmm_triggers_{args.run}_{args.sub}.npz"
    times, used = collect(sub_dir, [(args.vmm, args.ch)], args.max_captures,
                          args.source)
    t = times[(args.vmm, args.ch)]
    if not len(t):
        sys.exit(f"no hits on vmm {args.vmm} ch {args.ch} for "
                 f"{args.run}/{args.sub} -- try --scan / "
                 f"find_trigger_channel.py")
    n_raw = len(t)
    meta = write_npz(out, dedup(t), args.run, args.sub, src, used,
                     args.vmm, args.ch, n_raw)
    print(json.dumps({k: v for k, v in meta.items() if k != "captures"},
                     indent=1))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
