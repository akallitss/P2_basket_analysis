#!/usr/bin/env python3
"""Find which VMM channel actually carries the DREAM trigger in a sub_run.

VMM 0 channel 44 is the trigger fan-out for most of the campaign, but not for
all of it: run_25 has no hits on that channel at all, and in run_21/run_24 it
is busy with something that has no relation to the DREAM event stream. The
cabling is not documented anywhere, so read it off the data instead --
take the busiest channels and ask each one whether its hit times coincide
with DREAM triggers. The real trigger channel answers at >100 sigma; anything
else sits at the random-coincidence floor.

    ./find_trigger_channel.py run_25 meshscan_m00V --out trig.npz

Writes the winning channel's trigger times in the same format as
extract_vmm_triggers.py (and prints the full ranking, which is the useful
thing when nothing wins).
"""
import argparse
import json
import sys

import numpy as np

import extract_vmm_triggers as evt
import match_streams as ms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("sub")
    ap.add_argument("--out", default=None)
    ap.add_argument("--candidates", type=int, default=6)
    ap.add_argument("--scan-captures", type=int, default=2,
                    help="captures used to rank channels by occupancy")
    ap.add_argument("--min-sigma", type=float, default=25.0)
    args = ap.parse_args()

    sub_dir = evt.sub_dir_of(args.run, args.sub)
    src = evt.source_of(sub_dir)
    if src == "none":
        sys.exit(f"no VMM data for {args.run}/{args.sub}")

    # pass 1: who is busy at all
    counts = np.zeros((evt.N_VMM, evt.N_CH), np.int64)
    evt.collect(sub_dir, [], max_captures=args.scan_captures, counts=counts,
                quiet=True)
    tot = int(counts.sum())
    order = np.argsort(counts.ravel())[::-1][:args.candidates]
    cands = [divmod(int(k), evt.N_CH) for k in order]
    print(f"{args.run}/{args.sub} [{src}] {tot} hits in "
          f"{args.scan_captures} capture(s); candidates: "
          + ", ".join(f"vmm{v}ch{c} ({100*counts[v, c]/max(tot,1):.1f}%)"
                      for v, c in cands))

    # pass 2: their times over the whole sub_run, then ask DREAM
    times, used = evt.collect(sub_dir, cands, quiet=True)
    eid, td = ms.load_dream_triggers(args.run, args.sub)
    print(f"dream: {len(td)} events over {(td[-1]-td[0])/1e9:.1f} s")

    ranked = []
    for (v, c) in cands:
        t = evt.dedup(times[(v, c)])
        if len(t) < 1000:
            print(f"  vmm{v} ch{c:2d}: only {len(t)} triggers, skipped")
            continue
        try:
            _, lag, sig = ms.first_lock(t, td, verbose=False)
        except SystemExit as e:                  # no spill structure at all
            print(f"  vmm{v} ch{c:2d}: {len(t)} triggers, no lock ({e})")
            continue
        ranked.append((sig, v, c, lag, len(t)))
        print(f"  vmm{v} ch{c:2d}: {len(t):8d} triggers, "
              f"lock {sig:7.1f} sigma at {lag/1e9:+.6f} s")
    if not ranked:
        sys.exit("no candidate channel had enough hits")
    ranked.sort(reverse=True)
    sig, v, c, lag, n = ranked[0]
    verdict = "OK" if sig >= args.min_sigma else "NO LOCK ON ANY CHANNEL"
    print(f"=> vmm {v} ch {c} at {sig:.1f} sigma  [{verdict}]")

    if args.out and sig >= args.min_sigma:
        t = evt.dedup(times[(v, c)])
        meta = evt.write_npz(args.out, t, args.run, args.sub, src, used,
                             v, c, len(times[(v, c)]))
        print(json.dumps({k: x for k, x in meta.items() if k != "captures"},
                         indent=1))
        print(f"wrote {args.out}")
    sys.exit(0 if sig >= args.min_sigma else 1)


if __name__ == "__main__":
    main()
