#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tabulate_efficiency.py -- one row per sub_run of uRWELL-referenced VMM efficiency.

Reads the summaries urw_vmm_efficiency.py stored on EOS and writes a CSV plus a
markdown table, the same way tabulate_matches.py does for the matching.

    ./tabulate_efficiency.py --csv efficiency_table.csv --md EFFICIENCY_TABLE.md
"""
import argparse
import glob
import json
import os
import subprocess
import tempfile

BASE = "/eos/experiment/ntof/data/x17/p2_sps_july"
XRD = "root://eospublic.cern.ch"
STATIONS = ("P2_IN", "P2_MID", "P2_OUT")
FIELDS = ["run", "sub", "n_tracks", "n_tracks_timed"]
for s in STATIONS:
    FIELDS += [f"{s}_eff", f"{s}_lo", f"{s}_hi", f"{s}_n", f"{s}_accidental",
               f"{s}_latency_ns", f"{s}_res_dx_mm", f"{s}_rot_deg"]


def run_key(r):
    try:
        return int(r.split("_")[1])
    except (IndexError, ValueError):
        return 10 ** 6


def fetch(url, directory, dest):
    ls = subprocess.run(["xrdfs", url, "ls", directory],
                        capture_output=True, text=True, check=True)
    out = []
    for p in ls.stdout.split():
        if os.path.basename(p).startswith("vmm_eff_") and p.endswith(".json"):
            loc = f"{dest}/{os.path.basename(p)}"
            if subprocess.run(["xrdcp", "-f", "-s", f"{url}/{p}", loc],
                              capture_output=True).returncode == 0:
                out.append(loc)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=f"{BASE}/vmm/matching/efficiency")
    ap.add_argument("--csv", default="efficiency_table.csv")
    ap.add_argument("--md", default=None)
    ap.add_argument("--xrd", default=XRD)
    args = ap.parse_args()

    tmp = tempfile.TemporaryDirectory()
    files = (fetch(args.xrd, args.dir, tmp.name) if args.xrd
             else sorted(glob.glob(f"{args.dir}/vmm_eff_*.json")))
    rows = []
    for fn in files:
        with open(fn) as f:
            d = json.load(f)
        r = {"run": d["run"], "sub": d["sub"],
             "n_tracks": d.get("tracks", {}).get("n_good_track"),
             "n_tracks_timed": d.get("n_tracks_timed")}
        for st in d.get("stations", []):
            k, e = st["station"], st.get("efficiency")
            r[f"{k}_eff"] = None if not e else e["value"]
            r[f"{k}_lo"] = None if not e else e["lo"]
            r[f"{k}_hi"] = None if not e else e["hi"]
            r[f"{k}_n"] = None if not e else e["n"]
            r[f"{k}_accidental"] = None if not e else e["accidental_rate"]
            r[f"{k}_latency_ns"] = (st.get("latency") or {}).get("mu")
            r[f"{k}_res_dx_mm"] = (st.get("residual") or {}).get("dx_rms_mm")
            r[f"{k}_rot_deg"] = (st.get("frame") or {}).get(
                "affine_rotation_deg")
        rows.append(r)
    rows.sort(key=lambda r: (run_key(r["run"]), r["sub"]))
    if not rows:
        raise SystemExit(f"no vmm_eff_*.json under {args.dir}")

    with open(args.csv, "w") as f:
        f.write(",".join(FIELDS) + "\n")
        for r in rows:
            f.write(",".join(
                "" if r.get(k) is None else
                (f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]))
                for k in FIELDS) + "\n")

    lines = ["| run | sub_run | tracks | " + " | ".join(
        s.replace("P2_", "") for s in STATIONS) + " | latency ns | residual mm | frame ° |",
        "|" + "---|" * (3 + len(STATIONS) + 3)]
    for r in rows:
        def p(k):
            v = r.get(f"{k}_eff")
            return "—" if v is None else f"{100*v:.1f}%"
        lat = r.get("P2_MID_latency_ns")
        res = r.get("P2_MID_res_dx_mm")
        rot = r.get("P2_MID_rot_deg")
        lines.append(
            f"| {r['run']} | {r['sub']} | "
            f"{'' if r['n_tracks_timed'] is None else r['n_tracks_timed']} | "
            + " | ".join(p(s) for s in STATIONS) + " | "
            + ("—" if lat is None else f"{lat:.0f}") + " | "
            + ("—" if res is None else f"{res:.2f}") + " | "
            + ("—" if rot is None else f"{rot:.2f}") + " |")
    table = "\n".join(lines)
    if args.md:
        with open(args.md, "w") as f:
            f.write(table + "\n")
    print(table)
    best = {s: max((r for r in rows if r.get(f"{s}_eff") is not None),
                   key=lambda r: r[f"{s}_eff"], default=None)
            for s in STATIONS}
    print(f"\n{len(rows)} sub_runs measured")
    for s in STATIONS:
        b = best[s]
        if b:
            print(f"  best {s}: {100*b[f'{s}_eff']:.1f}% "
                  f"at {b['run']}/{b['sub']}")
    print(f"wrote {os.path.abspath(args.csv)}"
          + (f" and {os.path.abspath(args.md)}" if args.md else ""))


if __name__ == "__main__":
    main()
