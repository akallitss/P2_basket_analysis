#!/usr/bin/env python3
"""Tabulate the VMM<->DREAM match summaries of the whole campaign.

Reads every match_<run>_<sub>.json produced by match_streams.py (by default
the permanent store on EOS) and writes one CSV row per sub_run plus a
markdown table for the logbook.

    ./tabulate_matches.py --dir /eos/.../vmm/matching --csv matching_table.csv

Quality flag:
    OK      >=90% of DREAM events matched, residual rms < 50 ns
    PARTIAL >=50%
    FAIL    below that -- no usable lock (the residual rms then sits at the
            random-coincidence floor, ~1 us, and the fraction near the
            accidental rate)
    NO-VMM  the VMM side recorded nothing for this sub_run
"""
import argparse
import glob
import json
import os
import subprocess
import tempfile

FIELDS = ["run", "sub", "vmm_source", "trigger_vmm", "trigger_ch",
          "n_dream_events", "n_vmm_triggers",
          "n_matched", "match_frac_dream", "vmm_unmatched_frac", "drift_ppm",
          "offset_s", "residual_rms_ns", "residual_med_ns", "n_spills",
          "lock_sigma", "status"]


def run_key(r):
    try:
        return int(r.split("_")[1])
    except (IndexError, ValueError):
        return 10 ** 6


def status_of(d):
    if d.get("vmm_source") == "none":
        return "NO-VMM"
    f = d.get("match_frac_dream", 0.0)
    rms = d.get("residual_rms_ns") or 1e9
    if f >= 0.90 and rms < 50:
        return "OK"
    if f >= 0.50:
        return "PARTIAL"
    return "FAIL"


def fetch_over_xrootd(url, directory, dest):
    """Copy the summaries down with xrdcp.

    The EOS fuse mount serves a file written seconds ago through xrootd as an
    unstatable `-?????????` entry, so reading the sweep's own output back
    through /eos is unreliable; xrootd always sees it.
    """
    ls = subprocess.run(["xrdfs", url, "ls", directory],
                        capture_output=True, text=True, check=True)
    names = [p for p in ls.stdout.split()
             if os.path.basename(p).startswith("match_")
             and p.endswith(".json")]
    for p in names:
        subprocess.run(["xrdcp", "-f", "-s", f"{url}/{p}",
                        f"{dest}/{os.path.basename(p)}"], check=True)
    print(f"fetched {len(names)} summaries from {url}/{directory}")
    return dest


def rows_from(directory):
    rows = []
    for fn in sorted(glob.glob(f"{directory}/match_*.json")):
        with open(fn) as f:
            d = json.load(f)
        rows.append({
            "run": d["run"], "sub": d["sub"],
            "vmm_source": d.get("vmm_source", "?"),
            "trigger_vmm": d.get("trigger_vmm"),
            "trigger_ch": d.get("trigger_ch"),
            "n_dream_events": d["n_dream_events"],
            "n_vmm_triggers": d["n_vmm_triggers"],
            "n_matched": d["n_matched"],
            "match_frac_dream": d["match_frac_dream"],
            "vmm_unmatched_frac": d["vmm_unmatched_frac"],
            "drift_ppm": d["drift_ppm"],
            "offset_s": d["offset_seed_ns"] / 1e9,
            "residual_rms_ns": d["residual_rms_ns"],
            "residual_med_ns": d["residual_med_ns"],
            "n_spills": d["n_spills"],
            "lock_sigma": d.get("lock_sigma"),
            "status": status_of(d)})
    rows.sort(key=lambda r: (run_key(r["run"]), r["sub"]))
    return rows


def fmt(v, spec):
    return "" if v is None else format(v, spec)


def chan(r):
    if r["trigger_vmm"] is None:
        return "-"
    return f"{r['trigger_vmm']}:{r['trigger_ch']}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir",
                    default="/eos/experiment/ntof/data/x17/p2_sps_july/"
                            "vmm/matching")
    ap.add_argument("--csv", default="matching_table.csv")
    ap.add_argument("--md", default=None, help="also write a markdown table")
    ap.add_argument("--xrd", default="root://eospublic.cern.ch",
                    help="fetch the summaries with xrdcp instead of reading "
                         "--dir directly; '' to read the directory as-is")
    args = ap.parse_args()

    tmp = None
    src = args.dir
    if args.xrd:
        tmp = tempfile.TemporaryDirectory()
        src = fetch_over_xrootd(args.xrd, args.dir, tmp.name)
    rows = rows_from(src)
    if not rows:
        raise SystemExit(f"no match_*.json under {args.dir}")

    with open(args.csv, "w") as f:
        f.write(",".join(FIELDS) + "\n")
        for r in rows:
            f.write(",".join(
                "" if r[k] is None else
                (f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]))
                for k in FIELDS) + "\n")

    hdr = (f"| run | sub_run | src | trig ch | DREAM ev | VMM trig | "
           f"matched | VMM spare | rms ns | drift ppm | spills | status |")
    sep = "|" + "---|" * 12
    lines = [hdr, sep]
    for r in rows:
        lines.append(
            f"| {r['run']} | {r['sub']} | {r['vmm_source'][:5]} | "
            f"{chan(r)} | "
            f"{r['n_dream_events']} | {r['n_vmm_triggers']} | "
            f"{r['match_frac_dream']*100:.1f}% | "
            f"{r['vmm_unmatched_frac']*100:.1f}% | "
            f"{fmt(r['residual_rms_ns'], '.1f')} | "
            f"{fmt(r['drift_ppm'], '+.2f')} | {r['n_spills']} | "
            f"{r['status']} |")
    table = "\n".join(lines)
    if args.md:
        with open(args.md, "w") as f:
            f.write(table + "\n")
    print(table)

    n = len(rows)
    ok = [r for r in rows if r["status"] == "OK"]
    bad = [r for r in rows if r["status"] in ("FAIL", "NO-VMM")]
    ev = sum(r["n_dream_events"] for r in rows)
    mt = sum(r["n_matched"] for r in rows)
    print(f"\n{n} sub_runs: {len(ok)} OK, "
          f"{n - len(ok) - len(bad)} partial, {len(bad)} failed")
    print(f"DREAM events {ev}, matched {mt} ({100*mt/ev:.2f}%)")
    if ok:
        import statistics as st
        print("median residual rms (OK) "
              f"{st.median(r['residual_rms_ns'] for r in ok):.1f} ns; "
              "median clock drift "
              f"{st.median(r['drift_ppm'] for r in ok):+.2f} ppm")
    if bad:
        print("no lock: " + ", ".join(f"{r['run']}/{r['sub']}" for r in bad))
    print(f"wrote {os.path.abspath(args.csv)}"
          + (f" and {os.path.abspath(args.md)}" if args.md else ""))


if __name__ == "__main__":
    main()
