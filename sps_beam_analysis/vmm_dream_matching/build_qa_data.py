#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_qa_data.py -- one JSON describing the state of the whole SPS campaign.

Feeds the static QA page (sps-qa.html on the website), which has no server
behind it: everything the page can show has to be in this file. Run it on
lxplus, where the EOS tree is visible.

    python3 build_qa_data.py --out sps_qa_data.json

What it collects, per (run, sub_run):

  dream    events, chunk files, span, spill count           (from the match)
  vmm      capture source, trigger channel, trigger count   (from the match)
  match    status, matched fraction, residual rms, clock drift, and the
           PER-SPILL table with the reason each unmatched spill failed
  eff      uRWELL-track-referenced efficiency per station, if it has been run
  files    which data products exist on EOS

and per run, the configuration both DAQs recorded (gas, HV, VMM config name,
start time), so the views can group and colour by them.
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import tempfile

BASE = "/eos/experiment/ntof/data/x17/p2_sps_july"
XRD = "root://eospublic.cern.ch"


def run_key(name):
    try:
        return int(name.split("_")[1])
    except (IndexError, ValueError):
        return 10 ** 6


def fetch_json_dir(url, directory, prefix, dest):
    """xrdcp every <prefix>*.json down; the fuse mount is not reliable for
    files written through xrootd minutes earlier."""
    try:
        ls = subprocess.run(["xrdfs", url, "ls", directory],
                            capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return []
    names = [p for p in ls.stdout.split()
             if os.path.basename(p).startswith(prefix) and p.endswith(".json")]
    out = []
    for p in names:
        loc = f"{dest}/{os.path.basename(p)}"
        r = subprocess.run(["xrdcp", "-f", "-s", f"{url}/{p}", loc],
                           capture_output=True)
        if r.returncode == 0:
            out.append(loc)
    return out


def load_run_configs(base):
    """Per run: what the two DAQs recorded about the setup."""
    runs = {}
    for path in sorted(glob.glob(f"{base}/runs/run_*/run_config.json")):
        run = os.path.basename(os.path.dirname(path))
        try:
            with open(path) as f:
                d = json.load(f)
        except Exception:
            continue
        runs[run] = {"run": run, "dream_start": d.get("start_time"),
                     "gas": d.get("gas"), "beam": d.get("beam_type"),
                     "trigger": d.get("trigger")}
    for path in sorted(glob.glob(f"{base}/vmm/runs/run_*/run_config.json")):
        run = os.path.basename(os.path.dirname(path))
        try:
            with open(path) as f:
                d = json.load(f)
        except Exception:
            continue
        e = runs.setdefault(run, {"run": run})
        e["vmm_start"] = d.get("start_time")
        e["vmm_gas"] = d.get("gas")
        e["vmm_trigger"] = d.get("trigger")
        e["vmm_config"] = os.path.basename(
            (d.get("vmm_daq_info") or {}).get("config_file", "") or "")
        e["included_detectors"] = d.get("included_detectors")
        e["run_plan"] = d.get("run_plan")
        hv = d.get("hv_info") or {}
        if hv:
            e["hv"] = {k: v for k, v in list(hv.items())[:12]}
    return runs


def sub_run_files(base, run, sub):
    """Which products exist for this sub_run, and how big the raw data is."""
    d_dream = f"{base}/runs/{run}/{sub}"
    d_vmm = f"{base}/vmm/runs/{run}/{sub}"
    out = {}
    for tag, path, pattern in (
            ("combined_hits", d_dream + "/combined_hits_root", "*.root"),
            ("hits_root", d_dream + "/hits_root", "*.root"),
            ("decoded_root", d_dream + "/decoded_root", "*.root"),
            ("vmm_raw", d_vmm + "/raw_daq_data", "*.pcapng"),
            ("vmm_store", d_vmm + "/hits_store", "*")):
        try:
            out[tag] = len(glob.glob(f"{path}/{pattern}"))
        except OSError:
            out[tag] = 0
    out["recorded_events"] = int(os.path.exists(d_dream
                                                + "/recorded_events.npz"))
    out["hv_monitor"] = int(os.path.exists(d_dream + "/hv_monitor.csv"))
    return out


SPILL_KEYS = ("index", "status", "reason", "where", "n_dream", "n_matched",
              "n_vmm_in_window", "frac", "t0_s", "dur_s", "res_rms_ns")


def trim_spills(spills, keep_all_bad=True, max_good=400):
    """Keep every spill that is not fully matched, and thin the rest.

    The page's matching view is about the ones that failed, and a run with
    200 good spills does not need 200 rows to say so -- but it does need the
    time axis, so keep them, just with fewer fields.
    """
    out = []
    n_good = 0
    for sp in spills:
        row = {k: sp.get(k) for k in SPILL_KEYS}
        if sp.get("status") == "matched":
            n_good += 1
            if n_good > max_good:
                continue
            row.pop("reason", None)
        for k in ("frac", "t0_s", "dur_s", "res_rms_ns"):
            if isinstance(row.get(k), float):
                row[k] = round(row[k], 3)
        out.append(row)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=BASE)
    ap.add_argument("--match-dir", default=None)
    ap.add_argument("--eff-dir", default=None)
    ap.add_argument("--out", default="sps_qa_data.json")
    ap.add_argument("--no-xrd", action="store_true",
                    help="read the summaries through the fuse mount")
    args = ap.parse_args()
    match_dir = args.match_dir or f"{args.base}/vmm/matching"
    eff_dir = args.eff_dir or f"{args.base}/vmm/matching/efficiency"

    tmp = tempfile.TemporaryDirectory()
    if args.no_xrd:
        m_files = sorted(glob.glob(f"{match_dir}/match_*.json"))
        e_files = sorted(glob.glob(f"{eff_dir}/vmm_eff_*.json"))
    else:
        m_files = fetch_json_dir(XRD, match_dir, "match_", tmp.name)
        e_files = fetch_json_dir(XRD, eff_dir, "vmm_eff_", tmp.name)
    print(f"{len(m_files)} match summaries, {len(e_files)} efficiency results")

    eff = {}
    for fn in e_files:
        try:
            with open(fn) as f:
                d = json.load(f)
        except Exception:
            continue
        stations = {}
        for st in d.get("stations", []):
            e = st.get("efficiency")
            stations[st["station"]] = {
                "z_mm": st.get("z_mm"),
                "eff": None if not e else round(e["value"], 5),
                "lo": None if not e else round(e["lo"], 5),
                "hi": None if not e else round(e["hi"], 5),
                "n": None if not e else e["n"],
                "accidental": None if not e else round(e["accidental_rate"], 5),
                "vs_probe_r": None if not e else e.get("vs_probe_r"),
                "vs_window": None if not e else e.get("vs_window"),
                "latency_ns": (st.get("latency") or {}).get("mu"),
                "window_ns": st.get("window_ns"),
                "res_dx_mm": (st.get("residual") or {}).get("dx_rms_mm"),
                "res_dy_mm": (st.get("residual") or {}).get("dy_rms_mm"),
                "rotation_deg": (st.get("frame") or {}).get(
                    "affine_rotation_deg"),
                "per_pad": st.get("per_pad"),
                "error": st.get("error")}
        eff[(d["run"], d["sub"])] = {
            "stations": stations,
            "n_tracks": d.get("tracks", {}).get("n_good_track"),
            "n_tracks_timed": d.get("n_tracks_timed")}

    runs = load_run_configs(args.base)
    subs = []
    for fn in m_files:
        with open(fn) as f:
            m = json.load(f)
        run, sub = m["run"], m["sub"]
        spills = m.get("spills", [])
        row = {
            "run": run, "sub": sub,
            "dream": {"n_events": m["n_dream_events"],
                      "n_covered": m.get("n_dream_covered"),
                      "n_spills": m.get("n_spills")},
            "vmm": {"source": m.get("vmm_source"),
                    "n_triggers": m["n_vmm_triggers"],
                    "trigger_vmm": m.get("trigger_vmm"),
                    "trigger_ch": m.get("trigger_ch")},
            "match": {
                "n_matched": m["n_matched"],
                "frac": round(m["match_frac_dream"], 5),
                "frac_covered": (None if m.get("match_frac_covered") is None
                                 else round(m["match_frac_covered"], 5)),
                "vmm_spare": round(m["vmm_unmatched_frac"], 5),
                "rms_ns": m.get("residual_rms_ns"),
                "drift_ppm": m.get("drift_ppm"),
                "offset_s": round(m.get("offset_seed_ns", 0) / 1e9, 4),
                "lock_sigma": m.get("lock_sigma"),
                "n_spills": m.get("n_spills"),
                "n_spills_matched": m.get("n_spills_matched"),
                "n_spills_unmatched": m.get("n_spills_unmatched"),
                "n_dream_lost": m.get("n_dream_in_unmatched_spills"),
                "spills": trim_spills(spills)},
            "files": sub_run_files(args.base, run, sub),
            "eff": eff.get((run, sub), {}).get("stations"),
            "n_tracks": eff.get((run, sub), {}).get("n_tracks")}
        subs.append(row)

    # sub_runs whose DREAM data exists but which never went through matching
    seen = {(s["run"], s["sub"]) for s in subs}
    for d in sorted(glob.glob(f"{args.base}/runs/run_*/*/combined_hits_root")):
        sub = os.path.basename(os.path.dirname(d))
        run = os.path.basename(os.path.dirname(os.path.dirname(d)))
        if (run, sub) in seen:
            continue
        subs.append({"run": run, "sub": sub,
                     "dream": {"n_events": None}, "vmm": {}, "match": None,
                     "files": sub_run_files(args.base, run, sub),
                     "eff": None})
    subs.sort(key=lambda s: (run_key(s["run"]), s["sub"]))

    doc = {"base": args.base,
           "runs": [runs[k] for k in sorted(runs, key=run_key)],
           "subruns": subs}
    with open(args.out, "w") as f:
        json.dump(doc, f, separators=(",", ":"), default=float)
    n_eff = sum(1 for s in subs if s.get("eff"))
    print(f"{len(subs)} sub_runs, {len(doc['runs'])} runs, "
          f"{n_eff} with efficiency -> {args.out} "
          f"({os.path.getsize(args.out)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
