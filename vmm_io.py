#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/25/26 12:53 PM
Created in PyCharm
Created as vmm_io.py

@author: ak271430
"""
"""
vmm_io.py

Data loading and run management utilities for VMM config scan analysis.
"""
import os
import re
import uproot
import pandas as pd
from collections import defaultdict
from vmm_mapping import vmm_mapping

def get_root_file(run_dir, n_files=2, file_index=1):
    """
    Get a specific ROOT file from a run directory.
    Falls back to the last available file if file_index
    exceeds the number of files found.

    Returns
    -------
    str or None
        Full path to the ROOT file, or None if no files found.
    """
    files = list_root_files(run_dir, n=n_files)
    if not files:
        return None
    # Safe index — fall back to last file if index out of range
    idx = min(file_index, len(files) - 1)
    return os.path.join(run_dir, files[idx])


def load_hits_run(run_dir, n_files=1, branches=None):
    """
    Load and concatenate hits from up to n_files ROOT files in run_dir.

    Loads files in sorted order and merges them into a single DataFrame.
    Use this instead of get_root_file + load_hits_root when you want to
    merge statistics across multiple files per run.

    Parameters
    ----------
    run_dir : str
        Path to the run directory containing ROOT files.
    n_files : int
        Maximum number of files to load and merge (default 1).
        Set to a large number (e.g. 99) to load all available files.
    branches : list of str or None
        Branch names to load. None loads all branches.

    Returns
    -------
    pd.DataFrame or None
        Concatenated hit data, or None if no files found.
    """
    files = list_root_files(run_dir, n=n_files)
    if not files:
        return None
    dfs = [
        load_hits_root(os.path.join(run_dir, f), branches=branches)
        for f in files
    ]
    return pd.concat(dfs, ignore_index=True)

def load_hits_root(filename, branches=None, tree_name="hits"):
    """Load hits from ROOT file into a pandas DataFrame."""
    file = uproot.open(filename)

    if tree_name not in file:
        raise ValueError(
            f"Tree '{tree_name}' not found in ROOT file: {filename}"
        )

    tree = file[tree_name]

    if branches is not None:
        formatted = []
        for b in branches:
            if "/" not in b:
                formatted.append(f"{tree_name}/{b}")
            else:
                formatted.append(b)
        branches = formatted

    df = tree.arrays(branches, library="pd")
    df.columns = [c.replace(f"{tree_name}/", "") for c in df.columns]
    return df

def get_connected_channels(vmm_id):
    """
    Return the list of physically connected channels for a VMM.
    Returns None if all channels are active (detector VMMs).
    """
    for key, cfg in vmm_mapping.items():
        if vmm_id in cfg["vmm_ids"]:
            channels = cfg.get("channels", None)
            if channels is None:
                return None   # all channels
            return channels.get(vmm_id, None)
    return None

def load_run_table(csv_path):
    """Load run metadata CSV into a DataFrame."""
    return pd.read_csv(csv_path)


def filter_runs(df_run_scan, sg=None, sng=None):
    """Filter run table by gain and neighbor trigger settings."""
    df = df_run_scan.copy()
    if sg is not None:
        df = df[df["sg"] == sg]
    if sng is not None:
        df = df[df["sng"] == sng]
    return df["run_no"].tolist()


def get_run_dir(base_dir, run_no):
    """Return the directory path for a given run number."""
    return os.path.join(base_dir, f"run_{run_no}")


_SESSION_RE  = re.compile(r'^enp[^_]+_(\d{8}-\d{6})_')
_FILE_IDX_RE = re.compile(r'^enp[^_]+_\d{8}-\d{6}_(\d+)')

def _sort_key(fname):
    m = _FILE_IDX_RE.match(fname)
    return int(m.group(1)) if m else 0

def list_root_files(run_dir, n=None, min_size=10_000):
    """
    List ROOT files in a run directory, optionally limited to n.

    Applies two filters before returning:
    1. Skip files <= min_size bytes (empty ROOT files are ~10 KB).
    2. When multiple pcap sessions exist (different YYYYMMDD-HHMMSS
       prefix), keep only the session with the most files — the stray
       trigger-board pcap that precedes the main data capture is
       typically a single file and is silently dropped.
       A warning is printed whenever a session is skipped.
    """
    if not os.path.isdir(run_dir):
        return []

    candidates = sorted([
        f for f in os.listdir(run_dir)
        if f.startswith("enp") and f.endswith(".root")
    ], key=_sort_key)

    # Drop empty / near-empty files
    candidates = [
        f for f in candidates
        if os.path.getsize(os.path.join(run_dir, f)) > min_size
    ]

    if not candidates:
        return []

    # Group by pcap session prefix (YYYYMMDD-HHMMSS)
    sessions = defaultdict(list)
    for f in candidates:
        m = _SESSION_RE.match(f)
        sessions[m.group(1) if m else "_other"].append(f)

    if len(sessions) > 1:
        # Most files = main data session; break ties by latest timestamp
        main_key = max(sessions, key=lambda k: (len(sessions[k]), k))
        n_skipped = sum(len(v) for k, v in sessions.items()
                        if k != main_key)
        print(f"  NOTE {run_dir}: multiple pcap sessions — "
              f"using {main_key} ({len(sessions[main_key])} files), "
              f"skipping {n_skipped} file(s) from other session(s)")
        candidates = sorted(sessions[main_key], key=_sort_key)

    return candidates[:n] if n else candidates


def iter_hits_files(run_dir, n_files=1, branches=None, file_start=0):
    """
    Yield one DataFrame per ROOT file, processing files one at a time.

    Memory-efficient alternative to load_hits_run for large n_files.
    The caller should filter and delete each yielded DataFrame before
    the next iteration to keep peak memory low.

    Parameters
    ----------
    run_dir : str
        Path to the run directory containing ROOT files.
    n_files : int
        Maximum number of files to iterate over.
    branches : list of str or None
        Branch names to load. None loads all branches.
    file_start : int
        Index of the first file to read (0-based). Use 1 to skip the
        first file when it is known to have corrupted timestamps.

    Yields
    ------
    pd.DataFrame
        Hit data from one ROOT file.
    """
    files = list_root_files(run_dir)
    for fname in files[file_start:file_start + n_files]:
        yield load_hits_root(os.path.join(run_dir, fname),
                             branches=branches)


def get_run_groups(df_run_scan):
    """
    Group runs by configuration (sg, snt, sng).
    Respects 'valid' column if present — invalid runs are excluded.

    Returns
    -------
    dict with keys:
        'by_config'  : {(sg, snt, sng): [run_no, ...]}
        'sng0_runs'  : [run_no, ...]
        'sng1_runs'  : [run_no, ...]
        'pairs'      : [{'sg', 'snt', 'sng0', 'sng1'}, ...]
    """
    if "valid" in df_run_scan.columns:
        df = df_run_scan[df_run_scan["valid"] == 1].copy()
        n_excluded = len(df_run_scan) - len(df)
        if n_excluded > 0:
            excluded = df_run_scan[df_run_scan["valid"] == 0]
            print(f"Excluding {n_excluded} invalid run(s):")
            for _, row in excluded.iterrows():
                print(f"  run {int(row['run_no'])} "
                      f"(sg={row['sg']}, snt={row['snt']}, "
                      f"sng={row['sng']})")
    else:
        df = df_run_scan.copy()

    by_config = {}
    for _, row in df.iterrows():
        key = (row["sg"], row["snt"], row["sng"])
        by_config.setdefault(key, []).append(row["run_no"])

    sng0_runs = df[df["sng"] == 0]["run_no"].tolist()
    sng1_runs = df[df["sng"] == 1]["run_no"].tolist()

    pairs = []
    df0 = df[df["sng"] == 0]
    df1 = df[df["sng"] == 1]

    for _, r0 in df0.iterrows():
        match = df1[
            (df1["sg"]  == r0["sg"]) &
            (df1["snt"] == r0["snt"])
        ]
        if match.empty:
            print(f"  WARNING: no sng=1 match for "
                  f"run {int(r0['run_no'])} "
                  f"(sg={r0['sg']}, snt={r0['snt']})")
            continue
        pairs.append({
            "sg"  : r0["sg"],
            "snt" : r0["snt"],
            "sng0": int(r0["run_no"]),
            "sng1": int(match.iloc[0]["run_no"])
        })

    return {
        "by_config": by_config,
        "sng0_runs": sng0_runs,
        "sng1_runs": sng1_runs,
        "pairs"    : pairs
    }

def main():
    print('bonzo')


if __name__ == '__main__':
    main()
