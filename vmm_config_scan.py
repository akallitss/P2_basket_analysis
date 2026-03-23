#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 19/11/2025 17:37
Created in PyCharm
Created as vmm_config_scan.py

@author: akallits
"""
import os
import uproot
import pandas as pd
import matplotlib.pyplot as plt
from vmm_mapping import vmm_mapping
import numpy as np
from scipy.ndimage import gaussian_filter1d



# -----------------------------
# DEBUGGING CONTROL
# -----------------------------
DEBUG = False

def debug(msg):
    if DEBUG:
        print(msg)

# ------------------------------------
# PLOT CONTROL SWITCHES (put here)
# ------------------------------------
PLOTS = {
    "adc_hist_per_run": False,
    "adc_hist_separate_vmm": False,
    "mean_vs_peaking": True,
    "debug_samples": False,
    "plot_std_vs_peaking": True,
    "compare_full_vs_cut": True,
    "removed_fraction": False,
    "robust_stats": False,
# --- new QA investigations ---
    "qa_over_threshold_split": True,
    "qa_noise_pedestal_stability": True,
    "qa_mpv_estimation": True,
    "qa_channel_noise": True,
    "qa_adc16_artifact": True,
    "qa_signal_distributions": True,
    "qa_noise_quality_check": True,

}

# -----------------------------
# ROOT FILE LOADER
# -----------------------------
def load_hits_root(filename, branches=None, tree_name="hits"):
    """Load hits from ROOT into a pandas DataFrame."""
    debug(f"Opening ROOT file: {filename}")

    file = uproot.open(filename)

    if tree_name not in file:
        raise ValueError(f"Tree '{tree_name}' not found in ROOT file: {filename}")

    tree = file[tree_name]

    # Handle branch name formatting
    if branches is not None:
        formatted = []
        for b in branches:
            if "/" not in b:
                formatted.append(f"{tree_name}/{b}")
            else:
                formatted.append(b)
        branches = formatted

    df = tree.arrays(branches, library="pd")

    # Clean column names
    df.columns = [c.replace(f"{tree_name}/", "") for c in df.columns]

    debug(f"Loaded DataFrame with columns: {df.columns.tolist()}")
    return df


# -----------------------------
# RUN TABLE MANAGEMENT
# -----------------------------
def load_run_table(csv_path):
    debug(f"Loading run table: {csv_path}")
    return pd.read_csv(csv_path)


def filter_runs(df_run_scan, sg=None, sng=None):
    df = df_run_scan.copy()
    if sg is not None:
        df = df[df["sg"] == sg]
    if sng is not None:
        df = df[df["sng"] == sng]
    debug(f"Filtered runs: {df['run_no'].tolist()}")
    return df["run_no"].tolist()


# -----------------------------
# FILE SYSTEM HELPERS
# -----------------------------
def get_run_dir(base_dir, run_no):
    return os.path.join(base_dir, f"run_{run_no}")


def list_root_files(run_dir, n=None):
    if not os.path.isdir(run_dir):
        debug(f"Missing run dir: {run_dir}")
        return []
    files = sorted([f for f in os.listdir(run_dir) if f.startswith("enp") and f.endswith(".root")])
    return files[:n] if n else files


# -----------------------------
# ANALYSIS & PLOTTING FUNCTIONS
# -----------------------------
def plot_adc_histograms(df_hits, run_no, ax):
    for vmm_id in df_hits["vmm"].unique():
        adc_values = df_hits[df_hits["vmm"] == vmm_id]["adc"]
        ax.hist(adc_values, bins=50, alpha=0.4, label=f"Run {run_no} VMM {vmm_id}")


def compute_adc_stats(df_run_scan, vmm_ids, run_list, data_dir, adc_cut=100, use_over_threshold=True):
    """
    Compute ADC statistics per VMM per run.

    Parameters
    ----------
    df_run_scan : pd.DataFrame
        Run metadata (must include 'run_no' and 'snt').
    vmm_ids : list[int]
        List of VMM IDs to process.
    run_list : list[int]
        List of run numbers.
    data_dir : str
        Base directory with run data.
    adc_cut : int, optional
        Maximum ADC value to consider (used if use_over_threshold=False).
    use_over_threshold : bool, optional
        If True, select hits with over_threshold == 0 instead of ADC cut.

    Returns
    -------
    pd.DataFrame
        Statistics per VMM per run.
    """
    results = []

    for vmm_id in vmm_ids:
        for run_no in run_list:
            run_dir = get_run_dir(data_dir, run_no)
            root_files = list_root_files(run_dir, n=2)
            if not root_files:
                continue

            file_path = os.path.join(run_dir, root_files[1])
            df_hits = load_hits_root(file_path, branches=["adc", "vmm", "ch", "time", "over_threshold"])

            # --- Full ADCs for this VMM ---
            full_adc = df_hits[df_hits["vmm"] == vmm_id]["adc"].values
            if len(full_adc) == 0:
                continue

            # --- Select ADCs based on method ---
            if use_over_threshold:
                # Keep hits where over_threshold == 0
                if "over_threshold" not in df_hits.columns:
                    raise ValueError("Column 'over_threshold' not found in df_hits")
                adc_values = df_hits[(df_hits["vmm"] == vmm_id) & (df_hits["over_threshold"] == 0)]["adc"].values
            else:
                # Traditional ADC cut
                adc_values = full_adc[full_adc < adc_cut]

            if len(adc_values) == 0:
                continue

            # --- Robust statistics ---
            median_adc = np.median(adc_values)
            mad_adc = np.median(np.abs(adc_values - median_adc))
            robust_sigma = 1.4826 * mad_adc

            frac_removed = (len(full_adc) - len(adc_values)) / len(full_adc)

            # --- Standard statistics ---
            mean_adc = adc_values.mean()
            rms_adc = adc_values.std()
            num_hits = len(adc_values)
            rms_error = rms_adc / np.sqrt(2 * (num_hits - 1)) if num_hits > 1 else np.nan

            # --- Peaking time ---
            peaking_time = df_run_scan.loc[df_run_scan["run_no"] == run_no, "snt"].iloc[0]

            results.append({
                "run_no": run_no,
                "vmm_id": vmm_id,
                "peaking_time": peaking_time,
                "num_hits": num_hits,
                "mean_adc": mean_adc,
                "rms_adc": rms_adc,
                "median_adc": median_adc,
                "robust_sigma": robust_sigma,
                "frac_removed": frac_removed,
                "rms_error_adc": rms_error
            })

    return pd.DataFrame(results)

def compare_full_vs_cut(df_hits, vmm_id, run_no, adc_cut=100, use_over_threshold=True):
    """
    Compare full ADC distribution vs selected noise events.

    Parameters
    ----------
    df_hits : pd.DataFrame
        Hit data including 'adc', 'vmm', and optionally 'over_threshold'.
    vmm_id : int
        VMM to plot.
    run_no : int
        Run number (for title).
    adc_cut : int, optional
        ADC threshold for cut (ignored if use_over_threshold=True).
    use_over_threshold : bool, optional
        If True, select hits with over_threshold == 0 instead of ADC cut.
    """
    full_adc = df_hits[df_hits["vmm"] == vmm_id]["adc"].values

    if use_over_threshold:
        if "over_threshold" not in df_hits.columns:
            raise ValueError("Column 'over_threshold' not found in df_hits")
        cut_adc = df_hits[(df_hits["vmm"] == vmm_id) & (df_hits["over_threshold"] == 0)]["adc"].values
        label_cut = "over_threshold == 0"
    else:
        cut_adc = full_adc[full_adc < adc_cut]
        label_cut = f"ADC < {adc_cut}"

    plt.figure()
    plt.hist(full_adc, bins=80, alpha=0.3, label="Full", density=False)
    plt.hist(cut_adc, bins=80, alpha=0.6, label=label_cut, density=False)
    plt.title(f"Run {run_no} – VMM {vmm_id}")
    plt.xlabel("ADC")
    plt.ylabel("Numb of Hits")
    plt.legend()
    plt.show()


def compute_noise_baseline(df_hits, n_sigma=5,
                            adc_low_cut=20,
                            sigma_warn=10.0,
                            sigma_bad=13.0,
                            min_noise_hits=500):
    """
    Compute per-VMM noise pedestal parameters from over_threshold=0 hits.
    Returns a DataFrame with median, robust_sigma, and the n-sigma cut per VMM.

    Parameters
    ----------
    n_sigma : int
        Number of sigma for the noise cut above the median.
    adc_low_cut : int
        Remove hits at or below this ADC value (removes ADC=16 digital artifact).
    sigma_warn : float
        VMMs with robust_sigma above this are flagged as 'warn'.
    sigma_bad : float
        VMMs with robust_sigma above this are flagged as 'bad'.
    min_noise_hits : int
        Minimum number of noise hits required for a reliable estimate.
    """
    records = []

    for vmm_id in sorted(df_hits["vmm"].unique()):
        noise = df_hits[
            (df_hits["vmm"] == vmm_id) &
            (df_hits["over_threshold"] == 0) &
            (df_hits["adc"] > adc_low_cut)        # removes ADC=16 artifact
        ]["adc"].values

        if len(noise) < 100:
            continue

        median = np.median(noise)
        mad = np.median(np.abs(noise - median))
        robust_sigma = 1.4826 * mad
        cut = median + n_sigma * robust_sigma
        tail_frac = (noise > cut).mean() * 100

        # Three-tier quality flag
        if len(noise) < min_noise_hits:
            quality = "bad"
        elif robust_sigma >= sigma_bad:
            quality = "bad"
        elif robust_sigma >= sigma_warn:
            quality = "warn"
        else:
            quality = "ok"

        records.append({
            "vmm_id": vmm_id,
            "n_noise": len(noise),
            "median_adc": median,
            "robust_sigma": robust_sigma,
            "noise_cut": cut,
            "tail_frac_pct": tail_frac,
            "quality": quality
        })

    return pd.DataFrame(records)

def get_clean_signal(df_hits, vmm_id, exclude_trigger_vmms=(0, 1)):
    """
    Extract clean signal hits for a detector VMM:
    - over_threshold == 1
    - adc < 1023 (remove saturated hits)
    - excludes trigger VMMs
    """
    if vmm_id in exclude_trigger_vmms:
        return None

    signal = df_hits[
        (df_hits["vmm"] == vmm_id) &
        (df_hits["over_threshold"] == 1) &
        (df_hits["adc"] < 1023)
    ]["adc"].values

    return signal


def estimate_mpv(adc_values, adc_min, adc_max=800, bins=150, smooth_sigma=2):
    """
    Estimate MPV from clean Landau-like ADC distribution.
    adc_min should come from per-VMM noise cut, not hardcoded.
    """
    counts, edges = np.histogram(adc_values, bins=bins,
                                  range=(adc_min, adc_max))
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed = gaussian_filter1d(counts.astype(float), sigma=smooth_sigma)
    mpv = bin_centers[np.argmax(smoothed)]

    return mpv, counts, smoothed, bin_centers

def get_run_groups(df_run_scan):
    """
    Group runs by configuration (sg, snt, sng).
    Respects 'valid' column if present — invalid runs are excluded.
    """
    # Filter invalid runs if column exists
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

    # Pair sng=0 and sng=1 with same sg and snt
    pairs = []
    df0 = df[df["sng"] == 0]
    df1 = df[df["sng"] == 1]

    for _, r0 in df0.iterrows():
        match = df1[
            (df1["sg"]  == r0["sg"]) &
            (df1["snt"] == r0["snt"])
        ]
        if match.empty:
            print(f"  WARNING: no sng=1 match for run {int(r0['run_no'])} "
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

def qa_over_threshold_split(df_hits, run_no, vmm_ids):
    """
    Investigate the over_threshold flag split for a given run.
    Shows ADC distributions split by over_threshold flag per VMM.
    """
    for vmm_id in vmm_ids:
        df_vmm = df_hits[df_hits["vmm"] == vmm_id]
        n_total = len(df_vmm)
        if n_total == 0:
            continue

        n_above = (df_vmm["over_threshold"] == 1).sum()
        n_below = (df_vmm["over_threshold"] == 0).sum()

        print(f"\nRun {run_no} | VMM {vmm_id}")
        print(f"  Total hits      : {n_total}")
        print(f"  over_threshold=1: {n_above} ({100*n_above/n_total:.1f}%)")
        print(f"  over_threshold=0: {n_below} ({100*n_below/n_total:.1f}%)")

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        for flag, label, color in [(1, "over_threshold=1", "steelblue"),
                                    (0, "over_threshold=0", "orange")]:
            adc = df_vmm[df_vmm["over_threshold"] == flag]["adc"]
            axes[0].hist(adc, bins=100, alpha=0.6,
                         label=f"{label} (n={len(adc)})", color=color)
            axes[1].hist(adc, bins=100, alpha=0.6,
                         label=label, color=color)

        axes[0].set_title(f"Run {run_no} VMM {vmm_id} — full ADC range")
        axes[0].set_xlabel("ADC")
        axes[0].set_ylabel("Hits")
        axes[0].legend()

        axes[1].set_xlim(0, 200)
        axes[1].set_title(f"Run {run_no} VMM {vmm_id} — zoom low ADC")
        axes[1].set_xlabel("ADC")
        axes[1].set_ylabel("Hits")
        axes[1].legend()

        plt.tight_layout()
        plt.show()


def qa_noise_pedestal_stability(df_run_scan, data_dir,
                                 sng1_runs, vmm_groups):
    for group in vmm_groups:
        vmm_ids = group["vmm_ids"]
        label   = group["label"]

        fig, axes = plt.subplots(len(vmm_ids), 1,
                                  figsize=(10, 3.5*len(vmm_ids)))
        if len(vmm_ids) == 1:
            axes = [axes]

        for idx, vmm_id in enumerate(vmm_ids):
            for run_no in sng1_runs:
                run_dir = get_run_dir(data_dir, run_no)
                root_files = list_root_files(run_dir, n=2)
                if not root_files:
                    continue
                file_path = os.path.join(run_dir, root_files[1])
                df_hits = load_hits_root(
                    file_path, branches=["adc", "vmm", "over_threshold"]
                )

                noise = df_hits[
                    (df_hits["vmm"] == vmm_id) &
                    (df_hits["over_threshold"] == 0)
                ]["adc"].values

                if len(noise) < 100:
                    continue

                snt = df_run_scan.loc[
                    df_run_scan["run_no"] == run_no, "snt"
                ].iloc[0]

                axes[idx].hist(noise, bins=80, range=(0, 200),
                               alpha=0.5, histtype="step", linewidth=1.5,
                               density=True,
                               label=f"Run {run_no} (snt={snt:.0f})")

            axes[idx].set_title(f"VMM {vmm_id}")
            axes[idx].set_xlabel("ADC")
            axes[idx].set_ylabel("Density")
            axes[idx].legend(fontsize=9)

        plt.suptitle(
            f"Noise pedestal stability — {label} — sg=4.5, sng=1",
            y=1.01
        )
        plt.tight_layout()
        plt.show()


def qa_signal_distributions(df_hits, detector_vmms, run_no):
    """
    Plot signal ADC distributions per VMM for a given run.
    Also prints saturation statistics.
    """
    # Per VMM signal distributions
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()
    plot_idx = 0

    print(f"\nRun {run_no} — Saturation check per VMM:")
    print(f"{'VMM':>5} {'N_signal':>10} {'N_saturated':>12} "
          f"{'Sat%':>8} {'ADC_max':>10}")

    for vmm_id in sorted(detector_vmms):
        signal = df_hits[
            (df_hits["vmm"] == vmm_id) &
            (df_hits["over_threshold"] == 1)
        ]["adc"].values

        if len(signal) < 100:
            continue

        # Saturation stats
        adc_max = signal.max()
        n_sat = (signal == adc_max).sum()
        print(f"{vmm_id:>5} {len(signal):>10} {n_sat:>12} "
              f"{100*n_sat/len(signal):>8.2f}% {adc_max:>10}")

        if plot_idx >= len(axes):
            break

        axes[plot_idx].hist(signal, bins=100, color="steelblue", alpha=0.7)
        axes[plot_idx].set_xlim(100, 700)
        axes[plot_idx].set_xlabel("ADC")
        axes[plot_idx].set_ylabel("Hits")
        axes[plot_idx].set_title(f"VMM {vmm_id} (n={len(signal)})")
        plot_idx += 1

    plt.suptitle(f"Run {run_no} — Signal distribution per VMM", y=1.02)
    plt.tight_layout()
    plt.show()


def qa_noise_quality_check(df_run_scan, data_dir, runs_to_check):
    """
    Run noise quality check across a list of runs.
    Prints flagged VMMs and their sigma/noise_cut values.
    """
    for run_no in runs_to_check:
        run_dir = get_run_dir(data_dir, run_no)
        root_files = list_root_files(run_dir, n=2)
        if not root_files:
            continue
        file_path = os.path.join(run_dir, root_files[1])
        df_hits = load_hits_root(
            file_path, branches=["adc", "vmm", "over_threshold"]
        )

        df_noise = compute_noise_baseline(df_hits)
        snt = df_run_scan.loc[
            df_run_scan["run_no"] == run_no, "snt"
        ].iloc[0]

        bad = df_noise[df_noise["quality"] == "bad"]
        warn = df_noise[df_noise["quality"] == "warn"]

        if bad.empty and warn.empty:
            print(f"Run {run_no} (snt={snt:.0f}) — all VMMs OK")
        else:
            print(f"\nRun {run_no} (snt={snt:.0f}):")
            if not bad.empty:
                print("  BAD VMMs:")
                print(bad[["vmm_id", "robust_sigma",
                            "noise_cut", "quality"]].to_string(index=False))
            if not warn.empty:
                print("  WARN VMMs:")
                print(warn[["vmm_id", "robust_sigma",
                             "noise_cut", "quality"]].to_string(index=False))

def qa_mpv_estimation(df_hits, df_noise_baseline,
                       detector_vmms, run_no,
                       exclude_trigger_vmms=(0, 1)):
    """
    Visualize MPV estimation for each detector VMM.
    Shows smoothed signal distribution with noise cut and MPV marker.
    """
    for vmm_id in detector_vmms:
        signal_clean = get_clean_signal(df_hits, vmm_id,
                                         exclude_trigger_vmms)
        if signal_clean is None or len(signal_clean) < 100:
            continue

        # Get per-VMM noise cut as adc_min
        row = df_noise_baseline[df_noise_baseline["vmm_id"] == vmm_id]
        if row.empty:
            print(f"VMM {vmm_id} — no noise baseline, skipping")
            continue

        noise_cut = row["noise_cut"].iloc[0]
        noise_quality = row["quality"].iloc[0]

        mpv, counts, smoothed, bin_centers = estimate_mpv(
            signal_clean, adc_min=noise_cut
        )

        scale = counts.max() / smoothed.max()

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.bar(bin_centers, counts,
               width=bin_centers[1] - bin_centers[0],
               color="steelblue", alpha=0.5, label="Clean signal")
        ax.plot(bin_centers, smoothed * scale,
                color="black", linewidth=2, label="Smoothed")
        ax.axvline(noise_cut, color="orange", linestyle="--",
                   label=f"Noise cut={noise_cut:.0f} [{noise_quality}]")
        ax.axvline(mpv, color="red", linestyle="--",
                   label=f"MPV={mpv:.0f}")

        ax.set_xlabel("ADC")
        ax.set_ylabel("Hits")
        ax.set_title(f"VMM {vmm_id} Run {run_no} — MPV estimation")
        ax.legend()
        plt.tight_layout()
        plt.show()

        print(f"VMM {vmm_id} — noise_cut={noise_cut:.1f} "
              f"[{noise_quality}]  MPV={mpv:.1f} ADC  "
              f"N_signal={len(signal_clean)}")


def qa_channel_noise(df_hits, detector_vmms, run_no):
    """
    Per-channel noise investigation.
    Shows median ADC and robust sigma per channel for each detector VMM.
    """
    for vmm_id in detector_vmms:
        df_vmm_noise = df_hits[
            (df_hits["vmm"] == vmm_id) &
            (df_hits["over_threshold"] == 0)
        ]

        if df_vmm_noise.empty:
            continue

        # Compute per-channel stats
        ch_stats = []
        for ch_id, df_ch in df_vmm_noise.groupby("ch"):
            adc = df_ch["adc"].values
            if len(adc) < 10:
                continue
            median = np.median(adc)
            mad = np.median(np.abs(adc - median))
            ch_stats.append({
                "ch": ch_id,
                "n_hits": len(adc),
                "median_adc": median,
                "robust_sigma": 1.4826 * mad
            })

        if not ch_stats:
            continue

        df_ch_stats = pd.DataFrame(ch_stats)

        fig, axes = plt.subplots(1, 2, figsize=(14, 4))

        axes[0].bar(df_ch_stats["ch"], df_ch_stats["median_adc"],
                    color="steelblue", alpha=0.7)
        axes[0].axhline(df_ch_stats["median_adc"].median(),
                        color="red", linestyle="--",
                        label=f"VMM median="
                              f"{df_ch_stats['median_adc'].median():.1f}")
        axes[0].set_xlabel("Channel")
        axes[0].set_ylabel("Median ADC")
        axes[0].set_title(f"VMM {vmm_id} Run {run_no} — "
                          f"Noise median per channel")
        axes[0].legend()

        axes[1].bar(df_ch_stats["ch"], df_ch_stats["robust_sigma"],
                    color="orange", alpha=0.7)
        axes[1].axhline(df_ch_stats["robust_sigma"].median(),
                        color="red", linestyle="--",
                        label=f"VMM median sigma="
                              f"{df_ch_stats['robust_sigma'].median():.1f}")
        axes[1].set_xlabel("Channel")
        axes[1].set_ylabel("Robust sigma (ADC)")
        axes[1].set_title(f"VMM {vmm_id} Run {run_no} — "
                          f"Noise sigma per channel")
        axes[1].legend()

        plt.suptitle(f"VMM {vmm_id} — per channel noise — Run {run_no}")
        plt.tight_layout()
        plt.show()


def qa_adc16_artifact(df_hits, detector_vmms, run_no,
                       artifact_adc=16, low_adc_threshold=25):
    """
    Investigate the ADC=16 digital artifact.
    Shows which channels are affected and how uniformly.
    """
    print(f"\nRun {run_no} — ADC={artifact_adc} artifact investigation:")

    for vmm_id in detector_vmms:
        df_vmm_noise = df_hits[
            (df_hits["vmm"] == vmm_id) &
            (df_hits["over_threshold"] == 0)
        ]

        if df_vmm_noise.empty:
            continue

        total_noise = len(df_vmm_noise)
        artifact = df_vmm_noise[df_vmm_noise["adc"] == artifact_adc]
        n_artifact = len(artifact)

        if n_artifact == 0:
            print(f"  VMM {vmm_id} — no ADC={artifact_adc} hits found")
            continue

        ch_counts = artifact.groupby("ch").size().sort_values(ascending=False)

        print(f"\n  VMM {vmm_id}:")
        print(f"    Artifact hits     : {n_artifact} "
              f"({100*n_artifact/total_noise:.1f}% of noise hits)")
        print(f"    Channels affected : {len(ch_counts)} "
              f"out of {df_vmm_noise['ch'].nunique()}")
        print(f"    Hits/channel — "
              f"min={ch_counts.min()}  "
              f"max={ch_counts.max()}  "
              f"std={ch_counts.std():.1f}")

        # ADC value distribution below threshold
        low = df_vmm_noise[df_vmm_noise["adc"] < low_adc_threshold]["adc"].values
        vals, counts = np.unique(low, return_counts=True)
        print(f"    ADC distribution below {low_adc_threshold}:")
        for v, c in sorted(zip(vals, counts), key=lambda x: -x[1])[:5]:
            print(f"      ADC={v:3d} : {c:6d} hits")

        # Plot per-channel artifact counts
        fig, ax = plt.subplots(figsize=(12, 3))
        ax.bar(ch_counts.index, ch_counts.values,
               color="red", alpha=0.7)
        ax.set_xlabel("Channel")
        ax.set_ylabel(f"ADC={artifact_adc} hit count")
        ax.set_title(f"VMM {vmm_id} Run {run_no} — "
                     f"ADC={artifact_adc} artifact per channel")
        plt.tight_layout()
        plt.show()

def plot_removed_fraction(df_results, vmm_ids):
    plt.figure()
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.plot(df_vmm["peaking_time"], df_vmm["frac_removed"], "o-", label=f"VMM {vmm_id}")
    plt.xlabel("Peaking time (snt)")
    # plt.ylabel("Fraction removed by ADC<100 cut")
    plt.ylabel("Fraction removed by ADC over threshold cut")
    plt.title("Truncation bias vs peaking time")
    plt.legend()
    plt.show()



def plot_mean_vs_peaking(df_results, vmm_ids):
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.errorbar(df_vmm["peaking_time"], df_vmm["mean_adc"],
                     yerr=df_vmm["rms_adc"], fmt="o-", label=f"VMM {vmm_id}")

    plt.xlabel("Peaking Time (snt)")
    # plt.ylabel("Mean ADC (ADC < 100)")
    plt.ylabel("Mean ADC (ADC overthreshold cut)")
    plt.title("Mean ADC vs Peaking Time")
    plt.legend()
    plt.show()

def plot_std_vs_peaking(df_results, vmm_ids):
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.errorbar(df_vmm["peaking_time"], df_vmm["rms_adc"],
                     yerr=df_vmm["rms_error_adc"], fmt="o-", label=f"VMM {vmm_id}")

    plt.xlabel("Peaking Time (snt)")
    # plt.ylabel("Std ADC (ADC < 100)")
    plt.ylabel("Std ADC (ADC over threshold cut)")
    plt.title("Std ADC vs Peaking Time")
    plt.legend()
    plt.show()


def plot_adc_histograms_for_runs(run_list, data_dir):
    fig, ax = plt.subplots()

    for run_no in run_list:
        run_dir = get_run_dir(data_dir, run_no)
        root_files = list_root_files(run_dir, n=2)
        if not root_files:
            continue

        file_path = os.path.join(run_dir, root_files[1])
        df_hits = load_hits_root(file_path, branches=["adc", "vmm", "ch", "time", "over_threshold"])

        for vmm_id in df_hits["vmm"].unique():
            adc_values = df_hits[df_hits["vmm"] == vmm_id]["adc"]
            ax.hist(adc_values, bins=50, alpha=0.4,
                    label=f"Run {run_no} VMM {vmm_id}")

    ax.set_xlabel("ADC Value")
    ax.set_ylabel("Counts")
    ax.set_title("ADC Value Distribution by Run/VMM")
    ax.legend()
    plt.show()

def plot_adc_by_vmm(vmm_ids, run_list, df_run_scan, data_dir):
    for vmm_id in vmm_ids:
        fig, ax = plt.subplots()

        for run_no in run_list:
            run_dir = get_run_dir(data_dir, run_no)
            root_files = list_root_files(run_dir, n=2)
            if not root_files:
                continue

            file_path = os.path.join(run_dir, root_files[1])
            df_hits = load_hits_root(file_path, branches=["adc", "vmm", "ch", "time","over_threshold"])
            adc_values = df_hits[df_hits["vmm"] == vmm_id]["adc"]

            # string of parameters for legend
            param_row = df_run_scan[df_run_scan["run_no"] == run_no]
            param_str = ", ".join([f"{col}={param_row.iloc[0][col]}"
                                   for col in df_run_scan.columns if col != "run_no"])

            ax.hist(adc_values, bins=50, histtype='step',
                    linewidth=1.5, label=param_str, density=True)

        ax.set_xlabel("ADC")
        ax.set_ylabel("Normalized Counts")
        ax.set_title(f"ADC Distribution for VMM {vmm_id}")
        ax.legend()
        plt.show()


def plot_robust_vs_peaking(df_results, vmm_ids):
    # --- Individual plots per VMM ---
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue

        plt.figure()
        plt.plot(df_vmm["peaking_time"], df_vmm["median_adc"], "o-", label="Median ADC")
        plt.plot(df_vmm["peaking_time"], df_vmm["robust_sigma"], "o-", label="Robust σ (MAD)")
        plt.title(f"VMM {vmm_id} – Robust Noise Indicators")
        plt.xlabel("Peaking time (snt)")
        plt.ylabel("ADC / σ")
        plt.legend()
        plt.show()

    # --- Summary plot: all medians ---
    plt.figure()
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.plot(df_vmm["peaking_time"], df_vmm["median_adc"], "o-", label=f"VMM {vmm_id}")
    plt.xlabel("Peaking time (snt)")
    plt.ylabel("Median ADC (ADC < 100)")
    plt.title("Summary of Median ADC for all VMMs")
    plt.legend()
    plt.show()

    # --- Summary plot: all robust sigmas ---
    plt.figure()
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.plot(df_vmm["peaking_time"], df_vmm["robust_sigma"], "o-", label=f"VMM {vmm_id}")
    plt.xlabel("Peaking time (snt)")
    plt.ylabel("Robust σ (MAD) (ADC < 100)")
    plt.title("Summary of Robust σ for all VMMs")
    plt.legend()
    plt.show()


# -----------------------------
# MAIN LOGIC
# -----------------------------
def main():

    # ---- Load run metadata ----
    cnfg_dir = ("/drf/projets/clas12/P2/akallits/")
    df_run_scan = load_run_table(f"{cnfg_dir}vmm_config_scan.csv")

    # ---- Select runs ----
    # run_num_sg_st_sng = filter_runs(df_run_scan, sg=3.0, sng=1.0)
    run_num_sg_st_sng = filter_runs(df_run_scan, sg=4.5, sng=1.0)

    data_dir = "/drf/projets/clas12/cern_202511_p2_alinx/"
    # data_dir = "/media/akallits/EXTERNAL_USB/2312292/Extras/Physics/Post-Doc-Saclay/SPS_beam_test/data/VMM-alinx_data/5kHz-muons-config-scan"
    # data_dir = "/media/akallits/EXTERNAL_USB/2312292/Extras/Physics/Post-Doc-Saclay/SPS_beam_test/data/VMM-alinx_data/15kHz-muons-config-scan"

    # Derive VMM groups directly from vmm_mapping
    trigger_vmms = vmm_mapping["trigger"]["vmm_ids"]
    detector_vmms = [
        vid
        for key, cfg in vmm_mapping.items()
        if key != "trigger"
        for vid in cfg["vmm_ids"]
    ]

    # Groups for per-detector plotting — preserves physical meaning
    vmm_groups = [
        {"label": name, "vmm_ids": cfg["vmm_ids"]}
        for name, cfg in vmm_mapping.items()
        if name != "trigger"
    ]

    run_groups = get_run_groups(df_run_scan)

    print("Run pairs found:")
    print(f"{'sg':>6} {'snt':>6} {'sng0':>8} {'sng1':>8}")
    for p in run_groups["pairs"]:
        print(f"{p['sg']:>6} {p['snt']:>6} "
              f"{p['sng0']:>8} {p['sng1']:>8}")

    print(f"\nAll sng=0 runs: {run_groups['sng0_runs']}")
    print(f"All sng=1 runs: {run_groups['sng1_runs']}")

    input("Press ENTER to continue...")
    # ---- Get unique VMM IDs ----
    vmm_ids = sorted({vid for cfg in vmm_mapping.values() for vid in cfg["vmm_ids"]})

    # ---- Select runs ----
    run_num_sg_st_sng = filter_runs(df_run_scan, sg=4.5, sng=1.0)

    # ---- Compute legacy ADC stats ----
    df_results = compute_adc_stats(
        df_run_scan=df_run_scan,
        vmm_ids=vmm_ids,
        run_list=run_num_sg_st_sng,
        data_dir=data_dir,
        adc_cut=100
    )
    df_results.dropna(inplace=True)
    df_results.to_csv("vmm_adc_analysis.csv", index=False)

    # ---- QA investigations ----
    if PLOTS["qa_over_threshold_split"]:
        run_dir = get_run_dir(data_dir, 79)
        root_files = list_root_files(run_dir, n=2)
        df_hits = load_hits_root(
            os.path.join(run_dir, root_files[1]),
            branches=["adc", "vmm", "ch", "over_threshold"]
        )
        qa_over_threshold_split(df_hits, run_no=79,
                                vmm_ids=detector_vmms)
    if any([PLOTS["qa_channel_noise"],
            PLOTS["qa_adc16_artifact"],
            PLOTS["qa_mpv_estimation"]]):

        run_no_qa = 98  # run to investigate
        run_dir = get_run_dir(data_dir, run_no_qa)
        root_files = list_root_files(run_dir, n=2)
        file_path = os.path.join(run_dir, root_files[1])
        df_hits_qa = load_hits_root(
            file_path,
            branches=["adc", "vmm", "ch", "over_threshold"]
        )
        df_noise_qa = compute_noise_baseline(df_hits_qa)
        if PLOTS["qa_mpv_estimation"]:
            qa_mpv_estimation(df_hits_qa, df_noise_qa,
                              detector_vmms, run_no=run_no_qa)
        if PLOTS["qa_channel_noise"]:
            qa_channel_noise(df_hits_qa, detector_vmms, run_no=run_no_qa)

        if PLOTS["qa_adc16_artifact"]:
            qa_adc16_artifact(df_hits_qa, detector_vmms, run_no=run_no_qa)



    if PLOTS["qa_noise_pedestal_stability"]:
        qa_noise_pedestal_stability(df_run_scan, data_dir,
                                    run_num_sg_st_sng, vmm_groups)

    if PLOTS["qa_signal_distributions"]:
        run_dir = get_run_dir(data_dir, 79)
        root_files = list_root_files(run_dir, n=2)
        df_hits = load_hits_root(
            os.path.join(run_dir, root_files[1]),
            branches=["adc", "vmm", "ch", "over_threshold"]
        )
        qa_signal_distributions(df_hits, detector_vmms, run_no=79)

    if PLOTS["qa_noise_quality_check"]:
        qa_noise_quality_check(df_run_scan, data_dir,
                               runs_to_check=[79, 96, 98])

    # ---- Legacy plots ----

    if PLOTS["adc_hist_per_run"]:
        plot_adc_histograms_for_runs(run_num_sg_st_sng, data_dir)

    if PLOTS["adc_hist_separate_vmm"]:
        plot_adc_by_vmm(vmm_ids, run_num_sg_st_sng, df_run_scan, data_dir)

    if PLOTS["mean_vs_peaking"]:
        plot_mean_vs_peaking(df_results, vmm_ids)

    if PLOTS["plot_std_vs_peaking"]:
        plot_std_vs_peaking(df_results, vmm_ids)

    if PLOTS["compare_full_vs_cut"]:
        for vmm_id in vmm_ids:
            for run_no in run_num_sg_st_sng:
                run_dir = get_run_dir(data_dir, run_no)
                root_files = list_root_files(run_dir, n=2)
                if not root_files:
                    continue
                file_path = os.path.join(run_dir, root_files[1])
                df_hits = load_hits_root(
                    file_path,
                    branches=["adc", "vmm", "ch", "time", "over_threshold"]
                )
                compare_full_vs_cut(df_hits, vmm_id, run_no, adc_cut=100)


    if PLOTS["removed_fraction"]:
        plot_removed_fraction(df_results, vmm_ids)

    if PLOTS["robust_stats"]:
        plot_robust_vs_peaking(df_results, vmm_ids)


if __name__ == "__main__":
    main()


print("bonzo")
