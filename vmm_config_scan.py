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
    "compare_full_vs_cut": False,
    "removed_fraction": False,
    "robust_stats": False,
# --- new QA investigations ---
    "qa_over_threshold_split": False,
    "qa_noise_pedestal_stability": False,
    "qa_mpv_estimation": False,
    "qa_channel_noise": False,
    "qa_adc16_artifact": False,
    "qa_signal_distributions": False,
    "qa_noise_quality_check": False,
    "snr_vs_peaking": True,
    "snr_vs_gain":    True,
    "snr_heatmap":    True,
    "qa_noise_sigma_distribution": False,
    "qa_robust_vs_std_comparison": False,
    "qa_mpv_vs_median_comparison": True,


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


def compute_snr(data_dir, pairs, detector_vmms,
                exclude_trigger_vmms=(0, 1)):
    """
    Compute VMM-level SNR for each configuration pair.

    Strategy:
    - Noise sigma  : from sng=1 run (over_threshold=0, adc>20, MAD-based)
    - Signal MPV   : from sng=0 run (over_threshold=1, adc<1023)
    - SNR          : MPV / noise_sigma per VMM

    Parameters
    ----------
    data_dir : str
        Base directory containing run_* folders.
    pairs : list of dict
        Output of get_run_groups()["pairs"].
        Each dict has keys: sg, snt, sng0, sng1.
    detector_vmms : list of int
        VMM IDs to process (excludes trigger VMMs).
    exclude_trigger_vmms : tuple
        VMM IDs to skip entirely.

    Returns
    -------
    pd.DataFrame
        One row per (sg, snt, vmm_id) with SNR and diagnostics.
    """
    results = []

    for pair in pairs:
        sg       = pair["sg"]
        snt      = pair["snt"]
        run_sng0 = pair["sng0"]
        run_sng1 = pair["sng1"]

        print(f"\nProcessing sg={sg} snt={snt} "
              f"| noise=run {run_sng1} | signal=run {run_sng0}")

        # --- Load sng=1 run → noise baseline ---
        run_dir    = get_run_dir(data_dir, run_sng1)
        root_files = list_root_files(run_dir, n=2)
        if not root_files:
            print(f"  WARNING: no files found for run {run_sng1}, skipping")
            continue

        df_hits_noise = load_hits_root(
            os.path.join(run_dir, root_files[1]),
            branches=["adc", "vmm", "over_threshold"]
        )
        df_noise_baseline = compute_noise_baseline(df_hits_noise)

        # --- Load sng=0 run → signal ---
        run_dir    = get_run_dir(data_dir, run_sng0)
        root_files = list_root_files(run_dir, n=2)
        if not root_files:
            print(f"  WARNING: no files found for run {run_sng0}, skipping")
            continue

        df_hits_signal = load_hits_root(
            os.path.join(run_dir, root_files[1]),
            branches=["adc", "vmm", "over_threshold"]
        )

        # --- Per VMM ---
        for vmm_id in detector_vmms:

            # Noise
            noise_row = df_noise_baseline[
                df_noise_baseline["vmm_id"] == vmm_id
            ]
            if noise_row.empty:
                print(f"  VMM {vmm_id}: no noise baseline, skipping")
                continue

            noise_sigma   = noise_row["robust_sigma"].iloc[0]
            noise_cut     = noise_row["noise_cut"].iloc[0]
            noise_quality = noise_row["quality"].iloc[0]

            # Signal
            signal_clean = get_clean_signal(
                df_hits_signal, vmm_id, exclude_trigger_vmms
            )
            if signal_clean is None or len(signal_clean) < 100:
                print(f"  VMM {vmm_id}: insufficient signal hits, skipping")
                continue

            mpv, _, _, _ = estimate_mpv(
                signal_clean, adc_min=noise_cut
            )

            snr = mpv / noise_sigma if noise_sigma > 0 else np.nan

            print(f"  VMM {vmm_id}: noise_sigma={noise_sigma:.1f}  "
                  f"MPV={mpv:.1f}  SNR={snr:.1f}  [{noise_quality}]")

            results.append({
                "sg"           : sg,
                "snt"          : snt,
                "run_sng0"     : run_sng0,
                "run_sng1"     : run_sng1,
                "vmm_id"       : vmm_id,
                "noise_sigma"  : noise_sigma,
                "noise_cut"    : noise_cut,
                "noise_quality": noise_quality,
                "mpv"          : mpv,
                "snr"          : snr,
            })

    return pd.DataFrame(results)


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

def plot_snr_vs_peaking(df_snr):
    """
    SNR vs peaking time, one panel per sg.
    Skips sg values with only one peaking time available.
    """
    df = df_snr[df_snr["noise_quality"] == "ok"].copy()
    sg_values = sorted(df["sg"].unique())

    sg_to_plot = [
        sg for sg in sg_values
        if df[df["sg"] == sg]["snt"].nunique() > 1
    ]
    sg_skip = [sg for sg in sg_values if sg not in sg_to_plot]
    if sg_skip:
        print(f"plot_snr_vs_peaking: skipping sg={sg_skip} "
              f"— only one peaking time available")

    if not sg_to_plot:
        print("plot_snr_vs_peaking: nothing to plot")
        return

    fig, axes = plt.subplots(1, len(sg_to_plot),
                              figsize=(5*len(sg_to_plot), 5),
                              sharey=True)
    if len(sg_to_plot) == 1:
        axes = [axes]

    for ax, sg_val in zip(axes, sg_to_plot):
        df_sg = df[df["sg"] == sg_val]
        for vmm_id in sorted(df_sg["vmm_id"].unique()):
            df_vmm = df_sg[df_sg["vmm_id"] == vmm_id].sort_values("snt")
            if len(df_vmm) < 2:
                continue
            ax.plot(df_vmm["snt"], df_vmm["snr"],
                    "o-", label=f"VMM {vmm_id}")

        ax.set_title(f"sg={sg_val}")
        ax.set_xlabel("Peaking time (snt)")
        ax.set_ylabel("SNR (MPV / noise σ)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle("SNR vs Peaking Time\n"
                 "(quality=ok only)", y=1.02)
    plt.tight_layout()
    plt.show()


def plot_snr_vs_gain(df_snr):
    """
    SNR vs gain, one panel per snt.
    Skips snt values with only one gain setting available.
    """
    df = df_snr[df_snr["noise_quality"] == "ok"].copy()
    snt_values = sorted(df["snt"].unique())

    snt_to_plot = [
        snt for snt in snt_values
        if df[df["snt"] == snt]["sg"].nunique() > 1
    ]
    snt_skip = [snt for snt in snt_values if snt not in snt_to_plot]
    if snt_skip:
        print(f"plot_snr_vs_gain: skipping snt={snt_skip} "
              f"— only one gain setting available")

    if not snt_to_plot:
        print("plot_snr_vs_gain: nothing to plot")
        return

    fig, axes = plt.subplots(1, len(snt_to_plot),
                              figsize=(5*len(snt_to_plot), 5),
                              sharey=True)
    if len(snt_to_plot) == 1:
        axes = [axes]

    for ax, snt_val in zip(axes, snt_to_plot):
        df_snt = df[df["snt"] == snt_val]
        for vmm_id in sorted(df_snt["vmm_id"].unique()):
            df_vmm = df_snt[df_snt["vmm_id"] == vmm_id].sort_values("sg")
            if len(df_vmm) < 2:
                continue
            ax.plot(df_vmm["sg"], df_vmm["snr"],
                    "o-", label=f"VMM {vmm_id}")

        ax.set_title(f"snt={snt_val:.0f}")
        ax.set_xlabel("Gain (sg)")
        ax.set_ylabel("SNR (MPV / noise σ)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle("SNR vs Gain\n"
                 "(quality=ok only)", y=1.02)
    plt.tight_layout()
    plt.show()


def plot_snr_vs_peaking(df_snr):
    df = df_snr.copy()
    sg_values = sorted(df["sg"].unique())

    sg_to_plot = [
        sg for sg in sg_values
        if df[df["sg"] == sg]["snt"].nunique() > 1
    ]
    sg_skip = [sg for sg in sg_values if sg not in sg_to_plot]
    if sg_skip:
        print(f"plot_snr_vs_peaking: skipping sg={sg_skip} "
              f"— only one peaking time available")

    fig, axes = plt.subplots(1, len(sg_to_plot),
                              figsize=(6*len(sg_to_plot), 5),
                              sharey=True)
    if len(sg_to_plot) == 1:
        axes = [axes]

    # Quality marker map
    quality_style = {
        "ok"  : {"linestyle": "-",  "marker": "o", "alpha": 1.0},
        "warn": {"linestyle": "--", "marker": "s", "alpha": 0.7},
        "bad" : {"linestyle": ":",  "marker": "x", "alpha": 0.5},
    }

    for ax, sg_val in zip(axes, sg_to_plot):
        df_sg = df[df["sg"] == sg_val]
        for vmm_id in sorted(df_sg["vmm_id"].unique()):
            df_vmm = df_sg[
                df_sg["vmm_id"] == vmm_id
            ].sort_values("snt")
            if len(df_vmm) < 2:
                continue

            # Plot segments between consecutive points
            # coloring by worst quality in that segment
            line = ax.plot(df_vmm["snt"], df_vmm["snr"],
                           linestyle="-", marker="o",
                           label=f"VMM {vmm_id}")[0]
            color = line.get_color()

            # Overlay warn/bad points with different markers
            for _, row in df_vmm.iterrows():
                style = quality_style[row["noise_quality"]]
                ax.plot(row["snt"], row["snr"],
                        marker=style["marker"],
                        color=color,
                        alpha=style["alpha"],
                        markersize=8,
                        linestyle="")

        ax.set_title(f"sg={sg_val}")
        ax.set_xlabel("Peaking time (snt)")
        ax.set_ylabel("SNR (MPV / noise σ)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # Add quality legend
    from matplotlib.lines import Line2D
    quality_legend = [
        Line2D([0], [0], marker="o", color="gray",
               label="ok",   linestyle="-"),
        Line2D([0], [0], marker="s", color="gray",
               label="warn", linestyle="--", alpha=0.7),
        Line2D([0], [0], marker="x", color="gray",
               label="bad",  linestyle=":",  alpha=0.5),
    ]
    axes[-1].legend(handles=quality_legend,
                    title="Noise quality",
                    fontsize=8, loc="lower right")

    plt.suptitle("SNR vs Peaking Time (all data)", y=1.02)
    plt.tight_layout()
    plt.show()


def plot_snr_vs_gain(df_snr):
    df = df_snr.copy()
    snt_values = sorted(df["snt"].unique())

    snt_to_plot = [
        snt for snt in snt_values
        if df[df["snt"] == snt]["sg"].nunique() > 1
    ]
    snt_skip = [snt for snt in snt_values if snt not in snt_to_plot]
    if snt_skip:
        print(f"plot_snr_vs_gain: skipping snt={snt_skip} "
              f"— only one gain setting available")

    fig, axes = plt.subplots(1, len(snt_to_plot),
                              figsize=(6*len(snt_to_plot), 5),
                              sharey=True)
    if len(snt_to_plot) == 1:
        axes = [axes]

    quality_style = {
        "ok"  : {"linestyle": "-",  "marker": "o", "alpha": 1.0},
        "warn": {"linestyle": "--", "marker": "s", "alpha": 0.7},
        "bad" : {"linestyle": ":",  "marker": "x", "alpha": 0.5},
    }

    for ax, snt_val in zip(axes, snt_to_plot):
        df_snt = df[df["snt"] == snt_val]
        for vmm_id in sorted(df_snt["vmm_id"].unique()):
            df_vmm = df_snt[
                df_snt["vmm_id"] == vmm_id
            ].sort_values("sg")
            if len(df_vmm) < 2:
                continue

            line = ax.plot(df_vmm["sg"], df_vmm["snr"],
                           linestyle="-", marker="o",
                           label=f"VMM {vmm_id}")[0]
            color = line.get_color()

            for _, row in df_vmm.iterrows():
                style = quality_style[row["noise_quality"]]
                ax.plot(row["sg"], row["snr"],
                        marker=style["marker"],
                        color=color,
                        alpha=style["alpha"],
                        markersize=8,
                        linestyle="")

        ax.set_title(f"snt={snt_val:.0f}")
        ax.set_xlabel("Gain (sg)")
        ax.set_ylabel("SNR (MPV / noise σ)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    from matplotlib.lines import Line2D
    quality_legend = [
        Line2D([0], [0], marker="o", color="gray",
               label="ok",   linestyle="-"),
        Line2D([0], [0], marker="s", color="gray",
               label="warn", linestyle="--", alpha=0.7),
        Line2D([0], [0], marker="x", color="gray",
               label="bad",  linestyle=":",  alpha=0.5),
    ]
    axes[-1].legend(handles=quality_legend,
                    title="Noise quality",
                    fontsize=8, loc="lower right")

    plt.suptitle("SNR vs Gain (all data)", y=1.02)
    plt.tight_layout()
    plt.show()

def plot_snr_heatmap(df_snr):
    df = df_snr.copy()

    def config_label(row):
        return f"sg={row['sg']:.1f}\nsnt={row['snt']:.0f}"

    df["config"] = df.apply(config_label, axis=1)

    configs = (
        df[["sg", "snt", "config"]]
        .drop_duplicates()
        .sort_values(["sg", "snt"])["config"]
        .tolist()
    )
    vmm_ids = sorted(df["vmm_id"].unique())

    snr_matrix     = np.full((len(vmm_ids), len(configs)), np.nan)
    quality_matrix = np.full((len(vmm_ids), len(configs)), "",
                              dtype=object)

    for i, vmm_id in enumerate(vmm_ids):
        for j, cfg in enumerate(configs):
            row = df[
                (df["vmm_id"] == vmm_id) &
                (df["config"]  == cfg)
            ]
            if not row.empty:
                snr_matrix[i, j]     = row["snr"].iloc[0]
                quality_matrix[i, j] = row["noise_quality"].iloc[0]

    cell_w = 2.2
    cell_h = 1.0
    fig_w  = max(16, len(configs) * cell_w + 8)
    fig_h  = max(6,  len(vmm_ids) * cell_h + 3)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # More room on right for detector labels + colorbar
    fig.subplots_adjust(left=0.07, right=0.62, top=0.88, bottom=0.15)

    cmap = plt.cm.RdYlGn.copy()
    cmap.set_bad(color="lightgrey")

    masked = np.ma.masked_invalid(snr_matrix)
    im = ax.imshow(masked, aspect="auto", cmap=cmap,
                   vmin=np.nanmin(snr_matrix),
                   vmax=np.nanmax(snr_matrix))

    # Annotate cells
    for i in range(len(vmm_ids)):
        for j in range(len(configs)):
            val     = snr_matrix[i, j]
            quality = quality_matrix[i, j]
            if np.isnan(val):
                ax.text(j, i, "—", ha="center", va="center",
                        fontsize=9, color="gray")
            else:
                suffix = "" if quality == "ok" else f"\n({quality})"
                ax.text(j, i, f"{val:.1f}{suffix}",
                        ha="center", va="center",
                        fontsize=8, fontweight="bold")

    ax.set_xticks(range(len(configs)))
    ax.set_xticklabels(configs, fontsize=9)
    ax.set_yticks(range(len(vmm_ids)))
    ax.set_yticklabels([f"VMM {v}" for v in vmm_ids], fontsize=9)
    ax.set_title("SNR heatmap — VMM vs configuration\n"
                 "(warn/bad cells show quality label)",
                 fontsize=11, pad=12)

    # Detector labels — placed as text annotations
    # instead of twinx to avoid overlap with colorbar
    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    # Use figure coordinates to place detector labels
    # between the plot right edge and colorbar
    for i, vmm_id in enumerate(vmm_ids):
        label = vmm_to_detector.get(vmm_id, "")
        # Convert axes coords to figure coords
        x_fig = 0.635   # just right of plot area
        y_fig = ax.get_position().y0 + \
                ax.get_position().height * (1 - (i + 0.5) / len(vmm_ids))
        fig.text(x_fig, y_fig, label,
                 ha="left", va="center", fontsize=7,
                 clip_on=False)

    # Colorbar placed explicitly — further right of detector labels
    cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.73])
    fig.colorbar(im, cax=cbar_ax, label="SNR")

    plt.show()

def qa_noise_sigma_distribution(df_run_scan, data_dir, sng1_runs,
                                  sigma_warn=10.0, sigma_bad=13.0):
    """
    Plot distribution of robust_sigma across all VMMs and sng=1 runs.
    Shows where the warn/bad thresholds sit relative to the data.
    """
    all_sigmas = []

    for run_no in sng1_runs:
        run_dir = get_run_dir(data_dir, run_no)
        root_files = list_root_files(run_dir, n=2)
        if not root_files:
            continue
        df_hits = load_hits_root(
            os.path.join(run_dir, root_files[1]),
            branches=["adc", "vmm", "over_threshold"]
        )
        df_noise = compute_noise_baseline(df_hits)
        snt = df_run_scan.loc[
            df_run_scan["run_no"] == run_no, "snt"
        ].iloc[0]
        sg = df_run_scan.loc[
            df_run_scan["run_no"] == run_no, "sg"
        ].iloc[0]

        for _, row in df_noise.iterrows():
            all_sigmas.append({
                "run_no" : run_no,
                "sg"     : sg,
                "snt"    : snt,
                "vmm_id" : row["vmm_id"],
                "sigma"  : row["robust_sigma"],
                "quality": row["quality"]
            })

    df_sigmas = pd.DataFrame(all_sigmas)

    fig, ax = plt.subplots(figsize=(10, 4))
    for quality, color in [("ok", "green"),
                            ("warn", "orange"),
                            ("bad", "red")]:
        subset = df_sigmas[df_sigmas["quality"] == quality]["sigma"]
        if not subset.empty:
            ax.hist(subset, bins=30, alpha=0.7,
                    label=quality, color=color)

    ax.axvline(sigma_warn, color="orange", linestyle="--",
               label=f"warn threshold={sigma_warn}")
    ax.axvline(sigma_bad, color="red", linestyle="--",
               label=f"bad threshold={sigma_bad}")
    ax.set_xlabel("Robust sigma (ADC)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of noise sigma across all VMMs and runs")
    ax.legend()
    plt.tight_layout()
    plt.show()

    # Print borderline cases
    print("\nBorderline cases (sigma within 1 ADC of warn threshold):")
    borderline = df_sigmas[
        (df_sigmas["sigma"] >= sigma_warn - 1.0) &
        (df_sigmas["sigma"] <= sigma_warn + 1.0)
    ]
    if borderline.empty:
        print("  None")
    else:
        print(borderline[["run_no", "sg", "snt",
                           "vmm_id", "sigma",
                           "quality"]].to_string(index=False))

    return df_sigmas


def qa_robust_vs_std_comparison(df_run_scan, data_dir, sng1_runs):
    """
    Compare robust_sigma (MAD-based) vs standard deviation
    for all VMMs across all sng=1 runs.
    Validates choice of robust estimator for SNR computation.
    """
    comparison = []

    for run_no in sng1_runs:
        run_dir = get_run_dir(data_dir, run_no)
        root_files = list_root_files(run_dir, n=2)
        if not root_files:
            continue
        df_hits = load_hits_root(
            os.path.join(run_dir, root_files[1]),
            branches=["adc", "vmm", "over_threshold"]
        )
        snt = df_run_scan.loc[
            df_run_scan["run_no"] == run_no, "snt"
        ].iloc[0]
        sg = df_run_scan.loc[
            df_run_scan["run_no"] == run_no, "sg"
        ].iloc[0]

        for vmm_id in sorted(df_hits["vmm"].unique()):
            noise = df_hits[
                (df_hits["vmm"] == vmm_id) &
                (df_hits["over_threshold"] == 0) &
                (df_hits["adc"] > 20)
            ]["adc"].values

            if len(noise) < 100:
                continue

            median       = np.median(noise)
            mad          = np.median(np.abs(noise - median))
            robust_sigma = 1.4826 * mad
            std_sigma    = noise.std()

            comparison.append({
                "run_no"      : run_no,
                "sg"          : sg,
                "snt"         : snt,
                "vmm_id"      : vmm_id,
                "robust_sigma": robust_sigma,
                "std_sigma"   : std_sigma,
                "ratio"       : std_sigma / robust_sigma
                                if robust_sigma > 0 else np.nan
            })

    df_comp = pd.DataFrame(comparison)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sc0 = axes[0].scatter(df_comp["robust_sigma"],
                           df_comp["std_sigma"],
                           c=df_comp["snt"], cmap="viridis",
                           alpha=0.7)
    axes[0].plot([0, df_comp["robust_sigma"].max()],
                 [0, df_comp["robust_sigma"].max()],
                 "r--", label="1:1 line")
    axes[0].set_xlabel("Robust sigma (MAD)")
    axes[0].set_ylabel("Std sigma")
    axes[0].set_title("Robust σ vs Std σ per VMM per run")
    axes[0].legend()
    plt.colorbar(sc0, ax=axes[0], label="snt")

    sc1 = axes[1].scatter(df_comp["robust_sigma"],
                           df_comp["ratio"],
                           c=df_comp["snt"], cmap="viridis",
                           alpha=0.7)
    axes[1].axhline(1.0, color="red", linestyle="--",
                    label="ratio=1")
    axes[1].set_xlabel("Robust sigma (MAD)")
    axes[1].set_ylabel("Std / Robust ratio")
    axes[1].set_title("How much larger is Std vs Robust σ?")
    axes[1].legend()
    plt.colorbar(sc1, ax=axes[1], label="snt")

    plt.tight_layout()
    plt.show()

    print("\nSummary of std/robust ratio:")
    print(f"  Mean  : {df_comp['ratio'].mean():.3f}")
    print(f"  Median: {df_comp['ratio'].median():.3f}")
    print(f"  Min   : {df_comp['ratio'].min():.3f}")
    print(f"  Max   : {df_comp['ratio'].max():.3f}")

    inflated = df_comp[df_comp["ratio"] > 1.5]
    print(f"\nCases where std > 1.5 * robust_sigma "
          f"({len(inflated)} of {len(df_comp)}):")
    if inflated.empty:
        print("  None")
    else:
        print(inflated[["run_no", "sg", "snt", "vmm_id",
                         "robust_sigma", "std_sigma",
                         "ratio"]].to_string(index=False))

    return df_comp

def qa_mpv_vs_median_comparison(df_run_scan, data_dir, pairs,
                                  detector_vmms,
                                  exclude_trigger_vmms=(0, 1)):
    """
    Compare MPV-based vs median-based SNR across all config pairs.
    Validates that MPV is the correct signal estimator for Landau
    distributions — median drifts into the tail and gives
    configuration-dependent bias.

    Shows:
    - Scatter: SNR_mpv vs SNR_median colored by snt
    - Ratio SNR_median/SNR_mpv vs peaking time per VMM
    - Table of MPV, median, mean and both SNR estimates
    """
    records = []

    for pair in pairs:
        sg       = pair["sg"]
        snt      = pair["snt"]
        run_sng0 = pair["sng0"]
        run_sng1 = pair["sng1"]

        # Load signal run
        run_dir    = get_run_dir(data_dir, run_sng0)
        root_files = list_root_files(run_dir, n=2)
        if not root_files:
            continue
        df_hits_signal = load_hits_root(
            os.path.join(run_dir, root_files[1]),
            branches=["adc", "vmm", "over_threshold"]
        )

        # Load noise baseline from sng=1 run
        run_dir_n    = get_run_dir(data_dir, run_sng1)
        root_files_n = list_root_files(run_dir_n, n=2)
        if not root_files_n:
            continue
        df_hits_noise = load_hits_root(
            os.path.join(run_dir_n, root_files_n[1]),
            branches=["adc", "vmm", "over_threshold"]
        )
        df_noise = compute_noise_baseline(df_hits_noise)

        for vmm_id in detector_vmms:
            noise_row = df_noise[df_noise["vmm_id"] == vmm_id]
            if noise_row.empty:
                continue

            noise_sigma   = noise_row["robust_sigma"].iloc[0]
            noise_cut     = noise_row["noise_cut"].iloc[0]
            noise_quality = noise_row["quality"].iloc[0]

            signal = get_clean_signal(
                df_hits_signal, vmm_id, exclude_trigger_vmms
            )
            if signal is None or len(signal) < 100:
                continue

            mpv, _, _, _ = estimate_mpv(signal, adc_min=noise_cut)
            median_sig   = np.median(signal)
            mean_sig     = signal.mean()

            snr_mpv    = mpv        / noise_sigma
            snr_median = median_sig / noise_sigma
            snr_mean   = mean_sig   / noise_sigma

            records.append({
                "sg"           : sg,
                "snt"          : snt,
                "vmm_id"       : vmm_id,
                "noise_sigma"  : noise_sigma,
                "noise_quality": noise_quality,
                "mpv"          : mpv,
                "median_signal": median_sig,
                "mean_signal"  : mean_sig,
                "snr_mpv"      : snr_mpv,
                "snr_median"   : snr_median,
                "snr_mean"     : snr_mean,
                "ratio_med_mpv": snr_median / snr_mpv
            })

    df = pd.DataFrame(records)

    # --- Plot 1 — SNR_mpv vs SNR_median scatter ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    snt_values = sorted(df["snt"].unique())
    cmap       = plt.cm.viridis
    norm       = plt.Normalize(df["snt"].min(), df["snt"].max())

    for snt_val in snt_values:
        df_snt = df[df["snt"] == snt_val]
        color  = cmap(norm(snt_val))
        axes[0].scatter(df_snt["snr_mpv"], df_snt["snr_median"],
                        color=color, alpha=0.7,
                        label=f"snt={snt_val:.0f}")

    max_snr = max(df["snr_mpv"].max(), df["snr_median"].max()) * 1.1
    axes[0].plot([0, max_snr], [0, max_snr],
                 "r--", label="1:1 line")
    axes[0].set_xlabel("SNR (MPV-based)")
    axes[0].set_ylabel("SNR (Median-based)")
    axes[0].set_title("MPV vs Median SNR\n"
                      "points above 1:1 = median overestimates")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3)

    # --- Plot 2 — ratio SNR_median/SNR_mpv vs peaking time ---
    for vmm_id in sorted(df["vmm_id"].unique()):
        df_vmm = df[df["vmm_id"] == vmm_id].sort_values("snt")
        if df_vmm.empty:
            continue
        axes[1].plot(df_vmm["snt"], df_vmm["ratio_med_mpv"],
                     "o-", label=f"VMM {vmm_id}", alpha=0.8)

    axes[1].axhline(1.0, color="red", linestyle="--",
                    label="ratio=1 (no bias)")
    axes[1].set_xlabel("Peaking time (snt)")
    axes[1].set_ylabel("SNR_median / SNR_mpv")
    axes[1].set_title("Median/MPV ratio vs peaking time\n"
                      "increasing ratio = growing bias")
    axes[1].legend(fontsize=7)
    axes[1].grid(True, alpha=0.3)

    # --- Plot 3 — ratio vs sg ---
    for vmm_id in sorted(df["vmm_id"].unique()):
        df_vmm = df[df["vmm_id"] == vmm_id].sort_values("sg")
        if df_vmm.empty:
            continue
        axes[2].plot(df_vmm["sg"], df_vmm["ratio_med_mpv"],
                     "o-", label=f"VMM {vmm_id}", alpha=0.8)

    axes[2].axhline(1.0, color="red", linestyle="--",
                    label="ratio=1 (no bias)")
    axes[2].set_xlabel("Gain (sg)")
    axes[2].set_ylabel("SNR_median / SNR_mpv")
    axes[2].set_title("Median/MPV ratio vs gain\n"
                      "gain dependence of bias")
    axes[2].legend(fontsize=7)
    axes[2].grid(True, alpha=0.3)

    plt.suptitle("MPV vs Median as signal estimator\n"
                 "MPV is stable — median drifts into Landau tail",
                 y=1.02)
    plt.tight_layout()
    plt.show()

    # --- Summary table ---
    print("\nSNR estimator comparison summary:")
    print(f"{'sg':>5} {'snt':>6} {'VMM':>5} {'MPV':>8} "
          f"{'Median':>8} {'Mean':>8} "
          f"{'SNR_mpv':>10} {'SNR_med':>10} {'ratio':>8}")
    for _, row in df.sort_values(["sg", "snt", "vmm_id"]).iterrows():
        print(f"{row['sg']:>5.1f} {row['snt']:>6.0f} "
              f"{row['vmm_id']:>5.0f} "
              f"{row['mpv']:>8.1f} {row['median_signal']:>8.1f} "
              f"{row['mean_signal']:>8.1f} "
              f"{row['snr_mpv']:>10.1f} {row['snr_median']:>10.1f} "
              f"{row['ratio_med_mpv']:>8.2f}")

    return df
# -----------------------------
# MAIN LOGIC
# -----------------------------
def main():

    # ---- Load run metadata ----
    cnfg_dir = ("/drf/projets/clas12/P2/akallits/")
    df_run_scan = load_run_table(f"{cnfg_dir}vmm_config_scan.csv")

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

    df_snr = compute_snr(
        data_dir=data_dir,
        pairs=run_groups["pairs"],
        detector_vmms=detector_vmms
    )

    df_snr.to_csv("vmm_snr_results.csv", index=False)
    #
    # print("\n=== SNR Summary (noise quality=ok only) ===")
    # df_snr_clean = df_snr[df_snr["noise_quality"] == "ok"]
    # print(df_snr_clean[["sg", "snt", "vmm_id",
    #                     "noise_sigma", "mpv",
    #                     "snr"]].to_string(index=False))
    #
    # input("PRESS ENTER TO CONTINUE...")
    #
    # # Check MPV stability for VMM 10 across snt
    # print("VMM 10 MPV check:")
    # print(df_snr[df_snr["vmm_id"] == 10][
    #           ["sg", "snt", "mpv", "noise_sigma", "snr", "noise_quality"]
    #       ].to_string(index=False))
    #
    # # Same for VMM 9 since it's an outlier
    # print("\nVMM 9 MPV check:")
    # print(df_snr[df_snr["vmm_id"] == 9][
    #           ["sg", "snt", "mpv", "noise_sigma", "snr", "noise_quality"]
    #       ].to_string(index=False))
    #
    # input("PRESS ENTER TO CONTINUE...")
    #
    # # Check 1 — are signal runs being reused across pairs?
    # print("Signal run usage per pair:")
    # for p in run_groups["pairs"]:
    #     print(f"  sg={p['sg']} snt={p['snt']} → "
    #           f"signal=run {p['sng0']}  noise=run {p['sng1']}")
    #
    # # Check 2 — VMM 9 noise baseline in run 71 vs other sng=1 runs
    # print("\nVMM 9 noise_sigma across all sng=1 runs:")
    # for run_no in run_groups["sng1_runs"]:
    #     run_dir = get_run_dir(data_dir, run_no)
    #     root_files = list_root_files(run_dir, n=2)
    #     if not root_files:
    #         continue
    #     df_hits = load_hits_root(
    #         os.path.join(run_dir, root_files[1]),
    #         branches=["adc", "vmm", "over_threshold"]
    #     )
    #     df_noise = compute_noise_baseline(df_hits)
    #     row = df_noise[df_noise["vmm_id"] == 9]
    #     if row.empty:
    #         print(f"  run {run_no}: no data")
    #         continue
    #     snt = df_run_scan.loc[df_run_scan["run_no"] == run_no, "snt"].iloc[0]
    #     sg = df_run_scan.loc[df_run_scan["run_no"] == run_no, "sg"].iloc[0]
    #     print(f"  run {run_no} (sg={sg} snt={snt}): "
    #           f"n_noise={row['n_noise'].iloc[0]}  "
    #           f"sigma={row['robust_sigma'].iloc[0]:.2f}  "
    #           f"quality={row['quality'].iloc[0]}")
    #
    # input("PRESS ENTER TO CONTINUE...")


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
    if PLOTS["qa_mpv_vs_median_comparison"]:
        qa_mpv_vs_median_comparison(
            df_run_scan=df_run_scan,
            data_dir=data_dir,
            pairs=run_groups["pairs"],
            detector_vmms=detector_vmms
        )

    if PLOTS["qa_over_threshold_split"]:
        run_dir = get_run_dir(data_dir, 79)
        root_files = list_root_files(run_dir, n=2)
        df_hits = load_hits_root(
            os.path.join(run_dir, root_files[1]),
            branches=["adc", "vmm", "ch", "over_threshold"]
        )
        qa_over_threshold_split(df_hits, run_no=79, vmm_ids=detector_vmms)

    if any([PLOTS["qa_channel_noise"],
            PLOTS["qa_adc16_artifact"],
            PLOTS["qa_mpv_estimation"]]):

        run_no_qa = 98  # run to investigate
        run_dir = get_run_dir(data_dir, run_no_qa)
        root_files = list_root_files(run_dir, n=2)
        file_path = os.path.join(run_dir, root_files[1])
        df_hits_qa = load_hits_root(file_path, branches=["adc", "vmm", "ch", "over_threshold"])
        df_noise_qa = compute_noise_baseline(df_hits_qa)
        if PLOTS["qa_mpv_estimation"]:
            qa_mpv_estimation(df_hits_qa, df_noise_qa, detector_vmms, run_no=run_no_qa)
        if PLOTS["qa_channel_noise"]:
            qa_channel_noise(df_hits_qa, detector_vmms, run_no=run_no_qa)

        if PLOTS["qa_adc16_artifact"]:
            qa_adc16_artifact(df_hits_qa, detector_vmms, run_no=run_no_qa)



    if PLOTS["qa_noise_pedestal_stability"]:
        qa_noise_pedestal_stability(df_run_scan, data_dir, run_num_sg_st_sng, vmm_groups)

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

    if PLOTS["snr_vs_peaking"]:
        plot_snr_vs_peaking(df_snr)

    if PLOTS["snr_vs_gain"]:
        plot_snr_vs_gain(df_snr)

    if PLOTS["snr_heatmap"]:
        plot_snr_heatmap(df_snr)


    if PLOTS["qa_noise_sigma_distribution"]:
        qa_noise_sigma_distribution(
            df_run_scan, data_dir,
            run_groups["sng1_runs"]
        )

    if PLOTS["qa_robust_vs_std_comparison"]:
        qa_robust_vs_std_comparison(
            df_run_scan, data_dir,
            run_groups["sng1_runs"]
        )


if __name__ == "__main__":
    main()


print("bonzo")
