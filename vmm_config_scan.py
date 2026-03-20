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
    "robust_stats": False

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


def compute_noise_baseline(df_hits, n_sigma=5):
    """
    Compute per-VMM noise pedestal parameters from over_threshold=0 hits.
    Returns a DataFrame with median, robust_sigma, and the n-sigma cut per VMM.
    """
    records = []

    for vmm_id in sorted(df_hits["vmm"].unique()):
        noise = df_hits[
            (df_hits["vmm"] == vmm_id) &
            (df_hits["over_threshold"] == 0)
        ]["adc"].values

        if len(noise) < 100:
            continue

        median = np.median(noise)
        mad = np.median(np.abs(noise - median))
        robust_sigma = 1.4826 * mad
        cut = median + n_sigma * robust_sigma
        tail_frac = (noise > cut).mean() * 100

        records.append({
            "vmm_id": vmm_id,
            "n_noise": len(noise),
            "median_adc": median,
            "robust_sigma": robust_sigma,
            "noise_cut": cut,
            "tail_frac_pct": tail_frac
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

    # ---- Get unique VMM IDs ----
    vmm_ids = sorted({vid for cfg in vmm_mapping.values() for vid in cfg["vmm_ids"]})

    # ---- Compute stats ----
    df_results = compute_adc_stats(
        df_run_scan=df_run_scan,
        vmm_ids=vmm_ids,
        run_list=run_num_sg_st_sng,
        data_dir=data_dir,
        adc_cut=100
    )

    df_results.dropna(inplace=True)
    df_results.to_csv("vmm_adc_analysis.csv", index=False)

    # Step 1 - pick two controlled runs
    run_sng0 = 78  # sg=4.5, snt=200, sng=0
    run_sng1 = 79  # sg=4.5, snt=200, sng=1

    vmm_id = vmm_ids[5]  # just one VMM for now

    for run_no in [run_sng0, run_sng1]:
        run_dir = get_run_dir(data_dir, run_no)
        root_files = list_root_files(run_dir, n=2)
        file_path = os.path.join(run_dir, root_files[1])

        df_hits = load_hits_root(
            file_path,
            branches=["adc", "vmm", "ch", "over_threshold"]
        )

        df_vmm = df_hits[df_hits["vmm"] == vmm_id]
        n_total = len(df_vmm)
        n_above = (df_vmm["over_threshold"] == 1).sum()
        n_below = (df_vmm["over_threshold"] == 0).sum()

        print(f"\nRun {run_no} | VMM {vmm_id}")
        print(f"  Total hits      : {n_total}")
        print(f"  over_threshold=1: {n_above} ({100 * n_above / n_total:.1f}%)")
        print(f"  over_threshold=0: {n_below} ({100 * n_below / n_total:.1f}%)")

    run_dir = get_run_dir(data_dir, 79)
    root_files = list_root_files(run_dir, n=2)
    file_path = os.path.join(run_dir, root_files[1])

    df_hits = load_hits_root(
        file_path,
        branches=["adc", "vmm", "ch", "over_threshold"]
    )

    df_vmm = df_hits[df_hits["vmm"] == 11]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Left plot — full range
    for flag, label, color in [(1, "over_threshold=1", "steelblue"),
                               (0, "over_threshold=0", "orange")]:
        adc = df_vmm[df_vmm["over_threshold"] == flag]["adc"]
        axes[0].hist(adc, bins=100, alpha=0.6, label=f"{label} (n={len(adc)})", color=color)

    axes[0].set_xlabel("ADC")
    axes[0].set_ylabel("Hits")
    axes[0].set_title("Run 79 VMM 11 — full ADC range")
    axes[0].legend()

    # Right plot — zoom into low ADC to see noise region
    for flag, label, color in [(1, "over_threshold=1", "steelblue"),
                               (0, "over_threshold=0", "orange")]:
        adc = df_vmm[df_vmm["over_threshold"] == flag]["adc"]
        axes[1].hist(adc, bins=100, alpha=0.6, label=label, color=color)

    axes[1].set_xlim(0, 200)
    axes[1].set_xlabel("ADC")
    axes[1].set_ylabel("Hits")
    axes[1].set_title("Run 79 VMM 11 — zoom low ADC")
    axes[1].legend()

    plt.tight_layout()
    plt.show()

    df_vmm = df_hits[df_hits["vmm"] == 11]
    noise_hits = df_vmm[df_vmm["over_threshold"] == 0]["adc"].values

    # Basic stats on the noise population
    median = np.median(noise_hits)
    mad = np.median(np.abs(noise_hits - median))
    robust_sigma = 1.4826 * mad

    print(f"Noise pedestal:")
    print(f"  Median        : {median:.1f} ADC")
    print(f"  Robust sigma  : {robust_sigma:.1f} ADC")
    print(f"  5-sigma cut   : {median + 5 * robust_sigma:.1f} ADC")
    print(f"  Fraction of over_threshold=0 beyond 5-sigma: "
          f"{(noise_hits > median + 5 * robust_sigma).mean() * 100:.2f}%")


    input("Press ENTER to continue...")

    print("\n")
    print(f"\nRun 79 — Noise pedestal per VMM:")

    df_noise_baseline = compute_noise_baseline(df_hits, n_sigma=5)
    print(df_noise_baseline.to_string(index=False))

    input("Press ENTER to continue...")
    df_vmm = df_hits[df_hits["vmm"] == 11]
    signal = df_vmm[df_vmm["over_threshold"] == 1]["adc"].values

    print(f"Signal population — VMM 11, Run 79")
    print(f"  N hits         : {len(signal)}")
    print(f"  Mean ADC       : {signal.mean():.1f}")
    print(f"  Median ADC     : {np.median(signal):.1f}")
    print(f"  Std ADC        : {signal.std():.1f}")
    print(f"  Min / Max      : {signal.min()} / {signal.max()}")
    print(f"  Fraction at max (saturated): "
          f"{(signal == signal.max()).mean() * 100:.2f}%")

    input("Press ENTER to continue...")

    signal = df_vmm[df_vmm["over_threshold"] == 1]["adc"].values

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Left — full range
    axes[0].hist(signal, bins=100, color="steelblue", alpha=0.7)
    axes[0].axvline(signal.mean(), color="red", linestyle="--", label=f"Mean={signal.mean():.0f}")
    axes[0].axvline(np.median(signal), color="orange", linestyle="--", label=f"Median={np.median(signal):.0f}")
    axes[0].set_xlabel("ADC")
    axes[0].set_ylabel("Hits")
    axes[0].set_title("VMM 11 Run 79 — signal full range")
    axes[0].legend()

    # Right — zoom into the rising edge and peak region
    axes[1].hist(signal, bins=100, color="steelblue", alpha=0.7)
    axes[1].set_xlim(100, 600)
    axes[1].axvline(signal.mean(), color="red", linestyle="--", label=f"Mean={signal.mean():.0f}")
    axes[1].axvline(np.median(signal), color="orange", linestyle="--", label=f"Median={np.median(signal):.0f}")
    axes[1].set_xlabel("ADC")
    axes[1].set_ylabel("Hits")
    axes[1].set_title("VMM 11 Run 79 — signal zoom")
    axes[1].legend()

    plt.tight_layout()
    plt.show()

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    plot_idx = 0  # separate counter for axes slots

    for vmm_id in sorted(df_hits["vmm"].unique()):
        signal = df_hits[
            (df_hits["vmm"] == vmm_id) &
            (df_hits["over_threshold"] == 1)
            ]["adc"].values

        if len(signal) < 100:
            continue

        if plot_idx >= len(axes):  # safety guard
            break

        axes[plot_idx].hist(signal, bins=100, color="steelblue", alpha=0.7)
        axes[plot_idx].set_xlim(100, 700)
        axes[plot_idx].set_xlabel("ADC")
        axes[plot_idx].set_ylabel("Hits")
        axes[plot_idx].set_title(f"VMM {vmm_id} (n={len(signal)})")

        plot_idx += 1

    plt.suptitle("Run 79 — Signal distribution per VMM", y=1.02)
    plt.tight_layout()
    plt.show()

    print(f"\nRun 79 — Saturation check per VMM:")
    print(f"{'VMM':>5} {'N_signal':>10} {'N_saturated':>12} {'Sat%':>8} {'ADC_max':>10}")

    for vmm_id in sorted(df_hits["vmm"].unique()):
        signal = df_hits[
            (df_hits["vmm"] == vmm_id) &
            (df_hits["over_threshold"] == 1)
            ]["adc"].values

        if len(signal) < 100:
            continue

        adc_max = signal.max()
        n_sat = (signal == adc_max).sum()
        sat_frac = 100 * n_sat / len(signal)

        print(f"{vmm_id:>5} {len(signal):>10} {n_sat:>12} {sat_frac:>8.2f}% {adc_max:>10}")

        for vmm_id in [11, 13]:
            signal_clean = get_clean_signal(df_hits, vmm_id)
            signal_raw = df_hits[
                (df_hits["vmm"] == vmm_id) &
                (df_hits["over_threshold"] == 1)
                ]["adc"].values

            print(f"\nVMM {vmm_id}")
            print(f"  Raw   — N={len(signal_raw):>8},  Mean={signal_raw.mean():.1f},  "
                  f"Median={np.median(signal_raw):.1f}")
            print(f"  Clean — N={len(signal_clean):>8},  Mean={signal_clean.mean():.1f},  "
                  f"Median={np.median(signal_clean):.1f}")

    for vmm_id in [11, 13]:
        signal_clean = get_clean_signal(df_hits, vmm_id)

        # Get this VMM's noise cut from the baseline table
        noise_cut = df_noise_baseline.loc[
            df_noise_baseline["vmm_id"] == vmm_id, "noise_cut"
        ].iloc[0]

        print(f"\nVMM {vmm_id} — using noise_cut={noise_cut:.1f} as adc_min")

        mpv, counts, smoothed, bin_centers = estimate_mpv(
            signal_clean, adc_min=noise_cut
        )

        scale = counts.max() / smoothed.max()
        fig, ax = plt.subplots()
        ax.bar(bin_centers, counts, width=bin_centers[1] - bin_centers[0],
               color="steelblue", alpha=0.5, label="Clean signal")
        ax.plot(bin_centers, smoothed * scale, color="black",
                linewidth=2, label="Smoothed")
        ax.axvline(noise_cut, color="orange", linestyle="--",
                   label=f"Noise cut={noise_cut:.0f}")
        ax.axvline(mpv, color="red", linestyle="--",
                   label=f"MPV={mpv:.0f}")
        ax.set_xlabel("ADC")
        ax.set_ylabel("Hits")
        ax.set_title(f"VMM {vmm_id} Run 79 — MPV from noise cut")
        ax.legend()
        plt.tight_layout()
        plt.show()

        print(f"  MPV = {mpv:.1f} ADC")

        for run_no in [78, 79]:
            run_dir = get_run_dir(data_dir, run_no)
            root_files = list_root_files(run_dir, n=2)
            file_path = os.path.join(run_dir, root_files[1])

            df_hits = load_hits_root(
                file_path, branches=["adc", "vmm", "ch", "over_threshold"]
            )

            df_noise = compute_noise_baseline(df_hits)

            print(f"\nRun {run_no} — noise baseline table:")
            if df_noise.empty:
                print("  WARNING: no over_threshold=0 hits found — sng=0 run")
            else:
                print(df_noise.to_string(index=False))


    # ---- Plot result ----
    if PLOTS["adc_hist_per_run"]:
        plot_adc_histograms_for_runs(run_num_sg_st_sng, data_dir)

    if PLOTS["adc_hist_separate_vmm"]:
        plot_adc_by_vmm(vmm_ids, run_num_sg_st_sng, df_run_scan, data_dir)

    # plot_mean_vs_peaking(df_results, vmm_ids)
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
                df_hits = load_hits_root(file_path, branches=["adc", "vmm", "ch", "time", "over_threshold"])
                compare_full_vs_cut(df_hits, vmm_id, run_no, adc_cut=100)

    if PLOTS["removed_fraction"]:
        plot_removed_fraction(df_results, vmm_ids)

    if PLOTS["robust_stats"]:
        plot_robust_vs_peaking(df_results, vmm_ids)


if __name__ == "__main__":
    main()


print("bonzo")
