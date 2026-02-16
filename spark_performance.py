#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 14/01/2026 15:18
Created in PyCharm
Created as spark_performance.py

@author: akallits
"""


import os
import uproot
import numpy as np
# import pandas as pd
import matplotlib.pyplot as plt

def load_root_file(file_path, branches=None):
    """Load selected branches from ROOT file into pandas DataFrame."""
    with uproot.open(file_path) as file:
        tree = file["hits"]
        return tree.arrays(branches, library="pd")

def set_root_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 14,
        "axes.linewidth": 1.5,
        "axes.labelsize": 16,
        "axes.titlesize": 16,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 7,
        "ytick.major.size": 7,
        "xtick.major.width": 1.3,
        "ytick.major.width": 1.3,
        "xtick.minor.size": 4,
        "ytick.minor.size": 4,
        "xtick.minor.width": 1.0,
        "ytick.minor.width": 1.0,
        "legend.frameon": False,
        "figure.figsize": (7, 5),
    })


def find_sparks(hit_times, plot=False, bin_width = 0.01):
    """
    Find sparks from binning hit times into hit frequencies. Look for drops.
    :param hit_times:
    :param plot:
    :return:
    """

    hit_times = hit_times / 1e9  # Convert from ns to s
    t0 = hit_times[0]
    time_rel = hit_times - t0
    time_rel = time_rel[time_rel > 20]  # Ignore first 20 seconds
    # time_rel = time_rel[time_rel > .1]  # Ignore first 1 second
    t0 = time_rel[0]

    # Define bin width (window)

    bins = np.arange(t0, time_rel.max() + bin_width, bin_width)
    hist, bin_edges = np.histogram(time_rel, bins=bins)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    # Count hits per bin
    hist, bin_edges = np.histogram(time_rel, bins=bins)
    max_rate = hist.max()
    median = np.median(hist)
    rms = np.std(hist)
    # threshold = median - 2.1 * rms  # median - 1*sigma
    threshold = median - 3 * rms  # median - 3*sigma
    spark_mask = hist < threshold

    if plot:
        print("Number of bins:", len(hist))

        print("Maximum hit rate per bin:", max_rate)
        # print("Some histogram values:", hist[:20])

        print("Median hit rate per bin:", median)
        print("RMS of hit rate per bin:", rms)

        # threshold = 0.9 * max_rate  # 95% of max rate
        # threshold = 1  # median - 5*sigma
        print("Threshold for spark detection:", threshold)
        print(spark_mask)

    # Cluster consecutive bins below threshold and pick out the first bin of each cluster, then find the end of each
    # cluster as the first above threshold after a below threshold
    spark_starts = []
    spark_ends = []
    in_spark = False
    for i in range(len(spark_mask)):
        if spark_mask[i] and not in_spark:
            # Start of a new spark
            spark_starts.append(bin_edges[i])
            in_spark = True
        elif hist[i] >= median and in_spark:
            # End of the current spark
            spark_ends.append(bin_edges[i])
            in_spark = False

    # Handle case where spark goes till the end
    if in_spark:
        spark_starts.pop(-1)
    spark_durations = np.array([end - start for start, end in zip(spark_starts, spark_ends)])

    live_time = time_rel.max() - time_rel.min()

    ######## Plotting
    if plot:
        plt.figure()
        # plt.hist(
        #     time_rel,
        #     bins=bins,
        #     # range=(4405, 4450),
        #     histtype="step",
        #     linewidth=1.6,
        #     color="black",
        # )
        plt.step(bin_centers, hist, where="mid", color="black", linewidth=1.5, label="Hits per bin")
        plt.axhline(y=threshold, color="red", linestyle="--", linewidth=1.5, label=f"{int(threshold)}% threshold")
        plt.axhline(y=median, color="salmon", linestyle="--", linewidth=1.5, label="Median")
        plt.legend()

        for start, end in zip(spark_starts, spark_ends):
            plt.axvspan(start, end, color="red", alpha=0.3)
        plt.xlabel("Time [s]")
        plt.ylabel("Entries")
        plt.title("hits:time")
        # plt.yscale("log")  # very common in CERN timing plots
        plt.tight_layout()
        plt.show()
        print("Detected sparks:", len(spark_durations))
        print("Spark durations (s):", spark_durations)

    return spark_durations, live_time

def rate_distribution_diagnostic(hit_times, bin_width=0.01, plot=True, plot_rate=True):
    """
    Diagnostic (pre-spark):
      - occupancy per bin (entries/bin)
      - instantaneous hit rate per bin (Hz) = entries/bin_width
      - global hit rate = N_hits / live_time
    Uses the same 20 s cut as find_sparks().
    """

    hit_times = hit_times / 1e9  # ns → s
    t0 = hit_times[0]
    time_rel = hit_times - t0
    time_rel = time_rel[time_rel > 20]  # ignore first 20 s

    if len(time_rel) < 2:
        # Avoid crashes on tiny arrays
        return np.array([]), 0.0, 0, np.array([])

    bins = np.arange(time_rel.min(), time_rel.max() + bin_width, bin_width)
    hist, bin_edges = np.histogram(time_rel, bins=bins)

    live_time = time_rel.max() - time_rel.min()
    n_hits = len(time_rel)

    rate_per_bin = hist / bin_width  # Hz, instantaneous rate in each time bin
    #convert rate_per_bin in kHz
    rate_per_bin = rate_per_bin/1e3 #kHz

    global_hit_rate = n_hits / live_time if live_time > 0 else np.nan
    global_hit_rate_kHz = global_hit_rate/1e3

    if plot:
        plt.figure()
        plt.hist(hist, bins=100, histtype="step", linewidth=1.6)
        plt.xlabel("Entries per time bin")
        plt.ylabel("Number of bins")
        plt.yscale("log")
        plt.title("Distribution of time-bin occupancies")
        plt.tight_layout()
        plt.show()

        print("Rate distribution diagnostics (after 20 s cut):")
        print(f"  N hits        = {n_hits}")
        print(f"  Live time     = {live_time:.3f} s")
        print(f"  Global hit rate = {global_hit_rate:.3e} Hz")
        print(f"  Global hit rate in kHz = {global_hit_rate_kHz:.3e} kHz")
        print(f"  Median(occ)   = {np.median(hist):.2f}")
        print(f"  RMS(occ)      = {np.std(hist):.2f}")
        print(f"  Min/Max(occ)  = {hist.min()} / {hist.max()}")

    if plot and plot_rate:
        # New: instantaneous hit-rate distribution in Hz
        plt.figure()
        plt.hist(rate_per_bin, bins=100, histtype="step", linewidth=1.6)
        plt.xlabel("Instantaneous hit rate per bin [kHz]")
        plt.ylabel("Number of bins")
        plt.yscale("log")
        plt.title("Distribution of instantaneous hit rate")
        plt.tight_layout()
        plt.show()

    return hist, live_time, n_hits, rate_per_bin




def main():

    # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_no_protection_jan26"
    # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_no_protection_dec25"
    # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_AC_dec25"
    # data_dir = "/drf/projets/clas12/P2/spark_tests/sparks_test15" #without protection
    # data_dir = "/drf/projets/clas12/P2/spark_tests/sparks_test16" #only AC
    # data_dir = "/drf/projets/clas12/P2/spark_tests/sparks_test17" #without protection
    # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_AC_jan26"
    # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_DC_dec25"
    #data_dir = "/drf/projets/clas12/P2/spark_tests/sparks_test11" #AC-DC
    # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_R_dec25"
    data_dir = "/drf/projets/clas12/P2/spark_tests/sparks_test14" #AC-R


    set_root_style()


    # --- Hit-rate accumulators ---
    total_hits = 0
    total_hit_live_time = 0.0
    all_rate_per_bin = []

    # --- Spark accumulators ---
    total_sparks = 0
    total_live_time = 0.0
    all_spark_durations = []

    for file_name in sorted(os.listdir(data_dir)):
        if not (file_name.startswith("enp") and file_name.endswith(".root")):
            continue

        file_path = os.path.join(data_dir, file_name)
        print("Processing file:", file_path)

        df = load_root_file(file_path, branches=["time"])
        if df.empty:
            continue

        hit_times = np.array(df["time"])

        # ---- PRE-SPARK DIAGNOSTIC (per file) ----
        hist, hit_live_time, n_hits, rate_per_bin = rate_distribution_diagnostic(
            hit_times,
            bin_width=0.01,
            plot=True,
            plot_rate=True
        )

        total_hits += n_hits
        total_hit_live_time += hit_live_time
        if len(rate_per_bin) > 0:
            all_rate_per_bin.extend(rate_per_bin)

        # ---- SPARK FINDING (per file) ----
        spark_durations, live_time = find_sparks(
            hit_times,
            plot=True,
            bin_width=0.01
        )

        # accumulate spark totals
        total_sparks += len(spark_durations)
        total_live_time += live_time
        all_spark_durations.extend(spark_durations)

        # (optional) per-file spark rate print (correct place is inside loop)
        if live_time > 0:
            file_spark_rate = len(spark_durations) / live_time
            print(f"{file_name}: spark rate = {file_spark_rate:.3e} Hz")

    # ---- GLOBAL HIT RATE OVER ALL FILES (print once) ----
    if total_hit_live_time > 0:
        global_hit_rate_total = total_hits / total_hit_live_time
        global_hit_rate_total_kHz = global_hit_rate_total / 1e3
        global_hit_rate_err = np.sqrt(total_hits) / total_hit_live_time
        global_hit_rate_err_kHz = global_hit_rate_err / 1e3
        print("\n=== Hit-rate summary (after 20 s cut) ===")
        print(f"Total hits: {total_hits}")
        print(f"Total hit live time: {total_hit_live_time:.1f} s")
        print(f"Global hit rate: ({global_hit_rate_total:.3e} ± {global_hit_rate_err:.1e}) Hz")
        print(f"Global hit rate in kHz: ({global_hit_rate_total_kHz} ± {global_hit_rate_err:.1e}) kHz ")
    else:
        print("\nNo valid hit live time to compute global hit rate.")

    # ---- INSTANTANEOUS HIT RATE DISTRIBUTION OVER ALL FILES (plot once) ----
    if len(all_rate_per_bin) > 0:
        plt.figure()
        plt.hist(all_rate_per_bin, bins=120, histtype="step", linewidth=1.6)
        plt.xlabel("Instantaneous hit rate per bin [Hz]")
        plt.ylabel("Number of bins")
        plt.yscale("log")
        plt.title("Instantaneous hit rate distribution (all files)")
        plt.tight_layout()
        plt.show()

    # ---- SPARK RATE SUMMARY (global) ----
    if total_live_time > 0:
        sparking_rate = total_sparks / total_live_time
        rate_err = np.sqrt(total_sparks) / total_live_time
        print("\n=== Spark summary ===")
        print(f"Total sparks: {total_sparks}")
        print(f"Total live time: {total_live_time:.1f} s")
        print(f"Sparking rate: ({sparking_rate:.3e} ± {rate_err:.1e}) Hz")
    else:
        print("\nNo valid live time to compute spark rate.")

    # ---- Spark duration distribution ----
    if len(all_spark_durations) > 0:
        fig, ax = plt.subplots()
        ax.hist(all_spark_durations, bins=20, color="blue", alpha=0.7)
        mean = np.mean(all_spark_durations)
        rms = np.std(all_spark_durations)
        ax.axvline(mean, color="red", linestyle="--")
        ax.annotate(
            f"Sparks: {len(all_spark_durations)}\n"
            f"Mean: {mean:.2f} s\n"
            f"RMS: {rms:.2f} s\n",
            xy=(0.75, 0.95),
            xycoords="axes fraction",
            fontsize=14,
            va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.5),
        )
        ax.set_xlabel("Spark Duration [s]")
        ax.set_ylabel("Counts")
        ax.set_title("Distribution of Spark Durations")
        plt.tight_layout()
        plt.show()
    else:
        print("\nNo sparks found → spark duration histogram skipped.")


if __name__ == "__main__":
    main()
    print("bonzo")
