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
    # threshold = median - 1 * rms  # median - 5*sigma
    threshold = median - 3 * rms  # median - 5*sigma
    spark_mask = hist < threshold

    if plot:
        print("Number of bins:", len(hist))

        print("Maximum hit rate per bin:", max_rate)
        print("Some histogram values:", hist[:20])

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
        # plt.show()
        print("Detected sparks:", len(spark_durations))
        print("Spark durations (s):", spark_durations)

    return spark_durations


def main():
    # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_DC"
    # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_AC"
    data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_R"
    # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_no_protection"

    set_root_style()

    all_spark_durations = []
    for file_name in os.listdir(data_dir):
        if file_name.startswith("enp") and file_name.endswith(".root"):
            file_path = os.path.join(data_dir, file_name)
            print("Processing file:", file_path)
            df = load_root_file(file_path, branches=["time"])
            print(df)
            if df.empty:
                print("No data found in file.")
                continue
            spark_durations = find_sparks(np.array(df["time"]), plot=False)
            # if any(spark_durations > 0.5): #AC/DC only

            # find_sparks(np.array(df["time"]), plot=True, bin_width=0.01) default
            find_sparks(np.array(df["time"]), plot=False, bin_width=0.01)
            all_spark_durations.extend(spark_durations)

    fig, ax = plt.subplots()
    ax.hist(all_spark_durations, bins=20, color="blue", alpha=0.7)
    mean = np.mean(all_spark_durations)
    rms = np.std(all_spark_durations)
    ax.axvline(mean, color="red", linestyle="--")
    ax.annotate(f"Sparks: {len(all_spark_durations)}\nMean: {mean:.2f} s\nRMS: {rms:.2f} s",
                xy=(0.75, 0.95), xycoords='axes fraction',
                fontsize=14, va='top',
                bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.5))
    # Make an arrow pointing to the largest bin and label it "Double Sparks!!"
    counts, bin_edges = np.histogram(all_spark_durations, bins=20)
    max_bin_center = 0.5 * (bin_edges[-2] + bin_edges[-1])
    # ax.annotate("Double Sparks!!",
    #             xy=(max_bin_center, counts[-1]),
    #             xytext=(bin_edges[-1] * 0.75, counts[-1] + np.max(counts) * 0.4),
    #             arrowprops=dict(facecolor='black', arrowstyle="->", color="red"),
    #             fontsize=12, fontweight="bold", color="red"
    #             )
    ax.set_xlabel("Spark Duration [s]")
    ax.set_ylabel("Counts")
    ax.set_title("Distribution of Spark Durations")
    plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
    print("bonzo")
