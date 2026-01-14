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
import pandas as pd
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


data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_DC"

file_100 = [
    os.path.join(data_dir, f)
    for f in os.listdir(data_dir)
    if f.startswith("enp")
    and "_100_" in f
    and f.endswith(".root")
]

file_100 = file_100[0]
print(file_100)

df = load_root_file(file_100, branches=["time"])

# Convert to seconds
df["time_s"] = df["time"] / 1e9
t0 = df["time_s"].iloc[0]
time_rel = df["time_s"] - t0
# time_rel = df["time_s"].iloc[1:] - df["time_s"].iloc[0]
time_rel = time_rel[time_rel > 50]
# print(time_rel)
# input()
# time_s = df["time_s"]

# Define bin width (window) — you can tune
bin_width = 0.1  # seconds
bins = np.arange(0, time_rel.max() + bin_width, bin_width)
hist, bin_edges = np.histogram(time_rel, bins=bins)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])


# Count hits per bin
hist, bin_edges = np.histogram(time_rel, bins=bins)

# Mid-points of bins for plotting
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

print("Number of bins:", len(hist))
# print("First few counts:", hist[:10])

max_rate = hist.max()
print("Maximum hit rate per bin:", max_rate)

threshold = 0.9 * max_rate  # can tune to 0.9, 0.8 etc
print("Threshold for spark detection:", threshold)
spark_mask = hist < threshold
print(spark_mask)

######## Plotting
set_root_style()

plt.figure()
plt.hist(
    time_rel,
    bins=300,
    # range=(4405, 4450),
    histtype="step",
    linewidth=1.6,
    color="black",
)
# plt.step(bin_centers, hist, where="mid", color="black", linewidth=1.5, label="Hits per bin")
plt.axhline(y=threshold, color="red", linestyle="--", linewidth=1.5, label=f"{int(threshold)}% threshold")
plt.xlabel("Time [s]")
plt.ylabel("Entries")
plt.title("hits:time")
plt.yscale("log")  # very common in CERN timing plots
plt.tight_layout()
plt.show()





def main():
    pass

if __name__ == "__main__":
    main()
    print("bonzo")
