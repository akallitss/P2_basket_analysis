#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 14/11/2025 19:32
Created in PyCharm
Created as read_vmm_alinx_data.py

@author: akallits
"""

import os
import uproot
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.cluster import DBSCAN

from SquarePadDetector import SquareDetector


def cluster_events(df, time_threshold):
    """
    Create event IDs by grouping hits separated by a time gap > time_threshold.
    time_threshold in same units as df['time'] (e.g. ns).
    """
    df = df.sort_values("time").reset_index(drop=True)

    # Compute time differences between successive hits
    dt = df["time"].diff().fillna(0)

    # Start a new event when dt > threshold
    df["event_id"] = (dt > time_threshold).cumsum()

    return df


def cluster_spatial(df, radius, min_samples=1):
    """
    Spatial DBSCAN clustering per event, with ADC-weighted cluster centroids.
    Adds columns: spatial_cluster_id, cluster_x_mm, cluster_y_mm.
    """
    df = df.copy()
    cluster_ids = []
    cluster_x = []
    cluster_y = []

    for event_id, event_df in df.groupby("event_id"):
        coords = event_df[["X_mm", "Y_mm"]].values
        adcs = event_df["adc"].values

        # DBSCAN for spatial clustering
        db = DBSCAN(eps=radius, min_samples=min_samples).fit(coords)
        labels = db.labels_

        # Build global unique IDs
        global_ids = [event_id * 10000 + (l if l >= 0 else -1) for l in labels]

        # Store temporary results
        cluster_ids.extend(global_ids)

        # For centroid calculation:
        event_df_local = event_df.copy()
        event_df_local["local_cluster"] = labels

        # Compute centroid for each local cluster (excluding noise = -1)
        centroids = {}
        for c in np.unique(labels):
            if c == -1:
                continue  # skip noise hits

            mask = (event_df_local["local_cluster"] == c)
            adc_vals = adcs[mask]
            x_vals = coords[mask, 0]
            y_vals = coords[mask, 1]

            # ADC weighted centroid
            w = adc_vals.sum()
            cx = np.sum(x_vals * adc_vals) / w
            cy = np.sum(y_vals * adc_vals) / w

            centroids[c] = (cx, cy)

        # Assign centroid coordinates to each row
        for lab in labels:
            if lab == -1:
                cluster_x.append(np.nan)
                cluster_y.append(np.nan)
            else:
                cx, cy = centroids[lab]
                cluster_x.append(cx)
                cluster_y.append(cy)

    df["spatial_cluster_id"] = cluster_ids
    df["cluster_x_mm"] = cluster_x
    df["cluster_y_mm"] = cluster_y

    return df


def plot_cluster_2d_hist(df, x_col="cluster_x_mm", y_col="cluster_y_mm", cluster_col="spatial_cluster_id"):
    # Define 1 mm bins
    x_min, x_max = df[x_col].min(), df[x_col].max()
    y_min, y_max = df[y_col].min(), df[y_col].max()

    x_bins = np.arange(x_min, x_max + 1, 1)   # 1 mm
    y_bins = np.arange(y_min, y_max + 1, 1)

    # Plot: 2D histogram of hit density
    plt.figure(figsize=(8, 7))
    plt.hist2d(df[x_col], df[y_col], bins=[x_bins, y_bins])
    plt.xlabel("X [mm]")
    plt.ylabel("Y [mm]")
    plt.title("2D Histogram of Clusters (1 mm bins)")
    plt.colorbar(label="Counts per bin")
    plt.axis("equal")
    plt.tight_layout()

run_number = 'run_67'
print(run_number)
data_dir = f"/home/akallits/Documents/Saclay-PostDoc/SPS_beam_test/data/VMM-alinx_data/{run_number}"

files = os.listdir(data_dir)

enp_files = sorted([os.path.join(data_dir, f)
    for f in os.listdir(data_dir)
    if f.startswith("enp") and f.endswith(".root")])

# print("Selected files:")
# for f in enp_files:
#     print("  ", f)

file = uproot.open(enp_files[0])
tree = file["hits"]
branches = ["hits/id", "hits/det", "hits/plane", "hits/fec", "hits/vmm",
    "hits/readout_time", "hits/time", "hits/geo_id", "hits/ch",
    "hits/pos", "hits/bcid", "hits/tdc", "hits/adc",
    "hits/over_threshold", "hits/chip_time", "hits/event_counter" ]

df_hits = tree.arrays(library="pd")
print(df_hits.head())
#
df_padmap = pd.read_csv("p2_small_detector_map.csv")
# print(df_padmap.head())
## VMM to hybrid mapping
vmm_hybrid_mapping = {
    'trigger': [0, 1],
    'p2_large_1': [12, 13, 14, 15],
    'p2_small_1': [10, 11],
    'p2_small_3': [8, 9],
}

# VMM connector mapping and orientations
vmm_connector_mapping = {
    10: 0,
    11: 1
}

vmm_connector_orientations = {
    10: 'normal',
    11: 'normal',
}

vmm_connector_channel_mapping = {
    'normal': {x: x for x in range(64)},
    # 'reversed': {x: 63 - x for x in range(64)},  # Didn't seem to work
    'inverted': {x: x + 1 if x % 2 == 0 else x - 1 for x in range(64)},
    # 'reversed_inverted': {x: 62 - x if x % 2 == 0 else 64 - x for x in range(64)},  # Didn't check
}

# VMM mapping details combined
vmm_mapping = {
    'trigger': {
        'vmm_ids': [0, 1],
        'connector_ids': [10, 11],
        'orientation': None,
    },
    'p2_large_1': {
        'vmm_ids': [12, 13, 14, 15],
        'connector_ids': [0, 1, 2, 3],
        'orientation': ['normal', 'normal', 'normal', 'normal'],
    },
    'p2_small_1': {
        'vmm_ids': [10, 11],
        'connector_ids': [0, 1],
        'orientation': ['normal', 'normal', 'normal', 'normal'],

    },
    'p2_small_3': {
        'vmm_ids': [8, 9],
        'connector_ids': [8, 9],
        'orientation': ['normal', 'normal', 'normal', 'normal'],
    },
}




detector = 'trigger'
vmm_ids_trigger = vmm_hybrid_mapping[detector]
print(f'Detector: {detector} is associated with VMM IDs: {vmm_ids_trigger}')
# print(df_hits['vmm'])
df_trigger_hits = df_hits[df_hits["vmm"].isin(vmm_ids_trigger)]
print(df_trigger_hits.head())
print(f'Unique channels: {df_trigger_hits["ch"].unique()}')
print(f'Counts per unique channel:\n{df_trigger_hits["ch"].value_counts()}')

for vmm_num in df_trigger_hits['vmm'].unique():
    df_trigger_hits_ch = df_trigger_hits[df_trigger_hits['vmm'] == vmm_num]
    print(f'Counts per unique channel for VMM {vmm_num}:\n{df_trigger_hits_ch["ch"].value_counts()}')


#get the time of the trigger hits
for vmm in df_trigger_hits['vmm'].unique():
    df_trigger_hits_vmm = df_trigger_hits[df_trigger_hits['vmm'] == vmm]
    print(f'Processing VMM {vmm}')
    for ch in df_trigger_hits['ch'].unique():
        print(f' Processing channel {ch}')
        df_trigger_hits_ch = df_trigger_hits[df_trigger_hits['ch'] == ch]
        # df_trigger_hits_ch = df_trigger_hits_vmm[df_trigger_hits['ch'] == ch]
        trigger_times = df_trigger_hits_ch['time'].values
        # print("Trigger times:", trigger_times)
        # if ch == 40:
        #     print(df_trigger_hits_ch[['vmm', 'time']])

        #convert from ns to us
        trigger_times = trigger_times / 1000.0

        # plot histogram of trigger times
        # python
        # ROOT-like histogram: step outline + light fill, inward ticks, stats box, log y-scale

        # Save and tweak rcParams for a ROOT-like look
        # rc_backup = mpl.rcParams.copy()
        # mpl.rcParams.update({
        #     "figure.figsize": (10, 6),
        #     "axes.linewidth": 1.2,
        #     "xtick.direction": "in",
        #     "ytick.direction": "in",
        #     "xtick.top": True,
        #     "ytick.right": True,
        #     "font.size": 12,
        # })
        #
        # fig, ax = plt.subplots()
        #
        # # Draw step outline and a faint filled underlay (ROOT-esque)
        # counts, bins, _ = ax.hist(trigger_times, bins=100, histtype="step", color="black", linewidth=1.5)
        # ax.hist(trigger_times, bins=bins, histtype="stepfilled", color="C0", alpha=0.18)
        #
        # # Labels, title, grid and aspect
        # ax.set_title("Trigger Hit Times Distribution")
        # ax.set_xlabel("Time [s]")
        # ax.set_ylabel("Counts")
        # ax.grid(True, linestyle="--", alpha=0.3)
        #
        # # Use log scale for the y axis (common in ROOT plots)
        # ax.set_yscale("log")
        #
        # # Compute and show simple stats box (Entries, Mean, RMS)
        # entries = len(trigger_times)
        # mean = np.mean(trigger_times)
        # rms = np.sqrt(np.mean((trigger_times - mean) ** 2))
        # stats = f"Entries = {entries}\nMean = {mean:.3e} s\nRMS = {rms:.3e} s"
        # ax.text(0.98, 0.95, stats, transform=ax.transAxes, ha="right", va="top",
        #         bbox=dict(facecolor="white", edgecolor="black", pad=6))
        #
        # plt.tight_layout()
        #
        # # Restore rcParams if desired
        # mpl.rcParams.update(rc_backup)
        # plt.show()

        #find the time difference between each trigger hit and the next trigger hit
        time_differences = np.diff(trigger_times)
        # print("Time differences between trigger hits:", time_differences)
        #plot histogram of time differences
        # python
        # ROOT-like time-difference histogram (step + filled, inward ticks, stats box, log y-scale)
        rc_backup = mpl.rcParams.copy()
        mpl.rcParams.update({
            "figure.figsize": (10, 6),
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "font.size": 12,
        })

        fig, ax = plt.subplots()

        # Step outline + faint filled underlay (ROOT-esque)
        bins = np.linspace(-0.1, 1000, 50)
        counts, bins, _ = ax.hist(time_differences, bins=bins, histtype="step", color="black", linewidth=1.5)
        ax.hist(time_differences, bins=bins, histtype="stepfilled", color="C0", alpha=0.18)

        # Labels, title, grid
        ax.set_title(f"Time Differences Between Trigger Hits for Channel {ch}, on VMM {vmm}")
        ax.set_xlabel("Time Difference [us]")
        ax.set_ylabel("Counts")
        ax.grid(True, linestyle="--", alpha=0.3)

        # Log scale for the y axis (common in ROOT plots)
        ax.set_yscale("log")

        # Compute and show simple stats box (Entries, Mean, RMS) — handle empty arrays safely
        if len(time_differences) > 0:
            entries = len(time_differences)
            mean = np.nanmean(time_differences)
            rms = np.sqrt(np.nanmean((time_differences - mean) ** 2))
            stats = f"Entries = {entries}\nMean = {mean:.3e} us\nRMS = {rms:.3e} us"
        else:
            stats = "No entries"

        ax.text(0.98, 0.95, stats, transform=ax.transAxes, ha="right", va="top",
                bbox=dict(facecolor="white", edgecolor="black", pad=6))

        plt.tight_layout()
        mpl.rcParams.update(rc_backup)
plt.show()




detector = 'p2_small_1'  # Change this to analyze a different detector
vmm_ids = vmm_hybrid_mapping[detector]
print(f'Detector: {detector} is associated with VMM IDs: {vmm_ids}')
# print(df_hits['vmm'])
df_det_hits = df_hits[df_hits["vmm"].isin(vmm_ids)]
# print(df_det_hits.head())
# print(df_det_hits[["vmm", "ch", "adc"]])
print(df_hits['time'])
input()
# Invert: connector → vmm
connector_to_vmm = {v: k for k, v in vmm_connector_mapping.items()}

# Add real VMM number to padmap
df_padmap["vmm"] = df_padmap["connector"].map(connector_to_vmm)

print(df_padmap[["vmm", "channel"]])

# Correct channel numbering based on connector orientation
def correct_channel(row):
    vmm = row["vmm"]
    channel = row["channel"]
    orientation = vmm_connector_orientations[vmm]
    channel_map = vmm_connector_channel_mapping[orientation]
    return channel_map[channel]

df_padmap["channel_cor"] = df_padmap.apply(correct_channel, axis=1)

print(df_padmap[["vmm", "channel", "channel_cor"]])
# input()

# Merge hits with padmap using vmm + channel
df_with_xy = df_det_hits.merge(
    df_padmap,
    left_on=["vmm", "ch"],
    right_on=["vmm", "channel_cor"],
    how="left"
)

print(df_with_xy.head())

# Keep only hits with a valid pad match
df_xy = df_with_xy.dropna(subset=["X_mm", "Y_mm", "adc", "time"]).copy()

# print("For clustering")
# print(df_with_xy[["X_mm", "Y_mm", "adc", "time"]])
#

#
# df = df_xy.copy()
#
# print("Clustering hits into events and spatial clusters...")
# # 1) Cluster into events
# df = cluster_events(df, time_threshold=1e6)   # example: 1 microsecond
#
# # 2) Cluster spatially inside each event
# df = cluster_spatial(df, radius=20)            # radius in mm
#
# print(df[["adc", "X_mm", "Y_mm", "spatial_cluster_id", "cluster_x_mm", "cluster_y_mm"]].head(20))
#
# # Plot scatter plot of cluster centroids, one per cluster_id
# fig, ax = plt.subplots(figsize=(8, 7))
# sc = ax.scatter(
#     df["cluster_x_mm"],
#     df["cluster_y_mm"],
#     c=df["spatial_cluster_id"],
#     cmap="tab20",
#     s=50,
#     alpha=0.7,
#     edgecolors='k'
# )
# ax.set_xlabel("Cluster X [mm]")
# ax.set_ylabel("Cluster Y [mm]")
# ax.set_title("Spatial Clusters Colored by Cluster ID")
# ax.set_aspect("equal", adjustable="box")
# plt.colorbar(sc, label="Spatial Cluster ID")
# plt.tight_layout()
#
# # Plot histogram of cluster size -- number of hits per spatial cluster
# cluster_sizes = df.groupby("spatial_cluster_id").size()
# plt.figure(figsize=(8, 5))
# plt.hist(cluster_sizes, bins=np.arange(1, cluster_sizes.max() + 2) - 0.5, alpha=0.7)
# plt.xlabel("Number of Hits per Spatial Cluster")
# plt.ylabel("Counts")
# plt.title("Histogram of Spatial Cluster Sizes")
# plt.grid(True)
# plt.tight_layout()
#
# plot_cluster_2d_hist(df)
# plt.show()



# Count number of hits per pad coordinate
hit_counts = (
    df_xy.groupby(["X_mm", "Y_mm"])
         .size()
         .reset_index(name="hit_count")
)
# print(hit_counts.head())

# Merge counts back into df_xy so each hit knows its pad's hit count
df_xy = df_xy.merge(hit_counts, on=["X_mm", "Y_mm"], how="left")

detector_model = SquareDetector()
detector_model.read_mapping("p2_small_detector_map.csv")
detector_model.plot_hit_heatmap(df_xy, cmap="jet", global_coords=False, log_scale=True)
detector_model.plot_hit_heatmap(df_xy, cmap="jet", global_coords=False, area_norm=True, log_scale=True)
detector_model.plot_hit_heatmap(df_xy, cmap="jet", global_coords=False, adc_weighted=True, log_scale=True)
detector_model.plot_detector()
# plt.show()


####################### PLOTTING #########################

# #plot adc distribution from df_hits
# plt.figure(figsize=(10, 6))
# plt.hist(df_hits['adc'], bins=200, range=(0, 4096), alpha=0.7)
# plt.xlim(0, 1200)
# plt.title('ADC Distribution')
# plt.xlabel('ADC Value')
# plt.ylabel('Counts')
# plt.grid(True)
# # plt.show()
#
#
# #plot vmm distribution from df_hits
# plt.figure(figsize=(10, 6))
# plt.hist(df_hits['vmm'], bins=np.arange(df_hits['vmm'].min(), df_hits['vmm'].max() + 2) - 0.5, alpha=0.7)
# plt.title('VMM Distribution')
# plt.xlabel('VMM Value')
# plt.ylabel('Counts')
# plt.grid(True)
#
#
# # Plot channel vs adc
# plt.figure(figsize=(12, 6))
# plt.scatter(df_hits['ch'], df_hits['adc'], alpha=0.5, s=1)
# plt.title('Channel vs ADC')
# plt.xlabel('Channel')
# plt.ylabel('ADC Value')
# plt.grid(True)
# # plt.show()

for vmm in vmm_ids:
    df_vmm_det = df_hits[df_hits['vmm'] == vmm]
    plt.figure(figsize=(12, 6))
    plt.hist2d(df_vmm_det['ch'], df_vmm_det['adc'], bins=[50, 200], range=[[0, 65], [0, 1200]], cmap='viridis', cmin=1, norm=LogNorm())
    plt.title(f'2D Histogram of Channel vs ADC for VMM {vmm}')
    plt.xlabel('Channel')
    plt.ylabel('ADC Value')
    plt.grid(True)

# -------------/-----------------------------------------
# Scatter Plot
# ------------------------------------------------------
# plt.figure(figsize=(8, 6))
# sc = plt.scatter(
#     df_xy["X_mm"],
#     df_xy["Y_mm"],
#     c=df_xy["hit_count"],
#     marker='s',
#     s=200*df_xy["Size_mm"],                 # marker size
#     cmap="viridis"       # or "plasma", "turbo", "inferno"
# )
#
# plt.colorbar(sc, label="Number of hits at pad")
# plt.xlabel("X (mm)")
# plt.ylabel("Y (mm)")
# plt.title("Hit Occupancy Scatter Plot (Color = Hit Count)")
# plt.gca().set_aspect("equal", adjustable="box")
#
# plt.tight_layout()
# plt.show()

plt.show()

#mapping for Detector P2 on SPS beam test Nov 2025
#0,7,6,5,4 hybrids
#1,2,13,14,11,12,10,11,9,8 vmms
#trigger, large p2 card1, large p2 card2, small p2-1, small p2-3

