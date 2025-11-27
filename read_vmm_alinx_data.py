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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from sklearn.cluster import DBSCAN
from SquarePadDetector import SquareDetector

#######################
# DATA LOADING
#######################
def load_root_file(file_path, branches=None):
    """Load ROOT file into pandas DataFrame."""
    file = uproot.open(file_path)
    tree = file["hits"]
    return tree.arrays(library="pd")


def load_padmap(csv_path):
    """Load pad map CSV into DataFrame."""
    return pd.read_csv(csv_path)


def map_vmm_channels(df_padmap, connector_to_vmm, connector_orientations, channel_mappings):
    """Add corrected VMM channel numbers to padmap."""
    df_padmap["vmm"] = df_padmap["connector"].map(connector_to_vmm)

    def correct_channel(row):
        orientation = connector_orientations[row["vmm"]]
        mapping = channel_mappings[orientation]
        return mapping[row["channel"]]

    df_padmap["channel_cor"] = df_padmap.apply(correct_channel, axis=1)
    return df_padmap


def merge_hits_with_padmap(df_hits, df_padmap):
    """Merge hit data with padmap to get X/Y coordinates."""
    df_with_xy = df_hits.merge(
        df_padmap,
        left_on=["vmm", "ch"],
        right_on=["vmm", "channel_cor"],
        how="left"
    )
    return df_with_xy.dropna(subset=["X_mm", "Y_mm", "adc", "time"]).copy()


#######################
# CLUSTERING
#######################
def cluster_events(df, time_threshold):
    """Cluster hits into events by time gap."""
    df = df.sort_values("time").reset_index(drop=True)
    dt = df["time"].diff().fillna(0)
    df["event_id"] = (dt > time_threshold).cumsum()
    return df


def cluster_spatial(df, radius, min_samples=1):
    """Spatial DBSCAN clustering per event with ADC-weighted centroids."""
    df = df.copy()
    cluster_ids, cluster_x, cluster_y = [], [], []

    for event_id, event_df in df.groupby("event_id"):
        coords = event_df[["X_mm", "Y_mm"]].values
        adcs = event_df["adc"].values
        db = DBSCAN(eps=radius, min_samples=min_samples).fit(coords)
        labels = db.labels_
        global_ids = [event_id * 10000 + (l if l >= 0 else -1) for l in labels]
        cluster_ids.extend(global_ids)

        event_df_local = event_df.copy()
        event_df_local["local_cluster"] = labels

        centroids = {}
        for c in np.unique(labels):
            if c == -1:
                continue
            mask = (event_df_local["local_cluster"] == c)
            w = adcs[mask].sum()
            centroids[c] = (np.sum(coords[mask,0]*adcs[mask])/w,
                            np.sum(coords[mask,1]*adcs[mask])/w)

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


#######################
# PLOTTING FUNCTIONS
#######################
def plot_adc_distribution(df_hits, max_adc=1200, bins=200):
    plt.figure(figsize=(10,6))
    plt.hist(df_hits['adc'], bins=bins, range=(0,4096), alpha=0.7)
    plt.xlim(0, max_adc)
    plt.title('ADC Distribution')
    plt.xlabel('ADC Value')
    plt.ylabel('Counts')
    plt.grid(True)
    plt.show()


def plot_vmm_distribution(df_hits):
    plt.figure(figsize=(10,6))
    plt.hist(df_hits['vmm'], bins=np.arange(df_hits['vmm'].min(), df_hits['vmm'].max()+2)-0.5, alpha=0.7)
    plt.title('VMM Distribution')
    plt.xlabel('VMM Value')
    plt.ylabel('Counts')
    plt.grid(True)
    plt.show()


def plot_channel_vs_adc(df_hits):
    plt.figure(figsize=(12,6))
    plt.scatter(df_hits['ch'], df_hits['adc'], alpha=0.5, s=1)
    plt.title('Channel vs ADC')
    plt.xlabel('Channel')
    plt.ylabel('ADC Value')
    plt.grid(True)
    plt.show()


def plot_2d_hist_per_vmm(df_hits, vmm_ids, x='ch', y='adc', bins_x=50, bins_y=200, x_range=[0,65], y_range=[0,1200]):
    for vmm in vmm_ids:
        df_vmm = df_hits[df_hits['vmm'] == vmm]
        plt.figure(figsize=(12,6))
        plt.hist2d(df_vmm[x], df_vmm[y], bins=[bins_x,bins_y], range=[x_range,y_range],
                   cmap='viridis', cmin=1, norm=LogNorm())
        plt.title(f'2D Histogram of {x} vs {y} for VMM {vmm}')
        plt.xlabel(x)
        plt.ylabel(y)
        plt.grid(True)
        plt.show()


def plot_hit_scatter(df_xy):
    plt.figure(figsize=(8,6))
    sc = plt.scatter(
        df_xy["X_mm"],
        df_xy["Y_mm"],
        c=df_xy["hit_count"],
        marker='s',
        s=200*df_xy["Size_mm"],
        cmap="viridis"
    )
    plt.colorbar(sc, label="Number of hits at pad")
    plt.xlabel("X (mm)")
    plt.ylabel("Y (mm)")
    plt.title("Hit Occupancy Scatter Plot (Color = Hit Count)")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.show()


#######################
# TRIGGER TIME DIFFERENCES
#######################
def plot_trigger_time_differences(df_hits, vmm_ids):
    rc_backup = mpl.rcParams.copy()
    mpl.rcParams.update({
        "figure.figsize": (10,6),
        "axes.linewidth": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "font.size": 12,
    })

    for vmm in vmm_ids:
        df_vmm = df_hits[df_hits['vmm'] == vmm]
        for ch in df_vmm['ch'].unique():
            df_ch = df_vmm[df_vmm['ch'] == ch]
            trigger_times = df_ch['time'].values / 1000.0  # ns -> us
            if len(trigger_times) < 2:
                continue
            time_diff = np.diff(trigger_times)
            if len(time_diff) == 0:
                continue

            entries = len(time_diff)
            mean = np.nanmean(time_diff)
            rms = np.sqrt(np.nanmean((time_diff - mean)**2))

            fig, ax = plt.subplots()
            bins = np.linspace(-0.1, 1000, 50)
            ax.hist(time_diff, bins=bins, histtype="step", color="black", linewidth=1.5)
            ax.hist(time_diff, bins=bins, histtype="stepfilled", color="C0", alpha=0.18)
            ax.set_title(f"Time Differences Between Trigger Hits (VMM {vmm}, Ch {ch})")
            ax.set_xlabel("Time Difference [us]")
            ax.set_ylabel("Counts")
            ax.set_yscale("log")
            ax.grid(True, linestyle="--", alpha=0.3)
            stats = f"Entries = {entries}\nMean = {mean:.3e} us\nRMS = {rms:.3e} us"
            ax.text(0.98, 0.95, stats, transform=ax.transAxes, ha="right", va="top",
                    bbox=dict(facecolor="white", edgecolor="black", pad=6))
            plt.tight_layout()
            plt.show()

    mpl.rcParams.update(rc_backup)


#######################
# MAIN FUNCTION
#######################
def main():
   #######################
    # VMM MAPPINGS
    #######################
    vmm_hybrid_mapping = {
        'trigger': [0,1],
        'p2_large_1': [12,13,14,15],
        'p2_small_1': [10,11],
        'p2_small_3': [8,9],
    }

    vmm_connector_mapping = {10:0, 11:1}
    vmm_connector_orientations = {10:'normal', 11:'normal'}
    vmm_connector_channel_mapping = {
        'normal': {x:x for x in range(64)},
        'inverted': {x:x+1 if x%2==0 else x-1 for x in range(64)}
    }

   #######################
   # CONFIG OPTIONS
   #######################


    run_number = 'run_67'
    # detector = 'p2_small_1'
    detector = 'p2_small_3'
    # detector = 'p2_large_1'

    # PLOT FLAGS
    plot_adc = True
    plot_vmm = True
    plot_ch_vs_adc = True
    plot_2d_hist_all_vmm = True
    plot_hit_scatterplot = True
    plot_detector_heatmaps = True
    plot_trigger_times = True
    #######################
    # LOAD DATA
    #######################
    data_dir = f"/home/akallits/Documents/Saclay-PostDoc/SPS_beam_test/data/VMM-alinx_data/{run_number}"
    enp_files = sorted([os.path.join(data_dir, f) for f in os.listdir(data_dir)
                        if f.startswith("enp") and f.endswith(".root")])
    df_hits = load_root_file(enp_files[0])
    df_padmap = load_padmap("p2_small_detector_map.csv")

    df_padmap = map_vmm_channels(
        df_padmap,
        connector_to_vmm={v:k for k,v in vmm_connector_mapping.items()},
        connector_orientations=vmm_connector_orientations,
        channel_mappings=vmm_connector_channel_mapping
    )

    vmm_ids = vmm_hybrid_mapping[detector]
    df_det_hits = df_hits[df_hits["vmm"].isin(vmm_ids)]
    df_xy = merge_hits_with_padmap(df_det_hits, df_padmap)

    # Count hits per pad
    hit_counts = df_xy.groupby(["X_mm","Y_mm"]).size().reset_index(name="hit_count")
    df_xy = df_xy.merge(hit_counts, on=["X_mm","Y_mm"], how="left")

    #######################
    # PLOTTING
    #######################
    if plot_adc:
        plot_adc_distribution(df_hits)

    if plot_vmm:
        plot_vmm_distribution(df_hits)

    if plot_ch_vs_adc:
        plot_channel_vs_adc(df_hits)

    if plot_2d_hist_all_vmm:
        plot_2d_hist_per_vmm(df_hits, vmm_ids)

    if plot_hit_scatterplot:
        plot_hit_scatter(df_xy)

    if plot_detector_heatmaps:
        detector_model = SquareDetector()
        detector_model.read_mapping("p2_small_detector_map.csv")
        detector_model.plot_hit_heatmap(df_xy, cmap="jet", global_coords=False, log_scale=True)
        detector_model.plot_hit_heatmap(df_xy, cmap="jet", global_coords=False, area_norm=True, log_scale=True)
        detector_model.plot_hit_heatmap(df_xy, cmap="jet", global_coords=False, adc_weighted=True, log_scale=True)
        detector_model.plot_detector()

    if plot_trigger_times:
        trigger_vmm_ids = vmm_hybrid_mapping['trigger']
        plot_trigger_time_differences(df_hits, trigger_vmm_ids)

    if plot_trigger_times:
        trigger_vmm_ids = vmm_hybrid_mapping['p2_large_1']
        plot_trigger_time_differences(df_hits, trigger_vmm_ids)


if __name__ == "__main__":
    main()

#mapping for Detector P2 on SPS beam test Nov 2025
#0,7,6,5,4 hybrids
#1,2,13,14,11,12,10,11,9,8 vmms
#trigger, large p2 card1, large p2 card2, small p2-1, small p2-3

