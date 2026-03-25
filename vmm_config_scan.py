#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 19/11/2025 17:37
Created in PyCharm
Created as vmm_config_scan.py

@author: akallits
"""
import os
import pandas as pd
from vmm_mapping import vmm_mapping
import numpy as np
from vmm_io import (load_hits_root, load_run_table, filter_runs, get_run_dir, list_root_files, get_run_groups)
from vmm_noise import compute_noise_baseline, compute_channel_noise_sigma
from vmm_snr import compute_snr, compute_snr_per_channel, summarise_best_config, compute_adc_stats
from vmm_qa import (qa_over_threshold_split,
                    qa_noise_pedestal_stability,
                    qa_signal_distributions,
                    qa_noise_quality_check,
                    qa_mpv_estimation,
                    qa_channel_noise,
                    qa_adc16_artifact,
                    qa_noise_sigma_distribution,
                    qa_robust_vs_std_comparison,
                    qa_mpv_vs_median_comparison)

from vmm_plots import (plot_adc_histograms,
                       plot_adc_histograms_for_runs,
                       plot_adc_by_vmm,
                       compare_full_vs_cut,
                       plot_removed_fraction,
                       plot_mean_vs_peaking,
                       plot_std_vs_peaking,
                       plot_robust_vs_peaking,
                       plot_snr_vs_peaking,
                       plot_snr_vs_gain,
                       plot_snr_heatmap,
                       plot_snr_channel_heatmap_per_vmm,
                       plot_snr_channel_heatmap_all_configs,
                       plot_snr_channel_uniformity)
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
    "mean_vs_peaking": False,
    "debug_samples": False,
    "plot_std_vs_peaking": False,
    "compare_full_vs_cut": False,
    "removed_fraction": False,
    "robust_stats": False,
    "qa_over_threshold_split": False,
    "qa_noise_pedestal_stability": False,
    "qa_mpv_estimation": False,
    "qa_channel_noise": False,
    "qa_adc16_artifact": False,
    "qa_signal_distributions": False,
    "qa_noise_quality_check": False,
    "qa_noise_sigma_distribution": False,
    "qa_robust_vs_std_comparison": False,
    "qa_mpv_vs_median_comparison": False,
    "snr_vs_peaking": False,
    "snr_vs_gain":    False,
    "snr_heatmap":    False,
    "snr_channel_heatmap"     : False,
    "snr_channel_uniformity"  : False,
    "snr_channel_heatmap_per_vmm": False,

}

# -----------------------------
# ANALYSIS & PLOTTING FUNCTIONS
# -----------------------------


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


# -----------------------------
# MAIN LOGIC
# -----------------------------
def main():

    # ---- Load run metadata ----
    cnfg_dir = ("/drf/projets/clas12/P2/akallits/")
    df_run_scan = load_run_table(f"{cnfg_dir}vmm_config_scan.csv")

    data_dir = "/drf/projets/clas12/cern_202511_p2_alinx/"


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


    df_snr_ch = compute_snr_per_channel(
        data_dir=data_dir,
        pairs=run_groups["pairs"],
        detector_vmms=detector_vmms
    )

    df_snr_ch.to_csv("vmm_snr_per_channel.csv", index=False)

    # print(f"\nChannel-level SNR: {len(df_snr_ch)} entries")
    # print(f"SNR range : {df_snr_ch['snr_ch'].min():.1f} "
    #       f"— {df_snr_ch['snr_ch'].max():.1f}")
    # print(f"Mean SNR  : {df_snr_ch['snr_ch'].mean():.1f}")
    # print(f"Median SNR: {df_snr_ch['snr_ch'].median():.1f}")
    # print(f"\nPercentiles:")
    # for p in [1, 5, 25, 50, 75, 95, 99]:
    #     print(f"  {p:>3}th : "
    #           f"{np.percentile(df_snr_ch['snr_ch'], p):.1f}")



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

    df_snr_ch_uncut = compute_snr_per_channel(
        data_dir=data_dir,
        pairs=run_groups["pairs"],
        detector_vmms=detector_vmms,
        min_noise_hits=5,  # almost no cut
        min_signal_hits=10,  # almost no cut
        min_sigma=0.0,  # keep everything including stuck
        max_sigma=999.0,  # keep everything including noisy
        mpv_min=0.0,  # keep everything
        mpv_max=1023.0  # keep everything
    )

    df_snr_ch_uncut.to_csv("vmm_snr_per_channel_uncut.csv", index=False)

    # print(f"Uncut channel entries: {len(df_snr_ch_uncut)}")
    # print(f"SNR range: {df_snr_ch_uncut['snr_ch'].min():.1f} "
    #       f"— {df_snr_ch_uncut['snr_ch'].max():.1f}")

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

    # if PLOTS["snr_channel_heatmap_per_vmm"]:
    #     plot_snr_channel_heatmap_per_vmm(df_snr_ch, detector_vmms)
    if PLOTS["snr_channel_heatmap_per_vmm"]:
        # Version 1 — all channels shown, quality marked with X
        plot_snr_channel_heatmap_per_vmm(
            df_snr_ch_uncut,
            detector_vmms,
            show_quality_overlay=True
        )

        # Version 2 — only quality-passing channels shown
        plot_snr_channel_heatmap_per_vmm(
            df_snr_ch,
            detector_vmms,
            show_quality_overlay=False
        )

    if PLOTS["snr_channel_uniformity"]:
        plot_snr_channel_uniformity(df_snr_ch, detector_vmms)


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

    df_summary = summarise_best_config(df_snr, df_snr_ch, detector_vmms)
    df_summary.to_csv("vmm_snr_summary.csv", index=False)

if __name__ == "__main__":
    main()


print("bonzo")
