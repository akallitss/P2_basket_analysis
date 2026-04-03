#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/25/26 1:24 PM
Created in PyCharm
Created as vmm_config_scan_analysis.py

@author: ak271430
"""

"""
vmm_config_scan_analysis.py

Main entry point for VMM configuration scan analysis.

Orchestrates the full analysis pipeline:
1. Load run metadata and derive run groupings
2. Compute legacy ADC statistics (mean/std vs peaking time)
3. Compute VMM-level SNR per configuration pair
4. Compute channel-level SNR per configuration pair
5. Summarise best configuration per VMM
6. QA investigations (toggled via PLOTS dict)
7. Comparison plots (toggled via PLOTS dict)

All analysis logic lives in dedicated modules:
    vmm_io.py      — data loading and run management
    vmm_noise.py   — noise baseline estimation
    vmm_signal.py  — signal extraction and MPV estimation
    vmm_snr.py     — SNR computation and summary
    vmm_qa.py      — quality assurance investigations
    vmm_plots.py   — all plotting functions
    vmm_mapping.py — detector geometry and VMM groupings
"""
from vmm_mapping import vmm_mapping

from vmm_io import (load_run_table, get_run_groups,
                    get_run_dir, load_hits_run,
                    filter_runs)

from vmm_noise  import compute_noise_baseline

from vmm_snr import (compute_snr,
                     compute_snr_per_channel,
                     summarise_best_config)

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

from vmm_plots import (plot_adc_histograms_for_runs,
                       plot_adc_by_vmm,
                       compare_full_vs_cut,
                       plot_removed_fraction,
                       plot_mean_vs_peaking,
                       plot_std_vs_peaking,
                       plot_robust_vs_peaking,
                       plot_snr_vs_peaking,
                       plot_snr_vs_gain,
                       plot_snr_heatmap,
                       plot_adc_heatmap,
                       plot_snr_channel_heatmap_per_vmm,
                       plot_snr_channel_heatmap_all_configs,
                       plot_snr_channel_uniformity)

from vmm_snr import compute_adc_stats


# ─────────────────────────────────────────────
# PLOT / QA CONTROL
# ─────────────────────────────────────────────
PLOTS = {
    # Legacy ADC plots
    "adc_hist_per_run"          : False,
    "adc_hist_separate_vmm"     : False,
    "mean_vs_peaking"           : False,
    "plot_std_vs_peaking"       : False,
    "compare_full_vs_cut"       : False,
    "removed_fraction"          : False,
    "robust_stats"              : False,

    # QA investigations
    "qa_over_threshold_split"   : False,
    "qa_noise_pedestal_stability": True,
    "qa_signal_distributions"   : False,
    "qa_noise_quality_check"    : True,
    "qa_mpv_estimation"         : False,
    "qa_channel_noise"          : True,
    "qa_adc16_artifact"         : True,
    "qa_noise_sigma_distribution": False,
    "qa_robust_vs_std_comparison": False,
    "qa_mpv_vs_median_comparison": False,

    # SNR comparison plots
    "snr_vs_peaking"            : True,
    "snr_vs_gain"               : True,
    "snr_heatmap"               : True,
    "adc_heatmap"               : True,
    "snr_channel_heatmap_per_vmm": True,
    "snr_channel_heatmap_all_configs": True,
    "snr_channel_uniformity"    : True,
}


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():

    # ════════════════════════════════════════════════════════════
    # USER CONFIGURATION
    # Edit this section for each new dataset or config scan.
    # ════════════════════════════════════════════════════════════

    # --- Paths ---
    cnfg_dir = "/drf/projets/clas12/P2/akallits/"
    # cnfg_dir = "/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/VMM-alinx-data/"
    data_dir = "/drf/projets/clas12/cern_202511_p2_alinx/"
    # data_dir = "/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/VMM-alinx-data/5kHz-muons-config-scan/"

    # --- Legacy ADC stats filter ---
    # Which gain and neighbor setting to use for legacy plots
    # legacy_sg  = 4.5 sps beam test
    legacy_sg  = 1.0
    legacy_sng = 1.0

    # --- Files per run ---
    # How many ROOT files to load and merge per run directory.
    # n_root_files=1 uses only the first file.
    # Increase to merge more statistics (e.g. 2, 3, or 99 for all).
    n_root_files = 1

    # --- QA investigation runs ---
    # Run to use for over_threshold split and signal distribution QA.
    # Should be a well-behaved sng=1 run at your preferred config.
    # qa_run_signal = 79 beam test sps sng = 1
    qa_run_signal = 1

    # Run to use for channel noise and ADC=16 artifact QA.
    # Should be the most problematic run (short peaking time).
    # qa_run_noisy  = 98 beam test sps
    qa_run_noisy  = 2

    # Runs to check in noise quality check QA.
    # Include one run per peaking time to see the full picture.
    # qa_runs_quality_check = [79, 96, 98] beam test sps sng=1
    qa_runs_quality_check = [2,4]

    # ════════════════════════════════════════════════════════════
    # END USER CONFIGURATION
    # ════════════════════════════════════════════════════════════

    # ── Run metadata ───────────────────────────────────────────
    df_run_scan = load_run_table(f"{cnfg_dir}vmm_config_scan_lab.csv")
    run_groups  = get_run_groups(df_run_scan)

    sng1_runs = run_groups["sng1_runs"]
    pairs     = run_groups["pairs"]

    print("\nRun pairs found:")
    print(f"{'sg':>6} {'snt':>6} {'sng0':>8} {'sng1':>8}")
    for p in pairs:
        print(f"{p['sg']:>6} {p['snt']:>6} "
              f"{p['sng0']:>8} {p['sng1']:>8}")

    # ── VMM groupings from vmm_mapping ─────────────────────────
    detector_vmms = [
        vid
        for key, cfg in vmm_mapping.items()
        if key != "trigger"
        for vid in cfg["vmm_ids"]
    ]
    vmm_groups = [
        {"label": cfg.get("name", key), "vmm_ids": cfg["vmm_ids"]}
        for key, cfg in vmm_mapping.items()
        if key != "trigger"
    ]
    vmm_ids = sorted(
        {vid for cfg in vmm_mapping.values()
         for vid in cfg["vmm_ids"]}
    )

    # ── Legacy ADC stats ───────────────────────────────────────
    run_num_sg_st_sng = filter_runs(
        df_run_scan, sg=legacy_sg, sng=legacy_sng
    )

    df_results = compute_adc_stats(
        df_run_scan = df_run_scan,
        vmm_ids     = vmm_ids,
        run_list    = run_num_sg_st_sng,
        data_dir    = data_dir,
        adc_cut     = 100,
        n_files     = n_root_files
    )
    df_results.dropna(inplace=True)
    df_results.to_csv(f"{cnfg_dir}vmm_adc_analysis.csv", index=False)

    # ── VMM-level SNR ──────────────────────────────────────────
    df_snr = compute_snr(
        data_dir      = data_dir,
        pairs         = pairs,
        detector_vmms = detector_vmms,
        n_files       = n_root_files
    )
    df_snr.to_csv(f"{cnfg_dir}vmm_snr_results.csv", index=False)

    print("\n=== SNR Summary (noise quality=ok only) ===")
    df_snr_clean = df_snr[df_snr["noise_quality"] == "ok"]
    print(df_snr_clean[
        ["sg", "snt", "vmm_id",
         "noise_sigma", "mpv", "snr"]
    ].to_string(index=False))

    # ── Channel-level SNR ──────────────────────────────────────
    df_snr_ch = compute_snr_per_channel(
        data_dir      = data_dir,
        pairs         = pairs,
        detector_vmms = detector_vmms,
        n_files       = n_root_files
    )
    df_snr_ch.to_csv(f"{cnfg_dir}vmm_snr_per_channel.csv", index=False)

    print(f"\nChannel-level SNR: {len(df_snr_ch)} entries")
    print(f"SNR range : {df_snr_ch['snr_ch'].min():.1f} "
          f"— {df_snr_ch['snr_ch'].max():.1f}")
    print(f"Median SNR: {df_snr_ch['snr_ch'].median():.1f}")

    # Uncut version for overlay plots
    df_snr_ch_uncut = compute_snr_per_channel(
        data_dir       = data_dir,
        pairs          = pairs,
        detector_vmms  = detector_vmms,
        n_files        = n_root_files,
        min_noise_hits = 5,
        min_signal_hits   = 10,
        min_sigma         = 0.0,
        max_sigma         = 999.0,
        mpv_min           = 0.0,
        mpv_max           = 1023.0
    )

    # ── Best config summary ────────────────────────────────────
    df_summary = summarise_best_config(
        df_snr, df_snr_ch, detector_vmms
    )
    df_summary.to_csv(f"{cnfg_dir}vmm_snr_summary.csv", index=False)

    # ── QA investigations ──────────────────────────────────────
    if PLOTS["qa_over_threshold_split"]:
        run_dir = get_run_dir(data_dir, qa_run_signal)
        df_hits = load_hits_run(
            run_dir, n_files=n_root_files,
            branches=["adc", "vmm", "ch", "over_threshold"]
        )
        qa_over_threshold_split(
            df_hits, run_no=qa_run_signal,
            vmm_ids=detector_vmms
        )

    if PLOTS["qa_noise_pedestal_stability"]:
        qa_noise_pedestal_stability(
            df_run_scan, data_dir, sng1_runs, vmm_groups
        )

    if PLOTS["qa_signal_distributions"]:
        run_dir = get_run_dir(data_dir, qa_run_signal)
        df_hits = load_hits_run(
            run_dir, n_files=n_root_files,
            branches=["adc", "vmm", "ch", "over_threshold"]
        )
        qa_signal_distributions(
            df_hits, detector_vmms, qa_run_signal
        )

    if PLOTS["qa_noise_quality_check"]:
        qa_noise_quality_check(
            df_run_scan, data_dir,
            runs_to_check=qa_runs_quality_check
        )

    if any([PLOTS["qa_channel_noise"],
            PLOTS["qa_adc16_artifact"],
            PLOTS["qa_mpv_estimation"]]):

        run_dir    = get_run_dir(data_dir, qa_run_noisy)
        df_hits_qa = load_hits_run(
            run_dir, n_files=n_root_files,
            branches=["adc", "vmm", "ch", "over_threshold"]
        )
        df_noise_qa = compute_noise_baseline(df_hits_qa)

        if PLOTS["qa_channel_noise"]:
            qa_channel_noise(
                df_hits_qa, detector_vmms,
                run_no=qa_run_noisy
            )
        if PLOTS["qa_adc16_artifact"]:
            qa_adc16_artifact(
                df_hits_qa, detector_vmms,
                run_no=qa_run_noisy
            )
        if PLOTS["qa_mpv_estimation"]:
            qa_mpv_estimation(
                df_hits_qa, df_noise_qa,
                detector_vmms, run_no=qa_run_noisy
            )

    if PLOTS["qa_noise_sigma_distribution"]:
        qa_noise_sigma_distribution(
            df_run_scan, data_dir, sng1_runs
        )

    if PLOTS["qa_robust_vs_std_comparison"]:
        qa_robust_vs_std_comparison(
            df_run_scan, data_dir, sng1_runs
        )

    if PLOTS["qa_mpv_vs_median_comparison"]:
        qa_mpv_vs_median_comparison(
            df_run_scan, data_dir, pairs, detector_vmms
        )

    # ── Legacy plots ───────────────────────────────────────────
    if PLOTS["adc_hist_per_run"]:
        plot_adc_histograms_for_runs(
            run_num_sg_st_sng, data_dir
        )

    if PLOTS["adc_hist_separate_vmm"]:
        plot_adc_by_vmm(
            vmm_ids, run_num_sg_st_sng, df_run_scan, data_dir
        )

    if PLOTS["mean_vs_peaking"]:
        plot_mean_vs_peaking(df_results, vmm_ids)

    if PLOTS["plot_std_vs_peaking"]:
        plot_std_vs_peaking(df_results, vmm_ids)

    if PLOTS["compare_full_vs_cut"]:
        for vmm_id in vmm_ids:
            for run_no in run_num_sg_st_sng:
                run_dir = get_run_dir(data_dir, run_no)
                df_hits = load_hits_run(
                    run_dir, n_files=n_root_files,
                    branches=["adc", "vmm", "ch",
                              "time", "over_threshold"]
                )
                if df_hits is None:
                    continue
                compare_full_vs_cut(
                    df_hits, vmm_id, run_no, adc_cut=100
                )

    if PLOTS["removed_fraction"]:
        plot_removed_fraction(df_results, vmm_ids)

    if PLOTS["robust_stats"]:
        plot_robust_vs_peaking(df_results, vmm_ids)

    # ── SNR comparison plots ───────────────────────────────────
    if PLOTS["snr_vs_peaking"]:
        plot_snr_vs_peaking(df_snr)

    if PLOTS["snr_vs_gain"]:
        plot_snr_vs_gain(df_snr)

    if PLOTS["snr_heatmap"]:
        plot_snr_heatmap(df_snr)

    if PLOTS["adc_heatmap"]:
        plot_adc_heatmap(df_snr, metric="mpv")
        plot_adc_heatmap(df_snr, metric="noise_sigma")

    if PLOTS["snr_channel_heatmap_per_vmm"]:
        plot_snr_channel_heatmap_per_vmm(
            df_snr_ch_uncut,
            detector_vmms,
            show_quality_overlay=True
        )

    if PLOTS["snr_channel_heatmap_all_configs"]:
        plot_snr_channel_heatmap_all_configs(df_snr_ch)

    if PLOTS["snr_channel_uniformity"]:
        plot_snr_channel_uniformity(df_snr_ch, detector_vmms)

    print('bonzo')

if __name__ == "__main__":
    main()
