#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/25/26 1:24 PM
Created in PyCharm
Created as vmm_config_scan_analysis.py

@author: ak271430 Alexandra Kallitsopoulou (alexandra.kallitsopoulou@cea.fr)
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
from vmm_mapping import vmm_mapping  # default; overridden at runtime by dataset

from vmm_io import (load_run_table, get_run_groups,
                    get_run_dir, load_hits_run,
                    filter_runs)

from vmm_noise  import compute_noise_baseline

from vmm_snr import (compute_snr,
                     compute_snr_per_channel,
                     compute_tail_curves,
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
                    qa_mpv_vs_median_comparison,
                    qa_noise_run_diagnostic)

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
                       plot_snr_channel_uniformity,
                       plot_snr_method_comparison,
                       plot_all_methods_heatmap,
                       plot_tail_distributions,
                       plot_saturation_curves)

from vmm_snr import compute_adc_stats


# ─────────────────────────────────────────────
# DATASET PRESETS
# Add a new entry here when starting a new run campaign.
# Keys are the short names used to select a dataset below.
# ─────────────────────────────────────────────
DATASETS = {
    "15kHz": {
        "cnfg_dir"             : "/drf/projets/clas12/P2/akallits/",
        "data_dir"             : "/drf/projets/clas12/cern_202511_p2_alinx/",
        "run_table"            : "vmm_config_scan_15kHz.csv",
        "plot_subdir"          : "plots_15kHz/",
        "mapping_module"       : "vmm_mapping",
        "legacy_sg"            : 3.0,
        "legacy_sng"           : 1.0,
        "qa_run_signal"        : 149,
        "qa_run_noisy"         : 150,
        "qa_runs_quality_check": [149, 150],
        "qa_noise_diag_noise"  : 157,
        "qa_noise_diag_signal" : 156,
    },
    "5kHz": {
        "cnfg_dir"             : "/drf/projets/clas12/P2/akallits/",
        "data_dir"             : "/drf/projets/clas12/cern_202511_p2_alinx/",
        "run_table"            : "vmm_config_scan_5kHz.csv",
        "plot_subdir"          : "plots_5kHz/",
        "mapping_module"       : "vmm_mapping",
        "legacy_sg"            : 4.5,
        "legacy_sng"           : 1.0,
        "qa_run_signal"        : 79,
        "qa_run_noisy"         : 98,
        "qa_runs_quality_check": [79, 96, 98],
        "qa_noise_diag_noise"  : 79,
        "qa_noise_diag_signal" : 98,
    },
    "lab": {
        "cnfg_dir"             : "/drf/projets/clas12/P2/akallits/",
        "data_dir"             : "/drf/projets/clas12/P2/vmm_config_scan_lab/",
        "run_table"            : "vmm_config_scan_lab.csv",
        "plot_subdir"          : "plots_lab_test/",
        "mapping_module"       : "vmm_mapping_lab",
        "legacy_sg"            : 1.0,
        "legacy_sng"           : 1.0,
        "qa_run_signal"        : 1,
        "qa_run_noisy"         : 2,
        "qa_runs_quality_check": [2, 4],
        "qa_noise_diag_noise"  : 8,
        "qa_noise_diag_signal" : 7,
    },
}


# ─────────────────────────────────────────────
# PLOT / QA CONTROL
# ─────────────────────────────────────────────
PLOTS = {
    # Legacy ADC plots
    "adc_hist_per_run"          : True,
    "adc_hist_separate_vmm"     : True,
    "mean_vs_peaking"           : False,
    "plot_std_vs_peaking"       : False,
    "compare_full_vs_cut"       : False,
    "removed_fraction"          : False,
    "robust_stats"              : False,

    # QA investigations
    "qa_noise_run_diagnostic"   : True,
    "qa_over_threshold_split"   : True,
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
    "snr_method_comparison"     : True,

    # New distribution-based comparison plots
    "snr_all_methods_heatmap"   : True,   # 4-panel: MPV/mean/area-ratio/EER
    "tail_distributions"        : True,   # P(noise>x) vs P(signal≤x) per VMM
    "saturation_curves"         : True,   # P(signal>x) + saturation fraction
}


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    import time
    _t0 = time.time()

    def _hdr(title, width=57):
        """Print a timestamped section header."""
        t = time.time() - _t0
        print(f"\n{'─'*width}")
        print(f"  {title}  [{t:.0f}s]")
        print(f"{'─'*width}")

    def _done(msg=""):
        t = time.time() - _t0
        suffix = f"  {msg}" if msg else ""
        print(f"  ✓ done{suffix}  [{t:.0f}s]")

    # ════════════════════════════════════════════════════════════
    # USER CONFIGURATION
    # Normally only these three lines need changing between runs:
    # ════════════════════════════════════════════════════════════

    # --- Dataset ---
    # Pick one key from the DATASETS dict defined at the top of this file.
    # To add a new campaign, add an entry there — nothing else changes here.
    dataset = "5kHz"   # "15kHz" | "5kHz" | "lab"

    # --- Run mode ---
    # "analysis" : save all plots to plot_dir, no interactive windows
    # "debug"    : show each plot interactively (close window to advance),
    #              nothing saved — use the PLOTS dict to select which to step through
    # "both"     : save AND show interactively
    mode = "analysis"

    # --- Files per run ---
    # n_root_files : streaming analysis (compute_snr etc.) — safe at any N.
    #   Runs with fewer files automatically use all they have.
    n_root_files = 10
    # n_qa_files   : QA / legacy plots (loaded into RAM) — keep small.
    n_qa_files = 2

    # ════════════════════════════════════════════════════════════
    # END USER CONFIGURATION — nothing below this line needs editing
    # ════════════════════════════════════════════════════════════

    # Unpack dataset preset
    if dataset not in DATASETS:
        raise ValueError(
            f"Unknown dataset '{dataset}'. "
            f"Available: {list(DATASETS.keys())}"
        )
    _ds = DATASETS[dataset]
    cnfg_dir               = _ds["cnfg_dir"]
    data_dir               = _ds["data_dir"]
    legacy_sg              = _ds["legacy_sg"]
    legacy_sng             = _ds["legacy_sng"]
    qa_run_signal          = _ds["qa_run_signal"]
    qa_run_noisy           = _ds["qa_run_noisy"]
    qa_runs_quality_check  = _ds["qa_runs_quality_check"]
    qa_noise_diag_noise    = _ds["qa_noise_diag_noise"]
    qa_noise_diag_signal   = _ds["qa_noise_diag_signal"]
    plot_dir               = f"{cnfg_dir}{_ds['plot_subdir']}"

    # Load the right VMM mapping and patch all modules that use it
    import importlib as _il
    import vmm_plots as _vp, vmm_io as _vi, vmm_snr as _vs, vmm_qa as _vq
    vmm_mapping = _il.import_module(_ds["mapping_module"]).vmm_mapping
    _vp.vmm_mapping = vmm_mapping
    _vi.vmm_mapping = vmm_mapping
    _vs.vmm_mapping = vmm_mapping
    _vq.vmm_mapping = vmm_mapping

    # Derived from mode
    show_plots = mode in ("debug", "both")
    save_dir   = plot_dir if mode in ("analysis", "both") else None

    W = 57   # banner width

    print(f"\n{'='*W}")
    print(f"  Dataset      : {dataset}")
    print(f"  Mode         : {mode}")
    print(f"  Mapping      : {_ds['mapping_module']}")
    print(f"  data_dir     : {data_dir}")
    print(f"  plot_dir     : {plot_dir if save_dir else '(not saving)'}")
    print(f"  n_root_files : {n_root_files}   n_qa_files : {n_qa_files}")
    print(f"{'='*W}")

    # ── Run metadata ───────────────────────────────────────────
    _hdr("Loading run metadata")
    df_run_scan = load_run_table(f"{cnfg_dir}{_ds['run_table']}")
    run_groups  = get_run_groups(df_run_scan)

    sng1_runs = run_groups["sng1_runs"]
    pairs     = run_groups["pairs"]

    print(f"\n  {len(pairs)} config pair(s) found:")
    print(f"  {'sg':>6} {'snt':>6} {'sng0':>8} {'sng1':>8}")
    for p in pairs:
        print(f"  {p['sg']:>6} {p['snt']:>6} "
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
    print(f"\n  Detector VMMs : {detector_vmms}")

    # ── Legacy ADC stats ───────────────────────────────────────
    _hdr(f"[1/5] Legacy ADC stats  "
         f"(sg={legacy_sg}, sng={legacy_sng})")
    run_num_sg_st_sng = filter_runs(
        df_run_scan, sg=legacy_sg, sng=legacy_sng
    )
    print(f"  Runs: {run_num_sg_st_sng}")

    df_results = compute_adc_stats(
        df_run_scan = df_run_scan,
        vmm_ids     = vmm_ids,
        run_list    = run_num_sg_st_sng,
        data_dir    = data_dir,
        adc_cut     = 100,
        n_files     = n_root_files
    )
    df_results.dropna(inplace=True)
    df_results.to_csv(f"{cnfg_dir}vmm_adc_analysis_{dataset}.csv", index=False)
    _done(f"vmm_adc_analysis_{dataset}.csv  ({len(df_results)} rows)")

    # ── VMM-level SNR ──────────────────────────────────────────
    _hdr(f"[2/5] VMM-level SNR  "
         f"({len(pairs)} pairs, up to {n_root_files} files/run)")
    df_snr = compute_snr(
        data_dir      = data_dir,
        pairs         = pairs,
        detector_vmms = detector_vmms,
        n_files       = n_root_files
    )
    df_snr.to_csv(f"{cnfg_dir}vmm_snr_results_{dataset}.csv", index=False)
    _done(f"vmm_snr_results_{dataset}.csv  ({len(df_snr)} rows)")

    df_snr_clean = df_snr[df_snr["noise_quality"] == "ok"]
    print(f"\n  SNR summary (noise quality=ok only):")
    print(df_snr_clean[
        ["sg", "snt", "vmm_id", "noise_sigma", "mpv", "snr"]
    ].to_string(index=False))

    # ── Channel-level SNR ──────────────────────────────────────
    _hdr(f"[3/5] Channel-level SNR  "
         f"({len(pairs)} pairs, up to {n_root_files} files/run)")
    df_snr_ch = compute_snr_per_channel(
        data_dir      = data_dir,
        pairs         = pairs,
        detector_vmms = detector_vmms,
        n_files       = n_root_files
    )
    df_snr_ch.to_csv(f"{cnfg_dir}vmm_snr_per_channel_{dataset}.csv", index=False)
    if not df_snr_ch.empty:
        _done(f"{len(df_snr_ch)} channels  |  "
              f"SNR {df_snr_ch['snr_ch'].min():.1f}–"
              f"{df_snr_ch['snr_ch'].max():.1f}  |  "
              f"median {df_snr_ch['snr_ch'].median():.1f}")
    else:
        _done("0 channels (check noise quality flags)")

    print(f"\n  Re-running without quality cuts (for overlay plots)...")
    df_snr_ch_uncut = compute_snr_per_channel(
        data_dir        = data_dir,
        pairs           = pairs,
        detector_vmms   = detector_vmms,
        n_files         = n_root_files,
        min_noise_hits  = 5,
        min_signal_hits = 10,
        min_sigma       = 0.0,
        max_sigma       = 999.0,
        mpv_min         = 0.0,
        mpv_max         = 1023.0
    )
    print(f"  → {len(df_snr_ch_uncut)} channels (uncut)")

    # ── Best config summary ────────────────────────────────────
    _hdr("[4/5] Best configuration summary")
    df_summary = summarise_best_config(
        df_snr, df_snr_ch, detector_vmms
    )
    df_summary.to_csv(f"{cnfg_dir}vmm_snr_summary_{dataset}.csv", index=False)
    _done(f"vmm_snr_summary_{dataset}.csv")

    # ── QA investigations ──────────────────────────────────────
    n_qa_enabled = sum(
        1 for k in PLOTS if k.startswith("qa_") and PLOTS[k]
    )
    _hdr(f"[5/5] QA & plots  ({n_qa_enabled} QA + "
         f"{sum(1 for k,v in PLOTS.items() if not k.startswith('qa_') and v)} plots enabled)")

    if PLOTS["qa_noise_run_diagnostic"]:
        print(f"  → qa_noise_run_diagnostic  "
              f"(noise=run {qa_noise_diag_noise}, "
              f"signal=run {qa_noise_diag_signal})")
        qa_noise_run_diagnostic(
            data_dir,
            run_noise     = qa_noise_diag_noise,
            run_signal    = qa_noise_diag_signal,
            detector_vmms = detector_vmms,
            n_files       = n_qa_files,
            out_dir       = plot_dir,
            show          = show_plots
        )

    if PLOTS["qa_over_threshold_split"]:
        print(f"  → qa_over_threshold_split  (run {qa_run_signal})")
        run_dir = get_run_dir(data_dir, qa_run_signal)
        df_hits = load_hits_run(
            run_dir, n_files=n_qa_files,
            branches=["adc", "vmm", "ch", "over_threshold"]
        )
        qa_over_threshold_split(
            df_hits, run_no=qa_run_signal,
            vmm_ids=detector_vmms,
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["qa_noise_pedestal_stability"]:
        print(f"  → qa_noise_pedestal_stability")
        qa_noise_pedestal_stability(
            df_run_scan, data_dir, sng1_runs, vmm_groups,
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["qa_signal_distributions"]:
        print(f"  → qa_signal_distributions  (run {qa_run_signal})")
        run_dir = get_run_dir(data_dir, qa_run_signal)
        df_hits = load_hits_run(
            run_dir, n_files=n_qa_files,
            branches=["adc", "vmm", "ch", "over_threshold"]
        )
        qa_signal_distributions(
            df_hits, detector_vmms, qa_run_signal,
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["qa_noise_quality_check"]:
        print(f"  → qa_noise_quality_check  "
              f"(runs {qa_runs_quality_check})")
        qa_noise_quality_check(
            df_run_scan, data_dir,
            runs_to_check=qa_runs_quality_check
        )

    if any([PLOTS["qa_channel_noise"],
            PLOTS["qa_adc16_artifact"],
            PLOTS["qa_mpv_estimation"]]):

        run_dir    = get_run_dir(data_dir, qa_run_noisy)
        df_hits_qa = load_hits_run(
            run_dir, n_files=n_qa_files,
            branches=["adc", "vmm", "ch", "over_threshold"]
        )
        if df_hits_qa is None:
            print(f"  WARNING: no data for qa_run_noisy={qa_run_noisy}, "
                  f"skipping channel-noise / ADC16 / MPV-estimation QA")
        else:
            df_noise_qa = compute_noise_baseline(df_hits_qa)

            if PLOTS["qa_channel_noise"]:
                print(f"  → qa_channel_noise  (run {qa_run_noisy})")
                qa_channel_noise(
                    df_hits_qa, detector_vmms,
                    run_no=qa_run_noisy,
                    out_dir=save_dir, show=show_plots
                )
            if PLOTS["qa_adc16_artifact"]:
                print(f"  → qa_adc16_artifact  (run {qa_run_noisy})")
                qa_adc16_artifact(
                    df_hits_qa, detector_vmms,
                    run_no=qa_run_noisy,
                    out_dir=save_dir, show=show_plots
                )
            if PLOTS["qa_mpv_estimation"]:
                print(f"  → qa_mpv_estimation  (run {qa_run_noisy})")
                qa_mpv_estimation(
                    df_hits_qa, df_noise_qa,
                    detector_vmms, run_no=qa_run_noisy,
                    out_dir=save_dir, show=show_plots
                )

    if PLOTS["qa_noise_sigma_distribution"]:
        print(f"  → qa_noise_sigma_distribution")
        qa_noise_sigma_distribution(
            df_run_scan, data_dir, sng1_runs,
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["qa_robust_vs_std_comparison"]:
        print(f"  → qa_robust_vs_std_comparison")
        qa_robust_vs_std_comparison(
            df_run_scan, data_dir, sng1_runs,
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["qa_mpv_vs_median_comparison"]:
        print(f"  → qa_mpv_vs_median_comparison")
        qa_mpv_vs_median_comparison(
            df_run_scan, data_dir, pairs, detector_vmms,
            out_dir=save_dir, show=show_plots
        )

    # ── Legacy plots ───────────────────────────────────────────
    if PLOTS["adc_hist_per_run"]:
        print(f"  → adc_hist_per_run")
        plot_adc_histograms_for_runs(
            run_num_sg_st_sng, data_dir,
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["adc_hist_separate_vmm"]:
        print(f"  → adc_hist_separate_vmm")
        plot_adc_by_vmm(
            vmm_ids, run_num_sg_st_sng, df_run_scan, data_dir,
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["mean_vs_peaking"]:
        print(f"  → mean_vs_peaking")
        plot_mean_vs_peaking(df_results, vmm_ids,
                             out_dir=save_dir, show=show_plots)

    if PLOTS["plot_std_vs_peaking"]:
        print(f"  → std_vs_peaking")
        plot_std_vs_peaking(df_results, vmm_ids,
                            out_dir=save_dir, show=show_plots)

    if PLOTS["compare_full_vs_cut"]:
        print(f"  → compare_full_vs_cut")
        for vmm_id in vmm_ids:
            for run_no in run_num_sg_st_sng:
                run_dir = get_run_dir(data_dir, run_no)
                df_hits = load_hits_run(
                    run_dir, n_files=n_qa_files,
                    branches=["adc", "vmm", "ch",
                              "time", "over_threshold"]
                )
                if df_hits is None:
                    continue
                compare_full_vs_cut(
                    df_hits, vmm_id, run_no, adc_cut=100,
                    out_dir=save_dir, show=show_plots
                )

    if PLOTS["removed_fraction"]:
        print(f"  → removed_fraction")
        plot_removed_fraction(df_results, vmm_ids,
                              out_dir=save_dir, show=show_plots)

    if PLOTS["robust_stats"]:
        print(f"  → robust_stats")
        plot_robust_vs_peaking(df_results, vmm_ids,
                               out_dir=save_dir, show=show_plots)

    # ── SNR comparison plots ───────────────────────────────────
    if PLOTS["snr_vs_peaking"]:
        print(f"  → snr_vs_peaking  (MPV + mean)")
        plot_snr_vs_peaking(df_snr,
                            out_dir=save_dir, show=show_plots)
        plot_snr_vs_peaking(df_snr, snr_col="snr_mean",
                            out_dir=save_dir, show=show_plots)

    if PLOTS["snr_vs_gain"]:
        print(f"  → snr_vs_gain  (MPV + mean)")
        plot_snr_vs_gain(df_snr,
                         out_dir=save_dir, show=show_plots)
        plot_snr_vs_gain(df_snr, snr_col="snr_mean",
                         out_dir=save_dir, show=show_plots)

    if PLOTS["snr_heatmap"]:
        print(f"  → snr_heatmap  (MPV + mean)")
        plot_snr_heatmap(df_snr,
                         out_dir=save_dir, show=show_plots)
        plot_snr_heatmap(df_snr, snr_col="snr_mean",
                         out_dir=save_dir, show=show_plots)

    if PLOTS["adc_heatmap"]:
        print(f"  → adc_heatmap  (mpv, mean_signal, noise_sigma, mean_noise)")
        plot_adc_heatmap(df_snr, metric="mpv",
                         out_dir=save_dir, show=show_plots)
        plot_adc_heatmap(df_snr, metric="mean_signal",
                         out_dir=save_dir, show=show_plots)
        plot_adc_heatmap(df_snr, metric="noise_sigma",
                         out_dir=save_dir, show=show_plots)
        plot_adc_heatmap(df_snr, metric="mean_noise",
                         out_dir=save_dir, show=show_plots)

    if PLOTS["snr_channel_heatmap_per_vmm"]:
        print(f"  → snr_channel_heatmap_per_vmm  (MPV + mean)")
        plot_snr_channel_heatmap_per_vmm(
            df_snr_ch_uncut,
            detector_vmms,
            show_quality_overlay=True,
            out_dir=save_dir, show=show_plots
        )
        plot_snr_channel_heatmap_per_vmm(
            df_snr_ch_uncut,
            detector_vmms,
            snr_col="snr_mean_ch",
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["snr_channel_heatmap_all_configs"]:
        print(f"  → snr_channel_heatmap_all_configs  (MPV + mean)")
        plot_snr_channel_heatmap_all_configs(
            df_snr_ch,
            out_dir=save_dir, show=show_plots
        )
        plot_snr_channel_heatmap_all_configs(
            df_snr_ch, snr_col="snr_mean_ch",
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["snr_channel_uniformity"]:
        print(f"  → snr_channel_uniformity  (MPV + mean)")
        plot_snr_channel_uniformity(
            df_snr_ch, detector_vmms,
            out_dir=save_dir, show=show_plots
        )
        plot_snr_channel_uniformity(
            df_snr_ch, detector_vmms,
            snr_col="snr_mean_ch",
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["snr_method_comparison"]:
        print(f"  → snr_method_comparison")
        plot_snr_method_comparison(
            df_snr, df_snr_ch,
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["snr_all_methods_heatmap"]:
        print(f"  → snr_all_methods_heatmap  (4-panel)")
        plot_all_methods_heatmap(
            df_snr,
            out_dir=save_dir, show=show_plots
        )

    if PLOTS["tail_distributions"] or PLOTS["saturation_curves"]:
        print(f"  → tail curves  (loading histograms...)")
        tail_data = compute_tail_curves(
            data_dir      = data_dir,
            pairs         = pairs,
            detector_vmms = detector_vmms,
            n_files       = n_root_files
        )
        print(f"    {len(tail_data)} (config, VMM) entries loaded")

        if PLOTS["tail_distributions"]:
            print(f"  → tail_distributions")
            plot_tail_distributions(
                tail_data, detector_vmms,
                out_dir=save_dir, show=show_plots
            )

        if PLOTS["saturation_curves"]:
            print(f"  → saturation_curves")
            plot_saturation_curves(
                tail_data, detector_vmms,
                out_dir=save_dir, show=show_plots
            )

    print(f"\n{'='*W}")
    print(f"  Finished  [{time.time() - _t0:.0f}s total]")
    if save_dir:
        print(f"  Plots saved to: {save_dir}")
    print(f"{'='*W}\n")


if __name__ == "__main__":
    main()
