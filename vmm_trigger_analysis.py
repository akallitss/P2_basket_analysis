#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/30/26 1:18 PM
Created in PyCharm
Created as vmm_trigger_analysis.py
Trigger rate and efficiency analysis for VMM configuration scan.

Compares trigger rate (VMM 0/1) vs detector VMM hit rates
across different configurations (sg, snt).
@author: ak271430
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from vmm_mapping import vmm_mapping
from vmm_io import (load_run_table, get_run_groups,
                    get_run_dir, get_root_file,
                    load_hits_root,
                    iter_hits_files, list_root_files)
# ── Constants ──────────────────────────────────────────────
NS_PER_TICK  = 1.0          # 1 GHz clock
S_PER_TICK   = NS_PER_TICK * 1e-9


def _finish_fig(fig, stem, out_dir, show, rate_tag=""):
    """Save fig as PDF and PNG to out_dir, optionally display, then close."""
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        full = f"{stem}_{rate_tag}" if rate_tag else stem
        for ext in ("pdf", "png"):
            fig.savefig(os.path.join(out_dir, f"{full}.{ext}"),
                        bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)

# ── VMM groups from mapping ─────────────────────────────────
trigger_vmms  = vmm_mapping["trigger"]["vmm_ids"]
detector_vmms = [
    vid
    for key, cfg in vmm_mapping.items()
    if key != "trigger"
    for vid in cfg["vmm_ids"]
]
vmm_groups = [
    {"label": cfg.get("name", key),
     "vmm_ids": cfg["vmm_ids"]}
    for key, cfg in vmm_mapping.items()
    if key != "trigger"
]
trigger_ref_channels = {
    vmm_id: vmm_mapping["trigger"]["channels"][vmm_id][0]
    for vmm_id in trigger_vmms
}

print("Trigger reference channels:")
for vmm_id, ch_id in trigger_ref_channels.items():
    print(f"  VMM {vmm_id} → ch {ch_id}")

print("\nDetector VMM groups:")
for g in vmm_groups:
    print(f"  {g['label']}: VMMs {g['vmm_ids']}")


# ─────────────────────────────────────────────
# DATASET PRESETS
# Add a new entry here when starting a new run campaign.
# ─────────────────────────────────────────────
DATASETS = {
    "15kHz": {
        "cnfg_dir"        : "/drf/projets/clas12/P2/akallits/",
        "data_dir"        : "/drf/projets/clas12/cern_202511_p2_alinx/",
        "config_file"     : "vmm_config_scan_15kHz.csv",
        "plot_subdir"     : "plots_trigger_15kHz/",
        "test_run"        : 149,
        "root_file_index" : 2,
        "diagnostic_runs" : [
            (149, 3.0, 200),
            (156, 3.0, 100),
            (158, 3.0, 50),
            (152, 4.5, 200),
            (155, 6.0, 200),
        ],
    },
    "5kHz": {
        "cnfg_dir"        : "/drf/projets/clas12/P2/akallits/",
        "data_dir"        : "/drf/projets/clas12/cern_202511_p2_alinx/",
        "config_file"     : "vmm_config_scan_5kHz.csv",
        "plot_subdir"     : "plots_trigger_5kHz/",
        "test_run"        : 82,
        "root_file_index" : 2,
        "diagnostic_runs" : [
            (67, 3.0, 200),
            (70, 3.0, 100),
            (73, 3.0, 50),
            (97, 4.5, 50),
            (78, 4.5, 200),
            (82, 4.5, 100),
            (99, 6.0, 200),
        ],
    },
}


def main():
    # ════════════════════════════════════════════════════════════
    # USER CONFIGURATION — normally only these three lines change:
    # ════════════════════════════════════════════════════════════

    dataset = "15kHz"    # "15kHz" | "5kHz"
    mode    = "analysis" # "analysis": save, no window
                         # "debug"   : show window, nothing saved
                         # "both"    : save and show window
    n_files = 10         # files per run; fewer available → use all

    # ════════════════════════════════════════════════════════════
    # END USER CONFIGURATION
    # ════════════════════════════════════════════════════════════

    if dataset not in DATASETS:
        raise ValueError(
            f"Unknown dataset '{dataset}'. "
            f"Available: {list(DATASETS.keys())}"
        )
    _ds             = DATASETS[dataset]
    cnfg_dir        = _ds["cnfg_dir"]
    data_dir        = _ds["data_dir"]
    config_file     = _ds["config_file"]
    plot_dir        = f"{cnfg_dir}{_ds['plot_subdir']}"
    test_run        = _ds["test_run"]
    root_file_index = _ds["root_file_index"]
    diagnostic_runs = _ds["diagnostic_runs"]
    rate_tag        = config_file.replace("vmm_config_scan_", "").replace(".csv", "")

    show    = mode in ("debug", "both")
    out_dir = plot_dir if mode in ("analysis", "both") else None
    verbose = (mode == "debug")

    import time
    _t0 = time.time()

    def _hdr(title, W=57):
        t = time.time() - _t0
        print(f"\n{'─'*W}")
        print(f"  {title}  [{t:.0f}s]")
        print(f"{'─'*W}")

    def _done(msg=""):
        t = time.time() - _t0
        suffix = f"  {msg}" if msg else ""
        print(f"  ✓ done{suffix}  [{t:.0f}s]")

    W = 57
    print(f"\n{'='*W}")
    print(f"  Dataset  : {dataset}")
    print(f"  Mode     : {mode}")
    print(f"  data_dir : {data_dir}")
    print(f"  plot_dir : {plot_dir if out_dir else '(not saving)'}")
    print(f"  n_files  : {n_files}   test_run : {test_run}")
    print(f"{'='*W}")

    # ── Step 0: channel diagnostics (hits per channel) ──────
    all_vmm_ids  = trigger_vmms + detector_vmms
    diag_run_nos = [r for r, _, _ in diagnostic_runs]
    _hdr(f"Step 0: channel diagnostics  ({len(diag_run_nos)} runs)")
    ch_counts = collect_channel_hit_counts(
        data_dir, diag_run_nos, all_vmm_ids,
        n_files=n_files, file_start=1
    )
    ch_outliers = identify_outlier_channels(
        ch_counts,
        trigger_vmm_ids=trigger_vmms,
        dead_frac=0.1,
        noisy_factor=5.0
    )
    plot_hits_per_channel_all_vmms(
        ch_counts, ch_outliers,
        dead_frac=0.1, noisy_factor=5.0,
        out_dir=out_dir, show=show, rate_tag=rate_tag,
        verbose=verbose
    )
    _done()

    # ── Step 0b: per-channel on/off rate for diagnostic runs ──
    n_diag = len(diagnostic_runs)
    _hdr(f"Step 0b: per-channel spill rates  ({n_diag} runs)")
    per_ch_all_runs  = []
    per_ch_by_run    = {}
    for _i, (run_no, sg, snt) in enumerate(diagnostic_runs):
        if verbose:
            print(f"\n── Per-channel on/off rate: run {run_no} "
                  f"(sg={sg} snt={snt}) ──")
        else:
            filled = int(20 * (_i + 1) / n_diag)
            bar    = "=" * filled + "-" * (20 - filled)
            print(f"\r  [{bar}] {_i+1}/{n_diag}  run {run_no}",
                  end="", flush=True)
        per_ch = compute_per_channel_spill_rates(
            data_dir, run_no, all_vmm_ids,
            trigger_ref_channels=trigger_ref_channels,
            n_files=n_files, file_start=1,
            spill_threshold_khz=1.0,
            max_gap_s=2.0,
            verbose=verbose,
        )
        if per_ch is not None:
            per_ch_all_runs.append(per_ch)
            per_ch_by_run[run_no] = (sg, snt, per_ch)
    if not verbose:
        print()  # end progress-bar line
    _done()

    # ── Step 0c: build good-channel mask ───────────────────────
    _hdr("Step 0c: channel masking")
    if per_ch_all_runs:
        noisy_off   = identify_noisy_from_off_spill(
            per_ch_all_runs, noisy_factor=3.0, ch_outliers=ch_outliers
        )
        good_channels = build_good_channels(ch_outliers, noisy_off)
    else:
        noisy_off     = None
        good_channels = build_good_channels(ch_outliers)

    if verbose:
        print("\n── Good channels per detector VMM (global mask) ──")
        for vmm_id in detector_vmms:
            if vmm_id not in good_channels:
                continue
            chs = good_channels[vmm_id]
            print(f"  VMM {vmm_id}: {len(chs)} good channels  "
                  f"{sorted(chs.tolist())}")
    else:
        n_good_total = sum(len(v) for v in good_channels.values())
        print(f"  {n_good_total} good channels across "
              f"{len(good_channels)} detector VMMs")

    good_channels_by_config, noisy_off_by_config = build_good_channels_per_config(
        per_ch_by_run, ch_outliers, noisy_factor=3.0
    )
    if verbose:
        print("\n── Good channels per config (per-config mask) ──")
        for (sg, snt), gc in sorted(good_channels_by_config.items()):
            n_cfg_total = sum(len(v) for v in gc.values())
            print(f"  sg={sg} snt={snt:.0f}: {n_cfg_total} total good channels")

    for run_no, (sg, snt, per_ch) in per_ch_by_run.items():
        plot_rate_per_channel_on_off(
            per_ch, run_no, sg, snt,
            good_channels=good_channels_by_config.get((sg, snt), good_channels),
            ch_outliers=ch_outliers,
            noisy_off_spill=noisy_off_by_config.get((sg, snt), noisy_off),
            out_dir=out_dir, show=show,
            rate_tag=rate_tag
        )
    _done()

    # ── Steps 1-2: spill structure on test run ───────────────
    _hdr(f"Steps 1-2: spill structure  (run {test_run})")
    if verbose:
        print(f"  Loading run {test_run}...")
    df_hits = load_sorted_hits(
        data_dir, test_run,
        branches=["time", "vmm", "ch"],
        root_file_index=root_file_index
    )
    if df_hits is None or df_hits.empty:
        print(f"  No data for run {test_run} at "
              f"file_index={root_file_index}, "
              f"try a different root_file_index")
        return

    duration = compute_run_duration(df_hits)
    if verbose:
        print(f"  Run duration : {duration:.1f} s")
        print(f"  Total hits   : {len(df_hits)}")

    plot_trigger_rate_ms(
        df_hits,
        run_no=test_run,
        trigger_vmm=0,
        trigger_ch=trigger_ref_channels[0],
        bin_width_ms=1.0,
        out_dir=out_dir, show=show, rate_tag=rate_tag,
    )
    plot_spill_mask_diagnostic(
        df_hits,
        run_no=test_run,
        trigger_vmm=0,
        trigger_ch=trigger_ref_channels[0],
        spill_threshold_khz=1.0,
        bin_width_ms=1.0,
        max_gap_s=2.0,
        out_dir=out_dir, show=show, rate_tag=rate_tag,
        verbose=verbose,
    )
    _done()

    # ── Step 3: sanity check on test run ─────────────────────
    _hdr(f"Step 3: sanity check  (run {test_run})")
    df_run_scan = load_run_table(f"{cnfg_dir}{config_file}")

    df_spill_test = compute_spill_rates_all_runs(
        df_run_scan          = df_run_scan,
        data_dir             = data_dir,
        sng0_runs            = [test_run],
        trigger_ref_channels = trigger_ref_channels,
        detector_vmms        = detector_vmms,
        bin_width_s          = 0.001,
        n_files              = n_files,
        spill_threshold_khz  = 1.0,
        max_gap_s            = 2.0,
        good_channels        = good_channels,
        file_start           = 1,
        verbose              = verbose,
    )
    if verbose:
        print(f"\n  Spill rates on run {test_run}:")
        print(df_spill_test.to_string(index=False))
    _done()

    # ── Step 4: loop over all sng=0 configs ─────────────────
    run_groups = get_run_groups(df_run_scan)
    n_sng0     = len(run_groups["sng0_runs"])
    _hdr(f"Step 4: all configs  ({n_sng0} runs)")

    df_spill = compute_spill_rates_all_runs(
        df_run_scan              = df_run_scan,
        data_dir                 = data_dir,
        sng0_runs                = run_groups["sng0_runs"],
        trigger_ref_channels     = trigger_ref_channels,
        detector_vmms            = detector_vmms,
        bin_width_s              = 0.001,
        n_files                  = n_files,
        spill_threshold_khz      = 1.0,
        max_gap_s                = 2.0,
        good_channels            = good_channels,
        good_channels_per_config = good_channels_by_config,
        verbose                  = verbose,
    )

    csv_path = f"vmm_spill_rates_{rate_tag}.csv"
    df_spill.to_csv(csv_path, index=False)
    print(f"  Saved → {csv_path}")
    if verbose:
        print("\n  Spill rates — all configs:")
        print(df_spill.to_string(index=False))
    _done(f"{len(df_spill)} VMM-run records")

    # ── Step 5: summary plots ────────────────────────────────
    _hdr("Step 5: summary plots")
    trig_ref = list(trigger_ref_channels.keys())[0]

    plot_spill_rates_vs_config(
        df_spill        = df_spill,
        vmm_groups      = vmm_groups,
        trigger_ref_vmm = trig_ref,
        out_dir=out_dir, show=show, rate_tag=rate_tag,
    )
    plot_spill_on_all_vmms(
        df_spill        = df_spill,
        vmm_groups      = vmm_groups,
        trigger_ref_vmm = trig_ref,
        rate_col        = "rate_on_khz",
        out_dir=out_dir, show=show, rate_tag=rate_tag,
    )
    plot_spill_on_all_vmms(
        df_spill        = df_spill,
        vmm_groups      = vmm_groups,
        trigger_ref_vmm = trig_ref,
        rate_col        = "rate_off_khz",
        out_dir=out_dir, show=show, rate_tag=rate_tag,
    )
    _done()

def load_sorted_hits(data_dir, run_no, branches,
                      root_file_index=1,
                      connected_channels_only=False,
                      max_time_ticks=2e12):
    from vmm_io import get_connected_channels
    from functools import reduce
    import operator

    run_dir   = get_run_dir(data_dir, run_no)
    file_path = get_root_file(run_dir,
                               n_files=root_file_index + 1,
                               file_index=root_file_index)
    if file_path is None:
        return None

    df = load_hits_root(file_path, branches=branches)

    # Remove corrupted timestamps first — before anything else
    if "time" in df.columns:
        n_before  = len(df)
        df        = df[df["time"] < max_time_ticks].copy()
        n_corrupt = n_before - len(df)
        if n_corrupt > 0:
            print(f"  Removed {n_corrupt} hits with "
                  f"corrupted timestamps")

    # Filter connected channels
    if connected_channels_only and "ch" in df.columns:
        masks = []
        for vmm_id in df["vmm"].unique():
            channels = get_connected_channels(vmm_id)
            if channels is None:
                masks.append(df["vmm"] == vmm_id)
            else:
                masks.append(
                    (df["vmm"] == vmm_id) &
                    (df["ch"].isin(channels))
                )
        if masks:
            combined = reduce(operator.or_, masks)
            df       = df[combined]

    if df.empty:
        print(f"  Run {run_no} file {root_file_index}: "
              f"no hits after filtering — try a different file_index")
        return None

    df = df.sort_values("time").reset_index(drop=True)
    return df

def inspect_trigger_channels(data_dir, run_no,
                              root_file_index=1,
                              out_dir=None, show=True, rate_tag=""):
    """
    Plot hit count per channel for trigger VMMs.
    Connected channels show significantly more hits than
    floating/unconnected ones — used to validate channel mapping.
    """
    run_dir   = get_run_dir(data_dir, run_no)
    file_path = get_root_file(run_dir,
                               file_index=root_file_index)
    if file_path is None:
        print(f"No file found for run {run_no}")
        return

    df = load_hits_root(
        file_path,
        branches=["time", "vmm", "ch", "adc",
                  "over_threshold"]
    )

    for vmm_id in sorted(trigger_vmms):
        df_vmm = df[df["vmm"] == vmm_id]
        if df_vmm.empty:
            print(f"VMM {vmm_id}: no hits found")
            continue

        # Hit count per channel
        ch_counts = df_vmm.groupby("ch").size().reindex(
            range(64), fill_value=0
        )

        # Mean ADC per channel
        ch_adc = df_vmm.groupby("ch")["adc"].mean().reindex(
            range(64), fill_value=np.nan
        )

        fig, axes = plt.subplots(1, 2, figsize=(14, 4))

        # Hit count per channel
        axes[0].bar(ch_counts.index, ch_counts.values,
                    color="steelblue", alpha=0.7)
        axes[0].set_xlabel("Channel")
        axes[0].set_ylabel("Hit count")
        axes[0].set_title(f"VMM {vmm_id} — hits per channel")
        axes[0].set_yscale("log")

        # Annotate channels with significant hits
        threshold = ch_counts.max() * 0.01  # 1% of max
        active    = ch_counts[ch_counts > threshold]
        for ch, count in active.items():
            axes[0].annotate(
                f"ch{ch}\n{count}",
                xy=(ch, count),
                xytext=(ch, count * 2),
                fontsize=7,
                ha="center",
                color="darkred"
            )

        # Mean ADC per channel — distinguishes real signal
        # from noise (connected channels have higher ADC)
        axes[1].bar(
            ch_adc.dropna().index,
            ch_adc.dropna().values,
            color="orange", alpha=0.7
        )
        axes[1].set_xlabel("Channel")
        axes[1].set_ylabel("Mean ADC")
        axes[1].set_title(f"VMM {vmm_id} — mean ADC per channel")

        # Print summary
        print(f"\nVMM {vmm_id} — channels above 1% of max hits:")
        print(f"  Max hits on any channel: {ch_counts.max()}")
        print(f"  Threshold (1% of max)  : {threshold:.0f}")
        print(f"  Active channels        : "
              f"{sorted(active.index.tolist())}")
        print(f"  Hit counts             : "
              f"{active.to_dict()}")

        plt.suptitle(
            f"VMM {vmm_id} Run {run_no} — "
            f"trigger channel inspection"
        )
        plt.tight_layout()
        _finish_fig(fig, f"trigger_channels_vmm{vmm_id}_run{run_no}",
                    out_dir, show, rate_tag)

def inspect_inter_event_per_channel(df_hits, vmm_id, run_no,
                                     max_dt_us=1000,
                                     out_dir=None, show=True, rate_tag=""):
    """
    Plot inter-event time distribution separately per channel
    for a trigger VMM. Reveals which channel contributes the
    near-zero Δt spike.
    """
    df_vmm = df_hits[df_hits["vmm"] == vmm_id]
    channels = sorted(df_vmm["ch"].unique())

    n_ch  = len(channels)
    fig, axes = plt.subplots(n_ch, 2,
                              figsize=(12, 4 * n_ch))
    if n_ch == 1:
        axes = [axes]

    for idx, ch_id in enumerate(channels):
        t = df_vmm[df_vmm["ch"] == ch_id]["time"].values

        if len(t) < 2:
            continue

        dt_ns = np.diff(t) * NS_PER_TICK
        dt_us = dt_ns[(dt_ns > 0)] / 1000.0
        dt_us = dt_us[dt_us < max_dt_us]

        n_zero = (dt_us < 1).sum()
        n_total = len(dt_us)

        # Linear
        axes[idx][0].hist(dt_us, bins=200,
                          color="steelblue", alpha=0.7)
        axes[idx][0].set_xlabel("Δt (μs)")
        axes[idx][0].set_ylabel("Counts")
        axes[idx][0].set_title(
            f"VMM {vmm_id} ch {ch_id} — linear "
            f"(Δt<1μs: {n_zero}/{n_total} = "
            f"{100*n_zero/n_total:.1f}%)"
        )

        # Log
        axes[idx][1].hist(dt_us, bins=200,
                          color="steelblue", alpha=0.7)
        axes[idx][1].set_yscale("log")
        axes[idx][1].set_xlabel("Δt (μs)")
        axes[idx][1].set_ylabel("Counts (log)")
        axes[idx][1].set_title(
            f"VMM {vmm_id} ch {ch_id} — log scale"
        )

        print(f"VMM {vmm_id} ch {ch_id}: "
              f"n_hits={len(t)}  "
              f"Δt<1μs: {n_zero} ({100*n_zero/n_total:.1f}%)")

    plt.suptitle(
        f"VMM {vmm_id} Run {run_no} — "
        f"inter-event Δt per channel"
    )
    plt.tight_layout()
    _finish_fig(fig,
                f"inter_event_per_ch_vmm{vmm_id}_run{run_no}",
                out_dir, show, rate_tag)

def compute_per_channel_spill_rates(data_dir, run_no, all_vmm_ids,
                                     trigger_ref_channels,
                                     bin_width_s=0.001,
                                     n_files=5, file_start=1,
                                     spill_threshold_khz=1.0,
                                     min_spill_s=1.0,
                                     max_gap_s=0.0,
                                     max_time_ticks=2e12,
                                     verbose=True):
    """
    For every connected channel on each VMM compute the hit rate
    separately during spill-on and spill-off periods.

    Three-pass strategy (memory-safe):
      Pass 1 — time range only (time column).
      Pass 2 — trigger-channel histogram → spill mask.
      Pass 3 — per-channel hit classification using the mask.

    Returns
    -------
    dict {vmm_id: {
        'rate_on_hz' : np.ndarray shape (64,),
        'rate_off_hz': np.ndarray shape (64,),
        't_on_s'     : float,
        't_off_s'    : float,
    }}
    """
    from vmm_io import get_connected_channels

    trig_vmm = list(trigger_ref_channels.keys())[0]
    trig_ch  = trigger_ref_channels[trig_vmm]
    run_dir  = get_run_dir(data_dir, run_no)

    # ── Pass 1: time range ─────────────────────────────────
    t_start, t_end = _get_run_time_range(run_dir, n_files,
                                          file_start=file_start)
    if t_start is None:
        print(f"  Run {run_no}: no valid timestamps")
        return None

    bins   = np.arange(t_start, t_end + bin_width_s, bin_width_s)
    n_bins = len(bins) - 1

    # ── Pass 2: trigger histogram → spill mask ─────────────
    trig_counts = np.zeros(n_bins, dtype=np.int64)
    for df in iter_hits_files(run_dir, n_files=n_files,
                               branches=["time", "vmm", "ch"],
                               file_start=file_start):
        if df is None or df.empty:
            continue
        df = df[df["time"] < max_time_ticks].copy()
        t_s    = df["time"].values * S_PER_TICK
        m_trig = (df["vmm"] == trig_vmm) & (df["ch"] == trig_ch)
        c, _   = np.histogram(t_s[m_trig.values], bins=bins)
        trig_counts += c
        del df

    trig_rate_khz  = trig_counts / bin_width_s / 1e3
    min_spill_bins = max(1, int(min_spill_s / bin_width_s))
    max_gap_bins   = max(0, int(max_gap_s   / bin_width_s))
    on_mask, off_mask = compute_spill_masks(
        trig_rate_khz, spill_threshold_khz, min_spill_bins, max_gap_bins
    )
    t_on_s  = on_mask.sum()  * bin_width_s
    t_off_s = off_mask.sum() * bin_width_s

    if verbose:
        print(f"  Run {run_no}: spill-on {t_on_s:.0f}s  "
              f"spill-off {t_off_s:.0f}s  "
              f"trig_on {trig_rate_khz[on_mask].mean():.1f} kHz")

    # ── Pass 3: per-channel on/off counts ──────────────────
    ch_on  = {v: np.zeros(64, dtype=np.int64) for v in all_vmm_ids}
    ch_off = {v: np.zeros(64, dtype=np.int64) for v in all_vmm_ids}

    for df in iter_hits_files(run_dir, n_files=n_files,
                               branches=["time", "vmm", "ch"],
                               file_start=file_start):
        if df is None or df.empty:
            continue
        df  = df[df["time"] < max_time_ticks].copy()
        df.reset_index(drop=True, inplace=True)   # 0-based positions
        t_s = df["time"].values * S_PER_TICK

        # Bin index for each hit → spill classification
        bidx   = np.clip(np.searchsorted(bins[1:], t_s, side="left"),
                         0, n_bins - 1)
        is_on  = on_mask[bidx]

        for vmm_id in all_vmm_ids:
            connected = get_connected_channels(vmm_id)
            df_v  = df[df["vmm"] == vmm_id]
            if df_v.empty:
                continue

            ch_arr   = df_v["ch"].values
            on_arr   = is_on[df_v.index]   # index is 0-based after reset

            if connected is not None:
                valid = np.isin(ch_arr, connected)
            else:
                valid = (ch_arr >= 0) & (ch_arr < 64)

            ch_v  = ch_arr[valid]
            on_v  = on_arr[valid]
            off_v = ~on_v

            np.add.at(ch_on[vmm_id],  ch_v[on_v  & (ch_v < 64)], 1)
            np.add.at(ch_off[vmm_id], ch_v[off_v & (ch_v < 64)], 1)

        del df

    # ── Convert counts to rates ────────────────────────────
    result = {}
    for vmm_id in all_vmm_ids:
        result[vmm_id] = {
            "rate_on_hz" : ch_on[vmm_id]  / t_on_s  if t_on_s  > 0
                           else np.zeros(64),
            "rate_off_hz": ch_off[vmm_id] / t_off_s if t_off_s > 0
                           else np.zeros(64),
            "t_on_s"     : t_on_s,
            "t_off_s"    : t_off_s,
        }
    return result


def plot_rate_per_channel_on_off(per_ch_rates, run_no, sg, snt,
                                  good_channels=None,
                                  ch_outliers=None,
                                  noisy_off_spill=None,
                                  out_dir=None, show=True, rate_tag=""):
    """
    Per-channel spill-on vs spill-off rate for every detector VMM.

    One figure per detector group, two panels per VMM:
      Left  : spill-off rate per channel (Hz, log) — noise floor.
      Right : spill-on  rate per channel (Hz, log) — beam signal.

    Bar colours encode channel status:
      firebrick  (off) / steelblue (on) : good channel, used in analysis
      darkorange                         : masked as noisy (off-spill rate)
      lightgrey                          : masked as dead (low total hits)

    Parameters
    ----------
    good_channels : dict {vmm_id: array} or None
        Output of build_good_channels().
    ch_outliers : dict or None
        Output of identify_outlier_channels() — provides dead flags.
    noisy_off_spill : dict or None
        Output of identify_noisy_from_off_spill() — provides noisy flags.
    """
    from vmm_io import get_connected_channels

    for key, cfg in vmm_mapping.items():
        if key == "trigger":
            continue
        group_vmm_ids = [v for v in cfg["vmm_ids"] if v in per_ch_rates]
        if not group_vmm_ids:
            continue
        group_name = cfg.get("name", key)

        n   = len(group_vmm_ids)
        fig, axes = plt.subplots(n, 2, figsize=(18, 3.5 * n),
                                  squeeze=False)

        for row_idx, vmm_id in enumerate(group_vmm_ids):
            connected   = get_connected_channels(vmm_id)
            display_chs = np.array(connected if connected is not None
                                   else list(range(64)))

            r_off = per_ch_rates[vmm_id]["rate_off_hz"][display_chs]
            r_on  = per_ch_rates[vmm_id]["rate_on_hz"][display_chs]

            # ── Per-channel colour: good / noisy / dead ────────
            dead_set  = set()
            noisy_set = set()
            if ch_outliers is not None and vmm_id in ch_outliers:
                dead_set = set(ch_outliers[vmm_id]["dead"])
            if noisy_off_spill is not None and vmm_id in noisy_off_spill:
                noisy_set = set(noisy_off_spill[vmm_id]["noisy"])

            ax_off, ax_on = axes[row_idx]

            for ax, rates, good_color, panel_label in [
                (ax_off, r_off, "firebrick",  "Spill-OFF rate (noise floor)"),
                (ax_on,  r_on,  "steelblue",  "Spill-ON rate  (beam signal)"),
            ]:
                plot_r = np.where(rates == 0, 0.01, rates)

                bar_colors = []
                for ch in display_chs:
                    if ch in noisy_set:
                        bar_colors.append("darkorange")
                    elif ch in dead_set:
                        bar_colors.append("lightgrey")
                    else:
                        bar_colors.append(good_color)

                ax.bar(display_chs, plot_r, color=bar_colors,
                       alpha=0.85, width=0.8)

                # Median and threshold lines
                good_mask = np.array([ch not in dead_set and ch not in noisy_set
                                      for ch in display_chs])
                good_rates = rates[good_mask]
                med = (float(np.median(good_rates[good_rates > 0]))
                       if (good_rates > 0).any() else 0)
                if med > 0:
                    ax.axhline(med, color="k", ls=":", lw=1,
                               label=f"median (good) = {med:.2f} Hz")

                # On the spill-off panel draw the noisy threshold so it is
                # clear why some channels are masked and neighbours are not.
                is_off_panel = "OFF" in panel_label
                if is_off_panel and noisy_off_spill is not None and vmm_id in noisy_off_spill:
                    thr = noisy_off_spill[vmm_id].get("threshold_off_hz", 0)
                    if thr > 0:
                        ax.axhline(thr, color="darkorange", ls="--", lw=1.3,
                                   label=f"noisy thr = {thr:.2f} Hz")

                ax.set_yscale("log")
                ax.set_xlim(display_chs[0] - 0.5, display_chs[-1] + 0.5)
                ax.set_xlabel("Channel")
                ax.set_ylabel("Rate (Hz, log)")
                ax.set_title(f"VMM {vmm_id} — {panel_label}")
                ax.grid(True, alpha=0.3, axis="y")

                # Legend: colour patches for each category
                from matplotlib.patches import Patch
                import matplotlib.lines as mlines
                legend_items = [Patch(color=good_color, label="good channel")]
                if noisy_set:
                    legend_items.append(
                        Patch(color="darkorange", label="masked: noisy (off-spill)")
                    )
                if dead_set:
                    legend_items.append(
                        Patch(color="lightgrey", label="masked: dead")
                    )
                if med > 0:
                    legend_items.append(
                        mlines.Line2D([], [], color="k", ls=":",
                                      label=f"median (good) = {med:.2f} Hz")
                    )
                if is_off_panel and noisy_off_spill is not None and vmm_id in noisy_off_spill:
                    thr = noisy_off_spill[vmm_id].get("threshold_off_hz", 0)
                    if thr > 0:
                        legend_items.append(
                            mlines.Line2D([], [], color="darkorange", ls="--",
                                          label=f"noisy thr = {thr:.2f} Hz")
                        )
                ax.legend(handles=legend_items, fontsize=8, loc="upper right")

        fig.suptitle(
            f"{group_name} — per-channel rate  |  "
            f"run {run_no}  sg={sg}  snt={snt}",
            fontweight="bold"
        )
        plt.tight_layout()
        _finish_fig(fig,
                    f"rate_per_channel_{key}_run{run_no}",
                    out_dir, show, rate_tag)


def collect_channel_hit_counts(data_dir, run_nos, all_vmm_ids,
                               n_files=3, file_start=1,
                               max_time_ticks=2e12):
    """
    Aggregate hit counts per channel (0–63) for each VMM across runs.

    Returns
    -------
    dict {vmm_id: np.ndarray of shape (64,)}
    """
    counts = {v: np.zeros(64, dtype=np.int64) for v in all_vmm_ids}
    for run_no in run_nos:
        run_dir = get_run_dir(data_dir, run_no)
        for df in iter_hits_files(run_dir, n_files=n_files,
                                   branches=["time", "vmm", "ch"],
                                   file_start=file_start):
            if df is None or df.empty:
                continue
            df = df[df["time"] < max_time_ticks].copy()
            for vmm_id in all_vmm_ids:
                df_v = df[df["vmm"] == vmm_id]
                if df_v.empty:
                    continue
                ch_idx = df_v["ch"].values
                valid  = (ch_idx >= 0) & (ch_idx < 64)
                np.add.at(counts[vmm_id], ch_idx[valid], 1)
            del df
    return counts


def identify_outlier_channels(counts, trigger_vmm_ids=None,
                               dead_frac=0.1, noisy_factor=5.0):
    """
    Flag dead and noisy channels per detector VMM.

    Thresholds are computed only over physically connected channels
    (get_connected_channels). Trigger VMMs are skipped entirely.

    dead_frac    : channel below dead_frac × median(nonzero) → dead
    noisy_factor : channel above noisy_factor × median(nonzero) → noisy

    Returns
    -------
    dict {vmm_id: {'dead': [...], 'noisy': [...], 'median': float}}
    """
    from vmm_io import get_connected_channels

    trigger_vmm_ids = set(trigger_vmm_ids or [])
    result = {}
    for vmm_id, ch_counts in counts.items():
        if vmm_id in trigger_vmm_ids:
            result[vmm_id] = {"dead": [], "noisy": [], "median": np.nan}
            continue
        connected    = get_connected_channels(vmm_id)
        display_chs  = connected if connected is not None else list(range(64))
        subset       = ch_counts[np.array(display_chs)]
        nonzero      = subset[subset > 0]
        if len(nonzero) == 0:
            result[vmm_id] = {"dead": list(display_chs), "noisy": [],
                              "median": 0.0}
            continue
        median    = float(np.median(nonzero))
        dead_thr  = dead_frac    * median
        noisy_thr = noisy_factor * median
        dead  = [ch for ch in display_chs if ch_counts[ch] <  dead_thr]
        noisy = [ch for ch in display_chs if ch_counts[ch] >  noisy_thr]
        result[vmm_id] = {"dead": dead, "noisy": noisy, "median": median}
    return result


def identify_noisy_from_off_spill(per_ch_rates_list, noisy_factor=3.0,
                                   ch_outliers=None):
    """
    Flag noisy channels using spill-off (noise floor) rate.

    A channel is noisy if its off-spill rate exceeds noisy_factor × median
    off-spill rate among connected, non-dead channels. Self-triggering
    channels appear here but not (or much less) in the spill-on rate.

    Parameters
    ----------
    per_ch_rates_list : dict or list of dicts
        One or more outputs of compute_per_channel_spill_rates.
        When multiple dicts are given the per-channel max off-spill rate
        is used within this set.
    noisy_factor : float
        Threshold multiplier on the median off-spill rate. Default 3.0.
    ch_outliers : dict or None
        Output of identify_outlier_channels. Dead channels are excluded
        from the median computation so their near-zero off-spill rates
        do not pull the threshold down artificially.

    Returns
    -------
    dict {vmm_id: {'noisy': [ch, ...], 'median_off_hz': float,
                   'threshold_off_hz': float}}
    """
    from vmm_io import get_connected_channels

    if isinstance(per_ch_rates_list, dict):
        per_ch_rates_list = [per_ch_rates_list]

    # Merge across runs: take element-wise max off-spill rate per channel
    all_vmm_ids = set(v for pc in per_ch_rates_list for v in pc)
    merged = {}
    for vmm_id in all_vmm_ids:
        arrays = [pc[vmm_id]["rate_off_hz"]
                  for pc in per_ch_rates_list if vmm_id in pc]
        merged[vmm_id] = np.maximum.reduce(arrays)  # shape (64,)

    result = {}
    for vmm_id, r_off_all in merged.items():
        connected = get_connected_channels(vmm_id)
        if connected is None:
            connected = list(range(64))

        # Exclude dead channels from median — their near-zero off-spill rates
        # would otherwise pull the threshold down and cause inconsistent masking
        # (two channels with similar rates, one masked one not).
        dead = (set(ch_outliers[vmm_id]["dead"])
                if ch_outliers and vmm_id in ch_outliers else set())
        active = [ch for ch in connected if ch not in dead]

        r_off   = r_off_all[np.array(active)]
        nonzero = r_off[r_off > 0]
        if len(nonzero) == 0:
            result[vmm_id] = {"noisy": [], "median_off_hz": 0.0,
                              "threshold_off_hz": 0.0}
            continue
        median    = float(np.median(nonzero))
        threshold = noisy_factor * median
        noisy     = [ch for ch in active if r_off_all[ch] > threshold]
        result[vmm_id] = {"noisy": noisy, "median_off_hz": median,
                          "threshold_off_hz": threshold}
    return result


def build_good_channels(ch_outliers, noisy_off_spill=None):
    """
    Build the set of good channels per VMM: connected − dead − noisy.

    Parameters
    ----------
    ch_outliers : dict
        Output of identify_outlier_channels (dead flags from total hits).
    noisy_off_spill : dict or None
        Output of identify_noisy_from_off_spill.
        If None, falls back to the noisy flags already in ch_outliers.

    Returns
    -------
    dict {vmm_id: np.ndarray of good channel indices (sorted)}
    """
    from vmm_io import get_connected_channels

    good = {}
    for vmm_id, info in ch_outliers.items():
        connected = get_connected_channels(vmm_id)
        if connected is None:
            connected = list(range(64))

        bad = set(info["dead"])
        if noisy_off_spill is not None and vmm_id in noisy_off_spill:
            bad |= set(noisy_off_spill[vmm_id]["noisy"])
        else:
            bad |= set(info["noisy"])

        good[vmm_id] = np.array(sorted(ch for ch in connected
                                        if ch not in bad))
    return good


def build_good_channels_per_config(per_ch_by_run, ch_outliers,
                                    noisy_factor=3.0):
    """
    Build a per-(sg, snt) good-channel mask from diagnostic runs.

    Diagnostic runs sharing the same (sg, snt) are merged by element-wise
    max off-spill rate before thresholding, so a channel is only excluded
    for the gain/peaking-time setting where it is actually noisy — not
    globally across all configurations.

    Parameters
    ----------
    per_ch_by_run : dict {run_no: (sg, snt, per_ch_dict)}
        Output stored during the per-channel diagnostic loop in main().
    ch_outliers : dict
        Output of identify_outlier_channels() — dead flags from total hits.
    noisy_factor : float
        Threshold multiplier on the median off-spill rate. Default 3.0.

    Returns
    -------
    good_by_config  : dict {(sg, snt): good_channels_dict}
    noisy_by_config : dict {(sg, snt): noisy_off_spill_dict}
    """
    config_runs = {}
    for run_no, (sg, snt, per_ch) in per_ch_by_run.items():
        config_runs.setdefault((sg, snt), []).append(per_ch)

    good_by_config  = {}
    noisy_by_config = {}
    for (sg, snt), per_ch_list in config_runs.items():
        noisy_off = identify_noisy_from_off_spill(
            per_ch_list, noisy_factor=noisy_factor, ch_outliers=ch_outliers
        )
        good_by_config[(sg, snt)]  = build_good_channels(ch_outliers, noisy_off)
        noisy_by_config[(sg, snt)] = noisy_off
    return good_by_config, noisy_by_config


def plot_hits_per_channel_all_vmms(counts, outliers,
                                    dead_frac=0.1, noisy_factor=5.0,
                                    out_dir=None, show=True, rate_tag="",
                                    verbose=True):
    """
    Diagnostic plots restricted to physically connected channels only.

    Detector VMMs (all 64 ch connected) — two panels per VMM:
      Left  : hits per channel (bar, channel-number order).
      Right : rank plot — channels sorted lowest→highest hit count.
              Natural gaps between dead / normal / noisy clusters
              are immediately visible as vertical jumps.

    Trigger VMMs (few known channels) — single bar panel per VMM;
      channels not listed in vmm_mapping are excluded entirely.

    Colour coding (detector VMMs)
    -------------
    steelblue  : normal
    firebrick  : dead  (< dead_frac × median of nonzero channels)
    darkorange : noisy (> noisy_factor × median)
    """
    from vmm_io import get_connected_channels

    for key, cfg in vmm_mapping.items():
        group_vmm_ids = [v for v in cfg["vmm_ids"] if v in counts]
        if not group_vmm_ids:
            continue
        group_name     = cfg.get("name", key)
        is_trigger_grp = (key == "trigger")

        # Trigger groups: 1 column (bar only); detectors: 2 columns
        n_cols       = 1 if is_trigger_grp else 2
        width_ratios = None if is_trigger_grp else [2, 1]
        n            = len(group_vmm_ids)
        fig_h        = 3.5 * n
        fig_w        = 10 if is_trigger_grp else 18

        gs_kw = {"width_ratios": width_ratios} if width_ratios else {}
        fig, axes = plt.subplots(n, n_cols,
                                  figsize=(fig_w, fig_h),
                                  squeeze=False,
                                  gridspec_kw=gs_kw)

        for row_idx, vmm_id in enumerate(group_vmm_ids):
            # ── connected channels for this VMM ────────────────
            connected   = get_connected_channels(vmm_id)
            display_chs = np.array(connected if connected is not None
                                   else list(range(64)))
            ch_subset   = counts[vmm_id][display_chs].astype(float)

            info   = outliers[vmm_id]
            dead   = set(info["dead"])
            noisy  = set(info["noisy"])
            median = info["median"]

            bar_colors = []
            for ch in display_chs:
                if ch in noisy:
                    bar_colors.append("darkorange")
                elif ch in dead:
                    bar_colors.append("firebrick")
                else:
                    bar_colors.append("steelblue")

            plot_subset = np.where(ch_subset == 0, 0.5, ch_subset)

            # ── LEFT / ONLY panel: bar chart ───────────────────
            ax_bar = axes[row_idx, 0]

            if is_trigger_grp:
                # Sparse x-axis: actual channel numbers, not dense 0-63
                bar_width = min(3.0, 0.6 * (display_chs[-1] - display_chs[0])
                                / max(len(display_chs) - 1, 1))
                ax_bar.bar(display_chs, plot_subset,
                           color=bar_colors, alpha=0.8, width=bar_width)
                ax_bar.set_xlim(display_chs[0] - 5, display_chs[-1] + 5)
                ax_bar.set_xticks(display_chs)
                ax_bar.set_xticklabels([f"ch {c}" for c in display_chs],
                                       fontsize=9)
                ax_bar.set_title(
                    f"VMM {vmm_id}  (trigger — "
                    f"connected: {list(display_chs)})"
                )
            else:
                ax_bar.bar(display_chs, plot_subset,
                           color=bar_colors, alpha=0.8, width=0.8)
                ax_bar.set_xlim(-0.5, 63.5)
                ax_bar.set_xlabel("Channel")
                ax_bar.set_title(
                    f"VMM {vmm_id}  |  dead={sorted(dead)}  "
                    f"noisy={sorted(noisy)}"
                )

                if not np.isnan(median) and median > 0:
                    ax_bar.axhline(dead_frac * median, color="firebrick",
                                   ls="--", lw=1.2,
                                   label=f"dead < {dead_frac:.0%}·med")
                    ax_bar.axhline(noisy_factor * median, color="darkorange",
                                   ls="--", lw=1.2,
                                   label=f"noisy > {noisy_factor:.0f}×med")
                    ax_bar.axhline(median, color="limegreen", ls=":", lw=1,
                                   label=f"median = {median:.0f}")

            ax_bar.set_yscale("log")
            ax_bar.set_ylabel("Hits (log)")
            ax_bar.legend(fontsize=8, loc="upper right")
            ax_bar.grid(True, alpha=0.3, axis="y")

            # ── RIGHT panel: rank plot (detector VMMs only) ────
            if not is_trigger_grp:
                ax_rank = axes[row_idx, 1]

                sort_idx    = np.argsort(ch_subset)
                sorted_cnts = ch_subset[sort_idx]
                rank_colors = [bar_colors[i] for i in sort_idx]
                rank_plot   = np.where(sorted_cnts == 0, 0.5, sorted_cnts)
                n_disp      = len(display_chs)

                ax_rank.barh(np.arange(n_disp), rank_plot,
                             color=rank_colors, alpha=0.8, height=0.8)

                if not np.isnan(median) and median > 0:
                    ax_rank.axvline(dead_frac * median, color="firebrick",
                                    ls="--", lw=1.2)
                    ax_rank.axvline(noisy_factor * median, color="darkorange",
                                    ls="--", lw=1.2)
                    ax_rank.axvline(median, color="limegreen", ls=":", lw=1)

                ax_rank.set_xscale("log")
                ax_rank.set_xlabel("Hits (log)")
                ax_rank.set_ylabel("Channel rank (sorted)")
                ax_rank.set_title("Sorted by hit count")
                ax_rank.set_ylim(-0.5, n_disp - 0.5)
                ax_rank.grid(True, alpha=0.3, axis="x")

                # Annotate flagged channels by name on the rank plot
                for rank, i in enumerate(sort_idx):
                    ch = display_chs[i]
                    if ch in dead or ch in noisy:
                        ax_rank.text(rank_plot[rank] * 1.15, rank,
                                     f"ch{ch}", fontsize=6, va="center")

        fig.suptitle(f"{group_name} — hits per connected channel  "
                     f"[dead < {dead_frac:.0%}·med  |  "
                     f"noisy > {noisy_factor:.0f}×med]",
                     fontweight="bold")
        plt.tight_layout()
        _finish_fig(fig, f"hits_per_channel_{key}", out_dir, show, rate_tag)

        if verbose:
            print(f"\n{'='*50}")
            print(f"Group: {group_name}")
            for vmm_id in group_vmm_ids:
                info = outliers[vmm_id]
                connected = get_connected_channels(vmm_id)
                if is_trigger_grp:
                    print(f"  VMM {vmm_id}: connected = {connected}")
                else:
                    print(f"  VMM {vmm_id}: dead={info['dead']}  "
                          f"noisy={info['noisy']}  "
                          f"median={info['median']:.0f} hits")


def compute_run_duration(df_hits):
    """
    Compute run duration in seconds from timestamp range.
    Returns 0.0 if df_hits is None or empty.
    """
    if df_hits is None or df_hits.empty:
        return 0.0
    t = df_hits["time"].values
    return (t[-1] - t[0]) * S_PER_TICK


def compute_hit_rate(df_hits, vmm_id):
    """
    Compute hit rate in Hz for a single VMM.
    Uses total hits / time range of that VMM's hits.
    """
    t = df_hits[df_hits["vmm"] == vmm_id]["time"].values
    if len(t) < 2:
        return np.nan
    duration = (t[-1] - t[0]) * S_PER_TICK
    if duration <= 0:
        return np.nan
    return len(t) / duration


def compute_inter_event_times(df_hits, vmm_id,
                               min_dt_us=None,
                               max_dt_us=1000,
                               per_channel=True):
    """
    Compute inter-event time differences.

    Parameters
    ----------
    per_channel : bool
        If True, compute Δt within each channel separately
        then combine. This removes cross-channel coincidence
        spikes (e.g. two scintillator channels on same VMM
        firing for the same particle).
        If False, compute across all hits on the VMM.
    min_dt_us : float or None
        If set, apply a minimum Δt cut in microseconds.
    """
    df_vmm = df_hits[df_hits["vmm"] == vmm_id]
    if df_vmm.empty:
        return np.array([])

    if per_channel and "ch" in df_vmm.columns:
        # Compute Δt within each channel independently
        all_dt = []
        for ch_id, df_ch in df_vmm.groupby("ch"):
            t = df_ch["time"].values
            if len(t) < 2:
                continue
            dt_ns = np.diff(t) * NS_PER_TICK
            all_dt.append(dt_ns)

        if not all_dt:
            return np.array([])
        dt_ns = np.concatenate(all_dt)

    else:
        t     = df_vmm["time"].values
        dt_ns = np.diff(t) * NS_PER_TICK

    dt_us = dt_ns / 1000.0
    dt_us = dt_us[dt_us > 0]

    if min_dt_us is not None:
        dt_us = dt_us[dt_us > min_dt_us]

    dt_us = dt_us[dt_us < max_dt_us]
    return dt_us

def fit_exponential(dt_us, n_bins=200, auto_offset=True):
    """
    Fit exponential to inter-event time distribution.
    Supports shifted exponential: P(Δt) = λ·exp(-λ·(Δt - τ))
    where τ is the minimum observable Δt (delay module offset).

    Parameters
    ----------
    auto_offset : bool
        If True, automatically detect the offset τ as the
        minimum populated bin. Set False for direct trigger data.
    """
    if len(dt_us) < 50:
        return np.nan, np.nan, np.nan, np.nan, None, None

    counts, edges = np.histogram(dt_us, bins=n_bins)
    centers       = 0.5 * (edges[:-1] + edges[1:])
    bin_width     = edges[1] - edges[0]

    # Auto-detect offset: first bin with significant counts
    if auto_offset:
        threshold = counts.max() * 0.05  # 5% of peak
        first_populated = np.where(counts > threshold)[0]
        if len(first_populated) == 0:
            return np.nan, np.nan, np.nan, np.nan, None, None
        tau = centers[first_populated[0]]
    else:
        tau = 0.0

    print(f"  Detected offset τ = {tau:.1f} μs")

    # Shift the distribution
    dt_shifted = dt_us[dt_us > tau] - tau
    if len(dt_shifted) < 50:
        return np.nan, np.nan, np.nan, np.nan, None, None

    counts_s, edges_s = np.histogram(dt_shifted, bins=n_bins)
    centers_s         = 0.5 * (edges_s[:-1] + edges_s[1:])
    bin_width_s       = edges_s[1] - edges_s[0]

    mask = counts_s >= 5

    def exp_model(x, lam):
        A = len(dt_shifted) * bin_width_s * lam
        return A * np.exp(-lam * x)

    lam0 = 1.0 / dt_shifted.mean()

    try:
        popt, pcov = curve_fit(
            exp_model,
            centers_s[mask],
            counts_s[mask],
            p0=[lam0],
            sigma=np.sqrt(np.maximum(counts_s[mask], 1)),
            absolute_sigma=True
        )
        lam     = popt[0]
        lam_err = np.sqrt(pcov[0, 0])

        residuals = counts_s[mask] - exp_model(centers_s[mask], lam)
        chi2_ndf  = (np.sum(residuals**2 / np.maximum(counts_s[mask], 1))
                     / (mask.sum() - 1))

        rate_hz     = lam * 1e6
        rate_err_hz = lam_err * 1e6

        fit_x = np.linspace(0, centers_s[-1], 500)
        fit_y = exp_model(fit_x, lam)

        return rate_hz, rate_err_hz, chi2_ndf, tau, fit_x, fit_y

    except RuntimeError:
        return np.nan, np.nan, np.nan, tau, None, None


def fit_exponential_with_offset(dt_us, n_bins=200):
    """
    Fit exponential with minimum time offset for bunched beams.

    P(Δt) = λ · exp(-λ · (Δt - t_min))  for Δt > t_min

    This correctly handles beam structure where particles
    cannot arrive closer than t_min apart.

    Returns
    -------
    rate_hz : float
        Beam repetition rate in Hz (1/mean_spacing).
    t_min_us : float
        Minimum inter-event time in μs (bunch spacing floor).
    chi2_ndf : float
        Reduced chi-squared.
    """
    if len(dt_us) < 50:
        return np.nan, np.nan, np.nan, None, None

    counts, edges = np.histogram(dt_us, bins=n_bins)
    centers       = 0.5 * (edges[:-1] + edges[1:])
    bin_width     = edges[1] - edges[0]

    # Estimate t_min as the first bin with significant counts
    # (where distribution rises above noise floor)
    cumulative  = np.cumsum(counts)
    total       = cumulative[-1]
    t_min_guess = centers[np.argmax(cumulative > 0.01 * total)]

    mask = (counts >= 5) & (centers > t_min_guess)
    if mask.sum() < 10:
        return np.nan, np.nan, np.nan, None, None

    def exp_offset_model(x, lam, t_min):
        A = len(dt_us) * bin_width * lam
        return A * np.exp(-lam * (x - t_min))

    try:
        popt, pcov = curve_fit(
            exp_offset_model,
            centers[mask],
            counts[mask],
            p0=[1.0/dt_us.mean(), t_min_guess],
            sigma=np.sqrt(np.maximum(counts[mask], 1)),
            absolute_sigma=True,
            bounds=([0, 0], [np.inf, centers[-1]])
        )
        lam, t_min = popt
        lam_err    = np.sqrt(pcov[0, 0])

        residuals = counts[mask] - exp_offset_model(
            centers[mask], lam, t_min
        )
        chi2_ndf = (
            np.sum(residuals**2 / np.maximum(counts[mask], 1))
            / (mask.sum() - 2)
        )

        rate_hz     = lam * 1e6
        rate_err_hz = lam_err * 1e6

        fit_x = np.linspace(t_min, centers[-1], 500)
        fit_y = exp_offset_model(fit_x, lam, t_min)

        return rate_hz, t_min, chi2_ndf, fit_x, fit_y

    except RuntimeError:
        return np.nan, np.nan, np.nan, None, None


def compute_rates_all_runs(df_run_scan, data_dir,
                            sng0_runs, trigger_vmms,
                            detector_vmms,
                            trigger_ref_channels,
                            root_file_index=1):
    results = []

    for run_no in sng0_runs:
        # Use load_sorted_hits — applies timestamp cut automatically
        df_hits = load_sorted_hits(
            data_dir, run_no,
            branches=["time", "vmm", "ch"],
            root_file_index=root_file_index
        )
        if df_hits is None or df_hits.empty:
            print(f"  Run {run_no}: no data, skipping")
            continue

        t_all    = df_hits["time"].values
        duration = (t_all[-1] - t_all[0]) * S_PER_TICK

        if duration <= 0:
            print(f"  Run {run_no}: zero duration, skipping")
            continue

        sg  = df_run_scan.loc[
            df_run_scan["run_no"] == run_no, "sg"
        ].iloc[0]
        snt = df_run_scan.loc[
            df_run_scan["run_no"] == run_no, "snt"
        ].iloc[0]

        print(f"\nRun {run_no} (sg={sg} snt={snt}) "
              f"duration={duration:.1f}s")

        # Trigger rate — one reference channel per VMM
        for vmm_id, ch_id in trigger_ref_channels.items():
            n    = ((df_hits["vmm"] == vmm_id) &
                    (df_hits["ch"]  == ch_id)).sum()
            rate = n / duration

            print(f"  VMM {vmm_id} ch{ch_id} "
                  f"(trigger ref): {n} hits  {rate:.1f} Hz")

            results.append({
                "run_no"  : run_no,
                "sg"      : sg,
                "snt"     : snt,
                "vmm_id"  : vmm_id,
                "ch"      : ch_id,
                "type"    : "trigger",
                "n_hits"  : n,
                "duration": duration,
                "rate_hz" : rate
            })

        # Detector rates
        for vmm_id in detector_vmms:
            n    = (df_hits["vmm"] == vmm_id).sum()
            rate = n / duration

            print(f"  VMM {vmm_id} (detector): "
                  f"{n} hits  {rate:.1f} Hz")

            results.append({
                "run_no"  : run_no,
                "sg"      : sg,
                "snt"     : snt,
                "vmm_id"  : vmm_id,
                "ch"      : -1,
                "type"    : "detector",
                "n_hits"  : n,
                "duration": duration,
                "rate_hz" : rate
            })

    return pd.DataFrame(results)

def plot_inter_event_distribution(dt_us, vmm_id, run_no,
                                    rate_hz, chi2_ndf, tau,
                                    fit_x, fit_y,
                                    n_bins=200,
                                    out_dir=None, show=True, rate_tag=""):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))

    # Original distribution
    axes[0].hist(dt_us, bins=n_bins,
                  color="steelblue", alpha=0.7)
    if tau > 0:
        axes[0].axvline(tau, color="orange", linestyle="--",
                         label=f"Offset τ={tau:.0f} μs")
    axes[0].set_xlabel("Δt (μs)")
    axes[0].set_ylabel("Counts")
    axes[0].set_title(f"VMM {vmm_id} Run {run_no} — raw")
    axes[0].legend(fontsize=8)

    # Shifted distribution with fit — linear
    dt_shifted = dt_us[dt_us > tau] - tau
    axes[1].hist(dt_shifted, bins=n_bins,
                  color="steelblue", alpha=0.7, label="Shifted data")
    if fit_x is not None:
        axes[1].plot(fit_x, fit_y, color="red", linewidth=2,
                      label=f"λ={rate_hz:.0f}±"
                            f"{0:.0f} Hz\nχ²/ndf={chi2_ndf:.2f}")
    axes[1].set_xlabel("Δt - τ (μs)")
    axes[1].set_ylabel("Counts")
    axes[1].set_title("Shifted — linear")
    axes[1].legend(fontsize=8)

    # Shifted — log scale
    axes[2].hist(dt_shifted, bins=n_bins,
                  color="steelblue", alpha=0.7)
    if fit_x is not None:
        axes[2].plot(fit_x, fit_y, color="red", linewidth=2,
                      label=f"λ={rate_hz:.0f} Hz")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("Δt - τ (μs)")
    axes[2].set_ylabel("Counts (log)")
    axes[2].set_title("Shifted — log scale")
    axes[2].legend(fontsize=8)

    plt.suptitle(
        f"Inter-event time — VMM {vmm_id} Run {run_no}\n"
        f"Offset τ={tau:.0f} μs  |  "
        f"Fitted rate={rate_hz:.0f} Hz"
    )
    plt.tight_layout()
    _finish_fig(fig,
                f"inter_event_dist_vmm{vmm_id}_run{run_no}",
                out_dir, show, rate_tag)


def inspect_trigger_timing(df_hits, trigger_vmm=0, trigger_ch=48,
                            window_ns=600,
                            out_dir=None, show=True, rate_tag=""):
    """
    Diagnostic checks on the trigger timestamp stream before
    running coincidence efficiency.

    Checks
    ------
    1. Duplicate timestamps — same timestamp appearing more than once
       on the trigger channel (readout artifact / double-recording).
    2. Trigger pileup — consecutive trigger pairs closer than
       `window_ns`, meaning their coincidence windows would overlap
       and the same detector hit could be claimed by both triggers.

    Prints a summary and plots the inter-trigger time distribution.
    """
    t_trig = df_hits[
        (df_hits["vmm"] == trigger_vmm) &
        (df_hits["ch"]  == trigger_ch)
    ]["time"].values

    if len(t_trig) == 0:
        print(f"No hits on VMM {trigger_vmm} ch {trigger_ch}")
        return

    t_trig  = np.sort(t_trig)
    n_total = len(t_trig)

    # ── 1. Duplicate timestamps ─────────────────────────────
    _, counts    = np.unique(t_trig, return_counts=True)
    n_duplicates = int((counts > 1).sum())
    n_dup_hits   = int((counts[counts > 1] - 1).sum())

    print(f"\nTrigger VMM {trigger_vmm} ch {trigger_ch}:")
    print(f"  Total hits          : {n_total}")
    print(f"  Unique timestamps   : {len(counts)}")
    print(f"  Duplicate timestamps: {n_duplicates} "
          f"({100 * n_duplicates / len(counts):.2f}% of unique)")
    print(f"  Extra hits (dups)   : {n_dup_hits} "
          f"({100 * n_dup_hits / n_total:.2f}% of all hits)")

    # ── 2. Inter-trigger time vs window ────────────────────
    dt_trig_ns = np.diff(t_trig) * NS_PER_TICK
    n_pileup   = int((dt_trig_ns < window_ns).sum())
    pct_pileup = 100 * n_pileup / len(dt_trig_ns)

    print(f"\n  Coincidence window  : {window_ns} ns")
    print(f"  Consecutive trigger pairs closer than window: "
          f"{n_pileup} / {len(dt_trig_ns)} "
          f"({pct_pileup:.2f}%)")
    print(f"  Median inter-trigger Δt : "
          f"{np.median(dt_trig_ns):.0f} ns  "
          f"= {np.median(dt_trig_ns)/1e6:.2f} ms")
    print(f"  Min inter-trigger Δt   : "
          f"{dt_trig_ns.min():.0f} ns")

    # ── Plot inter-trigger time distribution ───────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Zoom: 0–5000 ns to see the pileup / delay-module region
    dt_zoom = dt_trig_ns[dt_trig_ns < 5000]
    axes[0].hist(dt_zoom, bins=200, color="steelblue", alpha=0.7)
    axes[0].axvline(window_ns, color="red", linestyle="--",
                    linewidth=1.5,
                    label=f"proposed window = {window_ns} ns")
    axes[0].set_xlabel("Inter-trigger Δt (ns)")
    axes[0].set_ylabel("Counts")
    axes[0].set_title("Inter-trigger Δt — zoom 0–5000 ns")
    axes[0].legend()

    # Full range in ms
    dt_ms = dt_trig_ns / 1e6
    axes[1].hist(dt_ms[dt_ms < 20], bins=200,
                 color="steelblue", alpha=0.7)
    axes[1].set_xlabel("Inter-trigger Δt (ms)")
    axes[1].set_ylabel("Counts")
    axes[1].set_title("Inter-trigger Δt — full range")

    fig.suptitle(
        f"Trigger timing — VMM {trigger_vmm} ch {trigger_ch}",
        fontweight="bold"
    )
    plt.tight_layout()
    _finish_fig(fig,
                f"trigger_timing_vmm{trigger_vmm}_ch{trigger_ch}",
                out_dir, show, rate_tag)

    return {
        "n_total"     : n_total,
        "n_duplicates": n_duplicates,
        "n_dup_hits"  : n_dup_hits,
        "n_pileup"    : n_pileup,
        "pct_pileup"  : pct_pileup,
    }


def compute_rate_vs_time(df_hits, vmm_id, bin_width_s=10.0,
                          channel=None):
    """
    Compute hit rate in kHz as a function of time within a run.

    Parameters
    ----------
    vmm_id : int
        VMM to compute rate for.
    bin_width_s : float
        Width of each time bin in seconds.
    channel : int or None
        If set, restrict to a single channel (e.g. trigger ref ch).

    Returns
    -------
    t_centers : np.ndarray
        Bin centres in seconds from run start.
    rates_khz : np.ndarray
        Hit rate in kHz for each bin.
    """
    mask = df_hits["vmm"] == vmm_id
    if channel is not None:
        mask &= df_hits["ch"] == channel

    t = df_hits[mask]["time"].values
    if len(t) < 2:
        return np.array([]), np.array([])

    t_s      = t * S_PER_TICK
    t_start  = t_s.min()
    t_end    = t_s.max()

    bins     = np.arange(t_start, t_end + bin_width_s, bin_width_s)
    counts, edges = np.histogram(t_s, bins=bins)
    t_centers     = 0.5 * (edges[:-1] + edges[1:]) - t_start
    rates_khz     = counts / bin_width_s / 1e3

    return t_centers, rates_khz


def plot_rate_vs_time(data_dir, df_run_scan, sng0_runs,
                       vmm_groups, trigger_ref_channels,
                       bin_width_s=10.0,
                       root_file_index=1,
                       out_dir=None, show=True, rate_tag=""):
    """
    Plot hit rate (kHz) vs time for each run, grouped by detector.

    One figure per detector group — each figure has:
      - Top panel  : trigger rate (VMM 0 ref channel) — beam monitor
      - One panel per VMM in the group

    Runs are shown as separate lines coloured by configuration (sg/snt).
    This reveals beam intensity variations within and across runs.

    Parameters
    ----------
    bin_width_s : float
        Time bin width in seconds. Use 5–30 s depending on statistics.
    """
    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    trig_vmm = list(trigger_ref_channels.keys())[0]
    trig_ch  = trigger_ref_channels[trig_vmm]

    # Colour each run by its configuration label
    configs      = df_run_scan[
        df_run_scan["run_no"].isin(sng0_runs)
    ][["run_no", "sg", "snt"]].drop_duplicates()
    config_labels = {
        row["run_no"]: f"sg={row['sg']:.1f} snt={row['snt']:.0f}"
        for _, row in configs.iterrows()
    }
    unique_cfgs  = configs[["sg", "snt"]].drop_duplicates()
    cmap         = plt.cm.tab10(
        np.linspace(0, 1, len(unique_cfgs))
    )
    cfg_colors   = {
        (r["sg"], r["snt"]): cmap[i]
        for i, (_, r) in enumerate(unique_cfgs.iterrows())
    }
    run_colors   = {
        row["run_no"]: cfg_colors[(row["sg"], row["snt"])]
        for _, row in configs.iterrows()
    }

    for group in vmm_groups:
        det_vmm_ids = group["vmm_ids"]
        group_label = group["label"]

        n_panels = 1 + len(det_vmm_ids)   # trigger + one per VMM
        fig, axes = plt.subplots(
            n_panels, 1,
            figsize=(14, 3 * n_panels),
            sharex=False
        )
        if n_panels == 1:
            axes = [axes]

        seen_labels = set()

        for run_no in sng0_runs:
            df_hits = load_sorted_hits(
                data_dir, run_no,
                branches=["time", "vmm", "ch"],
                root_file_index=root_file_index
            )
            if df_hits is None or df_hits.empty:
                continue

            color = run_colors.get(run_no, "gray")
            label = config_labels.get(run_no, f"run {run_no}")
            # Avoid duplicate legend entries for same config
            legend_label = label if label not in seen_labels else None
            seen_labels.add(label)

            # Top panel: trigger rate
            t_c, r_khz = compute_rate_vs_time(
                df_hits, trig_vmm,
                bin_width_s=bin_width_s,
                channel=trig_ch
            )
            if len(t_c):
                axes[0].plot(t_c, r_khz,
                             color=color, linewidth=1.2,
                             alpha=0.8, label=legend_label)

            # One panel per detector VMM
            for ax, vmm_id in zip(axes[1:], det_vmm_ids):
                t_c, r_khz = compute_rate_vs_time(
                    df_hits, vmm_id,
                    bin_width_s=bin_width_s
                )
                if len(t_c):
                    ax.plot(t_c, r_khz,
                            color=color, linewidth=1.2,
                            alpha=0.8, label=legend_label)

        # Formatting
        trig_name = vmm_to_detector.get(trig_vmm,
                                         f"VMM {trig_vmm}")
        axes[0].set_ylabel("Rate (kHz)")
        axes[0].set_title(
            f"Trigger — VMM {trig_vmm} ch {trig_ch}"
        )
        axes[0].legend(fontsize=9, loc="upper right")
        axes[0].grid(True, alpha=0.3)

        for ax, vmm_id in zip(axes[1:], det_vmm_ids):
            det_name = vmm_to_detector.get(vmm_id, "")
            ax.set_ylabel("Rate (kHz)")
            ax.set_title(f"VMM {vmm_id} — {det_name}")
            ax.grid(True, alpha=0.3)

        axes[-1].set_xlabel("Time from run start (s)")

        fig.suptitle(
            f"{group_label} — hit rate vs time\n"
            f"bin width = {bin_width_s:.0f} s  |  "
            f"colours = configurations",
            fontweight="bold"
        )
        plt.tight_layout()
        stem = f"rate_vs_time_{group_label.replace(' ', '_')}"
        _finish_fig(fig, stem, out_dir, show, rate_tag)


def plot_trigger_rate_ms(df_hits, run_no,
                          trigger_vmm=0, trigger_ch=48,
                          bin_width_ms=1.0,
                          out_dir=None, show=True, rate_tag=""):
    """
    Plot trigger count per ms as a function of time.

    With 1 ms bins on a bunched SPS beam the spill structure
    (flat top + ramp up/down + inter-spill gaps) becomes visible.

    With bin_width_ms=1: y-axis = triggers/ms = kHz directly.
    Narrow the bin to see intra-spill structure; widen it to
    see multi-spill patterns over the full run.

    Parameters
    ----------
    bin_width_ms : float
        Bin width in milliseconds. Default 1 ms.
    """
    t = df_hits[
        (df_hits["vmm"] == trigger_vmm) &
        (df_hits["ch"]  == trigger_ch)
    ]["time"].values

    if len(t) < 2:
        print(f"No trigger hits on VMM {trigger_vmm} ch {trigger_ch}")
        return

    t_ms     = t * NS_PER_TICK / 1e6          # ticks → ms
    t_start  = t_ms.min()
    t_ms_rel = t_ms - t_start                 # ms from run start

    bins          = np.arange(0, t_ms_rel.max() + bin_width_ms,
                               bin_width_ms)
    counts, edges = np.histogram(t_ms_rel, bins=bins)
    t_centers_ms  = 0.5 * (edges[:-1] + edges[1:])
    t_centers_s   = t_centers_ms / 1e3        # also keep seconds axis

    fig, axes = plt.subplots(2, 1, figsize=(14, 7))

    # Full run — shows all spills
    axes[0].plot(t_centers_s, counts / bin_width_ms,
                 color="steelblue", linewidth=0.8)
    axes[0].set_xlabel("Time from run start (s)")
    axes[0].set_ylabel("Triggers / ms")
    axes[0].set_title("Full run — spill structure")
    axes[0].grid(True, alpha=0.3)

    # Zoom into first 5 s — shows individual spill profile
    mask = t_centers_s <= 5.0
    if mask.sum() > 1:
        axes[1].plot(t_centers_s[mask],
                     (counts / bin_width_ms)[mask],
                     color="steelblue", linewidth=0.8)
        axes[1].set_xlabel("Time from run start (s)")
        axes[1].set_ylabel("Triggers / ms")
        axes[1].set_title("Zoom — first 5 s (spill profile detail)")
        axes[1].grid(True, alpha=0.3)

    fig.suptitle(
        f"Trigger rate — VMM {trigger_vmm} ch {trigger_ch}  |  "
        f"Run {run_no}  |  bin = {bin_width_ms:.1f} ms",
        fontweight="bold"
    )
    plt.tight_layout()
    _finish_fig(fig, f"trigger_rate_ms_run{run_no}", out_dir, show, rate_tag)


def plot_rate_overlay_with_trigger(df_hits, run_no, vmm_groups,
                                    trigger_vmm=0, trigger_ch=48,
                                    bin_width_ms=1.0,
                                    out_dir=None, show=True, rate_tag=""):
    """
    Overlay detector hit rate and trigger rate on the same time axis.

    One figure per detector group, one subplot per VMM.
    Each subplot has two y-axes:
      - Left  (blue)   : detector VMM rate in kHz
      - Right (orange) : trigger rate in kHz

    Both are binned at bin_width_ms so the spill structure is
    directly comparable — if detector spikes coincide with trigger
    spikes the VMM is seeing real beam signal.
    """
    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    bin_width_s = bin_width_ms / 1e3

    # Compute trigger rate once — shared across all panels
    t_trig_c, r_trig = compute_rate_vs_time(
        df_hits, trigger_vmm,
        bin_width_s=bin_width_s,
        channel=trigger_ch
    )

    if len(t_trig_c) == 0:
        print("No trigger hits found")
        return

    for group in vmm_groups:
        vmm_ids     = group["vmm_ids"]
        group_label = group["label"]

        n_vmm = len(vmm_ids)
        fig, axes = plt.subplots(
            n_vmm, 1,
            figsize=(14, 3.5 * n_vmm),
            sharex=True
        )
        if n_vmm == 1:
            axes = [axes]

        for ax, vmm_id in zip(axes, vmm_ids):
            det_name = vmm_to_detector.get(vmm_id, "")

            # Detector rate — left axis
            t_det_c, r_det = compute_rate_vs_time(
                df_hits, vmm_id,
                bin_width_s=bin_width_s
            )

            color_det  = "steelblue"
            color_trig = "darkorange"

            ax.plot(t_det_c, r_det,
                    color=color_det, linewidth=0.8,
                    label=f"VMM {vmm_id} — {det_name}")
            ax.set_ylabel("Detector rate (kHz)",
                          color=color_det)
            ax.tick_params(axis="y", labelcolor=color_det)

            # Trigger rate — right axis
            ax_r = ax.twinx()
            ax_r.plot(t_trig_c, r_trig,
                      color=color_trig, linewidth=0.8,
                      alpha=0.7,
                      label=f"Trigger VMM {trigger_vmm} ch {trigger_ch}")
            ax_r.set_ylabel("Trigger rate (kHz)",
                             color=color_trig)
            ax_r.tick_params(axis="y", labelcolor=color_trig)

            # Combined legend
            lines_l, labels_l = ax.get_legend_handles_labels()
            lines_r, labels_r = ax_r.get_legend_handles_labels()
            ax.legend(lines_l + lines_r, labels_l + labels_r,
                      fontsize=9, loc="upper right")
            ax.grid(True, alpha=0.3)

        axes[-1].set_xlabel("Time from run start (s)")

        fig.suptitle(
            f"{group_label} — detector vs trigger rate  |  "
            f"Run {run_no}  |  bin = {bin_width_ms:.1f} ms",
            fontweight="bold"
        )
        plt.tight_layout()
        stem = (f"rate_overlay_run{run_no}_"
                f"{group_label.replace(' ', '_')}")
        _finish_fig(fig, stem, out_dir, show, rate_tag)


def compute_spill_masks(rates_khz, threshold_khz,
                         min_spill_bins=0, max_gap_bins=0):
    """
    Derive spill-on / spill-off boolean masks from a binned rate array.

    Parameters
    ----------
    rates_khz : np.ndarray
        Binned hit or trigger rate in kHz.
    threshold_khz : float
        Bins above this value are classified as spill-on.
    min_spill_bins : int
        Minimum consecutive on-bins to keep an on-region. Shorter
        regions (e.g. startup artifacts) are reclassified as off.
        0 disables.
    max_gap_bins : int
        Maximum consecutive off-bins between two on-regions that will
        be reclassified as on (hole-filling). This handles bins that
        dip below threshold mid-spill due to SPS bunch microstructure
        at low beam intensity. 0 disables.

    Returns
    -------
    on_mask, off_mask : np.ndarray of bool
    """
    on_mask = rates_khz > threshold_khz

    # ── Hole-filling: bridge short off-gaps within spills ──────
    if max_gap_bins > 0:
        padded    = np.concatenate(([False], on_mask, [False]))
        on_starts = np.where(~padded[:-1] &  padded[1:])[0]  # F→T: on-region starts
        on_ends   = np.where( padded[:-1] & ~padded[1:])[0]  # T→F: on-region ends
        # Fill gaps between consecutive on-regions that are short enough.
        # on_mask[on_ends[i] : on_starts[i+1]] is the off-gap between region i and i+1.
        for i in range(len(on_starts) - 1):
            gap_len = on_starts[i + 1] - on_ends[i]
            if gap_len <= max_gap_bins:
                on_mask[on_ends[i] : on_starts[i + 1]] = True

    # ── Remove short on-regions (startup artifacts, noise spikes) ──
    if min_spill_bins > 1:
        padded = np.concatenate(([False], on_mask, [False]))
        starts = np.where(~padded[:-1] &  padded[1:])[0]
        ends   = np.where( padded[:-1] & ~padded[1:])[0]
        for s, e in zip(starts, ends):
            if (e - s) < min_spill_bins:
                on_mask[s:e] = False

    off_mask = ~on_mask
    return on_mask, off_mask


def plot_spill_mask_diagnostic(df_hits, run_no,
                               trigger_vmm=0, trigger_ch=48,
                               spill_threshold_khz=1.0,
                               bin_width_ms=1.0,
                               min_spill_s=1.0,
                               max_gap_s=0.0,
                               out_dir=None, show=True, rate_tag="",
                               verbose=True):
    """
    Step 1 validation: plot 1 ms trigger rate with spill-on/off
    regions shaded.

    Use this on a single run before running all configs to confirm
    the threshold correctly separates beam-on from beam-off.

    Green shading = spill-on (rate > threshold).
    Grey shading  = spill-off (noise floor).
    Red dashed line = threshold.
    """
    t = df_hits[
        (df_hits["vmm"] == trigger_vmm) &
        (df_hits["ch"]  == trigger_ch)
    ]["time"].values

    if len(t) < 2:
        print(f"No hits for VMM {trigger_vmm} ch {trigger_ch}")
        return

    bin_width_s = bin_width_ms * 1e-3
    t_s      = t * S_PER_TICK
    t_start  = t_s.min()
    bins     = np.arange(t_start, t_s.max() + bin_width_s, bin_width_s)
    counts, edges = np.histogram(t_s, bins=bins)
    t_centers_s   = 0.5 * (edges[:-1] + edges[1:]) - t_start
    rate_khz       = counts / bin_width_s / 1e3

    min_spill_bins = max(1, int(min_spill_s / bin_width_s))
    max_gap_bins   = max(0, int(max_gap_s   / bin_width_s))
    on_mask, off_mask = compute_spill_masks(
        rate_khz, spill_threshold_khz, min_spill_bins, max_gap_bins
    )
    n_on  = on_mask.sum()
    n_off = off_mask.sum()
    frac_on = 100 * n_on / len(rate_khz)

    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=False)

    for ax, xlim, title_suffix in zip(
        axes,
        [None, (0, min(10.0, t_centers_s.max()))],
        ["full run", "zoom first 10 s"],
    ):
        ax.plot(t_centers_s, rate_khz,
                color="steelblue", linewidth=0.7, zorder=3)
        ax.axhline(spill_threshold_khz, color="crimson",
                   linewidth=1.2, linestyle="--",
                   label=f"threshold = {spill_threshold_khz:.1f} kHz")

        # Shade contiguous spill-on / spill-off regions —
        # one axvspan per region, not per bin
        for mask_bool, color in [(on_mask, "limegreen"),
                                  (off_mask, "lightgrey")]:
            # Find starts and ends of each contiguous True block
            padded   = np.concatenate(([False], mask_bool, [False]))
            starts   = np.where(~padded[:-1] &  padded[1:])[0]
            ends     = np.where( padded[:-1] & ~padded[1:])[0]
            for s, e in zip(starts, ends):
                ax.axvspan(t_centers_s[s]   - bin_width_s / 2,
                           t_centers_s[e-1] + bin_width_s / 2,
                           alpha=0.25, color=color,
                           linewidth=0, zorder=1)

        if xlim is not None:
            ax.set_xlim(xlim)
        ax.set_ylabel("Trigger rate (kHz)")
        ax.set_title(f"Spill mask — {title_suffix}")
        ax.grid(True, alpha=0.3, zorder=2)
        ax.legend(loc="upper right")

    axes[-1].set_xlabel("Time from run start (s)")
    fig.suptitle(
        f"Run {run_no} — VMM {trigger_vmm} ch {trigger_ch}  |  "
        f"spill-on: {n_on} bins ({frac_on:.1f}%)  "
        f"spill-off: {n_off} bins",
        fontweight="bold"
    )
    plt.tight_layout()
    _finish_fig(fig, f"spill_mask_run{run_no}", out_dir, show, rate_tag)

    if verbose:
        print(f"  Spill-on  bins : {n_on}  ({frac_on:.1f}%)")
        print(f"  Spill-off bins : {n_off}  ({100-frac_on:.1f}%)")
        print(f"  Mean rate (on) : {rate_khz[on_mask].mean():.2f} kHz")
        print(f"  Mean rate (off): {rate_khz[off_mask].mean():.3f} kHz")


def _get_run_time_range(run_dir, n_files, max_time_ticks=2e12,
                        file_start=0):
    """
    Scan all files in a run to find the global timestamp range.
    Reads only the time column — one file at a time.
    Returns (t_start_s, t_end_s) or (None, None) if no data.
    """
    t_min = np.inf
    t_max = -np.inf

    for df in iter_hits_files(run_dir, n_files=n_files,
                               branches=["time"],
                               file_start=file_start):
        if df is None or df.empty:
            continue
        t = df["time"].values
        t = t[t < max_time_ticks]
        if len(t):
            t_min = min(t_min, t.min())
            t_max = max(t_max, t.max())
        del df
    if np.isinf(t_min):
        return None, None
    return t_min * S_PER_TICK, t_max * S_PER_TICK


def _accumulate_run_histograms(run_dir, trig_vmm, trig_ch,
                                detector_vmms, bins, n_files,
                                good_channels=None,
                                max_time_ticks=2e12, file_start=0):
    """
    Iterate over all files in a run one at a time, accumulating
    hit counts into pre-defined histogram bins.
    Each file's DataFrame is deleted immediately after histogramming
    to keep peak memory to a single file at a time.

    Parameters
    ----------
    good_channels : dict {vmm_id: array-like of channel indices} or None
        If provided, only hits on these channels are counted per VMM.
        Falls back to get_connected_channels when None.
    """
    from vmm_io import get_connected_channels

    n_bins      = len(bins) - 1
    trig_counts = np.zeros(n_bins, dtype=np.int64)
    det_counts  = {v: np.zeros(n_bins, dtype=np.int64)
                   for v in detector_vmms}

    for df in iter_hits_files(run_dir, n_files=n_files,
                               branches=["time", "vmm", "ch"],
                               file_start=file_start):
        if df is None or df.empty:
            continue
        df = df[df["time"] < max_time_ticks].copy()
        if df.empty:
            del df
            continue

        t_s = df["time"].values * S_PER_TICK

        # Trigger counts
        m_trig = (df["vmm"] == trig_vmm) & (df["ch"] == trig_ch)
        c, _   = np.histogram(t_s[m_trig.values], bins=bins)
        trig_counts += c

        # Detector counts — good channels only
        for vmm_id in detector_vmms:
            if good_channels is not None and vmm_id in good_channels:
                chs = good_channels[vmm_id]
            else:
                chs = get_connected_channels(vmm_id)
            m = df["vmm"] == vmm_id
            if chs is not None and len(chs) > 0:
                m &= df["ch"].isin(chs)
            c, _ = np.histogram(t_s[m.values], bins=bins)
            det_counts[vmm_id] += c

        del df

    return trig_counts, det_counts


def _accumulate_per_channel_on_off(run_dir, detector_vmms, good_channels,
                                   bins, on_mask, off_mask,
                                   n_files, file_start=0,
                                   max_time_ticks=2e12):
    """
    Third pass: per-channel hit counts split by spill-on / spill-off mask.
    Returns {vmm_id: {"channels": array, "on": array[n_ch], "off": array[n_ch]}}
    where on/off are total hit counts per good channel.
    """
    from vmm_io import get_connected_channels

    chs_per_vmm = {}
    ch_on  = {}
    ch_off = {}
    for vmm_id in detector_vmms:
        if good_channels is not None and vmm_id in good_channels:
            chs = np.array(good_channels[vmm_id])
        else:
            c = get_connected_channels(vmm_id)
            chs = np.arange(64) if c is None else np.array(c)
        chs_per_vmm[vmm_id] = chs
        ch_on[vmm_id]  = np.zeros(len(chs), dtype=np.int64)
        ch_off[vmm_id] = np.zeros(len(chs), dtype=np.int64)

    for df in iter_hits_files(run_dir, n_files=n_files,
                               branches=["time", "vmm", "ch"],
                               file_start=file_start):
        if df is None or df.empty:
            continue
        df = df[df["time"] < max_time_ticks]
        if df.empty:
            del df
            continue

        t_s     = df["time"].values * S_PER_TICK
        bin_idx = np.searchsorted(bins, t_s, side="right").astype(np.intp) - 1
        valid   = (bin_idx >= 0) & (bin_idx < len(on_mask))
        is_on   = np.zeros(len(df), dtype=bool)
        is_on[valid]  = on_mask[bin_idx[valid]]
        is_off  = np.zeros(len(df), dtype=bool)
        is_off[valid] = off_mask[bin_idx[valid]]

        vmm_arr = df["vmm"].values
        ch_arr  = df["ch"].values

        for vmm_id in detector_vmms:
            m_vmm = vmm_arr == vmm_id
            if not m_vmm.any():
                continue
            ch_vals  = ch_arr[m_vmm]
            is_on_v  = is_on[m_vmm]
            is_off_v = is_off[m_vmm]
            for i, ch in enumerate(chs_per_vmm[vmm_id]):
                m_ch = ch_vals == ch
                ch_on[vmm_id][i]  += (is_on_v  & m_ch).sum()
                ch_off[vmm_id][i] += (is_off_v & m_ch).sum()
        del df

    return {vmm_id: {"channels": chs_per_vmm[vmm_id],
                     "on":  ch_on[vmm_id],
                     "off": ch_off[vmm_id]}
            for vmm_id in detector_vmms}


def compute_spill_rates_all_runs(df_run_scan, data_dir,
                                  sng0_runs, trigger_ref_channels,
                                  detector_vmms,
                                  bin_width_s=0.001,
                                  n_files=99,
                                  spill_threshold_khz=1.0,
                                  min_spill_s=1.0,
                                  max_gap_s=0.0,
                                  good_channels=None,
                                  good_channels_per_config=None,
                                  normalize_per_channel=True,
                                  file_start=0,
                                  verbose=True):
    """
    For each run compute the mean hit rate during spill-on and spill-off
    periods for every detector VMM.

    Processes files one at a time to avoid loading a full run into memory.
    All VMMs share the same time bins so the trigger-derived spill mask
    applies consistently across detectors.

    Parameters
    ----------
    n_files : int
        Maximum number of ROOT files to read per run. Default 99
        loads all available files.
    file_start : int
        Index of the first file to read per run (0-based). Use 1 to
        skip file 0 when it has corrupted timestamps.
    spill_threshold_khz : float
        Trigger rate threshold (kHz) separating spill-on from spill-off.
    min_spill_s : float
        Minimum duration (seconds) of a contiguous above-threshold region
        to be counted as a real spill. Shorter bursts (e.g. startup
        artifacts at t≈0) are reclassified as spill-off. Default 1.0 s.
    max_gap_s : float
        Maximum gap (seconds) between two on-regions that will be filled
        in (hole-filling). Handles bins that dip below threshold mid-spill
        due to SPS bunch microstructure at low beam intensity. 0 disables.
    good_channels : dict {vmm_id: array-like} or None
        Global fallback good-channel mask (from build_good_channels()).
        Used when good_channels_per_config is None or has no entry for
        a given (sg, snt). If both are None, falls back to all connected
        channels.
    good_channels_per_config : dict {(sg, snt): good_channels_dict} or None
        Per-configuration good-channel masks built by
        build_good_channels_per_config(). When provided, each run uses
        the mask matching its own (sg, snt) so a channel noisy at one
        gain setting is not excluded at others.
    normalize_per_channel : bool
        If True (default), divide the VMM rate by the number of good
        channels so rates are expressed per channel — making VMMs with
        different numbers of connected pads directly comparable.

    Returns
    -------
    DataFrame with columns:
        run_no, sg, snt, vmm_id,
        rate_on_khz, rate_off_khz,         (per good channel if normalize_per_channel)
        rate_on_std_khz, rate_off_std_khz, (std across good channels)
        trig_rate_on_khz, n_on_bins, n_off_bins, n_good_channels
    """
    trig_vmm = list(trigger_ref_channels.keys())[0]
    trig_ch  = trigger_ref_channels[trig_vmm]

    records = []
    n_total = len(sng0_runs)

    for _i, run_no in enumerate(sng0_runs):
        row = df_run_scan[df_run_scan["run_no"] == run_no]
        if row.empty:
            continue
        sg  = row["sg"].iloc[0]
        snt = row["snt"].iloc[0]

        # Resolve per-config mask; fall back to global good_channels
        gc = (good_channels_per_config.get((sg, snt), good_channels)
              if good_channels_per_config is not None else good_channels)

        run_dir = get_run_dir(data_dir, run_no)
        n_avail = len(list_root_files(run_dir))
        if n_avail == 0:
            if not verbose:
                print()  # end progress-bar line before warning
            print(f"  Run {run_no}: no files found, skipping")
            continue
        n_load = min(n_files, n_avail)
        if verbose:
            print(f"  Run {run_no} (sg={sg}, snt={snt:.0f}): "
                  f"loading {n_load}/{n_avail} files...")
        else:
            filled = int(20 * (_i + 1) / max(n_total, 1))
            bar    = "=" * filled + "-" * (20 - filled)
            print(f"\r  [{bar}] {_i+1}/{n_total}  run {run_no}",
                  end="", flush=True)

        # Pass 1: get time range across all files (time column only)
        t_start, t_end = _get_run_time_range(run_dir, n_load,
                                              file_start=file_start)
        if t_start is None:
            print(f"    no valid timestamps, skipping")
            continue

        bins = np.arange(t_start, t_end + bin_width_s, bin_width_s)

        # Pass 2: accumulate histogram counts, one file at a time
        trig_counts, det_counts = _accumulate_run_histograms(
            run_dir, trig_vmm, trig_ch,
            detector_vmms, bins, n_load,
            good_channels=gc,
            file_start=file_start
        )

        trig_rate_khz  = trig_counts / bin_width_s / 1e3
        min_spill_bins = max(1, int(min_spill_s / bin_width_s))
        max_gap_bins   = max(0, int(max_gap_s   / bin_width_s))
        on_mask, off_mask = compute_spill_masks(
            trig_rate_khz, spill_threshold_khz, min_spill_bins, max_gap_bins
        )
        n_on  = on_mask.sum()
        n_off = off_mask.sum()

        if n_on == 0:
            print(f"    no spill-on bins "
                  f"(max trig rate={trig_rate_khz.max():.2f} kHz)")
            continue

        trig_rate_on = trig_rate_khz[on_mask].mean()

        # Per-channel counts for std-across-channels computation
        t_on_s  = n_on  * bin_width_s
        t_off_s = n_off * bin_width_s
        per_ch_counts = _accumulate_per_channel_on_off(
            run_dir, detector_vmms, gc,
            bins, on_mask, off_mask, n_load, file_start
        )

        for vmm_id in detector_vmms:
            det_rate_khz = det_counts[vmm_id] / bin_width_s / 1e3
            rate_on  = det_rate_khz[on_mask].mean()
            rate_off = (det_rate_khz[off_mask].mean()
                        if n_off > 0 else np.nan)

            if gc is not None and vmm_id in gc:
                n_good = max(1, len(gc[vmm_id]))
            else:
                from vmm_io import get_connected_channels
                chs = get_connected_channels(vmm_id)
                n_good = 64 if chs is None else len(chs)

            if normalize_per_channel:
                rate_on  = rate_on  / n_good
                rate_off = rate_off / n_good

            # Std across good channels
            pc = per_ch_counts[vmm_id]
            r_on_ch  = pc["on"]  / t_on_s  / 1e3 if t_on_s  > 0 else np.zeros(len(pc["on"]))
            r_off_ch = pc["off"] / t_off_s / 1e3 if t_off_s > 0 else np.zeros(len(pc["off"]))
            std_on  = float(r_on_ch.std())  if len(r_on_ch)  > 1 else 0.0
            std_off = float(r_off_ch.std()) if len(r_off_ch) > 1 else 0.0

            records.append({
                "run_no"            : run_no,
                "sg"                : sg,
                "snt"               : snt,
                "vmm_id"            : vmm_id,
                "rate_on_khz"       : rate_on,
                "rate_off_khz"      : rate_off,
                "rate_on_std_khz"   : std_on,
                "rate_off_std_khz"  : std_off,
                "trig_rate_on_khz"  : trig_rate_on,
                "n_on_bins"         : n_on,
                "n_off_bins"        : n_off,
                "n_good_channels"   : n_good,
            })

        if verbose:
            print(f"    on={n_on} bins  off={n_off} bins  "
                  f"trig_on={trig_rate_on:.1f} kHz")

    if not verbose and n_total > 0:
        print()  # end progress-bar line
    return pd.DataFrame(records)


def plot_spill_rates_vs_config(df_spill, vmm_groups,
                                trigger_ref_vmm=0,
                                out_dir=None, show=True, rate_tag=""):
    """
    Summary plot: mean VMM hit rate during spill-on and spill-off
    vs configuration (sg, snt).

    Layout: one figure per detector group.
      - Top panel   : trigger rate during spill-on — confirms beam
                      reference is stable across configurations.
      - One panel per VMM : spill-on (signal + noise) and spill-off
                      (noise floor) as separate lines.

    The x-axis lists (sg, snt) configurations sorted by gain then
    peaking time. Excess noise shows as elevated spill-off rate;
    signal rate should scale with gain and peak with optimal snt.

    Parameters
    ----------
    df_spill : DataFrame
        Output of compute_spill_rates_all_runs.
    vmm_groups : list of dict
        VMM groupings from vmm_mapping (detector groups only).
    trigger_ref_vmm : int
        VMM id of the trigger reference, shown in the top panel title.
    """
    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    # Configurations sorted by gain then peaking time
    configs = (
        df_spill[["sg", "snt"]]
        .drop_duplicates()
        .sort_values(["sg", "snt"])
        .reset_index(drop=True)
    )
    config_labels = [
        f"sg={r['sg']:.1f}\nsnt={r['snt']:.0f}"
        for _, r in configs.iterrows()
    ]
    x = np.arange(len(configs))

    for group in vmm_groups:
        vmm_ids     = group["vmm_ids"]
        group_label = group["label"]

        n_panels = 1 + len(vmm_ids)  # trigger panel + one per VMM
        fig, axes = plt.subplots(
            n_panels, 1,
            figsize=(max(10, 1.8 * len(configs)), 3.5 * n_panels),
            sharex=True
        )
        if n_panels == 1:
            axes = [axes]

        # ── Trigger rate panel ──────────────────────────────────
        trig_on = []
        for _, cfg_row in configs.iterrows():
            sub = df_spill[
                (df_spill["sg"]  == cfg_row["sg"]) &
                (df_spill["snt"] == cfg_row["snt"])
            ]
            trig_on.append(sub["trig_rate_on_khz"].mean()
                           if not sub.empty else np.nan)

        axes[0].plot(x, trig_on, "o-",
                     color="darkorange", linewidth=1.5,
                     label="Trigger (spill-on)")
        axes[0].set_ylabel("Rate (kHz)")
        axes[0].set_title(
            f"Trigger reference — VMM {trigger_ref_vmm}"
        )
        axes[0].legend(loc="upper right")
        axes[0].grid(True, alpha=0.3)
        axes[0].set_ylim(bottom=0)

        # ── Detector VMM panels ─────────────────────────────────
        for ax, vmm_id in zip(axes[1:], vmm_ids):
            det_name = vmm_to_detector.get(vmm_id, "")
            sub_vmm  = df_spill[df_spill["vmm_id"] == vmm_id]

            rate_on  = []
            rate_off = []
            for _, cfg_row in configs.iterrows():
                sub = sub_vmm[
                    (sub_vmm["sg"]  == cfg_row["sg"]) &
                    (sub_vmm["snt"] == cfg_row["snt"])
                ]
                if sub.empty:
                    rate_on.append(np.nan)
                    rate_off.append(np.nan)
                else:
                    rate_on.append(sub["rate_on_khz"].mean())
                    rate_off.append(sub["rate_off_khz"].mean())

            rate_on  = np.array(rate_on)
            rate_off = np.array(rate_off)

            ax.plot(x, rate_on,  "o-",  color="steelblue",
                    linewidth=1.5,
                    label="Spill-on  (signal + noise)")
            ax.plot(x, rate_off, "s--", color="firebrick",
                    linewidth=1.5,
                    label="Spill-off (noise floor)")
            ax.set_ylabel("Rate (kHz)")
            ax.set_title(f"VMM {vmm_id} — {det_name}")
            ax.legend(loc="upper right")
            ax.grid(True, alpha=0.3)
            ax.set_ylim(bottom=0)

        axes[-1].set_xticks(x)
        axes[-1].set_xticklabels(config_labels, fontsize=10)
        axes[-1].set_xlabel("Configuration  (sg, snt)")

        fig.suptitle(
            f"{group_label} — spill-on vs spill-off rate per configuration",
            fontweight="bold"
        )
        plt.tight_layout()
        stem = f"spill_rates_config_{group_label.replace(' ', '_')}"
        _finish_fig(fig, stem, out_dir, show, rate_tag)


def plot_spill_on_all_vmms(df_spill, vmm_groups, trigger_ref_vmm=0,
                            rate_col="rate_on_khz",
                            out_dir=None, show=True, rate_tag=""):
    """
    One figure per detector group with all VMMs on the same axes.

    rate_col : "rate_on_khz"  → spill-on summary
               "rate_off_khz" → spill-off (noise floor) summary
    Trigger reference (spill-on) is shown on a secondary axis.
    For rate_off plots the trigger line serves as a beam-stability
    reference rather than a noise comparison.
    """
    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    is_off    = rate_col == "rate_off_khz"
    tag       = "spill-off (noise floor)" if is_off else "spill-on"
    std_col   = rate_col.replace("_khz", "_std_khz")
    ylabel    = f"Avg. rate per channel — {tag} (kHz)"

    configs = (
        df_spill[["sg", "snt"]]
        .drop_duplicates()
        .sort_values(["sg", "snt"])
        .reset_index(drop=True)
    )
    config_labels = [
        f"sg={r['sg']:.1f}\nsnt={r['snt']:.0f}"
        for _, r in configs.iterrows()
    ]
    x = np.arange(len(configs))

    for group in vmm_groups:
        vmm_ids     = group["vmm_ids"]
        group_label = group["label"]

        fig, ax = plt.subplots(figsize=(max(10, 1.8 * len(configs)), 5))

        colors = plt.cm.tab10(np.linspace(0, 1, len(vmm_ids)))

        for color, vmm_id in zip(colors, vmm_ids):
            det_name = vmm_to_detector.get(vmm_id, "")
            sub_vmm  = df_spill[df_spill["vmm_id"] == vmm_id]
            has_std  = std_col in df_spill.columns

            rates = []
            stds  = []
            for _, cfg_row in configs.iterrows():
                sub = sub_vmm[
                    (sub_vmm["sg"]  == cfg_row["sg"]) &
                    (sub_vmm["snt"] == cfg_row["snt"])
                ]
                rates.append(sub[rate_col].mean() if not sub.empty else np.nan)
                stds.append(sub[std_col].mean()   if (has_std and not sub.empty) else 0.0)

            rates = np.array(rates)
            stds  = np.array(stds)

            ax.plot(x, rates, "o-", color=color, linewidth=1.5,
                    label=f"VMM {vmm_id} — {det_name}")
            if has_std:
                ax.fill_between(x, np.maximum(0, rates - stds), rates + stds,
                                alpha=0.15, color=color, linewidth=0)

        # Trigger reference on secondary y-axis
        trig_on = []
        for _, cfg_row in configs.iterrows():
            sub = df_spill[
                (df_spill["sg"]  == cfg_row["sg"]) &
                (df_spill["snt"] == cfg_row["snt"])
            ]
            trig_on.append(sub["trig_rate_on_khz"].mean()
                           if not sub.empty else np.nan)

        ax_r = ax.twinx()
        ax_r.plot(x, trig_on, "k--", linewidth=1.5,
                  alpha=0.5, label=f"Trigger VMM {trigger_ref_vmm}")
        ax_r.set_ylabel("Trigger rate (kHz)", color="black", alpha=0.6)
        ax_r.tick_params(axis="y", labelcolor="gray")
        ax_r.set_ylim(bottom=0)

        lines_l, labels_l = ax.get_legend_handles_labels()
        lines_r, labels_r = ax_r.get_legend_handles_labels()
        ax.legend(lines_l + lines_r, labels_l + labels_r,
                  fontsize=9, loc="upper left")

        ax.set_xticks(x)
        ax.set_xticklabels(config_labels, fontsize=10)
        ax.set_xlabel("Configuration  (sg, snt)")
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)

        fig.suptitle(
            f"{group_label} — avg. {tag} rate per channel vs configuration  "
            f"(shaded band = ±1σ across channels)",
            fontweight="bold"
        )
        plt.tight_layout()
        stem = (f"spill_on_vmms_{group_label.replace(' ', '_')}"
                f"_{rate_col}")
        _finish_fig(fig, stem, out_dir, show, rate_tag)


if __name__ == '__main__':
    main()
