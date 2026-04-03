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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from vmm_mapping import vmm_mapping
from vmm_io import (load_run_table, get_run_groups,
                    get_run_dir, get_root_file,
                    load_hits_root)
# ── Constants ──────────────────────────────────────────────
NS_PER_TICK  = 1.0          # 1 GHz clock
S_PER_TICK   = NS_PER_TICK * 1e-9

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


def main():
    # ── User configuration ──────────────────────────────────
    cnfg_dir = "/drf/projets/clas12/P2/akallits/"
    data_dir = "/drf/projets/clas12/cern_202511_p2_alinx/"

    # cnfg_dir = "/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/VMM-alinx-data/"
    # data_dir = "/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/VMM-alinx-data/5kHz-muons-config-scan/"


    root_file_index = 1

    # ── Load run metadata ───────────────────────────────────
    df_run_scan = load_run_table(f"{cnfg_dir}vmm_config_scan.csv")
    run_groups = get_run_groups(df_run_scan)

    for run_no in run_groups["sng0_runs"]:
        run_dir = get_run_dir(data_dir, run_no)
        file_path = get_root_file(run_dir, file_index=root_file_index)
        if file_path is None:
            continue

        df = load_hits_root(file_path, branches=["time", "vmm"])
        df = df.sort_values("time").reset_index(drop=True)
        t = df["time"].values

        sg = df_run_scan.loc[df_run_scan["run_no"] == run_no, "sg"].iloc[0]
        snt = df_run_scan.loc[df_run_scan["run_no"] == run_no, "snt"].iloc[0]

        print(f"\nRun {run_no} (sg={sg} snt={snt}):")
        print(f"  N hits          : {len(t)}")
        print(f"  t[0]            : {t[0]:.3e}")
        print(f"  t[-1]           : {t[-1]:.3e}")
        print(f"  t range (ticks) : {t[-1] - t[0]:.3e}")
        print(f"  t range (s)     : {(t[-1] - t[0]) * S_PER_TICK:.1f}")
        print(f"  t[0] (s)        : {t[0] * S_PER_TICK:.1f}")
        print(f"  t[-1] (s)       : {t[-1] * S_PER_TICK:.1f}")

        # Check for large jumps
        dt = np.diff(t)
        large_jumps = (dt > 1e11)
        if large_jumps.any():
            jump_idx = np.where(large_jumps)[0]
            print(f"  WARNING: {large_jumps.sum()} large jumps detected")
            for idx in jump_idx[:3]:
                print(f"    Jump at index {idx}: "
                      f"Δt={dt[idx]:.3e} ticks "
                      f"= {dt[idx] * S_PER_TICK:.1f}s")


    # input("Press enter to continue...")
    # ── Step 1: validate on one run ─────────────────────────
    # Use run 67 (sg=3, snt=200, sng=1) as the test case
    test_run = 67

    print(f"Loading run {test_run}...")
    df_hits = load_sorted_hits(
        data_dir, test_run,
        branches=["time", "vmm", "ch"],
        root_file_index=root_file_index
    )

    # Inspect each channel separately
    for vmm_id in trigger_vmms:
        inspect_inter_event_per_channel(df_hits, vmm_id,
                                        run_no=67)

    duration = compute_run_duration(df_hits)
    print(f"Run duration : {duration:.1f} s")

    # Hit rates for all VMMs
    print(f"\n{'VMM':>5} {'Type':>10} {'N_hits':>10} "
          f"{'Rate (Hz)':>12}")
    print("-" * 42)

    all_vmms = trigger_vmms + detector_vmms
    for vmm_id in sorted(all_vmms):
        n = (df_hits["vmm"] == vmm_id).sum()
        rate = compute_hit_rate(df_hits, vmm_id)
        vtype = "trigger" if vmm_id in trigger_vmms else "detector"
        print(f"{vmm_id:>5} {vtype:>10} {n:>10} {rate:>12.1f}")

    # For trigger rate — use one channel per physical trigger signal
    # VMM 0 ch48 and VMM 1 ch20 are the reference channels
    trigger_reference = {
        0: 48,  # VMM 0, channel 48
        1: 20,  # VMM 1, channel 0
    }

    print(f"\nFitting inter-event distributions...")
    for vmm_id, ref_ch in trigger_reference.items():
        # Get hits for this specific channel only
        df_ch = df_hits[
            (df_hits["vmm"] == vmm_id) &
            (df_hits["ch"] == ref_ch)
            ]
        t = df_ch["time"].values
        dt_ns = np.diff(t) * NS_PER_TICK
        dt_us = dt_ns[(dt_ns > 0)] / 1000.0
        dt_us = dt_us[dt_us < 1000]

        print(f"\nVMM {vmm_id} ch {ref_ch}: {len(dt_us)} Δt values")

        rate_hz, rate_err_hz, chi2_ndf, tau, fit_x, fit_y = \
            fit_exponential(dt_us)

        print(f"  Offset τ        : {tau:.1f} μs")
        print(f"  Fitted rate     : {rate_hz:.1f} ± {rate_err_hz:.1f} Hz")
        print(f"  Simple rate     : {len(t) / (t[-1] - t[0]) / S_PER_TICK:.1f} Hz")
        print(f"  χ²/ndf          : {chi2_ndf:.2f}")

        plot_inter_event_distribution(
            dt_us, f"{vmm_id} ch{ref_ch}", test_run,
            rate_hz, chi2_ndf, tau, fit_x, fit_y
        )

    # Detector VMMs — all channels combined per channel
    print(f"\nDetector VMM inter-event times...")
    for vmm_id in detector_vmms:
        dt_us = compute_inter_event_times(
            df_hits, vmm_id, max_dt_us=1000
        )
        if len(dt_us) < 50:
            continue

        print(f"\nVMM {vmm_id}: {len(dt_us)} Δt values")
        rate_hz, rate_err_hz, chi2_ndf, tau, fit_x, fit_y = \
            fit_exponential(dt_us)

        print(f"  Offset τ        : {tau:.1f} μs")
        print(f"  Fitted rate     : {rate_hz:.1f} ± {rate_err_hz:.1f} Hz")
        print(f"  Simple rate     : "
              f"{compute_hit_rate(df_hits, vmm_id):.1f} Hz")
        print(f"  χ²/ndf          : {chi2_ndf:.2f}")

        plot_inter_event_distribution(
            dt_us, vmm_id, test_run,
            rate_hz, chi2_ndf, tau, fit_x, fit_y
        )
    trigger_ref_channels = {0: 48, 1: 20}

    df_rates = compute_rates_all_runs(
        df_run_scan=df_run_scan,
        data_dir=data_dir,
        sng0_runs=run_groups["sng0_runs"],
        trigger_vmms=trigger_vmms,
        detector_vmms=detector_vmms,
        trigger_ref_channels=trigger_ref_channels,
        root_file_index=root_file_index
    )

    df_rates.to_csv("vmm_rates.csv", index=False)
    print("\n")
    print(df_rates.to_string(index=False))

    df_eff = compute_efficiency(df_rates)
    df_eff.to_csv("vmm_efficiency.csv", index=False)

    print("\n=== Efficiency summary ===")
    print(df_eff[["run_no", "sg", "snt", "vmm_id",
                  "det_rate", "trig_rate",
                  "efficiency"]].to_string(index=False))

    plot_rates_vs_config(df_rates, df_eff,
                         detector_vmms, vmm_groups)

    # ── Trigger stream QA on test run ──────────────────────────
    df_hits_dt = load_sorted_hits(
        data_dir, test_run,
        branches=["time", "vmm", "ch"],
        root_file_index=root_file_index
    )

    inspect_trigger_timing(
        df_hits_dt,
        trigger_vmm=0,
        trigger_ch=trigger_ref_channels[0],
        window_ns=600
    )

    # ── Triggers/ms vs time — spill structure ──────────────────
    plot_trigger_rate_ms(
        df_hits_dt,
        run_no=test_run,
        trigger_vmm=0,
        trigger_ch=trigger_ref_channels[0],
        bin_width_ms=1.0,
    )

    # ── Detector rate overlaid with trigger rate ────────────────
    plot_rate_overlay_with_trigger(
        df_hits_dt,
        run_no=test_run,
        vmm_groups=vmm_groups,
        trigger_vmm=0,
        trigger_ch=trigger_ref_channels[0],
        bin_width_ms=1.0,
    )

    # ── Rate vs time — all runs, 1 ms bins ─────────────────────
    plot_rate_vs_time(
        data_dir=data_dir,
        df_run_scan=df_run_scan,
        sng0_runs=run_groups["sng0_runs"],
        vmm_groups=vmm_groups,
        trigger_ref_channels=trigger_ref_channels,
        bin_width_s=0.001,
        root_file_index=root_file_index,
    )

    print("bonzo")

def load_sorted_hits(data_dir, run_no, branches,
                      root_file_index=1,
                      connected_channels_only=True,
                      max_time_ticks=2e12):
    from vmm_io import get_connected_channels
    from functools import reduce
    import operator

    run_dir   = get_run_dir(data_dir, run_no)
    file_path = get_root_file(run_dir,
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

    df = df.sort_values("time").reset_index(drop=True)
    return df

def inspect_trigger_channels(data_dir, run_no,
                              root_file_index=1):
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
        # plt.show()

def inspect_inter_event_per_channel(df_hits, vmm_id, run_no,
                                     max_dt_us=1000):
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
    plt.show()

def compute_run_duration(df_hits):
    """
    Compute run duration in seconds from timestamp range.
    """
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

def compute_efficiency(df_rates, trigger_ref_vmm=0):
    """
    Compute detector hit rate / trigger rate per run per VMM.
    Uses VMM 0 ch48 as the reference trigger.

    Returns
    -------
    pd.DataFrame
        One row per (run, detector VMM) with efficiency.
    """
    results = []

    for run_no in df_rates["run_no"].unique():
        df_run = df_rates[df_rates["run_no"] == run_no]

        # Get trigger rate for this run
        trig_row = df_run[
            (df_run["vmm_id"] == trigger_ref_vmm) &
            (df_run["type"]   == "trigger")
        ]
        if trig_row.empty:
            continue
        trig_rate = trig_row["rate_hz"].iloc[0]
        if trig_rate <= 0:
            continue

        sg  = df_run["sg"].iloc[0]
        snt = df_run["snt"].iloc[0]

        # Compute efficiency per detector VMM
        det_rows = df_run[df_run["type"] == "detector"]
        for _, row in det_rows.iterrows():
            efficiency = row["rate_hz"] / trig_rate

            results.append({
                "run_no"    : run_no,
                "sg"        : sg,
                "snt"       : snt,
                "vmm_id"    : row["vmm_id"],
                "det_rate"  : row["rate_hz"],
                "trig_rate" : trig_rate,
                "efficiency": efficiency
            })

    return pd.DataFrame(results)


def plot_rates_vs_config(df_rates, df_eff, detector_vmms,
                          vmm_groups):
    """
    Two plots:
    1. Trigger rate vs run (shows beam intensity variation)
    2. Detector efficiency (det_rate/trig_rate) vs configuration
       per VMM grouped by detector
    """
    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    # --- Plot 1: trigger rate per run ---
    trig = df_rates[
        (df_rates["type"]   == "trigger") &
        (df_rates["vmm_id"] == 0)
    ].sort_values("run_no")

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(trig["run_no"].astype(str),
           trig["rate_hz"],
           color="steelblue", alpha=0.7)
    for _, row in trig.iterrows():
        ax.text(str(row["run_no"]),
                row["rate_hz"] + 5,
                f"sg={row['sg']:.1f}\nsnt={row['snt']:.0f}",
                ha="center", va="bottom", fontsize=7)
    ax.set_xlabel("Run number")
    ax.set_ylabel("Trigger rate (Hz)")
    ax.set_title("Trigger rate per run (VMM 0 ch48)")
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    plt.show()

    # --- Plot 2: efficiency vs peaking time per VMM ---
    # One panel per detector group
    for group in vmm_groups:
        vmm_ids = group["vmm_ids"]
        label   = group["label"]

        fig, ax = plt.subplots(figsize=(10, 5))

        for vmm_id in vmm_ids:
            df_vmm = df_eff[
                df_eff["vmm_id"] == vmm_id
            ].sort_values(["sg", "snt"])

            if df_vmm.empty:
                continue

            # Config label for x axis
            x_labels = [
                f"sg={r['sg']:.1f}\nsnt={r['snt']:.0f}"
                for _, r in df_vmm.iterrows()
            ]
            x_pos = range(len(x_labels))

            ax.plot(x_pos, df_vmm["efficiency"],
                    "o-", label=f"VMM {vmm_id}",
                    linewidth=1.5, markersize=6)

        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels, fontsize=8)
        ax.set_xlabel("Configuration")
        ax.set_ylabel("Efficiency (detector rate / trigger rate)")
        ax.set_title(f"Hit efficiency per configuration — {label}")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()


def plot_inter_event_distribution(dt_us, vmm_id, run_no,
                                    rate_hz, chi2_ndf, tau,
                                    fit_x, fit_y,
                                    n_bins=200):
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
    plt.show()


def inspect_trigger_timing(df_hits, trigger_vmm=0, trigger_ch=48,
                            window_ns=600):
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
    plt.show()

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
                       root_file_index=1):
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
        plt.show()


def plot_trigger_rate_ms(df_hits, run_no,
                          trigger_vmm=0, trigger_ch=48,
                          bin_width_ms=1.0):
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
    plt.show()


def plot_rate_overlay_with_trigger(df_hits, run_no, vmm_groups,
                                    trigger_vmm=0, trigger_ch=48,
                                    bin_width_ms=1.0):
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
        plt.show()


if __name__ == '__main__':
    main()
