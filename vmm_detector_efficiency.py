#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vmm_detector_efficiency.py

Detector performance analysis using trigger–detector time coincidences.

Answers:
  - What is the timing offset between trigger and detector signals?
  - What is the true detector efficiency per VMM per configuration,
    corrected for accidental coincidences?

Requires:
  - sng=0 runs (signal + trigger hits present)
  - trigger reference channel defined in vmm_mapping
  - load_sorted_hits from vmm_trigger_analysis

@author: ak271430
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from vmm_mapping import vmm_mapping
from vmm_io import load_run_table, get_run_groups
from vmm_trigger_analysis import (load_sorted_hits,
                                   NS_PER_TICK,
                                   detector_vmms,
                                   vmm_groups,
                                   trigger_ref_channels)


# ─────────────────────────────────────────────────────────────
# TIMING DIAGNOSTIC
# ─────────────────────────────────────────────────────────────

def compute_trigger_detector_dt(df_hits,
                                  trigger_vmm=0,
                                  trigger_ch=48,
                                  search_window_ns=1000):
    """
    For each trigger hit, find all detector hits within
    ±search_window_ns and record Δt = t_det - t_trig (in ns).

    Uses np.searchsorted on a sorted detector time array —
    efficient for large datasets (O(N_trig × log N_det)).

    Parameters
    ----------
    search_window_ns : float
        Half-width of the search window in ns (= ticks at 1 GHz).
        Default 1000 ns = ±1 μs, wide enough to see any offset.

    Returns
    -------
    dict : vmm_id (int) → np.ndarray of Δt values in ns
    """
    t_trig = df_hits[
        (df_hits["vmm"] == trigger_vmm) &
        (df_hits["ch"]  == trigger_ch)
    ]["time"].values

    if len(t_trig) == 0:
        print(f"  No hits on trigger VMM {trigger_vmm} ch {trigger_ch}")
        return {}

    # Deduplicate trigger timestamps — identical timestamps from
    # the same physical pulse would inflate the denominator.
    n_before = len(t_trig)
    t_trig   = np.unique(t_trig)
    n_dedup  = n_before - len(t_trig)
    if n_dedup > 0:
        print(f"  Removed {n_dedup} duplicate trigger timestamps")

    window = search_window_ns  # 1 tick = 1 ns

    result = {}
    for vmm_id in detector_vmms:
        t_det = df_hits[df_hits["vmm"] == vmm_id]["time"].values
        if len(t_det) == 0:
            continue

        t_det_sorted = np.sort(t_det)
        all_dt = []

        for t_t in t_trig:
            lo = np.searchsorted(t_det_sorted, t_t - window)
            hi = np.searchsorted(t_det_sorted, t_t + window)
            if lo < hi:
                all_dt.append(
                    (t_det_sorted[lo:hi] - t_t) * NS_PER_TICK
                )

        result[vmm_id] = (np.concatenate(all_dt)
                          if all_dt else np.array([]))

    return result


def plot_trigger_detector_dt(dt_dict, run_no, vmm_groups,
                               search_window_ns=1000,
                               n_bins=200):
    """
    Plot Δt = t_detector - t_trigger distributions.

    One figure per detector group, one subplot per VMM.
    Each subplot shows linear (top) and log (bottom) scale
    so both the narrow coincidence peak and the flat
    accidental background are visible.

    A vertical dashed line marks Δt = 0.
    The peak position and width guide the choice of
    coincidence window in compute_time_correlated_efficiency().
    """
    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    for group in vmm_groups:
        vmm_ids     = [v for v in group["vmm_ids"] if v in dt_dict]
        group_label = group["label"]

        if not vmm_ids:
            continue

        n_vmm = len(vmm_ids)
        fig, axes = plt.subplots(
            2, n_vmm,
            figsize=(5 * n_vmm, 8),
            sharex=True
        )
        if n_vmm == 1:
            axes = np.array(axes).reshape(2, 1)

        for col, vmm_id in enumerate(vmm_ids):
            dt       = dt_dict[vmm_id]
            det_name = vmm_to_detector.get(vmm_id, "")
            title    = (f"VMM {vmm_id} — {det_name}\n"
                        f"n_pairs={len(dt):,}")

            for row_idx, (ax, yscale) in enumerate(
                zip(axes[:, col], ["linear", "log"])
            ):
                ax.hist(dt, bins=n_bins,
                        range=(-search_window_ns, search_window_ns),
                        color="steelblue", alpha=0.7)
                ax.axvline(0, color="red", linestyle="--",
                           linewidth=1.2,
                           label="Δt = 0" if row_idx == 0 else None)
                ax.set_yscale(yscale)
                ax.set_ylabel("Counts" if yscale == "linear"
                               else "Counts (log)")
                if row_idx == 0:
                    ax.set_title(title, fontsize=11)
                    ax.legend(fontsize=9)
                else:
                    ax.set_xlabel("Δt (ns)")

        fig.suptitle(
            f"{group_label} — trigger–detector Δt  |  "
            f"Run {run_no}  |  window ±{search_window_ns} ns",
            fontsize=13, fontweight="bold"
        )
        plt.tight_layout()
        plt.show()


# ─────────────────────────────────────────────────────────────
# COINCIDENCE EFFICIENCY
# ─────────────────────────────────────────────────────────────

def compute_time_correlated_efficiency(df_hits,
                                        trigger_vmm=0,
                                        trigger_ch=48,
                                        dt_min_ns=0,
                                        dt_max_ns=250,
                                        sideband_min_ns=400,
                                        sideband_max_ns=650):
    """
    Compute true coincidence efficiency per detector VMM.

    For each trigger hit, checks whether at least one detector hit
    falls in the signal window [dt_min_ns, dt_max_ns].
    Uses a sideband window of equal width to estimate and subtract
    the accidental coincidence rate.

    Definition
    ----------
    raw_efficiency  = n_triggers_with_≥1_signal_hit / n_triggers
    accidental_eff  = n_triggers_with_≥1_sideband_hit / n_triggers
    true_efficiency = raw_efficiency - accidental_eff

    Parameters
    ----------
    dt_min_ns, dt_max_ns : float
        Signal coincidence window in ns (Δt = t_det - t_trig).
    sideband_min_ns, sideband_max_ns : float
        Sideband window for accidental estimation — same width
        as signal window, in the flat background region.

    Returns
    -------
    pd.DataFrame
        One row per detector VMM with columns:
        vmm_id, n_triggers, n_matched, n_accidental,
        raw_efficiency, accidental_eff, true_efficiency
    """
    t_trig = df_hits[
        (df_hits["vmm"] == trigger_vmm) &
        (df_hits["ch"]  == trigger_ch)
    ]["time"].values

    if len(t_trig) == 0:
        print(f"No hits on trigger VMM {trigger_vmm} ch {trigger_ch}")
        return pd.DataFrame()

    t_trig     = np.sort(np.unique(t_trig))
    n_triggers = len(t_trig)

    sig_width  = dt_max_ns - dt_min_ns
    side_width = sideband_max_ns - sideband_min_ns
    if abs(sig_width - side_width) > 1:
        print(f"  WARNING: signal width ({sig_width} ns) != "
              f"sideband width ({side_width} ns) — "
              f"accidental subtraction will be scaled")
    scale = sig_width / side_width

    results = []

    for vmm_id in detector_vmms:
        t_det = df_hits[df_hits["vmm"] == vmm_id]["time"].values
        if len(t_det) == 0:
            continue

        t_det_sorted = np.sort(t_det)
        n_matched    = 0
        n_accidental = 0

        for t_t in t_trig:
            lo = np.searchsorted(t_det_sorted, t_t + dt_min_ns)
            hi = np.searchsorted(t_det_sorted, t_t + dt_max_ns)
            if hi > lo:
                n_matched += 1

            lo_s = np.searchsorted(t_det_sorted, t_t + sideband_min_ns)
            hi_s = np.searchsorted(t_det_sorted, t_t + sideband_max_ns)
            if hi_s > lo_s:
                n_accidental += 1

        raw_eff  = n_matched    / n_triggers
        acc_eff  = n_accidental / n_triggers * scale
        true_eff = raw_eff - acc_eff

        results.append({
            "vmm_id"         : vmm_id,
            "n_triggers"     : n_triggers,
            "n_matched"      : n_matched,
            "n_accidental"   : n_accidental,
            "raw_efficiency" : raw_eff,
            "accidental_eff" : acc_eff,
            "true_efficiency": true_eff,
        })

        print(f"  VMM {vmm_id}: "
              f"matched={n_matched}  "
              f"accidental={n_accidental}  "
              f"raw={raw_eff:.3f}  "
              f"acc={acc_eff:.4f}  "
              f"true={true_eff:.3f}")

    return pd.DataFrame(results)


def compute_time_correlated_efficiency_all_runs(
        df_run_scan, data_dir, sng0_runs,
        trigger_vmm=0, trigger_ch=48,
        dt_min_ns=0, dt_max_ns=250,
        sideband_min_ns=400, sideband_max_ns=650,
        root_file_index=1):
    """
    Run compute_time_correlated_efficiency() over all sng=0 runs.

    Returns
    -------
    pd.DataFrame
        One row per (run_no, vmm_id) with columns:
        run_no, sg, snt, vmm_id, n_triggers, n_matched,
        n_accidental, raw_efficiency, accidental_eff, true_efficiency
    """
    all_results = []

    for run_no in sng0_runs:
        sg  = df_run_scan.loc[
            df_run_scan["run_no"] == run_no, "sg"
        ].iloc[0]
        snt = df_run_scan.loc[
            df_run_scan["run_no"] == run_no, "snt"
        ].iloc[0]

        print(f"\nRun {run_no}  sg={sg}  snt={snt}")

        df_hits = load_sorted_hits(
            data_dir, run_no,
            branches=["time", "vmm", "ch"],
            root_file_index=root_file_index
        )
        if df_hits is None or df_hits.empty:
            print("  No data, skipping")
            continue

        df_eff = compute_time_correlated_efficiency(
            df_hits,
            trigger_vmm=trigger_vmm,
            trigger_ch=trigger_ch,
            dt_min_ns=dt_min_ns,
            dt_max_ns=dt_max_ns,
            sideband_min_ns=sideband_min_ns,
            sideband_max_ns=sideband_max_ns,
        )
        if df_eff.empty:
            continue

        df_eff["run_no"] = run_no
        df_eff["sg"]     = sg
        df_eff["snt"]    = snt
        all_results.append(df_eff)

    if not all_results:
        return pd.DataFrame()

    cols = ["run_no", "sg", "snt", "vmm_id",
            "n_triggers", "n_matched", "n_accidental",
            "raw_efficiency", "accidental_eff", "true_efficiency"]
    return pd.concat(all_results, ignore_index=True)[cols]


def plot_time_correlated_efficiency(df_coinc, vmm_groups):
    """
    Plot true coincidence efficiency vs configuration per VMM,
    grouped by detector.

    One figure per detector group, one line per VMM.
    X axis: configuration (sg / snt), sorted.
    Y axis: true_efficiency with binomial error bars.
    """
    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    configs = (
        df_coinc[["sg", "snt"]]
        .drop_duplicates()
        .sort_values(["sg", "snt"])
    )
    config_labels = [
        f"sg={r['sg']:.1f}\nsnt={r['snt']:.0f}"
        for _, r in configs.iterrows()
    ]
    config_keys = [
        (r["sg"], r["snt"])
        for _, r in configs.iterrows()
    ]

    for group in vmm_groups:
        vmm_ids     = group["vmm_ids"]
        group_label = group["label"]

        fig, ax = plt.subplots(
            figsize=(max(8, len(config_keys) * 1.5), 5)
        )

        for vmm_id in vmm_ids:
            df_vmm = df_coinc[df_coinc["vmm_id"] == vmm_id]
            if df_vmm.empty:
                continue

            x_pos, y, y_err = [], [], []
            for idx, (sg, snt) in enumerate(config_keys):
                row = df_vmm[
                    (df_vmm["sg"] == sg) & (df_vmm["snt"] == snt)
                ]
                if row.empty:
                    continue
                eff = row["true_efficiency"].iloc[0]
                n   = row["n_triggers"].iloc[0]
                x_pos.append(idx)
                y.append(eff)
                y_err.append(
                    np.sqrt(eff * (1 - eff) / n) if n > 0 else 0
                )

            det_name = vmm_to_detector.get(vmm_id, "")
            ax.errorbar(x_pos, y, yerr=y_err,
                        fmt="o-", linewidth=1.5, markersize=6,
                        label=f"VMM {vmm_id} — {det_name}",
                        capsize=3)

        ax.set_xticks(range(len(config_labels)))
        ax.set_xticklabels(config_labels)
        ax.set_xlabel("Configuration (sg / snt)")
        ax.set_ylabel("True coincidence efficiency")
        ax.set_title(
            f"{group_label}\n"
            f"Time-correlated efficiency  |  "
            f"signal [0, 250] ns  |  sideband [400, 650] ns"
        )
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)
        plt.tight_layout()
        plt.show()


# ─────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────

def main():
    # ── User configuration ──────────────────────────────────
    cnfg_dir = "/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/VMM-alinx-data/"
    data_dir = "/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/VMM-alinx-data/5kHz-muons-config-scan/"

    root_file_index = 1
    test_run        = 67       # sg=3, snt=200 — used for diagnostics

    trigger_vmm = 0
    trigger_ch  = trigger_ref_channels[trigger_vmm]

    # ── Load run metadata ───────────────────────────────────
    df_run_scan = load_run_table(f"{cnfg_dir}vmm_config_scan.csv")
    run_groups  = get_run_groups(df_run_scan)

    # ── Timing diagnostic on test run ──────────────────────
    print(f"\nLoading run {test_run} for timing diagnostic...")
    df_hits = load_sorted_hits(
        data_dir, test_run,
        branches=["time", "vmm", "ch"],
        root_file_index=root_file_index
    )

    print(f"\nComputing trigger–detector Δt on run {test_run}...")
    dt_dict = compute_trigger_detector_dt(
        df_hits,
        trigger_vmm=trigger_vmm,
        trigger_ch=trigger_ch,
        search_window_ns=1000
    )
    for vmm_id, dt in dt_dict.items():
        print(f"  VMM {vmm_id}: {len(dt):,} pairs in ±1000 ns window")

    plot_trigger_detector_dt(
        dt_dict, test_run, vmm_groups,
        search_window_ns=1000
    )

    # ── Single-run efficiency (verify window choice) ────────
    print(f"\nComputing time-correlated efficiency on run {test_run}...")
    df_eff_test = compute_time_correlated_efficiency(
        df_hits,
        trigger_vmm=trigger_vmm,
        trigger_ch=trigger_ch,
        dt_min_ns=0,
        dt_max_ns=250,
        sideband_min_ns=400,
        sideband_max_ns=650,
    )
    print("\n=== Time-correlated efficiency (test run) ===")
    print(df_eff_test[["vmm_id", "n_matched", "n_accidental",
                         "raw_efficiency", "accidental_eff",
                         "true_efficiency"]].to_string(index=False))

    # ── All runs ────────────────────────────────────────────
    print("\nComputing time-correlated efficiency across all runs...")
    df_coinc_all = compute_time_correlated_efficiency_all_runs(
        df_run_scan=df_run_scan,
        data_dir=data_dir,
        sng0_runs=run_groups["sng0_runs"],
        trigger_vmm=trigger_vmm,
        trigger_ch=trigger_ch,
        dt_min_ns=0,
        dt_max_ns=250,
        sideband_min_ns=400,
        sideband_max_ns=650,
        root_file_index=root_file_index,
    )

    df_coinc_all.to_csv(f"{cnfg_dir}vmm_coinc_efficiency.csv",
                         index=False)

    print("\n=== Time-correlated efficiency — all runs ===")
    print(df_coinc_all[["run_no", "sg", "snt", "vmm_id",
                          "n_triggers", "n_matched",
                          "true_efficiency"]].to_string(index=False))

    plot_time_correlated_efficiency(df_coinc_all, vmm_groups)


if __name__ == "__main__":
    main()
