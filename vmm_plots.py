#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/25/26 1:18 PM
Created in PyCharm
Created as vmm_plots.py

@author: ak271430
"""
"""
vmm_plots.py

Plotting functions for VMM config scan analysis.

Two categories:
- Legacy plots : ADC distributions, mean/std vs peaking time
                 (kept for backward compatibility)
- SNR plots    : configuration comparison at VMM and channel level
                 (main deliverables of the analysis)

All SNR comparison plots use a consistent color scale and
derive VMM groupings and detector names from vmm_mapping.py.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from vmm_mapping import vmm_mapping


# ─────────────────────────────────────────────
# LEGACY PLOTS
# ─────────────────────────────────────────────

def plot_adc_histograms(df_hits, run_no, ax):
    """Plot ADC histograms per VMM on a given axes."""
    for vmm_id in df_hits["vmm"].unique():
        adc_values = df_hits[df_hits["vmm"] == vmm_id]["adc"]
        ax.hist(adc_values, bins=50, alpha=0.4,
                label=f"Run {run_no} VMM {vmm_id}")


def plot_adc_histograms_for_runs(run_list, data_dir,
                                  get_run_dir_fn,
                                  list_root_files_fn,
                                  load_hits_root_fn):
    """Plot ADC distributions for all runs overlaid."""
    from vmm_io import get_run_dir, list_root_files, load_hits_root
    fig, ax = plt.subplots()

    for run_no in run_list:
        run_dir    = get_run_dir(data_dir, run_no)
        root_files = list_root_files(run_dir, n=2)
        if not root_files:
            continue
        file_path = os.path.join(run_dir, root_files[1])
        df_hits   = load_hits_root(
            file_path,
            branches=["adc", "vmm", "ch", "time", "over_threshold"]
        )
        for vmm_id in df_hits["vmm"].unique():
            adc_values = df_hits[df_hits["vmm"] == vmm_id]["adc"]
            ax.hist(adc_values, bins=50, alpha=0.4,
                    label=f"Run {run_no} VMM {vmm_id}")

    ax.set_xlabel("ADC Value")
    ax.set_ylabel("Counts")
    ax.set_title("ADC Value Distribution by Run/VMM")
    ax.legend()
    plt.show()


def compare_full_vs_cut(df_hits, vmm_id, run_no,
                         adc_cut=100, use_over_threshold=True):
    """Compare full ADC distribution vs selected noise events."""
    full_adc = df_hits[df_hits["vmm"] == vmm_id]["adc"].values

    if use_over_threshold:
        cut_adc   = df_hits[
            (df_hits["vmm"] == vmm_id) &
            (df_hits["over_threshold"] == 0)
        ]["adc"].values
        label_cut = "over_threshold == 0"
    else:
        cut_adc   = full_adc[full_adc < adc_cut]
        label_cut = f"ADC < {adc_cut}"

    plt.figure()
    plt.hist(full_adc, bins=80, alpha=0.3,
             label="Full", density=False)
    plt.hist(cut_adc, bins=80, alpha=0.6,
             label=label_cut, density=False)
    plt.title(f"Run {run_no} – VMM {vmm_id}")
    plt.xlabel("ADC")
    plt.ylabel("Numb of Hits")
    plt.legend()
    plt.show()


def plot_removed_fraction(df_results, vmm_ids):
    """Plot fraction of hits removed by ADC cut vs peaking time."""
    plt.figure()
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.plot(df_vmm["peaking_time"], df_vmm["frac_removed"],
                 "o-", label=f"VMM {vmm_id}")
    plt.xlabel("Peaking time (snt)")
    plt.ylabel("Fraction removed by ADC over threshold cut")
    plt.title("Truncation bias vs peaking time")
    plt.legend()
    plt.show()


def plot_mean_vs_peaking(df_results, vmm_ids):
    """Plot mean ADC vs peaking time per VMM."""
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.errorbar(df_vmm["peaking_time"], df_vmm["mean_adc"],
                     yerr=df_vmm["rms_adc"], fmt="o-",
                     label=f"VMM {vmm_id}")
    plt.xlabel("Peaking Time (snt)")
    plt.ylabel("Mean ADC (ADC overthreshold cut)")
    plt.title("Mean ADC vs Peaking Time")
    plt.legend()
    plt.show()


def plot_std_vs_peaking(df_results, vmm_ids):
    """Plot std ADC vs peaking time per VMM."""
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.errorbar(df_vmm["peaking_time"], df_vmm["rms_adc"],
                     yerr=df_vmm["rms_error_adc"], fmt="o-",
                     label=f"VMM {vmm_id}")
    plt.xlabel("Peaking Time (snt)")
    plt.ylabel("Std ADC (ADC over threshold cut)")
    plt.title("Std ADC vs Peaking Time")
    plt.legend()
    plt.show()


def plot_robust_vs_peaking(df_results, vmm_ids):
    """Plot robust noise indicators vs peaking time per VMM."""
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.figure()
        plt.plot(df_vmm["peaking_time"], df_vmm["median_adc"],
                 "o-", label="Median ADC")
        plt.plot(df_vmm["peaking_time"], df_vmm["robust_sigma"],
                 "o-", label="Robust σ (MAD)")
        plt.title(f"VMM {vmm_id} – Robust Noise Indicators")
        plt.xlabel("Peaking time (snt)")
        plt.ylabel("ADC / σ")
        plt.legend()
        plt.show()

    plt.figure()
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.plot(df_vmm["peaking_time"], df_vmm["median_adc"],
                 "o-", label=f"VMM {vmm_id}")
    plt.xlabel("Peaking time (snt)")
    plt.ylabel("Median ADC")
    plt.title("Summary of Median ADC for all VMMs")
    plt.legend()
    plt.show()

    plt.figure()
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.plot(df_vmm["peaking_time"], df_vmm["robust_sigma"],
                 "o-", label=f"VMM {vmm_id}")
    plt.xlabel("Peaking time (snt)")
    plt.ylabel("Robust σ (MAD)")
    plt.title("Summary of Robust σ for all VMMs")
    plt.legend()
    plt.show()


def plot_adc_by_vmm(vmm_ids, run_list, df_run_scan, data_dir):
    """Plot normalised ADC distributions per VMM across runs."""
    from vmm_io import get_run_dir, list_root_files, load_hits_root
    for vmm_id in vmm_ids:
        fig, ax = plt.subplots()
        for run_no in run_list:
            run_dir    = get_run_dir(data_dir, run_no)
            root_files = list_root_files(run_dir, n=2)
            if not root_files:
                continue
            file_path = os.path.join(run_dir, root_files[1])
            df_hits   = load_hits_root(
                file_path,
                branches=["adc", "vmm", "ch",
                          "time", "over_threshold"]
            )
            adc_values = df_hits[df_hits["vmm"] == vmm_id]["adc"]
            param_row  = df_run_scan[df_run_scan["run_no"] == run_no]
            param_str  = ", ".join([
                f"{col}={param_row.iloc[0][col]}"
                for col in df_run_scan.columns
                if col != "run_no"
            ])
            ax.hist(adc_values, bins=50, histtype="step",
                    linewidth=1.5, label=param_str, density=True)

        ax.set_xlabel("ADC")
        ax.set_ylabel("Normalized Counts")
        ax.set_title(f"ADC Distribution for VMM {vmm_id}")
        ax.legend()
        plt.show()


# ─────────────────────────────────────────────
# SNR COMPARISON PLOTS
# ─────────────────────────────────────────────

def _quality_overlay(ax, df_vmm, x_col, y_col):
    """
    Internal helper — overlay warn/bad markers on a line plot.
    Circles = ok, squares = warn, x = bad.
    """
    quality_style = {
        "ok"  : {"marker": "o", "alpha": 1.0},
        "warn": {"marker": "s", "alpha": 0.7},
        "bad" : {"marker": "x", "alpha": 0.5},
    }
    line = ax.plot(df_vmm[x_col], df_vmm[y_col],
                   linestyle="-", marker="o")[0]
    color = line.get_color()
    for _, row in df_vmm.iterrows():
        style = quality_style[row["noise_quality"]]
        ax.plot(row[x_col], row[y_col],
                marker=style["marker"],
                color=color,
                alpha=style["alpha"],
                markersize=8,
                linestyle="")
    return line


def _quality_legend(ax):
    """Add noise quality legend to an axes."""
    handles = [
        Line2D([0], [0], marker="o", color="gray",
               label="ok",   linestyle="-"),
        Line2D([0], [0], marker="s", color="gray",
               label="warn", linestyle="--", alpha=0.7),
        Line2D([0], [0], marker="x", color="gray",
               label="bad",  linestyle=":",  alpha=0.5),
    ]
    ax.legend(handles=handles, title="Noise quality",
              fontsize=8, loc="lower right")


def plot_snr_vs_peaking(df_snr):
    """
    SNR vs peaking time per sg.
    Skips sg values with only one peaking time.
    Marks warn/bad points with different markers.
    """
    df          = df_snr.copy()
    sg_values   = sorted(df["sg"].unique())
    sg_to_plot  = [
        sg for sg in sg_values
        if df[df["sg"] == sg]["snt"].nunique() > 1
    ]
    sg_skip = [sg for sg in sg_values if sg not in sg_to_plot]
    if sg_skip:
        print(f"plot_snr_vs_peaking: skipping sg={sg_skip} "
              f"— only one peaking time available")
    if not sg_to_plot:
        print("plot_snr_vs_peaking: nothing to plot")
        return

    fig, axes = plt.subplots(1, len(sg_to_plot),
                              figsize=(6*len(sg_to_plot), 5),
                              sharey=True)
    if len(sg_to_plot) == 1:
        axes = [axes]

    for ax, sg_val in zip(axes, sg_to_plot):
        df_sg = df[df["sg"] == sg_val]
        for vmm_id in sorted(df_sg["vmm_id"].unique()):
            df_vmm = df_sg[
                df_sg["vmm_id"] == vmm_id
            ].sort_values("snt")
            if len(df_vmm) < 2:
                continue
            line = _quality_overlay(ax, df_vmm, "snt", "snr")
            line.set_label(f"VMM {vmm_id}")

        ax.set_title(f"sg={sg_val}")
        ax.set_xlabel("Peaking time (snt)")
        ax.set_ylabel("SNR (MPV / noise σ)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        _quality_legend(ax)

    plt.suptitle("SNR vs Peaking Time (all data)", y=1.02)
    plt.tight_layout()
    plt.show()


def plot_snr_vs_gain(df_snr):
    """
    SNR vs gain per snt.
    Skips snt values with only one gain setting.
    Marks warn/bad points with different markers.
    """
    df           = df_snr.copy()
    snt_values   = sorted(df["snt"].unique())
    snt_to_plot  = [
        snt for snt in snt_values
        if df[df["snt"] == snt]["sg"].nunique() > 1
    ]
    snt_skip = [snt for snt in snt_values
                if snt not in snt_to_plot]
    if snt_skip:
        print(f"plot_snr_vs_gain: skipping snt={snt_skip} "
              f"— only one gain setting available")
    if not snt_to_plot:
        print("plot_snr_vs_gain: nothing to plot")
        return

    fig, axes = plt.subplots(1, len(snt_to_plot),
                              figsize=(6*len(snt_to_plot), 5),
                              sharey=True)
    if len(snt_to_plot) == 1:
        axes = [axes]

    for ax, snt_val in zip(axes, snt_to_plot):
        df_snt = df[df["snt"] == snt_val]
        for vmm_id in sorted(df_snt["vmm_id"].unique()):
            df_vmm = df_snt[
                df_snt["vmm_id"] == vmm_id
            ].sort_values("sg")
            if len(df_vmm) < 2:
                continue
            line = _quality_overlay(ax, df_vmm, "sg", "snr")
            line.set_label(f"VMM {vmm_id}")

        ax.set_title(f"snt={snt_val:.0f}")
        ax.set_xlabel("Gain (sg)")
        ax.set_ylabel("SNR (MPV / noise σ)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        _quality_legend(ax)

    plt.suptitle("SNR vs Gain (all data)", y=1.02)
    plt.tight_layout()
    plt.show()


def plot_snr_heatmap(df_snr):
    """
    Heatmap of SNR: rows=VMMs, columns=configurations.
    Shows all data — warn/bad cells annotated with quality label.
    Detector group labels on right axis from vmm_mapping.
    Fixed color scale across all cells for direct comparison.
    """
    df = df_snr.copy()

    def config_label(row):
        return f"sg={row['sg']:.1f}\nsnt={row['snt']:.0f}"

    df["config"] = df.apply(config_label, axis=1)

    configs = (
        df[["sg", "snt", "config"]]
        .drop_duplicates()
        .sort_values(["sg", "snt"])["config"]
        .tolist()
    )
    vmm_ids = sorted(df["vmm_id"].unique())

    snr_matrix     = np.full((len(vmm_ids), len(configs)), np.nan)
    quality_matrix = np.full((len(vmm_ids), len(configs)), "",
                              dtype=object)

    for i, vmm_id in enumerate(vmm_ids):
        for j, cfg in enumerate(configs):
            row = df[
                (df["vmm_id"] == vmm_id) &
                (df["config"] == cfg)
            ]
            if not row.empty:
                snr_matrix[i, j]     = row["snr"].iloc[0]
                quality_matrix[i, j] = row["noise_quality"].iloc[0]

    cell_w = 2.2
    cell_h = 1.0
    fig_w  = max(16, len(configs) * cell_w + 8)
    fig_h  = max(6,  len(vmm_ids) * cell_h + 3)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.07, right=0.62,
                        top=0.88, bottom=0.15)

    cmap = plt.cm.RdYlGn.copy()
    cmap.set_bad(color="lightgrey")

    masked = np.ma.masked_invalid(snr_matrix)
    im = ax.imshow(masked, aspect="auto", cmap=cmap,
                   vmin=np.nanmin(snr_matrix),
                   vmax=np.nanmax(snr_matrix))

    for i in range(len(vmm_ids)):
        for j in range(len(configs)):
            val     = snr_matrix[i, j]
            quality = quality_matrix[i, j]
            if np.isnan(val):
                ax.text(j, i, "—", ha="center", va="center",
                        fontsize=9, color="gray")
            else:
                suffix = ("" if quality == "ok"
                          else f"\n({quality})")
                ax.text(j, i, f"{val:.1f}{suffix}",
                        ha="center", va="center",
                        fontsize=8, fontweight="bold")

    ax.set_xticks(range(len(configs)))
    ax.set_xticklabels(configs, fontsize=9)
    ax.set_yticks(range(len(vmm_ids)))
    ax.set_yticklabels([f"VMM {v}" for v in vmm_ids], fontsize=9)
    ax.set_title("SNR heatmap — VMM vs configuration\n"
                 "(warn/bad cells show quality label)",
                 fontsize=11, pad=12)

    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }
    for i, vmm_id in enumerate(vmm_ids):
        label  = vmm_to_detector.get(vmm_id, "")
        x_fig  = 0.635
        y_fig  = (ax.get_position().y0 +
                  ax.get_position().height *
                  (1 - (i + 0.5) / len(vmm_ids)))
        fig.text(x_fig, y_fig, label,
                 ha="left", va="center",
                 fontsize=7, clip_on=False)

    cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.73])
    fig.colorbar(im, cax=cbar_ax, label="SNR")
    plt.show()


def plot_snr_channel_heatmap_per_vmm(df_snr_ch, detector_vmms,
                                      show_quality_overlay=False,
                                      min_noise_hits=50,
                                      min_sigma=2.0,
                                      max_sigma=20.0,
                                      mpv_min=100,
                                      mpv_max=300):
    """
    One heatmap per VMM showing SNR per channel per configuration.
    Rows = channels (0-63), columns = configurations (sg/snt).
    Consistent color scale across all VMMs for direct comparison.
    Gold border highlights the best configuration per VMM.

    Parameters
    ----------
    show_quality_overlay : bool
        If True show all channels with ✕ on low-quality ones.
        If False show only quality-passing channels (grey = excluded).
    """
    df = df_snr_ch.copy()

    configs = (
        df[["sg", "snt"]]
        .drop_duplicates()
        .sort_values(["sg", "snt"])
    )
    config_labels = [
        f"sg={row['sg']:.1f}\nsnt={row['snt']:.0f}"
        for _, row in configs.iterrows()
    ]
    config_keys = [
        (row["sg"], row["snt"])
        for _, row in configs.iterrows()
    ]

    vmin = 5
    vmax = 50
    cmap = plt.cm.RdYlGn.copy()
    cmap.set_bad(color="lightgrey")

    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    for vmm_id in sorted(detector_vmms):
        df_vmm = df[df["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue

        n_configs     = len(config_keys)
        matrix        = np.full((64, n_configs), np.nan)
        quality_mask  = np.zeros((64, n_configs), dtype=bool)

        for j, (sg, snt) in enumerate(config_keys):
            df_cfg = df_vmm[
                (df_vmm["sg"]  == sg) &
                (df_vmm["snt"] == snt)
            ]
            for _, row in df_cfg.iterrows():
                ch  = int(row["ch"])
                snr = row["snr_ch"]
                matrix[ch, j] = np.clip(snr, vmin, vmax)

                if show_quality_overlay:
                    is_bad = (
                        row["n_noise"]        < min_noise_hits or
                        row["noise_sigma_ch"] < min_sigma      or
                        row["noise_sigma_ch"] > max_sigma      or
                        row["mpv_ch"]         < mpv_min        or
                        row["mpv_ch"]         > mpv_max        or
                        row["n_signal"]       < 10
                    )
                    quality_mask[ch, j] = is_bad

        n_valid = (~np.isnan(matrix)).sum(axis=0)

        fig, ax = plt.subplots(figsize=(n_configs * 1.8 + 3, 12))
        fig.subplots_adjust(left=0.08, right=0.80,
                            top=0.92, bottom=0.12)

        masked = np.ma.masked_invalid(matrix)
        im = ax.imshow(masked, aspect="auto", cmap=cmap,
                       vmin=vmin, vmax=vmax)

        if show_quality_overlay:
            for ch in range(64):
                for j in range(n_configs):
                    if quality_mask[ch, j]:
                        ax.text(j, ch, "✕",
                                ha="center", va="center",
                                fontsize=6, color="black",
                                alpha=0.6)

        ax.set_xticks(range(n_configs))
        ax.set_xticklabels(config_labels, fontsize=9)
        ax.set_yticks(range(0, 64, 4))
        ax.set_yticklabels(range(0, 64, 4), fontsize=8)
        ax.set_ylabel("Channel", fontsize=10)

        detector_name = vmm_to_detector.get(vmm_id, "")
        quality_note  = (
            "✕ = low quality estimate"
            if show_quality_overlay
            else "grey = excluded by quality cuts"
        )
        ax.set_title(
            f"VMM {vmm_id} — {detector_name}\n"
            f"Per-channel SNR across configurations\n"
            f"({quality_note}  |  color scale: {vmin}–{vmax})",
            fontsize=11
        )

        for j in range(n_configs):
            col   = matrix[:, j]
            valid = col[~np.isnan(col)]
            if len(valid) > 0:
                ax.text(
                    j, 66,
                    f"med={np.median(valid):.0f}\nn={n_valid[j]}",
                    ha="center", va="top",
                    fontsize=7, color="black"
                )

        best_j = int(np.nanargmax(
            [np.nanmedian(matrix[:, j])
             for j in range(n_configs)]
        ))
        ax.axvline(best_j - 0.5, color="gold",
                   linewidth=2, linestyle="--")
        ax.axvline(best_j + 0.5, color="gold",
                   linewidth=2, linestyle="--")
        ax.text(best_j, -2, "best",
                ha="center", va="bottom",
                fontsize=8, color="goldenrod",
                fontweight="bold")

        cbar_ax = fig.add_axes([0.82, 0.12, 0.02, 0.80])
        fig.colorbar(im, cax=cbar_ax, label="SNR")
        plt.show()


def plot_snr_channel_heatmap_all_configs(df_snr_ch):
    """
    One heatmap per configuration showing SNR per channel per VMM.
    Rows = channels (0-63), columns = VMMs.
    Useful for spatial detector maps — shows which detector
    regions perform best at each configuration.
    """
    df = df_snr_ch.copy()

    configs = (
        df[["sg", "snt"]]
        .drop_duplicates()
        .sort_values(["sg", "snt"])
        .values.tolist()
    )
    vmm_ids = sorted(df["vmm_id"].unique())

    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    vmin = 5
    vmax = 50

    for sg, snt in configs:
        df_cfg = df[
            (df["sg"]  == sg) &
            (df["snt"] == snt)
        ]
        matrix = np.full((64, len(vmm_ids)), np.nan)
        for j, vmm_id in enumerate(vmm_ids):
            df_vmm = df_cfg[df_cfg["vmm_id"] == vmm_id]
            for _, row in df_vmm.iterrows():
                matrix[int(row["ch"]), j] = row["snr_ch"]

        fig, ax = plt.subplots(
            figsize=(len(vmm_ids) * 1.8 + 3, 14)
        )
        fig.subplots_adjust(left=0.07, right=0.80,
                            top=0.92, bottom=0.08)

        cmap = plt.cm.RdYlGn.copy()
        cmap.set_bad(color="lightgrey")

        masked = np.ma.masked_invalid(matrix)
        im = ax.imshow(masked, aspect="auto", cmap=cmap,
                       vmin=vmin, vmax=vmax)

        ax.set_xticks(range(len(vmm_ids)))
        ax.set_xticklabels(
            [f"VMM {v}" for v in vmm_ids], fontsize=9
        )
        ax.set_ylabel("Channel", fontsize=10)
        ax.set_yticks(range(0, 64, 4))
        ax.set_title(
            f"Per-channel SNR — sg={sg:.1f}  snt={snt:.0f}\n"
            f"(grey = excluded  |  color scale: {vmin}–{vmax})",
            fontsize=11
        )

        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(range(len(vmm_ids)))
        ax2.set_xticklabels(
            [vmm_to_detector.get(v, "") for v in vmm_ids],
            fontsize=7, rotation=15, ha="left"
        )

        for j, vmm_id in enumerate(vmm_ids):
            col   = matrix[:, j]
            valid = col[~np.isnan(col)]
            if len(valid) > 0:
                ax.text(j, 65,
                        f"med={np.median(valid):.0f}",
                        ha="center", va="top",
                        fontsize=7, color="black")

        cbar_ax = fig.add_axes([0.82, 0.08, 0.02, 0.84])
        fig.colorbar(im, cax=cbar_ax, label="SNR")
        plt.show()


def plot_snr_channel_uniformity(df_snr_ch, detector_vmms):
    """
    One figure per VMM showing SNR distribution across channels
    per configuration as box plots.
    Answers: which config gives best and most uniform channel SNR?
    """
    df = df_snr_ch.copy()
    df["config"] = df.apply(
        lambda r: f"sg={r['sg']:.1f}\nsnt={r['snt']:.0f}",
        axis=1
    )

    config_order = (
        df[["sg", "snt", "config"]]
        .drop_duplicates()
        .sort_values(["sg", "snt"])["config"]
        .tolist()
    )

    for vmm_id in sorted(detector_vmms):
        df_vmm = df[df["vmm_id"] == vmm_id]

        data   = []
        labels = []
        for cfg in config_order:
            vals = df_vmm[df_vmm["config"] == cfg]["snr_ch"].values
            if len(vals) > 0:
                data.append(vals)
                labels.append(cfg)

        if not data:
            continue

        fig, ax = plt.subplots(
            figsize=(len(data) * 1.8 + 2, 5)
        )

        bp = ax.boxplot(
            data,
            tick_labels=labels,
            patch_artist=True,
            medianprops={"color": "black", "linewidth": 2},
            flierprops={"marker": ".", "markersize": 4,
                        "alpha": 0.5}
        )

        colors = plt.cm.tab10(np.linspace(0, 1, len(data)))
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        for i, (vals, label) in enumerate(
            zip(data, labels), start=1
        ):
            ax.text(i, ax.get_ylim()[1] * 0.97,
                    f"n={len(vals)}\nmed={np.median(vals):.0f}",
                    ha="center", va="top", fontsize=7)

        ax.set_xlabel("Configuration")
        ax.set_ylabel("SNR per channel")
        ax.set_title(
            f"VMM {vmm_id} — channel SNR distribution "
            f"per configuration\n"
            f"(box = IQR, line = median, dots = outliers)"
        )
        ax.grid(True, alpha=0.3, axis="y")
        plt.tight_layout()
        plt.show()

def main():
    print('bonzo')


if __name__ == '__main__':
    main()
