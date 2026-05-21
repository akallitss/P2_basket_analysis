#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vmm_plots.py

Plotting functions for VMM config scan analysis.

Two categories:
- Legacy plots : ADC distributions, mean/std vs peaking time
- SNR plots    : configuration comparison at VMM and channel level

All plot functions accept out_dir and show parameters:
    out_dir : str or None — if set, saves PDF + PNG to that directory
    show    : bool        — if True, displays the figure interactively
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from vmm_mapping import vmm_mapping


# ─────────────────────────────────────────────
# GLOBAL PLOT STYLE
# ─────────────────────────────────────────────
plt.rcParams.update({
    "font.size"         : 13,
    "font.weight"       : "bold",
    "axes.titlesize"    : 15,
    "axes.titleweight"  : "bold",
    "axes.labelsize"    : 13,
    "axes.labelweight"  : "bold",
    "xtick.labelsize"   : 12,
    "ytick.labelsize"   : 12,
    "legend.fontsize"   : 12,
    "figure.titlesize"  : 15,
    "figure.titleweight": "bold",
})


# ─────────────────────────────────────────────
# FIGURE SAVE / SHOW HELPER
# ─────────────────────────────────────────────
def _finish_fig(fig, stem, out_dir, show):
    """Save fig as PDF and PNG to out_dir, optionally display, then close."""
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        for ext in ("pdf", "png"):
            fig.savefig(os.path.join(out_dir, f"{stem}.{ext}"),
                        bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


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
                                  get_run_dir_fn=None,
                                  list_root_files_fn=None,
                                  load_hits_root_fn=None,
                                  out_dir=None, show=True):
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
    _finish_fig(fig, "adc_hist_per_run", out_dir, show)


def compare_full_vs_cut(df_hits, vmm_id, run_no,
                         adc_cut=100, use_over_threshold=True,
                         out_dir=None, show=True):
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

    fig = plt.figure()
    plt.hist(full_adc, bins=80, alpha=0.3,
             label="Full", density=False)
    plt.hist(cut_adc, bins=80, alpha=0.6,
             label=label_cut, density=False)
    plt.title(f"Run {run_no} – VMM {vmm_id}")
    plt.xlabel("ADC")
    plt.ylabel("Numb of Hits")
    plt.legend()
    _finish_fig(fig, f"compare_full_vs_cut_vmm{vmm_id}_run{run_no}",
                out_dir, show)


def plot_removed_fraction(df_results, vmm_ids, out_dir=None, show=True):
    """Plot fraction of hits removed by ADC cut vs peaking time."""
    fig = plt.figure()
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
    _finish_fig(fig, "removed_fraction", out_dir, show)


def plot_mean_vs_peaking(df_results, vmm_ids, out_dir=None, show=True):
    """Plot mean ADC vs peaking time per VMM."""
    fig = plt.figure()
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
    _finish_fig(fig, "mean_vs_peaking", out_dir, show)


def plot_std_vs_peaking(df_results, vmm_ids, out_dir=None, show=True):
    """Plot std ADC vs peaking time per VMM."""
    fig = plt.figure()
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
    _finish_fig(fig, "std_vs_peaking", out_dir, show)


def plot_robust_vs_peaking(df_results, vmm_ids, out_dir=None, show=True):
    """Plot robust noise indicators vs peaking time per VMM."""
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        fig = plt.figure()
        plt.plot(df_vmm["peaking_time"], df_vmm["median_adc"],
                 "o-", label="Median ADC")
        plt.plot(df_vmm["peaking_time"], df_vmm["robust_sigma"],
                 "o-", label="Robust σ (MAD)")
        plt.title(f"VMM {vmm_id} – Robust Noise Indicators")
        plt.xlabel("Peaking time (snt)")
        plt.ylabel("ADC / σ")
        plt.legend()
        _finish_fig(fig, f"robust_indicators_vmm{vmm_id}", out_dir, show)

    fig = plt.figure()
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
    _finish_fig(fig, "robust_median_summary", out_dir, show)

    fig = plt.figure()
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
    _finish_fig(fig, "robust_sigma_summary", out_dir, show)


def plot_adc_by_vmm(vmm_ids, run_list, df_run_scan, data_dir,
                     out_dir=None, show=True):
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
        _finish_fig(fig, f"adc_by_vmm_vmm{vmm_id}", out_dir, show)


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
              fontsize=14, loc="lower right")


def plot_snr_vs_peaking(df_snr, snr_col="snr", out_dir=None, show=True):
    """
    SNR vs peaking time per sg.
    Skips sg values with only one peaking time.
    Marks warn/bad points with different markers.
    """
    method      = "mean" if "mean" in snr_col else "MPV"
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
            line = _quality_overlay(ax, df_vmm, "snt", snr_col)
            line.set_label(f"VMM {vmm_id}")

        ax.set_title(f"sg={sg_val}")
        ax.set_xlabel("Peaking time (snt)")
        ax.set_ylabel(f"SNR ({method} / noise σ)")
        ax.legend(fontsize=14)
        ax.grid(True, alpha=0.3)
        _quality_legend(ax)

    plt.suptitle(f"SNR vs Peaking Time — {method} method",
                 fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    _finish_fig(fig, f"snr_vs_peaking_{method.lower()}", out_dir, show)


def plot_snr_vs_gain(df_snr, snr_col="snr", out_dir=None, show=True):
    """
    SNR vs gain per snt.
    Skips snt values with only one gain setting.
    Marks warn/bad points with different markers.
    """
    method       = "mean" if "mean" in snr_col else "MPV"
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
            line = _quality_overlay(ax, df_vmm, "sg", snr_col)
            line.set_label(f"VMM {vmm_id}")

        ax.set_title(f"snt={snt_val:.0f}")
        ax.set_xlabel("Gain (sg)")
        ax.set_ylabel(f"SNR ({method} / noise σ)")
        ax.legend(fontsize=14)
        ax.grid(True, alpha=0.3)
        _quality_legend(ax)

    plt.suptitle(f"SNR vs Gain — {method} method", fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    _finish_fig(fig, f"snr_vs_gain_{method.lower()}", out_dir, show)


def plot_snr_heatmap(df_snr, snr_col="snr", out_dir=None, show=True):
    """
    Heatmap of SNR: rows=VMMs, columns=configurations.
    Shows all data — warn/bad cells annotated with quality label.
    Detector group labels on right axis from vmm_mapping.
    Fixed color scale across all cells for direct comparison.
    """
    method = "mean" if "mean" in snr_col else "MPV"
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
                snr_matrix[i, j]     = row[snr_col].iloc[0]
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
                        fontsize=12, color="gray")
            else:
                suffix = ("" if quality == "ok"
                          else f"\n({quality})")
                ax.text(j, i, f"{val:.1f}{suffix}",
                        ha="center", va="center",
                        fontsize=14, fontweight="bold")

    ax.set_xticks(range(len(configs)))
    ax.set_xticklabels(configs, fontsize=12)
    ax.set_yticks(range(len(vmm_ids)))
    ax.set_yticklabels([f"VMM {v}" for v in vmm_ids], fontsize=12)
    ax.set_title(f"SNR heatmap ({method} method) — VMM vs configuration\n"
                 "(warn/bad cells show quality label)",
                 fontsize=14, pad=12)

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
                 fontsize=13, clip_on=False)

    cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.73])
    fig.colorbar(im, cax=cbar_ax, label=f"SNR ({method})")
    _finish_fig(fig, f"snr_heatmap_{method.lower()}", out_dir, show)


def plot_adc_heatmap(df_snr, metric="mpv", out_dir=None, show=True):
    """
    Heatmap of ADC metric: rows=VMMs, columns=configurations.

    Parameters
    ----------
    metric : str
        'mpv'         — signal most probable value (ADC units)
        'noise_sigma' — noise width (ADC units)
    """
    if metric not in ("mpv", "noise_sigma", "mean_signal", "mean_noise"):
        raise ValueError(
            "metric must be 'mpv', 'noise_sigma', 'mean_signal', or 'mean_noise'"
        )

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

    val_matrix     = np.full((len(vmm_ids), len(configs)), np.nan)
    quality_matrix = np.full((len(vmm_ids), len(configs)), "",
                              dtype=object)

    for i, vmm_id in enumerate(vmm_ids):
        for j, cfg in enumerate(configs):
            row = df[
                (df["vmm_id"] == vmm_id) &
                (df["config"] == cfg)
            ]
            if not row.empty:
                val_matrix[i, j]     = row[metric].iloc[0]
                quality_matrix[i, j] = row["noise_quality"].iloc[0]

    cell_w = 2.2
    cell_h = 1.0
    fig_w  = max(16, len(configs) * cell_w + 8)
    fig_h  = max(6,  len(vmm_ids) * cell_h + 3)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.07, right=0.62,
                        top=0.88, bottom=0.15)

    cmap = (plt.cm.Blues.copy()   if metric in ("mpv", "mean_signal")
            else plt.cm.Oranges.copy())
    cmap.set_bad(color="lightgrey")

    masked = np.ma.masked_invalid(val_matrix)
    im = ax.imshow(masked, aspect="auto", cmap=cmap,
                   vmin=np.nanmin(val_matrix),
                   vmax=np.nanmax(val_matrix))

    for i in range(len(vmm_ids)):
        for j in range(len(configs)):
            val     = val_matrix[i, j]
            quality = quality_matrix[i, j]
            if np.isnan(val):
                ax.text(j, i, "—", ha="center", va="center",
                        fontsize=12, color="gray")
            else:
                suffix = ("" if quality == "ok"
                          else f"\n({quality})")
                ax.text(j, i, f"{val:.1f}{suffix}",
                        ha="center", va="center",
                        fontsize=14, fontweight="bold")

    ax.set_xticks(range(len(configs)))
    ax.set_xticklabels(configs, fontsize=12)
    ax.set_yticks(range(len(vmm_ids)))
    ax.set_yticklabels([f"VMM {v}" for v in vmm_ids], fontsize=12)

    label_str = {
        "mpv"         : "MPV [ADC]",
        "noise_sigma" : "Noise σ [ADC]",
        "mean_signal" : "Mean signal [ADC]",
        "mean_noise"  : "Mean noise [ADC]",
    }[metric]
    title_str = {
        "mpv"         : "ADC MPV heatmap — VMM vs configuration",
        "noise_sigma" : "Noise σ heatmap — VMM vs configuration",
        "mean_signal" : "Mean signal heatmap — VMM vs configuration",
        "mean_noise"  : "Mean noise heatmap — VMM vs configuration",
    }[metric]
    ax.set_title(f"{title_str}\n(warn/bad cells show quality label)",
                 fontsize=14, pad=12)

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
                 fontsize=13, clip_on=False)

    cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.73])
    fig.colorbar(im, cax=cbar_ax, label=label_str)
    _finish_fig(fig, f"adc_heatmap_{metric}", out_dir, show)


def plot_snr_channel_heatmap_per_vmm(df_snr_ch, detector_vmms,
                                      snr_col="snr_ch",
                                      show_quality_overlay=False,
                                      min_noise_hits=50,
                                      min_sigma=2.0,
                                      max_sigma=20.0,
                                      mpv_min=100,
                                      mpv_max=300,
                                      out_dir=None, show=True):
    """
    One heatmap per VMM showing SNR per channel per configuration.
    Rows = channels (0-63), columns = configurations (sg/snt).
    Consistent color scale across all VMMs for direct comparison.
    Gold border highlights the best configuration per VMM.
    """
    method = "mean" if "mean" in snr_col else "MPV"
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
                snr = row[snr_col]
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
                                fontsize=12, color="black",
                                alpha=0.6)

        ax.set_xticks(range(n_configs))
        ax.set_xticklabels(
            [f"{lbl}\nn={n}" for lbl, n in zip(config_labels, n_valid)],
            fontsize=11
        )
        ax.set_yticks(range(0, 64, 4))
        ax.set_yticklabels(range(0, 64, 4), fontsize=14)
        ax.set_ylabel("Channel", fontsize=13)

        detector_name = vmm_to_detector.get(vmm_id, "")
        quality_note  = (
            "✕ = low quality estimate"
            if show_quality_overlay
            else "grey = excluded by quality cuts"
        )
        ax.set_title(
            f"VMM {vmm_id} — {detector_name}\n"
            f"Per-channel SNR ({method}) across configurations\n"
            f"({quality_note}  |  color scale: {vmin}–{vmax})",
            fontsize=14
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
                fontsize=14, color="goldenrod",
                fontweight="bold")

        cbar_ax = fig.add_axes([0.82, 0.12, 0.02, 0.80])
        fig.colorbar(im, cax=cbar_ax, label=f"SNR ({method})")
        _finish_fig(fig,
                    f"snr_channel_heatmap_{method.lower()}_vmm{vmm_id}",
                    out_dir, show)


def plot_snr_channel_heatmap_all_configs(df_snr_ch, snr_col="snr_ch",
                                          out_dir=None, show=True):
    """
    One heatmap per configuration showing SNR per channel per VMM.
    Rows = channels (0-63), columns = VMMs.
    """
    method = "mean" if "mean" in snr_col else "MPV"
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
                matrix[int(row["ch"]), j] = row[snr_col]

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

        n_valid = [
            int((~np.isnan(matrix[:, j])).sum())
            for j in range(len(vmm_ids))
        ]
        ax.set_xticks(range(len(vmm_ids)))
        ax.set_xticklabels(
            [f"VMM {v}\nn={n}" for v, n in zip(vmm_ids, n_valid)],
            fontsize=11
        )
        ax.set_ylabel("Channel", fontsize=13)
        ax.set_yticks(range(0, 64, 4))
        ax.set_title(
            f"Per-channel SNR ({method}) — sg={sg:.1f}  snt={snt:.0f}\n"
            f"(grey = excluded  |  color scale: {vmin}–{vmax})",
            fontsize=14
        )

        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(range(len(vmm_ids)))
        ax2.set_xticklabels(
            [vmm_to_detector.get(v, "") for v in vmm_ids],
            fontsize=13, rotation=15, ha="left"
        )

        cbar_ax = fig.add_axes([0.82, 0.08, 0.02, 0.84])
        fig.colorbar(im, cax=cbar_ax, label=f"SNR ({method})")

        stem = (f"snr_channel_all_configs_{method.lower()}_"
                f"sg{sg:.1f}_snt{snt:.0f}".replace(".", "p"))
        _finish_fig(fig, stem, out_dir, show)


def plot_snr_channel_uniformity(df_snr_ch, detector_vmms,
                                 snr_col="snr_ch",
                                 out_dir=None, show=True):
    """
    One figure per VMM showing SNR distribution across channels
    per configuration as box plots.
    """
    method = "mean" if "mean" in snr_col else "MPV"
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
            vals = df_vmm[df_vmm["config"] == cfg][snr_col].values
            vals = vals[~np.isnan(vals)]
            if len(vals) > 0:
                data.append(vals)
                labels.append(f"{cfg}\nn={len(vals)}")

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

        ax.set_xlabel("Configuration")
        ax.set_ylabel(f"SNR ({method}) per channel")
        ax.set_title(
            f"VMM {vmm_id} — channel SNR ({method}) distribution "
            f"per configuration\n"
            f"(box = IQR, line = median, dots = outliers)"
        )
        ax.grid(True, alpha=0.3, axis="y")
        plt.tight_layout()
        _finish_fig(fig,
                    f"snr_channel_uniformity_{method.lower()}_vmm{vmm_id}",
                    out_dir, show)


def plot_snr_method_comparison(df_snr, df_snr_ch,
                                out_dir=None, show=True):
    """
    Two-panel comparison of MPV-based vs mean-charge SNR.

    Top panel   : VMM-level scatter — SNR_MPV vs SNR_mean,
                  one point per (VMM, config), colour = config.
    Bottom panel: channel-level scatter — snr_ch vs snr_mean_ch,
                  one point per (VMM, channel, config).
    Identity line y=x on both panels for visual alignment check.
    """
    config_labels = (
        df_snr[["sg", "snt"]]
        .drop_duplicates()
        .sort_values(["sg", "snt"])
        .apply(lambda r: f"sg={r['sg']:.1f} snt={r['snt']:.0f}", axis=1)
        .tolist()
    )
    palette = plt.cm.tab10(np.linspace(0, 1, len(config_labels)))
    color_map = dict(zip(config_labels, palette))

    def cfg_label(row):
        return f"sg={row['sg']:.1f} snt={row['snt']:.0f}"

    fig, (ax_vmm, ax_ch) = plt.subplots(1, 2, figsize=(14, 6))
    fig.subplots_adjust(wspace=0.35)

    # ── VMM level ──────────────────────────────────────────────
    df_v = df_snr.dropna(subset=["snr", "snr_mean"]).copy()
    df_v["cfg"] = df_v.apply(cfg_label, axis=1)

    for cfg in config_labels:
        sub = df_v[df_v["cfg"] == cfg]
        if sub.empty:
            continue
        ax_vmm.scatter(sub["snr"], sub["snr_mean"],
                       color=color_map[cfg], label=cfg,
                       s=60, alpha=0.8, zorder=3)

    lim_v = [0, max(df_v[["snr", "snr_mean"]].max()) * 1.1]
    ax_vmm.plot(lim_v, lim_v, "k--", lw=1, alpha=0.4, label="y = x")
    ax_vmm.set_xlim(lim_v)
    ax_vmm.set_ylim(lim_v)
    ax_vmm.set_xlabel("SNR  (MPV / noise σ)")
    ax_vmm.set_ylabel("SNR  (mean charge / noise σ)")
    ax_vmm.set_title("VMM-level SNR: MPV vs mean-charge")
    ax_vmm.legend(fontsize=10, title="Config")
    ax_vmm.grid(True, alpha=0.3)

    # ── Channel level ───────────────────────────────────────────
    df_c = df_snr_ch.dropna(subset=["snr_ch", "snr_mean_ch"]).copy()
    df_c["cfg"] = df_c.apply(cfg_label, axis=1)

    for cfg in config_labels:
        sub = df_c[df_c["cfg"] == cfg]
        if sub.empty:
            continue
        ax_ch.scatter(sub["snr_ch"], sub["snr_mean_ch"],
                      color=color_map[cfg], label=cfg,
                      s=8, alpha=0.3, zorder=3)

    all_vals = df_c[["snr_ch", "snr_mean_ch"]].values
    lim_c = [0, np.nanmax(all_vals) * 1.1]
    ax_ch.plot(lim_c, lim_c, "k--", lw=1, alpha=0.4, label="y = x")
    ax_ch.set_xlim(lim_c)
    ax_ch.set_ylim(lim_c)
    ax_ch.set_xlabel("SNR  (MPV / noise σ)")
    ax_ch.set_ylabel("SNR  (mean charge / noise σ)")
    ax_ch.set_title("Channel-level SNR: MPV vs mean-charge")
    ax_ch.legend(fontsize=10, title="Config",
                 markerscale=2)
    ax_ch.grid(True, alpha=0.3)

    plt.suptitle("MPV-based vs mean-charge SNR comparison",
                 fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    _finish_fig(fig, "snr_comparison_mpv_vs_mean", out_dir, show)


def main():
    pass


if __name__ == '__main__':
    main()