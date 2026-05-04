#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vmm_qa.py

Quality assurance and validation functions for VMM config scan analysis.

All qa_ functions accept out_dir and show parameters:
    out_dir : str or None — if set, saves PDF + PNG to that directory
    show    : bool        — if True, displays the figure interactively

Functions are grouped by investigation type:
- Over-threshold flag validation
- Noise pedestal stability and quality
- Signal distribution inspection
- ADC=16 digital artifact investigation
- Estimator validation (MPV vs median, robust sigma vs std)
- Noise run diagnostic (too few over_threshold=0 hits)
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from vmm_io     import get_run_dir, load_hits_run
from vmm_noise  import compute_noise_baseline
from vmm_signal import get_clean_signal, estimate_mpv
from vmm_mapping import vmm_mapping


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


def qa_over_threshold_split(df_hits, run_no, vmm_ids,
                              out_dir=None, show=True):
    """
    Investigate the over_threshold flag split for a given run.
    Shows ADC distributions split by over_threshold flag per VMM.
    """
    for vmm_id in vmm_ids:
        df_vmm  = df_hits[df_hits["vmm"] == vmm_id]
        n_total = len(df_vmm)
        if n_total == 0:
            continue

        n_above = (df_vmm["over_threshold"] == 1).sum()
        n_below = (df_vmm["over_threshold"] == 0).sum()

        print(f"\nRun {run_no} | VMM {vmm_id}")
        print(f"  Total hits      : {n_total}")
        print(f"  over_threshold=1: {n_above} "
              f"({100*n_above/n_total:.1f}%)")
        print(f"  over_threshold=0: {n_below} "
              f"({100*n_below/n_total:.1f}%)")

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        for flag, label, color in [
            (1, "over_threshold=1", "steelblue"),
            (0, "over_threshold=0", "orange")
        ]:
            adc = df_vmm[df_vmm["over_threshold"] == flag]["adc"]
            axes[0].hist(adc, bins=100, alpha=0.6,
                         label=f"{label} (n={len(adc)})",
                         color=color)
            axes[1].hist(adc, bins=100, alpha=0.6,
                         label=label, color=color)

        axes[0].set_title(f"Run {run_no} VMM {vmm_id} — full range")
        axes[0].set_xlabel("ADC")
        axes[0].set_ylabel("Hits")
        axes[0].legend()
        axes[1].set_xlim(0, 200)
        axes[1].set_title(f"Run {run_no} VMM {vmm_id} — zoom low ADC")
        axes[1].set_xlabel("ADC")
        axes[1].set_ylabel("Hits")
        axes[1].legend()

        plt.tight_layout()
        _finish_fig(fig, f"qa_ot_split_vmm{vmm_id}_run{run_no}",
                    out_dir, show)


def qa_noise_pedestal_stability(df_run_scan, data_dir,
                                 sng1_runs, vmm_groups, n_files=1,
                                 out_dir=None, show=True):
    """
    Check noise pedestal stability across sng=1 runs.
    Plots normalised noise distributions per VMM grouped by peaking time.
    """
    for group in vmm_groups:
        vmm_ids = group["vmm_ids"]
        label   = group["label"]

        fig, axes = plt.subplots(len(vmm_ids), 1,
                                  figsize=(10, 3.5*len(vmm_ids)))
        if len(vmm_ids) == 1:
            axes = [axes]

        for idx, vmm_id in enumerate(vmm_ids):
            for run_no in sng1_runs:
                run_dir = get_run_dir(data_dir, run_no)
                df_hits = load_hits_run(
                    run_dir, n_files=n_files,
                    branches=["adc", "vmm", "over_threshold"]
                )
                if df_hits is None:
                    continue
                noise = df_hits[
                    (df_hits["vmm"]            == vmm_id) &
                    (df_hits["over_threshold"] == 0)
                ]["adc"].values

                if len(noise) < 100:
                    continue

                snt_ser = df_run_scan.loc[df_run_scan["run_no"] == run_no, "snt"]
                snt_lbl = f"{snt_ser.iloc[0]:.0f}" if not snt_ser.empty else "?"

                axes[idx].hist(
                    noise, bins=80, range=(0, 200),
                    alpha=0.5, histtype="step", linewidth=1.5,
                    density=True,
                    label=f"Run {run_no} (snt={snt_lbl})"
                )

            axes[idx].set_title(f"VMM {vmm_id}")
            axes[idx].set_xlabel("ADC")
            axes[idx].set_ylabel("Density")
            axes[idx].legend(fontsize=9)

        plt.suptitle(
            f"Noise pedestal stability — {label} — sng=1"
        )
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        stem = f"qa_noise_pedestal_{label.replace(' ', '_')}"
        _finish_fig(fig, stem, out_dir, show)


def qa_signal_distributions(df_hits, detector_vmms, run_no,
                              out_dir=None, show=True):
    """
    Plot signal ADC distributions per VMM for a given run.
    Also prints saturation statistics per VMM.
    """
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes      = axes.flatten()
    plot_idx  = 0

    print(f"\nRun {run_no} — Saturation check per VMM:")
    print(f"{'VMM':>5} {'N_signal':>10} {'N_saturated':>12} "
          f"{'Sat%':>8} {'ADC_max':>10}")

    for vmm_id in sorted(detector_vmms):
        signal = df_hits[
            (df_hits["vmm"]            == vmm_id) &
            (df_hits["over_threshold"] == 1)
        ]["adc"].values

        if len(signal) < 100:
            continue

        adc_max = signal.max()
        n_sat   = (signal == adc_max).sum()
        print(f"{vmm_id:>5} {len(signal):>10} {n_sat:>12} "
              f"{100*n_sat/len(signal):>8.2f}% {adc_max:>10}")

        if plot_idx >= len(axes):
            break

        axes[plot_idx].hist(signal, bins=100,
                             color="steelblue", alpha=0.7)
        axes[plot_idx].set_xlim(100, 700)
        axes[plot_idx].set_xlabel("ADC")
        axes[plot_idx].set_ylabel("Hits")
        axes[plot_idx].set_title(f"VMM {vmm_id} (n={len(signal)})")
        plot_idx += 1

    plt.suptitle(f"Run {run_no} — Signal distribution per VMM")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    _finish_fig(fig, f"qa_signal_distributions_run{run_no}",
                out_dir, show)


def qa_noise_quality_check(df_run_scan, data_dir, runs_to_check,
                             n_files=1):
    """
    Run noise quality check across a list of runs.
    Prints flagged VMMs with their sigma and noise_cut values.
    """
    for run_no in runs_to_check:
        run_dir = get_run_dir(data_dir, run_no)
        df_hits = load_hits_run(
            run_dir, n_files=n_files,
            branches=["adc", "vmm", "over_threshold"]
        )
        if df_hits is None:
            continue
        df_noise = compute_noise_baseline(df_hits)
        snt_ser = df_run_scan.loc[df_run_scan["run_no"] == run_no, "snt"]
        if snt_ser.empty:
            print(f"Run {run_no} — not found in run scan table, skipping")
            continue
        snt = snt_ser.iloc[0]

        if "quality" not in df_noise.columns:
            print(f"Run {run_no} (snt={snt:.0f}) — "
                  f"no noise baseline (no over_threshold=0 hits)")
            continue

        bad  = df_noise[df_noise["quality"] == "bad"]
        warn = df_noise[df_noise["quality"] == "warn"]

        if bad.empty and warn.empty:
            print(f"Run {run_no} (snt={snt:.0f}) — all VMMs OK")
        else:
            print(f"\nRun {run_no} (snt={snt:.0f}):")
            if not bad.empty:
                print("  BAD VMMs:")
                print(bad[["vmm_id", "robust_sigma",
                            "noise_cut",
                            "quality"]].to_string(index=False))
            if not warn.empty:
                print("  WARN VMMs:")
                print(warn[["vmm_id", "robust_sigma",
                             "noise_cut",
                             "quality"]].to_string(index=False))


def qa_mpv_estimation(df_hits, df_noise_baseline,
                       detector_vmms, run_no,
                       exclude_trigger_vmms=(0, 1),
                       out_dir=None, show=True):
    """
    Visualise MPV estimation for each detector VMM.
    Shows smoothed signal distribution with noise cut and MPV marker.
    """
    for vmm_id in detector_vmms:
        signal_clean = get_clean_signal(
            df_hits, vmm_id, exclude_trigger_vmms
        )
        if signal_clean is None or len(signal_clean) < 100:
            continue

        row = df_noise_baseline[
            df_noise_baseline["vmm_id"] == vmm_id
        ]
        if row.empty:
            print(f"VMM {vmm_id} — no noise baseline, skipping")
            continue

        noise_cut     = row["noise_cut"].iloc[0]
        noise_quality = row["quality"].iloc[0]

        mpv, counts, smoothed, bin_centers = estimate_mpv(
            signal_clean, adc_min=noise_cut
        )
        scale = counts.max() / smoothed.max()

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.bar(bin_centers, counts,
               width=bin_centers[1] - bin_centers[0],
               color="steelblue", alpha=0.5,
               label="Clean signal")
        ax.plot(bin_centers, smoothed * scale,
                color="black", linewidth=2, label="Smoothed")
        ax.axvline(noise_cut, color="orange", linestyle="--",
                   label=f"Noise cut={noise_cut:.0f} [{noise_quality}]")
        ax.axvline(mpv, color="red", linestyle="--",
                   label=f"MPV={mpv:.0f}")
        ax.set_xlabel("ADC")
        ax.set_ylabel("Hits")
        ax.set_title(f"VMM {vmm_id} Run {run_no} — MPV estimation")
        ax.legend()
        plt.tight_layout()
        _finish_fig(fig, f"qa_mpv_vmm{vmm_id}_run{run_no}",
                    out_dir, show)

        print(f"VMM {vmm_id} — noise_cut={noise_cut:.1f} "
              f"[{noise_quality}]  MPV={mpv:.1f} ADC  "
              f"N_signal={len(signal_clean)}")


def qa_channel_noise(df_hits, detector_vmms, run_no,
                      out_dir=None, show=True):
    """
    Per-channel noise investigation.
    Shows median ADC and robust sigma per channel per VMM.
    """
    for vmm_id in detector_vmms:
        df_vmm_noise = df_hits[
            (df_hits["vmm"]            == vmm_id) &
            (df_hits["over_threshold"] == 0)
        ]
        if df_vmm_noise.empty:
            continue

        ch_stats = []
        for ch_id, df_ch in df_vmm_noise.groupby("ch"):
            adc = df_ch["adc"].values
            if len(adc) < 10:
                continue
            median = np.median(adc)
            mad    = np.median(np.abs(adc - median))
            ch_stats.append({
                "ch"          : ch_id,
                "n_hits"      : len(adc),
                "median_adc"  : median,
                "robust_sigma": 1.4826 * mad
            })

        if not ch_stats:
            continue

        df_ch_stats = pd.DataFrame(ch_stats)

        fig, axes = plt.subplots(1, 2, figsize=(14, 4))

        axes[0].bar(df_ch_stats["ch"], df_ch_stats["median_adc"],
                    color="steelblue", alpha=0.7)
        axes[0].axhline(
            df_ch_stats["median_adc"].median(),
            color="red", linestyle="--"
        )
        axes[0].set_xlabel("Channel")
        axes[0].set_ylabel("Median ADC")
        axes[0].set_title(
            f"VMM {vmm_id} Run {run_no} — Noise median per channel"
        )

        axes[1].bar(df_ch_stats["ch"], df_ch_stats["robust_sigma"],
                    color="orange", alpha=0.7)
        axes[1].axhline(
            df_ch_stats["robust_sigma"].median(),
            color="red", linestyle="--"
        )
        axes[1].set_xlabel("Channel")
        axes[1].set_ylabel("Robust sigma (ADC)")
        axes[1].set_title(
            f"VMM {vmm_id} Run {run_no} — Noise sigma per channel"
        )

        plt.suptitle(
            f"VMM {vmm_id} — per channel noise — Run {run_no}"
        )
        plt.tight_layout()
        _finish_fig(fig, f"qa_channel_noise_vmm{vmm_id}_run{run_no}",
                    out_dir, show)


def qa_adc16_artifact(df_hits, detector_vmms, run_no,
                       artifact_adc=16, low_adc_threshold=25,
                       out_dir=None, show=True):
    """
    Investigate the ADC=16 digital artifact.
    Shows which channels are affected and how uniformly.
    """
    print(f"\nRun {run_no} — ADC={artifact_adc} artifact "
          f"investigation:")

    for vmm_id in detector_vmms:
        df_vmm_noise = df_hits[
            (df_hits["vmm"]            == vmm_id) &
            (df_hits["over_threshold"] == 0)
        ]
        if df_vmm_noise.empty:
            continue

        total_noise = len(df_vmm_noise)
        artifact    = df_vmm_noise[
            df_vmm_noise["adc"] == artifact_adc
        ]
        n_artifact  = len(artifact)

        if n_artifact == 0:
            print(f"  VMM {vmm_id} — no ADC={artifact_adc} hits")
            continue

        ch_counts = (artifact.groupby("ch")
                              .size()
                              .sort_values(ascending=False))

        print(f"\n  VMM {vmm_id}:")
        print(f"    Artifact hits     : {n_artifact} "
              f"({100*n_artifact/total_noise:.1f}% of noise hits)")
        print(f"    Channels affected : {len(ch_counts)} "
              f"out of {df_vmm_noise['ch'].nunique()}")
        print(f"    Hits/channel — "
              f"min={ch_counts.min()}  "
              f"max={ch_counts.max()}  "
              f"std={ch_counts.std():.1f}")

        low = df_vmm_noise[
            df_vmm_noise["adc"] < low_adc_threshold
        ]["adc"].values
        vals, counts = np.unique(low, return_counts=True)
        print(f"    ADC distribution below {low_adc_threshold}:")
        for v, c in sorted(zip(vals, counts),
                            key=lambda x: -x[1])[:5]:
            print(f"      ADC={v:3d} : {c:6d} hits")

        fig, ax = plt.subplots(figsize=(12, 3))
        ax.bar(ch_counts.index, ch_counts.values,
               color="red", alpha=0.7)
        ax.set_xlabel("Channel")
        ax.set_ylabel(f"ADC={artifact_adc} hit count")
        ax.set_title(
            f"VMM {vmm_id} Run {run_no} — "
            f"ADC={artifact_adc} artifact per channel"
        )
        plt.tight_layout()
        _finish_fig(fig, f"qa_adc16_vmm{vmm_id}_run{run_no}",
                    out_dir, show)


def qa_noise_sigma_distribution(df_run_scan, data_dir,
                                  sng1_runs,
                                  sigma_warn=10.0,
                                  sigma_bad=13.0, n_files=1,
                                  out_dir=None, show=True):
    """
    Plot distribution of robust_sigma across all VMMs and sng=1 runs.
    Shows where warn/bad thresholds sit relative to the data.
    """
    all_sigmas = []

    for run_no in sng1_runs:
        run_dir = get_run_dir(data_dir, run_no)
        df_hits = load_hits_run(
            run_dir, n_files=n_files,
            branches=["adc", "vmm", "over_threshold"]
        )
        if df_hits is None:
            continue
        df_noise = compute_noise_baseline(df_hits)
        snt_ser = df_run_scan.loc[df_run_scan["run_no"] == run_no, "snt"]
        sg_ser  = df_run_scan.loc[df_run_scan["run_no"] == run_no, "sg"]
        if snt_ser.empty or sg_ser.empty:
            print(f"Run {run_no} — not found in run scan table, skipping")
            continue
        snt = snt_ser.iloc[0]
        sg  = sg_ser.iloc[0]

        if "vmm_id" not in df_noise.columns:
            continue

        for _, row in df_noise.iterrows():
            all_sigmas.append({
                "run_no" : run_no,
                "sg"     : sg,
                "snt"    : snt,
                "vmm_id" : row["vmm_id"],
                "sigma"  : row["robust_sigma"],
                "quality": row["quality"]
            })

    df_sigmas = pd.DataFrame(all_sigmas)

    fig, ax = plt.subplots(figsize=(10, 4))
    for quality, color in [("ok", "green"),
                            ("warn", "orange"),
                            ("bad", "red")]:
        subset = df_sigmas[
            df_sigmas["quality"] == quality
        ]["sigma"]
        if not subset.empty:
            ax.hist(subset, bins=30, alpha=0.7,
                    label=quality, color=color)

    ax.axvline(sigma_warn, color="orange", linestyle="--",
               label=f"warn threshold={sigma_warn}")
    ax.axvline(sigma_bad, color="red", linestyle="--",
               label=f"bad threshold={sigma_bad}")
    ax.set_xlabel("Robust sigma (ADC)")
    ax.set_ylabel("Count")
    ax.set_title(
        "Distribution of noise sigma across all VMMs and runs"
    )
    ax.legend()
    plt.tight_layout()
    _finish_fig(fig, "qa_noise_sigma_dist", out_dir, show)

    print("\nBorderline cases (sigma within 1 ADC of warn threshold):")
    borderline = df_sigmas[
        (df_sigmas["sigma"] >= sigma_warn - 1.0) &
        (df_sigmas["sigma"] <= sigma_warn + 1.0)
    ]
    if borderline.empty:
        print("  None")
    else:
        print(borderline[["run_no", "sg", "snt",
                           "vmm_id", "sigma",
                           "quality"]].to_string(index=False))

    return df_sigmas


def qa_robust_vs_std_comparison(df_run_scan, data_dir, sng1_runs,
                                  n_files=1, out_dir=None, show=True):
    """
    Compare robust_sigma (MAD-based) vs standard deviation
    across all VMMs and sng=1 runs.
    """
    comparison = []

    for run_no in sng1_runs:
        run_dir = get_run_dir(data_dir, run_no)
        df_hits = load_hits_run(
            run_dir, n_files=n_files,
            branches=["adc", "vmm", "over_threshold"]
        )
        if df_hits is None:
            continue
        snt_ser = df_run_scan.loc[df_run_scan["run_no"] == run_no, "snt"]
        sg_ser  = df_run_scan.loc[df_run_scan["run_no"] == run_no, "sg"]
        if snt_ser.empty or sg_ser.empty:
            print(f"Run {run_no} — not found in run scan table, skipping")
            continue
        snt = snt_ser.iloc[0]
        sg  = sg_ser.iloc[0]

        for vmm_id in sorted(df_hits["vmm"].unique()):
            noise = df_hits[
                (df_hits["vmm"]            == vmm_id) &
                (df_hits["over_threshold"] == 0) &
                (df_hits["adc"]            >  20)
            ]["adc"].values

            if len(noise) < 100:
                continue

            median       = np.median(noise)
            mad          = np.median(np.abs(noise - median))
            robust_sigma = 1.4826 * mad
            std_sigma    = noise.std()

            comparison.append({
                "run_no"      : run_no,
                "sg"          : sg,
                "snt"         : snt,
                "vmm_id"      : vmm_id,
                "robust_sigma": robust_sigma,
                "std_sigma"   : std_sigma,
                "ratio"       : (std_sigma / robust_sigma
                                 if robust_sigma > 0 else np.nan)
            })

    df_comp = pd.DataFrame(comparison)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sc0 = axes[0].scatter(
        df_comp["robust_sigma"], df_comp["std_sigma"],
        c=df_comp["snt"], cmap="viridis", alpha=0.7
    )
    axes[0].plot(
        [0, df_comp["robust_sigma"].max()],
        [0, df_comp["robust_sigma"].max()],
        "r--", label="1:1 line"
    )
    axes[0].set_xlabel("Robust sigma (MAD)")
    axes[0].set_ylabel("Std sigma")
    axes[0].set_title("Robust σ vs Std σ per VMM per run")
    axes[0].legend()
    plt.colorbar(sc0, ax=axes[0], label="snt")

    sc1 = axes[1].scatter(
        df_comp["robust_sigma"], df_comp["ratio"],
        c=df_comp["snt"], cmap="viridis", alpha=0.7
    )
    axes[1].axhline(1.0, color="red", linestyle="--",
                    label="ratio=1")
    axes[1].set_xlabel("Robust sigma (MAD)")
    axes[1].set_ylabel("Std / Robust ratio")
    axes[1].set_title("How much larger is Std vs Robust σ?")
    axes[1].legend()
    plt.colorbar(sc1, ax=axes[1], label="snt")

    plt.tight_layout()
    _finish_fig(fig, "qa_robust_vs_std", out_dir, show)

    print("\nSummary of std/robust ratio:")
    print(f"  Mean  : {df_comp['ratio'].mean():.3f}")
    print(f"  Median: {df_comp['ratio'].median():.3f}")
    print(f"  Min   : {df_comp['ratio'].min():.3f}")
    print(f"  Max   : {df_comp['ratio'].max():.3f}")

    inflated = df_comp[df_comp["ratio"] > 1.5]
    print(f"\nCases where std > 1.5 * robust_sigma "
          f"({len(inflated)} of {len(df_comp)}):")
    if inflated.empty:
        print("  None")
    else:
        print(inflated[["run_no", "sg", "snt", "vmm_id",
                         "robust_sigma", "std_sigma",
                         "ratio"]].to_string(index=False))

    return df_comp


def qa_mpv_vs_median_comparison(df_run_scan, data_dir, pairs,
                                  detector_vmms,
                                  exclude_trigger_vmms=(0, 1),
                                  n_files=1,
                                  out_dir=None, show=True):
    """
    Compare MPV-based vs median-based SNR across all config pairs.
    """
    records = []

    for pair in pairs:
        sg       = pair["sg"]
        snt      = pair["snt"]
        run_sng0 = pair["sng0"]
        run_sng1 = pair["sng1"]

        run_dir        = get_run_dir(data_dir, run_sng0)
        df_hits_signal = load_hits_run(
            run_dir, n_files=n_files,
            branches=["adc", "vmm", "over_threshold"]
        )
        if df_hits_signal is None:
            continue

        run_dir_n     = get_run_dir(data_dir, run_sng1)
        df_hits_noise = load_hits_run(
            run_dir_n, n_files=n_files,
            branches=["adc", "vmm", "over_threshold"]
        )
        if df_hits_noise is None:
            continue
        df_noise = compute_noise_baseline(df_hits_noise)

        if "vmm_id" not in df_noise.columns:
            continue

        for vmm_id in detector_vmms:
            noise_row = df_noise[df_noise["vmm_id"] == vmm_id]
            if noise_row.empty:
                continue

            noise_sigma   = noise_row["robust_sigma"].iloc[0]
            noise_cut     = noise_row["noise_cut"].iloc[0]
            noise_quality = noise_row["quality"].iloc[0]

            signal = get_clean_signal(
                df_hits_signal, vmm_id, exclude_trigger_vmms
            )
            if signal is None or len(signal) < 100:
                continue

            mpv, _, _, _ = estimate_mpv(signal, adc_min=noise_cut)
            median_sig   = np.median(signal)
            mean_sig     = signal.mean()

            snr_mpv    = mpv        / noise_sigma
            snr_median = median_sig / noise_sigma
            snr_mean   = mean_sig   / noise_sigma

            records.append({
                "sg"           : sg,
                "snt"          : snt,
                "vmm_id"       : vmm_id,
                "noise_sigma"  : noise_sigma,
                "noise_quality": noise_quality,
                "mpv"          : mpv,
                "median_signal": median_sig,
                "mean_signal"  : mean_sig,
                "snr_mpv"      : snr_mpv,
                "snr_median"   : snr_median,
                "snr_mean"     : snr_mean,
                "ratio_med_mpv": snr_median / snr_mpv
            })

    df = pd.DataFrame(records)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    snt_values = sorted(df["snt"].unique())
    cmap       = plt.cm.viridis
    norm       = plt.Normalize(df["snt"].min(), df["snt"].max())

    for snt_val in snt_values:
        df_snt = df[df["snt"] == snt_val]
        color  = cmap(norm(snt_val))
        axes[0].scatter(df_snt["snr_mpv"], df_snt["snr_median"],
                        color=color, alpha=0.7,
                        label=f"snt={snt_val:.0f}")

    max_snr = max(df["snr_mpv"].max(),
                  df["snr_median"].max()) * 1.1
    axes[0].plot([0, max_snr], [0, max_snr],
                 "r--", label="1:1 line")
    axes[0].set_xlabel("SNR (MPV-based)")
    axes[0].set_ylabel("SNR (Median-based)")
    axes[0].set_title("MPV vs Median SNR\n"
                      "points above 1:1 = median overestimates")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3)

    for vmm_id in sorted(df["vmm_id"].unique()):
        df_vmm = df[df["vmm_id"] == vmm_id].sort_values("snt")
        if df_vmm.empty:
            continue
        axes[1].plot(df_vmm["snt"], df_vmm["ratio_med_mpv"],
                     "o-", label=f"VMM {vmm_id}", alpha=0.8)

    axes[1].axhline(1.0, color="red", linestyle="--",
                    label="ratio=1 (no bias)")
    axes[1].set_xlabel("Peaking time (snt)")
    axes[1].set_ylabel("SNR_median / SNR_mpv")
    axes[1].set_title("Median/MPV ratio vs peaking time\n"
                      "increasing ratio = growing bias")
    axes[1].legend(fontsize=7)
    axes[1].grid(True, alpha=0.3)

    for vmm_id in sorted(df["vmm_id"].unique()):
        df_vmm = df[df["vmm_id"] == vmm_id].sort_values("sg")
        if df_vmm.empty:
            continue
        axes[2].plot(df_vmm["sg"], df_vmm["ratio_med_mpv"],
                     "o-", label=f"VMM {vmm_id}", alpha=0.8)

    axes[2].axhline(1.0, color="red", linestyle="--",
                    label="ratio=1 (no bias)")
    axes[2].set_xlabel("Gain (sg)")
    axes[2].set_ylabel("SNR_median / SNR_mpv")
    axes[2].set_title("Median/MPV ratio vs gain\n"
                      "gain dependence of bias")
    axes[2].legend(fontsize=7)
    axes[2].grid(True, alpha=0.3)

    plt.suptitle(
        "MPV vs Median as signal estimator\n"
        "MPV is stable — median drifts into Landau tail"
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    _finish_fig(fig, "qa_mpv_vs_median", out_dir, show)

    print("\nSNR estimator comparison summary:")
    print(f"{'sg':>5} {'snt':>6} {'VMM':>5} {'MPV':>8} "
          f"{'Median':>8} {'Mean':>8} "
          f"{'SNR_mpv':>10} {'SNR_med':>10} {'ratio':>8}")
    for _, row in df.sort_values(
        ["sg", "snt", "vmm_id"]
    ).iterrows():
        print(f"{row['sg']:>5.1f} {row['snt']:>6.0f} "
              f"{row['vmm_id']:>5.0f} "
              f"{row['mpv']:>8.1f} {row['median_signal']:>8.1f} "
              f"{row['mean_signal']:>8.1f} "
              f"{row['snr_mpv']:>10.1f} "
              f"{row['snr_median']:>10.1f} "
              f"{row['ratio_med_mpv']:>8.2f}")

    return df


def qa_noise_run_diagnostic(data_dir, run_noise, run_signal,
                             detector_vmms, n_files=1, adc_low_cut=20,
                             out_dir=None, show=True):
    """
    Diagnose why a noise run fails to build a baseline.

    Use when compute_snr warns:
        "no noise baseline for run X (too few over_threshold=0 hits)"

    Prints a per-VMM table and plots both runs' ADC distributions
    split by over_threshold flag with the adc_low_cut marked.
    """
    run_dir_n = get_run_dir(data_dir, run_noise)
    df_n = load_hits_run(run_dir_n, n_files=n_files,
                          branches=["adc", "vmm", "over_threshold"])

    run_dir_s = get_run_dir(data_dir, run_signal)
    df_s = load_hits_run(run_dir_s, n_files=n_files,
                          branches=["adc", "vmm", "over_threshold"])

    print(f"\n{'='*74}")
    print(f"Noise run diagnostic — noise=run {run_noise}  "
          f"signal=run {run_signal}")
    print(f"{'='*74}")
    hdr = (f"{'VMM':>5}  {'N_total':>8}  {'ot=0':>7}  "
           f"{'ot=0 & adc>'+str(adc_low_cut):>14}  "
           f"{'ot=1':>7}  {'%ot=0':>6}  {'%ot=0 (cut)':>11}")
    print(hdr)
    print("-" * 74)

    for vmm_id in sorted(detector_vmms):
        if df_n is None:
            print(f"{vmm_id:>5}  {'no data':>8}")
            continue

        df_vmm = df_n[df_n["vmm"] == vmm_id]
        n_total = len(df_vmm)
        if n_total == 0:
            print(f"{vmm_id:>5}  {'0':>8}  {'—':>7}  {'—':>14}  "
                  f"{'—':>7}  {'—':>6}  {'—':>11}")
            continue

        n_ot0     = int((df_vmm["over_threshold"] == 0).sum())
        n_ot0_cut = int(((df_vmm["over_threshold"] == 0) &
                         (df_vmm["adc"] > adc_low_cut)).sum())
        n_ot1     = int((df_vmm["over_threshold"] == 1).sum())
        pct_ot0   = 100.0 * n_ot0     / n_total
        pct_ot0c  = 100.0 * n_ot0_cut / n_total

        flag = " ← TOO FEW" if n_ot0_cut < 100 else ""

        print(f"{vmm_id:>5}  {n_total:>8}  {n_ot0:>7}  {n_ot0_cut:>14}  "
              f"{n_ot1:>7}  {pct_ot0:>5.1f}%  {pct_ot0c:>10.1f}%{flag}")

    print()

    vmm_list   = sorted(detector_vmms)
    chunk_size = 4

    for chunk_i, start in enumerate(range(0, len(vmm_list), chunk_size)):
        chunk = vmm_list[start:start + chunk_size]
        fig, axes = plt.subplots(len(chunk), 2,
                                  figsize=(14, 4 * len(chunk)))
        if len(chunk) == 1:
            axes = [axes]

        for row_i, vmm_id in enumerate(chunk):
            for col_i, (df_hits, run_no, run_label) in enumerate([
                (df_n, run_noise,  "noise run (sng=1)"),
                (df_s, run_signal, "signal run (sng=0)"),
            ]):
                ax = axes[row_i][col_i]

                if df_hits is None:
                    ax.text(0.5, 0.5, "No data",
                            transform=ax.transAxes,
                            ha="center", va="center", fontsize=13)
                    ax.set_title(f"VMM {vmm_id} — {run_label} "
                                 f"(run {run_no})")
                    continue

                df_vmm = df_hits[df_hits["vmm"] == vmm_id]
                if df_vmm.empty:
                    ax.text(0.5, 0.5, "No hits",
                            transform=ax.transAxes,
                            ha="center", va="center", fontsize=13)
                    ax.set_title(f"VMM {vmm_id} — {run_label} "
                                 f"(run {run_no})")
                    continue

                for flag, color, lbl in [
                    (0, "steelblue", "ot=0 (noise)"),
                    (1, "tomato",    "ot=1 (signal)"),
                ]:
                    adc = df_vmm[df_vmm["over_threshold"] == flag]["adc"]
                    ax.hist(adc, bins=100, range=(0, 300),
                            alpha=0.6, color=color,
                            label=f"{lbl} (n={len(adc)})")

                ax.axvline(adc_low_cut, color="black",
                           linestyle="--", linewidth=1.5,
                           label=f"adc_low_cut={adc_low_cut}")
                ax.set_xlabel("ADC")
                ax.set_ylabel("Hits")
                ax.set_title(f"VMM {vmm_id} — {run_label} "
                             f"(run {run_no})")
                ax.legend(fontsize=9)

        plt.suptitle(
            f"Noise run diagnostic — noise=run {run_noise}  "
            f"signal=run {run_signal}",
            fontweight="bold"
        )
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        stem = (f"qa_noise_diagnostic_"
                f"run{run_noise}_vs_{run_signal}_{chunk_i}")
        _finish_fig(fig, stem, out_dir, show)


def main():
    pass


if __name__ == '__main__':
    main()