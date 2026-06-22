#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vmm_detector_efficiency.py

Detector performance analysis using trigger–detector time coincidences.

Answers:
  - What is the timing offset between trigger and detector signals?
  - What is the true detector efficiency per VMM, corrected for
    accidental coincidences, using only on-spill trigger timestamps.

Workflow
--------
1. load_sorted_hits()               — load one ROOT file
2. filter_on_spill_triggers()       — gate trigger timestamps to beam-on only
3. compute_trigger_detector_dt()    — build Δt histogram over ±1000 ns
4. plot_trigger_detector_dt()       — raw Δt: linear + log, see peak vs background
5. fit_dt_peak()                    — fit Gaussian + flat bg → derive μ, σ, window
6. plot_dt_peak_fits()              — verify fit: green=signal window, orange=sideband
7. compute_time_correlated_efficiency() — efficiency with fitted window + sideband subtraction
8. compute_channel_efficiency()     — per-pad efficiency (vmm_id, ch)
9. plot_efficiency_map()            — 2D pad map coloured by true efficiency

@author: ak271430
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import os

from vmm_mapping import vmm_mapping
from vmm_io import (load_sorted_hits, compute_spill_masks,
                    NS_PER_TICK, S_PER_TICK,
                    iter_hits_files, get_run_dir)


def _save_fig(fig, out_dir, stem):
    """Save figure as PNG and PDF into out_dir if provided."""
    if out_dir is None:
        return
    os.makedirs(out_dir, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out_dir, f"{stem}.{ext}"),
                    bbox_inches="tight")


def load_detector_map(map_csv, detector_key):
    """
    Load pad-position CSV and resolve (connector, channel) → (vmm_id, ch).

    The CSV must have columns: connector, channel, X_mm, Y_mm, Size_mm.
    connector is a 0-based index into vmm_mapping[detector_key]['vmm_ids'].

    Returns
    -------
    pd.DataFrame with columns: vmm_id, ch, x_mm, y_mm, size_mm
    """
    cfg    = vmm_mapping[detector_key]
    vmm_ids = cfg["vmm_ids"]

    df = pd.read_csv(map_csv)
    df["vmm_id"] = df["connector"].map(lambda c: vmm_ids[c])
    df["ch"]     = df["channel"]
    df = df.rename(columns={"X_mm": "x_mm", "Y_mm": "y_mm",
                             "Size_mm": "size_mm"})
    return df[["vmm_id", "ch", "x_mm", "y_mm", "size_mm"]].copy()


# ── Derived from vmm_mapping ─────────────────────────────────
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
trigger_ref_channels = {
    vmm_id: vmm_mapping["trigger"]["channels"][vmm_id][0]
    for vmm_id in vmm_mapping["trigger"]["vmm_ids"]
}


# ─────────────────────────────────────────────────────────────
# SPILL FILTER
# ─────────────────────────────────────────────────────────────

def filter_on_spill_triggers(df_hits, trigger_vmm, trigger_ch,
                              spill_threshold_khz=1.0,
                              bin_width_ms=1.0,
                              min_spill_s=1.0,
                              max_gap_s=2.0):
    """
    Return sorted, deduplicated on-spill trigger timestamps.

    Bins the trigger hit rate into 1 ms bins, classifies bins as
    on/off spill via compute_spill_masks, then keeps only trigger
    timestamps that fall in on-spill bins.

    Parameters
    ----------
    spill_threshold_khz : float
        Rate threshold separating beam-on from beam-off bins.
    bin_width_ms : float
        Bin width in ms for rate computation (default 1 ms).
    min_spill_s : float
        Minimum on-region duration to retain (removes startup artifacts).
    max_gap_s : float
        Maximum off-gap to bridge within a spill (handles SPS bunch
        microstructure dips at low beam intensity). Default 2 s.

    Returns
    -------
    np.ndarray
        Sorted, deduplicated trigger timestamps (ticks) in on-spill bins.
    """
    t_trig_all = df_hits[
        (df_hits["vmm"] == trigger_vmm) &
        (df_hits["ch"]  == trigger_ch)
    ]["time"].values

    if len(t_trig_all) < 2:
        print(f"  Too few trigger hits on VMM {trigger_vmm} ch {trigger_ch}")
        return np.array([])

    bin_width_s   = bin_width_ms * 1e-3
    t_s           = t_trig_all * S_PER_TICK
    t_start       = t_s.min()
    edges         = np.arange(t_start, t_s.max() + bin_width_s, bin_width_s)
    counts, edges = np.histogram(t_s, bins=edges)
    rate_khz      = counts / bin_width_s / 1e3

    min_spill_bins = max(1, int(min_spill_s / bin_width_s))
    max_gap_bins   = max(0, int(max_gap_s   / bin_width_s))
    on_mask, _     = compute_spill_masks(rate_khz, spill_threshold_khz,
                                          min_spill_bins=min_spill_bins,
                                          max_gap_bins=max_gap_bins)

    bin_idx   = np.clip(np.digitize(t_s, edges) - 1, 0, len(on_mask) - 1)
    t_trig_on = t_trig_all[on_mask[bin_idx]]

    n_total = len(t_trig_all)
    n_on    = len(t_trig_on)
    print(f"  Spill filter: {n_on:,}/{n_total:,} trigger hits on-spill "
          f"({100 * n_on / n_total:.1f}%)")

    return np.sort(np.unique(t_trig_on))


# ─────────────────────────────────────────────────────────────
# TIMING DIAGNOSTIC
# ─────────────────────────────────────────────────────────────

def compute_trigger_detector_dt(df_hits, t_trig,
                                  search_window_ns=1000):
    """
    For each trigger timestamp, find all detector hits within
    ±search_window_ns and record Δt = t_det - t_trig (in ns).

    Uses np.searchsorted on a sorted detector time array —
    efficient for large datasets (O(N_trig × log N_det)).

    Parameters
    ----------
    t_trig : np.ndarray
        Trigger timestamps in ticks (sorted, deduplicated).
        Obtain from filter_on_spill_triggers().
    search_window_ns : float
        Half-width of the search window in ns.
        Default 1000 ns = ±1 μs, wide enough to see any offset.

    Returns
    -------
    dict : vmm_id (int) → np.ndarray of Δt values in ns
    """
    if len(t_trig) == 0:
        print("  No trigger timestamps provided")
        return {}

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
                               n_bins=200, out_dir=None):
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
        _save_fig(fig, out_dir, f"dt_raw_{group_label.replace(' ', '_')}")


# ─────────────────────────────────────────────────────────────
# PEAK FIT — coincidence window optimisation
# ─────────────────────────────────────────────────────────────

def _gaussian_plus_bg(x, amplitude, mu, sigma, bg):
    return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + bg


def fit_dt_peak(dt_dict, search_window_ns=1000, n_bins=200, n_sigma=3):
    """
    Fit Gaussian + flat background to each VMM's Δt distribution.

    Model: f(Δt) = A · exp(-(Δt - μ)² / 2σ²) + B

    The coincidence window is set to [μ - n_sigma·σ, μ + n_sigma·σ].
    The sideband is placed immediately after the signal window with
    the same width: [μ + n_sigma·σ + gap, μ + n_sigma·σ + gap + width],
    where gap = σ and width = 2·n_sigma·σ.

    Parameters
    ----------
    dt_dict : dict
        Output of compute_trigger_detector_dt(): vmm_id → Δt array (ns).
    search_window_ns : float
        Half-width of the Δt histogram range (must match the value used
        when calling compute_trigger_detector_dt).
    n_bins : int
        Number of histogram bins.
    n_sigma : float
        Half-width of the coincidence window in units of σ (default 3).

    Returns
    -------
    dict : vmm_id → {
        "mu"         : float  — peak centre (ns)
        "sigma"      : float  — peak width (ns)
        "amplitude"  : float  — Gaussian amplitude (counts/bin)
        "bg"         : float  — flat background (counts/bin)
        "dt_min"     : float  — window lower edge = μ - n_sigma·σ (ns)
        "dt_max"     : float  — window upper edge = μ + n_sigma·σ (ns)
        "sb_min"     : float  — sideband lower edge (ns)
        "sb_max"     : float  — sideband upper edge (ns)
        "bin_centers": array  — bin centres used for the fit
        "hist"       : array  — bin counts
        "success"    : bool
    }
    """
    bin_edges   = np.linspace(-search_window_ns, search_window_ns, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    results = {}

    for vmm_id, dt in dt_dict.items():
        if len(dt) < 20:
            print(f"  VMM {vmm_id}: too few Δt pairs ({len(dt)}) to fit")
            results[vmm_id] = {"success": False}
            continue

        hist, _ = np.histogram(dt, bins=bin_edges)
        hist     = hist.astype(float)

        # Initial estimates: background = low-percentile level,
        # peak centre = bin with max counts, peak width = 50 ns
        bg_init    = float(np.percentile(hist, 20))
        A_init     = float(hist.max() - bg_init)
        mu_init    = float(bin_centers[np.argmax(hist)])
        sigma_init = 50.0

        p0     = [A_init, mu_init, sigma_init, bg_init]
        bounds = (
            [0,    -search_window_ns / 2,  5,   0],
            [np.inf, search_window_ns / 2, 500, np.inf]
        )

        try:
            popt, pcov = curve_fit(_gaussian_plus_bg, bin_centers, hist,
                                   p0=p0, bounds=bounds, maxfev=5000)
            A_fit, mu_fit, sigma_fit, bg_fit = popt
            perr = np.sqrt(np.diag(pcov))

            win_half  = n_sigma * sigma_fit
            dt_min    = mu_fit - win_half
            dt_max    = mu_fit + win_half
            gap       = sigma_fit
            sb_min    = dt_max + gap
            sb_max    = sb_min + 2 * win_half

            print(f"  VMM {vmm_id}: "
                  f"μ={mu_fit:+.1f}±{perr[1]:.1f} ns  "
                  f"σ={sigma_fit:.1f}±{perr[2]:.1f} ns  "
                  f"signal=[{dt_min:.0f}, {dt_max:.0f}] ns  "
                  f"sideband=[{sb_min:.0f}, {sb_max:.0f}] ns")

            results[vmm_id] = {
                "mu"         : mu_fit,
                "sigma"      : sigma_fit,
                "amplitude"  : A_fit,
                "bg"         : bg_fit,
                "dt_min"     : dt_min,
                "dt_max"     : dt_max,
                "sb_min"     : sb_min,
                "sb_max"     : sb_max,
                "bin_centers": bin_centers,
                "hist"       : hist,
                "success"    : True,
            }

        except RuntimeError as exc:
            print(f"  VMM {vmm_id}: fit failed — {exc}")
            results[vmm_id] = {"success": False}

    return results


def plot_dt_peak_fits(peak_fits, run_no, vmm_groups, out_dir=None):
    """
    Plot the Δt histogram with the Gaussian+background fit overlaid.

    One figure per detector group, one subplot per VMM. Each subplot shows:
    - Blue bars  : Δt histogram (on-spill triggers only)
    - Red curve  : fitted Gaussian + flat background
    - Green band : signal coincidence window [μ - n_sigma·σ, μ + n_sigma·σ]
    - Orange band: sideband window (equal-width, for accidental subtraction)
    - Dashed lines: μ (black), window edges (green), sideband edges (orange)

    Use this immediately after fit_dt_peak() to visually verify that:
    1. The Gaussian describes the peak well (red curve follows blue bars).
    2. The signal window captures the peak and nothing more.
    3. The sideband sits in the flat accidental region.
    """
    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    for group in vmm_groups:
        vmm_ids     = [v for v in group["vmm_ids"] if v in peak_fits]
        group_label = group["label"]

        if not vmm_ids:
            continue

        n_vmm = len(vmm_ids)
        fig, axes = plt.subplots(1, n_vmm,
                                  figsize=(5 * n_vmm, 5),
                                  sharey=False)
        if n_vmm == 1:
            axes = [axes]

        for ax, vmm_id in zip(axes, vmm_ids):
            fit = peak_fits[vmm_id]
            det_name = vmm_to_detector.get(vmm_id, "")

            if not fit["success"]:
                ax.text(0.5, 0.5, "Fit failed",
                        ha="center", va="center",
                        transform=ax.transAxes, color="red")
                ax.set_title(f"VMM {vmm_id} — {det_name}")
                continue

            bc   = fit["bin_centers"]
            hist = fit["hist"]
            x_fine = np.linspace(bc[0], bc[-1], 500)
            y_fit  = _gaussian_plus_bg(x_fine,
                                        fit["amplitude"], fit["mu"],
                                        fit["sigma"],    fit["bg"])

            bw = bc[1] - bc[0]
            ax.bar(bc, hist, width=bw, color="steelblue",
                   alpha=0.6, label="Δt data")
            ax.plot(x_fine, y_fit, color="red",
                    linewidth=2, label="Gaussian + bg fit")

            # Signal window
            ax.axvspan(fit["dt_min"], fit["dt_max"],
                       color="green", alpha=0.15, label="signal window")
            for edge in (fit["dt_min"], fit["dt_max"]):
                ax.axvline(edge, color="green",
                           linestyle="--", linewidth=1.2)

            # Sideband window
            ax.axvspan(fit["sb_min"], fit["sb_max"],
                       color="orange", alpha=0.15, label="sideband")
            for edge in (fit["sb_min"], fit["sb_max"]):
                ax.axvline(edge, color="orange",
                           linestyle="--", linewidth=1.2)

            # Peak centre
            ax.axvline(fit["mu"], color="black",
                       linestyle=":", linewidth=1.5, label=f"μ={fit['mu']:+.1f} ns")

            info = (f"μ = {fit['mu']:+.1f} ns\n"
                    f"σ = {fit['sigma']:.1f} ns\n"
                    f"signal: [{fit['dt_min']:.0f}, {fit['dt_max']:.0f}] ns\n"
                    f"sideband: [{fit['sb_min']:.0f}, {fit['sb_max']:.0f}] ns")
            ax.text(0.97, 0.97, info,
                    ha="right", va="top",
                    transform=ax.transAxes,
                    fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.3",
                              facecolor="white", alpha=0.8))

            ax.set_xlabel("Δt (ns)")
            ax.set_ylabel("Counts")
            ax.set_title(f"VMM {vmm_id} — {det_name}")
            ax.legend(fontsize=8, loc="upper left")
            ax.grid(True, alpha=0.3)

        fig.suptitle(
            f"{group_label} — Δt peak fit  |  Run {run_no}",
            fontsize=13, fontweight="bold"
        )
        plt.tight_layout()
        _save_fig(fig, out_dir, f"dt_fit_{group_label.replace(' ', '_')}")


# ─────────────────────────────────────────────────────────────
# COINCIDENCE EFFICIENCY
# ─────────────────────────────────────────────────────────────

def compute_time_correlated_efficiency(df_hits, t_trig,
                                        dt_min_ns=0,
                                        dt_max_ns=250,
                                        sideband_min_ns=400,
                                        sideband_max_ns=650):
    """
    Compute true coincidence efficiency per detector VMM.

    For each on-spill trigger timestamp, checks whether at least one
    detector hit falls in the signal window [dt_min_ns, dt_max_ns].
    Uses a sideband window of equal width to estimate and subtract
    the accidental coincidence rate.

    Definition
    ----------
    raw_efficiency  = n_triggers_with_≥1_signal_hit / n_triggers
    accidental_eff  = n_triggers_with_≥1_sideband_hit / n_triggers
    true_efficiency = raw_efficiency - accidental_eff

    Parameters
    ----------
    t_trig : np.ndarray
        On-spill trigger timestamps in ticks (sorted, deduplicated).
        Obtain from filter_on_spill_triggers().
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
    if len(t_trig) == 0:
        print("  No trigger timestamps provided")
        return pd.DataFrame()

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


def plot_time_correlated_efficiency(df_coinc, vmm_groups, out_dir=None):
    """
    Plot true coincidence efficiency vs configuration per VMM,
    grouped by detector.

    One figure per detector group, one line per VMM.
    X axis: configuration (sg / snt), sorted.
    Y axis: true_efficiency with binomial error bars.

    df_coinc must have columns: sg, snt, vmm_id, true_efficiency, n_triggers.
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
        ax.set_title(f"{group_label} — time-correlated efficiency")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)
        plt.tight_layout()
        _save_fig(fig, out_dir, f"efficiency_{group_label.replace(' ', '_')}")


# ─────────────────────────────────────────────────────────────
# PAD EFFICIENCY MAP
# ─────────────────────────────────────────────────────────────

def compute_channel_efficiency(df_hits, t_trig,
                                dt_min_ns, dt_max_ns,
                                sb_min_ns, sb_max_ns):
    """
    Compute true coincidence efficiency per detector channel (vmm_id, ch).

    For each channel, counts how many on-spill triggers had at least one
    hit in the signal window [dt_min_ns, dt_max_ns], then subtracts the
    accidental rate estimated from the sideband [sb_min_ns, sb_max_ns].

    Parameters
    ----------
    df_hits : pd.DataFrame
        All hits (time, vmm, ch columns required), sorted by time.
    t_trig : np.ndarray
        On-spill trigger timestamps (ns), sorted.
    dt_min_ns, dt_max_ns : float
        Signal coincidence window relative to trigger (ns).
    sb_min_ns, sb_max_ns : float
        Sideband window for accidental estimation (ns).

    Returns
    -------
    pd.DataFrame
        One row per (vmm_id, ch) with columns:
        vmm_id, ch, n_triggers, n_signal, n_sideband,
        raw_efficiency, accidental_eff, true_efficiency
    """
    n_trig  = len(t_trig)
    t_hit   = df_hits["time"].values
    vmm_arr = df_hits["vmm"].values
    ch_arr  = df_hits["ch"].values
    results = []

    for vmm_id in np.unique(vmm_arr):
        vmm_mask = vmm_arr == vmm_id
        for ch in np.unique(ch_arr[vmm_mask]):
            t_ch = np.sort(t_hit[vmm_mask & (ch_arr == ch)])

            n_signal    = 0
            n_sideband  = 0
            for t in t_trig:
                i0 = np.searchsorted(t_ch, t + dt_min_ns)
                i1 = np.searchsorted(t_ch, t + dt_max_ns, side="right")
                if i1 > i0:
                    n_signal += 1

                j0 = np.searchsorted(t_ch, t + sb_min_ns)
                j1 = np.searchsorted(t_ch, t + sb_max_ns, side="right")
                if j1 > j0:
                    n_sideband += 1

            raw_eff = n_signal   / n_trig if n_trig > 0 else 0.0
            acc_eff = n_sideband / n_trig if n_trig > 0 else 0.0
            results.append({
                "vmm_id"         : vmm_id,
                "ch"             : ch,
                "n_triggers"     : n_trig,
                "n_signal"       : n_signal,
                "n_sideband"     : n_sideband,
                "raw_efficiency" : raw_eff,
                "accidental_eff" : acc_eff,
                "true_efficiency": raw_eff - acc_eff,
            })

    return pd.DataFrame(results)


def plot_efficiency_map(df_ch_eff, det_map, run_no,
                        detector_label="", out_dir=None):
    """
    Plot per-pad true coincidence efficiency as a 2D colour map.

    Each pad is drawn as a filled square at its physical (x_mm, y_mm)
    position, coloured by true_efficiency on a red–yellow–green scale.
    Pads with no data (not in df_ch_eff) are drawn in grey.

    Parameters
    ----------
    df_ch_eff : pd.DataFrame
        Output of compute_channel_efficiency(): must have vmm_id, ch,
        true_efficiency.
    det_map : pd.DataFrame
        Output of load_detector_map(): must have vmm_id, ch, x_mm,
        y_mm, size_mm.
    run_no : int
        Run number, used in the figure title and save filename.
    detector_label : str
        Human-readable detector name for title and filename.
    out_dir : str or None
        Directory for saving PNG/PDF. None → no save.
    """
    df = det_map.merge(
        df_ch_eff[["vmm_id", "ch", "true_efficiency"]],
        on=["vmm_id", "ch"],
        how="left",
    )

    cmap = plt.cm.RdYlGn
    norm = plt.Normalize(vmin=0, vmax=1)

    fig, ax = plt.subplots(figsize=(10, 8))

    for _, row in df.iterrows():
        half  = row["size_mm"] / 2.0
        eff   = row["true_efficiency"]
        color = cmap(norm(eff)) if pd.notna(eff) else "lightgrey"
        rect  = plt.Rectangle(
            (row["x_mm"] - half, row["y_mm"] - half),
            row["size_mm"], row["size_mm"],
            facecolor=color,
            edgecolor="black", linewidth=0.4,
        )
        ax.add_patch(rect)

        label = f"{eff:.2f}" if pd.notna(eff) else "—"
        ax.text(row["x_mm"], row["y_mm"], label,
                ha="center", va="center", fontsize=6, color="black")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="True coincidence efficiency", shrink=0.8)

    margin = det_map["size_mm"].max()
    ax.set_xlim(det_map["x_mm"].min() - margin,
                det_map["x_mm"].max() + margin)
    ax.set_ylim(det_map["y_mm"].min() - margin,
                det_map["y_mm"].max() + margin)
    ax.set_aspect("equal")
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title(
        f"{detector_label} — pad efficiency map  |  Run {run_no}",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    _save_fig(fig, out_dir,
              f"eff_map_{detector_label.replace(' ', '_')}_run{run_no}")


# ─────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────

def main():
    # ── User configuration ──────────────────────────────────
    # cnfg_dir = "/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/VMM-alinx-data/"
    # data_dir = "/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/VMM-alinx-data/5kHz-muons-config-scan/"

    cnfg_dir = "/drf/projets/clas12/P2/akallits/"
    data_dir = "/drf/projets/clas12/cern_202511_p2_alinx/"
    map_dir  = ("/local/home/ak271430/Documents/PostDocSaclay/"
                "data/det_mappings/")

    small_det_map_csv = os.path.join(map_dir, "p2_small_detector_map.csv")
    small_detectors   = ["p2_small_1", "p2_small_3"]

    # n_files: number of files for efficiency accumulation (Pass 2).
    # Memory safety: one file loaded at a time; only integer counts kept
    # between iterations. Peak memory = one ROOT file regardless of n_files.
    n_files       = 2
    diag_file_idx = 1       # file used for Δt diagnostic + peak fit (Pass 1)
    test_run      = 67      # sg=3, snt=200

    out_dir = os.path.join(cnfg_dir, "results", f"run_{test_run}")
    os.makedirs(out_dir, exist_ok=True)
    print(f"Results → {out_dir}")

    trigger_vmm         = 0
    trigger_ch          = trigger_ref_channels[trigger_vmm]

    spill_threshold_khz = 1.0
    bin_width_ms        = 1.0
    min_spill_s         = 1.0
    max_gap_s           = 2.0

    dt_min_ns       = 0
    dt_max_ns       = 250
    sideband_min_ns = 400
    sideband_max_ns = 650

    # ── Pass 1: diagnostic from one file ────────────────────
    # Load a single file to fit the coincidence peak and derive the
    # signal/sideband windows used in Pass 2.
    print(f"\nPass 1 — diagnostic (file index {diag_file_idx}, run {test_run})...")
    df_hits_diag = load_sorted_hits(
        data_dir, test_run,
        branches=["time", "vmm", "ch"],
        root_file_index=diag_file_idx,
    )
    if df_hits_diag is None or df_hits_diag.empty:
        print("No data loaded — check data_dir and test_run.")
        return

    t_trig_diag = filter_on_spill_triggers(
        df_hits_diag,
        trigger_vmm=trigger_vmm,
        trigger_ch=trigger_ch,
        spill_threshold_khz=spill_threshold_khz,
        bin_width_ms=bin_width_ms,
        min_spill_s=min_spill_s,
        max_gap_s=max_gap_s,
    )
    print(f"  {len(t_trig_diag):,} on-spill triggers")

    # Step 1: raw Δt histogram
    print(f"\nStep 1 — trigger–detector Δt (±1000 ns)...")
    dt_dict = compute_trigger_detector_dt(
        df_hits_diag, t_trig_diag, search_window_ns=1000
    )
    for vmm_id, dt in dt_dict.items():
        print(f"  VMM {vmm_id}: {len(dt):,} pairs in ±1000 ns window")
    plot_trigger_detector_dt(
        dt_dict, test_run, vmm_groups,
        search_window_ns=1000, out_dir=out_dir,
    )

    # Step 2: fit peak → derive coincidence window
    print(f"\nStep 2 — fitting Δt peak (Gaussian + flat background)...")
    peak_fits = fit_dt_peak(dt_dict, search_window_ns=1000,
                             n_bins=200, n_sigma=3)
    plot_dt_peak_fits(peak_fits, run_no=test_run, vmm_groups=vmm_groups,
                      out_dir=out_dir)

    del df_hits_diag, dt_dict, t_trig_diag   # free diagnostic data

    fitted = {v: r for v, r in peak_fits.items() if r.get("success")}
    if fitted:
        first = next(iter(fitted.values()))
        dt_min_ns       = first["dt_min"]
        dt_max_ns       = first["dt_max"]
        sideband_min_ns = first["sb_min"]
        sideband_max_ns = first["sb_max"]
        print(f"\n  Fitted window: "
              f"signal=[{dt_min_ns:.0f}, {dt_max_ns:.0f}] ns  "
              f"sideband=[{sideband_min_ns:.0f}, {sideband_max_ns:.0f}] ns")
    else:
        print(f"\n  Fits failed — using defaults: "
              f"signal=[{dt_min_ns}, {dt_max_ns}] ns  "
              f"sideband=[{sideband_min_ns}, {sideband_max_ns}] ns")

    # ── Pass 2: accumulate counts over n_files ───────────────
    # Files are streamed one at a time. After processing each file,
    # df_hits is deleted and only small count DataFrames are retained.
    # Efficiencies are recomputed from accumulated totals after the loop.
    print(f"\nPass 2 — accumulating over {n_files} file(s) "
          f"(starting at file index {diag_file_idx})...")
    run_dir = get_run_dir(data_dir, test_run)

    vmm_dfs      = []   # list of small DataFrames: [vmm_id, n_matched, n_accidental]
    ch_dfs       = []   # list of small DataFrames: [vmm_id, ch, n_signal, n_sideband]
    n_trig_total = 0

    for i, df_hits in enumerate(iter_hits_files(
            run_dir, n_files=n_files,
            branches=["time", "vmm", "ch"],
            file_start=diag_file_idx)):

        df_hits = (df_hits[df_hits["time"] < 2e12]
                   .sort_values("time")
                   .reset_index(drop=True))
        if df_hits.empty:
            del df_hits
            continue

        t_trig = filter_on_spill_triggers(
            df_hits,
            trigger_vmm=trigger_vmm,
            trigger_ch=trigger_ch,
            spill_threshold_khz=spill_threshold_khz,
            bin_width_ms=bin_width_ms,
            min_spill_s=min_spill_s,
            max_gap_s=max_gap_s,
        )
        n_trig = len(t_trig)
        print(f"  file {diag_file_idx + i}: {n_trig:,} on-spill triggers  "
              f"({len(df_hits):,} hits)")

        if n_trig == 0:
            del df_hits
            continue
        n_trig_total += n_trig

        # per-VMM counts (Step 3)
        df_vmm = compute_time_correlated_efficiency(
            df_hits, t_trig,
            dt_min_ns=dt_min_ns, dt_max_ns=dt_max_ns,
            sideband_min_ns=sideband_min_ns, sideband_max_ns=sideband_max_ns,
        )
        vmm_dfs.append(df_vmm[["vmm_id", "n_matched", "n_accidental"]].copy())

        # per-channel counts (Step 4)
        df_ch = compute_channel_efficiency(
            df_hits, t_trig,
            dt_min_ns=dt_min_ns, dt_max_ns=dt_max_ns,
            sb_min_ns=sideband_min_ns, sb_max_ns=sideband_max_ns,
        )
        ch_dfs.append(df_ch[["vmm_id", "ch", "n_signal", "n_sideband"]].copy())

        del df_hits, df_vmm, df_ch   # one file at a time

    print(f"\n  Total: {n_trig_total:,} on-spill triggers "
          f"across {len(vmm_dfs)} file(s)")

    if n_trig_total == 0 or not vmm_dfs:
        print("No on-spill triggers found — check spill threshold.")
        plt.show()
        return

    # ── Step 3: per-VMM summary (console only) ─────────────
    df_eff = (pd.concat(vmm_dfs)
                .groupby("vmm_id")[["n_matched", "n_accidental"]]
                .sum()
                .reset_index())
    df_eff["n_triggers"]      = n_trig_total
    df_eff["raw_efficiency"]  = df_eff["n_matched"]    / n_trig_total
    df_eff["accidental_eff"]  = df_eff["n_accidental"] / n_trig_total
    df_eff["true_efficiency"] = df_eff["raw_efficiency"] - df_eff["accidental_eff"]
    del vmm_dfs

    print("\n=== Time-correlated efficiency (per VMM) ===")
    print(df_eff[["vmm_id", "n_triggers", "n_matched", "n_accidental",
                   "raw_efficiency", "accidental_eff",
                   "true_efficiency"]].to_string(index=False))

    # ── Step 4: per-channel efficiency + map ────────────────
    df_ch_eff = (pd.concat(ch_dfs)
                   .groupby(["vmm_id", "ch"])[["n_signal", "n_sideband"]]
                   .sum()
                   .reset_index())
    df_ch_eff["n_triggers"]      = n_trig_total
    df_ch_eff["raw_efficiency"]  = df_ch_eff["n_signal"]   / n_trig_total
    df_ch_eff["accidental_eff"]  = df_ch_eff["n_sideband"] / n_trig_total
    df_ch_eff["true_efficiency"] = (df_ch_eff["raw_efficiency"]
                                    - df_ch_eff["accidental_eff"])
    del ch_dfs

    out_csv = os.path.join(out_dir, f"efficiency_per_channel_run{test_run}.csv")
    df_ch_eff.to_csv(out_csv, index=False)
    print(f"Saved → {out_csv}")

    print(f"\nStep 4 — per-pad efficiency map...")
    for det_key in small_detectors:
        det_label = vmm_mapping[det_key].get("name", det_key)
        det_vmms  = vmm_mapping[det_key]["vmm_ids"]
        df_det    = df_ch_eff[df_ch_eff["vmm_id"].isin(det_vmms)]
        if df_det.empty:
            print(f"  {det_label}: no hits — skipping map")
            continue
        det_map = load_detector_map(small_det_map_csv, det_key)
        print(f"  {det_label}: {len(df_det)} channels with hits")
        plot_efficiency_map(df_det, det_map, test_run,
                            detector_label=det_label, out_dir=out_dir)

    plt.show()


if __name__ == "__main__":
    main()

print("donzo")
