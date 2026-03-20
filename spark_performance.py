#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 14/01/2026 15:18
Created in PyCharm
Created as spark_performance.py

@author: akallits
"""

import os
from dataclasses import dataclass
from typing import Optional
from scipy.signal import find_peaks
from scipy.ndimage import uniform_filter1d


import uproot
# import ROOT
import numpy as np
import matplotlib.pyplot as plt


# ----------------------------
# Configuration
# ----------------------------
@dataclass
class Config:
    # Core analysis settings
    data_dir: str
    max_files: Optional[int] = None   # e.g. 3
    file_stride: int = 1              # e.g. 5 = take every 5th file
    bin_width: float = 0.01            # seconds for lab tests
    # bin_width: float = 0.001            # seconds for SPS Beam Tests
    # t_cut: float = 20.0                # seconds (ignore first t_cut) for lab tests
    t_cut: float = 0.0                # seconds (ignore first t_cut) for SPS Beam Tests
    fit_window_sigmas: float = 2.5     # window = mode ± fit_window_sigmas * sqrt(mode)

    # Mode control
    mode: str = "analysis"             # "debug" or "analysis"

    # Plot control
    plot_per_file: bool = True        # per-file plots (rate + hits:time)
    debug_fit: bool = True             # show gaussian-fit debug plot (window, points)
    plot_global: bool = False           # global plots (all-files rate distro + spark duration hist)

    # Spark definition behavior
    spark_threshold_type: str = "poisson"  # "poisson" or "rms"
    recovery_condition: str = "median"     # "median" or "threshold" (threshold = same as spark threshold)


# ----------------------------
# I/O + style
# ----------------------------
def load_root_file(file_path, branches=None):
    """Load selected branches from ROOT file into pandas DataFrame-like structure."""
    with uproot.open(file_path) as file:
        tree = file["hits"]
        return tree.arrays(branches, library="pd")


def set_root_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 14,
        "axes.linewidth": 1.5,
        "axes.labelsize": 16,
        "axes.titlesize": 16,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 7,
        "ytick.major.size": 7,
        "xtick.major.width": 1.3,
        "ytick.major.width": 1.3,
        "xtick.minor.size": 4,
        "ytick.minor.size": 4,
        "xtick.minor.width": 1.0,
        "ytick.minor.width": 1.0,
        "legend.frameon": False,
        "figure.figsize": (7, 5),
    })


# ----------------------------
# Shared preprocessing
# ----------------------------
def prepare_time_axis(hit_times_ns: np.ndarray, t_cut: float) -> np.ndarray:
    """
    Convert ns -> s and apply the initial time cut.
    Returns time_rel (seconds) after cut, relative to first timestamp.
    """
    hit_times_ns = np.asarray(hit_times_ns, dtype=float)
    if hit_times_ns.size < 2:
        return np.array([], dtype=float)

    hit_times_s = hit_times_ns / 1e9
    t0 = hit_times_s[0]
    time_rel = hit_times_s - t0
    time_rel = time_rel[time_rel > t_cut]
    return time_rel.astype(float)


# ----------------------------
# Gaussian main-peak fit (occupancy-value distribution)
# ----------------------------
def fit_main_peak_gaussian(occupancies, window_sigmas=2.5, min_points=6, debug=False):
    """
    Fit a Gaussian to the MAIN peak of the occupancy-value distribution.

    occupancies: array of counts-per-timebin (hist of time)
    window: [mode ± window_sigmas * sqrt(mode)] in occupancy units

    Returns:
      - if debug=False: (mu_fit, sigma_fit, A_fit, fit_ok)
      - if debug=True : (mu_fit, sigma_fit, A_fit, fit_ok, dbg_dict)
    """
    occ = np.asarray(occupancies, dtype=int)
    occ = occ[np.isfinite(occ)]
    occ = occ[occ >= 0]
    if occ.size < 20:
        out = (np.nan, np.nan, np.nan, False)
        return (*out, {}) if debug else out

    max_occ = int(occ.max())
    y_full = np.bincount(occ, minlength=max_occ + 1)
    x_full = np.arange(len(y_full))

    nonzero = y_full > 0
    x = x_full[nonzero]
    y = y_full[nonzero]
    if y.size < min_points:
        out = (np.nan, np.nan, np.nan, False)
        dbg = {"x_full": x_full, "y_full": y_full}
        return (*out, dbg) if debug else out

    mode = int(x[np.argmax(y)])

    #test for the SPS Beam Test Data
    # min_occ = 5  # still exclude noise floor
    # candidate_mask = x > min_occ
    # xc, yc = x[candidate_mask], y[candidate_mask]
    #
    # peaks, props = find_peaks(yc, height=yc.max() * 0.02, distance=50)
    #
    # if len(peaks) >= 2:
    #     # Sort by height, drop the tallest (spill peak), take the next tallest
    #     sorted_by_height = peaks[np.argsort(props["peak_heights"])[::-1]]
    #     signal_peak = sorted_by_height[1]  # second tallest = inter-spill baseline
    #     mode = int(xc[signal_peak])
    # elif len(peaks) == 1:
    #     mode = int(xc[peaks[0]])  # only one peak found, use it
    # else:
    #     mode = int(xc[np.argmax(yc)])  # fallback

    sigma_guess = float(np.sqrt(max(mode, 1)))  # Poisson-scale guess

    x_lo = max(0, int(np.floor(mode - window_sigmas * sigma_guess)))
    x_hi = int(np.ceil(mode + window_sigmas * sigma_guess))

    sel = (x >= x_lo) & (x <= x_hi)
    xw = x[sel].astype(float)
    yw = y[sel].astype(float)

    dbg = {
        "x_full": x_full, "y_full": y_full,
        "mode": mode, "sigma_guess": sigma_guess,
        "x_lo": x_lo, "x_hi": x_hi,
        "xw": xw, "yw": yw
    }

    if yw.size < min_points:
        out = (np.nan, np.nan, np.nan, False)
        return (*out, dbg) if debug else out

    try:
        from scipy.optimize import curve_fit

        def gaus(xx, A, mu, sig):
            return A * np.exp(-0.5 * ((xx - mu) / sig) ** 2)

        p0 = [float(yw.max()), float(mode), max(1.0, sigma_guess)]
        bounds = ([0.0, mode - 5 * sigma_guess, 0.2],
                  [np.inf, mode + 5 * sigma_guess, np.inf])

        popt, _ = curve_fit(gaus, xw, yw, p0=p0, bounds=bounds, maxfev=20000)
        A_fit, mu_fit, sig_fit = popt

        fit_ok = bool(np.isfinite(mu_fit) and np.isfinite(sig_fit) and sig_fit > 0)
        out = (float(mu_fit), float(sig_fit), float(A_fit), fit_ok)
        return (*out, dbg) if debug else out

    except Exception:
        occ_w = occ[(occ >= x_lo) & (occ <= x_hi)]
        if occ_w.size < 20:
            out = (np.nan, np.nan, np.nan, False)
            return (*out, dbg) if debug else out

        mu_fit = float(np.mean(occ_w))
        sig_fit = float(np.std(occ_w))
        A_fit = float(yw.max())
        fit_ok = bool(np.isfinite(mu_fit) and np.isfinite(sig_fit) and sig_fit > 0)
        out = (mu_fit, sig_fit, A_fit, fit_ok)
        return (*out, dbg) if debug else out


def plot_gaussian_fit_debug(mu_fit, sigma_fit, A_fit, fit_ok, dbg, window_sigmas=2.5):
    """Debug visualization: show occupancy-value histogram and highlight fit window and used points."""
    x_full = dbg["x_full"]
    y_full = dbg["y_full"]
    mode = dbg["mode"]
    sigma_guess = dbg["sigma_guess"]
    x_lo = dbg["x_lo"]
    x_hi = dbg["x_hi"]
    xw = dbg["xw"]
    yw = dbg["yw"]

    plt.figure()
    nz = y_full > 0
    plt.step(x_full[nz], y_full[nz], where="mid", linewidth=1.5, label="Occupancy-value histogram (nonzero)")

    plt.axvline(mode, linestyle="--", linewidth=1.2, label=f"Mode = {mode}")
    plt.axvline(x_lo, linestyle="--", linewidth=1.2, label=f"x_lo = {x_lo}")
    plt.axvline(x_hi, linestyle="--", linewidth=1.2, label=f"x_hi = {x_hi}")
    plt.axvspan(x_lo, x_hi, alpha=0.15,
                label=f"Fit window (±{window_sigmas}·σ_guess), σ_guess={sigma_guess:.2f}")

    plt.plot(xw, yw, marker="o", linestyle="none", markersize=4, label="Points used for fit")

    if fit_ok:
        xx = np.linspace(x_lo, x_hi, 400)
        yy = A_fit * np.exp(-0.5 * ((xx - mu_fit) / sigma_fit) ** 2)
        plt.plot(xx, yy, linewidth=1.5, label=f"Gaussian fit (μ={mu_fit:.2f}, σ={sigma_fit:.2f})")

    plt.xlabel("Occupancy (entries/bin)")
    plt.ylabel("Number of bins")
    plt.yscale("log")
    plt.title("Main-peak Gaussian fit debug")
    plt.legend()
    plt.tight_layout()
    plt.show()


# ----------------------------
# Threshold helpers
# ----------------------------
def compute_thresholds_from_occupancy(occ: np.ndarray, fit_window_sigmas: float):
    """
    Compute RMS threshold, Poisson threshold, and Gaussian-fit threshold (if fit_ok).
    Returns dict with values in occupancy units (entries/bin).
    """
    occ = np.asarray(occ, dtype=float)
    median = float(np.median(occ))
    rms = float(np.std(occ))
    thr_rms = median - 3.0 * rms

    sigma_poisson = float(np.sqrt(max(median, 1.0)))
    thr_pois = median - 3.0 * sigma_poisson

    mu_fit, sigma_fit, A_fit, fit_ok, dbg = fit_main_peak_gaussian(
        occ.astype(int), window_sigmas=fit_window_sigmas, debug=True
    )
    thr_fit = (mu_fit - 3.0 * sigma_fit) if fit_ok else np.nan

    return {
        "median": median,
        "rms": rms,
        "thr_rms": thr_rms,
        "thr_pois": thr_pois,
        "mu_fit": mu_fit,
        "sigma_fit": sigma_fit,
        "A_fit": A_fit,
        "fit_ok": fit_ok,
        "thr_fit": thr_fit,
        "dbg": dbg,
    }


def occ_to_khz(thr_occ: float, bin_width: float) -> float:
    """Convert a threshold in entries/bin to kHz."""
    return (thr_occ / bin_width) / 1e3


# ----------------------------
# Rate diagnostics
# ----------------------------
def rate_distribution_diagnostic(hit_times_ns: np.ndarray, cfg: Config):
    """
    Pre-spark diagnostic:
      - instantaneous hit rate per bin in kHz (distribution)
      - overlay thresholds in kHz: RMS, Poisson, Gaussian-fit
      - optional gaussian-fit debug plot

    Returns:
      live_time, n_hits, rate_per_bin_kHz
    """
    time_rel = prepare_time_axis(hit_times_ns, cfg.t_cut)
    if time_rel.size < 2:
        return 0.0, 0, np.array([])

    bins = np.arange(time_rel.min(), time_rel.max() + cfg.bin_width, cfg.bin_width)
    occ, _ = np.histogram(time_rel, bins=bins)

    live_time = float(time_rel.max() - time_rel.min())
    n_hits = int(time_rel.size)

    thr = compute_thresholds_from_occupancy(occ, cfg.fit_window_sigmas)

    # Optional debug fit window plot
    if cfg.mode == "debug" and cfg.debug_fit:
        plot_gaussian_fit_debug(
            thr["mu_fit"], thr["sigma_fit"], thr["A_fit"], thr["fit_ok"], thr["dbg"],
            window_sigmas=cfg.fit_window_sigmas
        )

    rate_per_bin_kHz = (occ / cfg.bin_width) / 1e3

    # Define bin width in kHz for the rate histogram
    rate_bin_width_kHz = 100e-6  # adjust to taste in seconds

    rate_edges = np.arange(
        rate_per_bin_kHz.min(),
        rate_per_bin_kHz.max() + rate_bin_width_kHz,
        rate_bin_width_kHz
    )

    if cfg.mode == "debug" and cfg.plot_per_file and rate_per_bin_kHz.size > 0:
        plt.figure()
        # plt.hist(rate_per_bin_kHz, bins=100, histtype="step", linewidth=1.6)
        plt.hist(rate_per_bin_kHz, bins=rate_edges, histtype="step", linewidth=1.6)
        # plt.axvline(occ_to_khz(thr["thr_rms"], cfg.bin_width), color="red", linestyle="--", linewidth=1.5,
        #             label=f"Thr = {occ_to_khz(thr['thr_rms'], cfg.bin_width):.2f} kHz (RMS)")
        plt.axvline(occ_to_khz(thr["thr_pois"], cfg.bin_width), color="blue", linestyle="--", linewidth=1.5,
                    label=f"Thr = {occ_to_khz(thr['thr_pois'], cfg.bin_width):.2f} kHz (Poisson)")
        if thr["fit_ok"]:
            plt.axvline(occ_to_khz(thr["thr_fit"], cfg.bin_width), color="green", linestyle="--", linewidth=1.5,
                        label=f"Thr = {occ_to_khz(thr['thr_fit'], cfg.bin_width):.2f} kHz (fit, σfit={thr['sigma_fit']:.2f})")

        plt.xlabel("Instantaneous hit rate per bin [kHz]")
        plt.ylabel("Number of bins")
        plt.yscale("log")
        plt.title("Distribution of instantaneous hit rate")
        plt.legend()
        plt.tight_layout()
        plt.show()

        global_hit_rate_Hz = (n_hits / live_time) if live_time > 0 else np.nan
        print("Rate diagnostics (after 20 s cut):")
        print(f"  N hits           = {n_hits}")
        print(f"  Live time        = {live_time:.3f} s")
        print(f"  Global hit rate  = {global_hit_rate_Hz:.3e} Hz")
        print(f"  Median(occ)      = {thr['median']:.2f}")
        print(f"  RMS(occ)         = {thr['rms']:.2f}")
        print(f"  Thr_RMS          = {thr['thr_rms']:.2f} entries/bin")
        print(f"  Thr_Poisson      = {thr['thr_pois']:.2f} entries/bin")
        if thr["fit_ok"]:
            print(f"  Fit: mu={thr['mu_fit']:.2f}, sigma={thr['sigma_fit']:.2f}")
            print(f"  Thr_Fit          = {thr['thr_fit']:.2f} entries/bin")
        else:
            print("  Fit: FAILED")

    return live_time, n_hits, rate_per_bin_kHz


# ----------------------------
# Spark finder
# ----------------------------
def find_sparks(hit_times_ns: np.ndarray, cfg: Config):
    """
    Find sparks by binning hit times into occupancies and identifying low-rate periods.
    Spark threshold uses cfg.spark_threshold_type.

    Returns:
      spark_durations [s], live_time [s], n_sparks
    """
    time_rel = prepare_time_axis(hit_times_ns, cfg.t_cut)
    if time_rel.size < 2:
        return np.array([]), 0.0, 0

    bins = np.arange(time_rel.min(), time_rel.max() + cfg.bin_width, cfg.bin_width)
    occ, bin_edges = np.histogram(time_rel, bins=bins)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    thr = compute_thresholds_from_occupancy(occ, cfg.fit_window_sigmas)

    # Choose spark threshold
    if cfg.spark_threshold_type.lower() == "poisson":
        spark_thr = thr["thr_pois"]
        spark_thr_label = f"Thr = {spark_thr:.1f} (median - 3·sqrt(median))"
    else:
        spark_thr = thr["thr_rms"]
        spark_thr_label = f"Thr = {spark_thr:.1f} (median - 3·RMS)"

    spark_mask = occ < spark_thr

    # Recovery condition
    if cfg.recovery_condition.lower() == "threshold":
        recovered = lambda v: v >= spark_thr
        recovery_label = "Recovery: above threshold"
    else:
        recovered = lambda v: v >= thr["median"]
        recovery_label = "Recovery: back to median"

    # Cluster consecutive below-threshold bins
    spark_starts, spark_ends = [], []
    in_spark = False
    for i in range(len(spark_mask)):
        if spark_mask[i] and not in_spark:
            spark_starts.append(bin_edges[i])
            in_spark = True
        elif recovered(occ[i]) and in_spark:
            spark_ends.append(bin_edges[i])
            in_spark = False

    # if a spark continues to the end, drop the incomplete one
    if in_spark and spark_starts:
        spark_starts.pop(-1)

    spark_durations = np.array([end - start for start, end in zip(spark_starts, spark_ends)], dtype=float)
    live_time = float(time_rel.max() - time_rel.min())

    # Per-file plot only in debug
    if cfg.mode == "debug" and cfg.plot_per_file:
        plt.figure()
        plt.step(bin_centers, occ, where="mid", color="black", linewidth=1.5, label="Hits per bin")

        # show all thresholds for reference
        # plt.axhline(thr["thr_rms"], color="red", linestyle="--", linewidth=1.5,
        #             label=f"Thr_RMS = {thr['thr_rms']:.1f}")
        plt.axhline(thr["thr_pois"], color="blue", linestyle="--", linewidth=1.5,
                    label=f"Thr_Poisson = {thr['thr_pois']:.1f}")
        if thr["fit_ok"]:
            plt.axhline(thr["thr_fit"], color="green", linestyle="--", linewidth=1.5,
                        label=f"Thr_Fit = {thr['thr_fit']:.1f} (σfit={thr['sigma_fit']:.2f})")

        # highlight the *active* spark threshold
        plt.axhline(spark_thr, color="black", linestyle=":", linewidth=1.2, label=f"USED: {spark_thr_label}")

        plt.axhline(thr["median"], color="salmon", linestyle="--", linewidth=1.5, label="Median")

        for start, end in zip(spark_starts, spark_ends):
            plt.axvspan(start, end, color="red", alpha=0.3)

        plt.xlabel("Time [s]")
        plt.ylabel("Entries / bin")
        plt.title(f"hits:time  ({recovery_label})")
        plt.legend()
        plt.tight_layout()
        plt.show()

        print(f"Detected sparks: {spark_durations.size}")
        if spark_durations.size:
            print("Spark durations (s):", spark_durations)

    return spark_durations, live_time, int(spark_durations.size)


# ----------------------------
# Per-file runner
# ----------------------------
def run_on_file(file_path: str, cfg: Config):
    df = load_root_file(file_path, branches=["time"])
    if df.empty:
        return {
            "n_hits": 0, "hit_live_time": 0.0, "rate_per_bin_kHz": np.array([]),
            "n_sparks": 0, "spark_live_time": 0.0, "spark_durations": np.array([])
        }

    hit_times = np.array(df["time"])

    # In analysis mode, we still compute diagnostics (needed for totals),
    # but plotting is suppressed.
    hit_live_time, n_hits, rate_per_bin_kHz = rate_distribution_diagnostic(hit_times, cfg)
    spark_durations, spark_live_time, n_sparks = find_sparks(hit_times, cfg)

    return {
        "n_hits": n_hits,
        "hit_live_time": hit_live_time,
        "rate_per_bin_kHz": rate_per_bin_kHz,
        "n_sparks": n_sparks,
        "spark_live_time": spark_live_time,
        "spark_durations": spark_durations
    }


# ----------------------------
# Main
# ----------------------------
def main():
    cfg = Config(
        # data_dir="/drf/projets/clas12/P2/spark_tests/sparks_test14", #AC-R
        # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_no_protection_jan26",
        # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_no_protection_dec25",
        # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_AC_dec25",
        # data_dir = "/drf/projets/clas12/P2/spark_tests/sparks_test15", #without protection
        # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_AC_jan26",
        # data_dir = "/drf/projets/clas12/P2/spark_tests/sparks_test16", #only AC
        # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_DC_dec25",
        # data_dir = "/drf/projets/clas12/P2/spark_tests/sparks_test11", #AC-DC
        # data_dir = "/mnt/data/P2_Basket_Analysis/spark_tests_data/sparks_AC_R_dec25",
        data_dir = "/drf/projets/clas12/cern_202511_p2_alinx_recovered/run_101",

        bin_width=0.01, #s
        # bin_width=0.01,
        # t_cut=20.0,
        t_cut=0.0,
        fit_window_sigmas=3,
        mode="debug",              # "debug" or "analysis"
        # mode="analysis",              # "debug" or "analysis"
        plot_per_file=True,        # per-file plots only if mode="debug"
        debug_fit=True,           # gaussian debug window plot only if mode="debug"
        # plot_global=True,          # global plots at the end
        spark_threshold_type="poisson",
        recovery_condition="median",
        max_files=1500,
        # file_stride=20,
    )

    set_root_style()

    # Accumulators
    total_hits = 0
    total_hit_live_time = 0.0
    all_rate_per_bin_kHz = []

    total_sparks = 0
    total_spark_live_time = 0.0
    all_spark_durations = []

    t_cursor_s = 0.0
    spark_rate_times_s = []
    spark_rate_kHz = []
    spark_rate_err_kHz = []

    files = sorted(
        f for f in os.listdir(cfg.data_dir)
        if f.startswith("enp") and f.endswith(".root")
    )

    # apply stride
    files = files[::cfg.file_stride]

    # apply max_files
    if cfg.max_files is not None:
        files = files[:cfg.max_files]

    for file_name in files:
        file_path = os.path.join(cfg.data_dir, file_name)
        print("Processing file:", file_path)

        out = run_on_file(file_path, cfg)
        # if out["n_sparks"] > 0:
        #     print(f"FOUND sparks in {file_name}: n_sparks={out['n_sparks']}, "
        #           f"live_time={out['spark_live_time']:.2f}s, "
        #           f"durations_sample={out['spark_durations'][:3]}")

        # --- spark-rate monitoring point per file ---
        if out["spark_live_time"] > 0:
            nsp = out["n_sparks"]
            T = out["spark_live_time"]

            rate_khz = (nsp / T) / 1e3
            err_khz = (np.sqrt(nsp) / T) / 1e3 if nsp > 0 else 0.0

            t_mid = t_cursor_s + 0.5 * T
            spark_rate_times_s.append(t_mid)
            spark_rate_kHz.append(rate_khz)
            spark_rate_err_kHz.append(err_khz)

            t_cursor_s += T

        total_hits += out["n_hits"]
        total_hit_live_time += out["hit_live_time"]
        if out["rate_per_bin_kHz"].size > 0:
            all_rate_per_bin_kHz.extend(out["rate_per_bin_kHz"])

        total_sparks += out["n_sparks"]
        total_spark_live_time += out["spark_live_time"]
        if out["spark_durations"].size > 0:
            all_spark_durations.extend(out["spark_durations"])

        if out["spark_live_time"] > 0:
            print(f"{file_name}: spark rate = {(out['n_sparks'] / out['spark_live_time'])/1e3:.3e} kHz")

    # Global hit-rate summary
    if total_hit_live_time > 0:
        global_hit_rate_Hz = total_hits / total_hit_live_time
        global_hit_rate_err_Hz = np.sqrt(total_hits) / total_hit_live_time
        print("\n=== Hit-rate summary (after 20 s cut) ===")
        print(f"Total hits: {total_hits}")
        print(f"Total hit live time: {total_hit_live_time:.1f} s")
        print(f"Global hit rate: ({global_hit_rate_Hz:.3e} ± {global_hit_rate_err_Hz:.1e}) Hz")
        print(f"Global hit rate: ({global_hit_rate_Hz/1e3:.3e} ± {global_hit_rate_err_Hz/1e3:.1e}) kHz")
    else:
        print("\nNo valid hit live time to compute global hit rate.")

    # Global instantaneous rate distribution
    if cfg.plot_global and len(all_rate_per_bin_kHz) > 0:
        plt.figure()
        plt.hist(all_rate_per_bin_kHz, bins=120, histtype="step", linewidth=1.6)
        plt.xlabel("Instantaneous hit rate per bin [kHz]")
        plt.ylabel("Number of bins")
        plt.yscale("log")
        plt.title("Instantaneous hit rate distribution (all files)")
        plt.tight_layout()
        plt.show()

    if cfg.plot_global and len(spark_rate_times_s) > 0:
        plt.figure()
        plt.errorbar(
            spark_rate_times_s,
            spark_rate_kHz,
            yerr=spark_rate_err_kHz,
            fmt="o",
            capsize=2,
            linewidth=1.2,
        )
        plt.xlabel("Run time (cumulative after t_cut) [s]")
        plt.ylabel("Spark rate per file [kHz]")
        plt.title("Spark-rate evolution vs run time")
        plt.tight_layout()
        plt.show()

    # Spark summary
    if total_spark_live_time > 0:
        sparking_rate = total_sparks / total_spark_live_time
        rate_err = np.sqrt(total_sparks) / total_spark_live_time
        print("\n=== Spark summary ===")
        print(f"Total sparks: {total_sparks}")
        print(f"Total live time: {total_spark_live_time:.1f} s")
        print(f"Sparking rate: ({sparking_rate:.3e} ± {rate_err:.1e}) Hz")
    else:
        print("\nNo valid live time to compute spark rate.")

    # Spark duration distribution
    if cfg.plot_global and len(all_spark_durations) > 0:
        all_spark_durations = np.asarray(all_spark_durations, dtype=float)
        fig, ax = plt.subplots()
        ax.hist(all_spark_durations, bins=20, color="blue", alpha=0.7)
        mean = float(np.mean(all_spark_durations))
        rms = float(np.std(all_spark_durations))
        # ax.axvline(mean, color="red", linestyle="--")
        ax.annotate(
            f"Sparks: {len(all_spark_durations)}\n"
            f"Mean: {mean:.2f} s\n"
            f"RMS: {rms:.2f} s\n",
            xy=(0.75, 0.95),
            xycoords="axes fraction",
            fontsize=14,
            va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.5),
        )
        ax.set_xlabel("Spark Duration [s]")
        ax.set_ylabel("Counts")
        ax.set_title("Distribution of Spark Durations")
        plt.tight_layout()
        plt.show()
    elif cfg.plot_global:
        print("\nNo sparks found → spark duration histogram skipped.")


if __name__ == "__main__":
    main()
    print("bonzo")
