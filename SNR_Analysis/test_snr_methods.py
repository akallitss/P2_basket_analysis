#!/usr/bin/env python3
"""
test_snr_methods.py

Smoke test for the three new SNR metrics and their plots.
No data files needed — uses synthetic histograms.

Run:
    python test_snr_methods.py
    python test_snr_methods.py --save /tmp/snr_test_plots/
"""
import argparse
import numpy as np
import pandas as pd

from vmm_snr import _compute_area_ratio, _compute_eer, _compute_p_saturation
from vmm_plots import (plot_all_methods_heatmap,
                       plot_tail_distributions,
                       plot_saturation_curves)

RNG = np.random.default_rng(42)


# ─────────────────────────────────────────────────────────────────
# Synthetic histogram builders
# ─────────────────────────────────────────────────────────────────

def _make_noise_hist(mean=80, sigma=5, n=50_000):
    """Gaussian noise centred near ADC=80."""
    vals = RNG.normal(mean, sigma, n).clip(21, 1023).astype(int)
    return np.bincount(vals, minlength=1024)[:1024].astype(np.int64)


def _make_signal_hist(mpv=180, width=30, tail=40, n=20_000, p_sat=0.0):
    """Landau-like signal: normal core + exponential tail."""
    core = RNG.normal(mpv, width, n) + RNG.exponential(tail, n)
    core = core.clip(0, 1021)
    hist = np.bincount(core.astype(int), minlength=1024)[:1024].astype(np.int64)
    # inject a controllable saturation fraction
    n_sat = int(n * p_sat)
    hist[1023] += n_sat
    return hist


def _make_tail_entry(sg, snt, vmm_id, noise_hist, signal_hist, noise_cut):
    """Build one tail_data dict from pre-made histograms."""
    adc      = np.arange(1024, dtype=np.float64)
    n_noise  = float(noise_hist.sum())
    n_signal = float(signal_hist.sum())

    noise_sf   = 1.0 - np.cumsum(noise_hist)  / n_noise
    signal_cdf = np.cumsum(signal_hist) / n_signal
    signal_sf  = 1.0 - signal_cdf

    from vmm_snr import _compute_eer, _compute_p_saturation
    eer_thr, eer_val = _compute_eer(signal_hist, noise_hist)
    p_sat            = _compute_p_saturation(signal_hist)

    return {
        "sg"           : sg,
        "snt"          : snt,
        "vmm_id"       : vmm_id,
        "adc"          : adc,
        "noise_sf"     : noise_sf,
        "signal_cdf"   : signal_cdf,
        "signal_sf"    : signal_sf,
        "noise_cut"    : noise_cut,
        "eer_threshold": eer_thr,
        "eer_value"    : eer_val,
        "p_saturation" : p_sat,
    }


# ─────────────────────────────────────────────────────────────────
# Helper metric tests
# ─────────────────────────────────────────────────────────────────

def test_helpers():
    print("─" * 55)
    print("Testing helper functions with synthetic histograms")
    print("─" * 55)

    noise_hist  = _make_noise_hist(mean=80, sigma=5)
    signal_hist = _make_signal_hist(mpv=180, width=30, tail=40, p_sat=0.02)
    noise_cut   = 105.0   # mean + 5*sigma

    # --- area_ratio ---
    ar = _compute_area_ratio(signal_hist, noise_hist, noise_cut)
    print(f"  area_ratio    = {ar:.1f}   (expect >> 1 for well-separated signal)")
    assert ar > 5, f"area_ratio suspiciously low: {ar:.2f}"

    # --- EER ---
    eer_thr, eer_val = _compute_eer(signal_hist, noise_hist)
    print(f"  eer_threshold = {eer_thr:.0f}    (expect between ~100 and ~180)")
    print(f"  eer_value     = {eer_val:.4f}  (expect << 0.5 for good separation)")
    assert 90 < eer_thr < 200,  f"EER threshold out of expected range: {eer_thr}"
    assert eer_val < 0.3,        f"EER value too high (poor separation?): {eer_val:.4f}"

    # --- p_saturation ---
    p_sat = _compute_p_saturation(signal_hist)
    print(f"  p_saturation  = {p_sat*100:.2f}%  (injected 2% sat — expect ~2%)")
    assert 0.01 < p_sat < 0.04, f"p_saturation outside expected range: {p_sat*100:.2f}%"

    # Edge case: empty histograms
    empty = np.zeros(1024, np.int64)
    # empty signal + real noise → ratio = 0, not NaN
    assert _compute_area_ratio(empty, noise_hist, noise_cut) == 0.0
    # both empty noise → ratio is NaN (division by zero)
    assert np.isnan(_compute_area_ratio(empty, empty, noise_cut))
    assert np.isnan(_compute_eer(empty, noise_hist)[0])
    assert np.isnan(_compute_p_saturation(empty))
    print("  Edge cases (empty histograms): OK")

    # Sanity check: heavily overlapping distributions → high EER value
    overlapping_noise   = _make_noise_hist(mean=150, sigma=20, n=50_000)
    overlapping_signal  = _make_signal_hist(mpv=160, width=20, tail=10, n=20_000)
    _, eer_val_bad = _compute_eer(overlapping_signal, overlapping_noise)
    _, eer_val_good = _compute_eer(signal_hist, noise_hist)
    print(f"  EER good separation = {eer_val_good:.4f}  (separated distributions)")
    print(f"  EER poor separation = {eer_val_bad:.4f}   (overlapping distributions)")
    assert eer_val_good < eer_val_bad, \
        "Well-separated distributions should have lower EER than overlapping ones"

    print("  All helper tests PASSED\n")


# ─────────────────────────────────────────────────────────────────
# Synthetic dataset: 2 VMMs × 4 configs (2 gains × 2 peaking times)
# ─────────────────────────────────────────────────────────────────

def build_synthetic_dataset():
    """Return (df_snr, tail_data, detector_vmms) with 4 configs × 2 VMMs."""
    configs = [
        {"sg": 1.0, "snt": 100, "noise_mean": 80,  "noise_sigma": 5,
         "mpv": 160, "width": 25, "p_sat": 0.005},
        {"sg": 1.0, "snt": 200, "noise_mean": 80,  "noise_sigma": 7,
         "mpv": 175, "width": 30, "p_sat": 0.008},
        {"sg": 4.5, "snt": 100, "noise_mean": 85,  "noise_sigma": 8,
         "mpv": 210, "width": 35, "p_sat": 0.03},
        {"sg": 4.5, "snt": 200, "noise_mean": 85,  "noise_sigma": 12,
         "mpv": 230, "width": 40, "p_sat": 0.05},
    ]
    vmm_ids      = [8, 12]
    detector_vmms = vmm_ids
    rows         = []
    tail_data    = []

    for vmm_id in vmm_ids:
        for cfg in configs:
            noise_hist  = _make_noise_hist(cfg["noise_mean"], cfg["noise_sigma"])
            signal_hist = _make_signal_hist(cfg["mpv"], cfg["width"],
                                            p_sat=cfg["p_sat"])
            noise_cut   = cfg["noise_mean"] + 5 * cfg["noise_sigma"]

            area_ratio          = _compute_area_ratio(signal_hist, noise_hist,
                                                       noise_cut)
            eer_thr, eer_val    = _compute_eer(signal_hist, noise_hist)
            p_sat               = _compute_p_saturation(signal_hist)

            # Rough MPV and noise_sigma for SNR columns
            adc       = np.arange(1024, dtype=float)
            mpv_idx   = int(np.argmax(signal_hist[int(noise_cut):]) + noise_cut)
            mpv       = float(mpv_idx)
            n_sig_sigma = cfg["noise_sigma"]
            snr       = mpv / n_sig_sigma if n_sig_sigma > 0 else np.nan
            mean_sig  = float((adc * signal_hist).sum()) / signal_hist.sum()
            snr_mean  = mean_sig / n_sig_sigma

            rows.append({
                "sg"           : cfg["sg"],
                "snt"          : cfg["snt"],
                "run_sng0"     : 0,
                "run_sng1"     : 1,
                "vmm_id"       : vmm_id,
                "noise_sigma"  : n_sig_sigma,
                "noise_cut"    : noise_cut,
                "noise_quality": "ok",
                "mpv"          : mpv,
                "snr"          : snr,
                "mean_signal"  : mean_sig,
                "mean_noise"   : cfg["noise_mean"],
                "snr_mean"     : snr_mean,
                "area_ratio"   : area_ratio,
                "eer_threshold": eer_thr,
                "eer_value"    : eer_val,
                "p_saturation" : p_sat,
            })

            tail_data.append(
                _make_tail_entry(cfg["sg"], cfg["snt"], vmm_id,
                                 noise_hist, signal_hist, noise_cut)
            )

    return pd.DataFrame(rows), tail_data, detector_vmms


# ─────────────────────────────────────────────────────────────────
# Plot tests
# ─────────────────────────────────────────────────────────────────

def test_plots(save_dir):
    print("─" * 55)
    print(f"Testing plot functions  (save_dir={save_dir or 'None (display only)'})")
    print("─" * 55)

    df_snr, tail_data, detector_vmms = build_synthetic_dataset()
    print(f"  Synthetic dataset: {len(df_snr)} rows, "
          f"{len(tail_data)} tail-curve entries")

    print("  plot_all_methods_heatmap ...", end=" ", flush=True)
    plot_all_methods_heatmap(df_snr, out_dir=save_dir, show=False)
    print("OK")

    print("  plot_tail_distributions  ...", end=" ", flush=True)
    plot_tail_distributions(tail_data, detector_vmms,
                            out_dir=save_dir, show=False)
    print("OK")

    print("  plot_saturation_curves   ...", end=" ", flush=True)
    plot_saturation_curves(tail_data, detector_vmms,
                           out_dir=save_dir, show=False)
    print("OK")

    if save_dir:
        import os
        saved = sorted(os.listdir(save_dir))
        print(f"\n  Saved {len(saved)} file(s) to {save_dir}:")
        for f in saved:
            print(f"    {f}")
    print()


# ─────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--save", metavar="DIR", default=None,
                        help="directory to save test plots (PDF + PNG)")
    args = parser.parse_args()

    test_helpers()
    test_plots(args.save)
    print("All tests passed.")
