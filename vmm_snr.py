#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vmm_snr.py

SNR computation for VMM config scan analysis.

Two levels of granularity:
- VMM level    : compute_snr()             one SNR per VMM per config
- Channel level: compute_snr_per_channel() one SNR per channel per config

Both use identical methodology:
- Noise sigma : robust MAD from sng=1 run (over_threshold=0, adc>20)
- Signal MPV  : smoothed histogram peak from sng=0 run
                (over_threshold=1, adc<1023)
- SNR         : MPV / noise_sigma

The noise baseline from the sng=1 run anchors both levels —
the VMM-level noise_cut is used as adc_min for the MPV search
at both VMM and channel level, ensuring consistent signal/noise
boundary definition across all granularities.
"""
import os
import numpy as np
import pandas as pd

from vmm_io     import get_run_dir, load_hits_run, iter_hits_files
from vmm_noise  import compute_noise_baseline, compute_noise_baseline_from_hists
from vmm_signal import get_clean_signal, estimate_mpv, estimate_mpv_from_hist
from vmm_mapping import vmm_mapping


def compute_adc_stats(df_run_scan, vmm_ids, run_list, data_dir,
                       adc_cut=100, n_files=1,
                       use_over_threshold=True):
    """
    Compute ADC statistics per VMM per run.

    Parameters
    ----------
    df_run_scan : pd.DataFrame
        Run metadata (must include 'run_no' and 'snt').
    vmm_ids : list[int]
        List of VMM IDs to process.
    run_list : list[int]
        List of run numbers.
    data_dir : str
        Base directory with run data.
    adc_cut : int, optional
        Maximum ADC value to consider (if use_over_threshold=False).
    n_files : int
        Number of ROOT files to load and merge per run directory.
    use_over_threshold : bool, optional
        If True, select hits with over_threshold==0 instead of ADC cut.

    Returns
    -------
    pd.DataFrame
        Statistics per VMM per run.
    """
    results = []

    # Outer loop over runs — load each run once, not once per VMM.
    for run_no in run_list:
        run_dir      = get_run_dir(data_dir, run_no)
        peaking_time = df_run_scan.loc[
            df_run_scan["run_no"] == run_no, "snt"
        ].iloc[0]

        # Accumulate per-VMM arrays via streaming (one file at a time).
        vmm_full_adc = {v: [] for v in vmm_ids}
        vmm_ot       = {v: [] for v in vmm_ids}

        for df in iter_hits_files(run_dir, n_files,
                                  branches=["adc", "vmm",
                                            "ch", "time",
                                            "over_threshold"]):
            for vmm_id in vmm_ids:
                mask = df["vmm"] == vmm_id
                vmm_full_adc[vmm_id].append(df[mask]["adc"].values)
                vmm_ot[vmm_id].append(
                    df[mask]["over_threshold"].values
                )
            del df

        for vmm_id in vmm_ids:
            if not vmm_full_adc[vmm_id]:
                continue
            full_adc = np.concatenate(vmm_full_adc[vmm_id])
            all_ot   = np.concatenate(vmm_ot[vmm_id])

            if len(full_adc) == 0:
                continue

            if use_over_threshold:
                adc_values = full_adc[all_ot == 0]
            else:
                adc_values = full_adc[full_adc < adc_cut]

            if len(adc_values) == 0:
                continue

            median_adc   = np.median(adc_values)
            mad_adc      = np.median(np.abs(adc_values - median_adc))
            robust_sigma = 1.4826 * mad_adc
            frac_removed = (len(full_adc) - len(adc_values)) / len(full_adc)
            mean_adc     = adc_values.mean()
            rms_adc      = adc_values.std()
            num_hits     = len(adc_values)
            rms_error    = (rms_adc / np.sqrt(2 * (num_hits - 1))
                            if num_hits > 1 else np.nan)

            results.append({
                "run_no"      : run_no,
                "vmm_id"      : vmm_id,
                "peaking_time": peaking_time,
                "num_hits"    : num_hits,
                "mean_adc"    : mean_adc,
                "rms_adc"     : rms_adc,
                "median_adc"  : median_adc,
                "robust_sigma": robust_sigma,
                "frac_removed": frac_removed,
                "rms_error_adc": rms_error
            })

        del vmm_full_adc, vmm_ot

    return pd.DataFrame(results)


def compute_snr(data_dir, pairs, detector_vmms,
                n_files=1,
                exclude_trigger_vmms=(0, 1)):
    """
    Compute VMM-level SNR for each configuration pair.

    Parameters
    ----------
    data_dir : str
        Base directory containing run_* folders.
    pairs : list of dict
        Output of get_run_groups()['pairs'].
        Each dict has keys: sg, snt, sng0, sng1.
    detector_vmms : list of int
        VMM IDs to process (excludes trigger VMMs).
    n_files : int
        Number of ROOT files to load and merge per run directory.
    exclude_trigger_vmms : tuple
        VMM IDs to skip entirely.

    Returns
    -------
    pd.DataFrame
        One row per (sg, snt, vmm_id) with columns:
        sg, snt, run_sng0, run_sng1, vmm_id,
        noise_sigma, noise_cut, noise_quality, mpv, snr
    """
    results = []

    for pair in pairs:
        sg       = pair["sg"]
        snt      = pair["snt"]
        run_sng0 = pair["sng0"]
        run_sng1 = pair["sng1"]

        print(f"\nProcessing sg={sg} snt={snt} "
              f"| noise=run {run_sng1} | signal=run {run_sng0}")

        # --- sng=1 run → noise histograms (one file at a time) ---
        # Accumulate ADC counts per VMM into 1024-bin integer arrays.
        # Peak memory = n_vmms * 1024 * 8 bytes, independent of file count.
        run_dir     = get_run_dir(data_dir, run_sng1)
        noise_hists = {v: np.zeros(1024, np.int64) for v in detector_vmms}
        n_noise_files = 0
        for df in iter_hits_files(run_dir, n_files,
                                   branches=["adc", "vmm",
                                             "over_threshold"]):
            df_det = df[df["vmm"].isin(detector_vmms)]
            for vmm_id in detector_vmms:
                adc = df_det.loc[
                    (df_det["vmm"]            == vmm_id) &
                    (df_det["over_threshold"] == 0) &
                    (df_det["adc"]            >  20), "adc"
                ].values.clip(0, 1023).astype(np.int32)
                if len(adc):
                    noise_hists[vmm_id] += np.bincount(
                        adc, minlength=1024
                    )[:1024]
            n_noise_files += 1
            del df, df_det

        if n_noise_files == 0:
            print(f"  WARNING: no files for run {run_sng1}, skipping")
            continue

        df_noise_baseline = compute_noise_baseline_from_hists(noise_hists)
        del noise_hists

        if "vmm_id" not in df_noise_baseline.columns:
            print(f"  WARNING: no noise baseline for run {run_sng1} "
                  f"(too few over_threshold=0 hits), skipping pair")
            continue

        # --- sng=0 run → signal histograms (one file at a time) ---
        run_dir      = get_run_dir(data_dir, run_sng0)
        signal_hists = {v: np.zeros(1024, np.int64) for v in detector_vmms}
        n_signal_files = 0
        for df in iter_hits_files(run_dir, n_files,
                                   branches=["adc", "vmm",
                                             "over_threshold"]):
            df_det = df[df["vmm"].isin(detector_vmms)]
            for vmm_id in detector_vmms:
                adc = df_det.loc[
                    (df_det["vmm"]            == vmm_id) &
                    (df_det["over_threshold"] == 1) &
                    (df_det["adc"]            <  1023), "adc"
                ].values.clip(0, 1023).astype(np.int32)
                if len(adc):
                    signal_hists[vmm_id] += np.bincount(
                        adc, minlength=1024
                    )[:1024]
            n_signal_files += 1
            del df, df_det

        if n_signal_files == 0:
            print(f"  WARNING: no files for run {run_sng0}, skipping")
            continue

        # --- Per VMM ---
        for vmm_id in detector_vmms:
            if vmm_id in exclude_trigger_vmms:
                continue

            noise_row = df_noise_baseline[
                df_noise_baseline["vmm_id"] == vmm_id
            ]
            if noise_row.empty:
                continue

            noise_sigma   = noise_row["robust_sigma"].iloc[0]
            noise_cut     = noise_row["noise_cut"].iloc[0]
            noise_quality = noise_row["quality"].iloc[0]

            if signal_hists[vmm_id].sum() < 100:
                continue

            mpv, _, _, _ = estimate_mpv_from_hist(
                signal_hists[vmm_id], adc_min=noise_cut
            )
            if np.isnan(mpv):
                continue
            snr = mpv / noise_sigma if noise_sigma > 0 else np.nan

            print(f"  VMM {vmm_id}: noise_sigma={noise_sigma:.1f}  "
                  f"MPV={mpv:.1f}  SNR={snr:.1f}  [{noise_quality}]")

            results.append({
                "sg"           : sg,
                "snt"          : snt,
                "run_sng0"     : run_sng0,
                "run_sng1"     : run_sng1,
                "vmm_id"       : vmm_id,
                "noise_sigma"  : noise_sigma,
                "noise_cut"    : noise_cut,
                "noise_quality": noise_quality,
                "mpv"          : mpv,
                "snr"          : snr,
            })

        del signal_hists, df_noise_baseline

    return pd.DataFrame(results)


def compute_snr_per_channel(data_dir, pairs, detector_vmms,
                             n_files=1,
                             exclude_trigger_vmms=(0, 1),
                             min_signal_hits=50,
                             min_noise_hits=50,
                             min_sigma=2.0,
                             max_sigma=20.0,
                             mpv_min=100,
                             mpv_max=300,
                             bins=80,
                             smooth_sigma=2):
    """
    Compute channel-level SNR for each configuration pair.

    Uses identical cuts and methodology as compute_snr() but
    applied per channel. The VMM-level noise_cut is used as
    adc_min for channel-level MPV search — ensuring consistency
    between VMM and channel level analyses.

    Quality cuts per channel:
    - min_noise_hits : minimum noise hits for stable MAD
    - min_sigma      : removes stuck channels (MAD floor artifact)
    - max_sigma      : removes floating/noisy channels
    - mpv_min/max    : removes MPV finder artifacts in sparse
                       histograms (too few hits → peak in tail)
    - min_signal_hits: minimum signal hits for reliable MPV

    Parameters
    ----------
    n_files : int
        Number of ROOT files to load and merge per run directory.
    bins : int
        Histogram bins for MPV — lower than VMM level to match
        reduced per-channel statistics.
    smooth_sigma : float
        Gaussian smoothing for MPV peak finder.

    Returns
    -------
    pd.DataFrame
        One row per (sg, snt, vmm_id, ch) with columns:
        sg, snt, run_sng0, run_sng1, vmm_id, ch,
        noise_sigma_ch, vmm_noise_cut, vmm_noise_quality,
        mpv_ch, snr_ch, n_signal, n_noise
    """
    results = []

    for pair in pairs:
        sg       = pair["sg"]
        snt      = pair["snt"]
        run_sng0 = pair["sng0"]
        run_sng1 = pair["sng1"]

        print(f"\nProcessing sg={sg} snt={snt} "
              f"| noise=run {run_sng1} | signal=run {run_sng0}")

        # --- sng=1 run → per-(vmm,ch) noise histograms ---
        run_dir        = get_run_dir(data_dir, run_sng1)
        noise_hists_ch = {}   # (vmm_id, ch_id) -> np.ndarray(1024)
        noise_hists_vmm = {v: np.zeros(1024, np.int64)
                           for v in detector_vmms}
        n_noise_files  = 0
        for df in iter_hits_files(run_dir, n_files,
                                   branches=["adc", "vmm", "ch",
                                             "over_threshold"]):
            df_det = df[
                df["vmm"].isin(detector_vmms) &
                (df["over_threshold"] == 0) &
                (df["adc"] > 20)
            ]
            if not df_det.empty:
                vmms = df_det["vmm"].values.astype(np.int32)
                chs  = df_det["ch"].values.astype(np.int32)
                adcs = df_det["adc"].values.clip(0, 1023).astype(np.int32)
                for vmm_id in np.unique(vmms):
                    if vmm_id not in detector_vmms:
                        continue
                    vm = vmms == vmm_id
                    noise_hists_vmm[vmm_id] += np.bincount(
                        adcs[vm], minlength=1024
                    )[:1024]
                    for ch_id in np.unique(chs[vm]):
                        key = (int(vmm_id), int(ch_id))
                        counts = np.bincount(
                            adcs[vm & (chs == ch_id)], minlength=1024
                        )[:1024]
                        if key not in noise_hists_ch:
                            noise_hists_ch[key] = np.zeros(1024, np.int64)
                        noise_hists_ch[key] += counts
            n_noise_files += 1
            del df, df_det

        if n_noise_files == 0:
            print(f"  WARNING: no files for run {run_sng1}, skipping")
            continue

        # VMM-level noise baseline (for noise_cut used in MPV search)
        df_noise_baseline = compute_noise_baseline_from_hists(
            noise_hists_vmm
        )
        del noise_hists_vmm

        if "vmm_id" not in df_noise_baseline.columns:
            print(f"  WARNING: no noise baseline for run {run_sng1} "
                  f"(too few over_threshold=0 hits), skipping pair")
            continue

        # --- sng=0 run → per-(vmm,ch) signal histograms ---
        run_dir         = get_run_dir(data_dir, run_sng0)
        signal_hists_ch = {}   # (vmm_id, ch_id) -> np.ndarray(1024)
        n_signal_files  = 0
        for df in iter_hits_files(run_dir, n_files,
                                   branches=["adc", "vmm", "ch",
                                             "over_threshold"]):
            df_det = df[
                df["vmm"].isin(detector_vmms) &
                (df["over_threshold"] == 1) &
                (df["adc"] < 1023)
            ]
            if not df_det.empty:
                vmms = df_det["vmm"].values.astype(np.int32)
                chs  = df_det["ch"].values.astype(np.int32)
                adcs = df_det["adc"].values.clip(0, 1023).astype(np.int32)
                for vmm_id in np.unique(vmms):
                    if vmm_id not in detector_vmms:
                        continue
                    vm = vmms == vmm_id
                    for ch_id in np.unique(chs[vm]):
                        key = (int(vmm_id), int(ch_id))
                        counts = np.bincount(
                            adcs[vm & (chs == ch_id)], minlength=1024
                        )[:1024]
                        if key not in signal_hists_ch:
                            signal_hists_ch[key] = np.zeros(1024, np.int64)
                        signal_hists_ch[key] += counts
            n_signal_files += 1
            del df, df_det

        if n_signal_files == 0:
            print(f"  WARNING: no files for run {run_sng0}, skipping")
            continue

        n_ok = n_low_n = n_stuck = n_noisy = n_mpv = n_sig = 0

        for vmm_id in detector_vmms:
            if vmm_id in exclude_trigger_vmms:
                continue

            noise_row = df_noise_baseline[
                df_noise_baseline["vmm_id"] == vmm_id
            ]
            if noise_row.empty:
                continue

            vmm_noise_cut     = noise_row["noise_cut"].iloc[0]
            vmm_noise_quality = noise_row["quality"].iloc[0]

            sig_keys = [k for k in signal_hists_ch
                        if k[0] == vmm_id]

            for key in sorted(sig_keys):
                ch_id     = key[1]
                sig_hist  = signal_hists_ch[key]
                noi_hist  = noise_hists_ch.get(key,
                                np.zeros(1024, np.int64))

                n_noise_ch = int(noi_hist.sum())
                if n_noise_ch < min_noise_hits:
                    n_low_n += 1
                    continue

                # Channel-level noise sigma from histogram
                adc_vals   = np.arange(1024, dtype=np.float64)
                cumsum_n   = np.cumsum(noi_hist)
                med_idx    = min(int(np.searchsorted(
                    cumsum_n, n_noise_ch / 2)), 1023)
                median_n   = adc_vals[med_idx]
                devs       = np.abs(adc_vals - median_n)
                order      = np.argsort(devs, kind="stable")
                cum_dev    = np.cumsum(noi_hist[order])
                mad_idx    = min(int(np.searchsorted(
                    cum_dev, n_noise_ch / 2)), 1023)
                mad_n      = float(devs[order[mad_idx]])
                noise_sigma = 1.4826 * mad_n

                if noise_sigma <= min_sigma:
                    n_stuck += 1
                    continue
                if noise_sigma >= max_sigma:
                    n_noisy += 1
                    continue

                n_signal_ch = int(sig_hist.sum())
                if n_signal_ch < min_signal_hits:
                    n_sig += 1
                    continue

                mpv, _, _, _ = estimate_mpv_from_hist(
                    sig_hist,
                    adc_min=vmm_noise_cut,
                    smooth_sigma=smooth_sigma
                )
                if np.isnan(mpv) or not (mpv_min <= mpv <= mpv_max):
                    n_mpv += 1
                    continue

                snr  = mpv / noise_sigma
                n_ok += 1

                results.append({
                    "sg"               : sg,
                    "snt"              : snt,
                    "run_sng0"         : run_sng0,
                    "run_sng1"         : run_sng1,
                    "vmm_id"           : vmm_id,
                    "ch"               : ch_id,
                    "noise_sigma_ch"   : noise_sigma,
                    "vmm_noise_cut"    : vmm_noise_cut,
                    "vmm_noise_quality": vmm_noise_quality,
                    "mpv_ch"           : mpv,
                    "snr_ch"           : snr,
                    "n_signal"         : n_signal_ch,
                    "n_noise"          : n_noise_ch,
                })

        print(f"  → {n_ok} channels ok | "
              f"low_n={n_low_n} stuck={n_stuck} "
              f"noisy={n_noisy} bad_mpv={n_mpv} "
              f"low_sig={n_sig}")

        del noise_hists_ch, signal_hists_ch, df_noise_baseline

    return pd.DataFrame(results)


def summarise_best_config(df_snr, df_snr_ch, detector_vmms):
    """
    Print and return a summary table of the best configuration
    per VMM at both VMM and channel level.

    Parameters
    ----------
    df_snr : pd.DataFrame
        VMM-level SNR results from compute_snr().
    df_snr_ch : pd.DataFrame
        Channel-level SNR results from compute_snr_per_channel().
    detector_vmms : list of int
        VMM IDs to include in summary.

    Returns
    -------
    pd.DataFrame
        One row per VMM with best config at both levels
        and an agreement flag.
    """
    vmm_to_detector = {
        vid: cfg.get("name", key)
        for key, cfg in vmm_mapping.items()
        for vid in cfg["vmm_ids"]
    }

    records = []

    print("\n" + "="*75)
    print("BEST CONFIGURATION SUMMARY PER VMM")
    print("="*75)
    print(f"{'VMM':>5} {'Detector':>22} "
          f"{'VMM_best':>14} {'VMM_SNR':>8} "
          f"{'CH_best':>14} {'CH_med':>7} {'Agreement':>10}")
    print("-"*75)

    for vmm_id in sorted(detector_vmms):

        df_v = df_snr[
            (df_snr["vmm_id"]        == vmm_id) &
            (df_snr["noise_quality"] == "ok")
        ]
        if df_v.empty:
            continue
        best_v  = df_v.loc[df_v["snr"].idxmax()]
        vmm_cfg = f"sg={best_v['sg']:.1f}/snt={best_v['snt']:.0f}"
        vmm_snr = best_v["snr"]

        df_c = df_snr_ch[df_snr_ch["vmm_id"] == vmm_id]
        if df_c.empty:
            continue
        ch_med     = df_c.groupby(["sg", "snt"])["snr_ch"].median()
        best_c_idx = ch_med.idxmax()
        best_c_snr = ch_med.max()
        ch_cfg     = (f"sg={best_c_idx[0]:.1f}"
                      f"/snt={best_c_idx[1]:.0f}")

        agreement = "✓" if vmm_cfg == ch_cfg else "✗"
        detector  = vmm_to_detector.get(vmm_id, "")

        print(f"{vmm_id:>5} {detector:>22} "
              f"{vmm_cfg:>14} {vmm_snr:>8.1f} "
              f"{ch_cfg:>14} {best_c_snr:>7.1f} "
              f"{agreement:>10}")

        records.append({
            "vmm_id"         : vmm_id,
            "detector"       : detector,
            "vmm_best_config": vmm_cfg,
            "vmm_snr"        : vmm_snr,
            "ch_best_config" : ch_cfg,
            "ch_median_snr"  : best_c_snr,
            "agreement"      : agreement == "✓"
        })

    df_summary = pd.DataFrame(records)

    print("="*75)
    print("\nOVERALL RECOMMENDATION:")

    vmm_votes = df_summary["vmm_best_config"].value_counts()
    ch_votes  = df_summary["ch_best_config"].value_counts()

    print(f"\n  VMM-level votes:")
    for cfg, n in vmm_votes.items():
        print(f"    {cfg} → {n} VMM(s)")

    print(f"\n  Channel-level votes:")
    for cfg, n in ch_votes.items():
        print(f"    {cfg} → {n} VMM(s)")

    overall_best = vmm_votes.idxmax()
    print(f"\n  → Recommended configuration: {overall_best}")
    print(f"    (wins {vmm_votes.max()} of {len(detector_vmms)} "
          f"VMMs at VMM level)")

    n_agree = df_summary["agreement"].sum()
    print(f"\n  VMM and channel level agree for "
          f"{n_agree}/{len(df_summary)} VMMs")
    print("="*75)

    return df_summary

def main():
    print('bonzo')


if __name__ == '__main__':
    main()
