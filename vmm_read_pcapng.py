#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vmm_read_pcapng.py

Read VMM3a SRS hit data from a raw pcapng capture file and identify
which VMM configuration (gain=sg, peaking time=snt) was active during
the run by matching the measured ADC MPV against the config scan
lookup table.

Binary decoding chain
---------------------
  pcapng Enhanced Packet Block
    → Ethernet II frame (14 bytes header)
    → IPv4 packet (20 bytes header)
    → UDP datagram (8 bytes header)
    → SRS payload (16-byte SRS header + N×6-byte hit/marker records)
    → VMM3a hit: vmm, channel, adc, over_threshold, bcid, tdc

SRS data format reference: ParserSRS.cpp / VMM3Parser.h from vmm-sdat
"""

import struct
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.ndimage import gaussian_filter1d


# ── Constants ────────────────────────────────────────────────────────────────

SRS_MAGIC     = 0x564D3300          # "VMM3" — top 24 bits of SRS dataId field
ETH_HDR_SIZE  = 14
IP_HDR_SIZE   = 20
UDP_HDR_SIZE  = 8
PKT_DATA_OFF  = ETH_HDR_SIZE + IP_HDR_SIZE + UDP_HDR_SIZE   # 42 bytes
SRS_HDR_SIZE  = 16                   # frameCounter + dataId + udpTS + offsetOverflow
HIT_SIZE      = 6                    # 4-byte data1 + 2-byte data2

TRIGGER_VMMS  = {0, 1}

VMM_TO_DETECTOR = {
    0:  "Trigger (fixed)",
    1:  "Trigger (fixed)",
    8:  "P2 Small Det 3 — X",
    9:  "P2 Small Det 3 — Y",
    10: "P2 Small Det 1 — X",
    11: "P2 Small Det 1 — Y",
    12: "P2 Large Det — A",
    13: "P2 Large Det — B",
    14: "P2 Large Det — C",
    15: "P2 Large Det — D",
}


# ── pcapng / SRS parsing ─────────────────────────────────────────────────────

def _gray2bin(n: int) -> int:
    mask = n >> 1
    while mask:
        n ^= mask
        mask >>= 1
    return n


def iter_srs_payloads(pcapng_path: Path):
    """
    Generator: stream UDP payloads from a pcapng file.

    Skips non-IPv4, non-UDP, and non-SRS packets.
    Yields raw bytes of each SRS UDP payload.
    """
    EPB_TYPE = 0x00000006

    with open(pcapng_path, "rb") as f:
        while True:
            hdr = f.read(8)
            if len(hdr) < 8:
                break
            btype, blen = struct.unpack("<II", hdr)
            body = f.read(blen - 12)
            f.read(4)                             # trailing block length

            if btype != EPB_TYPE:
                continue

            cap_len = struct.unpack("<I", body[12:16])[0]
            pkt = body[20: 20 + cap_len]

            if len(pkt) < PKT_DATA_OFF + SRS_HDR_SIZE:
                continue
            if struct.unpack(">H", pkt[12:14])[0] != 0x0800:  # not IPv4
                continue
            if pkt[23] != 17:                                  # not UDP
                continue

            udp_payload = pkt[PKT_DATA_OFF:]
            if len(udp_payload) < SRS_HDR_SIZE:
                continue

            data_id = struct.unpack(">I", udp_payload[4:8])[0]
            if (data_id & 0xFFFFFF00) != SRS_MAGIC:
                continue

            yield udp_payload


def accumulate_histograms(pcapng_path: Path) -> tuple[dict, dict]:
    """
    Single-pass streaming parse.

    Returns
    -------
    noise_hists  : {vmm_id: np.ndarray(1024, int64)}  — over_threshold=0 hits
    signal_hists : {vmm_id: np.ndarray(1024, int64)}  — over_threshold=1 hits
    """
    noise_hists  = {}
    signal_hists = {}

    for payload in iter_srs_payloads(pcapng_path):
        rec = payload[SRS_HDR_SIZE:]
        n_rec = len(rec) // HIT_SIZE

        for i in range(n_rec):
            d1, d2 = struct.unpack(">IH", rec[i * HIT_SIZE: i * HIT_SIZE + HIT_SIZE])

            if not ((d2 >> 15) & 1):       # bit 15 = 0 → marker, skip
                continue

            ot  = (d2 >> 14) & 1
            vmm = (d1 >> 22) & 0x1F
            adc = min((d1 >> 12) & 0x3FF, 1023)

            if vmm not in noise_hists:
                noise_hists[vmm]  = np.zeros(1024, dtype=np.int64)
                signal_hists[vmm] = np.zeros(1024, dtype=np.int64)

            if ot == 0:
                noise_hists[vmm][adc] += 1
            else:
                signal_hists[vmm][adc] += 1

    return noise_hists, signal_hists


# ── Noise and MPV estimation ─────────────────────────────────────────────────

def estimate_noise(hist: np.ndarray,
                   adc_low_cut: int = 20,
                   n_sigma: int = 5) -> dict | None:
    """
    MAD-based noise baseline from ADC histogram (over_threshold=0 hits).

    Removes the VMM3a ADC floor artifact (ADC ≤ adc_low_cut) before
    computing the median and robust sigma.
    """
    adc_vals = np.arange(1024, dtype=np.float64)
    h = hist.astype(float).copy()
    h[: adc_low_cut + 1] = 0
    total = int(h.sum())
    if total < 100:
        return None

    cumsum = np.cumsum(h)
    med_idx = min(int(np.searchsorted(cumsum, total / 2)), 1023)
    median  = adc_vals[med_idx]

    deviations = np.abs(adc_vals - median)
    order      = np.argsort(deviations, kind="stable")
    cumsum_dev = np.cumsum(h[order])
    mad_idx    = min(int(np.searchsorted(cumsum_dev, total / 2)), 1023)
    mad        = float(deviations[order[mad_idx]])
    robust_sigma = 1.4826 * mad
    noise_cut    = median + n_sigma * robust_sigma

    return {
        "n_noise"   : total,
        "median"    : median,
        "sigma"     : robust_sigma,
        "noise_cut" : noise_cut,
    }


def estimate_mpv(signal_hist: np.ndarray,
                 adc_min: float,
                 adc_max: float = 800,
                 smooth_sigma: float = 2) -> tuple:
    """
    Gaussian-smoothed histogram peak finder for Landau MPV.

    Returns (mpv, raw_counts, smoothed_counts, bin_centers).
    mpv is NaN if the histogram is empty in the search window.
    """
    adc_bins = np.arange(1024, dtype=np.float64)
    mask     = (adc_bins >= adc_min) & (adc_bins <= adc_max) & (adc_bins < 1023)

    counts   = signal_hist[mask].astype(float)
    centers  = adc_bins[mask]

    if counts.sum() == 0:
        return np.nan, counts, counts, centers

    smoothed = gaussian_filter1d(counts, sigma=smooth_sigma)
    mpv      = float(centers[np.argmax(smoothed)])
    return mpv, counts, smoothed, centers


# ── Configuration identification ─────────────────────────────────────────────

def load_snr_lookup(csv_path: Path) -> pd.DataFrame | None:
    if not csv_path.exists():
        return None
    df = pd.read_csv(csv_path)
    needed = {"sg", "snt", "vmm_id", "mpv"}
    if not needed.issubset(df.columns):
        return None
    return df


def identify_config(measured_mpv: dict,
                    snr_lookup: pd.DataFrame | None) -> tuple | None:
    """
    Identify the best-matching (sg, snt) by minimising the RMS of the
    relative MPV difference across all detector VMMs.

    Returns (sg, snt, rms_relative_error) or None.
    """
    if snr_lookup is None or not measured_mpv:
        return None

    configs = (
        snr_lookup[["sg", "snt"]]
        .drop_duplicates()
        .values.tolist()
    )
    best_config = None
    best_score  = float("inf")

    for sg, snt in configs:
        df_cfg = snr_lookup[
            (snr_lookup["sg"]  == sg) &
            (snr_lookup["snt"] == snt)
        ]
        sq_errors = []
        for vmm_id, mpv in measured_mpv.items():
            if vmm_id in TRIGGER_VMMS or np.isnan(mpv):
                continue
            row = df_cfg[df_cfg["vmm_id"] == vmm_id]
            if row.empty:
                continue
            ref = float(row["mpv"].iloc[0])
            sq_errors.append(((mpv - ref) / ref) ** 2)

        if not sq_errors:
            continue
        score = float(np.sqrt(np.mean(sq_errors)))
        if score < best_score:
            best_score  = score
            best_config = (sg, snt, score)

    return best_config


# ── Plots ────────────────────────────────────────────────────────────────────

def plot_adc_distributions(noise_hists, signal_hists,
                           noise_params, mpv_results,
                           title_suffix=""):
    """
    One panel per VMM: noise (grey) and signal (blue) ADC histograms,
    with vertical lines for the noise cut and the estimated MPV.
    """
    det_vmms = [v for v in sorted(set(noise_hists) | set(signal_hists))
                if v not in TRIGGER_VMMS]

    ncols = min(4, len(det_vmms))
    nrows = (len(det_vmms) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 5, nrows * 3.8),
                             squeeze=False)
    axes = axes.flatten()

    for ax, vmm_id in zip(axes, det_vmms):
        adc_bins = np.arange(1024)

        nh = noise_hists.get(vmm_id, np.zeros(1024, np.int64))
        sh = signal_hists.get(vmm_id, np.zeros(1024, np.int64))

        if nh.sum() > 0:
            ax.bar(adc_bins, nh, width=1, color="grey",
                   alpha=0.5, label="ot=0 (noise)")
        if sh.sum() > 0:
            ax.bar(adc_bins, sh, width=1, color="steelblue",
                   alpha=0.7, label="ot=1 (signal)")

        np_r = noise_params.get(vmm_id)
        mr   = mpv_results.get(vmm_id, {})
        mpv  = mr.get("mpv", np.nan)

        if np_r:
            ax.axvline(np_r["noise_cut"], color="orange", lw=1.5,
                       ls="--", label=f"noise cut={np_r['noise_cut']:.0f}")
        if not np.isnan(mpv):
            ax.axvline(mpv, color="red", lw=1.8,
                       label=f"MPV={mpv:.0f}")

        det_name = VMM_TO_DETECTOR.get(vmm_id, f"VMM {vmm_id}")
        ax.set_title(f"VMM {vmm_id}  —  {det_name}", fontsize=9)
        ax.set_xlabel("ADC", fontsize=9)
        ax.set_ylabel("Counts", fontsize=9)
        ax.set_xlim(0, 900)
        ax.set_yscale("log")
        ax.legend(fontsize=7, loc="upper right")

    for ax in axes[len(det_vmms):]:
        ax.set_visible(False)

    fig.suptitle(
        f"VMM3a ADC distributions — pcapng{title_suffix}",
        fontsize=13
    )
    plt.tight_layout()
    plt.show()


def plot_mpv_comparison(measured_mpv, snr_lookup, best_config):
    """
    Bar chart comparing measured MPV against the best-matching config
    and the other configs in the lookup table.
    """
    if snr_lookup is None:
        return

    det_vmms = [v for v in sorted(measured_mpv) if v not in TRIGGER_VMMS]
    configs   = (
        snr_lookup[["sg", "snt"]]
        .drop_duplicates()
        .sort_values(["sg", "snt"])
        .values.tolist()
    )

    x = np.arange(len(det_vmms))
    width = 0.8 / (len(configs) + 1)

    fig, ax = plt.subplots(figsize=(max(10, len(det_vmms) * 1.5), 5))

    # Measured values
    measured = [measured_mpv.get(v, np.nan) for v in det_vmms]
    ax.bar(x, measured, width=width * len(configs),
           color="black", alpha=0.85, label="Measured MPV", zorder=3)

    # Reference configs
    cmap = plt.cm.tab10(np.linspace(0, 0.9, len(configs)))
    for k, (sg, snt) in enumerate(configs):
        df_cfg = snr_lookup[
            (snr_lookup["sg"] == sg) & (snr_lookup["snt"] == snt)
        ]
        ref_mpv = [
            float(df_cfg[df_cfg["vmm_id"] == v]["mpv"].iloc[0])
            if not df_cfg[df_cfg["vmm_id"] == v].empty else np.nan
            for v in det_vmms
        ]
        lw = 3 if (best_config and sg == best_config[0] and snt == best_config[1]) else 1
        ls = "-" if (best_config and sg == best_config[0] and snt == best_config[1]) else "--"
        ax.plot(x, ref_mpv, marker="o", lw=lw, ls=ls, color=cmap[k],
                label=f"sg={sg:.1f} snt={snt:.0f}")

    ax.set_xticks(x)
    ax.set_xticklabels([f"VMM {v}" for v in det_vmms], fontsize=10)
    ax.set_ylabel("MPV (ADC)", fontsize=11)
    ax.set_title("Measured MPV vs config scan reference", fontsize=12)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.show()


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    pcapng_file = Path(
        "/drf/projets/clas12/P2/akallits/test_run/"
        "enp4s0f1_20251118-170050_3.pcapng"
    )
    snr_csv = Path(
        "/local/home/ak271430/Documents/PostDocSaclay/"
        "P2_basket_analysis/vmm_snr_results.csv"
    )

    print(f"File : {pcapng_file}")
    print(f"Size : {pcapng_file.stat().st_size / 1e6:.1f} MB")
    print("Parsing …")

    noise_hists, signal_hists = accumulate_histograms(pcapng_file)

    all_vmms = sorted(set(noise_hists) | set(signal_hists))
    print(f"\nVMMs found: {all_vmms}")
    print(f"\n{'VMM':>4}  {'detector':>24}  {'noise hits':>10}  {'signal hits':>11}")
    for vmm in all_vmms:
        nh = int(noise_hists.get(vmm,  np.zeros(1)).sum())
        sh = int(signal_hists.get(vmm, np.zeros(1)).sum())
        print(f"  {vmm:2d}  {VMM_TO_DETECTOR.get(vmm,'?'):>24}  {nh:10d}  {sh:11d}")

    # ── Noise baseline ────────────────────────────────────────────────────────
    noise_params = {}
    for vmm in all_vmms:
        if vmm in TRIGGER_VMMS:
            continue
        res = estimate_noise(noise_hists.get(vmm, np.zeros(1024, np.int64)))
        if res:
            noise_params[vmm] = res

    # ── MPV estimation ────────────────────────────────────────────────────────
    mpv_results  = {}
    measured_mpv = {}
    for vmm in all_vmms:
        if vmm in TRIGGER_VMMS:
            continue
        np_r    = noise_params.get(vmm)
        adc_min = np_r["noise_cut"] if np_r else 60.0
        hist    = signal_hists.get(vmm, np.zeros(1024, np.int64))
        mpv, counts, smoothed, centers = estimate_mpv(hist, adc_min)
        mpv_results[vmm]  = {"mpv": mpv, "counts": counts,
                              "smoothed": smoothed, "centers": centers}
        measured_mpv[vmm] = mpv

    # ── Summary table ─────────────────────────────────────────────────────────
    print(f"\n{'VMM':>4}  {'detector':>24}  "
          f"{'noise σ':>8}  {'cut':>6}  {'MPV':>7}  {'SNR':>6}")
    for vmm in sorted(noise_params):
        det  = VMM_TO_DETECTOR.get(vmm, "?")
        np_r = noise_params[vmm]
        mpv  = measured_mpv.get(vmm, np.nan)
        snr  = mpv / np_r["sigma"] if not np.isnan(mpv) else np.nan
        mpv_s = f"{mpv:.0f}" if not np.isnan(mpv) else "—"
        snr_s = f"{snr:.1f}" if not np.isnan(snr) else "—"
        print(f"  {vmm:2d}  {det:>24}  "
              f"{np_r['sigma']:>8.2f}  "
              f"{np_r['noise_cut']:>6.0f}  "
              f"{mpv_s:>7}  "
              f"{snr_s:>6}")

    # ── Config identification ─────────────────────────────────────────────────
    snr_lookup = load_snr_lookup(snr_csv)
    best       = identify_config(measured_mpv, snr_lookup)

    print()
    if best:
        sg, snt, score = best
        print("┌──────────────────────────────────────────────────────────┐")
        print(f"│  Best matching configuration:  sg={sg:.1f}   snt={snt:.0f} ns")
        print(f"│  RMS relative MPV error: {score:.3f}  ({score*100:.1f}%)")
        print("└──────────────────────────────────────────────────────────┘")

        # Show all configs ranked
        if snr_lookup is not None:
            print("\nAll configurations ranked by MPV agreement:")
            print(f"  {'sg':>5}  {'snt':>6}  {'RMS err %':>10}")
            configs = (
                snr_lookup[["sg", "snt"]]
                .drop_duplicates()
                .sort_values(["sg", "snt"])
                .values.tolist()
            )
            scores = []
            for s_sg, s_snt in configs:
                df_cfg = snr_lookup[
                    (snr_lookup["sg"]  == s_sg) &
                    (snr_lookup["snt"] == s_snt)
                ]
                errs = []
                for vmm_id, mpv in measured_mpv.items():
                    if vmm_id in TRIGGER_VMMS or np.isnan(mpv):
                        continue
                    row = df_cfg[df_cfg["vmm_id"] == vmm_id]
                    if row.empty: continue
                    ref = float(row["mpv"].iloc[0])
                    errs.append(((mpv - ref) / ref) ** 2)
                sc = float(np.sqrt(np.mean(errs))) if errs else np.inf
                scores.append((sc, s_sg, s_snt))
            for sc, s_sg, s_snt in sorted(scores):
                marker = "  ← best match" if (s_sg == sg and s_snt == snt) else ""
                print(f"  {s_sg:>5.1f}  {s_snt:>6.0f}  {sc*100:>10.1f}%{marker}")
    elif snr_lookup is None:
        print(f"Config lookup not found: {snr_csv}")
        print("Cannot identify configuration — raw ADC distributions plotted below.")
    else:
        print("Could not match any configuration.")

    # ── Plots ─────────────────────────────────────────────────────────────────
    plot_adc_distributions(noise_hists, signal_hists, noise_params, mpv_results)
    plot_mpv_comparison(measured_mpv, snr_lookup, best)


if __name__ == "__main__":
    main()
