#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/25/26 12:57 PM
Created in PyCharm
Created as vmm_noise.py

@author: ak271430
"""
"""
vmm_noise.py

Noise baseline estimation for VMM config scan analysis.
Provides per-VMM and per-channel noise characterisation
using robust MAD-based estimators.
"""
import numpy as np
import pandas as pd


def compute_noise_baseline(df_hits, n_sigma=5,
                            adc_low_cut=20,
                            sigma_warn=10.0,
                            sigma_bad=13.0,
                            min_noise_hits=500):
    """
    Compute per-VMM noise pedestal parameters from over_threshold=0 hits.

    Parameters
    ----------
    df_hits : pd.DataFrame
        Hit data with columns: adc, vmm, over_threshold.
    n_sigma : int
        Number of sigma for the noise cut above the median.
    adc_low_cut : int
        Remove hits at or below this ADC value.
        Removes ADC=16 digital artifact (VMM3a ADC floor code).
    sigma_warn : float
        VMMs with robust_sigma above this are flagged as 'warn'.
    sigma_bad : float
        VMMs with robust_sigma above this are flagged as 'bad'.
    min_noise_hits : int
        Minimum number of noise hits for a reliable estimate.

    Returns
    -------
    pd.DataFrame
        One row per VMM with columns:
        vmm_id, n_noise, median_adc, robust_sigma,
        noise_cut, tail_frac_pct, quality
    """
    records = []

    for vmm_id in sorted(df_hits["vmm"].unique()):
        noise = df_hits[
            (df_hits["vmm"]             == vmm_id) &
            (df_hits["over_threshold"]  == 0) &
            (df_hits["adc"]             >  adc_low_cut)
        ]["adc"].values

        if len(noise) < 100:
            continue

        median       = np.median(noise)
        mad          = np.median(np.abs(noise - median))
        robust_sigma = 1.4826 * mad
        cut          = median + n_sigma * robust_sigma
        tail_frac    = (noise > cut).mean() * 100

        if len(noise) < min_noise_hits:
            quality = "bad"
        elif robust_sigma >= sigma_bad:
            quality = "bad"
        elif robust_sigma >= sigma_warn:
            quality = "warn"
        else:
            quality = "ok"

        records.append({
            "vmm_id"      : vmm_id,
            "n_noise"     : len(noise),
            "median_adc"  : median,
            "robust_sigma": robust_sigma,
            "noise_cut"   : cut,
            "tail_frac_pct": tail_frac,
            "quality"     : quality
        })

    return pd.DataFrame(records)


def compute_channel_noise_sigma(noise_adc, use_dither=True,
                                 dither_scale=0.5):
    """
    Compute robust noise sigma for a single channel.

    Applies optional dithering to break integer ADC quantisation
    before computing MAD. This prevents sigma snapping to discrete
    multiples of 1.4826 — a known artifact when noise hits cluster
    tightly within a few ADC counts.

    Adding U(-0.5, 0.5) is mathematically equivalent to assuming
    the true analog value was uniformly distributed within each
    ADC bin, which is the correct physical interpretation of
    quantisation.

    Parameters
    ----------
    noise_adc : np.ndarray
        ADC values of noise hits for this channel.
    use_dither : bool
        Add uniform random noise to break integer clustering.
    dither_scale : float
        Half-width of dither noise in ADC counts.
        0.5 = standard half-bin dithering.

    Returns
    -------
    robust_sigma : float
    median : float
    """
    adc = noise_adc.astype(float)

    if use_dither:
        rng = np.random.default_rng(seed=42)
        adc = adc + rng.uniform(
            -dither_scale, dither_scale, size=len(adc)
        )

    median       = np.median(adc)
    mad          = np.median(np.abs(adc - median))
    robust_sigma = 1.4826 * mad

    return robust_sigma, median

def main():
    print('bonzo')


if __name__ == '__main__':
    main()
