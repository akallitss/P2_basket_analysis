#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/25/26 1:02 PM
Created in PyCharm
Created as vmm_signal.py

@author: ak271430
"""
"""
vmm_signal.py

Signal extraction and MPV estimation for VMM config scan analysis.
Provides clean signal selection and Most Probable Value (MPV)
estimation from Landau-distributed ADC distributions.

Key design decisions:
- Saturation cut at ADC=1023 (VMM3a 10-bit ADC hard ceiling)
- MPV estimated via Gaussian-smoothed histogram peak finder
- adc_min for MPV search derived from per-VMM noise_cut
  (ensures signal region starts exactly where noise ends)
- Trigger VMMs (0, 1) excluded from all signal analysis
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d


def get_clean_signal(df_hits, vmm_id,
                      exclude_trigger_vmms=(0, 1)):
    """
    Extract clean signal hits for a detector VMM.

    Applies:
    - over_threshold == 1  (hardware discriminator fired)
    - adc < 1023           (removes saturated hits)
    - excludes trigger VMMs

    Parameters
    ----------
    df_hits : pd.DataFrame
        Hit data with columns: adc, vmm, over_threshold.
    vmm_id : int
        VMM ID to extract signal for.
    exclude_trigger_vmms : tuple
        VMM IDs to skip entirely (trigger VMMs).

    Returns
    -------
    np.ndarray or None
        Clean signal ADC values, or None if VMM is excluded.
    """
    if vmm_id in exclude_trigger_vmms:
        return None

    signal = df_hits[
        (df_hits["vmm"]            == vmm_id) &
        (df_hits["over_threshold"] == 1) &
        (df_hits["adc"]            <  1023)
    ]["adc"].values

    return signal


def estimate_mpv(adc_values, adc_min, adc_max=800,
                  bins=150, smooth_sigma=2):
    """
    Estimate Most Probable Value from a Landau-like ADC distribution.

    Uses a Gaussian-smoothed histogram peak finder to locate the MPV
    robustly against bin-to-bin statistical fluctuations.

    Design notes:
    - adc_min should come from per-VMM noise_cut (not hardcoded)
      to ensure the signal region starts where noise ends
    - bins=150 works well at VMM level (millions of hits)
      use bins=80 at channel level (hundreds of hits)
    - smooth_sigma=2 provides stable peak finding without
      over-smoothing the Landau rising edge

    Parameters
    ----------
    adc_values : np.ndarray
        Clean signal ADC values (saturation already removed).
    adc_min : float
        Lower bound for MPV search — use per-VMM noise_cut.
    adc_max : float
        Upper bound — well below saturation.
    bins : int
        Histogram bins. Use 150 for VMM level, 80 for channel level.
    smooth_sigma : float
        Gaussian smoothing width in bins.

    Returns
    -------
    mpv : float
        Most Probable Value in ADC.
    counts : np.ndarray
        Raw histogram counts.
    smoothed : np.ndarray
        Gaussian-smoothed counts.
    bin_centers : np.ndarray
        Bin center positions.
    """
    counts, edges = np.histogram(
        adc_values, bins=bins, range=(adc_min, adc_max)
    )
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed    = gaussian_filter1d(counts.astype(float),
                                     sigma=smooth_sigma)
    mpv         = bin_centers[np.argmax(smoothed)]

    return mpv, counts, smoothed, bin_centers

def main():
    print('bonzo')


if __name__ == '__main__':
    main()
