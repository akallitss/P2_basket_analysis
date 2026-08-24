#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_vmm_adc_hist.py -- the track-matched VMM pulse-height spectrum, as an npz.

This is the missing input to mpgd2026/make_adc_comparison.py, whose panel (a)
has always taken --vmm-hist and has never been given one: the docstring there
says "it is not on this machine (it needs the autopsy parquet from EOS), so
without it panel (a) draws the VMM window from the measured quantiles" -- i.e.
two vertical lines instead of a distribution. With this npz the panel becomes a
real DREAM-vs-VMM overlay and the discriminator cut-off is visible rather than
asserted.

What "track-matched" means here, and why it has to be: the counts store holds
EVERY recorded hit and on run_46 that is dominated by noise (P2_IN alone is
1.14e9 hits in one narrow spike). Overlaying that on a DREAM cluster-charge
Landau would compare a noise spectrum with a signal spectrum. So the selection
is the same one the efficiency is defined by -- a uRWELL track landing on the
instrumented pads, inside a capture, with a hit within probe_r -- and the value
is that hit's ADC:

    fiducial : dpad_<st>  < fid_r      (nearest instrumented pad centre)
    live     : in_capture
    matched  : win_dmin_<st> < probe_r (nearest unmasked in-window hit)
    value    : win_adc_<st>            (that hit's per-pad ADC)

Run it where the autopsy parquet lives (lxplus/EOS):

    python3 make_vmm_adc_hist.py autopsy_run46_tracks.parquet \\
        --station P2_OUT -o vmm_adc_P2_OUT_run46.npz

@author: ak271430 Alexandra Kallitsopoulou
"""

import argparse
import json
import os

import numpy as np
import pandas as pd

ADC_MAX = 1024

# The VMM 10-bit ADC has strong differential non-linearity with a period of
# exactly 16 codes: every multiple of 16 is a coarse-bit boundary and collects
# 2-3x its neighbours. MEASURED on run_46: 13.5% (P2_OUT) and 14.9% (P2_MID) of
# all hits land exactly on a multiple of 16, against 6.2% if the codes were
# used evenly. A unit-bin mode therefore returns the tallest comb tooth, not
# the Landau peak -- which is how the autopsy came to quote "most probable 128"
# for P2_OUT (8x16) and 80 for P2_MID (5x16). Rebinning by the DNL period is
# what makes the peak an estimate of the signal instead of the ADC.
DNL_PERIOD = 16


def build(tracks, station, fid_r, probe_r, bins):
    need = [f"dpad_{station}", f"win_dmin_{station}", f"win_adc_{station}"]
    missing = [c for c in need if c not in tracks.columns]
    if missing:
        raise SystemExit(f"parquet has no {missing} -- is this the tracks "
                         f"file, and was {station} in the autopsy run?")

    fid = tracks[f"dpad_{station}"].to_numpy() < fid_r
    if "in_capture" in tracks.columns:
        fid &= tracks["in_capture"].to_numpy().astype(bool)
    hit = fid & (tracks[f"win_dmin_{station}"].to_numpy() < probe_r)

    adc = tracks[f"win_adc_{station}"].to_numpy()[hit].astype(float)
    # adc == 0 is the fill value for "no hit recorded", never a real pulse
    adc = adc[adc > 0]

    counts, edges = np.histogram(adc, bins=bins, range=(0, ADC_MAX))
    return counts, edges, adc, int(fid.sum()), int(hit.sum())


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tracks_parquet")
    ap.add_argument("--station", default="P2_OUT")
    ap.add_argument("--fid-r", type=float, default=9.0)
    ap.add_argument("--probe-r", type=float, default=15.0)
    # Unit ADC bins. The VMM ADC is 10-bit, so this is its native
    # resolution, and it is also what makes the mode here agree with the
    # autopsy's np.bincount().argmax(): at 4-wide bins P2_OUT peaks at 98,
    # at unit bins it peaks at 128, which is the published number.
    ap.add_argument("--bins", type=int, default=ADC_MAX)
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()

    tracks = pd.read_parquet(a.tracks_parquet)
    counts, edges, adc, n_fid, n_hit = build(
        tracks, a.station, a.fid_r, a.probe_r, a.bins)

    q = {f"p{p:02d}": float(np.percentile(adc, p)) for p in (1, 5, 25, 50, 95)}
    ctr = 0.5 * (edges[:-1] + edges[1:])
    mpv = float(ctr[int(np.argmax(counts))])
    lowest = float(adc.min())

    # DNL-robust peak: collapse the comb, then take the mode.
    mpv_r, on_teeth = mpv, None
    if len(counts) == ADC_MAX:
        n16 = (ADC_MAX // DNL_PERIOD) * DNL_PERIOD
        cb = counts[:n16].reshape(-1, DNL_PERIOD).sum(axis=1)
        mpv_r = float((np.argmax(cb) + 0.5) * DNL_PERIOD)
        on_teeth = float(counts[DNL_PERIOD::DNL_PERIOD].sum() / counts.sum())

    out = a.out or f"vmm_adc_{a.station}.npz"
    np.savez_compressed(
        out, counts=counts.astype(np.int64), edges=edges,
        meta=json.dumps(dict(
            station=a.station, source=os.path.abspath(a.tracks_parquet),
            fid_r=a.fid_r, probe_r=a.probe_r,
            n_fiducial_tracks=n_fid, n_matched=n_hit, n_in_hist=int(adc.size),
            lowest_adc=lowest, mpv_adc=mpv, mpv_adc_dnl_robust=mpv_r,
            frac_on_dnl_teeth=on_teeth, dnl_period=DNL_PERIOD,
            quantiles=q)))

    print(f"{a.station}: {n_fid:,} fiducial tracks, {n_hit:,} matched, "
          f"{adc.size:,} in the histogram")
    print(f"  lowest = {lowest:.0f}   5th pct = {q['p05']:.0f}   "
          f"median = {q['p50']:.0f}")
    print(f"  MPV: unit-bin mode {mpv:.0f} (DNL tooth if a multiple of "
          f"{DNL_PERIOD}), DNL-robust {mpv_r:.0f}"
          + (f"; {on_teeth:.1%} of hits on teeth" if on_teeth else ""))
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
