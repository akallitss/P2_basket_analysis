#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pad_adc.py -- per-pad pulse height vs per-pad VMM efficiency, run_46.

Task 1 of HANDOFF_VMM_PAD_ADC.md asks for a producer for the per-pad ADC /
efficiency correlation. This is the *untracked* route to it, and it needs no
pcapng decoding at all:

  * per-pad efficiency  (n, k) comes from the committed efficiency_v2 JSON,
    which already stores `per_pad` and reproduces the three headline numbers;
  * per-pad pulse height comes from `adc_vs_ch` in every capture's counts.npz,
    a (32 vmm, 64 ch, 128 adc-bin) additive histogram with bin = adc >> 3.
    vmm_reduce.py writes it whether or not --drop-columns was used, so it
    exists for all 48 captures of run_46 where the .npy columns do not.

The two are joined on (vmm, ch): `per_pad` is stored in (vmm, ch) sorted order,
which this script asserts rather than assumes.

WHAT THIS IS NOT.  `adc_vs_ch` counts every recorded hit on a channel, not the
in-window hit of a matched track (`win_adc_<station>` in the autopsy parquet).
It is a larger, independent sample -- 1.52e9 hits against 837k tracks -- but it
includes off-track beam. See compare_tracked.py for the tracked cross-check.

DNL.  The VMM3a ADC is ~2x too wide on every multiple of 16. adc_vs_ch bins by
8, so one DNL period is exactly 2 bins: rebinning by 2 gives 16-code bins each
containing exactly one full period, and the comb cancels. All quoted statistics
use that rebinning unless stated.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/dylan/PycharmProjects/DAQ_Control_VMM_Beam/vmm_qa")
import vmm_stations as VS                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

ADC_BIN = 8          # adc_vs_ch bin width, from vmm_reduce.ADC_CH_ADC_BINS
DNL_PERIOD = 16      # VMM3a ADC differential non-linearity period, in codes


def load_inputs(run="run_46", sub="cfg_gain4.5_peaktime200"):
    z = np.load(os.path.join(DATA, f"adc_vs_ch_{run}.npz"))
    eff = json.load(open(os.path.join(DATA, f"vmm_eff_{run}_{sub}.json")))
    return z, eff


def rebin(h, factor):
    """Sum adjacent bins. h is (..., nbins)."""
    n = (h.shape[-1] // factor) * factor
    return h[..., :n].reshape(*h.shape[:-1], n // factor, factor).sum(-1)


def quantile_from_hist(h, edges, q):
    """Interpolated quantile of a histogram. h (npads, nbins), edges (nbins+1)."""
    c = np.cumsum(h, axis=1, dtype=float)
    tot = c[:, -1:]
    out = np.full(len(h), np.nan)
    ok = tot[:, 0] > 0
    if not ok.any():
        return out
    f = np.where(tot > 0, c / np.where(tot > 0, tot, 1), 0.0)
    for i in np.flatnonzero(ok):
        j = int(np.searchsorted(f[i], q, side="left"))
        j = min(j, len(edges) - 2)
        below = f[i, j - 1] if j > 0 else 0.0
        frac = (q - below) / max(f[i, j] - below, 1e-12)
        out[i] = edges[j] + frac * (edges[j + 1] - edges[j])
    return out


def build(run="run_46", sub="cfg_gain4.5_peaktime200", min_tracks=500):
    z, eff = load_inputs(run, sub)
    H = z["adc_vs_ch"]                       # (32, 64, 128)
    tab = VS.build_pad_table()

    # pad_area is in the mapping CSV but not in build_pad_table's column list
    area = pd.read_csv(VS.MAP_CSV).set_index("channel_id")["pad_area"]
    tab["pad_area"] = tab["channel_id"].map(area)

    rows = []
    for st in eff["stations"]:
        name = st["station"]
        s = tab[tab["station"] == name].sort_values(["vmm", "ch"]).reset_index(drop=True)
        pp = st["per_pad"]
        x = np.asarray(pp["x"], float)
        y = np.asarray(pp["y"], float)
        # the JSON rounds x/y to 2 dp; assert the (vmm, ch) ordering holds
        assert len(x) == len(s), f"{name}: {len(x)} per_pad vs {len(s)} pads"
        dx = np.abs(s["pad_cx"].to_numpy() - x)
        dy = np.abs(s["pad_cy"].to_numpy() - y)
        assert np.nanmax(dx) < 0.01 and np.nanmax(dy) < 0.01, \
            f"{name}: per_pad is not in (vmm, ch) order (max dx={np.nanmax(dx)})"

        vmm = s["vmm"].to_numpy()
        ch = s["ch"].to_numpy()
        h = H[vmm, ch].astype(float)                 # (384, 128)

        mask = {int(v): set(c) for v, c in st["masking"]["mask"].items()}
        masked = np.array([c in mask.get(v, ()) for v, c in zip(vmm, ch)])

        # DNL-clean: rebin 8-code bins by 2 -> 16-code bins = one full period
        h16 = rebin(h, DNL_PERIOD // ADC_BIN)
        edges16 = np.arange(h16.shape[1] + 1) * DNL_PERIOD

        med = quantile_from_hist(h16, edges16, 0.50)
        p25 = quantile_from_hist(h16, edges16, 0.25)
        p75 = quantile_from_hist(h16, edges16, 0.75)
        # raw (8-code) median, to show the DNL is not driving the result
        edges8 = np.arange(h.shape[1] + 1) * ADC_BIN
        med8 = quantile_from_hist(h, edges8, 0.50)
        # MPV on the DNL-clean binning
        mpv = np.where(h16.sum(1) > 0,
                       edges16[np.argmax(h16, axis=1)] + DNL_PERIOD / 2.0,
                       np.nan)

        n = np.asarray(pp["n"], float)
        k = np.asarray(pp["k"], float)
        rows.append(pd.DataFrame({
            "station": name, "vmm": vmm, "ch": ch,
            "channel_id": s["channel_id"].to_numpy(),
            "x": s["pad_cx"].to_numpy(), "y": s["pad_cy"].to_numpy(),
            "radius": s["radius"].to_numpy(), "phi": s["phi"].to_numpy(),
            "pad_area": s["pad_area"].to_numpy(),
            "connector_N": s["connector_N"].to_numpy(),
            "half": s["half"].to_numpy(),
            "n_track": n, "k_track": k,
            "eff": np.where(n > 0, k / np.where(n > 0, n, 1), np.nan),
            "n_hits": h.sum(1),
            "adc_med": med, "adc_p25": p25, "adc_p75": p75,
            "adc_med_raw8": med8, "adc_mpv": mpv,
            "masked": masked,
        }))

    df = pd.concat(rows, ignore_index=True)
    df["illuminated"] = (df["n_track"] >= min_tracks) & (~df["masked"])
    # binomial error on the efficiency
    df["eff_err"] = np.sqrt(np.clip(df["eff"] * (1 - df["eff"]), 0, None)
                            / df["n_track"].clip(lower=1))
    return df, eff


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default="run_46")
    ap.add_argument("--sub", default="cfg_gain4.5_peaktime200")
    ap.add_argument("--min-tracks", type=int, default=500)
    ap.add_argument("-o", "--out", default=os.path.join(DATA, "pad_adc_run_46.csv"))
    a = ap.parse_args()

    df, _ = build(a.run, a.sub, a.min_tracks)
    df.to_csv(a.out, index=False)
    print(f"wrote {a.out}  ({len(df)} pads, "
          f"{int(df['illuminated'].sum())} illuminated)")
    for st, g in df.groupby("station"):
        i = g[g["illuminated"]]
        print(f"  {st:8s} illuminated={len(i):3d}  "
              f"eff={i['k_track'].sum()/i['n_track'].sum():.4f}  "
              f"median ADC {i['adc_med'].min():.0f}..{i['adc_med'].max():.0f}")


if __name__ == "__main__":
    main()
