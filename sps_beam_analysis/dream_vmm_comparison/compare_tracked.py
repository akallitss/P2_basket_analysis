#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""compare_tracked.py -- the tracked per-pad ADC (task 1) against the untracked
proxy, and the quintile table the handoff quotes.

Tracked   : win_adc_<station>, the ADC of the nearest in-window unmasked hit of
            a matched track, from the patched eff_autopsy_report.py.
Untracked : adc_vs_ch summed over every capture -- all recorded hits.
"""

import json
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
STATIONS = ["P2_IN", "P2_MID", "P2_OUT"]


def tracked(min_tracks=500, min_adc=50):
    d = json.load(open(os.path.join(DATA, "report_run_46_rate.json")))
    rows = []
    for st in STATIONS:
        pp = d["stations"][st]["per_pad"]
        g = pd.DataFrame({k: pp[k] for k in
                          ("channel_id", "x", "y", "vmm", "ch", "n", "k",
                           "adc_med", "adc_p25", "adc_p75", "adc_med_dnl",
                           "adc_n")})
        g["station"] = st
        rows.append(g)
    t = pd.concat(rows, ignore_index=True)
    t["eff"] = np.where(t.n > 0, t.k / t.n.clip(lower=1), np.nan)
    t["use"] = (t.n >= min_tracks) & (t.adc_n >= min_adc) & t.adc_med.notna()
    return t


def quintiles(g, col="adc_med_dnl"):
    g = g.sort_values(col)
    cw = np.cumsum(g.n.to_numpy())
    cut = np.searchsorted(cw, np.linspace(0, cw[-1], 6)[1:-1])
    out = []
    for p in np.split(np.arange(len(g)), cut):
        if not len(p):
            continue
        q = g.iloc[p]
        out.append((q[col].min(), q[col].max(), q.k.sum() / q.n.sum(),
                    int(q.n.sum()), len(q)))
    return out


def main():
    t = tracked()
    u = pd.read_csv(os.path.join(DATA, "pad_adc_run_46.csv"))
    u = u[["station", "vmm", "ch", "adc_med", "n_hits", "illuminated"]].rename(
        columns={"adc_med": "adc_med_untracked"})
    m = t.merge(u, on=["station", "vmm", "ch"], how="left")

    print("=" * 74)
    print("TRACKED per-pad ADC (task 1 producer), run_46")
    print("=" * 74)
    for st in STATIONS:
        g = m[(m.station == st) & m.use]
        r, p = stats.pearsonr(g.adc_med_dnl, g.eff)
        rho, _ = stats.spearmanr(g.adc_med_dnl, g.eff)
        gg = g.assign(a=g.adc_med_dnl - g.groupby("vmm").adc_med_dnl.transform("mean"),
                      e=g.eff - g.groupby("vmm").eff.transform("mean"))
        rw, pw = stats.pearsonr(gg.a, gg.e)
        print(f"\n--- {st}   pads={len(g)}   tracks={int(g.n.sum())}")
        print(f"    pooled     Pearson r = {r:+.3f} (p={p:.1e})   Spearman {rho:+.3f}")
        print(f"    within chip  Pearson r = {rw:+.3f} (p={pw:.1e})")
        # agreement of the two ADC estimates
        v = g.dropna(subset=["adc_med_untracked"])
        ra, _ = stats.pearsonr(v.adc_med_dnl, v.adc_med_untracked)
        off = (v.adc_med_untracked - v.adc_med_dnl)
        print(f"    tracked vs untracked median: r = {ra:+.3f}   "
              f"offset {off.median():+.1f} ADC (untracked − tracked), "
              f"spread {off.std():.1f}")
        print("      pad median ADC     efficiency    tracks   pads")
        for lo, hi, e, n, np_ in quintiles(g):
            print(f"      {lo:6.0f}-{hi:<6.0f}      {e:.3f}     {n:7d}  {np_:3d}")

    m.to_csv(os.path.join(DATA, "pad_adc_tracked_run_46.csv"), index=False)
    print(f"\nwrote {DATA}/pad_adc_tracked_run_46.csv")


if __name__ == "__main__":
    main()
