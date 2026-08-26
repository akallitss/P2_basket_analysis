#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eff_autopsy_report.py -- read what eff_autopsy.py dumped and answer the question.

Every number here is a decomposition of the SAME efficiency the production
script quotes, so the pieces add up to it by construction:

    fiducial track          the prediction lands on an instrumented pad and the
                            VMM was recording
    -> found                a pad hit within probe_r in the coincidence window
    -> missed, and then exactly one of
         no hit at all      the station produced nothing in +-10 us
         out of window      it produced something, but not in coincidence
         masked away        it produced a hit on the right pad, and the
                            hot-channel mask removed it
         wrong place        a hit in the window, further than probe_r

Then the same efficiency sliced by everything that could be causing it: pad,
VMM chip (each has its own threshold DAC), position in the spill, instantaneous
rate, time since the previous trigger, uRWELL track quality, pulse height, and
whether the OTHER two stations saw the same track.

    python3 eff_autopsy_report.py autopsy/autopsy_run_46_cfg_..._tracks.parquet \
        --json autopsy_report.json --figdir figs
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vmm_stations as VS                 # noqa: E402

STATIONS = ["P2_IN", "P2_MID", "P2_OUT"]
# per-VMM threshold DAC of the sub_run's configuration file (see the config
# txt next to the raw data); filled from --thresholds if given.
DEFAULT_SDT = {}


def clopper_pearson(k, n, cl=0.6827):
    from scipy.stats import beta
    if n == 0:
        return 0.0, 1.0
    a = (1 - cl) / 2
    lo = 0.0 if k == 0 else beta.ppf(a, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a, k + 1, n - k)
    return float(lo), float(hi)


def eff_of(mask_fid, mask_hit):
    n = int(mask_fid.sum())
    k = int((mask_fid & mask_hit).sum())
    lo, hi = clopper_pearson(k, n)
    return {"eff": (k / n if n else None), "k": k, "n": n, "lo": lo, "hi": hi}


def sliced(fid, hit, value, edges, label):
    """Efficiency in bins of `value` -- the generic 'does it depend on' test."""
    out = []
    b = np.digitize(value, edges)
    for i in range(1, len(edges)):
        m = fid & (b == i)
        r = eff_of(m, hit)
        r["lo_edge"], r["hi_edge"] = float(edges[i - 1]), float(edges[i])
        out.append(r)
    return {"variable": label, "bins": out}


# --------------------------------------------------------------------------- #
def channel_mask(meta, st, how, max_hz=20000.0):
    """{vmm: [ch]} under one of the three rules, from the stored occupancy."""
    occ = np.array(meta["occupancy"])
    vmms = VS.STATION_VMMS[st]
    if how == "ratio":                       # what the production run applied
        return {int(k): v for k, v in meta["stations"][st]["mask_prod"].items()}
    known = {v: c for v, c in VS.NOISY_CHANNELS.items() if v in vmms}
    if how == "none":
        return known
    # wall clock of the sub_run: the stored capture spans start at t=0
    # (a few hits per capture carry srs_timestamp = 0), so summing them is
    # 25x too long -- take the span the captures actually cover instead
    cov = meta.get("coverage") or []
    dur = ((max(hi for _, hi in cov) - min(lo for lo, _ in cov)) / 1e9
           if cov else 0.0)
    out = dict(known)
    if dur > 0:
        for v in vmms:
            hot = np.flatnonzero(occ[v] / dur > max_hz)
            if hot.size:
                out[int(v)] = sorted(set(out.get(int(v), []))
                                     | {int(x) for x in hot})
    return out


def reduce_hits(hits, n, keep):
    """Per track: distance to the nearest kept hit, and how many there were."""
    d = np.full(n, np.inf, np.float32)
    cnt = np.zeros(n, np.int32)
    sub = hits[keep]
    if len(sub):
        jj = sub["idx"].to_numpy(np.int64)
        dd = sub["dpred"].to_numpy(np.float32)
        np.add.at(cnt, jj, 1)
        ok = np.isfinite(dd)
        if ok.any():
            np.minimum.at(d, jj[ok], dd[ok])
    return d, cnt


DNL_PERIOD = 16      # VMM3a ADC differential non-linearity period, in codes


def pad_pulse_height(tracks, st, j, sel, n_pads):
    """Per-pad pulse height of the tracks that WERE detected on that pad.

    `win_adc_{st}` is the ADC of the nearest in-window unmasked hit, and is 0
    where the track produced none -- so the pulse height is defined on the k
    detected tracks, not on all n.  Returns arrays aligned with `ids`.

    Two medians are stored.  `adc_med` is the plain one.  `adc_med_dnl` is the
    same quantity read off a histogram rebinned to 16 codes: the VMM3a ADC is
    ~2x too wide on every multiple of 16, and a plain median lands on a comb
    tooth wherever the distribution is locally flat.  The two agree to a few
    codes on a well-populated pad and diverge on a sparse one, which is the
    useful signal -- quote `adc_med_dnl`, and use the pair as a warning flag.

    That 16-code histogram is kept as `adc_hist` (sparse: only pads that saw
    something), so the per-pad pulse-height DISTRIBUTION can be drawn, not just
    its median.  16 codes is the natural bin here for the same reason it fixes
    the median -- it is exactly one DNL period, so the comb averages out inside
    a bin instead of rippling across the plot.
    """
    adc = tracks[f"win_adc_{st}"].to_numpy().astype(np.float64)
    m = sel & (adc > 0)
    out = {k: np.full(n_pads, np.nan) for k in
           ("adc_med", "adc_p25", "adc_p75", "adc_med_dnl")}
    out["adc_n"] = np.zeros(n_pads, np.int64)
    if not m.any():
        d = {k: v.tolist() for k, v in out.items()}
        d.update(adc_hist={}, adc_hist_bin=DNL_PERIOD)
        return d

    df = pd.DataFrame({"j": j[m], "adc": adc[m]})
    g = df.groupby("j")["adc"]
    for name, q in (("adc_med", 0.5), ("adc_p25", 0.25), ("adc_p75", 0.75)):
        s = g.quantile(q)
        out[name][s.index.to_numpy()] = s.to_numpy()
    c = g.size()
    out["adc_n"][c.index.to_numpy()] = c.to_numpy()

    # DNL-clean median: interpolate the 50% point of a 16-code histogram
    nb = 1024 // DNL_PERIOD
    b = np.minimum((df["adc"].to_numpy() // DNL_PERIOD).astype(np.int64), nb - 1)
    h = np.zeros((n_pads, nb), np.int64)
    np.add.at(h, (df["j"].to_numpy(), b), 1)
    tot_h = h.sum(1)
    for i in np.flatnonzero(tot_h > 0):
        cum = np.cumsum(h[i]) / tot_h[i]
        k = int(np.searchsorted(cum, 0.5, side="left"))
        below = cum[k - 1] if k > 0 else 0.0
        frac = (0.5 - below) / max(cum[k] - below, 1e-12)
        out["adc_med_dnl"][i] = (k + frac) * DNL_PERIOD

    d = {k: v.tolist() for k, v in out.items()}
    d["adc_hist"] = {int(i): h[i].tolist() for i in np.flatnonzero(tot_h > 0)}
    d["adc_hist_bin"] = DNL_PERIOD
    return d


def analyse(tracks, hits, meta, st, args):
    pads = VS.build_pad_table()
    pads = pads[(pads["station"] == st) & pads["pad_cx"].notna()]
    sdt = args.sdt.get(st, {})

    px = tracks[f"px_{st}"].to_numpy()
    py = tracks[f"py_{st}"].to_numpy()
    dpad = tracks[f"dpad_{st}"].to_numpy()
    nearpad = tracks[f"nearpad_{st}"].to_numpy()
    fid = (dpad < args.fid_r) & tracks["in_capture"].to_numpy()

    # The mask is re-applied here rather than taken from the pass, so that the
    # three rules can be compared on one event sample instead of on three runs.
    mask = channel_mask(meta, st, args.mask, args.max_hz)
    if hits is not None and len(hits):
        flag = np.zeros(len(hits), bool)
        hv = hits["vmm"].to_numpy()
        hc = hits["ch"].to_numpy()
        for v, chans in mask.items():
            flag |= (hv == v) & np.isin(hc, chans)
        hits = hits.assign(masked=flag)
        n = len(tracks)
        dmin, _ = reduce_hits(hits, n, hits["in_win"].to_numpy() & ~flag)
        dmin_all, _ = reduce_hits(hits, n, hits["in_win"].to_numpy())
        d_any_new, _ = reduce_hits(hits, n, ~flag)
        tracks[f"win_dmin_{st}"] = dmin          # in place: the figures read
        tracks[f"winall_dmin_{st}"] = dmin_all    # the same frame afterwards
        tracks[f"any_dmin_{st}"] = d_any_new
    dmin = tracks[f"win_dmin_{st}"].to_numpy()
    dmin_all = tracks[f"winall_dmin_{st}"].to_numpy()
    hit = dmin < args.probe_r
    hit_all = dmin_all < args.probe_r
    n_any = tracks[f"nraw_any_{st}"].to_numpy()
    n_win = tracks[f"nraw_win_{st}"].to_numpy()

    res = {"station": st, "n_tracks": int(len(tracks)),
           "n_fiducial": int(fid.sum()),
           "mask_rule": args.mask,
           "mask": {int(k): list(v) for k, v in mask.items()},
           "efficiency": eff_of(fid, hit),
           "efficiency_no_mask": eff_of(fid, hit_all)}

    # ---- where the misses go ------------------------------------------------
    # "did the station produce anything" must be asked about the PLACE the
    # track went, not about the station as a whole: with a channel running at
    # 100 kHz almost every event has some hit somewhere in +-10 us, and a
    # budget built on that counts noise as a response.
    d_any = tracks[f"any_dmin_{st}"].to_numpy()
    miss = fid & ~hit
    masked_away = miss & hit_all
    out_of_time = miss & ~masked_away & (d_any < args.probe_r)
    nothing = miss & ~masked_away & ~out_of_time
    res["loss_budget"] = {
        "n_fiducial": int(fid.sum()),
        "n_miss": int(miss.sum()),
        "masked_away": int(masked_away.sum()),
        "hit_at_the_place_but_out_of_time": int(out_of_time.sum()),
        "nothing_within_probe_r_at_any_time": int(nothing.sum()),
        "of_which_station_totally_silent": int((nothing & (n_any == 0)).sum()),
        "n_hits_in_window_when_missing": float(n_win[miss].mean()),
        "comment": "out-of-time is an upper limit: a 10 us window around a "
                   "15 mm radius also catches accidentals"}

    # ---- the pad the track points at ---------------------------------------
    pad_vmm = dict(zip(pads["channel_id"], pads["vmm"]))
    pad_ch = dict(zip(pads["channel_id"], pads["ch"]))
    vmm_of_track = np.array([pad_vmm.get(int(p), -1) for p in nearpad])
    mask_prod = meta["stations"][st]["mask_prod"]
    masked_pad = np.array([
        int(pad_ch.get(int(p), -1)) in mask_prod.get(str(pad_vmm.get(int(p), -1)),
                                                     mask_prod.get(
                                                         pad_vmm.get(int(p), -1),
                                                         []))
        for p in nearpad])
    # kept under the ratio rule whatever --mask is set to: this IS the
    # diagnostic -- how much beam the old rule was deleting
    res["pointing_at_a_pad_the_ratio_mask_kills"] = {
        "n": int((fid & masked_pad).sum()),
        "frac_of_fiducial": float((fid & masked_pad).mean()),
        "efficiency_there": eff_of(fid & masked_pad, hit),
        "efficiency_there_no_mask": eff_of(fid & masked_pad, hit_all),
        "efficiency_elsewhere": eff_of(fid & ~masked_pad, hit)}

    # ---- per VMM chip -------------------------------------------------------
    occ = np.array(meta["occupancy"])
    per_vmm = []
    for v in VS.STATION_VMMS[st]:
        m = fid & (vmm_of_track == v)
        r = eff_of(m, hit)
        r.update(vmm=int(v), sdt=sdt.get(str(v), sdt.get(v)),
                 occupancy=int(occ[v].sum()),
                 loudest_channel=int(np.argmax(occ[v])),
                 loudest_hits=int(occ[v].max()),
                 n_dead_channels=int((occ[v] == 0).sum()),
                 eff_no_mask=eff_of(m, hit_all)["eff"])
        per_vmm.append(r)
    res["per_vmm"] = per_vmm

    # ---- per pad, for the map ----------------------------------------------
    ids = pads["channel_id"].to_numpy()
    order = {int(c): i for i, c in enumerate(ids)}
    j = np.array([order.get(int(p), -1) for p in nearpad])
    ok = fid & (j >= 0)
    tot = np.bincount(j[ok], minlength=len(ids))
    got = np.bincount(j[ok & hit], minlength=len(ids))
    res["per_pad"] = {"channel_id": ids.tolist(),
                      "x": pads["pad_cx"].round(2).tolist(),
                      "y": pads["pad_cy"].round(2).tolist(),
                      "vmm": pads["vmm"].tolist(), "ch": pads["ch"].tolist(),
                      "n": tot.tolist(), "k": got.tolist()}
    res["per_pad"].update(pad_pulse_height(tracks, st, j, ok & hit, len(ids)))
    live = tot >= args.min_pad
    e_pad = np.where(live, got / np.maximum(tot, 1), np.nan)
    res["pad_summary"] = {
        "n_pads_with_tracks": int(live.sum()),
        "n_pads_below_5pct": int(np.nansum(e_pad[live] < 0.05)),
        "frac_tracks_on_pads_below_5pct": float(
            tot[live][e_pad[live] < 0.05].sum() / max(tot[live].sum(), 1)),
        "median_pad_efficiency": float(np.nanmedian(e_pad[live])),
        "p90_pad_efficiency": float(np.nanpercentile(e_pad[live], 90))}

    # ---- does it depend on anything it should not? -------------------------
    res["vs"] = {}
    res["vs"]["t_in_spill"] = sliced(
        fid, hit, tracks["t_in_spill_ms"].to_numpy(),
        np.array([0, 250, 500, 1000, 1500, 2000, 3000, 4000, 6000, 1e9]),
        "time into the spill [ms]")
    res["vs"]["local_rate"] = sliced(
        fid, hit, tracks["n_trig_1ms"].to_numpy(),
        np.array([0, 1, 2, 3, 5, 8, 12, 20, 1e9]),
        "triggers within +-1 ms")
    res["vs"]["dt_prev_trigger"] = sliced(
        fid, hit, tracks["dt_prev_trig_us"].to_numpy(),
        np.array([0, 10, 25, 50, 100, 250, 500, 1000, 1e9]),
        "time since the previous trigger [us]")
    for ax in ("dx", "dy"):
        res["vs"][f"track_{ax}"] = sliced(
            fid, hit, tracks[ax].abs().to_numpy(),
            np.array([0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]),
            f"|front-back {ax}| [mm]")
    res["vs"]["n_hits_station"] = sliced(
        fid, hit, n_any, np.array([0, 1, 2, 3, 5, 10, 1e9]),
        "station hits in +-10 us")
    # Where inside its pad did the particle go? A threshold-limited readout
    # loses the edges and the corners first: the charge splits between two
    # pads and neither half clears the discriminator. A readout that is losing
    # whole events for some other reason is flat across the cell.
    res["vs"]["dist_to_pad_centre"] = sliced(
        fid, hit, dpad, np.arange(0, 9.5, 1.0),
        "distance from the track to the centre of its pad [mm]")

    # the same thing in two dimensions: every pad folded onto one cell
    cx = dict(zip(pads["channel_id"], pads["pad_cx"]))
    cy = dict(zip(pads["channel_id"], pads["pad_cy"]))
    ox = px - np.array([cx.get(int(p), np.nan) for p in nearpad])
    oy = py - np.array([cy.get(int(p), np.nan) for p in nearpad])
    good = fid & np.isfinite(ox) & np.isfinite(oy)
    b = np.linspace(-7, 7, 15)
    tot2, _, _ = np.histogram2d(ox[good], oy[good], bins=[b, b])
    got2, _, _ = np.histogram2d(ox[good & hit], oy[good & hit], bins=[b, b])
    res["intra_pad"] = {"edges": b.tolist(), "n": tot2.tolist(),
                        "k": got2.tolist()}

    # ---- amplitude ----------------------------------------------------------
    if hits is not None and len(hits):
        # when the track is missed, is there anything at that place at ANY
        # time? A flat dt means the answer is "only accidentals".
        near = hits[hits["dpred"] < args.probe_r]
        is_miss = miss[near["idx"].to_numpy()]
        h_miss, ed = np.histogram(near["dt"].to_numpy()[is_miss],
                                  bins=np.arange(-2000, 8001, 250))
        h_hit, _ = np.histogram(near["dt"].to_numpy()[~is_miss],
                                bins=np.arange(-2000, 8001, 250))
        res["dt_at_the_track"] = {
            "edges_ns": ed.tolist(), "missed_tracks": h_miss.tolist(),
            "found_tracks": h_hit.tolist()}

        on = hits[(hits["in_win"]) & (hits["dpred"] < args.probe_r)
                  & (~hits["masked"])]
        adc = on["adc"].to_numpy()
        res["adc"] = {
            "n": int(len(adc)),
            "percentiles": {str(p): float(np.percentile(adc, p))
                            for p in (1, 5, 10, 25, 50, 75, 90, 99)},
            "min": int(adc.min()) if len(adc) else None,
            "mode": int(np.bincount(adc.clip(0, 1023)).argmax()) if len(adc)
            else None}
        # efficiency if the threshold had been higher: cut the found hits
        thr_scan = {}
        idx = on["idx"].to_numpy()
        for t in (0, 50, 100, 150, 200, 250, 300, 400, 500):
            keep = np.zeros(len(tracks), bool)
            keep[idx[adc >= t]] = True
            thr_scan[str(t)] = eff_of(fid, keep)["eff"]
        res["adc"]["efficiency_vs_offline_threshold"] = thr_scan
        hh, he = np.histogram(adc, bins=np.arange(0, 1041, 10))
        res["adc"]["hist"] = {"edges": he.tolist(), "n": hh.tolist()}
        res["adc"]["hist_per_vmm"] = {
            str(int(v)): np.histogram(g["adc"].to_numpy(),
                                      bins=np.arange(0, 1041, 10))[0].tolist()
            for v, g in on.groupby("vmm") if len(g) > 200}
        res["adc"]["per_vmm_percentiles"] = {
            str(int(v)): {str(p): float(np.percentile(g["adc"], p))
                          for p in (1, 5, 50, 95)}
            for v, g in on.groupby("vmm") if len(g) > 200}
    return res


# --------------------------------------------------------------------------- #
def cross_station(tracks, args, res):
    """Are the three stations missing the same events?"""
    out = {}
    have = {}
    for st in STATIONS:
        if f"win_dmin_{st}" not in tracks:
            continue
        have[st] = (tracks[f"win_dmin_{st}"].to_numpy() < args.probe_r)
    fid = {st: ((tracks[f"dpad_{st}"].to_numpy() < args.fid_r)
                & tracks["in_capture"].to_numpy())
           for st in have}
    for a in have:
        for b in have:
            if a >= b:
                continue
            m = fid[a] & fid[b]
            pa, pb = have[a][m].mean(), have[b][m].mean()
            both = (have[a] & have[b])[m].mean()
            out[f"{a}|{b}"] = {
                "eff_a": float(pa), "eff_b": float(pb),
                "both": float(both), "independent_prediction": float(pa * pb),
                "eff_a_given_b": float((have[a] & have[b])[m].sum()
                                       / max(have[b][m].sum(), 1)),
                "eff_a_given_not_b": float((have[a] & ~have[b])[m].sum()
                                           / max((~have[b])[m].sum(), 1)),
                "n": int(m.sum())}
    res["cross_station"] = out


# --------------------------------------------------------------------------- #
def figures(tracks, meta, res, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for st in STATIONS:
        if st not in res["stations"]:
            continue
        r = res["stations"][st]
        px = tracks[f"px_{st}"].to_numpy()
        py = tracks[f"py_{st}"].to_numpy()
        fid = ((tracks[f"dpad_{st}"].to_numpy() < args.fid_r)
               & tracks["in_capture"].to_numpy())
        hit = tracks[f"win_dmin_{st}"].to_numpy() < args.probe_r
        pp = r["per_pad"]
        pn, pk = np.array(pp["n"]), np.array(pp["k"])
        live = pn >= args.min_pad
        e_pad = np.where(live, pk / np.maximum(pn, 1), np.nan)

        fig, ax = plt.subplots(2, 4, figsize=(25, 11))
        bs = args.bin
        xe = np.arange(px[fid].min(), px[fid].max() + bs, bs)
        ye = np.arange(py[fid].min(), py[fid].max() + bs, bs)
        tot, _, _ = np.histogram2d(px[fid], py[fid], bins=[xe, ye])
        ok, _, _ = np.histogram2d(px[fid & hit], py[fid & hit], bins=[xe, ye])
        with np.errstate(invalid="ignore", divide="ignore"):
            e = np.where(tot >= args.min_bin, ok / tot, np.nan)

        a = ax[0, 0]
        # same colour scale as the DREAM-era maps (urw_p2_efficiency.py), so
        # the two can be read against each other without rescaling by eye
        im = a.pcolormesh(xe, ye, e.T, cmap="viridis", vmin=args.eff_vmin,
                          vmax=1)
        fig.colorbar(im, ax=a, label="efficiency")
        a.set(xlabel="x in the P2 pad frame [mm]", ylabel="y [mm]",
              title=f"{st} efficiency map ({bs:g} mm bins, "
                    f">={args.min_bin} tracks)")
        a.set_aspect("equal")

        a = ax[0, 1]
        s = a.scatter(pp["x"], pp["y"], c=e_pad, s=26, marker="s",
                      cmap="viridis", vmin=0, vmax=1)
        fig.colorbar(s, ax=a, label="efficiency")
        dead = ~live
        a.scatter(np.array(pp["x"])[dead], np.array(pp["y"])[dead], s=26,
                  marker="x", color="0.6", label="no tracks")
        a.set(xlabel="pad x [mm]", ylabel="pad y [mm]",
              title=f"{st} per-pad efficiency ({int(live.sum())} pads)")
        a.legend(fontsize=8)
        a.set_aspect("equal")

        a = ax[1, 0]
        d = np.linspace(0, 60, 121)
        dm = tracks[f"win_dmin_{st}"].to_numpy()[fid]
        a.hist(dm[np.isfinite(dm)], bins=d, histtype="step", lw=1.6,
               label="in the coincidence window")
        dm2 = tracks[f"any_dmin_{st}"].to_numpy()[fid]
        a.hist(dm2[np.isfinite(dm2)], bins=d, histtype="step", lw=1.2,
               label="anywhere in +-10 us")
        a.axvline(args.probe_r, color="r", ls="--",
                  label=f"probe r = {args.probe_r:g} mm")
        lb = r["loss_budget"]
        a.set(xlabel="distance track -> nearest fired pad [mm]",
              ylabel="tracks", yscale="log",
              title=f"matching distance\n"
                    f"{lb['nothing_within_probe_r_at_any_time']} misses with "
                    f"nothing there at any time, "
                    f"{lb['hit_at_the_place_but_out_of_time']} out of time")
        a.legend(fontsize=8)
        a.grid(alpha=0.3)

        a = ax[1, 1]
        pv = r["per_vmm"]
        xs = np.arange(len(pv))
        a.bar(xs, [p["eff"] or 0 for p in pv], color="steelblue")
        a.errorbar(xs, [p["eff"] or 0 for p in pv],
                   yerr=[[(p["eff"] or 0) - p["lo"] for p in pv],
                         [p["hi"] - (p["eff"] or 0) for p in pv]],
                   fmt="none", ecolor="k", lw=1)
        for i, p in enumerate(pv):
            a.text(i, 0.02, f"sdt {p['sdt']}\n{p['occupancy']:.0e} hits",
                   ha="center", fontsize=7, color="w")
        a.set_xticks(xs)
        a.set_xticklabels([f"VMM {p['vmm']}" for p in pv], fontsize=8)
        a.set(ylabel="efficiency", ylim=(0, 1),
              title="efficiency per readout chip\n(each has its own threshold "
                    "DAC)")
        a.grid(alpha=0.3, axis="y")

        a = ax[0, 3]
        ip = r["intra_pad"]
        e2 = np.array(ip["k"], float) / np.maximum(np.array(ip["n"], float), 1)
        e2[np.array(ip["n"]) < args.min_bin] = np.nan
        im = a.pcolormesh(ip["edges"], ip["edges"], e2.T, cmap="viridis",
                          vmin=0, vmax=1)
        fig.colorbar(im, ax=a, label="efficiency")
        a.set(xlabel="track - pad centre, x [mm]", ylabel="y [mm]",
              title="every pad folded onto one cell\n"
                    "(a threshold loses the edges first)")
        a.set_aspect("equal")

        a = ax[1, 3]
        good = pn >= 200
        a.hist(e_pad[good], bins=np.linspace(0, 1, 41), histtype="stepfilled",
               color="steelblue", alpha=0.8,
               label=f"{int(good.sum())} pads with >=200 tracks")
        a.axvline(r["efficiency"]["eff"], color="k", lw=1.5,
                  label=f"station average {r['efficiency']['eff']:.3f}")
        if args.dream_eff:
            a.axvline(args.dream_eff, color="crimson", ls=":", lw=2,
                      label=f"DREAM readout {args.dream_eff:.3f}")
        a.set(xlabel="per-pad efficiency", ylabel="pads",
              title="how the same detector performs pad by pad\n"
                    "(one threshold DAC per chip, no per-channel trim)")
        a.legend(fontsize=7)
        a.grid(alpha=0.3)

        a = ax[1, 2]
        adc = r.get("adc")
        if adc and "hist" in adc:
            he = np.array(adc["hist"]["edges"])
            a.step(0.5 * (he[1:] + he[:-1]), adc["hist"]["n"], where="mid",
                   lw=1.6, color="k", label="all chips")
            for v, h in sorted(adc.get("hist_per_vmm", {}).items()):
                a.step(0.5 * (he[1:] + he[:-1]), h, where="mid", lw=1.0,
                       alpha=0.7, label=f"VMM {v}")
            a.axvline(adc["min"], color="r", ls="--",
                      label=f"lowest ADC seen = {adc['min']}")
            a.set(xlabel="VMM pulse height (PDO) of the hit on the track",
                  ylabel="hits", xlim=(0, 600),
                  title=f"pulse height, cut off at the discriminator\n"
                        f"most probable {adc['mode']}, "
                        f"5th percentile {adc['percentiles']['5']:.0f}")
            a.legend(fontsize=7)
            a.grid(alpha=0.3)

        a = ax[0, 2]
        t = sorted(int(k) for k in adc["efficiency_vs_offline_threshold"])
        v = [adc["efficiency_vs_offline_threshold"][str(k)] for k in t]
        a.plot(t, v, "o-")
        if args.dream_eff:
            a.axhline(args.dream_eff, color="crimson", ls=":",
                      label=f"DREAM readout, same detector: {args.dream_eff:.3f}")
            a.legend(fontsize=8)
        a.set(xlabel="extra ADC threshold, on top of the hardware one",
              ylabel="efficiency", ylim=(0, 1),
              title="what a HIGHER threshold would have cost\n"
                    "(the slope at zero says how close the hardware one is)")
        a.grid(alpha=0.3)
        fig.suptitle(
            f"{st}   {res['run']}/{res['sub']}   uRWELL-referenced VMM "
            f"efficiency {r['efficiency']['eff']:.4f} on "
            f"{r['efficiency']['n']} fiducial tracks", fontsize=13)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fn = f"{args.figdir}/autopsy_{res['run']}_{st}.png"
        fig.savefig(fn, dpi=95)
        plt.close(fig)
        print(f"  wrote {fn}")


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tracks")
    ap.add_argument("--hits-dir", default=None)
    ap.add_argument("--json", default="autopsy_report.json")
    ap.add_argument("--figdir", default=None)
    ap.add_argument("--thresholds", default=None,
                    help="json {station: {vmm: sdt}} from the config file")
    ap.add_argument("--probe-r", type=float, default=15.0)
    ap.add_argument("--fid-r", type=float, default=9.0)
    ap.add_argument("--mask", choices=["ratio", "rate", "none"],
                    default="rate",
                    help="hot-channel rule: 'ratio' is what the production "
                         "sweep applied (8x the chip median, which masks the "
                         "beam on a chip with a raised threshold), 'rate' is "
                         "an absolute 20 kHz cut, 'none' keeps only the "
                         "documented bad channels")
    ap.add_argument("--min-pad", type=int, default=30)
    ap.add_argument("--min-bin", type=int, default=25)
    ap.add_argument("--bin", type=float, default=4.0)
    ap.add_argument("--eff-vmin", type=float, default=0.5,
                    help="colour-scale floor; 0.5 matches the DREAM maps")
    ap.add_argument("--max-hz", type=float, default=20000.0)
    ap.add_argument("--dream-eff", type=float, default=None,
                    help="the DREAM-readout efficiency of the same station")
    args = ap.parse_args()
    args.sdt = json.load(open(args.thresholds)) if args.thresholds else {}

    base = args.tracks.replace("_tracks.parquet", "")
    tracks = pd.read_parquet(args.tracks)
    meta = json.load(open(base + ".json"))
    print(f"{len(tracks)} tracks, {meta['capture']['n_captures']} captures")

    res = {"run": meta["run"], "sub": meta["sub"],
           "n_tracks_timed": meta["n_tracks_timed"],
           "capture_coverage_frac": meta["capture_coverage_frac"],
           "track_quality": meta["tracks"], "stations": {}}
    for st in STATIONS:
        if f"px_{st}" not in tracks:
            continue
        fn = f"{base}_hits_{st}.parquet"
        hits = pd.read_parquet(fn) if os.path.exists(fn) else None
        r = analyse(tracks, hits, meta, st, args)
        r["latency"] = meta["stations"][st]["latency"]
        r["window_ns"] = meta["stations"][st]["window_ns"]
        r["frame"] = {k: meta["stations"][st]["frame"][k]
                      for k in ("affine_rotation_deg", "singular_values",
                                "det", "rigid_rmse_mm")}
        res["stations"][st] = r
        e = r["efficiency"]
        print(f"[{st}] efficiency {e['eff']:.4f} on {e['n']} tracks "
              f"(no mask {r['efficiency_no_mask']['eff']:.4f})")
        print(f"    losses: {r['loss_budget']}")
        for p in r["per_vmm"]:
            print(f"    VMM {p['vmm']:2d} sdt {str(p['sdt']):>4s} "
                  f"eff {0 if p['eff'] is None else p['eff']:.3f} "
                  f"on {p['n']:7d}  occ {p['occupancy']:.2e} "
                  f"(loudest ch {p['loudest_channel']}: "
                  f"{p['loudest_hits']:.1e})")
        del hits
    cross_station(tracks, args, res)
    print(f"cross-station: {json.dumps(res['cross_station'], indent=1)}")

    if args.figdir:
        os.makedirs(args.figdir, exist_ok=True)
        figures(tracks, meta, res, args)
    with open(args.json, "w") as f:
        json.dump(res, f, indent=1, default=float)
    print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
