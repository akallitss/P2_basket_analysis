#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
urw_vmm_efficiency.py -- VMM-read P2 efficiency referenced to uRWELL tracks.

The VMM campaign (run_21..run_71) reads the three P2 stations on the VMM/SRS
DAQ while the DREAM DAQ keeps only the two EIC uRWELL references on FEU 1. So
the stations have, until now, only had a *trigger-referenced* efficiency
(vmm_efficiency.py: is there a station hit in the trigger coincidence window),
which counts a hit anywhere in the detector and knows nothing about where the
particle went.

With the two streams matched event by event (match_streams.py) the real
measurement becomes possible, and it is the same one the DREAM-read campaign
had:

    efficiency = N(uRWELL track in the pads AND a VMM cluster within probe_r)
                 ----------------------------------------------------------
                            N(uRWELL track in the pads)

Chain, per sub_run:

  1. uRWELL tracks from the DREAM FEU-1 hits, front-back agreement cut
     (urw_p2_efficiency.build_tracks -- unchanged, this is the same reference
     as the DREAM-side measurement, which is what makes the two comparable).
  2. eventId -> t_vmm from the stored match, so every track carries the time it
     happened in the VMM clock.
  3. VMM hits of the station, kept only where they fall near one of those
     times; pads attached through vmm_stations.build_pad_table.
  4. The station's coincidence latency is fitted from the data (dt peak), and
     the signal window is +-n_sigma around it. An equal window offset well past
     the peak measures the accidental rate that survives it.
  5. uRWELL frame -> pad frame by free affine (fit_frame), then residuals and
     efficiency with Clopper-Pearson intervals, plus a per-pad map.

Run on lxplus:

    export SACLAY_MM_DIR=$PWD URW_DET_DIR=<eos>/config/detectors/
    python3 urw_vmm_efficiency.py run_32 meshscan_m00V --out out_eff
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extract_vmm_triggers as evt        # noqa: E402
import urw_p2_efficiency as UPE           # noqa: E402
import vmm_stations as VS                 # noqa: E402

BASE = evt.BASE
STATIONS = VS.STATIONS
# A single pad in the core of this beam takes a few kHz. Anything two decades
# above that is not the detector responding -- see the note in station_hits.
MAX_CHANNEL_HZ = 20000.0


# --------------------------------------------------------------------------- #
# inputs
# --------------------------------------------------------------------------- #
def load_match(run, sub, match_dir):
    """eventId -> t_vmm_ns for the matched events of this sub_run."""
    fn = f"{match_dir}/match_{run}_{sub}.npz"
    if not os.path.exists(fn):
        sys.exit(f"no match file {fn} -- run the matching sweep first")
    d = np.load(fn)
    m = d["matched"]
    if not m.any():
        sys.exit(f"{run}/{sub}: nothing matched, no track can be timed")
    return pd.DataFrame({"eventId": d["event_id"][m].astype(np.int64),
                         "t_vmm": d["t_vmm_ns"][m]})


def station_hits(sub_dir, station, t_want, pre_ns, post_ns, source="auto",
                 max_captures=None, quiet=True, max_hz=MAX_CHANNEL_HZ):
    """VMM hits of one station that land near one of the times in `t_want`.

    Streams the captures and throws away, per chunk, everything that is not
    within [-pre, +post] of a wanted trigger time: a sub_run has ~4e7 hits and
    only a few percent of them belong to a matched event, so this is what keeps
    the whole thing in memory. Returns (idx, dt, vmm, ch, adc) where idx is the
    index into t_want.
    """
    vmms = np.array(VS.STATION_VMMS[station])
    idx, dt, vv, cc, aa = [], [], [], [], []
    occ = np.zeros((32, 64), np.int64)
    span = {}                       # capture -> [t_first, t_last] of its hits
    for _name, cols in evt.iter_columns(sub_dir, max_captures, source, quiet):
        sel = np.isin(cols["vmm"], vmms)
        if not sel.any():
            continue
        np.add.at(occ, (np.clip(cols["vmm"][sel], 0, 31),
                        np.clip(cols["ch"][sel], 0, 63)), 1)
        t = evt.hit_times(cols, sel)
        if len(t):
            # quantiles, not min/max: a handful of hits per capture carry
            # srs_timestamp = 0, and taking the minimum stretched every
            # capture's span back to t = 0. Unioned, that made
            # covered_by_captures accept everything -- the cut silently did
            # nothing -- and made the live time come out 25x too long.
            q0, q1 = np.quantile(t, [0.002, 0.998])
            lo_hi = span.setdefault(_name, [float(q0), float(q1)])
            lo_hi[0] = min(lo_hi[0], float(q0))
            lo_hi[1] = max(lo_hi[1], float(q1))
        # nearest wanted time, then the window cut
        k = np.searchsorted(t_want, t)
        k = np.clip(k, 1, len(t_want) - 1)
        left = t - t_want[k - 1]
        right = t - t_want[k]
        use_left = np.abs(left) < np.abs(right)
        kk = np.where(use_left, k - 1, k)
        d = np.where(use_left, left, right)
        keep = (d >= -pre_ns) & (d <= post_ns)
        if not keep.any():
            continue
        idx.append(kk[keep])
        dt.append(d[keep])
        vv.append(cols["vmm"][sel][keep])
        cc.append(cols["ch"][sel][keep])
        aa.append(cols["adc"][sel][keep] if "adc" in cols
                  else np.ones(keep.sum(), np.int16))
    cov = sorted(span.values())
    if not idx:
        return (pd.DataFrame(columns=["idx", "dt", "vmm", "ch", "adc"]),
                {"coverage": cov})
    h = pd.DataFrame({
        "idx": np.concatenate(idx), "dt": np.concatenate(dt),
        "vmm": np.concatenate(vv).astype(np.int16),
        "ch": np.concatenate(cc).astype(np.int16),
        "adc": np.concatenate(aa).astype(np.int32)})
    # Hot channels have to be flagged on an ABSOLUTE rate, not on a ratio to
    # their own chip's median.
    #
    # The ratio test (the original: 8x the chip's median occupancy) flags the
    # beam spot on any chip whose threshold has been raised. VMM 9 of P2_MID
    # sits at sdt 256, which suppresses its noise; its median channel then
    # collapses, its illuminated pads stand 30-50x above it, and they get
    # masked. That threw away 64,630 hits arriving exactly at the +115 ns
    # trigger latency -- 7.7 points of P2_MID efficiency (0.5916 -> 0.6688),
    # deleted as noise. Flagging on the out-of-time occupancy does not fix it
    # either: the beam delivers particles at tens of kHz while the matched
    # trigger runs at ~1.4 kHz, so most of a beam pad's hits are out of time
    # with a matched trigger and it stays flagged.
    #
    # An absolute rate does separate them cleanly. A 12 x 12 mm pad in the
    # core of this beam takes a few kHz at most; the pathological channels are
    # at 10^5 Hz (VMM 4 ch 39 at 460 kHz, VMM 19 ch 20 at 116 kHz in run_46).
    # There are two decades of clear air between the two populations.
    #
    # And masking matters far less here than it did for the centroid analyses
    # this inherited the machinery from: the efficiency asks whether ANY pad
    # within probe_r fired, so a hot channel costs accidentals, and those are
    # measured in the sideband -- 5e-4 on P2_OUT with 131 kHz of noise left in.
    # Where a wrong mask costs 7.7 points and a missing one costs 0.05, the
    # cut belongs on the safe side.
    # wall clock of the sub_run, from the triggers themselves -- summing the
    # capture spans double-counts wherever two captures overlap
    dur = float(t_want[-1] - t_want[0]) / 1e9 if len(t_want) > 1 else 0.0
    mask = VS.merge_masks(VS.NOISY_CHANNELS,
                          _hot_from_rate(occ, vmms, dur, max_hz))
    h, n_drop = VS.apply_channel_mask(h, mask)
    return h.reset_index(drop=True), {"mask": {int(k): v for k, v in
                                               mask.items()},
                                      "n_masked_hits": int(n_drop),
                                      "mask_by_ratio_would_have_been": {
                                          int(k): v for k, v in
                                          _hot_from_occupancy(occ,
                                                              vmms).items()},
                                      "live_seconds": dur,
                                      "coverage": cov,
                                      "n_captures": len(cov)}


def _hot_from_rate(occ, vmms, seconds, max_hz=MAX_CHANNEL_HZ):
    """{vmm: [ch]} for channels whose rate is beyond anything the beam can do."""
    if seconds <= 0:
        return {}
    found = {}
    for v in vmms:
        hot = np.flatnonzero(occ[v] / seconds > max_hz)
        if hot.size:
            found[int(v)] = [int(x) for x in hot]
    return found


def _hot_from_occupancy(occ, vmms, ratio=VS.HOT_RATIO, min_hits=200):
    """{vmm: [ch]} for channels far above their own VMM's median occupancy.

    Kept only to record, in the summary, what the old ratio rule would have
    masked -- see the note in station_hits for why it is not used.
    """
    found = {}
    for v in vmms:
        c = occ[v].astype(float)
        if c.sum() < min_hits:
            continue
        nz = c[c > 0]
        med = float(np.median(nz)) if nz.size else 0.0
        if med <= 0:
            continue
        hot = np.flatnonzero(c > ratio * med)
        if hot.size:
            found[int(v)] = [int(x) for x in hot]
    return found


# --------------------------------------------------------------------------- #
# timing
# --------------------------------------------------------------------------- #
def fit_latency(dt, bin_ns=10.0, lo=-1000.0, hi=3000.0):
    """Peak position and width of the trigger-to-hit delay.

    A plain histogram peak plus a robust width: the shape is a drift-time
    distribution, not a Gaussian, so a Gaussian fit would understate the tail
    that the signal window has to contain.
    """
    h, ed = np.histogram(dt, bins=int((hi - lo) / bin_ns), range=(lo, hi))
    c = 0.5 * (ed[1:] + ed[:-1])
    if not h.sum():
        return dict(mu=0.0, sigma=50.0, peak=0, background=0.0, contrast=0.0)
    k = int(np.argmax(h))
    mu = float(c[k])
    # background from the far half of the range, per bin
    far = np.abs(c - mu) > 1500.0
    bg = float(h[far].mean()) if far.any() else 0.0
    # width: where the peak falls to a tenth above background
    thr = bg + 0.1 * (h[k] - bg)
    above = np.where(h > thr)[0]
    near = above[(above > k - 60) & (above < k + 60)] if len(above) else [k]
    w = float((c[near[-1]] - c[near[0]]) / 2) if len(near) > 1 else 25.0
    return dict(mu=mu, sigma=max(w, 12.5), peak=int(h[k]), background=bg,
                contrast=float(h[k] / bg) if bg > 0 else float("inf"))


# --------------------------------------------------------------------------- #
# pad hits in a time window
# --------------------------------------------------------------------------- #
def pads_in_window(hits, pads, lo, hi):
    """(idx, x, y, adc) of every mapped pad hit inside [lo, hi]."""
    h = hits[(hits["dt"] >= lo) & (hits["dt"] <= hi)]
    h = h.merge(pads[["vmm", "ch", "pad_cx", "pad_cy"]], on=["vmm", "ch"],
                how="left")
    h = h[h["pad_cx"].notna()]
    return (h["idx"].to_numpy(np.int64), h["pad_cx"].to_numpy(float),
            h["pad_cy"].to_numpy(float), h["adc"].to_numpy(float))


def nearest_pad(idx, x, y, px, py, n):
    """Per track: distance to the closest hit pad, and that pad's position.

    Vectorised: a groupby over ~3e5 events in python is minutes, this is
    milliseconds.
    """
    d = np.full(n, np.inf)
    bx = np.full(n, np.nan)
    by = np.full(n, np.nan)
    if not len(idx):
        return d, bx, by
    dd = np.hypot(x - px[idx], y - py[idx])
    np.minimum.at(d, idx, dd)
    best = dd == d[idx]                      # the winning hit of each event
    bx[idx[best]] = x[best]
    by[idx[best]] = y[best]
    return d, bx, by


def single_pad_events(idx, n):
    """Mask over hits selecting events whose window holds exactly one pad."""
    if not len(idx):
        return np.zeros(0, bool)
    cnt = np.bincount(idx, minlength=n)
    return cnt[idx] == 1


# --------------------------------------------------------------------------- #
# one station
# --------------------------------------------------------------------------- #
def covered_by_captures(t, coverage, margin_ns=0.0):
    """Which times fall inside a capture the decoder actually produced.

    Captures go missing -- run_32/meshscan_m00V has 20 pcapng files and 17
    entries in hits_store -- and a track that happened while the VMM wrote
    nothing is not an inefficiency, it is an absence. Without this cut the
    measurement is biased low by exactly the missing fraction.
    """
    ok = np.zeros(len(t), bool)
    for lo, hi in coverage:
        ok |= (t >= lo - margin_ns) & (t <= hi + margin_ns)
    return ok


def run_station(station, tracks, hits, pads, args, coverage=()):
    z = STATIONS[station]["z_mm"]
    pad_xy = pads.dropna(subset=["pad_cx"])[["pad_cx", "pad_cy"]].to_numpy()
    n = len(tracks)
    t_vmm = tracks["t_vmm"].to_numpy()
    in_cap = (covered_by_captures(t_vmm, coverage) if len(coverage)
              else np.ones(n, bool))
    out = {"station": station, "z_mm": z, "n_tracks": int(n),
           "n_tracks_in_captures": int(in_cap.sum()),
           "capture_coverage_frac": float(in_cap.mean()),
           "n_hits_near_trigger": int(len(hits))}
    if not len(hits):
        out["error"] = "no VMM hits near any matched trigger"
        return out, None

    lat = fit_latency(hits["dt"].to_numpy())
    out["latency"] = lat
    lo = lat["mu"] - args.n_sigma * lat["sigma"]
    hi = lat["mu"] + args.n_sigma * lat["sigma"]
    out["window_ns"] = [lo, hi]

    idx, hx, hy, hadc = pads_in_window(hits, pads, lo, hi)
    # the same window pushed past the peak: whatever it finds is accidental
    off = args.sideband_ns
    sidx, sx, sy, _ = pads_in_window(hits, pads, lo + off, hi + off)
    out["n_pad_hits_in_window"] = int(len(idx))
    out["n_events_with_pad"] = int(len(np.unique(idx))) if len(idx) else 0
    # ceiling: events with ANY mapped pad hit anywhere in the read window.
    # If the efficiency sits well below this, the coincidence window is the
    # thing costing it, not the detector.
    aidx, ax, ay, _ = pads_in_window(hits, pads, -1e9, 1e9)
    out["n_events_any_hit"] = int(len(np.unique(aidx))) if len(aidx) else 0

    # --- frame fit on unambiguous events -----------------------------------
    # One pad in the window means no competing cluster and no noise pad to
    # pick the wrong one, so this is the clean sample for geometry. The
    # efficiency below uses every event, ambiguous ones included.
    proj = UPE.project(tracks, z, args.dz)
    solo = single_pad_events(idx, n)
    if solo.sum() < args.min_events:
        out["error"] = (f"only {int(solo.sum())} single-pad events, need "
                        f"{args.min_events} to fit the frame")
        return out, None
    frame, aff, _ = UPE.fit_frame(proj[idx[solo]],
                                  np.column_stack([hx[solo], hy[solo]]),
                                  args.frame_win)
    out["frame"] = frame
    out["n_frame_fit"] = int(solo.sum())
    px, py = aff.apply(proj[:, 0], proj[:, 1])

    # --- residuals on that same clean sample -------------------------------
    dxs = hx[solo] - px[idx[solo]]
    dys = hy[solo] - py[idx[solo]]
    keep = (np.abs(dxs) < 40) & (np.abs(dys) < 40)
    out["residual"] = {
        "dx_rms_mm": float(np.std(dxs[keep])) if keep.any() else None,
        "dy_rms_mm": float(np.std(dys[keep])) if keep.any() else None,
        "dx_med_mm": float(np.median(dxs[keep])) if keep.any() else None,
        "dy_med_mm": float(np.median(dys[keep])) if keep.any() else None,
        "n": int(keep.sum())}

    # --- fiducial: the prediction has to land on the instrumented plane ----
    d2 = ((px[:, None] - pad_xy[None, :, 0]) ** 2
          + (py[:, None] - pad_xy[None, :, 1]) ** 2)
    dmin_pad = np.sqrt(d2.min(axis=1))
    near_pad = np.argmin(d2, axis=1)
    fid = (dmin_pad < args.fid_r) & in_cap
    out["n_fiducial"] = int(fid.sum())
    if not fid.any():
        out["error"] = "no track lands on the pad plane"
        return out, None

    # --- efficiency: is ANY pad hit close to where the track went? ---------
    # Not "is the leading cluster close": with a live neighbour firing in the
    # same window, taking the loudest pad would count a real hit as a miss.
    dmin, bx, by = nearest_pad(idx, hx, hy, px, py, n)
    sdmin, _, _ = nearest_pad(sidx, sx, sy, px, py, n)
    hit = fid & (dmin < args.probe_r)
    k, nn = int(hit.sum()), int(fid.sum())
    e = k / nn
    elo, ehi = UPE.clopper_pearson(k, nn)
    a_rate = float((fid & (sdmin < args.probe_r)).sum() / nn)
    out["efficiency"] = {
        "value": e, "lo": elo, "hi": ehi, "k": k, "n": nn,
        "accidental_rate": a_rate,
        "corrected": float((e - a_rate) / (1 - a_rate)) if a_rate < 1 else None,
        "vs_probe_r": {str(r): float((fid & (dmin < r)).sum() / nn)
                       for r in (5, 8, 10, 12, 15, 20, 25, 30, 40)}}
    # how much of the loss is the time window? Widen it and see.
    vs_win = {}
    for f in (1, 2, 3, 6, 10, 20):
        w0 = lat["mu"] - f * lat["sigma"]
        w1 = lat["mu"] + f * lat["sigma"]
        wi, wx, wy, _ = pads_in_window(hits, pads, w0, w1)
        wd, _, _ = nearest_pad(wi, wx, wy, px, py, n)
        vs_win[f"{f}sig"] = {
            "window_ns": [w0, w1],
            "eff": float((fid & (wd < args.probe_r)).sum() / nn)}
    out["efficiency"]["vs_window"] = vs_win

    # --- per-pad efficiency, for the map -----------------------------------
    npad = len(pad_xy)
    tot = np.bincount(near_pad[fid], minlength=npad)
    got = np.bincount(near_pad[fid & hit], minlength=npad)
    out["per_pad"] = {
        "x": [round(float(v), 2) for v in pad_xy[:, 0]],
        "y": [round(float(v), 2) for v in pad_xy[:, 1]],
        "n": tot.tolist(), "k": got.tolist()}

    det = pd.DataFrame({"px": px, "py": py, "dmin": dmin, "bx": bx, "by": by,
                        "fid": fid, "hit": hit,
                        "t_ns": tracks["t_ns"].to_numpy()})
    return out, det


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("sub")
    ap.add_argument("--base", default=BASE)
    ap.add_argument("--match-dir", default=f"{BASE}/vmm/matching")
    ap.add_argument("--out", default="out_vmm_eff")
    ap.add_argument("--stations", default="all")
    ap.add_argument("--max-chunks", type=int, default=0,
                    help="DREAM hit files to read (0 = all)")
    ap.add_argument("--max-captures", type=int, default=None,
                    help="VMM captures to read (default all)")
    ap.add_argument("--track-cut", type=float, default=3.0)
    ap.add_argument("--cluster-r", type=float, default=15.0)
    ap.add_argument("--probe-r", type=float, default=15.0)
    ap.add_argument("--fid-r", type=float, default=9.0)
    ap.add_argument("--frame-win", type=float, default=25.0)
    ap.add_argument("--n-sigma", type=float, default=3.0)
    ap.add_argument("--sideband-ns", type=float, default=4000.0)
    ap.add_argument("--pre-ns", type=float, default=1500.0)
    ap.add_argument("--post-ns", type=float, default=6000.0)
    ap.add_argument("--min-events", type=int, default=500)
    ap.add_argument("--save-detail", action="store_true")
    args = ap.parse_args()

    sub_dir_dream = f"{args.base}/runs/{args.run}/{args.sub}"
    sub_dir_vmm = evt.sub_dir_of(args.run, args.sub, args.base)
    run_json = f"{args.base}/runs/{args.run}/run_config.json"
    os.makedirs(args.out, exist_ok=True)

    print(f"[tracks] {args.run}/{args.sub}")
    tracks, tinfo, dz = UPE.build_tracks(run_json, sub_dir_dream,
                                         args.max_chunks or None,
                                         args.track_cut)
    args.dz = dz
    print(f"  {tinfo['n_good_track']} good tracks of {tinfo['n_4coord']} "
          f"4-coordinate events ({100*tinfo['good_track_frac']:.1f}%)")

    print("[match] attaching VMM times")
    mt = load_match(args.run, args.sub, args.match_dir)
    tracks = tracks.merge(mt, on="eventId", how="inner").sort_values("t_vmm")
    tracks = tracks.reset_index(drop=True)
    print(f"  {len(tracks)} tracks carry a matched VMM trigger time")
    if not len(tracks):
        sys.exit("no track is in the matched event set")

    t_want = tracks["t_vmm"].to_numpy()
    pads = VS.build_pad_table()
    want = (list(STATIONS) if args.stations == "all"
            else args.stations.split(","))

    summary = {"run": args.run, "sub": args.sub, "tracks": tinfo,
               "n_tracks_timed": int(len(tracks)),
               "args": {k: v for k, v in vars(args).items()},
               "stations": []}
    for st in want:
        print(f"[{st}] reading VMM hits near the matched triggers")
        h, minfo = station_hits(sub_dir_vmm, st, t_want, args.pre_ns,
                                args.post_ns,
                                max_captures=args.max_captures)
        print(f"  {len(h)} hits in [-{args.pre_ns:.0f}, +{args.post_ns:.0f}] ns"
              f" after masking {minfo.get('n_masked_hits', 0)} on "
              f"{len(minfo.get('mask', {}))} hot channel(s)")
        res, det = run_station(st, tracks, h,
                               pads[pads["station"] == st], args,
                               coverage=minfo.get("coverage", ()))
        res["masking"] = minfo
        summary["stations"].append(res)
        if "efficiency" in res:
            ef = res["efficiency"]
            print(f"  latency {res['latency']['mu']:+.0f} ns "
                  f"(contrast {res['latency']['contrast']:.1f}), "
                  f"window {res['window_ns'][0]:.0f}..{res['window_ns'][1]:.0f} ns")
            print(f"  residual rms {res['residual']['dx_rms_mm']:.2f} / "
                  f"{res['residual']['dy_rms_mm']:.2f} mm on "
                  f"{res['residual']['n']} single-pad events")
            print(f"  VMM captures cover "
                  f"{100*res['capture_coverage_frac']:.1f}% of the tracks")
            print(f"  EFFICIENCY {ef['value']:.4f} "
                  f"(+{ef['hi']-ef['value']:.4f}/-{ef['value']-ef['lo']:.4f}) "
                  f"on {ef['n']} fiducial tracks; "
                  f"accidental {ef['accidental_rate']:.4f}")
        else:
            print(f"  no efficiency: {res.get('error')}")
        if det is not None and args.save_detail:
            det.to_parquet(f"{args.out}/detail_{args.run}_{args.sub}_{st}.parquet")

    fn = f"{args.out}/vmm_eff_{args.run}_{args.sub}.json"
    with open(fn, "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print(f"wrote {fn}")


if __name__ == "__main__":
    main()
