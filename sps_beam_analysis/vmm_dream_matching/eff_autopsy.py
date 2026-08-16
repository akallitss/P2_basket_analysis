#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eff_autopsy.py -- why is the uRWELL-referenced VMM efficiency not 95%?

urw_vmm_efficiency.py gives 85% on P2_OUT at the best VMM configuration where
the DREAM readout of the same detectors gave 96-97%. That gap is either real
(and then it is the headline of the VMM campaign) or it is an artefact of the
measurement, and the only way to tell is to take one sub_run apart.

This script makes ONE streaming pass over the captures of one sub_run and
writes out everything the question needs, so that every hypothesis can then be
tested offline against the same event sample rather than against a re-run:

  autopsy_<run>_<sub>_tracks.parquet   one row per matched uRWELL track:
      the projected impact point in each station's pad frame, the nearest
      instrumented pad, the nearest FIRED pad (with and without the hot-channel
      mask), that hit's ADC and time, how many station hits the event had at
      all, the local trigger rate and the position in the spill.
  autopsy_<run>_<sub>_hits_<station>.parquet   every station hit near a matched
      trigger, un-masked and un-cut: track index, dt, vmm, ch, ADC,
      over-threshold flag, and the distance from the pad to the prediction.
  autopsy_<run>_<sub>.json   occupancy of all 32x64 channels, the masks, the
      per-capture coverage, latency and frame fits, and a wide (+-50 us) dt
      histogram per VMM.

Nothing is thrown away on the basis of a hypothesis being tested: the masks
are recorded as flags, not applied, and the time window is widened to +-8 us
so the sidebands and the second-trigger population are both in the file.

Memory: the pass streams chunk by chunk and keeps only hits within the widened
window of a matched trigger (~7e6 rows x 18 B ~ 130 MB) plus the track table;
it never holds a whole capture's hits.

    export SACLAY_MM_DIR=$PWD URW_DET_DIR=<eos>/config/detectors/
    python3 eff_autopsy.py run_46 cfg_gain4.5_peaktime200 --out autopsy
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extract_vmm_triggers as evt        # noqa: E402
import urw_p2_efficiency as UPE           # noqa: E402
import vmm_stations as VS                 # noqa: E402
import urw_vmm_efficiency as UVE          # noqa: E402

BASE = evt.BASE
STATIONS = VS.STATIONS
WIDE_NS = 50000.0          # half-range of the coarse dt histogram
WIDE_BIN = 250.0


# --------------------------------------------------------------------------- #
def collect(sub_dir, t_want, pre_ns, post_ns, source="auto", max_captures=None):
    """Every detector-VMM hit near a matched trigger, plus bookkeeping.

    One pass for all three stations: the production script pays the pcapng
    decode three times over, which is the only reason it is slow.
    """
    vmms = np.array(VS.DETECTOR_VMMS)
    occ = np.zeros((32, 64), np.int64)          # whole-capture occupancy
    nbin = int(2 * WIDE_NS / WIDE_BIN)
    wide = np.zeros((32, nbin), np.int64)       # dt histogram per VMM
    span = {}
    idx, dt, vv, cc, aa, ot = [], [], [], [], [], []
    n_hits_total = 0
    n_cap = 0
    t0 = time.time()
    for name, cols in evt.iter_columns(sub_dir, max_captures, source, False):
        n_hits_total += len(cols["vmm"])
        sel = np.isin(cols["vmm"], vmms)
        # Coverage is "the DAQ was recording", so it is taken from every hit of
        # the chunk, not from the station's: a capture in which this station
        # said nothing is exactly the case that must not look like coverage.
        # The SRS timestamp alone is enough (it is the coarse clock).
        if len(cols["srs_timestamp"]):
            t_lo = float(cols["srs_timestamp"].min()) * evt.SRS_TICK_NS
            t_hi = float(cols["srs_timestamp"].max()) * evt.SRS_TICK_NS
            lo_hi = span.setdefault(name, [t_lo, t_hi])
            lo_hi[0] = min(lo_hi[0], t_lo)
            lo_hi[1] = max(lo_hi[1], t_hi)
        if not sel.any():
            continue
        np.add.at(occ, (np.clip(cols["vmm"][sel], 0, 31),
                        np.clip(cols["ch"][sel], 0, 63)), 1)
        t = evt.hit_times(cols, sel)
        k = np.clip(np.searchsorted(t_want, t), 1, len(t_want) - 1)
        left = t - t_want[k - 1]
        right = t - t_want[k]
        use_left = np.abs(left) < np.abs(right)
        kk = np.where(use_left, k - 1, k)
        d = np.where(use_left, left, right)
        # coarse, wide dt histogram: the shape outside the signal window says
        # whether a station's hits are trigger-correlated at all
        w = np.floor((d + WIDE_NS) / WIDE_BIN).astype(np.int64)
        good = (w >= 0) & (w < nbin)
        np.add.at(wide, (np.clip(cols["vmm"][sel][good], 0, 31), w[good]), 1)
        keep = (d >= -pre_ns) & (d <= post_ns)
        if not keep.any():
            continue
        idx.append(kk[keep].astype(np.int32))
        dt.append(d[keep].astype(np.float32))
        vv.append(cols["vmm"][sel][keep].astype(np.int8))
        cc.append(cols["ch"][sel][keep].astype(np.int8))
        aa.append((cols["adc"][sel][keep] if "adc" in cols
                   else np.zeros(int(keep.sum()), np.int16)).astype(np.int16))
        ot.append((cols["over_threshold"][sel][keep] if "over_threshold" in cols
                   else np.ones(int(keep.sum()), bool)))
    n_cap = len(span)
    print(f"  {n_cap} captures, {n_hits_total} hits read in "
          f"{time.time()-t0:.0f} s")
    h = pd.DataFrame({
        "idx": np.concatenate(idx) if idx else np.zeros(0, np.int32),
        "dt": np.concatenate(dt) if dt else np.zeros(0, np.float32),
        "vmm": np.concatenate(vv) if vv else np.zeros(0, np.int8),
        "ch": np.concatenate(cc) if cc else np.zeros(0, np.int8),
        "adc": np.concatenate(aa) if aa else np.zeros(0, np.int16),
        "over_thr": np.concatenate(ot) if ot else np.zeros(0, bool)})
    cov = sorted(v for v in span.values() if v is not None)
    return h, occ, wide, cov, {"n_hits_total": int(n_hits_total),
                               "n_captures": int(n_cap)}


# --------------------------------------------------------------------------- #
def spill_structure(t, gap_ns=1.0e9):
    """Spill index and time-into-spill for a sorted array of trigger times."""
    if not len(t):
        return np.zeros(0, np.int32), np.zeros(0, np.float32)
    brk = np.concatenate([[0], np.flatnonzero(np.diff(t) > gap_ns) + 1,
                          [len(t)]])
    sid = np.zeros(len(t), np.int32)
    tin = np.zeros(len(t), np.float32)
    for i in range(len(brk) - 1):
        s, e = brk[i], brk[i + 1]
        sid[s:e] = i
        tin[s:e] = (t[s:e] - t[s]) / 1e6        # ms into the spill
    return sid, tin


def local_rate(t, win_ns=1.0e6):
    """Triggers within +-win of each trigger -- the instantaneous rate proxy."""
    lo = np.searchsorted(t, t - win_ns)
    hi = np.searchsorted(t, t + win_ns)
    return (hi - lo - 1).astype(np.int32)


# --------------------------------------------------------------------------- #
def station_frame(hits, pads, tracks, z, dz, args):
    """Latency, window and the uRWELL->pad affine, exactly as production."""
    lat = UVE.fit_latency(hits["dt"].to_numpy().astype(float))
    lo = lat["mu"] - args.n_sigma * lat["sigma"]
    hi = lat["mu"] + args.n_sigma * lat["sigma"]
    idx, hx, hy, _ = UVE.pads_in_window(hits, pads, lo, hi)
    proj = UPE.project(tracks, z, dz)
    solo = UVE.single_pad_events(idx, len(tracks))
    if solo.sum() < args.min_events:
        return lat, (lo, hi), None, None, proj
    frame, aff, _ = UPE.fit_frame(proj[idx[solo]],
                                  np.column_stack([hx[solo], hy[solo]]),
                                  args.frame_win)
    return lat, (lo, hi), frame, aff, proj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("sub")
    ap.add_argument("--base", default=BASE)
    ap.add_argument("--match-dir", default=f"{BASE}/vmm/matching")
    ap.add_argument("--out", default="autopsy")
    ap.add_argument("--max-chunks", type=int, default=0)
    ap.add_argument("--max-captures", type=int, default=None)
    ap.add_argument("--track-cut", type=float, default=3.0)
    ap.add_argument("--probe-r", type=float, default=15.0)
    ap.add_argument("--fid-r", type=float, default=9.0)
    ap.add_argument("--frame-win", type=float, default=25.0)
    ap.add_argument("--n-sigma", type=float, default=3.0)
    ap.add_argument("--pre-ns", type=float, default=2000.0)
    ap.add_argument("--post-ns", type=float, default=8000.0)
    ap.add_argument("--min-events", type=int, default=500)
    args = ap.parse_args()

    sub_dir_dream = f"{args.base}/runs/{args.run}/{args.sub}"
    sub_dir_vmm = evt.sub_dir_of(args.run, args.sub, args.base)
    run_json = f"{args.base}/runs/{args.run}/run_config.json"
    os.makedirs(args.out, exist_ok=True)
    tag = f"{args.run}_{args.sub}"

    print(f"[tracks] {args.run}/{args.sub}")
    tracks, tinfo, dz = UPE.build_tracks(run_json, sub_dir_dream,
                                         args.max_chunks or None,
                                         args.track_cut)
    print(f"  {tinfo['n_good_track']} good tracks "
          f"({100*tinfo['good_track_frac']:.1f}% of 4-coordinate events)")
    mt = UVE.load_match(args.run, args.sub, args.match_dir)
    tracks = tracks.merge(mt, on="eventId", how="inner").sort_values("t_vmm")
    tracks = tracks.reset_index(drop=True)
    n = len(tracks)
    print(f"  {n} tracks carry a matched VMM trigger time")

    t_want = tracks["t_vmm"].to_numpy()
    sid, tin = spill_structure(t_want)
    tracks["spill"] = sid
    tracks["t_in_spill_ms"] = tin
    tracks["n_trig_1ms"] = local_rate(t_want)
    tracks["dt_prev_trig_us"] = np.concatenate(
        [[np.inf], np.diff(t_want) / 1e3]).astype(np.float32)
    tracks["dt_next_trig_us"] = np.concatenate(
        [np.diff(t_want) / 1e3, [np.inf]]).astype(np.float32)

    print("[vmm] one pass over the captures for all three stations")
    hits, occ, wide, cov, cinfo = collect(sub_dir_vmm, t_want, args.pre_ns,
                                          args.post_ns,
                                          max_captures=args.max_captures)
    print(f"  {len(hits)} detector hits within "
          f"[-{args.pre_ns:.0f}, +{args.post_ns:.0f}] ns of a matched trigger")

    pads_all = VS.build_pad_table()
    in_cap = (UVE.covered_by_captures(t_want, cov) if cov
              else np.ones(n, bool))
    tracks["in_capture"] = in_cap

    summary = {"run": args.run, "sub": args.sub, "tracks": tinfo,
               "n_tracks_timed": int(n), "dz_mm": float(dz),
               "args": {k: v for k, v in vars(args).items()},
               "capture": cinfo, "coverage": cov,
               "capture_coverage_frac": float(in_cap.mean()),
               "occupancy": occ.tolist(),
               "wide_dt": {"bin_ns": WIDE_BIN, "half_range_ns": WIDE_NS,
                           "hist": wide.tolist()},
               "stations": {}}

    for st in STATIONS:
        z = STATIONS[st]["z_mm"]
        vmms = VS.STATION_VMMS[st]
        pads = pads_all[pads_all["station"] == st]
        h = hits[hits["vmm"].isin(vmms)].reset_index(drop=True)
        print(f"[{st}] {len(h)} hits near a trigger on VMMs {vmms}")

        # the two masks, as flags rather than as a cut
        auto = UVE._hot_from_occupancy(occ, np.array(vmms))
        prod = VS.merge_masks(VS.NOISY_CHANNELS, auto)
        def flag(mask):
            f = np.zeros(len(h), bool)
            for v, chans in mask.items():
                f |= (h["vmm"].to_numpy() == v) & np.isin(h["ch"].to_numpy(),
                                                          chans)
            return f
        h["masked"] = flag(prod)
        h["masked_auto"] = flag(auto)

        # geometry from the production recipe (masked hits excluded)
        lat, (lo, hi), frame, aff, proj = station_frame(
            h[~h["masked"]], pads, tracks, z, dz, args)
        if aff is None:
            print(f"  no frame fit for {st}")
            summary["stations"][st] = {"error": "no frame fit",
                                       "latency": lat}
            continue
        px, py = aff.apply(proj[:, 0], proj[:, 1])
        tracks[f"px_{st}"] = px.astype(np.float32)
        tracks[f"py_{st}"] = py.astype(np.float32)

        # the pad each hit belongs to, and its distance to the prediction
        hp = h.merge(pads[["vmm", "ch", "pad_cx", "pad_cy", "channel_id"]],
                     on=["vmm", "ch"], how="left")
        hx = hp["pad_cx"].to_numpy(float)
        hy = hp["pad_cy"].to_numpy(float)
        ii = hp["idx"].to_numpy(np.int64)
        hp["dpred"] = np.hypot(hx - px[ii], hy - py[ii]).astype(np.float32)
        hp["in_win"] = (hp["dt"] >= lo) & (hp["dt"] <= hi)

        # nearest instrumented pad to the prediction -- the fiducial
        pad_xy = pads.dropna(subset=["pad_cx"])[["pad_cx", "pad_cy"]].to_numpy()
        pad_id = pads.dropna(subset=["pad_cx"])["channel_id"].to_numpy()
        d2 = ((px[:, None] - pad_xy[None, :, 0]) ** 2
              + (py[:, None] - pad_xy[None, :, 1]) ** 2)
        near = np.argmin(d2, axis=1)
        tracks[f"dpad_{st}"] = np.sqrt(d2.min(axis=1)).astype(np.float32)
        tracks[f"nearpad_{st}"] = pad_id[near].astype(np.int32)
        del d2

        # per-track reductions, for several hit selections
        def reduce_to_tracks(sub, prefix):
            d = np.full(n, np.inf, np.float32)
            adc = np.zeros(n, np.int16)
            dtb = np.full(n, np.nan, np.float32)
            cnt = np.zeros(n, np.int32)
            if len(sub):
                jj = sub["idx"].to_numpy(np.int64)
                dd = sub["dpred"].to_numpy(np.float32)
                ok = np.isfinite(dd)
                np.add.at(cnt, jj, 1)
                if ok.any():
                    np.minimum.at(d, jj[ok], dd[ok])
                    win = dd[ok] == d[jj[ok]]
                    sel = np.flatnonzero(ok)[win]
                    adc[jj[sel]] = sub["adc"].to_numpy()[sel]
                    dtb[jj[sel]] = sub["dt"].to_numpy()[sel]
            tracks[f"{prefix}_dmin_{st}"] = d
            tracks[f"{prefix}_adc_{st}"] = adc
            tracks[f"{prefix}_dt_{st}"] = dtb
            tracks[f"{prefix}_n_{st}"] = cnt

        win = hp[hp["in_win"]]
        reduce_to_tracks(win[~win["masked"]], "win")
        reduce_to_tracks(win, "winall")            # mask ignored
        reduce_to_tracks(hp[~hp["masked"]], "any")  # whole +-us window
        # unmapped hits (no pad) still count as "the station responded"
        tracks[f"nraw_win_{st}"] = np.bincount(
            win["idx"].to_numpy(np.int64), minlength=n).astype(np.int32)
        tracks[f"nraw_any_{st}"] = np.bincount(
            hp["idx"].to_numpy(np.int64), minlength=n).astype(np.int32)

        eff = float(((tracks[f"dpad_{st}"] < args.fid_r) & in_cap
                     & (tracks[f"win_dmin_{st}"] < args.probe_r)).sum()
                    / max(((tracks[f"dpad_{st}"] < args.fid_r) & in_cap).sum(),
                          1))
        print(f"  latency {lat['mu']:+.0f} ns, window {lo:.0f}..{hi:.0f}, "
              f"rotation {frame['affine_rotation_deg']:.2f} deg, "
              f"efficiency {eff:.4f}")
        summary["stations"][st] = {
            "z_mm": z, "vmms": vmms, "latency": lat, "window_ns": [lo, hi],
            "frame": frame, "mask_auto": {int(k): v for k, v in auto.items()},
            "mask_prod": {int(k): v for k, v in prod.items()},
            "n_hits_near": int(len(h)),
            "n_hits_masked": int(h["masked"].sum()),
            "efficiency_reproduced": eff}

        cols = ["idx", "dt", "vmm", "ch", "adc", "over_thr", "masked",
                "masked_auto", "channel_id", "dpred", "in_win"]
        hp[cols].to_parquet(f"{args.out}/autopsy_{tag}_hits_{st}.parquet",
                            index=False)
        del hp, h, win

    tracks.to_parquet(f"{args.out}/autopsy_{tag}_tracks.parquet", index=False)
    with open(f"{args.out}/autopsy_{tag}.json", "w") as f:
        json.dump(summary, f, default=float)
    print(f"wrote {args.out}/autopsy_{tag}_*")


if __name__ == "__main__":
    main()
