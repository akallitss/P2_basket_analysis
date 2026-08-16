#!/usr/bin/env python3
"""Match the VMM trigger stream to the DREAM event stream for one sub_run.

VMM side  : vmm/runs/<run>/<sub>/hits_store/<capture>/{vmm,ch,offset,bcid,tdc,
            srs_timestamp}.npy  -> trigger-channel hits (VMM 0 ch 44),
            t = srs_timestamp*22.5 + (offset*4096 + bcid)*22.5 + TDC fine.
            NOTE the SRS timestamp tick is 22.5 ns (44.44 MHz, same clock as
            the BCID counter), NOT the 25 ns that vmm_decode.derive assumes —
            established by comparing SPS spill periods against the DREAM
            clock (16.00 s apparent vs 14.40 s true = exactly 10/9).
DREAM side: runs/<run>/<sub>/combined_hits_root/*.root ->
            unique (eventId, trigger_timestamp_ns)  (DREAM clock, 5 ns quanta).

Model: t_vmm = a * t_dream + b. a differs from 1 by only ~1 ppm but that
still accumulates to many ms over a run; b contains a trigger-path latency
of O(500 us) that is constant within a spill. Stages:
  1. seed (a, b) from paired spill starts (robust nearest-start pairing)
  2. walk spills in time order; per spill re-lock the lag with a residual-
     histogram scan around the previous spill's fit, then greedy nearest-
     neighbour match + robust linear refit at 20 us -> 5 us -> tol_us
Result on a good run: >99% of DREAM events matched, residual rms ~10 ns
(DREAM 5 ns quantisation + VMM 22.5 ns BCID clock). The ~30% of VMM
triggers with no DREAM partner are the DREAM busy-veto losses.
Outputs a JSON summary + npz of matched pairs (eventId <-> t_vmm).
"""
import argparse, glob, json, os, sys
import numpy as np

BASE = "/eos/experiment/ntof/data/x17/p2_sps_july"
CLOCK_PERIOD_NS = 22.5
# The SRS timestamp ticks on the same 44.44 MHz clock as the BCID counter:
# spill-period comparison against the DREAM 5-ns hardware clock gives
# 25/22.500011 exactly (the decoder's abs_time_ns wrongly assumes 25 ns).
SRS_TICK_NS = 22.5
TAC_SLOPE_NS = 60.0
TDC_RANGE = 255
TRIG_VMM, TRIG_CH = 0, 44
DEDUP_NS = 1500.0


def load_vmm_triggers(run, sub):
    store = f"{BASE}/vmm/runs/{run}/{sub}/hits_store"
    caps = sorted(glob.glob(store + "/*/"))
    if not caps:
        sys.exit(f"no hits_store captures under {store}")
    if not os.path.exists(caps[0] + "vmm.npy"):
        sys.exit(f"hits_store at {store} has no hit columns (dropped at "
                 "reduce time) — re-decode this run without --drop-columns")
    ts = []
    for cap in caps:
        try:
            vmm = np.load(cap + "vmm.npy")
            ch = np.load(cap + "ch.npy")
        except Exception as e:
            print(f"  skip {os.path.basename(cap.rstrip('/'))}: {e}")
            continue
        m = (vmm == TRIG_VMM) & (ch == TRIG_CH)
        if not m.any():
            continue
        off = np.load(cap + "offset.npy")[m].astype(np.float64)
        bcid = np.load(cap + "bcid.npy")[m].astype(np.float64)
        tdc = np.load(cap + "tdc.npy")[m].astype(np.float64)
        srs = np.load(cap + "srs_timestamp.npy")[m].astype(np.float64)
        t = (srs * SRS_TICK_NS + (off * 4096 + bcid) * CLOCK_PERIOD_NS
             + (CLOCK_PERIOD_NS - tdc * TAC_SLOPE_NS / TDC_RANGE))
        ts.append(t)
    t = np.sort(np.concatenate(ts))
    # dedup: a trigger pulse can fire several VMM samples; keep first of each burst
    keep = np.ones(len(t), bool)
    keep[1:] = np.diff(t) > DEDUP_NS
    return t[keep], len(t)


def load_dream_triggers(run, sub):
    import uproot
    files = sorted(glob.glob(f"{BASE}/runs/{run}/{sub}/combined_hits_root/*.root"))
    if not files:
        sys.exit(f"no combined_hits_root files for {run}/{sub}")
    eid, tts = [], []
    for fn in files:
        with uproot.open(fn) as f:
            t = f["hits"]
            a = t.arrays(["eventId", "trigger_timestamp_ns"], library="np")
        # one row per hit -> unique events
        e, idx = np.unique(a["eventId"], return_index=True)
        eid.append(e)
        tts.append(a["trigger_timestamp_ns"][idx])
    eid = np.concatenate(eid)
    tts = np.concatenate(tts).astype(np.float64)
    o = np.argsort(tts)
    return eid[o], tts[o]


def rate_profile(t, bin_ns, t0, t1):
    edges = np.arange(t0, t1 + bin_ns, bin_ns)
    h, _ = np.histogram(t, bins=edges)
    return h, edges


def spill_starts(t, gap_s=2.0, min_n=1000):
    """Start time of each real spill (noise blips with < min_n triggers dropped)."""
    g = np.where(np.diff(t) > gap_s * 1e9)[0]
    s = np.concatenate([[0], g + 1])
    e = np.concatenate([g + 1, [len(t)]])
    keep = (e - s) >= min_n
    return t[s[keep]]


def seed_from_spills(tv, td):
    """Fit t_vmm = a * t_dream + b on paired spill starts.

    The two DAQs run on different clock frequencies (SRS timestamp ticks are
    22.5 ns, not the 25 ns the decoder's abs_time_ns assumes), so a is far
    from 1 (~0.9). Pair spills by trying small index shifts and keeping the
    pairing whose linear fit has the smallest residual.
    """
    sv, sd = spill_starts(tv), spill_starts(td)
    # rough start: index-shift pairing on the first 20 spills only
    best = None
    for k in range(-4, 5):
        i0, j0 = max(0, k), max(0, -k)
        n = min(len(sv) - i0, len(sd) - j0, 20)
        if n < 3:
            continue
        x, y = sd[j0:j0 + n], sv[i0:i0 + n]
        A = np.polyfit(x, y, 1)
        rms = float(np.std(y - np.polyval(A, x)))
        if best is None or rms < best[0]:
            best = (rms, A[0], A[1], k)
    rms, a, b, k = best
    # refine: nearest-start pairing under the current map, outlier-clipped,
    # iterated — survives noise blips and spills missing on one side
    n_used = 0
    for _ in range(4):
        pred = a * sd + b
        idx = np.clip(np.searchsorted(sv, pred), 1, len(sv) - 1)
        near = np.where(pred - sv[idx - 1] < sv[idx] - pred, idx - 1, idx)
        res = sv[near] - pred
        mad = 1.4826 * np.median(np.abs(res - np.median(res)))
        keep = np.abs(res - np.median(res)) < max(3 * mad, 5e7)
        if keep.sum() < 3:
            break
        A = np.polyfit(sd[keep], sv[near[keep]], 1)
        a, b = A[0], A[1]
        rms = float(np.std(sv[near[keep]] - np.polyval(A, sd[keep])))
        n_used = int(keep.sum())
    print(f"[spill seed] {len(sv)} vmm / {len(sd)} dream spills, shift {k}, "
          f"{n_used} paired after clip: a = {a:.9f}, "
          f"spill-start residual rms {rms/1e3:.1f} us")
    return a, b


def lag_scan(tv, td, a, b, half=5e6, bin_ns=1000.0, sample=20, quiet=False):
    """Direct residual-histogram lag search: the trigger path has an
    O(500 us) latency, wandering by hundreds of us over a run, that the
    spill-edge seed cannot see. Returns the lag of the sharp peak to add
    to b."""
    x = a * td[::sample] + b
    lo = np.searchsorted(tv, x - half)
    hi = np.searchsorted(tv, x + half)
    res = np.concatenate([tv[l:h] - xi for xi, l, h in zip(x, lo, hi)])
    nb = int(2 * half / bin_ns)
    h, edges = np.histogram(res, bins=nb, range=(-half, half))
    k = np.argmax(h)
    lag = 0.5 * (edges[k] + edges[k + 1])
    if not quiet:
        print(f"[lag scan] peak at {lag/1e3:+.1f} us "
              f"({h[k]} entries vs mean {h.mean():.1f})")
    return lag


def greedy_match(tv, td_mapped, tol_ns):
    """For each mapped DREAM time, nearest VMM trigger within tol. Enforce
    one-to-one by keeping the better of colliding claims."""
    idx = np.searchsorted(tv, td_mapped)
    idx = np.clip(idx, 1, len(tv) - 1)
    left = tv[idx - 1]; right = tv[idx]
    use_left = (td_mapped - left) < (right - td_mapped)
    j = np.where(use_left, idx - 1, idx)
    res = td_mapped - tv[j]
    ok = np.abs(res) < tol_ns
    # one-to-one: among DREAM events claiming same VMM trigger keep min |res|
    order = np.lexsort((np.abs(res), j))
    j_s = j[order]
    first = np.ones(len(j_s), bool)
    first[1:] = np.diff(j_s) != 0
    keep = np.zeros(len(j), bool)
    keep[order[first]] = True
    ok &= keep
    return j, res, ok


def robust_linear(x, y, ok, n_iter=5, clip=5.0):
    a, b = 1.0, 0.0
    sel = ok.copy()
    for _ in range(n_iter):
        if sel.sum() < 10:
            break
        A = np.polyfit(x[sel], y[sel], 1)
        a, b = A[0], A[1]
        r = y - (a * x + b)
        s = 1.4826 * np.median(np.abs(r[sel] - np.median(r[sel])))
        s = max(s, 1.0)
        sel = ok & (np.abs(r - np.median(r[sel])) < clip * s)
    return a, b, sel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("sub")
    ap.add_argument("--out", default=".")
    ap.add_argument("--tol-us", type=float, default=50.0,
                    help="match tolerance for the final pass")
    args = ap.parse_args()

    print(f"[vmm] loading triggers {args.run}/{args.sub}")
    tv, n_raw = load_vmm_triggers(args.run, args.sub)
    print(f"  {n_raw} trigger hits -> {len(tv)} dedup'd triggers, "
          f"span {(tv[-1]-tv[0])/1e9:.1f} s")

    print("[dream] loading events")
    eid, td = load_dream_triggers(args.run, args.sub)
    print(f"  {len(td)} events, span {(td[-1]-td[0])/1e9:.1f} s, "
          f"t range [{td[0]:.0f}, {td[-1]:.0f}] ns")
    print(f"  vmm t range [{tv[0]:.3e}, {tv[-1]:.3e}] ns")

    # seed a (clock ratio) + b from spill-start pairing. The seed is only
    # good to O(10 ppm) / O(ms); the per-spill sequential tracker below does
    # the real locking.
    a, b = seed_from_spills(tv, td)

    # ---- per-spill sequential tracking: within a spill the clock relation
    # is linear to ~10 ns, but the seed's a-error accumulates to many ms over
    # a long run. So walk the spills in time order: each spill re-locks the
    # lag around the PREVIOUS spill's fit (extrapolation error ~us over the
    # 14.4 s SPS period), with a wide scan only for the first spill.
    gaps = np.where(np.diff(td) > 2e9)[0]
    starts = np.concatenate([[0], gaps + 1])
    ends = np.concatenate([gaps + 1, [len(td)]])
    j = np.zeros(len(td), np.int64)
    ok = np.zeros(len(td), bool)
    r = np.full(len(td), np.nan)
    spills = []
    a_cur, b_cur, first = a, b, True
    for s, e in zip(starts, ends):
        if (e - s) < 100:
            continue
        tds = td[s:e]
        half = 50e6 if first else 5e6
        a_s = a_cur
        b_s = b_cur + lag_scan(tv, tds, a_cur, b_cur, half=half,
                               bin_ns=500.0, sample=5, quiet=True)
        js = None
        for tol in (20e3, 5e3, args.tol_us * 1e3):
            tm = a_s * tds + b_s
            js, _, oks = greedy_match(tv, tm, tol)
            if oks.sum() < 50:
                break
            a_s, b_s, _ = robust_linear(tds, tv[js], oks)
        if js is None or oks.sum() < 50:
            continue
        if oks.mean() > 0.5:            # good lock: track from here
            a_cur, b_cur, first = a_s, b_s, False
        rr = tv[js] - (a_s * tds + b_s)
        j[s:e] = js; ok[s:e] = oks; r[s:e] = rr
        m = oks
        spills.append(dict(
            n_dream=int(e - s), n_matched=int(m.sum()),
            frac=float(m.mean()),
            t0_s=float((td[s] - td[0]) / 1e9),
            dur_s=float((tds[-1] - tds[0]) / 1e9),
            a_spill=float(a_s),
            b_minus_global_us=float((b_s - b) / 1e3),
            res_med_ns=float(np.median(rr[m])) if m.any() else None,
            res_rms_ns=float(np.std(rr[m])) if m.any() else None))
    fr = [sp["frac"] for sp in spills]
    print(f"[per-spill] {len(spills)} spills refit: match frac "
          f"min {min(fr):.3f} med {np.median(fr):.3f} max {max(fr):.3f}; "
          f"residual rms med "
          f"{np.median([sp['res_rms_ns'] for sp in spills]):.1f} ns")

    # how many VMM triggers went unmatched (expected: DREAM busy veto)
    used = np.zeros(len(tv), bool)
    used[j[ok]] = True
    summary = dict(
        run=args.run, sub=args.sub,
        n_vmm_trig_hits=int(n_raw), n_vmm_triggers=int(len(tv)),
        n_dream_events=int(len(td)),
        n_matched=int(ok.sum()),
        match_frac_dream=float(ok.sum() / len(td)),
        vmm_unmatched=int((~used).sum()),
        vmm_unmatched_frac=float((~used).mean()),
        clock_ratio_a=float(np.median([sp["a_spill"] for sp in spills])),
        drift_ppm=float((np.median([sp["a_spill"] for sp in spills]) - 1)
                        * 1e6),
        offset_seed_ns=float(b),
        residual_rms_ns=float(np.std(r[ok])),
        residual_med_ns=float(np.median(r[ok])),
        n_spills=len(spills), spills=spills)

    os.makedirs(args.out, exist_ok=True)
    tag = f"{args.run}_{args.sub}"
    with open(f"{args.out}/match_{tag}.json", "w") as f:
        json.dump(summary, f, indent=1)
    np.savez_compressed(
        f"{args.out}/match_{tag}.npz",
        event_id=eid[ok], t_dream_ns=td[ok], t_vmm_ns=tv[j[ok]],
        residual_ns=r[ok], a=a, b=b,
        tv_all=tv, td_all=td, matched_mask=ok, vmm_used=used)
    print(json.dumps({k: v for k, v in summary.items() if k != "spills"},
                     indent=1))
    print(f"wrote {args.out}/match_{tag}.json + .npz")


if __name__ == "__main__":
    main()
