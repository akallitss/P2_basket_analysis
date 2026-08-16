#!/usr/bin/env python3
"""Match the VMM trigger stream to the DREAM event stream for one sub_run.

VMM side  : trigger-channel hits (VMM 0 ch 44) via extract_vmm_triggers --
            the reduced hits_store columns where they exist, the raw pcapng
            otherwise (identical results),
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
  1. first lock: try pairing the densest DREAM spill with every VMM spill and
     keep the offset whose coincidence histogram spikes (spill-edge proximity
     alone picks the wrong spill on some runs)
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

import extract_vmm_triggers as evt

BASE = evt.BASE


def load_vmm_triggers(run, sub, npz=None, source="auto"):
    """Trigger times: from a pre-extracted npz, else straight from the data.

    extract_vmm_triggers reads the reduced column store when it has the hit
    columns and decodes the raw pcapng when it does not -- the two paths are
    bit-identical (verified on run_32/meshscan_m00V), so no sub_run of the
    campaign is out of reach.
    """
    if npz:
        d = np.load(npz)
        meta = json.loads(str(d["meta"]))
        return (d["t_ns"], meta["n_trigger_hits"], meta["source"],
                (meta.get("trigger_vmm", evt.TRIG_VMM),
                 meta.get("trigger_ch", evt.TRIG_CH)))
    sub_dir = evt.sub_dir_of(run, sub)
    src = evt.source_of(sub_dir)
    if src == "none":
        sys.exit(f"no VMM data for {run}/{sub} under {sub_dir}")
    times, _ = evt.collect(sub_dir, [(evt.TRIG_VMM, evt.TRIG_CH)],
                           source=source)
    t = times[(evt.TRIG_VMM, evt.TRIG_CH)]
    if not len(t):
        sys.exit(f"no hits on the trigger channel for {run}/{sub}")
    return evt.dedup(t), len(t), src, (evt.TRIG_VMM, evt.TRIG_CH)


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


def residual_peak(tv, td_sample, lag, half, bin_ns):
    """Peak of the pairwise-residual histogram of (tv - (td_sample + lag)).

    Returns (peak_lag_ns, peak_height, background_mean). Batched, so the
    window may be wide: at 1 us bins over +-0.6 s the true coincidences pile
    into a single bin while the random background stays at
    n_dream * n_vmm_in_window / n_bins.
    """
    nb = max(int(2 * half / bin_ns), 1)
    hist = np.zeros(nb, np.int64)
    lo = np.searchsorted(tv, td_sample + lag - half)
    hi = np.searchsorted(tv, td_sample + lag + half)
    for k in range(0, len(td_sample), 200):
        sl = slice(k, k + 200)
        parts = [tv[l:h] - (x + lag) for x, l, h
                 in zip(td_sample[sl], lo[sl], hi[sl])]
        if not parts:
            continue
        res = np.concatenate(parts)
        if len(res):
            hist += np.histogram(res, bins=nb, range=(-half, half))[0]
    j = int(np.argmax(hist))
    return lag + (j + 0.5) * bin_ns - half, int(hist[j]), float(hist.mean())


def first_lock(tv, td, verbose=True, max_lag_ns=120e9, min_sigma=25.0,
               pair_budget=2e7):
    """Find the absolute VMM-DREAM offset with no prior, by testing every
    plausible spill pairing.

    The spill-start linear seed is not reliable on its own: the first "spill"
    of a stream is often just the DAQ start, and one bad point can pull the
    fit onto a wrong spill pairing (run_57 was seeded 31 s away from the
    truth). Coincidence is a much stronger discriminator than spill-edge
    proximity, so pair one dense DREAM spill against each VMM spill in turn
    and keep the offset whose residual histogram actually spikes.
    """
    sv, sd = spill_starts(tv), spill_starts(td)
    g = np.where(np.diff(td) > 2e9)[0]
    s = np.concatenate([[0], g + 1])
    e = np.concatenate([g + 1, [len(td)]])
    n = e - s
    i = int(np.argmax(n))                     # the densest DREAM spill
    tds = td[s[i]:e[i]]
    # Cost per candidate is n_sample x (VMM triggers inside the +-0.6 s
    # window), and a junk trigger line can run at 200 kHz -- 100x a real one.
    # Scale the sample to keep that product bounded; 200 points still put a
    # true peak far above the background.
    half = 0.6e9
    rate = len(tv) / max(tv[-1] - tv[0], 1.0)         # triggers per ns
    n_sample = int(np.clip(pair_budget / max(2 * half * rate, 1.0),
                           200, 2000))
    step = max(1, len(tds) // n_sample)
    sample = tds[::step]

    def scan(lags):
        best = None
        for lag in lags:
            pk_lag, pk, mu = residual_peak(tv, sample, lag, half, 1e3)
            s = (pk - mu) / max(np.sqrt(max(mu, 1.0)), 1.0)
            if best is None or s > best[0]:
                best = (s, pk_lag, pk, mu)
        return best

    # the two DAQs are started together, so the offset is seconds, not
    # minutes: try the nearby spills first and only fall back to all of them
    # (which on a 200-spill run costs 20x more) if nothing convincing turns up
    all_lags = sv - tds[0]
    near = all_lags[np.abs(all_lags) < max_lag_ns]
    best = scan(near) if len(near) else None
    # the full fallback is for the case where the DAQs were started minutes
    # apart, which has not been seen (offsets measured so far: 0.8-6.4 s).
    # On a long run it is 200+ candidates, so only pay for it when the run is
    # short enough that it is cheap.
    if ((best is None or best[0] < min_sigma) and len(all_lags) > len(near)
            and len(all_lags) <= 40):
        cand = [b for b in (best, scan(all_lags)) if b is not None]
        best = max(cand, key=lambda b: b[0]) if cand else None
    if best is None:
        raise SystemExit("no VMM spill to lock on")
    sig, lag, pk, mu = best
    if verbose:
        print(f"[first lock] {len(sv)} vmm / {len(sd)} dream spills; densest "
              f"dream spill {i} ({len(tds)} ev, {len(sample)} sampled): "
              f"offset {lag/1e9:+.6f} s, peak {pk} vs background {mu:.1f} "
              f"({sig:.0f} sigma)")
    return 1.0, lag, sig


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
    ap.add_argument("--vmm-npz", default=None,
                    help="pre-extracted trigger times (extract_vmm_triggers)")
    ap.add_argument("--vmm-source", choices=["auto", "store", "pcapng"],
                    default="auto")
    args = ap.parse_args()

    print(f"[vmm] loading triggers {args.run}/{args.sub}")
    tv, n_raw, vmm_source, chan = load_vmm_triggers(
        args.run, args.sub, args.vmm_npz, args.vmm_source)
    print(f"  {n_raw} trigger hits on vmm {chan[0]} ch {chan[1]} -> "
          f"{len(tv)} dedup'd triggers [{vmm_source}], "
          f"span {(tv[-1]-tv[0])/1e9:.1f} s")

    print("[dream] loading events")
    eid, td = load_dream_triggers(args.run, args.sub)
    print(f"  {len(td)} events, span {(td[-1]-td[0])/1e9:.1f} s, "
          f"t range [{td[0]:.0f}, {td[-1]:.0f}] ns")
    print(f"  vmm t range [{tv[0]:.3e}, {tv[-1]:.3e}] ns")

    # absolute offset from the best spill pairing, judged by coincidence.
    # Good to O(us) at the spill it locked on; the per-spill tracker below
    # picks up the clock drift from there.
    a, b, lock_sigma = first_lock(tv, td)

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
    if spills:
        fr = [sp["frac"] for sp in spills]
        print(f"[per-spill] {len(spills)} spills refit: match frac "
              f"min {min(fr):.3f} med {np.median(fr):.3f} max {max(fr):.3f}; "
              f"residual rms med "
              f"{np.median([sp['res_rms_ns'] for sp in spills]):.1f} ns")
    else:
        print("[per-spill] NO spill locked — streams look mutually incoherent")

    # how many VMM triggers went unmatched (expected: DREAM busy veto)
    used = np.zeros(len(tv), bool)
    used[j[ok]] = True
    a_med = (float(np.median([sp["a_spill"] for sp in spills])) if spills
             else float("nan"))
    summary = dict(
        run=args.run, sub=args.sub, vmm_source=vmm_source,
        trigger_vmm=int(chan[0]), trigger_ch=int(chan[1]),
        n_vmm_trig_hits=int(n_raw), n_vmm_triggers=int(len(tv)),
        n_dream_events=int(len(td)),
        n_matched=int(ok.sum()),
        match_frac_dream=float(ok.sum() / len(td)),
        vmm_unmatched=int((~used).sum()),
        vmm_unmatched_frac=float((~used).mean()),
        clock_ratio_a=a_med,
        drift_ppm=(a_med - 1) * 1e6,
        offset_seed_ns=float(b), lock_sigma=float(lock_sigma),
        residual_rms_ns=float(np.std(r[ok])) if ok.any() else None,
        residual_med_ns=float(np.median(r[ok])) if ok.any() else None,
        n_spills=len(spills), spills=spills)

    os.makedirs(args.out, exist_ok=True)
    tag = f"{args.run}_{args.sub}"
    with open(f"{args.out}/match_{tag}.json", "w") as f:
        json.dump(summary, f, indent=1)
    # One row per DREAM event, matched or not (the unmatched ones are the
    # denominator of any efficiency built on this), plus which VMM triggers
    # were consumed. The VMM trigger times themselves live in the companion
    # vmm_triggers_<tag>.npz, indexed by vmm_index -- no need to duplicate
    # them here.
    vmm_index = np.where(ok, j, -1).astype(np.int64)
    np.savez_compressed(
        f"{args.out}/match_{tag}.npz",
        event_id=eid, t_dream_ns=td, matched=ok, vmm_index=vmm_index,
        t_vmm_ns=np.where(ok, tv[j], np.nan),
        residual_ns=r.astype(np.float32),
        vmm_used=used, n_vmm_triggers=len(tv), a_seed=a, b_seed=b)
    print(json.dumps({k: v for k, v in summary.items() if k != "spills"},
                     indent=1))
    print(f"wrote {args.out}/match_{tag}.json + .npz")


if __name__ == "__main__":
    main()
