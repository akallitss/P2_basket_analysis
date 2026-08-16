# VMM ↔ DREAM stream matching (SPS H4, run_NN campaign)

Event-level matching between the two DAQs of the July/August 2026 SPS test.
The DREAM DAQ reads the uRWELL trackers and records
`(eventId, trigger_timestamp_ns)`; the same external trigger is fanned into
the VMM3a/SRS DAQ on **VMM 0, channel 44**. DREAM vetoes triggers while it is
busy, so the VMM trigger stream is a superset (~1.35×) of the DREAM event
stream.

The matching is a **one-off**: it is run once over the whole campaign and the
result is stored permanently on EOS under

    /eos/experiment/ntof/data/x17/p2_sps_july/vmm/matching/
        vmm_triggers_<run>_<sub>.npz   every VMM trigger time of the sub_run
        match_<run>_<sub>.json         per-spill fit + match summary
        match_<run>_<sub>.npz          one row per DREAM event

`match_*.npz` holds `event_id`, `t_dream_ns`, `matched`, `t_vmm_ns`,
`residual_ns` and `vmm_index` (the index into the companion trigger file,
−1 when unmatched), plus `vmm_used` — which VMM triggers were consumed, i.e.
the complement is the DREAM busy-veto loss. **Unmatched DREAM events are kept**
in the table: they are the denominator of any efficiency built on this.

## Clock model

`t_vmm = a · t_dream + b`, refit per spill.

* **The SRS timestamp tick is 22.5 ns, not 25 ns.** `vmm_decode.derive()`
  computes `abs_time_ns = srs_timestamp * 25.0 + …`, which makes the SPS
  spill period read 16.00 s instead of the true 14.40 s (exactly 10/9). The
  SRS timestamp counts the same 44.44 MHz clock as the BCID counter. Any
  absolute-time use of `abs_time_ns` elsewhere inherits this 11% stretch.
* With the tick fixed, the residual clock drift is ≈ +0.7 ppm — but even
  1 ppm accumulates to ~ms over a run, hence per-spill locking.
* `b` contains a trigger-path latency of O(500 µs), constant within a spill.
  Markers arrive every 65536 ticks (1.47 ms) and the marker values are exact
  multiples, so within-spill VMM time is hardware-precise.

## The two scripts

**`extract_vmm_triggers.py <run> <sub>`** — pulls the trigger-channel times
out of the VMM data. It reads the reduced column store
(`hits_store/<capture>/*.npy`) when the online pass kept the hit columns, and
decodes the raw `raw_daq_data/*.pcapng` otherwise, via
`vmm_decode.iter_chunks` (put `P2_basket_online_analysis` on `PYTHONPATH`).
Most of the campaign was reduced with `--drop-columns`, so the pcapng path is
what makes the other half reachable — and it is cheap (~1 s per 60 MB
capture, since every chunk is filtered to channel 44 and dropped). The two
paths were verified bit-identical on run_32/meshscan_m00V.

**`find_trigger_channel.py <run> <sub>`** — the cabling is not the same all
campaign long, and nothing documents it, so read it off the data: take the
busiest channels and ask each whether its hit times coincide with DREAM
triggers. The real trigger channel answers at >100 σ; everything else sits at
the random-coincidence floor. `sweep_matching.sh` calls this automatically
whenever the default channel gives no usable lock.

**`match_streams.py <run> <sub>`** — the matcher:

1. *First lock.* Pair the densest DREAM spill against **every** VMM spill in
   turn and keep the offset whose coincidence histogram spikes. Judging the
   pairing by spill-edge proximity instead is not safe — the first "spill" of
   a stream is often just the DAQ start, and that one bad point seeded run_57
   31 s away from the truth (3% matched). Coincidence is decisive: the right
   offset stands out at 100–300 σ.
2. *Per-spill tracking.* Walk the spills in time order; each one re-locks the
   lag with a residual-histogram scan (±5 ms) centred on the previous spill's
   fit, then greedy one-to-one nearest matching and a robust linear refit at
   20 µs → 5 µs → `--tol-us`.

**Do not** try to find the offset with a global nearest-neighbour fit or an
FFT rate cross-correlation. Both lock onto the symmetric random-coincidence
background (median 0, no autocorrelation) or onto the spill envelope, and
bury the true peak; the FFT peak sits at a multiple of the spill period and
is only ~3 σ. Everything here is built on direct pairwise-residual histograms
for that reason.

## Running it

`sweep_matching.sh` does the whole campaign on lxplus, 6 at a time,
resumable (a sub_run whose json is already on EOS is skipped):

    ./sweep_matching.sh --jobs 6                # everything not yet done
    ./sweep_matching.sh --only run_57 --force   # redo one run
    ./sweep_matching.sh --dry-run               # just print the work list

It writes to `$TMPDIR` and pushes with `xrdcp` — the AFS home is 5 GB and a
single match npz is tens of MB, so nothing of size may land there.

`tabulate_matches.py` reads all the summaries back and writes
`matching_table.csv` (+ a markdown table) — see `MATCHING_TABLE.md` for the
committed result.

## Result

**110 sub_runs: 77 OK, 1 partial, 32 without a lock** — see
`MATCHING_TABLE.md` for the per-sub_run table and the reasons. 213.2 M DREAM
events matched, median residual rms **10.4 ns** (DREAM's 5 ns hardware
quantisation plus the VMM 22.5 ns BCID clock), median clock drift
**+0.67 ppm**, and 21–63% of VMM triggers left spare as the DREAM busy veto.
Every run from **run_31 onward** is matched; every sub_run without a lock is
**run_29 or earlier**, where the trigger line does not carry the DREAM
trigger at all.

A sub_run that fails to lock is unmistakable in the numbers: the match
fraction collapses to the accidental rate (a few %) and the residual rms
sits at the ~1 µs random-coincidence floor. Where it locks, it locks hard —
113 to 606 σ.

Two things that look like matching failures and are not: the two DAQs are
not stopped together (so judge on `match_frac_covered`, over the DREAM
events the VMM was recording, not on the raw fraction), and the start-of-run
offset between them ranges from −25 s to +88 s.
