# VMM ↔ DREAM stream matching (SPS H4, run_NN campaign)

Event-level matching between the two DAQs during the VMM campaign
(run_21..run_58+): the DREAM DAQ reads the uRWELL trackers and records
`(eventId, trigger_timestamp_ns)`; the VMM/SRS DAQ sees the same trigger
signal on VMM 0 channel 44. DREAM vetoes triggers while busy, so the VMM
trigger stream is a superset (~1.4×) of the DREAM event stream.

## Clock model

`t_vmm = a · t_dream + b`, per spill.

* **SRS tick is 22.5 ns, not 25 ns.** `vmm_decode.derive()` computes
  `abs_time_ns = srs_timestamp * 25.0 + …`, which makes the SPS spill
  period read 16.00 s instead of the true 14.40 s (exactly 10/9). The SRS
  timestamp counts the same 44.44 MHz clock as the BCID counter. Any
  absolute-time use of `abs_time_ns` elsewhere inherits this 11% stretch.
* With the tick fixed, the residual clock drift is ≈ +0.7 ppm — but even
  1 ppm accumulates to ~ms over a run, hence per-spill locking.
* `b` contains a trigger-path latency of O(500 µs), constant within a
  spill; markers arrive every 65536 ticks (1.47 ms) and the marker values
  are exact multiples, so within-spill VMM time is hardware-precise.

## Algorithm (match_streams.py)

1. Seed `(a, b)` from spill-start pairing (robust nearest-start,
   outlier-clipped — survives DAQs starting mid-spill and noise blips).
2. Walk spills in time order. For each spill: residual-histogram lag scan
   (±5 ms; ±50 ms until first lock) centred on the previous spill's fit,
   then greedy one-to-one nearest matching + robust linear refit at
   20 µs → 5 µs → 2 µs tolerance.

Run on lxplus (needs LCG for uproot):

    python3 match_streams.py run_32 meshscan_m00V --out out --tol-us 2

Outputs `match_<run>_<sub>.json` (per-spill summary) and `.npz` with the
`eventId ↔ t_vmm_ns` pairs — the input for uRWELL-track-referenced VMM
efficiency.

## Validation status (2026-08-16)

* run_32 (all meshscan sub_runs): **98.5–99.4% of DREAM events matched,
  residual rms ~10 ns**; unmatched VMM fraction ≈ 30% = DREAM busy veto.
* run_29 / run_31 (trigtest) and run_36: no event-level coherence found at
  all (both streams are internally µs-faithful — SPS revolution-period
  microstructure visible in each — but mutually incoherent even at
  ±30 ms). Cause not yet identified; possibly a different signal on ch 44
  during those periods.
* Data prerequisite: `hits_store/<capture>/*.npy` columns. Missing (need a
  re-decode of the raw pcapng without `--drop-columns`): runs ≤22, 24, 26,
  28, 37, 44, 57 (columns dropped), and everything from run_59 on (gas B).
