# 2026-07-24 — Mesh (gain) scan overnight + drift scan in progress (TB_July2026_H4)

Runs `beam_nominal_meshscan_1` (2026-07-23 18:17 → ~01:45, complete, QA done)
and `drift_scan_1` (started 12:28 today, **in progress** — points 450–750 of
450–900 recorded as of 17:20). **The mesh scan is clean: gain drops smoothly by
×3 over the 50 V scanned, saturation falls from 3.9 % to 0.4 %, and the
telescope-OR inefficiency (empty-trigger fraction) rises 10 % → 52 % — exactly
the efficiency-vs-gain curve this run was taken for.**

Beam much more intense than yesterday: **~4.6 kHz trigger rate averaged over a
run** (vs ~1200 Hz on 07-23), ~5.5 M triggers per 20-min sub-run.

Setup unchanged from 07-23 (5 detectors, external scintillator coincidence via
TCM, ZS, 16 samples × 60 ns, latency 32, Ar/iso 95/5) except where noted below.

---

## Run 3 — `beam_nominal_meshscan_1` (8 sub-runs × 20 min, overnight)

Two sub-runs at nominal HV (P2_IN 700/490, P2_MID/OUT 700/450 drift/mesh),
then the mesh stepped **down 10 V per sub-run with the drift stepped in
lockstep** (constant drift gap: 210 V for IN, 250 V for MID/OUT), from
in480/midout440 down to in430/midout390.

| sub-run | mesh IN / MID,OUT | P2_OUT med amp | P2_MID med amp | P2_IN med amp | P2_OUT sat | empty trig |
|---|---|---|---|---|---|---|
| nominal_00 | 490 / 450 | 710 | 512 | 657 | 3.9 % | 9.9 % |
| nominal_01 | 490 / 450 | 710 | 513 | 683 | 3.9 % | 10.2 % |
| meshscan_01 | 480 / 440 | 550 | 399 | 530 | 2.5 % | 15.5 % |
| meshscan_02 | 470 / 430 | 435 | 305 | 406 | 1.6 % | 23.1 % |
| meshscan_03 | 460 / 420 | 350 | 233 | 306 | 1.0 % | 39.8 % |
| meshscan_04 | 450 / 410 | 284 | 187 | 239 | 0.6 % | 41.8 % |
| meshscan_05 | 440 / 400 | 231 | 157 | 190 | 0.4 % | 51.5 % |
| meshscan_06 | 430 / 390 | — | — | — | — | — |

(amplitudes are median hit amplitude in ADC; "empty trig" = triggers with no
hit in **any** detector)

- **Gain slope: amplitude halves roughly every 30 V of mesh.** P2_OUT
  710 → 231 and P2_MID 512 → 157 over 50 V, i.e. a factor ~0.8 per 10 V step,
  consistent across all three P2 stations.
- **Repeatability is excellent**: nominal_00 vs nominal_01 agree to 0.1–0.4 %
  in amplitude — amplitude remains the reliable discriminant, as concluded from
  the latency scan.
- **P2_IN's problem is event count, not gain.** Its median amplitude (657–683
  ADC at nominal) is healthy and scales with HV exactly like MID/OUT, but its
  trigger share stays at ~0.5 % (vs 76 % / 55 % for OUT / MID). Whatever
  suppresses it removes events, it does not degrade the pulses it does see.
- P2_MID and P2_OUT: **0 dead channels in every sub-run**. The rising
  "dead-channel" counts for P2_IN and the uRWELLs at low gain are
  occupancy-starved channels (too few events to register), not real deads.
- The run-level QA verdict is FAIL, but that is the scan doing its job: the
  trigger-share and dead-channel checks trip on P2_IN (known issue) and on the
  intentionally low-gain points. Nothing new is wrong.
- **To do: `meshscan_06` has no telescope QA** (per-detector waveform QA
  finished 01:46 but `25_commissioning_qa` never ran on it) — re-run
  `run_beam_qa.sh beam_nominal_meshscan_1` to pick it up, and regenerate the
  trend panel so the scan curve includes the last point.

## Run 4 — `drift_scan_1` (10 points × 10 min, in progress)

Drift-field scan on P2_MID and P2_OUT: **drift 450 → 900 V in 50 V steps with
mesh fixed at 450 V**, i.e. drift-gap ΔV from 0 (zero field) to 450 V
(plateau). **P2_IN is off (HV 0) and out of the readout** — FEUs 1/4/5 only —
following yesterday's open item. uRWELLs unchanged (drift 600, resist 420).

Status as of 17:20: `drift_450` – `drift_750` recorded (~2.6 M triggers each,
~4.3 kHz), `drift_750` QA in progress; 800 / 850 / 900 remain (~1 h). Telescope
QA to be run once the scan completes.

**Two sub-runs have FEU dropouts — flag for the analysis:**

- `drift_450`: FEU 1 (both uRWELLs) missing 632 k of 2.61 M events (24 %).
- `drift_500`: FEU 1 missing 729 k (28 %), and **FEU 5 (P2_OUT) only recorded
  the last 369 k events** (event range 2.23–2.60 M) — this point is effectively
  lost for P2_OUT. **Consider retaking drift_500** after the scan finishes.
- `drift_550` onward is clean (0–22 missing events per FEU).

## Analysis

Same chain as 07-23: `run_beam_qa.sh <run_name>` on banco →
`TB_July2026_H4/analysis/telescope/<run>/...` (per sub-run rates / occupancy /
signal / HV plots + JSON/CSV, run-level trend panel
`trend_beam_nominal_meshscan_1.png`). Raw + decoded data under
`TB_July2026_H4/runs/<run>/`.
