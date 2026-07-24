# 2026-07-24 — Mesh (gain) scans + drift scan (TB_July2026_H4)

Runs `beam_nominal_meshscan_1` (2026-07-23 18:17 → ~01:45, complete, QA done),
`drift_scan_1` + `drift_scan_2` (12:28 → 15:10 today; split in two because the
HV max-voltage limit stopped scan 1 at the 800 V point — see below), and
`meshscan_fine_1` (18:38 → ~20:46 today, 5 V steps on MID/OUT), and
`p2in_check_1` (20:49 → 21:03, P2_IN at reduced HV). **The coarse mesh scan is
clean: gain drops smoothly by ×3 over the 50 V scanned, saturation falls from
3.9 % to 0.4 %, and the telescope-OR inefficiency (empty-trigger fraction)
rises 10 % → 52 % — exactly the efficiency-vs-gain curve this run was taken
for. The full drift curve 450–900 V is on disk, the fine mesh scan fills in the
5 V grid from 445 down to 390, and — best news of the day — P2_IN comes alive
at mesh 400 V: ~5× more events with a hit than at its 490 V nominal, so the
"dead" station was an HV working-point problem, not readout.**

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

## Run 4 — `drift_scan_1` + `drift_scan_2` (10 points × 10 min, complete)

Drift-field scan on P2_MID and P2_OUT: **drift 450 → 900 V in 50 V steps with
mesh fixed at 450 V**, i.e. drift-gap ΔV from 0 (zero field) to 450 V
(plateau). **P2_IN is off (HV 0) and out of the readout** — FEUs 1/4/5 only —
following yesterday's open item. uRWELLs unchanged (drift 600, resist 420).

`drift_scan_1` (12:28) took `drift_450` – `drift_750` (~2.6 M triggers each,
~4.3 kHz), then **stalled on the ramp to 800 V: the drift channels' max-voltage
limit on the HV crate was 750 V**, so with V0 = 800 the channels plateaued at
750.25 V measured and the run sat in ramp-wait from 13:43 until stopped
manually (~14:20). The `drift_800` sub-run of scan 1 contains only an
hv_monitor.csv and an empty raw dir — ignore it.

**Fix: raised the max-HV limit, then continued the remaining points as
`drift_scan_2`** (14:38, scan window overridden to 800–900). All three
sub-runs — `drift_800` / `drift_850` / `drift_900` — completed cleanly (~3 GB
raw each, same statistics as the other points); run finished normally at 15:10.

Per-sub-run waveform QA is still working through the backlog (~1 sub-run/h);
telescope QA to be run on both runs once it drains.

**Two scan-1 sub-runs have FEU dropouts — flag for the analysis:**

- `drift_450`: FEU 1 (both uRWELLs) missing 632 k of 2.61 M events (24 %).
- `drift_500`: FEU 1 missing 729 k (28 %), and **FEU 5 (P2_OUT) only recorded
  the last 369 k events** (event range 2.23–2.60 M) — this point is effectively
  lost for P2_OUT. **Consider retaking drift_500** after the scan finishes.
- `drift_550` onward is clean (0–22 missing events per FEU).

## Run 5 — `meshscan_fine_1` (12 points × 10 min, evening)

Fine mesh scan on **P2_MID and P2_OUT only** (P2_IN off and out of the readout,
FEUs 1/4/5, uRWELLs unchanged): **mesh 445 → 390 V in 5 V steps**, drift in
lockstep as usual (constant 250 V gap), i.e. the 5 V grid interleaving the
coarse scan's 10 V points. Started 18:38, last point (`meshscan_12_midout390`)
finished ~20:46; ~3 GB raw per sub-run — same statistics as the day's other
runs (~2.6 M triggers / point at ~4.5 kHz).

- **The `nominal_00` reference sub-run aborted after 3 s** (82 MB of raw
  instead of ~3 GB) — the scan itself then ran through untouched. So this run
  has **no nominal (450 V) anchor point of its own**: use `nominal_00/01` from
  `beam_nominal_meshscan_1` as the anchor, and the shared 440/430/420/410/400/
  390 points between the two scans as the coarse↔fine cross-normalization.
  Worth a quick look at why it died before the next run that starts with a
  nominal sub-run.
- QA: nothing processed yet — the qa_watcher is still draining the drift-scan
  backlog (~1 sub-run/h), so fine-scan waveform QA + telescope QA will follow
  overnight. Combined with the coarse scan this gives the mesh curve at
  5 V granularity over 390–450 V once QA is through.

## Run 6 — `p2in_check_1` (1 × 12 min): P2_IN responds at lower HV

Follow-up on the standing open item (P2_IN reading out but ~0.5 % trigger
share at nominal HV). Configuration: **P2_IN + the two uRWELL references only**
(FEUs 1/3; MID and OUT off), single sub-run `p2in_check_m400_d600` at **mesh
400 / drift 600** — both 90–100 V below the 490/700 nominal, keeping the drift
gap ≈ 200 V. 3.34 M triggers at ~4.6 kHz; HV rock stable (mesh imon 0.002 µA,
drift 0.012 µA); run finished normally.

**P2_IN is alive.** Quick look at the on-the-fly hit trees (raw ZS hits, no QA
selection — standard QA will follow once the watcher drains):

| raw-hits metric, FEU 3 (P2_IN) | nominal_00 (mesh 490) | p2in_check (mesh 400) |
|---|---|---|
| events with ≥1 hit / triggers | 0.79 M / 5.53 M = **14 %** | 2.47 M / 3.34 M = **74 %** |
| median hit amplitude | 61 ADC | 55 ADC |
| median hit significance | — | 11.7 σ |

**~5× more events with a P2_IN hit at 90 V *lower* mesh**, even though the
mesh-scan slope says 90 V should cost a factor ~8 in gain. That inverts the
gain logic and points squarely at the HV working point: at 490 V the mesh
channel showed imon spikes (07-23 entry) — the chamber was most likely
discharging/unstable above its max operating voltage, killing the response; at
400 V it is quiet and efficient. Consistent with P2_OUT's bench maximum being
420 V — 490 V was probably simply too high for this chamber.

**Next: scan P2_IN mesh upward from 400 V to find its real operating maximum
(watch imon for the onset of spikes), then bring it back into the telescope at
that point.** Note the 5 σ ZS threshold eats into efficiency at 55 ADC median
amplitude, so the operating point may want to be somewhat above 400 anyway.

## Analysis

Same chain as 07-23: `run_beam_qa.sh <run_name>` on banco →
`TB_July2026_H4/analysis/telescope/<run>/...` (per sub-run rates / occupancy /
signal / HV plots + JSON/CSV, run-level trend panel
`trend_beam_nominal_meshscan_1.png`). Raw + decoded data under
`TB_July2026_H4/runs/<run>/`.
