# 2026-07-23 — P2 telescope, first beam at SPS H4 (TB_July2026_H4)

Runs `beam_commissioning_1` and `latency_scan_1`. First beam on the 5-detector
telescope. **Both runs are good; the DAQ chain is validated end to end and the
running configuration needs no change.**

Setup: 5 detectors on 4 FEUs — EIC_uRWELL_front (z=0, FEU 1 ch 0-255), P2_IN
(z=320 mm, FEU 3), P2_MID (z=630 mm, FEU 4), P2_OUT (z=940 mm, FEU 5),
EIC_uRWELL_back (z=1370 mm, FEU 1 ch 256-511). External scintillator
coincidence into the TCM, zero suppression, 16 samples × 60 ns, Dream latency
32. Gas Ar/iso 95/5. HV: P2_IN 700/490, P2_MID 700/450, P2_OUT 700/450
(drift/mesh), uRWELL drift 600, resistive 420.

Beam: SPS slow extraction, **43 s spill period, ~5 s spills, ~11 kHz
instantaneous trigger rate in spill**, ~1200 Hz averaged over a run.

---

## Run 1 — `beam_commissioning_1` (2 sub-runs × ~160 s)

**Trigger distribution is correct.** Event-sync QA: **ALIGNED**, best shift 0
and median Δt = 0 µs on all four FEUs, both sub-runs. The TCM fans the external
trigger out correctly.

**P2_MID and P2_OUT are healthy.** 0–3 dead channels out of 512 each, clean
Landau spectra, full beam spot visible across connectors c_5/c_6.

| sub-run | triggers | empty | P2_OUT ev / share / amp | P2_MID ev / share / amp | P2_IN ev / amp |
|---|---|---|---|---|---|
| 00 | 191 458 | 10 % | 145 177 / 0.76 / 697 ADC | 103 664 / 0.54 / 510 ADC | 934 / 643 ADC |
| 01 | 192 346 | 40 % | 97 399 / 0.51 / 694 ADC | 69 351 / 0.36 / 511 ADC | 908 / 751 ADC |

- **Sub-run 00 caught the HV ramp** (200 V → nominal over its first ~55 s, only
  74 % of the monitored time at setpoint), so its first spill is at low gain.
  **Sub-run 01 is fully at nominal — use it for anything quantitative.**
- HV otherwise stable, imon < 1 µA, no sparks.
- **P2_IN reads out but is ~20× below its neighbours** (share 0.00). Its mesh
  (8:1, 490 V) does draw current with imon spikes, so it is biased — whatever
  limits it is downstream of the HV. **Main open item.**
- uRWELL references see far less (share 0.02 front / 0.08 back); their ~20 % of
  silent channels are most likely unconnected strips, not a fault.

## Run 2 — `latency_scan_1` (Dream latency 24 / 28 / 32 / 36 / 40, ~160 s each)

Event sync **ALIGNED at every point**. The peak sample moves exactly **1 sample
per latency unit and wraps modulo the 16-sample window**.

| latency | 24 | 28 | **32** | 36 | 40 |
|---|---|---|---|---|---|
| peak sample (of 16) | 15.6 ↩ | 3.4 | **7.1** | 11.1 | 15.1 |
| P2_OUT median amp [ADC] | 623 | 660 | **700** | 680 | 676 |
| P2_MID median amp [ADC] | 377 | 421 | **511** | 470 | 422 |
| P2_IN median amp [ADC] | 614 | 422 | **644** | 452 | 252 |

### → CONCLUSION: keep latency 32. No change, and a fine scan is not worth doing.

All three P2 stations peak in median amplitude at 32, and the peak sample lands
at 7.1 of 16 — about 7 samples of baseline before the pulse and 9 of tail after.
The neighbours at ±4 are already 8–10 % down, so the optimum is within ~±1 of 32
and fine-tuning would buy under 2 %. Latency 40 and 24 clip at the window edge;
28 cuts the leading edge.

**Method note — judge a latency scan on AMPLITUDE, not on trigger share.**
Three sub-runs sat at latency 32 with identical configuration
(`beam_commissioning_00/01`, `latency_032`):

- trigger share: 0.76 / 0.51 / 0.75 → **± 25 points**, dominated by beam intensity
- median amplitude: 697 / 694 / 700 → **± 0.4 %**

Share appears to peak at 36 (0.81 vs 0.75), but that difference is far inside
the run-to-run scatter. It is beam, not latency. Amplitude is the discriminant.

---

## Other results worth recording

- **Charge sharing is intrinsic, not a threshold artifact.** Only ~8 % of events
  have ≥2 pads (mean 1.15 pads/event), flat across all five latency points. Hits
  survive down to ~30 ADC (5th percentile 213), so the ZS threshold is not
  cutting neighbours off and lowering it (`PedRun Threshold` 5.00 σ) would not
  obviously buy cluster size. Position resolution here is pad-pitch limited.
- **Saturation** ~3.8 % of P2_OUT hits at latency 32 (8.4 % at 24). Acceptable;
  trimming mesh HV would reduce it at the cost of amplitude — that trade needs
  an HV scan, the natural next run.
- **Watch the mesh current**: imon reached 6.3 µA during `latency_032` and
  5.0 µA during `latency_028`, against < 1 µA in the commissioning run.
- `latency_028` caught one spill fewer than the others (2 spills, 130 543
  triggers, 808 Hz vs ~1200 Hz elsewhere) — lower statistics, not a fault.
- `latency_024` has 97 % empty triggers: with the window off the pulse, almost
  no trigger records a hit anywhere. Expected, and a good sanity check that the
  empty-trigger metric works.

## Operational note

`decode` **hung** on `beam_commissioning_01` FEU 03 — 36 min at 100 % CPU with
its output already written — blocking the whole processor pipeline, so the
latency scan did not process at all until it was killed. Re-decoding the same
file afterwards took **2.5 s**. Killing it leaves a truncated `.root` with no
keys, which reads back as "that detector has 0 events"; the sub-run must be
re-run through decode → analyze_waveforms → combine_feus_hits by hand. **If the
pipeline stalls, check `ps` for a long-running `decode`.**

## Analysis

`sps_beam_analysis/25_commissioning_qa.py` (+ `24_event_sync_qa.py`), run on
banco via `run_beam_qa.sh <run_name>`. Products under
`TB_July2026_H4/analysis/telescope/<run>/...` — per sub-run rates / occupancy /
signal / HV plots plus JSON and CSV, and a run-level trend panel which for a
scan run is the scan curve (`trend_latency_scan_1.png`).
