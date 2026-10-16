# Results: P2 gain–timing map and low-field timing

Answer to `HANDOFF_gap_timing_2026-10-15.md`. Status 2026-10-16: **A1, A2 and B done.**

## Where the code and the products are

- **Code:** `cosmic_bench_analysis/24_gain_timing.py`. This is stage 24 of the cosmic pipeline: run key in, products to the run's Analysis tree.
  - Command: `python3 24_gain_timing.py det1_long5`
- **Products:** `data/Cosmic_Bench/Analysis/det1/p2_det1_long_run_efficiency_7-19-26/long_run_det1_415_615/24_gain_timing/`. Every file has the suffix `_without_connectors_1_2_10_spark_vetoed`.
  - `summary*.txt`: all numbers
  - `slopes*.csv`
  - `pads_fe55_161639*.csv`: the per-pad table
  - `offsets_vs_fe55_gain*.png`, `maps_gain_vs_rise*.png`, `rise_residual_vs_channel*.png`
  - `inputs*.txt`: git commit, input files with their dates, and the cuts
- **Inputs it needs:** stage 13 and stage 19 for the same run key. Stage 13 was rerun on 2026-10-15 for `det1_long5`, the full run, using all 5 decoded chunks. Stage 19 for `det1_long5` already existed.

## Inputs, and corrections to the handoff

- **Cosmic timing.** det1 long run 7-19: drift 615 / mesh 415 V, **Ar/iso 95/5**, DREAM readout.
  - The handoff used the stage 13 table built from 3 of the 15 chunks: 000–002, about 3.8 h.
  - Only 5 chunks are decoded on the LaCie: 000–002, 008 and 009.
  - Stage 13 now uses **all 5 decoded chunks**: 000–002, 008 and 009, about 6.3 h.
  - The apparent FEU6 dropout at about 11.5 h was an M3 reconstruction artefact, not a detector fault. The raw waveforms confirm it: in chunk 009, FEU6 connectors 3/4/5 record the same signal rate as in the other chunks, before and after event 243,561. So the run is used whole (`det1_long5`), with no pre/post split.
  - Result: 54,133 hits pass the §1.2 cuts, on 687 pads with ≥ 40 hits.
  - Rise and frac30 need waveforms, so their coverage stays central. The peak-time cross-check from stage 19 uses the combined hits of all 15 chunks.
- **Fe55 gain map.** det1, drift 750 / mesh 420 V, VMM readout, taken 2026-09-09 after the July→September remount.
  - **Use run 161639** (reduced VMM gain), not `central_150306` or `150423`.
  - In `150423` (original VMM gain), 296 of the high-gain pads are **clipped by saturation**. That removes the range in gain needed for the fit, and the correlation collapses to R² ≤ 0.16. This is a selection artefact, not physics.
  - Run 163851 (720/390 V; its low-gain pads are censored by threshold) is the cross-check.
- **Independence.** Timing comes from DREAM, gain from the VMM: two separate electronics chains, two months apart, with a remount in between.
- **The gain map is a stable detector property.** The cosmic ln(amplitude) map from July follows ln G(Fe55, Sept) with R² = 0.69 within connectors.

## A1: timing offsets vs ln G(Fe55)

Model, fitted per pad: offset = a·ln G + c_connector. Errors are bootstraps over pads. R² is the variance explained within connectors. Matched pads: 488.

| quantity | slope (ns per e-fold of G) | R² | shuffle null (95 %) | within-pad slope vs ln A |
|---|---|---|---|---|
| rise = peak − frac30 | **−9.9 ± 0.3** | 0.77 | 0.009 | −1.3 |
| peak (t_parabola) | **−20.8 ± 0.7** | 0.69 | 0.010 | −1.7 |
| frac30 | **−10.9 ± 0.5** | 0.45 | 0.008 | −0.4 |

**Ratio of the frac30 slope to the peak slope: 0.53 ± 0.02.**

These agree with the 3-chunk numbers (−9.1 / −18.4 / −9.5).

Robustness, all consistent:

| variant | rise | peak | frac30 |
|---|---|---|---|
| Fe55 run 163851 | −9.0 | −18.4 | −9.3 |
| amp ≥ 150 ADC | −10.3 | −21.3 | −11.0 |
| equal-charge cut (amp ≥ 300·G_rel) | −10.2 | −21.8 | −11.3 |
| cosmic ln A as the regressor, no Fe55 | −12.2 | −23.8 | −11.6 |
| amp ≥ 500 ADC (high-gain pads only, shorter range in ln G) | −8.6 | −18.7 | −10.7 |
| stage 19 peak time, all 15 chunks, 632 matched pads | — | **−19.4**, R² 0.72 | — |

**Time walk and amplitude selection.** Within a pad the gaps are fixed. There the slope against amplitude is 7–30× smaller than the slope between pads, and the equal-charge cut leaves the slopes unchanged. Walk or a Landau-selection bias cannot produce the map.

**Electronics control.**
- Connector medians of rise span 157.0–161.8 ns. These are absorbed by the connector constants.
- After the gain term, the rise residual std is 1.5 ns, against 3.0 ns before the fit.
- A small residual trend of about ±0.9 ns remains with channel order inside a connector: blocks of 8 channels hold 14 % of the residual variance. Channel order follows the strip geometry, so this cannot be separated from position.
- No DREAM pulser or test-pulse run exists for det1.
- The gain is measured with different electronics (VMM).

**Spatial structure.**
- The ln G map is smooth (neighbour correlation 0.90).
- After subtracting a·ln G, the frac30 residual has no spatial structure (+0.06, against a null of 0.10).
- **The rise and peak residuals have weak but significant spatial structure: rise +0.20 and peak +0.14, against nulls of 0.10 and 0.09.** Both were at the null with 3 chunks, so this appeared only with more statistics.
- Its size is small: the residual std is 1.5 ns, against 3.0 ns of rise variation.
- Possible causes:
  - Fe55 and cosmics differ in HV (750/420 vs 615/415 V) and in date (Sept, after the remount, vs July).
  - Or there is a second, gain-independent contribution to the pulse shape.
- I have not localised it yet.

**Verdict.**
- **Detector origin, tied to the gain.** DREAM electronics are excluded as the main cause; at most about ±0.9 ns of channel-order structure is left.
- Higher gain means an earlier, faster pulse. The gain explains 77 % of the rise-time map.
- **Amplification gap vs drift gap: unclear.**
  - For: a single cause acting through the gain fits the amplification gap.
  - Against: frac30 moves **half as much as the peak** (ratio 0.53 ± 0.02), and rise and frac30 move by the same amount. The toy model's amplification-gap signature ("frac30 hardly at all") does not reproduce this.
  - The weak spatial structure in the rise residual leaves room for a small second contribution.
- **Request to the simulation session, at the cosmic-bench conditions.** These data are from the cosmic bench, not the beam, so the simulation must match:
  - gas Ar/iso **95/5**, not the beam Ar/CO2/iso
  - drift 615 / mesh 415 V for timing, 750 / 420 V for the Fe55 gain
  - 4 mm drift gap, DREAM shaping, 60 ns sampling

  At those conditions, please give:
  1. the ratio of the frac30 shift to the peak shift for an amplification-gap change. This is the discriminating number. Ion mobility, and so the ion-tail shape, differs between the gases.
  2. whether a peak shift of −21 ns per e-fold and a frac30 shift of −11 ns per e-fold come out at all.

## A2: implied amplification-gap variation

Fe55 run 161639, all 728 valid pads:

| measure | value |
|---|---|
| std of ln G | 0.269 |
| std/mean of G | 0.30 (the handoff quotes 0.32) |
| 5–95 % range of ln G | 0.88 |

**Provisional conversion.** I use the Ar/CO2/iso scaling from §3 of the handoff: ×1.5 per 15 µm, so d ln G/dg = 0.027 per µm. That is the **wrong gas** for this detector, which ran in Ar/iso 95/5 on the bench. With that scaling:

- σ_gap ≈ 10 µm
- 5–95 % range ≈ 33 µm

**Please redo the conversion with the simulation's Ar/iso 95/5 gain-vs-gap table at 420 V.**

These numbers are **upper limits on the gap part**. The Fe55 spread also contains VMM per-channel gain dispersion and any gain source other than the gap (mesh tension, local charging). Fe55 statistics are negligible at ~2,000 hits per pad.

For the timing side, the simulation needs to reproduce:
- a peak shift of −21 ns per e-fold, which is ±5.6 ns over the 1σ gain spread
- a frac30 shift of −11 ns per e-fold

## B: low field

**Drift foil (Alexandra):** the aluminium side **faces the gas**. A bare mylar surface charging up on the gas side is therefore not the mechanism.

### SPS (done 2026-10-16)

- **Code:** `sps_beam_analysis/32_low_field_timing.py`, commit `f85a439` and later. It has two modes:
  - `extract` runs on lxplus straight off EOS and keeps the leading hit per station for every trigger.
  - `analyse` runs locally.
- **Products:** EOS `p2_sps_july/analysis/low_field_timing/`, with `inputs.txt` and `CODE_VERSION.txt`. The per-trigger tables are in `extract/` (163 `.npz` files, with their own `CODE_VERSION.txt`).
  - `low_field_timing_vs_field.png` and `low_field_timing_within_subrun.png`
  - `low_field_timing_subruns.csv`, `low_field_timing_slices.csv`, `low_field_timing_trends.csv`
  - `low_field_timing_summary.txt`
- **Sample:** 8 runs, 166 sub-runs, all in Ar/CO2/iso 93/5/2.
- **Selection:**
  - leading hit per station with amplitude ≥ 60 ADC, not saturated or truncated, inside the 960 ns window;
  - a **P2_MID–P2_OUT time coincidence**, |Δt − median| < 150 ns. Accidentals are negligible: sidebands hold ≤ 10⁻³ of the peak, even at 60 ADC.
  - The width is the 68 % half-width of (MID − OUT)/√2. The time is the median of MID.
  - The 200 ADC cut of the simulation session was not used. At low field the median amplitude falls to about 150 ADC, so that cut keeps only the high tail.
- **Excluded:**
  - `mesh_drift_scan_up_1/drift2_11`: the decode hung, so P2_IN and P2_OUT are missing.
  - `drift_mesh_scan_1/drift_700`: coincidences only in the last sixth of the sub-run.

1. **No time dependence inside sub-runs.**
   - Each sub-run is 7–17 min, split into 6 slices, the first of which covers about the first minute after the HV step.
   - The width, the median time and the yield are flat. Over the 27 sub-runs at ≤ 400 V/cm, χ²/ndf against "no change" is **0.60 for the width and 0.94 for the median**.
   - The weighted mean change from the first to the last slice is −0.06 ns for the width and −0.23 ns for the median.
   - At 375–500 V/cm the sub-runs are flat to better than 1 ns and 2 %.
   - This rules out a charging-up with a time constant between about 1 and 10 minutes.
2. **No dependence on scan direction.**
   - **`drift_mesh_2d_2`** alternates direction row by row (mesh 330–450 V, 375–1125 V/cm). Each down row differs from the mean of its two up neighbours by **at most 2 ns, and mostly ≤ 0.8 ns**.
   - **The downward scan** (`mesh_drift_scan_up_1`, 07-26, 1062 → 25 V/cm) matches the upward scans of 07-27 and 07-28 at mesh 450 V within about 1 ns: for example 36.8 / 33.8 / 33.4 ns at 375 V/cm and 19.0 / 18.5 / 18.1 ns at 750 V/cm.
   - **At the one point with identical HV** (125 V/cm, mesh 450), up (07-25) and down (07-26) agree: width 43.2 ± 0.5 vs 45.0 ± 0.7 ns, median 278.0 vs 277.5 ns, relative yield 0.155 vs 0.153.
3. **A day-scale change, not a direction effect.**
   - 07-25 is better than every later day, but only at intermediate fields. At 375 V/cm the width is 28.4 ns against 33.4–36.8 ns, and the MID amplitude 293 ADC against 234–264 ADC.
   - At ≥ 875 V/cm all days agree: width 17.1–17.8 ns, MID amplitude 372–391 ADC. So it is not a change of gain.
   - P2_OUT shows the same extra low-field loss on the later days.
   - The later curves look like the 07-25 curve shifted by about +100 V/cm, which is about 40 V on the drift.
   - **This points to a slow change common to both chambers**, for example gas purity (electron attachment grows with drift time) or a small effective drift-field offset. It does not point to fast charging. I have not tested which; 07-24 cannot help, because P2_MID ran at low gain that day (amplitude about 97 ADC).
4. **What the low-field signal looks like** (mesh 450 V, 07-26 down scan):
   - The single-chamber relative yield, √(coincidences per trigger over the run's plateau), falls steeply below about 500 V/cm:

     | field (V/cm) | 1000 | 750 | 625 | 500 | 375 | 250 | 125 |
     |---|---|---|---|---|---|---|---|
     | 07-26 | 1.00 | 0.93 | 0.82 | 0.63 | 0.36 | 0.15 | 0.15 |

     07-25 gives 0.98, 0.93, 0.82, 0.57, 0.22 and 0.16 at 750–125 V/cm.
   - The median MID time is only **+11 ns** later than the plateau at 125 V/cm, +34 ns at about 250 V/cm, and **−35 ns at 0 V/cm**.
   - At zero drift field about 1,000 coincidences remain. They are narrow (17 ns) and early, so they must be charges produced at the mesh.
   - As the field falls, the surviving signal comes more and more from near the mesh. The yield loss is the electrons from the upper part of the gap that never arrive, as the handoff suspected.
   - Strong electron attachment, or a field strongly reduced over most of the gap, would both produce this. The flat time behaviour and the lack of a direction effect argue against a dynamic screening charge.
   - The relative yield is a proxy, not an efficiency. It assumes the muon fraction of the triggers is stable within a run, and the 60 ADC floor plus zero suppression bias it at the lowest fields.
5. **A P2_MID-only anomaly at 200–250 V/cm** on both days. Its amplitude halves (137–172 ADC against 260 at 162 V/cm), and its median arrives about 25 ns late, while P2_OUT stays smooth.
   - On 07-26 it spans `drift2_07`–`10`. The logbook reports sparklets on the MID mesh, worst in `drift2_08`.
   - On 07-25 it is the single 250 V/cm point.
   - Not understood. I flag it rather than explain it.
- **Caveat:** the `run_time.txt` of `drift_mesh_scan_1/drift_450` (0 V/cm) records a 3 s run, about 9 h after the rest. The hits in its file span the full sub-run (2.74 M triggers), so `run_time.txt` was probably overwritten by an aborted retake. Its timing within the day is uncertain.

**Verdict B (SPS):** the low-field degradation is **static on the scale of minutes and independent of scan direction**, so a charging-up of a dielectric surface during the scan is not supported. On top of that, there is a **slow, day-scale worsening at intermediate fields common to both chambers**. Suggested next test, simulation side: the attachment coefficient needed to reproduce the yield-vs-field curve above, and the matching drift-time cut-off.

### Cosmic (done 2026-10-16)

- **Code:** `cosmic_bench_analysis/25_drift_scan_time_slices.py`, commit `63f61b1`, run as `python3 25_drift_scan_time_slices.py det1_driftscan2`.
- **Products:** `data/Cosmic_Bench/Analysis/det1/p2_det1_drift_scan_7-19-26/drift_scan/25_drift_scan_time_slices/`, suffix `_without_connectors_1_2_10_spark_vetoed`:
  - `time_slices*.png`
  - `time_slices*.csv` and `time_trends*.csv`
  - `summary*.txt` and `inputs*.txt`
- **Setup:** det1, Ar/iso 95/5, mesh 415 V, drift 415 → 965 V (0–1375 V/cm). Each sub-run is 30 min with about 10k triggers, cut into 6 slices of 5 min.
- **Estimators**, the same as stage 16's `time_res_ns`: the leading pad (≥ 300 ADC), its `time_of_max`, and σ = (p84 − p16)/2. Hot pads are dropped and the spark veto is applied. The yield proxy is signal events per trigger.

1. **Flat inside every sub-run, at low field too:**

   | field (V/cm) | σ in the 6 slices (ns) | median in the 6 slices (ns) | yield proxy |
   |---|---|---|---|
   | 125 | 78.9, 81.0, 78.0, 81.5, 87.4, 84.9 (± 3–4 each) | 566–577 (± 4–6 each) | 0.32–0.36 |
   | 250 | 26.7–29.5 (± 1) | 460–467 | 0.36–0.39 |
   | 375 | 20.7–22.3 (± 0.8) | 446–450 | 0.37–0.41 |

   - The first-to-last changes are all below 2σ and point in opposite directions at different fields: 125 V/cm gives +6 ± 5 ns width and +11 ± 6 ns median; 375 V/cm gives −2.6 ± 1.3 ns median.
   - Over all 11 sub-runs, χ²/ndf against "no change" is 0.97 for σ and 1.54 for the median.
2. **The drift HV is at its setpoint within 0.5 V from the first slice** at every field. The only exception is the 0 V/cm point, where it ramps up from about 398 V during the first 5 min, so the scan ran upward from 0 V/cm.
3. **Cosmic and SPS behave differently at 125 V/cm.**
   - Cosmic (Ar/iso): the signal is **late (median +120 ns relative to 375 V/cm) and broad (σ ≈ 82 ns), and still present**. The yield proxy is 0.34, against 0.39 at 250–375 V/cm.
   - SPS (Ar/CO2/iso, MID–OUT coincidences): the signal at 125 V/cm is nearly gone (relative yield about 0.15) and only about 11 ns late.
   - The two yield measures differ (300 ADC on a single detector against a 60 ADC coincidence), so compare the trends, not the numbers.
   - At face value, in Ar/iso the electrons from the whole gap still arrive, slowly, while in the beam gas those from the upper gap are lost. That fits an attachment or loss mechanism specific to the beam gas or setup, rather than one common cause. It is not proven: the beam rate, the electronics and the selection all differ.

**Verdict B, overall:**
- In both setups the low-field timing is **constant over minutes to half an hour after an HV step**.
- At the SPS it is **independent of scan direction**.
- The drift HV is at its setpoint.
- So a dynamic charging-up of a surface is not supported by the data.
- What remains is static:
  - in Ar/iso (cosmic), late and broad signals that are still collected;
  - in the beam gas (SPS), a strong loss of the electrons from the upper gap, plus a day-scale worsening at intermediate fields.
- Suggested next step, simulation side: compare the cosmic 125 V/cm median shift (+120 ns) and width (82 ns) with Magboltz v_d and diffusion for Ar/iso 95/5. For the SPS, fit an attachment rate to the yield-vs-field curve.
