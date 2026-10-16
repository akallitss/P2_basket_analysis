# Results: P2 gain–timing map and low-field timing

Answer to `HANDOFF_gap_timing_2026-10-15.md`. Status 2026-10-15: **A1 and A2 done, B not started.**

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

**Drift foil (Alexandra):** the aluminium side **faces the gas**. A bare mylar surface charging up on the gas side is therefore not the mechanism. Charging would need another dielectric, for example at the frame edges.

Hysteresis and time-dependence analysis: not started yet.
- SPS products will go to EOS `p2_sps_july/analysis/low_field_timing/`, with a `CODE_VERSION.txt` and an `inputs.txt`.
- The cosmic part will go to the det1 drift-scan Analysis tree.
