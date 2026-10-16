# Handoff: P2 gain–timing map and low-field timing (experimental analysis)

Written 2026-10-15 by the simulation session (P2_Geant_Garfield), so the data
analysis can run in its own session. The simulation session carries on
separately. Please write results back to the place given in §6.

## 0. Why this matters

The simulation goal: find ways to raise the detector gain to make up for the VMM
threshold loss, while keeping time resolution ≲ 20 ns and e/γ separation. To
trust the simulation we validate each stage against data. Two features of the
data are not understood yet:

- **A.** Across the det1 pad plane, **high-gain regions arrive earlier**. This
  shows up in both the cosmic MIP timing map and the Fe55 gain map.
- **B.** **Below ~375 V/cm drift field the timing degrades far more than
  Magboltz predicts.** This happens in both gases and both setups (SPS
  Ar/CO2/iso 93/5/2, cosmic Ar/iso 95/5). At ≥ 500 V/cm the simulation agrees
  with the data within 1–2 ns.

Rule from Alexandra: **do not assume the bent-PCB explanation, prove or disprove it
with data.** Don't just agree with a hypothesis. Report what the data say,
including null results.

## 1. What the simulation session already found (please check, don't trust blindly)

Data used: `/local/home/ak271430/Documents/PostDocSaclay/data/Cosmic_Bench/Analysis/det1/p2_det1_long_run_efficiency_7-19-26/long_run_det1_415_615/`
(stages 13, 19, 20, 23; drift 615 / mesh 415 V, Ar/iso 95/5, 18 h).

1. **Stage 23's conclusion "time walk on a gain gradient" is not supported.**
   - From `13_timing/toa_hits_without_connectors_1_10_spark_vetoed.csv`: amp ≥ 300 ADC, not saturated, |t_parabola − median| < 150 ns. That leaves 31,956 hits on 358 pads with ≥ 40 hits.
   - The within-pad slope of `t_parabola` vs ln(amplitude) is only **−1.8 ns per e-fold**.
   - Walk therefore predicts a pad-to-pad std of 0.4 ns, but **7.8 ns** is observed.
   - Stage 23's r = −0.83 comes from only 11 points binned by distance from the apex.
2. **The pulse shape changes with gain.** Define rise = `t_parabola − t_frac30`, the per-pad median with connector medians subtracted.
   - Correlation with the per-pad ln(amplitude): **r = −0.77**. The map is reproducible: split by eventId parity, **r = 0.90**.
   - Per-pad amplitude quintiles, from low to high gain:

     | quantity | low gain → high gain |
     |---|---|
     | peak time | +7.4 → −6.5 ns |
     | frac30 time | +3.1 → −3.1 ns |
     | rise | +3.8 → −3.1 ns |

3. Interpretation from a toy model (`P2_Geant_Garfield/scripts/analysis/sps_drift_velocity_check.py`):
   - An **amplification-gap** change (ion transit ∝ g² at fixed V) moves the **peak** strongly and frac30 hardly at all.
   - A **drift-gap** change moves both. ±0.5 mm gives only ±5 ns.
   - **Leading hypothesis:** µm-scale amplification-gap non-uniformity (mesh/PCB non-planarity acting through ion collection), plus possibly a ~±0.3 mm drift-gap part.
   - The toy overpredicts the size by ~3×, because its ion model is crude. **Not proven.**
4. **Stage 16/23's "apparent v_d" is not a drift velocity.**
   - It is `vd_apparent = 0.8·gap / (p90−p10 of the peak time)`. It gives 57 µm/ns at 500 V/cm, against 35 from Magboltz Ar/iso.
   - The peak time does not span the full drift time, so the stage 23 v_d-scaling test (`23_drift_gap_bow/bow_vd_scaling_*.csv`) cannot reject or support a gap variation.
5. **Low field (B).**
   - Cosmic drift scan (`../p2_det1_drift_scan_7-19-26/drift_scan/16_drift_scan_efficiency/efficiency_vs_drift_without_connectors_1_2_10_spark_vetoed.csv`, column `time_res_ns`) vs the toy:

     | field (V/cm) | 125 | 250 | 375 | 500 | 750 | 1000 | 1250 |
     |---|---|---|---|---|---|---|---|
     | data (ns) | 82 | 28 | 21.7 | 20.9 | 23.1 | 25.5 | 26.9 |
     | sim (ns) | 21 | 19 | 19.7 | 20.7 | 23.1 | 24.4 | 25.7 |

   - SPS `drift_mesh_scan_1`, single-chamber width = 68 % half-width of (MID − OUT)/√2:

     | field (V/cm) | 125 | 250 | 375 | 500 | ≥ 625 |
     |---|---|---|---|---|---|
     | data (ns) | 45 | 35 | 25 | 20 | 16–17 |
     | sim (ns) | 48 | 20 | 15 | 14 | 14–16 |

   - At 125 V/cm the SPS median time is ~55 ns EARLIER than simulated, as if electrons from far up the gap were missing.
   - SPS HV monitor: vmon = setpoint and currents ~0, so a supply offset is ruled out.
   - In Ar/iso, Magboltz v_d is near its maximum at 125 V/cm, so a different gas composition can hardly explain a ~4× slowdown. Either the real field in the gap is much lower than (V_drift − V_mesh)/gap at low settings, or there is a loss mechanism.
   - **Open question for Alexandra:** does the Al of the 50 µm drift foil face the gas or the outside? A bare mylar surface facing the gas could charge up with ions and screen the low drift fields.

## 2. Test A1: rise-time map vs Fe55 gain map (the main test of the amplification-gap hypothesis)

- Fe55 gain map, same detector (det1, drift 750 / mesh 420 V):
  `/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/LAB_Measurements/source_tests/Fe55_data/position_scan/central_position/analysis/`
  (figure `central_150306_gain_map.png`; find the per-pad table behind it).
- Cosmic per-pad rise, peak and frac30 offsets: recompute from `toa_hits_*.csv` as in §1.2, or use the 1_2_10 connector set if a per-hit table exists for it.
- Do, pad by pad:
  1. Correlate rise, peak and frac30 offsets with **ln(Fe55 relative gain)**. Fe55 is an independent gain measure: point-like, no Landau or MIP threshold effects.
  2. Fit offset = a·ln(G_Fe55) + connector constants. Report the slope (ns per e-fold) for rise, peak and frac30 separately, and the variance explained.
  3. Mask out the Fe55 low-statistics and edge pads (grey or very few hits).
- Expectation if the hypothesis is right: the **rise** and **peak** offsets follow ln G strongly, and frac30 weakly.
  - If **frac30** follows ln G as strongly as the peak, a drift-gap / drift-field component is significant.
  - If none follow, look elsewhere (electronics: per-channel shaping, FEU or connector).
- Electronics control: check whether the rise-time map has connector/FEU structure (per-connector medians, a 64-channel pattern). Pulse-shape differences between DREAM channels must be excluded. If a pulser/test-pulse run exists, the per-channel rise there is the clean reference.

## 3. Test A2: size check

Convert the measured slope into an amplification-gap change. For Ar/iso at 415 V the simulated gain at 150 µm is 5.4e3. The simulated Ar/CO2/iso values are 1.2e4 / 1.8e4 / 2.7e4 / 4.1e4 at 165 / 150 / 135 / 120 µm, i.e. ~×1.5 per 15 µm. The simulation session can supply an Ar/iso gain-vs-gap table on request. Report the µm range of gap variation implied by the Fe55 gain spread (std/mean 0.32), so the simulation can include it.

## 4. Test B: low-field effect (cause of the slow/lossy behaviour below ~375 V/cm)

1. **Hysteresis / charging.** SPS has scans in both directions. Compare `drift_mesh_scan_1` (drift 450→900) with `mesh_drift_scan_up_1` and any downward drift scan at the same field, for timing width, median and efficiency. Also look at the time dependence WITHIN a low-field sub-run (first vs last minutes). Charging of a dielectric cathode surface predicts direction dependence and a time constant.
   - Products: `/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/mpgd26_workspace/products/analysis/`
   - Per-hit data on EOS: `/eos/experiment/ntof/data/x17/p2_sps_july/runs/<run>/<sub_run>/combined_hits_root/` (tree `hits`; FEU 3/4/5 = P2_IN/MID/OUT; `time` = 30 % leading edge, ftst-corrected).
   - Extraction script to reuse: `P2_Geant_Garfield/scripts/analysis/sps_hit_times_extract.py` (largest hit per event, 200–3500 ADC, not saturated or truncated).
2. **Cosmic.** Same within-run time dependence for the det1 drift-scan sub-runs at 465 / 515 V.
3. **Pitfall:** the `28_timing_qa` `mean_ns` has no amplitude cut, so noise hits pull it toward the window centre at low field (this faked a ~160 ns shift). Always use an amplitude cut and medians.
4. **Note:** SPS `drift_450` and `drift_500` were taken ~9 h after the rest of `drift_mesh_scan_1` (10:38 vs 01:25–03:00), so gas or charging conditions may differ.

## 5. Practicalities

- lxplus: plain non-interactive ssh fails (Kerberos gives "partial success" and a 2nd factor is required). Reuse an open master socket:
  `ssh -o BatchMode=yes -o ControlPath=$HOME/.ssh/master-%r@%h:%p akallits@lxplus.cern.ch`.
  If it is gone, Alexandra opens an interactive `ssh lxplus` first. `kinit` can't be typed through Claude's `!` prompt (no tty).
- Gas at SPS: take it from the corrected run_config only (Ar/CO2/iso 93/5/2 for runs ≤ 58). The cosmic bench is Ar/iso 95/5.
- Keep new figures minimal: white background, Arial, title + axes.

## 6. What to hand back to the simulation session

A short markdown at `P2_basket_analysis/RESULTS_gap_timing.md` with:
- A1: slopes (ns per e-fold of Fe55 gain) for rise / peak / frac30, variance explained, the electronics-control result, and a verdict: amplification gap / drift gap / electronics / unclear.
- A2: the implied amplification-gap variation in µm (range and std), or "not determinable".
- B: hysteresis / time-dependence result, and Alexandra's answer on the drift-foil Al orientation.
- Paths to the figures and tables.
