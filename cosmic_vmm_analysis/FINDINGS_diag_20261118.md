# Why do the P2 chambers lose ~15 % of the muons with the VMM read-out? Diagnostics of 18 Nov 2026

Lab cosmic stand (three P2 chambers IN / MID / OUT, VMM3a read-out, scintillator trigger on H0 ch 40), HV at the
operating point (IN 420/720 V, MID and OUT 430/730 V, mesh/drift), peaking time 200 ns. Runs 231-267 (16-18 Nov),
36 runs used (run_233 has a read-out fault and is left out), about 30 h of data, 79 574 on-time scintillator tags.
All processing on lxplus / HTCondor (never on the DAQ box).

## The question and what was already known

- The single-channel efficiency (counting estimator `channel_efficiency.py`, `CHANNEL_EFFICIENCY.md`) is about
  **0.85-0.90** with HV on, and it does **not** change with the threshold (17 Nov ladder, 14-20 mV above baseline,
  `DAQ_Control_VMM_Beam/docs/eff_ladder_20261117/README.md`) nor with the VMM gain (18 Nov gain scan, 4.5 / 6 / 9 mV/fC).
- The same chambers read out with Dream electronics during an HV scan reach **more than 90 %**, and the mesh is already
  at the highest voltage before sparking. So the missing efficiency is on the VMM side, not the detector gain.
- Question for this study: are the sparks (detector instabilities) contaminating the measurement, and where do the
  missed muons go? (user, 18 Nov: integrate the HV supply records + points 3, 4, 6 of the data list).

## What was done

| step | what | where |
|---|---|---|
| tag dump | `channel_efficiency.py dump --diag` (= `efficiency_tp.py --dump-diag`): for every tag its wall-clock time, every live probe hit within 15 mm of the prediction over +-2 us (time, ADC, distance, channel), chamber hits and saturated hits (ADC >= 1000) within +-100 us, time since the previous hit of the predicted channel and chip; per-second hit / saturated-hit rates per chamber. The fields of the normal dump are bit-identical (checked on chunk 26001). | 37 condor jobs, 1265 chunks |
| bad times | `spark_veto.py` on every run's own `cfg_*/hv_monitor.csv` (5 s sampling, DAQ-box clock = the clock of the capture files) and `tdk_lv_monitor.csv` | `veto/veto_run_N.json` |
| efficiency | `channel_efficiency.py measure` twice on the same tags: without and with `--bad-times`; regrouped per configuration with the pipeline's `regroup.py` | `results/all`, `results/veto` |
| diagnostics | `probe_diagnostics.py`: (1) HV flags vs chamber activity, (2) efficiency with/without veto, (3) ADC, (4) timing, (5) missed muons by cause | `results/diag/` |

Data: EOS `/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus/diag_20261118/`.
Copies of the summary tables and figures: `figures/diag_20261118/` (next to this file).
Every quantity is computed on on-time minus scaled off-time tags (fake tag pairs subtracted), scintillator window
140-310 ns, as in `channel_efficiency.py measure`.

## 1. Sparks and other HV / LV flags

Flags in `spark_veto.py`: `hv_spark` (current more than 1 uA above the run median), `hv_small` (0.15-1 uA),
`hv_vdev` (|vmon - v0| > 1.5 V), `hv_off`, `lv_glitch`; each flag vetoes [t - 6 s, t + 5 s].

- Per run: 2-16 real sparks, almost all on the **P2_MID mesh** (iseg NHR); many small current pulses on the
  **P2_IN and P2_OUT meshes** (CAEN N1471). The N1471 otherwise reads exactly 0.00 or 0.04 uA, so the 0.19-0.56 uA
  readings are real current pulses, not read-out steps. All flags together veto **9-25 % of the run time**.
- Do the flags come with chamber activity? Saturated hits per second within +-10 s of the flag, against random times
  (`diag_spark_flag_check.csv`, `.png`):

| flag class | flags | saturated hits / s, peak vs random | all hits / s, peak vs random |
|---|---|---|---|
| hv_spark IN mesh | 3 | x7.5 IN, x5 MID, x5.9 OUT | x1.5 IN, x2.3 OUT |
| hv_spark MID mesh | 220 | x2.4 MID only | x1.01 |
| hv_spark OUT mesh | 22 | x1.7-2.5 | x1.0-1.1 |
| hv_small IN mesh | 728 | x2.0 IN | x1.02 |
| hv_small OUT mesh | 473 | x2.0 OUT | x1.02 |
| hv_vdev MID / OUT mesh | 14 / 13 | x2-3.4 | ~x1 |
| lv_glitch hybrids | 162 | x1.1-2.3 | x1.03 |
| random times | 7200 | x1.03-1.06 | x1.00 |

  Only the three IN-mesh sparks produce a large burst. The others come with a small rise of saturated hits in the
  chamber of the electrode and no change of the hit rate.
- The data show many more discharge-like seconds than the HV log: 17 676 seconds whose saturated-hit count is above
  the run median + 10 MAD; only **22 %** of them fall inside the HV veto (about the time share of the veto itself),
  and only 35 % have an HV flag within 10 s. A 5 s sampled supply current misses most short discharges.
- **Effect on the efficiency: none.** Counting estimator, the same tags with and without the veto (the veto removes
  12-17 % of the muons): mean change **+0.001**, range -0.013 to +0.015, inside the errors (0.016-0.05) at every one of
  the 30 configuration x chamber points (`veto_by_config.csv` vs `all_by_config.csv`).

**Conclusion:** detector instabilities do not bias the efficiency (nor, therefore, what we conclude about the
thresholds). If they ever need removing, a data-side cut on bursts of saturated hits is the better tool; the HV
monitor at 5 s is too coarse.

## 2. Is the result limited by statistics?

- Chamber level, no: about 1000-3800 expected muons per configuration and chamber give +-0.015-0.03; the estimator
  is linear in the counts (no small-number bias), the toy closure gives 0.845-0.849 for 0.85, and repeated runs of
  the same file agree (four Km18 runs of 17 Nov: IN 0.88-0.92, OUT 0.92-0.93).
- What does bias numbers: (a) the "same pads in every configuration" table (`paired.csv`) shrinks to 417 pads once
  the 9 mV/fC files (many MID pads masked) are included, and MID then exceeds 1: use `by_config.csv` (all live pads);
  (b) chi2/ndf up to 7 (Km18 P2_IN): scale that error by sqrt(chi2/ndf); (c) the result depends on R_lo
  (+0.03 to +0.14 from 24 to 56 mm: tails of the track extrapolation); (d) per-chip values of 0.17-0.48 (chip 8) and
  1.0-1.4 (chips 3, 11, 12, 15) with errors of 0.05 are 4-8 sigma from 0.85 and above 1 is impossible: a systematic
  (channel-to-pad map, pad normalisation, or hits belonging to a neighbouring chip), not statistics.
- The no-veto pass of this study reproduces the pipeline exactly for the same runs (`all_by_config.csv` vs the
  pipeline's `results/current/by_config.csv`).

## 3. Pulse height of the hits that detect the muons

Largest ADC of the live hits within 15 mm and +-100 ns of the prediction, per detected muon (`diag_adc_spectra.png`,
quantiles in `diag_summary.txt` section 3, per chip in `diag_adc_per_chip.csv`):

| gain (mV/fC) | median ADC IN / MID / OUT | below 100 ADC | saturated (>= 1000) |
|---|---|---|---|
| 4.5 (all threshold files) | 209-222 / 166-177 / 141-148 | 11-26 % | 0.2-2.9 % |
| 6 | 279-285 / 230 / 182 | 6-14 % | 0.8-4.5 % |
| 9 | 384-395 / 313-331 / 256-272 | 0.2-2.8 % | 2.9-8.4 % |

- The spectrum scales with the gain (x1.8 from 4.5 to 9, expected x2).
- At 4.5 mV/fC the spectra have a sharp lower edge at 60-80 ADC: the threshold cuts into them (most clearly on P2_OUT,
  whose peak sits right at the edge). At 9 mV/fC the rising side of the spectrum is visible above the edge.
- Yet the efficiency at 9 mV/fC is not higher (0.79-0.85 vs 0.86-0.91 at 4.5, the Km18 reference). So the muons that
  are missed are not the small pulses below threshold. This is the argument the weekend 2D scan is designed to make
  quantitative (per-channel threshold in ADC from the noise hits, against the muon spectrum of the same channel).
- P2_OUT has the smallest pulses at every gain.

## 4. Timing of the hits

- The hits within 15 mm sit within about -75..+72 ns (5-95 %, worst files) of tag + latency; median shifts by up to -20 ns on MID in
  the bscan files (latency is measured per chunk).
- Efficiency vs the half-width of the time window (random hits subtracted, `diag_timing.png`): the +-100 ns window
  used for the efficiency is already on the plateau. Opening it to +-1.5 us adds **+1 to +3 points** at 4.5 mV/fC,
  up to **+8 points on P2_IN at 9 mV/fC** (late hits at high gain).
- Hits of the predicted pad's own channel within 10 us before the muon (dead time): ~0 % of the muons.

## 5. Where the missed muons go

Every muon (on - r off) without a live hit within 15 mm and +-100 ns, assigned to the FIRST matching cause
(`diag_summary.txt` section 5, `diag_misses_per_cfg_chip.csv`). This uses the simple nearest-hit efficiency
(0.32-0.82), which is lower than the counting estimator because it also counts tracking tails as misses; read the
classes, not the efficiency, from it.

| cause | fraction of all muons | comment |
|---|---|---|
| hv_bad_time (inside the HV/LV veto) | 2-8.5 % | the veto's time share x (1 - eff): no excess |
| discharge (a saturated hit in the probe chamber within +-100 us) | 0.6-9.4 % | rises with gain (9 mV/fC IN 8-9 %); random level about 1 % |
| out_of_time (hit within 15 mm but outside +-100 ns) | 0.1-1.1 % | |
| own_dead_time | ~0 | |
| next_to_masked (within 15 mm of a dead/masked pad) | 0.5-7.3 % | larger at 9 mV/fC on MID (many masks) |
| **chamber_silent (no live hit anywhere in the probe chamber within +-100 ns)** | **7.5-12.7 %** | **same size in every file and chamber; = the size of the real loss** |
| near_miss (hits in the window, nearest 15-40 mm) | 1-42 % | varies by FILE (MID bscan_hvon10 42 % vs Km18 2 %): prediction offsets, not losses; the counting estimator is immune |
| far_only (> 40 mm) | 0.7-2.3 % | |

Per chip, summed over the configurations, `chamber_silent` is **24 % for chip 2 (P2_IN) and 26 % for chip 14
(P2_OUT)**, against 5-13 % on all other chips.

**Reading:** about one muon in ten crosses a chamber that records nothing at all within +-100 ns, independent of
threshold and gain. That is a loss of the whole event in that chamber, not of a small pulse. Two explanations remain:
(a) acceptance: the prediction falls on a chamber edge or a region outside the read-out (chips 2 and 14 would be such
edges), or (b) read-out: the hits are lost between the chip and the capture (bursts of lost data, a hybrid not
sending for a while). The next analysis step (x-y map and time profile of the silent muons) and the pulser runs of
the weekend separate the two.

## Summary

| hypothesis for the missing ~10-15 % | verdict | evidence |
|---|---|---|
| detector gain / HV | ruled out | Dream > 90 % at the same HV (user); mesh at maximum |
| sparks / instabilities biasing the measurement | ruled out | veto changes the efficiency by +0.001 +- 0.007 |
| threshold too high (small pulses lost) | very unlikely | no threshold dependence; spectrum far above threshold at 9 mV/fC without an efficiency gain; to be made quantitative per channel |
| VMM gain | ruled out as a lever | 6 and 9 mV/fC are not better; 12/16 unusable (oscillation) |
| time window | minor | +1-3 points at +-1.5 us (4.5), up to +8 on IN at 9 |
| channel dead time | ruled out | ~0 |
| statistics | not at chamber level | errors 0.015-0.03; per-chip outliers are systematic |
| **whole chamber silent for ~10 % of the muons** | **open, dominant** | 7.5-12.7 % everywhere, 24-26 % on chips 2 and 14 |

Next steps: `HANDOFF_analysis_20261119.md` (this folder); data for the weekend:
`DAQ_Control_VMM_Beam/docs/HANDOFF_WEEKEND_RUNS_20261119.md`.
