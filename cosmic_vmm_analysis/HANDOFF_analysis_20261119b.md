# Handoff: VMM cosmic analysis, state at 19 Nov 2026 14:10 (continues HANDOFF_analysis_20261119.md)

Read `HANDOFF_analysis_20261119.md` first (rules, data locations, task list). This file says what was done on 19 Nov,
what was found, and the two decisions the user still has to take. Data taking runs in parallel until Mon 23 Nov 07:00
(`DAQ_Control_VMM_Beam/docs/HANDOFF_WEEKEND_RUNS_20261119.md`; run list with short names in the box file
`~/trim_work/out/weekend_20261119/runs.csv`).

`A=/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus`

## New code (main, one commit per file)

| commit | file | what |
|---|---|---|
| cb16da9 | `silent_muons.py` | chamber-silent muons: x-y maps, distance to the chamber outline and to chip boundaries, per chip, time structure (clustering, dispersion, prev_chip, chamber rate) |
| 22b4834 | `pulser_readout.py` | read-out efficiency from the weekend test pulses (lattice per pulsed channel, gaps, coincident losses, pulses around silent muons) - **its lattice result is wrong until the tick fix below** |
| 7c9a583 | `efficiency_tp.py` | `--dump-diag` adds `p*_chip_n1us`, `p*_chip_n100us` (hits on each of the probe chamber's 6 chips, columns = chips in ascending order, within +-1 us / +-100 us of tag + latency) and `adchist_h` (32 x 64 x 64, per-channel ADC histogram of all hits, 16-ADC bins). Old keys bit-identical (checked on chunk 26801). |
| b2f0da2 | `adc_threshold.py` | per-channel ADC edge (p01 and 50 % rising edge), own-pad muon ADC minus edge, loss fit - **the Moyal fit does not converge (values 1.0): replace it** |

EOS copies refreshed: `$A/weekend_20261119/pipeline/cva/` (used by the weekend loop: new weekend dumps carry the new
fields), `$A/diag_20261118/code/`, `$A/weekend_20261119/redump_test/cva/`.

## Data products made on 19 Nov

- Re-dumps with the new fields: 18 Nov runs 231-267 in `$A/diag_20261118/tags2/` (old `tags/` kept); weekend runs
  268-274 re-dumped in place in `$A/weekend_20261119/tags/`; later weekend runs dumped by the loop with the new code.
- Task 3: `$A/diag_20261118/results/all_thr`, `veto_thr` (18 Nov re-measured with the gain maps; build_thr log `results/build_thr_gain.log`).
- Task 1: `$A/diag_20261118/results/silent/` (summary, maps, csv); `$A/weekend_20261119/analysis/silent/` (weekend part, run of 09:00).
- Task 2: `$A/diag_20261118/results/adc/`, `$A/weekend_20261119/analysis/adc/` (`adc_edges.csv`, `adc_own_minus_edge.csv.gz`, spectra png).
- Task 5: `$A/weekend_20261119/analysis/pulser_run270/` (wrong lattice, see below).
- Condor scripts: lxplus `~/diag_20261118/{remeasure_thr,silent_job}.{sh,sub}`, `~/weekend_20261119/redump/{redump_run,adcthr,pulser,test}.{sh,sub}`.
  Weekend loop: cluster 13752682 (bigbird25), `~/weekend_20261119/loop.sh`, until Mon 18:00.
- Quick test scripts (not in git): lxplus `~/weekend_20261119/{tick_test,cmpdec,pgap,adcfine,adcnear}.py`.

## Findings

1. **Task 3 (threshold maps)**: rules for gain 3.0 / 6.0 / 9.0 in `weekend_20261119/pipeline/maps/config_maps.json`
   (maps `map_eq_g{3.0,6.0,9.0}.json` + their gscan_eq ext files). Realised dist_mV - bscan_target's target: median
   +0.5 mV for every gain file (as for 4.5); trim-clipped channels at 9.0 (341-375 of 1152) sit +12..18 mV (up to +45)
   above target. **In charge (dist_mV / gain): 9.0 points at 5-7 fC, 4.5 and 6.0 at ~4 fC** - the 18 Nov gain scan
   raised the threshold in charge together with the gain.
2. **Task 1 (silent muons)**, 18 Nov, 79 574 muons: silent fraction 14-15 %. P2_IN and P2_OUT strongly at the edges
   (0.36 / 0.29 at 0-5 mm inside the pad hull, 0.05-0.06 beyond 60 mm), a band ~20-30 mm wide along the lower long edge
   and the right edge, inside the pad area (pads exist, no signal) - candidate: area not under HV (the HV connectors were
   moved to sectors 2-3-4 on 26 Oct, overlap with the read-out area never confirmed) or outside the active area.
   Chips 2 / 14 high only because they sit at that edge (0.33 / 0.31 -> 0.18 / 0.06 beyond 20 mm). **P2_MID flat ~0.14
   everywhere** (not geometric). No time clustering (P(silent after silent) = P(silent after not silent) at every gap),
   per-run counts binomial (dispersion ~0.9), no chamber-rate dependence; prev_chip < 1 us: 0.24-0.46 (few muons).
   Per configuration: 9.0 files most silent (0.16-0.24), Km18 least (~0.12). Zero-silent chunks are binomial (~16 muons).
   **Pending user decision: tighter fiducial (> 20 mm inside the outline) on IN / OUT? The user said "show me first".**
   Still to do: whole-hybrid silence with `chip_n1us` / `chip_n100us` (dumps exist now).
3. **DECODER TIME BUG** (`coincidence_check.read_times`, used by every lxplus dump): `t = srs_h * 25.0 + chip`; the SRS
   marker tick is **22.5 ns** (frame = 4096 x 22.5 ns = 92.16 us). Effect: times stretched by 25/22.5 = 1.111 with
   fake ~256 us gaps (a 60 s capture spans 66.67 s; with 22.5 ns 60.01 s, the box p2-decode says 59.996 s live; pulse
   rate 9761 Hz -> 10845 Hz, box 10855 Hz). Hit counts agree with the box decoder (7 586 146 vs 7 592 028). Chip-to-chip
   timing is NOT affected (pulses of the 6 pulsed chips coincide within 200 ns in 99.8 % with either tick, median 14-23 ns)
   -> the efficiency (tag/probe within one frame) is fine. WRONG: wall clock (each file anchored at its first packet, error
   grows to ~6.7 s at the end of a 1-min file) -> spark-veto alignment (the 18 Nov "veto has no effect" result needs
   redoing), per-second rates, prev_chip / prev_chan, chip_n100us (mildly), the pulser lattice.
   **Pending user decision: fix to 22.5 and re-dump everything (18 Nov + weekend, ~45 condor jobs, ~1 h)?** If yes: one
   line in coincidence_check.py, refresh the three EOS code copies, re-dump (`redump_run.sh` pattern), re-run
   spark_veto/measure-with-veto, silent_muons, pulser_readout. Mid-weekend the weekend pipeline would mix old/new dumps
   unless the weekend tags are re-dumped too.
4. **Task 5 (pulser), first look**: run_270 (PRa), one file: pulses on all 6 chips coincide 99.8 % -> read-out loss
   ~0.1-0.2 %, far below the silent fraction -> silent muons are not read-out data loss (at least there). The
   pulser_readout.py result (eff 0.0999) is the tick bug (the lattice breaks at the fake gaps): re-run after the fix, or
   switch it to gap counting (round(gap / 92.16 us) - 1 missing) which is robust to it. Chip 9 of set B pulses ch 48
   (ch 16 gave no pulses in C10, changed 19 Nov 02:12; `maps/pulsed_channels.json` updated).
5. **Task 2 (ADC threshold)**: the p01 / half edges land on 48 / 64 because of the ADC DNL (codes 48, 64, 80 collect
   ~2x) and the 16-ADC bins; at full resolution (box `hits_store/*/counts.npz` 'adc', per VMM) the edge moves ~0.6 ADC per
   mV of threshold on an offset of ~40 ADC (VMM 3: gain 3.0 ~10 mV 48, Km18 53, B10 60, C10 73). Use p01. Fit-free result
   (p01 edge), fraction of detected muons within 0-32 ADC above the edge: **9.0: 3-7 %, 6.0: 4-11 %, 4.5: 7-19 %,
   3.0: 12-28 %** (P2_OUT highest); 5th percentile of ADC - edge 11-19 (3.0) -> 25-48 (9.0) ADC. The efficiency does
   not rise with the gain -> by the user's criterion **gain is not what limits the efficiency**. Caveat: the gain-3
   thresholds are at ~3-4 fC, close to REF in charge. To do: replace the Moyal fit by a fit-free extrapolation, or the
   9.0 shape scaled by the gain ratio; per chip.

## Next steps (suggested order)

1. Get the user's answers on the tick fix and the fiducial (above).
2. Tick fix + re-dump (if yes), then redo the veto result of 18 Nov and the pulser readout on all pulser runs
   (PRa/PRb/PCa/PCb, runs 270, 283, ...; runs.csv gives the short names).
3. Whole-hybrid silence (task 1 last item) from `chip_n1us`/`chip_n100us`: silent muons with the other chips of the
   hybrid / chamber active or not.
4. Task 4 (per-chip outliers), task 6 (weekend data as it comes), task 7 (FINDINGS for the weekend).
