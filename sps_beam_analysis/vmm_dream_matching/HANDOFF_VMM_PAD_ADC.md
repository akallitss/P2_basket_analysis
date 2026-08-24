# Handoff — VMM3a inefficiency vs per-pad median ADC (SPS H4, July 2026)

**Written 2026-08-21.** For whoever picks up the VMM efficiency problem.

The one-line version: **the VMM-read efficiency of the P2 telescope is limited by
the discriminator threshold, not by the chambers, and the sharpest remaining
handle on it is the per-pad median pulse height.** That handle is measured and
the correlation is real, but the code that produced the headline table was never
committed. Reproducing it is task 1 below.

---

## 1. The problem

At the best configuration of the campaign, the same three chambers read out by
DREAM and by the VMM3a give very different efficiencies:

| station | VMM (uRWELL-referenced) | DREAM, same chambers, same working point |
|---|---|---|
| P2_IN  | 0.155 | 0.9649 |
| P2_MID | 0.669 | 0.9706 |
| P2_OUT | **0.854** | **0.9604** |

Reference run for the VMM column: `run_46/cfg_gain4.5_peaktime200`
(2026-08-01 04:53, mesh 450 V / drift 750 V, gain 4.5 mV/fC, peaking time
200 ns, Ar/CO2/iC4H10 93/5/2, 48 captures, 1.52e9 hits, 837 338 matched uRWELL
tracks). DREAM reference: `highstat_eff_1/beam_commissioning_00` (25 Jul) and
two further runs, 0.955–0.966 over the three.

**P2_IN is a separate, unresolved problem** (three chips at raised thresholds,
VMM 4 carrying 1.14e9 hits from one screaming channel, a 20–30 mm shoulder in
the matching distance that says the channel->pad map is partly wrong). Its
number is a lower limit on the chamber. Do not fold it in with MID/OUT.

## 2. What is already established, and what it rests on

Excluded **by measurement**, each with a number in `EFFICIENCY_AUTOPSY.md`:
fiducial/active area (100% of tracks within 9 mm of an instrumented pad),
geometry and mapping (-59.98 deg rotation, 3.4 mm residual = pad/sqrt(12), and
*the VMM finds the same dead pad DREAM does* — pad 635, VMM 17 ch 4), the
coincidence window (3 sigma -> 20 sigma moves 0.8540 -> 0.8579), position
matching, time in spill, instantaneous rate, time since the previous trigger,
track quality, capture coverage, intra-pad position, cross-station common mode.

Established **as the cause** — five independent handles, none of which involves
the ADC:

1. threshold DAC alone (`run_47` vs `run_48`, same config file with the `sdt`
   lines commented out vs set, 36 min apart): P2_OUT 0.486 -> 0.822
2. amplifier gain 3.0 -> 4.5 mV/fC at fixed HV: 0.723 -> 0.822
3. within P2_MID, chips at `sdt` 256 read 0.53/0.31, chips at 224 read 0.70–0.78
4. the mesh scan and the gain step fit one Landau, V0 = 22.4 V (a physical bulk
   Micromegas value, and an output of the fit, not an input)
5. **per-pad efficiency correlates with per-pad median pulse height at +0.55**

and the punchline: **the best P2_OUT pads already reach 0.962 / 0.957 / 0.954,
i.e. the DREAM number.** Nothing is wrong with the chamber.

## 3. The specific question this handoff is about

Handle 5 is the one worth pushing, because it is per-pad rather than per-chip
and therefore has ~70x the statistics of the `sdt` comparison. Measured over the
69 illuminated P2_OUT pads (P2_MID: +0.52 over 54 pads):

| pad median pulse height [VMM ADC] | efficiency | tracks |
|---|---|---|
| 90–127  | 0.752 | 151 914 |
| 127–139 | 0.816 | 241 640 |
| 139–168 | 0.892 | 163 012 |
| 168–202 | 0.920 | 174 252 |
| 202–271 | 0.921 | 106 159 |

**Read the direction carefully.** If this were a selection effect it would run
the other way: a channel with a *high* threshold records only its big pulses and
would show a high median with a low efficiency. The correlation is positive, so
it is real signal-size variation across pads, and it dominates.

Open questions worth someone's time:

* Is the spread in per-pad median ADC **gain** (gas gain, mesh planarity, pad
  area) or **electronics** (per-channel baseline scatter under one global
  `sdt`)? The pad geometry is known exactly, so pad area and radius can be
  divided out — if the residual spread survives that, it is electronics.
* Does the residual correlate with position on the chamber (mesh support
  pillars, HV feed corner) or with channel index within a connector (cabling,
  baseline)? Those two hypotheses make different maps.
* What does the model predict for a per-channel trim? `gain_threshold_model.py`
  fits threshold/MPV = 0.68 and can be run per pad.

## 4. Task 1 — the table above is NOT reproducible from committed code

`eff_autopsy_report.py` stores, per pad, only:

```
res["per_pad"] = {"channel_id", "x", "y", "vmm", "ch", "n", "k"}
```

`n` = tracks on that pad, `k` = tracks with a hit. **There is no per-pad ADC.**
The quintile table and the +0.55 were computed ad hoc in an interactive session
and never committed. Before doing anything else, add per-pad ADC to that dict —
the input is already in hand, `autopsy_<tag>_tracks.parquet` carries
`win_adc_<station>` (the ADC of the nearest in-window unmasked hit) alongside
`nearpad_<station>`, so it is a groupby, not a new measurement.

Suggested: store per pad the median, the 25/75 percentiles and the count, then
recompute the correlation and the quintiles from the stored values so the number
in `EFFICIENCY_AUTOPSY.md` has a producer.

## 5. Task 2 — the ADC scale has a known defect, correct for it first

**The VMM3a 10-bit ADC has strong differential non-linearity with a period of
exactly 16 codes.** Every multiple of 16 is ~2x too wide and absorbs its
neighbours. Measured on run_46 (relative code width by phase, 16 parameters
fitted to ~700k hits):

```
phase:  0     1     2..15
width:  2.03  1.31  0.83–0.97     P2_OUT   13.5% of hits on multiples of 16
        2.10  1.00  0.86–0.95     P2_MID   14.9%          (6.2% if flat)
```

The pattern is identical on two independent stations, so it is a property of the
ASIC's ADC, not of a chip. Consequences:

* **A unit-bin mode returns a comb tooth, not the Landau peak.** This is how
  `EFFICIENCY_AUTOPSY.md` came to quote "most probable 128" for P2_OUT
  (128 = 8x16) and 80 for P2_MID (5x16). Corrected: **P2_OUT MPV = 104**
  (rebin-by-16) or **102** (periodic code-width correction) — two independent
  estimators agreeing. P2_MID ~ 88. **The 128 in that document is wrong and
  should be fixed.**
* Any *median* is much less affected than a mode, so the quintile table in
  section 3 is probably close to right — but recompute it on corrected codes
  before quoting it again, since the quintile *boundaries* sit near teeth.

**It does not explain the inefficiency**, and it is worth knowing why so nobody
re-opens it: the discriminator is an analog comparator upstream of the ADC, so a
code-width error cannot un-record a hit; our efficiency applies no ADC cut at
all (`urw_vmm_efficiency.py` unpacks `hadc` and never uses it); and the lowest
recorded pulse, ADC 47, is at phase 15, not on a boundary. All five handles in
section 2 are analog and untouched by digitisation.

Correction options, cheapest first: rebin to the DNL period (assumption-free,
costs resolution); fit the periodic code-width model in-situ from a spectrum
known to be smooth (done above, 16 parameters, hugely over-determined); or do it
properly with a code-density / histogram test (IEEE Std 1241) driven by the
**VMM's own internal test-pulse DAC** — which is the same bench run that would
give the per-channel threshold trim, so it buys both problems at once. No such
calibration run exists for this campaign (`qa_config.py: 'calibration': None`).

## 6. Data

Everything is on EOS; nothing needed here lives only on a laptop.

| what | where |
|---|---|
| VMM raw + reduced, per run/sub_run | `/eos/experiment/ntof/data/x17/p2_sps_july/vmm/runs/<run>/<sub_run>/` |
| — raw captures | `.../raw_daq_data/*.pcapng` |
| — reduced counts (additive histograms) | `.../hits_store/<capture>/{counts.npz,scalars.json}` |
| — full column store, only some sub_runs | `.../hits_store/<capture>/*.npy` |
| VMM<->DREAM stream matching | `/eos/experiment/ntof/data/x17/p2_sps_july/vmm/matching/` |
| VMM efficiency sweep (corrected) | `.../vmm/matching/efficiency_v2/` |
| rendered VMM QA plots, all 75 sub_runs | `/eos/experiment/ntof/data/x17/p2_sps_july/vmm/qa_plots/index.html` |
| DREAM analysis products | `/eos/experiment/ntof/data/x17/p2_sps_july/analysis/` |
| uRWELL detector config (needed by the autopsy) | `/eos/experiment/ntof/data/x17/p2_sps_july/config/detectors/` |

**Trap:** most sub_runs were reduced with `--drop-columns` and keep only
counts/scalars; the whole gas-B period (run_61 onward) has no `hits_store` at
all. Sub_runs that keep full columns: run_32/33/34/35/36, 38, 40–43, 49–52, 56.
Anything else must be decoded from pcapng first (`decode_to_store.py`, ~12 s per
capture). run_46 *does* keep columns, which is why it is the reference run.

## 7. Code

**On lxplus** (`akallits`, AFS home — read-only to others, ask for a copy):

| path | what |
|---|---|
| `~/p2_eff/eff_autopsy.py` | one streaming pass over a sub_run, all three stations; writes the tracks + hits parquet |
| `~/p2_eff/eff_autopsy_report.py` | slices the parquet every way; **this is the file task 1 edits** |
| `~/p2_eff/make_vmm_adc_hist.py` | track-matched ADC spectrum -> npz, with the DNL-robust MPV |
| `~/p2_eff/run_autopsy.sh`, `run_report.sh`, `run_adc_hist.sh` | LCG_110 + PYTHONPATH wrappers; use these, not bare python |
| `~/vmm_render/render_from_store.py` | rebuilds the QA plots from the counts store |

**Traps:** `run_autopsy.sh` writes to `${TMPDIR:-/tmp}/autopsy.$USER`, which is
**node-local and wiped** — copy the parquet somewhere durable or expect to
re-run (~6 min for run_46, of which 262 s is reading 1.52e9 hits over EOS).
`set -u` breaks the LCG `setup.sh`; the wrappers already work around it.

**On GitHub** — the code not on lxplus, and the authoritative copy of what is:

| repo | holds |
|---|---|
| `github.com/akallitss/P2_basket_analysis` | the analysis: `sps_beam_analysis/vmm_dream_matching/` (autopsy, matching, `urw_vmm_efficiency.py`, `gain_threshold_model.py`, `EFFICIENCY_AUTOPSY.md`, `TIMING.md`), `mpgd2026/` (figure producers incl. `make_adc_comparison.py`) |
| `github.com/akallitss/DAQ_Control_VMM_Beam` | the DAQ and online QA: `vmm_qa/vmm_stations.py` (cabling, pad map, masking), `vmm_reduce.py` (what `counts.npz` is), `render_from_store.py`, `qa_config.py` |

Mirrors of the analysis repo also exist at `drf-gitlab.cea.fr/p2/vmm-analysis`
and `.../vmm-analysis-extended`.

Two documents to read before touching anything:
`sps_beam_analysis/vmm_dream_matching/EFFICIENCY_AUTOPSY.md` (the full argument,
including everything that was excluded and how) and `.../COMPARISON_CLAIMS.md`
(which DREAM number is comparable to which VMM number, and why most pairs are
not).

## 8. Reproducing the reference numbers

```bash
ssh lxplus
~/p2_eff/run_autopsy.sh run_46 cfg_gain4.5_peaktime200     # ~6 min
#   -> $TMPDIR/autopsy.$USER/autopsy_run_46_cfg_gain4.5_peaktime200_{tracks,hits_*}.parquet
#   expect: P2_IN 0.1535  P2_MID 0.5916  P2_OUT 0.8540

~/p2_eff/run_report.sh $TMPDIR/autopsy.$USER/autopsy_run_46_*_tracks.parquet \
    --mask rate --figdir $TMPDIR/figs
#   --mask {ratio,rate,none} compares the hot-channel rules on one sample.
#   P2_MID is 0.5916 under `ratio` and 0.6688 under `rate`; `rate` is correct,
#   see task-3 note below.

python3 ~/p2_eff/make_vmm_adc_hist.py \
    $TMPDIR/autopsy.$USER/autopsy_run_46_*_tracks.parquet --station P2_OUT \
    -o vmm_adc_P2_OUT_run46.npz
#   expect: lowest 47, 5th pct 81, median 153, unit-bin mode 128 (a DNL tooth),
#           DNL-robust 104, 13.5% of hits on multiples of 16
```

The figure that puts DREAM and VMM pulse height on one axis is
`mpgd2026/make_adc_comparison.py --vmm-hist <that npz> --overlay`; it needs the
DREAM spectrum from `analysis/P2_OUT/highstat_eff_1/beam_commissioning_00/
20_beam_spectra/scan_row.json`.

## 9. Two traps that already cost time once

* **The auto hot-channel mask deletes beam.** The original rule masked any
  channel above 8x its chip's median occupancy. On a chip whose threshold has
  been raised, the noise is gone, the median collapses, and the *illuminated*
  pads stand 30–50x above it — so the rule masks the signal. It cost P2_MID 7.7
  points (0.5916 -> 0.6688) and made VMM 9 look dead when it was only masked.
  Now an absolute rate cut (20 kHz). Out-of-time occupancy does **not** fix it,
  because the beam rate far exceeds the matched-trigger rate.
* **A few hits per capture carry `srs_timestamp = 0`**, so (min,max) capture
  spans stretched back to t=0, the coverage cut accepted everything and the live
  time came out 25x too long. Quantiles now.

## 10. Suggested order of work

1. Add per-pad ADC to `eff_autopsy_report.py`; regenerate the quintile table and
   the correlation so they have a producer (section 4).
2. Apply the DNL correction and recompute; fix the "most probable 128" line in
   `EFFICIENCY_AUTOPSY.md` (section 5).
3. Divide out pad area/radius and re-test the correlation — gain or electronics?
4. Map the residual against position and against channel index within a
   connector; the two hypotheses predict different maps.
5. Cost the fix: run `gain_threshold_model.py` per pad to predict what a
   per-channel trim would recover, as input to the request for calibration
   bench time.

Contact: Alexandra Kallitsopoulou (`akallits`).
