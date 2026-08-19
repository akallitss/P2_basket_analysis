# DREAM vs VMM3a — what we can say at MPGD26, and what we cannot

Status register for the cross-readout comparison, written 2026-08-19 (talk 4 Sep).
The analysis is three days old and has already had two bugs corrected in it, so the
question is not "what did we measure" but "which of these survives a hostile question".

## The strategic point

**The load-bearing claims do not need DREAM at all.**

The argument "the deficit is the discriminator threshold, not the chamber" rests on
three VMM-internal comparisons — same detector, same day, same reference, same cuts,
minutes to hours apart:

| handle | comparison | P2_OUT |
|---|---|---|
| threshold DAC alone | run_47 `_deflt` vs run_48 `_opt`, same config file, 36 min apart | 0.486 → 0.822 |
| amplifier gain alone | run_41 (3.0 mV/fC) vs run_48 (4.5), same HV/gas | 0.723 → 0.822 |
| per-pad spread at one DAC | 69 illuminated pads, all six chips at sdt 224 | 0.59 → 0.96 |

None of these is exposed to the 25 Jul ↔ 1 Aug systematic, to the chamber swap, or to
the gas change. DREAM supplies only the *ceiling*, and even the ceiling is corroborated
internally: the best individual VMM pads read 0.962 / 0.957 / 0.954 against DREAM's
0.9604 on the same plane in the same run.

**Build the slide so that if the cross-date DREAM comparison is attacked, the conclusion
still stands.** Lead with the internal handles; use DREAM as the ceiling, not as the
argument.

> **Tiers updated 2026-08-19** by the V1-V8 checks below; where a check
> moved a claim, the tier entry says so.

## Tier 1 — say it plainly

* The two DAQ streams are matched. 77 of 110 sub_runs lock at 113–606 σ, residual rms
  10.4 ns; every failure is run_29 or earlier, before the trigger fan-out existed.
  Self-diagnosing: an unlocked sub_run collapses to the accidental rate.
* The VMM stations were measured against the *same* uRWELL tracks, fiducial and probe
  radius as the DREAM campaign, over 78 sub_runs.
* Geometry is independently confirmed: freely fitted uRWELL→pad transform gives a proper
  −59.98° rotation with unit singular values, residual 3.37 / 3.44 mm = 12 mm/√12.
* **The two readouts find the same dead pad** — pad 635 (VMM 17 ch 4) at (411, 237) mm,
  the same dark square at (410, 237) in the DREAM map. Two electronics chains, same
  defect, same coordinates. This is the single most persuasive slide element and it costs
  nothing to defend.
* The three internal handles above.
* The found-hit pulse-height spectrum is a Landau cut off on its low side (lowest ADC 47,
  5th percentile 81, MPV 128) — it starts at a wall, it does not rise from zero.
* The deficit is **not** rate, pile-up, spill time, coincidence window, clustering radius,
  track quality, capture coverage, charge sharing at pad boundaries, or a common-mode DAQ
  loss. Each excluded by measurement; the null results are flat and boring, which is the
  point.
* P2_OUT is the number to quote: the corrected hot-channel rule left it untouched
  (0.8540 both ways).

## Tier 2 — say it with the caveat attached, or not at all

* **0.854 vs 0.960 as a bare pair of numbers.** The DREAM reference is 25 Jul; run_46 is
  1 Aug. In between: a chamber swap (27 Jul) and, by our own backup slide, reference
  levels 3–5 pp lower and re-ordered on 27 Jul. Quote it as *"85 % against a 96 % ceiling
  measured on the same chambers the week before"*, never as a difference of two exact
  numbers. → task V2.
* **P2_MID = 0.669.** Moved +4.0 points on average (up to +12.6 on 42 of 78 sub_runs)
  when the hot-channel rule was corrected three days ago. It is the least stable number
  in the set. Keep it in the table, do not build a claim on it.
* **P2_IN = 0.155.** A lower bound on the readout, not a measurement of the detector:
  three of six chips at raised thresholds, residual 4.7 mm against 3.4 mm elsewhere, and
  a 20–30 mm shoulder in the matching distance = partly wrong channel→pad map.
  → **V6 done:** greyed and labelled on the figure.
* **The mesh scans never plateaued** at the highest HV used. So 85 % is the efficiency
  *at this working point*, not "the VMM's efficiency". Say so in the same breath.
* ~~**Cluster size 1.11–1.18 vs 1.05–1.08.**~~ → **V7 done: withdrawn.** The VMM
  efficiency code deliberately does not cluster at all. Retracted everywhere it appeared.

## Tier 3 — do not put on a slide as a measurement

* **V₀ = 22.4 V and "the discriminator sits at 0.68 × the MPV".** This is a Landau
  survival-model fit whose own documentation warns that its asymptote is an artefact of
  the Landau's unphysical low tail. Beautiful, and the right physics — but it is the most
  model-dependent statement we have, and a single sharp question about the fit range
  would put the whole slide in doubt. Show the *measured* mesh scan and the *measured*
  gain step; put the model in backup and introduce it as "a Landau survival model
  consistent with a bulk-Micromegas gain e-folding" if asked.
* **Clopper–Pearson intervals.** → **V4 done, and it is worse than feared.** At identical
  config, HV and gas the same measurement walks **5.6 points (P2_OUT) / 10.3 (P2_MID)**
  over 24 h, against a ±0.0004 binomial interval and ±0.005 of analysis systematics.
  Print no error bar on a VMM efficiency until that walk is explained.
* **Any VMM timing performance number.** σ is pinned at the 22.5 ns BCID quantisation for
  want of time-calibration runs, and the SRS tick fix (22.5 not 25 ns) means any absolute
  time elsewhere may carry an 11 % stretch. State the limitation, claim nothing.
* **A campaign-average VMM efficiency.** The 78 sub_runs sit at different HV, gain,
  peaking time and threshold. There is no meaningful average.

## V1-V8: run 2026-08-19, results

`verify_configs.py` reads the conditions off each run's own `run_config.json`;
`efficiency_conditions.csv` is `efficiency_table.csv` joined to them.
EOS was unreachable (Kerberos ticket expired 17 Aug), so two sub-checks are
parked — both are named below.

**V1 - gas. CONFIRMED, concern retired.** Every run of the comparison
(run_38 ... run_48, and run_54) is `Ar/CO2/Iso 93/5/2` at mesh 450 / drift 750
on all three stations, read from the VMM DAQ's own config, not from the run
number. run_46 is 1 Aug but sits *before* the CF4 changeover (first gas-B run
is run_61, 2 Aug 15:23). The headline comparison is gas-matched.

**V2 - the DREAM reference. CHANGED, and it changed two numbers.** The
reference in use was `highstat_eff_1` at **drift 700 V**, while run_46 ran at
**drift 750 V**. Taking every uRWELL-referenced DREAM measurement at run_46's
own working point instead:

| station | DREAM at mesh 450 / drift 750 | n | epochs |
|---|---|---|---|
| P2_OUT | **0.955 - 0.972** (median 0.960) | 18 | 25 + 28 Jul |
| P2_MID | 0.910 - 0.968 (median 0.915) | 18 | 25 + 28 Jul |
| P2_IN  | 0.943 | 1 | 28 Jul |

* **P2_OUT survives intact.** The quoted 0.960 lands on the median of 18
  independent measurements across two epochs. Quote the band.
* **P2_IN's reference was the wrong chamber.** 0.965 is the 25 Jul chamber,
  swapped out on 27 Jul; run_46 read the *replacement*, whose DREAM reference
  is **0.943**. Corrected in the figure and the deck.
* **P2_MID's DREAM reference is itself unstable** - 0.968 on 25 Jul against
  0.910-0.927 on 28 Jul, at identical HV. Confirms the Tier-2 flag: no claim
  may rest on P2_MID.

**V3 - the run_47/48 controlled experiment. PARTLY CONFIRMED, one check
parked.** Everything at run level is identical and verified: same gas, same
mesh 450 / drift 750, same gain, same peaking time, 36 minutes apart, and the
two chip-config filenames differ only by the `_deflt`/`_opt` suffix. Proving
they differ *only* in the per-chip `sdt` lines needs the two
`p2b-config-cern-ext_*.txt` files, which are on the DAQ machine and EOS, not
here. **Strengthened meanwhile:** the deflt/opt pair is not one experiment but
**five** on gas A, each 30-90 min apart - run_42/43, run_47/48, run_50/51,
run_52/53 - plus run_64/65 and run_66/67 on gas B. All six move the same way.
One config file being mislabelled would not do that.

**V4 - the error budget. This is the big one, and the answer is not the one we
expected.** Analysis systematics are small: probe radius 8/15/40 mm gives
0.850/0.854/0.857, the coincidence window 3 sigma -> 20 sigma gives
0.854 -> 0.858, and the hot-channel rule leaves P2_OUT untouched. Envelope
about +-0.005. **But repeatability is far worse than that.** At *identical*
chip config, HV and gas (`gain3.0_peaktime200`, mesh 450 / drift 750, gas A):

| | P2_IN | P2_MID | P2_OUT |
|---|---|---|---|
| within one run, minutes apart | 0.3 | 0.3 | 0.2 pts |
| 31 Jul 18:17 -> 1 Aug 18:30 | 3.6 | **10.3** | **5.6 pts** |

P2_OUT walks 0.734 -> 0.680 monotonically over 24 h at settings that never
changed. **What it is not:** the per-chip ADC percentiles are flat across the
whole window (P2_OUT p50 129/146/128/144/154/121 at the start against
132/154/127/143/154/119 at the end) and the live-channel counts are flat, so
it is neither a gain droop nor a threshold change. The masked-hit fraction
doubles over the same window (0.33 -> 0.78), which points at the analysis's
channel masking rather than the chamber - unresolved. Figure:
`mpgd2026/figs/vmm_stability.png`.

**Consequence: never print a Clopper-Pearson interval next to a VMM
efficiency.** +-0.0004 against a quantity that moves 5 points in a day is the
most attackable thing on any slide. Quote the number bare, or with the
repeatability spread.

**V5 - run_46 vs run_48. RESOLVED, and it is another threshold point.** They
are *different* chip config files (`..._peaktime200.txt` vs
`..._peaktime200_opt.txt`) at the same gain, peaking time, HV and gas, 72 min
apart: 0.854 against 0.822. run_46's P2_OUT chips are all at sdt 224 and the
autopsy reads run_48 at sdt 230 on five of six - higher threshold, lower
efficiency, the right direction. So this is *not* an unexplained
reproducibility term, and it does not enter V4's budget.

**V6 - P2_IN. DONE.** Greyed on `dream_vs_vmm.png` and labelled "LOWER BOUND,
not a measurement of the chamber". Owning it out loud costs ten seconds;
hiding it invites the question.

**V7 - cluster size. WITHDRAWN.** Not like-for-like, and not marginally so:
`urw_p2_efficiency.py` builds a charge-weighted leading cluster (leading pad
plus mapped pads within `cluster_r`), while `urw_vmm_efficiency.py`
deliberately does not cluster at all - its own comment reads *"Not 'is the
leading cluster close'"* - and counts any pad within the probe radius. The
1.11-1.18 vs 1.05-1.08 comparison is retracted in the report, the options deck
and the draft deck. Reinstate only if the VMM side is reclustered the same way.

**V8 - the contradiction between the two decks. DONE.** The 16 Aug options
deck's *"Why DREAM and VMM efficiencies must not share an axis"* slide and the
matching figure caption in the report both now say they are superseded, and
carry the three caveats that *do* survive.

### Still to run when Kerberos is back

```bash
kinit akallits@CERN.CH
```

1. **V3's chip-config diff** - fetch `p2b-config-cern-ext_gain4.5_peaktime200
   {_deflt,_opt}.txt` and confirm the only difference is the `sdt` lines.
2. **V4's cause** - the masked-hit fraction doubling. Re-run
   `urw_vmm_efficiency.py` on run_36 and run_58 with masking disabled; if the
   5.6-point walk disappears it is ours, and the number is recoverable.
3. **The VMM ADC histogram** for `adc_threshold_cost.png` panel (a):
   `python3 mpgd2026/make_adc_comparison.py --vmm-hist run46_P2_OUT_adc.npz`
   turns the drawn low edge into the measured spectrum.

## Figures produced

| file | what it shows |
|---|---|
| `mpgd2026/figs/adc_threshold_cost.png` | the framing figure - the DREAM Landau with the VMM's recording window on it, and per-pad efficiency vs per-pad pulse height |
| `mpgd2026/figs/dream_vs_vmm.png` | updated: drift-matched DREAM bands, P2_IN greyed |
| `mpgd2026/figs/vmm_stability.png` | the error budget - 5.6 points at fixed settings |

## Original task list (for the record)

## Tasks before the 26 Aug dry run

Priority ordered. Total ≈ 1.5 days.

**V1 — confirm the gas of run_41/46/47/48 from the VMM config, not the run number.**
30 min. The gas changed on 1 Aug and run_46 *is* 1 Aug; the autopsy asserts gas A. The
two DAQs number runs independently and drift apart after run_58, and 136 EOS gas labels
were already wrong once. Read it off the config file and write the answer down. A wrong
answer here invalidates the headline comparison, so it is the first thing to check.

**V2 — put a band on the DREAM reference.** ~1 h. We have two DREAM efficiency epochs
(25 Jul high-stat, 27 Jul) and know they differ by 3–5 pp. Quote the 25 Jul number with
the 27 Jul one as the spread, so the ceiling is a band, not a point.

**V3 — diff the run_47 and run_48 config files.** 30 min. This is the load-bearing
controlled experiment. Confirm they differ *only* in the per-chip `sdt` lines — no gain,
peaking, HV, latency or masking difference. If anything else moved, the cleanest argument
in the talk needs rewording.

**V4 — a systematic band on 0.854.** ~half a day, highest value. Vary what we already
know how to vary and quote the envelope: hot-channel rule (ratio vs absolute), probe
radius (8 / 15 / 40 mm → 0.850 / 0.854 / 0.857), coincidence window (3σ / 20σ → 0.854 /
0.858). Most of these numbers are already in the autopsy — this is assembly, not new
analysis. Output: **ε(P2_OUT) = 0.854 ± <syst>**, which is what makes the number quotable.

**V5 — resolve run_46 vs run_48.** 30 min, folds into V3. Same nominal gain and peaking
(4.5 mV/fC, 200 ns) yet 0.854 vs 0.822. run_46's P2_OUT chips are all at sdt 224;
the autopsy says run_48 sits at sdt 230 on five of six. If that is right it is another
threshold data point and strengthens the argument. If the configs match, we have a
3.2-point unexplained run-to-run reproducibility term that must go into V4's band.

**V6 — decide P2_IN's fate on the figure.** 15 min, editorial. Recommendation: keep the
point, grey it, label it "lower bound — channel map under revision", and say the sentence
out loud. Hiding it invites the question; owning it costs ten seconds.

**V7 — cluster size like-for-like, or cut.** 30 min.

**V8 — kill the contradiction between the two decks.** 15 min. The 16 Aug options deck
(`reports/mpgd26_sps_beam_2026-08/slides.html`) still carries a backup slide titled
*"Why DREAM and VMM efficiencies must not share an axis"*, describing the per-pad
acceptance-restricted efficiency as the missing step. That step is done, and the 17 Aug
draft deck now puts them on one axis. Both are published. Retire or rewrite the old
slide before anyone reads them side by side.

## The sentence to defend

> On the same three chambers, against the same external tracker and the same cuts, the
> DREAM readout gives 96 % and the VMM3a readout gives 85 % — and that gap is the VMM's
> discriminator threshold, not the chamber: the threshold DAC alone moves it from 49 % to
> 82 %, and the pads that happen to see the largest pulses already read 96 % in the very
> same run.

Everything in that sentence is Tier 1. V2 turned the 96 % into a band — 0.955-0.972 at
run_46's own working point, 18 measurements over two epochs — which is what "96 %" should
be read as. Say the *comparison*; do not attach an error bar to the 85 % (V4).

## What we are not claiming, and should say so once

The VMM was run with **one global threshold per chip and the per-channel trim DACs never
used**, and no time-calibration runs were taken. This is a statement about how we
configured the chip in a two-week beam test, not about the VMM3a. Said plainly and early,
it converts the whole result from "the VMM underperforms" into "here is what it costs to
skip the trim, measured" — which is a more useful result and a much easier one to defend.
