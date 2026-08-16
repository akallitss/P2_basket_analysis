# Why the VMM-read efficiency is not 96%

**run_46 / cfg_gain4.5_peaktime200, 2026-08-01, taken apart event by event.**

The uRWELL-referenced efficiency of the VMM-read P2 stations came out at 85%
on P2_OUT where the DREAM readout of the same detectors, with the same
reference and the same cuts, gives 96%. That is either the headline of the
campaign or an artefact of the measurement, so one sub_run was taken apart
completely.

**It is the VMM discriminator threshold, and nothing else.** The threshold sits
at roughly the most probable signal, so a few per cent of gain moves the
efficiency by ten points, and the pads that happen to see the largest pulses
already measure 96.2% — the DREAM number — on the same plane, in the same run.

---

## The run

`run_46/cfg_gain4.5_peaktime200`: mesh 450 V, drift 750 V on all three
stations — the nominal operating point, the same HV as `run_36/operating_*`
and `run_32/meshscan_m00V` — with the best VMM configuration of the campaign
(gain 4.5 mV/fC, peaking time 200 ns). It is the highest P2_OUT efficiency
measured anywhere in the campaign, so any deficit found here is a floor on the
problem, not a badly-chosen working point.

837,338 uRWELL tracks carry a matched VMM trigger time, 48 captures,
1.52 × 10⁹ VMM hits.

**Reference.** `highstat_eff_1/beam_commissioning_00`, 2026-07-25, same gas
(Ar/CO2/Iso 93/5/2), mesh 450 V, the same detectors read out by DREAM, the same
uRWELL tracks, the same analysis and the same cuts (fid_r 9 mm, probe_r 15 mm):

| station | DREAM | VMM (run_46) |
|---|---|---|
| P2_IN  | 0.9649 | 0.155 |
| P2_MID | 0.9706 | 0.669 |
| P2_OUT | 0.9604 | 0.854 |

P2_MID and P2_OUT are the same physical detectors on both dates. P2_IN is not
necessarily — a detector was swapped out on 27 July — so its two numbers are not
strictly the same object, and P2_IN is treated separately throughout.

---

## What it is not

Each of these is a measurement, not an argument.

**Not tracks counted outside the active area.** All 837,338 tracks land within
9 mm of an instrumented pad. The projected beam is a ~110 × 110 mm patch
covering 69–73 pads in the middle of the 384-pad instrumented region, nowhere
near its edge. The identical fiducial definition gives DREAM 96%.

**Not the geometry or the channel map** (on P2_MID and P2_OUT). The
uRWELL → pad transform is a proper −59.98° rotation with unit singular values,
and the residual is 3.37 / 3.44 mm = 12 mm/√12, the pad quantisation — i.e. the
pointing is as good as pads allow. Better: **the VMM finds the same dead pad as
DREAM**. Pad 635 (VMM 17 channel 4) at (411, 237) mm reads 29.8% here, and is
the same dark square at (410, 237) in the DREAM map. Two completely different
electronics chains, same defect, same coordinates.

**Not the coincidence window.** Widening it from 3σ to 20σ moves P2_OUT from
0.8540 to 0.8579. For the tracks that miss, the hits within 15 mm of the
prediction are flat across ±10 µs at the accidental level (~72 per 250 ns bin):
there is no out-of-time signal population being cut away.

**Not position matching or clustering.** The efficiency asks for *any* pad
within 15 mm, and it plateaus by 8 mm (0.8496 at 8 mm, 0.8567 at 40 mm). For
95.7% of the misses there is nothing within 15 mm at any time in ±10 µs — the
pad simply did not fire.

**Not rate, dead time or pile-up.**

| sliced by | P2_OUT efficiency |
|---|---|
| time into the spill, 0 → 6 s | 0.855 … 0.856 |
| triggers within ±1 ms, 0 → 12 | 0.853 … 0.856 |
| time since the previous trigger, <10 µs → >1 ms | 0.851 … 0.855 |

Flat everywhere. And VMM 19, whose channel 20 runs at 116 kHz, is the *second
most efficient* chip of P2_OUT at 86.2% — a screaming channel does not blind
the chip it sits on.

**Not the track sample.** Efficiency against the uRWELL front–back agreement
|dx| runs 0.83–0.88 with no trend, so tightening the track cut recovers nothing.

**Not missing captures.** Capture coverage is 100% of the tracks.

**Not charge sharing at pad boundaries.** Folding every pad onto one cell, the
efficiency is flat across it: 0.85 at the centre, 0.81 at 7.5 mm out.

**Not a common-mode DAQ loss.** P2_MID and P2_OUT miss almost independently —
both fire on 52.1% of tracks against 50.5% predicted by independence.

## What it is

**The discriminator sits at the peak of the signal.**

1. **The pulse height of the hits that are found is a Landau chopped off at its
   low side.** P2_OUT: lowest ADC seen 47, 5th percentile 81, most probable 128.
   The spectrum does not rise from zero, it starts at a wall.

2. **Efficiency is steeply threshold-dependent right there.** Adding 50 ADC
   counts offline costs 14 points (0.854 → 0.715); adding 100 costs 41 points.

3. **Efficiency responds to the amplifier gain, which cannot change a
   detector.** `run_41` (3.0 mV/fC) and `run_48` (4.5 mV/fC), same HV, same gas,
   sdt 230 on five of the six P2_OUT chips in both (the sixth, VMM 19, carries
   1.5% of the tracks), 4½ hours apart: P2_OUT **0.723 → 0.822**.

4. **Efficiency responds to the threshold DAC alone.** `run_47` and `run_48` are
   the *same configuration file* with every per-chip `sdt` line commented out
   or set — identical gain, peaking time, HV, gas, beam, 36 minutes apart:

   | | P2_IN | P2_MID | P2_OUT |
   |---|---|---|---|
   | `_deflt` (base thresholds) | 0.049 | 0.218 | 0.486 |
   | `_opt` (tuned per chip)    | 0.167 | 0.598 | 0.822 |

5. **Within one plane, at one instant, the chips separate by their threshold
   DAC.** P2_MID: VMM 9 and 12 at sdt 256 give 0.53 and 0.31; VMM 10, 11, 13 at
   sdt 224 give 0.78, 0.70, 0.72.

6. **Chips at equal DAC still differ.** All five illuminated P2_OUT chips are at
   sdt 224 and they read 0.897, 0.862, 0.844, 0.829, 0.763 — a 13-point spread
   at one threshold setting. Their pulse-height spectra differ correspondingly
   (medians 119 to 162 ADC), though not cleanly enough chip by chip to call it a
   trend on five points; the pad-level version below is the quantitative
   statement. The VMM's per-channel trim DACs were never used, so one global
   threshold per chip means the effective threshold moves with every channel's
   baseline.

7. **Pad by pad, efficiency tracks pulse height.** Over the 69 illuminated pads
   of P2_OUT the correlation between a pad's median pulse height and its
   efficiency is **+0.55** (P2_MID: +0.52 over 54 pads), and in quintiles:

   | pad median pulse height [ADC] | efficiency | tracks |
   |---|---|---|
   | 90–127  | 0.752 | 151,914 |
   | 127–139 | 0.816 | 241,640 |
   | 139–168 | 0.892 | 163,012 |
   | 168–202 | 0.920 | 174,252 |
   | 202–271 | 0.921 | 106,159 |

   Note the direction. If this were a selection effect it would run the other
   way: a channel with a *high* threshold records only its big pulses and so
   would show a high median with a low efficiency. The observed correlation is
   positive, so it is real signal-size variation, and it dominates.

8. **The best pads already reach DREAM.** Per-pad efficiency on P2_OUT spans
   0.59 (5th percentile) to 0.943 (95th), median 0.894, and the best individual
   pads are **0.962, 0.957, 0.954**. DREAM measured 0.9604. A detector cannot be
   96% efficient on some pads and 60% on others under one uniform beam; a
   readout with one global threshold per chip and no per-channel trim can.

9. **The mesh scan and the amplifier-gain step are one curve.** A Landau
   survival model fitted to the 11-point mesh scan (0.734 at nominal down to
   0.0078 at −100 V) and required to pass through the ×1.5 amplifier-gain point
   returns a gain e-folding of **V0 = 22.4 V** — what a bulk Micromegas actually
   has, and not an input to the fit. The discriminator comes out at 0.68 × the
   most probable signal. A Moyal cannot fit it; the scan decays as a power law
   below −40 V, which is the Landau's 1/x tail.

   | signal over threshold | equivalent mesh V | modelled efficiency |
   |---|---|---|
   | ×1 (now) | 0 | 0.743 |
   | ×1.5 | +9 | 0.822 *(measured: 0.822)* |
   | ×2 | +16 | 0.860 |
   | ×3 | +25 | 0.894 |
   | ×4 | +31 | 0.910 |

   Read the table where the scan measured it; the model's asymptote is an
   artefact of the Landau's unphysical low tail, and what the detector can do at
   this HV is not a model question anyway — DREAM measured it.

---

## Two bugs in my own analysis, found on the way

**The hot-channel mask was deleting beam.** It flagged any channel above 8× its
own chip's median occupancy. On a chip whose threshold has been raised the noise
disappears, the median collapses, and the illuminated pads stand 30–50× above
it — so they get masked. On P2_MID this removed 64,630 hits arriving *exactly at
the +115 ns trigger latency*; the pads it killed measure 61% efficiency with the
mask off and 1% with it on. VMM 9 was not dead, it was masked.

Flagging on the out-of-time occupancy does not fix it either: the beam delivers
particles at tens of kHz while the matched trigger runs at ~1.4 kHz, so most of
a beam pad's hits are out of time with a *matched* trigger. The rule is now an
absolute rate — a 12 × 12 mm pad in this beam takes a few kHz, the pathological
channels are at 10⁵ Hz, and there are two decades of clear air between them.

| station | ratio rule (as published) | absolute rate (corrected) |
|---|---|---|
| P2_IN  | 0.1535 | 0.1552 |
| P2_MID | 0.5916 | **0.6688** |
| P2_OUT | 0.8540 | 0.8540 |

P2_OUT is untouched, so the headline number stands. The rule also still does its
job where it matters: leaving VMM 4 channel 39 (460 kHz) in inflates P2_IN to
0.2253 by pure accidental coincidence, since at that rate the channel fires in
29% of the 630 ns windows all by itself.

**Capture spans started at t = 0.** A few hits per capture carry
`srs_timestamp = 0`, and taking each capture's span as (min, max) of its hit
times stretched every one of them back to zero. Unioned, that made the
"was the VMM recording" cut accept everything, and made the live time come out
25× too long. Both now use quantiles.

---

## The three stations

**P2_OUT — 0.854.** The clean case: no masking effect, no dead pads besides the
one DREAM also sees, uniform in time and rate, and the whole deficit in the
per-pad pulse-height spread.

**P2_MID — 0.669.** The same threshold-limited 0.70–0.78 as P2_OUT's chips on
its three healthy ones, minus two chips at sdt 256 (VMM 9 at 0.53, VMM 12 at
0.31, together 25% of the illuminated area), and its pulses are smaller than
P2_OUT's to begin with (most probable 80 ADC against 128) which puts it further
down the same curve. VMM 12 also has 27 of its 64 channels with zero hits in the
entire sub_run — worth a bench check.

**P2_IN — 0.155. A different problem, and not resolved here.** Three of its six
chips run at raised thresholds (sdt 256/288/300); VMM 4 carries 1.14 × 10⁹ hits
in the sub_run, channel 39 alone at ~460 kHz, and reads 4.9%. But even its two
chips at sdt 224 only reach 27.3% and 21.0%, far below what P2_MID and P2_OUT
manage at the same setting. Its residual is 4.66 / 4.46 mm against 3.37 / 3.44
elsewhere, and its matching-distance histogram has a shoulder at 20–30 mm — the
signature of a partly wrong channel → pad map, on the plane already known to
have mapping trouble (`STRATEGY_OVERRIDES` carries two hand-found exceptions for
P2_IN and none for the other two). Its number is a lower limit on the detector,
not a measurement of it.

---

## What to do about it

1. **Use the VMM's per-channel threshold trim.** Only one global `sdt` per chip
   was ever set. The per-pad efficiency spread of 0.59–0.96 at a single DAC
   value is what that costs.
2. **Get the global thresholds down.** sdt 224–230 is not a low threshold here:
   it is 0.68 × the most probable signal.
3. **Fix the noisy channels at source** (VMM 4 ch 39, VMM 19 ch 20) instead of
   raising the whole chip's threshold to hide them — that trade cost the entire
   chip on P2_IN.
4. **Bench-check VMM 12 of P2_MID** (27 dead channels) and re-derive the P2_IN
   channel → pad map, now that uRWELL tracks make a proper derivation possible.

## Reproducing it

```bash
# one streaming pass over the captures, all three stations, everything kept
python3 eff_autopsy.py run_46 cfg_gain4.5_peaktime200 --out autopsy

# every hypothesis tested against that one event sample
python3 eff_autopsy_report.py autopsy/autopsy_run_46_..._tracks.parquet \
    --thresholds thresholds_run_46.json --mask rate --dream-eff 0.960 \
    --json report.json --figdir figs

# the mesh scan and the gain step as one curve
python3 gain_threshold_model.py --station P2_OUT --observed-at-gain-step 0.822
```
