# AUDIT — det2 (P2_2) deformation from the cosmic-bench data, 27 Nov 2026

Question asked: the board is deformed in two ways (drift gap, amplification gap) that the hit-level data
could not separate. Can the waveforms isolate the drift gap, and what exactly is the deformation and what
does it do to the signal? Data: `det2_long_repos` (25 Nov, 430/730 V, 7.6 h, waveforms) and the 24 Nov
alignment / drift-scan / mesh-scan series (hits only). New stages 43-49, all outputs in `43_...`-`49_...`
folders of the run. Rules distilled from this audit are in `P2_basket_analysis/CLAUDE.md` section 7.

## 1. The findings in one page

1. **Yes, the waveforms separate the two.** A muon pulse is `box(T_d) (x) ion tail (x) shaper`, with
   `T_d = gap / v_d` the time the charge takes to arrive from the whole drift gap. The mean pulse length
   cannot tell the drift box from the ion tail (convolution commutes), which is why hit-level timing was
   ambiguous. Three things can: (a) the **pulse-to-pulse spread** of the arrival time is pure drift
   (where the primary clusters sit along the gap), the ion tail and the electronics are the same for every
   pulse; (b) pulses of two places with the same response differ only by the box, so the mean pulse of one is the
   mean pulse of the other with a longer box, and a cable delay would shift the start instead (it does not:
   the pulses start at the same time); (c) for an inclined track that crosses two pads, the pad holding the
   upper part of the track fires later by about `T_d / 2` (**micro-TPC**), with the sign given by the M3 slope.
2. **Drift gap.** It is not uniform: about 4.0 mm over most of the chamber (60 % of the pads within
   3.6-4.4 mm), **3.3 mm in the centre** (minimum 3.0 mm, zone centred at pad (378, 279) mm), and **6-7 mm in a
   ridge along the left side** (peak 7.2 mm, 8.2 mm with the dry-gas table; 71 pads above 5.5 mm). Five
   independent measurements agree (section 3). The ratios between zones are the solid result (centre 0.85 x,
   ridge 1.6 x the nominal zone, 1.8 x at its peak). The millimetre scale assumes the dry-gas drift velocity
   at 750 V/cm; if the gas holds 0.7 % of water, as the low-field data suggest (section 8), every gap is
   18 % larger (nominal zone 4.7 mm).
3. **Shape.** The thick zone is a ridge parallel to the left side of the fan, about 400 mm long end to end
   (240 mm above half height), peaking 60 mm inside the edge of the pad array, falling back both towards
   the edge and towards the interior (3.9 mm at 118 mm from the edge), and back to 3.6-3.8 mm about
   200 mm either side of its middle (section 4).
4. **Which surface.** The M3 tracks say the mid-height of the gap does not move where the gap opens:
   -0.14 +- 0.26 mm between the nominal and the thickest gap class (+3.0 mm of gap with the dry-gas
   table, +2.6 mm with the measured velocity). So **both faces move apart**, the window up by
   1.4 +- 0.26 mm and the board (mesh) down by 1.7 +- 0.26 mm (1.2 and 1.4 mm with the measured velocity). A window-only bulge and a board-only sag are each excluded at more than 5 sigma.
   (In the thin centre the same test cannot decide.)
5. **Amplification gap.** After removing what the drift gap does to the amplitude, the intrinsic gain
   still varies by 0.29 in ln G (robust sd), with large-scale structure: low at the top of the fan, high in
   the two lower corners. Where the gain is lower the ion tail is longer (-16 +- 5 ns per e-fold of gain,
   r = -0.5 to -0.6), which is what a larger amplification gap does and what a transparency or electronics
   gain difference would not do. The size matches the mesh-scan sensitivity. If it is all gap: rms 8 um,
   5-95 % range -16 to +12 um around 128 um. July Fe55 data show the same gain map (r = 0.81) and the
   same gain-rise relation on other electronics channels (section 5).
6. **Effects.** Amplitude follows the drift gap (982 / 1104 / 1452 / 1859 ADC for 3.4 / 4.0 / 4.8 / 6.3 mm),
   the time resolution follows the drift time in proportion (19.8 / 21.1 / 24.7 / 32.4 ns), and the lost
   muons fall with the gap: 3.8 / 3.4 / 2.4 / 0.9 % of the tracks without a hit; near a pillar, for vertical
   tracks, 23 / 21 / 15 / 4 % (section 6). The pillar shadow depends on the local drift gap, and at a
   given gap the loss follows the intrinsic gain (4.7 / 3.2 / 2.1 % for low / mid / high gain).
7. **Static.** Same pattern on 24 and 25 Nov (r = 0.97, the chamber was moved in between) and constant
   over the 7.6 h of the long run. No electronics delay is involved (zero-drift-field offsets 4.6 ns rms).

## 2. Why the hit-level data were ambiguous, and what the waveforms add

For a pulse on one pad, `SER` being the response to one electron (ion drift across the amplification gap,
then the DREAM shaper):

| quantity | drift gap +1 ns of T_d | ion tail +1 ns | cable delay +1 ns |
|---|---|---|---|
| peak time | +0.47 ns (drift scan) | about +0.5 ns | +1 |
| 30 % leading-edge time | +0.35 ns (drift scan) | small | +1 |
| rise (peak - 30 %) | +0.10 ns (drift scan) | large | 0 |
| pulse-to-pulse spread s68 | **proportional to T_d** | **0** | **0** |
| amplitude | more primary charge | gain falls if the gap is the cause | 0 |
| time difference of a pad pair under an inclined track | **T_d / 2** | 0 | cancels in the antisymmetrisation |

The peak-time map alone (stage 19/23) mixes all three columns. The forward model
(`p2_deform_model.py`: HEED clusters, Magboltz drift, ion box, shaper, 60 ns sampling) gives the same
pattern and was only used for response factors that the data could not give directly.

Model-free evidence that the late band of the stage-19 map is drift and not a delay:

- Quantiles of the 10 % leading-edge time, band (connector 9, y 150-350 mm) against connectors 2-7:
  1 % +8 ns, 16 % +13, 50 % +18, 84 % +34, 99 % +67 ns. The earliest pulses arrive at the same time; the
  distribution is stretched from a fixed start.
- The mean pulses (sum of waveforms on a 10 ns grid, using the trigger fine time) start together; the band's
  pulse is the reference pulse convolved with a longer box: fit of `box(T_ref) (x) M_band = box(T_band) (x)
  M_ref` with everything free gives T_ref = 125 ns, T_band = 189 ns, extra delay -1.4 ns, residual 0.09 %
  of the peak. 125 ns is what Magboltz gives for 4.0 mm at 750 V/cm (128 ns).
- On single-pad events the spread is proportional to the drift time (gap classes: spread ratios 0.94 / 1 /
  1.17 / 1.54 for T_d ratios 0.89 / 1 / 1.16 / 1.48).

## 3. The drift gap: five measurements

| # | method | data | what it needs | thin centre | nominal | thick ridge |
|---|---|---|---|---|---|---|
| 1 | mean-pulse fit against a reference pulse (stage 47) | 25 Nov waveforms | same SER in both places | 107-120 ns | 125 ns (ref) | 170-184 ns (region), 203 ns (ridge pads) |
| 2 | pulse-to-pulse spread, single-pad events (47) | 25 Nov waveforms | nothing about the SER | 113-123 ns | 125 ns (ref) | 159-182 ns (region), 212 ns (ridge pads) |
| 3 | pad-pair time difference, micro-TPC (48) | 25 Nov waveforms + M3 slopes | diffusion/threshold factor from the model | 113 +- 7 ns | **129 +- 3 ns (absolute)** | 174 +- 6 ns |
| 4 | regional peak time vs drift field, relative to the zero-field time (44) | 24 Nov hits | v_d(E) | 3.3 mm | 4.0 mm (assumed mean) | 5.8-6.1 mm |
| 5 | stretch of the amplitude-vs-drift-voltage curve (44) | 24 Nov hits | none (field scaling) | 0.81-0.88 x | 1 | 1.35-1.52 x |

Checks that the long drift time is a thicker gap and not slower drift in the same gap: the amplitude is
1.68 x larger and the Landau narrower in the ridge (ln-width 0.79 against 0.96); the amplitude curve is
stretched in voltage (method 5: the local field is lower by the same factor); inclined tracks (tan theta
0.24-0.32) fire two pads in 19.4 +- 1.5 % of the cases against 8.1 +- 0.4 % (forward model: 14.6 % for
7 mm, 7.8 % for 4 mm, 7.7 % for slow drift in 4 mm).

Final per-pad gap (stage 49): fine structure from method 1, large-scale level tied to method 2 (so that the
ion tail cannot leak in), converted with **the drift velocity the chamber itself shows** in the drift scan
(nominal-gap regions, ratio to the nominal field; only the anchor at 750 V/cm is from Magboltz and is
confirmed by method 3). With it the 24-region comparison with method 4 is 1 : 1 (slope 1.02, offset
-0.06 mm, rms 0.28 mm, r = 0.92) and with method 5 slope 0.89, rms 0.27 mm.

Measured v_d / dry Magboltz [um/ns]: 26.1 / 40.6 at 250 V/cm, 30.8 / 39.8 at 312, 33.8 / 38.2 at 375,
33.5 / 35.2 at 500, 32.7 / 32.9 at 625, 31.2 / 31.2 at 750 (anchor), 29.3 / 30.0 at 875, 27.5 / 27.9 at 1312.
The dry-gas table is right above 600 V/cm and too fast below (section 8); with it the ridge would come out
6.7 mm (class mean) and 8.2 mm (peak) instead of 6.3 and 7.2 mm.

Uncertainties: statistical +-0.1 mm per pad (sliding 18 mm). Drift-time scale: the nominal zone is 125 ns
(mean-pulse fit), 129 +- 3 ns (micro-TPC with the model's dilution factor; the raw plateau alone gives
>= 112 ns). Millimetre scale: drift time x v_d(750 V/cm), and v_d there is 31.2 um/ns in dry gas but
36.7 um/ns with 0.7 % water (section 8), so the common scale factor of all gaps is 1.0 (dry, nominal zone
3.96 mm, the design value) to 1.18 (wet, 4.7 mm). The data do not choose: the drift-scan fit trades gap
against water fraction (3.3 mm at 1.0 % to 4.8 mm at 0.5 %). Ratios between zones are not affected. In the
ridge the three time estimates differ by +-4 % (174 / 186 / 189 ns), i.e. +-0.3 mm. Lowering the amplitude
cut of the timing sample from 300 to 150 ADC moves the regional drift times by < 0.5 ns (mean-pulse fit)
and <= 3 ns (spread).

## 4. Shape and surfaces

- Ridge (stage 49, `ridge_profile.png`): axis at 61 deg in the pad frame through pad (222, 257) mm, i.e.
  parallel to the left side of the fan (connectors 9-10). Across it (central 140 mm): 3.9 mm at 117 mm from
  the edge of the pad array, 4.2 at 105, 4.7 at 94, 5.4 at 82, 6.3 at 71, **6.8 at 60**, 6.7 at 48, 6.5 at 36,
  5.8 mm at 23 mm (last live pads; connector 10 beyond is unplugged, so the last 23 mm are not seen; the
  single-pad drift times 203 -> 186 -> 168 ns confirm that it falls towards the edge). Along it: 3.8 / 4.3 /
  5.2 / 5.5 / 6.4 / 6.6 / 6.9 / 6.9 / 6.8 / 6.3 / 5.9 / 5.2 / 4.5 / 3.7 / 3.6 mm in 30 mm steps, i.e. an arch
  that is back to nominal about 200 mm either side of its middle (the low-radius corner beyond it is at
  4.2 mm).
- Thin zone: 130 pads more than 0.5 mm below the median, centred at (378, 279) mm.
- Orientation: the pad holding the higher part of the track fires later (D > 0 for every slope bin), so
  the window is on the upper side and the mesh below (M3 z up).
- Surfaces (stage 45, pad-level fit of the track height against the gap; slope -0.5 = only the board moves,
  +0.5 = only the window moves): boundary step on single-pad events -0.04 +- 0.10, centroid of single-pad
  clusters -0.02 +- 0.04 (response 0.54 in the toy), centroid x and y -0.05 +- 0.05 (gap = the preliminary
  dry-gas map of stage 47; with the measured velocity the slopes scale by 1.15). Delta-chi2 for
  -0.5 / 0 / +0.5: 19 / 0.2 / 27, 120 / 0.2 / 138, 73 / 0.8 / 106. Stable against a tilt (changes by <= 0.09)
  and against dropping the two apex regions.
  Estimators that use the amplitude ratio of two pads sit about 2 mm lower and are not used: the charge
  from the upper part of the gap diffuses more and leaks to the neighbour, so the mesh-side pad gets the larger
  share (in pad pairs with |slope| > 0.15 the larger pulse is on the mesh side in 57 % of the cases).
- The run_config describes this chamber as "bulked 25-6-26 with misaligned wall". Whether that is related
  is for the bench: the data only say that along the middle of the left side the two faces stand 3 mm
  further apart than they should, symmetrically, and that this does not change with time or handling.

## 5. The amplification gap

Amplitude decomposition per pad (stage 49): `A = primary charge (gap) x collection(E_local) x ballistic(T_d)
x G`. Collection vs field is the chamber amplitude curve of the drift scan (F(450)/F(750) = 0.98,
F(900) = 0.96, F(1200) = 0.84). After the correction ln G does not depend on the gap (r = +0.06; ln A
vs ln gap was +0.50) and the class means are flat (-0.05 / -0.01 / +0.10 / +0.05).

- Spread of ln A 0.39, of ln G 0.29 (robust sd, sliding 18 mm). Same map on the other day from the
  amplitude-curve collapse (g/s): r = 0.91.
- Ion tail: S = T_fit - T_spread (the part of the mean-pulse length that is not drift), 5-95 % -8 to
  +12 ns. S against ln G: -16 +- 5 ns per e-fold (48 regions, r = -0.48), -19 within connectors, -17 with the
  other day's gain (24 regions, r = -0.56), pad level r = -0.56. With 24 regions of this run alone it is
  -11 +- 7.
- Expected if the gain variation is the amplification gap: Magboltz gives d ln G / d gap = -0.036 per um at
  430 V (1 e-fold = 27 um = 21 % of the gap, so +42 % of ion transit time). The mesh scan gives the pulse
  response to the ion transit: the peak moves by -0.4 to -1.3 ns per 10 V at fixed amplitude (2.6 % of
  transit time), i.e. 7-20 ns of peak shift, 15-40 ns of box-equivalent, per e-fold. The measured 16 ns is
  inside. So the intrinsic gain variation is compatible with being all amplification gap.
- Fe55 (July, P2_OUT position, FEU 3, other channels; stage 46 + this audit): gain map vs November
  intrinsic gain r = 0.81 (sliding 18 mm, 469 pads of connectors 4-7); Fe55 rise vs Fe55 gain -12 to -15 ns
  per e-fold (point-like charge, no drift box); November muon rise with the drift part removed vs July
  Fe55 rise r = 0.72-0.74; ln(muon / Fe55) vs ln(gap): slope 1.16, r = 0.60 (primary charge proportional to
  the gap, with an independent gain).
  Not understood: the November gain contrast is 1.5 x the July one (slope of ln G vs ln G_Fe55), and in
  July 70-89 % of the photopeak events have one neighbour pad at 44-60 % of the amplitude.
- Electronics-sized terms: per-connector medians of S -5 to +7 ns, of the fitted delay -3 to +3 ns; the
  fitted delay map is smooth (-6 to +5 ns across the chamber), not connector-shaped (trigger light path is
  the candidate; not investigated).

## 6. What it does to the signal (stage 49, `effects.png`, `classes.csv`)

| drift-gap class | pads | gap | T_d | local field | time spread s68 | amplitude (median) | below 300 ADC | lost muons, all / 0-8 / 8-14 / > 14 deg | near a pillar, 0-8 deg |
|---|---|---|---|---|---|---|---|---|---|
| < 3.6 mm | 195 | 3.36 mm | 114 ns | 890 V/cm | 19.8 ns | 982 | 12.6 % | 3.78 / 5.2 / 3.7 / 2.4 % | 22.5 +- 2.7 % |
| 3.6-4.4 | 577 | 3.96 | 128 | 760 | 21.1 | 1104 | 10.9 % | 3.35 / 5.2 / 3.3 / 1.7 % | 21.3 +- 1.6 % |
| 4.4-5.5 | 112 | 4.81 | 148 | 630 | 24.7 | 1452 | 7.1 % | 2.42 / 4.4 / 1.8 / 1.2 % | 15.3 +- 3.0 % |
| > 5.5 | 71 | 6.32 | 189 | 480 | 32.4 | 1859 | 5.3 % | 0.88 / 2.0 / 0.8 / 0.1 % | 3.9 +- 1.9 % |

Far from the pillars (> 1.25 mm) the no-hit rate is 1.3-2.1 % in every class. Lattice: the measured 3.0 mm
one (stage 41, pad frame). Tracks at least 2 mm inside the live area, trimmed alignment, no-hit = no P2 hit
in the event. Time spread: single-pad pulses, 20 % leading edge, pad medians removed, includes the trigger.

- Both components matter for the lost muons. Logistic fit on 40,513 tracks: log-odds of no hit =
  a - (1.75 +- 0.20) ln(gap / 4 mm) - (1.13 +- 0.10) ln(G / median) - (4.8 +- 0.3) tan(theta). One e-fold of
  intrinsic gain multiplies the odds of losing a muon by 0.32, a gap 1.5 x thicker by 0.49. Within the
  nominal-gap class, by gain tercile (ln G = -0.36 / -0.04 / +0.30): 4.7 / 3.2 / 2.1 % lost (6.3 / 5.6 /
  3.7 % below 8 deg). Near a pillar (< 0.6 mm) the dependences are steeper (gap -2.4 +- 0.3, gain
  -1.3 +- 0.2); far from the pillars (> 1.2 mm) the gain term is the same (-1.1 +- 0.2) and the gap term
  weaker (-1.0 +- 0.4).
- The pillar shadow is set by the local drift gap: a thicker gap gives more primary charge, more lateral
  spread of an inclined track (gap x tan theta) and more transverse diffusion. This is the position
  dependence of the shadow depth that the 27 Nov audit left open.
- For the design: a 6.3 mm gap loses 0.9 % of the muons where a 4.0 mm gap loses 3.4 %, and pays for it
  with 32 ns instead of 21 ns of time spread.

## 7. Electronics, time and day stability

- Zero-drift-field point (dV = 0.5 V): regional time offsets 4.6 ns rms, within their errors; at
  625-1312 V/cm the regional offsets are 9.8 ns rms and scale with the drift time (stage 44).
- 24 Nov vs 25 Nov, 24 regions: peak-time offsets r = 0.97, slope 1.03; spreads r = 0.97; amplitude ratio
  0.90-1.06.
- Long run in four slices of 1.9 h: T_d of the ridge 190 / 191 / 195 / 191 ns, of the thin zone
  105.7 / 105.8 / 106.6 / 106.5 ns.

## 8. Drift velocity at low field

What the data say (stage 44, hits of the 24 Nov drift scan, mesh 430 V):

- Above about 500 V/cm the chamber follows dry-gas Magboltz for a 4 mm gap: (median peak time - zero-field
  time) / T_gas is 0.49 at 500 V/cm and 0.48 at 1312 V/cm.
- Below it the electrons are slower than dry gas, more and more so: in the nominal-gap regions v_d is 33.8
  instead of 38.2 um/ns at 375 V/cm, 30.8 instead of 39.8 at 312, 26.1 instead of 40.6 at 250; chamber-wide
  the implied drift time is about 370 ns at 125 V/cm and 780 ns at 62 V/cm against 108 and 174 ns. The
  velocity maximum sits near 400-500 V/cm instead of 200. A constant voltage offset does not describe it
  (the implied offset grows from 19 V at dV = 25 V to 54 V at 75 V).
- The same behaviour was seen on det1 in July and at the SPS in another gas (`RESULTS_gap_timing.md`, part B),
  so it is not specific to det2 or to this week.

Consequence for this audit: none on the method, a 7 % smaller ridge than the dry-gas table would give
(section 3), because stage 49 uses the measured velocity.

Cause (Magboltz scan, `gas_studies/water_scan/`, Garfield++ / Magboltz 11.19 on lxplus, 21 gases):

| gas model, chamber, d = 4 mm | chi2 / ndf (with 2 ns per-point syst.) | verdict |
|---|---|---|
| dry Ar/iC4H10 95/5 | 9223 / 21 | excluded, also misses the high-field slope |
| + water, fraction fitted | **0.70 %**: 167 / 20 (17.7 / 20) | fits the whole scan, 62-1312 V/cm |
| + air (N2/O2) | 7540 / 20 (3339 / 20) | excluded |
| more isobutane (up to 20 %) | 2308 / 20 (353 / 20) | disfavoured |

- Water at 0.7 % (0.5-1.2 % over data sets and variants) reproduces the drift-time curve with one gas and one
  gap; the nominal-gap regions alone give 0.78 %. With the gap free the fit trades it against the water
  fraction (3.3 mm at 1.0 %, 4.8 mm at 0.5 % with a free response scale), so the scan does not fix the
  absolute gap.
- Water moves the velocity maximum to ~450 V/cm and makes the gas 18 % faster at 750 V/cm (36.7 against
  31.1 um/ns). That is the scale uncertainty of section 3.
- The amplitude drop at low field (0.27 / 0.57 / 0.75 of the 750 V/cm value at 62 / 125 / 187 V/cm) is NOT
  attachment: explaining it needs eta = 5.9 / 2.2 / 1.3 per cm, which would remove the charge from the
  cathode side and bring the median peak time at 62 V/cm from 365 ns (measured) to 177 ns. The timing allows
  at most a quarter of that (about 0.3 % air). The longer pulses account for 0.64 / 0.86 / 0.96; the rest is
  a depth-independent loss, not understood (the low-field points were taken first, right after the mesh
  scan: repeating 62 and 125 V/cm at the end of a scan would test a gain transient).
- Caveats: Magboltz's water cross-sections are good to 10-20 % (maps one to one onto the fraction); the
  response model has one pad and no track angle; the zero-drift time wants a +6 to +17 ns shift. Not shown:
  where the water comes from. Permeation through the 50 um window gives the right order at low flow; a
  hygrometer on the exhaust and a flow scan (fraction proportional to 1 / flow) would settle it.

## 9. Other chambers

Stages 43, 47, 48, 49 were also run on `det1_long5` (P2_1, Saclay bulk, July, 415/615 V) and `det5_long1`
(P2_5, CERN bulk, September, 420/820 V), without a drift scan (dry-gas table, no collection correction).
They are a check that the method tells the two components apart, not a finished study of those chambers.

| | det2 (this audit) | det1 | det5 |
|---|---|---|---|
| drift time from the spread, 24 regions | 113-182 ns | 107-125 ns | 124-155 ns |
| drift time from the mean-pulse fit, 24 regions | 107-184 ns | 77-146 ns | 100-176 ns |
| micro-TPC, nominal zone (dry Magboltz for 4 mm) | 129 +- 3 (128) | 116 +- 4 (114) | see below (137) |
| pad-pair time difference at large slope | 55-59 ns | 45-50 ns | 80-86 ns |
| ion-tail term S against ln G | -16 +- 5 ns per e-fold, r = -0.48 | -35 +- 5, r = -0.71 | not usable (see below) |
| spread of intrinsic ln G | 0.29 | 0.40 | 0.18 |

- **det1 is the opposite case of det2.** Its drift time is nearly flat (spread 107-125 ns, gap 3.5-4.7 mm at
  5-95 %, and part of that follows the gain, so it is an upper limit), while its mean-pulse length varies by
  a factor 1.9 and follows the gain. So the det1 timing map is the amplification gap (rms 12 um if all
  gap), as the October study (`RESULTS_gap_timing.md`) suspected from the Fe55 correlation. For det1 the
  gap map needs `--s-radius 20`: the gain varies on a shorter scale than the default 45 mm and the
  mean-pulse fit follows the ion tail (it would say 2.5-5.7 mm).
  Not understood: in det1's highest-gain zone the micro-TPC gives 137 +- 6 ns (chi2 17/4) where the spread
  gives 108 ns.
- **det5 needs its own study.** Its gain is low (median 350 ADC, 40 % of the hits below the 300 ADC cut of
  the timing sample), so the mean pulse and the spread are selected on amplitude and the ion-tail term is
  not usable. The spread gives a dome of +-12 % (124-156 ns), smaller than the peak-time dome of stage 23
  suggested. The pad-pair time difference saturates at 80-86 ns against 55-59 ns in det2, which is a
  model-light number: a drift time of at least 165 ns, more likely 190-250 ns (the model fit gives 240-255 ns with a poor chi2), against 137 ns for 4 mm at
  1000 V/cm. Either det5's gap is 5-6.5 mm or its gas was different in September. Worth a check against
  the drawings before any det5 timing or drift-field number is used.


## 10. Open points

1. The mechanical cause of the ridge. Bench checks that would settle it: a straight edge or dial gauge
   across the left side at mid-length on both faces (expected: window +1.5 mm, back -1.8 mm relative to the
   centre, over a 100 mm wide strip); the same with the gas flow off and at a higher bubbler level (a
   pressure effect would change the ridge; the drift-time map of a 1 h run is enough to see 0.2 mm).
2. In the thin centre the window-or-board question is open (-0.6 mm of gap, +0.07 +- 0.25 mm of mid-height).
3. The last 23 mm before the left edge (connector 10) and connector 1 are not read out.
4. Gain contrast 1.5 x larger in November than in the July Fe55 map; the July neighbour-pad signal.
5. The smooth +-6 ns fitted-delay map (trigger light path?).
6. Amplification-gap numbers assume all the intrinsic gain variation is gap; the ion-tail slope supports it
   within a factor 1.5 but does not prove it. A Garfield++ ion-transit calculation for 128 um at 430 V
   would turn the 16 ns per e-fold into a number to compare.
7. The millimetre scale of the gap map: factor 1.0 (dry gas at 750 V/cm) to 1.18 (0.7 % water). A hygrometer
   on the exhaust, or one mechanical measurement of the gap anywhere in the live area, fixes it.
8. The amplitude loss at low drift field (not attachment, section 8).

## 11. Code and products

| stage | what | time |
|---|---|---|
| `43_deformation_waveforms.py` | one pass over the decoded files: waveform of every hit pad and of every pad within 22 mm of a track (`wf_table.pkl`, `wf.npy`, `regions.csv`) | 15-30 min |
| `p2_deform.py` | loader and pulse-shape features | |
| `p2_deform_model.py` | forward model (HEED, Magboltz, ion box, shaper, sampling; pad pairs) | |
| `44_deform_scan_regions.py` | hit-level drift/mesh scan by region: timing vs field, zero-field offsets, amplitude-curve collapse | 1 min |
| `45_deform_local_z.py` | track height (mid-gap) per region and per gap class, centroid and pad-boundary methods, toys | 4 min |
| `46_deform_fe55_vs_mip.py` | July Fe55 per pad vs November muons | few min |
| `47_deform_pulse_maps.py` | drift time per pad: mean-pulse fit, spread, variance | 6 min |
| `48_deform_microtpc.py` | pad-pair time difference vs slope, with the model | 2 min |
| `49_deform_decomposition.py` | gap, ion tail, intrinsic gain, amplification gap, surfaces, effects; needs `--scan-dir <44 output>` for the measured v_d and the collection curve | 1 min |
| `gas_studies/water_scan/` | Magboltz with water / air / other isobutane fractions against the drift scan | |

Order: 34 (alignment) -> 43 -> 47 -> 48 -> 45 -> 44 (scan run) -> 49. Outputs: `<run>/4x_.../summary.txt`
first; the per-pad table is `49_deform_decomposition/pads.csv` (gap, T_d, S, lnG, dgap_amp_um, f_prim,
f_field, f_ball, time spread, amplitude), figures `deformation_maps.png`, `consistency.png`, `effects.png`,
`ridge_profile.png`.
