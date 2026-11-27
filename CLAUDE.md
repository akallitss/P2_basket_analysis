# P2_basket_analysis — project rules for analysis sessions

These are lessons from the det2 cosmic-bench audit of 27 Nov 2026
(`cosmic_bench_analysis/AUDIT_det2_pillar_lattice_20261127.md`). They exist so the
same mistakes are not made again. Read them before any efficiency, pillar or
no-hit study.

## 1. Never trust a drawing for a geometry test: measure it from the data first

- The Saclay-bulked P2 chambers (det1-det4) have their small pillars on a
  **3.0 mm** triangular lattice, not the 4.0 mm of the Gerber
  (`Detector_Drawings/Version_Apr26/Insulation_masks/V1|V2`), and the mask film
  is shifted about (-2, +1) mm in the pad frame. Every bulk sheet and slide says
  "Mask_M2_V2, 4 mm"; the data say otherwise. Until a microscope check is on
  record, the Gerber is wrong for pillar positions.
- Before folding tracks on ANY assumed structure (pillar lattice, mesh, pad
  pattern, mask), run the model-free search: `41_pillar_lattice.py` (2D Fourier
  of the no-hit positions, all wavelengths and directions). A real periodic
  shadow shows up as a 60- or 90-degree set of peaks well above the
  sqrt(2 ln N_trials) noise ceiling. Only then fold.
- A negative result from a folded/offset scan ("no shadow at the Gerber
  positions, nor for any shift within a cell") only excludes THAT lattice. It
  says nothing about a lattice with a different pitch. Say so explicitly in the
  write-up instead of concluding "no pillar effect".
- Pillar lattices are fixed to the board: carry them between runs of the same
  chamber in the PAD frame (`41_pillar_lattice/lattice_pad.npz`), and
  cross-check the origin by the Fourier phase on the new run.

## 2. Cross-check a "mechanism is open" conclusion against the simplest geometry

- The no-hit loss had every signature of a sub-mm periodic dead structure
  (depends on lateral spread 4 mm * tan theta, uniform in position at the 0.7 mm
  scale, present in 87 % of cases as "charge never arrives"). The handoff listed
  pillars as excluded on the strength of the folded test above. When the
  signatures point at one mechanism and one test excludes it, re-examine the
  test's assumptions before looking for exotic mechanisms.
- A fit that prefers an unphysically wide edge (stage 31's big-pillar sigma 2.4
  mm against a 0.7 mm pointing) is a registration problem until proven
  otherwise. Fit each object separately (centre, radius, edge) before fitting
  them jointly at drawn positions.

## 3. The shadow depth is not geometry alone

- The effective dead radius of the pillars moves with the signal margin: 0.70 mm
  at low drift field or low gain, 0.23 mm at 800-1200 V/cm, 0.18-0.33 mm at the
  nominal 430/730 point depending on the run (gas quality, charging state).
  Never quote a pillar loss without the operating point and the run; use
  `42_pillar_shadow_vs_field.py` on the scan series.
- Run-to-run efficiency differences at identical voltages are first a shadow
  depth question, not a position or alignment question.
- Design extrapolations use loss = pi (r_drawn + delta)^2 / cell with delta from
  the measured margin table (audit section 5), never the drawn radius alone and
  never a model fitted on the wrong lattice (old stage 36 numbers are void).

## 4. Alignment and denominators (unchanged rules, still easy to get wrong)

- Pad->M3 alignment: always the trimmed fit (`fit_transform(..., trim_mm=15)`);
  the plain least squares is biased by second particles (scale 0.978, apparent
  pointing 3 mm instead of 0.7 mm). The det2 keys carry it; other keys still
  need switching and a full re-run before numbers are compared side by side.
- Guards: 2 mm at dead tiles and the outer edge with the trimmed alignment.
  The old 6 mm guard was compensating misalignment.
- Quote live-area efficiency AND the excluded fraction (det2: 322 of 1280 pads
  outside the denominator). A fixed-hardware number is a different number.
- Efficiency definition: closest CC-20 cluster within 20 mm of the track, with
  the event-mixing accidental subtraction (0.94 %). Seed thresholds remove the
  soft-edge pillar tracks first (50 ADC costs 0.8 pt) - do not add one.

## 5. Clustering: do not spend more sessions on it

- Tile adjacency and 20 mm links are identical on P2 (12 mm pads, no charge
  spreading); timing does not separate in-time second particles; neighbour-
  seeded thresholds find nothing (0.05 sub-threshold neighbours per hit).
- The only recoverable loss below the hit finder's 5-sigma threshold is a
  matched-filter second pass on the pad under the track: +0.4 pt at 3 sigma,
  +0.28 at 4 sigma, angle- and pillar-independent (Landau low tail). It is
  track-guided, so it is a bench-only number; keep it separate if ever used.

## 6. Bookkeeping

- Write the lattice/alignment/guard assumptions of every number into the
  summary file the stage produces; a handoff that says "no shadow" must say on
  which lattice.
- When a stage's conclusion is reversed, leave the stage, add the correction to
  the handoff banner and the audit, and point the CERN note at a correction
  block rather than silently re-deploying.
- Git: work on `main`, commit at wrap-up, push with
  `git push https://github.com/akallitss/P2_basket_analysis.git main`
  (the `origin` SSH remote has no key here; never the CEA GitLab remotes).
- Python: `P2_basket_analysis/.venv/bin/python`; `run_p2_pipeline.sh` by
  absolute path; `garfield_inputs/*.csv*` are git-ignored.

## 7. Deformation, drift gap and pulse timing (audit of 27 Nov 2026)

Source: `cosmic_bench_analysis/AUDIT_det2_deformation_20261127.md`. Read it before any timing-map,
gain-map or position-dependent efficiency study.

- The drift gap of det2 is **not 4 mm**: 3.3 mm in the centre, 4.0 mm over most of the area, 6-7 mm in a
  ridge along the left side (per-pad table: `<run>/49_deform_decomposition/pads.csv`). Anything that
  depends on the drift field, the primary charge or the lateral spread of a track (pillar shadow,
  efficiency, time resolution, amplitude) changes with it. Quote such numbers per drift-gap class, or say
  that they are chamber averages over an uneven gap.
- A map of peak time, rise time or pulse width is three things at once (drift time, ion tail, delay). Do
  not read it as one of them. Use the quantities that isolate one: the pulse-to-pulse time spread on
  single-pad events and the pad-pair time difference vs the M3 slope (drift only), the mean-pulse cross fit
  (drift if the response is the same), the zero-drift-field point and the early quantiles (delay).
  Stages 47, 48, 44.
- A muon amplitude map is not a gain map. Divide by the primary charge (proportional to the local gap) and
  the collection at the local field first (stage 49 `lnG`). The uncorrected map correlates with the gap
  at r = 0.5.
- Convert a drift time into a gap with the drift velocity measured in the chamber's own drift scan
  (`49_deform_decomposition.py --scan-dir <44 output>`), anchored at the nominal field. The dry Magboltz
  table for Ar/iC4H10 95/5 is right above 600 V/cm and 5-35 % too fast between 500 and 250 V/cm on the bench.
  The bench gas behaves like Ar/iC4H10 with ~0.7 % water (`gas_studies/water_scan/`): slow at low field and,
  18 % faster than dry gas at 750 V/cm. A drift scan alone cannot separate more water from a thinner gap (a
  valley: 0.9 % at 3.5 mm to 0.5 % at 5.0 mm), but a gap cannot replace the water. The drift time at the
  operating point from the pad pairs breaks the tie (det2: 0.55 % water, nominal zone 4.7 mm). Gap RATIOS
  between zones are solid; the millimetres of stage 49 are on the dry-gas scale, best factor 1.17 (1.05-1.25).
- Time spread from all hits contains pad sharing (a pad that sees only part of the gap fires early or late).
  The intrinsic number is on single-pad events.
- To say which surface is displaced, use the M3 track height per gap class with estimators that do not use
  the amplitude ratio of two pads (stage 45). Stage 35's tilt and bow numbers are void (corner-touching pads
  were counted as neighbours).
- Do not re-read `decoded_root` for a pulse-shape question: stage 43 writes every hit and track-near
  waveform of a run once (`43_deformation/wf_table.pkl`, `wf.npy`); `p2_deform.load(run_key)` returns it
  with the shape features.
