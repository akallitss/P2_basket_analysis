# AUDIT — det2 cosmic-bench analysis, 27 Nov 2026 (night): the lost muons are the pillars, on a 3 mm lattice

**Read this before the evening handoff.** It changes the conclusions of §3 (pillars "on hold", "no shadow")
and §4 (no-hit mechanism "open") of `HANDOFF_det2_20261127_evening.md`, and the CERN note §12.

## 1. The finding in one paragraph

The angle-dependent "no-hit" loss (5 % for vertical tracks, 1 % above 16°, falling with drift field and
with gain) is the shadow of the small insulation pillars. It was invisible to every direct test because the
tests used the Gerber lattice (0.8 mm pillars on a **4.0 mm** triangular lattice, Saclay masks V1/V2), while the
pillars physically on the board sit on a **3.0 mm** triangular lattice (row spacing 2.600 mm) with the same
orientation and ~2 mm shifted. A model-free 2D Fourier search of the vertical-track no-hit positions finds that
lattice at ~10 sigma in each of its three directions; folding the tracks onto it gives a 27 % no-hit rate within
0.25 mm of a lattice point for vertical tracks (1.5 % far away), shrinking with angle exactly as the angle trend
does. The same lattice, same phase, is found independently on the 11-24 det2 series and on det1's July data
(same Saclay film). The shadow depth is not fixed geometry: its effective dead radius goes from 0.70 mm at low
drift field or low gain to 0.23 mm at 800-1200 V/cm, and 0.18-0.33 mm at the nominal 430/730 point depending
on the run. The 1.7-2.2 points that stage 27 booked to "small pillars" from the angle model were right after
all; the 0.02 % of stage 36 and the "no shadow" of 38/`nohit_checks` were artefacts of the wrong lattice.

## 2. Evidence (det2_long_repos unless stated; outputs in `41_pillar_lattice/`, `42_pillar_shadow_vs_field_det2_initial2/`)

**Periodicity search** (11,812 vertical interior tracks, 544 no-hit). `F(k) = Σ_nohit e^{ik·x} − p0 Σ_all e^{ik·x}`
over 1.5-6 mm wavelengths and all directions (40,680 trials, noise max ≈ 4.6 σ):

| wavelength | direction | significance |
|---|---|---|
| 2.60 mm | 58.0° | 9.8 σ |
| 2.60 mm | 177.9° | 10.1 σ |
| 2.60 mm | 118.0° | 9.0 σ |
| 1.50 mm | 28 / 88 / 148° | 3.2-3.5 σ (second order) |
| 3.46 mm (Gerber row spacing) | 58 / 118 / 178° | < 1.2 σ (absent) |

Three directions 60° apart with the same wavelength = triangular lattice, pitch 2.60·2/√3 = **3.003 mm**.
The Gerber lattice (4.00 mm pitch, rows at 25/85/150°) has the same orientation (rows at 28/88/148°) but no
power in the data. Lattice origin refitted per quadrant: shifts ≤ 0.1 mm → one rigid lattice over the whole
detector. Hex cell 7.81 mm², so 0.8 mm pillars cover 6.45 % geometrically (3.6 % was designed).

**Cross-checks.** (a) 11-24 series, lattice carried over in the pad frame (it is fixed to the board): blind
re-search finds 2.59-2.60 mm at 29.5/89.5°; phase at the det2_long_repos origin −6/−14° (≤ 0.1 mm).
(b) det1_long5 (P2_1, July 2026, same Saclay film): 8.0 and 7.5 σ at 2.601/2.603 mm in the same directions,
origin within 0.4 mm of det2's; folded profile 12.6 % near vs 2.1 % far for vertical tracks. (c) det5_long1
(CERN bulk, 0.5 mm on 2 mm): no periodic shadow above 3.3 σ.

**Shadow profile, no-hit [%] vs distance to the nearest lattice point** (interior, ≥ 2 mm inside the live area):

| angle | 0-0.25 | 0.25-0.5 | 0.5-0.75 | 0.75-1.0 | 1.0-1.25 | 1.25-1.5 | 1.5-1.75 | n |
|---|---|---|---|---|---|---|---|---|
| 0-4° | 26.6 | 18.8 | 8.0 | 5.0 | 2.7 | 1.7 | 1.3 | 3,269 |
| 4-8° | 24.3 | 16.1 | 6.5 | 3.7 | 2.4 | 1.7 | 1.4 | 8,543 |
| 8-12° | 16.4 | 8.6 | 5.2 | 2.5 | 1.4 | 1.4 | 1.4 | 10,916 |
| 12-16° | 6.2 | 7.3 | 3.1 | 1.9 | 1.3 | 1.2 | 0.5 | 9,775 |
| 16-20° | 1.3 | 1.9 | 1.6 | 1.8 | 1.1 | 1.1 | 0.5 | 5,915 |

Hard dead disc ⊗ Gaussian pointing fitted to the vertical profile: r_eff 0.30 mm, pointing 0.40 mm, baseline
1.5 %, χ² 7.2/4. Vertical no-hit 4.82 % = 1.5 % floor + 3.3 points of shadow. Integrated over this run's angle
mix the shadow is **1.66 points** (model-free, far-baseline per class) to **2.2 points** (stage 36's electron
model refitted on the measured lattice; stage 27's angle method gave 2.18). The pointing from the shadow edge
(0.40 mm) is even better than the 0.72 mm of stage 34.

**The edge is soft.** Matched-cluster charge of HIT vertical tracks: median 600 ADC within 0.5 mm of a lattice
point, 1,115 at 0.75-1 mm, 1,410 beyond 1.25 mm; fraction below 300 ADC 33 % → 4 %. Tracks near a pillar lose
part of their charge, so whether they are lost depends on the margin above threshold, i.e. on gain and
collection.

**Shadow vs operating point** (11-24 series, vertical tracks, hard-disc r_eff with the 0.72 mm pointing):

| drift field [V/cm] (mesh 430) | < 200 | 200-400 | 400-600 | 600-800 | 800-1000 | 1000-1200 | 1200-1400 |
|---|---|---|---|---|---|---|---|
| no-hit near / far [%] | 71 / 32 | 38 / 6 | 18 / 4.5 | 18 / 1.8 | 9.0 / 0.5 | 7.5 / 1.5 | 19 / 1.2 |
| r_eff [mm] | 0.70 | 0.55 | 0.36 | 0.32 | **0.23** | 0.24 | 0.28 |
| vertical shadow [pts] | 9.4 | 7.2 | 2.8 | 2.8 | 1.4 | 1.3 | 2.5 |

| mesh [V] (drift − mesh = 300 V) | 340 | 350 | 360 | 370 | 380 | 390 | 400 | 410 | 420 |
|---|---|---|---|---|---|---|---|---|---|
| r_eff [mm] | 0.70 | 0.60 | 0.60 | 0.59 | 0.46 | 0.40 | 0.30 | 0.0 | 0.34 |

At 430/730: 11-24 alignment run r_eff 0.18 mm (shadow 1.6 pts), 11-25 long run 0.33 mm (3.3 pts). Same voltages,
different day: the run-to-run difference noted in the handoff is the shadow depth changing (gas quality /
charging state), not the geometry. The vertical excess re-rising above 1200 V/cm is where the field ratio
(430 V / 128 µm ÷ 1.3 kV/cm ≈ 26) starts to cost mesh transparency.

**Big pillars** (per-pillar disc fits, θ < 12°): centres offset from the Gerber by (−2.1 ± 0.4, +1.2 ± 0.4) mm
(individual: −3.1/+0.7, −1.7/+0.9, −1.6/+1.0, −1.8/+2.2), R_eff 3.45 / 4.48 / 4.08 / 4.67 mm against 3.08 drawn,
edge sigma 0.25-1.2 mm. The mask film is shifted ~2.4 mm on the board, and the big pillars carry a ~1 mm dead
ring. Stage 31's floating σ = 2.4 mm was this offset, not a gradual gain zone; the pointing is sub-mm.

## 3. What this changes in the earlier conclusions

| earlier statement | status |
|---|---|
| "Direct pillar tests: no shadow" (36, 38, nohit_checks) | wrong lattice (4 mm Gerber). On the real 3 mm lattice the shadow is 27 % at the pillar centre. |
| "Small-pillar loss 0.02 %, < 0.23 % at 95 % CL; 0.5 mm pillars gain ≤ 0.2 pt" (36) | invalid. Loss on this run 1.7-2.2 points; see §5 for the design numbers. |
| "27/31 book 2.18 % to small pillars — this is wrong" | 27/31 were right. Keep the row. |
| "No-hit mechanism open: lateral spread 4 mm·tanθ, sub-0.7 mm scale, drift-field dependent" | all three properties are the pillar shadow with a margin-dependent dead radius. |
| "Big-pillar fit floats σ ≈ 2.4 mm, maybe a low-gain zone" | mask film offset (−2.1, +1.2) mm; per-pillar edge is sharp. |
| "Run-to-run excess 1.6 vs 4 pts at the same voltages: position / mesh leak?" | shadow depth 0.18 → 0.33 mm; margin (gas, charging), not position. |
| CERN note §12 "pillars re-examined, no shadow"; pillars "on hold" | needs a correction block; the earlier §7/§8 pillar text was closer to the truth. |
| 13 % of no-hits have a 3-6 σ sub-threshold pulse, "angle independent" | consistent: tracks at the soft edge of a pillar. |

Unchanged and still valid: the trimmed alignment (z 246, scale 1.000), the 2 mm guards, the 0.7 mm pointing,
the accidental-corrected 95.9 % live-area efficiency, the hardware-defect list, the second-particle explanation
of "wrong cluster won".

## 4. Revised loss budget (det2_long_repos, 40,530 active tracks, live area only)

| loss | points | basis |
|---|---|---|
| small pillars (3 mm lattice, soft edge) | 1.7-2.2 | §2; 27's angle method 2.18 |
| big pillars (R_eff ≈ 4.2 mm, 4 inside the live area) | ≈ 0.25 | per-pillar fits |
| angle-independent no-hit floor (sub-threshold tail, far from pillars, > 16°) | ≈ 1.0 | far-baseline at large angle 0.3-1.0 %, 0.5 % floor after accidentals |
| distant / mis-matched hits (reference tails, second particles) | ≈ 0.6-0.9 | 27/33 (unchanged) |
| **total** | **4.1** | 100 − 95.9 |

Reminder for "detector performance": 322 of 1,280 pads are outside the denominator (connectors 1 and 10, the
dead pin row, 13 floating pads, shorted groups). The 95.9 % is a live-area number; the hardware list in the
morning handoff is what brings the whole fan to it.

## 5. Design: pillar diameter, pitch, operating point

Measured effective dead radius r_eff = r_drawn + δ, with δ set by the signal margin:

| regime | δ [mm] | r_eff for 0.8 mm pillars |
|---|---|---|
| good margin: 900-1200 V/cm, fresh gas (11-24 alignment run) | −0.25 … −0.17 | 0.15-0.23 |
| nominal 430/730 after hours of running (this long run) | −0.07 … −0.10 | 0.30-0.33 |
| poor margin: < 400 V/cm, or mesh ≤ 370 V, or 1.3 kV/cm | +0.15 … +0.30 | 0.55-0.70 |

Loss at normal incidence = π r_eff² / cell (triangular cell √3/2·p²); for this bench's cosmic angle mix multiply
by ≈ 0.5 (measured 1.66 / 3.30). Loss [%] at normal incidence:

| Ø [mm] \ pitch | 2.0 mm, good / nominal / poor | 3.0 mm, good / nominal / poor | 4.0 mm, good / nominal / poor |
|---|---|---|---|
| 0.3 | 0 / 0.3 / 14.5 | 0 / 0.15 / 6.5 | 0 / 0.08 / 3.6 |
| 0.4 | 0 / 1.1 / 18.4 | 0 / 0.5 / 8.2 | 0 / 0.3 / 4.6 |
| 0.5 | 0 / 2.3 / 22.7 | 0 / 1.0 / 10.1 | 0 / 0.6 / 5.7 |
| 0.6 | 0.2 / 4.0 / 27.4 | 0.1 / 1.8 / 12.2 | 0.06 / 1.0 / 6.9 |
| 0.8 (now) | 2.0 / 8.7 / 38.3 | 0.9 / 3.9 / 17.0 | 0.5 / 2.2 / 9.6 |

(δ held constant with diameter; this is the only assumption, and it is the conservative one for small pillars
because the dead ring is a field/charging effect at the pillar edge, not a fraction of the pillar.)

What it means:
1. **Operating point first.** On the present det2 (0.8 mm, 3 mm) the shadow is 3.3 points for vertical tracks at
   750 V/cm and 1.4 points at 800-1000 V/cm; drift 750-830 V at mesh 430 is worth ~1 point on the bench and
   ~2 points at normal incidence. The mesh scan says the same about gain: every 10 V of mesh below 420 V widens
   the dead radius by ~0.05 mm.
2. **Pitch before diameter.** Going from the physical 3 mm lattice to the drawn 4 mm halves the loss; a 2 mm
   pitch (V3B_final: 0.5 mm on 2 mm) at nominal margin gives 2.3 % at normal incidence — no better than today's
   0.8 mm on 3 mm — and 23 % if the margin is poor. Dense small pillars are the worst case when the gain or the
   drift field drops, because the dead ring does not shrink with the pillar.
3. **Diameter.** At nominal margin 0.5 mm on 3 mm is 1.0 % (normal incidence), 0.4 mm is 0.5 %; at good margin
   anything ≤ 0.5 mm is below 0.1 %. The 0.2-point ceiling of the old stage 36 is wrong; the gain from
   0.8 → 0.5 mm at the same pitch is 2-3 points at normal incidence, 1-1.5 on cosmics.
4. **Big pillars** cost ≈ 0.25 point with a ~1 mm dead ring; R_eff 4.2 mm, not 3.1.

## 6. What is still open

- **Physical check (bench, 30 min):** microscope on det2 (or any det1-det4 board / the leftover film): pillar
  pitch, diameter and the film offset. The data say 3.0 mm pitch, film shifted (−2, +1) mm in the pad frame.
  Every drawing and bulk sheet says Mask_M2_V2 at 4 mm; the printed film is not the Gerber in
  `Detector_Drawings/Version_Apr26`. Find which file was sent to the printer.
- **Why the margin moves the dead radius** (mechanism): the pillar edge is soft (charge loss, §2) and the
  r_eff grows when gain or collection drops. Candidates: field-line funnelling into the dielectric at low
  field ratio, positive charging of the pillar top. Garfield++ with the real geometry would settle it; the data
  only need the empirical table above.
- The all-angle no-hit at low drift field (32 % at < 200 V/cm far from any pillar) is a separate collection
  loss (attachment / slow drift) and is not in any handoff.
- det1 / det3 / det4: run `41_pillar_lattice.py` on their keys (needs `p2_study_cache` + `34`), then the loss
  budgets of those chambers get the same pillar row. det1_long5 already shows the lattice.
- CERN note §12 and the "pillars on hold" decision: correction block needed (not edited here).

## 7. Code (all on `main`)

| file | what |
|---|---|
| `41_pillar_lattice.py` | model-free lattice search, phase/origin, per-quadrant registration, shadow profiles, hard-disc fit, charge vs distance, per-big-pillar fits; writes `lattice_pad.npz` (board frame) |
| `42_pillar_shadow_vs_field.py` | the shadow on that lattice across a chamber's drift / mesh scans (r_eff and integrated shadow per point) |
| `36_pillar_extrapolation.py` | uses `41_pillar_lattice/lattice_pad.npz` when present instead of the Gerber small pillars |

Run order on a key: `p2_study_cache.py <key> --rebuild` → `35` → `34` → `41` (→ `36`); `42 <scan_key> <lattice_key>`
for a scan series. Outputs: `.../long_run_det2_430_730/41_pillar_lattice/`,
`.../p2_det2_long_run_430_repositioned_11-25-26/42_pillar_shadow_vs_field_det2_initial2/`.
Stages 38 and `nohit_checks/masktest|offscan*` still use the Gerber and should be read as "no shadow on the
4 mm lattice", which is now expected.

## 8. Clustering re-examined (same night)

Stage 34's matrix stands: tile adjacency and CC-20 clusters are identical (12.3 % multi-pad), time cuts do not
separate in-time second particles, seed thresholds only remove muons (50 ADC costs 0.8 pt, and those are the
soft-edge pillar tracks). The only thing below the hit finder's 5-sigma / 3-sample threshold is a matched-filter
second pass on the pad under the track (stage 39 tables, `pad_waveform_features.csv`):

| matched-filter cut | no-hit tracks with a pulse | noise control | net gain |
|---|---|---|---|
| 3 sigma | 14.0 % | 1.4 % | +0.40 pt |
| 4 sigma | 9.3 % | 0.6 % | +0.28 pt |
| 5 sigma | 4.2 % | 0.4 % | +0.12 pt |

The recoverable fraction is the same near and far from a pillar (13 / 14 / 15 %) and at all angles: it is the
Landau low tail, not the shadow. Pillar-lost tracks have 0.2 % of a normal pulse and cannot be recovered by any
clustering or threshold choice. Neighbour-seeded thresholds find 0.05 sub-threshold neighbour pads per hit track:
useless on 12 mm pads without charge spreading. Recommendation: keep the clustering; an optional track-guided
4-sigma second pass is a bench-only +0.3 pt, to be quoted separately if ever used. No code added.

## 9. How the mistake happened, and the rules that prevent it

1. Every pillar test (36, 38, `nohit_checks/masktest|offscan*`) folded the tracks on the Gerber lattice and
   scanned offsets within one cell. That excludes only that lattice and pitch; it was read as "no pillar
   effect". Rule: run the model-free periodicity search (`41_pillar_lattice.py`) BEFORE folding on any drawing.
2. The no-hit signatures (sub-mm scale, lateral-spread dependence, charge never arrives) all pointed at a dense
   periodic dead structure; the one contradicting test was trusted over three agreeing signatures. Rule: when
   one test contradicts the simplest mechanism, audit the test's assumptions first.
3. Stage 31's big-pillar edge of 2.4 mm against a 0.7 mm pointing was explained away as a "gain zone". It was
   the film offset. Rule: fit objects individually (centre + radius + edge) before fitting them jointly at
   drawn positions; an unphysical edge width is a registration error until proven otherwise.
4. Drawings, bulk sheets and slides all agreed on 4 mm; agreement among documents is not a measurement.
   Rule: a geometry used in a loss budget must be measured on the detector (data or microscope) and the
   measurement cited next to the number.
5. The run-to-run difference at identical voltages was filed under "position / mesh leak". Rule: a change in
   the vertical-track excess is first a shadow-depth (margin) question; check `42` before alignment.

These rules are also in `P2_basket_analysis/CLAUDE.md` so every session loads them.
