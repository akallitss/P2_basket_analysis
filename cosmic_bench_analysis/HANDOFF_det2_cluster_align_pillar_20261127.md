# HANDOFF — det2 clustering / 3D alignment / pillar extrapolation, 27 Nov 2026 (afternoon)

Branch `det2-cluster-align-pillar` (worktree `.claude/worktrees/det2-cluster-align-pillar`).
Run: `det2_long_repos` only. Outputs under `.../Analysis/det2/p2_det2_long_run_430_repositioned_11-25-26/long_run_det2_430_730/{34_study_cache,34_clustering_matching,35_alignment_3d,36_pillar_extrapolation}`.

Run order (each < 1 min; project venv `P2_basket_analysis/.venv/bin/python`):
`p2_study_cache.py <key> --rebuild` → `35_alignment_3d.py` (uses `34_study_cache/clusters_cc20.pkl`, made by the cache helper run in this session; regenerate with `p2_study_common.cc_clusters`) → `34_clustering_matching.py` → `36_pillar_extrapolation.py`.
Baseline check: the cache reproduces stage 06 exactly (36,382 active rays, 95.9 % / 96.8 %).

## 1. Alignment (question: 3D rotations / are border misses misalignment?)

- **No tilt.** The local crossing-height map in 40 mm cells is flat (χ² 96 → 94 with a plane, 76 ndf). The global linear model gives a tilt of 0.06° ± 0.05°, a bow of +0.10 ± 0.06 mm per (100 mm)². In-plane effect < 0.1 mm.
- **But the stage-06 2D transform is biased.** `fit_transform` is plain least squares on (track, max-cluster) pairs, and 4.7 % of those pairs are a second particle or noise ~100 mm away. The fit therefore shrinks the scale to **0.978** (true 1.000: 4.4 mm at ±200 mm), and the 03 z-scan plane is 240 mm where the slope fit gives **246.0 mm**. The residual map showed a radial inward pattern of up to 2 mm before the fix and < 0.5 mm after.
- Fix: `p2_align.fit_transform(..., trim_mm=15)` (opt-in, default unchanged): s = 1.0002, RMSE 39.7 → 4.7 mm.
- With the corrected alignment the live-area boundary (outer edge and dead tiles alike) is a sharp step: miss rate 48 % at −0.5 mm, 8 % at +0.5 mm, flat 3.9–4.4 % from +2 mm inside. With the old alignment the excess extended to about +5 mm, which is why `DEAD_GUARD_MM = 6` was needed. **The border losses were misalignment, not a tilt.**
- New denominator proposal: signed distance to the union of live tiles ≥ 2 mm. That gives 40,625 rays (+12 %) and **eff 95.9 ± 0.1 %**. Old: 36,382 rays, 95.9 %.

## 2. Pointing resolution and clustering (question 1)

- **Pointing σ = 0.70 / 0.75 mm** per pad axis (single-pad box ⊗ Gauss fit) and **0.63 mm** from the boundary step. The 3.0 mm quoted in the note came from the scale bias. The fitted box is ~0.9 mm narrower than the tile because single-pad events are selected away from the pad edges.
- Clustering: tile adjacency (pads whose tiles touch) and CC 20 mm give identical results (12.3 % multi-pad), so the clustering rule doesn't matter for P2.
- Timing doesn't separate second particles: all real pulses are within about ±60 ns, and distant pairs are as in-time as neighbours.
- A seed threshold only loses muons: 95.1 % at 50 ADC, 93.3 % at 100, 85.5 % at 300. Accidentals measured by event mixing (shift of 1000 events) are 0.94 % at centroid < 20 mm and 0.20 % at tile < 3σ. **Accidental-corrected efficiency 95.91 % (CC20 / centroid 20 mm, fiducial ≥ 2 mm)**.
- Tile < 3σ matching gives 95.56 % (corrected). The 0.35–0.5 % difference is tracks whose nearest fired pad is 2.2–20 mm away (median 5.6 mm). They are uniform over pads (no mapping swap), grow with angle (0.38 % at 0–10°, 1.38 % at 20–30°) and have muon-like amplitudes, so they are most likely reference-track tails and real hits. Recommended matching: tile distance ≤ 5σ (3.6 mm) or centroid ≤ 12 mm with accidental subtraction, or keep 20 mm and quote the corrected value.

## 3. Pillars (question: extrapolate the shadow to smaller pillars)

- The Gerber small pillars (r 0.40 mm, 4.0 mm pitch) are not one global lattice: rows are offset, about 10 per pad but not at the same place in each pad. Use the actual nearest pillar (KD-tree), never fold modulo the pitch.
- **Direct measurement with 0.7 mm pointing**: vertical tracks (< 6°) show 7.6 % no-hit within 0.4 mm of a pillar centre (n = 223) against 4.6 % ≥ 1.2 mm away. The local fit of the electron-level model (Heed clusters, Magboltz D_T 0.048 √cm so σ_T 0.31 mm over the gap, Polya 0.5 gain, a threshold, pointing smear) with a free baseline per angle class gives a dead ring of −0.08 mm (95 % CL −0.30 … +0.06), i.e. r_eff 0.32 mm. χ² is 42.4/34 against 43.5 with no pillar term at all.
- **Small-pillar loss in this run: 0.02 %, < 0.23 % at 95 % CL.** The 2.0 % "small-pillar shadow" in the note's loss budget is excluded: a pillar loss of 1 % gives Δχ² ≈ +140.
- Extrapolation (`loss_vs_diameter.csv`): at 0.5 mm diameter the loss is ≤ 0.02 % (cosmic) and ≤ 0.04 % (normal incidence) even at the 95 % CL upper edge. At the current 0.8 mm it is ≤ 0.21 % and ≤ 0.44 %. **Shrinking the pillars gains at most ~0.2 points.** Diffusion (0.3 mm) already makes 0.8 mm pillars mostly transparent.
- Big pillars (σ fixed to 0.72 mm): R_eff 3.40 ± 0.41 mm (Gerber 3.08), loss 0.14 %.

## 4. New open question: angle-dependent no-hit, uniform in position

No hit at all in the event: 4.9 % at 0–4°, 4.3 % at 4–8°, 2.9 % at 8–12°, 2.0 % at 12–16° and 0.95 % above 16°. It is the same ≥ 1.4 mm from any pillar, flat across the pad, round in (θx, θy) (not an M3 one-projection ridge), and independent of M3 χ² and of time. The low-amplitude tail grows the same way (< 100 ADC: 3.3 % at 0–4° against 1.6 % above 16°), while the median amplitude only moves 1150 → 1280. This is the largest remaining loss (~1.5–2 points on this run's angle mix).

Next: dump `decoded_root` waveforms of the pads under near-vertical no-hit tracks (with eventId and θ) to separate "sub-threshold pulse" from "nothing at all"; compare with the 11-24 alignment run and det5.

## 5. Not done / decisions for Alexandra

- The pipeline still uses the biased fit. Switching stages 06/11/12/16/17 to `trim_mm=15` and setting det2 `plane_z ≈ 246` changes every detector's maps. Re-run them all before quoting numbers side by side.
- Not merged to main. The CERN note still says pointing 3 mm and pillar shadow 2.0 %: both need correcting.

## 6. Implemented and validated on det2_long_repos only (27 Nov, later)

Alexandra's decisions: implement the corrected alignment and the new guard on ONE detector first, validate it,
hold the pillar theory, and put the findings in the CERN note.

- Config (`p2_qa_config.py`, default off for every other key): `fit_trim_mm` (FIT_TRIM_MM) and `edge_guard_mm`
  (EDGE_GUARD_MM). `det2_long_repos`: plane_z 240 → 246, fit_trim_mm 15, dead_guard 6 → 2, edge_guard 2.
- Code: `p2_align.outer_edge_distance` + edge guard in `mask_dead_cells`; `trim_mm=cfg.FIT_TRIM_MM` in 06/11/12/16;
  `ct.attrs['edge_guard']` set in 06/11/12/16/28. Stage 17 is det3-only (no cfg), so it was left alone.
- Guard choice: dead-border excess with the corrected alignment is +14.6/+5.2/+2.9/+1.7/+0.7 % in 0.5 mm bins and 0
  from 2.5 mm; the efficiency plateau starts at 1.5–2 mm → 2 mm (≈ 3σ of the 0.7 mm pointing), same at the outer edge.
- 37 step table (same data): A old 36,382 / 95.89 → B trimmed fit 37,892 / 95.38 → C z 246 37,827 / 95.40 →
  D guard 2 43,202 / 95.48 → E + edge guard 40,530 / 95.92; control (guards 2/2, old alignment) 39,103 / 95.09.
- Full pipeline re-run (`run_p2_pipeline.sh det2_long_repos`): 06 = 40,530 rays, 95.9 % / 96.8 %, s 1.000, RMSE 4.7 mm;
  12 = 0/7 flags, core σ 3.6/3.4 mm, z plateau peaks at 246; 03's own z-scan now finds z* = 246.0;
  27: only-distant 0.65 → 0.60 %, pad edges 0.14 → 0, big pillars 0.26 → 0.24 %; 33: det2 share 0.09–0.18 %.
  Old products are kept in `long_run_det2_430_730__v1_stage06_align_20261127/`.
- Still old/held: 27 books 2.18 % to "small pillars" from 31's angle model (χ²/ndf 42.6/10). 38 shows the angle
  trend is identical ≥ 1.4 mm from any pillar, so this attribution is wrong but left untouched (pillar theory on hold).
  Stage 31's big-pillar fit still floats σ = 2.4 mm (with σ fixed at 0.72: R_eff 3.40 mm): worth a look later.
- CERN note: callout at the top, correction markers in §7/§8/§11 and glossary, new §12 (alignment, pointing, step
  table, re-run table, new map, model-free no-hit vs angle), todo renumbered to §13. Edited in the cern-site working
  tree, NOT committed or deployed yet (waiting for her go-ahead).
