# HANDOFF — det2 (P2_2) cosmic bench, end of 27 Nov 2026: alignment fix, pointing, no-hit investigation

**Start here.** This supersedes the morning handoff (`HANDOFF_det2_efficiency_20261127.md` on `main`) for the
analysis. It also absorbs the running notes in `HANDOFF_det2_cluster_align_pillar_20261127.md` (same branch, more
detail per study). Bench/hardware state (§7) is unchanged since the morning.

## 0. Where things are

| What | Where | State |
|---|---|---|
| Analysis code | `P2_basket_analysis` branch **`det2-cluster-align-pillar`**, worktree `P2_basket_analysis/.claude/worktrees/det2-cluster-align-pillar` | commits on top of `main` (631c6b5), **NOT merged**. Pushed to GitHub `origin` (akallitss/P2_basket_analysis, branch `det2-cluster-align-pillar`). Per the new rule (work on `main`, no worktrees unless agreed), the next step is for Alexandra to merge it into `main` and continue on `main` |
| Python | `P2_basket_analysis/.venv/bin/python` (the worktree has no .venv; put that `bin` on PATH or call it explicitly) | |
| det2 outputs | `…/data/Cosmic_Bench/Analysis/det2/p2_det2_long_run_430_repositioned_11-25-26/long_run_det2_430_730/` (new) | re-run with the corrected config |
| "Before" outputs | same folder `…/long_run_det2_430_730__v1_stage06_align_20261127/` | backup of the old-alignment products |
| 11-24 scan study | `…/p2_det2_long_run_430_repositioned_11-25-26/40_nohit_vs_field_11-24/` | |
| CERN note | `cern-site` commit `702394f` (section 12 + correction markers) | **deployed and live** (verified). **Not pushed to GitHub** (`main` ahead 1): run `git -C ~/Documents/PostDocSaclay/cern-site push` |
| Raw data | LaCie: `…/data/Cosmic_Bench/<run>/` (`/mnt/cosmic_data` is not mounted on this PC) | decoded_root for the repositioned run is present |

Gotchas:
- Call `run_p2_pipeline.sh` by **absolute path** with the venv on PATH. It re-execs itself under systemd-run with
  an 8 GB cap.
- `garfield_inputs/*.csv*` are git-ignored. In a worktree, stage 36 reads the main checkout's copy.
- Stage 17 is det3-only (no `cfg`); do not add `cfg.*` arguments to it.

## 1. What changed in the pipeline (det2_long_repos ONLY)

Alexandra's decision: validate on one detector before generalising. Every other run key is byte-for-byte unchanged.

- `p2_qa_config.py`: new per-run fields `fit_trim_mm` (FIT_TRIM_MM) and `edge_guard_mm` (EDGE_GUARD_MM), both
  default `None` = historical behaviour. `det2_long_repos`: `plane_z` 240 → **246**, `fit_trim_mm` **15**,
  `dead_guard_mm` 6 → **2**, `edge_guard_mm` **2**.
- `p2_align.fit_transform(..., trim_mm=None)`: iterative trimmed Procrustes.
  `p2_align.outer_edge_distance`: rasterised signed distance to the outer boundary of the pad array.
  `mask_dead_cells` applies the edge guard if `ct.attrs['edge_guard']` is set.
- Stages 06/11/12/16 pass `trim_mm=cfg.FIT_TRIM_MM`. Stages 06/11/12/16/28 set `ct.attrs['edge_guard']`.

**Why.** The plain least-squares fit is biased: ~4.7 % of (track, max-cluster) pairs are a second particle or
noise ~10 cm away, which shrinks the scale to 0.978 (4.4 mm at ±200 mm). The trimmed fit gives 1.000, fit RMS
39.9 → 4.7 mm, and 03's z-scan now also lands on 246.0. **No tilt:** 0.06° ± 0.05° (stage 35).

**Step by step, same data (stage 37):**

| step | tracks | eff |
|---|---|---|
| A · old | 36,382 | 95.89 % |
| B · + trimmed fit | 37,892 | 95.38 % |
| C · + z 246 | 37,827 | 95.40 % |
| D · + dead guard 2 mm | 43,202 | 95.48 % |
| E · + edge guard 2 mm | 40,530 | 95.92 % |
| control · guards 2/2 with OLD alignment | 39,103 | 95.09 % |

The old number was right because two errors cancelled: the shrunk map hid out-of-edge tracks, and the 6 mm guard
hid the misalignment at dead pads.

**Full re-run (official):**
- 06: 40,530 tracks, **95.9 % / 96.8 %** (any hit).
- 12: 0/7 bias flags, core σ 3.6/3.4 mm (was 4.5).
- 27: only-distant 0.65 → 0.60 %, pad edges 0.14 → 0, big pillars 0.26 → 0.24 %.
- 33: det2's share 0.09–0.18 %.
- Accidentals (event mixing): 0.93 % → corrected eff 95.88 %.

## 2. Measured quantities to reuse

- **Pointing resolution at P2: σ = 0.70 / 0.75 mm** (single-pad box ⊗ Gauss), 0.63 mm (edge step). It grows from
  0.65 mm (vertical) to 0.9 mm (25°). The "3 mm" in the original note was the scale bias.
- Dead-border excess with the corrected alignment: +14.6 / +5.2 / +2.9 / +1.7 / +0.7 % in 0.5 mm steps, 0 from
  2.5 mm. The plateau starts at 1.5–2 mm → **2 mm guard (≈ 3σ)**.
- Clustering: tile adjacency ≡ CC 20 mm (12.3 % multi-pad); the rule doesn't matter for P2. Timing doesn't
  separate second particles. Seed thresholds only lose muons.
- Tile-within-3σ matching loses 0.35–0.5 % real hits (reference tails, grow with angle, uniform over pads, no
  mapping swap). Keep centroid < 20 mm and quote the accidental-corrected value.

## 3. Pillars: ON HOLD (Alexandra's decision); do not ship the old attribution

- Stage 31 / 27 still book **2.18 %** to "small pillars" from the old angle model (χ²/ndf 42.6/10).
  **This is wrong.** The no-hit vs angle trend (≈5.4 % vertical → 1.3 % above 16°) is identical ≥ 1.4 mm from any
  pillar (stage 38).
- Direct pillar tests with 0.7 mm pointing (Saclay V1/V2 and CERN masks; global and per-quadrant offset scans of
  the pattern over a full cell): no shadow. The electron-level model (stage 36) gives a small-pillar loss of
  0.02 % (< 0.23 % at 95 % CL), so going to 0.5 mm pillars gains ≤ 0.2 pt. **Not presented as a conclusion**:
  pillars are on hold; revisit only when asked.
- Saclay mask = staggered lattice (4 mm along rows, rows ~3.46 mm apart, offset 2 mm). Use a KD-tree on the real
  Gerber positions, never fold modulo the pitch.
- Stage 31's big-pillar fit still floats σ ≈ 2.4 mm (with σ fixed at 0.72: R_eff 3.40 mm). Possibly a gradual
  low-gain zone around big pillars; not investigated.

## 4. No-hit investigation (largest remaining loss): status

- **The charge never arrives** (stage 39, raw waveforms, all decoded files):
  - All 1,298 no-hit events are present in both P2 FEUs.
  - On the best pad within 3 mm: ~87 % pure noise; ~13 % a 3–6 σ sub-threshold pulse, angle-independent.
  - Stacked mean ≈ 3 ADC = 0.2 % of a hit.
  - The hit finder needs > 5 σ for ≥ 3 samples (`analyze_waveforms`, `thresholdSigma = 5`).
- **Real muons** (`nohit_checks/m3tests.py`): same top/bottom-pair kink, χ², NClus and trigger spacing as hit
  tracks. M3 planes are at z 24/144 and 1185/1302; P2 at 246.
- **Spatially random**, except 4 of the 5 big pillars (4 % of no-hits) (`spots*.py`).
- **Drift-field dependent** (stage 40, 11-24 scans). Vertical excess = no-hit 0–6° minus ≥ 12°:

  | drift field [V/cm] | < 200 | 200–400 | 400–600 | 600–800 | 800–1000 | > 1000 |
  |---|---|---|---|---|---|---|
  | excess [points] | ~11 | 8.8 | 4.9 | 2.1 | **0.5 ± 0.4** | 1.7–2.6 |

  On the mesh scan the excess is roughly flat (2–5 pts) while the overall level follows the gain.
- **Excluded:** sub-threshold, M3 fakes/kinks/pile-up, DAQ loss, path-length (1/cos θ) mechanisms (too weak),
  small pillars (any registration), mask lines, board vias (1 per pad, 0.1 % of area), fixed defects.
- **Open:** a mechanism depending on the lateral spread 4 mm·tan θ, below 0.7 mm scale, and on the drift field.
- Run-to-run at 430/730: 11-24 alignment run excess ≈ 1.6 pt (2.98 vs 1.36 %), 11-25 long run ≈ 4 pt
  (5.4 vs 1.3 %). Position changed; the mesh leak may have been developing.

## 5. Open items, in suggested order

1. **Validate the new config on a second det2 dataset**: `det2_initial2` (11-24 alignment run) with
   `plane_z ≈ 246.3`, `fit_trim_mm 15`, guards 2/2, and stage 37 on it. Do this before generalising to the
   other keys or detectors (Alexandra's explicit rule).
2. No-hit mechanism. Options offered, none chosen yet:
   - (a) stage 40 on det1/det5 drift scans (diagnostic only, no pipeline change);
   - (b) Garfield++ microscopic mesh-region simulation vs drift field;
   - (c) bench: efficiency at drift 780–830 V (mesh 430) on the next run, since the vertical excess is smallest at
     800–1000 V/cm. **Operating-point suggestion worth ~1 pt.**
3. Correct stage 27's loss table once the no-hit cause is understood (the "small pillars 2.18 %" row).
4. Only after the above: switch other run keys (det2 scans, det1/det3/det4/det5) to `fit_trim_mm`/guards, re-run
   everything, and compare. All of their maps carry the same ~2 % scale bias today.
5. Merge the branch into `main` (Alexandra), then work on `main` from now on. Push to GitHub origin via gh HTTPS, never to the CEA GitLab.
6. CERN note: section 12 describes the alignment and the pillar re-examination. The no-hit/drift-field findings
   (§4) are **not** in the note yet; add them when the cause is clearer or if asked.

## 6. New stages (all on the branch)

| stage | purpose |
|---|---|
| `p2_study_cache.py` / `p2_study_common.py` | per-run cache (hits, rays, geometry) + helpers; reproduces 06 exactly |
| `34_clustering_matching.py` | robust alignment, pointing resolution, clustering/matching/seed variants, accidentals |
| `35_alignment_3d.py` | tilt/bow fit, local z map, residual maps |
| `36_pillar_extrapolation.py` | electron-level pillar model and loss vs diameter (pillars on hold) |
| `37_alignment_validation.py` | step-by-step table A–E + control, profiles, residual maps |
| `38_nohit_vs_angle_pillar_distance.py` | model-free no-hit vs angle near/far from pillars |
| `39_nohit_waveforms.py` | raw waveforms under no-hit tracks (~8 min, reads all decoded files) |
| `40_nohit_angle_vs_field.py` | no-hit vs angle across the 11-24 drift/mesh scans |
| `nohit_checks/` | M3 quality, spatial clustering, mask and offset scans (see its README) |

Run order on a key: `p2_study_cache.py <key> --rebuild` → 35 → 34 → 37/38 → 39. Stage 40 is standalone (11-24 run).
Stages 35/37 need `34_study_cache/clusters_cc20.pkl`; 37 regenerates it if missing.

## 7. Bench / DAQ (unchanged since the morning; see `HANDOFF_det2_efficiency_20261127.md` §1–2)

- det2 HV all off; the mesh has a ≈ 10 MΩ resistive short. Inspect before any HV.
- Re-solder connectors 2/3/8; check the connector-5 pin row.
- `run_config.py` still points at the never-taken 3 h batch-2 run; `processor_config.py` `include_runs` points at
  it too. Repoint before the next run.
- No HV-trip / DAQ-stall watchdog yet (the 01:07 trip went unnoticed for 5 h).
