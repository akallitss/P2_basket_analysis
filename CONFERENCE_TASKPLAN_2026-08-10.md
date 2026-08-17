# Conference task plan — T–21 days

**Written 2026-08-10 (Monday). Conference 31 Aug (Mon), talk 4 Sep (Fri).**
Companion to `CONFERENCE_ROADMAP.md` (roadmap v2, 2026-07-12), which is now a month old and
whose WP structure needs re-weighting: **Act 3 happened.** This document is the execution plan
for the remaining 21 working days before the conference, 25 before the talk.

Sources for the status below: `sps_beam_analysis/PIPELINE_OVERVIEW_2026-07-26.md`,
`sps_beam_analysis/urw_reference/URW_TRACKING_HANDOFF_2026-07-25.md`,
`VMM_ONLINE_EFFICIENCY_PLAN.md` (2026-07-30), `sps_beam_analysis/logbook/`, git log to 2026-07-31.

---

## 1. What changed since the roadmap was written

The July 2026 H4 beam test (TB_July2026_H4, 23–31 July) ran and was analysed **during** the beam
time. The roadmap treated Act 3 as "depends on the beam-test data arriving in time"; it arrived,
and it is now the strongest material in the talk.

| Roadmap item | Status 2026-08-10 |
|---|---|
| WP6 3-plane analysis, inter-plane matching ("no code exists") | **Done and exceeded.** Stages 20–31 built; tag-and-probe (22) *plus* an absolute **uRWELL-referenced** efficiency (`urw_reference/`) that the roadmap never anticipated |
| 5b.1 2-of-3 tagging efficiency | Done — 92.5 / 96.3 / 92.5 % (IN/MID/OUT), stage 22 |
| Absolute efficiency (not in the plan at all) | **96–97 % ± 0.3 % syst** against uRWELL tracks, 1.15 M tracks/station |
| 5b.3 spatial residuals | Done — 3.39/3.47 mm rms = 12 mm/√12, i.e. binary pad response, no charge sharing |
| 5b.4 plane-to-plane time residuals | Done — σ = 15.5 (MID) / 18.4 (OUT) / 22.4 (IN) ns per station, walk-corrected, stable <0.5 ns over 4.5 h |
| 5a.2 bench-vs-beam ε(HV) | Beam half exists (mesh + drift curves, uRWELL-referenced and tag-probe). **Overlay not made** |
| 5a.1 DREAM vs VMM back-to-back | Data taken on the same detectors. VMM analysis has station scalars only — **not yet comparable** (see §3, T3) |
| 5b.2 3-plane event display | **Not made.** All ingredients exist |
| WP3 pion runs | **Still nothing.** `grep -ri pion` matches only the roadmap itself |
| WP2 efficiency vs rate (5/15 kHz) | Not run |
| WP1 SNR freeze | Unchanged since roadmap — polish only |
| Slides | **Do not exist.** No slide deck anywhere in the repo |

Two extra assets the roadmap did not know about:
- **Gas study → beam gas decision** (2026-08-01): Ar/CF₄/iC₄H₁₀ 88/10/2, 2.4× drift velocity,
  mesh 410→425 V for equal gain. Currently **uncommitted** (`gas_studies/decide.py` untracked,
  4 modified files).
- **Detector-lifetime autopsy** (`reports/det_lifetime_autopsy_2026-07`): det3 drift-foil HV
  decoupling, det1 mesh-contact arcing — a real methodology-value story for the QA slide.

**Consequence for the story:** the three-act structure survives, but the weight shifts. Act 3 is
no longer an outlook, it is the punchline — a P2-geometry telescope measured to 96–97 % absolute
efficiency with an independent reference, 16–22 ns timing, and a mapping error caught two ways.
Act 1 (Nov 2025 SNR/rate) should be compressed to make room.

---

## 2. Revised priority — what the talk needs

Ranked by *impact per remaining day*. P0 = the talk is weaker without it; P1 = strongly wanted;
P2 = cut candidates, explicitly listed in §5 so cutting is a decision, not an oversight.

**P0 — critical path**
- T1 Charging-up vs rate (blocks quoting 96–97 % honestly)
- T2 Bench ↔ beam ε(HV) overlay (validates the whole Act 2 methodology in one figure)
- T3 DREAM vs VMM like-for-like (the single justification of the dual-electronics strategy)
- T4 3-plane event display + residual/timing hero plots
- T7 Synthesis table (WP7)
- T8–T10 Slide build, dry run, freeze

**P1**
- T5 WP1 SNR freeze + presentation pass
- T6 WP2 efficiency vs rate, muons 5 vs 15 kHz
- T11 Timing narrative consolidation (bench 46 ns → beam 16–22 ns → Garfield 3–7 ns → gas)
- T12 QA/methodology slide from the existing reports

**P2 — cut candidates**
- T13 WP3 pion runs (timeboxed scoping only — see §5)
- T14 WP4 geometry comparison
- T15 Housekeeping (commit gas_studies, VMM `run_config_beam.py`)

---

## 3. Task detail

### T1 — Efficiency vs time within a sub-run: charging-up or rate? *(P0, ~1 day)*
`highstat_eff_1` loses 1.3 points monotonically over 2.5 h at constant 4.6 kHz on MID/OUT and
recovers after a 2.5 h pause — but sub-run 05 changes *both* the pause and the rate (3.2 vs
4.6 kHz), so the two explanations are not separated (handoff §13.4). Every track already carries
`t_ns`, so binning efficiency **within** a sub-run is nearly free and breaks the degeneracy.
- Add a `--time-bins` mode to `urw_reference/urw_p2_efficiency.py`; run on sub-runs 00–05.
- Deliverable: ε(t) across the full 4.5 h, one curve per station.
- **Why it is P0:** the headline number is 96–97 %; if it is drifting, the talk must say which
  number it is and why. Do this before any slide quotes an efficiency.

### T2 — Bench ↔ beam ε(HV) overlay *(P0, ~1 day)*
Wish-list item 5a.2, now doable entirely offline — both halves exist.
- Beam: `analysis/urw_referenced_efficiency/drift_mesh_scan_1/` mesh curves (MID 0.806→0.961,
  OUT 0.857→0.950 over 390→450 V; IN 0.314→0.917 over 370→430 V) and the drift curve
  (0.157 at zero field → 0.974 above 850 V).
- Bench: cosmic-bench stage 11 (`11_hv_scan_efficiency`) and stage 16 (`16_drift_scan_efficiency`)
  outputs for det1/det2/det4.
- One figure, ε vs mesh V, bench and beam on the same axes, with the caveats stated on the slide:
  different gas, different detectors of the same generation, cosmics vs MIPs, different reference.
- **Watch out:** the two are *not* the same detectors. The honest claim is "the bench reproduces
  the shape and the working point", not "the same number". Decide the wording early, it drives
  the figure.

### T3 — DREAM vs VMM on the same detectors *(P0, ~3–4 days, highest risk)*
The most-requested plot when two readouts coexist, and the direct justification of the
dual-electronics strategy. **Currently blocked by a denominator mismatch, not by physics.**
VMM efficiency (IN 0.11 / MID 0.53–0.58 / OUT 0.56–0.64) is a whole-station scalar over *every*
trigger, while only connectors c4–c6 of ~10 are instrumented; DREAM's 96 % is per-pad over the
illuminated overlap. They must never appear on the same axis as they stand.
- Path: the pad map is done (`vmm_stations.build_pad_table()`), and
  `vmm_efficiency.coincident_mask()` gives the in-time hit selection. What is missing is the
  **per-pad efficiency map** — restrict the denominator to triggers whose track crosses
  instrumented pads (VMM plan, Stage 3/4).
- Fallback if the acceptance restriction does not converge in ~2 days: drop the efficiency
  comparison and show what *is* already comparable on the same detector at the same HV —
  cluster size (VMM 1.11–1.18 pads vs DREAM 1.05–1.08), ADC/MPV spectra, timing σ (VMM pinned at
  the 22.5 ns BCID bin vs DREAM 15.5–22.4 ns), and the two hit maps side by side. That is still a
  real slide, and it makes the honest point that VMM timing is calibration-limited.
- **Known blocker to state on the slide, not to fix:** VMM σ is pinned at the BCID quantisation
  because there are no time-calibration runs (`qa_config.py: 'calibration': None`). This is an
  outlook item, not a defect.
- P2_IN's low VMM efficiency is **gain, not mapping** — it sat at mesh 440 V while MID/OUT were at
  450 V, still on the steep part of its turn-on. Say so explicitly or a reviewer will read it as a
  broken station.

### T4 — Hero plots from the beam test *(P0, ~2 days)*
1. **3-plane track event display** on the true pad geometry (wish-list 5b.2 — the picture that
   makes "this is the P2 tracker" credible). Ingredients: uRWELL track + three P2 leading clusters
   matched on `eventId`, drawn on the fan-pad tiles. Nothing new is needed physics-wise.
2. **Residual plot** — the 2D residual showing the pad footprint as a tilted filled square, with
   3.46 mm = 12 mm/√12 annotated. This is the spatial-resolution number for any Q² statement.
3. **Per-pad efficiency map** in the P2 frame at the working point, with the uRWELL-frame map
   beside it (the 5×5 mm reference hole appearing in all three stations is a nice methodology aside
   if there is room — otherwise backup).
4. **ε vs HV curves**, mesh and drift, uRWELL-referenced.

### T5 — WP1 SNR freeze *(P1, ~1 day)*
Unchanged from the roadmap: SNR vs (sg, snt) matrix per detector at 5 and 15 kHz; one channel-level
SNR map at best config; the headline "best config = (sg, snt), SNR = X large / Y small". Pipeline
is done — this is a presentation pass, not analysis.

### T6 — WP2 efficiency vs rate, muons *(P1, ~2 days)*
`vmm_detector_efficiency.py` over the matched config points of the `5kHz` and `15kHz` DATASETS
registries → ε per detector per rate; on-spill occupancy map at 15 kHz; ADC MPV vs rate.
**Memory note:** the registries already carry deliberate `n_files` settings (10 at 15 kHz,
3 at 5 kHz — "~105 MB/file, keep memory low"). Do not raise them without checking peak RSS.

### T7 — Synthesis table *(P0, ~1 day, do after T1–T6)*
One table: detector generation × electronics × environment → efficiency, SNR, timing, spatial
resolution, rate reach. Rows now genuinely available: Nov-2025 large/small + VMM + beam;
cosmic bench det1/det2/det4 + DREAM + cosmics; July-2026 IN/MID/OUT + DREAM + beam; same three
+ VMM + beam. This is the slide people photograph.

### T8 — Slides v1 *(P0, 17–21 Aug)*
Follow the roadmap §7 skeleton, re-weighted: Act 1 down to 2–3 slides, Act 3 up to 4–5.
Suggested 15-min structure: motivation (2) → Act 1 SNR + rate (3) → geometry verdict (1) →
Act 2 bench + QA methodology (3) → Act 3 telescope, event display, efficiency, timing,
DREAM/VMM (5) → synthesis + outlook (1).

### T9 — Internal dry run *(P0, by 26 Aug)*
Roadmap rule stands: only cosmetic plot changes after the dry run.

### T10 — Freeze + backup slides *(P0, 27–31 Aug)*
Backup set: spill mask / Δt fit / MAD noise methodology, mapping validation vs M3 **and** the
independent uRWELL confirmation of the flipped `c_5_top` ribbon, det1/det3 lifetime autopsy,
Garfield timing simulation, gas study.

### T11 — Timing narrative *(P1, ~0.5 day)*
Four numbers already exist and line up into one slide: bench 46 ns (drift-geometry limited) →
beam 15.5–22.4 ns per station at the working point → Garfield physics floor 3–7 ns → the 20 ns
P2 goal, now **met** at MID/OUT. Add the drift plateau (750–800 V, working point within 0.5 ns of
optimum) and the gas decision (Ar/CF₄/iC₄H₁₀ 88/10/2, 2.4× v_d) as the forward path. Cheap, and
it closes the loop opened by the July timing study.

### T12 — QA/methodology slide *(P1, ~0.5 day)*
From `reports/det_lifetime_autopsy_2026-07` and `reports/p2_efficiency_intermittency_2026-07`:
the bench *found* drift-foil HV decoupling and mesh-contact arcing before those detectors reached
a beam. Frame as methodology, never as failure.

---

## 4. Schedule

| Window | Work | Milestone |
|---|---|---|
| **Mon 10 – Sun 16 Aug** (T–21 → T–15) | T1, T2, T3 start, T5, T15 | Efficiency number settled; bench↔beam overlay exists; SNR frozen |
| **Mon 17 – Sun 23 Aug** (T–14 → T–8) | T3 finish or fall back, T4, T6, T11, T12, **T8 slides v1** | Every hero plot exists; deck assembled end to end |
| **Mon 24 – Sun 30 Aug** (T–7 → T–1) | T7 synthesis, **T9 dry run by Tue 26**, T10 backups | Content freeze after the dry run |
| **Mon 31 Aug – Fri 4 Sep** | Conference; talk Friday | Rehearse only |

Hard gate: **T3's go/no-go is Fri 21 Aug.** If the per-pad VMM efficiency is not converging by
then, take the fallback comparison and stop — do not spend the last week on it.

---

## 5. Explicit cut list

These are dropped unless the P0 work finishes early. Recording them so the cut is deliberate.

- **T13 WP3 pion runs.** The roadmap's highest-uncertainty item, still at zero code, on Nov-2025
  data. Its purpose was rate reach via local rate density — and the July telescope now delivers a
  rate story with real tracks. **Recommendation: timebox to half a day of scoping** (find the run
  numbers in the elog, confirm the data is where the registries expect it). If it is not trivially
  ready, it becomes one outlook line, not an analysis. Revisit for the paper, not the talk.
- **T14 WP4 geometry comparison.** Cheap *if* T5 and T6 land; otherwise a sentence on the Act-1
  verdict slide rather than a figure.
- **The `drift_mesh_2d_1` resume** and the other "nice-to-have additional runs" in the pipeline
  overview §6 — no beam before the conference, so these are outlook regardless.
- **VMM time calibration (Stage 5).** Blocked on calibration runs that do not exist. Outlook line.

## 6. Housekeeping (T15, ~1 h, do this week)

- `gas_studies/` has 4 modified files and untracked `decide.py` carrying the beam-gas decision —
  commit per-file, no Co-Authored-By trailer.
- `run_config_beam.py` on the DAQ side has the corrected VMM trigger channel (VMM 0 ch 44) applied
  in place but **not committed**, with a backup only in `/tmp` — that will not survive a reboot.

## 7. Open questions for you

1. **Talk length and title** — the roadmap assumed ≈15 min. Confirmed? The re-weighting in T8
   depends on it.
2. **Is the July beam-test material clearable for an external talk**, given there is no refereed
   BASKET paper and the specs come from collaboration internal notes? Worth asking the
   collaboration now rather than in the last week.
3. **Is any of the Nov-2025 pion data actually accessible** from here, or is it only on the
   `/drf/projets/clas12/` cluster path in the registries? That answer decides T13 immediately.
