# Conference Roadmap — P2 Basket Pad-Micromegas: from SPS Beam Test to Cosmic-Bench Characterization

*Status: **roadmap v3, 2026-08-10**. **Conference: 31 August 2026 — talk: 4 September 2026**
(T–21 days to the conference, T–25 to the talk).*

> **What changed from v2 (2026-07-12): Act 3 happened.** The July 2026 H4 beam test
> (TB_July2026_H4, 23–31 July) ran and was analysed during beam time. v2 treated the three-plane
> measurement as "depends on the beam-test data arriving in time" and as the talk's outlook; it is
> now the talk's strongest result, and it delivered more than was asked for — an **absolute**,
> uRWELL-referenced efficiency that v2 never anticipated. The work-package weighting, the
> timeline and the slide skeleton are all re-cut accordingly. §2 (literature anchors) is unchanged
> and still verified. The dated execution plan lives in `CONFERENCE_TASKPLAN_2026-08-10.md`.

---

## 1. The story to tell

One sentence: **"We are developing pad-Micromegas trackers for the P2 experiment; we optimized the
VMM3a front-end and validated two geometries at the SPS, characterized the new detector generation
on a cosmic bench with full external tracking, and then returned to the SPS with three detectors
mounted exactly as they will be installed in P2 — same plane distances, same electronics — and
measured three-plane track coincidences at 96–97 % efficiency."**

The three-act structure survives, with the weight moved to Act 3.

### Act 1 — Beam test at SPS (Nov 2025): geometry + electronics optimization
- Two detector geometries under test simultaneously with VMM3a readout:
  - **P2 Large detector** — fan-pad geometry, 256 pads (VMMs 12–15)
  - **P2 Small detectors 1 & 3** — compact pad geometry (VMMs 8–11)
- **SNR optimization**: full VMM3a shaping scan (gain `sg` = 3.0/4.5/6.0 mV/fC × peaking time
  `snt` = 50/100/200 ns), signal/noise run pairs via neighbour-triggering, MAD noise baseline + MPV
  signal → SNR per VMM and per channel. → *"We found the optimal operating point of the electronics."*
- **Rate performance**: identical config points taken at ~5 kHz and ~15 kHz muon rates, plus pion
  runs with localized, much higher rate density. → *"We probed the detector where it will actually
  have to live."*
- Trigger-based efficiency with accidental subtraction (Δt peak fit + sideband) and per-pad
  efficiency maps on the real Gerber geometry.

*Re-weighted for v3: compress to 2–3 slides. This act is now the setup for the punchline, not the
punchline. The pion sub-story is a cut candidate (§5).*

### Act 2 — Cosmic bench (2026): new detector generation, precision characterization
- New P2 detectors (det1–det4), 256-pad fan geometry, read out with DREAM FEUs behind the
  **M3 beam telescope** → track-projected, *geometry-resolved* measurements a beam trigger alone
  cannot give. (**DREAM was the only readout option available on the cosmic bench** — this is why
  the bench characterization is DREAM-based and why the beam test had to carry both electronics.)
  The bench delivers:
  - pad-level and sliding-window **efficiency maps** on the true pad footprint
  - **HV turn-on curves** (mesh and drift), per detector
  - **timing** from waveform TOA (~46 ns vs trigger, drift-geometry limited)
  - **stability**: HV spark monitoring/veto, per-pad spark flagging, pedestal QA, long-run
    duty-cycle tracking
- The bench methodology **found real detector pathologies before those detectors ever saw a beam**:
  det3's drift-foil HV decoupling and det1's mesh-contact arcing, both diagnosed to a cause and
  written up (`reports/det_lifetime_autopsy_2026-07`,
  `reports/p2_efficiency_intermittency_2026-07`). Frame as methodology value, never as failure.
- → *"Before the next beam test, every detector is fully mapped, efficiency-calibrated and
  stability-vetted on the bench."*

### Act 3 — Beam test at SPS (July 2026): the P2 telescope in its final configuration — **done**
No longer an outlook. This was a **dress rehearsal of the P2 tracker itself**, and it worked:

- **Three detectors in a row on the H4 beam, mounted as they will be installed in P2** — P2_IN,
  P2_MID, P2_OUT at z = 320, 630, 940 mm — same plane-to-plane distances, same mechanics, same
  electronics chain as the experiment.
- **Two external uRWELL reference planes** (EIC uRWELL front/back at z = 0 and 1370 mm, 1370 mm
  lever arm) bracket the telescope. This is the ingredient v2 did not plan for, and it changes the
  measurement from relative to **absolute**: real extrapolated tracks, ≲0.45 mm per reference
  plane, against which each P2 station is a device under test.
- **Two front-end electronics on the same detectors: DREAM and VMM3a**, closing the loop between
  the DREAM-based bench campaign and the Nov 2025 VMM config-scan optimum.
- The complementarity remains the selling point: **DREAM (full waveforms, external trigger) for
  understanding; VMM3a (self-triggered, high-rate) for operating.**

---

## 2. What P2 needs from these detectors (literature anchors for the motivation slide)

*Unchanged from v2. Verified against the published literature 2026-07-12.* Use these numbers on the
motivation slide; items flagged ⚠ have no public citation yet and must come from collaboration
documents.

**The P2 experiment** (design report: D. Becker et al., *Eur. Phys. J. A* 54 (2018) 208,
[arXiv:1802.04759](https://arxiv.org/abs/1802.04759)):
- Measures the weak mixing angle **sin²θ_W to 0.14 % relative precision** via the proton's weak
  charge, from the parity-violating asymmetry in elastic e-p scattering at MESA (Mainz).
- **A_PV = −39.94 ppb, target ΔA_PV = 0.56 ppb (1.4 %)** at **⟨Q²⟩ = 4.5×10⁻³ GeV²**; 155 MeV,
  150 µA polarized CW beam on a 60 cm LH₂ target; 10 000 h of running; total elastic-electron rate
  of order **0.1 THz** — this is the "high-rate environment" of the story's first sentence.
- BSM mass-scale sensitivity ~50 TeV.

**BASKET — the backward-angle Micromegas detector** (Mainz P2
[site](https://www.blogs.uni-mainz.de/fb08p2/the-p2-experiment/); design report Sec. 7 for the
physics case):
- **≈ 20 000 Micromegas channels in three detector planes** at large backward angles; the three
  planes reconstruct electron tracks in the 0.6 T solenoid field — *this is exactly the 3-plane
  coincidence configuration measured in July 2026 (Act 3)*.
- **Physics objectives:** determine the **strange-quark contributions to the nucleon
  electromagnetic form factors** and the **effective axial form factor of the proton** — the
  dominant hadronic uncertainties in the sin²θ_W extraction. The backward-angle program improves
  the uncertainty on **G_E^s ×4, G_M^s ×12, and G_A ×10** (design report, Sec. 7).
- Forward-tracker context (useful contrast on the slide): the *forward* Q² tracker is HV-MAPS, and
  the forward ⟨Q²⟩ must be known to **1 %** — tracking is what turns an integrating asymmetry
  measurement into physics.
- ⚠ **Not yet published for BASKET**: spatial-resolution spec, per-plane rate density, timing
  requirement, pad-vs-strip rationale, Saclay/IRFU role. No refereed BASKET paper exists as of
  July 2026 — cite the Mainz web pages + design report, and pull the specs from collaboration
  internal notes. (This also means the talk is presenting largely *unpublished* material — a plus,
  but see the clearance question in §8.)

**Electronics references for the comparison slides:**
- **DREAM**: 64 ch, CSA + shaper + 512-cell SCA analog memory, peaking time 50–900 ns (16
  settings), dead-timeless to 20 kHz external trigger — CLAS12 MVT paper, *NIM A* 957 (2020)
  163423. Full-waveform, externally-triggered → the "understanding" readout.
- **VMM3a**: 64 ch, self-triggered sparse readout, per-hit peak ADC + TDC, peaking times
  25/50/100/200 ns, gains 0.5–16 mV/fC → the "operating" readout. ⚠ canonical citation
  (G. Iakovidis, ATLAS NSW/RD51) to be confirmed before the talk.
- Context/outlook citation if asked about the future: SALSA MPGD readout ASIC,
  [arXiv:2501.10237](https://arxiv.org/abs/2501.10237).

---

## 3. What already exists (asset inventory, 2026-08-10)

| Analysis | Code | State |
|---|---|---|
| SNR config-scan pipeline (VMM + channel level) | `SNR_Analysis/vmm_config_scan_analysis.py` | **Done** for 5 & 15 kHz; report `vmm_analysis_report.html` — needs a presentation pass only |
| Trigger rate / spill structure / occupancy | `SNR_Analysis/vmm_trigger_analysis.py` | **Done** for muon datasets |
| Trigger-coincidence efficiency + per-pad maps | `vmm_detector_efficiency.py` (FanPadDetector, validated REVERSE ordering) | Works on single runs; the rate sweep across config points is **not run** |
| Pad geometry models (Gerber-true) | `FanPadDetector.py`, `SquarePadDetector.py`, `Detector_Mapping/` | **Done**, M3-cross-validated |
| Cosmic-bench QA pipeline (stages 02–18) | `cosmic_bench_analysis/` + `run_p2_pipeline.sh` | **Done** for det1–det4; PDF builders + LaTeX report exist |
| Bench HV / drift scan efficiency | stages 11, 16 + PDF builders | **Done** for det1/det2/det4 |
| Bench timing (waveform TOA) | stage 13 | **Done**, ~46 ns vs trigger |
| Detector-lifetime autopsies | `reports/det_lifetime_autopsy_2026-07`, `reports/p2_efficiency_intermittency_2026-07` | **Done** — det3 drift-foil decoupling, det1 mesh-contact arcing |
| Garfield++ / Magboltz simulation | `gas_studies/`, timing study report | **Done** — physics floor 3–7 ns; beam-gas decision Ar/CF₄/iC₄H₁₀ 88/10/2 (⚠ uncommitted) |
| **SPS beam-test stage battery (DREAM)** | `sps_beam_analysis/` stages 20–31 + condor | **Done** — spectra, alignment, tag-probe efficiency, beam profile, sync QA, commissioning QA, spark QA, timing QA, waveform timing, raw-stream efficiency, LaTeX summary |
| **uRWELL-referenced absolute efficiency** | `sps_beam_analysis/urw_reference/` | **Done** — the v2 plan had no equivalent of this at all |
| **VMM online analysis** | `vmm_qa/` (`vmm_stations.py`, `vmm_efficiency.py`, `vmm_beam_profile.py`) | Stages 0–2 done: trigger-referenced station efficiency, pad map, hit maps, wired into the live watcher |

### 3a. The Act 3 results in hand

| Measurement | Result |
|---|---|
| **Absolute efficiency** vs uRWELL tracks (`highstat_eff_1`, ~1.15 M tracks/station) | **P2_IN 96.5 %, P2_MID 97.1 %, P2_OUT 96.0 %**, ±0.3 % systematic (probe radius + track cut), statistical error negligible |
| 2-of-3 tag-and-probe, P2 planes only (same sub-run) | 92.5 / 96.3 / 92.5 % — 4 points lower and much less uniform on the outer stations, because the tag is smeared by the tagging planes' own 12 mm pads. *That gap is the clearest demonstration of what the external reference buys.* |
| **Spatial residuals** | 3.39 / 3.47 mm rms, identical at all three stations = 12 mm/√12 → binary pad response, **no charge sharing** (92–95 % single-pad clusters) |
| **Timing**, walk-corrected, per station | **P2_MID 15.5, P2_OUT 18.4, P2_IN 22.4 ns**, stable <0.5 ns over 4.5 h → **the 20 ns P2 goal is met at MID/OUT** |
| Timing vs drift V | Plateau at 750–800 V; the 700 V working point is within ~0.5 ns of optimum |
| Timing vs mesh V | Monotonic improvement, **no plateau reached at the ceiling** |
| **ε vs mesh** (constant drift field) | MID 0.806→0.961, OUT 0.857→0.950 over 390→450 V; IN 0.314→0.917 over 370→430 V — **none has reached a plateau** at the top of the scan, which is also the working point |
| **ε vs drift** (mesh 450 V) | 0.157 at zero drift field → 0.896 by 500 V → flattening above ~850 V at 0.974 |
| Frame / mechanical stability | uRWELL→P2 rotation −60°, proper (det +1), no shear; moves <0.02° over 4.5 h and across all 23 scan points |
| Mapping cross-check | The flipped `c_5_top` ribbon on P2_IN was found by the P2 planes against each other **and independently confirmed** through the uRWELL reference *and* through the VMM readout — three chains, same answer |

**Missing / to build:** DREAM↔VMM like-for-like comparison (blocked on the VMM acceptance
denominator — §4 WP-B), bench↔beam ε(HV) overlay, the charging-up-vs-rate separation, pion-run
analysis (still zero code), efficiency-vs-rate sweep across the muon datasets, cross-geometry
comparison plots, the cross-campaign summary figure, **and the slide deck, which does not exist**.

---

## 4. Work packages, re-cut for v3

v2's WP1–WP7 were ordered for a talk whose Act 3 had not happened. The packages below replace
them; the old numbering is given for cross-reference. Priority: **P0** = the talk is weaker
without it; **P1** = strongly wanted; **P2** = cut candidate (§5).

### WP-A — Settle the headline efficiency number *(P0, ~1 day — new in v3)*
`highstat_eff_1` loses 1.3 points monotonically over 2.5 h at constant 4.6 kHz on MID/OUT and
recovers after a 2.5 h pause — but the recovery sub-run changes **both** the pause and the rate,
so charging-up and rate loading are not separated. Every track already carries `t_ns`, so binning
efficiency *within* a sub-run is nearly free and breaks the degeneracy.
- Add a `--time-bins` mode to `urw_reference/urw_p2_efficiency.py`; run over sub-runs 00–05.
- **Blocks every slide that quotes 96–97 %.** Do it first.

### WP-B — DREAM vs VMM on the same detectors *(P0, ~3–4 days, highest risk — was wish-list 5a.1)*
The most-requested plot when two readouts coexist, and the direct justification of the
dual-electronics strategy. **Blocked by a denominator mismatch, not by physics:** VMM efficiency
(IN 0.11 / MID 0.53–0.58 / OUT 0.56–0.64) is a whole-station scalar over *every* trigger while
only connectors c4–c6 of ~10 are instrumented; DREAM's 96 % is per-pad over the illuminated
overlap. **They must never appear on the same axis as they stand.**
- Path: the pad map is done (`vmm_stations.build_pad_table()`) and `vmm_efficiency.coincident_mask()`
  gives the in-time selection; what is missing is the **per-pad efficiency map** with the
  denominator restricted to tracks crossing instrumented pads (`VMM_ONLINE_EFFICIENCY_PLAN.md`
  Stages 3–4).
- **Go/no-go Fri 21 Aug.** Fallback: compare what already *is* comparable at the same station and
  HV — cluster size (VMM 1.11–1.18 vs DREAM 1.05–1.08 pads), ADC/MPV spectra, timing σ, hit maps
  side by side. Still a real slide, and it makes the honest point that VMM timing is
  calibration-limited.
- **Corrected 2026-08-19 (`vmm_dream_matching/TIMING.md`):** VMM σ is *not* BCID-limited. The
  TDC fine time is already applied and toggling it changes σ by <1 ns; quantisation alone is
  22.5/√12 = 6.5 ns. The largest removable term is the **trigger channel itself** — a scintillator
  read through a VMM discriminator — at 33.4 ns against P2_MID's 26.6 ns intrinsic. Referencing to
  the DREAM trigger timestamp via the stream matching (10.4 ns rms) is the real improvement.
- Superseded, kept for the record: VMM σ is pinned at the 22.5 ns BCID quantisation for want
  of time-calibration runs (`qa_config.py: 'calibration': None`); and P2_IN reads low because of
  **gain** (mesh 440 V vs 450 V, still on the steep turn-on), not mapping — a reviewer will
  otherwise read it as a broken station.

### WP-C — Bench ↔ beam ε(HV) overlay *(P0, ~1 day — was wish-list 5a.2)*
Both halves now exist; nobody has overlaid them. Beam: the uRWELL-referenced mesh and drift curves
in `analysis/urw_referenced_efficiency/drift_mesh_scan_1/`. Bench: cosmic stages 11 and 16 for
det1/det2/det4. One figure, same axes.
- **Get the wording right before drawing it:** these are *not* the same detectors, and the gas
  differs. The defensible claim is "the bench reproduces the shape and predicts the working point",
  not "the same number".

### WP-D — Act 3 hero plots *(P0, ~2 days — was wish-list 5b.2/5b.3)*
1. **3-plane track event display** on the true fan-pad geometry — uRWELL track plus three matched
   P2 leading clusters. The picture that makes "this is the P2 tracker" instantly credible.
2. **2D residual plot** showing the pad footprint as a filled square tilted at −60°, annotated
   3.46 mm = 12 mm/√12.
3. **Per-pad efficiency map** in the P2 frame at the working point, beside the uRWELL-frame map.
4. **ε(HV)** curves, mesh and drift, uRWELL-referenced.

### WP-E — SNR result freeze *(P1, ~1 day — was WP1, unchanged and still mostly done)*
SNR vs (sg, snt) matrix per detector, 5 and 15 kHz side by side; one channel-level SNR map per
geometry at best config; the headline *best config = (sg, snt), SNR = X large / Y small*; example
signal-vs-noise ADC spectra at best/worst config. Presentation pass, not analysis.

### WP-F — Rate performance, muons *(P1, ~2 days — was WP2)*
`vmm_detector_efficiency.py` over the matched config points of the `5kHz` and `15kHz` registries →
ε per detector per rate; on-spill pad occupancy at 15 kHz; ADC MPV vs rate.
- **Memory safety:** the registries carry deliberate `n_files` values (10 at 15 kHz, 3 at 5 kHz —
  "~105 MB/file, keep memory low"). Verify and state peak RSS before raising either.

### WP-G — Timing narrative *(P1, ~0.5 day — new in v3)*
Four existing numbers form one slide: bench 46 ns (drift-geometry limited) → beam 15.5/18.4/22.4 ns
per station → Garfield physics floor 3–7 ns → the 20 ns P2 goal, **met at MID/OUT**. Close with the
drift plateau (750–800 V) and the gas decision (Ar/CF₄/iC₄H₁₀ 88/10/2, 2.4× v_d, mesh 410→425 V)
as the forward path. This closes the loop opened by the July timing study.

### WP-H — QA / methodology *(P1, ~0.5 day — was WP5.3)*
The bench found drift-foil HV decoupling and mesh-contact arcing before those detectors reached a
beam. One slide, framed as methodology.

### WP-I — Cross-campaign synthesis *(P0, ~1 day, do last — was WP7)*
One table: detector generation × electronics × environment → efficiency, SNR, timing, spatial
resolution, rate reach. Rows now genuinely available: Nov-2025 large/small + VMM + beam; cosmic
bench det1/det2/det4 + DREAM + cosmics; July-2026 IN/MID/OUT + DREAM + beam; the same three + VMM
+ beam. This is the slide people photograph. Follow with the outlook: what remains before
installation in P2.

### WP-J — Slides, dry run, freeze *(P0 — new in v3, and the real critical path)*
No deck exists. See the timeline in §6.

---

## 5. Wish list → scorecard, and the cut list

### 5a. What the v2 wish list asked of the beam test, and what came back

| v2 item | Outcome |
|---|---|
| 5a.1 DREAM vs VMM back-to-back | Data taken on the same detectors at the same HV. **Analysis blocked on the acceptance denominator** → WP-B |
| 5a.2 HV turn-on on beam, overlaid on the bench | Beam half **done** (mesh + drift, two independent methods). Overlay **not made** → WP-C |
| 5a.3 VMM at the Nov-2025 optimum vs a non-optimal config | Not taken as a dedicated comparison |
| 5a.4 Efficiency map on beam vs bench | Per-pad beam maps **done** at every scan point; the bench comparison rides on WP-C |
| 5a.5 Beam timing vs the bench 46 ns | **Done, and better than hoped**: 15.5–22.4 ns per station, plus σ vs drift and σ vs mesh curves |
| 5b.1 2-of-3 tagging efficiency | **Done** — 92.5 / 96.3 / 92.5 %, *and superseded* by the absolute uRWELL measurement |
| 5b.2 3-plane event display | Ingredients all present; **plot not made** → WP-D |
| 5b.3 Plane-to-plane spatial residuals | **Done** — 3.4 mm, = pad/√12 |
| 5b.4 Plane-to-plane time residuals | **Done** — per-station σ, plus the drift/mesh dependence |
| 5b.5 Track-level efficiency & fake rate vs rate | Partial: accidental slope measured (0.0006 per 10 mm probe radius); no dedicated rate scan |
| 5b.6 Angle scan | Not taken |
| 5b.7 Multi-track separation with pions | Not taken |

The v2 "minimum protected set" was 5a.1, 5b.1, 5b.2. Two of the three are in hand; 5a.1 is WP-B,
the one genuine gap.

### 5b. Cut list — dropped unless the P0 work finishes early

Recorded so that each cut is a decision, not an oversight.

- **Pion runs (v2 WP3).** The highest-uncertainty item in v2, still at **zero code**, on Nov-2025
  data — its purpose was rate reach via local rate density, and the July telescope now carries a
  rate story with real tracks. **Recommendation: timebox to half a day of scoping** (find the run
  numbers in the elog; confirm whether the data is reachable from here or only on the
  `/drf/projets/clas12/` cluster path). If it is not trivially ready, it becomes one outlook line
  and moves to the paper.
- **Geometry comparison (v2 WP4).** Cheap *if* WP-E and WP-F land; otherwise a sentence on the
  Act-1 verdict slide rather than a figure.
- **Resuming `drift_mesh_2d_1`** and the other "nice-to-have additional runs" — no beam before the
  conference, so these are outlook regardless.
- **VMM time calibration** (`VMM_ONLINE_EFFICIENCY_PLAN.md` Stage 5) — blocked on calibration runs
  that do not exist. Outlook line.

---

## 6. Timeline

| When | Work | Milestone |
|---|---|---|
| **Mon 10 – Sun 16 Aug** (T–21 → T–15) | WP-A, WP-C, WP-B starts, WP-E, housekeeping (§8) | Efficiency number settled; bench↔beam overlay exists; SNR frozen |
| **Mon 17 – Sun 23 Aug** (T–14 → T–8) | WP-B finish or fall back (**go/no-go Fri 21**), WP-D, WP-F, WP-G, WP-H, **slides v1** | Every hero plot exists; deck assembled end to end |
| **Mon 24 – Sun 30 Aug** (T–7 → T–1) | WP-I synthesis, **dry run by Tue 26**, backup slides | Content freeze after the dry run — cosmetic changes only |
| **Mon 31 Aug – Fri 4 Sep** | Conference; talk Friday | Rehearse only |

Hard gate: **WP-B's go/no-go is Fri 21 Aug.** If the per-pad VMM efficiency is not converging by
then, take the fallback and stop — the last week is not for analysis.

---

## 7. Slide skeleton (≈15 min talk — confirm the length, it drives this)

Re-weighted from v2: Act 1 down, Act 3 up.

1. Motivation & detector concept (P2 experiment, BASKET objectives from §2, pad-Micromegas
   concept) — 2 slides
2. **Act 1**: SPS Nov-2025 setup, two geometries, VMM3a readout + SNR scan matrix and best-config
   spectra — 2 slides
3. Rate performance (muons) and the geometry verdict — 1 slide
4. **Act 2**: cosmic bench + M3 telescope (note: DREAM-only readout available); pad-tile efficiency
   map, HV turn-on, timing — 2 slides
5. QA / stability methodology — what the bench caught before the beam — 1 slide
6. **Act 3**: the P2-geometry telescope at H4 — setup, uRWELL references, 3-plane event display
   — 1–2 slides
7. Act 3 results: absolute efficiency 96–97 %, residual 3.4 mm, per-pad maps, ε(HV) — 2 slides
8. Act 3 timing: 15.5–22.4 ns per station vs the 20 ns goal, vs Garfield, vs gas — 1 slide
9. DREAM vs VMM on the same detectors — 1 slide
10. Synthesis table + what remains before installation in P2 — 1 slide

**Backup slides:** spill mask / Δt fit / MAD noise methodology; mapping validation vs M3 *and* the
independent uRWELL and VMM confirmations of the flipped `c_5_top` ribbon; det1/det3 lifetime
autopsy; Garfield timing simulation; the gas study; the tag-probe-vs-uRWELL comparison as a
systematics discussion.

---

## 8. Practical notes and open questions

- Every plot destined for the talk gets a `--presentation` style pass (one consistent font/size, no
  debug titles, units on every axis); all pipelines save PNG+PDF — keep the PDFs for slides.
- The `DATASETS` / run-registry pattern is the mechanism for every new campaign entry — one dict
  entry each, nothing else changes.
- **Housekeeping, this week:** `gas_studies/` has four modified files plus untracked `decide.py`
  carrying the 2026-08-01 beam-gas decision — commit per-file. On the DAQ side,
  `run_config_beam.py` holds the corrected VMM trigger channel (VMM 0 ch 44) applied in place but
  **uncommitted**, with its only backup in `/tmp` — that will not survive a reboot.
- The `.tex`/`.pdf` mirrors of this document are still at **v2** and now diverge from it.

**Open questions that change the plan:**

1. **Talk length and title** — this document assumes ≈15 min. The §7 re-weighting depends on it.
2. **Clearance**: is the July beam-test material presentable externally, given that no refereed
   BASKET paper exists and the specs come from collaboration internal notes? Ask now, not in the
   last week.
3. **Is any Nov-2025 pion data reachable from here**, or only on the `/drf/projets/clas12/` cluster
   path in the registries? That answer settles the pion cut immediately.
4. **Which efficiency number goes on the slide** — this is WP-A's output, and it is the first thing
   an audience will write down.
