# Conference Roadmap — P2 Basket Pad-Micromegas: from SPS Beam Test to Cosmic-Bench Characterization

*Status: roadmap v2, 2026-07-12. **Conference: 31 August 2026 — talk: 4 September 2026** (≈7 weeks
from now to the conference).*

---

## 1. The story to tell

One sentence: **"We are developing pad-Micromegas trackers for the P2 experiment; we optimized the
VMM3a front-end and validated two geometries at the SPS, characterized the new detector generation
on a cosmic bench with full external tracking, and now return to the SPS with three detectors
mounted exactly as they will be installed in P2 — same plane distances, same electronics — to
measure three-plane track coincidences."**

The talk has a natural three-act structure:

### Act 1 — Beam test at SPS (Nov 2025): geometry + electronics optimization
- Two detector geometries under test simultaneously with VMM3a readout:
  - **P2 Large detector** — fan-pad geometry, 256 pads (VMMs 12–15)
  - **P2 Small detectors 1 & 3** — compact pad geometry (VMMs 8–11)
- **SNR optimization**: full VMM3a shaping scan (gain `sg` = 3.0/4.5/6.0 mV/fC × peaking time
  `snt` = 50/100/200 ns), signal/noise run pairs via neighbour-triggering, MAD noise baseline + MPV
  signal → SNR per VMM and per channel. → *"We found the optimal operating point of the electronics."*
- **Rate performance**: identical config points taken at ~5 kHz and ~15 kHz muon rates, plus **pion
  runs** with localized, much higher rate density. → *"We probed the detector where it will actually
  have to live."*
- Trigger-based efficiency with accidental subtraction (Δt peak fit + sideband) and per-pad
  efficiency maps on the real Gerber geometry.

### Act 2 — Cosmic bench (2026): new detector generation, precision characterization
- New P2 detectors (P2_1 “det1”, P2_2 “det2”), 256-pad fan geometry, read out with DREAM FEUs behind
  the **M3 beam telescope** → track-projected, *geometry-resolved* measurements a beam trigger alone
  cannot give. (**DREAM was the only readout option available on the cosmic bench** — this is why
  the bench characterization is DREAM-based and why the beam test must carry both electronics; see
  Act 3.) The bench delivers:
  - pad-level and sliding-window **efficiency maps** on the true pad footprint
  - **HV turn-on curve** (mesh 345–420 V, ε 0.40 → 0.66 measured so far)
  - **timing** from waveform TOA (~46 ns vs trigger, drift-geometry limited)
  - **stability**: HV spark monitoring/veto, per-pad spark flagging, pedestal QA, long-run duty-cycle
    tracking (this is how the det1 intermittent-gain problem was *found* — that is a selling point:
    the bench methodology catches real detector pathologies)
- → *"Before the next beam test, every detector is fully mapped, efficiency-calibrated and
  stability-vetted on the bench."*

### Act 3 — Beam test at SPS (summer 2026): the P2 telescope in its final configuration
This is no longer a generic "outlook" — it is a **dress rehearsal of the P2 tracker itself**:

- **Three detectors in a row on the SPS beam, mounted as they will be installed in P2**: same
  plane-to-plane distances, same mechanics, same electronics chain as the experiment. The
  measurement is the experiment's measurement: **coincidences of track segments among the three
  planes** — exactly how BASKET will reconstruct backward-scattered electron tracks in the P2
  solenoid field for the strange/axial form-factor program (see §2).
- **Two front-end electronics on the same detectors: DREAM and VMM3a.** On the cosmic bench, DREAM
  was the only available option, so the full bench characterization (efficiency maps, HV turn-on,
  timing, stability) exists in DREAM terms. The beam test is the *first opportunity to read the
  same, fully-mapped detectors with both electronics* — it closes the loop between the bench
  campaign (DREAM) and the Nov 2025 VMM config-scan optimum, and yields a like-for-like
  DREAM ↔ VMM comparison on identical hardware.
- The complementarity is the selling point: **DREAM (full waveforms, external trigger) for
  understanding; VMM3a (self-triggered, high-rate) for operating.** Same detectors, two readouts,
  and now the P2 geometry in between.
- Bench-validated detectors + SPS-optimized VMM configuration → measure the thing that matters:
  **three-plane tracking efficiency and stability vs particle rate at the optimal working point**,
  with known geometry and known per-pad response.

---

## 2. What P2 needs from these detectors (literature anchors for the motivation slide)

Verified against the published literature (2026-07-12). Use these numbers on the motivation slide;
items flagged ⚠ have no public citation yet and must come from collaboration documents.

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
  coincidence configuration of the summer 2026 beam test (Act 3)*.
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
  internal notes. (This also means the talk is presenting largely *unpublished* material — a plus.)

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

## 3. What already exists (asset inventory)

| Analysis | Code | State |
|---|---|---|
| SNR config-scan pipeline (VMM + channel level) | `SNR_Analysis/vmm_config_scan_analysis.py` + modules | **Done** for 5 kHz & 15 kHz; report `vmm_analysis_report.html` |
| Trigger rate / spill structure / occupancy | `SNR_Analysis/vmm_trigger_analysis.py` + `trigger_rate_analysis_report.md` | **Done** for muon datasets |
| Trigger-coincidence efficiency + per-pad maps | `vmm_detector_efficiency.py` (FanPadDetector, validated REVERSE ordering) | Works on single runs; needs a config-scan / rate sweep |
| Pad geometry models (Gerber-true) | `FanPadDetector.py`, `SquarePadDetector.py`, `Detector_Mapping/` | **Done**, M3-cross-validated |
| Cosmic-bench full QA pipeline (stages 02–13) | `cosmic_bench_analysis/` + `run_p2_pipeline.sh` | **Done** for det1 (4 long runs), det2 (1 run); PDF builders exist |
| HV-scan efficiency | stage 11 + `build_hv_scan_pdf.py` | Done for det1 345–420 V (365 V file to re-fetch) |
| Timing (waveform TOA) | stage 13 | Working on det2, **not committed yet** |
| Structural bias audit | stage 12 | Done (knob plateaus, estimator checks) — backs up systematic error bars |

**Missing / to build:** pion-run analysis (nothing in the codebase mentions pions), efficiency-vs-rate
sweep across the muon datasets, cross-geometry comparison plots, det2 HV scan, **inter-plane track
matching for the 3-plane beam test** (WP6 — no code exists for detector-to-detector coincidences yet;
the M3-projection logic in stages 03/05/06 is the closest starting point), and the final
cross-campaign summary figures.

---

## 4. Analysis roadmap (work packages, in priority order)

### WP1 — SNR result finalization *(mostly done — polish only)*
1. Freeze the SNR summary: **SNR vs (sg, snt)** matrix per detector, 5 kHz and 15 kHz side by side.
2. One channel-level SNR map per geometry at the best config (shows uniformity, not just averages).
3. Extract the headline number: *best config = (sg, snt) with SNR = X on large, Y on small.*

**Hero plots:** SNR heat matrix (sg × snt) per detector; example signal-vs-noise ADC spectra at
best/worst config (the "why this matters" plot); channel-level SNR pad map.

### WP2 — Rate performance, muons *(1–2 weeks of work; reuses existing pieces)*
1. Run `vmm_detector_efficiency.py` over the matched config points of the 5 kHz **and** 15 kHz
   run tables → efficiency per detector per rate.
2. Occupancy/hit-rate per pad on-spill (from trigger-analysis machinery) → rate-density maps.
3. Check for rate-dependent effects: Δt peak shape vs rate, ADC (gain) spectra on- vs off-spill,
   efficiency vs instantaneous rate within the spill.

**Hero plots:** **efficiency vs trigger rate** (one point per dataset, per geometry — this is the
money plot of Act 1); on-spill pad occupancy map at 15 kHz; ADC MPV vs rate (gain stability).

### WP3 — Pion runs *(new analysis — start early, highest uncertainty)*
1. Locate pion run numbers/conditions in the elog + run tables; add a `"pion"` entry to the
   `DATASETS` registries (the registry pattern means nothing else changes).
2. Same chain as WP2: spill mask → occupancy map → coincidence efficiency, but now the key variable
   is **local rate density** (hits/pad/s in the beam-spot core vs periphery).
3. Efficiency and MPV as a function of local pad rate — one curve combining muon + pion points
   extends the rate reach by an order of magnitude in a single figure.

**Hero plots:** pion beam-spot occupancy map (localized spot on the pad geometry — visually
striking); **efficiency (and/or MPV) vs local pad rate**, muon + pion points on one axis.

### WP4 — Geometry comparison *(synthesis of WP1–3, cheap once they exist)*
- Large fan-pad vs small pad detectors at identical (sg, snt, rate): SNR, efficiency, occupancy per
  unit area. State the verdict explicitly: which geometry carries forward, and why the new cosmic-
  bench detectors look the way they do. This is the hinge between Act 1 and Act 2.

**Hero plot:** two-panel same-scale comparison (SNR map / efficiency map, large vs small).

### WP5 — Cosmic-bench characterization of the new detectors *(pipeline exists — needs det2 completion)*
1. Finish det2: HV scan (stage 11 equivalent for det2), sliding efficiency map, timing (commit
   stage 13), longer statistics run.
2. Produce the *presentation* versions: pad-tile efficiency map at best HV, HV turn-on curve
   det1 + det2 overlaid, TOA resolution plot.
3. Stability narrative (optional but strong): duty-cycle/long-run monitoring plot showing how the
   bench catches gain instability — frame as methodology, not as a failure.

**Hero plots:** pad-tile efficiency map with M3 tracks (the "we see every pad" plot); HV turn-on
ε(HV); timing residual distribution with σ ≈ 46 ns; one spark/stability monitoring time series.

### WP6 — Three-plane beam test analysis *(new — depends on the beam-test data arriving in time)*
1. Extend the existing pipelines to the 3-plane setup: the `DATASETS` registry + per-detector
   config pattern already supports multiple FEUs/detectors; the new ingredient is **inter-plane
   track matching** (a 3-plane version of the M3-projection logic, but now the P2 detectors *are*
   the telescope).
2. First-priority results (see the wish list in §5): per-plane efficiency from 2-of-3 tagging,
   plane-to-plane timing, DREAM vs VMM on the same detector.
3. Even a *few hours* of good 3-plane coincidence data is enough for the conference: one event
   display of a 3-plane track + one 2-of-3 efficiency number per plane already makes Act 3 real.

### WP7 — Cross-campaign synthesis *(the conclusion slide — do last, 2–3 days)*
- One table/figure: detector generation × electronics × environment → efficiency, SNR, timing,
  rate reach. This is the single slide people will photograph.
- Explicit outlook: the 3-plane P2-geometry telescope with bench-validated detectors + VMM at the
  WP1 optimal config, and what remains before installation in P2.

---

## 5. Beam-test wish list — measurements & data that would be *nice to have* for the talk

Ordered by impact-per-beam-hour. Each item states what it buys for the conference.

### 5a. Comparison measurements (bench ↔ beam, DREAM ↔ VMM)

1. **Same detector, same HV, DREAM vs VMM back-to-back runs** (even one detector, one config):
   efficiency, ADC/MPV spectrum, cluster (pad-multiplicity) size, timing residual. *Buys:* the
   like-for-like electronics comparison slide — the single most requested plot when two readouts
   coexist in a project, and the direct justification of the dual-electronics strategy.
2. **HV turn-on curve on beam with DREAM** for at least one detector: overlay on the cosmic-bench
   ε(HV) from stage 11. *Buys:* "the bench predicts the beam" — validates the whole Act 2
   methodology in one figure.
3. **VMM at the Nov 2025 optimal (sg, snt) vs one non-optimal config** on the new detectors:
   *Buys:* closes the loop on Act 1 — the config-scan optimum demonstrably transfers to the new
   detector generation.
4. **Efficiency map on beam vs efficiency map on bench** for the same detector (beam spot scanned
   or wide beam): *Buys:* pad-level reproducibility, and cross-checks the dead zones (pillars,
   connector 10 class of problems) in an independent environment.
5. **Timing with beam tracks vs the bench 46 ns**: same TOA analysis (stage 13) on beam data,
   where track geometry is known and drift-path spread can be controlled. *Buys:* separates the
   intrinsic timing from the drift-geometry contribution identified in the timing study — and the
   Garfield physics-floor comparison (3–7 ns) becomes a "measured vs simulated" slide.

### 5b. Physics / telescope measurements (the P2-configuration measurements)

1. **Three-plane coincidence tracking with 2-out-of-3 efficiency tagging**: use track segments in
   two planes to predict the crossing point in the third → unbiased per-plane, per-pad efficiency
   *without any external trigger* — exactly the in-situ efficiency-monitoring scheme available in
   the real experiment. *Buys:* the headline Act 3 number and the method P2 itself will use.
2. **One 3-plane track event display** on the true pad geometry. *Buys:* the picture that makes
   the "this is the P2 tracker" claim instantly credible.
3. **Plane-to-plane spatial residuals** (predicted vs measured pad, in the fan-pad coordinates
   radius/φ): *Buys:* measured effective spatial resolution of the pad readout for tracking —
   the number that feeds any Q²-resolution statement.
4. **Plane-to-plane time residuals**: time coincidence window needed to match segments into
   tracks. *Buys:* the track-matching timing budget in the P2 rate environment — directly
   motivates (or relaxes) the 20 ns timing goal from the timing study.
5. **Track-level efficiency & fake rate vs beam rate** (muon rates, then pions for local density):
   coincidence rate, accidental/combinatorial track rate vs single-plane rates. *Buys:* the
   rate-reach money plot, now at *track* level rather than single-detector level — the quantity
   that actually limits P2 operation.
6. **Angle scan (if the setup allows tilting or off-axis positions)**: efficiency and cluster size
   vs incidence angle, since scattered electrons in P2 cross the planes at non-normal angles.
   *Buys:* realism of the acceptance; also feeds cluster-size systematics for tracking.
7. **Multi-track separation with pions/high intensity**: two tracks in the same coincidence
   window — can the three planes disentangle them? *Buys:* occupancy/ambiguity argument for the
   P2 rate environment.

*Fallback note:* items 5a.1, 5b.1 and 5b.2 are the minimum set worth protecting beam time for; all
other items degrade gracefully into "shown at the next opportunity" outlook material.

---

## 6. Suggested timeline (conference = 31 Aug 2026, talk = 4 Sep 2026)

| When | Milestone |
|---|---|
| **Now → end July (T–5 weeks)** | WP3 first pass (pion data located and read); WP2 efficiency-vs-rate sweep running; det2 HV scan taken on the bench; **beam-test preparation: run list + trigger plan for the §5 wish list agreed with the team** |
| **Early → mid August (T–4 → T–3 weeks)** | WP1 frozen; WP2/WP3 result plots v1; det2 characterization complete; beam-test data taken → WP6 first pass (event display + 2-of-3 efficiency) |
| **Mid August (T–3 → T–2 weeks)** | WP4 + WP5 presentation plots finalized (consistent style, labels, fonts); WP6 result plots; WP7 synthesis figure |
| **T–2 → T–1 weeks (17–24 Aug)** | Slides assembled; internal dry run; only cosmetic plot changes after this |
| **Last week (24–31 Aug)** | Freeze. Backup slides: methodology (spill mask, Δt fit, MAD noise), det1 stability study, mapping validation vs M3, Garfield timing simulation |

---

## 7. Slide skeleton (≈15 min talk)

1. Motivation & detector concept (P2 experiment, BASKET Micromegas objectives from §2,
   pad-Micromegas concept) — 1–2 slides
2. **Act 1**: SPS setup photo + two geometries + VMM3a readout — 1 slide
3. SNR optimization: scan matrix + best-config spectra — 2 slides
4. Rate performance: efficiency vs rate (muons), pion beam-spot + efficiency vs local rate — 2 slides
5. Geometry verdict — 1 slide
6. **Act 2**: cosmic bench + M3 telescope setup (note: DREAM-only readout available) — 1 slide
7. New detectors: pad-tile efficiency map, HV turn-on, timing — 2 slides
8. Stability & QA methodology (spark veto, pedestal, long-run monitoring) — 1 slide
9. **Act 3**: 3-plane P2-geometry beam test — setup photo, 3-plane track event display,
   2-of-3 efficiency, DREAM vs VMM comparison — 1–2 slides
10. Synthesis table + what remains before installation in P2 — 1 slide

---

## 8. Practical notes

- Every plot destined for the talk should get a `--presentation` style pass (single consistent
  font/size, no debug titles, units on every axis); both pipelines already save PNG+PDF, so keep
  the PDF versions for slides.
- The `DATASETS` registry pattern (SNR, trigger, cosmic pipelines) is the mechanism for every new
  campaign entry — pion runs, det2 HV scan — one dict entry each, nothing else changes.
- Open data questions to resolve early: pion run list/conditions (elog), the empty 365 V m3 file
  (re-fetch from `rays_daplxa`), whether det1 recovers enough to appear as a second bench detector.
