# Conference Roadmap — P2 Basket Pad-Micromegas: from SPS Beam Test to Cosmic-Bench Characterization

*Status: draft roadmap, 2026-07-11. Anchor the timeline to the actual conference date once fixed.*

---

## 1. The story to tell

One sentence: **"We are developing pad-Micromegas detectors for high-rate environments; we validated two
geometries and optimized the VMM3a front-end at the SPS, we are now characterizing the next detector
generation on a cosmic bench with full tracking, and we return to beam with VMM to close the loop."**

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
  cannot give:
  - pad-level and sliding-window **efficiency maps** on the true pad footprint
  - **HV turn-on curve** (mesh 345–420 V, ε 0.40 → 0.66 measured so far)
  - **timing** from waveform TOA (~46 ns vs trigger, drift-geometry limited)
  - **stability**: HV spark monitoring/veto, per-pad spark flagging, pedestal QA, long-run duty-cycle
    tracking (this is how the det1 intermittent-gain problem was *found* — that is a selling point:
    the bench methodology catches real detector pathologies)
- → *"Before the next beam test, every detector is fully mapped, efficiency-calibrated and
  stability-vetted on the bench."*

### Act 3 — Outlook: next beam test with VMM
- Bench-validated detectors + SPS-optimized VMM configuration → next campaign measures the thing
  that matters: **efficiency and stability vs particle rate at the optimal working point**, with
  known geometry and known electronics response.
- The cross-electronics angle is a strength, not a complication: **DREAM (full waveforms, external
  tracking) for understanding, VMM3a (self-triggered, high rate) for operating.** Same detectors,
  two complementary readouts, two beam environments.

---

## 2. What already exists (asset inventory)

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
sweep across the muon datasets, cross-geometry comparison plots, det2 HV scan, and the final
cross-campaign summary figures.

---

## 3. Analysis roadmap (work packages, in priority order)

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

### WP6 — Cross-campaign synthesis *(the conclusion slide — do last, 2–3 days)*
- One table/figure: detector generation × electronics × environment → efficiency, SNR, timing,
  rate reach. This is the single slide people will photograph.
- Explicit outlook: next beam test = bench-validated det2 (+det1 if recovered) + VMM at the WP1
  optimal config + rate scan informed by WP2/3.

---

## 4. Suggested timeline (back-count from the conference)

| When | Milestone |
|---|---|
| **Now → T–8 weeks** | WP3 first pass (pion data located and read); WP2 efficiency-vs-rate sweep running; det2 HV scan taken on the bench |
| **T–8 → T–5 weeks** | WP1 frozen; WP2/WP3 result plots v1; det2 characterization complete; commit stage 13 |
| **T–5 → T–3 weeks** | WP4 + WP5 presentation plots finalized (consistent style, labels, fonts); WP6 synthesis figure |
| **T–3 → T–1 weeks** | Slides assembled; internal dry run; only cosmetic plot changes after this |
| **T–1 week** | Freeze. Backup slides: methodology (spill mask, Δt fit, MAD noise), det1 stability study, mapping validation vs M3 |

---

## 5. Slide skeleton (≈15 min talk)

1. Motivation & detector concept (P2 basket pad-Micromegas, target environment) — 1–2 slides
2. **Act 1**: SPS setup photo + two geometries + VMM3a readout — 1 slide
3. SNR optimization: scan matrix + best-config spectra — 2 slides
4. Rate performance: efficiency vs rate (muons), pion beam-spot + efficiency vs local rate — 2 slides
5. Geometry verdict — 1 slide
6. **Act 2**: cosmic bench + M3 telescope setup — 1 slide
7. New detectors: pad-tile efficiency map, HV turn-on, timing — 2 slides
8. Stability & QA methodology (spark veto, pedestal, long-run monitoring) — 1 slide
9. **Act 3**: synthesis table + next beam test plan — 1–2 slides

---

## 6. Practical notes

- Every plot destined for the talk should get a `--presentation` style pass (single consistent
  font/size, no debug titles, units on every axis); both pipelines already save PNG+PDF, so keep
  the PDF versions for slides.
- The `DATASETS` registry pattern (SNR, trigger, cosmic pipelines) is the mechanism for every new
  campaign entry — pion runs, det2 HV scan — one dict entry each, nothing else changes.
- Open data questions to resolve early: pion run list/conditions (elog), the empty 365 V m3 file
  (re-fetch from `rays_daplxa`), whether det1 recovers enough to appear as a second bench detector.
