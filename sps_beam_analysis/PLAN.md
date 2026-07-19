# SPS beam-test analysis — adaptation plan

Adapting `cosmic_bench_analysis` for the P2 telescope at the SPS beam test.
Key differences from the bench:

| | Cosmic bench | SPS beam test |
|---|---|---|
| Trigger | M3 telescope (self/coincidence) or Fe55 self-trigger | **External scintillator trigger into the TCM** (`Sys DaqRun Trig Ext`, like Nov-2025 `Beam.cfg`) |
| Reference tracking | M3 (2× X/Y strip planes) | **None** (no M3). Candidates: mutual P2↔P2 tagging; BANCO ALPIDE telescope if its data can be synced (open question) |
| Source | cosmics (broad angle) / Fe55 (5.9 keV photopeak) | beam particles (~parallel, MIP-like → Landau, high in-spill rate) |
| Time structure | uniform | **SPS spill structure** (SFTPRO slow extraction, ~4.8 s spills) |

## 1. Layout: this directory (`sps_beam_analysis/`)

Sibling of `cosmic_bench_analysis`, NOT a copy of it. Shared pad-level
infrastructure is **imported** from `cosmic_bench_analysis` (path shim in
`sps_config.py`), so fixes land in one place:

- `p2_io.py` — streaming combined-hits access (unchanged)
- `p2_mapping.py` — Gerber pad map + FEU/connector wiring (unchanged; wiring
  comes from each run's `run_config.json` as now)
- `p2_channel_qa.py`, `p2_sparks.py` (HV spark veto), `p2_waveforms.py` — unchanged
- `p2_qa_config.py` internals reused where possible, but this directory gets its
  own **`sps_config.py`** (run registry, DATA/ANALYSIS roots, beam-specific
  knobs: spill gating, trigger config, telescope geometry z-positions from the
  beam-line survey — the run_config `det_center_coords` TODO-SPS)

Output tree mirrors the bench convention so `build_*_pdf.py` and the DAQ GUI
Analysis tab work as-is:
`<Analysis_SPS>/<detN>/<run>/<sub_run>/<stage>/...`

## 2. Carries over with config only (no code changes)

- `02_map_validation.py` — mapping sanity vs data
- `05_detector_deep_qa.py` — noise/pathology QA (detector-local)
- `09_pedestal_qa.py` — pedestal stability
- `13_timing_waveforms.py` — waveform shape/timing (sample window differs; knob)

## 3. Does NOT apply (M3-based — leave behind)

`03_m3_alignment`, `04_m3_reference_qa`, `06_efficiency_maps`,
`10_efficiency_map_sliding`, `12_validation` (M3-referenced parts),
`16_drift_scan_efficiency`, `m3_comparison`, `m3_signal_diagnostic`.

## 4. Adapt (new scripts, modeled on existing ones)

1. **`20_beam_spectra.py`** (from `18_fe55_spectra.py`): per-sub-run cluster
   charge spectrum, but **Landau MPV** fit instead of Gaussian photopeak;
   MPV / resolution / rate vs HV for scans. Same HV-settle cut + spark veto.
2. **`21_telescope_align.py`** (replaces `03_m3_alignment`): mutual alignment
   of the P2 detectors from cluster-position correlations of shared-trigger
   events (translation + rotation per detector pair; beam is ~parallel so no
   track fit needed — median residual minimisation).
3. **`22_tag_probe_efficiency.py`** (replaces `06`/`11` efficiency): with no
   external tracker, efficiency by tag-and-probe — event tagged by a clean
   cluster in the OTHER detector(s) (+ trigger implies beam particle), probe =
   detector under test, hit search in a radius around the tagged position.
   Honest caveat baked into the outputs: this is efficiency *relative to the
   tag selection*; geometric acceptance of the overlap region only. Per-pad
   efficiency maps like `06` where statistics allow.
   - If BANCO tracking becomes available: add an absolute-efficiency variant
     with real extrapolated tracks (separate script; needs event sync scheme).
4. **`23_beam_profile.py`** (new): beam-spot hit maps per detector, profile
   vs sub-run/HV, rate vs time (spill structure — fold with the DAQ
   `beam_monitor` CSVs once NXCALS access exists), in-spill vs out-of-spill
   occupancy, pile-up indicator (hits/event vs instantaneous rate).
5. **`24_event_sync_qa.py`** (new, small but critical): both FEUs see the same
   TCM trigger stream — verify event-ID/timestamp alignment between FEUs per
   sub-run (guards every correlation-based stage above against off-by-N
   event slips).
6. **`build_beam_pdf.py`** (from `build_hv_scan_pdf.py`): one summary PDF per
   run.

## 5. Suggested implementation order

1. `sps_config.py` + registry entry for a banco test run (can dry-run on the
   Fe55 data structure today — same combined-hits format)
2. `24_event_sync_qa` → `21_telescope_align` (foundations for everything)
3. `20_beam_spectra` (quick win, mostly a fit swap)
4. `22_tag_probe_efficiency`
5. `23_beam_profile` (+ beam-monitor folding when NXCALS is granted)
6. `build_beam_pdf`

Steps 1–3 are useful even before the beam (validate against Fe55/cosmic data);
4–5 need real beam geometry to be meaningful.

## 6. Open questions (answers change the plan)

- **BANCO telescope**: will its ALPIDE data be available + synchronizable with
  the DREAM event stream (shared trigger/busy)? If yes, absolute efficiency
  and resolution become possible and stage 22 gets a big brother.
- **Telescope composition**: just P2_OUT + P2_MID, or a third station
  (P2_IN?) — tag-and-probe becomes much stronger with ≥3 planes
  (majority tagging, unbiased probe).
- **Trigger**: which scintillators, coincidence logic, rate — needed for
  dead-time/pile-up interpretation in stage 23.
- Beam: particle type/momentum per period (Landau expectations, multiple
  scattering between stations for the alignment tolerances).
