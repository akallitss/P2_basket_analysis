# Handoff: VMM cosmic efficiency - calculation fixes and what the silent muons are (22 Nov 2026)

Single entry point for the next session. Follows `HANDOFF_analysis_20261120.md` (rules, funnel, gain scan, products,
re-run recipe: still valid) and supersedes it where they differ.

## 0. Rules (unchanged + new)

- All analysis on **lxplus** (HTCondor for decoding / loops); DAQ box read-only; never touch HV or the DAQ.
- Memory-cap every process: `systemd-run --user --scope -q -p MemoryMax=3G bash -c 'source pipeline/env.sh; ...'`.
- Run numbers / configurations: read the logbook, never guess. Pulser configs (`_tpA`, `_tpB`) excluded by name.
- One commit per file, no Co-Authored-By line (repo convention).
- **Efficiency definition still ON HOLD** (user). The funnel numbers below are the current best estimate, not a final
  definition.
- **Chambers on the stand** (user, 20 Nov): P2_IN, P2_MID, P2_OUT are all **Saclay-bulked (det1-4), insulation mask
  V1** (0.8 mm pillars on a 4 mm grid + 5 big 6.15 mm, **3.78 %** of the pad area). P2_IN is det3, NOT the CERN det5
  (det5 = CERN mask, 0.5 mm on 2 mm = 4.87 %, not on the stand). Pillar reference = the cosmic-bench (DREAM) analysis,
  `cosmic_bench_analysis/p2_qa_config.py: MASK_GBR_PATH`.
- User preference: use the **measured** pointing widths of this stand, not the ~5 mm SPS recollection.
- Do NOT put `--pointing` into the production tag dumps (`efficiency_tp.py` in the lxplus loop): the margin drops tags
  that the COUNT estimator (`channel_efficiency.py`) needs for its density.

## 1. What changed (20-22 Nov)

### 1a. Search radius scaled to the pointing resolution (calculation fix)
The prediction for P2_MID is an interpolation (tags IN+OUT, weights 0.5/0.5); for IN and OUT an extrapolation
(weights 2/-1): sigma_res = 1.22 s (MID) vs 2.45 s (IN/OUT). Measured (fit of the nearest-probe-hit distance, 2-D
Gaussian + chance; REF): **MID 4.7 mm, IN 8.5, OUT 8.4 mm** (ratio 1.8, model 2) -> s ~ 3.9 mm per chamber.
The old fixed 15 mm kept 99 % of real MID hits but only ~79 % for IN/OUT: the IN/OUT "far" excess was the pointing.
- Radius = 3.03 sigma_res (99 %): **IN 25.9, MID 14.4, OUT 25.5 mm**.
- Edge margin = 2 sigma_point inside the outline: **IN 15.6, MID 5.5, OUT 15.4 mm** (removes edge tracks whose true
  point may be outside; they were mostly found and partly "silent").
- `pointing.json`: `$W/analysis/muon_funnel_pointing/pointing.json` (fit on the first 7 REF runs).

### 1b. Unphysical tracks rejected (22 Nov)
The tag line must cross both scintillator paddles (the coincidence requires it). Line extrapolated to each paddle
face, error s*hypot(1-f, f) (+) 10 mm paddle-position systematic; rejected if it misses a paddle by > 3 sigma
(`--paddle-cut 3 --paddle-sys 10`). Removes 5-10 % of the tagged tracks, 75-82 % of them were FOUND -> efficiency
+0.1..+0.5 points, silent unchanged. Paddle positions = 4 Nov fit (`paddles_fit_20261104.json`), not re-measured.

### 1c. Current numbers (92 runs 268-371, REF = 14 runs; all cuts: pointing radius + margin + paddle cut)

| REF | fixed 15 mm (same runs) | **now: found** | silent | far | late | found / (1 - 3.78 %) |
|---|---|---|---|---|---|---|
| P2_IN | 71.9 % | **85.2 %** | 9.4 % | 4.9 % | 0.5 % | 88.6 % |
| P2_MID | 79.9 % | **80.2 %** | 13.1 % | 5.9 % | 0.8 % | 83.3 % |
| P2_OUT | 74.4 % | **87.8 %** | 8.4 % | 3.4 % | 0.4 % | 91.3 % |

Chance-found (window +4 us): ~0.01 %. Same improvement in every gain setting; MID moves < 0.5 points everywhere.

### 1d. The silent muons (P2_MID ~9 points beyond the pillars) - what they are NOT
- **Whole chamber silent**: 12.3 of MID's 13.0 % (REF dumps) have no hit on ANY chip within +-1 us.
- **Not the time window**: 140-310 ns vs -200..+1500 ns recovers < 1.5 points (`time_window_check.py`, runs 268-269).
- **Not a geometric hole**: silent flat vs tan(theta) (MID 14.0 / 12.1 / 12.4 / 11.6 % for 0-0.4); a pillar /
  dead-area hole would fall steeply by tan(theta) ~ 0.2 (4 mm drift gap).
- **Chamber-wide**: median pad silent IN 7.7 %, **MID 11.3 %**, OUT 6.5 % (REF, >= 15 muons/pad). Pad spread above
  statistical (chi2/ndf 2.3-3.3) -> some genuinely bad pads (IN/OUT tails), but MID's whole distribution is shifted.
  Small local effects: within 12 mm of a dead/masked pad 21-23 % silent (~0.3 points); big 6.15 mm pillars visible
  in MID (35 % silent within 8 mm, < 0.2 points). Worst pads repeat across chambers on the same channel numbers
  (ch 62 at (418,139) mm on MID and OUT; ch 42 at (330,115)): follow up.
- **Not fake tags**: silent flat vs tag pulse height at REF (all configs: rises, because high gain = larger pulses
  and more silent).
- Gain-scan extrapolation to 0 fC: MID silent -> 3-7 %, OUT -> ~1 % (linear fit over chips; a hint that MID's
  excess responds to threshold).
=> Working conclusion: **signal below threshold over the whole P2_MID** (lower gas gain or higher effective
threshold). Deciding tests: mesh-HV scan of MID with REF electronics (user decides HV), or DREAM on MID.

## 2. Code (branch `decoder-tick-225`, merged into main 22 Nov - see section 5)

| script | what |
|---|---|
| `pointing_resolution.py` | model + fit of sigma_res per probe and group -> `pointing.json` (radius, margin) |
| `muon_funnel.py` | `--pointing JSON` (per-probe radius, edge margin; `margin_cut`), `--margin-live`, `--paddle-cut K --paddle-sys MM --sigma-ch MM` (class `unphysical`, `unphys_cut`) |
| `muon_funnel_summary.py` | sums `margin_cut*`, `unphys_cut*`, keeps `probe_r`, `margin_mm` |
| `muon_funnel_job.sh` + `muon_funnel_pointing.sub`, `muon_funnel_physical.sub` | condor; `OUTDIR` / `EXTRA` from the environment; skips runs already done |
| `efficiency_tp.py` | `--pointing` (radius + margin, same as the funnel); bit-identical without it; deployed to `$W/pipeline/cva/` |
| `time_window_check.py` | hit-time spectrum of each probe around the scintillator for tagged muons |
| `silent_tests.py` | silent vs angle, per pad (spread, dead-pad / big-pillar distance), vs tag ADC, whole-chamber +-1 us; per config |
| `note/muon_funnel_page.py` | args: `SUMMARY NOTE_HTML GAIN_SCAN_JSON POINTING_JSON OLD_SUMMARY_JSON SILENT_TESTS_JSON`; sections 3b-3e, 5a, 5d (inline SVG illustrations) |

`$W = /eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus/weekend_20261119`

Outputs on EOS:
- `$W/analysis/muon_funnel_pointing/` - pointing.json, run*.json (radius + margin, no paddle cut), summary.json,
  gain_scan_out/, time_window_268_269.json, scripts copies (silent_tests.py, efficiency_tp.py, rerun_eff.py)
- `$W/analysis/muon_funnel_physical/` - **current**: run*.json (radius + margin + paddle cut), summary.json (92 runs),
  `summary_fixed15_same_runs.json` (fixed 15 mm, same runs, for before/after), gain_scan_out/, silent_tests.json
- `$W/analysis/muon_funnel/` - the fixed-15 mm funnel (kept up to date, 92 runs)

### Re-run for new runs (the loop keeps dumping)
```
ssh lxplus 'cd muon_funnel && condor_submit muon_funnel_physical.sub && condor_submit muon_funnel.sub'
# when done, build both summaries on the SAME run list (copy the common run*.json to /tmp dirs, as on 22 Nov):
#   muon_funnel_summary.py DIR OUT.json   for both; then
python3 analysis/gain_scan/gain_scan_eff.py analysis/muon_funnel_physical/summary.json results/current analysis/muon_funnel_physical/gain_scan_out
python3 analysis/muon_funnel_pointing/silent_tests.py --tags tags --chunks pseudo/chunks.json \
    --pointing analysis/muon_funnel_pointing/pointing.json --min-pad 15 --out analysis/muon_funnel_physical/silent_tests.json
# locally (copy the jsons to a scratch dir first):
python3 note/muon_funnel_page.py summary_phys.json NOTE_HTML gain_scan_phys.json pointing.json summary_old.json silent_tests.json
# then in cern-site: make hub; commit; scripts/deploy.sh; git push
```
Careful: summaries made at different times differ in run count (new runs arrive): always compare on the same runs.

## 3. Products
- CERN note (unlisted, SSO), updated 22 Nov, cern-site commits a9fb6df..aefb192, deployed and pushed:
  https://alexandra-kallitsopoulou.web.cern.ch/notes/unlisted/2026-11-19-where-do-the-missed-muons-go-the-vmm-cosmic-stand-step-by-step.html
  New: update box, 3b pointing resolution (+ before/after), 3c pillar ceiling (mask V1), 3d cuts table, 3e unphysical
  tracks, 5a distance histograms with fit and both radii, 5d why silent (illustrated), section 6 budget, open points.
  The PNG figures of section 8 (deck figures) are still from the fixed-15 mm funnel (MID changes < 0.5 points).
- The 20 Nov deck was NOT updated with the new radius / pillars / unphysical cut.

## 4. Open points (suggested order)
1. **P2_MID silent excess**: DREAM on MID (same pads; cosmic bench M3 hit/miss lists), and/or the mesh-HV scan (user's
   call); a MID <-> OUT position swap would separate chamber vs position in the stack.
2. **Bad pads repeating on the same channel numbers** across chambers (ch 62, 42, 37...): more statistics; check the
   channel's threshold / trim vs its silent share (calibration gap channels?).
3. **Pillar term at the VMM threshold**: SPS showed the pillar cost grows with threshold (2.4 -> 9.5 points); 3.78 % is
   a floor. Cosmics cannot resolve the 4 mm lattice (pads 12 mm, pointing 3-8 mm).
4. **Paddle positions**: refit on the weekend data (or measure) - then the paddle cut can be tightened (paddle-sys).
5. Refit `pointing.json` on all REF runs (now 14) - expect no change beyond 0.1 mm.
6. Update the deck figures (section 8 / 20 Nov deck) with the new cuts if they are shown again.
7. Older open points of HANDOFF_20261120 section 4 (gain-9 late hits, calibration gap, bscan FAR, crosstalk filter).
