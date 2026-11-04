> **Resolved 2026-11-04 (night): see `FINDINGS_paddle_efficiency.md`.** The cause is H2: the paddles' efficiency varies along their length (light attenuation against the threshold; PMT of 44 at the lower-left end, PMTs of 41 / 45 at the upper-right end). Paddle 44 (in every trigger) is 0.90 efficient on its lower-left half and 0.23 on its upper-right half; 41 is the mirror image; 45 is fine. Measured with `scint_dump.py` / `paddle_efficiency.py` / `paddle_study.py`; the MC with the measured maps reproduces all three tag distributions with no free parameter. H1 / H3 / H4 are ruled out; the paddle-position fit below is an artefact. The rest of this file is kept as the record of the question.

# Handoff: why does the det1 tag distribution disagree with the geometric expectation?

Written 2026-11-04 (evening). Branch `worktree-trigger-acceptance-mc` of P2_basket_analysis. Not merged, not pushed
(GitHub key refused). The user wants to understand **exactly** what is happening and believes the data in hand may be
enough to solve it.

## 1. The question in one paragraph

Lab cosmic stand, scintillator-referenced wk2 scan (runs 181-208). For each probe chamber P, a *tag* = the paddle
trigger plus a clean cluster in each of the two other chambers; the predicted point in P is the straight line through
the two clusters. The tag density (the kernel-radius row of `analysis_wk2/eff_definitions/maps_K*.png`) therefore
**cannot depend on P's own efficiency**. A toy MC with every chamber 100 % efficient reproduces the measured tag
density for det3 and det2. **It does not reproduce det1.** Geometry predicts det1's dense region upper-right
(x ≈ 450 mm), like the others. The data put it lower-left (x ≈ 250 mm), and the upper-right end of det1's read-out band
is nearly empty. Find the cause.

## 2. What is established (numbers you can rely on)

- **Stand (user, 2026-11-04):**
  - paddles 200 x 600 x 50 mm;
  - two top paddles exactly stacked (z -100..0), trigger = (top 0:41 OR top 0:45) AND bottom 0:44 (offline OR in `scintor_*.json`);
  - chambers det3 P2_IN / det2 P2_MID / det1 P2_OUT at z = 320 / 630 / 940 mm. These are the analysis-code values (SPS spacing), NOT measured on the stand; the user says "about 30 cm apart";
  - bottom paddle about 30 cm below det1, z ≈ 1250. The logbook's z = 970 (run_126) is outdated;
  - paddle x-y NOT recorded; the user says "roughly aligned with the read-out connectors".
- **Read-out area:** 384 pads per chamber = Gerber `connector_number` 3-5 (= `channel_id//128+1` 4-6). Band centre (381, 148) mm, long axis 21.5 deg to x. **Each VMM reads a strip running the whole length of the band**; VMM 2..7 (det3), 8..13 (det2), 14..19 (det1) sit side by side ACROSS the band. So no single chip can create an along-band gradient.
- **Data (all 10 runs pooled), along-band position of the tags relative to the band centre:**

  | probe | data median | data mean | MC median, paddles centred | MC, bottom paddle -100 mm along band |
  |---|---|---|---|---|
  | det3 | +6.6 | -4.1 | +3.7 | +9.4 |
  | det2 | -26.1 | -31.0 | +15.5 | -1.1 |
  | det1 | -83.6 | -68.5 | +1.6 | -42.4 |

  Across the band, all three match the MC well; the band width dominates there.
- **MC fit** (binned Poisson shape likelihood on the read-out footprint, 3 probes jointly; only the paddle x-y/angle free):
  - bottom paddle about (-73, -9) mm from the top pair, about 71 mm along the band;
  - chi2-like per probe: det3 1140/171, det2 1066/162, det1 2851/169 bins. The 1-D along-band fits per probe disagree: det2 wants -150..-200 mm, det1 about -50, det3 about -250 (weak). Paddle position alone cannot fit all three.
- **Weighting by tag-chamber efficiency changes almost nothing.** Model "geometry x tag-chamber efficiency": each MC tag weighted by the measured local efficiency (hit within 30 mm, 200 nearest tags) of the two tag chambers.
- **Absolute numbers:**
  - MC trigger rate 0.90 Hz at 1 /cm2/min vs 0.57 Hz measured (after veto);
  - only 0.44 of triggers cross all 3 read-out areas even at eff = 1 (data: 0.14 clean in all three);
  - implied clean-cluster efficiency per trigger crossing the area: det3 0.67, det2 0.53, det1 0.65.

## 3. Hypotheses, and what each predicts

Data minus MC (fitted paddles), along the band: roughly **det3 ≈ 0, det2 ≈ -30, det1 ≈ -55 mm**. The disagreement
**grows with depth and keeps one sign**. Use that to discriminate.

**H1. Chambers physically offset sideways on the stand** (a "shear": det2 displaced by s, det1 by about 2s along the
band, or any lateral offsets/rotations).
- This changes the *acceptance* (which tracks hit all three read-out bands), not the fit quality.
- A pure linear shear is **invisible to the straight-track residual alignment** done earlier. That alignment found det2 within 9 mm of the det3-det1 midpoint and rotations < 0.4 deg; it constrains only det2 - (det3 + det1)/2.
- `explore.py` (below) with det1 frame +80 mm along the band (det2 half) gave a shift of +35 / -16 / -72 relative to aligned. Right direction for det2/det1, wrong for det3, but combined with a paddle offset it may fit.
- **Most plausible mechanical explanation; test first.**

**H2. Paddle efficiency non-uniform along the 60 cm length** (light attenuation; PMT at one end; discriminator threshold).
- A bottom-paddle gradient imprints most strongly on the chamber nearest to it (det1), then det2, then det3. That is the observed depth ordering. A top-paddle gradient would do the opposite.
- Changes rates only, not overlap geometry.
- **Directly measurable** from data: the H0 VMM records the paddle pulses (ch 41, 45 top; ch 44 bottom; ch 40 hardware triple). Plot ch 44 ADC (and timing) against the track position extrapolated to z ≈ 1250 along the paddle. Attenuation along the paddle = H2.

**H3. A position error in one chamber's measured coordinates** (mapping/order errors, e.g. the P2_IN order overrides on c3/c5 top).
- A single-chamber error predicts shifts of *opposite* sign in the other two probes, as the linear algebra below shows. The data have the same sign in det2 and det1, so a single-chamber error is disfavoured as the main cause. (The note on the website currently names this as the likeliest cause; that claim should be revised.)
- Shifts of the predictions for an error δ (det3), ε (det2), η (det1):
  - det3 probe (2·x2 - x1): 2ε - η
  - det2 probe ((x3 + x1)/2): δ/2 + η/2
  - det1 probe (2·x2 - x3): 2ε - δ
- Solving for data 0 / -30 / -55 gives δ ≈ -2, ε ≈ -29, η ≈ -58: det3 unshifted, det2 -29, det1 -58. That is a linear shear, i.e. H1 in disguise if the frames' error is a physical offset. Local (non-uniform) mapping errors remain possible on top.

**H4. Wrong z positions.**
- Changes the extrapolation lever arm; mostly stretches rather than shifts. Only shifts through the mean slope of the selected tracks.
- Cheap to rule in or out by measuring the stand.

**H5. Clean-cluster efficiency varies along the band in the tag chambers** (not chip-level: e.g. radial gain gradient, noise hits joining the window).
- The cluster in `scint_efficiency.py` is the **unweighted mean of ALL live hits in ±150 ns**; it is "clean" if ≤ 6 hits all within 15 mm. Noise hits make it unclean or pull it.
- Would show as a tag deficit in all probes using that chamber as a tag.

## 4. What the existing products can and cannot do

- `analysis_wk2/eff_definitions/json/scintor_run*.json`: per probe, per tag, the **predicted** point (x, y, rounded 0.1 mm), `d` (distance to the nearest probe hit; 9999 = none), `any`, `live`, `w50`; plus probe pads (x, y, vmm, ch, live).
  - **No event index and no tag-cluster positions**, so the probes cannot be linked event by event.
  - The cheap test possible now: `d` and hit fraction against the along-band position for det1 (a position error would show as growing `d`).
- To decide H1/H2/H3, **re-run `scint_efficiency.py` with a per-trigger dump**:
  - trigger time;
  - per chamber: has / clean / cx / cy / number of hits / hit list (pad id + ADC);
  - the paddle channels' ADC and time (H0 ch 40, 41, 44, 45).

  With that you can:
  - fit 3-chamber tracks and get the slope distribution of the tagged sample. Compare to the MC: H1 vs H2 differ in the slope vs position correlation;
  - plot residuals against position per chamber (H3);
  - plot bottom-paddle ADC against extrapolated position (H2);
  - count tags by tag-cluster position directly (H5).
- MC side: `trigger_acceptance_mc.py` already supports chamber offsets in `Chamber(..., offset)`, but `Pool.pred` uses TRUE positions. The frame transform (measured = stand position - chamber offset) is done in `explore.py`. Add chamber dx/dy (and rotation) as fit parameters, plus a linear efficiency gradient along each paddle, and refit.

## 5. Where everything is

- **Code (this branch):**
  - `cosmic_vmm_analysis/trigger_acceptance_mc.py`: MC + paddle fit + figures + `geometry_3d.json`;
  - `cosmic_vmm_analysis/note/trigger_mc_note.py`: the CERN-site note generator.
- **Scratch scripts** copied to `cosmic_vmm_analysis/handoff_scratch/` (exploratory, paths hard-coded; run from cosmic_vmm_analysis with `python handoff_scratch/explore.py . `):
  - `explore.py` (median shifts for paddle / shear variants), `explore2.py` (2-D likelihood grid), `explore3.py` (1-D along-band likelihood per probe), `vmm_along.py` (VMM extent along / across the band), `run_mc.sh` (env + driver).
- **Environment:** `export VMM_QA_DIR=~/Documents/PostDocSaclay/DAQ_Control_VMM_Beam/vmm_qa`; python = `~/miniconda3/bin/python` (system python3 has no pandas).
  - Full MC: 4e7 muons, about 3 min, peak RSS 1.3 GB (2e6-muon chunks; the machine has about 5 GB free). Run as `run_mc.sh OUT 4e7`.
- **MC outputs:** `~/Documents/PostDocSaclay/data/LAB_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_wk2/trigger_acceptance_mc/`
  - `fit.json`, `geometry_3d.json`, `kernel_K+0.png` (data vs MC kernel radius), `projections_all.png` / `projections_K+0.png` (along / across band), `density_geometry*.png` (2-D data / MC).
- **Data jsons:** `.../analysis_wk2/eff_definitions/json/scintor_run{181..208}.json`.
- **Raw data:** `/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/runs/run_N/cfg_*/raw_daq_data/*.pcapng`.
- **lxplus:** code `~/cva3`, jsons `~/wk2/def`, driver `~/run_def.sh`; per-run spark veto `~/wk2/veto_run_N.json`.
  - TRAP: on lxplus987 `mkdir -p` / `os.makedirs` on the AFS home HANGS; write to /tmp, use plain mkdir.
  - ssh lxplus needs a 2FA login in a real terminal first (ControlMaster).
- **Map options used for every wk2 json:** `--conn-shift -1 --order "P2_IN:3=bot:rev|top:lin" --order "P2_IN:5=bot:rev|top:lin"`, triggers `--top 0:41,0:45 --bottom 0:44 --coinc-win 60`, `--bad-times veto_run_N.json`, `--half-win 150 --max-tag 6 --cluster-r 15`.
- **Note (live, unlisted):** https://alexandra-kallitsopoulou.web.cern.ch/notes/unlisted/2026-11-04-cosmic-stand-3d-trigger-geometry.html
  - Its section 5 names H3 (det3 positions) as the likeliest; update it once this is resolved.
  - Regenerate with `trigger_mc_note.py MC_DIR cern-site/notes/unlisted/<file>.html`, then `build_hub.py`, commit, `scripts/deploy.sh`.

## 6. Suggested order

1. **Ask the user to measure the stand:** paddle centres and orientation in the P2 pad frame, chamber heights, and the lateral position of each chamber (e.g. one corner of each frame from a plumb line). This is probably the fastest route to H1/H2/H4.
2. Cheap test on the existing jsons: det1 `d` (and `w50`) against the along-band position.
3. Per-trigger dump on lxplus (one or two K = 0 runs: 187, 196) and the three checks of section 4.
4. Extend the MC: chamber offsets/rotations + paddle efficiency gradient; refit all three probes. A good fit should bring det1 chi2 near det3/det2.
5. Consequence to report: if H1 holds, every chamber-only alignment or efficiency number that assumed a common footprint (read-out-area cut, margin, live-area) uses slightly wrong areas. The det1 efficiency in its upper-right is least trustworthy.

## 7. Side finding (not this question)

`cosmic_acceptance_mc.py` (2026-10-28) selects read-out pads with `channel_id//128+1 in --connectors`. The read-out
pads are channel_id//128+1 = 4..6 (Gerber connector_number 3..5), so its default `3,4,5` is probably one connector
off. Its expected-rate numbers in the Oct note should be rechecked.
