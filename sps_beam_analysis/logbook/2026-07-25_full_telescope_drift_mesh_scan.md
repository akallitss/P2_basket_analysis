# 2026-07-25 — Full-telescope drift + mesh scan, then the high-stat efficiency run (TB_July2026_H4)

Runs `drift_mesh_scan_1` (overnight, 23 sub-runs), `env_test_1` (10:54, config
shakedown) and `highstat_eff_1` (11:37 → 14:10, 5 × 30 min at the settled
working point, ~42 M triggers — the primary efficiency/alignment dataset so
far).

## Run 1 — `drift_mesh_scan_1` (01:03 → 05:07, 23 sub-runs × 10 min)

Finished normally. **The whole voltage program in a single run, with all five detectors
back in: P2_MID's gain collapse is fixed — median hit amplitude 287 ADC at
nominal, back on par with P2_OUT's 296, vs 41 at the same settings yesterday
midday — and P2_IN is in the telescope at its new working point (mesh 430 /
drift 630), recording hits in 90 % of triggers with a quiet mesh (imon
≤ 0.7 µA, vs 6 µA spikes at the old 490 V point).**

Fresh pedestals taken at 00:38 before the run. Four FEUs (1/3/4/5), processed
on the fly with the **fixed** analyze_waveforms, ~2.8 M triggers per sub-run
(~4.6 kHz). Program:

- `drift_450` … `drift_900` (10 pts): MID/OUT drift scan at mesh 450 — the
  repeat of yesterday's drift scan, now with the recovered MID and P2_IN in.
- `nominal_00`: MID/OUT 700/450, P2_IN 630/430.
- `meshscan_01..12`: MID/OUT mesh 445 → 390 in 5 V steps **and P2_IN mesh
  425 → 370 in lockstep** (drift gaps held at 250 V / 200 V) — the fine-scan
  repeat, plus P2_IN's first mesh curve.

## Quick numbers (direct hits-tree pass, no QA yet — see operational note)

"share" = events with ≥1 hit / triggers; amp = median hit amplitude [ADC].

| sub-run | P2_IN (mesh) | P2_MID | P2_OUT |
|---|---|---|---|
| nominal_00 | **0.90 / 97** (430) | **0.97 / 287** | 0.96 / 296 |
| drift_450 (zero drift field) | 0.90 / 95 (430) | 0.22 / 50 | 0.20 / 58 |
| meshscan_01 (midout 445) | 0.89 / 88 (425) | — | — |
| meshscan_06 (midout 420) | 0.73 / 54 (400) | — | — |
| meshscan_12 (midout 390) | 0.35 / 39 (370) | 0.84 / 71 | 0.89 / 92 |

- **P2_MID is recovered** — ×7 the collapsed amplitude, statistically on par
  with P2_OUT again. The gas line / mesh contact was the suspect; whatever the
  intervention was, record it here. All of yesterday's MID working-point
  conclusions (its fine-scan curve, its drift curve) are superseded by this
  run's data taken in the healthy state.
- **P2_IN's first mesh curve**: share 0.90 / 0.89 / 0.73 / 0.35 at mesh
  430 / 425 / 400 / 370 — still rising at 430, so its plateau is above 430,
  but 430 is already telescope-worthy and stable. The upward scan from 430
  toward the imon-spike onset (the real maximum) is still to be done.
- The drift-scan leg gives the transparency turn-on for MID/OUT in the healthy
  state; drift_450 (zero gap field) suppresses MID/OUT to ~0.2 share while
  P2_IN (fixed HV) stays at 0.90 — a clean control that the effect is the
  drift field, not the beam or DAQ.

## Data-quality flags

- `drift_700` FEU 5 (P2_OUT): **the decode-hang gotcha struck again** — chunk
  000 (971 MB) sits renamed as `.fdf.hang` in raw_daq_data and was skipped;
  only the last ~14 % of the sub-run (chunk 001) is processed. The raw file is
  intact, so the point is recoverable by re-decoding it by hand (last time the
  re-decode took 2.5 s). Meanwhile `nominal_00` is at the identical HV and
  covers the point.
- `recorded_events.npz` is not yet extracted for this run — run
  `extract_recorded_events.py` before any efficiency fits so the DAQ-overlap
  correction is available.

## Run 2 — `highstat_eff_1` (11:37 → 14:10, 5 × 30 min at the working point)

The high-statistics efficiency / alignment dataset in the healthy state: **all
three P2 stations at drift 700 / mesh 450** (uRWELLs 600 / 420 as always),
30-minute sub-runs `beam_commissioning_00..04`. A 6th sub-run was configured
but the run was stopped from the GUI at 14:10 after 04 completed — 5 sub-runs,
**~8.4 M triggers each at ~4.6 kHz → ~42 M triggers total**, ~13.5 GB raw per
sub-run (68 GB), zero decode hangs, all processed on the fly with the fixed
analyzer.

Quick hits-tree pass, first vs last sub-run (share = events with ≥1 hit /
triggers; amp = median hit amplitude):

| detector | share (00 → 04) | amp (00 → 04) |
|---|---|---|
| P2_IN (mesh **450**) | 0.94 → 0.94 | 180 → 177 |
| P2_MID | 0.98 → 0.97 | 305 → 295 |
| P2_OUT | 0.97 → 0.96 | 357 → 302 |
| uRWELLs (FEU 1) | 0.74 → 0.74 | 58 → 58 |

- **P2_IN held mesh 450 for the whole 2.5 h**: share 0.94, amp ~180 ADC, mesh
  current quiet (imon max 0.8–1.3 µA, no spikes). So 450 is still on the good
  side of its instability — the collapse edge sits somewhere in 450–490, and
  the working point no longer costs a factor in gain vs its neighbours.
- **Stability**: shares constant to ~1 pt; P2_IN and P2_MID amplitudes drift
  −2/−3 %, but **P2_OUT slid −15 % (357 → 302) over the run** — the same
  common-gain (P/T-type) drift seen overnight. For the efficiency analysis,
  calibrate amplitude per sub-run or use amplitude-free efficiency.
- **First dataset with the corrected uRWELL mapping**: before the run the
  detector config was fixed (backups kept next to run_config.json) —
  `det_type` → `urw_inter` (front) / `urw_strip` (back) and all uRWELL FEU
  orientations set to `inverted`. `env_test_1` (10:54, one sub-run) was the
  shakedown of that change.
- Caveat carried from the offline session: **P2_IN position reconstruction is
  still broken** (fires on 94 % of triggers here, but tag-probe ~2 % with a
  93 mm residual — channel→pad mapping/orientation wrong). This run is the
  dataset to validate the corrected P2_IN map against once it exists; until
  then P2_IN counts, don't position-analyse it.
- Sub-run 04 is missing its `.subrun_complete` marker (stop-button race, data
  and hits are complete on disk) — touch it or expect the fetch scripts to
  skip 04.

## Operational note

**qa_watcher is down** — last activity 18:50 on 07-24 (launching QA for
`drift_scan_2/drift_850`), not in the process list this morning, nothing in its
log about why. processor / backup / pedestal watchers and the Flask GUI are all
alive. Consequence: no waveform/telescope QA exists yet for `meshscan_fine_1`
(from `meshscan_06` on), `p2in_check_1`, `drift_scan_2` (850/900) or the
overnight run — the numbers above come from a direct pass over the hits trees.
Left as-is in case it was stopped deliberately; restart it (or run
`run_beam_qa.sh drift_mesh_scan_1`) when banco has headroom.

Also of note this morning: banco refused new SSH key-auth connections for ~40
minutes (~09:35–10:15, the known MaxStartups-under-load behaviour) while the
Flask GUI stayed reachable — GUI port 5001 answers `/get_current_run` and is a
usable fallback for checking run state when SSH is locked out.

## Analysis

Same chain and paths as the 07-24 entry. Products will land under
`TB_July2026_H4/analysis/telescope/drift_mesh_scan_1/` once QA runs; raw +
hits under `TB_July2026_H4/runs/drift_mesh_scan_1/`.

## P2_IN position reconstruction: FIXED (evening session)

The "broken pad mapping" (fires on 90-94 % of triggers but tag-probe ~2 % with
a ~90 mm residual that even a fitted rotation could not reduce) is a **single
flipped ribbon: P2_IN's `c_5_top` (dream conn 4) is cabled 'linear' while every
other connector half on all three stations is 'reverse'.**

Diagnosis was data-driven on `highstat_eff_1/beam_commissioning_00`: for every
P2_IN electronic channel, the median P2_MID single-cluster position in events
where that channel fires (P2_OUT run as the method control). In-beam channels
(n >= 5000 tags) give, per connector half, a 0.0 mm median residual for exactly
one strip-order hypothesis:

| half (dream conn) | reverse | linear |
|---|---|---|
| c_4_top (2), c_5_bot (3), c_6_bot (5), c_6_top (6) | **0.0** | 81-109 |
| **c_5_top (4)** | 138.0 | **0.0** |

(P2_OUT control: reverse wins at 0.0 mm on all seven illuminated halves.)

Fix: `p2_mapping.build_channel_table` gained `strategy_overrides=
{(connector_N, half): strategy}`; `sps_config.STRATEGY_OVERRIDES =
{'P2_IN': {(5, 'top'): 'linear'}}` is applied by default for every run of this
beam test. Validated: with the override the same in-beam channels give 0.0 mm
median AND max residual on all six illuminated halves.

Notes:
- The c_8 (vs MID/OUT's c_7) declaration on P2_IN dream conns 7/8 could NOT be
  tested conclusively: those channels carry <3 % of hits, are outside the beam
  spot, and their coincidence with P2_MID is diffuse (no localized partner —
  halo/shower hits). The position test weakly favours c_7 (same ordering as
  the P2_OUT control). Left as declared; needs a run with the beam (or a
  source) on that sector, or a cabling photo.
- Everything position-based involving P2_IN before this fix is invalid: its
  drift_mesh_scan_1 efficiency "curves" (~1.7 % flat), its 21 transforms, its
  20 per-pad maps. MID/OUT products are untouched (their maps were correct).
- Re-runs with the fixed map: 21+22 on highstat_eff_1 (in progress), then
  21+22 on drift_mesh_scan_1 and a P2_IN-specific 20 pass with `--fit-min`
  lowered (~60; its cluster MPV at mesh 430 sits below the default 200 ADC
  fit window, which is why all its drift-arm Landau fits failed).

Same session, analyzer fixes to the scan-axis handling: 22 now decides the
scan axis PER PROBE (fixed-HV control planes get the program's scanned value
as x and sub_run-named products — before, the drift-arm plot collapsed onto
mesh HV and fixed-HV probes overwrote one eff map per point); 20 gained
`--subruns-glob` (drift/mesh arms of a mixed run no longer overwrite each
other's `scan/` products) and the same per-det axis logic; its exponential
gain fit is now mesh-axis-only.

MID/OUT drift-arm gain curves (20, fixed axis): Landau MPV peaks at drift
650-700 (MID ~206, OUT ~320 ADC) and falls ~20 % by 900 (mesh transparency
rollover) while tag-probe efficiency still rises to 900 — so drift 700 is the
amplitude optimum, and the efficiency plateau above it is flat in exchange for
gain. Spark-sag contamination identified on four efficiency points (post-
discharge gain suppression, minutes-long, survives the imon veto): P2_MID @
mesh 445 (meshscan_01), P2_OUT @ drift_900, nominal_00 and mesh 440 — mesh
vmon sag episodes on 8:3 / 8:5 at exactly those sub_runs. Exclude from curve
fits.

Also: `highstat_eff_1` has a SIXTH sub-run — `beam_commissioning_05`,
16:11-16:42, full ~31 min at the identical working point (all P2 700/450,
uRWELLs 600/420), taken after the morning entry above was written. It lacks
`.subrun_complete` (touch it for the fetch scripts); analysis stages pick it
up automatically. ~50 M triggers total at the working point.

## highstat_eff_1 tag-probe efficiencies (fixed P2_IN map, 3-plane majority tags)

Per sub-run, DAQ-overlap-corrected, probe radius 27 mm (3 sigma of the 9 mm
alignment residual), spark-vetoed, ~5.5-7.8 M tags per point:

| sub-run | P2_IN | P2_MID | P2_OUT | note |
|---|---|---|---|---|
| 00 | 0.925 | 0.963 | 0.925 | |
| 01 | 0.922 | 0.957 | **0.843** | OUT spark episode, see below |
| 02 | 0.915 | 0.954 | 0.936 | OUT FEU 13.1 % DAQ loss, corrected |
| 03 | 0.913 | 0.944 | 0.907 | IN+MID FEUs ~10.7 % DAQ loss, corrected |
| 04 | 0.921 | 0.951 | 0.933 | |
| 05 | 0.917 | 0.962 | 0.920 | evening sub-run |

- **P2_IN with the corrected map: 0.913-0.925** (was ~0.017 with the broken
  one), alignment residual 9.0 mm (was 93), dy +3.2 mm, theta -0.22 deg. The
  telescope now majority-tags with all three planes.
- **P2_OUT sub-run 01 = spark-suppression case study at the working point**:
  28 mesh imon samples > 2 uA (vs 7 / 2 in neighbours), 24 of them in min
  0-11.5; the 2-min median amplitude sits 10-25 % low (360-456 ADC vs the
  490-520 plateau) exactly over that window, then recovers. The sub-run
  average 0.843 hides a suppressed first third + healthy remainder. Same
  post-discharge mechanism as the four flagged drift_mesh_scan_1 points.
- Method note for precision numbers: a gain-stability gate (drop time windows
  whose ~2-min median amplitude is < ~90 % of the sub-run plateau) would
  remove these transients cleanly; the imon spark veto alone cannot.

## drift_mesh_scan_1 rerun with the fixed P2_IN map (products regenerated)

- **P2_IN's first real mesh efficiency curve** (tag-probe, DAQ-corr., mesh
  370→430 + nominal): 0.32 / 0.45 / 0.59 / 0.70 / 0.79 / 0.85 / 0.88 at
  370/380/390/400/410/420/430 V — still rising at 430; with highstat's 0.92 at
  450, the plateau onset is around 440-450, right where MID/OUT sit.
- **P2_IN drift-arm control** (own HV fixed 630/430): flat 0.86-0.89 across
  drift 500-900, MPV flat at 63-65 ADC — clean confirmation the drift-scan
  response of MID/OUT is a real field effect.
- **Gain curves (20, fit-min fixed)**: P2_IN MPV 23→63 ADC over mesh 400→430;
  MID 124→238 over 430→450; OUT 179→333. MID's curve is cleanly exponential
  INCLUDING the spark-flagged 445 V and nominal points (x2 every ~21 V) —
  earlier "MPV 155 at 445" was a fit-window artifact (MPV below --fit-min
  200). OUT's nominal_00 MPV is ~3 % below its exponential. So the spark
  contamination shows up strongly in efficiency (2-4 pt dips) but only mildly
  in whole-sub-run MPV — a short suppressed window drags the efficiency
  integral much more than the fit's peak position. The gain-stability time
  gate remains the right tool; whole-sub-run MPVs are usable as-is.
- MID/OUT mesh efficiency curves regenerated with 3-plane tagging (P2_IN now
  tag-eligible): shapes unchanged, values shift <~1 pt.
