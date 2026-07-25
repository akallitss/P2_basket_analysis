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
