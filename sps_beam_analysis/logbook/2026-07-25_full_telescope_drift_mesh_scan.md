# 2026-07-25 — Overnight: full-telescope drift + fine mesh scan (TB_July2026_H4)

Run `drift_mesh_scan_1` (01:03 → 05:07, 23 sub-runs × 10 min, finished
normally). **The whole voltage program in a single run, with all five detectors
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
