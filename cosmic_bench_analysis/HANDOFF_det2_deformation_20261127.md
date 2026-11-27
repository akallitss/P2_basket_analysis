# HANDOFF — det2 deformation study, 27 Nov 2026 (night)

**Read `AUDIT_det2_deformation_20261127.md` first** (findings, evidence, numbers, open points) and
`P2_basket_analysis/CLAUDE.md` section 7 (rules). This file only says where things are and how to rerun.

## State

| what | where | state |
|---|---|---|
| code | `P2_basket_analysis` `main`: stages 43-49, `p2_deform.py`, `p2_deform_model.py`, `gas_studies/water_scan/` | committed |
| det2 products | `…/Analysis/det2/p2_det2_long_run_430_repositioned_11-25-26/long_run_det2_430_730/{43_deformation, 45_deform_local_z, 46_deform_fe55_vs_mip, 47_deform_pulse_maps, 48_deform_microtpc, 49_deform_decomposition}` | final re-run with the committed code |
| det2 scan products | `…/Analysis/det2/p2_det2_alignment_mesh_drift_scan_long_run_430_11-24-26/44_deform_scan_regions/` | done |
| det1 / det5 cross-check | same folder names under `det1/…/long_run_det1_415_615/` and `det5/…/long_run_det5_420_820/` (43, 47, 48, 49) | first pass only (audit section 9) |
| CERN note | not touched | the deformation results are not in the note |

## Rerun (det2)

```
cd cosmic_bench_analysis ; PY=../.venv/bin/python
$PY 43_deformation_waveforms.py det2_long_repos        # 15-30 min, reads decoded_root once; needs 34
$PY 47_deform_pulse_maps.py det2_long_repos --t-ref 125
$PY 48_deform_microtpc.py det2_long_repos
$PY 45_deform_local_z.py                               # reads 47's pad_gap_prelim.csv
$PY 44_deform_scan_regions.py                          # 24 Nov scans, hit level
$PY 49_deform_decomposition.py det2_long_repos --scan-dir <…/44_deform_scan_regions>
```

For another chamber: 34 (alignment) must exist for the run key; `47` takes `--ref-conn` (connectors used as
the reference zone) and `--amp-min`; `49` takes `--s-radius 20` when the gain varies on a short scale
(det1). Without `--scan-dir` stage 49 uses the dry-gas table and does not correct the collection vs field.
`43_deformation/wf.npy` + `wf_table.pkl` are 100-150 MB per run and are NOT in git (Analysis tree).

## What to look at first

`49_deform_decomposition/summary.txt`, `deformation_maps.png` (gap, drift time, spread, amplitude,
intrinsic gain, amplification gap, ion tail, delay), `consistency.png` (two days / methods, ion tail vs
gain, which surface), `effects.png` (lost muons, pillar shadow, amplitude, time spread vs gap),
`ridge_profile.png`; per-pad table `pads.csv`.

## Open, in suggested order

1. Bench check of the ridge on both faces (audit section 10.1), and a 1 h run at another gas overpressure.
2. det5: drift time 1.5 x det2's from the pad pairs. Check the drawings / gas of September before using any
   det5 drift-field or timing number; then redo 47 with `--amp-min 120` (low gain).
3. Put the deformation into the efficiency story: stage 27's loss table and the 27 Nov pillar audit quote
   chamber averages over an uneven gap; the shadow model of stage 36 should take the local gap.
4. Low-field drift velocity (audit section 8).
5. July Fe55 open points (gain contrast, neighbour-pad signal), amplification-gap number from Garfield++.
