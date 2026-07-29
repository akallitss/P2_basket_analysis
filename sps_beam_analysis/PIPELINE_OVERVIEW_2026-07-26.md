# P2 @ SPS H4 — analysis pipeline overview (TB_July2026_H4)

**Written 2026-07-26.** One document answering: what runs automatically, what the
numbered stages do, what each produces, where to look, what has been analysed so
far, and what is worth doing next. Companion to the dated logbook entries in
`logbook/`, which hold the per-run conclusions.

---

## 1. Data flow — the automatic layer

Everything below happens live during a run, driven from the DAQ repo
(`~/DAQ_Control_Dream_Beam`), without anyone touching it:

```
DREAM DAQ (FEUs) ──> raw FDF per FEU        runs/<run>/<sub_run>/raw_daq_data/
        │ processor_watcher (decode)   ──>  decoded_root/   (per-FEU nt trees, 16-sample waveforms)
        │ processor_watcher (combine)  ──>  combined_hits_root/  (analysed hits: time, amp, channel)
        │ qa_watcher                   ──>  analysis/<run>/<sub_run>/<det>/   hit & event maps
        │ hv_control monitor           ──>  runs/<run>/<sub_run>/hv_monitor.csv  (1 Hz V/I per channel)
        └ backup_watcher               ──>  EOS mirror (push-only)
```

- The **hits tree** (`combined_hits_root`) is enough for stages 20–28.
- The **decoded waveforms** (`decoded_root`) are needed only by stage 29 (and 30
  uses the raw decoded stream). Keep this in mind when pruning disk: on
  2026-07-25 `hits_root`/`raw_fdf` of drift_mesh_scan_1 were pruned but
  `decoded_root` survived — which is the only reason the timing-vs-drift study
  could still run locally on 07-26.
- `hv_monitor.csv` is the ground truth for what HV actually was (trips, sags),
  not the run config. The P2_MID trip during drift_mesh_2d_1 was diagnosed
  entirely from it.

## 2. The stage battery (`sps_beam_analysis/`, stages 20–30)

All stages share one convention:

```bash
export SPS_DATA_ROOT=/local/home/banco/P2_data/TB_July2026_H4/runs
export SPS_ANALYSIS_ROOT=/local/home/banco/P2_data/TB_July2026_H4/analysis
SPS_RUN=<run_name>  python3 <stage>.py live [--sub-run NAME]
```

No `--sub-run` = iterate every sub-run and (where meaningful) write scan-level
aggregations. `run_beam_qa.sh <run>` is the beam-time shortcut (25 + 24 on
every sub-run that has combined hits). Products land under
`analysis/<P2_IN|P2_MID|P2_OUT|telescope>/<run>/<sub_run|scan*>/<stage>/`,
browsable from the DAQ GUI (port 5001) → Analysis tab.

| Stage | What it does | Input | Key products |
|---|---|---|---|
| 20_beam_spectra | Cluster-charge spectra, Landau MPV vs HV/sub-run | hits | MPV-vs-HV curves (gain calibration of a scan) |
| 21_telescope_align | Mutual plane alignment from cluster-position correlations | hits | per-sub-run offsets; **gate for 22** |
| 22_tag_probe_efficiency | Tag-and-probe station efficiency, 3-plane majority tags | hits | per-pad `eff_map_*V.png/csv` per scan point (the GIF frames), scan-level curves |
| 23_beam_profile | Beam-spot + spill time structure | hits | profile maps |
| 24_event_sync_qa | FEU-to-FEU trigger alignment per station | hits | sync verdicts |
| 25_commissioning_qa | "Is the telescope alive": rates, occupancy, latency, HV verdicts | hits+HV | go/no-go tables (beam-time tool) |
| 26_hv_spark_qa | Spark counting/vetoing from hv_monitor + hits | hits+HV | spark rate vs HV; the **spark-vetoed** flag other stages consume |
| 27_pedestal_qa | Per-channel pedestal health | prg files | noisy/dead channel lists |
| 28_timing_qa | Hit-time distributions, time-walk, σ vs HV from the standard `time` branch | hits | `timing_vs_hv_<det>.png/csv` at scan level (light, fast) |
| 29_waveform_timing | TOA recomputed from raw waveforms (frac30/frac50/parabola/dCFD), fitted ftst clock-phase + walk corrections, pair Δt → single-station σ | **decoded waveforms** | per-point plots + `waveform_timing_summary.json`; ~15–19 min/sub-run |
| 30_raw_stream_efficiency | Efficiency per trigger from the raw ZS stream, bypassing clustering | decoded | cross-check of 22 |

**28 vs 29:** 28 is the quick QA (uses the one `time` number per hit — runs on
partial fetches, minutes per run). 29 is the real timing measurement (µs-level
corrections fitted from data). Trust 29 for physics numbers; use 28 for trends
while the beam is on.

## 3. Analysis status (as of 2026-07-26 evening)

| Run | 20 | 21 | 22 | 24 | 26 | 28 | 29 | Notes |
|---|---|---|---|---|---|---|---|---|
| beam_nominal_meshscan_1 (07-24) | – | ✓ | ✓ | – | ✓ | ✓+scan | 1 pt | most complete of the early runs |
| meshscan_fine_1 (07-24) | – | ✓ | ✓scan | – | ✓ | – | – | |
| drift_scan_1 / _2 (07-24) | – | ✓ | ✓scan | – | ✓ | – | – | |
| highstat_eff_1 (07-25) | – | ✓ | ✓ | ✓ | – | – | ✓ 2 pts | benchmark: σ = 15.5–22.9 ns/station at 700/450 |
| drift_mesh_scan_1 (07-25) | scan | ✓ all 23 | ✓ both halves | – | – | – | ✓ **both halves** | σ-vs-drift AND σ-vs-mesh curves complete 07-26 |
| low_mesh_scan_1 (07-26) | ✓* | ✓ | ✓ | ✓ | ✓ | ✓ | – | full pass done 07-26; *Landau MPV fit fails at every point — gain too low at gap 300 V, spectra themselves are saved |
| drift_mesh_2d_1 (07-26, 14/48, interrupted) | – | – | – | – | – | – | – | only qa_watcher maps; see caveats §4 |
| mesh_drift_scan_up_1 (07-26, complete 20:51, 37 pts) | **running** | **running** | **running** | **running** | **running** | **running** | – | hits-level pass launched 07-26 ~21:30; 29 not yet (37 pts ≈ 10 h — decide which points matter) |

## 4. Key results & important remarks

**Timing** (logbook 2026-07-25 entry, extended 07-26):
- Working point 700/450: **P2_MID 15.5, P2_OUT 18.4, P2_IN 22.4 ns** (triangulated,
  walk-corrected), stable <0.5 ns over 4.5 h.
- σ vs drift voltage (mesh 450): plateau at **750–800 V** (~16 ns MID); steep
  degradation below 600 V (drift-time spread dominates: ~150 ns at 500 V);
  slight worsening ≥850 V. **The working point is within ~0.5 ns of optimum.**
  Products: `analysis/telescope/drift_mesh_scan_1/scan/29_waveform_timing/`.
- σ vs mesh voltage (mesh half, gap constant): monotonic improvement, **no
  plateau reached by the ceiling** — MID 21.8→16.2 ns, OUT 24.3→19.2 ns over
  390→450 V; P2_IN flat ~32 ns below 385 V (amplitude-limited) then falling to
  27.2 ns at its 430 V point, still improving at its ceiling → argues for
  trying P2_IN at 445–450 V post-repair.
  Products: `.../scan/29_waveform_timing/timing_vs_mesh.{png,csv}`.
- Fitted ftst slopes differ per FEU (−9.6/−6.6/−5.5 ns/unit) — never assume a
  common clock phase; the per-station fit is load-bearing.

**Efficiency:** per-pad tag-probe maps exist for every mesh and drift point of
the 07-24/25 scans; animated versions (per detector + combined 3-panel) under
`analysis/*/drift_mesh_scan_1/scan_*/22_tag_probe_efficiency/*_anim.gif`.

**Data-quality caveats to carry into any future fit:**
- `drift_mesh_scan_1` drift_450 = zero drift field — starved stats, prompt
  conversions only; never fit it as a normal point. drift_700 has partial
  P2_OUT decoded data (403k vs ~2.3M hits) → MID–OUT pair absent there.
- `drift_mesh_2d_1`: dm_01_00/01/02 are beam-off (~91 MB each, misleadingly
  marked complete); P2_MID mesh tripped 4 min into dm_02_01 → its P2_MID data
  is dead for ~13/17 min. dm_02_00/01 exist **only on banco's disk** (EOS
  quota); treat as fragile until the quota is fixed.
- In stepped-together mesh scans, P2_IN runs at its own (lower) voltage — panel
  titles/filenames carry the true per-detector HV; frames advance by sub-run.
- The old P2_IN position map was wrong before the 07-25 evening fix; anything
  P2_IN-positional analysed before then was regenerated — do not resurrect old
  copies from EOS.
- EOS backups have been failing since 07-25 (salsachip quota). Until fixed,
  banco's disk is the only copy of the newest runs — mind `space_manager` /
  `prune_active_run.py` (it refuses to delete anything not verified on EOS,
  which currently means it deletes nothing new).

## 5. Potential improvements

1. **Fold the scan-level aggregation into stage 29** (28 already has one). The
   timing-vs-drift table/plot was assembled by an ad-hoc script on 07-26; a
   `--scan` mode writing `timing_vs_hv` like 28's would make it one command.
2. **Auto-run the deep pass at run end.** The watchers already know when a run
   finishes; a hook that fires 21→22→26→28 (hits-level, cheap) would have kept
   §3 green without anyone remembering. 29 stays manual (hours, needs babysitting).
3. **Triangulated per-station σ at every scan point** (drop the equal-resolution
   assumption) — needs all three pairs, currently blocked wherever a station is
   parked (drift scans) or partial; worth it on mesh scans where all three step.
4. **Backfill P2_OUT at drift_700** by fetching its raw FDF from EOS
   (`fetch_run_partial.sh`) and re-decoding — closes the one hole in the
   MID–OUT pair curve.
5. **Run 30_raw_stream_efficiency on one scan** and compare with 22 — the
   clustering-independent cross-check exists but has barely been used.
6. **Pedestal coverage:** pedestals in use cover FEUs 1/4/5 only (07-25 note) —
   retake full-telescope pedestals before the next physics run, and run 27 on
   them as a matter of routine.

## 6. Nice-to-have additional runs

- **Resume/finish `drift_mesh_2d_1`** (34 sub-runs, ~11.5 h) once P2_MID is
  cleared — the 2D efficiency+timing surface is the natural summary plot of the
  campaign. Resume procedure is staged (see DAQ repo,
  `config/json_run_configs/drift_mesh_2d_1_resume.json`); decide first whether
  to retake dm_01_00/01/02 (beam-off) and dm_02_01 (P2_MID dead).
- **P2_IN fine mesh scan 430→450 V** (its ceiling): its efficiency plateau
  starts ~440–450 and its timing (22–23 ns) is amplitude-limited — establish
  whether 445–450 buys timing before the next long run. Now possible at full
  quality with the repaired 8:1 channel.
- **Post-repair high-stat benchmark**: highstat_eff_1's P2_IN predates the
  connector repair and map fix era — a 1–2 h run at the working point with all
  three healthy stations would replace it as the reference.
- **A short 750–800 V drift confirmation run** only if the working point is
  ever revisited — the 07-26 curve says 700 V is already near-optimal, so this
  is low priority.

---
*Batch launched 2026-07-26 evening: low_mesh_scan_1 full pass, then stage 29 on
the drift_mesh_scan_1 mesh half. Update §3 when it lands; mesh_drift_scan_up_1
gets the same pass after it ends (~21:00).*
