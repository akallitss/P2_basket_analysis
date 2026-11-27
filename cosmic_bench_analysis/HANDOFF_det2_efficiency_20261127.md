# HANDOFF — det2 (P2_2) cosmic-bench efficiency investigation, 25-27 Nov 2026

Start here for the next session on det2 / the DREAM cosmic bench. Everything below is
committed; the code is on `main` after the 27 Nov end-of-day sweep (branch was
`det2-p1-scans-11-24-26`).

Public write-up (pedagogical, all figures): `cern-site/notes/2026-11-26-chasing-det2s-lost-muons-a-cosmic-bench-efficiency-investigation.html`
-> https://alexandra-kallitsopoulou.web.cern.ch/notes/2026-11-26-chasing-det2s-lost-muons-a-cosmic-bench-efficiency-investigation.html
(**committed + pushed to GitHub but NOT deployed yet**: run `scripts/deploy.sh` in cern-site from a
real terminal, lxplus needs the password).

---

## 1. Bench / hardware state (rays machine, 27 Nov morning)

| Item | State |
|---|---|
| det2 HV | **all channels OFF** (switched off by hand 26 Nov 11:20) |
| det2 mesh | **resistive short ≈10 MΩ**: tripped 26 Nov 01:07:36 (10 µA limit, sagged to 110 V); on 26 Nov 11:14 it drew current ∝ V (0.27 µA at 2.75 V, 2 µA at 21 V) and tripped at 110 V. Needs bench inspection before any new HV. |
| DAQ tmux | daq_control / dream_daq / hv_control / flask_server restarted 26 Nov 11:22 (idle) |
| processor tmux | running, idle; its `include_runs` points at `p2_det2_efficiency_batch2_430_11-26-26` which was **never taken**. Kill it or repoint before the next run. |
| run_config.py (rays + repo) | points at the 3 h batch-2 run (`run_time 180`, 430/730). Do not start it until the mesh holds 430 V with ~0 µA. |
| Disk /mnt/cosmic_data | 155 GB free (83 %); all dream_run scratch removed after byte-size checks (manifest `P2/DELETED_manifest_dream_run_20261126_1130.txt`); mesh-off segments 007-010 of the repositioned run deleted |
| Data fetched | `p2_det2_long_run_430_repositioned_11-25-26` on the LaCie (light tier + decoded_root, 6.7 GB) |

Runs of the series (all det2 on P1, FEU 6/7, M3 trigger FEU 1):
- `p2_det2_alignment_mesh_drift_scan_long_run_430_11-24-26` — alignment, mesh scan 340-420, drift scan 430-955 (DAQ hung from drift 980 V)
- `p2_det2_long_run_430_11-25-26` — 3 h, **det2 had moved out of the M3 acceptance**: unusable, not analysed
- `p2_det2_long_run_430_repositioned_11-25-26` — the long run; good until the mesh trip (`t_max_h = 7.58`, trip at 7.614 h)

## 2. Hardware defects found (for the bench, all confirmed on 11-24 AND 11-25 data)

| Defect | Channels | Evidence | Bench action |
|---|---|---|---|
| Connector 9 top half read **reversed** | FEU 7 ch 448-511 | muon over ch c fires ~959-c | label/flip; analysis fixed by `strategy_overrides={(9,'top'):'linear'}` |
| **Dead pin row** | FEU 6 every even ch 384-446 (connector 5 bot) | ~10 % firing, no muon response | inspect connector 5 (whole row open) |
| **Dead-but-noisy pads** (floating inputs) | FEU 6 odd 217-227, even 246-254 (conn 3 top); FEU 7 254, 281 | fire 1.4-4.7x median on noise, never on their muon | re-solder connector 3 |
| **Shorted groups** (identical waveforms, ratio 1.00±0.04) | FEU6 48-51, 209/211, 236/238, 241/243; FEU7 271-273, 338-340, 351/352, 353-355 | each group contains footprint neighbours 0.5 mm apart in a pin row; on conn 3 next to the open pins | re-solder board connectors 2, 3, 8 (solder bridges) |
| Weak pads (kept in eff) | FEU 6 452, 456; FEU 7 261; candidate FEU 7 29 | ~27 % response on own pad | genuine low gain |

## 3. Analysis results (det2_long_repos, 36,382 tracks, M3 chi2 < 5)

**Efficiency 95.9 % (any pad fired 96.8 %), 0/7 validation flags.** Timing 18.7 ns (multi-pad, invvar-weighted).

Journey (same data): 82.1 first look → 88.8 (conn-9 flip + dead pin row) → 89.3 (+13 noisy pads) → 92.4
(+ 20 mm clustering, shorted channels, dead-pad Voronoi cells) → 94.4 (+ 6 mm guard) → **95.9** (efficiency =
a cluster within 20 mm, not the max-charge cluster).

Remaining loss (% of tracks): only-distant hits 0.65 (≈0.5 M3 mis-pointing, 0.09-0.19 P2, 3 methods) ·
resolution tail 0.27 (mostly M3) · small-pillar shadow 2.0 · below threshold ≈0.8 · big pillars 0.26 · pad edges 0.14.
With a perfect reference ≈96.6 %.

11-24 series re-analysed with all fixes: alignment run 97.3 %; mesh scan 64 % (340 V) → 97.0 % (420 V);
drift scan 96-97.8 % from 705 V (96.2 % at 730). Old PDFs kept as `*_v1_flipped-map.pdf`; new `*_v2.pdf`.

Key physics findings:
- **Pointing** (M3 + alignment at the P2 plane): ≈3.0 mm/axis (core residual 4.5 mm minus pad term 3.4 mm);
  residual grows to 5.9 mm above 15° (projection to the mesh vs mid-gap charge).
- **Pillars are effectively bigger than drawn**: small r_eff 0.48 mm (Gerber 0.40), big R_eff 4.26 mm (Gerber
  3.0) — a dead ring at the pillar edge. Small-pillar loss is only visible vs **track angle** (no-hit 6 % at 0-2°,
  0.8 % above 20°); a distance-to-pillar test is blind to it (3 mm pointing ≫ 0.8 mm pillar).
- **"Wrong cluster won" (1.32 % under the old definition) = a second real in-time particle** (2.7 % of triggers
  have one; muon-like amplitude and time; ~14 cm away). Not crosstalk, not pointing. Hence the definition switch.
- **Readout crosstalk**: same-shape 0.85 % copy on Dream ch ±1/±2 (0.15 % ±3, 0 far); fires only for pulses
  >2500 ADC; harmless after clustering.
- Inter-pad gaps (Gerber) median 0.07 mm: negligible.
- Interior misses drop in the last 2 h before the trip (5.8 → 4.5 %): possibly gain change as the leak grew.

## 4. Code changes (P2_basket_analysis/cosmic_bench_analysis)

Config (`p2_qa_config.py`, all 4 det2 keys `det2_initial2/meshscan2/driftscan2/long_repos`):
`strategy_overrides`, `dead_channels`, `shorted_channels` (groups → `SHORTED_GROUPS`), `dead_guard_mm=6`,
`cluster_link_mm=20` (default for every run key), `t_max_h=7.58` (long_repos).

Behaviour changes that affect EVERY run when it is next reprocessed:
- **Clustered P2 position** (`p2_io.event_centroids`, link `CLUSTER_LINK_MM`, default 20 mm).
- **Efficiency = closest cluster within MATCH_R** (`p2_align.nearest_cluster`, stages 06/11/12/16/28).
- **Dead-pad cells** (Voronoi of dropped pads incl. dead connectors) + optional guard leave the active area
  (`p2_align.mask_dead_cells`, stages 06/10/11/12/16/28). Dead connectors 1/10 cells are now excluded too.
- Stage 10 masks at `ACTIVE_R` (was 15 mm).
- `MASK_GBR_PATH` resolves the main checkout from inside a worktree (pillar overlay was silently skipped).
→ det1/det3/det4/det5 numbers will shift (upwards) when re-run; re-run before quoting them next to det2.

New stages (in `run_p2_pipeline.sh` unless noted): 26 crosstalk waveforms (manual, needs decoded_root),
27 loss breakdown, 28 eff vs M3 chi2, 29 crosstalk/shorts illustration (manual, after 26), 30 miss anatomy,
31 pillar shadow vs angle (+ big-pillar blurred-disc fit), 32 mis-reco anatomy (max-cluster view),
33 reference-vs-P2 blame. All have pages in `build_final_pdf.py`.

Gotchas:
- Call `run_p2_pipeline.sh` by **absolute path** (its systemd-run re-exec fails on a relative `$0`, silently).
- In a worktree there is no `.venv`: put `P2_basket_analysis/.venv/bin` on PATH.
- Pedestal `hv_monitor.csv` column order differs from run CSVs: parse by header name.

## 5. Open items / next steps

1. **Bench**: inspect det2 mesh leak (burn mark? humidity? cable? try another HV channel), re-solder board
   connectors 2/3/8, check connector 5 row; only then retake the 3 h batch 2 (config ready).
2. Deploy the CERN note (`cd cern-site && scripts/deploy.sh`, password).
3. Re-run the other detectors' keys (det1/det5…) with the new clustering + dead-cell + any-cluster definitions
   before comparing them with det2; det5's earlier "pillars 1.09 pp" estimate was position-based and likely low.
4. Optional: sub-threshold-pulse study vs angle to separate pillar shadow from gain tail more cleanly; project
   the track to mid-gap (z + 2 mm) to shrink the inclined-track residual; FEU 7 ch 29 weak-pad check.
5. DAQ: no alarm on HV trips — the 01:07 trip went unnoticed for 5 h; a trip/stall watchdog in daq_control is
   still not implemented (also the 11-24 RunCtrl hang). Consider it before the next long run.
