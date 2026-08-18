# Roadmap §5 → figures, producers, products

Where every wish-list item lives. Figure paths are relative to
`conference/figures/`; producers and products are relative to
`~/Documents/PostDocSaclay/`. Re-gather with `python3 gather_figures.py`.

Three trees matter:

| what | where |
|---|---|
| **figure producers** | `P2_basket_analysis/mpgd2026/make_talk_figs.py` (cross-campaign), `P2_basket_analysis/sps_beam_analysis/mpgd26_figs/fig_*.py` (beam), `P2_basket_analysis/cosmic_bench_analysis/22_conference_figures.py` (bench), `nTof_x17/mpgd26/make_coincidence.py` (3D) |
| **stage products** | EOS `/eos/experiment/ntof/data/x17/p2_sps_july/analysis/`, mirrored into the figure workspace `data/SPS_Beam_Test/mpgd26_workspace/products/analysis/` |
| **uRWELL-referenced results** | LaCie `…/TB_July2026_H4/analysis/urw_referenced_efficiency/<run>/urw_p2_efficiency_<run>.csv` |

---

## 5a — comparison measurements (bench ↔ beam, DREAM ↔ VMM)

| item | figure | producer | product behind it |
|---|---|---|---|
| **5a.1** DREAM vs VMM, same chambers | `4_act3_beam/s23_dream_vs_vmm`, `s24_vmm_threshold` | `make_talk_figs.py` → `fig_dream_vs_vmm`, `fig_vmm_threshold` | `sps_beam_analysis/vmm_dream_matching/efficiency_table.csv`, `thresholds_run_46.json`, `model_P2_OUT.json` |
| **5a.2** ε(HV) bench ↔ beam | `4_act3_beam/s21_bench_beam_mesh` (gain) and `4_act3_beam/bench_beam_drift` (transport) | `make_talk_figs.py` → `fig_bench_beam_mesh`, `fig_bench_beam_drift` | bench: `data/Cosmic_Bench/Analysis/det*/…/11_hv_scan_efficiency/`, `…/16_drift_scan_efficiency/`; beam: `urw_referenced_efficiency/drift_mesh_scan_1/` |
| **5a.3** VMM config optimum | `5_optional/s35_vmm_config_scan`, `2_act1_nov2025/s08_snr_matrix` | `mpgd26_figs/fig_vmm.py`; `make_talk_figs.py` → `fig_snr_matrix` | `reports/mpgd26_sps_beam_2026-08/data/vmm_config_scan_table.csv`; `data/SPS_Beam_Test/VMM-alinx-data/vmm_snr_results.csv` |
| **5a.4** efficiency map bench vs beam | `4_act3_beam/bench_beam_maps` | `make_talk_figs.py` → `fig_bench_beam_maps` | bench `06_efficiency/efficiency_map_sliding_*.npz`; beam `22_tag_probe_efficiency/eff_map_*.csv` |
| **5a.5** timing beam vs bench | `4_act3_beam/s22_timing_campaigns`, `s22_timing_vs_drift_magboltz`, `5_optional/s27_timing_ladder` | `make_talk_figs.py` → `fig_timing_campaigns`; `mpgd26_figs/fig_timing.py` | `report_data/dream_timing_scans.csv`, `dream_timing_persubrun.csv`, `data/gas_timing_model.csv` |

Extra, produced 2026-08-17/18 and not on the wish list:
`5_optional/eff_2d_curves_*`, `timing_2d_curves_*`,
`timing_2d_curves_vs_drift_*`, `timing_drift_choice_table_*`.

---

## 5b — physics / telescope measurements

| item | state | figure | producer | product |
|---|---|---|---|---|
| **5b.1** 3-plane 2-of-3 tagging | **done**, and superseded by the absolute uRWELL measurement | `5_optional/s32_eff_methods_comparison` (the three methods side by side), `4_act3_beam/s17_urw_panels_P2_MID`, `4_act3_beam/s20_padmap_working_point` | stage `sps_beam_analysis/22_tag_probe_efficiency.py`; figures from `mpgd26_figs/fig_efficiency.py`, `fig_padmap.py` | `analysis/<det>/<run>/scan/22_tag_probe_efficiency/tag_probe_efficiency_spark_vetoed.csv` and the per-sub_run `eff_map_<det>_*.csv` |
| **5b.2** 3-plane event display | **done** | `4_act3_beam/s18_coincidence_side`, `s19_coincidence_hero`, `9_spare/coincidence_beam`, `4_act3_beam/s19_event_display` | `nTof_x17/mpgd26/make_coincidence.py`, event chosen by `tools/extract_coincident_event.py` | `nTof_x17/mpgd26/data/coincident_events.json` (25 events, subrun 23) |
| **5b.3** plane-to-plane spatial residuals | **done** — rms 3.37 / 3.44 mm = 12 mm/√12 | `4_act3_beam/s17_urw_panels_P2_MID` (panels 1–2: 1D and 2D residual, the tilted pad square) | `sps_beam_analysis/urw_reference/urw_p2_efficiency.py` | `urw_referenced_efficiency/highstat_eff_1/` — `res_x_mm`, `res_y_mm` in the CSV, panels in `per_subrun/` |
| **5b.4** plane-to-plane time residuals | **done** | `5_optional/s27_timing_ladder`, `s27_timing_vs_mesh_bothruns`, `5_optional/timing_2d_curves_vs_drift_*` | stage `sps_beam_analysis/29_waveform_timing.py`; figures from `mpgd26_figs/fig_timing.py` | `analysis/telescope/<run>/<sub>/29_waveform_timing/waveform_timing_summary.json` → `report_data/dream_timing_persubrun.csv` (`kind='pair'` rows are the plane-to-plane σ) |
| **5b.5** track efficiency & fake rate vs rate | **done, muons only** — within sub-runs the slope is +0.6 ± 0.4 (MID) and +0.5 ± 0.3 (OUT) points/kHz, i.e. no rate dependence up to ~700 Hz of reference tracks ≈ 2.2 kHz of trigger; fake-match 0.05–0.11 % per 10 mm | `4_act3_beam/rate_performance` | `mpgd26_figs/fig_rate.py` | per-time-bin `efficiency_vs_time` in `urw_timebins/<run>/<sub>/urw_p2_efficiency_*.json`; `accidental_per_10mm` in the CSVs |
| **5b.6** angle scan | **not taken** — outlook line | — | — | — |
| **5b.7** multi-track separation, pions | **not taken** — outlook line | — | — | — |
| *new* uRWELL track vs P2-only track | **measured** — 508 coincidences (80.9 % of tracked events); the P2-only track localises to 3.3 mm at the middle plane (= 12/√12) but **cannot measure angle**: one pad step over the 620 mm lever arm is 19 mrad against a 1.2 mrad beam divergence, so 73 % of coincidences fire the same pad on all three planes | `4_act3_beam/p2_standalone_tracking` | extraction `nTof_x17/mpgd26/tools/extract_track_pairs.py` (lxplus, `~/mpgd26_tracks/run_pairs.sh`); figure `mpgd26_figs/fig_p2_tracking.py` | `nTof_x17/mpgd26/data/track_pairs.csv` (508 rows) and `track_pairs.json` (40 display events) |
| *new* 3D display of both tracks | **pending** — extraction done, scene not rendered | — | `nTof_x17/mpgd26/make_coincidence.py`, to be extended with the P2-only line | `data/track_pairs.json` |

---

## The other things a reader will ask for

| topic | where |
|---|---|
| chamber ↔ station history (which detector, which days) | `sps_beam_analysis/chamber_history.py` |
| how the beam figures are rebuilt | `sps_beam_analysis/mpgd26_figs/README.md` |
| the uRWELL condor pipeline | `sps_beam_analysis/condor/urw/README.md` |
| the DREAM stage pipeline | `sps_beam_analysis/PIPELINE_OVERVIEW_2026-07-26.md`, `condor/README_CONDOR.md` |
| bench pipeline and its conference figures | `cosmic_bench_analysis/22_conference_figures.py`, `run_p2_pipeline.sh` |
| campaign logbook notes | LaCie `…/TB_July2026_H4/analysis/2026-07-*.md` (also on EOS `analysis/`) |
| the full written report | `reports/mpgd26_sps_beam_2026-08/index.html` |
| the deck | `cern-site/notes/2026-08-17-mpgd26-talk-draft-deck.html` |
