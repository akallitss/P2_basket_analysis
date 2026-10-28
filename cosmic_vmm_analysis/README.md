# Cosmic VMM analysis (Saclay three-chamber stand, October 2026)

Analysis of the VMM threshold campaign on the cosmic stand: P2_IN (det3, VMM 2-7), P2_MID (det2, VMM 8-13) and
P2_OUT (det1, VMM 14-19) lying flat, at the SPS spacing (z = 320 / 630 / 940 mm), read with Samtec cables (2 m)
instead of the Hitachi ribbons used at SPS. No external trigger.

The DAQ side (run control, config generator `quick_scripts/make_cosmic_sdt_configs.py`, mask lists, logbook,
observations log `docs/LAB_COSMICS_OBSERVATIONS.md`) lives in `DAQ_Control_VMM_Beam` (branch `lab-cosmics-oct26`),
which links this directory as `DAQ_Control_VMM_Beam/cosmic_vmm_analysis` (symlink to the sibling checkout).

Data: `/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/runs/run_N/cfg_*/raw_daq_data/*.pcapng`.
Heavy jobs run on lxplus; the scripts that need pandas want the LCG view
(`source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh`). The pad map comes from
`DAQ_Control_VMM_Beam/vmm_qa/vmm_stations.py` (the SPS map) through `vmm_qa_path.py`: set `VMM_QA_DIR`, or keep
the two repos side by side, or copy `vmm_stations.py` + `P2_BASKET_mapping.csv` into `./vmm_qa` (lxplus).
Outputs are kept outside the repo, under `data/Lab_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_*`.

## Scripts

| script | what |
|---|---|
| `coincidence_check.py` | shared pcapng reader `read_times(path, state, with_adc=False)` (absolute hit time) + chamber-pair time differences |
| `noise_scan_summary.py` | per-channel hit rate of every point of a threshold scan |
| `compare_hv_noise.py` | HV on vs HV off per-channel rates |
| `neighbor_baseline.py` | per-channel ADC baseline from neighbour read-out (sng 1) |
| `baseline_vs_peak.py` | does the baseline predict the noise peak position |
| `threshold_table.py` | per-chip noise-onset thresholds and per-channel baselines tables |
| `adc_spectra.py` | ADC spectra per chip |
| `burst_anatomy.py` | what the bursts of VMM 4/7/9 look like |
| `event_multiplicity.py` | time groups, three-chamber "muon-like" proxy |
| `tag_and_probe.py` | chamber-level tag-and-probe (no geometry) -- superseded, tags are noise-dominated |
| `tag_and_probe_pads.py` | geometric tag-and-probe with the pad map (method of `sps_beam_analysis/22_tag_probe_efficiency.py`) |
| `muon_signal_check.py` | ADC by event class and chip-pair time coincidences (map independent) |
| `dump_candidates.py` | keep the hits of small multi-chamber groups (input of the mapping / orientation studies) |
| `orientation_check.py` | channel order per connector (adjacency, smoothness, straight tracks) and rotation/mirror scan between chambers |
| `straight_track_rate.py` | straight three-chamber tracks above chance per run (group and chip-preserving nulls); `--conn-shift`, `--order`, `--veto` |
| `event_display.py` | event displays of random straight / unselected groups, det2 residual against chance |
| `residual_diagnostics.py` | what widens the det2 residual: pad size, z spacing, rotation, per connector half |
| `half_order_scan.py` | channel order per connector from straight tracks, chosen on half the groups, scored on the other half |
| `inventory_runs.py` | one row per sub-run: settings, masks, HV declaration, connectors, sampled hit rate per VMM |
| `run_metadata_corrections.json` | facts that override the run_config.json texts (HV connectors NOT moved on 26 Oct) |
| `vmm_qa_path.py` | puts `vmm_stations` on `sys.path` |

## The CERN-site note

`note/` holds the generator of the unlisted note "Measuring the VMM noise threshold with cosmic muons: what we did and
what we saw" (cern-site `notes/unlisted/2026-10-22-...html`): `note/build_note.sh` rebuilds it from the analysis outputs,
then `scripts/build_hub.py` + `scripts/deploy.sh` in cern-site publish it. Figures are inline SVG written by the parts.

## Map for the lab runs

Found from the data (`straight_track_rate.py` connector-shift scan, `half_order_scan.py` with held-out halves). The
run_config.json texts say c4-c6 for runs 1-111, which is NOT what the data show:

| runs | hybrids H1-H3 / H4-H6 / H7-H9 on | options |
|---|---|---|
| 1-111 | c5-c7 | `--conn-shift 1 --order 'P2_IN:5=bot:rev\|top:lin'` |
| 112- | c3-c5 | `--conn-shift -1 --order 'P2_IN:3=bot:rev\|top:lin' --order 'P2_IN:5=bot:rev\|top:lin'` |

The linear top halves belong to det3 (P2_IN) connectors c3 and c5, whichever hybrid reads them (c5 top is linear when
read by VMM 3 in runs 1-111 and by VMM 7 from run_112); all other halves keep the SPS order. On the c4-c6 map the old
runs show nothing (about 10 straight tracks/h against 511 +- 10/h at shift +1, pooled 76-91): the 25 Oct "no muon
signal" was this map error. Always scan the connector shift before trusting a null result.

Discharges: check repeated pad triplets before trusting a track count. In run_121 det3 VMM 2 ch 27 and VMM 6 ch 44
fire saturated every few seconds and induce small hits on fixed det2/det1 channels: 131 identical fake "tracks" that
pass both chance nulls (`--veto 2:27,6:44`). Run_124 has none. Runs 76-91: det2 VMM 11 ch 35 saturated with pickup on
VMM 4 ch 9 and VMM 18 ch 10 (`--veto 11:35,18:10,4:9`).

Tag-and-probe trap: noise firing two chambers at once passes the off-time subtraction and has a cosmic-like tag angle
distribution (run_76 on the wrong map: 98 % "purity", 2 % efficiency).

## Status (2026-10-28)

Muons are seen in run_121 and run_124 (HV P2_IN 420/720, P2_MID/OUT 430/730, c3-c5): straight tracks above the
chip-preserving chance level 604 +- 43 /h (run_121, VMM 4/7/9 parked, discharge vetoed) and 813 +- 49 /h (run_124, all
chips live). det2 residual core 5.7 / 7.8 mm against 4.1 / 4.3 mm from the pads alone; det2 ~9 mm from half way in z,
rotations < 0.4 deg. Outputs: `data/Lab_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_run121`, `analysis_run124`,
`analysis_inventory`. Next: efficiency per chamber, run by run.

## Status (2026-10-25)

No measurable muon signal at any threshold setting (DAQ repo `docs/LAB_COSMICS_OBSERVATIONS.md` items 61-64): the
SPS map is the best wiring found (no alternative beats it on held-out events, no mirror), the time coincidences
between chambers are shared noise (chip-level), the bursts are whole-chamber events, and the low-threshold
"tracks" are a fixed pattern of hot channels. Next: HV / gas at the lab vs SPS, and a source test.
