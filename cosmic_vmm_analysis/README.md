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
| `vmm_qa_path.py` | puts `vmm_stations` on `sys.path` |

## Status (2026-10-25)

The geometric tag-and-probe does not find tracks yet: the time coincidences between chambers are real (narrow
peaks) but spread over all chip pairs, and the SPS channel order inside the connector halves does not hold for
every half with the Samtec cables (the top halves of c5 / c6 do follow it; the P2_IN ribbon overrides needed with
the Hitachi cables do not). Next: per-half channel-order scan on clean events, then the chamber-to-chamber
rotation + mirror scan (as `cosmic_bench_analysis/03_m3_alignment.py`), then the efficiency.
