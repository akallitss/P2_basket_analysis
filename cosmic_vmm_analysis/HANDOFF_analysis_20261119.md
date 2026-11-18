# Handoff: VMM cosmic efficiency analysis, to pick up on 19 Nov 2026

Read first: `FINDINGS_diag_20261118.md` (this folder) for what is known. The data taking that runs in parallel
from 19 Nov to Monday 23 Nov is in `DAQ_Control_VMM_Beam/docs/HANDOFF_WEEKEND_RUNS_20261119.md`.

## Rules

- **All analysis on lxplus** (interactive for light work, HTCondor for anything that decodes). Never on the DAQ box
  vmm_daplxa: it takes data, has 7 GB of RAM; there only read logs and monitor files.
- `ssh lxplus` works only while the user's ControlMaster is open; if you get "Permission denied (publickey,
  gssapi-with-mic,keyboard-interactive)", ask the user to run `! ssh -fN lxplus`.
- `condor_q` only shows the schedd of the node you land on: use `condor_q -name bigbird25.cern.ch` (or the schedd
  named in the job's `.log`). A job's stdout arrives only when it ends.
- On EOS, `mkdir -p` of deep paths fails: `cd` into an existing directory and `mkdir` there.
- Environment: `source /eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus/eff_pipeline/env.sh`
  (LCG_106: pandas, scipy; sets VMM_QA_DIR, EFF_RUNS_DIR, EFF_ANA).
- Code lives in P2_basket_analysis `cosmic_vmm_analysis/` on main. The EOS copies (`diag_20261118/code/`,
  `weekend_20261119/pipeline/cva/`) must be refreshed from git after every change
  (`tar c *.py | ssh lxplus 'cd <copy> && tar x'`); commit on main, one commit per file, no Co-Authored-By line.

## Where the data are

`A=/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus`

| path | content |
|---|---|
| `$A/diag_20261118/pseudo/chunks.json` | 1265 chunks of runs 231-267 (pid = run x 100 + k), configuration per chunk |
| `$A/diag_20261118/tags/tags_run<pid>.npz` | `--dump-diag` tag dumps (keys `p3_*`, `p2_*`, `p1_*` = IN, MID, OUT as probe; `rates_*`); `eff_run<pid>.json` (dead/masked lists, run-level after fix_dead) |
| `$A/diag_20261118/veto/veto_run_<run>.json` | spark_veto intervals and flags |
| `$A/diag_20261118/results/{all,veto,diag}` | counting estimator without / with veto; probe_diagnostics outputs |
| `$A/results/current/` | the incremental pipeline of 18 Nov (no --diag); its loop job ended 19 Nov 04:00 |
| `$A/weekend_20261119/` | the weekend pipeline (set up by the data session; runs >= 268, with --diag, veto, diagnostics) |

Per-tag keys: `px, py` (prediction, mm), `cx, cy, vmm, ch` (predicted pad), `on`, `dt_scint`, `slope`, `wall`,
`t_tag` (ns), `n_any`, `d_win`, `n100us`, `sat100us`, `prev_chip`, `prev_chan`; flattened per hit `fired_nr_i, _dt,
_adc, _d, _key` (key = vmm x 64 + ch), plus the standard `fired_i/x/y`, `side_i/x/y`. `r_scale` per chunk: weight
on = 1, off = -r_scale; scintillator window 140-310 ns on `dt_scint`. `probe_diagnostics.load_all()` builds the
pandas tables (tags, near hits) and is the easiest starting point.

## Tasks, in order

### 1. The "chamber silent" muons (dominant loss, ~10 %; 24-26 % on chips 2 and 14)

Goal: acceptance (edge / region not read out) or read-out (lost data)?
- x-y map of silent / all muons per chamber (weights on - r off), with the pads (`channel_efficiency.pad_geometry()`);
  distance of the prediction to the chamber's outer boundary and to the chip boundaries. If the silent muons pile up
  at the edge of chips 2 / 14 -> acceptance or geometry (check the fiducial: nearest pad live within `--fid-r` 9 mm).
- Time structure: silent fraction per chunk and per run; are silent muons clustered in time (seconds, minutes)?
  versus the per-second chamber rate (`rates_hits_p*`) and versus `prev_chip` (a chip that sent nothing for a long
  time before the muon points to a stalled read-out).
- Is the whole hybrid silent? Needs per-chip hit counts around the tag, not in the dump yet: add to
  `efficiency_tp.diag_fields` the hits per chip of the probe chamber within +-1 us and +-100 us (6 small int arrays
  per tag); re-dump one run first (one condor job, ~4 min), then all.
- Cross-check: does the per-chip counting efficiency follow the per-chip silent fraction?

### 2. Per-channel threshold in ADC and the muon spectrum above it (what the weekend 2D scan is for)

The user's criterion: **if at high gain the low tail of the muon ADC distribution sits above the threshold, gain is
not what limits the efficiency.** Needed:
- Per channel and run, the threshold in ADC = the lower edge of the ADC distribution of ALL hits of that channel
  (noise piles up just above threshold). Add to the dump a per-channel ADC histogram of all hits (32 x 64 x 64 bins
  of 16 ADC, int32, small) and take the edge (first percentile, or the 50 % point of the rising edge: ask the user
  which). Remember the ADC DNL (every 16th code ~2x wide, memory reference_vmm_adc_dnl): use 16-ADC bins aligned on it.
- For every detected muon, own-pad ADC minus its channel's edge; per chip and configuration the distribution of that
  distance. Fit the part above 0 (Landau x Gaussian, or the 9 mV/fC shape scaled by the gain ratio) and integrate
  below 0 = expected fraction lost under threshold.
- On the weekend gain x threshold grid: at 9 mV/fC the expected sub-threshold loss should be ~0 while the measured
  efficiency stays ~0.85 -> gain is not limiting. At 4.5 (and 3.0 if taken) compare the estimated loss with the
  (absent) efficiency difference.
- Threshold in charge: dist_mV / gain (mV/fC), once task 3 gives dist_mV for the gain files.

### 3. Threshold maps for the gain files

`build_thr.py` gives dist_mV = NaN for configurations without a rule in `maps/config_maps.json`. The per-gain maps
are on the DAQ box, `~/trim_work/out/gain_20261118/hvon_g6.0/map_eq.json` and `hvon_g9.0/map_eq.json` (equalised
start maps of the HV-on scans; eq files `p2b-config-cern-ext_gain{6.0,9.0}_peaktime200_opt_cosmic_gscan_eq.txt`).
The data session copies them to `$A/weekend_20261119/pipeline/maps/`. Add the rules (regex on the configuration
name -> map + reference ext file, as for the bscan files) and re-run build_thr / measure / regroup on both the
18 Nov and the weekend tags. Check: realised dist_mV of an R-target file vs the threshold bscan_target intended
(+0.55 mV median for the 4.5 files).

### 4. Per-chip outliers (chip 8 at 0.2-0.5; chips 3, 11, 12, 15 above 1)

Systematic, not statistics. Offsets of the fired pads (`fired_x/y`) from the prediction (`px/py`) per predicted
chip: a mapping error shows as a displaced or doubled peak; a chip whose hits are attributed to a neighbour shows an
excess near the chip boundary. Compare with `orientation_check.py` and the connector orders (`--order` options).

### 5. Pulser runs of the weekend (read-out efficiency)

The weekend files `*_tpA*` / `*_tpB*` pulse one channel on six chips (set A: chips 2, 5, 8, 11, 14, 17; set B: 3, 6,
9, 12, 15, 18; the exact channels are in the `st` lines of the ext files, copied next to every run) at 10.85 kHz, one
pulse per 4096-BC frame of 92.16 us. Write `pulser_readout.py`:
- per pulsed channel and chunk: pulses received vs expected (live time x 10850.3 Hz) and the gaps between consecutive
  pulses in SRS time (a missing pulse = a gap of 2, 3, ... frames); losses per chip, per hybrid, per FEC; random or
  in bursts; simultaneous on all chips of a hybrid?
- for the silent muons on pulsed chips: was the pulse of that chip in the muon's frame received? Pulses missing
  exactly when the chamber is silent = a read-out loss, quantified.
- exclude the pulsed channels and their two neighbours (crosstalk ~3.5 %) from the efficiency of those runs.

### 6. Weekend data, as it comes

The weekend pipeline updates `$A/weekend_20261119/results/current/` every 30 min. Check its log
(`$A/weekend_20261119/logs/`) once a day; re-run tasks 1-5 on the growing data set. The run plan and the meaning of
each configuration name: `DAQ_Control_VMM_Beam/docs/HANDOFF_WEEKEND_RUNS_20261119.md`.

### 7. Write-up

A FINDINGS for the weekend data next to `FINDINGS_diag_20261118.md`; a CERN-site explainer only if the user asks
(no analogies, schematics of the real thing, light theme only).

## Questions to bring to the user

- Which ADC-edge definition for "threshold in ADC" (first percentile vs 50 % of the rising edge)?
- If task 1 shows the silent muons are geometric (edge), remove them by a tighter fiducial?
