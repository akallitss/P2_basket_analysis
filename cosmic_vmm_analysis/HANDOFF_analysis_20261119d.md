# Handoff: VMM cosmic efficiency analysis - overview at 19 Nov 2026 ~20:30

This is the single entry point for the next session. It merges the two analysis sessions of 19 Nov:
- session 1bd4ab0b ("analysis session"): `HANDOFF_analysis_20261119.md` (rules, data, tasks), `..._20261119b.md`
  (morning findings), `..._20261119c.md` (evening: tick fix + re-dump, fiducials);
- session 08f0123f ("dead-time session", this file): dead time, sparks, pulser read-out, reconstruction-free and
  connector-matched efficiencies, comparison with the beam test and the Dream M3 bench.

Read this file first; go to a/b/c for details (paths, commands, per-configuration tables). Where this file and c
disagree, this file is the later state (one correction: chip 8, finding 4 below).

Data taking continues until **Mon 23 Nov 07:00** (supervisor on the box, `DAQ_Control_VMM_Beam/docs/HANDOFF_WEEKEND_RUNS_20261119.md`);
the lxplus loop (condor 13752682, bigbird25) re-measures every 30 min until **Mon 18:00**.

## 1. The question and where it stands

**Question (user):** the VMM-read P2 chambers reach ~0.85-0.90 on the cosmic stand, Dream reaches ~0.93-0.98 on the
same chambers. Is the gap the VMM electronics, or how we compute the efficiency?

**Answer so far: it is the electronics settings, not our efficiency calculation.** The chain of evidence:

1. Read-out losses, dead time and sparks are too small (section 4, findings 1-3): together < 0.5 % of the muons.
2. Reconstruction-free: of the muons through P2_IN and P2_OUT (scintillator-confirmed), **P2_MID records no hit at
   all in 10-15 %** (any channel, any position; finding 5). No reconstruction choice can recover those.
3. **Same connectors, cosmic VMM = beam VMM, both ~8-10 points below Dream** (finding 6): P2_OUT connector 5
   VMM 0.86-0.88 here vs 0.872 at the SPS beam (uRWELL tracks); Dream 0.94-0.96 on the same pads (beam and M3 bench).
4. The beam autopsy (`sps_beam_analysis/vmm_dream_matching/EFFICIENCY_AUTOPSY.md`, run_46) already traced the beam gap
   to the **VMM discriminator threshold at ~0.7x the most probable signal**, plus pad-to-pad signal-size variation
   (`sps_beam_analysis/dream_vmm_comparison/README.md`: one threshold cut into Dream's per-pad spectra reproduces the
   VMM efficiency pad by pad, r = 0.90).
5. Why raising the VMM gain did not help on cosmics: the per-gain calibration sets every channel at a fixed **noise
   rate** (R = 2/10/50 Hz, `bscan_target.py`), i.e. a fixed number of noise sigmas; VMM noise scales with the gain,
   so the threshold in charge stays ~4 fC (4.5, 6.0) or rises (9.0: 5-7 fC, trims clip on 341-375 of 1152 channels).
   The gain scan could not lower the threshold in charge by construction. (At the beam the DAC was fixed in mV, so
   gain raised signal over threshold: 0.72 -> 0.82.)

**User's view (19 Nov evening):** the next real step is a measurement of the electronics transfer function (per
channel ADC/fC and threshold in fC with the test pulse) - **not now**. Until then the analysis question is the
efficiency calculation itself; keep it decoupled from the reconstruction.

## 2. Rules (from a, still valid)

- All analysis on **lxplus** (interactive for light work, HTCondor for decoding). Never compute on the DAQ box
  vmm_daplxa (7 GB RAM, takes data); there only read logs / monitor files.
- Memory-cap every process: `systemd-run --user --scope -q -p MemoryMax=3G bash -c 'source .../env.sh; ...'` (the env
  must be sourced INSIDE the scope, and PYTHONPATH prepended, not overwritten).
- `ssh lxplus` needs the user's ControlMaster (password + OTP). On "Permission denied" ask the user for `! ssh lxplus true`.
- On EOS `mkdir -p` of deep paths fails: mkdir one level at a time. `cp -p` fails on EOS (use plain cp).
- Code: P2_basket_analysis `cosmic_vmm_analysis/`; **one commit per file, no Co-Authored-By line**; EOS copies are
  plain copies (no git), refresh by scp/tar after a commit, check against git first.
- Never mask / skip chips or make a fiducial the default without asking (memory feedback-vmm-study-keep-all-chips).
- The user decides runs and HV; never touch HV or the DAQ.

## 3. Code state

Branch **`decoder-tick-225`** of P2_basket_analysis (worktree `.claude/worktrees/decoder-tick-225`, pushed to GitHub,
**not merged into main**) now contains everything of both sessions:

| commit(s) | what |
|---|---|
| 8b9f493 / 67145d6 | `coincidence_check.read_times`: SRS marker tick **22.5 ns** (was 25 ns; same change in both sessions) |
| 466a107 | `pulser_readout.py`: lattice from the **nominal** frame period (92.16 us, exact) and the **peak** of the residuals, re-centred; the old median-based fit walked off the pulses on noisy / ns-rounded channels (chips 8 and 11 at 0.03) |
| 85c6efd | `channel_efficiency measure --fiducial` (pad-level) - **biased, do not use for results** (c, finding 4) |
| 3bbcf0b | `fiducial_muons.py` - muon-level fiducial A10/A20/A30 (the valid one) |
| 215e46e | `HANDOFF_analysis_20261119c.md` |
| 9f0519c | merge of `tick225-fiducial` into `decoder-tick-225` |
| ab39cb7 .. a0c6cee | this session's scripts (one commit each, table below) |
| (this file) | `HANDOFF_analysis_20261119d.md` |

`tick225-fiducial` (other session) is fully contained; merge **`decoder-tick-225` into main with the user** (after the
weekend, together with decision 3 of c: lxplus git clone + VERSION file for the EOS copies).

This session's scripts (all run from `$W` with `pipeline/env.sh` and `PYTHONPATH=<cva copy>:$PYTHONPATH`):

| script | what | outputs |
|---|---|---|
| `spark_gap_chunk.py` + `deadtime_sparks_plots.py` | per muon: time since the previous chip hit, since/until sparks (>= 1/3/10 saturated hits in 1 us), hit-free interval around the muon vs random times | `$W/analysis/deadtime_sparks/results/fig1-3, summary.txt` |
| `spark_profile_chunk.py` + `profile_pulser_plots.py` | hit rate after each spark (same chamber, others, test pulses) 300 ns - 1 s | `.../deadtime_sparks/results/fig4_spark_profile.png` |
| `pulser_summary_chunk.py` (+ `profile_pulser_plots.py`) | test pulses vs the real pcapng duration per file | `.../deadtime_sparks/results/fig5_pulser.png, pulser_summary.csv` |
| `file_clock.py` | decoded span vs packet span per capture file (the tick check) | stdout |
| `scint_any_hit.py` | any chamber hit per scintillator coincidence, alone / given the other two, several windows, random-corrected | `$W/analysis/pulser_tick225/scint_any_hit.log` |
| `conn_eff.py` | per predicted connector: hit within 15 mm / +-150 ns (beam definition) and any-chamber-hit | `.../pulser_tick225/conn_eff.log` |
| `window_check.py` | efficiency vs probe half-window, from the dumps | stdout (numbers in finding 7) |

NB: the `deadtime_sparks` outputs were made **14:32-14:53 with the old 25 ns clock** for both tags and decoding
(consistent; log-scale timing results hold, rates per second ~10 % low). Everything in `pulser_tick225/` is 22.5 ns.

EOS code copies (`A=/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus`, `W=$A/weekend_20261119`):
- `$W/pipeline/cva/` (weekend loop + `post_t225.sh`): 22.5 ns `coincidence_check.py`, `channel_efficiency.py`,
  `fiducial_muons.py` (from c), and since 20:15 the fixed **`pulser_readout.py`** (backup
  `pulser_readout.py.bak-20261119-prelattice`). The other new scripts are not there (they live in the repo and in
  `$W/analysis/pulser_tick225/` / `deadtime_sparks/`).
- `$W/analysis/pulser_tick225/cva/` = `$W/pipeline/cva` + fixed `pulser_readout.py` (used for the fixed pulser runs).
- `$A/diag_20261118/code/`, `$W/redump_test/cva/`: 22.5 ns (c); `pulser_readout.py` there still the old lattice.
- `$A/eff_pipeline/cva/`: stale (25 ns), unused.

LaCie copies (rsync `--no-o --no-g --no-p`): `.../vmm_threshold_scan/analysis_weekend_20261119/analysis/{deadtime_sparks,pulser_tick225}/`.

## 4. Findings (both sessions, consolidated)

1. **Read-out is not losing data.** Pulser run_270 (set A, ch 31 of chips 2 5 8 11 14 17), 22.5 ns + fixed lattice,
   49 files: **99.82-99.87 % of test pulses received on all six chips** (`pulser_tick225/results_run270.log`);
   missing frames are independent per chip (53 of 31.9 M frames lost on >= 5 chips at once). Around the silent muons
   both pulses arrived 98.7-100 %. Against the packet clock (`deadtime_sparks/results/pulser_summary.csv`) 99.84-99.91 %.
   **Run_283 (set A at gain 9, `gscan_hvon10_hvm_tpA`, 49 chunks, done 20:45, `pulser_tick225/results_run283b.log`):
   chips 2, 14, 17 at 99.86 % as in run_270, but chip 5 0.970, chip 11 0.988 and **chip 8 0.877** (3.9 M of 31.9 M
   pulses missing, only 0.47 M of them in runs of >= 5 frames). Around silent muons the pulses of chip 8 were
   received 0.81, chip 5 0.94. PRELIMINARY: mostly isolated missing pulses -> candidate = channel busy with noise hits
   at gain 9 (chips 8 and 11 were the noisy ones in the 18 Nov coarse probe); check the hit rate of 5:31, 8:31,
   11:31 in run_283 vs run_270 and the gaps (pulser_gaps.csv) before calling it a read-out loss. At gain 9 this is
   part of the lower MID efficiency (0.70).
2. **Dead time is real but small.** Silent fraction vs time since the previous hit on the predicted chip, against
   that chip's own rate in the run: excess only below ~10 us (0.58 at < 100 ns, 0.24 at 1-3 us vs ~0.15) = **83 of
   28 480 muons, 0.29 %** (2 % of the silent ones). The < 100 ns bin partly counts the muon's own early hit
   (prev_chip is measured from tag + latency - win). (c finding 6 sees the same with the fixed clock.)
3. **Sparks do not blind the chamber.** After a spark (>= 3 saturated hits in 1 us; 0.2-1 /s per chamber) the chamber
   is **flooded** (10-100x its rate) for ~1 ms, then normal; other chambers rise too. Only after big sparks
   (>= 10) in P2_IN do the test pulses drop to ~0.6 for 1-10 ms (~0.03 % of the time). Muons 1-100 ms after big
   sparks: silent 0.24-0.28 vs 0.15 (~230 muons) - ~0.1 % of all muons. The HV-monitor spark veto changes eps by
   <= 0.02 (c finding 2). Silent muons do not cluster in time; silent and detected muons sit in equal hit-free
   intervals (no "read-out goes dark" signature).
4. **Correction to c finding 3: chip 8 ch 31 is NOT losing pulses.** Its 2.9 % came from the old lattice fit (noise
   pairs on that channel biased the fitted period by ~5 ns -> the lattice drifted off). With the fix chip 8
   receives 99.86 % (as the others). Chip 8's low efficiency (0.2-0.5, task 4) still needs an explanation - but not
   read-out. Do not chase 8:31 pulses on the box.
5. **Reconstruction-free efficiency** (`scint_any_hit.py`, 15 files per run, window 140-310 ns after the
   scintillator): P(probe has any hit | scint and the other two chambers hit), random-corrected:

   | run (config) | IN | MID | OUT | | P(hit given scint alone) IN / MID / OUT |
   |---|---|---|---|---|---|
   | 268 REF g4.5 | 0.75 | **0.85** | 0.73 | | 0.56 / 0.56 / 0.54 |
   | 269 REF g4.5 | 0.72 | **0.90** | 0.70 | | 0.50 / 0.55 / 0.50 |
   | 275 B10 g4.5 | 0.76 | **0.87** | 0.73 | | 0.56 / 0.55 / 0.51 |
   | 278 g6.0 | 0.74 | **0.87** | 0.74 | | 0.53 / 0.55 / 0.54 |
   | 281 g9.0 | 0.70 | **0.70** | 0.67 | | 0.49 / 0.42 / 0.47 |

   MID is bracketed (cleanest). IN/OUT include extrapolation acceptance. "Scint alone" ~0.5 = about half of the
   paddle coincidences do not cross a chamber's read-out area (paddle acceptance) - never compare that column with
   a Dream "per trigger" number. Chamber hits arrive 165-265 ns after the scintillator (5-95 %): widening the window
   to 2 us adds <= 1 %.
6. **Connector-matched comparison** (user: P2_IN = det3, P2_MID = det2, P2_OUT = det1; the beam read connectors 4,5,6
   with chips 2-7 / 8-13 / 14-19 in that order, the cosmic stand reads 3,4,5 (`--conn-shift -1`), so the same pads
   sit on chips one hybrid later: cosmic c3 = chips 2,3/8,9/14,15, c4 = 4,5/10,11/16,17, c5 = 6,7/12,13/18,19).
   "any" = any chamber hit; VMM from `conn_eff.py` (REF, 370-750 muons per cell, +-0.015); Dream M3 bench from
   `cosmic_bench_analysis` stage 06 `ray_hit_miss_list*.csv` (`has_any`, active area; `pad_idx // 128` -> the kept
   connectors in ascending order, verified); beam from `sps_beam_analysis/dream_vmm_comparison/data/compare_dream_vmm_P2_OUT.csv`:

   | | c3 | c4 | c5 |
   |---|---|---|---|
   | P2_OUT Dream M3 bench (det1 19 Jul, mesh 415) | 0.984 | 0.968 | 0.943 |
   | P2_OUT VMM cosmic (mesh 430) | 0.888 | 0.875 | 0.860 |
   | P2_OUT beam Dream (eff_nominal_1) | - | (no tracks) | **0.960** |
   | P2_OUT beam VMM (run_46) | - | 0.844 | **0.872** |
   | P2_MID Dream M3 bench (det2 9 Jul) | 0.923 | 0.938 | 0.942 |
   | P2_MID VMM cosmic | 0.886 | 0.876 | 0.869 |
   | P2_IN Dream M3 bench (det3 16 Jul, short run) | 0.891 | 0.890 | 0.902 |
   | P2_IN VMM cosmic | 0.879 | 0.912 | 0.798 |

   Cosmic VMM = beam VMM on the shared connector (chips 16,17 read c5 at the beam: 0.872; c4 here: 0.875). Gap to
   Dream: OUT 8-10, MID 4-7 points; IN close on c3/c4, c5 (chips 6,7) ~10 low. The Dream bench numbers depend
   strongly on which connectors are excluded (det1 85 % -> 97 % by also dropping connector 2): always quote the
   connector set.
7. **The +-100 ns probe window of the tag-and-probe costs 1-4 %** (live pad within 15 mm: IN 0.657 -> 0.700, MID
   0.731 -> 0.754, OUT 0.658 -> 0.670 from +-100 ns to +-2 us; side-band not subtracted). Since chamber hits are tight
   with respect to the scintillator (finding 5), this is **tag-time jitter** (tag = earliest hit of a clean cluster),
   not late VMM hits. Fix candidates: a wider window, or the tag time from the scintillator.
8. **Per-chip / per-region eps of the counting estimator are not absolute** (both sessions independently): pad-level
   fiducial gives eps > 1 (c finding 4), and per-chip eps of the middle connector are 0.95-1.16 (`by_chip_config.csv`).
   The estimator takes expected crossings from the density of the smeared predictions; only the full chamber
   balances. Use the muon-level A20 (`fiducial_muons.py`) or the any-hit / within-R numbers for regions.
9. **Muon-level fiducial** (c finding 5): OUT +4-5 points with a 10-20 mm cut and flat inside; MID flat; **IN rises
   towards the centre** (0.69-0.73 at 10-40 mm, 0.80 at >= 60 mm) - a gradient, not only an edge.
10. **Tick fix** (c finding 1): efficiencies unchanged to <= 0.002; wall clock, rates, gaps and the pulser lattice
    were wrong before. The SPS matching code had found the same 10/9 (`vmm_dream_matching/README.md`).
11. **How tags handle several clusters** (efficiency_tp): a tag chamber's "cluster" is a TIME group (gap 100 ns), so
    a muon plus any other hit on the chamber within ~100 ns, > 6 pads, a dead/masked pad, or a stack burst (> 36
    hits) removes the muon from the sample (not counted as a miss). The probe has no cluster requirement. Muons are
    therefore tested in quiet moments only; no bias seen yet, but untested.

## 5. Next steps (suggested order)

1. **Gain-9 pulser losses** (section 4.1, run_283: chip 8 0.88, chip 5 0.97, chip 11 0.99): noise-busy or real? Then set-B runs (tpB, chips
   3 6 9 12 15 18) as they appear (`post_t225.sh w_pulser <run>` now uses the fixed `pulser_readout.py`).
2. **Matched-definition table for the write-up**: for every configuration of the weekend grid, MID any-hit given IN
   and OUT (`scint_any_hit.py` on all files, not 15) and per connector (`conn_eff.py`) next to the counting eps,
   with the Dream bench numbers on connectors 3-5. Add the side-band subtraction to `conn_eff.py` /
   `window_check.py` (shift the window by 4 us, as scint_any_hit does).
3. **Tag time from the scintillator** instead of the earliest cluster hit (finding 7) in `efficiency_tp.py` (or as a
   variant), and re-measure: does the 1-4 % window loss disappear?
4. **Multi-cluster tags** (finding 11): count the scintillator-confirmed muons rejected because the tag chamber had
   more than one spatial cluster; cluster in space inside each time group and keep exactly-one-muon-like-cluster tags;
   compare the efficiency of the recovered muons.
5. Chip 8 (eps 0.2-0.5, not read-out): offsets of fired pads vs prediction per chip (task 4 of a), its threshold
   (`thr.csv`), dead/masked neighbours.
6. P2_IN gradient (c next step 2): x-y map of A20.
7. Remaining c/b items: whole-hybrid silence (`chip_n1us`/`chip_n100us`), ADC-threshold fit replacement, FINDINGS
   for the weekend.
8. After Mon 18:00: merge `decoder-tick-225` into main with the user; lxplus git clone + VERSION file.

**Not now (user):** the electronics transfer function (test-pulse amplitude scan per channel -> ADC/fC and threshold
in fC per channel and gain), HV changes, new runs.

## 6. Watch points

- The loop re-measures only when inputs change; force with `rm $W/results/current/.inputs_sig`. Edit
  `$W/pipeline/run_new.sh` only between passes (no `$W/state/lock.d`), replace with mv (c).
- `condor_q -name bigbird25.cern.ch` for the loop; condor `arguments = $(step) $(arg)` without quotes.
- Background jobs started with `nohup setsid` on an lxplus node can vanish (run_pulser_fixed.sh waiter); prefer
  condor for anything that must survive, and check outputs, not processes.
- The local `DAQ_Control_VMM_Beam/vmm_qa/vmm_stations.py` lacks `MAP_CSV`: run mapping code on lxplus with the EOS
  pipeline's `vmm_qa` (or set `VMM_QA_DIR` to it).
