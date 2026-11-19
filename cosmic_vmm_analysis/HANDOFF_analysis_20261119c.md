# Handoff: VMM cosmic analysis, state at 19 Nov 2026 ~20:15 (continues HANDOFF_analysis_20261119b.md)

Read `HANDOFF_analysis_20261119.md` (rules, data locations, task list) and `HANDOFF_analysis_20261119b.md` (morning
findings) first. This file covers the afternoon/evening of 19 Nov: the user's two decisions from 20261119b were taken
and carried out, plus one new estimator. Data taking continues until Mon 23 Nov 07:00; the lxplus loop until Mon 18:00.
This session (1bd4ab0b) was asked by the user to run the analysis until the end of the weekend scan.

`A=/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus`, `D=$A/diag_20261118`,
`W=$A/weekend_20261119`. Condor scripts of this session: lxplus `~/weekend_20261119/tick225/` (logs in `logs/`).

## User decisions (19 Nov)

1. **Decoder SRS tick 25 -> 22.5 ns: fix and re-dump everything.** Done (below).
2. **Fiducial: "several, compare"** (5/10/15/20 mm, all edges, IN and OUT), always next to the unchanged full-chamber
   result. A 20 mm pad cut drops 170/384 pads per chamber, incl. ~60/64 of chips 2, 7, 14, 19 (the top/bottom bands) -
   do not make any fiducial the default without asking (memory feedback-vmm-study-keep-all-chips).
   After the pad-level version proved biased (below) the user said: try the fiducial on the muons and compare. Done.
3. After the weekend (not before the loop ends Mon 18:00): git clone of P2_basket_analysis on lxplus, pipeline scripts
   into the repo, EOS code copies refreshed from a commit with a VERSION file (see "Code" below for why).

## Code (P2_basket_analysis branch `tick225-fiducial`, pushed to GitHub, NOT merged)

Worktree `P2_basket_analysis/.claude/worktrees/tick225-fiducial`. Branch made from local main 7a66ba1 (main was 5
commits ahead of origin; the push carried those commits too, GitHub main untouched).

| commit | file | what |
|---|---|---|
| 67145d6 | `coincidence_check.py` | `SRS_TICK_NS = 22.5` in `read_times` (was `srs_h * 25.0`) |
| 85c6efd | `channel_efficiency.py` | `measure --fiducial CHAMBER:MM[,...]`: pads closer than MM to their chamber's pad outline count as not `clean` (so measure, regroup.py, checks.py all drop them); `d_edge`, `fid` in pads.csv. **Biased - do not use for results** (kept as a tool; see findings 4) |
| 3bbcf0b | `fiducial_muons.py` | new: chamber-level efficiency A10/A20/A30 per muon with the fiducial cut on the muon's predicted point; per cut and per band of distance to the outline |
| (this file) | `HANDOFF_analysis_20261119c.md` | |

Pushing from the DAQ-analysis box: `origin` is SSH and no GitHub key is loaded there; use
`git -c credential.helper= -c credential.helper='!gh auth git-credential' push https://github.com/akallitss/P2_basket_analysis.git tick225-fiducial`.

**EOS code is plain copies, not git** (no `.git` anywhere under `$A`). Deployed by scp, checksum-checked against main
first: `$W/pipeline/cva/`, `$D/code/`, `$W/redump_test/cva/` have the new `coincidence_check.py` +
`channel_efficiency.py`; `fiducial_muons.py` only in `$W/pipeline/cva/`. `$A/eff_pipeline/cva/` is NOT updated
(stale, still 25 ns; not used this weekend). The pipeline scripts (`$W/pipeline/run_new.sh`, `env.sh`, `scripts/*.py`)
exist only on EOS.

`$W/pipeline/run_new.sh` (loop pass) changes, backups next to it:
- `run_new.sh.bak-20261119-prefid` = before any change today;
- `run_new.sh.bak-20261119-padfid` = with the pad-level fid block (writes `results/current/fid/fid{5,10,15,20}`);
- current: that block replaced by `fiducial_muons.py tags_ready pseudo $W/fid_muons --scripts $S` ->
  `results/current/fid_muons/` + `fid_muons.log`. Replaced atomically (mv) between passes.
  NB: a running bash reads its script while executing - edit run_new.sh only between passes (no `$W/state/lock.d`)
  and replace it with mv, never edit in place.

## What was run

| step | where | output |
|---|---|---|
| re-dump 22.5 ns, all chunks with the 19 Nov --diag fields | condor 13753921, 51 jobs: 18 Nov runs 231-267 -> `$D/tags2`, weekend 268-281 in place in `$W/tags` (weekend loop paused by holding `$W/state/lock.d` meanwhile; released 15:36) | 51/51 runs complete |
| `$D/tags_use` | all 2470 links now -> `../tags2/` (were `../tags/`, the 25 ns dumps, kept) | |
| old outputs kept | `$D/results/{all_thr,veto_thr,silent}_tick25`, `silent_tick25.log`; `$W/analysis/{silent,pulser_run270}_tick25`, `silent_tick25.log` | |
| re-analysis | condor 13754128 (`post_t225.sh` step arg): d_measure (= ~/diag_20261118/remeasure_thr.sh -> `$D/results/all_thr`, `veto_thr`), d_silent (`$D/results/silent`), w_silent (`$W/analysis/silent`), w_pulser 270 (`$W/analysis/pulser_run270`), fid D (`$D/results/fid/fid{0,5,10,15,20}`, pad-level, biased - **still running at 20:15**, informational only) | |
| weekend re-measure | loop pass 15:40 after removing `results/current/.inputs_sig`: runs 268-282 on the new dumps (+ the pad-level fid set); later passes add new runs | `$W/results/current/` |
| muon fiducial | condor 13754351 (`fidmu.sh D|W`) | `$D/results/fid_muons/`, `$W/analysis/fid_muons/` (+ `.log`); and every loop pass `$W/results/current/fid_muons/` |

Check of the fix: chunk 26001 (one 1-min file) wall span 66.67 s -> 60.00 s (ratio 0.9000); chunk 25510 (two files)
120.5 s. Tags per chunk +0.2-1 % (likely pairs across an SRS marker that were mis-timed before).

LaCie copy (user's disk, rsync `--no-o --no-g`, excl. tags/pseudo/redump_test):
`/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/LAB_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_weekend_20261119/`
(refreshed ~19:30; `analysis/silent_20mm/` = local redraw of the silent map in 20 mm squares with pads coloured by chip,
`redraw.py` next to it; it runs with the p2-mapping venv and `VMM_QA_DIR` = the EOS pipeline's `vmm_qa` copy, because
the local DAQ repo's `vmm_stations.py` has no `MAP_CSV`).

## Findings

1. **Tick fix: efficiencies unchanged** (as expected: tag/probe within one frame). 18 Nov per configuration `all_thr`
   vs `all_thr_tick25` identical to <= 0.002.
2. **Spark veto, now with a correct clock: still no effect.** 18 Nov: the veto removes 14.5 % of the
   scintillator-confirmed on-time tags per chamber and changes eps by <= 0.02 per configuration (within errors;
   `$D/results/veto_thr` vs `all_thr`). Sparks are not the loss.
3. **Pulser run_270 (PRa, set A: ch31 on chips 2,5,8,11,14,17), fixed lattice:** read-out 99.87 % on chips 2, 5, 11,
   14, 17 (`$W/analysis/pulser_run270.log`). Around the chamber-silent muons both lattice pulses were received
   ~100 % -> **silent muons are not read-out data loss.**
   **Chip 8 ch31 gets only 2.9 % of its pulses** (918 549 of 31.9 M). Not in the dead/masked list of run_270; in the
   config (`/local/p2/p2testbench/TestBenchCERN/config_ext/p2b-config-cern-ext_gain4.5_peaktime200_opt_cosmic_trimmed_Km18_v1_hvm_tpA.txt`
   on the DAQ box, chip N = Hyb N//2 Vmm N%2, chip 8 = Hyb 4 Vmm 0) it has the same `st 1`, `sdp10 40`, `sng 0` as the
   others, trim 25 (sd 19), global sdt 244 (all chips 226-268). Chip 8 is also the chip with eps 0.2-0.5.
   **Open:** a second set-A pulser run now exists in the weekend data (`gain9.0_gscan_hvon10_hvm_tpA`) - run
   `post_t225.sh w_pulser <run>` on it (find the run number in `$W/pseudo/chunks.json`) to see if chip 8 ch31 fails
   again; then the baseline / threshold of 8:31 (vmm_monitor_tool) and its hit rate in the box QA.
4. **The pad-level fiducial is biased: eps > 1.** Weekend `results/current/fid/`: P2_OUT 1.03-1.13 at 15-20 mm, P2_IN
   ~1.0; MID unchanged (not cut). Cause: the counting estimator (`channel_efficiency.py`) takes a pad's expected
   crossings from the density of the PREDICTED positions around it; the tracking smears them, so central pads are
   over-estimated and edge pads under-estimated, and only the whole chamber balances. **The same bias explains the
   per-chip eps > 1 of the central chips 3/11/12/15 (task 4)** - per-chip / per-region eps from this estimator are not
   absolute. Full-chamber eps is fine. Do not quote any fid/fidN numbers.
5. **Muon-level fiducial (`fiducial_muons.py`), the valid version.** A20 = P(any pad within 20 mm of the predicted
   point), muon window minus sideband, scint-confirmed tags within 6 mm of a clean live pad. Its absolute level is
   below the counting eps (containment: a 20 mm circle does not catch every muon's hit), so compare cuts, not levels.
   All configurations together:

   | cut (mm inside outline) | 18 Nov IN | MID | OUT | weekend 268-284 IN | MID | OUT |
   |---|---|---|---|---|---|---|
   | 0 | 0.710 | 0.779 | 0.726 | 0.738 | 0.787 | 0.764 |
   | 10 | 0.734 | 0.782 | 0.760 | 0.763 | 0.788 | 0.805 |
   | 20 | 0.740 | 0.782 | 0.765 | 0.769 | 0.787 | 0.815 |
   | 30 | 0.755 | 0.780 | 0.758 | 0.780 | 0.787 | 0.811 |

   Per band of distance (18 Nov / weekend): 0-5 mm IN 0.53/0.54, OUT 0.54/0.59, MID 0.71/0.73; 5-10 mm IN 0.64/0.62,
   OUT 0.65/0.66. **OUT** recovers fully by ~10-20 mm and is flat inside (+4-5 points from the cut). **MID** is flat
   everywhere (control, no edge effect). **IN keeps rising towards the centre** (18 Nov 0.69-0.73 at 10-40 mm, 0.76 at
   40-60, 0.80 at >= 60 mm): not only an edge band - a gradient across the chamber (HV / gas gap / mesh? or the
   chips 2-3 vs 5-6 difference). The cut moves every configuration by about the same amount, so the configuration
   ranking does not change. Same per-configuration A20 table in the logs (`fid_muons.log`, section "per configuration
   and cut").
6. **Silent muons (fixed clock):** 14-15 % per chamber on 18 Nov (79 574 muons) and on the weekend runs 268-281
   (30 851); edge band as before; P(silent | previous hit on the predicted chip < 1 us) still ~0.46 (97 / 222 muons in
   that bin), 0.19-0.27 at 1-10 us, ~0.12-0.18 beyond - a small dead-time-like effect (see also memory
   vmm-deadtime-sparks-study-nov19: dead time 0.3 %, sparks ~0.1 %, written by another session).
7. **Single-run spread (user's question in the morning, "how did D10 reach 0.95?"):** D10 P2_IN 0.951 +- 0.036 is one
   run, ~1.4 sigma above REF (0.884 +- 0.032); the table of 18 values scatters both ways (A50 IN 0.76). Muon-level A20
   (weekend, cut 0) puts the gain-3 files (D50/D10/D2) at or slightly above REF: IN 0.79-0.82 vs REF 0.78-0.79, OUT
   0.79-0.82 vs 0.81-0.82 - one run each; wait for the repeat passes P2/P3/P4 before calling it.

## Next steps (suggested order)

1. Pulser: `w_pulser` on the second set-A run (gain9.0_gscan_hvon10_hvm_tpA) and every set-B run (PRb/PCb) as they
   arrive -> is chip 8 ch31 always at ~3 %? Then look at 8:31 (baseline, threshold, trim) on the box (read-only).
2. P2_IN gradient: A20 x-y map with the muon-level estimator (fiducial_profile.csv is 1-D; extend `fiducial_muons.py`
   with a 20 mm x-y map, or reuse `silent_muons.py` map code with A20 instead of n_any == 0).
3. Per configuration over the repeat passes (runs.csv on the box gives the short names): D10/D50/D2 vs REF with both
   the counting eps (full chamber) and the muon-level A20 (cut 0 and 20 mm).
4. Remaining items of 20261119b: whole-hybrid silence from `chip_n1us`/`chip_n100us`; ADC threshold fit replacement
   (task 2); per-chip outliers (task 4: now explained for the central chips by finding 4, chip 8 still open);
   FINDINGS for the weekend (task 7).
5. Mon 23 Nov after 18:00: the lxplus git clone / VERSION file (user decision 3); merge `tick225-fiducial` into main
   with the user.

## Watch points

- lxplus ssh from the box goes through a ControlMaster opened by the user with password + OTP (Kerberos alone is not
  enough). If ssh fails with "Permission denied (publickey,keyboard-interactive)", ask the user to run
  `! ssh lxplus true`. Condor jobs and the loop keep running meanwhile.
- The loop runs every 30 min (cluster 13752682 on bigbird25, `condor_q -name bigbird25.cern.ch`); a pass re-measures
  only when the dumped runs or maps change. To force a re-measure: `rm $W/results/current/.inputs_sig`.
- Condor submit files: `arguments = $(step) $(arg)` - quotes around a queue variable are rejected.
