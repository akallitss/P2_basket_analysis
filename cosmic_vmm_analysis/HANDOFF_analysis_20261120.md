# Handoff: VMM cosmic efficiency - where the inefficiency comes from (20 Nov 2026, ~16:00)

Single entry point for the next session. It follows `HANDOFF_analysis_20261119d.md` (still valid for the rules, the
read-out / dead-time / spark results and the EOS layout) and supersedes it where they differ.

## 0. Rules (unchanged, from the user and memory)

- All analysis on **lxplus** (HTCondor for anything that decodes or loops over the tag dumps). The DAQ box `vmm_daplxa`
  is read-only for us: logs, logbook, small report files; never compute there. Never touch HV or the DAQ.
- Memory-cap every process: `systemd-run --user --scope -q -p MemoryMax=3G bash -c 'source pipeline/env.sh; ...'`.
- EOS: `mkdir -p` of deep paths fails, mkdir one level at a time; `cp -p` fails.
- **Run numbers and configurations: read the logbook, never guess** (user, 20 Nov):
  `/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/logbook/LOGBOOK.md` (current copy of the box's).
- Pulser runs (`_tpA`, `_tpB`) are excluded **by configuration name**, never by run number.
- Commits: one per file; the user's memory says no Co-Authored-By line. GitHub push fails (ssh key refused) since 19 Nov:
  the branch exists only locally until the user fixes access.
- **The efficiency definition is ON HOLD** (user): do not adopt one; the funnel is a diagnostic.
- `ssh lxplus` needs the user's ControlMaster; on "Permission denied" ask for `! ssh lxplus true`.

## 1. The question and the answer so far

VMM-read P2 chambers: ~80 % found (P2_MID, REF) on the cosmic stand; DREAM 93-94 % on the same pads. Is it the
counting, and can the VMM settings (gain, thresholds) close the gap?

1. **Counting is not the problem.** Read-out 0.14, dead time 0.29, sparks 0.12 points (19 Nov). Following every
   scintillator coincidence (muon funnel), P2_MID records **no hit at all** for 13-14 % of the muons that certainly
   crossed it (REF: 13.8 % of 6 068). Tag-line back-extrapolation to the (fitted) paddles: removing the 13 % of lines
   that miss moves P2_MID 79.6 -> 80.6 %, silent unchanged.
2. **More VMM gain does not help, and nothing is mixed up** (20 Nov tests):
   - the gain is real: muon pulse height / gain-3 = 1.8-2.2 (gain 6), 2.8-3.2 (gain 9) (`muon_adc_gain.py`);
   - effect 1: the same channel needs 1.19x (gain 6) / 1.47x (gain 9) more charge for the same noise rate
     (`trim_reports.py` on the bscan_target reports) -> threshold/signal gets worse with gain;
   - effect 2: per-chip spread of the needed thresholds 10 / 20 / 39 mV (gain 3/6/9) vs trim range ~34 mV -> 360-708
     channels clipped at gain 9; **but clipped channels are not more silent** (`silent_vs_clipped.py`): clipping is
     not the cause of the gain-9 loss;
   - P2_MID silent vs threshold in charge: gains 3/4.5/6 at 3.2-4.4 fC -> 12.4-14.5 %; gain 9 at 6.7-8 fC -> 22-25 %;
     slope 2.9 +- 0.2 points/fC overall, 1.3 +- 0.8 inside the reachable range; best ~87-88 % "any hit".
3. **What the silent muons are** (P2_MID): only 3-5 points look like a threshold tail (spectrum of found muons extended
   below the edge, `threshold_spectrum.py`; edge at 12-19 % of the median muon signal); charge sharing is not the cause
   (93-97 % single-pad); at gain 9 ~4 points are LATE hits (100 ns-1 us after the window, `silent_late_hits.py`);
   **8-10 points (16 at gain 9) leave no signal above threshold at all and are not explained.**
4. Crosstalk: rejected tag groups (6-11 % of active muon windows) are neighbour-channel crosstalk (single pad, small ADC,
   in time, 60-65 % 1-3 channels from a muon channel; consecutive channels of a chip are never adjacent pads, median
   149 mm). Rate vs muon pad ADC: ~1 % for medium signals, 14-29 % (+-1) and 27-35 % (+-2!) for saturated pads.

User's reading (20 Nov): "if the electronics cannot do better than ~90 %, we have to push the detector gain".
Our nuance: detector gain recovers the threshold-tail part (3-5 points) for sure; the no-signal part (8-10 points) is
unknown - it may not respond to gain at all. **The mesh-HV scan with REF electronics is the direct test (user decides).**

## 2. Code (branch `decoder-tick-225`, worktree `.claude/worktrees/decoder-tick-225`, NOT pushed, not merged)

All in `cosmic_vmm_analysis/`, run from `$W` (below) with `source pipeline/env.sh` and
`PYTHONPATH=$W/pipeline/cva:$PYTHONPATH` (the 22.5 ns coincidence_check etc.).

| script | what | output on EOS |
|---|---|---|
| `muon_funnel.py` (+ `_summary.py`, `_job.sh`, `.sub`) | every scintillator coincidence -> one class (no tag / tag rejected / outside / dead / silent / far / late / found), per chip, rejected-group pieces, examples, `--paddles` check, crosstalk rate | `$W/analysis/muon_funnel/run*.json`, `summary.json` |
| `gain_scan_eff.py` | funnel classes per gain-scan setting and chip vs threshold in fC; fits | `$W/analysis/gain_scan/out/gain_scan_eff.*`, `gain_scan_chip.csv` |
| `muon_adc_gain.py` | real gain from muon pad ADC; ADC edge | `.../out/muon_adc_gain.csv` |
| `trim_reports.py` | calibration reports: no-crossing channels, clipped, targets in fC, per-chip spread, same-channel ratios | stdout (local reports) |
| `threshold_spectrum.py` | found-muon max-pad spectrum (DNL corrected), below-edge extrapolation, cluster size | `.../out/threshold_spectrum.json` |
| `silent_vs_clipped.py`, `silent_late_hits.py`, `silent_vs_pad_position.py` | silent-muon tests | `.../out/*.txt` |
| `tag_rejects.py` (+ `_summary`, job) | rejected tag groups (19 Nov) | `$W/analysis/tag_rejects/` |
| `wrapup_losses.py` | read-out / dead time / sparks figures + definitions table (19 Nov) | `$W/analysis/wrapup_20261119/out/` |
| `note/muon_funnel_page.py` | builds the CERN-site note | (cern-site) |
| `note/gain_deck_figs.py` | deck figures (lxplus) | `$W/analysis/gain_scan/deck_figs/` |
| `note/deck_update_20261120.py` | inserts the new section into the 18 Nov deck (uses build.py helpers) | LaCie deck |

`$W = /eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus/weekend_20261119`.
Threshold reports copied to `$W/analysis/gain_scan/reports/` from the box
(`~/trim_work/out/gain_20261118/hvon_g*_R*_report.json`, `~/trim_work/out/bscan_20261116/hvon*_v1_report.json`).
Paddle fit: `$W/analysis/muon_funnel/paddles_fit_20261104.json` (copy of the 4 Nov trigger-acceptance `fit.json`).

**Re-running for new runs** (the weekend continues until Mon 23 Nov 07:00; the lxplus loop dumps new runs every 30 min):
```
scp cosmic_vmm_analysis/muon_funnel_job.sh lxplus:muon_funnel/      # (already there)
ssh lxplus 'cd muon_funnel && condor_submit muon_funnel.sub'          # skips done / partly dumped runs
cd $W; python3 analysis/muon_funnel/muon_funnel_summary.py analysis/muon_funnel analysis/muon_funnel/summary.json
python3 analysis/gain_scan/gain_scan_eff.py analysis/muon_funnel/summary.json results/current analysis/gain_scan/out
python3 analysis/gain_scan/{muon_adc_gain,threshold_spectrum}.py tags pseudo/chunks.json analysis/gain_scan/out/...
python3 analysis/gain_scan/gain_deck_figs.py analysis/gain_scan/out analysis/gain_scan/reports \
        analysis/muon_funnel/summary.json analysis/gain_scan/deck_figs
```
Then locally: the page (`python3 note/muon_funnel_page.py SUMMARY NOTE_HTML GAIN_SCAN_JSON`, copy the 7 deck figures to
`notes/unlisted/<slug>_files/`) and the deck (`<deck venv python> note/deck_update_20261120.py BUILD_PY OLD_PPTX
FIGDIR SUMMARY GAIN_SCAN_JSON OUT_PPTX`; venv with python-pptx: `~/.claude/jobs/6a168f88/tmp/venv` - copy it if that
job dir is cleaned; template inputs are the 18 Nov pptx and `.build_P2_Detector_Performance/build.py`).

## 3. Products

- CERN-site note (unlisted, behind CERN SSO):
  https://alexandra-kallitsopoulou.web.cern.ch/notes/unlisted/2026-11-19-where-do-the-missed-muons-go-the-vmm-cosmic-stand-step-by-step.html
  (sections: set-up and words, funnel, efficiency step, event displays per class, crosstalk, paddles, gain scan,
  four tests, open points, references). cern-site repo: committed, deployed with `scripts/deploy.sh`; not pushed.
- Deck: LaCie `.../1.General_Info/presentations/WeeklyMeetings/P2_Detector_Performance_Gain_2026-11-20.pptx` + `.pdf`
  (70 slides: the 18 Nov deck + new section after its summary, 24 slides incl. three step-by-step builds, + 4 backups).
  Build inputs: `.build_P2_Detector_Performance/update_20261120/`.
- `WRAPUP_20261119.md` (19 Nov tests and the efficiency-definitions table).

## 4. Open points (suggested order)

1. **The no-signal muons** (8-10 points of P2_MID): map them (x-y, chip, pad, time of day, run), compare with DREAM on
   the same pads (does DREAM lose the same places? M3 bench `ray_hit_miss_list*.csv`), check pillars / known dead
   regions, check whether they cluster in time (HV / spark state). The silent-muon dumps (`n_any == 0`) have px, py,
   vmm, ch, wall time.
2. **Mesh-HV scan** proposal for the user: P2_MID +10/+20 V with REF electronics, silent vs HV (user decides HV).
3. **Gain-9 late hits**: time of the late hit vs its ADC (saturated pulses?); is the VMM time stamp of large pulses
   delayed? `fired_nr_dt`, `fired_nr_adc` in the dumps.
4. **Calibration gap**: 6-38 % of the channels had no measured noise-rate crossing in the bscan_target scans and got the
   chip median; at gain 6 R=2 and gain 9 R=10 their muons are more often silent (28-33 %).
5. **bscan settings**: 18-45 % FAR with chance hits < 1 % -> tag positions wrong in those settings; unexplained.
6. Crosstalk: filter single neighbour-channel pads before the tag cluster cuts (+15 % tags); why +-2 > +-1.
7. Paddle positions: refit on the weekend data or measure; then the paddle check can become a cut.
8. After Mon 23 Nov 18:00: merge `decoder-tick-225` into main with the user (HANDOFF d step 8).

## 5. Watch points

- `n_any` (dumps, +-100 ns around tag + latency, live hits) and the funnel's silent (any hit 140-310 ns after the
  scintillator) differ slightly; the "parts" figure mixes them (stated in the notes).
- `threshold_spectrum.py` edge = 1st percentile of the found muons' max-pad ADC; baseline fixed at 40 ADC (per-channel
  baselines vary); the below-edge extrapolations are smooth assumptions - a steep rise towards zero (few-electron
  events) would hide more.
- Funnel REF numbers change as more REF runs arrive (now 268, 269, 284, 299, 301, 303, 305); rebuild page and deck
  from the same summary.json.
- The funnel `run285.json` (tpB) of the first pass is kept in `$W/analysis/muon_funnel/v1_20261119/` - not used.
