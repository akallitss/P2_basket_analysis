# Handoff — det5 at P2_MID, Kapton check, night runs 23→24 Nov 2026 (state 24 Nov 10:40)

Read this first. Branch `note-dream-vmm-same-pads` of P2_basket_analysis (worktree
`.claude/worktrees/note-dream-vmm-same-pads`), based on `decoder-tick-225`; **not pushed** (it carries the unpushed
`decoder-tick-225` history — pushing is the user's call: `git push cea note-dream-vmm-same-pads`).

## 0. Hardware state right now

| what | state |
|---|---|
| Chambers | P2_IN = det3, **P2_MID = det5** (CERN bulk, since run_389), P2_OUT = det1 |
| HV | **OFF** (guard ramped down 08:07, all channels < 4 V) |
| LV | **OFF** (end script, 08:08, all three TDK supplies 0 V / 0 A) |
| Before the next run | LV on (GUI power "on"), then a normal run start re-applies the chip config (VMMs lose it at power off); HV up with `hv_lab.py ramp-up`; `config/lab_hv_hands_off` is still set (runs never set HV) |
| Pruner / backup | running, untouched: tmux `vmm_prune` (`prune_backed_up.py --low 70 --target 90`), `vmm_backup_watcher` (EOS copy). Disk 78 GB free |
| Box Kerberos | renewed by the night script; renewable until 26 Nov 14:13 — nothing renews it now |
| lxplus condor (bigbird25) | night loops `13770497` dump (until 24 Nov 12:00), `13770498` measure (until 14:00); `condor_q -name bigbird25.cern.ch` |

## 1. What changed on 23 Nov

- **det5 replaced det2 at P2_MID** (user). Same connectors c3–c5, hybrids H4–H6 = VMM 8–13, same HV cabling, det5 at
  430/730 V like det2. det5 has the CERN insulation mask (0.5 mm pillars on 2 mm, **4.87 %** dead area; det1–4 have
  Saclay mask V1, 3.78 %). Recorded in `cosmic_vmm_analysis/run_metadata_corrections.json` (runs ≥ 389) and the lab
  logbook. **`run_config.json` of runs ≥ 389 still says MID = det2** (DAQ detector table not edited).
- **Kapton on MID c4** (user): cable to the bottom connector of the adaptation card, side B, ~3 contacts "after pin 38"
  taped; the photo suggests the tape may cover ~pins 37–46. Photos:
  `~/Documents/PostDocSaclay/analysis_local/det5_kapton_20261123/*.jpg`. **Not masked in the chip configs** (user).
  Card pins 37–46 = **VMM 10 and VMM 11 channels 30, 31, 28, 29, 26, 27, 24, 25, 22, 23** (same channel numbers on
  both halves; VMM 10 = bot = strips 1–64, VMM 11 = top = strips 65–128; from `P2_BASKET_mapping.csv` mec8_pin and
  `vmm_stations.build_pad_table`, SPS c5 → lab c4 = connector_number 3).

## 2. Kapton check: what was tried, what is left

1. run_389 (REF, 5 min) and run_390 (bscan_hvon50, 3 min), det5 HV only: per-channel rates 0.3 / 1.2 Hz — too quiet
   to separate channels. Odd but not conclusive: VMM 11 ch 26/27 (pins 41/42) and VMM 10 ch 23 (pin 46) at ~20 Hz.
2. HV-off noise σ scan of chips 10–11 (23 Nov 17:17, `~/trim_work/out/det5_kapton_20261123/sigscan`) and of all
   chips (night, `out/night_20261123/fullscan`) vs the 16 Nov `fastscan_20261116/full2` (det2): **no channel on pins
   37–46 drops** (0.86–0.97 × 16 Nov, chip medians 0.91–0.95 on MID). The only clear drop: **VMM 10 ch 43**
   (card pin 24, strip 21, channel_id 404): σ 2.85 → 1.45–1.69 mV, noise rate 38 Hz → 0 — in both scans.
   **Caveat:** σ is dominated by the 2 m Samtec cable, so a pad cut at the DETECTOR-side card barely changes σ. The
   σ test only excludes a cut between cable and hybrid.
3. **Left to do (decisive):** muon signals with HV on. From the night REF runs 391–406: per channel of VMM 10/11 the
   high-ADC hit count (QA `hits_store/*/counts.npz` → `adc_vs_ch` (32, 64, 128) on the box, small files; or the
   per-pad found share of `results/current/pads.csv` on lxplus). Disconnected pads: no muon Landau tail / found ≈ 0.
   Check pins 37–46 AND VMM 10 ch 43.
   Scripts: `analysis_local/det5_kapton_20261123/` (`kapton.py` per-channel rates from hits_store, `sigcmp*.py` σ
   comparison, outputs `sigcmp_out.txt`, `sigcmp_all_out.txt`).

## 3. Night 23→24 Nov (all done, no failures)

Supervisor `~/trim_work/night_20261123.sh` (log `out/night_20261123/supervisor.log`, `runs.csv` there), end script
`~/trim_work/end_20261124.sh` (log `end.log`). Sequence: HV-off σ scan of chips 2–19 (17:32–17:51) → ramp-up
(IN 420/720, MID 430/730, OUT 430/730; currents ≤ 0.15 µA, no trips all night) → 16 runs of 49 min:

| runs | config |
|---|---|
| 391, 392, 393, 395, 396, 398, 400, 401, 402, 404, 405 | REF `trimmed_Km18_v1_hvm` (11 runs, ~9 h) |
| 394, 403 | D50 `gain3.0 gscan_hvon50` |
| 397, 406 | B50 `gain6.0 gscan_hvon50` |
| 399 | PRa `trimmed_Km18_v1_hvm_tpA` (REF + test pulses on ch 31 of chips 2,5,8,11,14,17) |

Plus 389 (REF 5 min) and 390 (bscan_hvon50 3 min) from 23 Nov afternoon with only MID HV on (ignore for efficiency).
Deadline 08:00 → last run ended 07:24; guard ramp-down 08:07; LV off 08:08.

## 4. First results (lxplus pipeline, counting estimator, 24 Nov 08:25)

EOS `…/vmm_cosmic_run_oct26/analysis_lxplus/night_20261123/results/current/` (runs 389–406) vs the det2 weekend
`…/weekend_20261119/results/current/` (runs 268–388, final). `eps` at the default radius:

| config | P2_IN det3 night / weekend | **P2_MID det5 night / det2 weekend** | P2_OUT det1 night / weekend |
|---|---|---|---|
| REF | 0.884 / 0.887 | **0.842 ± 0.016 / 0.850 ± 0.010** | 0.919 / 0.913 |
| D50 | 0.921 / 0.926 | 0.860 / 0.865 | 0.916 / 0.898 |
| B50 | 0.803 / 0.785 | 0.850 / 0.846 | 0.889 / 0.875 |

- IN and OUT reproduce the weekend → stand, trigger and electronics stable.
- **det5 at MID is as efficient as det2 was** (within errors), although det5 has ~1.1 pt more pillar area. So the
  P2_MID deficit seen all along (note section 5d/9: "MID silent ~13 %, uniform, not a weak det2 area") does **not**
  follow the chamber: it stays with the MID position / its electronics (H4–H6, VMM 8–13) / the stack geometry. This
  is the chamber-swap test of the note's open point 2, done with det5 instead of a MID↔OUT swap.
- Not yet run on the night data: muon funnel (`muon_funnel.py`, condor `.sub` in `$W/analysis/muon_funnel*`),
  silent tests, per-pad maps. Do these next to confirm with the funnel's FOUND/SILENT numbers (det2 REF: MID found
  79.7 %, silent 13.5 %).

Electronics reference (HV-off σ, night vs 16 Nov, per-chip median ratio): IN chips 2–7 0.99–1.01, OUT chips 14–19
0.97–1.00, **MID chips 8–13 0.91–0.95** (det5 quieter than det2 at the inputs). A few channels of chip 6 have no new fit
(ch 1, 23, 31; chip 6 had only 48 fitted channels) — chip 6 is the noisy P2_IN chip of §6.

## 5. Other results of this session (22–23 Nov)

- **Missed-muons note, section 9 "The same pads with DREAM"** — deployed (unlisted)
  https://alexandra-kallitsopoulou.web.cern.ch/notes/unlisted/2026-11-19-where-do-the-missed-muons-go-the-vmm-cosmic-stand-step-by-step.html
  (cern-site commit c53558c, not pushed to GitHub). Code: `cosmic_vmm_analysis/dream_pad_export.py`,
  `note/dream_same_pads_section.py`, hook in `note/muon_funnel_page.py` (`SAME_PADS_HTML`). Inputs:
  LaCie `data/LAB_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_same_pads_20261122/`.
  Findings: DREAM bench on the stand pads c3–c5: det1 98 % / det2 94 % / det3 89 % any hit; VMM gap uniform on det2,
  smaller uniform on det1, det3 loss = one chip (VMM 6, c5). Beam det3 on 24 common pads: DREAM 97.3, VMM 89.0,
  bench 92.4, stand 78.4 %.
- **Chip 6 autopsy (not in the note yet):** VMM 6 (hybrid 3 bot, P2_IN c5): DAC settings ordinary (sdt 241, trims
  12–29, none masked) but 27–29 mV above baseline on all 57 channels vs ~16 mV on chips 3–5 (REF 6.2 vs 3.6 fC); the
  noisiest chip in every noise scan (gain 3 R10: 13.8 vs 9–11 mV); loss uniform over channels; detected muons normal in
  ADC and timing; even at equal nominal fC ~2× more silent; its silent muons are independent of chip activity (unlike
  other chips). Never pulsed in the pulser runs. Proposed: short REF run with chip 6 sdt −10…−15 DAC; swap the c5
  cable / hybrid 3. Scripts + outputs: LaCie `…/analysis_same_pads_20261122/chip6/`.
- Side observation, not followed up: on chips other than 6, silent muons coincide with the predicted chip being
  empty for ±100 µs far more than chance (60–90 % vs 14–39 %).

## 6. Open points, in order

1. Kapton: identify the cut channels from the night muon data (§2.3); tell the user which pads of det5 are blind.
2. Funnel + silent tests on runs 391–406 (REF), side by side with the det2 weekend; update the note (new section:
   "det5 at P2_MID") — the MID deficit does not follow the chamber.
3. Chip 6 tests (§5) — needs the user's go for data taking.
4. DAQ `run_config` detector table: MID = det2 → det5 (user's call; corrections file covers the analysis).
5. lxplus: the night loops end today 12:00 / 14:00; for new runs either extend (`~/night_20261123/loop_*.sh`, END
   lines) or make a new area the same way (copy `pipeline/`, set `EFF_ANA` and `EFF_MIN_NEW` in `pipeline/env.sh`).
6. After the analysis: merge `decoder-tick-225` + this branch with the user (end-of-day skill).

## 7. Access notes

- Local ssh agent lost its key on 23 Nov ~17:30 (`ssh-add -l` empty): box via `vmm_daplxa` needs
  `! ssh-add ~/.ssh/id_ed25519`; lxplus needs `! ssh -fN lxplus` (2FA). Analysis on lxplus, not on the box (user rule).
- `fast_ratescan.py` refuses HV > 5 V; `baseline_scan.py --hv on` needs all six HV channels at target.
- LaCie (`~/Documents/PostDocSaclay/data` symlink) was unmounted on 23 Nov evening; 23 Nov Kapton files are on the
  internal disk `~/Documents/PostDocSaclay/analysis_local/det5_kapton_20261123/`.
