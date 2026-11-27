# Handoff — det5 at P2_MID, MID↔OUT read-out swap, new DAQ PC, DREAM vs stand (state 27 Nov 2026 12:10)

Read this first. Follows `HANDOFF_det5_night_20261124.md` (det5 install, Kapton, night 23→24) and
`HANDOFF_analysis_20261122.md` (funnel cuts, products, re-run recipe): both still valid where not superseded here.

Branch `note-dream-vmm-same-pads` of P2_basket_analysis (worktree `.claude/worktrees/note-dream-vmm-same-pads`),
based on `decoder-tick-225`; **not pushed** (carries the unpushed `decoder-tick-225` history — pushing / merging is the
user's call, e.g. via the end-of-day skill). Working tree clean at `7e9ae9d`.

## 0. Hardware and machines right now

| what | state |
|---|---|
| Chambers | P2_IN = det3, **P2_MID = det5** (CERN bulk, since run_389), P2_OUT = det1 |
| Read-out | **swapped since 24 Nov 12:14:51**: det5 (MID) read by H7-H9 = VMM 14-19, det1 (OUT) by H4-H6 = VMM 8-13 (adaptation cards moved with their cables, connector for connector); P2_IN unchanged (H1-H3, VMM 2-7) |
| Kapton | on the card now on det5 c4 (→ H8, VMM 16/17); it now cuts c4 top strip 97 (VMM 17 ch 31, 0 hits); VMM 17 ch 27/23 are live |
| HV / LV | **both OFF** since 25 Nov 07:42 (end of the night script) |
| DAQ PC | **new: dedippcy245** (`ssh vmm-DAQ-daplxa`). DAQ services running idle in tmux (flask, daq, processor, qa, backup watcher, rate guard, logbook…); last run **run_424**, next run_425 |
| Old DAQ PC | dedippce185 (`vmm_daplxa`) retired; full backups: EOS `vmm_cosmic_run_oct26/machine_backup_20261123/` + `machine_backup_20261124/` (MANIFEST in each) and the Transcend (other session, verified) |
| lxplus | swap-area loops ended 25 Nov; nothing running for this work |

### New-PC facts (details: DAQ project memory `new-daq-box-dedippcy245`, `status-daq-box-migration`)
- Paths changed: `/local/p2/data` (not p2data), `/local/p2/testbench/TestBenchCERN` (not p2testbench), data NIC `enp2s0`
  (capture files `enp2s0_*.pcapng`). Fix code, never symlink (user rule).
- Login on site needs the **CEA Kerberos ticket** (gssapi) — not the ssh key; if refused, ask the user for `kinit` or
  `ssh-add`. lxplus: `! ssh -fN lxplus` (2FA) when the master is gone.
- Start services from a KCM shell: `export KRB5CCNAME=KCM:; cd /local/p2/DAQ_Control_VMM_Beam; ./start_servers.sh`,
  then GUI `POST localhost:5002/start_processor`, `/start_qa`, `/start_backup`, and the rate guard
  (`bash_scripts/start_tmux.sh vmm_rate_guard ".venv/bin/python3 …/rate_guard.py --max-mb-min 100" 2000`).
- LV via GUI `POST /power/run {"action":"on"|"off"|"measure"}`. Needs pyvisa for SYSTEM python3: installed
  `--user --break-system-packages` as a stopgap; proper fix `sudo apt install python3-pyvisa python3-pyvisa-py`
  (user/IT). **After "on" (it power-cycles the TDKs) wait ≥ 2 min** before a run: the ALINX needs > 30 s to boot
  ("Fec does not respond to ping" / warm reset "1 failed hybrid" otherwise). Arduino step fails harmlessly (no
  /dev/ttyACM0, same on the old box).
- CERN ticket (KCM) renewable until **29 Nov 18:20** → re-kinit before then for the EOS backup.
- HV tool: `python3 /local/p2/DAQ_Control_VMM_Beam/vmm_monitor_tool/hv_lab.py state|ramp-up|ramp-down|guard|check-on`
  (logs `~/trim_work/out/hv_lab/`).

### Run scripts (copies in `~/Documents/PostDocSaclay/analysis_local/swap_mid_out_20261124/`)
- `night_20261124.sh` (new PC; on the box `~/trim_work/`): HV up + guard (trip → re-ramp once, then detector off),
  SEQ of 49-min runs until DEADLINE, then **HV down first, then LV off** (user rule), ALERT files on problems.
  `STOP` flag between runs; reusable for the next night: copy, set W/DEADLINE/SEQ/notes.
- `swap_20261124.sh` (old PC, LV power-cycle at start), `machine_backup_20261124.sh`.

## 1. Run list since 23 Nov

| runs | what | use |
|---|---|---|
| 389, 390 | det5 MID only HV, Kapton check | no |
| 391-406 | night 23→24, det5 at MID, default read-out (old PC) | yes (391-406 REF/D50/B50, 399 pulser) |
| 407 | swap, REF, **stopped by hand at 24 of 49 min** (colleague, PC change) — no `.subrun_complete` | yes (`EFF_FORCE_COMPLETE="407"`) |
| 408 | started as HV went off | **exclude** |
| 409 | new PC test, 5 min REF, **HV off** | **exclude** |
| 410-424 | night 24→25, swap, new PC: 11 REF, D50 413/422, B50 416, pulser 418; no trips | yes |

Corrections file: `run_metadata_corrections.json` (runs ≥ 389 det5) — **runs ≥ 407 swap not yet recorded there**
(the box setup entry `bdba478` and run_config.json of runs ≥ 407 are correct, so only add a note if needed).

## 2. Results

### 2a. Swap (funnel, current cuts, REF) — the read-out is cleared
| | before swap (391-406) | after swap (407, 410-424) |
|---|---|---|
| det5 @ MID silent | 17.3 ± 0.4 % (VMM 8-13) | **17.0 ± 0.4 %** (VMM 14-19) |
| det1 @ OUT silent | 8.0 ± 0.3 % (VMM 14-19) | **8.0 ± 0.3 %** (VMM 8-13) |
VMM 8-13 silent on det5 19/19/25/14/11/8 % → on det1 7/9/10/7/7/7 %. det2 at MID (weekend) was 13.1 %.
→ Not H4-H6 / VMM 8-13. Left: chamber vs MID position. **Only MID is on the iseg NHR HV supply** (IN/OUT on CAEN
N1471) → proposed next test.

### 2b. New PC / P2_IN (`in_timing_check.py`)
Muon hit timing unchanged (core median −3 → −8 ns, width 36 → 44 ns, window ±100). P2_IN latency +26 ns = swap
(tags now from other hybrids) — harmless. More random hits near predictions (IN ×3.8, MID ×3, OUT ×2.7): chip rates
changed — VMM 4 50 → 170 M/run, VMM 5 124 → 13 M (same H2), VMM 8/9 (now det1 c3) 1/2 → 13/18 M, VMM 15 (det5 c3)
4 → 38 M. They inflate IN late/far (mostly chip 4). **Open: VMM 2 silent 9 → 17 % at unchanged rate.**
The counting estimator (results/current chamber_level) is NOT usable on the swap area (noise-dependent tags; B50
collapses) — use the funnel.

### 2c. det5 DREAM bench vs stand (`det5_dream_vmm_pads.py`)
- **Bench fault:** DREAM **FEU 6 / DREAM connector 7, even channels dead** since Sep (det5 9 Sep on its c4 bottom:
  32 pads; det2 bench 24 Nov on its c5: 31/32; alive in July for det1-4). Bench det5 c4 = 75 % because of it.
  → tell the bench people / fix before the next DREAM run.
- Stand-dead before swap (VMM 9 ch 35, 11 ch 27, 11 ch 52): fine on DREAM, alive after swap → old read-out path /
  old Kapton, not det5 pads.
- c3 top strips 126/128 dead on DREAM and before the swap (VMM 9 ch 2 low, ch 0 zero), fire after the swap (VMM 15,
  noisy; edge pads, too few muons to test) → connector-end contact? c4 top strip 97 weak everywhere.
- The stand's silent excess is not in DREAM-weak pads; per-pad spread = statistics; the stand tests only ~160-220 of
  the 384 c3-c5 pads (≥ 15 muons, paddle footprint).
- det2 bench run 24 Nov (`p2_det2_alignment_mesh_drift_scan_long_run_430_11-24-26`): short (9.6k rays), alignment fit
  RMSE 54 mm → "within" 82 % unreliable; any hit 93.9 % (= July 94 %).

## 3. Code (branch, one commit per file)

| file | what |
|---|---|
| `stand_map.py` | chip → chamber per layout, env `VMM_STAND` (`default` / `swap_mid_out_20261124`); used by `coincidence_check.CH` (→ `efficiency_tp`), `channel_efficiency.DETV`, `probe_diagnostics.DETV`. ~20 other scripts still hard-code 8-13 = MID (only these three matter for the pipeline + funnel) |
| `muon_funnel_job.sh` | `W=` (area) and `S=` (scripts/paddles/pointing area) from env |
| `in_timing_check.py` | probe hit timing + latency per chip, two areas |
| `det5_dream_vmm_pads.py` | det5 pad by pad DREAM vs stand (inputs: stage-02 csv, `dream_pad_export.py` output, `silent_tests.py` JSONs night/swap, Mapper xy JSON) |
| `note/det5_swap_section.py` + hook in `note/muon_funnel_page.py` (`SWAP_HTML`) | note section 10 |

Not in the repo (EOS only): lxplus pipeline area `analysis_lxplus/swap_20261124/` = copy of the night pipeline with
the 4 new/changed cva files, `env.sh` `EFF_MIN_NEW=407`, `VMM_STAND=swap_mid_out_20261124`,
`EFF_FORCE_COMPLETE="407"` + a 3-line patch in its `scripts/find_complete_runs.py` (accept listed runs without
`.subrun_complete`). Loops `~/swap_20261124/loop_{dump,measure}.{sh,sub}` (END dates to edit, submit with
`-name bigbird25.cern.ch`). Funnel: `~/muon_funnel/muon_funnel_swap.sub` → `swap_20261124/analysis/muon_funnel_physical/`
(run408/409 moved to `excluded/`), `summary.json` built with `muon_funnel_summary.py`. silent_tests JSONs only in
lxplus `/tmp/ak_silent_tests_{night_20261123,swap_20261124}.json` (copies in the job tmp — regenerate if needed:
`silent_tests.py --tags tags --chunks pseudo/chunks.json --pointing <weekend>/analysis/muon_funnel_pointing/pointing.json --min-pad 15`
with `PYTHONPATH=<area>/pipeline/cva` and `VMM_STAND` for the swap area).

## 4. Products
- CERN note (unlisted) updated 25 Nov with **section 10** "det5 at P2_MID and the read-out swap" (what changed,
  Kapton, chamber table, per-connector bars, D50/B50, new-PC timing, conclusion), update box, open point 2 rewritten,
  open points 11 / references 12. cern-site `28f1453`, deployed and pushed.
  https://alexandra-kallitsopoulou.web.cern.ch/notes/unlisted/2026-11-19-where-do-the-missed-muons-go-the-vmm-cosmic-stand-step-by-step.html
  Rebuild section 10 in place: `note/det5_swap_section.py --weekend W --night N --swap S --timing T --out frag --note <html>`
  (idempotent `<!-- swap:* -->` markers), then cern-site `python3 scripts/build_hub.py`, commit, `scripts/deploy.sh`, push.
- The DREAM-vs-stand det5 findings (2c) are NOT in the note yet (offered as 10g).

## 5. Open points, in order
1. **MID position vs chamber**: swap the HV channels (MID on a CAEN N1471 channel, OUT on the NHR) — no mechanics,
   2-3 REF runs (~2.5 h incl. ramp) answer it; else a physical det5 ↔ det1 swap. User decides the HV change.
   Remember: a new layout needs a setup entry on the box (`p2_mapping.edit.add_entry`, as for `bdba478`) if the
   read-out moves; an HV-channel swap needs the setup's `hv.channels` changed (hv_lab and the DAQ read targets from it).
2. **Bench FEU 6 / connector 7 even row** dead — report / fix (cable/adapter/FEU).
3. **VMM 2 silent 9 → 17 %** after the PC change; and the large per-chip rate changes (VMM 4/5, 8/9, 15) — check
   noise scans on the new PC vs the old (HV-off sigma scan like `fastscan`, needs LV on, HV off).
4. Note 10g with the DREAM-vs-stand det5 results (2c), if the user wants it.
5. `run_metadata_corrections.json`: add the runs ≥ 407 swap / 408-409 exclusion entry.
6. Analysis scripts hard-coding `enp4s0f1` (vmm_read_pcapng.py, recover_root_files.py, cable_length_study.py) and the
   ~20 scripts with fixed 8-13 = MID ranges: switch to `stand_map` / interface-agnostic when touched.
7. Proper `apt install python3-pyvisa python3-pyvisa-py` on the new PC (sudo); CERN KCM re-kinit before 29 Nov.
8. Older open points of HANDOFF_20261122 §4 (efficiency definition on hold, pillar term, paddle refit, crosstalk…).
9. Merge/push `decoder-tick-225` + this branch with the user (end-of-day skill).
