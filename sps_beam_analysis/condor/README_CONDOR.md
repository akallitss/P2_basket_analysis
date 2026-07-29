# Running the SPS P2 analysis on lxplus HTCondor, from the EOS backup

**Written 2026-07-28.** Why this exists: banco's disk filled up and the run data
moved to EOS. The analysis stages read a local directory tree, and a serial pass
over the whole campaign is measured in days. This runs one job per sub_run on
HTCondor, straight from the EOS copy, and brings the (small) products back to
the directory the DAQ GUI's Analysis tab already reads.

**The GUI keeps working, unchanged.** `flask_app/app.py` serves
`…/TB_July2026_H4/analysis` from local disk. `merge_and_pull.sh` fills that same
directory. No Flask code, config or restart is involved.

---

## 1. The three commands

```bash
cd ~/P2_basket_analysis/sps_beam_analysis/condor

./submit.sh --group rec  --source ntof,salsachip   # FIRST: seconds per job
./submit.sh --group hits --source ntof,salsachip   # 257 jobs, ~15 min wall
ssh -o GSSAPIDelegateCredentials=yes akallits@lxplus.cern.ch 'condor_q'
./merge_and_pull.sh                                # products home + curves
```

`--group wave` runs stage 29 (waveform timing) instead — same sub_runs, but each
reads `decoded_root` and takes ~18 min, i.e. ~50 h serial compressed into one
batch. `--run NAME[,NAME…]` restricts to particular runs; `--dry-run` builds the
joblist and stops.

### The three EOS sources

The campaign moved destination twice in three days while quotas were being
sorted out, so it is spread across three locations and **the early runs exist
only on salsachip**:

| `--source` | location | holds |
|---|---|---|
| `ntof` | `eospublic:/eos/experiment/ntof/data/x17/p2_sps_july` | 169 sub_runs, everything from 07-26 on |
| `salsachip` | `eosproject:/eos/project/s/salsachip/…/P2_SPS_Dream_Data` | 88 sub_runs — `beam_commissioning_1`, `beam_nominal_meshscan_1`, `drift_mesh_scan_1`, `drift_scan_1/2/final`, `env_test_1`, `highstat_eff_1`, `latency_scan_1`, `meshscan_fine_1`, `p2in_check_1`, `run_1` |
| `user` | `eosuser:/eos/user/a/akallits/P2_SPS_backup_temp` | only `drift_mesh_2d_2` + `eff_nominal_1`, both already on nTOF — nothing unique |

`--source ntof,salsachip` gives **257 jobs**, the whole campaign. Order matters:
the first source listed that actually holds a run's data wins it, so runs
present in both (`drift_mesh_2d_1`, `low_mesh_scan_1`, `mesh_drift_scan_up_1`)
come from nTOF.

**Reads may come from any source; writes ALWAYS go to nTOF.** salsachip has been
over quota since 2026-07-25 (`[3021]` on every write) and the user space was a
stopgap, so both are read-only here. `job.sh` keeps `EOS_URL/EOS_BASE` (input,
per-source) separate from `OUT_URL/OUT_BASE` (output, always nTOF); the `rec`
group therefore writes `recorded_events.npz` into the nTOF tree even for a
salsachip run, and the `hits` group looks there first.

## 2. What runs where

```
banco                     lxplus (AFS)              worker node          EOS
-----                     ------------              -----------          ---
submit.sh
 ├ tar the code  ───────> p2code.tgz  ───────────>  unpacked to scratch
 ├ make_joblist.py  <──────────────────────────────────────────────  xrdfs ls -R
 │   (169 lines)   ─────> joblist.txt
 └ condor_submit  ──────> analysis.sub ─┐
                                        └────────>  job.sh
                                                     ├ stage inputs   <── runs/
                                                     ├ run stages
                                                     └ push products  ──> analysis/
merge_and_pull.sh  <───────────────────────────────────────────────────── analysis/
 ├ --scan-only merges (local, reads only scan_row.json)
 └ push merged curves ─────────────────────────────────────────────────> analysis/
```

`job.sh` rebuilds `<DATA_ROOT>/<run>/<sub_run>/…` in the worker's scratch dir
exactly as the DAQ writes it, then points `SPS_DATA_ROOT` at it. **The stages
have no idea EOS is involved** — no code in them needed changing for this.

## 3. Stage groups

| group | stages | input | measured |
|---|---|---|---|
| `probe` | 21, 28 | `combined_hits_root` | 205 s, shakedown only |
| `rec` | `extract_recorded_events` | `decoded_root` **read remotely** | seconds |
| `hits` | 21 → 22 → 20 → 23 → 26 → 28 | `combined_hits_root` (0.4–2.1 GB) | **6–19 min, up to 4.7 GB peak** |
| `wave` | 29 | `decoded_root` (~370 MB) | ~18 min (banco figure) |

**Run `rec` before `hits`.** Stage 22 needs `<sub_run>/recorded_events.npz` to
correct for DAQ overlap: with zero suppression, a trigger the probe FEU never
recorded looks exactly like one where it recorded nothing, so without the npz
the efficiency is biased low (it is the difference between `eff` and
`eff_corr`). `rec` writes that file into the nTOF tree (never beside the data, which for
salsachip is unwritable), and a later `hits` job stages it automatically — it
looks in the nTOF tree first, then falls back to beside-the-data, which is how
the runs that already carry an npz from banco keep working.

Measured on three sub_runs so far, DAQ loss came out at **0.0 %** every time, so
`eff_corr == eff`. That looks like a genuine null for these runs rather than a
missing correction — but it is measured now, not assumed, and runs whose FEUs
dropped triggers (decoder hangs) are where a non-zero value should show up.

`rec` is nearly free because it needs only the `eventId` branch: LCG's uproot
reads that one column over xrootd (measured **2.4 M eventIds in 0.5 s**), so
the group stages nothing at all rather than moving 370 MB per sub_run — 63 GB
of transfer avoided across the campaign.

21 must precede 22 (it writes the alignment 22 consumes) — they are in the same
job, so no DAGMan is needed. Across sub_runs there are no dependencies at all.

## 4. The one thing that needed new code: scan-level products

Stages 20, 22, 26 and 28 build a `rows` list across every sub_run and then draw
a scan-level curve from it. In one process that is free. Split across 169
machines, `rows` never exists anywhere.

So each stage now **persists its row** as `scan_row.json` beside that sub_run's
products, and takes `--scan-only`, which skips the expensive pass entirely and
rebuilds the curve from those files — a few kB each. `merge_and_pull.sh` runs
that on banco after pulling.

Verified 2026-07-28 on `eff_drift_ab_1` stage 26: 12 independent per-sub_run
runs plus `--scan-only` reproduce the serial whole-run CSV **byte-for-byte**,
both stations.

Two run-wide properties had to stop being inferred from "the sub_runs I can
see", because a worker sees exactly one:

- **the scan axis** — `sps_config.scan_axis()` now falls back to
  `run_config.json` when given fewer than two sub_runs. Without this, every
  point of a *drift* scan was labelled with the (constant) *mesh* voltage.
  This was a real bug, caught by the byte-for-byte comparison above.
- **the product directory** — `make_joblist.py` emits a `prod_sub` column
  (`scan` for a multi-point run), which `job.sh` passes to stages 20/22 as
  `--prod-sub`, so a scan's per-point eff maps file together under `scan/`
  where a serial pass would have put them (and where the GIF builder looks).

Stage 20's channel-space mode has no scan-level aggregation and is skipped by
`--scan-only`; it is per-sub_run only by design.

## 5. Things that will bite

- **`GSSAPIDelegateCredentials=yes` on every ssh/scp.** Without it you reach
  lxplus but hold no AFS token, and cannot read your own `.bashrc`. The error
  looks nothing like a Kerberos problem. All scripts here set it.
- **AFS home is ~88 % full (4.6 of 5 GB).** Nothing but scripts and logs may
  live there. `products.tgz` comes back *empty* whenever the EOS upload
  succeeded, precisely so a 169-job sweep does not deposit ~340 MB on AFS.
- **`set -u` breaks the LCG `setup.sh`** (it reads `COMPILER` without setting
  it). `job.sh` drops the flag across that one source and restores it.
- **A missing `transfer_output_files` puts the job on hold** instead of
  returning its logs — the worst outcome for a job that died early. `job.sh`
  therefore creates `products.tgz` in its first second and repacks it from an
  `EXIT` trap, however it exits.
- **`xrdcp -r` fails on an empty directory.** `out_dir()` creates a tree for
  every detector before knowing whether it has products, so the uRWELL
  references (no pad map) leave empty trees that failed the recursive upload
  and dragged the whole job into the per-file fallback. `job.sh` now skips
  dirs with no files. Both `job.sh` and `merge_and_pull.sh` prefer one
  recursive `xrdcp -r` per station dir (3–4 round trips instead of ~50); the
  per-file loop remains in `job.sh` as a fallback and measured 50 files in 6 s,
  so it is a slower path, not a disastrous one.
- **Retries are for infrastructure only.** `job.sh` exits 1 for a staging/setup
  failure (worth retrying) and 2 when a stage failed (deterministic —
  `retry_until` in `analysis.sub` stops those). Triage a sweep with:
  ```bash
  ssh … 'grep -h "STAGE FAILURES" ~/p2_condor/logs/*.out'
  ```
- **Sub_runs vary ~5x in size, and memory tracks input at ~2.2x.** A
  `low_mesh_scan_1` point stages 385 MB and peaks at 813 MB;
  `highstat_eff_1/beam_commissioning_00` (3 stations x 8 M events) stages
  2.1 GB and peaks at **4720 MB** -- which exceeded the 4000 MB originally
  requested and survived only because this pool allocates 6000 MB regardless.
  `request_memory` is now 8000 and `request_disk` 8 GB. Size the requests from
  the biggest sub_runs, never the typical ones.
- **`make_joblist.py` reports what it skips.** A sub_run on EOS with no
  `combined_hits_root` is printed, not silently dropped — currently
  `drift_mesh_2d_1/dm_02_02_m420_d530`, whose backup is incomplete.

## 6. Verified on 2026-07-28

- Kerberos reaches the worker (`MY.SendCredential = true`) with a ccache good
  for both reading the nTOF space and **writing back** to it. Products land as
  `dneff:za` via the inherited `sys.owner.auth`, so they bill to the nTOF quota.
- LCG_110 `x86_64-el9-gcc13-opt`: python 3.13.11, uproot 5.7.1, numpy 2.4.4,
  pandas 2.2.3, scipy 1.17.1, matplotlib 3.11.0 — reads our files as-is.
- A condor-produced `alignment.json` matches banco's to ~13 significant figures
  (float rounding from a different BLAS; same 1 520 221 pairs, same RMSE).
- Idle → running ≈ 4 min with ~7300 idle jobs in the pool.
