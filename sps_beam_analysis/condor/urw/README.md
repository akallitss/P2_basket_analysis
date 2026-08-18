# uRWELL-referenced efficiency on HTCondor

The copy of record of the `~/p2_condor_urw` pipeline on lxplus, which runs
`urw_reference/urw_p2_efficiency.py` one sub_run per job straight off the nTOF
EOS tree and uploads to `analysis/urw_timebins/<run>/<sub_run>/`.  It lived only
on lxplus until 2026-08-17; if that AFS home is lost, this is what rebuilds it.

```bash
# on lxplus, from ~/p2_condor_urw
condor_submit urwtb_hvscan.sub          # queue … from urwtb_hvscan.txt
# each line:  <run> <sub_run> <time_bins>      (0 = one number per sub_run)
```

Then, locally:

```bash
cd sps_beam_analysis/urw_reference
python3 merge_timebin_products.py --run low_mesh_scan_1   # -> one scan CSV
```

`urwtb_hvscan.txt` is the 2026-08-17 job list: the three HV scans that had never
been put through the uRWELL reference — `low_mesh_scan_1` (MID/OUT mesh 385 → 330
V), `p2in_hvrange_1` and `p2in_hvrange_2` (the replacement P2_IN chamber, 200 →
450 V).  29 of 30 landed; `low_mesh_scan_1/nominal_00` fails inside the stage on
an empty array and is redundant anyway (390/690 is also the first point of
`drift_mesh_scan_1`).

## The one correction the job applies

`p2in_hvrange_2/run_config.json` — and only that run, out of 74 — records
`det_type: "P2"` for **both uRWELL reference planes**, so the tracking cannot
build its geometry and every job dies on a missing `P2.json`.  `urwtb_job.sh`
repairs the *staged* copy (`urw_inter` for the front plane, `urw_strip` for the
back) before running; EOS is never modified.  Anything else that reads that run
directly will hit the same wall.
