# MPGD26 beam-report figures

The 26 figures under `reports/mpgd26_sps_beam_2026-08/figs/`. They were built on
2026-08-16 in a scratch tree that no longer exists, which is why for a day
nothing in the repo could re-make them; the scripts were recovered and are now
here.

```bash
./fetch_from_eos.sh          # workspace from the nTOF EOS copy (~90 MB)
./fetch_from_eos.sh --hv     # + HV monitor traces (136 MB) for fig_stability
python3 aggregate.py         # products -> report_data/*.csv
python3 fig_efficiency.py    # and fig_timing / fig_vmm / fig_padmap /
                             # fig_charging / fig_stability
```

`paths.py` names the workspace (`$MPGD26_WORKSPACE`) and the output directory
(`$MPGD26_FIGS`, default = the report's own `figs/`). Figures are written where
they belong and only then copied into `conference/` by its gather script.

**EOS, not the LaCie.** The LaCie copy of `analysis/` predates the condor
re-processing: its stage-29 products still use the pre-rename column names
(`mesh_in`/`mesh_midout` instead of `mesh_v_<det>`) and it holds 1 timing run
where EOS has 16. `fetch_from_eos.sh` replaces the tree, symlink included.

| script | figures |
|---|---|
| `fig_efficiency.py` | ε(mesh), ε(drift), method comparison, the 2D surface as boxes **and** as curves, the P2_IN chamber swap |
| `fig_timing.py` | σ(mesh), σ(drift) + Magboltz, the 2D surface both ways, the correction ladder |
| `fig_vmm.py` | gas A/B turn-on, drift latency, config scan, the superseded DREAM caveat |
| `fig_padmap.py` | per-pad efficiency at the working point |
| `fig_stability.py` | campaign HV timeline, sparks, efficiency drift (needs `--hv`) |
| `fig_charging.py` | efficiency vs time within sub-runs |
| `make_hv_setpoints.py` | the (run, sub_run, det) → mesh/drift table every voltage axis joins against |

`p2style.finish()` writes PDF beside every PNG.

ε(mesh) spans 330–450 V because `low_mesh_scan_1` and `p2in_hvrange_1/2` were
put through the uRWELL reference on 2026-08-17 (`../condor/urw/`). The three
runs are joined per station **only where it is the same chamber** — the
replacement P2_IN is its own series, never continued from the CERN-built
chamber it replaced.
