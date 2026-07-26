# `urw_reference` — the EIC uRWELL reference telescope

Tracking with the two EIC uRWELL planes that bracket the P2 BASKET stations on
the TB_July2026_H4 beamline, and the P2 residual/efficiency measurement that
falls out of it.

```
beam ---------------------------------------------------------------> z
       uRWELL_front    P2_IN      P2_MID      P2_OUT     uRWELL_back
         z = 0        z = 320    z = 630     z = 940      z = 1370   [mm]
         FEU 1         FEU 3      FEU 4       FEU 5        FEU 1
```

**Why this exists.** Stage 22 (`../22_tag_probe_efficiency.py`) can only measure
a P2 station against the *other* P2 stations, so it is limited by their 12 mm
pads and is honest only about efficiency relative to the tag selection. The
uRWELLs are ~1 mm-pitch strip detectors, so once mapped correctly they point at
every P2 plane to about 1 mm — three times better than a P2 pad — and turn the
measurement into a real tracking-referenced efficiency that needs no other P2
plane to have fired. On `highstat_eff_1/beam_commissioning_00` the two methods
give 96.5 / 97.1 / 96.0 % (here) versus 92.5 / 96.3 / 92.5 % (stage 22), with a
*tighter* probe radius.

## Start here

* **`URW_TRACKING_HANDOFF_2026-07-25.md`** — the full document, written for
  someone who has never seen this data, this DAQ or these detectors. Strip maps,
  file formats, the channel → strip wiring, clustering, alignment, tracks, and
  §13 on the P2 measurement. Read §6 before plotting anything.
* `ORDERING.md` — the connector-ordering investigation. Read the box at the top:
  its body predates the 2026-07-26 correction.

## Scripts

| script | what it does |
|---|---|
| `urw_lib.py` | geometry + clustering. `VIEW_MODE_DEFAULT` is the authoritative channel → strip wiring |
| `align_and_track.py` | front↔back alignment and two-point tracks (handoff §8, §9) |
| `loop_subruns.sh` | drives `align_and_track.py` over every sub-run of a run |
| `urw_p2_efficiency.py` | **the measurement** — P2 residuals and efficiency vs uRWELL tracks (§13) |
| `plot_hv_curves.py` | efficiency vs mesh / drift voltage from that stage's CSV |
| `record_mapping_alignment.py` | freezes the mapping + alignment to the analysis dir |
| `write_analysis_readme.py` | generates a README for an output tree |
| `urw_qa.py` | per-sub-run uRWELL QA: occupancy, spectra, cluster size, beam spot, timing |
| `explore*.py`, `check_*.py` | the diagnostics behind the findings; each docstring says what it established. `check_*.py` are **superseded**, kept for the record |

## Running

```bash
unset PYTHONPATH                      # the ISEG HV SDK shadows uproot/ROOT
PY=/local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python
cd $(git rev-parse --show-toplevel)/sps_beam_analysis/urw_reference

$PY align_and_track.py --min-amp 0 --plot align.png
$PY urw_p2_efficiency.py --run highstat_eff_1 \
    --out $SPS_ANALYSIS_ROOT/urw_referenced_efficiency/highstat_eff_1
$PY plot_hv_curves.py --csv <that dir>/urw_p2_efficiency_<run>.csv --out <that dir>
```

The scripts put `sps_beam_analysis` on `sys.path` and call
`sps_config.setup_paths()`, so the shared core (`p2_io`, `p2_mapping`,
`p2_sparks`, `sps_cluster`) is the same code the other stages use — there is no
second copy. The one external dependency is `DetectorConfigLoader` /
`DreamDetector` from `~/dylan/saclay_micromegas`, which `urw_lib` imports by
absolute path.

**Results live under `$SPS_ANALYSIS_ROOT`, not here** — see
`analysis/urw_referenced_efficiency/`, which carries its own READMEs plus the
frozen `MAPPING_AND_ALIGNMENT.md`.

## Two things to know before you change anything

1. **The channel → strip wiring is measured, not assumed.** All four uRWELL
   views are wired differently from the map file, and the back needs a reversal
   of channel order *inside* each connector rather than a swap of the two
   connectors. Getting it wrong is not obvious — the wrong choice left the back
   pointing at 4.4 mm instead of 0.9 mm and looked plausible. Handoff §6.2.
2. **Most `hits_root/` directories have been deleted** (all 23 sub-runs of
   `drift_mesh_scan_1`). `urw_lib.feu_hit_files` falls back to
   `combined_hits_root/`, which holds all four FEUs in one tree, so `iter_hits`
   must be given `feu=` or P2 channels get mapped onto uRWELL strips silently.
   Handoff §5.3.
