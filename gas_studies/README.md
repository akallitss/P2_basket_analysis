# gas_studies — HV / gain equivalence between gases (Magboltz + Garfield++)

Map operating high voltages between gas mixtures for the **P2 Micromegas**, so
that a detector calibrated on one gas can be run at equivalent **gain** on
another. Built for the SPS beam test where we switch from the Saclay lab gas
**Ar/iC₄H₁₀ 95/5** to the NSW gas **Ar/CO₂/iC₄H₁₀ 95/3/2**, and designed so a
**new gas** can be characterised in minutes at the beam.

Detector geometry (fixed):
- amplification gap **d_amp = 150 µm = 0.015 cm**
- drift/conversion gap **d_drift = 3 mm = 0.30 cm**

## What it produces

1. **Gain-equivalence HV mapping.** Gas gain in the parallel-plate amplification
   gap is `G(V_mesh) = exp[(α − η)(E) · d_amp]`, with `E = V_mesh / d_amp` and
   the Townsend `α` / attachment `η` coefficients from Magboltz. For every
   reference mesh voltage we invert `G` of the candidate gas to find the mesh
   voltage giving the **same gain** → `figs/hv_equivalence.png`,
   `results/hv_equivalence.csv`, and headline numbers in
   `results/gain_summary.txt` (e.g. "415 V in Ar/iso ≈ X V in NSW gas").
2. **Drift velocity vs the mesh→drift HV difference** across the 3 mm gap,
   both gases overlaid → `figs/drift_velocity_vs_dV.png`
   (x-axis is `ΔV = E · d_drift`, i.e. the HV difference you dial between mesh
   and drift electrode).
3. Diagnostics: `figs/gain_vs_vmesh.png`, `figs/townsend_vs_E.png`.

## Layout

| file | role |
|------|------|
| `gases.py`       | **single source of truth** — gas registry (name, Magboltz composition, role) + geometry constants |
| `p2_gas_scan.cpp`| Garfield++/Magboltz scanner: reads a composition on argv, scans amp (Townsend/attachment) + drift (v_d) fields → one CSV per gas |
| `CMakeLists.txt` | builds the scanner against an LCG Garfield++ view |
| `run_lxplus.sh`  | ssh to lxplus, build, run the scan for each gas, fetch CSVs into `results/` |
| `analyze.py`     | reads `results/*.csv` → the plots and mapping tables above |
| `results/`       | Magboltz CSVs (committed) + `hv_equivalence.csv`, `gain_summary.txt` |
| `figs/`          | output plots |

## Run it

```bash
cd gas_studies
./run_lxplus.sh                 # build + scan every gas in gases.py, on lxplus
python3 analyze.py              # build the plots + mapping (uses the working point below)
```

Reference working point defaults to **mesh 415 V** (det1 working point). Change
with `python3 analyze.py --ref-vmesh 415 --working-points 375 400 415 450`.

**Runtime.** The low-field drift points are fast (~15 s each); the high-field
amplification points are the slow part — Magboltz's steady-state Townsend method
takes a few minutes per point. The default grid is **6 amplification points**
(≈15–25 min/gas). To trade accuracy for speed at the beam, edit `ampV` /
`ncoll_amp` in `p2_gas_scan.cpp`; the gain curve is smooth so even 5 points
interpolate well.

### Adding a NEW gas at the beam

1. Append an entry to `GASES` in `gases.py`, e.g.
   ```python
   "ar_co2_93_7": {
       "label": "Ar/CO2 93/7",
       "magboltz": [("ar", 93.0), ("co2", 7.0)],
       "role": "candidate",
   },
   ```
   (Magboltz names: `ar`, `co2`, `ch4`, `ic4h10`, `cf4`, `c2h6`, `n2`, `xe`, …;
   up to 6 components.)
2. Scan just that gas and re-plot:
   ```bash
   ./run_lxplus.sh ar_co2_93_7
   python3 analyze.py
   ```

## Environment / how the remote build works

- `run_lxplus.sh` uses `ssh lxplus` in `BatchMode` (needs a working
  Kerberos ticket / key — run `kinit`/`ssh lxplus` once if it prompts).
- Garfield++ comes from an **LCG view** (default
  `LCG_105/x86_64-el9-gcc13-opt`); override with `LCG_VIEW=…`. Only Magboltz is
  used (no Heed gas-file database needed).
- lxplus is load-balanced (node-local `/tmp`, shared AFS `$HOME`). Sources are
  staged under `~/p2_gas_scan` on AFS; the build + run happen in a node-local
  `/tmp` dir (AFS quota is tight) and only the tiny CSVs are copied back.

## Method & caveats

- Gain uses the **parallel-plate / Rose–Korff** approximation
  `G = exp[(α−η) d_amp]` at uniform field `E = V_mesh/d_amp`. This is the right
  tool for *relative* HV↔gain mapping between gases; absolute gain also carries
  a mesh-transparency / field-line factor that largely cancels in the ratio.
- `α`, `η`, `v_d` are from Magboltz at 20 °C, 760 Torr (`gases.py` /
  `p2_gas_scan.cpp` if you need other T/p).
- Consistent with the timing-study transport in
  `cosmic_bench_analysis/garfield_inputs/` (same RunMagboltz call).
