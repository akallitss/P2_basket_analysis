# gas_studies/water_scan - which gas makes det2 drift so slowly at low field?

det2 (P2_2) drift-field scan of 24-25 Nov 2026 (mesh 430 V, dV = 25 ... 525 V, nominal
gap 4.0 mm, Ar/iC4H10 95/5 nominal, 50 um aluminised-mylar cathode window). Stage
`cosmic_bench_analysis/44_deform_scan_regions.py` found the timing consistent with dry-gas
Magboltz above ~400 V/cm but much slower below. This folder tests gas compositions with
Magboltz against that curve and gives a drift-velocity table for converting drift times
into gap thicknesses.

## Result (27 Nov 2026)

**Water vapour at 0.7 % (range 0.5-1.2 %) reproduces the drift-time curve over the full scan,
62-1312 V/cm, with one gas and one gap. Dry gas, extra air and extra isobutane do not.** The
amplitude drop at low field is NOT explained by attachment along the drift path. The timing
excludes that, and water has no attachment in Magboltz. That drop needs a separate,
depth-independent cause (see "Open").

Fits of y(dV) = median t_max - t0 over all 21 drift points (dV = 25 ... 525 V). Details in
`results/summary.txt`, `fits_family_profiles.csv`, `fits_per_gas.csv`.

| gas model (chamber) | d = 4 mm, c from MC | + 2 ns syst. | d free (+2 ns) | verdict |
|---|---|---|---|---|
| dry Ar/iC4H10 95/5 | chi2 9223/21 | - | s and d run away (d = 9.9 mm, chi2 969/19) | excluded; also misses the high-field slope (pulls -5...-10 at 500-1300 V/cm) |
| + H2O (fraction fitted) | **0.70 %**, chi2 167/20 | **0.68 %, chi2 17.7/20** | 0.88 %, d = 3.57 mm, 9.4/19; with s free 0.52 %, d = 4.8 mm, 7.8/18 | **fits** |
| + air (N2/O2 79/21) | 0.85 %, chi2 7540/20 | 1.0 %, 3339/20 | d -> 9-11 mm, dt0 -70 ns, >= 289/19 | excluded |
| iC4H10 fraction (5-20 %) | runs to 15-20 %, 2308/20 | 20 % (edge), 353/20 | 8 % at d = 7 mm, dt0 -27 ns (5 sigma on prior), 74/19 | disfavoured (unphysical d, t0) |
| H2O + air (Magboltz eta) | 0.75 % + 0.5 % air: 337/21; + 1 % air: 1173/21 | - | - | excluded by the timing (see below) |

The good-gap regions give the same answer: reg24 regions 5, 11, 17, 20, 22, 23: 0.78 %, chi2 54/20
(18.5/20 with 2 ns syst.); reg8 regions 2, 3, 6, 7: 0.72 %, 101/20. The conclusion does not depend
on the gap non-uniformity.

**Water fraction: 0.70 % (reference: chamber, d = 4 mm, response factor from the MC); 0.52-1.22 %
over data sets and fit variants.** The statistical error is +-0.02 %. The spread comes from the
d-fraction correlation: free d trades a higher water fraction for a smaller gap. The zero-drift t0
also wants a +6...+17 ns shift, in the opposite direction to the expected walk bias.

Data-implied drift velocity (4 mm / T, `vd_models_vs_data.png`): it rises slowly to a maximum of
4.0 cm/us near 450 V/cm, then falls to 3.2 cm/us at 1300 V/cm. That is the 0.7 % H2O curve. Dry
gas peaks at 4.1 cm/us at 220 V/cm and falls to 2.8 cm/us. **Water also makes the gas ~18 %
faster at 750 V/cm** (3.67 vs 3.11 cm/us). The "agreement with dry Magboltz at high field" was
partly hidden in the response factor c = 0.47 assumed earlier. The MC gives c = 0.51-0.56.

### Response factor c

The cluster MC of the leading-pad peak time gives M(T)/T = 0.51 (T = 60 ns) ... 0.56 (350 ns) ...
0.51 (900 ns). Across all response variants (`response_systematics.csv`) the implied T at 62 V/cm
is 681-706 ns and at 125 V/cm 307-315 ns, i.e. +-2 %. The response is not a free parameter at
the level that matters. With s free, the water fits give s = 0.85-0.99 and the same fraction
within 0.06 %. Implied full drift times (chamber, d = 4 mm): 684, 309, 189, 141, 121, 111 ns at
62, 125, 187, 250, 312, 500 V/cm. Dry Magboltz gives 165, 107, 98, 98, 101, 114 ns.

### Attachment and the low-field amplitude

- Median amplitude relative to 750 V/cm: 0.27 / 0.57 / 0.75 / 0.83 at 62 / 125 / 187 / 250 V/cm.
  Ballistic deficit of the long drift times (MC) explains only 0.64 / 0.86 / 0.96 / 0.99.
- If the missing charge were attachment along the drift (exp(-eta z)), the amplitudes would need
  eta*d = 2.36 / 0.89 / 0.51 / 0.38 (eta = 5.9 / 2.2 / 1.3 / 1.0 per cm;
  `attachment_required.csv`). Magboltz gives that size and shape for H2O 0.6-0.75 % + 1 % air.
- **But the timing rules it out.** With that attachment the far-side charge is lost and the
  median peak time at 62 V/cm drops from 365 to 177 ns. A depth-independent loss of the same size
  keeps it at 361 ns (data 365 +- 7) (`loss_mechanism_check.csv`). Scaling the required
  attachment by k with the water fraction refitted, the chamber timing allows k <~ 0.25, the reg24
  set k <~ 0.35 (Delta chi2 = 4, 2 ns syst.; `attachment_scale_scan.csv`). That means eta*d(62 V/cm)
  <~ 0.6, roughly <~ 0.3 % air in the Magboltz O2 convention.
- So in the best model (water, eta = 0) the signal survives at low field as observed. The event
  yield drop (0.37 / 0.42 vs 0.48 per trigger at 62 / 125 V/cm) follows the amplitude drop, which
  comes from a mechanism not described by gas transport.
- Integral/amplitude (6.30 at 62 V/cm) lies between the no-loss MC (7.65) and the attachment MC
  (6.41). Its threshold-defined integration span probably truncates long, low pulses, so it is
  not used to decide.

### Drift time -> gap (`T_vs_gap_dry_and_best.csv`, `gap_from_T_dry_vs_best.csv`)

T = d / v_d(dV/d), full drift time (not c*T). Band = envelope over 0.52-1.22 % H2O.

| dV | T [ns] | d dry [mm] | d 0.70 % H2O [mm] | band [mm] | change |
|---|---|---|---|---|---|
| 300 V | 107 | 3.15 | **3.89** | 3.68-4.27 | +0.74 mm (+24 %) |
| 300 V | 125 | 3.85 | **4.87** | 4.65-5.09 | +1.02 mm (+27 %) |
| 300 V | 170 | 5.98 | **6.94** | 6.44-6.94 | +0.97 mm (+16 %) |
| 300 V | 190 | 7.01 | **7.65** | 6.92-7.77 | +0.64 mm (+9 %) |
| 200 V | 107 / 125 / 170 / 190 | 3.67 / 4.59 / 6.82 / 7.74 | 4.37 / 5.04 / 6.32 / 6.75 | 4.13-4.38 / 4.58-5.11 / 5.43-6.62 / 5.73-7.14 | +19 / +10 / -7 / -13 % |
| 400 V | 107 / 125 / 170 / 190 | 2.98 / 3.56 / 5.27 / 6.18 | 3.42 / 4.25 / 6.67 / 7.67 | 3.29-3.83 / 4.05-4.81 / 6.39-6.94 / 7.45-7.76 | +15 / +19 / +27 / +24 % |

T(d) at dV = 300 V (ns), d = 3 / 4 / 5 / 6 / 7 / 8 mm: dry 103 / 129 / 150 / 170 / 190 / 209;
0.70 % H2O 90 / 109 / 128 / 147 / 172 / 200 (0.52 %: 93 / 113 / 131 / 150 / 171 / 196;
1.22 %: 83 / 101 / 124 / 153 / 194 / 246). In wet gas a gap of 4 mm at dV = 300 V gives T = 109 ns,
not 129 ns. Above ~6.5 mm the field at dV = 200-300 V drops below the shifted v_d maximum. The
conversion then becomes steep and fraction-dependent, so large gaps are better measured at
dV = 400 V (or more).

### Plausibility of the source

0.7 % H2O is ~60 % of the vapour pressure of lab air at 20 C / 50 % RH. With a typical PET
water-vapour transmission of order 1-few g/(m^2 day) for 50 um at room temperature (an
order-of-magnitude literature value, not measured here), 0.12 m^2 lets in ~0.1-0.4 l of vapour
per day. At 1-2 l/h that gives 10^-4-10^-3. Reaching 7x10^-3 needs a lower flow (~0.1-0.5 l/h),
a more permeable window or an additional source (outgassing of the chamber materials). So the
order of magnitude is plausible but not demonstrated. O2 and N2 permeate PET much less than water
per unit partial pressure, so permeation predicts water without significant air. That fits the
timing limit on attachment. A leak would bring ~100x more air than water and is excluded by the
timing.

### Open

- **The low-field amplitude loss** (x0.42 at 62 V/cm beyond ballistic deficit, x0.67 at 125 V/cm)
  is depth-independent. Candidates: a gain transient (the low-field points were the first taken,
  16:49-18:21, right after the mesh scan had run the mesh at 340-420 V), or a field effect at the
  mesh at very low drift/amplification field ratio. Test: repeat 62 / 125 V/cm at the end of a
  scan, or interleave high and low dV.
- **Not tested**: CO2 or other molecular admixtures with large low-energy inelastic
  cross-sections (CO2 at the per-cent level gives a similar shift of the v_d maximum). Water is
  favoured because <1 % suffices and a source exists. A hygrometer / O2 meter on the exhaust and
  a flow scan would settle it: the fraction should scale ~1/flow if it is permeation.
- sigma68(t_max): water reproduces its rise (dry gas misses it by a factor 1.5-7 at 250-62 V/cm).
  The data are ~20-40 % above the water model at 125-300 V/cm (chi2 852/20 with sigma0 free). The
  MC underestimates the spread (track angles, partial pad coverage, gap non-uniformity within
  the chamber).

### Caveats

- Magboltz H2O ("H2O anisotropic 2017"): the source states that the measured data in water vapour
  are not accurate and that its set reproduces v_d and diffusion to ~2 %. In mixtures with argon
  the low-energy rotational/vibrational cross-sections dominate. A 10-20 % error on them maps
  roughly to the same relative error on the fraction.
- Magboltz O2 3-body attachment uses the full gas density with the O2 third-body efficiency for
  every partner (T3B = 1). Ar is a much weaker third body, H2O and hydrocarbons much stronger, so
  "% air" from attachment is uncertain by a factor of a few. The timing limit is on eta*d itself,
  not on the air fraction.
- Temperature / pressure fixed at 293.15 K / 760 Torr. v_d depends on E/N: a 2 % density change
  shifts the effective field by 2 %. That is negligible against the factor-4 low-field effect but
  is part of the +-few % in the conversion.
- The response MC is one leading pad with full-depth ionisation (cosmic angles and edge sharing
  not modelled). The SER comes from the stage-47 mean pulse with T_ref = 125 ns (an input of
  that stage). Varying T_ref 100-150 ns changes the implied T by <1 %.
- The chamber is a mixture of regions whose gap differs by up to +50 %. The chamber "d" is an
  effective value. The good-gap subsets were selected on a gap estimate that itself used dry
  Magboltz at high field. The selection is relative, so it holds, but those absolute gaps scale
  with v_wet / v_dry (+15-20 % at 625-1300 V/cm).
- Low-field statistics: 62 V/cm has 1985 events (A >= 300), median error 6.6 ns. The 2 ns
  systematic per point covers the overnight sequence (gas or charging drifts are confounded with
  field). Without it, chi2/ndf = 8 for the best water fit.
- dt0 (zero-drift t0) floats with a 5 ns prior. The fits want +6 ns (d = 4) to +17 ns (d free).
- One combination (0.75 % H2O + 2 % air) stalled at its first field point in Magboltz and was
  dropped. The 1 % air combinations already bracket it.

## What was run

| step | file | where / how |
|---|---|---|
| Magboltz scan | `water_scan.cpp`, `CMakeLists.txt`, `gaslist.txt`, `gaslist_combos.txt`, `run_lxplus.sh` | lxplus955 + lxplus999, LCG_109 view (`x86_64-el9-gcc13-opt`), Garfield++ `78fe1bd3` with **Magboltz 11.19** (H2O = "H2O anisotropic 2017" set with rotational dipole levels; O2 2016 with 3-body attachment). 293.15 K, 760 Torr, thermal motion on, 1e7 collisions per point, 25 fields 20-1500 V/cm. The local build (`~/Software/garfieldpp`, Garfield++ git `22af4690`, same Magboltz 11.19) gives the same v_d at 100 V/cm (0.003398 vs 0.003400 cm/ns in `garfield_inputs/gas_table.csv`) but is ~3x slower, so it was only used for the check. |
| pulse-response MC | `pulse_response.py` -> `results/pulse_response_table.csv` | local. Single-electron response (SER) from the chamber mean pulse (stage 47, `mean_pulse_ref.csv` = box(125 ns) (x) SER) by a CR-RC^n + undershoot fit through the box (rms residual 0.15 % of the peak). Cluster MC of one leading pad (2.8 clusters/mm, Ar-like cluster sizes, Polya gain, longitudinal diffusion, 60 ns sampling with random phase, parabola peak like `WaveformAnalyzer`, 32-sample window, A >= 300 ADC cut, binomial attachment loss exp(-eta z)). Tabulates, for full drift time T and attenuation eta*d, the median peak time minus the T = 0 peak time M(T), sigma68, median amplitude and integral/amplitude. |
| response systematics | `sens_response.py` -> `results/response_systematics.csv` | SER from T_ref = 100 / 150 ns or exact deconvolution, cluster density 2.0 / 3.6 per mm, Polya theta 0.5 / 3, diffusion x0.5 / x2, amplitude cut 150 / 600 ADC. |
| fits | `analyze_water.py` (model, data, fit), `make_results.py` (all tables and figures) | local, `.venv/bin/python`. |

### Data observable and model

- `y(dV)` = median `time_of_max` of the leading pad (A >= 300 ADC, not saturated or truncated)
  minus the zero-drift reference t0 (393.09 +- 0.54 ns, chamber; per region for the region sets).
  Errors: bootstrap from `region_stats.csv`. Data sets: `chamber` (all pads), `reg24_good`
  (reg24 regions whose high-field gap estimate is within 5 % of 4.0 mm: 5, 11, 17, 20, 22, 23),
  `reg8_good` (reg8 regions 2, 3, 6, 7).
- Model: `E = dV/d`, `T = d / v_d(E)`, `att = eta(E) d`, `y = s * M(T, att) + dt0`.
  `s = 1` means the response factor c(T) = M/T comes from the MC ("c fixed"), `s` free is an extra
  empirical factor ("c free"). dt0 has a Gaussian prior of 5 ns (the zero-drift t0 is
  measured on small pulses, median 82 ADC, so walk can shift it by a few ns). d = 4 mm fixed or free.
- Contaminant fractions are fitted continuously: ln v_d and eta are interpolated linearly in the
  fraction between Magboltz grid points. The fraction interval is Delta chi2 = max(1, chi2/ndf).
  Fits are repeated with a 2 ns systematic added in quadrature to every point ("sys 2").

## Files in `results/`

| file | content |
|---|---|
| `magboltz/<gas>.csv` (+ `magboltz_logs.tar.gz`) | Magboltz output per gas: E, v_d (cm/ns) and error (%), D_L, D_T (sqrt(cm)), alpha, eta (1/cm) |
| `magboltz_all_gases.csv` | all gases in one table |
| `pulse_response_table.csv` | MC response M, S68, amplitude, ioa vs (T, eta d) |
| `response_systematics.csv` | c(T) and data-implied T for each response-model variant |
| `data_implied_drift_time.csv` | data y, implied full drift time T and v_d = 4 mm / T, dry-gas T |
| `fits_per_gas.csv` | fit of every Magboltz gas, 4 variants, dV >= 25 and dV >= 200 V |
| `fits_family_profiles.csv` | best fraction per family (H2O, air, iC4H10), variant, data set |
| `attachment_required.csv` | eta d needed by the amplitude drop, given T from the timing |
| `shape_amplitude_check.csv` | amplitude, ioa, sigma68 predicted by the timing-fit models |
| `vd_table_dry_and_best.csv` | E -> v_d for dry gas, best water fit and its +-1 sigma fractions |
| `T_vs_gap_dry_and_best.csv` | T = d / v_d(dV/d) for d = 2.5 ... 9 mm at dV = 200, 300, 400 V |
| `gap_from_T_dry_vs_best.csv` | gap d for T = 107, 125, 170, 190 ns, dry vs best gas |
| `summary.txt` | the printed fit summary |
| `*.png` | `fit_data_vs_models`, `chi2_profiles`, `vd_models_vs_data`, `amplitude_shape_check`, `T_to_gap_conversion` |

## Reproduce

```
cd gas_studies/water_scan
./run_lxplus.sh launch                       # Magboltz, gaslist.txt (detached on lxplus)
SKIP_BUILD=1 GASLIST=$PWD/gaslist_combos.txt NPAR=4 ./run_lxplus.sh launch
./run_lxplus.sh status ; ./run_lxplus.sh fetch
../../.venv/bin/python pulse_response.py results/pulse_response_table.csv
../../.venv/bin/python sens_response.py
../../.venv/bin/python make_results.py
```
