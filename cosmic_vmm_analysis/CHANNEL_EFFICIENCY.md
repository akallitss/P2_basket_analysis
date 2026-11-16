# Single-channel efficiency on the VMM cosmic bench

How efficient is ONE read-out channel (one pad) at a given threshold? This document explains the measurement in plain
words, gives every command to run it on new data, lists the checks that tell you whether the answer can be trusted,
and ends with a glossary. Code: `channel_efficiency.py` (this directory). Explainer with animations (CERN site):
<https://alexandra-kallitsopoulou.web.cern.ch/notes/2026-11-17-single-channel-efficiency-on-the-vmm-cosmic-bench.html>.

## 1. The measurement in five sentences

1. A muon seen by two chambers predicts where it crossed the third one (the *probe*): a straight line through two
   points. Only muons that also fired the scintillators are used (*scintillator-confirmed tags*), because two chambers
   firing at the same time without the scintillators are often correlated noise, not muons.
2. The prediction is not exact: it is wrong by sigma ≈ 3 mm when the probe sits between the two others (P2_MID) and
   ≈ 8 mm when it sits outside (P2_IN, P2_OUT). The pads are 12 mm wide, so "did the predicted pad fire?" often looks
   at the wrong pad: that number measures the tracking, not the pad (it gave 0.1-0.3).
3. Counting fixes it: the prediction error moves muons around, it does not create or remove them. The number of
   muons that really crossed a pad = (local muon density) x (pad area), and the density is counted from the predicted
   points in a 40 mm circle around the pad.
4. Count, for growing circles of radius R around the pad, how often the pad fired (random hits from an off-time
   window removed). Divided by the expected crossings, this curve rises while R is still smaller than the prediction
   error and levels off at the efficiency once all the pad's own muons are inside the circle.
5. Beyond that radius (R_lo, default 40 mm) the curve only grows slowly, as b R^2 (hits that come with a muon that
   crossed elsewhere); the efficiency is the value at R_lo minus that slow growth. No model of the prediction error is
   needed - which matters, because the predicted points are not spread smoothly: the tag clusters sit on pad centres
   and the three chambers share the same pad lattice, so the predictions sit on (P2_IN/OUT, ~55 %) or half-way
   between (P2_MID) pad centres. A Gaussian resolution model fails on that (tried: chi2/ndf 40-60).

Formula (per group of pads; *E* = expected crossings, from the tag density counted in pad cells; *K(R)* = pad firings
from tags predicted within *R*, random hits from the off-time window removed):

    K(R) = eps * E + b * R^2     for R >= R_lo,     eps = (K(R_lo) - b R_lo^2) / E,   b from the increments beyond R_lo

The density is counted in pad cells (each tag belongs to its nearest live pad; tags in the cells within 40 mm / the
area of those cells): exact for both smooth and lattice-like predictions, and it follows edges and masked pads.

## 2. Running it

All commands are complete. `$A` is the analysis directory; run on lxplus for the decoding step (the raw data are on
EOS and a run needs a few GB of RAM).

**Step 0 - environment (lxplus).** Loads Python with numpy/pandas/scipy (the "LCG view") and tells the code where
the pad map is.

    source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
    export VMM_QA_DIR=~/tpc/vmm_qa          # folder with vmm_stations.py + P2_BASKET_mapping.csv
    cd ~/cva_own                            # a copy of this cosmic_vmm_analysis folder

**Step 1 - decode the runs into per-muon records** (`dump`). For each run it reads every capture file, builds the
tags, and writes `tags_runN.npz` (one record per tag: predicted point, scintillator time, every pad that fired within
80 mm, in the muon window and in an off-time window) plus `eff_runN.json` (run summary, list of dead/masked pads).

    python3 channel_efficiency.py dump --runs-dir /eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/runs \
        --out ~/wk2tags80 --jobs 3 181 184 187 190 193 196 199 202 205 208 211 214 218 221 223 225

`--jobs 3`: three runs at the same time. `--mem 10G`: memory cap per run (default; the job is stopped instead of
crashing the machine). `--scint-ch 0:40`: the scintillator coincidence channel (H0 VMM 0, channel 40). `--conn-shift`,
`--order`, `--dead-vmm`: the read-out map, defaults = the one used since run_112. About 1 min per run.

**Step 2 - threshold of each chip in each run** (`sdt-table`, on the DAQ box where the configuration files are).

    python3 channel_efficiency.py sdt-table /local/p2/p2testbench/TestBenchCERN/config_ext \
        p2b-config-cern-ext_gain4.5_peaktime200_opt_cosmic_safe_m9c16to21_m4c36to41_h0open_adwk2_K0_i1.txt -o run_sdt.json

(one file name per configuration used; the names are in each run's `eff_runN.json` title.)

**Step 3 - threshold and baseline of every channel in mV** (`thr-table`). Needs one monitor map (baseline and
threshold per channel, measured by `vmm_monitor_tool` with one configuration) and the mV-per-sdt calibration of each
chip; every other run is shifted by its own sdt. A run with its own measured map (e.g. trimmed files) is given with
`--run-map RUN=map.json` and used as is.

    python3 channel_efficiency.py thr-table ~/wk2tags80 --map map_Km18_i5_unmasked.json --sdtcal sdtcal.json \
        --run-sdt run_sdt.json -o thr.csv

**Step 4 - the efficiency** (`measure`).

    python3 channel_efficiency.py measure ~/wk2tags80 thr.csv out_eff --scint 140:310
    python3 channel_efficiency_report.py out_eff ~/wk2tags80 out_eff/figs      # verification figures + numbers.json

`--scint 140:310`: keep tags whose time is 140-310 ns after a scintillator coincidence (the printed histogram must show
the muon peak inside). `--rlo 40`: radius from which all of a pad's own muons are counted; the output also gives the
result for R_lo = 24 ... 56 mm (`eps_rloX` columns), which must agree once R_lo is larger than the prediction error.
`--clean 40`: leave out pads with a dead or masked pad
within 40 mm (the density must be smooth). `--bins`: threshold - baseline bins in mV. `--plateau 18 40`: range of the
summary line. `--per-channel`: one fit per channel (needs much more data than one week).

Outputs: `efficiency.csv` (one row per group: kind = run / distance / plateau / blind / channel; eps, err, b,
chi2/ndf of the b R^2 line beyond R_lo, pads, eps_rloX), `pads.csv`, `curves.npz` (all counts, to redo any grouping), `efficiency.png` (left: the counting
curve and the fit; right: efficiency against threshold - baseline).

**Closure test (always after changing the code)**: `python3 test_channel_efficiency.py` (3 s) and
`python3 channel_efficiency.py toy --sigmas 3 8 12` (simulated muons with a known efficiency).

## 3. When to trust the number (checks, in the order to look at them)

| check | where | good |
|---|---|---|
| scintillator peak inside the `--scint` window | printed histogram | peak at +150..+300 ns, nothing outside |
| R_lo scan | `eps_rlo24 ... eps_rlo56` | flat from R_lo = 40 on; a steady rise = scattering tails: quote the 48/56 value, the spread as systematic |
| line beyond R_lo | `chi2ndf` | ≲ 2 (P2_MID 0.6); 4-6 on P2_IN/OUT = the same tails |
| counting curve | `efficiency.png` left | rises, then nearly flat (a slow straight rise in R^2) |
| blind pads (baseline above threshold) | `kind == blind` | 0 within the error |
| the three chambers | `kind == plateau` | similar efficiency |
| repeated runs at the same threshold | `kind == run` | agree within the errors |

Known limits: (1) tags that are not muons lower the efficiency by exactly their fraction, which is why only
scintillator-confirmed tags are used; (2) the quoted errors are conservative (toy spread ≈ 0.6 of them); (2b) the tails of the
prediction error (multiple scattering, amplified in the extrapolated chambers) are not all inside 40 mm: +0.02..0.03
on P2_IN/OUT at R_lo = 56;
(3) pads near masked pads are dropped, so low thresholds (many masks) keep fewer pads; (4) the baseline map should be
measured close in time to the runs.

## 4. Results of the wk2 scan (3-5 Nov 2026)

Runs 181-225 (16 long runs, 12 thresholds K = -16 ... +64), 26 746 scintillator-confirmed tags (10 % of all tags),
outputs `data/Lab_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_wk2/channel_eff/` (`out/`, `figs/`).

| chamber | efficiency, pads 18-40 mV above baseline | R_lo 24 / 40 / 56 mm | chi2/ndf |
|---|---|---|---|
| P2_IN  | 0.75 ± 0.02 | 0.71 / 0.75 / 0.77 | 5.9 |
| P2_MID | 0.75 ± 0.03 | 0.74 / 0.75 / 0.75 | 0.6 |
| P2_OUT | 0.79 ± 0.02 | 0.76 / 0.79 / 0.80 | 3.8 |

Against threshold - baseline (all chambers alike within errors): ~0.75-0.8 from 20 to ~45 mV, ~0.6-0.7 at 45-60 mV,
~0.45-0.65 above 60 mV. Blind pads (baseline above threshold, VMM 4/7/9 blocks): 0.000 ± 0.002 (470 pad-runs). The
three repeated runs at K = 0 (187, 196, 218; three days apart) agree within errors in every chamber.

The same pads by other methods: naive "did the predicted pad fire" 0.17 / 0.40 / 0.20 (scintillator tags) or
0.09 / 0.28 / 0.12 (all tags); counting with all tags (no scintillator) 0.44 / 0.57 / 0.51 - the self-made tags are
~30-40 % non-muons. Independent check: the chamber-level scintillator efficiency (any hit within 30 mm, live pads,
`analysis_wk2/eff_definitions/scan.txt`) at K = 0: 0.74 / 0.68-0.71 / 0.80-0.82.

Closure (toy, true 0.85): 0.85 ± 0.01 for prediction errors of 3, 8 and 12 mm, with smooth or lattice predictions and
with random hits; 15 % tails at 20 mm: 0.83-0.84 at R_lo = 40, 0.85 at R_lo = 56; 30 % fake tags: 0.60 (= 0.85 x 0.7).

## 5. Glossary

- **tag**: a muon candidate made of a clean cluster in each of two chambers at the right time difference.
- **probe chamber**: the third chamber, whose pads are being measured.
- **on-time / off-time tags**: tag pairs with the muon time difference, or a clearly wrong one (pure chance); the
  off-time ones, scaled, are subtracted.
- **sideband window**: the same tags, but pad hits looked for 4 µs later: only random hits are there.
- **scintillator-confirmed**: the tag time is 140-310 ns after a hit of the scintillator coincidence channel.
- **prediction error (resolution) sigma**: how far the predicted point is from where the muon really crossed.
- **containment C(R)**: fraction of the muons crossing a pad whose prediction lands within R of the pad centre.
- **threshold - baseline (distance)**: how far above its own baseline a channel's threshold sits, in mV.
- **sdt**: the chip's global threshold DAC; one count ≈ 0.75 mV.
- **monitor map**: baseline and threshold of every channel measured on the chip (vmm_monitor_tool).
- **blind pad**: baseline above the threshold, the channel cannot fire on signals.
- **LCG view**: a ready-made Python installation on CVMFS (lxplus).
- **npz / csv / json**: numpy archive / table / structured text output files.
