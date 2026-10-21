# VMM threshold theory and threshold/masking trade-off

Code behind the unlisted note
`cern-site/notes/unlisted/2026-10-21-testing-the-vmm-threshold-theory-step-by-step.html`
(https://alexandra-kallitsopoulou.web.cern.ch/notes/unlisted/2026-10-21-testing-the-vmm-threshold-theory-step-by-step.html).
The figure scripts write straight into that note's `_files/` folder, so running one overwrites the published PNGs
(commit + `scripts/deploy.sh` in cern-site to publish).

Question: is VMM efficiency (~85 %) below DREAM (95.6 %) because the VMM threshold is too high for the gain,
and how far can it come down before noisy channels have to be masked?

## Scripts

| script | what | note |
|---|---|---|
| `gain_step_test.py`, `mesh_test.py` | does the run_46 threshold model predict the gain step, VMM 19 sdt step, mesh scan | Part 1 |
| `note_figs.py` | Part 1 figures | needs `eff/` |
| `note_figs2.py` | Part 2: sdt → mV → ADC → fC, noise floor, trim (figs 6–12) | needs `vmm_adc_all.npz`, `Threshold_093415.csv` |
| `sum_vmm_adc.py` | sums `counts.npz` of every EOS sub-run → `vmm_adc_all.npz` | run on lxplus |
| `build_thrscan_summary.py` | bench scan ROOT trees → `thrscan_summary.json` | needs `thrscan/qa/` |
| `tradeoff.py` | core: per-channel noise onsets (`onsets(gain, R)`), chip baselines `BASE`, efficiency model `eff_d` | Part 3 |
| `tradeoff3.py` | trimmed scenario | Part 3 J |
| `note_figs3.py` | figs 13–16 (noise maps, bench vs beam, P2_OUT trade-off, station + trim) | |
| `note_figs3b.py` | figs 15b/15c (P2_IN, P2_MID); **no main guard, importing it redraws them** | |
| `note_figs3c.py` | figs 17/18: best efficiency and masks per VMM id. `GAIN=3.0 python3 note_figs3c.py` → `_g3` versions | |
| `oos/oos_bin.py` | **lxplus**: decode every pcapng of a sub-run, hits per (vmm, ch) per 0.1 s → `run46_tbin.npz` (~5 min for run_46) | |
| `oos/oos_load.py` | spills from VMM 0 ch 44, in/out-of-spill rates per channel (`load()`) | |
| `oos/oos_adc.py` | **lxplus**, second pass: out-of-spill ADC spectra and ±150 ns coincidences → `run46_adc.npz` | |
| `note_figs3d.py` | fig 14b: bench vs run_46 between-spill hits | |

lxplus recipe (oos): copy `oos_bin.py`/`oos_adc.py`, `../vmm_dream_matching/extract_vmm_triggers.py` and
`P2_basket_online_analysis/vmm_decode.py` into `~/p2_oos/`, `source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh`,
`export PYTHONPATH=$HOME/p2_oos:$PYTHONPATH` (append, LCG needs its own), run detached with `nohup`.
`oos_adc.py` needs `run46_masks.npz` (bins written from `oos_load.load()`).

## Data (gitignored, kept next to the scripts on this machine)

`thrscan/` (bench scan configs, logs, qa ROOT trees, 125 MB), `eff/` (EOS `vmm/matching/efficiency_v2` JSONs),
`vmm_adc_all.npz`, `oos/*.npz`, the vmm_daplxa calibration CSVs. Committed: `thrscan_summary.json`, `runcfg/`
(per-run per-chip sdt), the `*_numbers.json` / `part3_allchips*.json` outputs.

## Results so far (2026-10-22)

- No trim, 4.5 mV/fC, P2_OUT: best sdt 14:220 15:≤210 16:220 17:≤210 18:≤210 19:230 with 3 masks → model 91.2 % vs July 83.7 % (measured 85.3 %).
  With trim 92–94 %. 3.0 mV/fC is worse everywhere (P2_OUT ≤ 87.7 %).
- All 18 chips: 7 masks of 1152 outside the walking-burst chips VMM 4/7/9 (whose hot channel hops; fixed masks won't hold).
- Between spills in run_46 (fig 14b): the bench picks the right noisy channels, but they are ~100× hotter in the beam hall
  → beam noise onsets are 5–10 DAC higher than on the bench; P2_OUT chip-average optimum 89.9 → 87.7 / 85.2 % (July 80.4 %).
  Separate real-particle floor on every pad between spills (OUT 7.7 / MID 1.5 / IN 0.4 Hz), absent on the bench.
  VMM 10 ch 22 has 0 hits in run_46.

## Open / next

1. In-beam threshold scan using only between-spill hits (a few min per sdt), down to sdt ~170, steps 2–3.
2. vmmsc pedestal + threshold calibration on the P2 hybrids (per-channel baselines, trims).
3. Same out-of-spill check on run_48/52/66 (July "opt" settings) — `oos_bin.py run_52 <sub> ...`.
4. What is the between-spill particle floor (activation near P2_OUT?) — check its time dependence across the day.
