# VMM Configuration Scan SNR Analysis

Signal-to-Noise Ratio analysis for VMM3a ASIC configuration scans on P2 BASKET detector data from the CERN SPS beam test.

## Overview

The pipeline sweeps two VMM shaping parameters — **gain** (`sg`) and **peaking time** (`snt`) — and finds the configuration
that maximises SNR. Each configuration is represented by a pair of runs:


| Run type | VMM setting | Selects hits with |
|----------|-------------|-------------------|
| Signal run (`sng=0`) | Neighbor triggering **off** | `over_threshold == 1` (hardware discriminator fired) |
| Noise run (`sng=1`)  | Neighbor triggering **on**  | `over_threshold == 0` (sub-threshold, noise-like) |

SNR is computed at two levels of granularity:

- **VMM level** — one SNR value per VMM per configuration pair
- **Channel level** — one SNR value per (VMM, channel) per configuration pair

Both levels use identical methodology and share a common noise baseline.

## Repository Structure

```
vmm_config_scan_analysis.py  — main entry point; orchestrates the full pipeline
vmm_snr.py                   — SNR computation (VMM-level and channel-level)
vmm_noise.py                 — noise baseline estimation (MAD-based)
vmm_signal.py                — signal extraction and MPV estimation
vmm_io.py                    — data loading and run management (ROOT files)
vmm_plots.py                 — all plotting functions
vmm_mapping.py               — detector geometry for CERN/SPS data
vmm_mapping_lab.py           — detector geometry for lab data
vmm_qa.py                    — quality assurance investigations
test_snr_methods.py          — synthetic smoke test (no data files needed)
```

## Detector Layout

| Detector | VMM IDs (CERN) | Channels |
|----------|----------------|----------|
| P2 Large Detector | 12, 13, 14, 15 | 0–63 per VMM |
| P2 Small Detector 1 | 10 (conn. 0), 11 (conn. 1) | ch 0–31, ch 14–63 |
| P2 Small Detector 3 | 8 (conn. 0), 9 (conn. 1) | ch 0–31, ch 14–63 |
| Trigger (excluded) | 0, 1 | — |

Lab data uses `vmm_mapping_lab.py` (different VMM IDs, fewer detectors).

## Running the Analysis

### Switching datasets

Set **one line** at the top of the user configuration block:

```python
dataset = "15kHz"   # "15kHz" | "5kHz" | "lab"
```

All paths, run-table CSV, QA run numbers, plot subdirectory, and VMM mapping are
loaded automatically from the `DATASETS` registry at the top of the file.
To add a new campaign, add one entry there — nothing else changes.

### Run mode

```python
mode = "analysis"   # save all plots to plot_dir, no interactive windows
mode = "debug"      # show each plot interactively (close window to advance),
                    # nothing saved — use the PLOTS dict to select which to step through
mode = "both"       # save AND show interactively
```

### Files per run

```python
n_root_files = 10   # streaming analysis — safe at any N, memory is O(n_vmms × 1024)
n_qa_files   = 2    # QA / legacy plots — loaded fully into RAM, keep small
```

Runs with fewer files than `n_root_files` automatically use all available.

### Enabling / disabling plots

Toggle entries in the `PLOTS` dict at module level. Useful in `debug` mode to
step through only the plots you care about.

### Running

```bash
python vmm_config_scan_analysis.py
```

The run-table CSV (e.g. `vmm_config_scan_15kHz.csv`) must exist in `cnfg_dir`.

## Analysis Methodology

### 1. Noise Estimation (from `sng=1` run)

For each VMM (or channel), noise hits are selected with:
- `over_threshold == 0`
- `adc > 20` (removes the VMM3a ADC=16 digital floor artifact)

The noise width is estimated via the **Median Absolute Deviation (MAD)**:

```
robust_sigma = 1.4826 × MAD
noise_cut    = median + 5 × robust_sigma
```

A histogram-based implementation is used for memory efficiency — peak memory scales
as `O(n_vmms × 1024 bytes)` regardless of the number of ROOT files loaded.

Noise quality flags:

| Flag | Condition |
|------|-----------|
| `ok`   | `robust_sigma < 10` and `n_noise ≥ 500` |
| `warn` | `10 ≤ robust_sigma < 13` |
| `bad`  | `robust_sigma ≥ 13` or `n_noise < 500` |

### 2. Signal Extraction and MPV Estimation (from `sng=0` run)

Signal hits are selected with:
- `over_threshold == 1`
- `adc < 1023` (removes saturated hits at the VMM3a 10-bit ceiling)

The **Most Probable Value (MPV)** of the Landau-like charge distribution is found
via a Gaussian-smoothed histogram peak finder:
- Search window: `[noise_cut, 800]` ADC
- Gaussian smoothing: `sigma = 2` bins
- Resolution: 150 bins at VMM level, 80 bins at channel level

### 3. SNR Metrics

Four metrics are computed at both VMM and channel level:

| Metric | Definition | Interpretation |
|--------|-----------|----------------|
| `snr` | MPV / noise_sigma | primary; robust peak position relative to noise width |
| `snr_mean` | mean_signal / noise_sigma | secondary cross-check using mean charge |
| `area_ratio` | signal counts / noise counts above `noise_cut` | how much signal dominates noise above threshold |
| `eer_value` | P at the crossing of P(noise>x) and P(signal≤x) | equal-error-rate; lower = better separation |

The **EER (equal-error-rate) threshold** is the ADC value where the fraction of
noise hits above it equals the fraction of signal hits below it — the optimal
discrimination threshold. In log scale the tails are nearly linear, making
extrapolation stable. The saturation probability `p_saturation` reports what
fraction of signal hits are lost to ADC=1023 clipping.

### 4. Channel-Level Quality Cuts

| Cut | Parameter | Default | Purpose |
|-----|-----------|---------|---------|
| Minimum noise hits | `min_noise_hits` | 50 | stable MAD estimate |
| Minimum sigma | `min_sigma` | 2.0 ADC | removes stuck channels (MAD floor) |
| Maximum sigma | `max_sigma` | 20.0 ADC | removes floating/noisy channels |
| MPV range | `mpv_min`–`mpv_max` | 100–300 ADC | removes peak-finder artifacts |
| Minimum signal hits | `min_signal_hits` | 50 | reliable MPV estimation |

### 5. Best Configuration Selection

`summarise_best_config()` compares all configurations per VMM:
- **VMM level**: highest `snr` (excluding bad noise quality)
- **Channel level**: highest median `snr_ch` across channels
- Prints a side-by-side table with an agreement flag (`✓` / `✗`)
- Reports an overall recommended configuration (plurality vote)

## Output Files

CSV files written to `cnfg_dir`, named with the dataset tag:

| File | Contents |
|------|----------|
| `vmm_snr_results_{dataset}.csv` | VMM-level SNR per configuration pair |
| `vmm_snr_per_channel_{dataset}.csv` | Channel-level SNR (quality-cut channels only) |
| `vmm_snr_per_channel_uncut_{dataset}.csv` | Channel-level SNR (all channels) |
| `vmm_snr_summary_{dataset}.csv` | Best configuration per VMM with agreement flag |
| `vmm_adc_analysis_{dataset}.csv` | Legacy ADC statistics |

Plots are saved as both PDF and PNG in `plot_dir`.

## Plots Reference

### SNR comparison plots

| Plot key | Function | Description |
|----------|----------|-------------|
| `snr_heatmap` | `plot_snr_heatmap` | SNR (MPV or mean) — rows=VMMs, cols=configs |
| `adc_heatmap` | `plot_adc_heatmap` | MPV / noise_sigma / mean_signal heatmaps |
| `snr_vs_peaking` | `plot_snr_vs_peaking` | SNR vs peaking time per VMM |
| `snr_vs_gain` | `plot_snr_vs_gain` | SNR vs gain per VMM |
| `snr_all_methods_heatmap` | `plot_all_methods_heatmap` | 4-panel: MPV-SNR, mean-SNR, area-ratio, EER side by side |
| `snr_method_comparison` | `plot_snr_method_comparison` | MPV-based vs mean-charge SNR scatter |

### Distribution-based plots

| Plot key | Function | Description |
|----------|----------|-------------|
| `tail_distributions` | `plot_tail_distributions` | Per VMM, log-scale: P(noise>x) dashed, P(signal≤x) solid, all configs overlaid; EER crossing marked with dot |
| `saturation_curves` | `plot_saturation_curves` | Per VMM, P(signal>x) for each config; red line at ADC=1022 shows saturation fraction |

### Channel-level plots

| Plot key | Function | Description |
|----------|----------|-------------|
| `snr_channel_heatmap_per_vmm` | `plot_snr_channel_heatmap_per_vmm` | Per-channel SNR heatmap per VMM |
| `snr_channel_heatmap_all_configs` | `plot_snr_channel_heatmap_all_configs` | Per-channel SNR across all configs |
| `snr_channel_uniformity` | `plot_snr_channel_uniformity` | Channel-to-channel SNR spread |

## Key Output Columns

### VMM-level (`vmm_snr_results_{dataset}.csv`)

| Column | Description |
|--------|-------------|
| `sg` | Shaping gain |
| `snt` | Shaping time (peaking time, ns) |
| `vmm_id` | VMM identifier |
| `noise_sigma` | Robust noise width (1.4826 × MAD) |
| `noise_cut` | ADC threshold separating noise from signal |
| `noise_quality` | `ok` / `warn` / `bad` |
| `mpv` | Most Probable Value of signal distribution (ADC) |
| `snr` | MPV / noise_sigma |
| `snr_mean` | mean_signal / noise_sigma |
| `area_ratio` | signal counts / noise counts above noise_cut |
| `eer_threshold` | ADC at the equal-error-rate crossing |
| `eer_value` | Probability at the EER crossing (lower = better) |
| `p_saturation` | Fraction of signal hits at ADC saturation (1023) |

### Channel-level (`vmm_snr_per_channel_{dataset}.csv`)

| Column | Description |
|--------|-------------|
| `vmm_id`, `ch` | VMM and channel identifiers |
| `noise_sigma_ch` | Channel-level robust noise width |
| `vmm_noise_cut` | VMM-level noise cut (lower bound for MPV search) |
| `mpv_ch` | Channel-level MPV (ADC) |
| `snr_ch` | Channel-level MPV / noise_sigma |
| `snr_mean_ch` | Channel-level mean_signal / noise_sigma |
| `area_ratio_ch` | Channel-level area ratio |
| `eer_threshold_ch` | Channel-level EER threshold |
| `eer_value_ch` | Channel-level EER value |
| `p_saturation_ch` | Channel-level saturation fraction |
| `n_signal`, `n_noise` | Hit counts used in the estimates |

## Testing Without Data

A standalone smoke test runs all new metrics and plot functions on synthetic
histograms — no ROOT files needed:

```bash
# Assertions only
python test_snr_methods.py

# Save plots to inspect visually
python test_snr_methods.py --save /tmp/snr_test_plots/
```

## Dependencies

```
numpy
pandas
scipy
uproot      (ROOT file I/O)
matplotlib
```
