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
vmm_mapping.py               — detector geometry and VMM groupings
vmm_qa.py                    — quality assurance investigations
```

## Detector Layout

| Detector | VMM IDs | Channels |
|----------|---------|----------|
| P2 Large Detector | 12, 13, 14, 15 | 0–63 per VMM |
| P2 Small Detector 1 | 10 (conn. 0), 11 (conn. 1) | ch 0–31, ch 14–63 |
| P2 Small Detector 3 | 8 (conn. 0), 9 (conn. 1) | ch 0–31, ch 14–63 |
| Trigger (excluded) | 0, 1 | — |

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

The factor 1.4826 makes MAD a consistent estimator of the Gaussian sigma. A histogram-based implementation is used for memory efficiency — peak memory scales as `O(n_vmms × 1024 bytes)` regardless of the number of ROOT files loaded.

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

The **Most Probable Value (MPV)** of the Landau-like charge distribution is found via a Gaussian-smoothed histogram peak finder:
- Search window: `[noise_cut, 800]` ADC (lower bound anchored to the per-VMM noise_cut)
- Gaussian smoothing: `sigma = 2` bins
- Resolution: 150 bins at VMM level, 80 bins at channel level (fewer statistics)

Two SNR estimators are computed:

| Estimator | Definition |
|-----------|------------|
| `snr`      | `MPV / noise_sigma` (primary) |
| `snr_mean` | `mean_signal / noise_sigma` (secondary cross-check) |

### 3. Channel-Level Quality Cuts

In addition to the VMM-level selection, each channel must pass:

| Cut | Parameter | Default | Purpose |
|-----|-----------|---------|---------|
| Minimum noise hits | `min_noise_hits` | 50 | stable MAD estimate |
| Minimum sigma | `min_sigma` | 2.0 ADC | removes stuck channels (MAD floor) |
| Maximum sigma | `max_sigma` | 20.0 ADC | removes floating/noisy channels |
| MPV range | `mpv_min`–`mpv_max` | 100–300 ADC | removes peak-finder artifacts |
| Minimum signal hits | `min_signal_hits` | 50 | reliable MPV estimation |

### 4. Best Configuration Selection

`summarise_best_config()` compares all configurations per VMM:
- **VMM level**: selects the configuration with highest `snr` (excluding bad noise quality)
- **Channel level**: selects the configuration with highest median `snr_ch` across channels
- Prints a side-by-side comparison table with an agreement flag (`✓` / `✗`)
- Reports an overall recommended configuration (plurality vote across VMMs)

## Running the Analysis

Edit the **User Configuration** block at the top of `vmm_config_scan_analysis.py`:

```python
# Paths
cnfg_dir = "/path/to/config/and/output/"
data_dir = "/path/to/run_*/data/"

# Number of ROOT files to merge per run (increase for more statistics)
n_root_files = 5

# Run numbers for QA investigations
qa_run_signal = 149   # well-behaved sng=1 run at preferred config
qa_run_noisy  = 150   # most problematic run (shortest peaking time)
qa_runs_quality_check = [149, 150]

# Plot output
plot_dir   = f"{cnfg_dir}plots/"
show_plots = False    # True to display figures interactively
```

Enable/disable individual plots and QA checks via the `PLOTS` dictionary:

```python
PLOTS = {
    "snr_heatmap"               : True,   # SNR vs (sg, snt) heatmap
    "snr_channel_heatmap_per_vmm": True,  # per-channel SNR heatmap
    "qa_noise_run_diagnostic"   : True,   # noise run sanity check
    ...
}
```

Run:

```bash
python vmm_config_scan_analysis.py
```

The run table CSV (e.g. `vmm_config_scan_15kHz.csv`) must exist in `cnfg_dir`. It maps run numbers to their VMM settings (`sg`, `snt`, `sng`).

## Output Files

CSV files written to the working directory:

| File | Contents |
|------|----------|
| `vmm_snr_results.csv` | VMM-level SNR per configuration pair |
| `vmm_snr_per_channel.csv` | Channel-level SNR (quality-cut channels only) |
| `vmm_snr_per_channel_uncut.csv` | Channel-level SNR (all channels, no quality cuts) |
| `vmm_snr_summary.csv` | Best configuration per VMM with agreement flag |
| `vmm_adc_analysis.csv` | Legacy ADC statistics (mean, RMS, robust sigma) |

Plots are saved as both PDF and PNG in `plot_dir`.

## Key Output Columns

### VMM-level (`vmm_snr_results.csv`)

| Column | Description |
|--------|-------------|
| `sg` | Shaping gain setting |
| `snt` | Shaping time (peaking time, ns) |
| `vmm_id` | VMM identifier |
| `noise_sigma` | Robust noise width (1.4826 × MAD) |
| `noise_cut` | ADC threshold separating noise from signal |
| `noise_quality` | `ok` / `warn` / `bad` |
| `mpv` | Most Probable Value of signal distribution (ADC) |
| `snr` | MPV / noise_sigma |
| `snr_mean` | mean_signal / noise_sigma |

### Channel-level (`vmm_snr_per_channel.csv`)

| Column | Description |
|--------|-------------|
| `vmm_id`, `ch` | VMM and channel identifiers |
| `noise_sigma_ch` | Channel-level robust noise width |
| `vmm_noise_cut` | VMM-level noise cut (used as MPV search lower bound) |
| `mpv_ch` | Channel-level MPV (ADC) |
| `snr_ch` | Channel-level MPV / noise_sigma |
| `snr_mean_ch` | Channel-level mean_signal / noise_sigma |
| `n_signal`, `n_noise` | Hit counts used in the estimates |

## Dependencies

```
numpy
pandas
scipy
uproot      (ROOT file I/O)
matplotlib
```
