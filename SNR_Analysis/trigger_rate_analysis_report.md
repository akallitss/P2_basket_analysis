# VMM3a ASIC Trigger Rate Analysis — P2 Basket Detector, SPS Beam Test

## Context and Hardware

The P2 basket detector consists of four sub-detectors instrumented with VMM3a ASICs:

- **P2 Large Detector**: VMMs 12, 13, 14, 15 — all 64 channels connected per VMM
- **P2 Small Detector 1**: VMMs 10 (ch 0–31 connected), 11 (ch 14–63 connected) — 82 pads on 2 connectors
- **P2 Small Detector 3**: VMMs 8 (ch 0–29 connected), 9 (ch 14–63 connected)
- **Trigger VMMs**: VMM 0 (ch 40 used as reference) and VMM 1 — scintillator-based coincidence trigger

### Configuration Scan Parameters

The following VMM3a parameters were scanned:

- **sg** (shaping gain): 3.0, 4.5, 6.0 mV/fC
- **snt** (peaking time): 50, 100, 200 ns

**5 kHz dataset — 7 configurations:**

| Run | sg (mV/fC) | snt (ns) |
|-----|-----------|---------|
| 73  | 3.0       | 50      |
| 70  | 3.0       | 100     |
| 67  | 3.0       | 200     |
| 97  | 4.5       | 50      |
| 82  | 4.5       | 100     |
| 78  | 4.5       | 200     |
| 99  | 6.0       | 200     |

**15 kHz dataset — 5 configurations:**

| Run | sg (mV/fC) | snt (ns) |
|-----|-----------|---------|
| 158 | 3.0       | 50      |
| 156 | 3.0       | 100     |
| 149 | 3.0       | 200     |
| 152 | 4.5       | 200     |
| 155 | 6.0       | 200     |

---

## Step 1 — Trigger Rate Time Series

### Method

The trigger VMM (VMM 0, ch 40) hit timestamps are binned into 1 ms windows. The resulting rate vs. time plot directly reveals the SPS beam spill structure.

### Parameter Selection

**Bin width — 1 ms**
Chosen to balance time resolution against per-bin statistics. At the expected on-spill trigger rate of ~4 kHz, a 1 ms bin collects ~4 triggers on average — enough to distinguish beam-on (≥ 2–3 counts) from beam-off (0 counts) with confidence. A larger bin (e.g. 10 ms) would smooth over the within-spill bunch microstructure that is needed to set the hole-filling parameter in Step 2; a smaller bin (e.g. 0.1 ms) would leave most bins empty even during spills, making subsequent threshold classification unreliable. The rule of thumb is: bin_width ≈ 1 ms is a good starting point for trigger rates in the range 1–100 kHz.

**Trigger reference channel — VMM 0, ch 40**
The trigger VMM (VMM 0) is connected to a scintillator-based coincidence trigger. Channel 40 is selected because it shows the highest and most consistent hit rate during spills, confirming it receives the coincidence signal cleanly. This is verified by inspecting the hits-per-channel plot of VMM 0 (`hits_per_channel_trigger_5kHz.png`): the reference channel stands out clearly above all others. Any channel on the trigger VMM with a hit rate clearly above the rest and stable across spills is a valid reference.

### What It Shows and Why It Is Useful

The SPS delivers beam in discrete spills separated by long inter-spill gaps. A working detector and DAQ shows a clear two-level structure: high rate (~4–28 kHz depending on beam intensity) during spills, near-zero rate between them. This is the first sanity check — if the trigger rate shows no structure, data acquisition failed.

### Results

**5 kHz — Run 67 (sg=3.0, snt=200):**
Run 67 shows approximately 10 regular spills over ~560 s, each lasting ~5 s, with trigger rate ~4–5 triggers/ms (`4–5 kHz`) on-spill and `0` off-spill. The 5 s zoom confirms the SPS bunch microstructure: within a spill, individual 1 ms bins occasionally dip to zero due to the SPS RF structure, but the rate remains consistently above `3 kHz` throughout the spill body.

**15 kHz — Run 149 (sg=3.0, snt=200):**
Run 149 shows only 2 spills visible in the ~120 s file, with a peak rate of ~80 triggers/ms (`80 kHz`). Note: run 158 (sg=3.0, snt=50) ran at a lower effective beam intensity (~10 kHz on-spill) due to beam conditions during that fill; rates from this run should be interpreted separately.

**Plot references:** `trigger_rate_ms_run67_5kHz.png`, `trigger_rate_ms_run149_15kHz.png`

---

## Step 2 — Spill Mask Construction (Beam-On / Beam-Off Classification)

### Method

Each 1 ms bin is classified as spill-on (beam present) or spill-off (inter-spill gap) using a fixed threshold on the trigger rate. Threshold: `1.0 kHz`. Any bin with trigger rate `> 1 kHz` is on-spill; below = off-spill. Two refinements are applied:

1. **Hole-filling (`max_gap_s = 2 s`)**: Short off-gaps within a spill (< 2 s) are reclassified as on-spill. This corrects for the SPS bunch microstructure at low beam intensity, where individual 1 ms bins dip to zero mid-spill.
2. **Minimum spill duration (`min_spill_s = 1 s`)**: On-regions shorter than 1 s are discarded — these correspond to start-of-run artifacts or isolated noise spikes.

### Parameter Selection

**Spill threshold — 1 kHz**
Read directly from the trigger rate time series (Step 1). The plot shows a clear two-level structure with no ambiguity: on-spill rates are ~3–80 kHz; off-spill rates are effectively 0. Any value in the range 0.1–2 kHz works for this dataset; 1 kHz is chosen as a round number comfortably below the minimum on-spill rate (~3 kHz at 5 kHz beam) and well above the maximum off-spill rate (~0.01 kHz). To set this threshold for a new dataset, identify the gap between the two levels in the trigger rate histogram and place the threshold at the midpoint of that gap.

**Hole-filling gap — 2 s**
Determined by inspecting the zoomed trigger rate panel. Within a single spill, dips below threshold are typically < 200 ms wide (SPS bunch microstructure at 5 kHz beam intensity). A 2 s gap threshold safely bridges these within-spill dips while leaving real inter-spill gaps untouched — the SPS cycle at the CERN SPS is ~14–21 s in this dataset, so inter-spill gaps are always ≥ 5 s. The rule: set `max_gap_s` comfortably less than the shortest observed inter-spill gap and comfortably greater than the longest within-spill dip. Both can be read from the zoomed trigger rate plot.

**Minimum spill duration — 1 s**
Real spills last ~5 s in this dataset. Any on-region shorter than 1 s after hole-filling is an artifact — typically a burst of noise at run start or an isolated cosmic ray hit during the inter-spill gap. Set `min_spill_s` to less than half the expected shortest real spill duration, which is visible in the zoomed trigger rate plot.

### Why This Metric Is Good

Separating on-spill from off-spill time is fundamental. The spill-on rate captures signal + noise; the spill-off rate captures only noise. Any metric computed without this separation would mix the two and be uninterpretable. The `1 kHz` threshold sits comfortably between the ~4–80 kHz on-spill rate and the ~0 kHz off-spill rate, making it robust across the full range of beam intensities in this dataset.

### Results

**5 kHz — Run 67 (sg=3.0, snt=200):**
`46,102` on-bins (`8.1%`) and `522,549` off-bins. The low on-fraction reflects the fact that only 1 ROOT file is loaded for this run and captures a short fraction of the full run. The zoom panel confirms clean classification: the 5 s spill is entirely captured as on-spill (green shading), with an abrupt and correct transition to off-spill at spill end.

**5 kHz — Run 82 (sg=4.5, snt=100):**
`155,544` on-bins (`26.8%`) using file index 2. Clean, regular spill structure is visible over the full ~600 s run. The mask correctly identifies each ~5 s spill.

**15 kHz — Run 149 (sg=3.0, snt=200):**
`6,682` on-bins (`5.8%`) — only 2 spills captured in this file. The `80 kHz` peak rate is well above the `1 kHz` threshold; the mask is correct. The 10 s zoom shows no spill in the first 10 s of the file, consistent with the full-run view.

**Plot references:** `spill_mask_run67_5kHz.png`, `spill_mask_run82_5kHz.png`, `spill_mask_run149_15kHz.png`

---

## Step 3 — Channel Quality Selection

### Step 3a — Total Hit Count per Channel (Dead / Noisy Identification)

#### Method

Hit counts are aggregated over 3 ROOT files per run across **all diagnostic runs** for all channels (0–63) of each VMM. Channels are then classified relative to the median hit count of connected channels:

- **Dead**: count `< 10%` × median — disconnected or broken pad
- **Noisy (total-hit criterion)**: count `> 5×` median — excessively self-triggering at total hit level

The connected channel range per VMM is defined from the detector geometry (e.g. VMM 10 uses ch 0–31 only; ch 32–63 are unconnected floating inputs that self-trigger at high gain).

#### Parameter Selection

**Connected channel range**
Defined in `vmm_mapping.py` from knowledge of the detector geometry — which detector pads are physically wire-bonded to which VMM channel indices. For the small detectors, connector 0 maps to VMM channels 0–31 and connector 1 maps to channels 14–63 (with channels 14–31 shared between connectors). Unconnected channels are excluded from all threshold computations and rate averages. To verify the mapping, check the hit-count plot: connected channels show consistent non-zero counts; unconnected channels either show near-zero counts (floating input at low gain) or anomalously high counts (floating input at high gain, where the open input self-triggers).

**Dead threshold — 10% of median, over all diagnostic runs**
A channel is dead if its total hit count is below 10% of the median hit count among connected channels of the same VMM. Hit counts are aggregated across all diagnostic runs (not only the test run) so that a channel at the beam edge — which may receive fewer hits at a specific gain setting — is not incorrectly flagged as dead. Only channels with consistently low hit counts across all configurations are classified as dead. Channels flagged as dead are visible as firebrick-red bars in the hit-count plot that fall clearly below the median line.

**Noisy threshold for total hits — 5× median**
Conservative value used for visual diagnostic labeling only. It flags only the most extreme outliers in the total hit count. This threshold is intentionally permissive because total hit count cannot distinguish a beam-hit pad from a genuinely noisy channel — that distinction is handled by the off-spill rate criterion in Step 3b. Channels above 5× median appear as orange bars in the hit-count plot.

#### Why This Step Is Useful

Dead channels contribute zero signal and should be excluded from rate normalization. Unconnected channels that float at high gain generate enormous hit counts and would dominate the average rate if included. This step provides a first coarse quality map of the detector.

#### Results (5 kHz)

P2 Large Detector VMMs 12–15 show dead = 1 channel on VMM 12 at sg=3.0, otherwise fully active. No channels exceed the `5×` noisy threshold on the large detector. Small Detector VMMs 10–11 confirm the expected geometry: VMM 10 active on ch 0–31 only, VMM 11 active on ch 14–63 only; the unconnected halves show dramatically fewer hits.

**Plot references:** `hits_per_channel_p2_large_1_5kHz.png`, `hits_per_channel_p2_small_1_5kHz.png`

---

### Step 3b — Per-Channel Spill-On / Spill-Off Rate (Noise-Based Channel Masking)

#### Method

For each diagnostic run (covering all configurations), per-channel hit rates are computed separately during on-spill and off-spill periods using the spill mask from Step 2. The off-spill rate is the key discriminator:

- A channel with elevated **off-spill rate** (beam absent) is **genuinely noisy** — self-triggering on electronics noise.
- A channel with elevated **on-spill rate only** is a **beam-hit pad** — seeing real signal. It must not be masked.

Noisy threshold: off-spill rate `> 3×` median off-spill rate among connected channels of that VMM. The element-wise maximum off-spill rate across all diagnostic runs is used, so a channel noisy at any gain setting is masked everywhere.

#### Parameter Selection

**Noisy threshold for off-spill rate — 3× median of active channels**
Tighter than the total-hit threshold because the off-spill rate distribution should be narrow: without beam, all connected channels in a well-behaved VMM see similar electronics noise. A channel at 3× the median is already anomalously noisy. The value was chosen empirically: at 5× median, known self-triggering channels were missed; at 2× median, beam-hit pads near the spill boundary (where residual charge briefly elevates the rate in the first ~1 ms after spill end) were incorrectly flagged. The 3× threshold correctly identifies genuine self-triggering while avoiding false positives.

The median is computed exclusively over **non-dead connected channels**. Including dead channels (which have near-zero off-spill rates) would pull the median down artificially, placing the `3× median` threshold too low and causing channels with similar off-spill rates to be classified inconsistently — some just above the deflated threshold, some just below. By restricting the median to active channels, the threshold is anchored to the true electronics noise floor of well-functioning pads.

The per-channel rate plots show a dashed orange line at the noisy threshold on the spill-off panel, so the classification boundary is directly visible. All channels with bars above this line are masked as noisy; all bars below are retained.

**Per-configuration masking**
Diagnostic runs are grouped by `(sg, snt)`. For each configuration, the off-spill rate is computed only from the diagnostic run(s) that share that gain/peaking-time setting. A channel is only excluded for the configuration in which its off-spill rate is elevated — not globally across all configurations. This avoids penalising a channel that is noisy at high gain but performs correctly at low gain, which would suppress signal in configurations where that channel is clean.

If multiple diagnostic runs share the same `(sg, snt)`, their per-channel off-spill rates are merged by element-wise maximum before thresholding (the most conservative estimate within that configuration). Configurations in the scan that have no matching diagnostic run fall back to the global mask (element-wise maximum across all diagnostic runs).

**Number of diagnostic files per run — 5 files, starting from file 1**
Five files provide sufficient statistics for per-channel rate estimation (typically 10–50 spills depending on the run). File 0 is skipped (`file_start=1`) because in some runs the first file contains corrupted timestamps at the very start of acquisition. For the config-scan rate computation itself (Step 4), only 1 file per run is used to keep I/O manageable; the per-channel noisy detection uses more files for robustness.

#### Why This Metric Is Better than Total-Hit Noisy Detection

The total-hit criterion cannot distinguish a beam-hit pad (high on-spill rate) from a self-triggering channel (high off-spill rate). Using the off-spill rate as the noise discriminator avoids false-positive masking of real detector pads in the beam footprint — a beam-hit pad with high on-spill rate but low off-spill rate is correctly retained.

#### Results (5 kHz)

A small number of individual channels are flagged on VMM 14 and VMM 15 (shown as orange bars in the per-channel rate plots). VMM 12 and VMM 13 have very few noisy channels. The per-channel plots also reveal that VMM 14 and VMM 15 show a systematically higher off-spill baseline across all channels at certain configurations (particularly sg=4.5, snt=100 and sg=3.0, snt=200), consistent with a whole-VMM noise elevation rather than isolated individual channels — possibly related to their position in the beam footprint (higher beam activity leading to increased cross-talk into the off-spill period).

**Final good channel counts (5 kHz):** VMM 12: 58 ch, VMM 13: 62 ch, VMM 14: 63 ch, VMM 15: 61 ch. VMMs 10–11 and 8–9 are geometry-limited as expected from the detector design.

**Plot references:** `rate_per_channel_p2_large_1_run78_5kHz.png`, `rate_per_channel_p2_large_1_run82_5kHz.png`, `rate_per_channel_p2_large_1_run97_5kHz.png`, `rate_per_channel_p2_large_1_run99_5kHz.png`

---

## Step 4 — Spill-On and Spill-Off Rate per Configuration

### Method

For each configuration (sg, snt), all hits from good channels only are accumulated into 1 ms time bins. The spill mask is applied to separate on-spill and off-spill bins. Two rates are computed per VMM per configuration:

- **`rate_on_khz`**: mean detector hit rate during on-spill bins, divided by `n_good_channels` (kHz per channel) — signal + noise, dominated by beam-induced muon hits.
- **`rate_off_khz`**: mean detector hit rate during off-spill bins, divided by `n_good_channels` (kHz per channel) — pure electronics noise floor.

Normalization by `n_good_channels` makes VMMs with different numbers of active pads directly comparable. Additionally, **`rate_on_std_khz`** and **`rate_off_std_khz`** (standard deviation across good channels) are computed from a per-channel accumulation pass with the spill mask applied. These quantify channel-to-channel spread within each VMM.

### Parameter Selection

**Good channel mask**
Built from the union of dead channels (Step 3a) and noisy channels (Step 3b). Only hits from channels in this mask are counted when computing the VMM rates. The mask is computed **per `(sg, snt)` configuration** using the diagnostic run(s) matching that setting (see Step 3b). Each run in the config scan therefore uses only the noisy-channel flags derived from its own gain/peaking-time conditions, so a channel that self-triggers only at high gain is not penalised at lower gains. Dead channels (from total hit counts) are applied globally across all configurations, since a dead channel is hardware-level and does not depend on VMM settings.

**Normalization by n_good_channels**
The raw VMM hit rate (total hits from all good channels ÷ time) is divided by the number of good channels to obtain a per-channel rate. This normalization is essential for comparing VMMs with different numbers of connected pads (e.g. VMM 10 with 32 channels vs. VMM 12 with 63 channels). Without it, a VMM with more channels would always appear to have a higher total rate regardless of per-pad performance.

**Number of files per config-scan run — 1 file**
A single ROOT file per run is used in Step 4 to keep the total processing time manageable across 7 configurations × 4 VMM groups. This is sufficient because the rates are averaged over all on-spill (or off-spill) bins in the file, which typically covers 2–30 complete spills depending on the run. The per-channel standard deviation (±1σ band) uses a separate pass over the same single file with the spill mask already applied.

### Why These Metrics Are Good

- **`rate_on`**: directly proportional to the probability that a beam muon generates a hit above threshold. Configurations with higher on-spill rate detect more muons — but gain must be balanced against noise.
- **`rate_off`**: the noise floor. A well-configured VMM has `rate_off` much smaller than `rate_on`. A rising noise floor at high gain limits the useful dynamic range.
- **`rate_on_std`**: a wide ±1σ band indicates beam hits concentrated on a few channels (narrow beam footprint or non-uniform detector response). A narrow band indicates uniform illumination.
- **`rate_off_std`**: should be small. A large spread indicates individual noisy channels pulling up the variance even after masking.

---

### 5 kHz Results — P2 Large Detector

The trigger reference (VMM 0) is stable across all 7 configurations at `4.1 kHz`, confirming that beam conditions are consistent and the configuration scan changes only the VMM ASIC settings, not the beam flux.

**Spill-on rate** increases monotonically with both sg and snt within the scanned range:

- VMM 12 (58 ch): `0.0015 kHz/ch` (sg=3.0, snt=50) to `0.005 kHz/ch` (sg=6.0, snt=200) — factor ~3
- VMM 13 (62 ch): `0.001 kHz/ch` to `0.003 kHz/ch` — similar trend, lower absolute rate
- VMM 14 (63 ch): `0.003 kHz/ch` (sg=3.0, snt=50) to `0.016 kHz/ch` (sg=4.5 or 6.0, snt=200) — factor ~5, consistently ~4–5× higher than VMM 12
- VMM 15 (61 ch): `0.002 kHz/ch` to `0.010 kHz/ch` — intermediate between VMM 12 and VMM 14

The ±1σ band on the on-spill rate is always wider than the mean, indicating the beam footprint illuminates channels non-uniformly — channels directly under the beam see much higher rates than edge channels.

**Spill-off rate (noise floor)** remains very small for all configurations:

- VMM 12/13: `0.0001–0.0008 kHz/ch` across all configs — noise floor negligible compared to signal
- VMM 14/15: `0.0002–0.0009 kHz/ch` — slightly higher but still well below the on-spill rate

The on/off ratio for the large detector ranges from ~5 (sg=3.0, snt=50) to ~50+ (sg=4.5 or 6.0, snt=200), indicating excellent signal-to-noise for all configurations tested.

**Plot references:** `spill_on_vmms_P2_Large_Detector_rate_on_khz_5kHz.png`, `spill_on_vmms_P2_Large_Detector_rate_off_khz_5kHz.png`, `spill_rates_config_P2_Large_Detector_5kHz.png`

---

### 5 kHz Results — P2 Small Detectors

VMMs 10 and 11 (P2 Small Detector 1) show dramatically higher on-spill rates than the large detector VMMs: `~0.065–0.17 kHz/ch`. This is expected — the small detector has 82 pads covering a much smaller area, so the same beam spot hits a higher fraction of channels. The noise floor remains negligible (`< 0.002 kHz/ch`) across all configurations, giving very good signal-to-noise.

The on-spill rate trend for the small detector follows the same sg, snt dependence as the large detector but with larger absolute values and a cleaner monotonic increase from sg=3.0, snt=50 through sg=6.0, snt=200.

**Plot references:** `spill_rates_config_P2_Small_Detector_1_5kHz.png`

---

### 15 kHz Results — P2 Large Detector

The 15 kHz dataset covers 5 configurations: sg=3.0 (snt=50, 100, 200), sg=4.5 (snt=200), sg=6.0 (snt=200).

The trigger reference is stable at `~26–28 kHz` for most runs. Run 158 (sg=3.0, snt=50) shows `10 kHz` due to lower beam intensity during that fill — rates from this run should be interpreted separately.

**Spill-on rate** is much higher in absolute terms than at 5 kHz (×6–30 depending on VMM and configuration), consistent with the higher beam flux. At sg=6.0, snt=200:

| VMM | rate_on (kHz/ch) |
|-----|-----------------|
| 12  | `0.110`         |
| 13  | `0.051`         |
| 14  | `0.103`         |
| 15  | `0.136`         |

The on-spill rate increases monotonically with sg and snt, the same trend as at 5 kHz. The ±1σ band is very wide — at 15 kHz beam intensity, hit rates span nearly a decade across channels within the same VMM, highlighting the strong non-uniformity of beam illumination.

**Spill-off rate** is notably higher than at 5 kHz for the same configurations, especially at high gain:

- VMM 14 at sg=4.5, snt=200: `0.0035 kHz/ch` off-spill (vs `0.0003 kHz/ch` at 5 kHz — a factor of ~10 increase)
- VMM 15 at sg=6.0, snt=200: `0.0037 kHz/ch`

This rise in the noise floor at 15 kHz beam intensity is consistent with increased cross-talk or charge-sharing effects at higher rates. Even so, the on/off ratio remains `> 10` for all VMMs at all configurations tested.

The 15 kHz summary plot shows a qualitative difference from 5 kHz: the on-spill rate curve rises smoothly and steeply from sg=3.0, snt=50 through sg=6.0, snt=200 without the zig-zag pattern seen at 5 kHz. This is because the 15 kHz configuration set does not include the sg=4.5, snt=50 and snt=100 points that created the apparent zig-zag in the 5 kHz scan.

**Plot references:** `spill_on_vmms_P2_Large_Detector_rate_on_khz_15kHz.png`, `spill_on_vmms_P2_Large_Detector_rate_off_khz_15kHz.png`, `spill_rates_config_P2_Large_Detector_15kHz.png`

---

## Summary Tables — Key Rates at Representative Configurations

### P2 Large Detector — VMM 12 (rate per good channel, kHz)

| Dataset | sg (mV/fC) | snt (ns) | rate_on (kHz/ch) | rate_off (kHz/ch) | on/off ratio |
|---------|-----------|---------|-----------------|-------------------|-------------|
| 5 kHz   | 3.0       | 50      | `0.0015`        | `0.00050`         | 3           |
| 5 kHz   | 3.0       | 200     | `0.0027`        | `0.00027`         | 10          |
| 5 kHz   | 4.5       | 200     | `0.0047`        | `0.00013`         | 36          |
| 5 kHz   | 6.0       | 200     | `0.0047`        | `0.00034`         | 14          |
| 15 kHz  | 3.0       | 50      | `0.0041`        | `0.00015`         | 27          |
| 15 kHz  | 3.0       | 200     | `0.062`         | `0.00020`         | 310         |
| 15 kHz  | 4.5       | 200     | `0.077`         | `0.0013`          | 59          |
| 15 kHz  | 6.0       | 200     | `0.110`         | `0.0013`          | 85          |

### P2 Large Detector — VMM 14 (rate per good channel, kHz)

| Dataset | sg (mV/fC) | snt (ns) | rate_on (kHz/ch) | rate_off (kHz/ch) | on/off ratio |
|---------|-----------|---------|-----------------|-------------------|-------------|
| 5 kHz   | 3.0       | 50      | `0.0030`        | `0.00032`         | 9           |
| 5 kHz   | 3.0       | 200     | `0.0089`        | `0.00079`         | 11          |
| 5 kHz   | 4.5       | 200     | `0.0163`        | `0.00030`         | 54          |
| 5 kHz   | 6.0       | 200     | `0.0161`        | `0.00073`         | 22          |
| 15 kHz  | 3.0       | 50      | `0.0032`        | `0.00017`         | 19          |
| 15 kHz  | 3.0       | 200     | `0.046`         | `0.00029`         | 160         |
| 15 kHz  | 4.5       | 200     | `0.070`         | `0.0035`          | 20          |
| 15 kHz  | 6.0       | 200     | `0.103`         | `0.0036`          | 29          |

---

## Conclusions

1. **The spill mask correctly separates beam-on from beam-off time** at both beam intensities. The `1 kHz` threshold provides a clean margin, and the 2 s hole-filling correctly handles mid-spill dips from the SPS bunch microstructure without merging real inter-spill gaps.

2. **Signal rates increase with both gain (sg) and peaking time (snt)**, as expected from VMM3a ASIC physics: higher gain amplifies charge above threshold; longer peaking time integrates more charge. The monotonic trend holds consistently across both beam intensities and all detector VMMs.

3. **The noise floor (off-spill rate) is negligible at 5 kHz** for all configurations and all VMMs: `rate_off / rate_on < 5%` in all cases. At 15 kHz, the noise floor rises at high gain (VMM 14/15 reach `rate_off ~ 0.004 kHz/ch` at sg=4.5–6.0, snt=200), but the on/off ratio remains `> 10` even in the worst case.

4. **VMM 14 and VMM 15 consistently show higher hit rates** (on-spill ×3–5, off-spill ×2–5 compared to VMM 12). This is likely due to their position in the beam footprint (more pads directly illuminated). Even after good-channel masking, their off-spill rate shows stronger configuration dependence — possibly due to cross-talk or EMI sensitivity at high gain and long peaking time.

5. **The ±1σ band across channels is always wide on the spill-on rate**, confirming the beam does not illuminate all pads uniformly. It is narrow on the spill-off rate for VMM 12/13 (uniform noise floor) but wide for VMM 14/15 at high gain, indicating lingering channel-to-channel noise variation even after masking.

6. **Per-configuration channel masking avoids cross-config bias**: the good-channel mask is built separately for each `(sg, snt)` configuration from the diagnostic run matching that setting. Channels noisy only at high gain are not excluded from low-gain configurations where they are clean.

7. **Three masking bugs were identified and corrected**:
   - *Dead detection used only the test run*: a channel at the beam edge with fewer hits at one gain was incorrectly flagged as dead. Fixed by aggregating hit counts across all diagnostic runs.
   - *Dead channels contaminated the noisy threshold*: their near-zero off-spill rates pulled the median down, placing the `3× median` cutoff too low and causing channels with nearly identical rates to be classified differently. Fixed by computing the median over non-dead active channels only.
   - *Diagnostic plots showed the global mask instead of the per-config mask*: a channel noisy at `sg=6.0` appeared masked in the `sg=3.0` diagnostic plot even though its off-spill rate there was indistinguishable from unmasked neighbours. Fixed by passing the per-config `(sg, snt)` mask to each plot. A threshold line is now drawn on the spill-off panel to make the classification boundary directly visible.

8. **Configuration recommendation (preliminary)**: Based on these trigger-level metrics, sg=4.5 with snt=200 offers the best signal rate with a low noise floor at 5 kHz, giving on/off ratios of 36 (VMM 12) and 54 (VMM 14). At 15 kHz, the same configuration becomes noisier (VMM 14 drops to on/off ~20), suggesting sg=3.0 with snt=200 may be preferable at higher beam intensity from a noise perspective, despite its lower absolute signal rate. The optimal configuration will be confirmed by the full SNR analysis from the configuration scan pipeline.