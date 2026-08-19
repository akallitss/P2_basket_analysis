# VMM timing: what exists, what it says, and what would improve it

## Yes, it exists — and it ran on the whole campaign

The measurement remembered from the beam is `vmm_efficiency.py`
(P2_basket_online_analysis), and it is exactly as described:

* **trigger-referenced** — the external trigger is digitised on hybrid 0, VMM 0
  channel 44, in the *same* stream as the detectors, so no external reference
  is needed;
* **BCID-phase coincidence inside one SRS marker** — the firmware leaves the
  `offset` (BCID rollover) counter stuck, and markers arrive every 1.6384 ms
  while BCID wraps every 92.16 µs, so absolute time is ambiguous by an unknown
  multiple of 92.16 µs and is not even monotonic. Hits are therefore paired
  only with triggers carrying the **same `srs_timestamp`**, and compared by
  **BCID phase wrapped into ±46.08 µs**. This lifts the coincidence peak from
  1.01× the accidental background to 8.6–17.1×;
* **accidental-subtracted** — a Gaussian peak on a flat background gives the
  signal window, and an equal-width sideband offset past the peak is
  subtracted. ~17.8 BCID cycles fit in one marker interval, so the true partner
  is one of ~18 candidates; that is what the sideband removes;
* **efficiency relative to the trigger acceptance** — the scintillator defines
  the denominator, uncorrected for its geometric overlap with each station, so
  the absolute value is a lower bound and only station-to-station comparison is
  meaningful.

**It does not need reproducing: it already ran on every capture.** Each
capture's `scalars.json` carries `mu_ns`, `sigma_ns`, `contrast`, the raw /
accidental / corrected efficiency and the binomial interval, per station.

`vmm_trigger_timing.py` gathers all of it and joins each capture to its run's
own `run_config.json`:

* **10 122 station-captures over 45 runs**, 9 321 (92 %) past the fit-quality
  cut (`contrast > 5`, |µ| < 600 ns — a failed fit is unmistakable, µ runs to
  the edge of the search window and the peak sits at the accidental floor).
* Products: `vmm_trigger_timing.csv` (per capture), `vmm_trigger_timing_by_hv.csv`
  (median per station × gas × mesh × drift).

## Different gases: already measured, on all three stations

Median σ [ns] at mesh 450, by drift voltage — **this is the gas comparison,
and it is already in hand**:

| drift [V] | P2_MID gas A | P2_MID gas B | P2_OUT gas A | P2_OUT gas B |
|---|---|---|---|---|
| 600 | 37.6 | **33.1** | 90.2 | **69.5** |
| 650 | 26.5 | **24.6** | 50.6 | **44.4** |
| 700 | 23.0 | – | 33.4 | – |
| 750 | 22.1 | **19.4** | 26.4 | **25.7** |
| 800 | 21.9 | – | 24.1 | – |
| 850 | 22.3 | – | 23.4 | – |

gas A = Ar/CO₂/iC₄H₁₀ 93/5/2, gas B = Ar/CF₄/iC₄H₁₀ 88/10/2.

* **Gas B is better at every point where both were taken**, and the margin is
  largest at low field (P2_OUT 90.2 → 69.5 at drift 600) — which is where a
  faster drift velocity should show, and is the shape the Magboltz study
  predicted.
* **P2_MID on gas B reaches 19.4 ns — under the 20 ns P2 goal**, measured on
  the VMM readout.
* **σ plateaus above drift 750** (21.9 / 22.3 at 800 / 850 on gas A), so the
  working point is on the flat.
* Mean latency drops on gas B at every drift point (P2_MID 171.2 → 157.4 ns at
  600, 108.4 → 103.4 at 750) — the drift is measurably faster, as designed.
* **Gas B is limited to 3 drift points** (600, 650, 750) by the known holes:
  run_59 recorded nothing, run_62 stopped after 2 of 6.
* **P2_IN is excluded** — raised thresholds, few captures, unstable fits.

Figure: `mpgd2026/figs/vmm_timing.png` panel (a).

## Can it be improved? Yes — but not the way the deck says

The roadmap, the report and the draft deck all say σ is *"pinned at the 22.5 ns
BCID quantisation for want of time-calibration runs"*. **That is wrong on the
mechanism**, and `vmm_timing_budget.py` is the measurement.

**1. The TDC fine time is already applied.** `vmm_decode.derive()` computes
`t_fine = 22.5 − tdc × 60/255` ns with the nominal TAC slope, so the time is
*not* quantised at 22.5 ns. And the observed TDC spread (1st–99th percentile,
100 counts ≈ 23.5 ns) maps onto almost exactly one clock period, so the nominal
slope is about right. Toggling the fine time off changes nothing measurable:

| | BCID only | BCID + TDC |
|---|---|---|
| P2_MID − trigger | 43.8 ns | 42.7 ns |
| P2_OUT − trigger | 63.6 ns | 62.9 ns |

which is exactly what should happen: **BCID quantisation alone is
22.5/√12 = 6.5 ns rms**, and 6.5 ns in quadrature with 43 ns is invisible.
The clock was never the limit.

**2. Per-channel t0 and time walk are real but subdominant.** Both corrections
are free offline (the VMM ships a per-hit ADC), and both were measured:

| station | measured | after per-channel t0 | after t0 + slewing |
|---|---|---|---|
| P2_MID | 42.7 | 41.6 | 41.4 ns |
| P2_OUT | 62.9 | 60.4 | 60.3 ns |
| P2_IN | 69.2 | 67.7 | 67.9 ns |

The channel-to-channel t0 spread is 12.9–16.7 ns rms (range up to 104 ns) and
the walk is ~19 ns across the ADC deciles — genuine, but 13 ns in quadrature
with 43 ns buys 2–4 %. Worth doing, not worth a campaign.

**3. The trigger channel is the largest removable term.** Station-minus-station
coincidence contains no trigger term at all, which over-determines the system:

```
sigma(P2_MID - TRIG)  = 42.7 ns        sMID^2 + sTRIG^2
sigma(P2_OUT - TRIG)  = 62.9 ns        sOUT^2 + sTRIG^2
sigma(P2_MID - P2_OUT)= 59.5 ns        sMID^2 + sOUT^2
-->  trigger channel 33.4 ns | P2_MID intrinsic 26.6 | P2_OUT intrinsic 53.2
```

**The trigger is a scintillator read through a VMM discriminator** — its own
threshold, its own walk, no calibration — and it costs more than the chamber
does. Removing it takes P2_MID from 42.7 to 26.6 ns, a 38 % improvement, with
no change to the detector and no new runs.

### So: what the uRWELL reference buys

1. **A better time reference — the biggest single win.** The stream matching
   (`match_streams.py`) locks the two DAQs to a **10.4 ns residual rms**.
   Referencing VMM hits to the *DREAM* trigger timestamp instead of the VMM's
   own trigger channel replaces a 33.4 ns term with ~10 ns.
2. **A track-defined denominator.** Efficiency stops being "relative to the
   trigger acceptance, a lower bound" and becomes absolute — already done in
   `urw_vmm_efficiency.py`.
3. **No sideband subtraction.** Requiring a *spatial* coincidence with the
   track kills the accidental background that the ~18-fold BCID ambiguity
   creates, instead of subtracting it statistically.
4. **σ as a function of track position, angle and cluster size**, which a
   scintillator cannot give — and the ability to select single-pad on-track
   hits, removing the wrong-pad contamination.

### Caveats on the budget numbers

The widths in §1–3 are an **rms inside ±120 ns of the peak**, not the Gaussian
core σ that `vmm_efficiency.py` fits, so they sit above the campaign table —
use them for the *decomposition*, not as absolutes. And the only sub_run with
hit-level columns on this machine is **run_33/driftscan_gap150V**: a low drift
field on the *bench* gas (Ar/iC₄H₁₀ 95/5), which is a worst case for the drift
term. **Repeat at the nominal working point when EOS is reachable** — the
trigger term should be unchanged, but the intrinsic terms will drop.

## Reproducing

```bash
python3 vmm_trigger_timing.py      # campaign table, from the stored products
python3 vmm_timing_budget.py       # the decomposition, from hit columns
python3 ../../mpgd2026/make_vmm_timing_figs.py -o ../../mpgd2026/figs
```
