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


## Best-timing working point per detector per gas (2026-08-19)

From `vmm_timing_by_subrun.csv` (median over captures, >=8 captures, contrast
>=8), scanning **every** chip configuration and HV point of the campaign.

**The selection needs an efficiency requirement, and this is not a detail.**
Sorted on sigma alone the "best" points are all mesh 350-400 V with efficiency
0.000 and sigma 4.8-8.8 ns -- *below the 6.5 ns quantisation floor*, so they
cannot be a real coincidence. They are dead-detector runs where the fit latches
onto a narrow accidental structure. Requiring a real signal (efficiency > 0.30):

| station | gas | best sigma | gain | peaking | mesh / drift | sub_run |
|---|---|---|---|---|---|---|
| P2_MID | A | **21.7 ns** | 4.5 | 200 ns | 450 / 750 | `run_48/cfg_gain4.5_peaktime200_opt` |
| P2_MID | B | **18.9 ns** | 4.5 | 200 ns | 450 / 750 | `run_66/cfg_gain4.5_peaktime200_opt` |
| P2_OUT | A | **23.5 ns** | 3.0 | 200 ns | 450 / **850** | `run_57/driftscan_gap400V` |
| P2_OUT | B | **24.4 ns** | 4.5 | 200 ns | 450 / 750 | `run_66/cfg_gain4.5_peaktime200_opt` |

Two things worth saying out loud:

* **`run_66` is the best-timing point for both stations on gas B**, so one
  sub_run gives the whole gas-B column.
* **P2_OUT's best timing is at drift 850, not 750** -- it is still improving
  where P2_MID has plateaued. Its own drift optimum has not been reached.

### Shorter shaping does NOT improve the timing

At the nominal point on gas A, sweeping the peaking time:

| peaking [ns] | P2_MID sigma | P2_MID eff | P2_OUT sigma | P2_OUT eff |
|---|---|---|---|---|
| 25 | 28.9 | 0.216 | 35.3 | 0.334 |
| 50 | 23.7 | 0.362 | 34.7 | 0.568 |
| 100 | 22.1 | 0.560 | 29.0 | 0.656 |
| 200 | 22.8 | 0.375 | 28.3 | 0.603 |

Timing gets *worse* towards short shaping, which is backwards for a shaper and
is the threshold story again: at short peaking the pulse is smaller relative to
a fixed discriminator, so more of the sample sits on the slow part of the
leading edge and time-walks. There is no timing-versus-efficiency trade-off to
optimise here -- 100-200 ns wins on both axes.

## The fit machinery, and how it was validated

`vmm_timing_peaks.py` rebuilds the coincidence from hit columns and fits a
Gaussian on a flat background (the flat term is the ~18-fold BCID ambiguity,
which the campaign subtracts by sideband). It is a reimplementation, so it was
checked against the campaign fit on the one sub_run available locally:

| | this fit | campaign `scalars.json` |
|---|---|---|
| P2_MID, run_33 drift 600 | 37.3 ns | 37.6 ns |
| P2_OUT, run_33 drift 600 | 86.0 ns | 90.2 ns |

Figure `mpgd2026/figs/vmm_timing_peaks_validation.png` is that check, and shows
what the distributions look like: P2_MID a clean Gaussian on a flat pedestal at
peak/bg 19, P2_OUT broad and visibly non-Gaussian at that low drift field.

## Still to run -- needs a Kerberos ticket

**EOS was not reachable on 2026-08-19**: the ticket expired 17 Aug 05:56 and the
renewable window had lapsed (`kinit -R` returns "Ticket expired while renewing"),
so a fresh password login is needed. Everything below is written and validated,
and runs unattended once that is done:

```bash
kinit akallits@CERN.CH
./run_timing_nominal.sh
```

It stages four sub_runs from `/eos/.../vmm/runs/`, produces
`timing_peaks_{midA,midB,outA,nomA}.npz`, re-runs the timing budget at the
**nominal** working point into `TIMING_BUDGET_nominal.txt`, and draws
`mpgd2026/figs/vmm_timing_peaks.png`.

**What to expect from the nominal budget, and what would change the story.**
The existing budget is from run_33 at drift 600 on the bench gas, where the
drift term is large: trigger 33.4 ns against P2_MID's 26.6 ns intrinsic. At
drift 750 the campaign fit gives 22 ns *total* for P2_MID -- already below the
33.4 ns trigger term measured at low field. Those cannot both be right, so one
of two things is true, and the run settles it:

* the trigger term is **smaller** at the nominal point than at run_33 (the
  trigger channel's own walk depends on its pulse height, which is not the same
  run), or
* the rms-inside-+-120 ns estimator used for the budget sits well above the
  fitted core sigma, and the two must be compared like for like.

The budget script should therefore be re-run **and** compared against a core-
sigma estimator on the same sample before the 33.4 ns number goes on a slide.
It is currently in `TIMING.md` and in the deck as the headline of the timing
outlook; if the nominal run does not reproduce it, that claim comes back out.

## Reproducing

```bash
python3 vmm_trigger_timing.py      # campaign table, from the stored products
python3 vmm_timing_budget.py       # the decomposition, from hit columns
python3 ../../mpgd2026/make_vmm_timing_figs.py -o ../../mpgd2026/figs
```
