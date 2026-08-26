# p2_sps — VMM3a and DREAM on the same uRWELL pads

P2 basket, SPS H4, July 2026. The two readouts of the P2 stations were run over
the same detectors; this package puts them side by side **pad by pad** and
answers why the VMM is less efficient than DREAM on P2_OUT.

## The result

The pad-to-pad gain variation belongs to the **chamber**: a factor 3.9 roll-off
across the beam spot, one smooth gradient, measured identically by both
readouts (per-pad r = +0.94, plane fit same direction to 5°). What the VMM adds
is a **discriminator level that sits inside the Landau**. One fitted level,
162 DREAM ADC, cut into DREAM's own per-pad spectra reproduces the VMM's
efficiency pad by pad — 53 pads, one free number, r = +0.90, 85.5 % predicted
against 85.3 % measured. On an average pad that level is 0.69× the most probable
pulse; on a weak pad it is 1.19×, so there it eats the peak, not the tail.

The same cut also explains why the VMM's copy of the gain map looks *flatter*
than DREAM's: relative pad-to-pad rms goes 0.423 (DREAM, everything kept) →
0.355 (cut at the threshold) → 0.260 (VMM), and the last step is an additive
+43-count offset in the VMM ADC. Nothing is left over.

Notes (unlisted, world-readable by URL):

- [Why the VMM loses the weak pads](https://dylan-neff.web.cern.ch/notes/p2out-vmm-threshold-weak-pads.html) — the full case, with the caveats
- [VMM vs DREAM pad gain](https://dylan-neff.web.cern.ch/notes/p2out-vmm-vs-dream-pad-gain.html) — the comparison it rests on
- [Three slides](https://dylan-neff.web.cern.ch/notes/p2out-vmm-three-slides.html) — the MPGD2026 sequence, with speaker notes
- Figures: <https://dylan-neff.web.cern.ch/p2_sps/>

## The second question: how well do the three of them track?

`p2_selftrack.py` asks what `urw_p2_efficiency.py` cannot — not "did the
station fire?" but "where does the P2-only track say the particle went?", with
the uRWELL standing in for truth. All three stations are put in one frame at
once (each station's uRWELL &rarr; pad affine inverted), a straight line is
fitted through the three clusters, and the result is compared with the
reference track.

The answer, on `eff_nominal_1`: detection is fine, geometry is pad-limited.
Every station's residual is a 12 mm box, σ_core ≈ 3.65 mm — pitch/√12, with no
charge sharing — and a two-pad cluster is not better than a single pad. The
P2-only track points to ~4.3 mm and ~0.75 mrad in the core, on a heavy tail.

**And P2 cannot measure any of that on its own.** Checking P2_MID against the
line through P2_IN and P2_OUT returns 0.18 mm, 23× better than the 4.3 mm the
reference measures on the same events, because **70 % of three-station tracks
have all three stations reporting the identical pad**. The self-consistency
residual is the same rounding error three times over. That is the tracking
version of the 4-point efficiency gap in the handoff's §13.4, and it is much
bigger.

Report: `report_track.html` —
[published](https://dylan-neff.web.cern.ch/notes/p2-tracking-vs-reference.html).
Analysis chain in the table below.

The committed `data/p2_selftrack_*.npz` carries every histogram and map whole
but a **thinned** per-track sample (`track_stats.py --trim`); the full sample is
on EOS at `analysis/selftrack/<run>/`. Widths recomputed from the thinned
sample agree with the full one to 2 %.

## The third question: the dead areas, and the bulk pillars

Two things at once, off one new extraction stage (`p2_pillars.py`, which reuses
`p2_selftrack.build_joined` so the tracks, alignment, fiducial and vetoes are
literally the same): **cut the dead areas out of each chamber and re-measure the
stack**, and **can the reference resolve the bulk support pillars?**

The pillars, yes — comfortably. The reference points to **0.21–0.26 mm** at the
P2 planes (measured two independent ways: the edge of the pad-residual box, and
the ratio of the pillar lattice's first and second harmonics, in which the
pillar's size and deadness cancel). At that resolution a 3 mm period keeps 80 %
of its contrast, and the lattice shows up in the efficiency map as six spots
60° apart at **15–19 σ**. Folding every pillar in the beam spot onto one cell
images the average one: efficiency drops to 0.39 (IN) / 0.58 (OUT) at the
centre and the leading pad loses 40 % of its charge.

Three results worth carrying:

- **The as-built lattice is triangular, d = 2.995 ± 0.004 mm, rows along the
  board x axis, effective ⌀ 0.75 mm in charge** — the same on all three
  chambers (each station's fitted lattice read on its neighbours' maps recovers
  ~90 % of the amplitude). That matches the *symmetry, orientation and diameter*
  of the archived April artwork (`P2_BASKET-Mask_M2_V*.gbr`, triangular
  ⌀0.80 mm) but **not its 4.000 mm spacing**, and nothing about the CERN May
  bulk order (`P2_Mask2.gbr`, square 2.000 mm, ⌀0.500 mm). Neither file in the
  repositories is the mask this bulk was made from.
- **The pillars are most of the inefficiency**: 9.5 pp of P2_IN's 10.7,
  3.5 pp of P2_OUT's 4.5. The *charge* footprint is the same on all three, the
  *efficiency* footprint is not — a pillar only costs a hit when the charge it
  removes takes the pad under threshold, so the weakest station pays the most.
- **There is almost nothing to mask** in the `eff_nominal_1` beam spot: ~75 live
  pads per station, no dead channels, no weak connector block, one ⌀6.15 mm
  support pillar. Masking it moves every tracking number by < 0.1 %. The
  documented P2_IN connector-11 problem is under the `highstat_eff_1` spot,
  which is processed the same way.

The ⌀6.15 mm support pillar is worth looking at on its own: binned at 0.2 mm it
is a round hole sitting exactly on its gerber position, efficiency 1.6–2.3 %
inside against 91–98 % a few millimetres away. That is the calibration for
everything above — the reference images a millimetre-scale feature at its true
place, and a pillar is dead, not merely attenuating.

Report: `report_pillar.html`. **Only the derived products are committed**
(`data/pillar_numbers_*.json`, `data/pillar_arrays_*.npz`); the raw
`p2_pillars_*.npz` are 45–78 MB and live on EOS under `analysis/pillars/`.

## The runs

The comparison only means anything at a **matched operating point**. Every VMM
`cfg_*` run sits at mesh 450 V / drift 750 V; the DREAM run at the identical
point is **`eff_nominal_1`** (27 July, 17 sub_runs, 1.87 M fiducial tracks) —
*not* `highstat_eff_1`, which is at drift 700. VMM side: `run_46 /
cfg_gain4.5_peaktime200`, 0.84 M good tracks. All six P2_OUT chips are at
`sdt = 224`, which is what makes "one threshold" a meaningful statement.

Pad 635 is dead in both readouts and is excluded from the threshold fit.

## Chain

Everything downstream of `data/` runs locally from what is committed here; only
the two extraction steps need lxplus and the raw data.

| Step | Script | Writes |
|---|---|---|
| DREAM per-pad efficiency + spectra (lxplus, LCG_110) | `urw_p2_padadc.py` | `data/dream_padadc_*` |
| P2 self-tracking vs the reference (lxplus, LCG_110) | `p2_selftrack.py` | `data/p2_selftrack_*` |
| Fine maps, pad ledger, box-edge histograms (lxplus) | `p2_pillars.py` | EOS `analysis/pillars/p2_pillars_*` |
| Bulk-pillar geometry from the mask gerber | `pillar_geom.py --build` | `data/p2_pillars_gerber.npz` |
| Lattice fit, folds, masks, resolution | `pillar_stats.py`, `pillar_numbers.py` | `data/pillar_numbers_*`, `data/pillar_arrays_*` |
| VMM per-pad ADC, tracked and untracked | `pad_adc.py`, `compare_tracked.py` | `data/pad_adc_*` |
| Join the two on the pad map | `compare_dream_vmm.py` | `data/compare_dream_vmm_P2_OUT.csv` |
| Fit the one threshold, cost the fixes | `threshold_model.py` | `data/threshold_model_P2_OUT.json` |
| Aggregate the self-tracking counts | `track_stats.py` | (in memory) |
| Figures | `figures.py`, `figures_p2out.py`, `figures_dv.py`, `figures_slide.py`, `figures_deck.py`, `figures_track.py`, `figures_pillar.py` | `figures/*.png` |
| Reports | `make_report.py`, `make_report_dv.py`, `make_report_slide.py`, `make_deck_mockup.py`, `make_report_track.py`, `make_report_pillar.py` | `report*.html`, `deck_mockup.html` |
| Notes for the site | `make_note.py` | `*_note.html` (PNGs inlined as data: URIs) |

`figures/` and `*_note.html` are **not committed** — they are regenerated from
`data/` by the scripts above, which is the convention in this repository. Run
them with the project venv (`../../.venv/bin/python` from here).

## The MPGD2026 slides

`figures_deck.py` builds the three-slide sequence, 16:9 at 160 dpi:

1. `deck_1_deficit.png` — same chamber, same beam, same working point: DREAM
   95.6 %, VMM 85.3 %, and the loss is a *place*, not a scale factor.
2. `deck_2_gainmap.png` — that place is the low-gain corner, and both readouts
   measure the same gain map.
3. `deck_3_threshold.png` — the ridgeline: on a log axis a gain factor is a
   sideways shift, and the Landau slides onto a fixed discriminator line.

`figures_slide.py` builds the backup slides — `slide_proof.png` (predicted vs
measured, and the spread waterfall) and `slide_fix.png` (efficiency against
signal-over-threshold).

## Gotchas

- `urw_p2_padadc.py` hardcodes the banco data path — pass `--data-root`
  explicitly, and make sure `REPO_ROOT/Detector_Mapping/P2_BASKET/` exists.
- `extract_vmm_triggers.iter_columns` yields no `adc` on the hits_store branch
  and `urw_vmm_efficiency.station_hits` silently substitutes `np.ones(...)`. Any
  sub_run with a full column store gets **fake ADC = 1**. `run_46` escapes only
  because its capture 00001 lacks npy and falls back to pcapng.
- Pad **area** is not a driver of the gain spread: across the illuminated pads
  area spans 0.16 % while pulse height spans a factor 2.9, so a regression on
  area is a regression on position in disguise. An earlier version of the
  published note said otherwise and was corrected.
- P2_MID and P2_IN are **not** in the VMM/DREAM comparison. Their chips sit at
  higher `sdt`, and P2_MID has 27 of 64 channels dead on VMM 12. They *are* in
  the tracking study, which is DREAM-only.
- **Three sub-runs of `eff_nominal_1` are missing P2 FEUs from
  `combined_hits_root`** — `eff_nominal_10` and `_15` have no FEU 4 (P2_MID),
  and `_14` has neither FEU 4 nor FEU 5 (P2_OUT). Not a detector
  inefficiency: the FEU is absent from the merged file. `p2_selftrack.py`
  refuses those sub-runs and says so; 14 of 17 are used.
- **The pillar lattice is fitted, never assumed.** `hex_scan` scans spacing and
  orientation together rather than three independent wavevectors: fitted
  separately they disagree by parts per thousand in |G|, which over a 64 mm
  window is a radian of phase and puts the folded origin on an interstitial
  site. For the same reason the fold's origin is taken from the *minimum of the
  folded cell*, not from the harmonic phases, and `fold_hex` must carry the
  fitted θ — a fifth of a degree is 0.2 mm of drift across the window, which is
  the width of the thing being folded. Both of those were bugs first.
