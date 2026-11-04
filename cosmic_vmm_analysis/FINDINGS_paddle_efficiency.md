# Why the tag distributions are not the paddle geometry: the paddles are not uniformly efficient

Written 2026-11-04 (night). Resolves `HANDOFF_det1_tag_residual.md`. Data: wk2 scan runs 181-208, lab stand
(top paddles 41 / 45 on hybrid 0 at z = -100..0, det3 / det2 / det1 at z = 320 / 630 / 940, bottom paddle 44 at
z ≈ 1250). Band coordinates: *along* = the read-out band axis (21.5° to x) from the pad centroid, positive toward the
wide upper-right end of the wedge; *across* = its normal.

## 1. Answer

The det1 tag density (and the det2 one, and part of det3's) is shaped by the **position-dependent efficiency of the
scintillator paddles**, not by their geometry, and not by the chambers:

| paddle | plane | efficiency, lower-left half | efficiency, upper-right half | pulse delay along the band |
|---|---|---|---|---|
| 45 (top) | z = -50 | 0.75 | 0.92 | -4.7 ns/m (earlier toward upper right) |
| 41 (top) | z = -50 | 0.15 | 0.87 | same as 45 (t41 - t45 flat) |
| 44 (bottom, in every trigger) | z = 1250 | **0.90** | **0.23** | **+11 ns/m relative to 45** |

(`paddle_study/paddle_efficiency.json`, chamber-referenced, 12 649 straight tracks of runs 184, 187, 196, 202, 205.)

- **Paddle 44 is in AND in every trigger.** Its efficiency drops from 0.9 to about 0.25 between along = -150 and
  -30 mm at its plane, then stays at 0.2-0.3 up to its upper-right end. Det1, 310 mm above it, sees this step almost
  unsmeared: that is the sharp two-level structure of the det1 tag density at along ≈ 0 (the step moves by about
  -90 → 0 mm because the tracks are tilted, see below), and the "nearly empty upper-right end of det1".
- **Paddle 41 is the mirror image**: 0.95 at its upper-right end, 0.15 at the lower-left. Tracks that fire 41 are
  therefore tilted (they enter at the upper right, where 41 works, and exit through 44 at the lower left, where 44
  works): mean along-band slope -0.22 for "41 only" triggers and -0.15 for "both", against +0.01 for "45 only".
  This is the apparent "slope" of the track population the earlier paddle-position fit tried to absorb.
- **Paddle 45 is fine** (0.9, mildly attenuated toward its lower-left end). It fires in 93 % of the triggers; 41 in
  only 52 %; 7-8 % of the triggers are "41 only" (a thin rim just outside 45's edges).
- The pulse timing says where the PMTs are: 44's pulse comes later toward the upper right, 45's (and 41's) earlier,
  11 ns/m relative to each other, i.e. about 5.5 ns/m each = light propagation in the scintillator. So the PMT of 44
  is at the **lower-left** end and the PMTs of 41 and 45 at the **upper-right** end; each paddle is efficient near
  its own PMT and loses efficiency at the far end: light attenuation against a discriminator threshold, with 45
  having the best light yield / lowest threshold. (The ADC of the four paddle channels is constant, 513-554: the
  VMM sees logic pulses, so pulse height cannot be checked directly.)
- The hardware triple (ch 40) is exactly 41 AND 45 AND 44 (99.9 % identical to "both" in the offline trigger). On
  chamber-only straight tracks: 45 fires in 0.63, 44 in 0.40, 41 in 0.38, ch 40 in 0.17; the offline trigger
  (41 OR 45) AND 44 in 0.31. (These integrate over the whole wedge, including its ends outside the paddles.)

## 2. Closure

`paddle_study.py` feeds the measured 2-D efficiency maps of the three paddles into the geometric MC of
`trigger_acceptance_mc.py` (all chambers 100 % efficient, **no free parameter**: no paddle position fit):

| probe | data mean along | MC, measured paddle maps | MC, uniform 600 x 200 paddles centred on the band |
|---|---|---|---|
| det3 | -4.1 | -3.1 (chi2 1100 / 48 bins) | +1.9 (2133) |
| det2 | -31.0 | -29.7 (106 / 32) | +7.8 (1713) |
| det1 | -68.5 | -65.7 (913 / 48) | +3.2 (10 482) |

`paddle_study/closure.png`: the two-level det1 distribution, the det2 shape and the across-band widths of all
three probes are reproduced. What remains (det3: data below MC at along -100..0 and above at +50..+120; the
narrower across-band peak of det3 / det1; the tails beyond the wedge ends) is the chambers' own efficiency
non-uniformity and the mis-predicted tags, which this MC does not model.

## 3. Why the earlier hypotheses were wrong

- H1 (sheared chamber stack), H3 (det3 map), H4 (z): the straight-track residuals (det2 offsets 0-2 mm in every run)
  were never the issue; the "disagreement growing with depth" is simply the distance to paddle 44.
- The paddle-position fit of `trigger_acceptance_mc.py` could not work: a position moves an edge, it cannot halve the
  density over half the paddle. Its fitted -73 mm bottom-paddle offset is an artefact.
- H2 (paddle efficiency along the length) was the right one; the data in hand (the paddle channels in the capture
  files) were enough to measure it, in two independent ways (triggered sample split by top paddle; chamber-only
  straight tracks as the reference).

## 4. Consequences and what to do on the stand

1. **Efficiency numbers of the chambers are not affected**: the tag is defined by the two other chambers, the
   paddles only select which tracks exist. But the statistics per unit area is 3-4 times lower on the upper-right
   half of every chamber, worst for det1; the kernel-radius maps and the "empty upper-right of det1" are the
   paddle, not the chamber.
2. **Fix the trigger**: lower the discriminator thresholds (or raise the HV) of 44 and 41 until their efficiency
   along the whole length matches 45; or at least swap 44 end-for-end so that its good half is under the half
   where 41 is good (then the triple coincidence covers the same region). Re-measure with
   `paddle_efficiency.py` (no trigger needed: it uses chamber-only tracks) before and after.
3. The hardware triple ch 40 (17 % of the straight tracks) is far too restrictive; the offline (41 OR 45) AND 44 is
   31 %. 45 AND 44 alone would be better than the triple; 45 alone with the chambers as the second constraint is the
   most efficient reference.
4. In `trigger_acceptance_mc.py`, the paddle-position fit should be replaced by the measured maps (as in
   `paddle_study.closure`), and the note on the CERN site (section 5, "det3 positions likeliest") revised.

## 5. Files

- `scint_dump.py`: per-trigger dump of the scintor selection (paddle pulses ADC / time, clean clusters, hits).
- `paddle_efficiency.py`: chamber-only straight tracks (det2 seed, det3 and det1 within ±150 ns, clean clusters,
  det2 within 25 mm of the line), paddle pulse looked up at t(det2) - 242 ns ± 120 ns; sideband 5 µs earlier for
  the chance level (0.000).
- `paddle_study.py`: the figures and jsons above (`analysis_wk2/paddle_study/`): `trigger_composition.json`,
  `paddle_images.png`, `paddle_efficiency.png/.json`, `closure.png/.json`.
- Dumps: `analysis_wk2/dumps/dump_run{181..208}.npz` (all ten runs), `peff_run{184,187,196,202,205}.npz`, with
  their `.json` / `.log`. On lxplus the same under `/tmp/akallits_dump/` (drivers `~/run_dump.sh RUN...`,
  `~/run_peff.sh RUN...`, code in `~/cva3/`).
- Rerun locally: `python paddle_study.py $A/paddle_study --dumps $A/dumps --jsons "$A/eff_definitions/json/scintor_run*.json"`
  with `A=.../analysis_wk2` and `VMM_QA_DIR` set (65 s, 0.75 GB at `--n 3e7`).
