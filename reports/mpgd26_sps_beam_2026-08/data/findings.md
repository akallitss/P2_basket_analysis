# Where do we lose efficiency — P2_IN / P2_MID / P2_OUT (SPS July 2026)

Working point: mesh 450 V / drift 700 V, run `highstat_eff_1`, sub_run `beam_commissioning_00`
(~1.15 M uRWELL tracks per station per sub_run; 6 sub_runs, ~6.1-6.3 M tracks total).
All efficiencies below are uRWELL-referenced absolute values unless labelled "tag-probe".
Cold-pad / uniformity numbers use the per-pad tag-probe map, pads with n_tag >= 200 only.

## 1. Loss budget (absolute percentage points, working-point sub_run)

| mechanism                                | P2_IN | P2_MID | P2_OUT |
|------------------------------------------|------:|-------:|-------:|
| total inefficiency (1 - eff)             |  3.51 |   2.94 |   3.96 |
| — miss with NO P2 cluster in event       |  2.92 |   2.42 |   3.33 |
| — miss with far cluster (>15 mm)         |  0.59 |   0.52 |   0.63 |
| cold-pad loss (tag-probe map)            |  2.92 |   0.07 |   0.07 |
| uniform floor (median healthy-pad ineff.)|  3.95 |   3.60 |   7.76 |
| spark-veto deadtime (live-time, not eff) |  0.00 |   0.00 |   3.13 |
| efficiency drift over run (max-min)      |  0.40 |   1.31 |   1.40 |

Working-point efficiencies: P2_IN 96.49 %, P2_MID 97.06 %, P2_OUT 96.04 %
(track-weighted aggregate over all 6 sub_runs: 96.34 / 96.31 / 95.16 %).

## 2. Miss composition: "empty event" dominates

83.3 / 82.2 / 84.1 % of all misses (IN/MID/OUT) have no P2 cluster anywhere in the
event; only 16-18 % have a cluster somewhere but >15 mm from the track. So the
inefficiency is overwhelmingly a "detector did not respond" effect (gain/threshold),
not mis-reconstruction or misalignment. The eff-vs-probe-radius curves are flat from
r = 8 mm on (e.g. P2_IN 96.05 % at 8 mm -> 96.66 % at 40 mm), and accidentals are
negligible (~7e-4 per 10 mm window).

## 3. Cold regions

### P2_IN — connector 11 partially cold + connector 14 nearly dead (real hardware issue)
- Connector 11 (channel_id 704-767): beam-illuminated block mean eff 0.72 vs 0.955
  (MID) and 0.922 (OUT) on the same channels. Channels 759, 764, 761, 762, 766, 765
  sit in the beam core (n_tag 13k-104k) at eff 0.28-0.48 and alone account for
  ~2.8 pp of the 2.92 pp cold-pad loss.
- Connector 14 (896-959): 6 beam-fringe pads at eff 0.019 — effectively dead
  (channels 896-906; small n_tag so only ~0.1 pp weight in this beam spot).
- This is a channel_id-block effect (contiguous channels within //64 groups), i.e.
  connector/FEU-input, not a gas/mesh region. It is unique to P2_IN.
- Note: the uRWELL-referenced P2_IN eff (96.5 %) is only 0.6 pp below P2_MID even
  though the tag-probe cold loss is 2.9 pp — with the 15 mm probe window part of the
  weak-connector tracks are still recovered by neighbouring pads; in a strict per-pad
  sense the conn-11 region is the dominant localized loss of the whole system.

### P2_MID / P2_OUT — pillar footprints, negligible loss
Cold pads in MID and OUT are the SAME channel_ids in both stations and form two
tight spatial clusters:
- cluster A at ~(325, 193) mm: channels 557/562/575/579/598/618/643/657/659/678/681/697/701/717/736/740/755
- cluster B at ~(347, 300) mm: channels 852/861/863/877/878/879/883/886/887/888/889/890/891/892

Identical positions on all three detectors (562/579/618/755 are cold on P2_IN too)
= the known support-pillar dead spots (2 of the 4 big pillars fall inside this beam
spot). The pillars also suppress tagging there (n_tag 200-500 vs 10k-120k on
neighbours, because the tag stations' own pillars kill the tag), so their weight is
tiny: 0.07 pp for both MID and OUT, ~0.1 % of tags.

## 4. Uniformity (per-pad tag-probe eff, n_tag >= 200)

| station | pads | p5    | p25   | p50   | p75   | p95   |
|---------|-----:|------:|------:|------:|------:|------:|
| P2_IN   | 195  | 0.286 | 0.912 | 0.957 | 0.965 | 0.973 |
| P2_MID  | 242  | 0.335 | 0.940 | 0.961 | 0.968 | 0.977 |
| P2_OUT  | 247  | 0.331 | 0.888 | 0.918 | 0.933 | 0.946 |

(the low p5 reflects the pillar/connector pads above). P2_OUT's whole distribution
is shifted ~4 pp low — a uniform gain/threshold deficit, not localized damage.

## 5. Tag-probe vs uRWELL-referenced scalars (why 92.5/96.3/92.5 vs 96-97)

| station | tag-probe scalar | map tag-weighted mean | uRW-referenced |
|---------|-----------------:|----------------------:|---------------:|
| P2_IN   | 0.9251           | 0.9263                | 0.9649         |
| P2_MID  | 0.9628           | 0.9643                | 0.9706         |
| P2_OUT  | 0.9253           | 0.9266                | 0.9604         |

Cross-check: the n_tag-weighted mean of the per-pad map reproduces the station
scalar to 0.1-0.2 pp — internally consistent.

The tag-probe numbers are biased LOW because the tag position is built from the
other two P2 stations' 12 mm pad centroids: the extrapolated/interpolated prediction
is smeared by several mm, so tracks near pad borders and near the beam edge fail
the probe match. The bias is smallest for P2_MID (-0.7 pp), which is *interpolated*
between IN and OUT, and largest for the outer stations (-4.0 / -3.5 pp), which
require *extrapolation*. The uRWELL telescope (sub-mm, rigid-frame RMSE ~4.9 mm
residual, flat eff-vs-r above 8 mm) is the absolute reference; use 96-97 %.
Consequently the tag-probe "uniform floor" (3.6-7.8 pp) overstates the true uniform
inefficiency; the uRW numbers bound it at ~2.4-3.3 pp (the no-hit fraction).

## 6. Spark veto / live time (highstat_eff_1)

- P2_IN: 3 sparks total in 6 sub_runs (live_fraction 0.993-1.000); mean deadtime 0.46 %.
- P2_MID: 1 spark (sub_run 03, live 0.9934); mean deadtime 0.16 %.
- P2_OUT: 32 sparks; per-sub_run deadtime 3.1 / 10.9 / 0.9 / 5.5 / 0.8 / 3.4 %
  (mean 4.1 %, worst sub_run 01 at 10.9 %). Mean imon 0.07-0.13 uA vs 0.03-0.05 uA
  for IN/MID.
- This is removed from the efficiency denominator, so it costs DAQ live-time, not
  detector efficiency — but P2_OUT is the sparky station and loses ~4 % of beam time.

## 7. Efficiency drift over the run

Eff decreases monotonically over sub_runs 00 -> 04 and partially recovers in 05:
P2_IN 96.49 -> 96.21 -> 96.60 (span 0.40 pp), P2_MID 97.06 -> 95.75 -> 96.72
(1.31 pp), P2_OUT 96.04 -> 94.63 -> 95.67 (1.40 pp). A ~1.3-1.4 pp rate/charging
drift on MID and OUT; P2_IN is 3x more stable.

## Bottom line
- True working-point efficiency 96-97 % everywhere; the deficit is ~83 % "no
  response at all", uniformly spread — i.e. a gain/threshold floor, not geometry.
- The one real localized hardware problem is P2_IN connector 11 (weak, eff ~0.72
  in-beam) with connector 14 nearly dead on the fringe — worth a cabling/FEU check.
- Pillar dead spots are visible but cost only ~0.07 pp in this beam spot.
- P2_OUT: overall ~1 pp lower, 4 pp lower per-pad distribution, plus it sparks
  (~4 % live-time loss) and drifts (1.4 pp) — the station to watch at HV.
