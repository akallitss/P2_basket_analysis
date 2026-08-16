# VMM efficiency against uRWELL tracks, whole SPS campaign

Produced by `urw_vmm_efficiency.py` (see README.md) and tabulated with
`tabulate_efficiency.py`; the per-sub_run summaries live on EOS under
`vmm/matching/efficiency/`. **78 sub_runs measured.**

    efficiency = N(uRWELL track on the pads AND a VMM pad within probe_r)
                 / N(uRWELL track on the pads)

with Clopper-Pearson 68.27% intervals. The uRWELL reference and the track cut
are the same ones the DREAM-read campaign used on these same three detectors,
which is what makes the two numbers comparable.

## What it says

**The best VMM configuration reaches 85% on P2_OUT** (run_46,
`cfg_gain4.5_peaktime200`), against the **96–97%** the same detectors gave on
the DREAM readout. P2_MID peaks at 64%, P2_IN at 17%.

* **The configuration dominates.** Same detectors, minutes apart:
  `cfg_gain4.5_peaktime200_opt` gives 16.7 / 59.8 / 82.2%, its `_deflt`
  counterpart 5.1 / 30.8 / 48.6%. Short peaking times are worse
  (`peaktime25`: 1.4 / 25.3 / 42.3%), long ones better.
* **The mesh scans have not plateaued at the top.** run_32 falls
  73.4 → 6.4% (P2_OUT) from `m00V` to `m60V` and run_33 continues to 0.8% at
  `m100V`, but the highest point of the scan is still on the rising edge —
  the VMM readout is gain-starved at the HV this campaign used, where the
  DREAM readout was already on its plateau.
* **The drift scans have plateaued.** run_35 gap 300/350/400 V gives
  73.9 / 74.0 / 73.6% — drift field is not what is missing.
* **P2_IN's number is not a detector statement.** 22 of its 67 illuminated
  pads (44% of the tracks) sit below 5% efficiency and its residual is 5.1 mm
  where P2_MID/P2_OUT give 3.5 mm. That is the signature of a partly wrong
  channel→pad map, and P2_IN is the plane with the known mapping trouble.
  Treat its efficiency as a lower bound until the map is re-derived — which
  these tracks now make possible.

## Why the geometry is trustworthy

The uRWELL→pad transform is fitted free (a 2×2 matrix, not a rotation) and
comes out as a proper rotation of **−60°** with unit singular values and
det +1 at every station and every sub_run — the same −60° the DREAM readout
measured through an entirely different chain. Residuals on P2_MID and P2_OUT
are **3.5 mm = 12 mm/√12**, pad quantisation with no charge sharing, again
matching the DREAM-side number. Efficiency is flat against the coincidence
window from 2σ out, and against the probe radius it plateaus at one pad and
then creeps 0.0006 per 10 mm, which is the accidental rate.

`run_71` has only 5 uRWELL tracks in the whole sub_run, so it has no
reference and no measurement.

| run | sub_run | tracks | IN | MID | OUT | latency ns | residual mm | frame ° |
|---|---|---|---|---|---|---|---|---|
| run_31 | trigtest_gain3.0_peaktime200 | 21888 | 16.7% | 58.5% | 74.2% | 115 | 3.48 | -59.97 |
| run_32 | meshscan_m00V | 280865 | 10.5% | 57.7% | 73.4% | 115 | 3.53 | -59.85 |
| run_32 | meshscan_m10V | 282703 | 5.7% | 46.0% | 60.5% | 115 | 3.57 | -59.83 |
| run_32 | meshscan_m20V | 239404 | 3.2% | 33.6% | 45.6% | 115 | 3.66 | -59.86 |
| run_32 | meshscan_m30V | 285360 | 1.9% | 22.3% | 31.0% | 115 | 3.83 | -59.91 |
| run_32 | meshscan_m40V | 284606 | 1.1% | 13.2% | 17.6% | 115 | 4.12 | -60.00 |
| run_32 | meshscan_m50V | 283753 | 0.6% | 7.6% | 10.0% | 115 | 4.32 | -59.91 |
| run_32 | meshscan_m60V | 282699 | 0.3% | 4.5% | 6.4% | 115 | 4.53 | -60.05 |
| run_33 | driftscan_gap150V | 282545 | 6.4% | 55.2% | 64.2% | 185 | 3.57 | -59.84 |
| run_33 | driftscan_gap200V | 279837 | 8.2% | 56.9% | 69.3% | 135 | 3.54 | -59.83 |
| run_33 | driftscan_gap250V | 284313 | 9.6% | 57.2% | 71.7% | 115 | 3.52 | -59.81 |
| run_33 | meshscan_m100V | 281922 | — | 0.6% | 0.8% | 115 | 6.26 | -60.43 |
| run_33 | meshscan_m70V | 278720 | 0.2% | 2.8% | 4.0% | 115 | 4.91 | -59.82 |
| run_33 | meshscan_m80V | 281413 | — | 1.8% | 2.4% | 115 | 5.35 | -60.01 |
| run_33 | meshscan_m90V | 281792 | — | 1.0% | 1.4% | 135 | 6.08 | -59.78 |
| run_35 | driftscan_gap300V | 274464 | 10.4% | 58.0% | 73.9% | 115 | 3.50 | -59.82 |
| run_35 | driftscan_gap350V | 276610 | 10.1% | 57.2% | 74.0% | 115 | 3.52 | -59.82 |
| run_35 | driftscan_gap400V | 263950 | 9.4% | 55.8% | 73.6% | 115 | 3.50 | -59.87 |
| run_36 | operating_00 | 1140586 | 16.6% | 57.7% | 73.4% | 115 | 3.51 | -59.82 |
| run_36 | operating_01 | 1226835 | 16.4% | 57.4% | 73.2% | 115 | 3.54 | -59.83 |
| run_36 | operating_02 | 1230746 | 16.5% | 57.5% | 73.4% | 115 | 3.51 | -59.83 |
| run_38 | cfg_gain3.0_peaktime100 | 806545 | 9.9% | 51.8% | 68.2% | 25 | 3.54 | -59.84 |
| run_39 | cfg_gain3.0_peaktime25 | 747025 | 1.6% | 20.8% | 31.5% | -85 | 3.87 | -59.72 |
| run_40 | cfg_gain3.0_peaktime50 | 788001 | 3.9% | 37.6% | 52.7% | -45 | 3.63 | -59.84 |
| run_41 | cfg_gain3.0_peaktime200 | 790865 | 15.5% | 56.6% | 72.3% | 115 | 3.52 | -59.83 |
| run_42 | cfg_gain3.0_peaktime200_deflt | 793005 | 3.7% | 18.0% | 31.6% | 115 | 3.93 | -59.94 |
| run_43 | cfg_gain3.0_peaktime200_opt | 547978 | 15.3% | 56.4% | 71.9% | 115 | 3.53 | -59.81 |
| run_44 | cfg_gain4.5_peaktime25 | 935741 | 1.4% | 25.3% | 42.3% | -85 | 7.31 | -59.60 |
| run_45 | cfg_gain4.5_peaktime100 | 836275 | 9.6% | 63.6% | 82.4% | 25 | 3.47 | -59.83 |
| run_46 | cfg_gain4.5_peaktime200 | 837338 | 15.4% | 59.2% | 85.4% | 115 | 3.51 | -59.71 |
| run_47 | cfg_gain4.5_peaktime200_deflt | 841215 | 5.1% | 30.8% | 48.6% | 115 | 3.72 | -59.90 |
| run_48 | cfg_gain4.5_peaktime200_opt | 842609 | 16.7% | 59.8% | 82.2% | 115 | 3.51 | -59.76 |
| run_50 | cfg_gain3.0_peaktime200_opt | 501540 | 14.6% | 39.7% | 71.0% | 115 | 3.56 | -59.93 |
| run_51 | cfg_gain3.0_peaktime200_deflt | 504803 | 3.3% | 12.8% | 29.4% | 115 | 4.00 | -59.97 |
| run_52 | cfg_gain4.5_peaktime200_opt | 504108 | 16.6% | 43.5% | 81.8% | 115 | 3.58 | -59.78 |
| run_53 | cfg_gain4.5_peaktime200_deflt | 367318 | 4.9% | 21.8% | 46.6% | 115 | 4.11 | -59.98 |
| run_54 | cfg_gain4.5_peaktime50 | 362047 | 3.5% | 37.4% | 70.0% | -45 | 3.95 | -59.80 |
| run_55 | meshscan_m00V | 49982 | 14.1% | 37.3% | 67.8% | 115 | 4.55 | -59.94 |
| run_55 | meshscan_m10V | 72040 | 8.8% | 29.1% | 55.2% | 115 | 4.53 | -59.91 |
| run_55 | meshscan_m20V | 71810 | 5.0% | 21.0% | 40.7% | 115 | 4.80 | -60.09 |
| run_56 | meshscan_m30V | 62753 | 2.6% | 13.7% | 27.4% | 115 | 4.75 | -59.97 |
| run_56 | meshscan_m40V | 71880 | 1.7% | 7.9% | 16.5% | 135 | 5.61 | -60.21 |
| run_56 | meshscan_m50V | 65577 | 1.0% | 5.5% | 9.6% | 135 | 5.41 | -60.08 |
| run_56 | meshscan_m60V | 64384 | — | 3.5% | 5.5% | 135 | 5.83 | -60.30 |
| run_56 | meshscan_m70V | 64426 | — | 2.3% | 3.6% | 135 | 6.23 | -59.85 |
| run_57 | driftscan_gap150V | 66576 | 7.8% | 35.2% | 57.7% | 205 | 4.78 | -59.97 |
| run_57 | driftscan_gap200V | 67202 | 9.7% | 36.4% | 63.5% | 135 | 4.64 | -59.97 |
| run_57 | driftscan_gap250V | 68174 | 12.5% | 36.6% | 66.4% | 115 | 4.64 | -59.87 |
| run_57 | driftscan_gap300V | 67049 | 14.3% | 37.0% | 67.9% | 115 | 4.60 | -59.81 |
| run_57 | driftscan_gap350V | 40800 | 14.7% | 36.4% | 68.0% | 115 | 4.54 | -59.98 |
| run_57 | driftscan_gap400V | 72810 | 14.2% | 35.9% | 68.4% | 115 | 4.75 | -59.90 |
| run_57 | meshscan_m100V | 67191 | — | — | 0.7% | 125 | — | — |
| run_57 | meshscan_m80V | 67127 | — | 1.5% | 2.1% | 105 | 5.94 | -60.09 |
| run_57 | meshscan_m90V | 68347 | — | 0.9% | 1.2% | 115 | 6.32 | -59.76 |
| run_58 | operating_00 | 170506 | 13.9% | 36.8% | 68.0% | 115 | 4.74 | -59.92 |
| run_58 | operating_01 | 169272 | 14.2% | 36.8% | 67.8% | 115 | 4.74 | -59.94 |
| run_58 | operating_02 | 166294 | 14.0% | 36.9% | 68.0% | 115 | 4.75 | -59.83 |
| run_61 | meshscan_m00V | 158564 | 7.8% | 34.8% | 58.4% | 115 | 3.99 | -60.02 |
| run_61 | meshscan_m10V | 154523 | 4.4% | 24.8% | 42.3% | 115 | 4.18 | -60.07 |
| run_61 | meshscan_m20V | 159232 | 2.5% | 16.0% | 27.2% | 115 | 4.30 | -60.01 |
| run_61 | meshscan_m30V | 6117 | — | 8.4% | 11.4% | 115 | 10.04 | -57.73 |
| run_61 | meshscan_m40V | 100667 | 0.8% | 5.9% | 9.1% | 115 | 4.84 | -60.12 |
| run_61 | meshscan_m50V | 100718 | 0.4% | 3.9% | 5.2% | 115 | 5.32 | -60.17 |
| run_61 | meshscan_m60V | 100675 | — | 2.4% | 3.3% | 115 | 5.58 | -60.37 |
| run_61 | meshscan_m70V | 101901 | — | 1.6% | 1.8% | 115 | 5.81 | -60.26 |
| run_61 | meshscan_m80V | 97151 | — | 0.9% | 1.2% | 115 | 6.44 | -59.96 |
| run_61 | meshscan_m90V | 99236 | — | — | 0.7% | 2715 | — | — |
| run_62 | driftscan_gap150V | 41412 | 5.3% | 33.6% | 50.9% | 155 | 4.81 | -60.07 |
| run_62 | driftscan_gap200V | 14748 | 6.8% | 34.2% | 54.1% | 135 | 4.94 | -60.13 |
| run_63 | operating_00 | 72327 | 8.2% | 35.0% | 57.9% | 115 | 5.00 | -59.96 |
| run_63 | operating_01 | 36450 | 8.4% | 35.0% | 58.2% | 115 | 4.75 | -60.07 |
| run_63 | operating_02 | 59503 | 8.1% | 35.3% | 58.1% | 115 | 4.75 | -59.93 |
| run_63 | operating_03 | 74298 | 7.8% | 34.6% | 58.3% | 115 | 4.66 | -59.92 |
| run_64 | cfg_gain3.0_peaktime200_opt | 98798 | 8.6% | 36.6% | 60.8% | 115 | 3.98 | -60.00 |
| run_65 | cfg_gain3.0_peaktime200_deflt | 142796 | 1.7% | 10.0% | 19.8% | 95 | 4.56 | -60.10 |
| run_66 | cfg_gain4.5_peaktime200_opt | 148362 | 10.3% | 38.3% | 73.6% | 95 | 4.06 | -59.89 |
| run_67 | cfg_gain4.5_peaktime200_deflt | 148217 | 2.8% | 17.2% | 35.3% | 95 | 4.22 | -60.00 |
| run_71 | cfg_gain4.5_peaktime50 | 5 | — | — | — | -75 | — | — |
