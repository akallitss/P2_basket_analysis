# VMM efficiency against uRWELL tracks, whole SPS campaign

Produced by `urw_vmm_efficiency.py` (see README.md) and tabulated with
`tabulate_efficiency.py`; the per-sub_run summaries live on EOS under
`vmm/matching/efficiency_v2/`. **78 sub_runs measured.**

    efficiency = N(uRWELL track on the pads AND a VMM pad within probe_r)
                 / N(uRWELL track on the pads)

with Clopper-Pearson 68.27% intervals. The uRWELL reference and the track cut
are the same ones the DREAM-read campaign used on these same three detectors,
which is what makes the two numbers comparable.

> **Re-measured 2026-08-16** with the corrected hot-channel rule (commit
> b3c6b90). The old rule flagged channels at 8x their own chip's median
> occupancy, which masks the beam on any chip whose threshold has been raised.
> P2_MID moves **+4.0 points on average** (up to +12.6, on 42 of the 78
> sub_runs), P2_IN +0.4, P2_OUT +0.07. `EFFICIENCY_AUTOPSY.md` has the story.

## What it says

**The best VMM configuration reaches 85.4% on P2_OUT and 66.9% on P2_MID**
(run_46, `cfg_gain4.5_peaktime200`), against the **96-97%** the same detectors
gave on the DREAM readout of 25 July.

**That gap is the VMM discriminator threshold, not the detector.**
`EFFICIENCY_AUTOPSY.md` takes run_46 apart and excludes everything else by
measurement; the short version is that the threshold sits at ~0.7x the most
probable signal, and the P2_OUT pads that happen to see the largest pulses
already read 0.962 -- the DREAM number -- in this very run.

* **The configuration dominates, because it is the threshold.** Same
  detectors, 36 minutes apart, the same file with the per-chip `sdt` lines
  commented out or set: `cfg_gain4.5_peaktime200_opt` gives 16.9 / 61.5 /
  82.2%, its `_deflt` counterpart 5.3 / 32.5 / 48.6%. Every `_opt`/`_deflt`
  pair in the campaign does the same thing (runs 42/43, 47/48, 50/51, 52/53,
  64/65, 66/67).
* **So does the amplifier gain**, which cannot change a detector: at fixed HV
  and peaking time, 3.0 -> 4.5 mV/fC gives 72.3 -> 82.2% on P2_OUT (run_41 vs
  run_48). Peaking time works the same way: 25 / 50 / 100 / 200 ns gives
  31.5 / 52.7 / 68.2 / 72.3% at gain 3.0.
* **The mesh scans have not plateaued at the top.** run_32 falls
  73.4 -> 6.4% (P2_OUT) from `m00V` to `m60V` and run_33 continues to 0.8% at
  `m100V`, but the highest point of the scan is still on the rising edge --
  the VMM readout is gain-starved at the HV this campaign used, where the
  DREAM readout was already on its plateau. Fitted as a Landau swept past a
  fixed threshold, the scan and the amplifier-gain step are one curve with a
  gain e-folding of 22.4 V, the physical bulk-Micromegas value.
* **The drift scans have plateaued.** run_35 gap 300/350/400 V gives
  73.9 / 74.0 / 73.6% -- drift field is not what is missing.
* **P2_IN's number is not a detector statement.** Three of its six chips run
  at raised thresholds (sdt 256/288/300), VMM 4 carries 1.14e9 hits in one
  sub_run, and its residual is 4.7 mm where P2_MID and P2_OUT give 3.4 mm,
  with a shoulder at 20-30 mm in the matching distance -- the signature of a
  partly wrong channel->pad map, on the plane with the known mapping trouble.
  Treat its efficiency as a lower bound until the map is re-derived, which
  these tracks now make possible.
* **Read P2_IN's low-gain points against its accidental rate.** Below
  `m40V` its accidental column reaches 3-5%, which is the whole of what the
  efficiency column shows there: those points are consistent with zero. On
  P2_MID and P2_OUT the accidental never exceeds 0.14%.

## Why the geometry is trustworthy

The uRWELL->pad transform is fitted free (a 2x2 matrix, not a rotation) and
comes out as a proper rotation of **-60 deg** with unit singular values and
det +1 at every station and every sub_run -- the same -60 deg the DREAM readout
measured through an entirely different chain. Residuals on P2_MID and P2_OUT
are **3.4 mm = 12 mm/sqrt(12)**, pad quantisation with no charge sharing, again
matching the DREAM-side number. And the VMM finds **the same dead pad as
DREAM** -- pad 635 at (411, 237) mm -- through completely different
electronics. Efficiency is flat against the coincidence window from 2 sigma
out, and against the probe radius it plateaus at one pad and then creeps
0.0006 per 10 mm, which is the accidental rate.

`run_71` has only 5 uRWELL tracks in the whole sub_run, so it has no
reference and no measurement.

| run | sub_run | tracks | IN | MID | OUT | latency ns | residual mm | frame ° |
|---|---|---|---|---|---|---|---|---|
| run_31 | trigtest_gain3.0_peaktime200 | 21888 | 16.7% | 58.5% | 74.2% | 115 | 3.48 | -59.97 |
| run_32 | meshscan_m00V | 280865 | 11.4% | 57.7% | 73.4% | 115 | 3.53 | -59.85 |
| run_32 | meshscan_m10V | 282703 | 7.0% | 46.0% | 60.5% | 115 | 3.57 | -59.83 |
| run_32 | meshscan_m20V | 239404 | 4.3% | 33.7% | 45.6% | 115 | 3.67 | -59.86 |
| run_32 | meshscan_m30V | 285360 | 2.7% | 22.3% | 31.0% | 115 | 3.84 | -59.90 |
| run_32 | meshscan_m40V | 284606 | 3.3% | 13.3% | 19.1% | 115 | 4.15 | -59.99 |
| run_32 | meshscan_m50V | 283753 | 3.6% | 7.6% | 11.1% | 115 | 4.34 | -59.88 |
| run_32 | meshscan_m60V | 282699 | 3.8% | 4.5% | 6.4% | 115 | 4.58 | -60.05 |
| run_33 | driftscan_gap150V | 282545 | 7.3% | 55.2% | 64.2% | 185 | 3.57 | -59.84 |
| run_33 | driftscan_gap200V | 279837 | 9.1% | 56.9% | 69.3% | 135 | 3.54 | -59.83 |
| run_33 | driftscan_gap250V | 284313 | 10.5% | 57.2% | 71.7% | 115 | 3.52 | -59.81 |
| run_33 | meshscan_m100V | 281922 | 4.6% | 0.6% | 0.8% | 115 | 6.34 | -60.36 |
| run_33 | meshscan_m70V | 278720 | 5.6% | 2.8% | 4.0% | 115 | 4.94 | -59.82 |
| run_33 | meshscan_m80V | 281413 | 5.3% | 1.8% | 2.4% | 115 | 5.39 | -59.99 |
| run_33 | meshscan_m90V | 281792 | 4.4% | 1.0% | 1.4% | 135 | 6.12 | -59.76 |
| run_35 | driftscan_gap300V | 274464 | 11.3% | 58.0% | 73.9% | 115 | 3.50 | -59.82 |
| run_35 | driftscan_gap350V | 276610 | 11.0% | 57.2% | 74.0% | 115 | 3.52 | -59.82 |
| run_35 | driftscan_gap400V | 263950 | 10.3% | 55.8% | 73.6% | 115 | 3.50 | -59.87 |
| run_36 | operating_00 | 1140586 | 17.5% | 57.7% | 73.4% | 115 | 3.51 | -59.82 |
| run_36 | operating_01 | 1226835 | 17.3% | 57.4% | 73.2% | 115 | 3.54 | -59.83 |
| run_36 | operating_02 | 1230746 | 17.4% | 57.5% | 73.4% | 115 | 3.51 | -59.83 |
| run_38 | cfg_gain3.0_peaktime100 | 806545 | 10.7% | 51.8% | 68.2% | 25 | 3.54 | -59.84 |
| run_39 | cfg_gain3.0_peaktime25 | 747025 | 1.7% | 21.6% | 31.5% | -85 | 3.99 | -59.87 |
| run_40 | cfg_gain3.0_peaktime50 | 788001 | 3.9% | 38.6% | 52.7% | -45 | 3.61 | -59.87 |
| run_41 | cfg_gain3.0_peaktime200 | 790865 | 15.5% | 56.6% | 72.3% | 115 | 3.52 | -59.83 |
| run_42 | cfg_gain3.0_peaktime200_deflt | 793005 | 3.7% | 19.2% | 31.6% | 115 | 3.88 | -59.96 |
| run_43 | cfg_gain3.0_peaktime200_opt | 547978 | 15.3% | 56.4% | 71.9% | 115 | 3.53 | -59.81 |
| run_44 | cfg_gain4.5_peaktime25 | 935741 | 1.1% | 32.9% | 45.2% | -85 | 7.22 | -59.77 |
| run_45 | cfg_gain4.5_peaktime100 | 836275 | 9.6% | 64.0% | 82.4% | 25 | 3.48 | -59.78 |
| run_46 | cfg_gain4.5_peaktime200 | 837338 | 15.5% | 66.9% | 85.4% | 115 | 3.47 | -59.77 |
| run_47 | cfg_gain4.5_peaktime200_deflt | 841215 | 5.3% | 32.5% | 48.6% | 115 | 3.70 | -59.93 |
| run_48 | cfg_gain4.5_peaktime200_opt | 842609 | 16.9% | 61.5% | 82.2% | 115 | 3.50 | -59.78 |
| run_50 | cfg_gain3.0_peaktime200_opt | 501540 | 14.6% | 49.9% | 71.0% | 115 | 3.55 | -59.94 |
| run_51 | cfg_gain3.0_peaktime200_deflt | 504803 | 3.3% | 15.7% | 29.4% | 115 | 3.94 | -60.11 |
| run_52 | cfg_gain4.5_peaktime200_opt | 504108 | 16.8% | 54.6% | 81.8% | 115 | 3.57 | -59.87 |
| run_53 | cfg_gain4.5_peaktime200_deflt | 367318 | 5.1% | 27.7% | 46.6% | 115 | 4.22 | -60.10 |
| run_54 | cfg_gain4.5_peaktime50 | 362047 | 3.5% | 47.0% | 70.1% | -45 | 3.99 | -59.85 |
| run_55 | meshscan_m00V | 49982 | 14.1% | 47.4% | 67.8% | 115 | 4.54 | -59.92 |
| run_55 | meshscan_m10V | 72040 | 8.8% | 37.6% | 55.2% | 115 | 4.63 | -59.97 |
| run_55 | meshscan_m20V | 71810 | 5.0% | 26.9% | 40.7% | 115 | 4.87 | -60.11 |
| run_56 | meshscan_m30V | 62753 | 2.9% | 18.0% | 27.4% | 115 | 4.84 | -60.07 |
| run_56 | meshscan_m40V | 71880 | 1.8% | 10.3% | 16.5% | 115 | 5.45 | -60.22 |
| run_56 | meshscan_m50V | 65577 | 1.0% | 6.2% | 9.6% | 135 | 5.52 | -60.09 |
| run_56 | meshscan_m60V | 64384 | — | 3.5% | 5.5% | 135 | 5.85 | -60.33 |
| run_56 | meshscan_m70V | 64426 | — | 2.3% | 3.6% | 135 | 6.22 | -59.79 |
| run_57 | driftscan_gap150V | 66576 | 7.8% | 45.1% | 57.7% | 205 | 4.78 | -59.99 |
| run_57 | driftscan_gap200V | 67202 | 9.7% | 47.1% | 63.5% | 135 | 4.62 | -59.96 |
| run_57 | driftscan_gap250V | 68174 | 12.5% | 47.5% | 66.4% | 115 | 4.64 | -59.90 |
| run_57 | driftscan_gap300V | 67049 | 14.3% | 48.0% | 67.9% | 115 | 4.64 | -59.84 |
| run_57 | driftscan_gap350V | 40800 | 14.7% | 47.4% | 68.0% | 115 | 4.55 | -59.95 |
| run_57 | driftscan_gap400V | 72810 | 14.2% | 46.6% | 68.4% | 115 | 4.75 | -59.90 |
| run_57 | meshscan_m100V | 67191 | — | — | 0.7% | 125 | — | — |
| run_57 | meshscan_m80V | 67127 | — | 1.5% | 2.1% | 105 | 5.94 | -60.09 |
| run_57 | meshscan_m90V | 68347 | — | 0.9% | 1.2% | 115 | 6.32 | -59.76 |
| run_58 | operating_00 | 170506 | 13.9% | 47.4% | 68.0% | 115 | 4.75 | -59.92 |
| run_58 | operating_01 | 169272 | 14.2% | 47.6% | 67.8% | 115 | 4.72 | -59.91 |
| run_58 | operating_02 | 166294 | 14.0% | 47.7% | 68.0% | 115 | 4.78 | -59.86 |
| run_61 | meshscan_m00V | 158564 | 7.8% | 43.9% | 58.4% | 115 | 4.01 | -60.00 |
| run_61 | meshscan_m10V | 154523 | 4.4% | 31.5% | 42.3% | 115 | 4.20 | -60.09 |
| run_61 | meshscan_m20V | 159232 | 2.5% | 20.5% | 27.2% | 115 | 4.37 | -60.06 |
| run_61 | meshscan_m30V | 6117 | — | 10.6% | 11.4% | 115 | 9.58 | -60.85 |
| run_61 | meshscan_m40V | 100667 | 0.8% | 7.1% | 9.1% | 115 | 4.77 | -60.11 |
| run_61 | meshscan_m50V | 100718 | 0.4% | 4.1% | 5.2% | 115 | 5.24 | -60.32 |
| run_61 | meshscan_m60V | 100675 | — | 2.4% | 3.3% | 115 | 5.57 | -60.38 |
| run_61 | meshscan_m70V | 101901 | — | 1.6% | 1.8% | 115 | 5.90 | -60.23 |
| run_61 | meshscan_m80V | 97151 | — | 0.9% | 1.2% | 115 | 6.44 | -59.96 |
| run_61 | meshscan_m90V | 99236 | — | — | 0.7% | 2715 | — | — |
| run_62 | driftscan_gap150V | 41412 | 5.3% | 42.3% | 50.9% | 155 | 4.94 | -60.07 |
| run_62 | driftscan_gap200V | 14748 | 6.8% | 43.4% | 54.1% | 135 | 5.01 | -60.09 |
| run_63 | operating_00 | 72327 | 8.2% | 43.8% | 57.9% | 105 | 5.04 | -59.98 |
| run_63 | operating_01 | 36450 | 8.4% | 44.0% | 58.2% | 105 | 4.86 | -60.09 |
| run_63 | operating_02 | 59503 | 8.1% | 44.1% | 58.1% | 115 | 4.82 | -59.96 |
| run_63 | operating_03 | 74298 | 7.8% | 43.5% | 58.3% | 115 | 4.81 | -59.93 |
| run_64 | cfg_gain3.0_peaktime200_opt | 98798 | 8.6% | 46.1% | 60.8% | 105 | 3.96 | -60.01 |
| run_65 | cfg_gain3.0_peaktime200_deflt | 142796 | 1.7% | 12.3% | 19.8% | 95 | 4.53 | -60.22 |
| run_66 | cfg_gain4.5_peaktime200_opt | 148362 | 10.5% | 50.9% | 73.6% | 95 | 4.06 | -59.89 |
| run_67 | cfg_gain4.5_peaktime200_deflt | 148217 | 3.0% | 23.4% | 35.3% | 95 | 4.32 | -60.12 |
| run_71 | cfg_gain4.5_peaktime50 | 5 | — | — | — | -75 | — | — |
