# BCIDs and read-out-stream positions of the beam tracks (SPS H4, VMM)

Written 2026-11-05. The same study as `FINDINGS_bcid_stream.md` (cosmic stand), on the SPS beam data. Same FEC,
firmware and cabling: P2_IN / P2_MID / P2_OUT = VMM 2-7 / 8-13 / 14-19 (det5 / det1 / det3), the beam trigger
(scintillators) on VMM 0 ch 44.

- **run_46** `cfg_gain4.5_peaktime200` (1 Aug, mesh 450 / drift 750, Ar/CO2/iC4H10 93/5/2, the reference run of the
  VMM efficiency study): 48 captures, 35 min, 9.74 M triggers, **856 953 straight tracks**.
- **run_66** `cfg_gain4.5_peaktime200_opt` (3 Aug, Ar/CF4/iC4H10 88/10/2): 41 captures, 30 min, 4.23 M triggers,
  **129 011 straight tracks** (fewer: P2_IN efficiency is lower on this run).

Track: a trigger with no other trigger within ±1 µs, and in each station a clean cluster (≤ 6 hits within 15 mm,
hot channels excluded) within ±150 ns of trigger + the station's latency; P2_MID within 25 mm of the P2_IN-P2_OUT
midpoint. The latencies, fitted on the first capture, are P2_IN 225 / P2_MID 115 / P2_OUT 125 ns (run_46) and
205 / 105 / 115 ns (run_66), as in `vmm_dream_matching/efficiency_table.csv`. Code: `bcid_stream_sps.py` (lxplus,
`~/run_bcid_sps.sh RUN SUB`) + `bcid_stream_plots.py --setup sps`. Outputs:
`data/SPS_Beam_Test/bcid_stream/` (`bcid_sps_run_{46,66}.npz/.json`; `run_46/`, `run_66/`: `bcid_tracks.png`,
`bcid_marker.png`, `bcid_summary.json`; each run json has 12 stream excerpts under `examples`). Numbers below are
run_46; run_66 agrees except where noted.

## 0. The stream at the SPS

- Same structure as in the lab: markers every 65 536 ticks (1.47 ms), in rounds of all 20 VMMs with the same
  timestamp. The **`offset` field counts 0-15 uniformly**, as in the lab: the statement in
  `vmm_dream_matching/TIMING.md` that it is stuck at a constant in this firmware does not hold for these data (and
  the marker period is 1.47 ms, not 1.6384 ms, with the 22.5 ns tick). Full hit times are therefore unambiguous
  within a marker period. Only 0.03 % of the track hits have offset -1.
- What is different is the load: 33 M hit words per 45 s capture, so markers are only 1.8 % of the words (lab:
  92 %), and about 600 hit words lie between a VMM's marker and a track hit (lab: 7-18). Most of it is one channel:
  VMM 4 ch 39 (P2_IN) carries 25 M of the 33 M hit words.

## 1. Time to the SRS marker

Uniform over 0-1.47 ms (median 737 µs) for all three stations, offset uniform 0-15, the same as for all hits. Word
position in the UDP frame uniform. As in the lab, this carries no information.

## 2. BCID differences between the three stations

All three stations' hits use the **same marker** (99.995 %); 99.99 % of the tracks lie in one UDP frame.

| pair | ΔBCID 68 % | ΔBCID 90 % | \|ΔBCID\| ≤ 1 | \|ΔBCID\| ≤ 3 | Δt (BCID + TDC) 68 % | sd |
|---|---|---|---|---|---|---|
| P2_MID − P2_IN | −7 … −2 | −8 … 0 | 12 % | 37 % | −152 … −41 ns | 56 ns |
| P2_OUT − P2_IN | −7 … −1 | −9 … +1 | 16 % | 42 % | −150 … −23 ns | 66 ns |
| P2_OUT − P2_MID | −1 … +2 | −2 … +4 | 64 % | 93 % | −27 … +45 ns | 42 ns |

- P2_MID and P2_OUT behave as the lab chambers: ±1.5 BC around zero (93 % within ±3 BC).
- **P2_IN comes about 5 BC (110 ns) later** than the other two: its latency is 225 ns against 115 / 125 ns, and
  its spread is wider (68 % width 5 BC instead of 3). That is a property of the P2_IN chamber (its timing is also
  the worst in `TIMING.md`), not of the read-out: relative to its own latency, P2_IN is within ±2.5 BC.
- **Trigger → stations:** the stations' BCIDs are later than the trigger's by P2_MID 5 [4, 7], P2_OUT 6 [4, 9],
  P2_IN 10 [6, 13] BC (median [5-95 %]); 99.99 % share the trigger's marker.

## 3. Where the track hits sit in the hit list

- **One block, as in the lab.** Span from the first to the last track hit: 2 hit words (median; 2-3 for 95 %, 4 at
  99 %). 92.6 % of the tracks have nothing else in between.
- **The trigger is read 2-4 words before P2_MID** (median −3; 5-95 %: −4 … −2), i.e. just before the block: the
  order is trigger, then the three stations.
- Non-track hits inside the span: 6.6 % of the tracks have another over-threshold chamber hit inside, but 6.3 %
  is a hot channel (mostly VMM 4 ch 39); without the hot channels it is about 0.3 %. Neighbour hits 0.003 %, another
  hybrid-0 (trigger) hit 0.9 %.
- **Station order in the block:** P2_MID-P2_OUT-P2_IN 47 %, P2_IN-P2_MID-P2_OUT 35 %, P2_MID-P2_IN-P2_OUT 10 %,
  P2_OUT-P2_IN-P2_MID 5 %, the rest 3 % (run_66: 43 / 39 / 10 / 5 %). In the lab the VMM-ID order
  (det3-det2-det1, the same VMM slots as IN-MID-OUT) dominated at 66 %. Here P2_IN, whose hit comes 5 BC later,
  ends up last almost half the time, so the order inside the block follows the hit time at least partly. Not
  understood in detail.
- **Chance reference** (the first hit of each station 5 µs later): span median 48 words, 9-142 for 68 % (run_66:
  65): at beam rates the list is dense, and the 3-word track block stands out much more clearly than in the lab.

## 4. Time order of the stream

At the SPS the stream is almost time-ordered: 0.19 % of the chamber hits are read after a hit more than 1 µs later,
0.08 % > 20 µs (lab: 17 % and 1.3 %). Within one VMM, 2-8 % of consecutive hits go backwards in time (lab:
14-27 %), with p99 of the backward step 0.2-2 µs. A flat 0.07 % tail up to > 100 µs, the same for the trigger
hits, is probably mis-assigned markers (not investigated). For event building at the SPS a look-back of a few µs is
enough; for the trigger-referenced tracks a window of a few hit words in stream order already contains them.

## 5. Same as the lab / different

| | lab (cosmics, 21 885 tracks) | SPS run_46 (856 953 tracks) |
|---|---|---|
| same marker for all hits of a track | 99.995 % | 99.995 % |
| in one UDP frame | 99.97 % | 99.99 % |
| span in hit list (median / p99) | 2 / 4 | 2 / 4 |
| nothing else inside the span | 91.8 % | 92.6 % |
| ΔBCID 68 % (chambers at the same latency) | ±2 | −1 … +2 (MID-OUT) |
| reference (trigger / paddle) read before the middle chamber | 2-3 words | 3 words |
| reference → chamber BCID | 10-11 BC (paddles) | 5-6 BC (MID / OUT), 10 BC (IN) |
| order = VMM-ID order | 66 % | 35 % |
| hits read > 1 µs late | 17 % | 0.19 % |
| markers / all words | 92 % | 1.8 % |
