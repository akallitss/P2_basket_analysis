# BCIDs and read-out-stream positions of the reconstructed muon tracks

Written 2026-11-05. Data: wk2 scan runs 181-208 (ten full runs, 30 h), lab stand, det3 / det2 / det1 = VMM 2-7 /
8-13 / 14-19, paddles on hybrid 0 (VMM 0). Tracks: the chamber-only straight tracks of `paddle_efficiency.py`
(det2 seed, det3 and det1 within ±150 ns, clean clusters, det2 within 25 mm of the det3-det1 midpoint), spark-vetoed:
**21 885 tracks**. Code: `bcid_stream.py` (lxplus, `~/run_bcid.sh RUN`) + `bcid_stream_plots.py`. Outputs:
`data/Lab_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_wk2/bcid_stream/` (`bcid_run*.npz/.json`,
`bcid_tracks.png`, `bcid_marker.png`, `bcid_summary.json`; each run json has 12 stream excerpts under `examples`).

## 0. The stream and a clock correction

- The FEC sends 1491-word UDP frames. About 92 % of the words are SRS markers: every VMM gets a marker every
  **65 536 ticks**, all 20 VMMs in one round with the same timestamp (99 % of the rounds complete). Between two
  rounds come the hits of that period, with `offset` = 0..15 (4096-BC overflows since the marker). Only 0.05 % of
  the hits have offset -1 (read after the next marker).
- **The marker tick is 22.5 ns, not 25 ns.** Against the pcap wall clock the marker timestamp advances by 2.670e9
  ticks per 60.08 s (22.4998 ns per tick; runs 121, 130, 196 agree). This also follows from the format: the hits
  of one marker period span exactly 16 × 4096 BC = 65 536 BC. `coincidence_check.read_times` and the QA compute
  `srs * 25 ns`. That does not change time differences inside one marker period (coincidences, latencies,
  efficiencies), but it stretches absolute times by 11 % and leaves a 164 µs gap at each marker. Rates per second
  of span come out about 10 % low, and the per-file `WallClock` mapping used by the spark veto drifts by up to
  about 6 s at the end of each 60 s file. `bcid_stream.py` uses 22.5 ns (`--srs-tick`). Fixing `read_times` is a
  separate change: it affects every earlier rate number.

## 1. Time to the SRS marker

| chamber (first cluster hit) | t − marker | offset | words since its VMM's marker | word in frame |
|---|---|---|---|---|
| det3 | uniform 0 – 1.47 ms (median 734 µs) | uniform 0 – 15 | 18 (14 – 22) | uniform 0 – 1490 |
| det2 | same | same | 12 (10 – 16) | same |
| det1 | same | same | 7 (3 – 11) | same |

The time to the marker is just where in the 1.47 ms marker period the muon came: uniform, the same as for all hits
(`bcid_marker.png`). It carries no information. The "words since marker" differs by chamber only because the marker
round is written in VMM order (VMM 0 … 19): the hits come right after the round, and det3's marker is further back.

## 2. BCID differences between the three chambers

All three chambers' hits of a track use the **same marker** (99.995 %). The BCIDs are therefore directly comparable:

| pair | ΔBCID 68 % | ΔBCID 90 % | \|ΔBCID\| ≤ 1 | \|ΔBCID\| ≤ 3 | Δt (BCID + TDC) 68 % |
|---|---|---|---|---|---|
| det2 − det3 | −2 … +2 | −3 … +3 | 64 % | 95.5 % | −40 … +34 ns |
| det1 − det3 | −2 … +2 | −3 … +3 | 65 % | 95.6 % | −36 … +35 ns |
| det1 − det2 | −1 … +2 | −2 … +3 | 68 % | 97.0 % | −30 … +36 ns |

So a muon gives the same BCID to within ±2 BC (±45 ns, the drift-time spread), and to within ±4 BC (99 %). The
paddles come 10-11 BC (225-250 ns) **before** the chambers (44: 10 [8, 12], 45: 11 [9, 13], 41: 10 [9, 12],
5-95 %): this is the chamber latency, and they share the same marker as well.

## 3. Where the track hits sit in the hit list

- **The hits of a muon come out as one block.** The span from the first to the last track hit in the list is 2 hit
  words (median; 2-3 for 90 %, 4 at 99 %). With typically one hit per chamber, that means the three are
  **consecutive** hit words: 91.8 % of the tracks have no other hit in between at all.
- 99.97 % of the tracks lie in one UDP frame (the rest straddle two frames).
- Non-track hits inside the track span: another over-threshold chamber hit in 0.4 % of the tracks, a neighbour hit
  (flag 0) in 0.04 %, a hybrid-0 paddle hit in 7.8 %. The paddle hits are the trigger pulses of the same muon:
  they are read 2-3 words **before** det2 (5-95 %: −4 … +3).
- Order in the block: det3 → det2 → det1 in 66 % of the tracks, i.e. the VMM-ID order (hybrid 0 first, then
  VMM 2 … 19). The other orders (det2-det1-det3 14 %, det3-det1-det2 7 %, det2-det3-det1 7 %, …) are not
  understood yet. A guess to check: the VMMs are read in turns and the block straddles one turn.
- Reference with unrelated hits (the first hit of each chamber 5 µs later): span 6 (median), 2-24 for 90 %, 48 at
  99 %. So "nearby in the list" is not a strong muon signature by itself: the list is sparse, about 1.2 kHz of hits
  between 680 marker rounds per second.

## 4. The stream is not time-ordered, but muon blocks are

For all over-threshold chamber hits, "lateness" = (latest hit time already read) − (own time): 17 % of the hits are
read after a hit more than 1 µs later, 9 % > 5 µs, 1.3 % > 20 µs, 0.2 % > 90 µs. Within one VMM, 14-27 % of
consecutive hits go backwards in time (p99 about 1.2-2.6 µs). An online event builder that groups hits by time
therefore needs to look back about 20 µs (99 %), or sort per marker period (1.47 ms) before grouping. A muon's hits,
by contrast, come out together and in time order: for them a window of a few hit words in stream order is enough.

## 5. Per run

| run | straight tracks | per hour | in one frame | span median | lateness p99 |
|---|---|---|---|---|---|
| 181 K+64 | 883 | 294 | 100 % | 2 | 16 µs |
| 184 K+24 | 1689 | 563 | 99.9 % | 2 | 20 µs |
| 187 K0 | 2540 | 847 | 100 % | 2 | 24 µs |
| 190 K−8 | 2096 | 699 | 100 % | 2 | 29 µs |
| 193 K−12 | 1851 | 617 | 100 % | 2 | 34 µs |
| 196 K0 | 2479 | 828 | 99.96 % | 2 | 23 µs |
| 199 K+40 | 1998 | 666 | 100 % | 2 | 19 µs |
| 202 K+12 | 2850 | 950 | 100 % | 2 | 22 µs |
| 205 K−4 | 3046 | 1015 | 100 % | 2 | 20 µs |
| 208 K−10 | 2453 | 818 | 99.8 % | 2 | 54 µs |

(Per hour on the corrected 22.5 ns clock. The stream structure does not depend on the threshold. The lateness tail
is longest at the lowest thresholds, K−12 and K−10, where the noise rate is highest.)
