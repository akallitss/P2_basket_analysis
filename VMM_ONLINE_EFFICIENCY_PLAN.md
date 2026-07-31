# VMM online analysis — efficiency, maps, timing-vs-HV

**Written 2026-07-30.** Stages 0 and 1 are **done and running**; the rest is the
roadmap. Companions: `sps_beam_analysis/PIPELINE_OVERVIEW_2026-07-26.md` (Dream
side), `README_vmm_hybrid_pcapng_monitoring.md` (data format).

---

## 1. The setup

Same three stations as the Dream readout, same detectors, same z
(`run_config_beam.P2_VMM_CABLING`):

| Station | z (mm) | c4 (hybrid, bot, top) | c5 | c6 |
|---|---|---|---|---|
| P2_IN | 320 | (1, 2, 3) | (2, 4, 5) | (3, 6, 7) |
| P2_MID | 630 | (4, 8, 9) | (5, 10, 11) | (6, 12, 13) |
| P2_OUT | 940 | (7, 14, 15) | (8, 16, 17) | (9, 18, 19) |

18 VMMs = 3 stations × 3 connectors × 2 VMMs (128 ch per connector). Hybrid 0
(VMMs 0/1) is the external trigger digitiser. uRWELL references are read by
Dream, not VMM.

It is a telescope, so the SPS methods port. And because the trigger is digitised
**in the same stream**, there are two independent efficiency handles where Dream
had one — trigger-referenced (Stage 1, done) and tag-and-probe (Stage 4). They
have different systematics, so each validates the other.

## 2. What was built (Stages 0 + 1)

Two new modules, deployed to the live tree
`/local/p2/DAQ_Control_VMM_Beam/vmm_qa/` (which is the deployment the running
`qa_watcher` uses — **not** `~/DAQ_Control_VMM_Beam`, which is a stale checkout):

- **`vmm_stations.py`** — cabling, trigger, hot-channel masking, pcapng parser,
  marker-cadence diagnostics.
- **`vmm_efficiency.py`** — trigger-referenced station efficiency, plots, JSON.

Both run in the DAQ venv (no ROOT, **no scipy** — numpy fallbacks throughout,
with scipy used automatically offline). Nothing that the live QA runs was
modified, so the running watcher was never at risk.

### 2a. Trigger: VMM 0 ch 44, not VMM 1 ch 60

VMM 0 ch 44 holds **100.0%** of that VMM's hits; VMM 1 never appears at all.
`run_config_beam.py` recorded `{'vmm': 1, 'channel': 60}` from July notes, which
predate the 07-29 recabling. Corrected in place, with the old value kept as
`trigger_channel_july_note` for provenance. Backup at
`/tmp/run_config_beam.py.bak_ak`; **not committed** (that file has your own
uncommitted work in it).

Re-verify after any hybrid-0 change: `vmm_efficiency.py <file> --find-trigger`.

### 2b. Hot-channel masking, in analysis

Seeded with the known ones (VMM 7 ch 62, VMM 4 ch 58) plus an automatic flag at
8× the VMM's own median occupancy. On the test file it masked 6 102 hits from
`{4: [58, 61], 7: [37, 41, 43, 60, 62]}`. Trigger hits are never masked.

### 2c. The timing was aliased — and that is now solved

`offset`, the BCID rollover counter, is stuck at a constant, and SRS markers
arrive only every **1.6384 ms** while BCID wraps every **92.16 µs** — so
**17.8 wraps per marker interval**, and `abs_time_ns` is not even monotonic.
Coincidence-hunting on it directly returns noise: σ ≈ 300 ns, i.e. the entire
search window, with 17–55% accidentals.

What *is* unambiguous is the BCID phase and which marker interval a hit belongs
to. Markers are broadcast (2 799 of 2 808 values shared across VMMs), so the
working method is:

1. pair a trigger only with station hits carrying the **same `srs_timestamp`**;
2. compare **BCID phases wrapped into ±46.08 µs**;
3. subtract the accidental shape by **event mixing** (correlate each interval's
   triggers with the *next* interval's hits) — necessary because BCID phases are
   not uniform, so a flat-background fit is badly biased;
4. fit the **dominant** peak locally (a global width straddles the secondary
   peak of §2e and returns 300 ns for a 21 ns peak).

Effect: peak/background went 1.01 → **8.9–15.4**, accidentals 17–55% → **2–6%**.

### 2d. First efficiency numbers — and why they are NOT comparable to Dream's 96%

400 packets (487 k hits) of `run_24/nominal_05`, after trigger de-duplication
(§2e), which cut the denominator from 171 437 to 73 341:

| Station | µ (ns) | σ (ns) | peak/bg | raw | accidental | **efficiency** |
|---|---|---|---|---|---|---|
| P2_IN | +179.6 | 40.3 | 15.1 | 0.135 | 0.025 | **0.110** |
| P2_MID | +108.7 | 22.0 | 28.8 | 0.595 | 0.061 | **0.534** |
| P2_OUT | +108.7 | 23.4 | 24.8 | 0.165 | 0.025 | **0.139** |

MID and OUT share a latency to 0.0 ns; IN sits 71 ns later and is 2× wider.

**These are not the same quantity as the Dream 96%, and must not be compared to
it.** Dream's number comes from `22_tag_probe_efficiency`, which is *per pad*,
computed *only over the geometric overlap where the tagging planes actually
illuminate the probe* — its own docstring says "regions the tag never
illuminates are simply absent from the denominator". The number above is a
whole-station scalar whose denominator is **every trigger**, with no spatial
information at all. On stations where only 3 of ~10 connectors are instrumented
(c4–c6), most triggered particles land on pads that are not read out and are
counted as inefficient.

Evidence this is a denominator effect and not a sick detector:

- Requiring the other two stations to fire (temporal tag-and-probe, still no
  spatial matching) moves the numbers to IN 0.155 / MID 0.758 / OUT 0.180 —
  same ordering, still diluted by the probe's un-instrumented area.
- Only **1.7 %** of triggers fire all three stations, consistent with three
  partially-overlapping small acceptances.
- The coincidence peaks themselves are healthy: σ at the quantisation floor,
  peak/background 15–29, clean flat baseline.

**Unresolved:** why MID (0.53) is ~4× IN and OUT, when per-VMM occupancy shows
all three roughly centred on their instrumented region. Needs the pad map to
answer — see below.

**Consequence for planning: the pad map (Stage 2) is the critical path, not a
nice-to-have.** No number comparable to Dream's 96 % can be produced without it,
because the denominator has to be restricted to tracks that actually cross
instrumented pads.

### 2e. The trigger re-fires — found and fixed

The trigger channel's inter-arrival distribution peaks at **500–600 ns**
(p50 = 626 ns), and every station showed a matching secondary coincidence peak
~500 ns before the main one at 36–40 % of its height. Same particle, counted
several times. `TRIGGER_DEADTIME_NS = 1500` now collapses repeats into one
trigger; this raised efficiency by ~30 % relative and nearly doubled
peak/background (MID 15.4 → 28.8).

The underlying cause — discriminator ringing, a wide trigger pulse crossing
threshold in successive BCIDs, or something else — is still worth understanding
on the DAQ side, since the dead time is a workaround.

Confirmed at the same time that the trigger is genuinely beam: its rate shows
textbook SPS spill structure (~6 s on, ~10 s off), 1.4 kHz average / ~50 kHz
in-spill. Station cluster size is ~1 pad (mean 1.11–1.18).

### 2f. σ is pinned at the BCID bin

σ ≈ 22–23 ns is exactly the 22.5 ns BCID quantisation, so the TDC fine time is
contributing nothing. With per-channel time calibration σ should fall well below
one BCID. Same gap that gates Stage 5.

### 2g. Also observed

Corrupt VMM ids (id 27, not in the cabling) at 0.07% of hits — now counted and
reported rather than silently entering the histograms.

---

## 2A. Stage 2 (mapping) and Stage 4 wiring — done 2026-07-30

### 2A.1 A marker-parser bug was hiding a third of the detector

`vmm_pcapng_qa.py` inherits vmm-sdat's `ParserSRS` logic, which gates VMM marker
words on `vmmid_marker < 16` — upstream a FEC hosts at most 16 VMMs and ids ≥ 16
encode TRG trigger words instead. **This FEC hosts 20** (hybrids 0–9), so in SRS
mode every marker for VMMs 16–19 was thrown away. Those four had
`srs_timestamp = 0` on every hit, hence meaningless `abs_time_ns` in all online
QA plots and **exactly zero** trigger coincidences — 4 of P2_OUT's 6 VMMs, 115 613
hits, 37 % of all detector hits.

Fixed in both `vmm_stations.py` and `vmm_pcapng_qa.py` (SRS mode takes the id at
face value; TRG behaviour untouched). P2_OUT went from **0.139 → 0.560**, and it
now agrees with P2_MID (0.534) as two identical stations should.

### 2A.2 The pad map, and an independent confirmation of the flipped ribbon

`vmm_stations.build_pad_table()` implements
`(vmm, ch)` → station/connector/half → strip → `channel_id` → `pad_cx/pad_cy`,
reusing the M3-validated logic from `cosmic_bench_analysis/p2_mapping.py`. The
map CSV is deployed next to it so the beam path stays self-contained.

Both uncertain links were resolved from the data, using **trigger-coincident hits
only** (noise smears the spot and halves the discrimination):

- `vmm_map_scan.py` — global scan of 3 orderings × 2 half-assignments.
  `reverse` + bottom-VMM-reads-`bot` wins: RMS 76.9 mm vs 115–135 mm for the
  alternatives, and 59/61 mm for MID/OUT individually.
- `vmm_map_perhalf.py` — per-connector-half refinement. Every half is `reverse`
  **except `('P2_IN', 5, 'top') → 'linear'`** (59.7 mm vs 103.8 mm).

That override is the **same flipped `c_5_top` ribbon the Dream side found**
(commit `659dfbb`), recovered independently through a different readout chain —
strong evidence the whole mapping is right. Cross-check: P2_MID and P2_OUT beam
centroids now agree to ~1 mm ((424, 217) vs (425, 216)); P2_IN sits at (372, 193).

### 2A.2b Compactness was the wrong metric — speckle is the right one

Alexandra spotted from the hit maps that P2_OUT looked right while MID and IN
did not: **holes next to bright pads, and stripes**. The compactness scan had
missed this entirely, because permuting channels *inside* a 64-strip half leaves
the hits in the same part of the fan — the RMS barely moves. What a scrambled
ordering actually produces is high-frequency texture.

`vmm_stations.map_roughness()` measures that: each pad's count against the
median of its 6 nearest neighbours, over the illuminated region.
`vmm_map_rough.py` minimises it per half. Result:

| Station | roughness | verdict |
|---|---|---|
| P2_OUT | 0.057 | correct, no ordering change helps |
| P2_MID | 0.105 | no ordering change helps — see below |
| P2_IN | 0.239 → 0.184 | two halves fixed, **still wrong** |

It found a second override, `('P2_IN', 6, 'top') → 'linear'`, that the
compactness scan could not see.

**P2_MID's holes are not a mapping error.** Only **345 of 384** of its channels
ever fire, against 378/384 for OUT and 381/384 for IN — ~39 genuinely dead
channels, which is exactly what the metric says no ordering can fix.

**P2_IN is still visibly striped after the fix.** Roughness 0.184 vs OUT's
0.057, and the texture is still there by eye. The honest conclusion is that the
true permutation for those halves is **not in {linear, reverse, pairswap}** —
that set came from the DREAM/K59V adapter and the VMM hybrid need not share it.
Next step is a wider search (e.g. fit the permutation directly by maximising
smoothness over all 64 channel positions, or bit-level patterns like
`ch ^ k` / `ch >> n`), not another guess from the same three.

### 2A.2b-ii Neighbour-based hot-pad masking

With the pad map in hand, a much better hot-channel test became possible:
`vmm_stations.auto_hot_pads()` flags a pad holding more than 6× the median of
its **6 spatial neighbours**. Occupancy varies strongly across the fan, so the
old per-VMM median hid channels that were hot relative to the pads actually next
to them. Calibrated on run_25: the worst pad on the healthy P2_OUT is 4.7×, so
6× flags nothing there while catching the real outliers on IN (13–30×) and MID
(8–986×, including pads whose neighbours are entirely dead).

It caught P2_IN **VMM 5 ch 12** (1 046 hits vs neighbour median 80) which the
per-VMM test had missed entirely, plus several on P2_MID.

**Masking made the efficiency go DOWN, not up:**

| | before mask | after mask |
|---|---|---|
| run_24 P2_IN | 0.1108 | **0.1055** |
| run_25 P2_IN | 0.1244 | **0.1168** |
| run_24/25 P2_MID | 0.5400 / 0.5790 | 0.5401 / 0.5748 |
| run_24/25 P2_OUT | 0.5674 / 0.6397 | 0.5676 / 0.6397 |

The noisy pads had been contributing chance hits inside the coincidence window
and **inflating** P2_IN's efficiency. Removing them gives a cleaner and slightly
lower number, and better signal-to-background (P2_IN peak/bg 12.4 → 16.7 on
run_25, P2_MID 28.9 → 31.0). So noisy channels are **not** the cause of P2_IN's
deficit — which leaves the gain explanation of §2A.2c standing.

Also fixed: the hit maps were being drawn from **unmasked** hits while the
efficiency used masked ones, so the two products disagreed about what the
detector looked like. `vmm_beam_profile.make()` now applies the identical mask.

### 2A.2c Why P2_IN reads low — gain, not mapping

Note first that **the mapping cannot affect the efficiency at all**: the
trigger-referenced measurement uses only `(vmm, ch)` and time, never a pad
position. Re-running with the corrected map reproduced the efficiencies to four
decimals, as it must.

P2_IN reads low because its signals are small (run_25 `meshscan_m00V`):

| Station | hits | ADC p50 | ADC p90 | live channels |
|---|---|---|---|---|
| P2_IN | 43,906 | 101 | 199 | 381/384 |
| P2_MID | 173,885 | 112 | 246 | 345/384 |
| P2_OUT | 241,593 | 137 | 398 | 378/384 |

Nearly all its channels are alive, so this is not a connection or cabling fault
— the pulse-height spectrum is simply compressed towards threshold, ~2× lower
than P2_OUT. With a fixed threshold, hit yield falls much faster than gain,
which is the 4–5× hit deficit.

P2_IN runs at **mesh 440 V** against 450 V for the others, and the run_25 mesh
curve shows it still on the **steep rising part** of its turn-on at 440 V
(0.121) while MID/OUT have reached 0.58/0.64 at 450 V. So the Dream result
(`P2in_hvrange_2`, excellent efficiency) is not contradicted — Dream measured it
where it was properly powered. **The test is to raise P2_IN's mesh toward
450–460 V**, which is the same conclusion the 2026-07-26 logbook reached from
timing ("argues for trying P2_IN at 445–450 V post-repair").

### 2A.3 Hit maps and beam profiles

`vmm_beam_profile.py` renders the actual fan geometry (rotated pad tiles from
`pad_cx/pad_cy/pad_w/pad_h/pad_angle`), per capture file:

- `<base>_hitmap.png` — pad occupancy per station. Top row all hits, bottom row
  **trigger-coincident only**, i.e. the beam. Dead pads show as gaps.
- `<base>_beam_profile.png` — x and y projections of the in-time hits with
  centroid and RMS; the numbers also go into `events.json`.

Only the instrumented connectors (c4–c6) are drawn — the rest of each chamber is
absent by construction, which is exactly the acceptance limit behind §2d.

Beam centroids across the telescope (full file):

| Station | z (mm) | x (mm) | y (mm) |
|---|---|---|---|
| P2_IN | 320 | 414.3 ± 58.0 | 216.2 ± 36.4 |
| P2_MID | 630 | 422.3 ± 46.2 | 216.2 ± 36.9 |
| P2_OUT | 940 | 420.7 ± 50.1 | 211.8 ± 37.6 |

Agreement to ~8 mm in x and ~4 mm in y over a 620 mm lever arm — a beam
essentially parallel to the telescope axis. This is an **independent check on the
mapping**: a wrong wiring would not put three separately-mapped stations on the
same beam line.

### 2A.4 Wired into the watcher

`vmm_pcapng_qa.py` now runs `vmm_efficiency.analyse()` on the DataFrame it has
**already parsed** — no second parse — writes the two PNGs into the same
`--out-dir`, and folds the scalars into `events.json` (`n_triggers_dedup`,
`n_masked_hits`, and per-station `efficiency`/`mu_ns`/`sigma_ns`/`contrast`), so
sub-run aggregation is a cheap JSON scan. `qa_watcher` spawns a fresh process per
file, so **it picks this up automatically on the next capture** — no restart.

The whole stage is wrapped in a broad `try/except`: if efficiency fails, QA still
produces all its normal plots. New flags `--no-efficiency` and `--eff-window`.

**Cost, measured on a full file** (1 749 packets, 1.70 M hits):

| | wall | peak RSS |
|---|---|---|
| baseline (`--no-efficiency`) | 65.5 s | 1.36 GB |
| with efficiency | 74.1 s | 1.87 GB |

So **+8.6 s (+13 %) and +0.5 GB**. Two honest caveats: QA was *already* slower
than the 44.4 s capture rotation before this change (65.5 s), so any backlog
pre-exists and is not caused by the efficiency stage; and the machine has only
7 GB total with ~4 GB free, so 1.87 GB peak is comfortable but not negligible
against the watcher's 80 % memory kill.

Full-file numbers are stable against the 400-packet sample (0.111 / 0.540 /
0.567 vs 0.110 / 0.534 / 0.560).

## 3. What ports from the SPS pipeline

| SPS stage | Ports? | Notes |
|---|---|---|
| 20_beam_spectra (Landau MPV vs HV) | Yes | needs clustering + mapping |
| 21_telescope_align | Yes | 3 planes, known z |
| 22_tag_probe_efficiency | Yes | the main port; cross-checks Stage 1 |
| 23_beam_profile | Yes | |
| 24_event_sync_qa | Adapt | FEU sync → VMM/hybrid marker sync |
| 25_commissioning_qa | Yes | most beam-useful item |
| 26_hv_spark_qa | Yes | needs the VMM-side HV csv |
| 27_pedestal_qa | Adapt | VMM pedestals, not Dream prg files |
| 28_timing_qa | Yes | BCID+TDC replaces the Dream `time` branch |
| 29_waveform_timing | No | VMM ships no 16-sample waveforms |
| 30_raw_stream_efficiency | Superseded | Stage 1 is the better cross-check here |

---

## 4. Still to do

<details>
<summary><b>Next up — explain the secondary peak (blocks quoting any efficiency)</b></summary>

- Histogram the trigger channel's own hit-to-hit spacing: if VMM 0 ch 44 fires
  twice ~500 ns apart, that is the whole story and the fix is to de-duplicate
  triggers within a dead time.
- Check whether the secondary peak's share tracks beam intensity (accidental) or
  is constant (instrumental).
- Compare the ADC spectrum of hits in the main vs secondary peak — a re-fire
  should look different from a real particle.
- If it is a trigger re-fire, add a configurable dead time to
  `vs.trigger_times()` and re-run; efficiency will rise.
</details>

<details>
<summary><b>Stage 1 completion — from one file to the whole campaign</b></summary>

- Run over full files (drop `--max-packets`) and check stability file to file.
- Wire into `qa_watcher`: call `vmm_efficiency.analyse()` from
  `vmm_pcapng_qa.py` using the DataFrame it has **already parsed** (no second
  parse), writing the two PNGs into the same `--out-dir`. They then appear in
  the GUI gallery automatically — `/list_pngs` is filename-driven, so no GUI
  work is needed.
- Extend `events.json` with `eff_<station>`, `mu_ns`, `sigma_ns`, `n_masked` so
  aggregation is a cheap JSON scan.
- `vmm_qa/aggregate_subrun.py`, fired by the watcher at sub-run end → 
  `analysis/<run>/<sub_run>/_summary/`, giving **efficiency vs HV** across a
  scan. This is the same gap as SPS improvement #2.
- Backfill the whole campaign on EOS/condor (you offered — this is the step
  where it pays off; one job per pcapng, results merged from the JSONs).
- Refactor: make `vmm_pcapng_qa.py` import `parse_pcapng()` from
  `vmm_stations.py` instead of keeping its own copy. Deliberately deferred —
  that script is what the live watcher runs. Do it between runs.
</details>

<details>
<summary><b>Stage 2 — pad mapping ✅ DONE (see §2A.2)</b></summary>

`reverse` + bottom-VMM-reads-`bot`, with one override
(`('P2_IN', 5, 'top') → 'linear'`). Independently confirms the Dream-side
flipped-ribbon fix. Re-run `vmm_map_perhalf.py` after any recabling.

Remaining: the per-pad **efficiency map** is now possible
(`vmm_efficiency.coincident_mask()` gives the in-time hit selection) but is not
yet plotted — that is the natural next product.
</details>

<details>
<summary><b>Stage 3 — clustering</b></summary>

- Event building on marker interval + BCID phase (reuse Stage 1's grouping),
  then spatial clustering (mirror `sps_cluster.py`: leading pad + pads within a
  radius).
- Products: cluster size; **cluster charge + Landau MPV** (= SPS 20, the online
  gain monitor — arguably the most useful thing during an HV scan); beam profile
  in mm (= SPS 23).
- Also gives the **per-pad efficiency map** for Stage 1, which currently
  produces only a per-station scalar.
</details>

<details>
<summary><b>Stage 4 — alignment + tag-and-probe (the SPS 21 → 22 port)</b></summary>

- Per-pad efficiency maps and efficiency-vs-HV curves.
- Cross-check against Stage 1. Where the two disagree, the disagreement is
  itself the measurement — different systematics, same detector.
- Carry over SPS 22's caveat verbatim: efficiency relative to the tag selection
  and its acceptance, not absolute.
</details>

<details>
<summary><b>Stage 5 — timing resolution vs drift/mesh (blocked)</b></summary>

- **Blocked: there are no time-calibration runs.** `qa_config.py` has
  `'calibration': None`, and `load_calibration()` reads `timewalk_a/b/c/d` but
  explicitly does not apply them.
- Evidence it matters: σ is pinned at the 22.5 ns BCID bin (§2e), so the TDC
  fine time is currently contributing nothing.
- SPS 29 found per-FEU clock phases of −9.6/−6.6/−5.5 ns per unit and called the
  per-station fit load-bearing; the VMM analogue is per-VMM/per-channel offsets.
- Needs a calibration campaign first, then the σ-vs-drift / σ-vs-mesh curves
  follow the SPS 29 shape.
</details>

<details>
<summary><b>DAQ-side fixes that would sharpen everything (worth raising)</b></summary>

- **Restore the `offset` field** (BCID rollover counter), or **raise the SRS
  marker rate above ~11 kHz** (period < 92.16 µs). Either removes the 17.8-fold
  ambiguity entirely, which would cut the accidental background by ~18× and let
  the analysis use plain `abs_time_ns`.
- Fix `hit_valid` in `vmm_pcapng_qa.py`: it currently flags 100% of hits as
  anomalous, because it inherits the ESS SRS spec's valid offset range while
  this firmware pins `offset` at a constant. As a QA flag it is pure noise.
- Consider adding **scipy** to the DAQ venv — exact Clopper–Pearson and proper
  least-squares fitting instead of the numpy fallbacks.
</details>

## 5. How to run it

```bash
cd /local/p2/DAQ_Control_VMM_Beam/vmm_qa
nice -n 19 ../.venv/bin/python vmm_efficiency.py <file.pcapng> \
     --out-dir <dir> --json            # add --max-packets N for a quick look
nice -n 19 ../.venv/bin/python vmm_efficiency.py <file.pcapng> --find-trigger
```

Products: `<base>_trigger_dt.png` (Δt per station, signal + sideband windows,
event-mixing subtracted), `<base>_trigger_efficiency.png`, `efficiency.json`.
