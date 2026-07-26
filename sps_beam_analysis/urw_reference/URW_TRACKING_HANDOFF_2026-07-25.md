# uRWELL tracking handoff — TB_July2026_H4

**Written 2026-07-25. Revised 2026-07-26. For someone who has never seen this
data, this DAQ, or these detectors.**

**The goal:** reconstruct straight tracks between the two EIC uRWELL reference
detectors, then use them to measure the three P2 BASKET stations sitting
between them. The uRWELLs are 1370 mm apart. The beam particles come in nearly
straight, so a translation-only alignment is enough — no rotations are
expected, and the data confirms that (§8).

> **What changed on 2026-07-26.** The back uRWELL's channel → strip wiring was
> wrong (§6.2): the fix takes it from 4.4 mm pointing to ~0.9 mm. Every back
> number in §8 and §9 moved as a result, and the "4.5 mrad beam divergence" it
> implied was an artefact. The P2 work (§13) is new and is the point of the
> exercise.

Everything below has been run on this machine. Numbers quoted are real outputs
you should be able to reproduce, not estimates. Where something is unverified I
say so explicitly.

---

## 0. Ten-minute quickstart

```bash
# 1. The environment. BOTH lines matter.
unset PYTHONPATH                                  # see §2, this is not optional
PY=/local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python

# 2. Run the worked example: align the two uRWELLs and build tracks.
cd /local/home/banco/P2_basket_analysis/sps_beam_analysis/urw_reference
$PY align_and_track.py --min-amp 0 --plot align.png

# 3. The actual measurement: P2 residuals and efficiency against those tracks.
#    (~25 min over all six sub-runs; add --max-chunks 1 for a 3-minute look.)
$PY urw_p2_efficiency.py --run highstat_eff_1 --out out_eff
```

Step 2 prints the alignment constants and the track angles and writes
`align.png`; §9 explains every number. Step 3 writes one eight-panel PNG per
station per sub-run plus a summary; §13 explains those. Then read §4 (the strip
maps) and §6 (the channel → strip wiring), because those two sections are where
all the traps are.

---

## 1. The setup

Five detector planes on the H4 beamline at the SPS, all perpendicular to the
beam. The beam travels along **z**. The two uRWELLs are the *references* — the
things we trust to define a track. The three P2 BASKET stations in the middle
are the *devices under test*.

```
     beam  ------------------------------------------------------------>  z
             |            |             |             |               |
             |            |             |             |               |
        uRWELL_front    P2_IN        P2_MID        P2_OUT       uRWELL_back
          z =    0      z = 320      z = 630      z = 940       z = 1370   [mm]
          "urw_inter"                                            "urw_strip"
          FEU 1                                                  FEU 1
          conn 1,2,3,4    FEU 3        FEU 4        FEU 5         conn 5,6,7,8
```

| plane | `det_type` | z [mm] | FEU | FEU Dream connectors |
|---|---|---|---|---|
| `EIC_uRWELL_front` | `urw_inter` | 0 | 1 | 1, 2, 3, 4 |
| `P2_IN` | `P2` | 320 | 3 | 4, 5, 6, 8 |
| `P2_MID` | `P2` | 630 | 4 | 4, 5, 6, 7 |
| `P2_OUT` | `P2` | 940 | 5 | 4, 5, 6, 7 |
| `EIC_uRWELL_back` | `urw_strip` | 1370 | 1 | 5, 6, 7, 8 |

**The single most useful fact in this document:** *both uRWELLs are read out by
the same FEU (FEU 1)*. Front and back hits are therefore already in the same
file and already share an event number. For uRWELL-to-uRWELL tracking you never
have to synchronise anything. You only need the cross-FEU merged files (§5.4)
when you bring the P2 stations in.

Run conditions (from `run_config.json`): gas Ar/Iso 95/5, trigger "SPS external
scintillator coincidence via TCM", both uRWELLs held at a fixed operating point
(drift 600 V, resistive 420 V) for the entire campaign.

---

## 2. Environment — read this or nothing will work

```bash
unset PYTHONPATH
PY=/local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python
```

The `banco` login shell puts an ISEG high-voltage SDK on `PYTHONPATH`. That SDK
ships its own partial copies of modules that shadow ROOT and uproot, and the
failures are confusing (import errors deep inside uproot, or a `uproot` that
loads but cannot read a tree). `unset PYTHONPATH` in every shell before you run
anything. If a script suddenly breaks, check this first.

Use that interpreter, not the system `python3`. It is the only one with uproot,
awkward, pandas, matplotlib and ROOT all working together.

---

## 3. What a uRWELL is, in one paragraph

A micro-pattern gas detector. A charged particle crosses a thin gas gap, ionises
the gas, the ionisation is amplified in a micro-well structure, and the charge
lands on a plane of readout **strips**. The strips are one-dimensional: a strip
tells you the coordinate *across* it, and nothing about the coordinate *along*
it. So each detector has two independent sets of strips — one set running
vertically that measures **x**, one set running horizontally that measures **y**
— stacked in the same plane. One particle therefore produces a signal on a few
adjacent x-strips *and* a few adjacent y-strips, and the crossing point of the
two gives you a 2D point.

Both of our uRWELLs have a ~127 × 127 mm active area and 128 strips per view,
i.e. 256 strips per detector, 512 channels for the pair — which is exactly the
8 × 64 channels of one FEU.

The two are **not** identical. They are test structures with deliberately
different strip geometry, which is the subject of the next section.

---

## 4. The strip maps — the confusing part

The strip geometry lives in two plain CSV files:

```
/local/home/banco/P2_data/TB_July2026_H4/config/detectors/inter_map.txt   <- front (urw_inter)
/local/home/banco/P2_data/TB_July2026_H4/config/detectors/strip_map.txt   <- back  (urw_strip)
```

Each has 256 rows (one per strip) and these columns:

```
connector,connectorChannel,stripNb,axis,pitch(mm),interpitch(mm),neighbours(:separated),xGerber,yGerber
0,0,0,y,1,0.9:0.75,1,-113.875,-13.45
```

| column | meaning |
|---|---|
| `connector` | which 64-channel Dream connector, **0-based**, 0–3 |
| `connectorChannel` | 0–63 within that connector |
| `stripNb` | strip index within the view, 0–127 |
| `axis` | **the direction the strip RUNS** — see the warning below |
| `pitch(mm)` | centre-to-centre spacing to the next strip |
| `interpitch(mm)` | gap between strips (a design parameter of the test structure) |
| `neighbours` | adjacent `stripNb`s |
| `xGerber`, `yGerber` | strip position in the PCB design frame, mm |

### 4.1 The `axis` trap

> **`axis` is the direction the strip runs, NOT the coordinate it measures.**
> A strip with `axis='y'` runs along y, so it measures **x**.
> A strip with `axis='x'` runs along x, so it measures **y**.

This inverts once and it is very easy to get backwards. Picture it:

```
        axis = 'y'  strips                     axis = 'x' strips
        (run vertically, measure X)            (run horizontally, measure Y)

    y                                      y
    ^  | | | | | | | | | | | |             ^  ═══════════════════════
    |  | | | | | | | | | | | |             |  ═══════════════════════
    |  | | | | | | | | | | | |             |  ═══════════════════════
    |  | | | | | | | | | | | |             |  ═══════════════════════
    |  | | | | | | | | | | | |             |  ═══════════════════════
    +---------------------------> x        +---------------------------> x
       ^                                        the strip's position
       the strip's position                     varies along Y
       varies along X                           -> it measures y
       -> it measures x
```

And the corresponding column you read for the position:

* `axis='y'` → the strip measures x → use **`xGerber`**
* `axis='x'` → the strip measures y → use **`yGerber`**

(The other Gerber column is constant for all strips in that view — it is just
the centre of the strip, which carries no information.)

The library does this for you:

```python
meas   = 'x' if row['axis'] == 'y' else 'y'
gerber = row['xs_gerber'] if meas == 'x' else row['ys_gerber']
```

### 4.2 Positions are shifted to a local frame

The `xGerber`/`yGerber` values are in the PCB design frame and run from about
−113.9 to +13.6 mm. The detector loader shifts each view so its minimum is 0.
Everything downstream — my code, the numbers in this document — is in this
**local frame, 0 to ~127 mm, origin at the low-coordinate corner of the active
area**. Convert to the global beamline frame only at the very end, and only if
you need it (`UrwGeometry.to_global`).

### 4.3 The FRONT detector (`urw_inter`) — uniform pitch, two interpitch halves

Every strip is 1.0 mm pitch, both views. What varies is the *interpitch* (the
gap between strips), and only in the x-measuring view, which is split into two
halves:

```
  FRONT (urw_inter), local coordinates, 127.1 x 127.0 mm
  ┌───────────────────────────────┬───────────────────────────────┐
  │                               │                               │
  │   x-strips 0-63               │   x-strips 64-127             │  y-strips:
  │   pitch      1.0 mm           │   pitch      1.0 mm           │  all 128,
  │   interpitch 0.9 : 0.75       │   interpitch 0.67 : 0.5       │  pitch 1.0,
  │                               │                               │  interpitch 0.1
  │   map connector 0             │   map connector 1             │  uniform
  │                               │                               │
  └───────────────────────────────┴───────────────────────────────┘
  x = 0                        x = 63.0                     x = 127.1
```

So: **uniform 1 mm pitch everywhere; the left and right halves differ only in
inter-strip gap.** For tracking, the front detector is the simple one.

### 4.4 The BACK detector (`urw_strip`) — three pitch zones per view

This is the one that surprises people. *Both* views are divided into three
zones with **different strip pitch**, and the zone boundaries fall at the same
strip numbers in x and y. The result is a 3 × 3 patchwork:

```
  BACK (urw_strip), local coordinates, 127.4 x 127.3 mm

  y=127.3 ┌─────────────────────┬──────────────┬──────┐
          │                     │              │      │   y-strips 96-127
          │  x:1.0   y:0.5      │ x:1.5  y:0.5 │ 0.5  │   pitch 0.5 mm
          │                     │              │ 0.5  │   (a narrow 15.5 mm band)
  y=111.8 ├─────────────────────┼──────────────┼──────┤
          │                     │              │      │
          │                     │              │      │   y-strips 64-95
          │  x:1.0   y:1.5      │ x:1.5  y:1.5 │ x:0.5│   pitch 1.5 mm
          │                     │              │ y:1.5│
          │                     │              │      │
  y=64.3  ├─────────────────────┼──────────────┼──────┤
          │                     │              │      │
          │                     │              │      │
          │   ***  BEAM  ***    │              │      │   y-strips 0-63
          │                     │              │      │   pitch 1.0 mm
          │  x:1.0   y:1.0      │ x:1.5  y:1.0 │ x:0.5│
          │                     │              │ y:1.0│
  y=0     └─────────────────────┴──────────────┴──────┘
          x=0                x=63.0        x=110.6  x=127.4

            x-strips 0-63       x-strips      x-strips 96-127
            pitch 1.0 mm        64-95         pitch 0.5 mm
                                pitch 1.5 mm  (narrow 15.5 mm band)
```

Exact zone table:

| view | `stripNb` | pitch | interpitch | local coordinate range | map connector |
|---|---|---|---|---|---|
| back x | 0–63 | 1.0 | 0.75 | 0 → 63.00 | 0 |
| back x | 64–95 | 1.5 | 1.12 | 64.06 → 110.56 | 1 |
| back x | 96–127 | 0.5 | 0.37 | 111.94 → 127.44 | 1 |
| back y | 0–63 | 1.0 | 0.1 | 0 → 63.00 | 2 |
| back y | 64–95 | 1.5 | 0.1 | 64.25 → 110.75 | 3 |
| back y | 96–127 | 0.5 | 0.1 | 111.75 → 127.25 | 3 |

Two consequences you must not forget:

1. **A cluster's width in *millimetres* is not its width in *strips*.** Always
   cluster in position space with a pitch-aware tolerance, never "channel ± 1".
   The library uses `gap <= 1.05 * max(pitch_i, pitch_j)`, which correctly
   merges across a zone boundary without over-merging inside the 1.5 mm zone.
2. **The spatial resolution is not uniform across the plane.** If you quote a
   residual for the back detector, say which zone it came from. The 0.5 mm band
   is only 15.5 mm wide and is at the far edge — it gets very few tracks (~3 000
   clusters against ~50 000 in the 1.0 mm zone in the run I checked).

The beam lands mostly in the 1.0 mm × 1.0 mm corner (beam centre is around
x ≈ 61, y ≈ 51 mm), spilling into the 1.5 mm zones. That is good news: most of
your tracks are in the uniform region.

---

## 5. Where the data is and how it is structured

### 5.1 Directory layout

```
/local/home/banco/P2_data/TB_July2026_H4/
├── config/detectors/           <- strip maps + detector json (inter_map.txt, strip_map.txt, ...)
├── pedestals/                  <- pedestal runs
└── runs/
    ├── drift_mesh_scan_1/      <- 23 sub-runs. USE THIS ONE (fully processed, verified)
    │   ├── run_config.json     <- geometry, cabling, HV, z positions. Read by all the tools.
    │   ├── nominal_00/         <- one sub-run
    │   │   ├── raw_daq_data/       *.fdf     raw DAQ, one file per FEU per chunk
    │   │   ├── decoded_root/       *.root    unpacked waveforms (zero-suppressed)
    │   │   ├── hits_root/          *_hits.root      <- FITTED HITS, one file per FEU
    │   │   ├── combined_hits_root/ *_feu-combined_hits.root  <- ALL FEUS MERGED
    │   │   ├── hv_monitor.csv
    │   │   ├── recorded_events.npz   <- trigger bookkeeping, see §5.5
    │   │   └── .subrun_complete      <- present only when the sub-run finished
    │   ├── drift_450/ ... drift_900/          (10 sub-runs, P2 drift scan)
    │   └── meshscan_01_midout445/ ...         (12 sub-runs, P2 mesh scan)
    ├── env_test_1/
    └── highstat_eff_1/         <- newest, high statistics; sub-runs beam_commissioning_00..04
```

**Which run should you use?** Start with `drift_mesh_scan_1/nominal_00`. It is
complete, fully processed, and every number in this document comes from it. The
scan sub-runs vary *P2* voltages — **the uRWELLs were at a fixed operating point
for the whole campaign**, so as far as the uRWELLs are concerned all 23 sub-runs
are the same conditions and can be summed for statistics. `highstat_eff_1` is
the high-statistics run and is the one to use for a final efficiency number.

### 5.2 File naming

```
EicP2Bt_nominal_00_datrun_260725_02H50_000_01_hits.root
   |        |         |       |      |    |   |
   |        |         |       |      |    |   └─ FEU number (01, 03, 04, 05)
   |        |         |       |      |    └───── chunk number (000, 001, ...)
   |        |         |       |      └────────── start time 02H50
   |        |         |       └───────────────── date 2026-07-25
   |        |         └───────────────────────── "datrun" = data (vs "pedthr" = pedestal)
   |        └─────────────────────────────────── sub-run name
   └──────────────────────────────────────────── campaign prefix
```

A sub-run is split into **chunks** of a few million hits. There are several
chunk files per FEU. **Read all of them** — a single chunk is only part of the
sub-run.

### 5.3 The `hits_root` files — one FEU each

Produced by a C++ matched-filter analyser that fits each channel's waveform.
One entry = **one strip that fired in one event** (not one event).

> ⚠️ **These files no longer exist everywhere.** As of 2026-07-26 the per-FEU
> `hits_root/` directories have been deleted from **all 23 sub-runs of
> `drift_mesh_scan_1`** (and one sub-run each of `drift_mesh_2d_1` and
> `low_mesh_scan_1`); `highstat_eff_1` still has all six. `urw_lib.feu_hit_files`
> therefore falls back to `combined_hits_root/` (§5.4) when it finds nothing,
> and those files hold **all four FEUs in one tree**. When that fallback fires
> you *must* pass `feu=` to `iter_hits` as well — otherwise the P2 channels get
> mapped onto uRWELL strips and you get plausible-looking garbage with no error
> message. Every script in `urw_analysis/` now does this; if you write your own,
> do not forget it. The two sources agree exactly: `feu == 1` of the combined
> file is byte-for-byte the old `*_01_hits.root` content, and re-running the
> alignment through the fallback reproduces the per-FEU numbers to 4 decimals.

Tree `hits`:

| branch | type | meaning |
|---|---|---|
| `eventId` | uint64 | trigger number (see §5.5 — this is the join key) |
| `trigger_timestamp_ns` | uint64 | wall-clock timestamp of the trigger |
| `channel` | uint16 | **global channel 0–511** within the FEU |
| `amplitude` | float | fitted pulse amplitude — your "charge" |
| `time_of_max` | float | time of the pulse maximum, ns. In-time signal ≈ 380 ns |
| `max_sample` | float | interpolated sample index at the maximum |
| `significance` | float | amplitude / noise. **Hard-thresholded at 5 by the analyser** |
| `saturated` | bool | ADC saturated |
| `integral`, `time_over_threshold`, `local_baseline`, `local_max`, `sample`, `left_sample`, `right_sample`, `time`, `trunc_left`, `trunc_right` | float/bool | other fit outputs, not needed for tracking |

There is also a small `pedestals` tree (512 entries: `channel`, `mean`, `rms`,
`rms_gate`) — the per-channel noise, useful for spotting dead or hot channels.

> **Gotcha — `significance` is already cut.** Its minimum value in the file is
> exactly 5.00. Adding `--min-signif 5` does nothing. Do not mistake that for
> "my cut had no effect on the data".

> **Gotcha — ROOT cycles.** These files contain several *cycles* of the same
> tree (`hits;7`, `hits;8`) from reprocessing. `uproot.open(f)['hits']` picks the
> highest cycle, which is the newest and correct one. **Never add cycles
> together** — you would double-count. If you use bare ROOT/PyROOT, be explicit.

### 5.4 The `combined_hits_root` files — all FEUs, event-synchronised

Same tree, same branches, **plus one**:

| branch | type | meaning |
|---|---|---|
| `feu` | int32 | which FEU this hit came from: 1, 3, 4 or 5 |

This is the file to use when you want uRWELL **and** P2 information in the same
event. Map `feu` → detector using the table in §1: `feu == 1` is both uRWELLs,
`feu == 4` is P2_MID, and so on.

Structure of one such file (`nominal_00`, chunk 000):

```
16,765,346 entries total
  feu 1: 6,210,307 hits    <- both uRWELLs
  feu 3: 2,805,044 hits    <- P2_IN
  feu 4: 3,542,457 hits    <- P2_MID
  feu 5: 4,207,538 hits    <- P2_OUT
```

> **Gotcha — the tree is ordered by FEU, not by event.** All of FEU 1 comes
> first, then all of FEU 3, and so on. If you read only the first N entries as a
> quick test you will get *only FEU 1* and conclude the P2 stations are empty.
> Iterate over the whole tree, or filter on `feu` after reading everything.

For **uRWELL-only tracking you do not need these files at all** — `feu == 1` in
the combined file is byte-for-byte the same 6,210,307 hits as the
`*_01_hits.root` file, which is a third of the size. Use `hits_root` for §7–§9
and switch to `combined_hits_root` when you add the P2 planes.

### 5.5 `eventId` — the join key

This is what lets you say "the front hit and the back hit are the same particle".

Verified properties on `nominal_00`:

* **Continuous across chunks.** Chunk 000 covers eventId 1 … 2,395,958; chunk
  001 continues 2,395,959 … 2,795,239. It does **not** restart per file, so you
  can concatenate chunks and merge on `eventId` safely.
* **Common across FEUs.** All four FEUs use the same numbering, 1 … 2,795,239.
  So merging uRWELL and P2 hits is a plain join on `eventId`.
* **Not every event appears in every file.** An event only shows up if some
  channel exceeded threshold. In `nominal_00`, FEU 1 has hits in 2,081,959 of
  the 2,795,239 triggers (74 %).

`recorded_events.npz` in each sub-run directory is the bookkeeping:

```python
d = np.load('recorded_events.npz', allow_pickle=True)
d['feu1_range']    # array([      1, 2795239])
d['feu1_n']        # 2795239
d['feu1_missing']  # array([], dtype=...)   <- events the FEU dropped
```

**Check `feuN_missing` before trusting a sub-run.** It is empty for `nominal_00`,
but it is *not* empty everywhere — FEU 1 dropped ~632 000 of 2.61 M events in
the `drift_450` sub-run. Missing events are a readout problem, not a physics
one, but they will skew any efficiency you compute.

---

## 6. From channel number to strip position — including the trap

### 6.1 The channel numbering

Each FEU has 8 Dream connectors of 64 channels. The global channel in the file is

```
channel = (FEU connector number - 1) * 64 + connectorChannel        # 0 .. 511
```

which lays the two detectors out like this:

```
 channel:  0        64       128      192      256      320      384      448      511
           |--------|--------|--------|--------|--------|--------|--------|--------|
 FEU conn:     1        2        3        4        5        6        7        8
 detector: <-------- uRWELL FRONT -------->  <--------- uRWELL BACK --------->
 view:     <-- front x --> <-- front y -->   <-- back x --->  <--- back y --->
```

So: front x = channels 0–127, front y = 128–255, back x = 256–383,
back y = 384–511. Each **view** is 128 channels = **two** Dream connectors.

Note the map file's `connector` column is **0-based** (0–3) while the FEU
connectors in `run_config.json` are **1-based** (1–8). Map connector 0 → FEU
connector 1 for the front, → FEU connector 5 for the back.

### 6.2 THE CHANNEL → STRIP WIRING — read this before you plot anything

> **None of the four views is wired the way the strip map file says.** If you
> build the mapping naively from the map you will get wrong positions, and the
> wrongness is subtle enough to look plausible.

A view is read out on **two** 64-channel Dream connectors, so there are two
independent binary choices — which connector carries the low strips, and which
way the channel order runs inside a connector — hence four candidate wirings:

| mode | meaning |
|---|---|
| `AB` | map connector 0 → the lower FEU connector, channel order as-is |
| `BA` | the two connectors interchanged |
| `AB_rev` | as `AB`, but channel order **reversed inside each connector** |
| `BA_rev` | as `BA`, but reversed inside each connector |

The measured wiring (`urw_lib.VIEW_MODE_DEFAULT`):

| view | mode |
|---|---|
| front x | `BA` |
| front y | `AB` |
| back x | `AB_rev` |
| back y | `AB_rev` |

This is equivalent to: **every connector's channel order is reversed** — which
is exactly what `dream_feu_orientation: "inverted"` in `run_config.json` says —
**and in addition front y has its connector pair interchanged.** (For the front,
whose map is uniform, `BA`/`AB` and `AB_rev`/`BA_rev` are the same thing up to a
mirror of the axis, so the table above and "all inverted, front y also swapped"
describe the same detector. See §6.3.)

> ⚠️ **Corrected 2026-07-26.** Until then the code applied `BA` to back x and
> back y ("the connectors are interchanged"). That is the wrong member of the
> pair: it left the back pointing at **4.4 mm** instead of ~0.9 mm, which made
> the back look useless for tracking and is why §8/§9 below concluded the beam
> had a 4.4 mrad divergence. Redo any back-detector number produced before that
> date. The front was and is correct.

Why the connector order matters, and why it is easy to miss — the beam sits near
the middle of the plane, so getting the two halves the wrong way round **cuts the
beam spot in half and throws the two pieces to opposite edges**:

```
   WRONG (as the map says)                  RIGHT (connectors swapped)

   counts                                   counts
     |█                           █|          |            ████
     |██                         ██|          |          ████████
     |███                       ███|          |         ██████████
     |████                     ████|          |        ████████████
     +-----------------------------+          +-----------------------------+
     0          x [mm]          127           0          x [mm]          127
     "two blobs at the edges,                 "one blob in the middle"
      a hole in the middle"
```

The *connector order* was fixed first, with a test that needs no external
reference: for the FRONT, whose map is uniform (strip position increases
monotonically with channel index across the whole 0–127 of the view), swapping
the connectors is exactly a cyclic shift by half the plane, so the right choice
is the one that produces **one compact blob**. Fraction of the profile in a
single contiguous blob, as-is → connectors interchanged:

| view | as-is | interchanged |
|---|---|---|
| front x | 0.39 | **0.98** |
| front y | **0.998** | 0.62 |

(`python check_view_transform.py` — but read its header: it is superseded for the
back, see below.)

The *within-connector* direction is invisible to that test, because for a
uniform map reversing the channel order inside both connectors is the same as
mirroring the whole axis. The BACK's map is **not** uniform — three pitch zones
per view, §4.4 — so there the four modes are genuinely different, and the test
above silently picked the wrong one. The way to settle it is to use the front,
now known good, as the reference: the beam is parallel to well under a mrad, so
the back must reproduce the front position up to a constant. Width of
`back − front` per candidate (`python explore4_back_map.py`):

| view | `AB` | `BA` | `AB_rev` | `BA_rev` |
|---|---|---|---|---|
| back x | 45.6 mm | 5.42 mm | **0.80 mm** | 47.0 mm |
| back y | 35.2 mm | 6.40 mm | **0.88 mm** | 38.8 mm |

The full write-up is `~/P2_basket_analysis/sps_beam_analysis/urw_reference/ORDERING.md` (read the box at the
top — the body of that file predates this correction). **You do not have to do
anything about this** — `urw_lib.py` applies it automatically via
`VIEW_MODE_DEFAULT`. Just do not "fix" it by removing that.

### 6.3 Two conventions that are choices, not measurements

Two things in `urw_lib.py` are *labelling*, and it is worth knowing which:

1. **The mirror ambiguity.** No amount of uRWELL data can distinguish a view's
   wiring from its mirror image (`BA` ↔ `AB_rev` for a uniform map): both put the
   beam in the same place, one just counts strips the other way. What breaks it
   is the P2 pad map, which is an absolute physical frame: the uRWELL → P2
   transform has to be a **proper rotation** (both detectors are seen from the
   same side along the same beam, so a reflection is impossible). With the table
   in §6.2 it comes out as a −60° rotation with det = +1. The mirror partner
   would need det = −1. That is the only reason the front's entry is `BA`/`AB`
   rather than `AB_rev`/`BA_rev`.
2. **`AXIS_FLIP_DEFAULT`.** With the correct wiring the back reads
   *anti-parallel* to the front on both views (fitted front→back slope ≈ −1).
   That is physical, but it makes every downstream "the slope should be +1" check
   awkward, so the back's two axes are mirrored in software. This changes no
   strip → position assignment, only the sign of the axis label.

Finally, `run_config.json`'s `dream_feu_orientation` (now `"inverted"` on all
eight uRWELL connectors) records the plug orientation — the permutation *within*
a 64-channel connector. As of 2026-07-26 that field and the measured wiring in
§6.2 **agree**; no analysis code reads it, but it is no longer independent of the
mapping the way an earlier version of this handoff claimed.

---

## 7. From hits to points — clustering

One particle fires several adjacent strips. Turn them into one position:

1. Take all hits in one event, in one view.
2. Sort by strip position.
3. Start a new cluster whenever the gap to the previous strip exceeds
   `1.05 × max(pitch_i, pitch_j)` — pitch-aware, so it works across the back
   detector's zone boundaries.
4. Cluster position = **amplitude-weighted centroid**;
   cluster charge = sum of amplitudes.
5. Per event and view, keep the **highest-charge** cluster as *the* point.

`urw_lib.cluster_hits()` does all of this vectorised; `urw_lib.leading_points()`
does step 5 and gives you one `(x, y)` row per event.

**Do not apply an extra amplitude cut for tracking.** The analyser's
significance > 5 is already the real threshold. I measured this directly on
`nominal_00`:

| cut | events with all four coordinates | fit spread | events on the ridge |
|---|---|---|---|
| `--min-amp 100` | 7 274 | 4.99 / 5.03 mm | 86 % |
| `--min-amp 0` | **266 797** | **4.55 / 4.83 mm** | **93 %** |

An amplitude cut of 100 throws away 97 % of your four-fold coincidences *and
makes the correlation worse*. Use `--min-amp 0`.

One caveat, unexplained and worth someone's attention: the mean cluster size is
only **1.1–1.2 strips** in every zone, including the 1.0 mm pitch region. That
is low for a minimum-ionising particle — you would normally expect 2–4. The
likely causes are the analyser's per-strip threshold or the fixed
drift 600 / resist 420 operating point, but I have not tested either. It does
not block tracking; it does mean the position resolution is closer to
pitch/√12 (binary) than to a proper charge-weighted centroid.

---

## 8. Alignment — translation only

The two detectors are mounted square to the beam, so the transformation from
front to back is expected to be a pure shift. **You should verify that rather
than assume it**, and the check is free: fit a straight line of back position
against front position and look at the slope.

* slope = **1** → pure translation. Alignment is one number per axis, the offset.
* slope ≠ 1 → a scale error (usually a wrong pitch in some zone) or a rotation.

Measured on `drift_mesh_scan_1/nominal_00`, 266 797 matched events,
`--min-amp 0`, robust fit with 2.5σ clipping:

```
x:  back = +0.9996 * front  - 1.001 mm     spread 0.66 mm    81.0 % of events used
y:  back = +1.0418 * front  - 4.246 mm     spread 0.60 mm    82.7 % of events used
```

So the alignment constants are **dx = −1.00 mm, dy = −4.25 mm** — a couple of
millimetres of mechanical offset, exactly as expected, and no rotation.

> **The y slope of 1.042 is not a scale error — it is a divergent beam.** An
> earlier version of this document treated it as an open question about the
> strip map. It is not: the same anisotropic magnification is measured
> independently at the three P2 stations, where it grows *linearly with z*
> (stretch 1.007 → 1.013 → 1.023 at z = 320 → 630 → 940 mm) with no shear
> (off-diagonal < 0.0013), and extrapolating it to dz = 1370 mm reproduces the
> slope above to within 0.16 %. That is a beam diverging from a virtual source
> **~40 m upstream in y** and hundreds of m in x — focused in one plane and not
> the other, which is normal for an SPS extraction line. A wrong pitch in a strip
> map would give the same wrong scale at every z. Reproduce with
> `python explore6_divergence.py` (front-only projection, one chunk: 41.7 m in y,
> 224 m in x) or read it off the frozen record, which uses the full track model
> over all six sub-runs (38.8 m in y, 320 m in x — §13.7). The y numbers agree to
> 7 %; x is barely constrained either way, because x really is parallel.

**What the "spread" of 0.6 mm is, and what it is not.** It is *not* the detector
resolution, though it is now close enough to be a useful bound. With only two
planes a straight line through two points is exact by construction — there is no
residual to measure and no redundancy. The spread is the width of the
front-to-back correlation: the two planes' resolutions in quadrature, plus the
beam's genuine angular spread over the 1370 mm lever arm. Taking it as an upper
bound gives ≲0.45 mm per plane. **You still cannot measure uRWELL resolution or
uRWELL efficiency from the two uRWELLs alone** — that is what the three P2 planes
in the middle are for.

> Historical note: before the back's wiring was corrected (§6.2) these spreads
> read 4.55 and 4.83 mm and the x slope read 1.045. Every one of those numbers
> was the broken back map, not physics.

---

## 9. Tracks

With alignment in hand, a track is two points and a straight line:

```python
xb_aligned = xb - dx
slope_x    = (xb_aligned - xf) / dz          # dz = 1370 mm
x_at(z)    = xf + slope_x * z                # z measured from the front plane
```

Measured track angles on `nominal_00`, after alignment:

```
x:  median +0.04 mrad    sigma (from IQR) 0.57 mrad
y:  median +1.64 mrad    sigma (from IQR) 1.22 mrad
```

Well under 2 mrad, centred on zero — the beam really is parallel, which is what
"the muons come in pretty straight" means quantitatively. Use the IQR-based
sigma, not the plain rms (7.5–8.3 mrad): the distribution has tails from wrong
pairings and the rms is dominated by them.

These are *not* pure divergence. `align_and_track.py` removes only the offset,
not the scale, so the y number still contains the magnification of §8 (0.042 ×
a 24 mm beam σ over 1370 mm ≈ 0.7 mrad) on top of the two planes' resolutions
(0.6 mm / 1370 mm ≈ 0.44 mrad). Both axes are consistent with a genuine
divergence well under 1 mrad — which is why, in practice, the **front plane
alone predicts a P2 position about as well as a two-point track does** (§13).

> Before the back's wiring was corrected (§6.2) these read 4.33 and 4.59 mrad.
> That apparent divergence was almost entirely back-plane noise. The practical
> consequence is large: with a genuinely parallel beam the *front alone* points
> at a P2 plane about as well as a two-point track does, so a broken back plane
> does not stop you — see §13.

Extrapolated to the three P2 planes (front local frame, mm):

```
P2_IN  (z=320):  x  64.7 ± 20.2    y  56.2 ± 24.4
P2_MID (z=630):  x  64.7 ± 20.0    y  56.6 ± 24.3
P2_OUT (z=940):  x  64.7 ± 20.1    y  57.0 ± 24.4
```

The beam barely moves across the telescope, as it should.

**The next step, once you trust the above,** is to merge in a P2 station from
`combined_hits_root`, predict its position from the uRWELL track, and histogram
(measured − predicted). **That is now done — see §13.**

---

## 10. The code

Everything lives in `/local/home/banco/P2_basket_analysis/sps_beam_analysis/urw_reference/`.

| file | what it does |
|---|---|
| `urw_lib.py` | the library: geometry (`UrwGeometry`), streaming hit reader (`iter_hits`), clustering (`cluster_hits`), one point per event (`leading_points`) |
| `align_and_track.py` | **start here** — alignment + tracks + plots. §8 and §9 are its output |
| `urw_p2_efficiency.py` | **the deliverable** — P2 residuals and efficiency against uRWELL tracks. §13 is its output |
| `urw_qa.py` | per-sub-run QA: occupancy, spectra, cluster size, beam spot, timing |
| `explore4_back_map.py` | scores the four candidate wirings per view (§6.2) |
| `explore5_back_zones.py` | back residual vs position, split by pitch zone |
| `explore6_divergence.py` | the beam-divergence measurement (§8) |
| `check_view_transform.py` | the old anchor-free connector-order test — **superseded**, see its header |
| `check_connector_order.py` | the same superseded verdict via the P2 pad map |
| `ORDERING.md` | write-up of the connector-order finding; read the box at the top first |

Minimal example — build points for one detector yourself:

```python
import sys, os
sys.path.insert(0, '/local/home/banco/P2_basket_analysis/sps_beam_analysis/urw_reference')
import urw_lib as U
import pandas as pd

run_dir = '/local/home/banco/P2_data/TB_July2026_H4/runs/drift_mesh_scan_1'
sub_dir = os.path.join(run_dir, 'nominal_00')
run_json = os.path.join(run_dir, 'run_config.json')

# geometry: reads run_config.json + the strip map, applies the measured wiring
geo = U.UrwGeometry('EIC_uRWELL_front', run_json, sub_run_name='nominal_00')
print(geo)                 # active size, z position, which FEU connectors
print(geo.zone_table())    # the pitch zones of §4

# stream the hits (never splits an event across chunks) and cluster
points = []
for chunk in U.iter_hits(U.feu_hit_files(sub_dir, feu=geo.feu_num)):
    clusters = U.cluster_hits(chunk, geo, min_amp=0)
    points.append(U.leading_points(clusters))
points = pd.concat(points, ignore_index=True)

print(points.head())       # eventId, x, y, x_size, y_size, x_charge, ...
```

`zone_table()` is the quickest way to see §4 and §6 at once. For the back
detector it prints:

```
 zone view  pitch interpitch  n_strips feu_connectors  pos_min  pos_max
    0    y    1.0        0.1        64            [7]  64.2500 127.2500
    1    y    0.5        0.1        32            [8]   0.0000  15.5000
    2    y    1.5        0.1        32            [8]  16.5000  63.0000
    3    x    1.0       0.75        64            [5]  64.4375 127.4375
    4    x    0.5       0.37        32            [6]   0.0000  15.5000
    5    x    1.5       1.12        32            [6]  16.8750  63.3750
```

There are the three pitch zones per view from §4.4. `feu_connectors` and the
positions are reported *after* the wiring of §6.2 and the axis flip of §6.3, so
this is what the analysis actually uses — the back's connectors are **not**
interchanged (connector 7 carries y's high positions, 8 the low ones, in map
order), and the 1.0 mm zone is the upper half of each view.

Useful arrays on `geo`, all indexed by global channel 0–511:
`geo.view` (`'x'`/`'y'`/`''`), `geo.pos` (local mm), `geo.pitch`, `geo.zone`,
`geo.mapped` (bool).

To do the same for the back detector, change the name to
`'EIC_uRWELL_back'` — it is the same FEU, so `feu_hit_files` returns the same
files, and the geometry object handles the different channel range and map.

---

## 11. Checklist of traps

In rough order of how much time each will cost you:

1. `unset PYTHONPATH`, and use the venv interpreter. (§2)
2. `axis` is the direction the strip **runs**, so `axis='y'` measures **x**. (§4.1)
3. All four views are wired differently from the map file, and the back needs a
   **within-connector reversal**, not a connector swap. Handled by `urw_lib`; do
   not undo it. (§6.2)
4. The back detector has **three pitch zones per view**. Never cluster by channel index. (§4.4)
5. **Do not apply an amplitude cut** — you lose 97 % of coincidences for nothing. (§7)
6. `significance` is already cut at 5 by the analyser; cutting at 5 again is a no-op. (§5.3)
7. ROOT files have multiple **cycles**; take the highest, never the sum. (§5.3)
8. The combined file is **ordered by FEU**; truncating the read gives you only FEU 1. (§5.4)
9. Read **all chunk files** for a sub-run, not just chunk 000. (§5.2)
9b. Most `hits_root/` directories have been **deleted**; the fallback to the
    combined files needs `feu=` on `iter_hits` or P2 hits become uRWELL
    strips, silently. (§5.3)
10. Check `recorded_events.npz` → `feuN_missing` before trusting a sub-run. (§5.5)
11. Two planes give **no residual and no efficiency** for the uRWELLs themselves. (§8)

---

## 12. Status — what is done, what is not

**Verified and reproducible:**

* Geometry, strip maps and the channel → strip wiring of all four views, with
  the back's within-connector reversal settled against the front (§6.2).
* Clustering, occupancy, amplitude spectra, cluster size, timing.
* Front-to-back alignment: dx = −1.00 mm, dy = −4.25 mm, no rotation; spread
  0.6–0.7 mm, so ≲0.45 mm per plane.
* Two-point tracks: the beam is parallel to 0.6–1.2 mrad.
* The uRWELL → P2 frame relation: a proper rotation of −60°, det +1 (§13).
* P2 residuals and uRWELL-referenced P2 efficiency over all six sub-runs of
  `highstat_eff_1` (§13).

**Not started:**

* Running over `drift_mesh_scan_1`'s 23 sub-runs — a mesh/drift scan, so the
  interesting product there is efficiency *versus HV*, which
  `urw_p2_efficiency.py` will produce as-is (it already loops sub-runs and
  reads each one's HV from `run_config.json`).
* uRWELL efficiency and resolution proper, which needs the P2 stations used the
  other way round (P2 tags, uRWELL probes). The pieces all exist.
* A resolution estimate per pitch zone of the back detector.

**Open questions:**

* Mean cluster size of 1.1–1.2 strips is low for a MIP (§7). Cause unknown —
  analyser threshold or gain. This is the likely reason the uRWELL resolution
  sits near 0.45 mm rather than pitch/√12 = 0.29 mm.
* A small ~5 × 5 mm patch of the front uRWELL, at local (63, 70) mm, gives a
  wrong position: it shows up as a hole in the efficiency map of **all three**
  P2 stations at the same place (§13). Cause not yet identified.
* The absolute lab direction of x and y is still not determined by this data.
  The mirror ambiguity of §6.3 is fixed *relative to the P2 pad frame*; tying
  that to the lab needs the P2 frame's own orientation or a survey.

---

## 13. Referencing the P2 stations — residuals and efficiency

This is what the telescope is for. `urw_p2_efficiency.py` does the whole thing:

```bash
unset PYTHONPATH
cd ~/P2_basket_analysis/sps_beam_analysis/urw_reference
/local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python urw_p2_efficiency.py \
    --run highstat_eff_1 --out out_eff
```

It builds uRWELL tracks from the FEU-1 hits files, streams each P2 station's
`combined_hits_root` into one leading cluster per event, matches on `eventId`,
fits the uRWELL frame onto the P2 pad frame, and writes per-station plots, a
JSON and a CSV.

### 13.1 The two frames

The P2 pad map lives in the BASKET PCB's own frame (a fan of ~12 × 11.6 mm pads
described by radius and φ), which is rotated with respect to the uRWELL. Fitting
a **free 2×2 matrix** rather than assuming a rotation is deliberate: the two
detectors are seen from the same side along the same beam, so the answer *has*
to be a proper rotation, and getting one is a genuine test that both strip maps
and the pad map are right. Measured on `highstat_eff_1/beam_commissioning_00`:

| station | rotation | det | singular values | shear |
|---|---|---|---|---|
| P2_IN | −59.68° | +1.005 | 1.0072 / 0.9974 | −0.0012 |
| P2_MID | −59.74° | +1.012 | 1.0133 / 0.9991 | −0.0011 |
| P2_OUT | −60.01° | +1.023 | 1.0230 / 0.9995 | −0.0008 |

A rotation of −60° with no reflection and no shear. The 0.5–2 % departure from
orthogonal is the beam divergence of §8, and it is the reason the applied
transform is the affine rather than the rigid one: it would otherwise leak into
the residuals as a ±1 mm trend across the plane.

### 13.2 Residuals

```
P2_IN    dx rms 3.39 mm    dy rms 3.47 mm
P2_MID   dx rms 3.37 mm    dy rms 3.44 mm
P2_OUT   dx rms 3.39 mm    dy rms 3.45 mm
```

All three agree, as they must — same PCB. The pads are 12.05 × 11.60 mm, and
12 mm / √12 = **3.46 mm**, so the P2 position resolution here is exactly what a
uniformly illuminated pad gives with **no charge sharing at all** — consistent
with 92–95 % of clusters being a single pad and a mean cluster size of 1.05–1.08.
The 2D residual plot shows the pad footprint directly, as a filled square tilted
at the −60° of §13.1. The uRWELL contributes ~1 mm to this in quadrature, i.e.
essentially nothing — which is the whole point of using it as the reference.

> This is a statement about the *readout as configured*, not a fundamental limit:
> a resistive pad detector that shared charge over two or three pads would do
> considerably better than 3.4 mm. Whether that is threshold, gain or the
> analyser's `significance > 5` cut is worth following up.

### 13.3 Efficiency

```
efficiency = N(track in the pads AND a P2 cluster within --probe-r)
             ------------------------------------------------------
                        N(track in the pads)
```

with 68.27 % Clopper-Pearson intervals. Dropped from the denominator, and each
reported in the JSON: triggers the probe's FEU never recorded
(`recorded_events.npz`), events inside that station's HV spark window or before
its mesh settled, and tracks whose projection falls outside the pad plane. Dead
connectors are deliberately **not** removed — they stay in the denominator and
appear as holes in the map.

On `beam_commissioning_00`, ~1.15 M tracks per station:

| station | efficiency | misses with no P2 hit | misses out of range |
|---|---|---|---|
| P2_IN | 0.9649 ± 0.0002 | 33 664 | 6 745 |
| P2_MID | 0.9706 ± 0.0002 | 27 945 | 6 040 |
| P2_OUT | 0.9604 ± 0.0002 | 37 295 | 7 037 |

The errors are statistical only. The systematic that matters is the choice of
`--probe-r`, and the stage measures it for you by scanning it (in the JSON as
`efficiency.vs_probe_r`, and printed):

```
r [mm]   5       8       10      12      15      18      20      25      30      40
eff      0.5854  0.9624  0.9655  0.9663  0.9670  0.9673  0.9674  0.9677  0.9681  0.9686
```

It plateaus by ~10 mm — one pad — and then creeps up by only **0.0006 per
10 mm**, which is the accidental-match rate. Tightening to 10 mm costs 0.15 %,
loosening to 40 mm gains 0.16 %.

The other knob is `--track-cut`, which decides how good a uRWELL track has to be.
Varying it over a factor of ~7 (73 % → 91 % of four-coordinate events kept):

| `--track-cut` | tracks kept | P2_IN | P2_MID | P2_OUT |
|---|---|---|---|---|
| 1.5 mm | 73.4 % | 0.9682 | 0.9751 | 0.9685 |
| 3.0 mm (default) | 83.3 % | 0.9670 | 0.9749 | 0.9685 |
| 10 mm | 91.0 % | 0.9642 | 0.9726 | 0.9660 |

≤0.4 % across the whole range, and drifting the right way — a looser cut admits
tracks that point worse, so more of them miss. So read the efficiencies as
**96–97 % with a ~0.3 % systematic**; the statistical error is negligible next
to it. Note the working point: all six sub-runs sit at P2 drift 700 V / mesh
450 V, which is the *top* of the mesh range scanned in `drift_mesh_scan_1`
(445 → 390 V), so this is the highest-gain point taken, and there is no evidence
here that the efficiency plateau was reached.

### 13.4 What the six sub-runs show

All six sit at the same working point, so `summary_highstat_eff_1.png` is a
stability measurement rather than a scan, and two things come out of it.

The **frame is mechanically stable**: the fitted uRWELL → P2 rotation moves by
less than 0.02° per station across the ~4.5 hours the run spans, and the residual
widths sit in 3.35–3.47 mm throughout. Nothing moved.

The **efficiency drifts down and then recovers**:

| sub-run | rate | P2_IN | P2_MID | P2_OUT |
|---|---|---|---|---|
| 00 | 4629 Hz | 0.9649 | 0.9706 | 0.9604 |
| 01 | 4642 Hz | 0.9633 | 0.9643 | 0.9513 |
| 02 | 4621 Hz | 0.9626 | 0.9608 | 0.9493 |
| 03 | 4641 Hz | 0.9623 | 0.9589 | 0.9470 |
| 04 | 4619 Hz | 0.9621 | 0.9575 | 0.9463 |
| 05 | 3232 Hz | 0.9660 | 0.9672 | 0.9567 |

Sub-runs 00–04 are consecutive 30-minute blocks at a constant ~4.6 kHz, and
P2_MID and P2_OUT lose about **1.3 points** monotonically over the two and a half
hours. Sub-run 05 comes after a ~2.5 hour gap and at 30 % lower rate, and both
recover almost fully. P2_IN is nearly flat throughout (0.3 points).

That shape — monotonic loss under continuous irradiation, recovery after a
rest — is what charging-up of a resistive layer looks like. It is **not proven
here**, because sub-run 05 changes two things at once (the pause *and* the rate),
so the two explanations are not separated. Doing so needs either a rate scan at
fixed time or a long run split finely in time; the second is nearly free, since
every track already carries `t_ns` and the efficiency could be binned within a
sub-run. Worth doing before quoting an efficiency for these detectors.

**Compare with the P2's own tag-and-probe** (stage 22 of
`P2_basket_analysis/sps_beam_analysis`, same sub-run), which has to use the other
two P2 stations as the reference: 92.5 % / 96.3 % / 92.5 %. It is 4 points lower
on the outer stations and much less uniform, because its tag position is smeared
by the tagging planes' own 12 mm pads even with a 27 mm probe radius. That gap —
not the absolute number — is the clearest demonstration of what the uRWELL
reference buys.

### 13.5 A free cross-check of the P2 mapping

The P2 team found on 2026-07-25 that P2_IN's `c_5_top` ribbon is mounted flipped
relative to P2_MID and P2_OUT, and encoded it as
`sps_config.STRATEGY_OVERRIDES = {'P2_IN': {(5, 'top'): 'linear'}}`. That was
derived from the P2 stations comparing themselves with each other. The uRWELL
reference confirms it independently, and the sensitivity is brutal — same chunk,
same tracks, only the override toggled:

| | residual rms | efficiency |
|---|---|---|
| with the override | 3.39 / 3.47 mm | 0.9670 |
| without it | 5.24 / 6.69 mm | 0.5959 |

It also shows what a mapping error looks like in these plots, which is worth
knowing: the residual widens but does **not** blow up (the correctly-mapped 7/8
of the channels still dominate it), the frame fit still returns a plausible
−60.5° rotation, and what really gives it away is the **efficiency** collapsing
and the accidental slope going from 0.0007 to 0.0272 per 10 mm. If you ever see
an efficiency far below 90 % here, suspect the mapping before the detector.

### 13.6 Reading the plots

Each station gets an eight-panel PNG. The two efficiency maps are the point:
one binned in the **P2 pad frame** and one in the **uRWELL frame**. A defect of
the probe sits still in the P2 frame; a defect of the reference sits still in the
uRWELL frame and appears in *all three* stations at the same place. That is how
the ~5 × 5 mm hole at uRWELL local (63, 70) mm was identified as a front-uRWELL
problem rather than three coincidentally dead pads.

The track-density panel (uRWELL frame) shows the reference's own dead strips as
white lines. They cancel in the efficiency ratio — the efficiency map is flat
across them — which is a useful sanity check that the ratio is being formed
correctly.

### 13.7 Where the mapping and alignment are frozen

Everything the measurement rests on that is *not* in the raw data — the
channel → strip wiring, the front→back alignment, the uRWELL → P2 transforms,
the beam optics — is written out by `record_mapping_alignment.py` into

```
/local/home/banco/P2_data/TB_July2026_H4/analysis/urw_referenced_efficiency/
    mapping_urwell.csv         512 rows: detector, channel -> view, position, pitch, zone
    mapping_alignment.json     the same numbers, machine readable
    MAPPING_AND_ALIGNMENT.md   the same numbers, with the evidence for each
```

`mapping_urwell.csv` is the one to reach for if you are not running this code:
join on `(detector, channel)` and you have a position in mm, with the wiring of
§6.2 and the axis convention of §6.3 already applied. The `view` column is the
coordinate the strip **measures**, so the `axis` trap of §4.1 is already undone.

The code remains authoritative — regenerate the record after changing
`urw_lib.VIEW_MODE_DEFAULT`, do not hand-edit it.

### 13.8 Efficiency versus HV

`drift_mesh_scan_1` scans both electrodes, and `plot_hv_curves.py` turns the
stage's CSV into curves in the style of stage 22:

```bash
$PY urw_p2_efficiency.py --run drift_mesh_scan_1 --out <analysis>/urw_referenced_efficiency/drift_mesh_scan_1
$PY plot_hv_curves.py --csv <that dir>/urw_p2_efficiency_drift_mesh_scan_1.csv --out <that dir>
```

**Grouping the points is the part to get right.** A mesh scan here does *not*
hold the drift fixed: the drift tracks the mesh so the drift *field* stays
constant (P2_MID and P2_OUT at drift = mesh + 250 V, P2_IN at mesh + 200 V),
while the separate `drift_*` sub-runs hold the mesh at 450 V and walk the drift
from 450 to 900 V. So the invariant that defines a scan family is
**drift − mesh**, not drift. Points are grouped by that, families are matched
across stations by rank, and a genuine 2D scan therefore comes out as one figure
per drift setting. Repeated HV settings — P2_IN sits at mesh 430 V through the
whole drift scan, and `nominal_00` repeats `drift_700` exactly — are combined by
adding their counts, not by overplotting.

Two traps that this run walked straight into:

* `drift_450` sets P2_MID and P2_OUT to drift = mesh = 450 V, i.e. **zero drift
  field**. Efficiency there is ~15 %, and that is real, not a bug.
* `sps_cluster.settle_t_min` derives the HV-settle cut by matching
  `hv_monitor.csv` timestamps against a DAQ start parsed out of the chunk
  filename. On this run that lands ~9.6 h off and would silently discard every
  event. The stage now refuses a settle cut longer than the sub-run's own span
  of trigger timestamps and says so. Worth checking whether stages 21/22 are
  affected the same way.

**Results on `drift_mesh_scan_1`** (all 23 sub-runs, in
`analysis/urw_referenced_efficiency/drift_mesh_scan_1/`):

* Mesh curve, at constant drift field: P2_MID 0.806 → 0.961 and P2_OUT 0.857 →
  0.950 over 390 → 450 V; P2_IN 0.314 → 0.917 over 370 → 430 V, which captures
  the whole turn-on. **None of the three has reached a plateau** at the top of
  the scan, which is also the `highstat_eff_1` working point.
* Drift curve, mesh held at 450 V: 0.157 at drift = mesh (zero field), 0.896 by
  500 V, then a slow rise flattening above ~850 V at 0.974 for both stations —
  +1.1 (P2_MID) and +2.1 (P2_OUT) points over the 700 V nominal.
* P2_IN, held at fixed HV through the 100-minute drift scan, stayed within
  0.9146–0.9179. A 0.3-point spread is a useful bound on how much of any curve
  could be time drift rather than voltage.
* The frame fit is unchanged across all 23 points: rotation within 0.02° and
  residual widths 3.35–3.56 mm, the latter improving slightly with mesh HV as
  clusters grow.
