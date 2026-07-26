# uRWELL Dream-connector ordering — finding, 2026-07-25

> ## SUPERSEDED IN PART, 2026-07-26 — read this box first
>
> The verdict table below is **right for the front and wrong for the back**, and
> the "orientation vs pair order" argument near the end is **wrong**.
>
> Every test in this document is blind to a reversal of channel order *inside* a
> connector, because for a **uniform** map that reversal is degenerate with a
> mirror of the whole view. The front's map is uniform, so the front is fine.
> The **back's is not** — `strip_map.txt` gives it three pitch zones per view
> (1.0 / 1.5 / 0.5 mm) — and there the degeneracy breaks and the tests below
> picked the wrong member of the pair.
>
> With the front established as an external reference good to <1 mm (it points
> at all three P2 stations, `explore6_divergence.py`), `explore4_back_map.py`
> scores all four candidate wirings per view by the width of `back − front`:
>
> | view | `AB` | `BA` (this doc's answer) | `AB_rev` | `BA_rev` |
> |---|---|---|---|---|
> | back x | 45.6 mm | 5.42 mm | **0.80 mm** | 47.0 mm |
> | back y | 35.2 mm | 6.40 mm | **0.88 mm** | 38.8 mm |
>
> So the back's connectors are **not** interchanged; the channel order inside
> each is **reversed**. Using `BA` left the back pointing at 4.4 mm — useless
> for tracking, and it is why the back looked like it contributed nothing to a
> two-point track. Fixed, it points at ~0.7–1.2 mm, like the front.
>
> This also means the run config's `dream_feu_orientation: inverted` on all eight
> connectors (set 2026-07-25) **is** consistent and physical: every connector
> carries a within-connector reversal, and on top of that exactly one view
> (front y) has its connector pair interchanged. The claim below that a uniform
> orientation could not explain the pattern only considered the pair order.
>
> Live encoding: `urw_lib.VIEW_MODE_DEFAULT` (+ `AXIS_FLIP_DEFAULT`, a pure
> labelling choice so both planes measure in the same direction).
> `urw_qa.py --raw-view-mode` reproduces the untouched map order.

Follow-on to `~/P2_basket_analysis/sps_beam_analysis/urw_reference/URW_TRACKING_HANDOFF_2026-07-25.md`. The strip maps are
deployed and bound correctly (handoff §2/§5 all reproduce), but the
**channel → strip assignment is still wrong for three of the four views**: the
two 64-channel Dream connectors of a view are interchanged relative to the map's
connector column. Everything below was measured on this machine, on
`drift_mesh_scan_1/nominal_00`.

## Verdict

| view | FEU 1 conn | map order as-is | correction |
|---|---|---|---|
| front x | 1, 2 | wrong | **swap** |
| front y | 3, 4 | correct | keep |
| back x | 5, 6 | wrong | **swap** |
| back y | 7, 8 | wrong | **swap** |

*(2026-07-26: this table's back rows are wrong — see the box at the top. The
live encoding is `VIEW_MODE_DEFAULT`; `CONNECTOR_SWAP_DEFAULT` no longer exists,
and the flag is now `urw_qa.py --raw-view-mode`.)*

## Why the handoff did not catch it

Handoff §4.2 verified the channel *arithmetic*
(`connector = channel//64 + 1`, offset by `starting_connector`) against the data,
and that is right. What it could not verify is which physical cable sits on which
FEU connector, because the check used raw occupancy — and §7 already noted that
raw uRWELL occupancy "is near-uniform across the plane … raw hits are not a beam
spot". A near-uniform quantity is exactly the one that cannot detect a
permutation of two halves of the plane. Once hits are cut at `amplitude > 100`
the occupancy is strongly structured, and the permutation is obvious.

The run config records a related quantity, `dream_feu_orientation`, but **no
analysis code reads that field** — neither `Detector_Classes` nor
`sps_beam_analysis`/`p2_mapping`. It is metadata only. See
"Orientation vs pair order" below for why it is not the same degree of freedom.

## Evidence

### 1. Anchor against the P2 pad map — `check_connector_order.py`

`combined_hits_root/*_feu-combined_hits.root` is event-synchronised across all
four FEUs. For every trigger take the leading uRWELL channel in a view and the
leading P2 pad, whose position comes from the independently validated
`P2_BASKET_mapping.csv`. The median P2 coordinate vs uRWELL channel index must be
a single continuous monotone relation — one beam, same tracks. With the map
order as-is it breaks into two disjoint branches with a jump at the connector
boundary; `(idx + 64) % 128` repairs it.

Spearman rho, `amplitude > 150`:

| view | anchor axis | natural | swapped | verdict |
|---|---|---|---|---|
| front_x | P2_MID pad_cy | +0.323 | **−0.665** | swap |
| front_y | P2_MID pad_cx | **+0.812** | −0.344 | keep |
| back_x | P2_MID pad_cy | +0.335 | **−0.682** | swap |
| back_y | P2_MID pad_cx | −0.357 | **+0.761** | swap |

`--anchor P2_OUT` reproduces this exactly (+0.812/−0.334, −0.661, −0.672,
+0.764). `--anchor P2_IN` also says swap for x and back_y, but its front_y
numbers are +0.152 vs −0.205 — both consistent with noise, so it does not
discriminate for that view. P2_IN has ~3× fewer events, is the repaired station,
and is cabled on connectors 4,5,6,8; it is the weakest of the three anchors and
its front_y "verdict" should be ignored.

Plot: `diag_order_test.png` (top row natural, bottom row swapped).

### 2. Front vs back in position space

Independent of the P2 anchor. The two uRWELLs are 1370 mm apart on the same
beam line and are mapped by *different* files (`inter_map.txt` /
`strip_map.txt`), so their agreement is a real cross-check.

* map order as-is: both x and y break into two disjoint diagonal islands
  (`r = +0.78 / −0.65`, the y slope coming out negative).
* with the swaps: one continuous band over the whole plane. Robust
  (2.5σ-clipped) fit on 10 904 matched events, `amplitude > 100`:

  | | slope | offset | residual rms | on the ridge |
  |---|---|---|---|---|
  | x | +1.036 | +2.89 mm | 5.05 mm | 86 % |
  | y | +0.996 | +6.81 mm | 5.09 mm | 86 % |

  Slope 1 to within a few percent between two independently mapped detectors
  also validates the back map's 1.0 / 1.5 / 0.5 mm pitch zones — a wrong pitch
  in any zone would show up as a kinked slope. The few-mm offsets are the real
  relative misalignment; the 5 mm residual is beam divergence over 1370 mm plus
  scattering plus cluster resolution, and is the starting point for a proper
  alignment.

Plot: `diag_front_back_position.png`.

### 3. What it looked like before

`out/` (as-is) vs `out_fixed/` (corrected), same cuts:

* as-is: front x profile is a symmetric V with a hard zero at x ≈ 60 mm; the
  back populates all four corners of the plane with the middle empty.
* corrected: single contiguous blob in both detectors. Back core fit
  x = 61.3 ± 19.8 mm, y = 50.9 ± 15.9 mm; front x = 57.2 ± 21.9 mm.

### 4. Anchor-free per-view test — `check_view_transform.py`

Added 2026-07-25 after the cabling was re-checked at the beam. This needs no P2
and no front↔back matching, so it tests each of the four views on its own.

Both map files are uniform: the global index (map connector)·64 + channel runs
monotonically with the Gerber coordinate over the full 0–127 of every view, in
`inter_map.txt` and `strip_map.txt` alike. Interchanging a view's two connectors
is therefore exactly a cyclic shift of the position array by half the plane. The
beam straddles the middle of both detectors, so the wrong choice **splits the
spot in two and throws the pieces to opposite edges**.

Leading-cluster profile, `amplitude > 100`, `nominal_00`:

| view | σ as-is | one-blob fraction as-is | σ interchanged | one-blob fraction | verdict |
|---|---|---|---|---|---|
| front x | 46.2 mm | 0.39 | **21.9 mm** | **0.98** | swap |
| front y | **28.4 mm** | **0.998** | 38.5 mm | 0.62 | keep |
| back x | 52.8 mm | 0.48 | **15.2 mm** | **0.97** | swap |
| back y | 50.0 mm | 0.38 | **15.3 mm** | **0.95** | swap |

The raw profiles make it plain — as-is, front x / back x / back y each have a
hard zero straddling the plane centre with the entries piled at both edges,
which is the wraparound signature; front y has no such gap and is already one
contiguous blob. Three independent methods now agree on all four views.

## Orientation vs pair order

`dream_feu_orientation` was corrected on 2026-07-25 to `inverted` on all eight
uRWELL Dream connectors (previously `x1/x2 normal` with `y1/y2 inverted` on the
front and `rotated` on the back — a bookkeeping error). Generator
`run_config_beam.py::_urwell`, propagated to every `run_config.json` carrying
these detectors.

**That change does not affect the table above, and must not be read as fixing
it.** The two are different degrees of freedom:

* *orientation* is how a flat cable is plugged — a permutation **within** a
  64-channel connector;
* *pair order* is which of a view's two detector connectors lands on the lower
  FEU Dream connector — a permutation **between** the two 64-channel blocks.
  That is what `dream_feus` (`x1`/`x2`, `y1`/`y2`) claims and what the data
  above contradicts for three views.

A uniform orientation cannot produce a non-uniform pair order. Writing it out:
let `I` be the map as-written, `S` the pair interchange, `W` a within-connector
reversal, so `W∘S = R`, the full 128-channel reversal. The tests above
distinguish `{I, R}` from `{S, SR}` (a global reflection of an axis leaves both
the rank correlation and the profile compactness unchanged). Measured, front y
is in `{I, R}` and the other three are in `{S, SR}`. A uniform cabling rule
applied to a uniform map gives one class for all four views: `W` everywhere
lands all four in `{S, SR}`, and `W∘S` everywhere lands all four in `{I, R}`.
Neither matches. So the front y pair is genuinely cabled in the opposite order
from the other three, whatever the plug orientation is.

The practical rule: `CONNECTOR_SWAP_DEFAULT` in `urw_lib.py` stays as it is, and
anyone wiring `dream_feu_orientation` into the mapping code must apply it as a
within-connector permutation only — never as a uniform view-level correction, or
front y breaks.

**2026-07-26: the paragraph above is wrong**, and the box at the top of this file
says why. The step it skips is that `{I, R}` and `{S, SR}` are classes of *two*
elements each, so "front y is in `{I, R}`, the others are in `{S, SR}`" does not
pin any of the four. A uniform within-connector reversal `W` on all eight
connectors, plus a pair interchange on front y alone, reproduces the measured
classes exactly — `W ∈ {S, SR}` for the three, `W∘S ∈ {I, R}` for front y. That
is the wiring the data actually has, and it is what the run config says.
The live table is `urw_lib.VIEW_MODE_DEFAULT`, not `CONNECTOR_SWAP_DEFAULT`
(removed).

## Consequences

* Any uRWELL geometry result produced before this — positions, profiles,
  residuals, efficiency, alignment — is folded and must be redone. Nothing on
  disk is affected yet: the only uRWELL product in `analysis/` is
  `23_beam_profile` output for the front, and handoff §6 stages 21/22/26/28 never
  enrolled the uRWELLs at all.
* Handoff §6 step 2 (teach `sps_config.channel_table()` a strip branch) must
  apply the swap, or it will inherit the fold.
* Worth checking whether the P2 stations have the same class of problem. Their
  pad map is validated in position space, so probably not, but the same
  connector-order degree of freedom exists there and `dream_feu_orientation` is
  equally unread for them.

## Open

* The swap is degenerate with "reverse the channel order within each connector"
  — the two differ by an overall reflection of the axis. The P2 anchor fixes
  the sign *relative to the P2 basket frame*; establishing the absolute lab
  direction needs the P2 frame's own orientation, which the anchor test does not
  supply.
* `active_size` still comes from the map extent (handoff §4.5) and is unchanged
  by the swap.
