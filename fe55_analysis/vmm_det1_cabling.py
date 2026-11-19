#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vmm_det1_cabling.py — (vmm, ch) -> pad for ONE P2 basket detector read out by
nine VMM hybrids, one per connector: the Fe55 bench.

The cabling itself is a p2-mapping setup, `setups/fe55_det1_2026.yaml`, and the
whole chain (vmm, ch) -> connector, half -> strip -> channel_id -> pad geometry
is p2-mapping's (`p2_mapping.table.build_table`), so a pad number means the same
thing here as at SPS and on the cosmic stand. Live time comes from p2-decode's
pcapng walker. What is left in this module is the bench-specific part:

  * the hybrids that are powered but NOT plugged into det1 (connectors 1 and 9),
    kept in the setup file under `disconnected:` and used to check the VMM
    numbering from data (`verify_cabling`);
  * variants of the cabling for the mapping scan and the self-test (another
    channel order, another disconnected list), built as modified p2-mapping
    setups rather than with arithmetic of their own;
  * the per-pad helpers the Fe55 tools share (pad keys, polygons, live time).

Chain
-----
    (vmm, ch)
      -> hybrid, connector_N (1..9), half           setups/fe55_det1_2026.yaml
      -> strip  = half_base + f(within)             p2_mapping.routing (order)
      -> sector = connector_N - 1
      -> channel_id = sector*128 + (strip - 1)
      -> pad_cx, pad_cy, pad_area, ...              p2-mapping p2_basket_pads.csv

CABLING AS OF 2026-09-09 (Alexandra, Fe55 bench, det1, drift 750 V / mesh 420 V)
  hybrid 1 -> connector 1 -> VMM 2 (bot) / 3 (top)
  hybrid 2 -> connector 2 -> VMM 4 (bot) / 5 (top)
  ...
  hybrid 9 -> connector 9 -> VMM 18 (bot) / 19 (top)
  connectors 1 and 9 DISCONNECTED  -> sectors 0 and 8 read out by nobody
  connector 10 (sector 9) is not cabled in this setup at all

So the live pads are sectors 1-7 (connectors 2-8), VMMs 4-17: 7*128 = 896 of
the detector's 1280 pads.

Moved from P2_basket_online_analysis (tag legacy-final), where it carried its
own copy of the cabling arithmetic and of P2_BASKET_mapping.csv.

@author: ak271430 Alexandra Kallitsopoulou
"""

import dataclasses
import os

import numpy as np
import pandas as pd

from p2_decode.pcapng import walk_blocks
from p2_mapping import routing
from p2_mapping.setups import load_setup
from p2_mapping.table import build_table, load_pads

__all__ = ["build_pad_table", "full_pad_map", "attach_pads", "pad_polygons",
           "vmm_half_map", "expected_vmms", "verify_cabling", "live_time",
           "capture_window", "SETUP_FILE", "CABLING_WHEN"]

HERE = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------------------------------------
# The setup
# --------------------------------------------------------------------------
SETUP_FILE = os.path.join(HERE, "setups", "fe55_det1_2026.yaml")

# Which cabling entry of the setup the tools use. Fixed rather than taken from
# each capture's time: the self-test writes synthetic captures with arbitrary
# timestamps, and a re-cabling of the bench is a new dated entry in the setup
# file plus a new value here (or `when=` on the calls), not a silent change.
CABLING_WHEN = "2026-09-09T00:00:00"

SETUP = load_setup(SETUP_FILE)
DETECTOR = "det1"
_DESIGN = routing.design(SETUP.design)

STRIPS_PER_CONN = _DESIGN["strips_per_connector"]   # 128 = 2 halves
CH_PER_HALF = _DESIGN["channels_per_half"]          # 64 = one VMM
N_CONNECTORS = _DESIGN["n_connectors"]              # 10
N_PADS = N_CONNECTORS * STRIPS_PER_CONN             # 1280

# Kept for callers that ask where the geometry came from; it is p2-mapping's.
MAP_CSV = str(routing.pads_path(SETUP.design))


def _hybrid_rows(when=CABLING_WHEN):
    """{hybrid: row} for every hybrid on the FEC, plugged in or not, and the
    set of hybrids that are not plugged in. Rows are p2-mapping cabling rows."""
    live = {int(r["hybrid"]): dict(r) for r in SETUP.cabling_at(when).hybrids}
    dead = {}
    for r in SETUP.extra.get("disconnected", ()):
        r = {"detector": DETECTOR, **r}
        dead[int(r["hybrid"])] = r
    clash = set(live) & set(dead)
    if clash:
        raise ValueError(f"{SETUP_FILE}: hybrids {sorted(clash)} are both cabled "
                         "and listed under disconnected")
    return {**live, **dead}, set(dead)


_ROWS, _DEAD = _hybrid_rows()

# Hybrids present, in connector order. Hybrid n sits on connector n.
HYBRIDS = tuple(sorted(_ROWS))
# Connectors with no hybrid plugged in. Their channels are omitted from the
# table entirely, so their hits (there should be none) never reach a pad and
# their pads never enter a gain/rate map as "dead pads".
DISCONNECTED = tuple(sorted(int(_ROWS[h]["connector"]) for h in _DEAD))

# Global VMM id of the BOTTOM chip of hybrid 1: hybrid n carries VMM
# base + 2(n-1) / base + 2(n-1) + 1. MEASURED 2026-09-09 (see the setup file):
# the initial guess base = 0 described the LOCAL chip numbering on the hybrid,
# shifted every connector by one sector, and was caught by the two disconnected
# connectors -- base = 0 predicts VMM 0,1 and 16,17 silent, base = 2 predicts
# VMM 2,3 and 18,19 silent, and the data shows the latter.
VMM_BASE = int(_ROWS[HYBRIDS[0]]["bottom_vmm"]) - 2 * (HYBRIDS[0] - 1)
for _h, _r in _ROWS.items():
    if (_r["bottom_vmm"], _r["top_vmm"]) != (VMM_BASE + 2 * (_h - 1),
                                             VMM_BASE + 2 * (_h - 1) + 1):
        raise ValueError(f"{SETUP_FILE}: hybrid {_h} breaks the 2n/2n+1 VMM "
                         "convention that VMM_BASE and the self-test assume")

# Which detector half the hybrid's BOTTOM VMM reads ('bot' = strips 1-64), and
# the within-half order every half uses unless its row says otherwise: the
# p2_basket design's defaults. 'reverse' is validated on the SPS beam and for
# the DREAM/K59V readout; it is NOT independently validated for this bench
# cabling -- run `vmm_fe55_map.py <file> --scan-mapping` to confirm.
BOTTOM_VMM_HALF = _DESIGN["bottom_vmm_half"]
WITHIN_STRATEGY = _DESIGN["default_order"]

# Per-(connector_N, half) exceptions: the rows' own bottom_order / top_order.
STRATEGY_OVERRIDES = {}
for _r in _ROWS.values():
    _top = "top" if BOTTOM_VMM_HALF == "bot" else "bot"
    for _k, _half in (("bottom_order", BOTTOM_VMM_HALF), ("top_order", _top)):
        if _r.get(_k):
            STRATEGY_OVERRIDES[(int(_r["connector"]), _half)] = _r[_k]

# Fe55 is a 5.9 keV photon: the conversion is a point on the scale of a 12 mm
# pad, so charge is contained in one pad and the single-pad ADC spectrum IS the
# Fe55 spectrum. Recorded here because the number drives the gain estimator.
PAD_PITCH_MM = 12.0


def hybrid_vmms(hybrid, base=VMM_BASE):
    """(bottom_vmm, top_vmm) global ids for a hybrid."""
    return base + 2 * (hybrid - 1), base + 2 * (hybrid - 1) + 1


def _setup_variant(disconnected=DISCONNECTED, base=VMM_BASE, hybrids=HYBRIDS,
                   when=CABLING_WHEN):
    """The p2-mapping setup with `hybrids` minus `disconnected` plugged in.

    The real cabling is the default; other lists are how the self-test and the
    command line describe a bench that is cabled differently. `base` shifts the
    VMM ids of every row, for testing the numbering convention."""
    drop = {int(c) for c in (disconnected or ())}
    pool = _ROWS if when == CABLING_WHEN else _hybrid_rows(when)[0]
    rows = []
    for h in hybrids:
        r = dict(pool[h]) if h in pool else {"detector": DETECTOR,
                                               "connector": h, "hybrid": h}
        if int(r["connector"]) in drop:
            continue
        r["bottom_vmm"], r["top_vmm"] = hybrid_vmms(h, base)
        rows.append(r)
    entry = dataclasses.replace(SETUP.cabling_at(when), hybrids=tuple(rows))
    return dataclasses.replace(SETUP, cabling=(entry,))


def vmm_half_map(disconnected=DISCONNECTED, bottom_half=BOTTOM_VMM_HALF,
                 base=VMM_BASE, hybrids=HYBRIDS):
    """{vmm_id: (connector_N, half)} for every instrumented VMM."""
    top_half = "top" if bottom_half == "bot" else "bot"
    out = {}
    for r in _setup_variant(disconnected, base, hybrids).cabling[0].hybrids:
        n = int(r["connector"])
        out[int(r["bottom_vmm"])] = (n, bottom_half)
        out[int(r["top_vmm"])] = (n, top_half)
    return out


def expected_vmms(**kw):
    """Sorted list of VMM ids that should appear in the data."""
    return sorted(vmm_half_map(**kw))


# --------------------------------------------------------------------------
# Pad map
# --------------------------------------------------------------------------
# p2-mapping's table carries pad_cx .. delta_phi; the Fe55 rates also need the
# pad area, so these come straight from the same geometry table.
_EXTRA_GEOM = ("pad_area", "mec8_connector")
_TABLE_COLS = ["vmm", "ch", "hybrid", "connector_N", "half", "sector", "strip",
               "channel_id", "pad_cx", "pad_cy", "radius", "phi", "delta_phi",
               "pad_w", "pad_h", "pad_angle", "pad_area", "mec8_connector",
               "mapped"]


def load_pad_map(map_csv=None):
    """The pad geometry indexed by channel_id (0..1279): p2-mapping's, or a
    CSV of the same layout if one is given."""
    if map_csv:
        m = pd.read_csv(map_csv).set_index("channel_id").sort_index()
    else:
        m = load_pads(SETUP.design).copy()
    if len(m) != N_PADS:
        print(f"[det1_cabling] WARNING: map has {len(m)} rows, expected {N_PADS}")
    return m


def build_pad_table(map_csv=None, strategy=WITHIN_STRATEGY,
                    bottom_half=BOTTOM_VMM_HALF, disconnected=DISCONNECTED,
                    overrides=None, base=VMM_BASE, hybrids=HYBRIDS,
                    when=CABLING_WHEN):
    """One row per instrumented (vmm, ch), joined to pad geometry.

    Columns: vmm, ch, hybrid, connector_N, half, sector, strip, channel_id,
             pad_cx, pad_cy, radius, phi, delta_phi, pad_w, pad_h, pad_angle,
             pad_area, mec8_connector, mapped.

    `strategy` is the within-half order of every half, `overrides` a
    {(connector_N, half): order} on top of it (default: the setup rows' own
    bottom_order / top_order).
    """
    over = dict(STRATEGY_OVERRIDES if overrides is None else overrides)
    setup = _setup_variant(disconnected, base, hybrids, when)
    tab = build_table(setup, when, order=strategy, bottom_half=bottom_half,
                      orders={(DETECTOR, n, h): o for (n, h), o in over.items()})
    pad = load_pad_map(map_csv)
    if map_csv:                         # geometry from the given CSV, not p2-mapping's
        tab = tab.drop(columns=[c for c in pad.columns if c in tab.columns])
        have = [c for c in pad.columns if c in _TABLE_COLS]
    else:
        have = [c for c in _EXTRA_GEOM if c in pad.columns]
    tab = tab.merge(pad[have], left_on="channel_id", right_index=True, how="left")
    tab["mapped"] = tab["pad_cx"].notna()
    tab = tab[[c for c in _TABLE_COLS if c in tab.columns]]

    # A scrambled strategy would collide two channels onto one pad; a silent
    # duplicate would double-count that pad's rate. Cheap, so always checked.
    dup = tab.loc[tab["mapped"], "channel_id"].duplicated().sum()
    if dup:
        raise ValueError(f"cabling produces {dup} duplicated channel_id -- "
                         f"strategy={strategy!r} overrides={over!r}")

    tab.attrs.update(strategy=strategy, bottom_half=bottom_half,
                     overrides=dict(over), base=base,
                     disconnected=sorted(int(c) for c in (disconnected or ())),
                     hybrids=list(hybrids), setup=SETUP.name,
                     setup_sha256=SETUP.sha256, cabling_when=str(when))
    return tab


def full_pad_map(map_csv=None, **kw):
    """All 1280 pads with an `instrumented` flag, for drawing the whole fan.

    The uninstrumented parts of the detector are real estate that exists and is
    simply not read out. Drawing them grey rather than omitting them is what
    makes a rate map honest about how much of the chamber the source covers.
    """
    pad = load_pad_map(map_csv).reset_index()
    tab = build_pad_table(map_csv=map_csv, **kw)
    live = tab.loc[tab["mapped"], ["channel_id", "vmm", "ch", "connector_N",
                                   "half", "hybrid"]]
    out = pad.merge(live, on="channel_id", how="left", suffixes=("", "_tab"))
    out["instrumented"] = out["vmm"].notna()
    return out


def attach_pads(hits, tab):
    """Left-join hits to the pad table on (vmm, ch).

    Only for small hit samples: the join materialises every pad column on every
    hit (~1.5 GB on a 5 M hit capture, as vmm_stations.auto_hot_pads documents).
    The per-pad products in vmm_fe55_map.py histogram on (vmm, ch) instead and
    attach geometry per PAD, which is 64x32 rows rather than millions.
    """
    return hits.merge(tab, on=["vmm", "ch"], how="left", copy=False)


def pad_key(vmm, ch):
    """Flat key vmm*64 + ch, the index used by every per-pad histogram here."""
    return np.asarray(vmm, dtype=np.int64) * 64 + np.asarray(ch, dtype=np.int64)


def counts_per_pad(vmm, ch, tab, weights=None):
    """Total (or weight-summed) hits per row of `tab`, aligned to tab's order.

    One bincount over the 32*64 key space; never joins geometry onto hits.
    """
    n_key = 32 * 64
    k = pad_key(vmm, ch)
    ok = (k >= 0) & (k < n_key)
    hist = np.bincount(k[ok], weights=None if weights is None
                       else np.asarray(weights, dtype=np.float64)[ok],
                       minlength=n_key)
    return hist[pad_key(tab["vmm"].to_numpy(), tab["ch"].to_numpy())]


def pad_polygons(tab):
    """(N,4,2) pad corners from centre + size + rotation.

    The P2 pads are fan tiles, each with its own angle; drawing them as rotated
    rectangles reproduces the fan instead of a square grid. Same construction as
    vmm_beam_profile.pad_polygons -- kept here so this module stands alone.
    """
    cx = tab["pad_cx"].to_numpy(float)
    cy = tab["pad_cy"].to_numpy(float)
    w = np.nan_to_num(tab["pad_w"].to_numpy(float), nan=PAD_PITCH_MM) \
        if "pad_w" in tab else np.full(len(tab), PAD_PITCH_MM)
    h = np.nan_to_num(tab["pad_h"].to_numpy(float), nan=PAD_PITCH_MM) \
        if "pad_h" in tab else np.full(len(tab), PAD_PITCH_MM)
    ang = np.nan_to_num(np.deg2rad(tab["pad_angle"].to_numpy(float)), nan=0.0) \
        if "pad_angle" in tab else np.zeros(len(tab))

    dx = np.array([-0.5, 0.5, 0.5, -0.5])
    dy = np.array([-0.5, -0.5, 0.5, 0.5])
    ca, sa = np.cos(ang)[:, None], np.sin(ang)[:, None]
    ox = dx[None, :] * w[:, None]
    oy = dy[None, :] * h[:, None]
    return np.stack([cx[:, None] + ox * ca - oy * sa,
                     cy[:, None] + ox * sa + oy * ca], axis=-1)


# --------------------------------------------------------------------------
# Does the data agree with the cabling above?
# --------------------------------------------------------------------------
def verify_cabling(vmm_ids, tab=None, disconnected=DISCONNECTED, base=VMM_BASE,
                   hybrids=HYBRIDS, quiet_frac=0.002):
    """Compare the VMM ids present in the data against the declared cabling.

    `vmm_ids` is a {vmm: n_hits} mapping (or any Series/dict-like).

    The two disconnected connectors are a free calibration of the whole
    numbering convention: each candidate convention predicts a DIFFERENT pair of
    silent VMM ids, so one file settles it. Anything that fires on a connector
    declared disconnected, or is silent on one declared live, is reported.

    Returns a dict; nothing is raised, because at the bench a surprising answer
    is information, not a reason to stop.
    """
    counts = {int(k): int(v) for k, v in dict(vmm_ids).items()}
    total = sum(counts.values()) or 1
    live = vmm_half_map(disconnected, BOTTOM_VMM_HALF, base, hybrids)
    dead = {}
    for n in hybrids:
        if n in set(disconnected or ()):
            for v in hybrid_vmms(n, base):
                dead[int(v)] = int(n)

    silent = [v for v in sorted(live) if counts.get(v, 0) / total < quiet_frac]
    loud_dead = {v: counts[v] for v in sorted(dead)
                 if counts.get(v, 0) / total >= quiet_frac}
    unknown = {v: c for v, c in sorted(counts.items())
               if v not in live and v not in dead and c / total >= quiet_frac}

    return {
        "expected_vmms": sorted(live),
        "declared_disconnected_connectors": sorted(disconnected or ()),
        "silent_but_expected": silent,
        "firing_but_declared_disconnected": loud_dead,
        "unexpected_vmm_ids": unknown,
        "n_hits_total": int(total),
        "consistent": not (silent or loud_dead or unknown),
    }


def format_cabling_report(rep):
    """One-screen summary of verify_cabling(), for the console and the JSON."""
    lines = []
    ok = "OK" if rep["consistent"] else "MISMATCH"
    lines.append(f"cabling check: {ok}   ({rep['n_hits_total']:,} hits)")
    lines.append(f"  expected VMMs: {rep['expected_vmms']}")
    if rep["silent_but_expected"]:
        lines.append(f"  SILENT but expected live : {rep['silent_but_expected']}"
                     "  <- dead hybrid, unpowered, or wrong VMM_BASE")
    if rep["firing_but_declared_disconnected"]:
        lines.append("  FIRING but declared disconnected : "
                     f"{rep['firing_but_declared_disconnected']}"
                     "  <- DISCONNECTED list or numbering convention is wrong")
    if rep["unexpected_vmm_ids"]:
        lines.append(f"  unexpected VMM ids : {rep['unexpected_vmm_ids']}")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# Live time
# --------------------------------------------------------------------------
# Every rate in Hz comes from the pcapng's own capture clock: the Enhanced
# Packet Block timestamps, which dumpcap writes from the system clock and whose
# resolution is declared in the Interface Description Block. Getting
# `if_tsresol` wrong is a silent factor-1000 on every rate, so it is parsed
# rather than assumed (the DAQ's dumpcap writes nanoseconds, not the pcapng
# default of microseconds).
#
# The SRS marker span is kept only as a CROSS-CHECK, never as the live time.
# Measured 2026-09-09 on data/hybrid_performance_data, the two agree to ~10 %
# (packet clock 19.9 s vs SRS 22.0 s at the nominal 25 ns tick on the same
# file), which is close enough to catch a gross error and far too loose to
# normalise a rate map with. It also cannot be trusted in principle:
# vmm_stations documents that `offset` is stuck in this firmware, so there is
# no BCID rollover counter and chip time is only unambiguous within 92.16 us
# of a marker. The wall clock has no such caveat.
def capture_window(path):
    """(t_first, t_last, n_packets, tsresol) from a pcapng's packet clock.

    Unix epoch seconds, from p2-decode's block walker (`walk_blocks`), which
    reads only block headers and honours each interface's if_tsresol. Simple
    Packet Blocks carry no timestamp and are counted only. `tsresol` is
    returned as None: p2-decode converts to ns per interface and does not
    expose it, and nothing downstream uses it beyond the per-file record.
    """
    walked = walk_blocks(path, times=True)
    t_ns = walked[4]
    del walked
    n_pkt = int(t_ns.size)
    t_ns = t_ns[np.isfinite(t_ns)]
    if t_ns.size == 0:
        return None, None, n_pkt, None
    return float(t_ns[0]) / 1e9, float(t_ns[-1]) / 1e9, n_pkt, None


def live_time(paths, capture_duration_s=None, srs_span_s=None):
    """Live time in seconds over one or more captures, plus how it was got.

    The number reported is the UNION of the files' [t_first, t_last] intervals.
    Union, not a sum of spans: captures handed in twice, or overlapping ring
    segments, must not each be credited with their own live time. Duplicated
    files are common enough at a bench (a copied directory, a re-listed glob)
    that silently doubling the denominator would halve every rate.

    A first-to-last-packet span always UNDERCOUNTS by the dead interval before
    the first packet and after the last -- at a few Hz of Fe55 with packets
    arriving sparsely, that is not negligible. So when the DAQ's configured
    capture length is known (run_config_beam.CAPTURE_DURATION_S), the nominal
    n_files * capture_duration_s is reported alongside and preferred when it is
    larger, with the ratio kept so the caller can see how much was inferred.

    `srs_span_s` is an optional cross-check from the SRS marker span; it is
    never used for the value, only compared (see the note above).

    Returns dict(live_s, span_s, nominal_s, t_start, t_end, n_files,
                 n_packets, method, warnings, per_file).
    """
    paths = [str(p) for p in ([paths] if isinstance(paths, (str, bytes))
                              else list(paths))]
    per, warnings = [], []
    for p in paths:
        try:
            t0, t1, n, res = capture_window(p)
        except Exception as e:                      # a half-written file
            per.append({"path": p, "error": f"{type(e).__name__}: {e}"})
            continue
        if t0 is None:
            per.append({"path": p, "error": "no packets"})
            continue
        per.append({"path": p, "t_first": t0, "t_last": t1,
                    "span_s": t1 - t0, "n_packets": n, "tsresol": res})

    good = [d for d in per if "t_first" in d]
    if not good:
        return {"live_s": 0.0, "span_s": 0.0, "nominal_s": None,
                "t_start": None, "t_end": None, "n_files": 0, "n_packets": 0,
                "method": "no usable capture",
                "warnings": ["no capture yielded a packet timestamp"],
                "per_file": per}

    good.sort(key=lambda d: d["t_first"])

    # Union of the covered intervals.
    merged = []
    for d in good:
        if merged and d["t_first"] <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], d["t_last"])
        else:
            merged.append([d["t_first"], d["t_last"]])
    span = float(sum(b - a for a, b in merged))

    n_distinct = len(merged)
    if n_distinct < len(good):
        warnings.append(
            f"{len(good) - n_distinct} of {len(good)} captures overlap in time "
            "(duplicated or ring-overlapping files); counted once")

    live, method = span, f"packet clock, union of {n_distinct} interval(s)"
    nominal = None
    if capture_duration_s:
        nominal = float(capture_duration_s) * n_distinct
        if nominal > span:
            live, method = nominal, (
                f"configured capture duration {capture_duration_s:.1f} s x "
                f"{n_distinct} (packet-clock union was {span:.2f} s)")

    if srs_span_s and span > 0:
        ratio = srs_span_s / span
        if not 0.5 < ratio < 2.0:
            warnings.append(
                f"SRS marker span {srs_span_s:.1f} s disagrees with the packet "
                f"clock {span:.1f} s by {ratio:.1f}x -- rates use the packet "
                "clock; check the capture is not truncated")

    return {"live_s": float(live), "span_s": span, "nominal_s": nominal,
            "t_start": float(merged[0][0]), "t_end": float(merged[-1][1]),
            "n_files": n_distinct, "n_packets": sum(d["n_packets"] for d in good),
            "method": method, "warnings": warnings, "per_file": per}


if __name__ == "__main__":
    tab = build_pad_table()
    print(f"instrumented (vmm,ch): {len(tab)}   mapped pads: {tab.mapped.sum()}")
    print(f"VMMs: {sorted(tab.vmm.unique())}")
    print(f"connectors: {sorted(tab.connector_N.unique())}  "
          f"sectors: {sorted(tab.sector.unique())}")
    full = full_pad_map()
    print(f"detector pads: {len(full)}  instrumented: {int(full.instrumented.sum())}")
    live = tab[tab.mapped]
    print(f"active area: x {live.pad_cx.min():.0f}..{live.pad_cx.max():.0f} mm, "
          f"y {live.pad_cy.min():.0f}..{live.pad_cy.max():.0f} mm, "
          f"pad area {live.pad_area.min():.1f}..{live.pad_area.max():.1f} mm^2")
