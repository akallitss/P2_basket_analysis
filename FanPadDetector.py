
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FanPadDetector.py

Geometry model for the P2 BASKET fan-pad Micromegas detector.

Each pad is drawn as its true Gerber rectangle: centred on the measured pad
centroid (pad_cx, pad_cy), sized to the measured pad extent (pad_w tangential,
pad_h radial ≈ 12.0 × 11.6 mm) and rotated to the per-sector pad_angle.  This
reproduces the physical layout including the real ~0.4 mm inter-pad gap.  These
columns are produced by P2_mapping.parse_pad_regions from the F_Cu region fills;
if they are absent the model falls back to the strip endpoint (x, y) with a
synthetic square of side pad_size_mm.

Channel mapping is embedded at load time:
  (vmm_id, ch) → pad index
with support for four connector orientations — 'normal', 'flipped'
(end-to-end pin reverse), 'back' (back-side insertion = adjacent pair swap),
and 'flipped_back' (both) — to handle unknown FEB pin-order direction and
connector seating.  See pin_to_ch().

Usage
-----
    from FanPadDetector import FanPadDetector

    det = FanPadDetector.from_mapping_csv(
        csv_path    = "Detector_Mapping/P2_BASKET/P2_BASKET_mapping.csv",
        mec8_to_vmm = {(6, 1): 12, (6, 0): 13, (7, 1): 14, (7, 0): 15},
        fpc_connectors = [6, 7],
        orientation = "normal",   # or "flipped"
        pad_size_mm = 14.0,       # pad pitch → pads tile edge-to-edge
    )
    fig, axes = det.plot_hit_heatmap(df_hits_per_channel)
    fig, axes = det.plot_efficiency_map(df_ch_eff)

@author: ak271430
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ─────────────────────────────────────────────────────────────
# CONNECTOR ORIENTATION
# ─────────────────────────────────────────────────────────────

# A flat ribbon / FPC connection has two independent geometric toggles:
#   • end-to-end reverse   ("flipped")  — cable plugged with pin order reversed
#   • back-side insertion  ("back")     — 2-row connector mirrored about its long
#                                          axis, swapping adjacent pin pairs
#                                          (1↔2, 3↔4, …  i.e. ch ^ 1)
# Both are XOR masks on the base channel ch0 = mec8_pin − 3, so they commute and
# compose into four orientations.
ORIENTATIONS = ("normal", "flipped", "back", "flipped_back")


def pin_to_ch(mec8_pin, orientation="normal"):
    """
    Map a MEC8 pin (3–66) to a VMM channel (0–63) for a given connector
    orientation.  Returns a channel that may be outside 0–63 for pins 67/68
    or invalid combinations; callers must range-check.

    orientation : {'normal', 'flipped', 'back', 'flipped_back'}
        'normal'        → ch = pin − 3
        'flipped'       → end-to-end reverse        (ch = 63 − ch0 = 66 − pin)
        'back'          → back-side plug, pair swap  (ch = ch0 ^ 1)
        'flipped_back'  → both                       (ch = 62 ^ ch0)
    """
    ch = mec8_pin - 3                       # base
    if "flipped" in orientation:            # end-to-end reverse
        ch = 63 - ch
    if "back" in orientation and 0 <= ch <= 63:   # back-side pair swap
        ch = ch ^ 1
    return ch


# ─────────────────────────────────────────────────────────────
# STRIP-BASED CHANNEL ORDERING  (validated against M3 tracks)
# ─────────────────────────────────────────────────────────────

# The cosmic-bench DREAM readout was validated against M3 telescope tracks
# (cosmic_bench_analysis/p2_mapping.py, 03_m3_alignment.py): the order of the
# 64 channels inside a connector half is 'reverse'. Unlike pin_to_ch, these
# strategies treat the VMM channel as the within-half index and read the
# geometry straight from the CSV's own `strip` column, so they are immune to
# the 2-pin (ground) discontinuity in the MEC8 pin numbering — that gap makes
# the validated order impossible to reproduce with any single ch = f(mec8_pin)
# orientation.
STRIP_STRATEGIES = ("reverse", "linear", "pairswap")


def strip_to_ch(strip, strategy="reverse"):
    """
    Map a physical strip (1–128) to a VMM channel (0–63), treating the VMM
    channel as the within-connector-half index.

    base = 1 for the bottom half (strips 1–64), 65 for the top half
    (strips 65–128); within = strip − base (0–63).

    strategy : {'reverse', 'linear', 'pairswap'}
        'reverse'  → ch = 63 − within   ** VALIDATED against M3 tracks **
        'linear'   → ch = within
        'pairswap' → ch = within ^ 1
    """
    base   = 1 if strip <= 64 else 65
    within = strip - base
    if strategy == "reverse":
        return 63 - within
    if strategy == "linear":
        return within
    if strategy == "pairswap":
        return within ^ 1
    raise ValueError(f"unknown strip strategy {strategy!r}")


def _fanpad_from_row(row, columns, vmm_id, ch, apex, fallback_size):
    """
    Build a FanPad from a mapping-CSV row.

    Uses the true Gerber pad geometry (pad_cx/pad_cy centroid, pad_w/pad_h
    extent) when those columns are present; otherwise falls back to the strip
    endpoint (x, y) with a synthetic square of side ``fallback_size``.
    """
    angle = (np.deg2rad(float(row["pad_angle"]))
             if "pad_angle" in columns and pd.notna(row["pad_angle"]) else None)
    via_x = (float(row["via_x"])
             if "via_x" in columns and pd.notna(row["via_x"]) else None)
    via_y = (float(row["via_y"])
             if "via_y" in columns and pd.notna(row["via_y"]) else None)

    has_geo = (all(c in columns for c in ("pad_cx", "pad_cy", "pad_w", "pad_h"))
               and pd.notna(row["pad_cx"]))
    if has_geo:
        cx, cy = float(row["pad_cx"]), float(row["pad_cy"])
        w,  h  = float(row["pad_w"]),  float(row["pad_h"])
    else:
        cx, cy = float(row["x"]), float(row["y"])
        w = h  = fallback_size

    return FanPad(cx, cy, w, h, angle, vmm_id, ch,
                  via_x=via_x, via_y=via_y, apex=apex)


# ─────────────────────────────────────────────────────────────
# PAD PRIMITIVE
# ─────────────────────────────────────────────────────────────

class FanPad:
    """
    Single detector pad drawn as its true Gerber rectangle.

    The pad centre (``pad_cx, pad_cy``), tangential/radial extents
    (``pad_w``/``pad_h``) and orientation (``pad_angle``) all come from the
    F_Cu region fill measured directly from the Gerber (see
    P2_mapping.parse_pad_regions).  Drawing each pad at its real centroid with
    its real size reproduces the physical layout — including the genuine ~0.4 mm
    inter-pad gap — rather than tiling synthetic squares.

    Parameters
    ----------
    cx, cy : float
        True pad centroid (mm) from the Gerber region fill.
    w, h : float
        Pad extent (mm) along the pad_angle axis (w, tangential) and
        perpendicular to it (h, radial).
    angle : float or None
        Pad orientation in radians (per-sector pad_angle).
    vmm_id, ch : int
        VMM id and channel that read this pad.
    via_x, via_y : float or None
        Via position; used only as an orientation fallback.
    apex : (float, float) or None
        Fan apex; used only as an orientation fallback.
    """
    __slots__ = ("cx", "cy", "w", "h", "angle",
                 "via_x", "via_y", "apex", "vmm_id", "ch")

    def __init__(self, cx, cy, w, h, angle, vmm_id, ch,
                 via_x=None, via_y=None, apex=None):
        self.cx     = cx
        self.cy     = cy
        self.w      = w
        self.h      = h
        self.angle  = angle
        self.vmm_id = vmm_id
        self.ch     = ch
        self.via_x  = via_x
        self.via_y  = via_y
        self.apex   = apex

    @property
    def center(self):
        return (self.cx, self.cy)

    @property
    def rotation(self):
        """
        Angle (rad) used to orient the pad rectangle.

        Priority:
          1. ``angle`` — the per-sector pad_angle from the mapping CSV.
          2. apex→centre — smooth radial fallback.
          3. via→centre — last-resort per-pad fallback.
        """
        if self.angle is not None:
            return self.angle
        if self.apex is not None:
            return np.arctan2(self.cy - self.apex[1],
                              self.cx - self.apex[0])
        if self.via_x is not None:
            return np.arctan2(self.cy - self.via_y,
                              self.cx - self.via_x)
        return 0.0

    def vertices(self):
        """4 corners of the pad rectangle (pad_w × pad_h) in detector coords."""
        hw  = 0.5 * self.w
        hh  = 0.5 * self.h
        rot = self.rotation
        corners = np.array([
            [-hw, -hh], [hw, -hh], [hw, hh], [-hw, hh],
        ], dtype=float)
        c, s = np.cos(rot), np.sin(rot)
        R = np.array([[c, -s], [s, c]])
        return corners @ R.T + np.array([self.cx, self.cy])


# ─────────────────────────────────────────────────────────────
# DETECTOR MODEL
# ─────────────────────────────────────────────────────────────

class FanPadDetector:
    """
    P2 BASKET fan-pad Micromegas detector model.

    Attributes
    ----------
    pads : list of FanPad
        All active pads, in load order.
    channel_map : dict
        (vmm_id, ch) → pad index in self.pads.
    apex : (float, float)
        Fan apex position (mm) used for polar coordinates.
    orientation : str
        'normal' or 'flipped' — VMM channel assignment direction.
    """

    def __init__(self):
        self.pads        = []
        self.channel_map = {}
        self.apex        = (0.0, 0.0)
        self.orientation = "normal"

    # ── Construction ─────────────────────────────────────────

    @classmethod
    def from_mapping_csv(cls, csv_path, mec8_to_vmm, fpc_connectors,
                         orientation="normal", pad_size_mm=14.0,
                         half_width_mm=None, strategy=None):
        """
        Build detector from P2_BASKET_mapping.csv.

        The channel→pad ordering is chosen by one of two independent families:

        * ``strategy`` in {'reverse', 'linear', 'pairswap'} — the VMM channel is
          the within-connector-half index and the geometry is read from the
          CSV's ``strip`` column (see strip_to_ch). ``'reverse'`` is the order
          validated against M3 tracks in the cosmic-bench DREAM readout; prefer
          it. Immune to the MEC8 2-pin numbering gap.
        * ``strategy=None`` (default) — legacy behaviour: the channel is derived
          from ``mec8_pin`` via ``pin_to_ch(pin, orientation)``. Kept so the
          ``orientation`` diagnostic (normal/flipped/back/flipped_back) still
          works, but note no ``ch = f(mec8_pin)`` formula can reproduce the
          validated order exactly.

        Parameters
        ----------
        csv_path : str
            Path to CSV generated by P2_mapping.py (must have via_x, via_y).
        mec8_to_vmm : dict
            {(fpc_connector, mec8_connector): vmm_id}
        fpc_connectors : list of int
            FPC connector numbers that were cabled (subset of 0–9).
        orientation : {'normal', 'flipped', 'back', 'flipped_back'}
            Pin-based connector orientation (used only when strategy is None;
            see pin_to_ch).
        pad_size_mm : float
            Full side length of each square pad (mm). Default 14 mm (slightly
            above the ~11-12 mm pitch) so neighbouring pads tile edge-to-edge.
        half_width_mm : float or None
            Deprecated alias kept for backward compatibility; if given,
            pad_size_mm = 2 × half_width_mm.
        strategy : {'reverse', 'linear', 'pairswap'} or None
            Strip-based within-half ordering. Overrides ``orientation`` when set.
        """
        if half_width_mm is not None:
            pad_size_mm = 2.0 * half_width_mm

        use_strip = strategy in STRIP_STRATEGIES
        if strategy is not None and not use_strip:
            raise ValueError(f"unknown strip strategy {strategy!r}; "
                             f"choose one of {STRIP_STRATEGIES} or None")

        df = pd.read_csv(csv_path)
        df = df[df["connector_number"].isin(fpc_connectors)].copy()
        subset = ["x", "y", "via_x", "via_y", "mec8_connector", "mec8_pin"]
        if use_strip:
            subset.append("strip")
        df = df.dropna(subset=subset)

        det = cls()
        det.orientation = strategy if use_strip else orientation

        # Recover apex from polar coordinates if available
        if "phi" in df.columns and "radius" in df.columns:
            row0 = df.iloc[0]
            det.apex = (
                row0["x"] - row0["radius"] * np.cos(row0["phi"]),
                row0["y"] - row0["radius"] * np.sin(row0["phi"]),
            )

        for _, row in df.iterrows():
            fpc  = int(row["connector_number"])
            mec8 = int(row["mec8_connector"])

            vmm_id = mec8_to_vmm.get((fpc, mec8))
            if vmm_id is None:
                continue

            if use_strip:
                ch = strip_to_ch(int(row["strip"]), strategy)
            else:
                ch = pin_to_ch(int(row["mec8_pin"]), orientation)

            if not (0 <= ch <= 63):
                continue

            pad = _fanpad_from_row(row, df.columns, vmm_id, ch,
                                   det.apex, pad_size_mm)
            idx = len(det.pads)
            det.pads.append(pad)
            det.channel_map[(vmm_id, ch)] = idx

        print(f"FanPadDetector: loaded {len(det.pads)} pads "
              f"({'strategy' if use_strip else 'orientation'}={det.orientation})")
        return det

    @classmethod
    def from_mapping_csv_full(cls, csv_path, pad_size_mm=14.0,
                              fpc_connectors=None, half_width_mm=None):
        """
        Load the full detector geometry for display and geometry testing.

        No VMM wiring required.  Uses connector_number as pseudo-vmm_id and
        (strip − 1) as ch so the full 1280-pad fan can be visualised.

        Parameters
        ----------
        csv_path : str
            Path to CSV generated by P2_mapping.py (must have via_x, via_y).
        pad_size_mm : float
            Full side length of each square pad (mm). Default 14 mm (slightly
            above the ~11-12 mm pitch) so neighbouring pads tile edge-to-edge.
        fpc_connectors : list of int or None
            Subset of connectors to load (default: all present in the CSV).
        half_width_mm : float or None
            Deprecated alias; if given, pad_size_mm = 2 × half_width_mm.
        """
        if half_width_mm is not None:
            pad_size_mm = 2.0 * half_width_mm
        df = pd.read_csv(csv_path)
        if fpc_connectors is not None:
            df = df[df["connector_number"].isin(fpc_connectors)].copy()
        df = df.dropna(subset=["x", "y", "via_x", "via_y"])

        det = cls()
        det.orientation = "geometry_only"

        if "phi" in df.columns and "radius" in df.columns:
            row0 = df.iloc[0]
            det.apex = (
                row0["x"] - row0["radius"] * np.cos(row0["phi"]),
                row0["y"] - row0["radius"] * np.sin(row0["phi"]),
            )

        for _, row in df.iterrows():
            fpc   = int(row["connector_number"])
            strip = int(row["strip"]) if "strip" in df.columns else 1
            ch    = strip - 1          # 0-indexed strip within connector (0–127)

            pad = _fanpad_from_row(row, df.columns, fpc, ch,
                                   det.apex, pad_size_mm)
            idx = len(det.pads)
            det.pads.append(pad)
            det.channel_map[(fpc, ch)] = idx

        n_conn = df["connector_number"].nunique()
        print(f"FanPadDetector: loaded {len(det.pads)} pads "
              f"({n_conn} FPC connectors, geometry_only)")
        return det

    # ── Plotting helpers ──────────────────────────────────────

    def _draw_pads(self, ax, values, cmap, norm, missing_color="lightgrey"):
        """
        Draw all pads as rotated rectangles coloured by values[i].

        Parameters
        ----------
        values : array-like, length = len(self.pads)
            One value per pad; np.nan → missing_color.
        """
        for pad, val in zip(self.pads, values):
            verts = pad.vertices()
            color = (cmap(norm(val))
                     if np.isfinite(val) and val >= 0
                     else missing_color)
            poly = plt.Polygon(verts,
                               facecolor=color,
                               edgecolor="black",
                               linewidth=0.3,
                               zorder=2)
            ax.add_patch(poly)

        # Auto-set limits
        all_x = [p.cx for p in self.pads]
        all_y = [p.cy for p in self.pads]
        margin = 20
        ax.set_xlim(min(all_x) - margin, max(all_x) + margin)
        ax.set_ylim(min(all_y) - margin, max(all_y) + margin)
        ax.set_aspect("equal")
        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")

    # ── Public plot methods ───────────────────────────────────

    def plot_hit_heatmap(self, df_ch_stats, run_no=None,
                         detector_label="P2 Large Detector",
                         out_dir=None):
        """
        2×2 plot of hit count (left) and total ADC (right) per pad, with a
        linear colour scale on the top row and a log scale on the bottom row.

        Parameters
        ----------
        df_ch_stats : pd.DataFrame
            Must have columns vmm_id, ch, n_hits, adc_sum.
            Output of compute_pad_hit_stats() accumulated over all files.
        run_no : int or None
        detector_label : str
        out_dir : str or None
        """
        from matplotlib.colors import LogNorm

        hit_map = {}
        adc_map = {}
        for _, row in df_ch_stats.iterrows():
            key = (int(row["vmm_id"]), int(row["ch"]))
            idx = self.channel_map.get(key)
            if idx is not None:
                hit_map[idx] = row["n_hits"]
                adc_map[idx] = row["adc_sum"]

        n = len(self.pads)
        hit_vals = np.array([hit_map.get(i, np.nan) for i in range(n)])
        adc_vals = np.array([adc_map.get(i, np.nan) for i in range(n)])

        fig, axes = plt.subplots(2, 2, figsize=(18, 17))

        for col, (vals, title) in enumerate([
            (hit_vals, "Hit count"),
            (adc_vals, "Total ADC (charge)"),
        ]):
            valid = vals[np.isfinite(vals)]
            vmax  = valid.max() if len(valid) else 1.0

            # ── linear (top row) ──
            ax_lin = axes[0, col]
            norm   = plt.Normalize(vmin=0, vmax=vmax)
            self._draw_pads(ax_lin, vals, plt.cm.viridis, norm)
            sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
            sm.set_array([])
            fig.colorbar(sm, ax=ax_lin, shrink=0.8, label=title)
            ax_lin.set_title(f"{title} — linear", fontsize=11)

            # ── log (bottom row); non-positive pads shown as missing ──
            ax_log   = axes[1, col]
            log_vals = np.where(np.isfinite(vals) & (vals > 0), vals, np.nan)
            positive = log_vals[np.isfinite(log_vals)]
            vmin_log = positive.min() if len(positive) else 1.0
            log_norm = LogNorm(vmin=vmin_log, vmax=max(vmax, vmin_log * 10))
            self._draw_pads(ax_log, log_vals, plt.cm.viridis, log_norm)
            sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=log_norm)
            sm.set_array([])
            fig.colorbar(sm, ax=ax_log, shrink=0.8, label=title)
            ax_log.set_title(f"{title} — log", fontsize=11)

        run_str = f"  |  Run {run_no}" if run_no is not None else ""
        fig.suptitle(f"{detector_label} — pad hit heatmaps{run_str}",
                     fontsize=13, fontweight="bold")
        plt.tight_layout()

        if out_dir is not None:
            import os
            os.makedirs(out_dir, exist_ok=True)
            stem = f"hit_map_{detector_label.replace(' ', '_')}_run{run_no}"
            for ext in ("png", "pdf"):
                fig.savefig(os.path.join(out_dir, f"{stem}.{ext}"),
                            bbox_inches="tight")
        return fig, axes

    def plot_efficiency_map(self, df_ch_eff, run_no=None,
                            detector_label="P2 Large Detector",
                            out_dir=None):
        """
        Per-pad hit-rate map coloured by true_efficiency.

        Parameters
        ----------
        df_ch_eff : pd.DataFrame
            Must have columns vmm_id, ch, true_efficiency.
        """
        eff_map = {}
        for _, row in df_ch_eff.iterrows():
            key = (int(row["vmm_id"]), int(row["ch"]))
            idx = self.channel_map.get(key)
            if idx is not None:
                val = max(0.0, float(row["true_efficiency"]))
                eff_map[idx] = val

        n = len(self.pads)
        eff_vals = np.array([eff_map.get(i, np.nan) for i in range(n)])

        valid = eff_vals[np.isfinite(eff_vals)]
        vmax  = valid.max() if len(valid) else 1.0
        cmap  = plt.cm.RdYlGn
        norm  = plt.Normalize(vmin=0, vmax=vmax)

        fig, ax = plt.subplots(figsize=(10, 9))
        self._draw_pads(ax, eff_vals, cmap, norm)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, shrink=0.8)
        cb.set_label("Hit rate per trigger (relative, normalised to data max)")
        cb.formatter = plt.FuncFormatter(lambda x, _: f"{x*100:.3f}%")
        cb.update_ticks()

        run_str = f"  |  Run {run_no}" if run_no is not None else ""
        ax.set_title(
            f"{detector_label} — pad hit-rate map{run_str}\n"
            f"(max = {vmax*100:.3f}%)",
            fontsize=12, fontweight="bold",
        )
        plt.tight_layout()

        if out_dir is not None:
            import os
            os.makedirs(out_dir, exist_ok=True)
            stem = f"eff_map_{detector_label.replace(' ', '_')}_run{run_no}"
            for ext in ("png", "pdf"):
                fig.savefig(os.path.join(out_dir, f"{stem}.{ext}"),
                            bbox_inches="tight")
        return fig, ax

    def plot_geometry(self, ax=None, color="lightblue", label_pads=False):
        """Quick geometry check: draw all pads with uniform colour."""
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 9))
        dummy = np.zeros(len(self.pads))
        cmap = plt.cm.Blues
        norm = plt.Normalize(vmin=0, vmax=1)
        self._draw_pads(ax, dummy, cmap, norm, missing_color=color)
        if label_pads:
            for pad in self.pads:
                cx, cy = pad.center
                ax.text(cx, cy, f"({pad.vmm_id},{pad.ch})",
                        fontsize=3, ha="center", va="center")
        ax.set_title(f"FanPadDetector geometry — {len(self.pads)} pads  "
                     f"(orientation={self.orientation})",
                     fontsize=11)
        return ax
