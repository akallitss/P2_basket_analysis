#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 15/11/2025 19:05
Created in PyCharm
Created as SquarePadDetector.py

@author: akallits
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class SquarePad:
    def __init__(self, half_width, x=0, y=0, rotation=0):
        self.half_width = half_width
        self.x = x
        self.y = y
        self.rotation = rotation   # per-pad rotation (radians)

    def __repr__(self):
        return (f"SquarePad(hw={self.half_width}, "
                f"x={self.x:.2f}, y={self.y:.2f}, rot={self.rotation:.3f})")


class SquareDetector:
    def __init__(self, x=0, y=0, rotation=0):
        # No forced first pad — user may add any pad they want
        self.square_pads = []
        self.x = x
        self.y = y
        self.rotation = rotation   # global rotation

    # --------------------------
    #  POSITIONING / GEOMETRY
    # --------------------------

    def set_rotation(self, rotation):
        self.rotation = rotation

    def set_center(self, x, y):
        self.x = x
        self.y = y

    # --------------------------
    #   ADDING PADS
    # --------------------------

    def add_pad_absolute(self, x, y, half_width, rotation=0):
        """Add a pad specifying full geometry directly."""
        self.square_pads.append(SquarePad(half_width, x, y, rotation))

    def add_pad_relative(self, ref_index, dx, dy, half_width, rotation=0):
        """
        Add pad relative to an existing pad:
            (dx, dy) is in same coordinate system as reference pad.
        """
        ref = self.square_pads[ref_index]
        new_x = ref.x + dx
        new_y = ref.y + dy

        self.square_pads.append(SquarePad(half_width, new_x, new_y, rotation))

    # --------------------------
    #  COORDINATES
    # --------------------------

    def get_pad_center(self, pad_index):
        """Global coordinates (after detector rotation + translation)."""
        pad = self.square_pads[pad_index]

        # first rotate pad local position by detector rotation
        c = np.cos(self.rotation)
        s = np.sin(self.rotation)
        rotated_x = c * pad.x - s * pad.y
        rotated_y = s * pad.x + c * pad.y

        return rotated_x + self.x, rotated_y + self.y

    def get_pad_global_rotation(self, pad_index):
        """Pad rotation including global detector rotation."""
        return self.rotation + self.square_pads[pad_index].rotation

    # --------------------------
    #  PLOTTING
    # --------------------------

    def plot_detector(self, global_coords=False, ax_in=None,
                      zorder=10, pad_color='lightgreen', pad_alpha=0.5):

        if ax_in is None:
            fig, ax = plt.subplots()
        else:
            ax = ax_in

        for i, pad in enumerate(self.square_pads):

            if global_coords:
                x, y = self.get_pad_center(i)
                rotation = self.get_pad_global_rotation(i)
            else:
                x, y = pad.x, pad.y
                rotation = pad.rotation

            verts = square_vertices(x, y, pad.half_width, rotation)

            square = plt.Polygon(
                verts, edgecolor='black', facecolor=pad_color,
                alpha=pad_alpha, zorder=zorder
            )

            ax.add_patch(square)
            ax.scatter([x], [y], s=5, color='black')

        if ax_in is None:
            ax.set_aspect('equal')
            # plt.show()

    def __repr__(self):
        text = "SquareDetector:\n"
        for i, p in enumerate(self.square_pads):
            text += f"  {i}: {p}\n"
        return text

# ---------------------------------------------------------------
# Mapping file
# ---------------------------------------------------------------

    def write_mapping(self, filename="detector_map.csv"):
        # Write CSV mapping file for the detector pads

        data = []

        for i, pad in enumerate(self.square_pads):
            size = 2 * pad.half_width
            rotation_deg = getattr(pad, "rotation", 0) * 180 / np.pi

            data.append({
                "PadID": i,
                "X_mm": round(pad.x, 3),
                "Y_mm": round(pad.y, 3),
                "Size_mm": round(size, 3),
                "Rotation_deg": round(rotation_deg, 3)
            })

        df = pd.DataFrame(data)
        df.to_csv(filename, index=False)

        print(f"\n✔ CSV mapping saved to {filename}\n")

    def read_mapping(self, filename):
        # Read CSV mapping file and populate the detector pads

        df = pd.read_csv(filename)

        self.square_pads = []  # Clear existing pads

        for _, row in df.iterrows():
            pad_id = int(row["PadID"])
            x = float(row["X_mm"])
            y = float(row["Y_mm"])
            size = float(row["Size_mm"])
            rotation_deg = float(row.get("Rotation_deg", 0))
            rotation_rad = rotation_deg * np.pi / 180

            half_width = size / 2

            self.add_pad_absolute(x, y, half_width, rotation_rad)

        print(f"\n✔ Loaded {len(self.square_pads)} pads from {filename}\n")

    def plot_hit_heatmap(self, df, ax=None, cmap='viridis', area_norm=False, adc_weighted=False, global_coords=False):
        """
        Plot pads colored by hit_count using a colormap.

        df must contain columns: PadID, X_mm, Y_mm, hit_count
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 8))

        if area_norm:
            # Normalize hit counts by pad area
            area_map = {i: (4 * pad.half_width ** 2) for i, pad in enumerate(self.square_pads)}
            df = df.copy()
            df["hit_count"] = df.apply(lambda row: row["hit_count"] / area_map.get(row["PadID"], 1), axis=1)
        if adc_weighted:
            # Weight hit counts by ADC value if available
            if "adc" in df.columns:
                df = df.copy()
                df["hit_count"] = df["hit_count"] * df["adc"]

        # Merge df by PadID (ensure correct ordering)
        hit_map = df.groupby("PadID")["hit_count"].sum().to_dict()
        vmax = max(hit_map.values())

        for pad_id, pad in enumerate(self.square_pads):

            # Default to 0 hits if missing
            hits = hit_map.get(pad_id, 0)

            # Normalized 0-1 value for colormap
            color_val = hits / vmax

            # Use global or local coordinates
            if global_coords:
                x, y = self.get_pad_center(pad_id)
                rot = self.get_pad_global_rotation(pad_id)
            else:
                x, y = pad.x, pad.y
                rot = pad.rotation

            verts = square_vertices(x, y, pad.half_width, rot)

            square = plt.Polygon(
                verts,
                edgecolor='black',
                facecolor=plt.get_cmap(cmap)(color_val),
                alpha=0.9
            )
            ax.add_patch(square)

        # Scatter pad centers
        xs = df["X_mm"]
        ys = df["Y_mm"]
        sc = ax.scatter(xs, ys, c=df["hit_count"], cmap=cmap, s=5)

        plt.colorbar(sc, ax=ax, label="Hit count")
        ax.set_aspect("equal")
        ax.set_xlabel("X[mm]")
        ax.set_ylabel("Y[mm]")
        ax.set_title("Hit heatmap per pad")
        return ax


# -----------------------------------------
# helper for square vertices
# -----------------------------------------

def square_vertices(x, y, half_width, rotation=0):
    corners = np.array([
        [-half_width, -half_width],
        [half_width, -half_width],
        [half_width,  half_width],
        [-half_width,  half_width],
    ])

    # local pad rotation
    if rotation != 0:
        c = np.cos(rotation)
        s = np.sin(rotation)
        R = np.array([[c, -s], [s, c]])
        corners = corners @ R.T

    # translate to position
    corners += np.array([x, y])
    return corners


