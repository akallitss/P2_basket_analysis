#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 15/11/2025 19:05
Created in PyCharm
Created as SquarePadDetector.py

@author: akallits
"""

import numpy as np
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
