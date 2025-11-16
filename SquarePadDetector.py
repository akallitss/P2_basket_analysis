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
    def __init__(self, half_width, x=0, y=0):
        self.half_width = half_width
        self.x = x
        self.y = y

    def __repr__(self):
        return f"SquarePad(half_width={self.half_width}, x={self.x:.2f}, y={self.y:.2f})"


class SquareDetector:
    def __init__(self, pad_half_width, x=0, y=0, rotation=0):
        self.square_pads = []
        self.pad_half_width = pad_half_width
        self.x = x
        self.y = y
        self.rotation = rotation

        self.add_pad()

    def set_rotation(self, rotation):
        self.rotation = rotation

    def set_center(self, x, y):
        self.x = x
        self.y = y

    def add_pad(self, reference_pad_index=None, dx_index=1, dy_index=0):
        """Add a square pad relative to an existing pad. Default is one pad width right."""
        if reference_pad_index is None:
            if not self.square_pads:
                self.square_pads.append(SquarePad(self.pad_half_width, 0, 0))
            else:
                raise ValueError("No reference pad index given, but first pad already exists.")
            return

        ref_pad = self.square_pads[reference_pad_index]
        dx = dx_index * 2 * self.pad_half_width
        dy = dy_index * 2 * self.pad_half_width
        new_x = ref_pad.x + dx
        new_y = ref_pad.y + dy

        if any(np.isclose(pad.x, new_x) and np.isclose(pad.y, new_y) for pad in self.square_pads):
            raise ValueError("Pad to be added already exists.")

        self.square_pads.append(SquarePad(self.pad_half_width, new_x, new_y))

    def get_pad_center(self, pad_index):
        """Get the center of a pad with global rotation and translation."""
        pad = self.square_pads[pad_index]
        cos_theta = np.cos(self.rotation)
        sin_theta = np.sin(self.rotation)
        rotated_x = cos_theta * pad.x - sin_theta * pad.y
        rotated_y = sin_theta * pad.x + cos_theta * pad.y
        return rotated_x + self.x, rotated_y + self.y

    def plot_detector(self, global_coords=False, ax_in=None, zorder=10, pad_colors='lightgreen', pad_alpha=0.5):
        if ax_in is None:
            fig, ax = plt.subplots()
        else:
            ax = ax_in

        for i, pad in enumerate(self.square_pads):
            if global_coords:
                x, y = self.get_pad_center(i)
                rotation = self.rotation
            else:
                x, y = pad.x, pad.y
                rotation = 0
            if isinstance(pad_colors, str):
                pad_color = pad_colors
            elif isinstance(pad_colors, dict):
                pad_color = pad_colors.get(i, 'none')
            square = plt.Polygon(square_vertices(x, y, self.pad_half_width, rotation), edgecolor='black',
                                 facecolor=pad_color, alpha=pad_alpha, zorder=zorder)
            ax.add_patch(square)
            ax.scatter(x, y, s=1, color='black', zorder=zorder)

        if ax_in is None:
            ax.set_aspect('equal')
            margin = self.pad_half_width * 2
            xs = [p.x for p in self.square_pads]
            ys = [p.y for p in self.square_pads]
            ax.set_xlim(min(xs) - margin, max(xs) + margin)
            ax.set_ylim(min(ys) - margin, max(ys) + margin)
            # plt.show()

    def __repr__(self):
        return f"SquareDetector with {len(self.square_pads)} pads:\n" + "\n".join(map(str, self.square_pads))


def square_vertices(x, y, half_width, rotation=0):
    """Return the 4 corners of a square centered at (x, y), rotated by rotation (rad)."""
    corners = np.array([
        [-half_width, -half_width],
        [half_width, -half_width],
        [half_width, half_width],
        [-half_width, half_width],
    ], dtype=float)

    if rotation != 0:
        rotation_matrix = np.array([
            [np.cos(rotation), -np.sin(rotation)],
            [np.sin(rotation),  np.cos(rotation)]
        ])
        corners = corners @ rotation_matrix.T

    corners += np.array([x, y])
    return corners
