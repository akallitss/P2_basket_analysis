#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 15/11/2025 19:06
Created in PyCharm
Created as SquarePadDetectorTest.py

@author: akallits
"""
import matplotlib.pyplot as plt

from SquarePadDetector import SquareDetector

# def main():
#     detector = SquareDetector(6.2, 100, 25)
#     detector.add_pad(0, 1, 0)
#     detector.add_pad(0, -1, 0)
#     detector.add_pad(0, 1, 1)
#     detector.add_pad(0, 1, -1)
#     detector.add_pad(0, -1, 1)
#     detector.add_pad(0, -1, -1)
#     detector.plot_detector(global_coords=False)
#     detector.plot_detector(global_coords=True)
#     print('bonzo')

# ---------------------------------------------------------------
# Build LARGE pads (8 columns × 4 rows)
# Using:
#   pad size = 12.4 mm (half-width = 6.2 mm)
#   pitch     = 12.5 mm
# ---------------------------------------------------------------

def build_large_pads():
    half_width = 6.2
    pitch = 12.5

    det = SquareDetector()

    # First pad at (0, 0)
    det.add_pad_absolute(0, 0, half_width)

    # Now place the others
    for c in range(8):
        for r in range(4):
            if r == 0 and c == 0:
                continue  # pad already added

            x = c * pitch
            y = r * pitch

            det.add_pad_absolute(x, y, half_width)

    return det


# ---------------------------------------------------------------
# Build SMALL pads (10 columns × 5 rows)
# Using:
#   pad size = 9.9 mm (half-width = 4.95 mm)
#   pitch     = 10.0 mm
# ---------------------------------------------------------------

def build_small_pads(y_offset):
    half_width = 4.95
    pitch = 10.0

    det = SquareDetector()

    # First pad
    det.add_pad_absolute(0, y_offset, half_width)

    for c in range(10):
        for r in range(5):
            if r == 0 and c == 0:
                continue

            x = c * pitch
            y = y_offset + r * pitch

            det.add_pad_absolute(x, y, half_width)

    return det


# ---------------------------------------------------------------
# WRITE PAD MAPPING
# ---------------------------------------------------------------

def write_mapping_to_txt(large, small, filename="p2_small_pad_mapping.txt"):
    with open(filename, "w") as f:
        f.write("PadID\tX(mm)\tY(mm)\tSize(mm)\tType\n")

        # Large pads
        for i, pad in enumerate(large.square_pads):
            f.write(
                f"{i}\t{pad.x:.2f}\t{pad.y:.2f}\t{2*pad.half_width:.2f}\tlarge\n"
            )

        # Small pads (continue numbering)
        offset = len(large.square_pads)

        for i, pad in enumerate(small.square_pads):
            f.write(
                f"{offset + i}\t{pad.x:.2f}\t{pad.y:.2f}\t{2*pad.half_width:.2f}\tsmall\n"
            )

    print(f"\n✔ Pad mapping written to {filename}\n")


# ---------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------

def main():
    large = build_large_pads()

    # small pads begin after 4 large-pad rows
    y_offset = 4 * 12.5

    small = build_small_pads(y_offset)

    # Plot
    large.plot_detector(global_coords=False)
    small.plot_detector(global_coords=False)

    # Export mapping
    write_mapping_to_txt(large, small, "p2_small_pad_mapping.txt")

    plt.show()


if __name__ == '__main__':
    main()

    print("bonzo")


#3 15.400 mm pad width
#interspace 0.100 mm
#16channels with large pads(18 pins not connected)
#4 13.900 mm pad width small
#interspace 0.100 mm
#25channels with small pads (first 8pins not connected)
#10x10 total active area
