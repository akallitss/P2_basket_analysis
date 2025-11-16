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

def build_large_pads():
    large = SquareDetector(6.2)   # half width 12.4 / 2
    pitch = 12.5                  # 12.4 + 0.1

    # First pad stays at (0,0)
    idx = 0

    # 4 rows × 8 columns
    for c in range(8):
        for r in range(4):
            if r == 0 and c == 0:
                continue  # already placed first pad
            large.add_pad(0, (pitch * c) / (2*large.pad_half_width),
                             (pitch * r) / (2*large.pad_half_width))

    return large


def build_small_pads(y_offset):
    small = SquareDetector(4.95)  # half width 9.9 / 2
    pitch = 10.0                  # 9.9 + 0.1

    # First pad in small region
    small.square_pads[0].y = y_offset

    # 5 rows × 10 columns
    for r in range(5):
        for c in range(10):
            if r == 0 and c == 0:
                continue
            small.add_pad(0, (pitch * c) / (2*small.pad_half_width),
                             (pitch * r) / (2*small.pad_half_width))

    return small

# ---------------------------------------------------------------
# NEW: Create mapping text file directly from your pad objects
# ---------------------------------------------------------------
def write_mapping_to_txt(large, small, filename="p2_small_pad_mapping.txt"):
    with open(filename, "w") as f:
        f.write("PadID\tX(mm)\tY(mm)\tSize(mm)\tType\n")

        # Large pads
        for i, pad in enumerate(large.square_pads):
            f.write(
                f"{i}\t{pad.x:.2f}\t{pad.y:.2f}\t{2*pad.half_width:.2f}\tlarge\n"
            )

        # Small pads (ID continues from large)
        offset = len(large.square_pads)
        for i, pad in enumerate(small.square_pads):
            f.write(
                f"{offset + i}\t{pad.x:.2f}\t{pad.y:.2f}\t{2*pad.half_width:.2f}\tsmall\n"
            )

    print(f"\n Pad mapping written to {filename}\n")


# ---------------------------------------------------------------


def main():
    # Build large region
    large = build_large_pads()

    # Compute where small pads should begin in y:
    y_offset = 4 * 12.5   # 4 rows * pitch

    small = build_small_pads(y_offset)

    large.plot_detector(global_coords=False)
    small.plot_detector(global_coords=False)

    write_mapping_to_txt(large, small, "p2_small_pad_mapping.txt")
    plt.show()

    print("bonzo")


if __name__ == '__main__':
    main()


#3 15.400 mm pad width
#interspace 0.100 mm
#16channels with large pads(18 pins not connected)
#4 13.900 mm pad width small
#interspace 0.100 mm
#25channels with small pads (first 8pins not connected)
#10x10 total active area
