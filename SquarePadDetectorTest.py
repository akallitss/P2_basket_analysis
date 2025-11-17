#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 15/11/2025 19:06
Created in PyCharm
Created as SquarePadDetectorTest.py

@author: akallits
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from SquarePadDetector import SquareDetector

# ---------------------------------------------------------------
# Build full detector: large + small pads in same object
# ---------------------------------------------------------------

def build_detector():

    det = SquareDetector()   # SINGLE detector

    # -----------------------------
    # Large pad block
    # -----------------------------
    large_half = 6.2
    large_pitch = 12.5

    for c in range(8):
        for r in range(4):
            x = c * large_pitch + large_half
            y = r * large_pitch + large_half
            det.add_pad_absolute(x, y, large_half)

    # -----------------------------
    # Small pad block (offset below large)
    # -----------------------------
    small_half = 4.95
    small_pitch = 10.0

    y_offset = 4 * large_pitch     # same as before

    for r in range(5):
        for c in range(10):
            x = c * small_pitch + small_half
            y = y_offset + r * small_pitch + small_half
            det.add_pad_absolute(x, y, small_half)

    return det

# ---------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------

def main():

    detector = build_detector()

    detector.plot_detector(global_coords=False)

    detector.write_mapping("p2_small_detector_map.csv")

    plt.show()


if __name__ == "__main__":
    main()

    print("bonzo")


#3 15.400 mm pad width
#interspace 0.100 mm
#16channels with large pads(18 pins not connected)
#4 13.900 mm pad width small
#interspace 0.100 mm
#25channels with small pads (first 8pins not connected)
#10x10 total active area
