#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 20/11/2025 12:21
Created in PyCharm
Created as plot_memory_log.py

@author: akallits
"""

import re
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

def read_memory_log(log_path):
    """
    Reads a system monitor log and returns a DataFrame with:
      - timestamp
      - mem_available (GB)
      - mem_total (GB)
      - mem_used (GB)
    """
    pattern = re.compile(
        r"^(.*?), Mem:\s*([0-9.]+)/([0-9.]+)\s*GB available"
    )

    timestamps = []
    mem_available = []
    mem_total = []

    with open(log_path, "r") as f:
        for line in f:
            match = pattern.search(line)
            if match:
                ts_str, avail_str, total_str = match.groups()
                timestamps.append(datetime.strptime(ts_str, "%Y-%m-%d %H:%M:%S"))
                mem_available.append(float(avail_str))
                mem_total.append(float(total_str))

    df = pd.DataFrame({
        "timestamp": timestamps,
        "mem_available": mem_available,
        "mem_total": mem_total,
    })

    df["mem_used"] = df["mem_total"] - df["mem_available"]
    return df


def plot_memory(df):
    plt.figure(figsize=(12, 5))
    plt.plot(df["timestamp"], df["mem_used"])
    plt.xlabel("Time")
    plt.ylabel("Memory Used (GB)")
    plt.title("Memory Usage Over Time")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    log_file = ("/home/akallits/Documents/Saclay-PostDoc/SPS_beam_test/system_monitor_daq/system_usage.log")   # ← Change this to your log file path
    df = read_memory_log(log_file)
    plot_memory(df)
