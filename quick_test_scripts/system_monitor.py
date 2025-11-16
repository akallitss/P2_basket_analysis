#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 15/11/2025 13:43
Created in PyCharm
Created as system_monitor.py

@author: akallits
"""


import psutil
import time
from datetime import datetime

#128.141.41.210:/mnt/data /mnt/data nfs defaults 0 0 #command on the /etc/fstab file to mount the remote directory
LOG_FILE = "/mnt/data/beam_sps_25/P2_logs/system_usage.log" #link between two computers the /mnt/data directories
# LOG_FILE = "system_usage.log"
INTERVAL = 2  # seconds between logs

def log_system_status():
    while True:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Memory
        mem = psutil.virtual_memory()
        mem_available_gb = mem.available / (1024**3)
        mem_total_gb = mem.total / (1024**3)

        # Disk storage ("/" for Linux root partition)
        disk = psutil.disk_usage("/")
        disk_free_gb = disk.free / (1024**3)
        disk_total_gb = disk.total / (1024**3)

        # Write to file
        with open(LOG_FILE, "a") as f:
            f.write(
                f"{timestamp}, "
                f"Mem: {mem_available_gb:.2f}/{mem_total_gb:.2f} GB available, "
                f"Disk: {disk_free_gb:.2f}/{disk_total_gb:.2f} GB free\n"
            )

        time.sleep(INTERVAL)


if __name__ == "__main__":
    print(f"Logging to {LOG_FILE} every {INTERVAL}s... (Ctrl+C to stop)")
    log_system_status()
