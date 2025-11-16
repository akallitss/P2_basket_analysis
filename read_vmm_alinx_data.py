#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 14/11/2025 19:32
Created in PyCharm
Created as read_vmm_alinx_data.py

@author: akallits
"""


import os
import uproot
from matplotlib.colors import LogNorm
import awkward as ak
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

run_number = 'run_67'
data_dir = f"/home/akallits/Documents/Saclay-PostDoc/SPS_beam_test/data/VMM-alinx_data/{run_number}"

files = os.listdir(data_dir)
print(run_number)

enp_files = sorted([
    os.path.join(data_dir, f)
    for f in os.listdir(data_dir)
    if f.startswith("enp") and f.endswith(".root")
])

print("Selected files:")
for f in enp_files:
    print("  ", f)

file = uproot.open(enp_files[0])
tree = file["hits"]
print(tree.keys())
print(tree.show())
# input()
branches = [
    "hits/id", "hits/det", "hits/plane", "hits/fec", "hits/vmm",
    "hits/readout_time", "hits/time", "hits/geo_id", "hits/ch",
    "hits/pos", "hits/bcid", "hits/tdc", "hits/adc",
    "hits/over_threshold", "hits/chip_time", "hits/event_counter"
]


df_hits = tree.arrays(library="pd")
print(df_hits.head())

#0,7,6,5,4 hybrids
#1,2,13,14,11,12,10,11,9,8 vmms
#trigger, large p2 card1, large p2 card2, small p2-1, small p2-3

#make a dictionary of the detectors under test and the vmm ids associated with them
vmm_hybrid_mapping = {
    'trigger': [0, 1],
    'p2_large_1': [12, 13, 14, 15],
    'p2_small_1': [10, 11],
    'p2_small_3': [8, 9],
}

detector = 'p2_small_1'  # Change this to analyze a different detector
vmm_ids = vmm_hybrid_mapping[detector]
print(f'Detector: {detector} is associated with VMM IDs: {vmm_ids}')

# print(df_hits[df_hits['vmm'].isin(vmm_ids)].head())
# input()


####################### PLOTTING #########################

#plot adc distribution from df_hits
plt.figure(figsize=(10, 6))
plt.hist(df_hits['adc'], bins=200, range=(0, 4096), alpha=0.7)
plt.xlim(0, 1200)
plt.title('ADC Distribution')
plt.xlabel('ADC Value')
plt.ylabel('Counts')
plt.grid(True)
# plt.show()


#plot vmm distribution from df_hits
plt.figure(figsize=(10, 6))
plt.hist(df_hits['vmm'], bins=np.arange(df_hits['vmm'].min(), df_hits['vmm'].max() + 2) - 0.5, alpha=0.7)
plt.title('VMM Distribution')
plt.xlabel('VMM Value')
plt.ylabel('Counts')
plt.grid(True)


# Plot channel vs adc
plt.figure(figsize=(12, 6))
plt.scatter(df_hits['ch'], df_hits['adc'], alpha=0.5, s=1)
plt.title('Channel vs ADC')
plt.xlabel('Channel')
plt.ylabel('ADC Value')
plt.grid(True)
# plt.show()

#plot 2D histogram of channel vs adc
# vmms = df_hits['vmm'].unique()
# for vmm in vmms:
#     df_vmm = df_hits[df_hits['vmm'] == vmm]
#     plt.figure(figsize=(12, 6))
#     plt.hist2d(df_vmm['ch'], df_vmm['adc'], bins=[50, 200], range=[[0, 65], [0, 1200]], cmap='viridis', cmin=1, norm=LogNorm())
#     plt.title(f'2D Histogram of Channel vs ADC for VMM {vmm}')
#     plt.xlabel('Channel')
#     plt.ylabel('ADC Value')
#     plt.grid(True)

for vmm in vmm_ids:
    df_vmm_det = df_hits[df_hits['vmm'] == vmm]
    plt.figure(figsize=(12, 6))
    plt.hist2d(df_vmm_det['ch'], df_vmm_det['adc'], bins=[50, 200], range=[[0, 65], [0, 1200]], cmap='viridis', cmin=1, norm=LogNorm())
    plt.title(f'2D Histogram of Channel vs ADC for VMM {vmm}')
    plt.xlabel('Channel')
    plt.ylabel('ADC Value')
    plt.grid(True)

plt.show()

