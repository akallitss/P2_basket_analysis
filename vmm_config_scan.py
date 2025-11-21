#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 19/11/2025 17:37
Created in PyCharm
Created as vmm_config_scan.py

@author: akallits
"""
import os
import uproot
import pandas as pd
import matplotlib.pyplot as plt
from  vmm_mapping import vmm_mapping


def load_hits_root(filename, branches=None, tree_name="hits"):
    """
    Load a ROOT file containing a hits tree and return a pandas DataFrame.

    Parameters
    ----------
    filename : str
        Path to the ROOT file.
    branches : list of str, optional
        List of branch names to read. If None, all branches are loaded.
        Branch names should be given without the tree prefix (e.g. "adc", not "hits/adc").
    tree_name : str
        Name of the tree inside the ROOT file.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing the selected branches.
    """

    # Open the ROOT file
    file = uproot.open(filename)

    if tree_name not in file:
        raise ValueError(f"Tree '{tree_name}' not found in ROOT file: {filename}")

    tree = file[tree_name]

    # If user gives branch names WITHOUT prefix — clean way
    if branches is not None:
        formatted = []
        for b in branches:
            if "/" not in b:  # user passed "adc"
                formatted.append(f"{tree_name}/{b}")
            else:  # rare case user passed full path
                formatted.append(b)
        branches = formatted

    df = tree.arrays(branches, library="pd")

    # Remove "hits/" prefix → make column names clean
    df.columns = [col.replace(f"{tree_name}/", "") for col in df.columns]

    return df



def main():

    df_run_scan = pd.read_csv("vmm_config_scan .csv")
    print(df_run_scan)

    #select the first row and print it
    first_row = df_run_scan.iloc[0]
    print(first_row)
    #print the values of each column in the first row
    for col in df_run_scan.columns:
        print(f"{col}: {first_row[col]}")

    #access a specific value
    gain = first_row["sg"]
    print(f"Gain: {gain}")

    #find all run numbers where gain is 3.0
    runs_with_gain_3 = df_run_scan[df_run_scan["sg"] == 3.0]["run_no"]
    print("Runs with gain 3.0:")
    print(runs_with_gain_3.tolist())

    #filter runs with sng = 0.0
    runs_with_sng_0 = df_run_scan[(df_run_scan["sng"] == 0.0) & (df_run_scan["sg"] == 3.0)]["run_no"]
    print("Runs with sng 0.0 and gain 3.0:")
    print(runs_with_sng_0.tolist())
    # input()

    #filter runs with sng = 1.0
    runs_with_sng_1 = df_run_scan[(df_run_scan["sng"] == 1.0) & (df_run_scan["sg"] == 3.0)]["run_no"]
    print("Runs with sng 1.0 and gain 3.0:")
    print(runs_with_sng_1.tolist())

    #path to run files
    data_dir = "/media/akallits/EXTERNAL_USB/2312292/Extras/Physics/Post-Doc-Saclay/SPS_beam_test/data/VMM-alinx_data/5kHz-muons-config-scan"

    #check if the run numbers exist as directories in the data directory
    for run_no in df_run_scan["run_no"]:
        run_dir = os.path.join(data_dir, f"run_{run_no}")
        if os.path.isdir(run_dir):
            print(f"Run directory exists: {run_dir}")
        else:
            print(f"Run directory does not exist: {run_dir}")

    #get the first five enp*root files from runs with gain 3.0
    # for run_no in runs_with_gain_3:
    for run_no in runs_with_sng_0:
        run_dir = os.path.join(data_dir, f"run_{run_no}")
        if os.path.isdir(run_dir):
            root_files = [f for f in os.listdir(run_dir) if f.startswith("enp") and f.endswith(".root")]
            first_five_files = root_files[:2]
            print(f"First five enp*root files in run {run_no}:")
            for file in first_five_files:
                print(file)
        else:
            print(f"Run directory does not exist: {run_dir}")



    #load hits from the first file of the runs with gain 3.0
    # for run_no in runs_with_sng_0:
    for run_no in runs_with_sng_1:
    # for run_no in runs_with_gain_3:
        run_dir = os.path.join(data_dir, f"run_{run_no}")
        if os.path.isdir(run_dir):
            root_files = [f for f in os.listdir(run_dir) if f.startswith("enp") and f.endswith(".root")]
            if root_files:
                first_file = root_files[1]
                file_path = os.path.join(run_dir, first_file)
                print(f"Loading hits from file: {file_path}")
                df_hits = load_hits_root(file_path, branches=["adc", "vmm", "ch", "time"])
                print(df_hits.head())
            else:
                print(f"No enp*root files found in run {run_no}")
        else:
            print(f"Run directory does not exist: {run_dir}")


    #for each run with gain 3.0, print the vmm ids used
    # for run_no in runs_with_gain_3:
    for run_no in runs_with_sng_1:
    # for run_no in runs_with_sng_0:
        print(f"Run {run_no} VMM IDs:")
        for config_name, config in vmm_mapping.items():
            print(f"  Configuration: {config_name} - VMM IDs: {config['vmm_ids']}")

    #for each run with gain 3.0, plot a histogram of adc values for each vmm id
    fig, ax = plt.subplots()
    # for run_no in runs_with_sng_0:
    for run_no in runs_with_sng_1:
    # for run_no in runs_with_gain_3:
        run_dir = os.path.join(data_dir, f"run_{run_no}")
        if os.path.isdir(run_dir):
            root_files = [f for f in os.listdir(run_dir) if f.startswith("enp") and f.endswith(".root")]
            if root_files:
                first_file = root_files[1]
                file_path = os.path.join(run_dir, first_file)
                df_hits = load_hits_root(file_path, branches=["hits/adc", "hits/vmm", "hits/ch", "hits/time"])
                for vmm_id in df_hits["vmm"].unique():
                    adc_values = df_hits[df_hits["vmm"] == vmm_id]["adc"]
                    ax.hist(adc_values, bins=50, alpha=0.5, label=f"Run {run_no} VMM {vmm_id}")
    ax.set_xlabel("ADC Value")
    ax.set_ylabel("Counts")
    ax.set_title("ADC Value Distribution for Runs with Gain 3.0")
    ax.legend()
    # plt.show()

    #get all unique vmm ids from vmm mapping
    vmm_ids = set()
    for config in vmm_mapping.values():
        vmm_ids.update(config["vmm_ids"])
    vmm_ids = sorted(list(vmm_ids))
    print(f"Unique VMM IDs: {vmm_ids}")
    # input()

    #plot on seperate figures for each vmm id the adc value distribution from runs with gain 3.0
    for vmm_id in vmm_ids:
        # For vmm_id, get corresponding config name
        name = None
        for config_name, config in vmm_mapping.items():
            if vmm_id in config["vmm_ids"]:
                name = config_name
                break
        print(f"Plotting ADC distribution for VMM {vmm_id} ({name})")

        fig, ax = plt.subplots()
        # for run_no in runs_with_gain_3:
        # for run_no in runs_with_sng_0:
        for run_no in runs_with_sng_1:
            run_dir = os.path.join(data_dir, f"run_{run_no}")

            # Get parameter string for run for legend
            param_row = df_run_scan[df_run_scan["run_no"] == run_no]
            param_str = ", ".join([f"{col}={param_row.iloc[0][col]}" for col in df_run_scan.columns if col != "run_no"])

            if os.path.isdir(run_dir):
                root_files = [f for f in os.listdir(run_dir) if f.startswith("enp") and f.endswith(".root")]
                if root_files:
                    first_file = root_files[1]
                    file_path = os.path.join(run_dir, first_file)
                    df_hits = load_hits_root(file_path, branches=["hits/adc", "hits/vmm", "hits/ch", "hits/time"])
                    adc_values = df_hits[df_hits["vmm"] == vmm_id]["adc"]
                    # ax.hist(adc_values, bins=50, histtype='step', linewidth=1.5, label=f"Run {run_no}", density=True)
                    ax.hist(adc_values, bins=50, histtype='step', linewidth=1.5, label=param_str, density=True)
        ax.set_xlabel("ADC Value")
        ax.set_ylabel("Counts")
        ax.set_title(f"ADC Value Distribution for {name} VMM {vmm_id} (Gain 3.0 Runs)")
        ax.legend()
    # plt.show()


    #set a cut on adc < 100 and save the mean adc value and the rms for each vmm id and each run with gain 3.0 and sng 1.0
    adc_cut = 100
    results = []
    for vmm_id in vmm_ids:
        for run_no in runs_with_sng_1:
            run_dir = os.path.join(data_dir, f"run_{run_no}")
            if os.path.isdir(run_dir):
                root_files = [f for f in os.listdir(run_dir) if f.startswith("enp") and f.endswith(".root")]
                if root_files:
                    first_file = root_files[1]
                    file_path = os.path.join(run_dir, first_file)
                    df_hits = load_hits_root(file_path, branches=["hits/adc", "hits/vmm", "hits/ch", "hits/time"])
                    adc_values = df_hits[(df_hits["vmm"]==vmm_id) & (df_hits["adc"] < adc_cut)]["adc"]
                    mean_adc = adc_values.mean()
                    rms_adc = adc_values.std()
                    results.append({
                        "run_no": run_no,
                        "vmm_id": vmm_id,
                        "peaking_time": df_run_scan[df_run_scan["run_no"] == run_no].iloc[0]["snt"],
                        "mean_adc": mean_adc,
                        "rms_adc": rms_adc
                    })
    df_results = pd.DataFrame(results)
    print("Mean and RMS ADC values for each VMM ID and run (ADC < 100):")
    print(df_results)
    df_results.to_csv("vmm_adc_analysis.csv", index=False)

    #filter the nan values from the dataframe
    df_results = df_results.dropna()
    print(df_results)

    #plot mean adc vs peaking time for each vmm id
    for vmm_id in vmm_ids:
        df_vmm = df_results[df_results["vmm_id"] == vmm_id]
        if df_vmm.empty:
            continue
        plt.errorbar(df_vmm["peaking_time"], df_vmm["mean_adc"], yerr=df_vmm["rms_adc"], fmt='o-', label=f"VMM {vmm_id}")
    plt.xlabel("Peaking Time (snt)")
    plt.ylabel("Mean ADC Value (ADC < 100)")
    plt.title("Mean ADC vs Peaking Time for Each VMM ID (Gain 3.0, sng 1.0)")
    plt.legend()
    # plt.savefig("mean_adc_vs_peaking_time.png")
    plt.show()

    pass




if __name__ == "__main__":
    main()

    print("bonzo")
