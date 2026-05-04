#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vmm_combine_results.py

Combines direct-trigger and lower-rate VMM config scan results.

Strategy
--------
Each source contributes only the configurations that were actually measured
in that beam condition.  Rows outside the declared config lists are dropped
before merging so stale or mis-assigned configs never pollute the output.

A 'source' column is added to every output so every row is traceable:
    'direct_trigger'  — from direct-trigger data
    'lower_rate'      — from lower-rate data

Outputs (written to OUT_DIR)
-------
    vmm_snr_results_combined.csv      — VMM-level SNR
    vmm_snr_per_channel_combined.csv  — per-channel SNR
    vmm_snr_summary_combined.csv      — best-config summary per VMM
    snr_heatmap_combined.pdf/.png     — SNR heatmap with source annotation
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# ─────────────────────────────────────────────
# PATHS  — edit if files move
# ─────────────────────────────────────────────

# Direct-trigger results (primary)
DT_DIR = "/drf/projets/clas12/P2/akallits/"
DT_SNR    = DT_DIR + "vmm_snr_results_15kHz.csv"
DT_SNR_CH = DT_DIR + "vmm_snr_per_channel_15kHz.csv"

# Lower-rate results (fill-in)
LR_DIR = "/drf/projets/clas12/P2/akallits/"
LR_SNR    = LR_DIR + "vmm_snr_results_5kHz.csv"
LR_SNR_CH = LR_DIR + "vmm_snr_per_channel_5kHz.csv"

# lab results (fill-in)
LAB_DIR    = "/drf/projets/clas12/P2/akallits/"
LAB_SNR    = LAB_DIR + "vmm_snr_results_lab.csv"
LAB_SNR_CH = LAB_DIR + "vmm_snr_per_channel_lab.csv"


# Output directory
OUT_DIR = DT_DIR

# ─────────────────────────────────────────────
# CONFIG ASSIGNMENT
# Each tuple is (sg, snt).  Only these configs are kept from each source.
# ─────────────────────────────────────────────

DT_CONFIGS = [
    (3.0,  50.0),
    (3.0, 100.0),
    (3.0, 200.0),
    (4.5, 200.0),
    (6.0, 200.0),
]

LR_CONFIGS = [
    (4.5,  50.0),
    (4.5, 100.0),
]

LAB_CONFIGS = [
    (1.0,  50.0),
    (1.0, 100.0),
    (1.0, 200.0),
]

# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────

def _filter_configs(df, configs, label):
    """Keep only rows whose (sg, snt) is in configs; warn about missing ones."""
    mask = pd.Series(
        [tuple(r) in configs for r in df[["sg", "snt"]].values.tolist()],
        index=df.index
    )
    filtered = df[mask].copy()
    present = set(map(tuple, filtered[["sg", "snt"]].drop_duplicates().values.tolist()))
    for cfg in configs:
        if cfg not in present:
            print(f"  WARNING [{label}]: config sg={cfg[0]:.1f}/snt={cfg[1]:.0f} "
                  f"not found in source — will be absent from combined output")
    return filtered


def _combine(*dfs):
    """
    Concatenate any number of source-tagged DataFrames.
    Each has already been filtered to its config list so there is no overlap.
    """
    combined = pd.concat(list(dfs), ignore_index=True)
    combined = combined.sort_values(["sg", "snt", "vmm_id"]).reset_index(drop=True)
    return combined


def _summarise(df_snr, df_snr_ch):
    """
    Recompute best-configuration summary from combined VMM-level SNR
    and combined channel-level SNR.
    """
    records = []

    for vmm_id in sorted(df_snr["vmm_id"].unique()):
        # VMM-level best: highest SNR among quality='ok' rows
        df_v = df_snr[
            (df_snr["vmm_id"]        == vmm_id) &
            (df_snr["noise_quality"] == "ok")
        ]
        if df_v.empty:
            continue
        best_v  = df_v.loc[df_v["snr"].idxmax()]
        vmm_cfg = f"sg={best_v['sg']:.1f}/snt={best_v['snt']:.0f}"
        vmm_snr = best_v["snr"]
        vmm_src = best_v["source"]

        # Channel-level best: config with highest median channel SNR
        df_c = df_snr_ch[df_snr_ch["vmm_id"] == vmm_id]
        if df_c.empty:
            ch_cfg = ""
            ch_snr = np.nan
            ch_src = ""
        else:
            ch_med     = df_c.groupby(["sg", "snt"])["snr_ch"].median()
            best_c_idx = ch_med.idxmax()
            ch_snr     = ch_med.max()
            ch_cfg     = (f"sg={best_c_idx[0]:.1f}"
                          f"/snt={best_c_idx[1]:.0f}")
            ch_src = df_snr_ch.loc[
                (df_snr_ch["vmm_id"] == vmm_id) &
                (df_snr_ch["sg"]     == best_c_idx[0]) &
                (df_snr_ch["snt"]    == best_c_idx[1]),
                "source"
            ].iloc[0]

        records.append({
            "vmm_id"         : vmm_id,
            "vmm_best_config": vmm_cfg,
            "vmm_snr"        : round(vmm_snr, 2),
            "vmm_source"     : vmm_src,
            "ch_best_config" : ch_cfg,
            "ch_median_snr"  : round(ch_snr, 2),
            "ch_source"      : ch_src,
            "agreement"      : vmm_cfg == ch_cfg,
        })

    return pd.DataFrame(records)


def _plot_snr_heatmap(df_snr, out_dir):
    """
    SNR heatmap: rows = VMM IDs, columns = (sg, snt) configs.
    Cell colour = SNR value; hatching = noise_quality != 'ok';
    border colour = source (blue=DT, orange=LR).
    Saved as PDF and PNG.
    """
    vmm_ids = sorted(df_snr["vmm_id"].unique())
    configs  = sorted(df_snr[["sg", "snt"]].drop_duplicates().values.tolist())
    col_labels = [f"sg={c[0]:.1f}\nsnt={c[1]:.0f}" for c in configs]

    snr_mat  = np.full((len(vmm_ids), len(configs)), np.nan)
    qual_mat = np.full((len(vmm_ids), len(configs)), "", dtype=object)
    src_mat  = np.full((len(vmm_ids), len(configs)), "", dtype=object)

    for i, vmm in enumerate(vmm_ids):
        for j, (sg, snt) in enumerate(configs):
            row = df_snr[
                (df_snr["vmm_id"] == vmm) &
                (df_snr["sg"]     == sg)  &
                (df_snr["snt"]    == snt)
            ]
            if not row.empty:
                snr_mat[i, j]  = row["snr"].values[0]
                qual_mat[i, j] = row["noise_quality"].values[0]
                src_mat[i, j]  = row["source"].values[0]

    vmin = np.nanmin(snr_mat)
    vmax = np.nanmax(snr_mat)
    cmap = plt.cm.viridis
    norm = Normalize(vmin=vmin, vmax=vmax)

    fig, ax = plt.subplots(figsize=(max(8, len(configs) * 1.4), max(5, len(vmm_ids) * 0.7)))

    for i in range(len(vmm_ids)):
        for j in range(len(configs)):
            val = snr_mat[i, j]
            if np.isnan(val):
                color = "#e0e0e0"
                txt   = "—"
            else:
                color = cmap(norm(val))
                txt   = f"{val:.1f}"

            # cell background
            ax.add_patch(plt.Rectangle(
                (j - 0.5, i - 0.5), 1, 1,
                facecolor=color, edgecolor="none"
            ))

            # source border
            src = src_mat[i, j]
            if src == "direct_trigger":
                border_color = "#1f77b4"
            elif src == "lower_rate":
                border_color = "#ff7f0e"
            elif src == "lab":
                border_color = "#2ca02c"
            else:
                border_color = "none"
            ax.add_patch(plt.Rectangle(
                (j - 0.5, i - 0.5), 1, 1,
                facecolor="none",
                edgecolor=border_color, linewidth=2.5
            ))

            # hatch for non-ok quality
            if qual_mat[i, j] not in ("ok", ""):
                ax.add_patch(plt.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    facecolor="none",
                    edgecolor="red", linewidth=0, hatch="//", alpha=0.5
                ))

            # SNR text
            brightness = np.mean(cmap(norm(val))[:3]) if not np.isnan(val) else 0.9
            txt_color  = "white" if brightness < 0.55 else "black"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=9, color=txt_color, fontweight="bold")

    ax.set_xlim(-0.5, len(configs)  - 0.5)
    ax.set_ylim(-0.5, len(vmm_ids) - 0.5)
    ax.set_xticks(range(len(configs)))
    ax.set_xticklabels(col_labels, fontsize=9)
    ax.set_yticks(range(len(vmm_ids)))
    ax.set_yticklabels([f"VMM {v}" for v in vmm_ids], fontsize=9)
    ax.set_xlabel("Configuration (sg / snt)", fontsize=11)
    ax.set_ylabel("VMM ID", fontsize=11)
    ax.set_title("Combined SNR Heatmap\n"
                 "(border: blue=direct-trigger, orange=lower-rate, green=lab  |  hatch=quality≠ok)",
                 fontsize=11)

    # colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("SNR", fontsize=10)

    # legend patches
    legend_handles = [
        mpatches.Patch(edgecolor="#1f77b4", facecolor="none", linewidth=2.5,
                       label="direct trigger"),
        mpatches.Patch(edgecolor="#ff7f0e", facecolor="none", linewidth=2.5,
                       label="lower rate"),
        mpatches.Patch(edgecolor="#2ca02c", facecolor="none", linewidth=2.5,
                       label="lab"),
        mpatches.Patch(facecolor="none", edgecolor="red", hatch="//", alpha=0.5,
                       label="quality ≠ ok"),
    ]
    ax.legend(handles=legend_handles, loc="upper left",
              bbox_to_anchor=(1.12, 1.0), fontsize=9, frameon=True)

    fig.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out_dir, f"snr_heatmap_combined.{ext}"),
                    bbox_inches="tight")
    plt.close(fig)
    print(f"  snr_heatmap_combined.pdf/png  →  {out_dir}")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    # ── Load ──────────────────────────────────────────────────────
    dt_snr    = pd.read_csv(DT_SNR)
    dt_snr_ch = pd.read_csv(DT_SNR_CH)
    lr_snr    = pd.read_csv(LR_SNR)
    lr_snr_ch = pd.read_csv(LR_SNR_CH)
    lab_snr = pd.read_csv(LAB_SNR)
    lab_snr_ch = pd.read_csv(LAB_SNR_CH)

    # ── Filter to declared configs ────────────────────────────────
    print("Filtering direct-trigger to declared DT configs …")
    dt_snr    = _filter_configs(dt_snr,    set(DT_CONFIGS), "DT VMM-level")
    dt_snr_ch = _filter_configs(dt_snr_ch, set(DT_CONFIGS), "DT channel-level")

    print("Filtering lower-rate to declared LR configs …")
    lr_snr    = _filter_configs(lr_snr,    set(LR_CONFIGS), "LR VMM-level")
    lr_snr_ch = _filter_configs(lr_snr_ch, set(LR_CONFIGS), "LR channel-level")

    print("Filtering lower-rate to declared LAB configs …")
    lab_snr = _filter_configs(lab_snr, set(LAB_CONFIGS), "LAB VMM-level")
    lab_snr_ch = _filter_configs(lab_snr_ch, set(LAB_CONFIGS), "LAB channel-level")

    dt_snr["source"]    = "direct_trigger"
    dt_snr_ch["source"] = "direct_trigger"
    lr_snr["source"]    = "lower_rate"
    lr_snr_ch["source"] = "lower_rate"

    lab_snr["source"] = "lab"
    lab_snr_ch["source"] = "lab"

    print(f"\nDirect trigger — VMM-level    : {len(dt_snr):>5} rows  "
          f"| configs: {sorted(dt_snr[['sg','snt']].drop_duplicates().values.tolist())}")
    print(f"Direct trigger — channel-level: {len(dt_snr_ch):>5} rows")
    print(f"Lower rate     — VMM-level    : {len(lr_snr):>5} rows  "
          f"| configs: {sorted(lr_snr[['sg','snt']].drop_duplicates().values.tolist())}")
    print(f"Lower rate     — channel-level: {len(lr_snr_ch):>5} rows")

    print(f"Lab            — VMM-level    : {len(lab_snr):>5} rows  "
          f"| configs: {sorted(lab_snr[['sg','snt']].drop_duplicates().values.tolist())}")
    print(f"Lab            — channel-level: {len(lab_snr_ch):>5} rows")

    # ── Combine ───────────────────────────────────────────────────
    combined_snr    = _combine(dt_snr,    lr_snr,   lab_snr)
    combined_snr_ch = _combine(dt_snr_ch, lr_snr_ch, lab_snr_ch)

    n_dt  = (combined_snr["source"] == "direct_trigger").sum()
    n_lr  = (combined_snr["source"] == "lower_rate").sum()
    n_lab = (combined_snr["source"] == "lab").sum()
    print(f"\nCombined VMM-level SNR   : {len(combined_snr)} rows "
          f"({n_dt} direct_trigger  +  {n_lr} lower_rate  +  {n_lab} lab)")

    n_dt_ch  = (combined_snr_ch["source"] == "direct_trigger").sum()
    n_lr_ch  = (combined_snr_ch["source"] == "lower_rate").sum()
    n_lab_ch = (combined_snr_ch["source"] == "lab").sum()
    print(f"Combined channel-level SNR: {len(combined_snr_ch)} rows "
          f"({n_dt_ch} direct_trigger  +  {n_lr_ch} lower_rate  +  {n_lab_ch} lab)")

    # ── Summary ───────────────────────────────────────────────────
    summary = _summarise(combined_snr, combined_snr_ch)

    print(f"\n{'='*70}")
    print("COMBINED BEST-CONFIGURATION SUMMARY")
    print(f"{'='*70}")
    print(summary.to_string(index=False))

    votes = summary["vmm_best_config"].value_counts()
    print(f"\nVMM-level votes:")
    for cfg, n in votes.items():
        print(f"  {cfg} → {n} VMM(s)")
    print(f"\nRecommended: {votes.idxmax()}  "
          f"({votes.max()} of {len(summary)} VMMs)")

    n_agree = summary["agreement"].sum()
    print(f"VMM and channel level agree for {n_agree}/{len(summary)} VMMs")

    # ── Save CSVs ─────────────────────────────────────────────────
    os.makedirs(OUT_DIR, exist_ok=True)

    out_snr    = OUT_DIR + "vmm_snr_results_combined.csv"
    out_snr_ch = OUT_DIR + "vmm_snr_per_channel_combined.csv"
    out_sum    = OUT_DIR + "vmm_snr_summary_combined.csv"

    combined_snr.to_csv(out_snr,       index=False)
    combined_snr_ch.to_csv(out_snr_ch, index=False)
    summary.to_csv(out_sum,            index=False)

    print(f"\nSaved CSVs:")
    print(f"  {out_snr}")
    print(f"  {out_snr_ch}")
    print(f"  {out_sum}")

    # ── Plot ──────────────────────────────────────────────────────
    _plot_snr_heatmap(combined_snr, OUT_DIR)


if __name__ == "__main__":
    main()
