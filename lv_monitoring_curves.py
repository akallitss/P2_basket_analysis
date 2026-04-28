#!/usr/bin/env python3
"""
Low Voltage Monitoring Curves
==============================
Reads tti1, tti2, and enxa triplet log files from the lv_monitoring directory,
lets you pick one or more triplets, and saves voltage/current figures.

Layout (time-based overview)
------
  Col 0            Col 1
  ┌────────────────┬────────────────┐
  │ TTI1 Voltage   │ TTI2 Voltage   │  row 0
  ├────────────────┼────────────────┤
  │ TTI1 Current   │ TTI2 Current   │  row 1
  ├────────────────┴────────────────┤
  │    Backend Connectivity (ENXA)  │  row 2
  └─────────────────────────────────┘

Run (interactive):   python lv_monitoring_curves.py
Run (headless):      python lv_monitoring_curves.py --list
                     python lv_monitoring_curves.py --triplet 3 5 7
                     python lv_monitoring_curves.py --date 20260420 20260421
"""
import argparse
import os
import re
import sys
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timedelta

import matplotlib
_headless = not os.environ.get("DISPLAY") and sys.platform != "darwin"
if _headless:
    matplotlib.use("Agg")   # avoids display probe on headless machines
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import pandas as pd

_SHOW = not _headless

# ── Configuration ─────────────────────────────────────────────────────────────

DATA_DIR = Path(
    "/local/home/ak271430/Documents/PostDocSaclay/data/TestPagure/lv_monitoring"
    # "/local/p2/testbench/TestBenchPagure/monitoring"
)

# One colour per channel (up to 3 channels per TTI unit)
CHANNEL_COLORS = ["tab:blue", "tab:orange", "tab:green"]

# ── Irradiation parameters ─────────────────────────────────────────────────────
IRRAD_START_HHMM     = "12:45"  # HH:MM on the first day of irradiation
DOSE_RATE_KRAD_PER_H = 1.0      # krad per hour


# =============================================================================
# STEP 1 — Discover files and group them into triplets
# =============================================================================

FILE_PATTERN = re.compile(r"^(tti1|tti2|enxa[^_]+)_(\d{8}-\d{6})_mon\.log$")


def find_triplets(data_dir: Path) -> dict:
    """
    Scan *data_dir* and return a dict keyed by timestamp string.
    Each value is {'tti1': Path, 'tti2': Path, 'enxa': Path}.
    Only complete triplets are included.
    """
    groups = defaultdict(dict)
    for f in sorted(data_dir.glob("*_mon.log")):
        m = FILE_PATTERN.match(f.name)
        if not m:
            continue
        prefix, ts = m.group(1), m.group(2)
        key = "enxa" if prefix.startswith("enxa") else prefix
        groups[ts][key] = f
    return {
        ts: files for ts, files in groups.items()
        if all(k in files for k in ("tti1", "tti2", "enxa"))
    }


def _print_triplet_menu(timestamps: list) -> None:
    print("\nAvailable triplets:")
    print(f"  {'#':>3}   {'Date':10}  {'Time':8}")
    print(f"  {'─'*3}   {'─'*10}  {'─'*8}")
    for i, ts in enumerate(timestamps, start=1):
        dt = datetime.strptime(ts, "%Y%m%d-%H%M%S")
        print(f"  [{i:>2}]  {dt.strftime('%Y-%m-%d')}  {dt.strftime('%H:%M:%S')}")


# =============================================================================
# STEP 2 — Triplet selection (interactive or CLI)
# =============================================================================

def select_triplets_interactive(triplets: dict) -> list:
    """
    Print a numbered menu and return a list of (timestamp, files) tuples.
    Accepts: a single number, comma-separated numbers (e.g. "3,5,7"), or "all".
    """
    timestamps = sorted(triplets.keys())
    _print_triplet_menu(timestamps)
    n = len(timestamps)

    while True:
        try:
            raw = input(f"\nSelect triplet(s) [1–{n}, comma-separated, or 'all']: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nAborted.")
            sys.exit(0)

        if raw.lower() == "all":
            indices = list(range(1, n + 1))
        else:
            parts = [p.strip() for p in raw.split(",")]
            try:
                indices = [int(p) for p in parts if p]
            except ValueError:
                print("  Please enter numbers, comma-separated, or 'all'.")
                continue

        if all(1 <= idx <= n for idx in indices):
            return [(timestamps[idx - 1], triplets[timestamps[idx - 1]]) for idx in indices]
        print(f"  All numbers must be between 1 and {n}.")


def select_triplets_by_number(triplets: dict, numbers: list) -> list:
    timestamps = sorted(triplets.keys())
    result = []
    for n in numbers:
        if not 1 <= n <= len(timestamps):
            print(f"ERROR: --triplet {n} out of range (1–{len(timestamps)})")
            sys.exit(1)
        ts = timestamps[n - 1]
        result.append((ts, triplets[ts]))
    return result


def select_triplets_by_date(triplets: dict, dates: list) -> list:
    """dates is a list of YYYYMMDD strings."""
    result = []
    for ts, files in sorted(triplets.items()):
        if any(ts.startswith(d) for d in dates):
            result.append((ts, files))
    if not result:
        print(f"ERROR: no triplets found for dates: {', '.join(dates)}")
        sys.exit(1)
    return result


# =============================================================================
# STEP 3 — Parse TTI log files
# =============================================================================

TTI_LINE_RE = re.compile(
    r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d+)"
    r"\s+TTI Channel (\d+)"
    r"\s+([\d.]+) V"
    r"\s*;\s*([\d.]+) A"
)


def parse_tti_file(filepath: Path) -> pd.DataFrame:
    records = []
    with open(filepath) as fh:
        for line in fh:
            m = TTI_LINE_RE.match(line.strip())
            if m:
                records.append({
                    "timestamp": datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S.%f"),
                    "channel":   int(m.group(2)),
                    "voltage_V": float(m.group(3)),
                    "current_A": float(m.group(4)),
                })
    return pd.DataFrame(records)


# =============================================================================
# STEP 4 — Parse ENXA connectivity log
# =============================================================================

ENXA_TS_RE = re.compile(r"^(\d{8}-\d{6})\s+Capturing on")


def parse_enxa_file(filepath: Path) -> pd.DataFrame:
    records = []
    current_ts = None
    with open(filepath) as fh:
        for line in fh:
            ts_match = ENXA_TS_RE.match(line)
            if ts_match:
                current_ts = datetime.strptime(ts_match.group(1), "%Y%m%d-%H%M%S")
            elif current_ts is not None and "→" in line:
                records.append({"timestamp": current_ts, "reachable": "ThurlbyThand" in line})
                current_ts = None
    return pd.DataFrame(records)


def load_and_merge(selected: list) -> tuple:
    """
    Parse all selected (timestamp, files) triplets and return
    (df1, df2, df_enxa) merged and sorted by timestamp.
    """
    dfs1, dfs2, dfs_enxa = [], [], []
    for ts, files in selected:
        print(f"  parsing {ts} …")
        dfs1.append(parse_tti_file(files["tti1"]))
        dfs2.append(parse_tti_file(files["tti2"]))
        dfs_enxa.append(parse_enxa_file(files["enxa"]))

    def merge(dfs):
        combined = pd.concat([d for d in dfs if not d.empty], ignore_index=True)
        return combined.sort_values("timestamp").drop_duplicates().reset_index(drop=True)

    df1     = merge(dfs1)
    df2     = merge(dfs2)
    df_enxa = merge(dfs_enxa)

    for label, df in [("TTI1", df1), ("TTI2", df2)]:
        if df.empty:
            print(f"  WARNING: no data for {label}")
        else:
            print(f"  {label}: {len(df)} readings, channels {sorted(df['channel'].unique())}")
    if df_enxa.empty:
        print("  ENXA: no connectivity data found")
    else:
        n_ok = df_enxa["reachable"].sum()
        print(f"  ENXA: {len(df_enxa)} probes, {n_ok} reachable / {len(df_enxa) - n_ok} unreachable")

    return df1, df2, df_enxa


def make_label(selected: list) -> str:
    """Build a short label from the first and last selected timestamps."""
    timestamps = [ts for ts, _ in selected]
    if len(timestamps) == 1:
        return timestamps[0]
    first = datetime.strptime(timestamps[0],  "%Y%m%d-%H%M%S").strftime("%Y%m%d")
    last  = datetime.strptime(timestamps[-1], "%Y%m%d-%H%M%S").strftime("%Y%m%d")
    return first if first == last else f"{first}_{last}"


# =============================================================================
# STEP 5 — Plot and save
# =============================================================================

MIN_ABS_GAP = 0.5


def _split_channels(df: pd.DataFrame, quantity: str) -> tuple:
    channels = sorted(df["channel"].unique())
    if len(channels) <= 1:
        return channels, []
    medians = {ch: df[df["channel"] == ch][quantity].median() for ch in channels}
    sorted_ch = sorted(channels, key=lambda c: medians[c])
    sorted_meds = [medians[c] for c in sorted_ch]
    total_range = sorted_meds[-1] - sorted_meds[0]
    gaps = [sorted_meds[i + 1] - sorted_meds[i] for i in range(len(sorted_meds) - 1)]
    max_gap = max(gaps)
    split_idx = gaps.index(max_gap)
    if total_range > 0 and max_gap > MIN_ABS_GAP and max_gap / total_range > 0.40:
        return sorted_ch[: split_idx + 1], sorted_ch[split_idx + 1:]
    return sorted_ch, []


def _plot_channels(ax, df, quantity, unit_label, title, ch_colors, all_channels) -> list:
    low_chs, high_chs = _split_channels(df, quantity)
    ax2 = ax.twinx() if high_chs else None
    twin_axes = [ax2] if ax2 else []
    legend_lines, legend_labels = [], []

    for ch in low_chs + high_chs:
        color = ch_colors[list(all_channels).index(ch) % len(ch_colors)]
        ch_data = df[df["channel"] == ch].sort_values("timestamp")
        on_twin = ch in high_chs
        suffix = " (right)" if on_twin else ""
        (line,) = (ax2 if on_twin else ax).plot(
            ch_data["timestamp"], ch_data[quantity],
            color=color, linewidth=1.4, marker=".", markersize=3,
            linestyle="--" if on_twin else "-",
            label=f"CH {ch}{suffix}",
        )
        legend_lines.append(line)
        legend_labels.append(f"CH {ch}{suffix}")

    ax.set_title(title, fontsize=11)
    ax.set_ylabel(unit_label)
    ax.grid(True, alpha=0.3)
    if ax2 is not None:
        ax2.set_ylabel(unit_label + "  (right scale)", color="dimgrey", fontsize=8)
        ax2.tick_params(axis="y", labelcolor="dimgrey")
    ax.legend(legend_lines, legend_labels, fontsize=8, loc="best")
    return twin_axes


def _shade_connectivity(axes, df_enxa, all_times) -> None:
    if df_enxa.empty:
        return
    t_end_global = (
        all_times.max() if not all_times.empty
        else df_enxa["timestamp"].iloc[-1] + timedelta(minutes=2)
    )
    for i in range(len(df_enxa)):
        row = df_enxa.iloc[i]
        t_end = df_enxa["timestamp"].iloc[i + 1] if i + 1 < len(df_enxa) else t_end_global
        if not row["reachable"]:
            for ax in axes:
                ax.axvspan(row["timestamp"], t_end, color="lightsalmon", alpha=0.25, zorder=0)


def _fmt_xaxis(axes, all_times) -> None:
    """Use date+time format when data spans more than one day."""
    if all_times.empty:
        return
    span = all_times.max() - all_times.min()
    if span > timedelta(hours=20):
        fmt = mdates.DateFormatter("%b %d\n%H:%M")
    else:
        fmt = mdates.DateFormatter("%H:%M")
    for ax in axes:
        ax.xaxis.set_major_formatter(fmt)
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())


# ── 5d. Time-based overview ───────────────────────────────────────────────────

def plot_triplet(label: str, df1: pd.DataFrame, df2: pd.DataFrame, df_enxa: pd.DataFrame) -> None:
    """Full TTI1+TTI2 overview vs time. Saves lv_monitoring_<label>.png"""
    fig = plt.figure(figsize=(15, 12))
    gs = fig.add_gridspec(3, 2, height_ratios=[2, 2, 1], hspace=0.45, wspace=0.30)
    ax_v1  = fig.add_subplot(gs[0, 0])
    ax_v2  = fig.add_subplot(gs[0, 1])
    ax_i1  = fig.add_subplot(gs[1, 0], sharex=ax_v1)
    ax_i2  = fig.add_subplot(gs[1, 1], sharex=ax_v2)
    ax_net = fig.add_subplot(gs[2, :])

    fig.suptitle(f"LV Monitoring  —  {label.replace('_', ' → ')}", fontsize=14, fontweight="bold", y=0.98)

    all_times = pd.concat([
        df1["timestamp"] if not df1.empty else pd.Series(dtype="object"),
        df2["timestamp"] if not df2.empty else pd.Series(dtype="object"),
    ])

    extra_axes = []
    for df, ax_v, ax_i, title in [(df1, ax_v1, ax_i1, "TTI 1"), (df2, ax_v2, ax_i2, "TTI 2")]:
        if df.empty:
            for ax in (ax_v, ax_i):
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes, color="grey")
            ax_v.set_title(f"{title} — Voltage", fontsize=11)
            ax_i.set_title(f"{title} — Current", fontsize=11)
            continue
        all_chs = sorted(df["channel"].unique())
        extra_axes.extend(
            _plot_channels(ax_v, df, "voltage_V", "Voltage (V)", f"{title} — Voltage", CHANNEL_COLORS, all_chs)
            + _plot_channels(ax_i, df, "current_A", "Current (A)", f"{title} — Current", CHANNEL_COLORS, all_chs)
        )

    _shade_connectivity([ax_v1, ax_i1, ax_v2, ax_i2] + extra_axes, df_enxa, all_times)

    ax_net.set_title("Backend Connectivity  (ENXA ARP probe)", fontsize=11)
    if not df_enxa.empty:
        status = df_enxa["reachable"].astype(int)
        ax_net.step(df_enxa["timestamp"], status, where="post", color="black", linewidth=1.2)
        ax_net.fill_between(df_enxa["timestamp"], status, step="post", alpha=0.25,
                            color=["lightsalmon" if not r else "lightgreen" for r in df_enxa["reachable"]])
        ax_net.set_yticks([0, 1])
        ax_net.set_yticklabels(["Unreachable", "Reachable"], fontsize=9)
        ax_net.set_ylim(-0.2, 1.4)
        ax_net.legend(handles=[
            mpatches.Patch(color="lightgreen",  alpha=0.5, label="Reachable"),
            mpatches.Patch(color="lightsalmon", alpha=0.5, label="Unreachable"),
        ], fontsize=8, loc="upper right")
    else:
        ax_net.text(0.5, 0.5, "No ENXA data available", ha="center", va="center",
                    transform=ax_net.transAxes, color="grey")
    ax_net.grid(True, alpha=0.3)

    _fmt_xaxis([ax_i1, ax_i2, ax_net], all_times)
    fig.autofmt_xdate(rotation=30, ha="right")

    out_dir = DATA_DIR / label
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"lv_monitoring_{label}.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\nPlot saved to:\n  {out_path}")
    if _SHOW:
        plt.show()


# =============================================================================
# STEP 6 — TTI1 vs cumulative radiation dose
# =============================================================================

def _dose_df(df1: pd.DataFrame, irrad_start: datetime) -> pd.DataFrame:
    """Filter to post-irradiation and add dose_krad column."""
    df = df1[df1["timestamp"] >= irrad_start].copy()
    if not df.empty:
        df["dose_krad"] = (
            (df["timestamp"] - irrad_start).dt.total_seconds() / 3600.0
            * DOSE_RATE_KRAD_PER_H
        )
    return df


def plot_tti1_vs_dose(label: str, df1: pd.DataFrame, irrad_start: datetime) -> None:
    """TTI1 voltage + current vs dose. Saves lv_monitoring_<label>_dose.png"""
    df = _dose_df(df1, irrad_start)
    if df.empty:
        print("  No TTI1 data at or after irradiation start — skipping dose plot.")
        return

    fig, (ax_v, ax_i) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle(
        f"TTI1 vs Radiation Dose  —  {label.replace('_', ' → ')}"
        f"   (start {irrad_start.strftime('%Y-%m-%d %H:%M')}, {DOSE_RATE_KRAD_PER_H:.1f} krad/h)",
        fontsize=13, fontweight="bold",
    )
    all_chs = sorted(df["channel"].unique())

    for quantity, ax, ylabel in [("voltage_V", ax_v, "Voltage (V)"), ("current_A", ax_i, "Current (A)")]:
        low_chs, high_chs = _split_channels(df, quantity)
        ax2 = ax.twinx() if high_chs else None
        legend_lines, legend_labels = [], []
        for ch in low_chs + high_chs:
            color = CHANNEL_COLORS[all_chs.index(ch) % len(CHANNEL_COLORS)]
            ch_data = df[df["channel"] == ch].sort_values("dose_krad")
            on_twin = ch in high_chs
            suffix = " (right)" if on_twin else ""
            (line,) = (ax2 if on_twin else ax).plot(
                ch_data["dose_krad"], ch_data[quantity],
                color=color, linewidth=1.4, marker=".", markersize=4,
                linestyle="--" if on_twin else "-", label=f"CH {ch}{suffix}",
            )
            legend_lines.append(line)
            legend_labels.append(f"CH {ch}{suffix}")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        if ax2 is not None:
            ax2.set_ylabel(ylabel + "  (right scale)", color="dimgrey", fontsize=8)
            ax2.tick_params(axis="y", labelcolor="dimgrey")
        ax.legend(legend_lines, legend_labels, fontsize=8, loc="best")

    ax_i.set_xlabel("Cumulative Dose (krad)")
    out_dir = DATA_DIR / label
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"lv_monitoring_{label}_dose.png"
    fig.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Dose plot saved to:\n  {out_path}")
    if _SHOW:
        plt.show()


# =============================================================================
# STEP 7 — TTI1 all-channel current vs dose (single panel)
# =============================================================================

def plot_tti1_current_vs_dose(label: str, df1: pd.DataFrame, irrad_start: datetime) -> None:
    """Single-panel TTI1 current vs dose. Saves lv_monitoring_<label>_current_dose.png"""
    df = _dose_df(df1, irrad_start)
    if df.empty:
        print("  No TTI1 data at or after irradiation start — skipping current-vs-dose plot.")
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(
        f"TTI1 — Current vs Radiation Dose  —  {label.replace('_', ' → ')}"
        f"   (start {irrad_start.strftime('%Y-%m-%d %H:%M')}, {DOSE_RATE_KRAD_PER_H:.1f} krad/h)",
        fontsize=13, fontweight="bold",
    )
    all_chs = sorted(df["channel"].unique())
    low_chs, high_chs = _split_channels(df, "current_A")
    ax2 = ax.twinx() if high_chs else None

    legend_lines, legend_labels = [], []
    for ch in low_chs + high_chs:
        color = CHANNEL_COLORS[all_chs.index(ch) % len(CHANNEL_COLORS)]
        ch_data = df[df["channel"] == ch].sort_values("dose_krad")
        on_twin = ch in high_chs
        suffix = " (right)" if on_twin else ""
        (line,) = (ax2 if on_twin else ax).plot(
            ch_data["dose_krad"], ch_data["current_A"],
            color=color, linewidth=1.6, marker=".", markersize=5,
            linestyle="--" if on_twin else "-", label=f"CH {ch}{suffix}",
        )
        legend_lines.append(line)
        legend_labels.append(f"CH {ch}{suffix}")

    ax.set_xlabel("Cumulative Dose (krad)")
    ax.set_ylabel("Current (A)")
    ax.grid(True, alpha=0.3)
    if ax2 is not None:
        ax2.set_ylabel("Current (A)  (right scale)", color="dimgrey", fontsize=8)
        ax2.tick_params(axis="y", labelcolor="dimgrey")
    ax.legend(legend_lines, legend_labels, fontsize=9, loc="best")

    out_dir = DATA_DIR / label
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"lv_monitoring_{label}_current_dose.png"
    fig.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Current-vs-dose plot saved to:\n  {out_path}")
    if _SHOW:
        plt.show()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LV monitoring curves")
    parser.add_argument(
        "--triplet", "-t", type=int, nargs="+", metavar="N",
        help="select one or more triplet numbers (1-based); omit for interactive menu",
    )
    parser.add_argument(
        "--date", "-d", nargs="+", metavar="YYYYMMDD",
        help="select all triplets whose timestamp starts with these dates",
    )
    parser.add_argument(
        "--list", "-l", action="store_true",
        help="list available triplets and exit",
    )
    args = parser.parse_args()

    triplets = find_triplets(DATA_DIR)
    if not triplets:
        print(f"ERROR: no complete triplets found in\n  {DATA_DIR}")
        sys.exit(1)

    timestamps = sorted(triplets.keys())

    if args.list:
        _print_triplet_menu(timestamps)
        sys.exit(0)

    # ── Select triplets ───────────────────────────────────────────────────────
    if args.triplet:
        selected = select_triplets_by_number(triplets, args.triplet)
    elif args.date:
        selected = select_triplets_by_date(triplets, args.date)
    elif not _SHOW:
        # Headless with no selection argument: just print the list and exit
        _print_triplet_menu(timestamps)
        print("\nRe-run with --triplet N [N ...] or --date YYYYMMDD [...]")
        sys.exit(0)
    else:
        selected = select_triplets_interactive(triplets)

    # ── Load and merge ────────────────────────────────────────────────────────
    print(f"\nLoading {len(selected)} triplet(s) …")
    df1, df2, df_enxa = load_and_merge(selected)

    label = make_label(selected)

    # Irradiation start is anchored to the earliest selected date
    first_date = datetime.strptime(selected[0][0], "%Y%m%d-%H%M%S").date()
    irrad_start = datetime.combine(first_date, datetime.strptime(IRRAD_START_HHMM, "%H:%M").time())

    # ── Plot ──────────────────────────────────────────────────────────────────
    plot_triplet(label, df1, df2, df_enxa)
    plot_tti1_vs_dose(label, df1, irrad_start)
    plot_tti1_current_vs_dose(label, df1, irrad_start)
