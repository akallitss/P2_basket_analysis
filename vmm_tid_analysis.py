#!/usr/bin/env python3
"""
vmm_tid_analysis.py

Total Ionizing Dose (TID) analysis — VMM3a chips, no detector connected.

Two-stage workflow
------------------
  extract  Parse pcapng snapshots → per-file summary CSVs in --processed-dir.
           Runs in parallel (--n-workers).  Resume-safe: already-done files
           are skipped.  Raw hit data is freed immediately after each file.

  plot     Load all summary CSVs from --processed-dir and produce plots.
           No pcapng re-reading required.

  both     (default) Run extract then plot in one call.

Physical observables (no detector = dark noise only)
----------------------------------------------------
  pedestal_median  median ADC of OT=0 hits    →  Vth drift
  noise_sigma      MAD-based ENC proxy         →  leakage / ENC increase
  dark_rate_hz     total hits / capture time   →  spontaneous firing rate
  ot_fraction      fraction of OT=1 hits       →  threshold-crossing rate
  n_active_ch      channels with ≥ 1 hit       →  dead-channel onset

File naming (timestamp auto-extracted, no schedule CSV needed)
-------------------------------------------------------------
  {prefix}_{YYYYMMDD}-{HHMMSS}_{idx}.pcapng  or  .pcap
  e.g.  enxa0cec8769e15_20260421-190835_114.pcapng

Usage
-----
  # Extract only (slow, parallelised, resumable)
  python vmm_tid_analysis.py /path/to/pcapng_dir \\
      --mode extract \\
      --irrad-start "2026-04-20 12:45" --dose-rate 1.0 \\
      --processed-dir processed/ --n-workers 8

  # Plot only (fast, loads CSVs)
  python vmm_tid_analysis.py /path/to/pcapng_dir \\
      --mode plot \\
      --irrad-start "2026-04-20 12:45" --irrad-end "2026-04-22 18:00" \\
      --dose-rate 1.0 --processed-dir processed/ --out-dir results/tid

  # Both in one call — Pagure testbench paths
  python vmm_tid_analysis.py /local/p2/data/pagure_202604 \\
      --run-range 1 8 \\
      --irrad-start "2026-04-20 12:45" --irrad-end "2026-04-22 18:00" \\
      --dose-rate 1.0 \\
      --lv /local/p2/testbench/TestBenchPagure/monitoring \\
      --out-dir results/tid

  # Exclude NoIrad folders (alternative to --run-range when folder names vary)
  python vmm_tid_analysis.py /local/p2/data/pagure_202604 \\
      --exclude NoIrad --mode extract --n-workers 8 --processed-dir processed/
"""

# ── Testbench default paths ────────────────────────────────────────────────────
DATA_DIR_DEFAULT   = "/local/p2/data/pagure_202604"
LV_DIR_DEFAULT     = "/local/p2/testbench/TestBenchPagure/monitoring"

import argparse
import array
import multiprocessing as mp
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
_headless = not os.environ.get("DISPLAY") and sys.platform != "darwin"
if _headless:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SHOW = not _headless

# ── Defaults ─────────────────────────────────────────────────────────────────
IRRAD_START_DEFAULT = datetime(2026, 4, 20, 12, 45)
DOSE_RATE_DEFAULT   = 1.0       # krad / h

# ── Noise quality thresholds (shared with vmm_noise.py conventions) ───────────
ADC_LOW_CUT    = 20    # remove ADC ≤ 20  (VMM3a digital floor artefact)
MIN_NOISE_HITS = 50    # minimum OT=0 hits for a reliable MAD estimate

# ── File naming regex ─────────────────────────────────────────────────────────
_FNAME_RE = re.compile(r".*_(\d{8})-(\d{6})_\d+\.pca(?:p|png)$")

# ── TTI log line regex ────────────────────────────────────────────────────────
_TTI_RE = re.compile(
    r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d+)"
    r"\s+TTI Channel (\d+)\s+([\d.]+) V\s*;\s*([\d.]+) A"
)


# =============================================================================
# 1. PCAPNG PARSER
# =============================================================================

def _parse_vm3_block(block, frame_counter, fec_id,
                      fec_buf, vmm_buf, time_buf, ch_buf, adc_buf, ot_buf):
    if len(block) < 22 or block[4:7] != b'VM3':
        return
    for i in range(0, len(block) - 22, 6):
        d1 = "{:032b}".format(int(block[i+16:i+20].hex(), 16))
        d2 = "{:016b}".format(int(block[i+20:i+22].hex(), 16))
        if d2[0] == '1':
            fec_buf.append(fec_id)
            vmm_buf.append(int(d1[5:10],  2))
            time_buf.append(frame_counter)
            ch_buf.append(int(d2[2:8],   2))
            adc_buf.append(int(d1[10:20], 2))
            ot_buf.append(int(d2[1],      2))


def _detect_src_ips(reader, n_probe=500):
    from scapy.all import UDP, IP
    found = {}
    for i, pkt in enumerate(reader):
        if UDP in pkt and IP in pkt:
            payload = bytes(pkt[UDP].payload)
            if len(payload) > 6 and payload[4:7] == b'VM3':
                ip = pkt[IP].src
                found[ip] = found.get(ip, 0) + 1
        if i >= n_probe:
            break
    return found


def parse_pcapng(filepath, src_ip=None):
    """
    Parse one pcapng/pcap file.

    Returns
    -------
    hits       : pd.DataFrame (fec, vmm, time, ch, adc, over_threshold)
    duration_s : float  wall-clock capture duration from packet timestamps
    """
    from scapy.all import PcapReader, UDP, IP

    filepath = str(filepath)
    fec_buf  = array.array('B')
    vmm_buf  = array.array('B')
    time_buf = array.array('I')
    ch_buf   = array.array('B')
    adc_buf  = array.array('H')
    ot_buf   = array.array('B')
    t_first = t_last = None

    with PcapReader(filepath) as reader:
        if src_ip is None:
            detected = _detect_src_ips(reader, n_probe=500)
            if not detected:
                return pd.DataFrame(), 0.0
            src_ips = set(detected.keys())
        else:
            src_ips = {src_ip}

    with PcapReader(filepath) as reader:
        for pkt in reader:
            if not (UDP in pkt and IP in pkt):
                continue
            if pkt[IP].src not in src_ips:
                continue
            t = float(pkt.time)
            if t_first is None:
                t_first = t
            t_last = t
            payload = bytes(pkt[UDP].payload)
            fc      = int(payload[0:4].hex(), 16)
            fec_id  = int(pkt[IP].src.split('.')[-1])
            _parse_vm3_block(payload, fc, fec_id,
                              fec_buf, vmm_buf, time_buf,
                              ch_buf, adc_buf, ot_buf)

    if not fec_buf:
        return pd.DataFrame(), 0.0

    hits = pd.DataFrame({
        'fec'           : np.frombuffer(fec_buf,  dtype=np.uint8).copy(),
        'vmm'           : np.frombuffer(vmm_buf,  dtype=np.uint8).copy(),
        'time'          : np.frombuffer(time_buf, dtype=np.uint32).copy(),
        'ch'            : np.frombuffer(ch_buf,   dtype=np.uint8).copy(),
        'adc'           : np.frombuffer(adc_buf,  dtype=np.uint16).copy(),
        'over_threshold': np.frombuffer(ot_buf,   dtype=np.uint8).astype(bool).copy(),
    })
    duration_s = (t_last - t_first) if t_first is not None else 0.0
    return hits, duration_s, t_first


# =============================================================================
# 2. METRICS FROM HITS
# =============================================================================

def compute_vmm_metrics(hits, duration_s):
    """Per-VMM summary metrics from a single snapshot's hit DataFrame."""
    records = []
    for vmm_id in sorted(hits["vmm"].unique()):
        df_v = hits[hits["vmm"] == vmm_id]
        n    = len(df_v)
        if n == 0:
            continue

        # Use ALL hits (not just OT=0) — under irradiation ot_fraction→1 so
        # there are no OT=0 hits to filter on; the full ADC distribution
        # still encodes the pedestal position and noise.
        ped = df_v.loc[df_v["adc"] > ADC_LOW_CUT, "adc"].values

        if len(ped) >= MIN_NOISE_HITS:
            med   = float(np.median(ped))
            sigma = 1.4826 * float(np.median(np.abs(ped - med)))
        else:
            med = sigma = np.nan

        records.append({
            "vmm"            : int(vmm_id),
            "n_hits"         : n,
            "pedestal_median": med,
            "noise_sigma"    : sigma,
            "dark_rate_hz"   : n / duration_s if duration_s > 0 else np.nan,
            "ot_fraction"    : float(df_v["over_threshold"].sum()) / n,
            "n_active_ch"    : int(df_v["ch"].nunique()),
        })
    return pd.DataFrame(records)


def compute_channel_metrics(hits):
    """Per-channel noise metrics from a single snapshot's hit DataFrame."""
    records = []
    for vmm_id in sorted(hits["vmm"].unique()):
        df_v = hits[hits["vmm"] == vmm_id]
        for ch_id in range(64):
            ped = df_v.loc[
                (df_v["ch"] == ch_id) &
                (df_v["adc"] > ADC_LOW_CUT), "adc"
            ].values
            n_total = int((df_v["ch"] == ch_id).sum())

            if len(ped) >= MIN_NOISE_HITS:
                med   = float(np.median(ped))
                sigma = 1.4826 * float(np.median(np.abs(ped - med)))
            else:
                med = sigma = np.nan

            records.append({
                "vmm"            : int(vmm_id),
                "ch"             : ch_id,
                "noise_sigma"    : sigma,
                "pedestal_median": med,
                "n_hits"         : n_total,
            })
    return pd.DataFrame(records)


# =============================================================================
# 3. FILE DISCOVERY AND DOSE ASSIGNMENT
# =============================================================================

def find_pcapng_files(data_dir, run_dirs=None, exclude_pattern=None):
    """
    Find pcap/pcapng files under data_dir whose names contain a timestamp
    in the standard format.

    Parameters
    ----------
    run_dirs : list of str or None
        If given, search only these immediate subdirectories of data_dir.
        e.g. ['Run_1', 'Run_2', ..., 'Run_8'].
        None = search all subdirectories (subject to exclude_pattern).
    exclude_pattern : str or None
        Skip any subdirectory whose name contains this string
        (case-insensitive).  e.g. 'NoIrad'.
        Ignored when run_dirs is explicitly provided.

    Returns
    -------
    pd.DataFrame : filepath, filename, capture_time, run_dir
    """
    data_dir = Path(data_dir)

    if run_dirs is not None:
        search_dirs = []
        for rd in run_dirs:
            d = data_dir / rd
            if d.is_dir():
                search_dirs.append(d)
            else:
                print(f"  WARNING: requested run directory not found: {d}")
    else:
        search_dirs = []
        for d in sorted(data_dir.iterdir()):
            if not d.is_dir():
                continue
            if exclude_pattern and exclude_pattern.lower() in d.name.lower():
                print(f"  Skipping (excluded): {d.name}")
                continue
            search_dirs.append(d)
        # also pick up files sitting directly in data_dir
        search_dirs.append(data_dir)

    records = []
    for search_dir in search_dirs:
        run_label = search_dir.name if search_dir != data_dir else "."
        n_before  = len(records)
        for ext in ("*.pcapng", "*.pcap"):
            for fp in sorted(search_dir.rglob(ext)):
                m = _FNAME_RE.match(fp.name)
                if not m:
                    continue
                ts = datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M%S")
                records.append({"filepath": fp, "filename": fp.name,
                                 "capture_time": ts, "run_dir": run_label})
        n_found = len(records) - n_before
        if n_found:
            print(f"  {run_label}: {n_found} file(s)")

    if not records:
        print(f"WARNING: no timestamped pcap/pcapng files found under {data_dir}")
    return (pd.DataFrame(records)
              .sort_values("capture_time")
              .reset_index(drop=True))


def assign_dose(df_files, irrad_start, dose_rate_krad_h, irrad_end=None):
    """
    Add dose_krad, phase, and anneal_h to the file list.

    phase: 'pre' before irrad_start, 'irr' during, 'anneal' after irrad_end.
    Annealing runs inherit the EOI dose on the dose axis.
    """
    df = df_files.copy()
    dt_h = (df["capture_time"] - irrad_start).dt.total_seconds() / 3600.0
    df["dose_krad"] = np.clip(dt_h * dose_rate_krad_h, 0.0, None)
    df["phase"]     = "irr"
    df.loc[df["capture_time"] < irrad_start, ["phase", "dose_krad"]] = ["pre", 0.0]

    df["anneal_h"] = np.nan
    if irrad_end is not None:
        eoi_dose = ((irrad_end - irrad_start).total_seconds() / 3600.0
                    * dose_rate_krad_h)
        ann = df["capture_time"] > irrad_end
        df.loc[ann, "phase"]     = "anneal"
        df.loc[ann, "dose_krad"] = eoi_dose
        df.loc[ann, "anneal_h"]  = (
            (df.loc[ann, "capture_time"] - irrad_end)
            .dt.total_seconds() / 3600.0
        )
    return df


# =============================================================================
# 4. WORKER  (top-level so it is picklable for multiprocessing)
# =============================================================================

def _worker(task):
    """
    Process one pcapng file: parse → compute metrics → save CSVs.

    task : dict with keys
        filepath, dose_krad, phase, anneal_h, capture_time, filename,
        processed_dir, src_ip, skip_channels
    Returns : (vmm_csv_path or None, error_msg or None)
    """
    fp          = Path(task["filepath"])
    proc_dir    = Path(task["processed_dir"])
    stem        = fp.stem
    vmm_out     = proc_dir / f"{stem}_vmm.csv"
    ch_out      = proc_dir / f"{stem}_ch.csv"

    # Resume: skip if already done
    if vmm_out.exists():
        return str(vmm_out), None

    try:
        hits, duration_s, t_first = parse_pcapng(fp, src_ip=task["src_ip"])
        if hits.empty:
            return None, f"no hits: {fp.name}"

        # ── Compute dose from actual packet timestamp (finer than filename) ──
        irrad_start_ts = task.get("irrad_start_ts")
        dose_rate      = task.get("dose_rate", 1.0)
        irrad_end_ts   = task.get("irrad_end_ts")

        if t_first is not None and irrad_start_ts is not None:
            dt_h = (t_first - irrad_start_ts) / 3600.0
            if irrad_end_ts is not None and t_first > irrad_end_ts:
                dose_krad  = (irrad_end_ts - irrad_start_ts) / 3600.0 * dose_rate
                phase      = "anneal"
                anneal_h   = (t_first - irrad_end_ts) / 3600.0
            elif dt_h < 0:
                dose_krad  = 0.0
                phase      = "pre"
                anneal_h   = float("nan")
            else:
                dose_krad  = dt_h * dose_rate
                phase      = "irr"
                anneal_h   = float("nan")
            capture_time = datetime.fromtimestamp(t_first).strftime("%Y-%m-%d %H:%M:%S")
        else:
            dose_krad    = task["dose_krad"]
            phase        = task["phase"]
            anneal_h     = task["anneal_h"]
            capture_time = str(task["capture_time"])

        # ── VMM metrics ──────────────────────────────────────────────────
        df_vmm = compute_vmm_metrics(hits, duration_s)
        df_vmm["dose_krad"]    = dose_krad
        df_vmm["phase"]        = phase
        df_vmm["anneal_h"]     = anneal_h
        df_vmm["capture_time"] = capture_time
        df_vmm["filename"]     = task["filename"]
        df_vmm.to_csv(vmm_out, index=False)

        # ── Channel metrics (optional) ───────────────────────────────────
        if not task["skip_channels"]:
            df_ch = compute_channel_metrics(hits)
            df_ch["dose_krad"]    = dose_krad
            df_ch["phase"]        = phase
            df_ch["capture_time"] = capture_time
            df_ch.to_csv(ch_out, index=False)
            del df_ch

        del hits  # free memory before returning
        return str(vmm_out), None

    except Exception as exc:
        return None, f"{fp.name}: {exc}"


# =============================================================================
# 5. EXTRACTION PIPELINE
# =============================================================================

def run_extract(df_files, processed_dir, src_ip=None,
                 skip_channels=False, n_workers=1,
                 irrad_start=None, dose_rate=1.0, irrad_end=None,
                 stride=1, checkpoint_every=20):
    """
    Parse all pcapng files in df_files, save per-file summary CSVs.

    Parameters
    ----------
    stride         : int  process every Nth file (1 = all, 10 = every 10th)
    checkpoint_every : int  print progress every N completed files
    n_workers      : int  parallel worker processes (1 = sequential)
    """
    proc_dir = Path(processed_dir)
    proc_dir.mkdir(parents=True, exist_ok=True)

    df_sel = df_files.iloc[::stride].reset_index(drop=True)
    n      = len(df_sel)
    print(f"\nExtraction: {n} files  (stride={stride}, workers={n_workers})")

    irrad_start_ts = irrad_start.timestamp() if irrad_start is not None else None
    irrad_end_ts   = irrad_end.timestamp()   if irrad_end   is not None else None

    tasks = [
        {
            "filepath"       : str(row["filepath"]),
            "dose_krad"      : row["dose_krad"],
            "phase"          : row["phase"],
            "anneal_h"       : row["anneal_h"],
            "capture_time"   : row["capture_time"],
            "filename"       : row["filename"],
            "processed_dir"  : str(proc_dir),
            "src_ip"         : src_ip,
            "skip_channels"  : skip_channels,
            "irrad_start_ts" : irrad_start_ts,
            "dose_rate"      : dose_rate,
            "irrad_end_ts"   : irrad_end_ts,
        }
        for _, row in df_sel.iterrows()
    ]

    n_done = n_skip = n_err = 0

    if n_workers > 1:
        ctx = mp.get_context("fork" if sys.platform == "linux" else "spawn")
        with ctx.Pool(processes=n_workers) as pool:
            for i, (out, err) in enumerate(pool.imap_unordered(_worker, tasks), 1):
                if err and out is None:
                    n_err += 1
                    print(f"  ERR [{i}/{n}] {err}")
                elif out is None:
                    n_skip += 1
                else:
                    n_done += 1
                if i % checkpoint_every == 0:
                    print(f"  [{i}/{n}]  done={n_done}  skip={n_skip}  err={n_err}")
    else:
        for i, task in enumerate(tasks, 1):
            out, err = _worker(task)
            if err and out is None:
                n_err += 1
                print(f"  ERR [{i}/{n}] {err}")
            elif out is None:
                n_skip += 1
            else:
                n_done += 1
            if i % checkpoint_every == 0:
                print(f"  [{i}/{n}]  done={n_done}  skip={n_skip}  err={n_err}")

    print(f"Extraction complete: done={n_done}  "
          f"skipped(resume)={n_skip}  errors={n_err}")


# =============================================================================
# 6. AGGREGATION — load processed CSVs
# =============================================================================

def load_metrics(processed_dir):
    """
    Load all *_vmm.csv and *_ch.csv files from processed_dir.

    Returns
    -------
    df_vmm : pd.DataFrame
    df_ch  : pd.DataFrame  (empty if no channel CSVs found)
    """
    proc_dir  = Path(processed_dir)

    vmm_files = sorted(proc_dir.glob("*_vmm.csv"))
    ch_files  = sorted(proc_dir.glob("*_ch.csv"))

    if not vmm_files:
        print(f"WARNING: no *_vmm.csv files in {proc_dir}")
        return pd.DataFrame(), pd.DataFrame()

    df_vmm = pd.concat([pd.read_csv(f) for f in vmm_files],
                        ignore_index=True)
    df_vmm["capture_time"] = pd.to_datetime(df_vmm["capture_time"])
    df_vmm = df_vmm.sort_values("dose_krad").reset_index(drop=True)

    df_ch = pd.DataFrame()
    if ch_files:
        df_ch = pd.concat([pd.read_csv(f) for f in ch_files],
                           ignore_index=True)
        df_ch["capture_time"] = pd.to_datetime(df_ch["capture_time"])
        df_ch = df_ch.sort_values("dose_krad").reset_index(drop=True)

    print(f"Loaded: {len(vmm_files)} VMM CSVs → {len(df_vmm)} rows  |  "
          f"{len(ch_files)} channel CSVs → {len(df_ch)} rows")
    return df_vmm, df_ch


# =============================================================================
# 7. PLOT HELPERS
# =============================================================================

def _vmm_label(vmm_id):
    try:
        from vmm_mapping import vmm_mapping
        for key, cfg in vmm_mapping.items():
            if int(vmm_id) in cfg["vmm_ids"]:
                return f"VMM {vmm_id} — {cfg.get('name', key)}"
    except Exception:
        pass
    return f"VMM {vmm_id}"


def _add_eoi(ax, eoi_dose):
    if eoi_dose is None:
        return
    ax.axvline(eoi_dose, color="crimson", lw=1.2, ls="--",
               zorder=5, label="End of Irradiation")
    xmin, xmax = ax.get_xlim()
    if xmax > eoi_dose:
        ax.axvspan(eoi_dose, xmax, color="gold", alpha=0.10, zorder=0)


def _draw_vmm_lines(ax, df, col, vmm_ids, colors):
    for vmm_id, color in zip(vmm_ids, colors):
        sub = (df[df["vmm"] == vmm_id]
                 .dropna(subset=["dose_krad", col])
                 .sort_values("dose_krad"))
        if sub.empty:
            continue

        # Aggregate multiple files sharing the same dose step → median ± std
        agg = (sub.groupby(["dose_krad", "phase"])[col]
                   .agg(median="median", std="std")
                   .reset_index()
                   .sort_values("dose_krad"))

        pre = agg[agg["phase"] == "pre"]
        irr = agg[agg["phase"] == "irr"]
        ann = agg[agg["phase"] == "anneal"]

        if not pre.empty:
            ax.scatter(pre["dose_krad"], pre["median"],
                       color=color, marker="s", s=55, zorder=5)

        if not irr.empty:
            ax.plot(irr["dose_krad"], irr["median"], color=color,
                    lw=1.5, marker="o", ms=5, zorder=4,
                    label=_vmm_label(vmm_id))
            band = irr[irr["std"].notna() & (irr["std"] > 0)]
            if not band.empty:
                ax.fill_between(band["dose_krad"],
                                band["median"] - band["std"],
                                band["median"] + band["std"],
                                color=color, alpha=0.15, zorder=2)

        if not ann.empty:
            bridge = pd.concat([irr.tail(1), ann]) if not irr.empty else ann
            ax.plot(bridge["dose_krad"], bridge["median"], color=color,
                    lw=1.3, marker="^", ms=5, ls="--", zorder=4)


def _save(fig, out_path):
    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {out_path}")
    if _SHOW:
        plt.show()
    plt.close(fig)


def _metric_plot(df, col, ylabel, title, eoi_dose=None, out_path=None):
    vmm_ids = sorted(df["vmm"].unique())
    colors  = plt.cm.tab10(np.linspace(0, 1, max(len(vmm_ids), 1)))
    fig, ax = plt.subplots(figsize=(10, 5))
    _draw_vmm_lines(ax, df, col, vmm_ids, colors)
    _add_eoi(ax, eoi_dose)
    ax.set_xlabel("Cumulative Dose (krad)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    _save(fig, out_path)


# =============================================================================
# 8. INDIVIDUAL METRIC PLOTS
# =============================================================================

def plot_pedestal_vs_dose(df, eoi_dose=None, out_path=None):
    _metric_plot(df, "pedestal_median",
                 "Pedestal  (median ADC)",
                 "Pedestal Drift vs. Cumulative Dose  [Vth shift]",
                 eoi_dose, out_path)


def plot_noise_vs_dose(df, eoi_dose=None, out_path=None):
    _metric_plot(df, "noise_sigma",
                 "Noise σ  (ADC counts ≈ ENC)",
                 "Equivalent Noise Charge vs. Cumulative Dose",
                 eoi_dose, out_path)


def plot_dark_rate_vs_dose(df, eoi_dose=None, out_path=None):
    _metric_plot(df, "dark_rate_hz",
                 "Dark count rate  (Hz)",
                 "Dark Count Rate vs. Cumulative Dose  [leakage current proxy]",
                 eoi_dose, out_path)


def plot_ot_fraction_vs_dose(df, eoi_dose=None, out_path=None):
    _metric_plot(df, "ot_fraction",
                 "OT fraction",
                 "Over-Threshold Fraction vs. Cumulative Dose  [spontaneous crossings]",
                 eoi_dose, out_path)


def plot_active_channels_vs_dose(df, eoi_dose=None, out_path=None):
    _metric_plot(df, "n_active_ch",
                 "Active channels  (channels with ≥ 1 hit)",
                 "Active Channel Count vs. Cumulative Dose  [dead-channel onset]",
                 eoi_dose, out_path)


# =============================================================================
# 9. CHANNEL NOISE HEATMAP  (channel × dose step)
# =============================================================================

def plot_channel_heatmap(df_ch, vmm_id, metric="noise_sigma",
                          eoi_dose=None, out_path=None):
    """
    2-D heatmap: channel (y-axis) × dose step (x-axis).
    Grey = missing / too few hits.  Crimson dashed line = EOI.
    """
    sub = df_ch[df_ch["vmm"] == vmm_id].copy()
    if sub.empty:
        return

    dose_steps = sorted(sub["dose_krad"].unique())
    grid = np.full((64, len(dose_steps)), np.nan)
    for col_i, dose in enumerate(dose_steps):
        step = sub[sub["dose_krad"] == dose]
        # Multiple files share the same dose step → take median per channel
        for ch_id in range(64):
            vals = step.loc[step["ch"] == ch_id, metric].dropna()
            if len(vals) > 0:
                grid[ch_id, col_i] = float(vals.median())

    fig, ax = plt.subplots(figsize=(max(8, len(dose_steps) * 0.55), 7))
    cmap = plt.cm.plasma.copy()
    cmap.set_bad("lightgrey")
    vmax = np.nanpercentile(grid, 95) or 1.0
    vmin = np.nanmin(grid[grid > 0]) if np.any(grid > 0) else 0.0

    im = ax.imshow(grid, aspect="auto", origin="lower",
                   cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=[-0.5, len(dose_steps) - 0.5, -0.5, 63.5])
    plt.colorbar(im, ax=ax,
                 label=("Noise σ (ADC)" if metric == "noise_sigma"
                         else "Pedestal median (ADC)"))

    ax.set_xticks(range(len(dose_steps)))
    ax.set_xticklabels([f"{d:.2f}" for d in dose_steps],
                        rotation=45, ha="right", fontsize=7)
    ax.set_xlabel("Cumulative Dose (krad)")
    ax.set_ylabel("Channel")
    ax.set_title(f"VMM {vmm_id} — {metric} per channel vs. dose",
                  fontweight="bold")

    if eoi_dose is not None and eoi_dose in dose_steps:
        ax.axvline(dose_steps.index(eoi_dose), color="crimson",
                   lw=1.5, ls="--", label="EOI")
        ax.legend(fontsize=8)
    fig.tight_layout()
    _save(fig, out_path)


# =============================================================================
# 10. LV CURRENT
# =============================================================================

def load_lv_current(lv_log_dir, irrad_start, dose_rate_krad_h, irrad_end=None):
    records = []
    for fpath in sorted(Path(lv_log_dir).glob("tti1_*_mon.log")):
        with open(fpath) as fh:
            for line in fh:
                m = _TTI_RE.match(line.strip())
                if m:
                    records.append({
                        "timestamp": datetime.strptime(m.group(1),
                                                       "%Y-%m-%d %H:%M:%S.%f"),
                        "channel":   int(m.group(2)),
                        "current_A": float(m.group(4)),
                    })
    if not records:
        return pd.DataFrame()

    df = (pd.DataFrame(records)
            .sort_values("timestamp")
            .drop_duplicates()
            .reset_index(drop=True))
    df = df[df["timestamp"] >= irrad_start].copy()
    if df.empty:
        return df

    eoi_dose = None
    if irrad_end is not None:
        eoi_dose = ((irrad_end - irrad_start).total_seconds() / 3600.0
                    * dose_rate_krad_h)

    dt_h = (df["timestamp"] - irrad_start).dt.total_seconds() / 3600.0
    df["dose_krad"] = np.clip(dt_h * dose_rate_krad_h, 0.0, None)
    if eoi_dose is not None:
        df.loc[df["timestamp"] > irrad_end, "dose_krad"] = eoi_dose
    return df[["timestamp", "dose_krad", "channel", "current_A"]]


def plot_lv_current_vs_dose(df_lv, eoi_dose=None, out_path=None):
    if df_lv.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    channels = sorted(df_lv["channel"].unique())
    colors   = plt.cm.tab10(np.linspace(0, 1, max(len(channels), 1)))
    for ch, color in zip(channels, colors):
        s = df_lv[df_lv["channel"] == ch].sort_values("dose_krad")
        ax.plot(s["dose_krad"], s["current_A"] * 1e3,
                color=color, lw=1.0, label=f"CH {ch}")
    _add_eoi(ax, eoi_dose)
    ax.set_xlabel("Cumulative Dose (krad)")
    ax.set_ylabel("Supply Current  (mA)")
    ax.set_title("LV Supply Current vs. Cumulative Dose", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    _save(fig, out_path)


# =============================================================================
# 11. DASHBOARD
# =============================================================================

def plot_tid_dashboard(df_vmm, df_lv, eoi_dose=None, label="", out_path=None):
    """
    2×3 panel dashboard.
      Row 0: pedestal drift | noise σ | LV current
      Row 1: dark rate      | OT fraction | active channels
    """
    vmm_ids = sorted(df_vmm["vmm"].unique()) if not df_vmm.empty else []
    colors  = plt.cm.tab10(np.linspace(0, 1, max(len(vmm_ids), 1)))

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    suptitle = "VMM3a  —  TID Pedestal & Noise Monitor"
    if label:
        suptitle += f"\n{label}"
    fig.suptitle(suptitle, fontsize=14, fontweight="bold", y=1.01)

    metric_panels = [
        (axes[0, 0], "pedestal_median", "Pedestal (ADC)",  "Pedestal Drift"),
        (axes[0, 1], "noise_sigma",     "Noise σ (ADC)",   "ENC vs. Dose"),
        (axes[1, 0], "dark_rate_hz",    "Dark rate (Hz)",  "Dark Count Rate"),
        (axes[1, 1], "ot_fraction",     "OT fraction",     "OT Fraction"),
        (axes[1, 2], "n_active_ch",     "Active channels", "Active Channels"),
    ]
    for ax, col, ylabel, ptitle in metric_panels:
        if not df_vmm.empty and col in df_vmm.columns:
            _draw_vmm_lines(ax, df_vmm, col, vmm_ids, colors)
        else:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes, color="grey")
        _add_eoi(ax, eoi_dose)
        ax.set_xlabel("Dose (krad)", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(ptitle, fontsize=10)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    ax_lv = axes[0, 2]
    if not df_lv.empty:
        ch_list = sorted(df_lv["channel"].unique())
        lv_colors = plt.cm.tab10(np.linspace(0, 1, max(len(ch_list), 1)))
        for ch, color in zip(ch_list, lv_colors):
            s = df_lv[df_lv["channel"] == ch].sort_values("dose_krad")
            ax_lv.plot(s["dose_krad"], s["current_A"] * 1e3,
                       color=color, lw=0.9, label=f"CH {ch}")
        _add_eoi(ax_lv, eoi_dose)
    else:
        ax_lv.text(0.5, 0.5, "No LV data", ha="center", va="center",
                   transform=ax_lv.transAxes, color="grey")
    ax_lv.set_xlabel("Dose (krad)", fontsize=9)
    ax_lv.set_ylabel("Current (mA)", fontsize=9)
    ax_lv.set_title("LV Supply Current", fontsize=10)
    ax_lv.legend(fontsize=7)
    ax_lv.grid(True, alpha=0.3)

    phase_handles = [
        plt.Line2D([0],[0], marker="s", color="grey", ls="none", ms=6,
                   label="Pre-irradiation"),
        plt.Line2D([0],[0], marker="o", color="grey", ls="-",  lw=1.4, ms=4,
                   label="Irradiation"),
        plt.Line2D([0],[0], marker="^", color="grey", ls="--", lw=1.2, ms=5,
                   label="Annealing"),
    ]
    if eoi_dose is not None:
        phase_handles += [
            plt.Line2D([0],[0], color="crimson", ls="--", lw=1.2, label="EOI"),
            plt.Rectangle((0,0),1,1, facecolor="gold", alpha=0.3,
                           label="Annealing region"),
        ]
    fig.legend(handles=phase_handles, loc="lower center",
               ncol=len(phase_handles), fontsize=9,
               bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout()
    _save(fig, out_path)


# =============================================================================
# 12. DEGRADATION SUMMARY
# =============================================================================

def compute_degradation_summary(df_vmm, pre_dose_krad=0.5):
    """
    % change at EOI and % recovery after annealing vs. pre-irradiation baseline.
    """
    metrics = ["pedestal_median", "noise_sigma", "dark_rate_hz", "ot_fraction"]
    records = []

    for vmm_id in sorted(df_vmm["vmm"].unique()):
        sub = df_vmm[df_vmm["vmm"] == vmm_id].sort_values("dose_krad")

        for metric in metrics:
            if metric not in sub.columns:
                continue
            s = sub.dropna(subset=[metric])
            if s.empty:
                continue

            pre      = s.loc[s["dose_krad"] < pre_dose_krad, metric]
            baseline = pre.mean() if not pre.empty else np.nan
            if np.isnan(baseline):
                continue

            irr     = s[s["phase"] == "irr"]
            eoi_val = (irr.loc[irr["dose_krad"].idxmax(), metric]
                       if not irr.empty else np.nan)
            ann     = s[s["phase"] == "anneal"]
            ann_val = ann[metric].mean() if not ann.empty else np.nan

            change   = 100 * (eoi_val - baseline) / baseline if baseline else np.nan
            recovery = np.nan
            if not np.isnan(ann_val) and (baseline - eoi_val) != 0:
                recovery = 100 * (ann_val - eoi_val) / (baseline - eoi_val)

            records.append({
                "vmm"         : vmm_id,
                "metric"      : metric,
                "baseline"    : round(baseline, 3),
                "eoi_value"   : round(eoi_val, 3) if not np.isnan(eoi_val) else np.nan,
                "change_pct"  : round(change,  1) if not np.isnan(change)  else np.nan,
                "anneal_value": round(ann_val, 3) if not np.isnan(ann_val) else np.nan,
                "recovery_pct": round(recovery,1) if not np.isnan(recovery)else np.nan,
            })

    return pd.DataFrame(records)


# =============================================================================
# 13. MAIN
# =============================================================================

def main():
    mp.freeze_support()   # needed for Windows frozen executables

    parser = argparse.ArgumentParser(
        description="VMM3a TID analysis  (extract → processed CSVs → plots)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("data_dir", metavar="DATA_DIR",
                        help="root directory containing pcap/pcapng snapshots")
    parser.add_argument("--mode", choices=["extract", "plot", "both"],
                        default="both",
                        help="extract: parse files only  |  "
                             "plot: load CSVs and plot  |  "
                             "both: do both  (default)")
    parser.add_argument("--irrad-start", metavar="DATETIME",
                        default=IRRAD_START_DEFAULT.strftime("%Y-%m-%d %H:%M"),
                        help="irradiation start  YYYY-MM-DD HH:MM  "
                             "(default: %(default)s)")
    parser.add_argument("--irrad-end", metavar="DATETIME", default=None,
                        help="end of irradiation  YYYY-MM-DD HH:MM  "
                             "(enables annealing phase labelling)")
    parser.add_argument("--dose-rate", type=float, default=DOSE_RATE_DEFAULT,
                        metavar="KRAD_PER_H",
                        help="dose rate in krad/h  (default: %(default)s)")
    parser.add_argument("--processed-dir", metavar="DIR", default="processed",
                        help="directory for per-file summary CSVs  "
                             "(default: ./processed)")
    parser.add_argument("--out-dir", metavar="DIR", default=".",
                        help="output directory for plots  (default: .)")
    parser.add_argument("--lv", metavar="DIR", default=LV_DIR_DEFAULT,
                        help="directory with tti1_*_mon.log files  "
                             "(default: %(default)s)")
    parser.add_argument("--src-ip", metavar="IP", default=None,
                        help="force FEC source IP (auto-detected if omitted)")
    parser.add_argument("--n-workers", type=int, default=1,
                        metavar="N",
                        help="parallel worker processes for extraction  "
                             "(default: 1 = sequential)")
    parser.add_argument("--stride", type=int, default=1,
                        metavar="N",
                        help="use every Nth file  (1=all, 60=one per minute → hourly)")
    parser.add_argument("--run-range", nargs=2, type=int, metavar=("FIRST", "LAST"),
                        default=None,
                        help="include only Run_FIRST through Run_LAST subdirectories "
                             "(e.g. --run-range 1 8 → Run_1 … Run_8)")
    parser.add_argument("--exclude", metavar="PATTERN", default=None,
                        help="exclude subdirectories whose name contains PATTERN "
                             "(e.g. --exclude NoIrad)")
    parser.add_argument("--min-dose", type=float, default=None, metavar="KRAD",
                        help="skip files with estimated dose below this value (krad)")
    parser.add_argument("--max-dose", type=float, default=None, metavar="KRAD",
                        help="skip files with estimated dose above this value (krad)")
    parser.add_argument("--no-channels", action="store_true",
                        help="skip per-channel metrics (faster for large datasets)")
    parser.add_argument("--heatmap-vmms", nargs="+", type=int, metavar="VMM",
                        default=None,
                        help="VMM IDs to generate channel heatmaps for "
                             "(default: all)")
    args = parser.parse_args()

    irrad_start = datetime.strptime(args.irrad_start, "%Y-%m-%d %H:%M")
    irrad_end   = (datetime.strptime(args.irrad_end, "%Y-%m-%d %H:%M")
                   if args.irrad_end else None)
    dose_rate   = args.dose_rate
    out_dir     = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    eoi_dose = None
    if irrad_end is not None:
        eoi_dose = ((irrad_end - irrad_start).total_seconds() / 3600.0
                    * dose_rate)
        print(f"EOI dose: {eoi_dose:.2f} krad")

    # ── File discovery ────────────────────────────────────────────────────────
    run_dirs = None
    if args.run_range:
        first, last = args.run_range
        run_dirs = [f"Run_{i}" for i in range(first, last + 1)]

    print(f"\nScanning {args.data_dir} ...")
    df_files = find_pcapng_files(args.data_dir,
                                 run_dirs=run_dirs,
                                 exclude_pattern=args.exclude)
    if df_files.empty:
        print("No files found — exiting.")
        sys.exit(1)

    df_files = assign_dose(df_files, irrad_start, dose_rate, irrad_end)

    if args.min_dose is not None:
        df_files = df_files[df_files["dose_krad"] >= args.min_dose].reset_index(drop=True)
    if args.max_dose is not None:
        df_files = df_files[df_files["dose_krad"] <= args.max_dose].reset_index(drop=True)
    if df_files.empty:
        print("No files in the requested dose range — exiting.")
        sys.exit(1)

    n_total  = len(df_files)
    for ph in ("pre", "irr", "anneal"):
        sub = df_files[df_files["phase"] == ph]
        if not sub.empty:
            print(f"  {ph:8s}: {len(sub):4d} files  "
                  f"{sub['dose_krad'].min():.2f} – {sub['dose_krad'].max():.2f} krad")

    # ── Extract ───────────────────────────────────────────────────────────────
    if args.mode in ("extract", "both"):
        run_extract(
            df_files,
            processed_dir   = args.processed_dir,
            src_ip          = args.src_ip,
            skip_channels   = args.no_channels,
            n_workers       = args.n_workers,
            stride          = args.stride,
            checkpoint_every= 20,
            irrad_start     = irrad_start,
            dose_rate       = dose_rate,
            irrad_end       = irrad_end,
        )

    # ── Plot ──────────────────────────────────────────────────────────────────
    if args.mode in ("plot", "both"):
        df_vmm, df_ch = load_metrics(args.processed_dir)
        if df_vmm.empty:
            print("No processed metrics found — run with --mode extract first.")
            sys.exit(1)

        df_lv = pd.DataFrame()
        if args.lv:
            df_lv = load_lv_current(args.lv, irrad_start, dose_rate, irrad_end)
            print(f"LV: {len(df_lv)} readings")

        print("\nGenerating plots ...")

        plot_pedestal_vs_dose(df_vmm, eoi_dose,
            out_path=str(out_dir / "tid_pedestal_vs_dose.png"))
        plot_noise_vs_dose(df_vmm, eoi_dose,
            out_path=str(out_dir / "tid_noise_vs_dose.png"))
        plot_dark_rate_vs_dose(df_vmm, eoi_dose,
            out_path=str(out_dir / "tid_dark_rate_vs_dose.png"))
        plot_ot_fraction_vs_dose(df_vmm, eoi_dose,
            out_path=str(out_dir / "tid_ot_fraction_vs_dose.png"))
        plot_active_channels_vs_dose(df_vmm, eoi_dose,
            out_path=str(out_dir / "tid_active_channels_vs_dose.png"))

        if not df_lv.empty:
            plot_lv_current_vs_dose(df_lv, eoi_dose,
                out_path=str(out_dir / "tid_lv_current_vs_dose.png"))

        if not df_ch.empty:
            hmap_vmms = (args.heatmap_vmms
                         if args.heatmap_vmms
                         else sorted(df_ch["vmm"].unique()))
            for vmm_id in hmap_vmms:
                plot_channel_heatmap(
                    df_ch, vmm_id, metric="noise_sigma",
                    eoi_dose=eoi_dose,
                    out_path=str(out_dir / f"tid_ch_noise_heatmap_vmm{vmm_id}.png"))
                plot_channel_heatmap(
                    df_ch, vmm_id, metric="pedestal_median",
                    eoi_dose=eoi_dose,
                    out_path=str(out_dir / f"tid_ch_pedestal_heatmap_vmm{vmm_id}.png"))

        lbl = (f"{irrad_start.strftime('%Y-%m-%d')}  |  {dose_rate:.1f} krad/h"
               + (f"  |  EOI = {eoi_dose:.1f} krad" if eoi_dose else ""))
        plot_tid_dashboard(df_vmm, df_lv, eoi_dose=eoi_dose, label=lbl,
                            out_path=str(out_dir / "tid_dashboard.png"))

        df_deg = compute_degradation_summary(df_vmm)
        if not df_deg.empty:
            deg_csv = out_dir / "tid_degradation_summary.csv"
            df_deg.to_csv(deg_csv, index=False)
            print(f"\nDegradation summary → {deg_csv}")
            print(df_deg.to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()
