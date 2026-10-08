#!/usr/bin/env python3
"""
Signal loss vs off-detector cable length: external-pulser runs (June 2026) vs Fe55 on
the P2 detector (Sep/Oct 2026), and what the "edge-channel breakdown" seen on
the external-pulser runs really is.

Two campaigns, one observable per channel each:

  pulser   hybrid_performance/ (vmm_daplxa), hybrid 4118, VMM id 7, Zedboard
           pulser at 1 kHz into every third channel at once (3k, 3k+1, 3k+2),
           cable = 50 cm injection + X m transmission. Read from the converted
           `hits` trees (the June pcapng are owner-only). Per channel:
             * in-time amplitude: hits within +-500 ns of a pulse, the pulse
               clock taken from clean mid-connector channels of the same run
             * in-time efficiency: in-time hits per pulse
             * out-of-time rate: everything else, i.e. self-triggering
           The pulsed residue is detected from the data: the labels 3k+1/3k+2
           are swapped in most June runs but not in all of them.

  fe55     source_tests/ (vmm_daplxa), det on 5 hybrids (VMM 4-13), pcapng
           decoded with vmm_decode. Per channel: Fe55 photopeak position
           (windowed median of the over-threshold spectrum), and the
           over-threshold rate in the no-source reference runs. Per-file
           photopeaks are kept for the gain-stability check.

Edge channels: the VMM channels at the connector borders, 0, 31, 32, 63.

    python3 cable_length/cable_length_study.py [--refresh] [--jobs N]

Caches per-run tables in cable_length/cache/, writes figures and
summary.txt to cable_length/figs/ (both regenerable, gitignored).
"""
import argparse
import datetime as dt
import glob
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = ("/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/"
        "LAB_Measurements/cable_length")
PULSER_DIR = os.path.join(DATA, "pulser_hybrid_performance")
# 7 Oct runs: the curated source_tests/cable_length/ set (byte-identical copies
# of the source_tests/ files where they overlap); 29 Sep refs: source_tests/
FE55_DIRS = (os.path.join(DATA, "fe55_cable_length"), os.path.join(DATA, "fe55_source_tests"))
ONLINE_REPO = "/local/home/ak271430/Documents/PostDocSaclay/P2_basket_online_analysis"
sys.path.insert(0, ONLINE_REPO)
from vmm_decode import iter_chunks              # noqa: E402
from vmm_det1_cabling import capture_window     # noqa: E402
from vmm_fe55_map import windowed_median_peak   # noqa: E402

CACHE = os.path.join(HERE, "cache")
FIGS = os.path.join(HERE, "figs")

EDGE = (0, 31, 32, 63)
PULSER_VMM = 7
FE55_VMMS = tuple(range(4, 14))
INTIME_NS = 500.0
PULSE_GAP_NS = 1000.0         # reference hits closer than this are one pulse
MIN_PULSES = 5000             # 10 s at 1 kHz; fewer = the pulser did not run
# mid-connector channels: each VMM's 64 channels leave on two connector blocks
# (ch 0-31, 32-63); keep >= 3 channels away from either end, i.e. ch 3-27 and 36-60
MID = [c for c in range(3, 61) if not 28 <= c <= 35]

# --- pulser runs: (cable, total length m) -> directory pattern ------------
# total = 0.5 m injection + transmission. Logbook: June 5 and June 8 morning
# transmission cables are Samtec; the 2.5 m one is not named in the logbook
# but was taken in the same series (assumed Samtec, flagged in the plots).
PULSER_SETS = {
    ("Samtec", 1.0): "hybrid_4118_test_SNR_SNG_AC_2cables_50cm_50cm_channels_{g}_vmm1_{v}_2026*",
    ("Samtec", 1.5): "hybrid_4118_test_SNR_SNG_AC_2cables_50cm_1meter_channels_{g}_vmm1_{v}_2026*",
    ("Samtec", 2.0): "hybrid_4118_test_SNR_SNG_AC_2cables_50cm_1meter5_channels_{g}_vmm1_{v}_2026*",
    ("Samtec", 2.5): "hybrid_4118_test_SNR_SNG_AC_2cables_50cm_2meter_channels_{g}_vmm1_{v}_2026*",
    ("Samtec", 3.0): "hybrid_4118_test_SNR_SNG_AC_2cables_50cm_2meter5_channels_{g}_vmm1_{v}_2026*",
    ("Hitachi", 1.5): "hybrid_hybrid_4118_test_SNR_SNG_AC_hitachi_50cm_1meter_channels_{g}_vmm1_{v}_2026*",
    ("Hitachi", 2.5): "hybrid_hybrid_4118_test_SNR_SNG_AC_hitachi_50cm_2meter_channels_{g}_vmm1_{v}_2026*",
}
ASSUMED_CABLE = {("Samtec", 3.0)}
GROUPS = ("3k", "3k+1", "3k+2")
VOLTS = ("1V8", "2V5", "3V3")

# --- fe55 runs (7 Oct is the clean set; 29 Sep saturates, refs only) -----
FE55_RUNS = [
    # label, kind, cable, length, gain mV/fC, stem
    ("Hitachi 2.0 m (A)", "fe55", "Hitachi", 2.0, 1, "enp4s0f1_enp4s0f1_fe55_hitachi_2.0m_m410d750_sng_1mVfC_20261007-114534"),
    ("Hitachi 2.0 m (B)", "fe55", "Hitachi", 2.0, 1, "enp4s0f1_enp4s0f1_fe55_hitachi_2.0m_m410d750_sng_1mVfC_20261007-125430"),
    ("Samtec 1.5 m", "fe55", "Samtec", 1.5, 1, "enp4s0f1_enp4s0f1_fe55_samtec_1.5m_m410d750_sng_1mVfC_20261007-151716"),
    ("Samtec 2.0 m", "fe55", "Samtec", 2.0, 1, "enp4s0f1_enp4s0f1_fe55_samtec_2.0m_m410d750_sng_1mVfC_20261007-162811"),
    ("Samtec 2.5 m", "fe55", "Samtec", 2.5, 1, "enp4s0f1_enp4s0f1_fe55_samtec_2.5m_m410d750_sng_1mVfC_20261007-164656"),
    ("Hitachi 2.0 m", "fe55", "Hitachi", 2.0, 3, "enp4s0f1_enp4s0f1_fe55_hitachi_2.0m_m410d750_sng_3mVfC_20261007-135852"),
    ("Samtec 1.5 m", "fe55", "Samtec", 1.5, 3, "enp4s0f1_enp4s0f1_fe55_samtec_1.5m_m410d750_sng_3mVfC_20261007-144517"),
    ("Samtec 2.0 m", "fe55", "Samtec", 2.0, 3, "enp4s0f1_enp4s0f1_fe55_samtec_2.0m_m410d750_sng_3mVfC_20261007-155347"),
    ("Samtec 2.5 m", "fe55", "Samtec", 2.5, 3, "enp4s0f1_enp4s0f1_fe55_samtec_2.5m_m410d750_sng_3mVfC_20261007-171625"),
    ("Hitachi 2.0 m (r3)", "ref", "Hitachi", 2.0, 1, "enp4s0f1_enp4s0f1_ref3_hitachi_2.0m_m410d750_sng_1mVfC_20261007-113537"),
    ("Hitachi 2.0 m (r4)", "ref", "Hitachi", 2.0, 1, "enp4s0f1_enp4s0f1_ref4_hitachi_2.0m_m410d750_sng_1mVfC_20261007-134400"),
    ("Samtec 1.5 m", "ref", "Samtec", 1.5, 1, "enp4s0f1_enp4s0f1_ref4_samtec_1.5m_m410d750_sng_1mVfC_20261007-150833"),
    ("Samtec 2.0 m", "ref", "Samtec", 2.0, 1, "enp4s0f1_enp4s0f1_ref4_samtec_2.0m_m410d750_sng_1mVfC_20261007-161257"),
    ("Samtec 2.5 m", "ref", "Samtec", 2.5, 1, "enp4s0f1_enp4s0f1_ref4_samtec_2.5m_m410d750_sng_1mVfC_20261007-170053"),
    ("Hitachi 2.0 m", "ref", "Hitachi", 2.0, 3, "enp4s0f1_enp4s0f1_ref4_hitachi_2.0m_m410d750_sng_3mVfC_20261007-135238"),
    ("Samtec 1.5 m", "ref", "Samtec", 1.5, 3, "enp4s0f1_enp4s0f1_ref4_samtec_1.5m_m410d750_sng_3mVfC_20261007-145714"),
    ("Samtec 2.0 m", "ref", "Samtec", 2.0, 3, "enp4s0f1_enp4s0f1_ref4_samtec_2.0m_m410d750_sng_3mVfC_20261007-160431"),
    ("Samtec 2.5 m", "ref", "Samtec", 2.5, 3, "enp4s0f1_enp4s0f1_ref4_samtec_2.5m_m410d750_sng_3mVfC_20261007-170954"),
    ("29 Sep 1.5 m", "ref", "?", 1.5, 3, "enp4s0f1_ref_1.5m_5hyb_20260929-143010"),
    ("29 Sep 2.0 m", "ref", "?", 2.0, 3, "enp4s0f1_ref_2.0m_5hyb_20260929-145448"),
    ("29 Sep 2.5 m", "ref", "?", 2.5, 3, "enp4s0f1_ref_2.5m_5hyb_20260929-161201"),
]

# --- style: palette slots from the dataviz reference instance --------------
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
CABLE_COLOR = {"Samtec": "#2a78d6", "Hitachi": "#eb6834"}
# ordinal blue ramp (steps 250..650) keyed on total length
LEN_COLOR = {1.0: "#86b6ef", 1.5: "#5598e7", 2.0: "#2a78d6", 2.5: "#1c5cab", 3.0: "#104281"}
EDGE_SHADE = "#f0efec"
plt.rcParams.update({
    "font.size": 12, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "figure.dpi": 110,
    "savefig.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
})


# ======================================================================
# pulser
# ======================================================================
def _pulser_run(d):
    """Per-channel in-time amplitude / efficiency / out-of-time rate of one run."""
    f = [x for x in glob.glob(os.path.join(d, "*.root")) if "validation" not in x][0]
    a = uproot.open(f)["hits"].arrays(["hits/vmm", "hits/ch", "hits/adc", "hits/time"],
                                      library="np")
    v, c, adc, t = (a[k] for k in ("hits/vmm", "hits/ch", "hits/adc", "hits/time"))
    live_s = (t.max() - t.min()) * 1e-9
    m = v == PULSER_VMM
    c, adc, t = c[m], adc[m].astype(float), t[m]
    # the pulsed residue is the one whose mid-connector channels sit far above the rest
    meds = [np.median(adc[(c % 3 == r) & np.isin(c, MID)]) for r in range(3)]
    res = int(np.argmax(meds))
    ref = np.sort(t[(c % 3 == res) & np.isin(c, MID)])
    pulses = ref[np.r_[True, np.diff(ref) > PULSE_GAP_NS]]
    rows = []
    # neighbour-trigger readout of the unpulsed channels at pulse time: baseline + noise + crosstalk
    h_nb = np.zeros((64, 1024), np.int32)
    for ch in range(64):
        sel = c == ch
        tt, aa = t[sel], adc[sel]
        if len(tt):
            i = np.clip(np.searchsorted(pulses, tt), 1, len(pulses) - 1)
            it = np.minimum(abs(tt - pulses[i - 1]), abs(tt - pulses[i])) < INTIME_NS
        else:
            it = np.zeros(0, bool)
        pulsed = ch % 3 == res
        if not pulsed and it.any():
            h_nb[ch] = np.bincount(np.clip(aa[it].astype(int), 0, 1023), minlength=1024)
        rows.append(dict(
            ch=ch, pulsed=pulsed, n_hits=len(tt),
            amp=np.median(aa[it]) if pulsed and it.any() else np.nan,
            eff=it.sum() / len(pulses) if pulsed else np.nan,
            # off-pulse hits on a pulsed channel = self-triggering. On an
            # unpulsed channel in-time hits are neighbour-trigger readout.
            noise_hz=(~it).sum() / live_s,
            noise_adc=np.median(aa[~it]) if (~it).any() else np.nan))
    return pd.DataFrame(rows), res, len(pulses), live_s, h_nb


def pulser_table(refresh=False):
    out = os.path.join(CACHE, "pulser_channels.csv")
    out_nb = os.path.join(CACHE, "pulser_neighbour_hist.npz")
    if os.path.exists(out) and os.path.exists(out_nb) and not refresh:
        return pd.read_csv(out)
    frames, nb = [], {}
    for (cable, length), pat in PULSER_SETS.items():
        for volt in VOLTS:
            for g in GROUPS:
                ds = [d for d in sorted(glob.glob(os.path.join(PULSER_DIR, pat.format(g=g, v=volt))))
                      if "twisted" not in d]
                if not ds:
                    continue
                d = ds[-1]       # latest take of that configuration
                df, res, npul, live, h_nb = _pulser_run(d)
                if npul >= MIN_PULSES:   # pooled over the three groups of a configuration
                    key = f"{cable}|{length}|{volt}"
                    nb[key] = nb.get(key, 0) + h_nb
                df = df.assign(cable=cable, length=length, volt=volt, group=g,
                               residue=res, n_pulses=npul, live_s=live,
                               ok=npul >= MIN_PULSES, run=os.path.basename(d))
                frames.append(df)
                print(f"pulser {cable:7s} {length} m {volt} {g:4s} -> residue {res}  "
                      f"{npul} pulses  {os.path.basename(d)}")
    tab = pd.concat(frames, ignore_index=True)
    tab.to_csv(out, index=False)
    np.savez_compressed(out_nb, **nb)
    return tab


def pulser_neighbour_hist():
    """{(cable, length, volt): 64 x 1024 ADC histogram of in-time neighbour hits}."""
    z = np.load(os.path.join(CACHE, "pulser_neighbour_hist.npz"))
    return {(k.split("|")[0], float(k.split("|")[1]), k.split("|")[2]): z[k] for k in z.files}


def hist_level_width(h, nmin=30):
    """Median and half 16-84 % width of ADC histograms along the last axis (NaN below nmin)."""
    n = h.sum(-1)
    cdf = np.cumsum(h, -1) / np.maximum(n, 1)[..., None]
    q = lambda p: (cdf < p).sum(-1).astype(float)
    med, wid = q(0.5), (q(0.84) - q(0.16)) / 2
    bad = n < nmin
    med[bad] = np.nan
    wid[bad] = np.nan
    return med, wid, n


def hist_mean(h, nmin=30):
    """Mean ADC of histograms along the last axis (NaN below nmin)."""
    n = h.sum(-1)
    m = (h * np.arange(h.shape[-1])).sum(-1) / np.maximum(n, 1)
    return np.where(n >= nmin, m, np.nan)


def pulser_per_channel(tab, volt):
    """One row per (cable, length, ch) with the measurement from the run that pulsed it.

    Noise is the out-of-time rate on the pulsed run, so every channel's
    self-triggering is measured while it is actually being pulsed.
    """
    p = tab[(tab.volt == volt) & tab.pulsed & tab.ok]
    return p.groupby(["cable", "length", "ch"], as_index=False).first()


# ======================================================================
# fe55
# ======================================================================
def _fe55_files(stem):
    fs = []
    for d in FE55_DIRS:
        fs = glob.glob(os.path.join(d, f"{stem}_*.pcapng"))
        if fs:
            break
    return sorted(fs, key=lambda p: int(re.search(r"_(\d+)\.pcapng$", p).group(1)))


def _fe55_file(path):
    """ADC histograms of VMM 4-13 for one capture: over-threshold hits, and
    neighbour-trigger hits (not over threshold = baseline + noise + crosstalk)."""
    h = np.zeros((len(FE55_VMMS), 64, 1024), np.int32)
    h_nb = np.zeros_like(h)
    for ch in iter_chunks(path):
        v = ch["vmm"].astype(int) - FE55_VMMS[0]
        c = ch["ch"].astype(int)
        a = np.clip(ch["adc"].astype(int), 0, 1023)
        ok = (v >= 0) & (v < len(FE55_VMMS))
        ot = ch["over_threshold"]
        np.add.at(h, (v[ok & ot], c[ok & ot], a[ok & ot]), 1)
        np.add.at(h_nb, (v[ok & ~ot], c[ok & ~ot], a[ok & ~ot]), 1)
    t0, t1, _, _ = capture_window(path)
    return h, h_nb, t0, t1


def _peaks(h, nmin):
    return np.array([[windowed_median_peak(h[v, c].astype(float))
                      if h[v, c, 60:1015].sum() > nmin else np.nan
                      for c in range(64)] for v in range(h.shape[0])])


def _fe55_run(stem, kind):
    out = os.path.join(CACHE, f"{stem}.npz")
    if os.path.exists(out) and "h_nb" in np.load(out).files:
        return out
    hsum, nbsum, per_file_pk, t_start, t_live = None, None, [], [], 0.0
    for f in _fe55_files(stem):
        h, h_nb, t0, t1 = _fe55_file(f)
        hsum = h if hsum is None else hsum + h
        nbsum = h_nb if nbsum is None else nbsum + h_nb
        t_start.append(t0)
        t_live += t1 - t0
        if kind == "fe55":
            per_file_pk.append(_peaks(h, nmin=80))
    np.savez_compressed(out, h=hsum, h_nb=nbsum, t_start=np.array(t_start), live_s=t_live,
                        file_peaks=np.array(per_file_pk) if per_file_pk else np.zeros(0))
    return out


def fe55_cache(refresh=False, jobs=6):
    if refresh:
        for _, _, _, _, _, stem in FE55_RUNS:
            p = os.path.join(CACHE, f"{stem}.npz")
            if os.path.exists(p):
                os.remove(p)
    with ProcessPoolExecutor(jobs) as ex:
        list(ex.map(_fe55_run, [r[5] for r in FE55_RUNS], [r[1] for r in FE55_RUNS]))
    data = {}
    for label, kind, cable, length, gain, stem in FE55_RUNS:
        z = np.load(os.path.join(CACHE, f"{stem}.npz"))
        data[stem] = dict(label=label, kind=kind, cable=cable, length=length, gain=gain,
                          h=z["h"], h_nb=z["h_nb"], live_s=float(z["live_s"]), t_start=z["t_start"],
                          file_peaks=z["file_peaks"])
        if kind == "fe55":
            data[stem]["peaks"] = _peaks(z["h"], nmin=300)
    return data


def _run(data, kind, cable, length, gain, which=0):
    hits = [s for s, d in data.items() if d["kind"] == kind and d["cable"] == cable
            and d["length"] == length and d["gain"] == gain]
    return hits[which] if len(hits) > which else None


def ratio_stats(a, b):
    """Per-channel photopeak ratio a/b: mid-connector median+spread, and edge/mid-connector per channel."""
    r = a / b
    mid = r[:, MID].ravel()
    mid = mid[np.isfinite(mid)]
    med = np.median(mid)
    spread = 1.4826 * np.median(abs(mid - med))
    per_ch = np.nanmedian(r, axis=0) / med
    n_vmm = np.isfinite(r).sum(axis=0)
    # uncertainty on the per-channel median across VMMs: MAD / sqrt(n)
    err = 1.4826 * np.nanmedian(abs(r / med - per_ch), axis=0) / np.sqrt(np.maximum(n_vmm, 1))
    return med, spread, per_ch, err, n_vmm


# ======================================================================
# figures
# ======================================================================
def _edge_bands(ax):
    for c in EDGE:
        ax.axvspan(c - 0.5, c + 0.5, color=EDGE_SHADE, zorder=0, lw=0)
    ax.set_xlim(-1, 64)
    ax.set_xticks([0, 8, 16, 24, 31.5, 40, 48, 56, 63],
                  ["0", "8", "16", "24", "31|32", "40", "48", "56", "63"])
    ax.tick_params(axis="x", labelsize=9)


def _save(fig, name):
    p = os.path.join(FIGS, name)
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print("wrote", p)


def fig_pulser_amplitude(pc, volt):
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    for ax, cable in zip(axs, ("Samtec", "Hitachi")):
        sub = pc[pc.cable == cable]
        for length in sorted(sub.length.unique()):
            s = sub[sub.length == length].sort_values("ch")
            mid = np.nanmedian(s.set_index("ch").amp.reindex(MID))
            lab = f"{length:.1f} m" + (" (cable assumed)" if (cable, length) in ASSUMED_CABLE else "")
            ax.plot(s.ch, s.amp / mid, "-o", ms=3.5, lw=1.5, color=LEN_COLOR[length], label=lab)
        _edge_bands(ax)
        ax.set_title(f"{cable} transmission cable", color=INK, fontsize=12, loc="left")
        ax.set_xlabel("VMM channel")
    axs[0].set_ylabel("in-time pulse amplitude /\nmid-connector median")
    axs[0].legend(title="total cable", fontsize=9, title_fontsize=9, loc="lower center", ncol=2)
    axs[1].legend(title="total cable", fontsize=9, title_fontsize=9, loc="lower center")
    axs[0].set_ylim(0.6, 1.15)
    fig.suptitle(f"Pulser, {volt}: the edge dip is the same at every length "
                 "(grey bands = connector edges)", color=INK, fontsize=12, x=0.01, ha="left")
    _save(fig, f"pulser_amplitude_vs_channel_{volt}.png")


def fig_pulser_noise(pc, volt):
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    for ax, cable in zip(axs, ("Samtec", "Hitachi")):
        sub = pc[pc.cable == cable]
        for length in sorted(sub.length.unique()):
            s = sub[sub.length == length].sort_values("ch")
            lab = f"{length:.1f} m" + (" (cable assumed)" if (cable, length) in ASSUMED_CABLE else "")
            ax.plot(s.ch, np.maximum(s.noise_hz, 1), "-o", ms=3.5, lw=1.5,
                    color=LEN_COLOR[length], label=lab)
        _edge_bands(ax)
        ax.set_yscale("log")
        ax.set_title(f"{cable} transmission cable", color=INK, fontsize=12, loc="left")
        ax.set_xlabel("VMM channel")
    axs[0].set_ylabel("out-of-time (self-trigger) rate [Hz]")
    axs[0].set_ylim(0.8, 1e5)
    axs[0].legend(title="total cable", fontsize=9, title_fontsize=9, loc="upper center", ncol=2)
    axs[1].legend(title="total cable", fontsize=9, title_fontsize=9, loc="upper center")
    fig.suptitle(f"Pulser, {volt}: edge channels self-trigger once the cable is long "
                 "(floor at 1 Hz = none seen)", color=INK, fontsize=12, x=0.01, ha="left")
    _save(fig, f"pulser_noise_vs_channel_{volt}.png")


def fig_pulser_efficiency(pc, volt):
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for cable, mk in (("Samtec", "o"), ("Hitachi", "s")):
        sub = pc[pc.cable == cable]
        for length in sorted(sub.length.unique()):
            s = sub[sub.length == length].sort_values("ch")
            ax.plot(s.ch, s.eff, "-" + mk, ms=3.5, lw=1.2, color=LEN_COLOR[length],
                    mfc="white" if cable == "Hitachi" else LEN_COLOR[length],
                    label=f"{cable} {length:.1f} m")
    _edge_bands(ax)
    ax.set_ylim(0.6, 1.05)
    ax.set_xlabel("VMM channel")
    ax.set_ylabel("in-time hits per pulse")
    ax.legend(fontsize=8, ncol=2, loc="lower center")
    ax.set_title(f"Pulser, {volt}: pulses lost on the channels that self-trigger",
                 color=INK, fontsize=12, loc="left")
    _save(fig, f"pulser_efficiency_vs_channel_{volt}.png")


def fig_fe55_edge_ratio(data, gain=1):
    """Top: raw per-channel photopeak ratio (the cable loss is visible).
    Bottom: the same divided by the mid-connector ratio (is any channel worse than typical?)."""
    s15, s20, s25 = (_run(data, "fe55", "Samtec", L, gain) for L in (1.5, 2.0, 2.5))
    hA, hB = _run(data, "fe55", "Hitachi", 2.0, gain, 0), _run(data, "fe55", "Hitachi", 2.0, gain, 1)
    pairs = [(s20, s15, "Samtec 2.0 / 1.5 m", LEN_COLOR[2.0]),
             (s25, s15, "Samtec 2.5 / 1.5 m", LEN_COLOR[3.0]),
             (hB, hA, "Hitachi 2.0 m, run B / run A (same cable)", MUTED)]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7.6), sharex=True,
                                   gridspec_kw=dict(height_ratios=(1.15, 1)))
    for i, (a, b, lab, col) in enumerate(pairs):
        if a is None or b is None:
            continue
        r = data[a]["peaks"] / data[b]["peaks"]
        med_all = np.nanmedian(r)
        med, spread, per_ch, err, _ = ratio_stats(data[a]["peaks"], data[b]["peaks"])
        x = np.arange(64) + (i - 1) * 0.22
        ax1.errorbar(x, np.nanmedian(r, axis=0), yerr=err * med, fmt="o", ms=3.5, lw=1, color=col,
                     label=f"{lab}: mid-connector ratio {med:.3f} (all channels {med_all:.3f})")
        ax1.axhline(med, color=col, lw=1, ls="--")
        ax2.errorbar(x, per_ch, yerr=err, fmt="o", ms=3.5, lw=1, color=col)
    for ax in (ax1, ax2):
        _edge_bands(ax)
        ax.axhline(1, color="#c3c2b7", lw=1)
    ax1.set_ylim(0.76, 1.04)
    ax1.set_ylabel("photopeak ratio\n(median over VMM 4-13)")
    ax1.legend(fontsize=8.5, loc="lower center")
    ax1.set_title(f"Fe55 on the detector, {gain} mV/fC. Top: raw ratio, the cable loss "
                  "(dashed = mid-connector ratio, ch 3-27 + 36-60).\nBottom: divided by the mid-connector ratio",
                  color=INK, fontsize=11, loc="left")
    ax2.set_ylim(0.95, 1.05)
    ax2.set_xlabel("VMM channel")
    ax2.set_ylabel("ratio / mid-connector ratio")
    _save(fig, f"fe55_edge_ratio_vs_channel_{gain}mVfC.png")


def fig_noise_compare(data, volt="3V3"):
    """Neighbour-trigger hits, external pulser vs detector Fe55 runs (Samtec): level (median ADC) on top,
    spread (half 16-84 % width) below; per channel and vs length."""
    nb = pulser_neighbour_hist()
    edge = list(EDGE_ZONE)
    lengths_b = [L for L in (1.0, 1.5, 2.0, 2.5, 3.0) if ("Samtec", L, volt) in nb]
    det = {L: _run(data, "fe55", "Samtec", L, 3) for L in (1.5, 2.0, 2.5)}
    fig, axs = plt.subplots(2, 3, figsize=(20, 9), sharex="col",
                            gridspec_kw=dict(width_ratios=(1.25, 1.25, 1)))
    for row, (idx, ylab, ylim) in enumerate(((0, "median ADC of\nneighbour-trigger hits", (0, 130)),
                                              (1, "spread of neighbour-trigger hits\n[ADC, half 16-84 % width]", (0, 50)))):
        ax = axs[row, 0]
        for L in lengths_b:
            v = hist_level_width(nb[("Samtec", L, volt)])[idx]
            ax.plot(np.arange(64), v, "-o", ms=3, lw=1.3, color=LEN_COLOR[L], label=f"{L:.1f} m")
        ax = axs[row, 1]
        for L, s in det.items():
            if s is None:
                continue
            v = np.nanmedian(hist_level_width(data[s]["h_nb"])[idx], axis=0)
            ax.plot(np.arange(64), v, "-s", ms=3, lw=1.3, color=LEN_COLOR[L], label=f"{L:.1f} m")
        for ax in axs[row, :2]:
            _edge_bands(ax)
            ax.set_ylim(*ylim)
            ax.legend(title="total cable", fontsize=9, title_fontsize=9, loc="upper center", ncol=3)
        axs[row, 0].set_ylabel(ylab)
        # vs length
        ax = axs[row, 2]
        for v_, ls in ((volt, "-"), ("1V8", "--")) if volt != "1V8" else ((volt, "-"),):
            xb, bb, eb = [], [], []
            for L in lengths_b:
                if ("Samtec", L, v_) not in nb:
                    continue
                q = hist_level_width(nb[("Samtec", L, v_)])[idx]
                xb.append(L); bb.append(np.nanmedian(q[MID]))
                eb.append(np.nanmedian(q[edge]) if idx == 0 else np.nanmax(q[edge]))
            ax.plot(xb, bb, ls + "o", ms=7, lw=2, color=CABLE_COLOR["Samtec"], label=f"external pulser {v_}, mid-connector")
            ax.plot(xb, eb, ls + "o", ms=7, lw=1.2, color=CABLE_COLOR["Samtec"], mfc="white",
                    label=f"external pulser {v_}, edge " + ("median" if idx == 0 else "worst"))
            if idx == 0:
                mb = [np.nanmedian(hist_mean(nb[("Samtec", L, v_)])[MID]) for L in xb]
                ax.plot(xb, mb, ls + "^", ms=6, lw=1, color=MUTED, label=f"external pulser {v_}, mid-connector mean")
        for gain, mk in ((3, "s"), (1, "D")):
            xd, bd, ed, md = [], [], [], []
            for L in (1.5, 2.0, 2.5):
                s = _run(data, "fe55", "Samtec", L, gain)
                if s is None:
                    continue
                q = hist_level_width(data[s]["h_nb"])[idx]
                qc = np.nanmedian(q, axis=0)
                xd.append(L); bd.append(np.nanmedian(q[:, MID]))
                ed.append(np.nanmedian(qc[edge]) if idx == 0 else np.nanmax(qc[edge]))
                md.append(np.nanmedian(hist_mean(data[s]["h_nb"])[:, MID]))
            ax.plot(xd, bd, "-" + mk, ms=7, lw=1.5, color=INK2, label=f"detector Fe55 {gain} mV/fC, mid-connector")
            ax.plot(xd, ed, ":" + mk, ms=7, lw=1.5, color=INK2, mfc="white",
                    label=f"detector Fe55 {gain} mV/fC, edge " + ("median" if idx == 0 else "worst"))
            if idx == 0:
                ax.plot(xd, md, ":^", ms=6, lw=1, color=MUTED, mfc="white",
                        label=f"detector Fe55 {gain} mV/fC, mid-connector mean")
        ax.set_ylim(*ylim)
        ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0,
                  handlelength=2.2, handletextpad=0.5)
    axs[0, 0].set_title(f"(a) external pulser runs, Samtec, {volt} (3 mV/fC)", color=INK, fontsize=12, loc="left")
    axs[0, 1].set_title("(b) detector, Samtec, Fe55 runs 3 mV/fC", color=INK, fontsize=12, loc="left")
    axs[0, 2].set_title("(c) vs length", color=INK, fontsize=12, loc="left")
    for ax in axs[1, :2]:
        ax.set_xlabel("VMM channel")
    axs[1, 2].set_xlabel("total Samtec cable length [m]")
    fig.suptitle("Neighbour-trigger hits (baseline + noise + crosstalk). Top: level. Bottom: spread. "
                 "External pulser runs 3-5x noisier than the detector, worse with length", color=INK, fontsize=13,
                 x=0.01, ha="left")
    fig.tight_layout()
    _save(fig, f"noise_pulser_vs_detector_{volt}.png")


def source_vs_nosource(data):
    """Detector neighbour-trigger hits, Fe55 run vs no-source reference run, Samtec.
    Hits of all mid-connector channels of VMM 4-13 pooled into one ADC histogram per
    run (the reference runs have too few hits per channel for a per-channel width)."""
    out = []
    for gain in (1, 3):
        for L in (1.5, 2.0, 2.5):
            for kind in ("fe55", "ref"):
                s = _run(data, kind, "Samtec", L, gain)
                if s is None:
                    continue
                h = data[s]["h_nb"][:, MID]
                pooled = h.sum((0, 1))
                med, wid, n = (a[0] for a in hist_level_width(pooled[None], nmin=1))
                out.append(dict(gain=gain, L=L, kind=kind, stem=s, live_s=data[s]["live_s"], hist=pooled,
                                median=float(med), width=float(wid), mean=float(hist_mean(pooled[None], 1)[0]),
                                n=int(n), n_per_ch=float(np.median(h.sum(-1)))))
    return out


def fig_source_vs_nosource(sv):
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.8), sharey=True)
    for ax, gain in zip(axs, (1, 3)):
        for r in (r for r in sv if r["gain"] == gain):
            x = np.arange(r["hist"].size)
            y = r["hist"] / max(r["n"], 1)
            src = r["kind"] == "fe55"
            ax.step(x, y, where="mid", lw=1.8 if src else 1.2, ls="-" if src else "--", color=LEN_COLOR[r["L"]],
                    label=f"{r['L']:.1f} m, {'Fe55' if src else 'no source'}: median {r['median']:.0f}, "
                          f"{r['n_per_ch']:.0f} hits/ch")
        ax.set_yscale("log")
        ax.set_xlim(0, 200)
        ax.set_ylim(1e-5, 1)
        ax.set_xlabel("ADC of neighbour-trigger hits (mid-connector ch, VMM 4-13 pooled)")
        ax.set_title(f"({'ab'[gain == 3]}) {gain} mV/fC", color=INK, fontsize=12, loc="left")
        ax.legend(fontsize=8, loc="upper right")
    axs[0].set_ylabel("fraction of hits")
    fig.suptitle("Detector, Samtec: neighbour-trigger hits with the Fe55 source (solid) and without (dashed). "
                 "Same median; the no-source runs have 150-400x fewer hits", color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout()
    _save(fig, "detector_neighbour_hits_source_vs_nosource.png")


def fig_fe55_noise(data):
    fig, ax = plt.subplots(figsize=(10, 4.4))
    refs = [s for s, d in data.items() if d["kind"] == "ref" and d["gain"] == 1]
    for s in refs:
        d = data[s]
        rate = d["h"].sum(axis=2) / d["live_s"]
        col = CABLE_COLOR[d["cable"]] if d["cable"] in CABLE_COLOR else MUTED
        ax.plot(np.arange(64), np.median(rate, axis=0), "-o", ms=3, lw=1.2, color=col,
                alpha=0.4 + 0.2 * (d["length"] - 1.5) / 0.5 if d["cable"] == "Samtec" else 0.9,
                label=d["label"])
    _edge_bands(ax)
    ax.set_yscale("log")
    ax.set_ylim(0.01, 1e5)
    ax.set_xlabel("VMM channel")
    ax.set_ylabel("over-threshold rate, no source [Hz]\n(median over VMM 4-13)")
    ax.legend(fontsize=8, ncol=2, loc="upper center")
    ax.set_title("Fe55 setup, reference runs (1 mV/fC): edge channels as quiet as mid-connector ones; "
                 "same y-range as the pulser plot", color=INK, fontsize=11, loc="left")
    _save(fig, "fe55_noise_vs_channel_1mVfC.png")


def fig_mid_vs_length(pc3, data):
    fig, ax = plt.subplots(figsize=(8.5, 5))
    # pulser: mid-connector in-time amplitude, indexed to Samtec 1.5 m
    mid = (pc3[pc3.ch.isin(MID)].groupby(["cable", "length"]).amp.median())
    ref = mid[("Samtec", 1.5)]
    for cable in ("Samtec", "Hitachi"):
        b = mid[cable]
        ax.plot(b.index, b.values / ref, "-o", ms=8, lw=2, color=CABLE_COLOR[cable],
                label=f"pulser, {cable}")
    # fe55: mid-connector photopeak ratio, indexed to Samtec 1.5 m (1 mV/fC)
    s15 = _run(data, "fe55", "Samtec", 1.5, 1)
    pts = {}
    for s, d in data.items():
        if d["kind"] != "fe55" or d["gain"] != 1:
            continue
        med, *_ = ratio_stats(d["peaks"], data[s15]["peaks"])
        pts.setdefault(d["cable"], []).append((d["length"], med, d["label"]))
    for cable, p in pts.items():
        p = sorted(p)
        x = np.array([q[0] for q in p], float)
        y = np.array([q[1] for q in p])
        if cable == "Hitachi":   # two takes of the same cable: nudge apart, label
            x = x + np.array([-0.03, 0.03])[:len(x)]
            for xi, yi, q in zip(x, y, p):
                ax.annotate(q[2][-3:], (xi, yi), textcoords="offset points",
                            xytext=(-14 if xi < 2 else 6, -3), fontsize=8, color=INK2)
        ax.errorbar(x, y, yerr=0.04, fmt="s--" if cable == "Samtec" else "s",
                    ms=8, lw=1.5, capsize=3, color=CABLE_COLOR[cable], mfc="white",
                    label=f"Fe55, {cable} (±4 % run-to-run)")
    ax.axhline(1, color="#c3c2b7", lw=1)
    ax.set_xlabel("total cable length [m]")
    ax.set_ylabel("mid-connector signal / Samtec 1.5 m")
    ax.legend(fontsize=9, loc="lower left")
    ax.set_title("Mid-connector signal vs cable length, both campaigns", color=INK, fontsize=12, loc="left")
    ax.annotate("pulser 3.0 m: cable type\nnot in logbook (assumed Samtec)",
                xy=(3.0, mid[("Samtec", 3.0)] / ref), xytext=(2.6, 0.93),
                fontsize=8, color=INK2, arrowprops=dict(arrowstyle="-", color=MUTED))
    _save(fig, "mid_connector_signal_vs_length.png")


def fig_gain_stability(data):
    s15 = _run(data, "fe55", "Samtec", 1.5, 1)
    ref = data[s15]["peaks"]
    fig, ax = plt.subplots(figsize=(10, 4.2))
    for s, d in data.items():
        if d["kind"] != "fe55" or d["gain"] != 1 or not len(d["file_peaks"]):
            continue
        y = [np.nanmedian(fp / ref) for fp in d["file_peaks"]]
        x = [dt.datetime.fromtimestamp(t) for t in d["t_start"]]
        col = LEN_COLOR[d["length"] + (0.5 if d["cable"] == "Samtec" else 0)] \
            if d["cable"] == "Samtec" else CABLE_COLOR["Hitachi"]
        ax.plot(x, y, "o", ms=5, color=col)
        ax.annotate(d["label"], (x[len(x) // 2], max(y)), textcoords="offset points",
                    xytext=(0, 8), ha="center", fontsize=9, color=INK2)
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
    ax.set_ylabel("photopeak / Samtec 1.5 m\n(median over channels, per 1-min file)")
    ax.set_xlabel("7 Oct 2026")
    ax.set_ylim(0.82, 1.13)
    ax.set_title("Fe55 1 mV/fC: flat within a run, steps between runs "
                 "(Hitachi A→B: same cable, -3.6 %)", color=INK, fontsize=12, loc="left")
    _save(fig, "fe55_gain_stability.png")


# ======================================================================
# Samtec only: the cable that will actually be used (Hitachi unavailable)
# ======================================================================
EDGE_ZONE = (0, 1, 29, 30, 31, 32, 33, 34, 62, 63)   # where self-triggering shows up


def samtec_numbers(tab, data):
    """Everything the Samtec conclusions rest on, in one dict."""
    out = {"pulser": {}, "fe55": {}}
    for volt in VOLTS:
        pc = pulser_per_channel(tab, volt)
        pc = pc[pc.cable == "Samtec"]
        rows = {}
        for length, s in pc.groupby("length"):
            s = s.set_index("ch")
            rows[length] = dict(
                mid=np.nanmedian(s.amp.reindex(MID)),
                edge_amp={c: s.amp.get(c, np.nan) / np.nanmedian(s.amp.reindex(MID)) for c in EDGE},
                noise_max=np.nanmax(s.noise_hz.reindex(EDGE_ZONE)),
                # eff == 0 is a channel with no hits at all (dead / not cabled), not a loss
                eff_min=np.nanmin(s.eff.reindex(EDGE_ZONE).replace(0, np.nan)),
                dead=[c for c in range(64) if s.n_hits.get(c, 0) == 0])
        L = np.array(sorted(rows))
        rel = np.array([rows[x]["mid"] for x in L]) / rows[1.5]["mid"]
        slope = np.polyfit(L, rel, 1)[0]
        # edge drift: edge amp/mid-connector at each length relative to the shortest cable
        drift = {x: np.nanmedian([rows[x]["edge_amp"][c] / rows[L[0]]["edge_amp"][c]
                                  for c in EDGE]) for x in L}
        out["pulser"][volt] = dict(rows=rows, L=L, rel=rel, slope=slope, edge_drift=drift)
    for gain in (1, 3):
        s15 = _run(data, "fe55", "Samtec", 1.5, gain)
        pts = {1.5: (1.0, np.ones(len(EDGE)))}
        for L in (2.0, 2.5):
            s = _run(data, "fe55", "Samtec", L, gain)
            if s:
                med, _, per_ch, _, _ = ratio_stats(data[s]["peaks"], data[s15]["peaks"])
                pts[L] = (med, per_ch[list(EDGE)])
        L = np.array(sorted(pts))
        rel = np.array([pts[x][0] for x in L])
        refs = {d["length"]: np.nanmax(np.median(d["h"].sum(axis=2) / d["live_s"], axis=0)[list(EDGE_ZONE)])
                for d in data.values() if d["kind"] == "ref" and d["cable"] == "Samtec" and d["gain"] == gain}
        out["fe55"][gain] = dict(L=L, rel=rel, slope=np.polyfit(L, rel, 1)[0] if len(L) > 1 else np.nan,
                                 edge=dict((x, pts[x][1]) for x in L), noise_max=refs)
    return out


def fig_samtec_summary(sn):
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    # (a) mid-connector signal vs length
    ax = axs[0]
    p = sn["pulser"]["3V3"]
    ax.plot(p["L"], p["rel"], "-o", ms=8, lw=2, color=CABLE_COLOR["Samtec"],
            label=f"pulser 3.3 V: {100 * p['slope']:+.1f} %/m")
    for gain, mk in ((1, "s"), (3, "D")):
        f = sn["fe55"][gain]
        ax.errorbar(f["L"], f["rel"], yerr=np.where(f["L"] == 1.5, 0, 0.04), fmt=mk + "--", ms=8,
                    lw=1.5, capsize=3, color=CABLE_COLOR["Samtec"], mfc="white",
                    label=f"Fe55 {gain} mV/fC: {100 * f['slope']:+.1f} %/m")
    ax.axhline(1, color="#c3c2b7", lw=1)
    ax.set_xlabel("total Samtec cable length [m]")
    ax.set_ylabel("mid-connector signal / 1.5 m")
    ax.set_title("(a) mid-connector channels: signal vs length", color=INK, fontsize=12, loc="left")
    ax.legend(fontsize=9, loc="lower left")
    # (b) edge channels relative to mid-connector, relative to the shortest cable
    ax = axs[1]
    d = p["edge_drift"]
    ax.plot(list(d), list(d.values()), "-o", ms=8, lw=2, color=CABLE_COLOR["Samtec"],
            label="pulser 3.3 V (vs 1.0 m)")
    f = sn["fe55"][1]
    ax.plot(f["L"], [np.nanmedian(f["edge"][x]) for x in f["L"]], "s--", ms=8, lw=1.5,
            color=CABLE_COLOR["Samtec"], mfc="white", label="Fe55 1 mV/fC (vs 1.5 m)")
    ax.axhline(1, color="#c3c2b7", lw=1)
    ax.set_ylim(0.9, 1.1)
    ax.set_xlabel("total Samtec cable length [m]")
    ax.set_ylabel("edge-channel signal / mid-connector,\nrelative to the shortest cable")
    ax.set_title("(b) edge ch 0/31/32/63 vs mid-connector", color=INK, fontsize=12, loc="left")
    ax.legend(fontsize=9, loc="upper left")
    # (c) self-triggering on the edge channels
    ax = axs[2]
    for volt, col in (("3V3", CABLE_COLOR["Samtec"]), ("1V8", LEN_COLOR[3.0])):
        rows = sn["pulser"][volt]["rows"]
        L = sorted(rows)
        y = np.array([rows[x]["noise_max"] for x in L])
        seen = y > 0
        ax.plot(np.array(L)[seen], y[seen], "-o", ms=8, lw=2, color=col, label=f"external pulser, {volt}")
        ax.plot(np.array(L)[~seen], np.full((~seen).sum(), 0.1), "v", ms=8, color=col, mfc="white",
                label="pulser: none in 10 s" if volt == "3V3" else None)
    nz = sn["fe55"][1]["noise_max"]
    ax.plot(sorted(nz), [nz[x] for x in sorted(nz)], "s--", ms=8, lw=1.5, color=CABLE_COLOR["Samtec"],
            mfc="white", label="on detector (Fe55 ref runs)")
    ax.set_yscale("log")
    ax.set_ylim(0.05, 3e4)
    ax.set_xlim(0.9, 3.4)
    ax.set_xlabel("total Samtec cable length [m]")
    ax.set_ylabel("highest self-trigger rate among\nedge channels 0/1, 29-34, 62/63 [Hz]")
    ax.set_title("(c) edge self-triggering: external pulser only", color=INK, fontsize=12, loc="left")
    ax.legend(fontsize=9, loc="upper left")
    for x in sorted(sn["pulser"]["1V8"]["rows"]):
        e = sn["pulser"]["1V8"]["rows"][x]["eff_min"]
        if e < 0.99 and x > 1.5:
            ax.annotate(f"1.8 V: {e:.2f} of pulses kept", (x, sn["pulser"]["1V8"]["rows"][x]["noise_max"]),
                        textcoords="offset points", xytext=(8, -14), ha="left", fontsize=8, color=INK2)
    fig.suptitle("Samtec cables (the ones P2 will use): what changes with length",
                 color=INK, fontsize=13, x=0.01, ha="left")
    fig.tight_layout()
    _save(fig, "samtec_summary.png")


# ======================================================================
# summary
# ======================================================================
def summary(tab, data):
    L = []
    w = L.append
    w("Cable-length study: external-pulser runs (June) vs Fe55 on detector (Oct)")
    w("=" * 68)
    w("\nPulser group labels -> pulsed VMM residue (ch % 3), per run family:")
    g = tab[tab.ok].groupby(["cable", "length", "group"]).residue.agg(lambda s: sorted(set(s)))
    for k, v in g.items():
        w(f"  {k[0]:7s} {k[1]:.1f} m  {k[2]:5s} -> {v}")
    bad = tab[~tab.ok].drop_duplicates("run")
    w("\nPulser runs dropped (no pulses recorded):")
    for _, r in bad.iterrows():
        w(f"  {r.run}  ({r.n_pulses} pulses)")
    for volt in VOLTS:
        pc = pulser_per_channel(tab, volt)
        w(f"\nPulser {volt}: mid-connector in-time ADC | edge amp/mid-connector | self-trigger rate [Hz] on 0,1,30,31,62,63")
        for (cable, length), s in pc.groupby(["cable", "length"]):
            s = s.set_index("ch")
            b = np.nanmedian(s.amp.reindex(MID))
            amp = " ".join(f"{c}:{s.amp.get(c, np.nan) / b:.2f}" for c in EDGE)
            nz = " ".join(f"{s.noise_hz.get(c, np.nan):.0f}" for c in (0, 1, 30, 31, 62, 63))
            eff = " ".join(f"{s.eff.get(c, np.nan):.2f}" for c in (0, 1, 62, 63))
            w(f"  {cable:7s} {length:.1f} m  mid-connector {b:5.0f} | {amp} | noise {nz} | eff 0/1/62/63 {eff}")
    for gain in (1, 3):
        w(f"\nFe55 {gain} mV/fC photopeak ratios (mid-connector median ± per-channel spread; edge ch / mid-connector):")
        s15 = _run(data, "fe55", "Samtec", 1.5, gain)
        for s, d in data.items():
            if d["kind"] != "fe55" or d["gain"] != gain or s == s15:
                continue
            med, spread, per_ch, err, n = ratio_stats(d["peaks"], data[s15]["peaks"])
            e = " ".join(f"{c}:{per_ch[c]:.3f}" for c in EDGE)
            w(f"  {d['label']:18s} / Samtec 1.5 m: {med:.3f} ± {spread:.3f} | {e}")
        hA, hB = _run(data, "fe55", "Hitachi", 2.0, gain, 0), _run(data, "fe55", "Hitachi", 2.0, gain, 1)
        if hA and hB:
            med, *_ = ratio_stats(data[hB]["peaks"], data[hA]["peaks"])
            w(f"  repeatability Hitachi 2.0 m B/A: {med:.3f}")
        h = data[s15]["h"]
        w(f"  saturated fraction (ADC>=1020), Samtec 1.5 m: {h[..., 1020:].sum() / h.sum():.3f}")
    w("\nFe55 reference runs, over-threshold rate [Hz], median over VMMs (mid-connector | 0 31 32 63):")
    for s, d in data.items():
        if d["kind"] != "ref":
            continue
        r = d["h"].sum(axis=2) / d["live_s"]
        w(f"  {d['label']:20s} {d['gain']} mV/fC  {np.median(r[:, MID]):5.2f} | "
          + " ".join(f"{np.median(r[:, c]):5.2f}" for c in EDGE))
    sn = samtec_numbers(tab, data)
    w("\nSAMTEC ONLY (Hitachi will not be available)")
    for volt, p in sn["pulser"].items():
        w(f"  pulser {volt}: mid-connector / 1.5 m " + " ".join(f"{x:.1f}m:{r:.3f}" for x, r in zip(p["L"], p["rel"]))
          + f"  -> {100 * p['slope']:+.1f} %/m")
        w("     edge amp/mid-connector vs 1.0 m: " + " ".join(f"{x:.1f}m:{v:.3f}" for x, v in p["edge_drift"].items()))
        w("     max self-trigger on edge zone [Hz]: " + " ".join(
            f"{x:.1f}m:{r['noise_max']:.0f}" for x, r in p["rows"].items())
          + " | min eff: " + " ".join(f"{x:.1f}m:{r['eff_min']:.2f}" for x, r in p["rows"].items())
          + " | dead ch: " + " ".join(f"{x:.1f}m:{r['dead']}" for x, r in p["rows"].items() if r["dead"]))
    for gain, f in sn["fe55"].items():
        w(f"  Fe55 {gain} mV/fC: mid-connector / 1.5 m " + " ".join(f"{x:.1f}m:{r:.3f}" for x, r in zip(f["L"], f["rel"]))
          + f"  -> {100 * f['slope']:+.1f} %/m (±4 % per point)")
        w("     edge/mid-connector (0,31,32,63): " + " | ".join(
            f"{x:.1f}m: " + " ".join(f"{v:.3f}" for v in f["edge"][x]) for x in f["L"]))
        w("     max ref-run rate on edge zone [Hz]: " + " ".join(f"{x:.1f}m:{v:.2f}" for x, v in sorted(f["noise_max"].items())))
    fig_samtec_summary(sn)

    w("\nNEIGHBOUR-TRIGGER HITS, Samtec (level = median ADC, spread = half 16-84 % width)")
    w("  external-pulser runs: 1.0/1.5 m on 5 June, 2.0/2.5/3.0 m on 8 June (length and day are confounded)")
    nb = pulser_neighbour_hist()
    for Lb in (1.0, 1.5, 2.0, 2.5, 3.0):
        parts = []
        for volt in VOLTS:
            if ("Samtec", Lb, volt) in nb:
                m, wd, _ = hist_level_width(nb[("Samtec", Lb, volt)])
                mn = hist_mean(nb[("Samtec", Lb, volt)])
                parts.append(f"{volt}: median {np.nanmedian(m[MID]):.0f} mean {np.nanmedian(mn[MID]):.0f} "
                             f"edge median {np.nanmedian(m[list(EDGE_ZONE)]):.0f} spread {np.nanmedian(wd[MID]):.1f} "
                             f"worst edge {np.nanmax(wd[list(EDGE_ZONE)]):.0f}")
        w(f"  external pulser {Lb:.1f} m | " + " | ".join(parts))
    for gain in (1, 3):
        for Lf in (1.5, 2.0, 2.5):
            s = _run(data, "fe55", "Samtec", Lf, gain)
            if s is None:
                continue
            m, wd, _ = hist_level_width(data[s]["h_nb"])
            mn = hist_mean(data[s]["h_nb"])
            w(f"  detector Fe55 {gain} mV/fC {Lf:.1f} m | median {np.nanmedian(m[:, MID]):.0f} "
              f"mean {np.nanmedian(mn[:, MID]):.0f} edge median {np.nanmedian(np.nanmedian(m, axis=0)[list(EDGE_ZONE)]):.0f} "
              f"spread {np.nanmedian(wd[:, MID]):.1f} worst edge {np.nanmax(np.nanmedian(wd, axis=0)[list(EDGE_ZONE)]):.1f}")
    w("  detector, Fe55 run vs no-source run (mid-connector ch of VMM 4-13 pooled):")
    sv = source_vs_nosource(data)
    for r in sv:
        w(f"    {r['gain']} mV/fC {r['L']:.1f} m {'Fe55     ' if r['kind'] == 'fe55' else 'no source'} "
          f"{r['stem']:28s} live {r['live_s']:4.0f} s | hits/ch {r['n_per_ch']:6.0f} | "
          f"median {r['median']:.0f} mean {r['mean']:.1f} spread {r['width']:.1f}")
    fig_source_vs_nosource(sv)
    w("  mid-connector ratio vs all-channel ratio (Fe55, Samtec / 1.5 m):")
    for gain in (1, 3):
        s15 = _run(data, "fe55", "Samtec", 1.5, gain)
        for Lf in (2.0, 2.5):
            s = _run(data, "fe55", "Samtec", Lf, gain)
            if s:
                r = data[s]["peaks"] / data[s15]["peaks"]
                w(f"    {gain} mV/fC {Lf:.1f} m: mid-connector {ratio_stats(data[s]['peaks'], data[s15]['peaks'])[0]:.3f}"
                  f"  all channels {np.nanmedian(r):.3f}")
    txt = "\n".join(L)
    with open(os.path.join(FIGS, "summary.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--refresh", action="store_true", help="rebuild the caches")
    ap.add_argument("--jobs", type=int, default=6, help="parallel pcapng decoders")
    args = ap.parse_args()
    os.makedirs(CACHE, exist_ok=True)
    os.makedirs(FIGS, exist_ok=True)

    tab = pulser_table(args.refresh)
    data = fe55_cache(args.refresh, args.jobs)

    for volt in ("3V3", "1V8"):
        pc = pulser_per_channel(tab, volt)
        fig_pulser_amplitude(pc, volt)
        fig_pulser_noise(pc, volt)
        fig_pulser_efficiency(pc, volt)
    fig_fe55_edge_ratio(data, 1)
    fig_fe55_edge_ratio(data, 3)
    fig_noise_compare(data, "3V3")
    fig_noise_compare(data, "1V8")
    fig_fe55_noise(data)
    fig_mid_vs_length(pulser_per_channel(tab, "3V3"), data)
    fig_gain_stability(data)
    summary(tab, data)


if __name__ == "__main__":
    main()
