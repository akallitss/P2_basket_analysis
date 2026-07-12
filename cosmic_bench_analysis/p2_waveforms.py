#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_waveforms.py

Waveform access + time-of-arrival extraction for the P2 pad detector, working
directly on the per-FEU decoded_root files (nt tree written by
mm_dream_reconstruction/decoder: eventId, timestamp, ftst + parallel
sample/channel/amplitude vectors; non-ZS = full 32-sample waveforms for all
512 DREAM channels of the FEU).

Processing chain (mirrors mm_dream_reconstruction/waveform_analysis):
  1. pedestal mean per channel      : median of raw samples over events
                                      (signals are rare, median is immune)
  2. common-noise subtraction (CNS) : per event, per 64-channel block, per
                                      sample, subtract the block median
  3. noise RMS per channel          : robust MAD of the CNS-subtracted samples
  4. hit selection                  : peak > n_sig_peak * rms and >=
                                      min_samples_above samples > n_sig_samp*rms

Fine-timestamp correction (same convention as WaveformAnalyzer.cpp):
  ftst is the trigger phase within the coarse DREAM clock, in units of 10 ns
  (100 MHz). The corrected hit time is
      t_ns = t_samples * time_per_sample + ftst * 10
  i.e. the ftst shift is ADDED to the sample-clock time to re-reference it to
  the trigger. 13_timing_waveforms.py validates the sign empirically (mean
  t vs ftst must be flat for the correct convention).

Timing algorithms (each maps one baseline-subtracted waveform -> time in
SAMPLE units; NaN if it cannot be computed). All leading-edge crossings are
linearly interpolated between samples.
"""

import glob
import os
import re

import numpy as np
import pandas as pd

N_CH = 512          # channels per DREAM FEU
N_SAMP = 32         # samples per waveform (run_config n_samples_per_waveform)
TIME_PER_FTST = 10.0   # ns per ftst unit (DREAM 100 MHz clock)


# --------------------------------------------------------------------------- #
# decoded_root access
# --------------------------------------------------------------------------- #
_FEU_RE = re.compile(r'_(\d{3})_(\d{2})\.root$')


def list_decoded_files(decoded_dir, feus=None):
    """Decoded files as a DataFrame(path, chunk, feu), sorted by (chunk, feu).
    `feus`: optional iterable restricting to those FEU numbers."""
    rows = []
    for fp in sorted(glob.glob(os.path.join(decoded_dir, '*.root'))):
        m = _FEU_RE.search(os.path.basename(fp))
        if not m:
            continue
        chunk, feu = int(m.group(1)), int(m.group(2))
        if feus is not None and feu not in set(int(f) for f in feus):
            continue
        rows.append((fp, chunk, feu))
    return pd.DataFrame(rows, columns=['path', 'chunk', 'feu'])


def iter_decoded_events(path, step=2000):
    """Yield (eventId[n], timestamp[n], ftst[n], wf[n, 512, 32]) chunks from a
    decoded nt tree. Dense fill via index arithmetic — no ordering assumption
    on the (sample, channel, amplitude) triplets."""
    import uproot
    with uproot.open(f'{path}:nt') as t:
        for arrs in t.iterate(['eventId', 'timestamp', 'ftst', 'sample',
                               'channel', 'amplitude'],
                              step_size=step, library='ak'):
            import awkward as ak
            ev = np.asarray(arrs['eventId'])
            ts = np.asarray(arrs['timestamp'])
            ft = np.asarray(arrs['ftst'])
            n = len(ev)
            counts = ak.num(arrs['channel'])
            ch = np.asarray(ak.flatten(arrs['channel']), dtype=np.int64)
            sa = np.asarray(ak.flatten(arrs['sample']), dtype=np.int64)
            am = np.asarray(ak.flatten(arrs['amplitude']), dtype=np.float32)
            iev = np.repeat(np.arange(n, dtype=np.int64), np.asarray(counts))
            wf = np.zeros((n, N_CH, N_SAMP), dtype=np.float32)
            ok = (ch >= 0) & (ch < N_CH) & (sa >= 0) & (sa < N_SAMP)
            wf[iev[ok], ch[ok], sa[ok]] = am[ok]
            yield ev, ts, ft, wf


def pedestal_from_data(path, n_events=400):
    """Per-channel pedestal mean: median of raw samples over the first
    `n_events` events (cosmic occupancy is tiny, so the median never sits on
    a signal). Returns (512,) float32."""
    for ev, ts, ft, wf in iter_decoded_events(path, step=n_events):
        flat = wf.transpose(1, 0, 2).reshape(N_CH, -1)   # (512, n*32)
        return np.median(flat, axis=1).astype(np.float32)
    return np.zeros(N_CH, dtype=np.float32)


def subtract_common_noise(wf):
    """CNS in-place equivalent: per event, per 64-channel block, per sample,
    subtract the median over the block channels. wf: (n, 512, 32)."""
    n = wf.shape[0]
    blocks = wf.reshape(n, N_CH // 64, 64, N_SAMP)
    med = np.median(blocks, axis=2, keepdims=True)
    return (blocks - med).reshape(n, N_CH, N_SAMP)


def noise_rms_from_data(path, pedestal, n_events=400):
    """Per-channel robust noise RMS (MAD * 1.4826) of the CNS-subtracted,
    pedestal-subtracted samples."""
    for ev, ts, ft, wf in iter_decoded_events(path, step=n_events):
        w = subtract_common_noise(wf - pedestal[None, :, None])
        flat = w.transpose(1, 0, 2).reshape(N_CH, -1)
        mad = np.median(np.abs(flat - np.median(flat, axis=1, keepdims=True)),
                        axis=1)
        return (1.4826 * mad).astype(np.float32)
    return np.full(N_CH, 3.0, dtype=np.float32)


def extract_hits(path, pedestal, rms, n_sig_peak=5.0, n_sig_samp=3.0,
                 min_samples_above=3, step=2000, max_events=None,
                 adc_max=4095.0):
    """Scan a decoded file and return the hit waveforms.

    Returns (df, waves):
      df    : DataFrame(eventId, timestamp, ftst, channel, amp_max, saturated)
      waves : (n_hits, 32) float32, pedestal- and CNS-subtracted
    """
    dfs, wlist = [], []
    seen = 0
    for ev, ts, ft, wf in iter_decoded_events(path, step=step):
        w = subtract_common_noise(wf - pedestal[None, :, None])
        peak = w.max(axis=2)                                  # (n, 512)
        nabove = (w > (n_sig_samp * rms)[None, :, None]).sum(axis=2)
        hit = (peak > n_sig_peak * rms[None, :]) & (nabove >= min_samples_above)
        iev, ich = np.nonzero(hit)
        if len(iev):
            raw_max = wf[iev, ich, :].max(axis=1)
            dfs.append(pd.DataFrame({
                'eventId': ev[iev], 'timestamp': ts[iev], 'ftst': ft[iev],
                'channel': ich, 'amp_max': peak[iev, ich],
                'noise_rms': rms[ich],
                'saturated': raw_max >= adc_max - 4}))
            wlist.append(w[iev, ich, :].copy())
        seen += len(ev)
        if max_events is not None and seen >= max_events:
            break
    if not dfs:
        return (pd.DataFrame(columns=['eventId', 'timestamp', 'ftst', 'channel',
                                      'amp_max', 'noise_rms', 'saturated']),
                np.empty((0, N_SAMP), dtype=np.float32))
    return pd.concat(dfs, ignore_index=True), np.vstack(wlist)


# --------------------------------------------------------------------------- #
# per-waveform helpers (vectorized over a (n, 32) matrix)
# --------------------------------------------------------------------------- #
def _peak_parabola(w):
    """Sub-sample peak (time, amplitude) via 3-point parabola around argmax."""
    n, m = w.shape
    imax = w.argmax(axis=1)
    t = imax.astype(float)
    a = w[np.arange(n), imax]
    inner = (imax > 0) & (imax < m - 1)
    i = imax[inner]
    y1 = w[inner, i - 1]; y2 = w[inner, i]; y3 = w[inner, i + 1]
    den = y1 - 2 * y2 + y3
    off = np.zeros(inner.sum())
    good = np.abs(den) > 1e-9
    off[good] = 0.5 * (y1 - y3)[good] / den[good]
    off = np.clip(off, -1, 1)
    t[inner] = i + off
    a[inner] = y2 - 0.25 * (y1 - y3) * off
    return t, a


def _leading_edge_cross(w, level, imax):
    """Interpolated crossing time of `level` walking left from imax.
    w: (n, 32); level, imax: (n,). NaN when the edge never drops below level."""
    n, m = w.shape
    t = np.full(n, np.nan)
    for k in range(n):          # short walks (<32), plain loop is fine
        i = int(imax[k])
        lv = level[k]
        while i > 0 and w[k, i - 1] > lv:
            i -= 1
        if i == 0 and w[k, 0] > lv:
            continue            # pulse starts before the window
        y0, y1 = w[k, i - 1] if i > 0 else 0.0, w[k, i]
        if i == 0:
            t[k] = 0.0
        elif abs(y1 - y0) > 1e-9:
            t[k] = (i - 1) + (lv - y0) / (y1 - y0)
        else:
            t[k] = i - 0.5
    return t


# --------------------------------------------------------------------------- #
# timing algorithms: (waves, df) -> time in SAMPLE units (n,), NaN = failed
# --------------------------------------------------------------------------- #
def t_frac(w, frac=0.3):
    """Leading-edge crossing of frac * (parabola peak amplitude) — the
    mm_dream_reconstruction reference uses frac=0.3."""
    tp, ap = _peak_parabola(w)
    imax = w.argmax(axis=1)
    return _leading_edge_cross(w, frac * ap, imax)


def t_thr(w, noise_rms, n_sig=5.0):
    """Fixed-threshold leading edge at n_sig * channel noise RMS."""
    imax = w.argmax(axis=1)
    return _leading_edge_cross(w, n_sig * np.asarray(noise_rms), imax)


def t_parabola(w):
    """Sub-sample time of the pulse maximum (parabola)."""
    tp, _ = _peak_parabola(w)
    return tp


def t_centroid(w, n_sig=3.0, noise_rms=None):
    """Amplitude-weighted mean sample over the above-threshold pulse window."""
    thr = (n_sig * np.asarray(noise_rms))[:, None] if noise_rms is not None \
        else 0.2 * w.max(axis=1, keepdims=True)
    ww = np.where(w > thr, w, 0.0)
    s = ww.sum(axis=1)
    t = (ww * np.arange(w.shape[1])[None, :]).sum(axis=1)
    out = np.full(len(w), np.nan)
    ok = s > 0
    out[ok] = t[ok] / s[ok]
    return out


def t_dcfd(w, frac=0.5, delay=2):
    """Digital constant-fraction: zero crossing of frac*w[t] - w[t-delay],
    interpolated, taken on the leading edge (last crossing before the max)."""
    n, m = w.shape
    c = frac * w[:, delay:] - w[:, :-delay]        # c[t] for t = delay..m-1
    out = np.full(n, np.nan)
    imax = w.argmax(axis=1)
    for k in range(n):
        hi = int(imax[k]) - delay                  # index into c
        found = np.nan
        for i in range(max(hi, 0), 0, -1):
            if c[k, i] >= 0 > c[k, i - 1]:
                y0, y1 = c[k, i - 1], c[k, i]
                found = (i - 1) + (0 - y0) / (y1 - y0) + delay
                break
        out[k] = found
    return out


def t_sigmoid(w, max_fits=None):
    """Half-max time from a logistic fit  A / (1 + exp(-(t-t0)/tau))  to the
    leading edge (start of pulse .. one sample past the max)."""
    from scipy.optimize import curve_fit

    def f(t, A, t0, tau):
        return A / (1.0 + np.exp(-(t - t0) / tau))

    n, m = w.shape
    out = np.full(n, np.nan)
    imax = w.argmax(axis=1)
    t50 = t_frac(w, 0.5)
    idx = np.arange(n) if max_fits is None else np.arange(min(n, max_fits))
    tt = np.arange(m, dtype=float)
    for k in idx:
        i1 = int(imax[k])
        i0 = i1
        lv = 0.05 * w[k, i1]
        while i0 > 0 and w[k, i0 - 1] > lv:
            i0 -= 1
        lo, hi = max(i0 - 1, 0), min(i1 + 2, m)
        if hi - lo < 4 or not np.isfinite(t50[k]):
            continue
        try:
            p, _ = curve_fit(f, tt[lo:hi], w[k, lo:hi],
                             p0=[w[k, i1], t50[k], 0.7],
                             maxfev=200)
            if lo - 2 <= p[1] <= hi + 2 and 0.05 < p[2] < 10:
                out[k] = p[1]
        except Exception:
            pass
    return out


def build_template(w, amp, n_top=500, upsample=10, sat=None):
    """Average normalized pulse shape from the highest-amplitude clean
    waveforms, aligned at their 50% leading-edge crossing. Returns
    (t_up, template) on an upsampled grid (sample units)."""
    order = np.argsort(amp)[::-1]
    if sat is not None:
        order = order[~np.asarray(sat)[order]]
    sel = order[:n_top]
    t50 = t_frac(w[sel], 0.5)
    ok = np.isfinite(t50)
    sel, t50 = sel[ok], t50[ok]
    m = w.shape[1]
    t_up = np.arange(-m, m, 1.0 / upsample)
    acc = np.zeros_like(t_up)
    cnt = np.zeros_like(t_up)
    tt = np.arange(m, dtype=float)
    for k, s in enumerate(sel):
        a = w[s].max()
        if a <= 0:
            continue
        shifted = tt - t50[k]
        vals = np.interp(t_up, shifted, w[s] / a, left=np.nan, right=np.nan)
        good = np.isfinite(vals)
        acc[good] += vals[good]
        cnt[good] += 1
    tem = np.where(cnt > 0, acc / np.maximum(cnt, 1), 0.0)
    return t_up, tem


def t_template(w, t_up, template, upsample=10, search=(-8.0, 8.0)):
    """Best time shift of the normalized template by maximising the dot
    product on an upsampled grid (parabola-refined). Time returned is the
    shift of the template's 50%-crossing reference, in sample units."""
    m = w.shape[1]
    tt = np.arange(m, dtype=float)
    shifts = np.arange(search[0], search[1] + 1e-9, 1.0 / upsample)
    # matrix of template values sampled at (tt - shift): (n_shifts, m)
    tem_at = np.vstack([np.interp(tt - s, t_up, template, left=0, right=0)
                        for s in shifts])
    norms = (tem_at ** 2).sum(axis=1)
    out = np.full(len(w), np.nan)
    for k in range(len(w)):
        score = tem_at @ w[k] / np.sqrt(np.maximum(norms, 1e-9))
        j = int(score.argmax())
        if 0 < j < len(shifts) - 1:
            y1, y2, y3 = score[j - 1], score[j], score[j + 1]
            den = y1 - 2 * y2 + y3
            off = 0.5 * (y1 - y3) / den if abs(den) > 1e-12 else 0.0
            out[k] = shifts[j] + np.clip(off, -1, 1) / upsample
        else:
            out[k] = shifts[j]
    return out


ALGORITHMS = {
    'frac20':   lambda w, df: t_frac(w, 0.20),
    'frac30':   lambda w, df: t_frac(w, 0.30),   # reference implementation
    'frac50':   lambda w, df: t_frac(w, 0.50),
    'thr5sig':  lambda w, df: t_thr(w, df['noise_rms'].to_numpy(), 5.0),
    'parabola': lambda w, df: t_parabola(w),
    'centroid': lambda w, df: t_centroid(w, 3.0, df['noise_rms'].to_numpy()),
    'dcfd':     lambda w, df: t_dcfd(w, 0.5, 2),
    'sigmoid':  lambda w, df: t_sigmoid(w),
    # 'template' is added by the stage after building the template
}


# --------------------------------------------------------------------------- #
# statistics helpers
# --------------------------------------------------------------------------- #
def robust_sigma(x, clip=2.5, iters=5):
    """Gaussian-core sigma by iterative clipping (insensitive to tails)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 10:
        return np.nan
    mu, sig = np.median(x), 1.4826 * np.median(np.abs(x - np.median(x)))
    for _ in range(iters):
        m = np.abs(x - mu) < clip * max(sig, 1e-9)
        if m.sum() < 10:
            break
        mu, sig = x[m].mean(), x[m].std()
    return float(sig)


def q68_half_width(x):
    """Half-width of the central 68.3% interval (quantile-based sigma)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 10:
        return np.nan
    lo, hi = np.percentile(x, [15.865, 84.135])
    return float((hi - lo) / 2)


def slewing_correction(t_ns, amp, n_bins=25):
    """Amplitude-walk correction: median t per log10(amp) bin, interpolated
    and subtracted. Returns (t_corrected, bin_centers_amp, bin_median_t)."""
    ok = np.isfinite(t_ns) & (amp > 0)
    la = np.log10(np.where(amp > 0, amp, np.nan))
    edges = np.nanpercentile(la[ok], np.linspace(0, 100, n_bins + 1))
    edges = np.unique(edges)
    if len(edges) < 3:
        return t_ns - np.nanmedian(t_ns), np.array([]), np.array([])
    ib = np.clip(np.digitize(la, edges) - 1, 0, len(edges) - 2)
    med = np.full(len(edges) - 1, np.nan)
    for b in range(len(edges) - 1):
        m = ok & (ib == b)
        if m.sum() >= 20:
            med[b] = np.median(t_ns[m])
    ctr = 0.5 * (edges[:-1] + edges[1:])
    good = np.isfinite(med)
    if good.sum() < 2:
        return t_ns - np.nanmedian(t_ns), ctr, med
    corr = np.interp(la, ctr[good], med[good])
    return t_ns - corr, 10 ** ctr, med
