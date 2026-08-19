#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vmm_timing_peaks.py -- the actual dt distributions behind the timing numbers,
with Gaussian-on-flat-background fits, for the best-timing working point of
each detector and each gas.

The campaign table (vmm_trigger_timing.csv) carries only the FITTED mu and
sigma per capture.  To show the distribution itself -- which is what a fit is
worth looking at -- the coincidence has to be rebuilt from hit columns.  This
script does that for a named sub_run and stores the histograms, so the figure
can then be drawn without touching the data again.

Method, identical to vmm_efficiency.py's:
  hits are paired with a trigger only inside the SAME srs_timestamp (markers),
  BCID phases are wrapped into +-46.08 us, and each station hit takes its
  nearest trigger.  Fit = Gaussian + flat background over the search window;
  the flat term is the ~18-fold BCID ambiguity, not a systematic.

Products:  timing_peaks_<tag>.npz   histograms + fit results per station

Usage
-----
  # local store (a directory of capture directories holding *.npy columns)
  python3 vmm_timing_peaks.py --store /path/to/_subrun_store --tag run33_gap150

  # a sub_run staged from EOS
  python3 vmm_timing_peaks.py --store $TMPDIR/run_66/cfg_gain4.5_peaktime200_opt \
      --tag run66_gasB --captures 40
"""

import argparse
import glob
import json
import os

import numpy as np

CLOCK = 22.5
TAC_NS, TDC_RANGE = 60.0, 255
PERIOD = 4096 * CLOCK
STATIONS = {'P2_IN': range(2, 8), 'P2_MID': range(8, 14),
            'P2_OUT': range(14, 20)}
TRG_VMM, TRG_CH = 0, 44
MAX_PAIRS_PER_MARKER = 200_000

try:
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


def wrap(d):
    return (d + PERIOD / 2.0) % PERIOD - PERIOD / 2.0


def gauss_bg(x, a, mu, sigma, bg):
    return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + bg


def hot_channels(c, max_rate_hz):
    """Channels above an absolute rate, which is the only rule that survives.

    Flagging at N x the chip's own median occupancy deletes beam: on a chip
    whose threshold has been raised the noise disappears, the median collapses
    and the illuminated pads stand 30-50x above it (EFFICIENCY_AUTOPSY.md).
    A 12 x 12 mm pad in this beam takes a few kHz and the pathological channels
    sit at 1e5 Hz, so an absolute cut has two decades of clear air.
    """
    srs = c['srs_timestamp'].astype(np.float64)
    lo, hi = np.quantile(srs, [0.001, 0.999])     # quantiles: a few hits read 0
    live_s = max((hi - lo) * CLOCK * 1e-9, 1e-3)
    key = c['vmm'].astype(np.int64) * 64 + c['ch'].astype(np.int64)
    u, n = np.unique(key, return_counts=True)
    return set(u[(n / live_s) > max_rate_hz].tolist()), live_s


def dt_of(store, station, captures, win, trg_vmm, trg_ch, max_rate_hz=0.0):
    """dt for every station hit against its nearest trigger in the same marker."""
    caps = sorted(glob.glob(os.path.join(store, '*')))
    caps = [c for c in caps if os.path.isdir(c)][:captures]
    vmms = np.fromiter(STATIONS[station], int)
    out = []
    for cap in caps:
        try:
            c = {k: np.load(os.path.join(cap, f'{k}.npy')) for k in
                 ('vmm', 'ch', 'bcid', 'tdc', 'srs_timestamp')}
        except FileNotFoundError:
            continue
        t = (c['bcid'].astype(np.float64) * CLOCK
             + CLOCK - c['tdc'].astype(np.float64) * TAC_NS / TDC_RANGE)
        is_t = (c['vmm'] == trg_vmm) & (c['ch'] == trg_ch)
        is_d = np.isin(c['vmm'], vmms)
        if max_rate_hz > 0:
            hot, _ = hot_channels(c, max_rate_hz)
            if hot:
                key = c['vmm'].astype(np.int64) * 64 + c['ch'].astype(np.int64)
                is_d &= ~np.isin(key, list(hot))
        order = np.argsort(c['srs_timestamp'], kind='stable')
        srs = c['srs_timestamp'][order]
        bounds = np.flatnonzero(np.diff(srs)) + 1
        for lo, hi in zip(np.r_[0, bounds], np.r_[bounds, srs.size]):
            idx = order[lo:hi]
            ti, di = idx[is_t[idx]], idx[is_d[idx]]
            if ti.size == 0 or di.size == 0:
                continue
            if ti.size * di.size > MAX_PAIRS_PER_MARKER:
                continue
            d = wrap(t[di][:, None] - t[ti][None, :])
            j = np.abs(d).argmin(axis=1)
            dd = d[np.arange(d.shape[0]), j]
            k = np.abs(dd) < win
            if k.any():
                out.append(dd[k].astype(np.float32))
    return np.concatenate(out) if out else np.zeros(0, np.float32)


def fit_peak(dt, win, nbins):
    """Gaussian on a flat background.  Returns the histogram and the fit."""
    h, e = np.histogram(dt, bins=nbins, range=(-win, win))
    ctr = 0.5 * (e[:-1] + e[1:])
    if h.sum() < 500:
        return dict(counts=h, edges=e, ok=False)
    bg0 = float(np.median(h[np.abs(ctr) > win * 0.6]))
    k = int((h - bg0).argmax())
    mu0 = float(ctr[k])
    a0 = float(h[k] - bg0)
    # a narrow seed: the core within +-100 ns of the crude peak
    core = np.abs(ctr - mu0) < 100
    w = (h - bg0).clip(0)[core]
    s0 = float(np.sqrt(np.average((ctr[core] - mu0) ** 2,
                                  weights=w if w.sum() else None))) or 20.0
    p = [a0, mu0, s0, bg0]
    ok = False
    if HAVE_SCIPY and a0 > 0:
        try:
            fitr = np.abs(ctr - mu0) < max(6 * s0, 150)
            p, _ = curve_fit(gauss_bg, ctr[fitr], h[fitr], p0=p,
                             sigma=np.sqrt(np.maximum(h[fitr], 1)),
                             maxfev=20000)
            ok = True
        except (RuntimeError, ValueError):
            pass
    a, mu, sg, bg = p
    sg = abs(float(sg))
    return dict(counts=h, edges=e, ok=ok, amp=float(a), mu=float(mu),
                sigma=sg, bg=float(bg), n=int(h.sum()),
                peak_over_bg=float(a / bg) if bg > 0 else float('nan'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--store', required=True,
                    help='directory of capture directories holding *.npy columns')
    ap.add_argument('--tag', required=True, help='name for the output npz')
    ap.add_argument('--captures', type=int, default=12)
    ap.add_argument('--win', type=float, default=500.0)
    ap.add_argument('--nbins', type=int, default=200)
    ap.add_argument('--trigger-vmm', type=int, default=TRG_VMM)
    ap.add_argument('--trigger-ch', type=int, default=TRG_CH)
    ap.add_argument('--label', default='', help='free text stored with the fits')
    ap.add_argument('--max-rate-khz', type=float, default=20.0,
                    help='drop station channels above this absolute rate; '
                         '0 disables. A beam pad takes a few kHz, the '
                         'pathological channels 100+.')
    a = ap.parse_args()

    res, meta = {}, dict(tag=a.tag, store=a.store, label=a.label,
                         win=a.win, nbins=a.nbins, captures=a.captures,
                         max_rate_khz=a.max_rate_khz)
    for st in STATIONS:
        dt = dt_of(a.store, st, a.captures, a.win,
                   a.trigger_vmm, a.trigger_ch,
                   max_rate_hz=a.max_rate_khz * 1e3)
        f = fit_peak(dt, a.win, a.nbins)
        res[f'{st}_counts'] = f['counts']
        res[f'{st}_edges'] = f['edges']
        meta[st] = {k: v for k, v in f.items()
                    if k not in ('counts', 'edges')}
        if f.get('ok'):
            print(f'  {st:8s} mu {f["mu"]:7.1f} ns   sigma {f["sigma"]:6.2f} ns'
                  f'   peak/bg {f["peak_over_bg"]:6.1f}   {f["n"]} pairs')
        else:
            print(f'  {st:8s} no usable peak ({f.get("n", 0)} pairs)')

    out = f'timing_peaks_{a.tag}.npz'
    np.savez_compressed(out, meta=json.dumps(meta), **res)
    print(f'-> {out}')


if __name__ == '__main__':
    main()
