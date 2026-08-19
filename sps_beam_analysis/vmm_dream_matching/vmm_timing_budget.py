#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vmm_timing_budget.py -- where the VMM coincidence width comes from, and which
handles would actually remove it.

The campaign quotes the VMM time width as "pinned at the 22.5 ns BCID
quantisation for want of calibration runs".  That is not what the data says,
and this script is the measurement.  Three tests, all on hit-level columns, no
new runs needed:

  1. TDC fine time on/off.  vmm_decode.derive() already applies
     t_fine = 22.5 - tdc * 60/255 with the nominal TAC slope, so the time is
     NOT quantised at 22.5 ns.  Toggling it measures what it is worth.

  2. Per-channel t0 and time walk.  The VMM has no per-channel time
     calibration and ships a per-hit ADC, so both corrections are free
     offline.  This measures the spread each one removes.

  3. Station-minus-station coincidence, which contains no trigger term at all:

        sigma(MID - TRIG)^2 = sMID^2 + sTRIG^2
        sigma(MID - OUT)^2  = sMID^2 + sOUT^2

     Three such pairs over-determine the system, so the trigger channel's own
     contribution can be solved for rather than assumed.

Caveats, and they matter: the width here is an rms inside +-120 ns of the peak,
not the Gaussian core sigma that `vmm_efficiency.py` fits, so it sits above the
campaign numbers -- use it for the DECOMPOSITION, not as an absolute.  And the
only sub_run with hit columns on this machine is a low-drift-field point on the
bench gas, which is a worst case for the drift term; repeat at the nominal
working point when EOS is reachable.

Usage:  python3 vmm_timing_budget.py [--store DIR] [--captures N]
"""

import argparse
import glob
import os

import numpy as np

STORE = ('/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/'
         'SPS_Beam_Test/TB_July26_H4_vmm/_subrun_store')
CLOCK = 22.5            # ns per BCID count (44.44 MHz)
TAC_NS, TDC_RANGE = 60.0, 255
PERIOD = 4096 * CLOCK   # 92160 ns, one BCID wrap
STATIONS = {'P2_IN': range(2, 8), 'P2_MID': range(8, 14),
            'P2_OUT': range(14, 20)}
TRG_VMM, TRG_CH = 0, 44
WIN, CORE = 400.0, 120.0
MR = 0.0                # hot-channel rate veto, set from the CLI in main()
MAX_PAIRS_PER_MARKER = 200_000     # skip pathological markers, bounds memory


def wrap(d):
    return (d + PERIOD / 2.0) % PERIOD - PERIOD / 2.0


def load(cap, cols):
    return {k: np.load(os.path.join(cap, f'{k}.npy')) for k in cols}


def times(c, use_fine=True):
    t = c['bcid'].astype(np.float64) * CLOCK
    if use_fine:
        t = t + (CLOCK - c['tdc'].astype(np.float64) * TAC_NS / TDC_RANGE)
    return t


def hot_mask(c, max_rate_hz):
    """Absolute-rate hot-channel veto, the same rule the efficiency uses.

    A ratio-to-chip-median rule deletes beam on any chip whose threshold has
    been raised (EFFICIENCY_AUTOPSY.md); an absolute cut has two decades of
    clear air between a beam pad (few kHz) and a pathological channel (1e5 Hz).
    """
    if max_rate_hz <= 0:
        return None
    srs = c['srs_timestamp'].astype(np.float64)
    lo, hi = np.quantile(srs, [0.001, 0.999])
    live_s = max((hi - lo) * CLOCK * 1e-9, 1e-3)
    key = c['vmm'].astype(np.int64) * 64 + c['ch'].astype(np.int64)
    u, n = np.unique(key, return_counts=True)
    hot = u[(n / live_s) > max_rate_hz]
    return np.isin(key, hot) if hot.size else None


def picker(c, name):
    if name == 'TRIG':
        return (c['vmm'] == TRG_VMM) & (c['ch'] == TRG_CH)
    return np.isin(c['vmm'], np.fromiter(STATIONS[name], int))


def pair(caps, a, b, use_fine=True, want=(), max_rate_hz=0.0):
    """Nearest-in-time partner for every `a` hit, within one marker interval."""
    need = ['vmm', 'ch', 'bcid', 'tdc', 'srs_timestamp'] + list(want)
    out = {k: [] for k in ('dt',) + tuple(want)}
    for cap in caps:
        try:
            c = load(cap, need)
        except FileNotFoundError:
            return None
        t = times(c, use_fine)
        ia, ib = picker(c, a), picker(c, b)
        hot = hot_mask(c, max_rate_hz)
        if hot is not None:
            ia = ia & ~hot
        order = np.argsort(c['srs_timestamp'], kind='stable')
        srs = c['srs_timestamp'][order]
        bounds = np.flatnonzero(np.diff(srs)) + 1
        for lo, hi in zip(np.r_[0, bounds], np.r_[bounds, srs.size]):
            idx = order[lo:hi]
            ai, bi = idx[ia[idx]], idx[ib[idx]]
            if ai.size == 0 or bi.size == 0:
                continue
            if ai.size * bi.size > MAX_PAIRS_PER_MARKER:
                continue
            d = wrap(t[ai][:, None] - t[bi][None, :])
            j = np.abs(d).argmin(axis=1)
            dd = d[np.arange(d.shape[0]), j]
            keep = np.abs(dd) < WIN
            if not keep.any():
                continue
            out['dt'].append(dd[keep].astype(np.float32))
            for k in want:
                out[k].append(c[k][ai][keep])
    if not out['dt']:
        return None
    res = {k: np.concatenate(v) for k, v in out.items() if v}
    m = np.median(res['dt'])
    core = np.abs(res['dt'] - m) < CORE
    return {k: v[core] for k, v in res.items()}


def rms(x):
    return float(np.sqrt(np.mean((x - np.mean(x)) ** 2)))


def core_sigma(x, win=CORE, nbins=120):
    """Gaussian-on-flat-background sigma -- the estimator the campaign fits.

    The rms above includes the tails, so it sits well above this; the two must
    never be mixed inside one decomposition.
    """
    if x.size < 2000:
        return float('nan')
    m = float(np.median(x))
    h, e = np.histogram(x, bins=nbins, range=(m - win, m + win))
    ctr = 0.5 * (e[:-1] + e[1:])
    bg0 = float(np.median(h[np.abs(ctr - m) > win * 0.75]))
    k = int((h - bg0).argmax())
    a0, mu0 = float(h[k] - bg0), float(ctr[k])
    w = (h - bg0).clip(0)
    s0 = float(np.sqrt(np.average((ctr - mu0) ** 2, weights=w))) if w.sum() \
        else 20.0
    try:
        from scipy.optimize import curve_fit

        def f(x_, a, mu, sg, bg):
            return a * np.exp(-0.5 * ((x_ - mu) / sg) ** 2) + bg

        p_, _ = curve_fit(f, ctr, h, p0=[a0, mu0, s0, bg0],
                          sigma=np.sqrt(np.maximum(h, 1)), maxfev=20000)
        return abs(float(p_[2]))
    except Exception:
        return float('nan')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--store', default=STORE)
    ap.add_argument('--captures', type=int, default=6)
    ap.add_argument('--max-rate-khz', type=float, default=20.0,
                    help='drop channels above this absolute rate; 0 disables')
    a = ap.parse_args()
    caps = sorted(glob.glob(os.path.join(a.store, '*')))[:a.captures]
    caps = [c for c in caps if os.path.isdir(c)]
    global MR
    MR = a.max_rate_khz * 1e3
    print(f'{len(caps)} captures from {a.store}'
          f'   (hot-channel veto at {a.max_rate_khz:g} kHz)\n')

    print('== 1. what the TDC fine time is worth ==')
    for st in STATIONS:
        line = f'   {st:8s}'
        for fine, lbl in ((False, 'BCID only'), (True, 'BCID+TDC')):
            r = pair(caps, st, 'TRIG', use_fine=fine, max_rate_hz=MR)
            line += f'   {lbl} {rms(r["dt"]):5.1f} ns' if r else f'   {lbl}  --'
        print(line)
    print('   -> the clock is not the limit: quantisation alone is '
          f'{CLOCK / np.sqrt(12):.1f} ns rms\n')

    print('== 2. per-channel t0 and time walk ==')
    for st in STATIONS:
        r = pair(caps, st, 'TRIG', want=('adc', 'vmm', 'ch'),
                 max_rate_hz=MR)
        if r is None:
            print(f'   {st:8s} skipped -- no adc column in this store')
            continue
        if r['dt'].size < 5000:
            continue
        d = r['dt']
        base = rms(d)
        key = r['vmm'].astype(np.int64) * 64 + r['ch'].astype(np.int64)
        uk, inv = np.unique(key, return_inverse=True)
        cnt = np.bincount(inv, minlength=len(uk))
        med = np.array([np.median(d[inv == i]) if cnt[i] >= 200 else np.nan
                        for i in range(len(uk))])
        ok = np.isfinite(med)
        d1 = d - np.where(np.isfinite(med[inv]), med[inv], np.median(d))
        q = np.unique(np.quantile(r['adc'], np.linspace(0, 1, 11)))
        ib = np.clip(np.digitize(r['adc'], q) - 1, 0, len(q) - 2)
        wm = np.array([np.median(d1[ib == i]) if (ib == i).sum() >= 200
                       else np.nan for i in range(len(q) - 1)])
        d2 = d1 - np.where(np.isfinite(wm[ib]), wm[ib], 0.0)
        print(f'   {st:8s} {base:5.1f} -> t0 {rms(d1):5.1f} -> +walk '
              f'{rms(d2):5.1f} ns   (channel t0 spread {rms(med[ok]):.1f} ns '
              f'rms, walk {np.nanmax(wm) - np.nanmin(wm):.0f} ns over the '
              f'ADC deciles)')
    print()

    print('== 3. how much of it is the trigger channel ==')
    print(f'   {"pair":26s} {"rms":>8s} {"core fit":>9s}')
    S = {}
    for x, y in (('P2_MID', 'TRIG'), ('P2_OUT', 'TRIG'), ('P2_IN', 'TRIG'),
                 ('P2_MID', 'P2_OUT')):
        r = pair(caps, x, y, max_rate_hz=MR)
        if r is None or r['dt'].size <= 1000:
            S[(x, y)] = (float('nan'), float('nan'))
        else:
            S[(x, y)] = (rms(r['dt']), core_sigma(r['dt']))
        print(f'   sigma({x:6s} - {y:6s}) = {S[(x, y)][0]:8.1f} '
              f'{S[(x, y)][1]:9.1f} ns')

    for i, name in ((0, 'rms'), (1, 'core fit')):
        s = {k: v[i] for k, v in S.items()}
        t2 = 0.5 * (s[('P2_MID', 'TRIG')] ** 2 + s[('P2_OUT', 'TRIG')] ** 2
                    - s[('P2_MID', 'P2_OUT')] ** 2)
        if not np.isfinite(t2) or t2 <= 0:
            continue
        print(f'\n   --- solved with the {name} estimator ---')
        print(f'   trigger channel      {np.sqrt(t2):5.1f} ns')
        for st in ('P2_MID', 'P2_OUT'):
            v = s[(st, 'TRIG')] ** 2 - t2
            print(f'   {st} intrinsic     {np.sqrt(v):5.1f} ns' if v > 0
                  else f'   {st} intrinsic    unphysical')
    print('\n   -> the trigger channel is a scintillator read through a VMM '
          'discriminator.\n      Referencing to the DREAM trigger timestamp '
          'instead (stream matching,\n      10.4 ns residual rms) replaces '
          'it -- see README.md.')


if __name__ == '__main__':
    main()
