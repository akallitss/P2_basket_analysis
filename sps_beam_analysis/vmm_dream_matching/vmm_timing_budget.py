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


def picker(c, name):
    if name == 'TRIG':
        return (c['vmm'] == TRG_VMM) & (c['ch'] == TRG_CH)
    return np.isin(c['vmm'], np.fromiter(STATIONS[name], int))


def pair(caps, a, b, use_fine=True, want=()):
    """Nearest-in-time partner for every `a` hit, within one marker interval."""
    need = ['vmm', 'ch', 'bcid', 'tdc', 'srs_timestamp'] + list(want)
    out = {k: [] for k in ('dt',) + tuple(want)}
    for cap in caps:
        c = load(cap, need)
        t = times(c, use_fine)
        ia, ib = picker(c, a), picker(c, b)
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--store', default=STORE)
    ap.add_argument('--captures', type=int, default=6)
    a = ap.parse_args()
    caps = sorted(glob.glob(os.path.join(a.store, '*')))[:a.captures]
    print(f'{len(caps)} captures from {a.store}\n')

    print('== 1. what the TDC fine time is worth ==')
    for st in STATIONS:
        line = f'   {st:8s}'
        for fine, lbl in ((False, 'BCID only'), (True, 'BCID+TDC')):
            r = pair(caps, st, 'TRIG', use_fine=fine)
            line += f'   {lbl} {rms(r["dt"]):5.1f} ns' if r else f'   {lbl}  --'
        print(line)
    print('   -> the clock is not the limit: quantisation alone is '
          f'{CLOCK / np.sqrt(12):.1f} ns rms\n')

    print('== 2. per-channel t0 and time walk ==')
    for st in STATIONS:
        r = pair(caps, st, 'TRIG', want=('adc', 'vmm', 'ch'))
        if r is None or r['dt'].size < 5000:
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
    s = {}
    for x, y in (('P2_MID', 'TRIG'), ('P2_OUT', 'TRIG'), ('P2_IN', 'TRIG'),
                 ('P2_MID', 'P2_OUT')):
        r = pair(caps, x, y)
        s[(x, y)] = rms(r['dt']) if r is not None and r['dt'].size > 1000 \
            else float('nan')
        print(f'   sigma({x:6s} - {y:6s}) = {s[(x, y)]:5.1f} ns')
    t2 = 0.5 * (s[('P2_MID', 'TRIG')] ** 2 + s[('P2_OUT', 'TRIG')] ** 2
                - s[('P2_MID', 'P2_OUT')] ** 2)
    if t2 > 0:
        print(f'\n   trigger channel      {np.sqrt(t2):5.1f} ns')
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
