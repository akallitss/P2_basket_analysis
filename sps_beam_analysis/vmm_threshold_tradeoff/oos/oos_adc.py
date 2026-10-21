"""Run on lxplus, second pass: ADC spectra per (vmm, ch) in and out of spill, and how many out-of-spill
hits have a partner hit on another channel of the same station within +-150 ns (particles cluster, noise doesn't).
    python3 oos_adc.py run_46 cfg_gain4.5_peaktime200 run46_masks.npz out.npz
"""
import sys, time
import numpy as np
import extract_vmm_triggers as evt

run, sub, mfn, out = sys.argv[1:5]
m = np.load(mfn); bs = float(m['bin_s'])
OUT, INN = set(m['out'].tolist()), set(m['inn'].tolist())
lo = int(m['lo']); L = int(max(m['out'].max(), m['inn'].max()) - lo + 10)
cls = np.zeros(L + 10**6, np.int8)              # 1 out, 2 in, 0 other; generous length, clipped below
cls[np.array(sorted(OUT)) - lo] = 1; cls[np.array(sorted(INN)) - lo] = 2
H = np.zeros((2, 20, 64, 64), np.int64)          # [out/in, vmm, ch, adc//16]
coinc = np.zeros((2, 20, 64), np.int64)          # out-of-spill hits with / without a partner
WIN = 150.0
ST = np.repeat(np.arange(4), [2, 6, 6, 6])        # station id per vmm: trig, IN, MID, OUT
t0 = time.time(); caps = []
for name, cols in evt.iter_columns(evt.sub_dir_of(run, sub), None, source='pcapng'):
    v = cols['vmm'].astype(np.int64); c = cols['ch'].astype(np.int64)
    ok = (v >= 0) & (v < 20) & (c >= 0) & (c < 64)
    v, c = v[ok], c[ok]; a = cols['adc'][ok].astype(np.int64)
    t = evt.hit_times(cols, ok)
    k = np.floor(t / (bs * 1e9)).astype(np.int64) - lo
    good = (k >= 0) & (k < len(cls))
    cl = np.zeros(len(k), np.int8); cl[good] = cls[k[good]]
    for j, sel in ((0, cl == 1), (1, cl == 2)):
        np.add.at(H[j], (v[sel], c[sel], np.clip(a[sel] // 16, 0, 63)), 1)
    o = (cl == 1) & (v >= 2)
    if o.any():
        st = ST[v[o]]; to = t[o]; vo, co = v[o], c[o]
        for s in (1, 2, 3):
            q = st == s
            if q.sum() < 2: continue
            tt = to[q]; ordr = np.argsort(tt); tt = tt[ordr]; idq = (vo[q] * 64 + co[q])[ordr]
            dn = np.r_[np.inf, np.diff(tt)]; up = np.r_[np.diff(tt), np.inf]
            # partner must be a different channel: approximate with nearest neighbours in time
            pn = np.r_[-1, idq[:-1]]; pu = np.r_[idq[1:], -1]
            has = ((dn < WIN) & (pn != idq)) | ((up < WIN) & (pu != idq))
            np.add.at(coinc[0], (idq[has] // 64, idq[has] % 64), 1)
            np.add.at(coinc[1], (idq[~has] // 64, idq[~has] % 64), 1)
    if name not in caps:
        caps.append(name); print(len(caps), name, f'{time.time()-t0:.0f} s', flush=True)
np.savez_compressed(out, H=H, coinc=coinc, captures=np.array(caps))
print('DONE', len(caps), f'{time.time()-t0:.0f} s')
