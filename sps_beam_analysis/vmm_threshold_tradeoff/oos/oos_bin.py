"""Run on lxplus: every hit of one VMM sub_run binned in time (0.1 s) per (vmm, ch).

Output is sparse (bin, vmm*64+ch, count) so spills can be found afterwards from the
trigger channel and in-/out-of-spill counts taken per channel.
    python3 oos_bin.py run_46 cfg_gain4.5_peaktime200 out.npz [--captures N]
"""
import sys, time, argparse
import numpy as np
import extract_vmm_triggers as evt

ap = argparse.ArgumentParser()
ap.add_argument('run'); ap.add_argument('sub'); ap.add_argument('out')
ap.add_argument('--captures', type=int, default=None)
ap.add_argument('--bin-s', type=float, default=0.1)
a = ap.parse_args()
NI = 20 * 64
parts, caps, nhit, nbad = [], [], 0, 0
t0 = time.time()
for name, cols in evt.iter_columns(evt.sub_dir_of(a.run, a.sub), a.captures, source='pcapng'):
    v = cols['vmm'].astype(np.int64); c = cols['ch'].astype(np.int64)
    ok = (v >= 0) & (v < 20) & (c >= 0) & (c < 64)
    t = evt.hit_times(cols, ok)
    k = np.floor(t / (a.bin_s * 1e9)).astype(np.int64)
    key, n = np.unique(k * NI + v[ok] * 64 + c[ok], return_counts=True)
    parts.append(np.stack([key, n]))
    nhit += int(ok.sum()); nbad += int((~ok).sum())
    if name not in caps:
        caps.append(name); print(f'{len(caps):3d} {name} {nhit/1e6:.0f} M hits {time.time()-t0:.0f} s', flush=True)
        # merge to keep memory flat
        allk = np.concatenate(parts, axis=1); u, inv = np.unique(allk[0], return_inverse=True)
        parts = [np.stack([u, np.bincount(inv, weights=allk[1]).astype(np.int64)])]
allk = np.concatenate(parts, axis=1); u, inv = np.unique(allk[0], return_inverse=True)
cnt = np.bincount(inv, weights=allk[1]).astype(np.int64)
np.savez_compressed(a.out, tbin=u // NI, idx=u % NI, n=cnt, bin_s=a.bin_s, captures=np.array(caps),
                    nhit=nhit, nbad=nbad)
print('DONE', len(caps), 'captures', nhit, 'hits', nbad, 'bad', f'{time.time()-t0:.0f} s')
