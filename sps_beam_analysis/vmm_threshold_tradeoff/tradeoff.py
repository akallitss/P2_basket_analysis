"""Threshold / masking trade-off from the 29-30 Jul bench threshold scan + DREAM-spectrum model."""
import json, sys, os
import numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/'
sys.path.insert(0, REPO + 'dream_vmm_comparison')
import threshold_model as tm
SCR = os.path.dirname(os.path.abspath(__file__)) + '/'
ADC_PER_DAC, EDGE_PER_T = 0.9, 0.196
TRIM_DAC = 34 / 0.801            # full trim range in DAC counts (~42)
# chip baselines in DAC counts: mean of sdt-(edge-38)/0.9 over runs 48/52/66 (gain 4.5)
BASE = {2: 208, 3: 182, 5: 212, 8: 197, 10: 210, 11: 204, 12: 210, 13: 188,
        14: 201, 15: 180, 16: 199, 17: 187, 18: 182, 19: 204}
HOP = {4, 7, 9}                  # one runaway channel that hops between steps; not threshold noise
rows = json.load(open(SCR + 'thrscan_summary.json'))


def onsets(gain, R=10.0):
    """per chip: array(64) of first scanned sdt from which the channel stays below R Hz.
    value 210 means 'quiet over the whole scan' (true onset <= 210, unmeasured)."""
    rr = sorted([r for r in rows if r['gain'] == gain], key=lambda r: r['sdt'])
    sdts = np.array([r['sdt'] for r in rr])
    out = {}
    for v in range(2, 20):
        rate = np.array([np.array(r['per'][str(v)]['ch']) / r['span_s'] for r in rr])  # (nsdt, 64)
        noisy = rate > R
        on = np.full(64, sdts[0])
        for c in range(64):
            idx = np.where(noisy[:, c])[0]
            if len(idx):
                on[c] = sdts[idx.max() + 1] if idx.max() + 1 < len(sdts) else 999
        out[v] = on
    return out, sdts


g, H, b = tm.load(); S = tm.Spectra(H, b)
w = g.n_track_v.values.astype(float); ed = g.eff_d.values
vm = g.vmm.values.astype(int); ch = g.ch.values.astype(int)


def pad_eff(d_dac):
    """efficiency of each pad if its channel sits d_dac DAC counts above its baseline"""
    T = np.maximum(d_dac, 0) * ADC_PER_DAC / EDGE_PER_T
    # frac_above takes a scalar; evaluate per unique T
    out = np.empty(len(g))
    for t in np.unique(T):
        m = T == t
        out[m] = (ed * S.frac_above(t))[m]
    return out


def eff_d(d, gain):
    """pad efficiencies when the channel sits d DAC above its baseline (signal scales with gain)"""
    T = max(d, 0) * ADC_PER_DAC / EDGE_PER_T * 4.5 / gain
    return ed * S.frac_above(min(T, (S.cum.shape[1] - 1) * S.binw))


def station(on, d, trim):
    """all P2_OUT chips at d DAC above their own baseline.
    noisy channels: raised to their onset with trim (mask if beyond trim range), masked without."""
    dch = np.empty(len(g)); masked = np.zeros(len(g), bool)
    nmask = nraise = 0; unknown = 0
    for v in range(14, 20):
        S_t = BASE[v] + d
        if S_t < 210:
            unknown += 64
        need = on[v] > S_t
        if trim:
            mk = need & (on[v] - S_t > TRIM_DAC); rs = need & ~mk
        else:
            mk = need; rs = np.zeros(64, bool)
        nmask += mk.sum(); nraise += rs.sum()
        sel = vm == v
        thr = np.where(rs, on[v], S_t)
        dch[sel] = thr[ch[sel]] - BASE[v]
        masked[sel] = mk[ch[sel]]
    e = pad_eff(dch); e[masked] = 0
    return np.average(e, weights=w), nmask, nraise, unknown


if __name__ == '__main__':
    for gain in [4.5, 3.0]:
        on, sd = onsets(gain)
        print(f'\n=== gain {gain}: per-chip noisy-channel onsets (10 Hz), base, quiet-at-210 count')
        for v in range(2, 20):
            o = on[v]; nz = o > 210
            print(f'VMM {v:2d} base {BASE.get(v,"?")!s:>4}  quiet@210: {(~nz).sum():2d}  '
                  f'noisy chans: ' + ' '.join(f'{c}:{o[c]}' for c in np.where(nz)[0]))
    on, _ = onsets(4.5)
    print('\nP2_OUT station, gain 4.5 (eff, masked, raised, chans-below-scan)')
    for d in range(5, 50, 2):
        a = station(on, d, False); bb = station(on, d, True)
        print(f'd={d:2d} DAC ({d*0.9:4.1f} ADC): no-trim eff {a[0]:.3f} mask {a[1]:2d} | '
              f'trim eff {bb[0]:.3f} mask {bb[1]} raise {bb[2]:2d} | unmeasured chans {a[3]}')
    # where the July config sat
    cfg = {14: 230, 15: 230, 16: 230, 17: 230, 18: 230, 19: 265}
    dch = np.array([cfg[v] - BASE[v] for v in vm], float)
    print('July opt config model eff:', np.average(pad_eff(dch), weights=w), ' measured eff_v:', np.average(g.eff_v, weights=w))
