import numpy as np
from tradeoff import eff_d, onsets, BASE, g, w, vm, ch, TRIM_DAC

def trimmed(gain, d, R=10.0, pessimistic=False):
    """clean channels at d DAC above baseline; noisy channels raised to their own quiet sdt."""
    on, _ = onsets(gain, R)
    e_all = np.empty(len(g)); nraise = nmask = 0; uni = []
    for v in range(14, 20):
        St = BASE[v] + d
        req = np.where(on[v] > 210, on[v], St if not pessimistic else max(St, 210))
        req = np.maximum(req, St)
        mk = (req - St > TRIM_DAC) | (on[v] >= 999)
        nraise += ((req > St) & ~mk).sum(); nmask += mk.sum()
        de = {r: eff_d(r - BASE[v], gain) for r in np.unique(req)}
        sel = np.where(vm == v)[0]
        for i in sel:
            e_all[i] = 0 if mk[ch[i]] else de[req[ch[i]]][i]
        uni.append(np.mean([0 if mk[c] else np.average(de[req[c]], weights=w) for c in range(64)]))
    return np.average(e_all, weights=w), np.mean(uni), nraise, nmask

if __name__ == '__main__':
    for gain in [4.5, 3.0]:
        print('gain', gain)
        for d in [5, 8, 10, 12, 15, 18, 20, 23, 25, 30]:
            o = trimmed(gain, d); p = trimmed(gain, d, pessimistic=True)
            print(f'  d={d:2d} ({d*0.9:4.1f} ADC, {d*0.9/(1.12*gain):.1f} fC): optimistic beam {o[0]:.3f} uni {o[1]:.3f} raise {o[2]:3d} mask {o[3]} | pessimistic beam {p[0]:.3f} uni {p[1]:.3f} raise {p[2]}')
