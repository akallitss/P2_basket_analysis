"""run_46 per-channel hit rates in spill and between spills, from oos_bin.py's time-binned counts."""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
GUARD_S = 0.5            # keep this far from any spill bin on both sides
TRIG = 0 * 64 + 44       # VMM 0 ch 44 = DREAM trigger fan-in, fires only in spill
SDT46 = {int(k): v for k, v in json.load(open(HERE + '../runcfg/sdt.json'))['run_46'].items()}


def load(fn=HERE + 'run46_tbin.npz'):
    z = np.load(fn)
    tb, idx, n, bs = z['tbin'], z['idx'], z['n'], float(z['bin_s'])
    lo, hi = np.percentile(np.repeat(tb, 1), [0.01, 99.99])       # a handful of junk timestamps
    m = (tb >= lo) & (tb <= hi)
    tb, idx, n = tb[m] - int(lo), idx[m], n[m]
    L = int(tb.max()) + 1
    M = np.zeros((L, 20 * 64), np.int64)
    np.add.at(M, (tb, idx), n)
    trig = M[:, TRIG]
    alive = M.sum(1) > 0                                          # DAQ was writing in that bin
    spill = trig > 0
    g = int(round(GUARD_S / bs))
    near = np.convolve(spill, np.ones(2 * g + 1), 'same') > 0
    out = alive & ~near
    inn = alive & spill & (trig > 0.5 * np.median(trig[spill]))   # full-intensity spill bins only
    return dict(M=M, bs=bs, trig=trig, alive=alive, out=out, inn=inn,
                r_out=M[out].sum(0).reshape(20, 64) / (out.sum() * bs),
                r_in=M[inn].sum(0).reshape(20, 64) / (inn.sum() * bs),
                n_out=M[out].sum(0).reshape(20, 64), t_out=out.sum() * bs, t_in=inn.sum() * bs,
                n_tot=M.sum(0).reshape(20, 64))


if __name__ == '__main__':
    d = load()
    print(f"bins {len(d['trig'])}  alive {d['alive'].sum()*d['bs']:.0f} s  in-spill {d['t_in']:.0f} s  "
          f"out-of-spill {d['t_out']:.0f} s  spill bins {int((d['trig']>0).sum())}")
    for v in range(2, 20):
        ro, ri = d['r_out'][v], d['r_in'][v]
        top = np.argsort(ro)[::-1][:3]
        print(f"VMM {v:2d} sdt {SDT46[v]}  out: chip {ro.sum():9.0f} Hz  median {np.median(ro):7.2f}  "
              f">10Hz {int((ro>10).sum()):2d}  top " + ', '.join(f'ch{c}:{ro[c]:.0f}' for c in top) +
              f"   in-out chip {ri.sum()-ro.sum():8.0f} Hz")
