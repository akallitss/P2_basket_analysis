"""Does the run_46 threshold model predict the VMM runs where only the
amplifier gain (or sdt, or mesh HV) changed?  Per pad, DREAM spectra fixed."""
import json, os, sys, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'dream_vmm_comparison'))
import threshold_model as tm
S0 = os.path.dirname(os.path.abspath(__file__)) + '/eff/'

g, H, binw = tm.load()                       # 53 used pads, DREAM eff_nominal_1
S = tm.Spectra(H, binw)
full = pd.read_csv(tm.DATA + '/compare_dream_vmm_P2_OUT.csv')
idx = full.index[full['use']].to_numpy()      # rows of the 384-pad per_pad arrays
print('spectrum range', H.shape, binw, 'max ADC', H.shape[1]*binw)

def vmm(run_sub):
    f = glob.glob(S0 + f'vmm_eff_{run_sub}*.json')[0]
    d = json.load(open(f))
    s = [x for x in d['stations'] if x['station'] == 'P2_OUT'][0]['per_pad']
    n = np.array(s['n'])[idx].astype(float); k = np.array(s['k'])[idx].astype(float)
    return n, k

effd = g.eff_d.to_numpy()
Ts = np.arange(4.0, 3000.0, 1.0)
F = np.array([S.frac_above(T) for T in Ts])          # (nT, npad)

def fit(n, k, mask):
    e = k / np.maximum(n, 1)
    chi = ((n * (effd * F - e) ** 2)[:, mask]).sum(1)
    i = np.argmin(chi); T = Ts[i]
    pred = effd * F[i]
    ok = chi <= chi.min() * 1.21
    return T, Ts[ok].min(), Ts[ok].max(), pred, e

def summ(n, k, pred, mask):
    e = k / np.maximum(n, 1)
    w = n[mask]
    return dict(meas=float((k[mask]).sum() / w.sum()),
                pred=float((pred[mask] * w).sum() / w.sum()),
                r=float(np.corrcoef(pred[mask], e[mask])[0, 1]),
                rms=float(np.sqrt(np.average((pred[mask] - e[mask]) ** 2, weights=w))))

not19 = (g.vmm.to_numpy() != 19)
print('pads used', len(g), 'excluding VMM19:', not19.sum())

out = {}
print('\n=== 1. amplifier gain step at fixed sdt 230 (VMM19 pads dropped: sdt 265 at 4.5) ===')
for lo, hi in [('run_43_cfg_gain3.0_peaktime200_opt', 'run_48_cfg_gain4.5_peaktime200_opt'),
               ('run_50_cfg_gain3.0_peaktime200_opt', 'run_52_cfg_gain4.5_peaktime200_opt')]:
    nh, kh = vmm(hi); nl, kl = vmm(lo)
    Th, a, b, ph, _ = fit(nh, kh, not19)
    Tl, c, d_, pl_free, _ = fit(nl, kl, not19)
    i = np.argmin(abs(Ts - 1.5 * Th)); pl = effd * F[i]
    sh = summ(nh, kh, ph, not19); sl = summ(nl, kl, pl, not19); slf = summ(nl, kl, pl_free, not19)
    print(f'{hi[:6]} gain4.5: T={Th:.0f} [{a:.0f},{b:.0f}]  meas {sh["meas"]:.3f} fit {sh["pred"]:.3f} r={sh["r"]:.2f}')
    print(f'{lo[:6]} gain3.0: PREDICTED (T=1.5x{Th:.0f}={1.5*Th:.0f}, no fit): {sl["pred"]:.3f}  measured {sl["meas"]:.3f}  per-pad r={sl["r"]:.2f} rms={sl["rms"]:.3f}')
    print(f'          free fit T={Tl:.0f} [{c:.0f},{d_:.0f}] -> T30/T45 = {Tl/Th:.2f}  (expect 1.50; band {c/b:.2f}-{d_/a:.2f}) r={slf["r"]:.2f}')
    out[lo+'/'+hi] = dict(T45=Th, T30_free=Tl, ratio=Tl/Th, pred30=sl['pred'], meas30=sl['meas'], r30=sl['r'], meas45=sh['meas'])

print('\n=== 2. sdt step 224 -> 230 at gain 4.5 (VMM19 dropped) ===')
for rs in ['run_46_cfg_gain4.5_peaktime200', 'run_48_cfg_gain4.5_peaktime200_opt', 'run_52_cfg_gain4.5_peaktime200_opt', 'run_45_cfg_gain4.5_peaktime100']:
    n, k = vmm(rs); T, a, b, p, _ = fit(n, k, not19); s = summ(n, k, p, not19)
    print(f'{rs:40s} T={T:.0f} [{a:.0f},{b:.0f}] meas {s["meas"]:.3f} r={s["r"]:.2f}')

print('\n=== 2b. VMM19 alone (sdt 224 in run_46, 265 in run_48/52) ===')
m19 = ~not19
for rs in ['run_46_cfg_gain4.5_peaktime200', 'run_48_cfg_gain4.5_peaktime200_opt', 'run_52_cfg_gain4.5_peaktime200_opt']:
    n, k = vmm(rs); T, a, b, p, _ = fit(n, k, m19); s = summ(n, k, p, m19)
    print(f'{rs:40s} VMM19 pads={m19.sum()} T={T:.0f} [{a:.0f},{b:.0f}] meas {s["meas"]:.3f}')
    nn, kk = n[not19], k[not19]
    print(f'{"":40s} other chips meas {kk.sum()/nn.sum():.3f}')

print('\n=== 3. mesh scans, gain 3.0 sdt 230 (run_55/56) and run_32/33 ===')
allm = np.ones(len(g), bool)
rows = []
for f in sorted(glob.glob(S0 + 'vmm_eff_run_5[56]_meshscan*.json')) + sorted(glob.glob(S0 + 'vmm_eff_run_3[23]_meshscan*.json')):
    rs = f.split('vmm_eff_')[1][:-5]
    dv = -float(rs.split('_m')[-1][:-1])
    n, k = vmm(rs)
    if n.sum() == 0: continue
    T, a, b, p, e = fit(n, k, allm); s = summ(n, k, p, allm)
    rows.append(dict(run=rs[:6], dv=dv, T=T, lo=a, hi=b, meas=s['meas'], pred=s['pred'], r=s['r']))
m = pd.DataFrame(rows).sort_values(['run', 'dv'], ascending=[True, False]); print(m.to_string(index=False))
for run, q in m.groupby('run'):
    q = q[q['T'] < 2500]
    sl = np.polyfit(q.dv, np.log(q['T']), 1)
    print(f'{run}: log T vs dV slope -> V0 = {-1/sl[0]:.1f} V  (T at dV=0: {np.exp(sl[1]):.0f})  using {len(q)} points')
json.dump(dict(gain=out, mesh=m.to_dict('records')), open(S0 + '../gain_step_test.json', 'w'), indent=1, default=float)

print('\n=== 4. per-chip threshold at one sdt (run_46, all at 224) ===')
n, k = vmm('run_46_cfg_gain4.5_peaktime200')
for v in sorted(g.vmm.unique()):
    mk = g.vmm.to_numpy() == v
    if mk.sum() < 2: continue
    T, a, b, p, _ = fit(n, k, mk)
    print(f'VMM {v}: pads {mk.sum():2d}  T={T:.0f} [{a:.0f},{b:.0f}]  eff {k[mk].sum()/n[mk].sum():.3f}')
Tp, sT, lc, hc = tm.fit_threshold_per_pad(g, S)
ch = g.vmm.to_numpy()
grand = Tp.mean(); betw = sum((ch == v).sum() * (Tp[ch == v].mean() - grand) ** 2 for v in np.unique(ch))
tot = ((Tp - grand) ** 2).sum()
print(f'per-pad T: fraction of variance between chips = {betw/tot:.2f}')
