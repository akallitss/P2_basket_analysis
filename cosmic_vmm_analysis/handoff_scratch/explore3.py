import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'explore2.py')).read().split("ref = F.nll(base)")[0])
effmap = {nm: T.local_efficiency(runs, nm, pads[nm], 30.0) for nm in T.NAMES}
bins = np.arange(-360, 361, 15)
D = {nm: np.histogram((np.c_[data[nm][0], data[nm][1]] - centre) @ u, bins)[0] for nm in T.NAMES}


def nll1(p):
    m = P.trigger(p)
    out = {}
    for nm in T.NAMES:
        ok, xy = P.pred[nm]
        s = m & ok
        h = np.histogram((xy[s] - centre) @ u, bins, weights=P.weight[nm][s])[0]
        mu = D[nm].sum() * (h + 0.5) / (h + 0.5).sum()
        out[nm] = float((mu - D[nm] * np.log(mu)).sum())
    return out


for lab, em in (('geometry', None), ('geometry x tag eff', effmap)):
    P.set_tag_efficiency(em)
    r0 = nll1(base)
    print(lab, ': along-band 1D dNLL per probe (det3, det2, det1, sum) vs bottom-paddle shift along the band')
    for al in range(0, -401, -50):
        p = base.copy()
        p[3:5] += al * u
        r = nll1(p)
        d = [r[nm] - r0[nm] for nm in T.NAMES]
        print(f'{al:5d}', ' '.join(f'{x:7.0f}' for x in d), f'{sum(d):8.0f}')
