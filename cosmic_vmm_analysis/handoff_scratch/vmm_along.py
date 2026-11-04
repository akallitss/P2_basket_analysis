import json, numpy as np
J = '/local/home/ak271430/Documents/PostDocSaclay/data/LAB_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_wk2/eff_definitions/json/scintor_run187.json'
d = json.load(open(J))
c = np.array([381., 148.])
a = np.radians(21.5)
u = np.array([np.cos(a), np.sin(a)])
v = np.array([-np.sin(a), np.cos(a)])
for k in ('probe_det3 P2_IN', 'probe_det2 P2_MID'):
    p = d[k]['pads']
    xy = np.c_[p['x'], p['y']] - c
    vm = np.array(p['vmm'])
    print(k)
    for m in np.unique(vm):
        s = vm == m
        al, ac = xy[s] @ u, xy[s] @ v
        print(f'  VMM {m:2d}: along {al.min():6.0f} .. {al.max():6.0f} (median {np.median(al):5.0f}), across {ac.min():5.0f} .. {ac.max():5.0f}')
