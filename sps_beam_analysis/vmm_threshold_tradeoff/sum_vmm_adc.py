import numpy as np, glob, os, sys
B = '/eos/experiment/ntof/data/x17/p2_sps_july/analysis/vmm/'
out = {}
for d in sorted(glob.glob(B + 'run_*/*/')):
    run, sub = d.rstrip('/').split('/')[-2:]
    if not run.startswith('run_') or int(run[4:]) < 31: continue
    caps = sorted(glob.glob(d + '*/counts.npz'))
    if not caps: continue
    adc = np.zeros((32, 1024), np.uint64); ot0 = np.zeros((32, 1024), np.uint64); avc = np.zeros((32, 64, 128), np.uint64)
    for c in caps:
        try:
            z = np.load(c)
            adc += z['adc']; ot0 += z['adc_ot0']; avc += z['adc_vs_ch']
        except Exception as e:
            print('skip', c, e, file=sys.stderr)
    k = f'{run}__{sub}'
    out[k + '__adc'] = adc; out[k + '__ot0'] = ot0; out[k + '__avc'] = avc
    print(k, len(caps), int(adc[14:20].sum()))
np.savez_compressed('/tmp/akallits/vmm_adc_all.npz', **out)   # runs on lxplus; copy the result next to these scripts
