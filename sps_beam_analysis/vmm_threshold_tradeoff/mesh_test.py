"""Station-level: cut DREAM's own spectrum AT EACH MESH VOLTAGE at one fixed
threshold (fitted at 450 V only) and compare with the VMM mesh scan."""
import json, glob, numpy as np, pandas as pd
A = '/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/SPS_beam_test/TB_July2026_H4/analysis/det3/'
spec = {}
for run in ['beam_nominal_meshscan_1', 'meshscan_fine_1', 'low_mesh_scan_1']:
    for f in glob.glob(A + run + '/*/20_beam_spectra/scan_row.json'):
        d = json.load(open(f)); s = d['spectrum']
        spec[(run, s['hv'], f.split('/')[-3])] = (np.array(s['edges']), np.array(s['counts'], float))
te = pd.read_csv(A + 'beam_nominal_meshscan_1/scan/22_tag_probe_efficiency/tag_probe_efficiency_spark_vetoed.csv')
effd = {int(r.hv): r.eff for r in te.itertuples() if r.sub_run != 'nominal_01'}

def frac_above(e, c, T):
    cum = np.concatenate([[0], np.cumsum(c)]) / c.sum()
    return 1 - np.interp(T, e, cum)

def quant(e, c, q):
    cum = np.concatenate([[0], np.cumsum(c)]) / c.sum()
    return np.interp(q, cum, e)

c = pd.read_csv('/local/home/ak271430/Documents/PostDocSaclay/P2_basket_analysis/sps_beam_analysis/vmm_dream_matching/efficiency_conditions.csv')
run = 'beam_nominal_meshscan_1'
e450, c450 = spec[(run, 450, 'nominal_00')]
print('DREAM own quantiles vs mesh (beam_nominal_meshscan_1): gain proxy')
for hv in [450, 440, 430, 420, 410, 400, 390]:
    k = [x for x in spec if x[0] == run and x[1] == hv and x[2] != 'nominal_01'][0]
    e, cc = spec[k]
    print(hv, ' q50=%.0f q75=%.0f q90=%.0f  DREAM eff %.3f' % (quant(e, cc, .5), quant(e, cc, .75), quant(e, cc, .9), effd[hv]))
for vrun in ['run_32', 'run_55']:
    v = c[(c.run.isin([vrun, 'run_56'] if vrun == 'run_55' else [vrun, 'run_33'])) & c['sub'].str.startswith('mesh')]
    v = v.set_index('mesh').eff_OUT
    target = v[450] / effd[450]
    Ts = np.arange(30, 3000, 1.0)
    T = Ts[np.argmin(abs(np.array([frac_above(e450, c450, t) for t in Ts]) - target))]
    print(f'\n{vrun}: T fitted at 450 V only = {T:.0f} DREAM ADC (this run\'s own scale)')
    print(' mesh  VMM meas   predicted (DREAM spectrum at that HV, same T)')
    for hv in [450, 440, 430, 420, 410, 400, 390]:
        k = [x for x in spec if x[0] == run and x[1] == hv and x[2] != 'nominal_01'][0]
        e, cc = spec[k]
        print(f' {hv}   {v.get(hv, np.nan):.3f}    {effd[hv]*frac_above(e, cc, T):.3f}')
