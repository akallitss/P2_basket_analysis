"""Per-chip, per-channel hit counts of every step of the 29-30 Jul bench threshold scan -> thrscan_summary.json.

Input: thrscan/qa/<step>/*.root, the qa_plots hit trees copied from
vmm_daplxa:/local/p2/p2testbench/TestBenchCERN/config_scan/ (the pcapng files there are unreadable).
Span = 0.1-99.9 % of abs_time_ns, times 0.9 for the 22.5 ns SRS tick that the tree assumes is 25 ns.
"""
import glob, json, os, re
import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
rows = []
for f in sorted(glob.glob(HERE + 'thrscan/qa/*/*.root')):
    m = re.search(r'gain([\d.]+)_.*scan_(\d+)_', f); g = float(m.group(1)); sdt = int(m.group(2))
    a = uproot.open(f)['hits'].arrays(['vmm', 'ch', 'adc', 'abs_time_ns'], library='np')
    span = (np.percentile(a['abs_time_ns'], 99.9) - np.percentile(a['abs_time_ns'], 0.1)) * 0.9e-9
    per = {}
    for v in range(2, 20):
        m_ = a['vmm'] == v
        per[v] = dict(n=int(m_.sum()), ch=np.bincount(a['ch'][m_], minlength=64).tolist(),
                      adc_med=float(np.median(a['adc'][m_])) if m_.any() else None,
                      adc_p01=float(np.percentile(a['adc'][m_], 1)) if m_.sum() > 20 else None)
    rows.append(dict(gain=g, sdt=sdt, span_s=span, n=int(len(a['vmm'])), per=per))
json.dump(rows, open(HERE + 'thrscan_summary.json', 'w'))
print(len(rows), 'scan steps')
