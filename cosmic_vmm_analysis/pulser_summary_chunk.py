"""pulser chunk: pulses received per pulsed chip against the real duration (pcapng packet timestamps), per file;
and from the tag dump whether the chamber's pulsed chips had hits within +-100 us of every ON muon (silent / detected).
Usage (from the weekend dir $W, env.sh sourced, PYTHONPATH = the cva copy): pulser_summary_chunk.py PID OUT_DIR
  -> OUT_DIR/ps_run<PID>.json; plots: profile_pulser_plots.py PROFILE_DIR OUT_DIR(ps_*.json) RESULTS_DIR.
"""
import glob, json, struct, sys, numpy as np
from coincidence_check import read_times
P = 4096 * 22.5
CHIPS = (2, 5, 8, 11, 14, 17)
COL = {3: (0, 3), 2: (0, 3), 1: (0, 3)}      # columns of the pulsed chips in chip_n100us (chamber chips ascending)

def pcap_span(path):
    data = open(path, 'rb').read(); off = 0; first = last = None
    while off + 12 <= len(data):
        bt, bl = struct.unpack_from('<II', data, off)
        if bl < 12: break
        if bt == 6:
            hi, lo = struct.unpack_from('<II', data, off + 12); ts = ((hi << 32) | lo) * 1e-9
            first = ts if first is None else first; last = ts
        off += bl
    return last - first

pid, out = sys.argv[1], sys.argv[2]
d = sorted(glob.glob(f'pseudo/run_{pid}/cfg_*'))[0]
rows = []
for fn in sorted(glob.glob(f'{d}/raw_daq_data/*.pcapng')):
    recs = read_times(fn, {}, with_adc=True)
    V = np.concatenate([r[0] for r in recs]); T = np.concatenate([r[1] for r in recs]); C = np.concatenate([r[2] for r in recs])
    span = pcap_span(fn)
    for v in CHIPS:
        t = np.sort(T[(V == v) & (C == 31)])
        keep, last = 0, -np.inf
        for x in t:
            if x - last > 0.5 * P: keep += 1; last = x
        rows.append(dict(pid=int(pid), file=fn.split('/')[-1], chip=v, real_s=span, decoded_s=float(T.max() - T.min()) * 1e-9,
                         pulses=keep, expected=span / (P * 1e-9)))
Z = np.load(f'tags/tags_run{pid}.npz')
mu = []
for p in (3, 2, 1):
    dts = Z[f'p{p}_dt_scint']; idx = np.nonzero((dts >= 140) & (dts <= 310))[0]
    c100 = Z[f'p{p}_chip_n100us'][idx]
    for i, k in enumerate(idx):
        mu.append(dict(pid=int(pid), det=p, on=bool(Z[f'p{p}_on'][k]), silent=bool(Z[f'p{p}_n_any'][k] == 0),
                       pulsedA=int(c100[i, COL[p][0]]), pulsedB=int(c100[i, COL[p][1]]), r=float(Z[f'p{p}_r_scale'])))
json.dump(dict(files=rows, muons=mu), open(f'{out}/ps_run{pid}.json', 'w'))
print(pid, 'done', sum(r['pulses'] for r in rows if r['chip'] == 5), round(sum(r['expected'] for r in rows if r['chip'] == 5)))
