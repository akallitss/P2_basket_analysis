"""per capture file: decoded hit-time span vs the pcapng packet timestamps span, and the pulses received.
Usage (from $W, env.sh, PYTHONPATH = cva copy): file_clock.py PID [PID ...]   (prints one line per capture file;
with the 22.5 ns tick the decoded span equals the packet span to ~0.02 %, with the old 25 ns tick it is 1.111x).
"""
import glob, struct, sys, numpy as np
from coincidence_check import read_times
P = 4096 * 22.5

def pcap_ts(path):
    """first and last Enhanced Packet Block timestamps (s), assuming the default microsecond resolution."""
    first = last = None; n = 0
    with open(path, 'rb') as f:
        data = f.read()
    off = 0
    while off + 12 <= len(data):
        bt, bl = struct.unpack_from('<II', data, off)
        if bl < 12: break
        if bt == 6:
            hi, lo = struct.unpack_from('<II', data, off + 12)
            ts = ((hi << 32) | lo) * 1e-6
            first = ts if first is None else first; last = ts; n += 1
        off += bl
    return first, last, n

for pid in sys.argv[1:]:
    d = sorted(glob.glob(f'pseudo/run_{pid}/cfg_*'))[0]
    for fn in sorted(glob.glob(f'{d}/raw_daq_data/*.pcapng')):
        state = {}
        recs = read_times(fn, state, with_adc=True)
        V = np.concatenate([r[0] for r in recs]); T = np.concatenate([r[1] for r in recs]); C = np.concatenate([r[2] for r in recs])
        t0, t1, n = pcap_ts(fn)
        tp = np.sort(T[(V == 5) & (C == 31)])
        print(f'{pid} {fn.split("/")[-1]}: packets {n}, packet-time span {t1 - t0:.3f} s, decoded hit span '
              f'{(T.max() - T.min()) * 1e-9:.3f} s, ratio {(T.max() - T.min()) * 1e-9 / (t1 - t0):.4f}; chip-5 pulses {len(tp)}'
              f' = {len(tp) * P * 1e-9:.3f} s of pulsing; monotonic share of T {np.mean(np.diff(T) >= 0):.4f}')
