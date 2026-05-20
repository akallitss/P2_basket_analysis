#!/usr/bin/env python3
"""
P2 Basket Micromegas readout PCB – channel mapping from KiCad Gerber files.

Signal path:
  RotRect SMDPad (connector signal pin, F_Cu) — carries TO.P (J<n>, pin, pinname)
                                                      and TO.N (/Sector<n>/Sig<m>)
    → F_Cu trace → ViaPad (0.8 mm, carries same TO.N) → B_Cu trace → strip endpoint
    → FPC connector → K59V adapter card (Fx2Mec8.net) → MEC8 connector + pin → VMM

Strategy (net-name guided, exact):
  1. Read F_Cu RotRect connector pads with TO.P / TO.N attributes
     → 1280 signal pads, each tagged with connector, pin, sector, strip
  2. Read F_Cu via pads (0.8 mm circles) with TO.N → net → via position
  3. Read B_Cu conductor polylines; for each long polyline find the endpoint
     nearest to a known via → net → strip endpoint
  4. Join connector pad ↔ via ↔ strip endpoint by net name (exact, no geometry guessing)
  5. (Optional) parse K59V adapter netlist → strip → MEC8 connector + pin

Output CSV:
  pad_number, x, y, connector_number, connector_pin, pin_name,
  sector, strip, channel_id, pin_x, pin_y, pad_angle,
  radius, phi, delta_phi,            ← polar coords relative to fan apex
  mec8_connector, mec8_pin          ← only when --k59v netlist is provided
"""

import re, csv, argparse
import numpy as np
from scipy.spatial import KDTree


# ── F_Cu attribute parser ──────────────────────────────────────────────────────

def parse_fcu_attributes(filepath):
    """
    Parse F_Cu Gerber for connector pad and via pad attributes.

    Reads TO.P (ref, pin number, pin name) and TO.N (net name) attributes
    attached to each pad flash.

    Returns
    -------
    signal_pads : list of dicts
        One entry per RotRect signal pad whose net name contains 'Sig'.
        Keys: ref, pin, pinname, net, sector, strip, x, y
    net_to_via : dict
        {net_name: (x_via, y_via)}
        One entry per signal net, from the 0.8 mm circular via-pad flash.
    """
    scale = 1e-6
    aperture_defs   = {}
    current_net     = None
    current_ref     = None
    current_pin     = None
    current_pinname = None
    current_ap      = None

    signal_pads = []
    net_to_via  = {}

    with open(filepath) as fh:
        for raw in fh:
            line = raw.strip()

            # Aperture definition  %ADD<n><type>,<params>*%
            m = re.match(r'%ADD(\d+)(\w+),(.+?)\*%', line)
            if m:
                aperture_defs[m.group(1)] = (m.group(2), m.group(3))
                continue

            # Net name attribute  %TO.N,<net>*%
            m = re.search(r'%TO\.N,(.+?)\*%', line)
            if m:
                current_net = m.group(1)
                continue

            # Pad reference attribute  %TO.P,<ref>,<pin>,<pinname>*%
            m = re.search(r'%TO\.P,([^,]+),([^,]+),([^*]+)\*%', line)
            if m:
                current_ref     = m.group(1)
                current_pin     = m.group(2)
                current_pinname = m.group(3)
                continue

            # Clear attributes  %TD*%
            if '%TD*%' in line:
                current_net = current_ref = current_pin = current_pinname = None
                continue

            # Aperture select  D<n>*
            m = re.match(r'^D(\d+)\*$', line)
            if m:
                current_ap = m.group(1)
                continue

            # Pad flash  X<x>Y<y>D03*
            m = re.match(r'X(-?\d+)Y(-?\d+)D03\*', line)
            if m and current_ap:
                x = int(m.group(1)) * scale
                y = int(m.group(2)) * scale
                atype, asize = aperture_defs.get(current_ap, ('?', '0'))

                # RotRect connector signal pad
                if (atype == 'RotRect'
                        and current_net and 'Sig' in current_net
                        and current_ref and current_ref.startswith('J')
                        and current_pinname and current_pinname.startswith('S')):
                    mn = re.search(r'/Sector(\d+)/Sig(\d+)', current_net)
                    parts = asize.split('X')
                    angle = float(parts[2]) if len(parts) == 3 else 0.0
                    if mn:
                        signal_pads.append({
                            'ref':     current_ref,
                            'pin':     int(current_pin),
                            'pinname': current_pinname,
                            'net':     current_net,
                            'sector':  int(mn.group(1)),
                            'strip':   int(mn.group(2)),
                            'angle':   round(angle, 3),
                            'x':       x,
                            'y':       y,
                        })

                # C,0.8 mm via pad with signal net (first flash per net wins)
                elif (atype == 'C'
                        and current_net and 'Sig' in current_net
                        and abs(float(asize) - 0.8) < 0.05
                        and current_net not in net_to_via):
                    net_to_via[current_net] = (x, y)

    return signal_pads, net_to_via


# ── B_Cu polyline parser ───────────────────────────────────────────────────────

def parse_bcu_polylines(filepath, trace_max_diam=0.2):
    """
    Parse conductor polylines from a KiCad B_Cu Gerber file.

    Groups consecutive D01 draw commands (Conductor/EtchedComponent aperture)
    into polylines, broken at D02 moves or aperture changes.

    Returns list of polylines, each a list of [x, y] points (mm).
    """
    aperture_func = {}
    aperture_def  = {}
    current_func  = None
    current_ap    = None
    scale = 1e-6

    polylines = []
    cur_poly  = []
    prev_x = prev_y = None

    TRACE_FUNCS = ('EtchedComponent', 'Conductor')

    def flush():
        if len(cur_poly) >= 2:
            polylines.append(list(cur_poly))
        cur_poly.clear()

    with open(filepath) as fh:
        for raw in fh:
            line = raw.strip()
            line = re.sub(r'\*%?$', '', line)
            if line.startswith('%'):
                line = line[1:]

            m = re.match(r'TA\.AperFunction,(.+)$', line)
            if m:
                current_func = m.group(1)
                continue
            if line == 'TD':
                current_func = None
                continue

            m = re.match(r'ADD(\d+)([A-Za-z]+),(.+)$', line)
            if m:
                aid, atype, params = m.group(1), m.group(2), m.group(3)
                aperture_def[aid] = (atype, params)
                if current_func:
                    aperture_func[aid] = current_func
                continue

            m = re.match(r'^D(\d+)$', line)
            if m:
                new_ap = m.group(1)
                if new_ap != current_ap:
                    flush()
                current_ap = new_ap
                continue

            m = re.match(r'X(-?\d+)Y(-?\d+)D(\d+)$', line)
            if m:
                xi, yi, op = int(m.group(1)), int(m.group(2)), m.group(3)
                nx, ny = xi * scale, yi * scale

                if op == '02':
                    flush()
                    prev_x, prev_y = nx, ny
                    cur_poly.clear()
                elif op == '01' and current_ap:
                    atype, params = aperture_def.get(current_ap, ('?', ''))
                    func = aperture_func.get(current_ap, '')
                    is_signal = (func in TRACE_FUNCS
                                 and atype == 'C'
                                 and float(params) < trace_max_diam)
                    if is_signal and prev_x is not None:
                        if not cur_poly:
                            cur_poly.append([prev_x, prev_y])
                        cur_poly.append([nx, ny])
                    else:
                        flush()
                    prev_x, prev_y = nx, ny
                elif op == '03' and current_ap:
                    flush()
                    prev_x, prev_y = nx, ny

    flush()
    return polylines


# ── Via-guided strip endpoint finder ──────────────────────────────────────────

def find_strip_endpoints(bcu_path, net_to_via, min_span_mm=5.0, via_threshold=0.5):
    """
    For each signal net, find the far endpoint of the B_Cu trace whose
    near endpoint sits at the via position.

    Parameters
    ----------
    bcu_path      : path to P2_BASKET-B_Cu.gbr
    net_to_via    : {net_name: (x_via, y_via)} from parse_fcu_attributes
    min_span_mm   : minimum head-to-tail span to consider a polyline as a signal trace
    via_threshold : max allowed distance (mm) from polyline endpoint to via centre

    Returns
    -------
    net_to_strip : dict {net_name: (x_strip, y_strip)}
    """
    polys_b = parse_bcu_polylines(bcu_path)
    print(f"INFO: {len(polys_b)} B_Cu polylines parsed")

    via_nets = list(net_to_via.keys())
    via_pts  = np.array([net_to_via[n] for n in via_nets])
    via_tree = KDTree(via_pts)

    # net → (x_strip, y_strip, d_best) — keep closest-to-via endpoint per net
    net_to_strip_raw = {}

    for poly in polys_b:
        head = np.array(poly[0])
        tail = np.array(poly[-1])
        if np.linalg.norm(head - tail) < min_span_mm:
            continue

        dh, vi_h = via_tree.query(head)
        dt, vi_t = via_tree.query(tail)

        if dh < dt and dh < via_threshold:
            net = via_nets[vi_h]
            if net not in net_to_strip_raw or dh < net_to_strip_raw[net][2]:
                net_to_strip_raw[net] = (float(tail[0]), float(tail[1]), dh)
        elif dt <= dh and dt < via_threshold:
            net = via_nets[vi_t]
            if net not in net_to_strip_raw or dt < net_to_strip_raw[net][2]:
                net_to_strip_raw[net] = (float(head[0]), float(head[1]), dt)

    matched = len(net_to_strip_raw)
    print(f"INFO: {matched}/{len(net_to_via)} signal nets matched to strip endpoint")

    if matched < len(net_to_via):
        missing = [n for n in net_to_via if n not in net_to_strip_raw]
        print(f"  Unmatched: {missing[:10]}{'…' if len(missing) > 10 else ''}")

    dists = [v[2] for v in net_to_strip_raw.values()]
    if dists:
        print(f"INFO: Via-to-endpoint dist — "
              f"min={min(dists):.4f}  median={np.median(dists):.4f}  "
              f"max={max(dists):.4f} mm")

    return {net: (x, y) for net, (x, y, _) in net_to_strip_raw.items()}


# ── K59V adapter card netlist parser ──────────────────────────────────────────

def parse_k59v_netlist(filepath):
    """
    Parse the Fx2Mec8 KiCad netlist (KiCad format E) to extract the
    strip → MEC8 connector + pin mapping.

    The netlist lists 128 signal nets (/Sig1 … /Sig128).  Each net has two
    nodes:
      - J2  (FPC socket):  pinfunction "S<strip>_<fpc_pin>" encodes the strip
      - JM_1 or JM_2 (MEC8 connectors): pin gives the MEC8 pin number

    Returns
    -------
    strip_to_mec8 : dict  {strip_number: (mec8_connector_idx, mec8_pin)}
        strip_number      : int 1–128
        mec8_connector_idx: 0 for JM_1, 1 for JM_2
        mec8_pin          : int (MEC8 physical pin number)
    """
    text = open(filepath).read()

    # Extract all signal net blocks
    net_blocks = re.findall(
        r'\(net\s+\(code "[^"]+"\)\s+\(name "(/Sig\d+)"\)(.*?)\n\t\t\)',
        text, re.DOTALL)

    strip_to_mec8 = {}
    for name, body in net_blocks:
        nodes = re.findall(
            r'\(node\s+\(ref "([^"]+)"\)\s+\(pin "([^"]+)"\)'
            r'\s+\(pinfunction "([^"]+)"\)',
            body)
        strip = mec8_idx = mec8_pin = None
        for ref, pin, pfunc in nodes:
            if ref == 'J2':
                m = re.match(r'S(\d+)_\d+', pfunc)
                if m:
                    strip = int(m.group(1))
            elif ref.startswith('JM_'):
                mec8_idx = int(ref.split('_')[1]) - 1   # JM_1 → 0, JM_2 → 1
                mec8_pin = int(pin)
        if strip is not None and mec8_idx is not None:
            strip_to_mec8[strip] = (mec8_idx, mec8_pin)

    print(f"INFO: K59V netlist — {len(strip_to_mec8)} strips mapped to MEC8 pins")
    jm_counts = [0, 0]
    for idx, _ in strip_to_mec8.values():
        jm_counts[idx] += 1
    print(f"  JM_1 (connector 0): {jm_counts[0]} strips")
    print(f"  JM_2 (connector 1): {jm_counts[1]} strips")
    return strip_to_mec8


# ── Fan-apex finder ───────────────────────────────────────────────────────────

def apex_from_strips(rows):
    """
    Find the convergence point (apex) of the radial strip fan using least-squares
    line intersection.

    Each strip defines a line through its connector-pad centre (pin_x, pin_y)
    and strip endpoint (x, y).  The apex is the point minimising the sum of
    squared perpendicular distances to all these lines.

    Parameters
    ----------
    rows : list of dicts with keys 'x', 'y', 'pin_x', 'pin_y'

    Returns
    -------
    apex : np.ndarray shape (2,)  — (x_apex, y_apex) in mm
    """
    valid = [r for r in rows if r['x'] != '' and r['x'] is not None]
    M = np.zeros((2, 2))
    b = np.zeros(2)
    for r in valid:
        p  = np.array([float(r['pin_x']), float(r['pin_y'])])
        q  = np.array([float(r['x']),     float(r['y'])])
        d  = q - p
        dn = np.linalg.norm(d)
        if dn < 1e-6:
            continue
        d = d / dn
        # projection matrix onto line normal: I - d d^T
        P = np.eye(2) - np.outer(d, d)
        M += P
        b += P @ p
    apex, _, _, _ = np.linalg.lstsq(M, b, rcond=None)
    return apex


def add_polar_coords(rows, apex):
    """
    Append radius, phi (radians), delta_phi (radians) to each row in-place.
    delta_phi = phi[strip_n] - phi[strip_{n-1}] within each connector
    (sorted by strip number); strip 1 gets delta_phi = 0.
    """
    for r in rows:
        if r['x'] == '' or r['x'] is None:
            r['radius'] = r['phi'] = r['delta_phi'] = ''
            continue
        dx = float(r['x']) - apex[0]
        dy = float(r['y']) - apex[1]
        r['radius'] = round(np.hypot(dx, dy), 4)
        r['phi']    = round(np.arctan2(dy, dx), 8)

    # delta_phi per connector, sorted by strip number
    from itertools import groupby
    key = lambda r: r['connector_number']
    for conn, grp in groupby(sorted(rows, key=key), key=key):
        strips = sorted([r for r in grp if r['phi'] != ''],
                        key=lambda r: r['strip'])
        for i, r in enumerate(strips):
            if i == 0:
                r['delta_phi'] = 0.0
            else:
                r['delta_phi'] = round(r['phi'] - strips[i - 1]['phi'], 8)


# ── Main pipeline ──────────────────────────────────────────────────────────────

def build_mapping(fcu_path, bcu_path, k59v_netlist=None,
                  min_span_mm=5.0, via_threshold=0.5):
    """
    Full pipeline: Gerber files (+ optional K59V netlist) → list of mapping row dicts.

    Row keys (always present):
      pad_number, x, y, connector_number, connector_pin, pin_name,
      sector, strip, channel_id, pin_x, pin_y, pad_angle

    Additional keys when k59v_netlist is provided:
      mec8_connector  — 0 (JM_1) or 1 (JM_2)
      mec8_pin        — physical MEC8 pin number (3–68)

    channel_id = sector * 128 + (strip - 1), unique 0-based index 0–1279.
    """
    print("INFO: Parsing F_Cu pad attributes …")
    signal_pads, net_to_via = parse_fcu_attributes(fcu_path)
    print(f"  → {len(signal_pads)} signal connector pads, "
          f"{len(net_to_via)} signal via positions")

    print("INFO: Matching B_Cu traces to via positions …")
    net_to_strip = find_strip_endpoints(bcu_path, net_to_via,
                                        min_span_mm=min_span_mm,
                                        via_threshold=via_threshold)

    strip_to_mec8 = {}
    if k59v_netlist:
        print("INFO: Parsing K59V adapter netlist …")
        strip_to_mec8 = parse_k59v_netlist(k59v_netlist)

    # Sort: connector J0…J9, then physical pin number ascending
    signal_pads.sort(key=lambda p: (int(p['ref'][1:]), p['pin']))

    rows = []
    for pad_number, pad in enumerate(signal_pads):
        net = pad['net']
        strip_pt = net_to_strip.get(net)
        sx = round(strip_pt[0], 4) if strip_pt else ''
        sy = round(strip_pt[1], 4) if strip_pt else ''

        row = {
            'pad_number':       pad_number,
            'x':                sx,
            'y':                sy,
            'connector_number': int(pad['ref'][1:]),
            'connector_pin':    pad['pin'],
            'pin_name':         pad['pinname'],
            'sector':           pad['sector'],
            'strip':            pad['strip'],
            'channel_id':       pad['sector'] * 128 + (pad['strip'] - 1),
            'pin_x':            round(pad['x'], 4),
            'pin_y':            round(pad['y'], 4),
            'pad_angle':        pad['angle'],
        }

        if strip_to_mec8:
            mec8 = strip_to_mec8.get(pad['strip'])
            row['mec8_connector'] = mec8[0] if mec8 else ''
            row['mec8_pin']       = mec8[1] if mec8 else ''

        rows.append(row)

    matched = sum(1 for r in rows if r['x'] != '')
    print(f"INFO: {matched}/{len(rows)} channels have valid strip coordinates")

    print("INFO: Computing fan apex and polar coordinates …")
    apex = apex_from_strips(rows)
    print(f"  Apex: ({apex[0]:.3f}, {apex[1]:.3f}) mm")
    add_polar_coords(rows, apex)

    return rows


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='P2 Basket Micromegas channel mapping from KiCad Gerbers')
    parser.add_argument('--fcu', required=True,
                        help='Path to P2_BASKET-F_Cu.gbr')
    parser.add_argument('--bcu', required=True,
                        help='Path to P2_BASKET-B_Cu.gbr')
    parser.add_argument('--k59v', default=None,
                        help='Path to Fx2Mec8.net (K59V adapter card netlist); '
                             'adds mec8_connector and mec8_pin columns')
    parser.add_argument('--out', default='P2_BASKET_mapping.csv',
                        help='Output CSV (default: P2_BASKET_mapping.csv)')
    parser.add_argument('--min_span_mm', type=float, default=5.0,
                        help='Min B_Cu trace span (mm) to consider as signal track '
                             '(filters copper-fill fragments; default 5)')
    parser.add_argument('--via_threshold', type=float, default=0.5,
                        help='Max distance (mm) from trace endpoint to via centre '
                             '(default 0.5)')
    args = parser.parse_args()

    rows = build_mapping(
        fcu_path=args.fcu,
        bcu_path=args.bcu,
        k59v_netlist=args.k59v,
        min_span_mm=args.min_span_mm,
        via_threshold=args.via_threshold,
    )

    fieldnames = ['pad_number', 'x', 'y',
                  'connector_number', 'connector_pin', 'pin_name',
                  'sector', 'strip', 'channel_id',
                  'pin_x', 'pin_y', 'pad_angle',
                  'radius', 'phi', 'delta_phi']
    if args.k59v:
        fieldnames += ['mec8_connector', 'mec8_pin']

    with open(args.out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"INFO: Written {len(rows)} rows → {args.out}")


if __name__ == '__main__':
    main()