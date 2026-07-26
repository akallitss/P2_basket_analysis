#!/usr/bin/env python3
"""
Geometry + clustering for the two EIC uRWELL references of TB_July2026_H4.

Geometry comes from the strip maps deployed 2026-07-25 (see
~/P2_basket_analysis/sps_beam_analysis/urw_reference/URW_TRACKING_HANDOFF_2026-07-25.md §4) via DetectorConfigLoader /
DreamDetector.  Hits come from the C++ matched-filter analyzer output
(`<sub_run>/hits_root/*_01_hits.root`), NOT from DreamData: the decoded_root
files here are zero-suppressed (vector<uint16_t> sample/channel/amplitude) and
carry no `_array` flag, so DreamData.read_data() finds no files and its
split_det_data() reshape does not apply.

Conventions (handoff §4):
  * a map row's `axis` is the direction the strip RUNS, so axis='y' strips
    measure x and axis='x' strips measure y;
  * global channel = (feu_connector - 1) * 64 + connectorChannel;
  * strip position is the Gerber coordinate, already shifted so min -> 0 by the
    loader, i.e. local detector coordinates in [0, active_size].
"""

import os
import sys
import glob

import numpy as np
import pandas as pd

sys.path.insert(0, '/local/home/banco/dylan/saclay_micromegas')
from Detector_Classes.DetectorConfigLoader import DetectorConfigLoader  # noqa: E402
from Detector_Classes.DreamDetector import DreamDetector  # noqa: E402

DET_DIR_DEFAULT = '/local/home/banco/P2_data/TB_July2026_H4/config/detectors/'
URW_DETS = ('EIC_uRWELL_front', 'EIC_uRWELL_back')
N_CHAN = 512

# How each view's 128 strip-map entries are laid onto its 128 Dream channels.
# A view is read out on two 64-channel connectors, so there are two binary
# choices and four candidate wirings:
#
#   'AB'      map connector 0 -> the lower FEU connector, channel order as-is
#   'BA'      the two connectors interchanged
#   'AB_rev'  as AB, but channel order REVERSED inside each connector
#   'BA_rev'  as BA, but channel order reversed inside each connector
#
# 'BA' and 'AB_rev' (and likewise 'AB'/'BA_rev') differ only by a mirror of the
# whole view WHEN the map is uniform, which is why the July 2026 tests - all of
# which are mirror-blind - could pin the pair order but never the within-
# connector reversal (ORDERING.md, "Orientation vs pair order").  The back's
# map is NOT uniform (three pitch zones per view: 1.0 / 1.5 / 0.5 mm), so for
# the back the four candidates are genuinely distinct and the degeneracy breaks.
#
# Measured 2026-07-26 on highstat_eff_1/beam_commissioning_00 by
# explore4_back_map.py.  The front points at the three P2 stations to better
# than ~1 mm and the beam is parallel (the best front->back interpolation
# fraction for predicting a P2 plane is 0.000, i.e. the "track angle" carries no
# information), so the back must reproduce the front position up to a constant.
# Width of (back - front) per candidate:
#
#   view      AB       BA       AB_rev    BA_rev
#   back x   45.6 mm   5.42 mm   0.80 mm  47.0 mm
#   back y   35.2 mm   6.40 mm   0.88 mm  38.8 mm
#
# So the back is 'AB_rev' on both views: its connectors are NOT interchanged,
# the channel order inside each is reversed.  Using 'BA' (this file until
# 2026-07-26) left the back pointing at 4.4 mm and made it useless for tracking.
#
# The front is 'BA'/'AB' here.  Its own data cannot separate that from the
# mirror partner 'AB_rev'/'BA_rev' (they agree to 0.1 mm), but the choice below
# is the one that maps the front into the P2 BASKET pad frame by a PROPER
# rotation (det +0.99, -60 deg); the mirror partner would need a reflection,
# which two detectors viewed from the same side along the beam cannot need.
#
# This makes the run config's `dream_feu_orientation: inverted`, corrected on
# all eight uRWELL connectors 2026-07-25, a real and consistent statement: every
# connector carries a within-connector reversal, and on top of that exactly one
# view (front y) has its connector pair interchanged.  The note previously here
# - that a uniform 'inverted' could not explain the measured pattern - was
# wrong because it only considered the pair order.
VIEW_MODE_DEFAULT = {
    'EIC_uRWELL_front': {'x': 'BA', 'y': 'AB'},
    'EIC_uRWELL_back': {'x': 'AB_rev', 'y': 'AB_rev'},
}

# Pure labelling convention applied AFTER the wiring above: mirror a view's
# local coordinate (pos -> min + max - pos) so that both uRWELLs measure in the
# same direction.  With the correct wiring the back reads anti-parallel to the
# front on both views (fitted front->back slope -0.993 in x, -1.029 in y), which
# is physical but makes every downstream "slope should be +1" check awkward.
# Flipping the back here costs nothing - the strip -> position assignment is
# untouched, only the sign of the axis label - and leaves front -> back at +1.
AXIS_FLIP_DEFAULT = {
    'EIC_uRWELL_back': {'x': True, 'y': True},
}


class UrwGeometry:
    """Channel -> (view, local position, pitch, zone) lookup for one uRWELL.

    Attributes
    ----------
    view   : (512,) 'x', 'y' or '' for channels not belonging to this detector
    pos    : (512,) local coordinate of the strip [mm], NaN if unmapped
    pitch  : (512,) strip pitch of the zone the strip lives in [mm]
    zone   : (512,) det_map group index, -1 if unmapped
    """

    def __init__(self, name, run_json, det_dir=DET_DIR_DEFAULT, sub_run_name=None,
                 view_mode=None, axis_flip=None):
        loader = DetectorConfigLoader(run_json, det_dir)
        cfg = loader.get_det_config(name, sub_run_name=sub_run_name)
        if 'det_map' not in cfg:
            raise RuntimeError(f'{name}: no det_map bound (det_type={cfg.get("det_type")!r}) '
                               f'- check the strip map deployment')
        self.name = name
        self.det_type = cfg['det_type']
        self.det = DreamDetector(config=cfg)
        self.feu_num = self.det.feu_num
        self.feu_connectors = list(self.det.feu_connectors)
        self.active_size = np.asarray(self.det.active_size, float)
        self.center = np.asarray(self.det.center, float)

        self.view = np.full(N_CHAN, '', dtype='<U1')
        self.pos = np.full(N_CHAN, np.nan)
        self.pitch = np.full(N_CHAN, np.nan)
        self.interpitch = np.full(N_CHAN, np.nan)
        self.zone = np.full(N_CHAN, -1, dtype=np.int16)

        for gi, row in self.det.det_map.iterrows():
            # axis is the direction the strip runs -> it measures the other axis
            meas = 'x' if row['axis'] == 'y' else 'y'
            gerber = row['xs_gerber'] if meas == 'x' else row['ys_gerber']
            conns = np.asarray(row['connectors'], int)
            chans = np.asarray(row['channels'], int)
            gch = (conns - 1) * 64 + chans
            self.view[gch] = meas
            self.pos[gch] = gerber
            self.pitch[gch] = float(row['pitch(mm)'])
            self.zone[gch] = gi
            ip = str(row['interpitch(mm)'])
            self.interpitch[gch] = float(ip.split(':')[0])

        self.view_mode = dict(VIEW_MODE_DEFAULT.get(name, {'x': 'AB', 'y': 'AB'})
                              if view_mode is None else view_mode)
        for view, mode in self.view_mode.items():
            self._apply_view_mode(view, mode)

        self.axis_flip = dict(AXIS_FLIP_DEFAULT.get(name, {'x': False, 'y': False})
                              if axis_flip is None else axis_flip)
        for view, do_flip in self.axis_flip.items():
            if do_flip:
                m = self.view == view
                self.pos[m] = np.nanmin(self.pos[m]) + np.nanmax(self.pos[m]) - self.pos[m]

        self.mapped = self.zone >= 0
        if self.mapped.sum() != 256:
            raise RuntimeError(f'{name}: {self.mapped.sum()} mapped channels, expected 256')

    def _view_connectors(self, view):
        conns = sorted(c for c in self.feu_connectors
                       if (self.view[(c - 1) * 64:c * 64] == view).any())
        if len(conns) != 2:
            raise RuntimeError(f'{self.name}: view {view!r} spans {conns}, '
                               f'expected exactly 2 connectors')
        return conns

    def _apply_view_mode(self, view, mode):
        """Re-lay a view's 128 map entries onto its 128 channels (VIEW_MODE_DEFAULT).

        The loader assigns the map's connector column to FEU connectors in
        order; `mode` says how the cabling actually differs from that - which
        of the two connectors carries the low strips ('AB' vs 'BA') and whether
        the channel order runs backwards inside each connector ('_rev').
        """
        if mode not in ('AB', 'BA', 'AB_rev', 'BA_rev'):
            raise ValueError(f'{self.name}: unknown view mode {mode!r}')
        if mode == 'AB':
            return
        lo, hi = self._view_connectors(view)
        a = slice((lo - 1) * 64, lo * 64)
        b = slice((hi - 1) * 64, hi * 64)
        for arr in (self.pos, self.pitch, self.interpitch, self.zone):
            blk = [arr[a].copy(), arr[b].copy()]
            if mode.endswith('_rev'):
                blk = [x[::-1] for x in blk]
            if mode.startswith('BA'):
                blk = blk[::-1]
            arr[a], arr[b] = blk[0], blk[1]

    def zone_table(self):
        rows = []
        for gi, row in self.det.det_map.iterrows():
            m = self.zone == gi
            meas = 'x' if row['axis'] == 'y' else 'y'
            # report the FEU connectors the zone actually lands on AFTER the
            # wiring of VIEW_MODE_DEFAULT - row['connectors'] is the map's own
            # column and is misleading for any view that is not plain 'AB'.
            feu_conns = sorted(set((np.flatnonzero(m) // 64 + 1).tolist()))
            rows.append(dict(zone=gi, view=meas, pitch=float(row['pitch(mm)']),
                             interpitch=str(row['interpitch(mm)']), n_strips=int(m.sum()),
                             feu_connectors=feu_conns,
                             pos_min=float(np.nanmin(self.pos[m])),
                             pos_max=float(np.nanmax(self.pos[m]))))
        return pd.DataFrame(rows)

    def to_global(self, x_local, y_local):
        """Local (x, y) [mm] -> global (x, y, z) [mm] via Detector.convert_coords_to_global."""
        x_local = np.atleast_1d(np.asarray(x_local, float))
        y_local = np.atleast_1d(np.asarray(y_local, float))
        out = np.empty((len(x_local), 3))
        for i, (xi, yi) in enumerate(zip(x_local, y_local)):
            out[i] = self.det.convert_coords_to_global(np.array([xi, yi, 0.0]))
        return out

    def __repr__(self):
        return (f'<UrwGeometry {self.name} ({self.det_type}) feu={self.feu_num} '
                f'conn={self.feu_connectors} active={np.round(self.active_size[:2], 3)} '
                f'z={self.center[2]}>')


# ---------------------------------------------------------------- hits I/O ---

HIT_BRANCHES = ['eventId', 'channel', 'amplitude', 'time_of_max',
                'significance', 'saturated', 'trigger_timestamp_ns']


def feu_hit_files(sub_run_dir, feu=1):
    """Hits files of a sub-run holding FEU <feu>, sorted by chunk number.

    Two layouts exist.  The original one is a file per FEU, whose name ends in
    the FEU id, so FEU 1 is `*_01_hits.root` (handoff §5.3).  Those files have
    since been deleted from this campaign - no `hits_root/` directory survives
    anywhere under P2_data/TB_July2026_H4/runs - so we fall back to the event
    synchronised `combined_hits_root/` files, which hold every FEU in one tree
    with a `feu` branch to tell them apart (handoff §5.4).  For FEU 1 the two
    are the same hits: the handoff records feu==1 of the combined file as
    byte-for-byte the `*_01_hits.root` content.

    When the fallback is used the caller must filter on `feu`; `iter_hits`
    does that for you if you pass its `feu=` argument.
    """
    pat = os.path.join(sub_run_dir, 'hits_root', f'*_{feu:02d}_hits.root')
    files = sorted(glob.glob(pat))
    if files:
        return files
    return sorted(glob.glob(os.path.join(sub_run_dir, 'combined_hits_root',
                                         '*_feu-combined_hits.root')))


def iter_hits(files, branches=HIT_BRANCHES, step=2_000_000, max_hits=None,
              progress=True, feu=None):
    """Stream hits, never splitting an eventId across yields.

    uproot's `f['hits']` resolves to the highest cycle, which is the newest
    reprocessing - do not sum cycles (handoff §7).

    `feu` selects one FEU when the files are `combined_hits_root` ones, which
    hold all four FEUs in a single tree ordered by FEU, not by event (handoff
    §5.4).  The filter is applied before anything else, so the eventIds seen
    downstream are monotonic again and `max_hits` counts hits of the FEU you
    asked for rather than hits of all four.  It is ignored for per-FEU files,
    which have no `feu` branch.
    """
    import uproot
    carry = None
    n_read = 0
    for path in files:
        tree = uproot.open(path)['hits']
        want = list(branches)
        use_feu = feu is not None and 'feu' in tree
        if use_feu and 'feu' not in want:
            want.append('feu')
        if progress:
            print(f'    {os.path.basename(path)}: {tree.num_entries} hits', flush=True)
        for chunk in tree.iterate(want, library='np', step_size=step):
            if max_hits is not None and n_read >= max_hits:
                carry = None
                break
            if use_feu:
                sel = chunk['feu'] == feu
                chunk = {k: v[sel] for k, v in chunk.items() if k in branches}
            if carry is not None:
                chunk = {k: np.concatenate([carry[k], chunk[k]]) for k in chunk}
                carry = None
            ev = chunk['eventId']
            if not len(ev):
                # the whole step fell outside the requested FEU's block
                continue
            n_read += len(ev)
            # hold back the trailing (possibly incomplete) event
            last = ev[-1]
            keep = ev != last
            if keep.all():
                yield chunk
            else:
                carry = {k: v[~keep] for k, v in chunk.items()}
                if keep.any():
                    yield {k: v[keep] for k, v in chunk.items()}
        if carry is not None:
            yield carry
            carry = None


# ------------------------------------------------------------- clustering ---

def cluster_hits(chunk, geo, min_amp=0.0, min_signif=0.0, drop_saturated=False,
                 gap_factor=1.05):
    """Amplitude-weighted 1D clusters, per event and per view.

    Two strips join the same cluster when their position gap is
    <= gap_factor * max(pitch_i, pitch_j).  With gap_factor=1.05 this merges
    across the pitch-zone boundaries (gaps 1.00-1.375 mm, handoff §4) the way
    the C++ Clusterizer1D does, while still resolving neighbouring strips inside
    the 1.5 mm zone - a flat 2.0 mm threshold cannot do both.

    Returns a DataFrame with one row per cluster:
      eventId, view, pos, size, charge, max_amp, time, zone, zone_mixed, t_ns
    """
    ev = chunk['eventId'].astype(np.int64)
    ch = chunk['channel'].astype(np.int32)
    amp = chunk['amplitude'].astype(np.float64)

    # amplitude is a fitted quantity and can come out <= 0; those would give a
    # zero-charge cluster and a NaN centroid, so drop them unconditionally.
    sel = geo.mapped[ch] & (amp > 0)
    if min_amp > 0:
        sel &= amp >= min_amp
    if min_signif > 0:
        sel &= chunk['significance'] >= min_signif
    if drop_saturated:
        sel &= ~chunk['saturated'].astype(bool)
    if not sel.any():
        return _empty_clusters()

    ev, ch, amp = ev[sel], ch[sel], amp[sel]
    tmax = chunk['time_of_max'][sel].astype(np.float64)
    tns = chunk['trigger_timestamp_ns'][sel].astype(np.float64)
    pos = geo.pos[ch]
    pitch = geo.pitch[ch]
    zone = geo.zone[ch]
    vcode = (geo.view[ch] == 'y').astype(np.int8)  # 0 = x view, 1 = y view

    order = np.lexsort((pos, vcode, ev))
    ev, ch, amp, tmax, tns = ev[order], ch[order], amp[order], tmax[order], tns[order]
    pos, pitch, zone, vcode = pos[order], pitch[order], zone[order], vcode[order]

    new = np.ones(len(ev), bool)
    if len(ev) > 1:
        gap = np.diff(pos)
        tol = gap_factor * np.maximum(pitch[:-1], pitch[1:])
        new[1:] = (np.diff(ev) != 0) | (np.diff(vcode) != 0) | (gap > tol)
    starts = np.flatnonzero(new)

    q = np.add.reduceat(amp, starts)
    centroid = np.add.reduceat(amp * pos, starts) / q
    size = np.diff(np.append(starts, len(ev)))
    tmean = np.add.reduceat(amp * tmax, starts) / q
    zmin = np.minimum.reduceat(zone, starts)
    zmax = np.maximum.reduceat(zone, starts)
    amax = np.maximum.reduceat(amp, starts)

    return pd.DataFrame(dict(
        eventId=ev[starts], view=np.where(vcode[starts] == 1, 'y', 'x'),
        pos=centroid, size=size.astype(np.int32), charge=q, max_amp=amax,
        time=tmean, zone=zmin.astype(np.int16), zone_mixed=(zmin != zmax),
        t_ns=tns[starts]))


def _empty_clusters():
    return pd.DataFrame(dict(eventId=np.array([], np.int64), view=np.array([], '<U1'),
                             pos=np.array([], float), size=np.array([], np.int32),
                             charge=np.array([], float), max_amp=np.array([], float),
                             time=np.array([], float), zone=np.array([], np.int16),
                             zone_mixed=np.array([], bool), t_ns=np.array([], float)))


def leading_points(clusters):
    """One (x, y) point per event from the highest-charge cluster in each view."""
    if not len(clusters):
        return pd.DataFrame(columns=['eventId', 'x', 'y', 'x_size', 'y_size',
                                     'x_charge', 'y_charge', 'x_time', 'y_time',
                                     'n_x', 'n_y', 't_ns'])
    c = clusters.sort_values('charge', ascending=False)
    lead = c.groupby(['eventId', 'view'], sort=False).first()
    n = c.groupby(['eventId', 'view'], sort=False).size().rename('n')
    lead = lead.join(n).reset_index()
    wide = lead.pivot(index='eventId', columns='view')
    out = pd.DataFrame(index=wide.index)
    for view in ('x', 'y'):
        for col, name in (('pos', ''), ('size', '_size'), ('charge', '_charge'),
                          ('time', '_time'), ('n', 'n_')):
            key = (col, view)
            tgt = f'n_{view}' if name == 'n_' else f'{view}{name}'
            out[tgt] = wide[key] if key in wide.columns else np.nan
    out['t_ns'] = wide[('t_ns', 'x')] if ('t_ns', 'x') in wide.columns else np.nan
    out['t_ns'] = out['t_ns'].fillna(wide[('t_ns', 'y')] if ('t_ns', 'y') in wide.columns else np.nan)
    return out.reset_index()
