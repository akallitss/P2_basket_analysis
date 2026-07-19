#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sps_config.py

Central configuration + run registry for the P2 (BASKET) SPS beam-test
analysis. Sibling of cosmic_bench_analysis/p2_qa_config.py, NOT a copy: the
pad-level infrastructure (p2_io, p2_mapping, p2_channel_qa, p2_sparks,
p2_waveforms) is *imported* from ../cosmic_bench_analysis via a sys.path shim
(setup_paths()), so a fix in the shared core lands in one place.

Difference from the bench registry
-----------------------------------
On the bench every registry entry is ONE detector-under-test referenced to the
M3 telescope. At the SPS there is no external tracker: the P2 telescope is a set
of N P2 stations (P2_IN / P2_MID / P2_OUT) that tag-and-probe each other. A
registry entry here therefore describes a WHOLE RUN, and the per-detector
information (name, FEUs, mesh spark channel, wiring, survey z) is read straight
from the run's run_config.json -- the number of stations is discovered from the
`included_detectors` wiring, never hardcoded. Telescope-wide stages (24 event
sync, 21 alignment) iterate every station; single-detector stages (20 spectra,
23 profile) and the probe of 22 write under that station's own folder.

Output layout (mirrors the bench so build_*_pdf.py + the DAQ GUI work as-is):
    <ANALYSIS_ROOT>/<det_tag>/<run>/<sub_run>/<stage>/...
det_tag is the station name (e.g. 'P2_OUT'); telescope-wide products use the
special tag 'telescope'.

Run registry
------------
Runs are registered in RUNS, keyed by a short name. Scripts pick a run with
`cfg = config_from_argv()` (first CLI arg = run key) or `get_config('<key>')`.
"""

import json
import os
import re
import sys

# --------------------------------------------------------------------------- #
# Path shim: pull the shared pad-level core from cosmic_bench_analysis
# --------------------------------------------------------------------------- #
_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(_HERE)                       # .../P2_basket_analysis
CORE_DIR = os.path.join(REPO_ROOT, 'cosmic_bench_analysis')


def setup_paths() -> None:
    """Put cosmic_bench_analysis/ (shared core) and this dir on sys.path.
    Call this at the top of every stage before importing p2_io / p2_mapping /
    p2_sparks etc. so the imported modules are the shared bench originals."""
    for p in (CORE_DIR, _HERE):
        if p not in sys.path:
            sys.path.insert(0, p)


setup_paths()

# --------------------------------------------------------------------------- #
# Paths — SITE-AWARE so the same code runs on the laptop and on the banco DAQ
# machine during the beam test (no file copying). Resolution order per root:
#   1. explicit env var  (SPS_DATA_ROOT / SPS_ANALYSIS_ROOT / SPS_COSMIC_BENCH_ROOT)
#   2. banco auto-detect (the DAQ machine: /local/home/banco exists)
#   3. laptop defaults
# On banco the beam DAQ writes runs under <base_data_dir>/runs/, so point
# SPS_DATA_ROOT at that 'runs' dir (e.g. /local/home/banco/P2_data/<beamtest>/runs)
# and select a run with the 'live' registry entry + SPS_RUN=<run_name>.
# --------------------------------------------------------------------------- #
_ON_BANCO = os.path.isdir('/local/home/banco')

def _root(env, banco_default, laptop_default):
    v = os.environ.get(env)
    if v:
        return v
    return banco_default if _ON_BANCO else laptop_default

# Where the SPS beam-test run data lives (mirrors the DAQ layout). Registry
# entries may override data_root (the laptop dry-run fixtures do).
DATA_ROOT = _root(
    'SPS_DATA_ROOT',
    '/local/home/banco/P2_data',
    '/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/SPS_Beam_Test')
# The cosmic-bench tree (Fe55 dry-run stand-in data on the laptop).
COSMIC_BENCH_ROOT = _root(
    'SPS_COSMIC_BENCH_ROOT',
    '/local/home/banco/P2_data',
    '/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/Cosmic_Bench')

# Analysis products (plots, csv, pdf).
ANALYSIS_ROOT = _root(
    'SPS_ANALYSIS_ROOT',
    '/local/home/banco/P2_data/Analysis',
    '/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/Analysis')

# Gerber-derived pad mapping + insulation-mask pillar Gerber (shared with the
# bench: same P2 BASKET PCB).
MAP_CSV_PATH = os.path.join(
    REPO_ROOT, 'Detector_Mapping', 'P2_BASKET', 'P2_BASKET_mapping.csv')
MASK_GBR_PATH = os.path.join(
    os.path.dirname(REPO_ROOT), 'Detector_Drawings', 'Version_Apr26',
    'Insulation_masks', 'V2', 'P2_BASKET-Mask_M2_V2.gbr')

DEFAULT_RUN = 'fe55_telescope'

# The special det_tag under which telescope-wide products (event sync,
# alignment) are written.
TELESCOPE_TAG = 'telescope'


def _slug(name: str) -> str:
    """Filesystem-safe station tag from a detector name ('P2_OUT' -> 'P2_OUT')."""
    return re.sub(r'[^A-Za-z0-9_.-]', '_', str(name))


# --------------------------------------------------------------------------- #
# Per-station info, built entirely from the run's run_config.json
# --------------------------------------------------------------------------- #
class DetInfo:
    """One P2 telescope station, derived from a run_config `detectors` entry.

    Everything (FEUs, mesh spark channel, survey z, wiring) comes from the
    config so the code adapts to N stations without hardcoding counts or FEU
    numbers.
    """

    def __init__(self, det_dict, det_tags=None):
        self.name = det_dict['name']
        # tag priority: registry det_tags override > run_config 'det_tag' > slug
        self.det_tag = (det_tags or {}).get(
            self.name, det_dict.get('det_tag') or _slug(self.name))
        self.det_type = det_dict.get('det_type', 'P2')
        dream = det_dict.get('dream_feus', {})
        self.feus = sorted({int(v[0]) for v in dream.values()})
        self.dream_feus = dream
        # mesh HV channel -> spark channel string 'card:chan' (p2_sparks fmt).
        mesh = det_dict.get('hv_channels', {}).get('mesh')
        self.mesh_card, self.mesh_chan = (
            (int(mesh[0]), int(mesh[1])) if mesh else (None, None))
        self.spark_channel = (f'{self.mesh_card}:{self.mesh_chan}'
                              if mesh else None)
        drift = det_dict.get('hv_channels', {}).get('drift')
        self.drift_channel = f'{int(drift[0])}:{int(drift[1])}' if drift else None
        c = det_dict.get('det_center_coords', {}) or {}
        self.x, self.y, self.z = (float(c.get('x', 0.0)),
                                  float(c.get('y', 0.0)),
                                  float(c.get('z', 0.0)))

    def __repr__(self):
        return (f'<DetInfo {self.name} tag={self.det_tag} FEUs={self.feus} '
                f'mesh={self.spark_channel} z={self.z}>')


# --------------------------------------------------------------------------- #
# Run config object
# --------------------------------------------------------------------------- #
class RunConfig:
    """A whole SPS beam-test run (all P2 stations of the telescope)."""

    def __init__(self, key, run, data_root=None, analysis_root=ANALYSIS_ROOT,
                 map_csv=MAP_CSV_PATH, mask_gbr=MASK_GBR_PATH,
                 ref_plane=None, min_tag=1, strategy='reverse',
                 det_tags=None, drop_connectors=None,
                 spark_imon_thr=2.0, spark_guard_before=2.0,
                 spark_guard_after=10.0, burst_npads=0,
                 min_amp=0.0, cluster_r=15.0, has_geometry=True, note=''):
        self.KEY = key
        self.RUN = run
        self.DATA_ROOT = data_root or DATA_ROOT
        self.ANALYSIS_ROOT = analysis_root
        self.MAP_CSV_PATH = map_csv
        self.MASK_GBR_PATH = mask_gbr
        self.STRATEGY = strategy
        # Reference plane name for alignment; None -> first station (highest z
        # or first in config). Others align to it.
        self.REF_PLANE = ref_plane
        # tag-and-probe: an event is TAGGED for a probe plane when >= MIN_TAG of
        # the OTHER planes carry a clean single cluster (use a majority, i.e.
        # >= ceil(n_other/2)+? , when >=3 planes exist -- resolved at run time).
        self.MIN_TAG = int(min_tag)
        # Optional {det_name: det_tag} override (e.g. {'P2_OUT':'det2'}).
        self._det_tags = det_tags or {}
        # Optional {det_name: (connectors,)} dead-connector drops per station.
        self._drop_connectors = drop_connectors or {}
        # HV spark veto knobs (shared cfg surface with p2_sparks.SparkVeto).
        self.SPARK_IMON_THR = float(spark_imon_thr)
        self.SPARK_GUARD_BEFORE = float(spark_guard_before)
        self.SPARK_GUARD_AFTER = float(spark_guard_after)
        # Beam self-/scintillator-trigger events are naturally high-multiplicity;
        # the cosmic burst veto (default 20) would eat real events, so default 0.
        self.BURST_NPADS = int(burst_npads)
        self.MIN_AMP = float(min_amp)
        self.CLUSTER_R = float(cluster_r)
        # False when the (feu, channel) -> pad wiring / pad-position map is not
        # available for this detector (e.g. a large multi-FEU detector with no
        # dream_feus table): geometry stages fall back to channel space.
        self.HAS_GEOMETRY = bool(has_geometry)
        self.NOTE = note
        self._cfg_cache = None
        self._det_cache = None

    # -- run-level paths ---------------------------------------------------- #
    @property
    def run_dir(self):
        return os.path.join(self.DATA_ROOT, self.RUN)

    @property
    def run_config_path(self):
        return os.path.join(self.run_dir, 'run_config.json')

    def subrun_dir(self, sub_run):
        return os.path.join(self.run_dir, sub_run)

    def combined_hits_dir(self, sub_run):
        return os.path.join(self.subrun_dir(sub_run), 'combined_hits_root')

    def hv_monitor_csv(self, sub_run):
        return os.path.join(self.subrun_dir(sub_run), 'hv_monitor.csv')

    # -- run_config.json ---------------------------------------------------- #
    def _load_run_config(self):
        if self._cfg_cache is None:
            with open(self.run_config_path) as fh:
                self._cfg_cache = json.load(fh)
        return self._cfg_cache

    def detectors(self):
        """List of DetInfo for every included station, in config order."""
        if self._det_cache is None:
            cfg = self._load_run_config()
            included = cfg.get('included_detectors')
            dets = cfg.get('detectors', [])
            if included:
                order = {n: i for i, n in enumerate(included)}
                dets = [d for d in dets if d.get('name') in order]
                dets.sort(key=lambda d: order[d['name']])
            self._det_cache = [DetInfo(d, self._det_tags) for d in dets]
        return self._det_cache

    def det_by_tag(self, tag):
        for d in self.detectors():
            if d.det_tag == tag:
                return d
        raise KeyError(f'No station with det_tag={tag!r} in {self.RUN}')

    def ref_det(self):
        """The reference station (REF_PLANE by name, else the first station)."""
        dets = self.detectors()
        if self.REF_PLANE:
            for d in dets:
                if d.name == self.REF_PLANE:
                    return d
        return dets[0]

    def all_feus(self):
        feus = set()
        for d in self.detectors():
            feus.update(d.feus)
        return sorted(feus)

    def drop_connectors_for(self, det):
        return tuple(self._drop_connectors.get(det.name, ()))

    # -- sub_runs ----------------------------------------------------------- #
    def find_subruns(self):
        """On-disk sub_run directories that carry combined-hits ROOT files,
        in directory-name order (which is the HV/time order for scan runs)."""
        import glob
        out = []
        if not os.path.isdir(self.run_dir):
            return out
        for name in sorted(os.listdir(self.run_dir)):
            d = self.subrun_dir(name)
            if not os.path.isdir(d):
                continue
            if glob.glob(os.path.join(d, 'combined_hits_root', '*.root')):
                out.append(name)
        return out

    def subrun_mesh_hv(self, sub_run, det):
        """Mesh HV [V] this station was set to in `sub_run`, from run_config
        `sub_runs` (hvs[card][chan]); None if not tabulated."""
        if det.mesh_card is None:
            return None
        cfg = self._load_run_config()
        for s in cfg.get('sub_runs', []):
            if s.get('sub_run_name') == sub_run:
                hvs = s.get('hvs', {})
                v = hvs.get(str(det.mesh_card), {}).get(str(det.mesh_chan))
                return int(v) if v is not None else None
        return None

    # -- channel tables ----------------------------------------------------- #
    def channel_table(self, det):
        """Pad channel table for one station (its FEUs only), via the shared
        p2_mapping. Cached on the DetInfo."""
        import p2_mapping as pmap
        if getattr(det, '_ct', None) is None:
            det._ct = pmap.build_channel_table(
                self.run_config_path, self.MAP_CSV_PATH,
                det_type=det.det_type, det_name=det.name,
                strategy=self.STRATEGY,
                drop_connectors=self.drop_connectors_for(det))
        return det._ct

    # -- output ------------------------------------------------------------- #
    def out_dir(self, det_tag, sub_run, *stage_parts):
        d = os.path.join(self.ANALYSIS_ROOT, det_tag, self.RUN, sub_run,
                         *stage_parts)
        os.makedirs(d, exist_ok=True)
        return d

    SPARK_SUFFIX = '_spark_vetoed'

    def product_suffix(self, veto_sparks=False):
        return self.SPARK_SUFFIX if veto_sparks else ''

    def alignment_json(self, sub_run):
        """Canonical path of the per-sub_run alignment JSON that 21 writes and
        22 consumes (telescope-wide product)."""
        d = self.out_dir(TELESCOPE_TAG, sub_run, '21_telescope_align')
        return os.path.join(d, 'alignment.json')

    def __repr__(self):
        return (f'<SPS RunConfig {self.KEY}: {self.RUN} '
                f'-> {os.path.join(self.ANALYSIS_ROOT, "<det>", self.RUN)}>')


# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #
RUNS = {
    # Dry-run stand-in for the SPS P2 telescope: the 7-18-26 banco Fe55 mesh-HV
    # scan of the two telescope stations P2_OUT (FEU 3) and P2_MID (FEU 4).
    # Same combined-hits format and two-station wiring the beam telescope will
    # have, so every stage can be exercised end-to-end today. The Fe55 source
    # sat near P2_OUT, so P2_MID sees far fewer events -> the two-plane
    # correlation stages (21/22) are statistics-limited here (expected; real
    # 3-plane beam data is what makes them meaningful).
    'fe55_telescope': RunConfig(
        key='fe55_telescope',
        run='p2_fe55_det2_det3_mesh_scan_7-18-26',
        data_root=COSMIC_BENCH_ROOT,
        ref_plane='P2_OUT',        # highest-stats plane = alignment reference
        det_tags={'P2_OUT': 'det2', 'P2_MID': 'det3'},
        min_tag=1,                 # only one other plane exists (2 stations)
        burst_npads=0,             # self-trigger source: high multiplicity
        spark_imon_thr=2.0,
        note='Fe55 telescope-scan dry-run (2 stations: P2_OUT FEU3, P2_MID FEU4)'),

    # November 2025 SPS test beam (TbSPS25), run_19: ONE large P2 detector
    # 'P2_SPS25' read out across 5 FEUs (cfg Feu 1-5 = RunCtrl Ids 69/70/101/
    # 102/103), external TCM/scintillator trigger, ZS, 16 samples, 150/80 GeV
    # mu/pi beam. There is NO run_config.json from the DAQ and the (feu,
    # channel) -> pad wiring for this detector is unknown, so HAS_GEOMETRY is
    # False: the geometry stages (20/23) run in channel space, and stage 24
    # (event sync) -- which only needs the FEU list -- runs fully. A minimal
    # run_config.json enumerating the 5 FEUs is written into the run dir so
    # RunConfig.feus = [1, 2, 3, 4, 5].
    'run19_sps25': RunConfig(
        key='run19_sps25',
        run='run_19',
        data_root='/local/home/ak271430/Documents/PostDocSaclay/data/'
                  'nov25_run19_test',
        has_geometry=True,
        burst_npads=0,
        note='Nov-2025 SPS TbSPS25 run_19. Large P2 detector P2_SPS25 = physical '
             'connectors 1-4 on FEU 5 (full 512 ch = channel_id 0-511, sectors '
             '0-3), covered by P2_BASKET_mapping.csv. FEUs 1-4 read other, small '
             'P2 detectors (not modelled). Wiring from Alexandra 2026-07-19.'),
}


# 'live' — the beam-test workhorse on banco: select any run the DAQ has written
# under DATA_ROOT by name, e.g.  SPS_DATA_ROOT=/local/home/banco/P2_data/<bt>/runs
# SPS_RUN=run_5 python 24_event_sync_qa.py live . Detector wiring comes from that
# run's DAQ-written run_config.json, so no per-run code change is needed. Only
# added when SPS_RUN is set, so importing the module is always safe.
_live_run = os.environ.get('SPS_RUN')
if _live_run:
    try:
        RUNS['live'] = RunConfig(
            key='live', run=_live_run, data_root=DATA_ROOT, has_geometry=True,
            note=f'Live beam-test run {_live_run!r} under {DATA_ROOT} (SPS_RUN env)')
    except Exception as _e:
        print(f'[sps_config] could not register live run {_live_run!r}: {_e}')


def get_config(key=None) -> RunConfig:
    key = key or DEFAULT_RUN
    if key not in RUNS:
        raise KeyError(f'Unknown run key {key!r}. Known: {sorted(RUNS)}')
    return RUNS[key]


def config_from_argv() -> RunConfig:
    """First CLI arg (if present and not an option) selects the run."""
    key = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith('-') \
        else DEFAULT_RUN
    return get_config(key)


if __name__ == '__main__':
    for k in RUNS:
        c = get_config(k)
        print(c)
        print('  run_config :', c.run_config_path,
              '(exists)' if os.path.isfile(c.run_config_path) else '(MISSING)')
        try:
            for d in c.detectors():
                print('   ', d, '  hits:',
                      'ok' if os.path.isdir(c.combined_hits_dir(
                          c.find_subruns()[0])) else '—')
            print('  sub_runs   :', len(c.find_subruns()),
                  '->', c.find_subruns()[:3], '...')
        except Exception as e:  # noqa: BLE001
            print('  [could not read run_config detectors]:', e)
