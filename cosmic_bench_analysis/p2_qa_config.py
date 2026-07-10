#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_qa_config.py

Central configuration + run registry for the P2 (BASKET) cosmic-bench QA.
Structured after nTof_x17/mx_june_cosmic_qa/qa_config.py, but adapted for the
P2 detector, which is a *metallic, non-resistive* Micromegas read out on
**pads** (not resistive strips). There is therefore no strip/micro-TPC analysis
here; the spatial unit is the pad, positioned from the Gerber-derived mapping
(Detector_Mapping/P2_BASKET/P2_BASKET_mapping.csv).

Run registry
------------
Runs are registered in RUNS, keyed by a short name. Scripts pick a run with
`cfg = config_from_argv()` (first CLI arg = run key, default DEFAULT_RUN) or
`get_config('<key>')`.

Output layout
-------------
Every QA product is written under an Analysis/ tree that lives next to the data
(kept separate from the code), organised **by detector number** so each detector
under test gets its own folder:

    <Cosmic_Bench>/Analysis/<detN>/<run>/<sub_run>/<stage>/...

The detector number (detN) is taken from the run name (e.g. 'p2_det1_...' -> 'det1').
"""

import json
import os
import re
import sys

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(_HERE)  # .../P2_basket_analysis

# Where the fetched bench data lives (pull_run mirrors the DAQ layout here).
DATA_ROOT = '/local/home/ak271430/Documents/PostDocSaclay/data/Cosmic_Bench'

# Gerber-derived pad mapping (channel_id -> pad geometry).
MAP_CSV_PATH = os.path.join(
    REPO_ROOT, 'Detector_Mapping', 'P2_BASKET', 'P2_BASKET_mapping.csv')

DEFAULT_RUN = 'det1_long'

# --------------------------------------------------------------------------- #
# M3 reference-tracking consumer cut (tracking v2)
# --------------------------------------------------------------------------- #
# Recommended recipe for a "good" M3 ray after the tracking-v2 reconstruction:
#     NClusX >= M3_MIN_NCLUS  &&  NClusY >= M3_MIN_NCLUS  &&  Chi2X,Chi2Y < M3_CHI2_CUT
# The v2 reco changed the chi2 scale, so the old thresholds (1.5 / 20) are stale.
# The NClus requirement is essential: a 2-point-per-coordinate fit is exact, so it
# carries a denormal-tiny (not exactly zero) chi2 and slips through a naive chi2
# cut -- NClus>=3 forces a genuine 3-4 layer fit. These match the M3RefTracking
# defaults; we pass them explicitly so the recipe does not silently ride on the
# reader's defaults.
M3_CHI2_CUT = 5.0
M3_MIN_NCLUS = 3


def setup_paths() -> None:
    """Put the repo root and cosmic_bench_analysis/ on sys.path for imports."""
    for p in (REPO_ROOT, _HERE):
        if p not in sys.path:
            sys.path.insert(0, p)


def det_tag_from_run(run_name: str) -> str:
    """'p2_det1_long_run_6-30-26' -> 'det1'.  Falls back to 'detX'."""
    m = re.search(r'det(\d+)', run_name)
    return f'det{m.group(1)}' if m else 'detX'


# --------------------------------------------------------------------------- #
# Config object
# --------------------------------------------------------------------------- #
class _Config:
    def __init__(self, key, run, sub_run, det_name='P2_1',
                 det_type='P2', ref_det_type='m3', data_root=DATA_ROOT,
                 map_csv=MAP_CSV_PATH, dead_connectors=(),
                 spark_channel='1:0', spark_imon_thr=2.0,
                 spark_guard_before=2.0, spark_guard_after=10.0,
                 burst_npads=20, det_tag=None, match_r=20.0, plane_z=None):
        self.KEY = key
        self.DATA_ROOT = data_root
        self.RUN = run
        self.SUB_RUN = sub_run
        self.DET_NAME = det_name
        self.DET_TYPE = det_type            # detector-under-test det_type in run_config
        self.REF_DET_TYPE = ref_det_type    # reference tracker det_type (m3 telescope)
        self.MAP_CSV_PATH = map_csv
        # Explicit det_tag override for runs whose name contains several
        # detector tags (e.g. 'p2_det1_det2_...'), where the regex would
        # always pick the first one.
        self.DET_TAG = det_tag or det_tag_from_run(run)
        # Physical connectors (1..10) that are disconnected/dead on this
        # detector. They are dropped from the pad map so they do not show up as
        # spurious dead regions or bias the efficiency (a track pointing at a
        # dead connector is not a real 'miss').
        self.DEAD_CONNECTORS = tuple(dead_connectors)
        # --- HV spark flagging (hv_monitor.csv) ---------------------------- #
        # The mesh HV channel discharges (sparks) as brief imon spikes. Events
        # taken during a spark + its recovery are vetoed from every stage.
        self.SPARK_CHANNEL = spark_channel       # CAEN 'board:chan' of the mesh
        self.SPARK_IMON_THR = spark_imon_thr     # imon >= this [uA] is a spark
        self.SPARK_GUARD_BEFORE = spark_guard_before   # veto pad before [s]
        self.SPARK_GUARD_AFTER = spark_guard_after     # veto pad after [s]
        # P2-internal discharges the mesh-current veto misses show up as burst
        # events where tens of pads fire at once (real muon clusters are 1-3
        # pads, median 1). Events with >= this many pads are vetoed together
        # with the HV sparks. Set to 0/None to disable.
        self.BURST_NPADS = burst_npads
        # Track-hit match radius [mm] used by the efficiency stages (06/11/12).
        # Pick it on the plateau of 12_validation's eff-vs-R knob test.
        self.MATCH_R = float(match_r)
        # Measured detector plane z [mm] (03_m3_alignment z-scan). When set it
        # overrides the run_config det_center z in det_plane_z() — the fitted
        # plane wins over the nominal one if they disagree.
        self.PLANE_Z = plane_z

        # Analysis/ tree, keyed by detector -> run -> sub_run.
        self.OUT_BASE = os.path.join(data_root, 'Analysis',
                                     self.DET_TAG, run, sub_run)

    # -- input paths (mirror the DAQ layout under DATA_ROOT) ---------------- #
    @property
    def run_dir(self):
        return os.path.join(self.DATA_ROOT, self.RUN)

    @property
    def sub_run_dir(self):
        return os.path.join(self.run_dir, self.SUB_RUN)

    @property
    def run_config_path(self):
        return os.path.join(self.run_dir, 'run_config.json')

    def subrun_dir(self, sub_run):
        """Path to a specific sub_run under this run (for multi-subrun scans)."""
        return os.path.join(self.run_dir, sub_run)

    @property
    def combined_hits_dir(self):
        return os.path.join(self.sub_run_dir, 'combined_hits_root')

    @property
    def hits_root_dir(self):
        """Per-FEU hits_root files (pre-combination). These carry the
        'pedestals' TTree (per-channel mean/rms) that the combined files drop."""
        return os.path.join(self.sub_run_dir, 'hits_root')

    @property
    def m3_tracking_dir(self):
        return os.path.join(self.sub_run_dir, 'm3_tracking_root')

    @property
    def hv_monitor_csv(self):
        return os.path.join(self.sub_run_dir, 'hv_monitor.csv')

    def det_plane_z(self, default=232.0):
        """Detector plane z [mm] for M3 projection. Priority: the measured
        PLANE_Z from the 03 z-scan (if registered), else the run_config
        det_center z of DET_NAME, else `default`. This keeps every stage on
        the fitted plane when it disagrees with the nominal mounting height
        (P2_1 at p1_z ~232, P2_2 nominal 702 / fitted 712)."""
        if self.PLANE_Z is not None:
            return float(self.PLANE_Z)
        try:
            with open(self.run_config_path) as f:
                for d in json.load(f)['detectors']:
                    if d.get('name') == self.DET_NAME:
                        return float(d['det_center_coords']['z'])
        except (OSError, KeyError, ValueError, TypeError):
            pass
        return default

    # -- output helper ------------------------------------------------------ #
    def out_dir(self, *parts):
        d = os.path.join(self.OUT_BASE, *parts)
        os.makedirs(d, exist_ok=True)
        return d

    @property
    def dead_suffix(self):
        """Filename tag noting the dropped dead connectors, e.g.
        '_without_connector_10' (empty string if none are dropped). Appended to
        product filenames so it is obvious a plot excludes a dead connector."""
        if not self.DEAD_CONNECTORS:
            return ''
        conns = '_'.join(str(c) for c in sorted(self.DEAD_CONNECTORS))
        word = 'connector' if len(self.DEAD_CONNECTORS) == 1 else 'connectors'
        return f'_without_{word}_{conns}'

    SPARK_SUFFIX = '_spark_vetoed'

    def product_suffix(self, veto_sparks=False):
        """Filename tag for a stage's products: always notes dropped dead
        connectors, and appends the spark-veto tag when the veto is applied.
        Keeps spark-vetoed products from overwriting the un-vetoed ones."""
        return self.dead_suffix + (self.SPARK_SUFFIX if veto_sparks else '')

    def __repr__(self):
        return (f'<P2 cfg {self.KEY}: {self.DET_TAG} {self.RUN}/{self.SUB_RUN} '
                f'-> {self.OUT_BASE}>')


# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #
RUNS = {
    # First P2 detector-1 efficiency run (fetched 2026-07-01).
    'det1_long': _Config(
        'det1_long',
        run='p2_det1_long_run_6-30-26',
        sub_run='efficiency_long_run_p2_det1',
        det_name='P2_1',
        dead_connectors=(10,)),   # connector 10 disconnected on P2 det1

    # Second P2 detector-1 long efficiency run (fetched 2026-07-06), taken at a
    # single working point (mesh 440 V, drift 600 V) -- NOT an HV scan. Same
    # detector as det1_long; connectors 1 and 10 were both disconnected.
    'det1_long2': _Config(
        'det1_long2',
        run='p2_det1_long_run_7-4-26',
        sub_run='mesh_440V_drift_600V',
        det_name='P2_1',
        dead_connectors=(1, 10)),

    # Third P2 detector-1 long efficiency run (fetched 2026-07-08), same single
    # working point as det1_long2 (mesh 440 V, drift 600 V). Connectors 1 and
    # 10 still disconnected.
    'det1_long3': _Config(
        'det1_long3',
        run='p2_det1_long_run_7-7-26',
        sub_run='mesh_440V_drift_600V',
        det_name='P2_1',
        dead_connectors=(1, 10)),

    # Combined P2_1 + P2_2 long run 7-9-26 (10 h at mesh 430 V / drift 600 V,
    # drift gap 170 V), first run with detector 2 installed at p2_z. Data can
    # be fetched mid-run; products grow as the on-the-fly processor catches up.
    # One registry entry per detector, same run/sub_run:
    'det1_long4': _Config(
        'det1_long4',
        run='p2_det1_det2_long_run_mesh_scan_7-9-26',
        sub_run='long_run_mesh_430V_drift_600V',
        det_name='P2_1',
        det_tag='det1',
        spark_channel='1:0',       # P2_1 mesh = CAEN card 1 ch 0
        dead_connectors=(1, 10)),  # same disconnections as det1_long2/3
    'det2_long1': _Config(
        'det2_long1',
        run='p2_det1_det2_long_run_mesh_scan_7-9-26',
        sub_run='long_run_mesh_430V_drift_600V',
        det_name='P2_2',
        det_tag='det2',
        spark_channel='1:2',       # P2_2 mesh = CAEN card 1 ch 2
        # Connector 1 is physically disconnected from the detector; connectors
        # 8-10 sit on FEU 8, which is excluded from the run (crashes the DAQ).
        # Live readout = connectors 2-7 (deduced from track-hit correlation,
        # hit-level fit: scale 0.92, median residual 11 mm).
        dead_connectors=(1, 8, 9, 10),
        match_r=40.0,      # eff-vs-R plateau (12_validation knob test)
        plane_z=712.0),    # 03 z-scan fit; nominal p2_z would be 702

    # Mesh-HV scan tacked onto the 7-9-26 combined run: scan_mesh_<NNN>V
    # sub_runs, mesh 430->395 V in 5 V steps with the drift stepped in tandem
    # (drift gap fixed at 170 V). The long_run_mesh_430V sub_run is picked up
    # by find_subruns too and acts as a high-stats 430 V point.
    'det1_hvscan2': _Config(
        'det1_hvscan2',
        run='p2_det1_det2_long_run_mesh_scan_7-9-26',
        sub_run='hv_scan',
        det_name='P2_1',
        det_tag='det1',
        spark_channel='1:0',
        dead_connectors=(1, 10)),
    'det2_hvscan1': _Config(
        'det2_hvscan1',
        run='p2_det1_det2_long_run_mesh_scan_7-9-26',
        sub_run='hv_scan',
        det_name='P2_2',
        det_tag='det2',
        spark_channel='1:2',
        dead_connectors=(1, 8, 9, 10),   # see det2_long1
        match_r=40.0,
        plane_z=712.0),

    # Mesh-HV scan (mesh 345->420 V in 5 V steps, drift = mesh + 180 V), one
    # sub_run per HV point (mesh_<NNN>V_drift_<MMM>V).
    'det1_hvscan': _Config(
        'det1_hvscan',
        run='p2_det1_mesh_hv_scan_7-2-26',
        sub_run='hv_scan',
        det_name='P2_1',
        # connectors 1 and 10 both disconnected during the HV scan
        dead_connectors=(1, 10)),
}


def get_config(key=None) -> _Config:
    key = key or DEFAULT_RUN
    if key not in RUNS:
        raise KeyError(f'Unknown run key {key!r}. Known: {sorted(RUNS)}')
    return RUNS[key]


def config_from_argv() -> _Config:
    """First CLI arg (if present and a known key) selects the run."""
    key = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith('-') \
        else DEFAULT_RUN
    return get_config(key)


if __name__ == '__main__':
    for k in RUNS:
        c = get_config(k)
        print(c)
        print('  run_config :', c.run_config_path,
              '(exists)' if os.path.isfile(c.run_config_path) else '(MISSING)')
        print('  hits_dir   :', c.combined_hits_dir,
              '(exists)' if os.path.isdir(c.combined_hits_dir) else '(MISSING)')
        print('  map_csv    :', c.MAP_CSV_PATH,
              '(exists)' if os.path.isfile(c.MAP_CSV_PATH) else '(MISSING)')
