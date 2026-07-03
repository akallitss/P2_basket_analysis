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
                 spark_guard_before=2.0, spark_guard_after=10.0):
        self.KEY = key
        self.DATA_ROOT = data_root
        self.RUN = run
        self.SUB_RUN = sub_run
        self.DET_NAME = det_name
        self.DET_TYPE = det_type            # detector-under-test det_type in run_config
        self.REF_DET_TYPE = ref_det_type    # reference tracker det_type (m3 telescope)
        self.MAP_CSV_PATH = map_csv
        self.DET_TAG = det_tag_from_run(run)
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

    # Mesh-HV scan (mesh 345->420 V in 5 V steps, drift = mesh + 180 V), one
    # sub_run per HV point (mesh_<NNN>V_drift_<MMM>V).
    'det1_hvscan': _Config(
        'det1_hvscan',
        run='p2_det1_mesh_hv_scan_7-2-26',
        sub_run='hv_scan',
        det_name='P2_1',
        dead_connectors=(10,)),
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
