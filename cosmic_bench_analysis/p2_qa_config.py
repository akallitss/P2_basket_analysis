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
# Input run data (combined_hits_root / m3_tracking_root / hv_monitor.csv ...)
# lives on the LaCie external drive; the internal disk was chronically full.
DATA_ROOT = ('/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/'
             'data/Cosmic_Bench')
# Analysis products (plots, csv, pdf) stay on the internal disk as before.
ANALYSIS_ROOT = ('/local/home/ak271430/Documents/PostDocSaclay/'
                 'data/Cosmic_Bench/Analysis')

# Gerber-derived pad mapping (channel_id -> pad geometry).
MAP_CSV_PATH = os.path.join(
    REPO_ROOT, 'Detector_Mapping', 'P2_BASKET', 'P2_BASKET_mapping.csv')

# Insulation-mask Gerber (same KiCad board frame as the pad map): the exact
# positions of the mesh-support pillars — ~11.7k small (~0.5 mm) on a ~4 mm
# grid plus 5 big (~3.3 mm) ones. Overlaid on hitmaps/efficiency maps so
# pillar-shadow dead spots can be told apart from real gain defects.
MASK_GBR_PATH = os.path.join(
    os.path.dirname(REPO_ROOT), 'Detector_Drawings', 'Version_Apr26',
    'Insulation_masks', 'V2', 'P2_BASKET-Mask_M2_V2.gbr')

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

# Filename tag for products made with a non-default chi2 cut (cut-sensitivity
# checks). Empty for the default M3_CHI2_CUT so the standard products keep
# their names; '_chi2' otherwise so variant products never overwrite them.
CHI2_SUFFIX = '_chi2'


def chi2_tag(chi2_cut) -> str:
    return '' if float(chi2_cut) == M3_CHI2_CUT else CHI2_SUFFIX


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
                 map_csv=MAP_CSV_PATH, mask_gbr=MASK_GBR_PATH,
                 dead_connectors=(),
                 spark_channel='1:0', spark_imon_thr=2.0,
                 spark_guard_before=2.0, spark_guard_after=10.0,
                 burst_npads=20, det_tag=None, match_r=20.0, plane_z=None,
                 t_max_h=None, min_amp=0.0, out_tag=None, noisy_pads=(),
                 hot_pad_ratio=5.0):
        self.KEY = key
        self.DATA_ROOT = data_root
        self.RUN = run
        self.SUB_RUN = sub_run
        self.DET_NAME = det_name
        self.DET_TYPE = det_type            # detector-under-test det_type in run_config
        self.REF_DET_TYPE = ref_det_type    # reference tracker det_type (m3 telescope)
        self.MAP_CSV_PATH = map_csv
        # Insulation-mask Gerber with the mesh-support pillar positions (pad
        # frame). None/missing file -> pillar overlays are simply skipped.
        self.MASK_GBR_PATH = mask_gbr
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
        # --- per-run data-quality cuts (p2_io applies them at read time) ---- #
        # Analysis time window [h since first trigger]. Events after t_max_h
        # are dropped from BOTH the hits and the M3 rays (a detector that
        # tripped mid-run must not dilute occupancy or efficiency). None = all.
        self.T_MAX_H = t_max_h
        # Minimum hit amplitude [ADC]. Kills a stale-pedestal noise floor
        # (hits at ~threshold on every pad) before centroids / efficiency /
        # burst counting. 0 = keep everything.
        self.MIN_AMP = float(min_amp)
        # channel_ids of pads with a KNOWN noisy/oscillating channel whose
        # hits pass MIN_AMP (e.g. det3 pad 510 = FEU6 ch321 firing at
        # 120-150 ADC, 55% of the 7-16 initial run's signal band). Excluded
        # from centroids/efficiency; QA stages (02/05) still show them so
        # the pathology stays visible.
        self.NOISY_PADS = tuple(noisy_pads)
        # Automatic hot-pad cut (p2_io.auto_hot_pads): a pad carrying more
        # than hot_pad_ratio x the median occupancy of the fired pads (and
        # >= 30 hits) is constantly firing -- not spark-like, so the HV/burst
        # vetoes never catch it, but it kills hitmap colour scales and biases
        # centroids (e.g. det1 7-19 pad 1089 at 12x median). Flagged per
        # sub_run at load time and dropped like NOISY_PADS. 0 disables.
        self.HOT_PAD_RATIO = float(hot_pad_ratio)

        # Analysis/ tree (ANALYSIS_ROOT, internal disk), keyed by detector ->
        # run -> sub_run. out_tag adds a suffix directory so a windowed variant
        # of the same sub_run (e.g. the pre-discharge hours only) never
        # overwrites the full-run products.
        self.ANALYSIS_ROOT = ANALYSIS_ROOT
        self.OUT_BASE = os.path.join(ANALYSIS_ROOT, self.DET_TAG, run,
                                     sub_run + (f'_{out_tag}' if out_tag else ''))

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
    def decoded_root_dir(self):
        # Raw per-FEU decoded waveforms (nt tree: eventId/timestamp/ftst +
        # sample/channel/amplitude vectors) — needed by the waveform-timing
        # stage (13). Bulky: fetch only the FEUs/chunks needed from rays.
        return os.path.join(self.sub_run_dir, 'decoded_root')

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

    # det4 (P2_4, bulked 8-7-26 Alex+Enzo, first with the cleanest evac line)
    # at p1_z, FEUs 3+4 (connectors 2-9), mesh 1:0 / drift 1:1. 12 h long run
    # at mesh 480 V / drift 700 V, started 2026-07-15 22:42; the same run also
    # carries drift_scan_det4_* (drift 700->500, 50 V steps) and
    # mesh_scan_det4_* (mesh 480->405, 5 V steps, drift in tandem) sub_runs.
    'det4_long1': _Config(
        'det4_long1',
        run='p2_det4_long_run_drift_mesh_scan_7-15-26',
        sub_run='long_run_det4_480_700',
        det_name='P2_4',
        det_tag='det4',
        spark_channel='1:0',       # P2_4 mesh = CAEN card 1 ch 0
        dead_connectors=(1, 10),   # physically disconnected
        # Mesh 1:0 imon climbed from 2.2 h and CAEN tripped it at 5.72 h
        # (04:25 7-16); everything after is mesh-off noise. The run also has
        # a ~680 hits/event noise floor at ~18 ADC on every pad (stale
        # pedestals: do_pedestal_threshold_run=False), so only hits above
        # 100 ADC carry signal information.
        t_max_h=5.7,
        min_amp=100.0),

    # Stable window ONLY: mesh imon is quiet (~0.02 uA) until 2.2 h, then the
    # detector discharges continuously until the 5.72 h trip. This variant
    # keeps just the healthy hours (products land in <sub_run>_stable0-2.2h/).
    # NOTE: the m3 ray files for chunks 000-001 are empty on the DAQ, so there
    # are NO reference tracks before 2.56 h -> no alignment/efficiency here
    # unless the rays are reprocessed from the raw .fdf on rays_daplxa.
    'det4_long1_stable': _Config(
        'det4_long1_stable',
        run='p2_det4_long_run_drift_mesh_scan_7-15-26',
        sub_run='long_run_det4_480_700',
        det_name='P2_4',
        det_tag='det4',
        spark_channel='1:0',
        dead_connectors=(1, 10),
        t_max_h=2.2,
        min_amp=100.0,
        out_tag='stable0-2.2h'),

    # det3+det4 simultaneous drift scan (7-16-26): one sub_run per point,
    # named drift_scan_det4_<mesh>_<drift>_det3_<mesh>_<drift>, 30 min each.
    # det4 mesh 450 V / drift 450..900; det3 mesh 420 V / drift 420..900.
    # The DAQ again ran with do_pedestal_threshold_run=False, so the ~18 ADC
    # noise floor covers every pad -> min_amp cut required (as det4_long1).
    # sub_run here is only the products directory (stage 16 loops sub_runs).
    # Initial (working-point) run of the 7-16 evening det3-mesh-scan session:
    # det3 at its drift-scan plateau (mesh 420 / drift 820), det4 at mesh 430 /
    # drift 830 (mesh raised vs the flat drift scan at 450... note 430<450 —
    # probing below the discharge threshold seen at 480 V). Same stale-pedestal
    # noise floor (do_pedestal_threshold_run=False) -> min_amp=100.
    'det3_initial1': _Config(
        'det3_initial1',
        run='p2_det3_mesh_scan_det4_initial_7-16-26',
        sub_run='initial_run_det3_420_820_det4_430_830',
        det_name='P2_3',
        det_tag='det3',
        spark_channel='1:2',
        dead_connectors=(1, 8, 9, 10),
        match_r=40.0,
        min_amp=0.0,      # reprocessed w/ real pedestals (thr ~28 ADC)
        noisy_pads=(510,)),   # FEU6 ch321: 120-150 ADC oscillation, 55% of hits
    # det3 mesh scan of the same 7-16 evening run (mesh 345..420 in 5 V steps,
    # drift in tandem = mesh + 400). sub_run is only the products dir; stage 16
    # --scan mesh loops the mesh_scan_det3_* sub_runs. Same noisy pad 510.
    'det3_meshscan1': _Config(
        'det3_meshscan1',
        run='p2_det3_mesh_scan_det4_initial_7-16-26',
        sub_run='mesh_scan',
        det_name='P2_3',
        det_tag='det3',
        spark_channel='1:2',
        dead_connectors=(1, 8, 9, 10),
        match_r=40.0,
        min_amp=0.0,      # reprocessed w/ real pedestals (thr ~28 ADC)
        noisy_pads=(510,)),

    # 10h run at det3's working point (mesh 420 / drift 820), 7-17 05:56-15:57,
    # right after the descending mesh scan. Reprocessed with the 19H40 pedestal
    # (FEU-regex fix) -> min_amp 0. Same noisy pad 510 masked from centroids.
    'det3_final1': _Config(
        'det3_final1',
        run='p2_det3_mesh_scan_det4_initial_7-16-26',
        sub_run='final_run_det3_420_820',
        det_name='P2_3',
        det_tag='det3',
        spark_channel='1:2',
        dead_connectors=(1, 8, 9, 10),
        match_r=40.0,
        min_amp=0.0,
        noisy_pads=(510,)),

    # det1 (P2_1) repaired and re-installed at p2_z (det3's old slot) for the
    # 7-19-26 run: FEUs 6+7 with connectors 2-9 cabled (1 and 10 unbonded),
    # mesh on CAEN 1:2 / drift 1:3. 2 h initial run at mesh 430 / drift 740
    # (gap 310 V), then a descending mesh scan 430->355 V in 5 V steps with
    # the drift in tandem (drift = mesh + 310), 30 min per point. A fresh
    # pedestal run was taken at 00:21, 2 min before the run -> min_amp 0.
    'det1_initial1': _Config(
        'det1_initial1',
        run='p2_det1_long_run_mesh_scan_7-19-26',
        sub_run='initial_run_det1_430_740',
        det_name='P2_1',
        det_tag='det1',
        spark_channel='1:2',
        dead_connectors=(1, 10),
        match_r=40.0,               # far plane, like det2/det3 (z ~702)
        min_amp=0.0),
    # Mesh scan of the same run; sub_run is only the products dir, stage 16
    # --scan mesh (and stage 11) loop the mesh_scan_det1_* sub_runs.
    'det1_meshscan1': _Config(
        'det1_meshscan1',
        run='p2_det1_long_run_mesh_scan_7-19-26',
        sub_run='mesh_scan',
        det_name='P2_1',
        det_tag='det1',
        spark_channel='1:2',
        dead_connectors=(1, 10),
        match_r=40.0,
        min_amp=0.0),

    'det4_initial1': _Config(
        'det4_initial1',
        run='p2_det3_mesh_scan_det4_initial_7-16-26',
        sub_run='initial_run_det3_420_820_det4_430_830',
        det_name='P2_4',
        det_tag='det4',
        spark_channel='1:0',
        dead_connectors=(1, 10),
        min_amp=0.0),     # reprocessed w/ real pedestals (thr ~28 ADC)

    # det1 DRIFT scan (7-19-26 afternoon), taken right after the mesh scan
    # settled on 415 V as the working point: mesh FIXED at 415 V, drift stepped
    # 415->965 V in 50 V steps (drift gap 0->550 V), 30 min per point. Probes
    # the drift-field response at fixed gain. sub_runs named
    # drift_scan_det1_415_<drift>; sub_run here is only the products dir (stage
    # 16 --scan drift loops them). Same 00:21 pedestal as det1_meshscan1 ->
    # min_amp 0. Fetched mid-run (processor lags the DAQ ~5-10 min/point).
    'det1_driftscan2': _Config(
        'det1_driftscan2',
        run='p2_det1_drift_scan_7-19-26',
        sub_run='drift_scan',
        det_name='P2_1',
        det_tag='det1',
        spark_channel='1:2',
        dead_connectors=(1, 10),
        match_r=40.0,
        min_amp=0.0),

    # det1 24 h efficiency long run (7-19-26 19:24) at the working point chosen
    # from today's mesh + drift scans: mesh 415 V / drift 615 V (drift gap
    # 200 V, E ~ 667 V/cm) -- the timing-optimal AND efficiency-plateau point.
    # Single sub_run long_run_det1_415_615. Same 00:21 pedestal -> min_amp 0.
    # Fetched mid-run (products grow as the on-the-fly processor catches up).
    'det1_long5': _Config(
        'det1_long5',
        run='p2_det1_long_run_efficiency_7-19-26',
        sub_run='long_run_det1_415_615',
        det_name='P2_1',
        det_tag='det1',
        spark_channel='1:2',
        dead_connectors=(1, 10),
        match_r=40.0,
        min_amp=0.0),

    'det3_driftscan1': _Config(
        'det3_driftscan1',
        run='p2_det3_det4_drift_scan_7-16-26',
        sub_run='drift_scan',
        det_name='P2_3',
        det_tag='det3',
        spark_channel='1:2',        # P2_3 mesh = CAEN card 1 ch 2
        # c_2..c_7 wired (FEU6: conns 2-5 both halves, FEU7: conns 6-7)
        dead_connectors=(1, 8, 9, 10),
        match_r=40.0,               # far plane, like det2 (z ~702)
        min_amp=0.0),     # reprocessed w/ real pedestals (thr ~28 ADC)
    'det4_driftscan1': _Config(
        'det4_driftscan1',
        run='p2_det3_det4_drift_scan_7-16-26',
        sub_run='drift_scan',
        det_name='P2_4',
        det_tag='det4',
        spark_channel='1:0',        # P2_4 mesh = CAEN card 1 ch 0
        dead_connectors=(1, 10),
        min_amp=0.0),     # reprocessed w/ real pedestals (thr ~28 ADC)

    # Fe55 mesh-HV scan of the SPS-telescope detectors on the banco bench
    # (7-18-26, banco_daplxa:/local/home/banco/P2_data/Fe55/runs/run_1):
    # det2 = P2_OUT (bulked 25-6-26, misaligned wall) and det3 = P2_MID
    # (bulked 2-7-26). NO M3 telescope here: the DREAM DAQ self-triggers on
    # the Fe55 source (TCM multiplicity), so none of the tracking/efficiency
    # stages apply -- the analysis is 18_fe55_spectra.py (per-event charge
    # spectrum + photopeak fit -> gain / resolution / rate vs mesh HV).
    # One sub_run per HV point, fe55_<NN>_mesh_out<V>_mid<V> (out = P2_OUT
    # mesh 420->365, mid = P2_MID mesh 510->455, both in 5 V steps with the
    # drifts stepped in tandem 700->645), 5 min per point. FEU 3 = P2_OUT,
    # FEU 4 = P2_MID, connectors 4-7 (both halves) each -- the channel table
    # is built from the run_config wiring, so no dead-connector drop needed.
    # Real pedestal/threshold run this time (do_pedestal_threshold_run=True).
    # burst_npads=0: in self-trigger mode with the source the per-event pad
    # multiplicity is naturally high (~12 hits/event on FEU 3 at 420 V), so
    # the cosmic burst veto would eat real events; HV spark veto still applies.
    'det2_fe55scan1': _Config(
        'det2_fe55scan1',
        run='p2_fe55_det2_det3_mesh_scan_7-18-26',
        sub_run='fe55_scan',
        det_name='P2_OUT',
        det_tag='det2',
        spark_channel='8:1',       # P2_OUT mesh = CAEN card 8 ch 1
        burst_npads=0),
    'det3_fe55scan1': _Config(
        'det3_fe55scan1',
        run='p2_fe55_det2_det3_mesh_scan_7-18-26',
        sub_run='fe55_scan',
        det_name='P2_MID',
        det_tag='det3',
        spark_channel='8:3',       # P2_MID mesh = CAEN card 8 ch 3
        burst_npads=0),

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
