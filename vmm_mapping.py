#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 19/11/2025 19:17
Created in PyCharm
Created as vmm_mapping.py

@author: akallits
"""


# VMM mapping details combined
vmm_mapping = {
    'trigger': {
    'vmm_ids'      : [0, 1],
    'connector_ids': [10, 11],
    'orientation'  : None,
    'channels'     : {
        0: [40, 48],
        1: [0, 20, 40, 60]
    }
    },
    'p2_large_1': {
        'name': 'P2 Large Detector',
        'vmm_ids': [12, 13, 14, 15],
        # K59V adapter: each FPC connector → 2 MEC8 connectors → 2 VMMs
        # (fpc_connector_number, mec8_connector) → vmm_id
        # physical FPC connectors 4 and 5 were cabled in run 67;
        # FPC 4 → VMMs 12,13 and FPC 5 → VMMs 14,15 (A-style wiring).
        # Connector orientation 'flipped_back': end-to-end pin reverse AND
        # back-side insertion (adjacent channel pairs swapped, ch ^ 1).
        # Chosen by visual inspection of the 16-panel orientation scan
        # (A-flipped_back shows the most concentrated beam profile). NB: run 67
        # is a wide muon beam, so the spread metric only weakly distinguishes
        # orientations (all 16 within 107–132 mm) — revisit with a focused beam.
        'mec8_to_vmm': {
            (4, 1): 12,
            (4, 0): 13,
            (5, 1): 14,
            (5, 0): 15,
        },
        'fpc_connectors': [4, 5],
        # Channel→pad ordering. 'strategy' (strip-based) takes precedence over
        # 'orientation' (pin-based). 'reverse' is the within-half order
        # validated against M3 tracks in the cosmic-bench DREAM readout
        # (cosmic_bench_analysis/p2_mapping.py, session 2026-07-01): it reads
        # the geometry from the CSV's own `strip` column, so it is not subject
        # to the MEC8 2-pin numbering gap that no ch=f(mec8_pin) orientation
        # can reproduce. 'orientation' is retained only for the pin-based
        # 16-panel comparison diagnostic (weak/inconclusive on the wide run-67
        # muon beam).
        'strategy': 'reverse',
        'orientation': 'flipped_back',
    },
    'p2_small_1': {
        'name': 'P2 Small Detector 1',
        'vmm_ids': [10, 11],
        'connector_ids': [0, 1],
        'orientation': ['normal', 'normal', 'normal', 'normal'],
        # Connected channels derived from p2_small_detector_map.csv:
        #   connector 0 (VMM 10): pads 0-31  → VMM ch 0-31
        #   connector 1 (VMM 11): pads 32-81 → VMM ch 14-63
        'channels': {
            10: list(range(32)),
            11: list(range(14, 64)),
        },
    },
    'p2_small_3': {
        'name': 'P2 Small Detector 3',
        'vmm_ids': [8, 9],
        'connector_ids': [0, 1],
        'orientation': ['normal', 'normal', 'normal', 'normal'],
        # Same detector type as p2_small_1, same channel mapping:
        #   connector 0 (VMM 8):  VMM ch 0-31
        #   connector 1 (VMM 9):  VMM ch 14-63
        'channels': {
            8: list(range(32)),
            9: list(range(14, 64)),
        },
    },
}
