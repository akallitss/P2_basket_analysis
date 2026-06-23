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
        # FPC 4 → VMMs 14,15 and FPC 5 → VMMs 12,13 (C-style wiring);
        # flat ribbon cables plugged with pins reversed: ch = 66 − mec8_pin.
        # Determined from 8-combination orientation scan: C-flipped gives
        # tightest beam-spot cluster (top-5 spread 21 mm, centre ≈ 348,241 mm).
        'mec8_to_vmm': {
            (5, 1): 12,
            (5, 0): 13,
            (4, 1): 14,
            (4, 0): 15,
        },
        'fpc_connectors': [4, 5],
        'orientation': 'flipped',
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
