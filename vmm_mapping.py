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
        'connector_ids': [0, 1, 2, 3],
        'orientation': ['normal', 'normal', 'normal', 'normal'],
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
