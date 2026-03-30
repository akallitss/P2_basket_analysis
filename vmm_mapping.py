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

    },
    'p2_small_3': {
        'name': 'P2 Small Detector 3',
        'vmm_ids': [8, 9],
        'connector_ids': [0, 1],
        'orientation': ['normal', 'normal', 'normal', 'normal'],
    },
}
