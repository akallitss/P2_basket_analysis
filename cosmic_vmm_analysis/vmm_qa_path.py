# -*- coding: utf-8 -*-
"""Puts DAQ_Control_VMM_Beam/vmm_qa on sys.path, so vmm_stations (the SPS channel->pad map) can be imported.

Looked up in order: $VMM_QA_DIR; the sibling checkout ../../DAQ_Control_VMM_Beam/vmm_qa; ./vmm_qa next to these
scripts (the copy used on lxplus).
"""
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
for _d in (os.environ.get('VMM_QA_DIR', ''),
           os.path.join(_here, '..', '..', 'DAQ_Control_VMM_Beam', 'vmm_qa'),
           os.path.join(_here, 'vmm_qa')):
    if _d and os.path.isfile(os.path.join(_d, 'vmm_stations.py')):
        sys.path.insert(0, os.path.abspath(_d))
        break
else:
    raise ImportError('vmm_stations.py not found: set VMM_QA_DIR to DAQ_Control_VMM_Beam/vmm_qa')
