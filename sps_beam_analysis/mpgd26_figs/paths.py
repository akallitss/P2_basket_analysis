"""Where the MPGD26 beam figures read from and write to.

The figures were first made in a scratch tree that no longer exists, which is
why nothing in the repo could re-make them.  This module is the fix: one place
naming the workspace, so every `fig_*.py` is re-runnable from a checkout.

    S    workspace root   $MPGD26_WORKSPACE, or the default below
    A    stage products   S/products/analysis      (usually a symlink to the
                          campaign analysis tree on the LaCie)
    URW  uRWELL-referenced efficiency products
    RD   report_data/     tidy CSVs written by aggregate.py
    OUT  figures          the report's own figs/ directory -- figures are
                          produced in their canonical location and only then
                          copied to conference/ by its gather script

`setup_workspace.sh` builds S; `aggregate.py` fills RD; each `fig_*.py`
writes into OUT.
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))

S = os.environ.get(
    'MPGD26_WORKSPACE',
    os.path.join(os.path.expanduser('~'), 'Documents', 'PostDocSaclay', 'data',
                 'SPS_Beam_Test', 'mpgd26_workspace'))

A = os.path.join(S, 'products', 'analysis')
URW = os.path.join(S, 'products', 'urw_local', 'urw_referenced_efficiency')
RD = os.path.join(S, 'report_data')
OUT = os.environ.get(
    'MPGD26_FIGS',
    os.path.join(REPO, 'reports', 'mpgd26_sps_beam_2026-08', 'figs'))

os.makedirs(RD, exist_ok=True)
os.makedirs(OUT, exist_ok=True)
