#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""check_radii.py -- guard: every analysis radius comes from the run config.

The match radius R and the active-area radius are physics cuts, tuned once
(Analysis/DET1_EFFICIENCY.md) and then required to be identical at every point
in the chain. They have escaped twice already:

  * `--active-r` was a hard-coded argparse default of 30 mm in 06/11/12/16.
    Fixing only 06 and 12 left the mesh and drift scans computing efficiency
    against a denominator ~3.6x larger than the long runs -- a silent 10-point
    disagreement between two stages of the same pipeline.
  * 17_lifetime_autopsy.py carried its own module-level MATCH_R / ACTIVE_R.

Both were invisible because nothing crashed; the numbers were just wrong.
This script fails loudly instead.

    python3 check_radii.py          # exit 0 = all clean
"""

import ast
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# argparse options that carry a physics radius. Each MUST default to None so
# the stage falls back to the per-run config value.
RADIUS_OPTS = {'--r', '--active-r', '--match-r'}

# Options that are display/smoothing knobs, not cuts -- a literal default is
# fine for these.
EXEMPT_OPTS = {'--kernel', '--grid', '--min', '--chi2-cut', '--fit-fiducial',
               '--time-window-h'}

# Module-level names that would shadow the config.
SHADOW_NAMES = {'MATCH_R', 'ACTIVE_R'}


def stage_files():
    pats = ('[0-9][0-9]_*.py', 'build_*.py', 'p2_*.py', 'm3_*.py', 'zz_*.py')
    out = []
    for p in pats:
        out += glob.glob(os.path.join(HERE, p))
    return sorted(set(out))


def check(path):
    """Return a list of problem strings for one file."""
    problems = []
    src = open(path).read()
    name = os.path.basename(path)
    try:
        tree = ast.parse(src)
    except SyntaxError as e:
        return [f'{name}: does not parse ({e})']

    # 1. argparse radius options must default to None
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and
                isinstance(node.func, ast.Attribute) and
                node.func.attr == 'add_argument'):
            continue
        if not node.args or not isinstance(node.args[0], ast.Constant):
            continue
        opt = node.args[0].value
        if opt not in RADIUS_OPTS or opt in EXEMPT_OPTS:
            continue
        default = next((k.value for k in node.keywords if k.arg == 'default'), None)
        if default is None:                      # no default= at all -> None
            continue
        if isinstance(default, ast.Constant) and default.value is None:
            continue
        lit = getattr(default, 'value', ast.dump(default))
        problems.append(
            f"{name}: '{opt}' defaults to {lit!r}, not None -- it will "
            f"override the run config silently")

    # 2. the fallback must actually be wired to the config
    for opt, cfg_attr in (('--r', 'MATCH_R'), ('--active-r', 'ACTIVE_R')):
        if f"'{opt}'" not in src:
            continue
        if not re.search(r'cfg\.' + cfg_attr, src):
            problems.append(
                f'{name}: accepts {opt} but never reads cfg.{cfg_attr} -- '
                f'no fallback to the run config')

    # 3. no module-level constant shadowing the config value with a literal
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for t in node.targets:
            if isinstance(t, ast.Name) and t.id in SHADOW_NAMES:
                if isinstance(node.value, ast.Constant):
                    problems.append(
                        f'{name}: module-level {t.id} = {node.value.value!r} '
                        f'shadows the run config with a literal')
    return problems


def main():
    sys.path.insert(0, HERE)
    all_problems = []
    checked = 0
    for f in stage_files():
        if os.path.basename(f) in ('check_radii.py', 'p2_qa_config.py'):
            continue
        checked += 1
        all_problems += check(f)

    print(f'checked {checked} stage files')
    if all_problems:
        print(f'\n{len(all_problems)} PROBLEM(S):')
        for p in all_problems:
            print('  ✗ ' + p)
        return 1
    print('  ✓ every radius option defaults to None and falls back to the run config')
    print('  ✓ no module-level MATCH_R / ACTIVE_R literals')
    return 0


if __name__ == '__main__':
    sys.exit(main())
