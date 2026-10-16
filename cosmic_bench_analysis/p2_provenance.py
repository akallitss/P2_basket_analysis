"""
Provenance record for a stage's products: which code and which inputs made
them. One `inputs<sfx>.txt` per stage directory, written next to the figures,
so a product can be traced without hashing files after the fact.

    import p2_provenance as prov
    prov.write(out_dir, sfx, stage=__file__, inputs={'toa table': path, ...},
               settings={'amp_min': 300, ...})
"""
import datetime as dt
import os
import socket
import subprocess
import sys


def _git(repo, *args):
    try:
        # rstrip only: `status --porcelain` lines start with a status column
        return subprocess.run(['git', '-C', repo, *args], capture_output=True,
                              text=True, check=True).stdout.rstrip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ''


def code_version(path):
    """Commit, nearest tag and uncommitted files of the repo holding `path`."""
    repo = os.path.dirname(os.path.abspath(path))
    top = _git(repo, 'rev-parse', '--show-toplevel')
    if not top:
        return ['code:      not in a git repository']
    dirty = _git(top, 'status', '--porcelain', '--', 'cosmic_bench_analysis',
                 'sps_beam_analysis')
    lines = [f'repo:      {os.path.basename(top)}',
             f'commit:    {_git(top, "rev-parse", "HEAD")}',
             f'describe:  {_git(top, "describe", "--tags", "--always")}',
             f'committed: {_git(top, "log", "-1", "--format=%ci")}']
    if dirty:
        lines.append('uncommitted (M = modified, ?? = untracked):')
        lines += ['  ' + l for l in dirty.splitlines()]
    else:
        lines.append('uncommitted: none')
    return lines


def _describe_input(p):
    if isinstance(p, str) and os.path.exists(p):
        t = dt.datetime.fromtimestamp(os.path.getmtime(p)).strftime('%Y-%m-%d %H:%M')
        return f'{p}  (modified {t})'
    return str(p)


def write(out_dir, sfx, stage, inputs, settings=None, notes=None):
    """Write <out_dir>/inputs<sfx>.txt and return its path."""
    now = dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    lines = [f'stage:     {os.path.basename(stage)}',
             f'command:   {" ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])}',
             f'run at:    {now} on {socket.gethostname()}', '']
    lines += code_version(stage) + ['', 'inputs:']
    lines += [f'  {k}: {_describe_input(v)}' for k, v in inputs.items()]
    if settings:
        lines += ['', 'settings:'] + [f'  {k}: {v}' for k, v in settings.items()]
    if notes:
        lines += ['', 'notes:'] + [f'  {n}' for n in notes]
    path = os.path.join(out_dir, f'inputs{sfx}.txt')
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return path
