#!/usr/bin/env bash
# code_version.sh -- the provenance stamp that travels with a condor code tarball.
#
#   code_version.sh <repo> <patch_file> <path>...
#
# Prints a CODE_VERSION section for <repo> to stdout: commit, nearest tag,
# packing time, and every uncommitted or untracked file under the packed
# <path>s. When tracked files there differ from the commit, their diff is
# written to <patch_file> (removed otherwise), so a product can always be
# traced to the exact code even when it was packed from a dirty tree.
set -euo pipefail
REPO=$1; PATCH=$(realpath -m "$2"); shift 2
cd "$REPO"
dirty=$(git status --porcelain -- "$@")
echo "repo:      $(basename "$(git rev-parse --show-toplevel)")"
echo "commit:    $(git rev-parse HEAD)"
echo "describe:  $(git describe --tags --always)"
echo "committed: $(git log -1 --format=%ci)"
echo "packed:    $(date -u +%Y-%m-%dT%H:%M:%SZ) by $(whoami)@$(hostname -s)"
if [ $# -le 5 ]; then echo "paths:     $*"; else echo "paths:     $# files"; fi
if [ -n "$dirty" ]; then
    echo "uncommitted (M = modified, in $(basename "$PATCH"); ?? = untracked):"
    echo "$dirty" | sed 's/^/  /'
else
    echo "uncommitted: none"
fi
if git diff --quiet HEAD -- "$@"; then
    rm -f "$PATCH"
else
    git diff HEAD -- "$@" > "$PATCH"
fi
