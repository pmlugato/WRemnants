#!/bin/bash
# Source setup.sh and exec the given command, preserving argument boundaries.
#
# Why this exists rather than `run_with_singularity.sh bash -c "..."`: the
# container's runscript re-parses what it is given, so a command passed as one
# string goes through a SECOND round of shell word-splitting inside the
# container. Anything quoted containing a space is then silently torn apart --
# `--unblind 'sig_JPsi.*'` survives, but a quoted `python -c` snippet or any
# argument with an embedded space does not.
#
# The CI wrappers (scripts/ci/setup_and_run{,_python}.sh) use unquoted `$@`,
# which splits for the same reason. `exec "$@"` here keeps argv exactly as the
# caller wrote it.
#
# Usage:
#   scripts/ci/run_with_singularity.sh scripts/rabbit/run_in_container.sh \
#       python scripts/histmakers/mz_dilepton.py --axes a b --postfix 'x y'

cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
. ./setup.sh >/dev/null 2>&1
exec "$@"
