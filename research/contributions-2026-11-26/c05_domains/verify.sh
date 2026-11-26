#!/bin/sh
set -eu
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
PYTHON=${PYTHON:-python3}
mkdir -p results/tmp
TMPDIR="$(pwd)/results/tmp"
export TMPDIR
"$PYTHON" check_inputs.py
"$PYTHON" build_solver.py
"$PYTHON" run.py
"$PYTHON" analyze.py
"$PYTHON" write_report.py
"$PYTHON" seal.py
