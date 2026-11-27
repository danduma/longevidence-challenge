#!/bin/sh
set -eu
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONOPTIMIZE=0 PYTHONDONTWRITEBYTECODE=1
exec "${PYTHON:-python3}" -B analysis.py "$@"
