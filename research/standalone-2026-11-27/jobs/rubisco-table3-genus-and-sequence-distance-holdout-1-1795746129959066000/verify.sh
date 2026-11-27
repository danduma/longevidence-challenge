#!/bin/sh
set -eu
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
if [ -x .venv/bin/python ]; then
    exec .venv/bin/python -B verify.py
fi
exec python3 -B verify.py
