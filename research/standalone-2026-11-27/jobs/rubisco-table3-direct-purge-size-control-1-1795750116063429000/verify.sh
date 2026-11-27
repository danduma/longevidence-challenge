#!/bin/sh
set -eu
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
if [ -z "${PYTHON:-}" ] && [ -x .venv/bin/python ]; then
    PYTHON=.venv/bin/python
fi
exec "${PYTHON:-python3}" -B verify.py
