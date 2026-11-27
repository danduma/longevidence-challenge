#!/usr/bin/env bash
# Refit from cached feature matrices using an installed or locally provisioned environment.
set -eu
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
if [ -d runtime/verify-deps ] && command -v python3.14 >/dev/null; then
    export PYTHONPATH="$PWD/runtime/verify-deps"
    exec python3.14 -B verify.py "$@"
fi
exec python3 -B verify.py "$@"
