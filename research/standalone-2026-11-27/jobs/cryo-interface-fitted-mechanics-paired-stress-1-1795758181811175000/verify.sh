#!/bin/sh
set -eu
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONOPTIMIZE=0 PYTHONDONTWRITEBYTECODE=1
python=${PYTHON:-python3}
if ! "$python" -B -c 'import numpy, scipy' 2>/dev/null; then
    if [ ! -x .venv/bin/python ]; then
        "$python" -m venv .venv
    fi
    python=.venv/bin/python
    if ! "$python" -B -c 'import numpy, scipy' 2>/dev/null; then
        mkdir -p .scratch/tmp
        TMPDIR="$PWD/.scratch/tmp" PIP_NO_CACHE_DIR=1 "$python" -m pip \
            --disable-pip-version-check install -r requirements.txt
    fi
fi
exec "$python" -B verify.py
