#!/bin/sh
set -eu
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONOPTIMIZE=0
python=${PYTHON:-python3}
if ! "$python" -B -c 'import numpy, scipy' 2>/dev/null; then
    mkdir -p .scratch/tmp
    export TMPDIR="$PWD/.scratch/tmp" PIP_NO_CACHE_DIR=1
    if [ ! -x .venv/bin/python ]; then
        "$python" -m venv .venv
    fi
    python=.venv/bin/python
    if ! "$python" -B -c 'import numpy, scipy' 2>/dev/null; then
        "$python" -m pip \
            --disable-pip-version-check install -r requirements.txt
    fi
fi
"$python" -B verify.py
"$python" -B check_claim.py
