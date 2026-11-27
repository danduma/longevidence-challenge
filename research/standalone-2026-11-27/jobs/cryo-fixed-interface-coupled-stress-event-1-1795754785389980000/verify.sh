#!/bin/sh
set -eu
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONOPTIMIZE=0
if python3 -B -c 'import numpy, scipy' 2>/dev/null; then
    python=python3
else
    python=.venv/bin/python
    if [ ! -x "$python" ]; then
        python3 -m venv .venv
    fi
    if ! "$python" -B -c 'import numpy, scipy' 2>/dev/null; then
        mkdir -p .tmp
        TMPDIR="$PWD/.tmp" PIP_NO_CACHE_DIR=1 "$python" -m pip --disable-pip-version-check install -r requirements.txt
    fi
fi
"$python" -B verify.py
