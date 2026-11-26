#!/bin/sh
set -eu
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
C01_PYTHON="${PYTHON:-python3}"
"$C01_PYTHON" run.py --jobs 2
"$C01_PYTHON" checks.py
"$C01_PYTHON" analyze.py
"$C01_PYTHON" verify.py
"$C01_PYTHON" write_report.py
