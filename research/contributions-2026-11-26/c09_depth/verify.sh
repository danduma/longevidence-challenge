#!/bin/sh
set -eu
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
"${C09_PYTHON:-python3}" -B analyze.py
"${C09_PYTHON:-python3}" -B verify_expected.py
