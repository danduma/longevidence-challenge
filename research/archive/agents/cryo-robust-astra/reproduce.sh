#!/bin/sh
set -eu
cd "$(dirname "$0")"
PYTHON=${PYTHON:-../../.venv/bin/python}
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="$PWD/.mplcache"
test "$(git -C upstream rev-parse HEAD)" = 1d498ee5c583e0112a762e60892d576ec2e27f57
shasum -a 256 -c outputs/preregistration.sha256
clang++ -O3 -std=c++17 -shared -fPIC thermal.cpp -o thermal.so
"$PYTHON" -u investigate.py >outputs/run.log
"$PYTHON" verify.py >outputs/independent.log
"$PYTHON" analytic_limits.py >outputs/analytic.log
"$PYTHON" refinement.py >outputs/refinement.log
"$PYTHON" boundary.py >outputs/boundary.log
"$PYTHON" split_fine.py >outputs/split_fine.log
"$PYTHON" local_tolerance.py >outputs/local_tolerance.log
"$PYTHON" assay_bound.py >outputs/assay_bound.json
"$PYTHON" figure.py
"$PYTHON" audit.py >outputs/final_checks.log
