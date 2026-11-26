#!/bin/sh
set -eu
cd "$(dirname "$0")"
PYTHON=${PROTEASE_PYTHON:-/path/to/longevidence-checkout/research/millennium/.venv/bin/python}
"$PYTHON" -m unittest test_reconstruct.py
"$PYTHON" reconstruct.py > results/reconstruction.log 2>&1
"$PYTHON" verify_sources.py > results/source_validation.log 2>&1
"$PYTHON" exclusions.py > results/exclusions.log 2>&1
"$PYTHON" evaluate.py > results/evaluate.log 2>&1
"$PYTHON" impact.py > results/impact.log 2>&1
"$PYTHON" uncertainty.py > results/uncertainty.log 2>&1
"$PYTHON" manual_parent_sensitivity.py
"$PYTHON" provenance.py
"$PYTHON" final_checks.py
