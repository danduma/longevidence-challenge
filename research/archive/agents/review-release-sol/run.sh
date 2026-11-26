#!/bin/zsh
set -euo pipefail
cd "${0:A:h}"
PYTHON_AUDIT=/path/to/longevidence-checkout/research/millennium/.venv/bin/python
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
"$PYTHON_AUDIT" inspect_release.py > outputs_initial.txt
"$PYTHON_AUDIT" audit.py > audit_run.log
"$PYTHON_AUDIT" baseline.py > baseline_run.log
"$PYTHON_AUDIT" supplemental_checks.py > verification.log
# Optional network evidence refresh: "$PYTHON_AUDIT" fetch_evidence.py
# Pin verification uses preserved API evidence:
"$PYTHON_AUDIT" verify_release.py >> verification.log
