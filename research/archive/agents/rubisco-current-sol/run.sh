#!/bin/sh
set -eu
cd "$(dirname "$0")"
PYTHON=/path/to/longevidence-checkout/research/millennium/.venv/bin/python
"$PYTHON" scripts/analyze_controls.py
"$PYTHON" scripts/analyze_library.py
"$PYTHON" scripts/audit_reused.py
"$PYTHON" scripts/build_ledger.py
"$PYTHON" scripts/build_frontier.py
"$PYTHON" scripts/control_sensitivity.py
"$PYTHON" scripts/plot_results.py
"$PYTHON" scripts/verify.py
"$PYTHON" scripts/make_manifest.py
