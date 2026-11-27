#!/bin/sh
set -eu
cd "$(dirname "$0")"
python3 -m venv .venv
PIP_NO_CACHE_DIR=1 .venv/bin/python -m pip --use-deprecated=legacy-certs --disable-pip-version-check install -r requirements.txt
.venv/bin/python -B verify.py
