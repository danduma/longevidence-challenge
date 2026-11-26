#!/bin/sh
set -eu
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
exec "${C07_PYTHON:-python3}" -B verify.py
