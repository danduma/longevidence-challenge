#!/bin/sh
set -eu
cd "$(dirname "$0")"
export UV_CACHE_DIR="$PWD/.runtime/cache"
export UV_LINK_MODE=copy
export UV_PYTHON_DOWNLOADS=never
if [ ! -x .runtime/venv/bin/python ]; then
    uv venv --python 3.14 .runtime/venv
fi
uv pip sync --python .runtime/venv/bin/python requirements.txt
exec .runtime/venv/bin/python -B verify.py "$@"
