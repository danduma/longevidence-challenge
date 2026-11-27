#!/bin/sh
set -eu
cd "$(dirname "$0")"
export UV_CACHE_DIR="$PWD/.runtime/cache"
export UV_LINK_MODE=copy
export UV_PYTHON_DOWNLOADS=never
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
if [ ! -x .runtime/venv/bin/python ]; then
    uv venv --python 3.14 .runtime/venv
fi
uv pip sync --python .runtime/venv/bin/python requirements.txt
case "${1:-main}" in
    main) script=verify.py ;;
    margin) script=margin.py ;;
    *) echo "Usage: sh verify.sh [main|margin]" >&2; exit 2 ;;
esac
exec .runtime/venv/bin/python -B "$script"
