"""Run the independent verifier in a temporary, job-local dependency tree."""
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parent


def main():
    with TemporaryDirectory(prefix=".verify-", dir=ROOT) as temporary:
        packages = Path(temporary) / "packages"
        env = dict(os.environ, TMPDIR=temporary, TEMP=temporary, TMP=temporary,
                   OPENBLAS_NUM_THREADS="1", PYTHONNOUSERSITE="1",
                   PYTHONPATH=str(packages), PYTHONDONTWRITEBYTECODE="1",
                   PIP_CONFIG_FILE=os.devnull)
        pip = [sys.executable, "-B", "-m", "pip", "--isolated",
               "--disable-pip-version-check", "install", "--no-cache-dir",
               "--only-binary=:all:", "--use-deprecated=legacy-certs",
               "--target", str(packages), "-r", "requirements.txt"]
        subprocess.run(pip, cwd=ROOT, env=env, check=True)
        return subprocess.run(
            [sys.executable, "-B", "verify.py"], cwd=ROOT, env=env
        ).returncode


if __name__ == "__main__":
    raise SystemExit(main())
