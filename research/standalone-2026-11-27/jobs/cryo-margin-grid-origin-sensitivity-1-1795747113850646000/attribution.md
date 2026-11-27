# Attribution and reproducibility

Marc Davici's Vitrify files in `sources/upstream` are unmodified snapshots at commit `1d498ee5c583e0112a762e60892d576ec2e27f57` under the included Apache-2.0 license. The NumPy/SciPy solver in `forward.py` and uniform candidates in `candidates.json` are reused from the coordinator's preceding uniform-surround comparison. The solver independently implements the pinned finite-volume and Maxwell equilibrium equations. The translation adapter is in `analysis.py`.

Joshi and Rabin's 2023 paper snapshot retains author, DOI and Creative Commons attribution metadata. The 2022 article and forum files contain attributed excerpts and reading notes. GitHub API snapshots were retrieved for the preceding comparison on the same date. `protocol.json` records the phase set and metrics fixed before the first translated run.

Run `sh verify.sh` with uv and an existing Python 3.14 installation, or `python3 -B verify.py` in an environment with the packages in `requirements.txt`. The shell script keeps its environment and package cache in this job. A run uses two workers by default. Set `CRYO_WORKERS=1` to reduce peak memory.
