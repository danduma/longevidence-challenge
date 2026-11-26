"""Fresh coordinator execution of the separate reviewer's thermal implementation.

Reuses the documented benchmark geometry/parameters; integrates independently
of the Astra solver. All new outputs stay in verification/. No model calls.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
JOB = HERE.parent / "agents/review-cryo-sol"
OUT = HERE / "cryo-independent-reproduction"
OUT.mkdir(exist_ok=True)
subprocess.run(["c++", "-O3", "-std=c++17", str(JOB / "independent.cpp"), "-o", str(OUT / "independent")], check=True)
spec = importlib.util.spec_from_file_location("independent_cryo_review", JOB / "run_review.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.ROOT = OUT
module.OUT = OUT
rows = []
for mode, mean in [("clamped", 1.83), ("extrapolated", 1.47)]:
    for design in ["ramp", "flat"]:
        rows.append(module.run(.5, design, mean=mean, mode=mode))
ramp, flat, nominal_ramp, nominal_flat = rows
assert not ramp["pass"] and flat["pass"]
assert flat["window_rate"] < ramp["window_rate"]
assert nominal_flat["min_rate"] > nominal_ramp["min_rate"]
assert nominal_flat["window_rate"] < nominal_ramp["window_rate"]
assert all(r["missing"] == 0 and abs(r["energy_rel_error"]) < 1e-6 for r in rows)
result = {
    "command": "verification/cryo_independent_reproduction.py with shared scientific Python",
    "source_sha256": {n: hashlib.sha256((JOB / n).read_bytes()).hexdigest()
                      for n in ["independent.cpp", "run_review.py"]},
    "checks": "Four exact-horizon independent RK2 solves preserve conditional thermal feasibility and both schedule-ranking reversals; all crossings complete and energy balances pass.",
    "limitations": "Same prescribed geometry, properties and criteria; no independent physical or biological validation.",
    "results": rows,
}
(OUT / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
print("PASS: independent clamped feasibility and extrapolated/clamped rate-screen reversals.")
