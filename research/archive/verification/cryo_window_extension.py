"""Coordinator sensitivity to rate-screen temperature support, not new physics.

Uses the reviewed worker's flux equations with parameterized screen bounds.
Original source bytes and outputs remain untouched. The broader −123 to 0 C
interval follows Joshi/Rabin's glass-to-end guidance (PMC8823202); the
−100 to −40 C interval is Gao's highlighted concern region. Neither screen
models ice kinetics or establishes tissue survival.
"""
from pathlib import Path
import ctypes
import hashlib
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
JOB = HERE.parent / "agents/cryo-robust-astra"
OUT = HERE / "cryo-window-extension"
OUT.mkdir(exist_ok=True)
source = (JOB / "thermal.cpp").read_text()
assert source.count("int chunk,double *out") == 1
assert source.count("t>=-100 && t<=-40") == 1
modified = source.replace("int chunk,double *out", "int chunk,double lower,double upper,double *out")
modified = modified.replace("t>=-100 && t<=-40", "t>=lower && t<=upper")
(OUT / "thermal_window.cpp").write_text(modified)
subprocess.run(["clang++", "-O3", "-std=c++17", "-shared", "-fPIC", str(OUT / "thermal_window.cpp"), "-o", str(OUT / "thermal_window.so")], check=True)
sys.path.insert(0, str(JOB))
os.chdir(JOB)
from investigate import setup

f = ctypes.CDLL(str(OUT / "thermal_window.so")).solve
ip = np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS")
dp = np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")
I, D = ctypes.c_int, ctypes.c_double
f.argtypes = [I, I, ip, ip, dp, D, D, I, dp, I, D, I, D, D, dp, dp]
f.restype = None
rows = []
for design in ["margin", "equal"]:
    for scenario in ["nominal", "high", "deepLow", "deepHigh"]:
        for sar in ["extrapolated", "clamped"]:
            p, faces, labels, iron, dt, duty = setup(.5, design, scenario, 1, 70)
            for lo, hi in [(-150., 0.), (-123., 0.), (-123., -35.), (-100., -40.)]:
                values = np.zeros(15)
                end = np.zeros(len(iron))
                f(len(iron), len(faces), faces, labels, iron, p.spacing_m, dt, len(duty), duty,
                  int(sar == "clamped"), 0., int(np.ceil(len(duty) / 120)), lo, hi, values, end)
                rate = values[4]
                rows.append(dict(design=design, scenario=scenario, sar=sar, mm=.5,
                                 lower=lo, upper=hi, rate=rate, min_time=values[5], min_temp=values[6],
                                 cold=values[8], hot=values[9],
                                 passes=bool(rate>=50 and values[8]>=-35 and values[9]<=0)))
            pd.DataFrame(rows).to_csv(OUT / "windows.csv", index=False)
df = pd.DataFrame(rows)
# The extended bounds must reduce exactly to the existing narrow-screen data.
original = pd.read_csv(JOB / "outputs/sweep.csv").query("mm==.5 and dt_divisor==1")
merged = df.query("lower==-100").merge(original, on=["design", "scenario", "sar", "mm"])
assert len(merged) == 16
assert np.max(np.abs(merged.rate - merged.window_rate)) < 1e-9
wide = df.query("lower==-150").merge(original, on=["design", "scenario", "sar", "mm"])
assert np.max(np.abs(wide.rate - wide.min_rate)) < 1e-9
comparison = df.pivot(index=["scenario", "sar", "lower", "upper"], columns="design", values="rate")
comparison["margin_minus_equal"] = comparison.margin - comparison.equal
comparison.to_csv(OUT / "rankings.csv")
(OUT / "provenance.json").write_text(json.dumps({
    "base_source_sha256": hashlib.sha256(source.encode()).hexdigest(),
    "modified_source_sha256": hashlib.sha256(modified.encode()).hexdigest(),
    "scope": "Same numerical equations; temperature-support sensitivity only. Not an independent physics solver.",
    "checks": "All 16 narrow-window and all 16 all-time outputs match the worker sweep within 1e-9."
}, indent=2) + "\n")
print(comparison.to_string())
