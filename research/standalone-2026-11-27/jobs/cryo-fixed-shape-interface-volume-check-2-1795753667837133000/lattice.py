"""Fixed schedules extracted from the pinned release and predecessor."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "sources/upstream/src"))

from vitrify import geometry

release = ROOT / "sources/upstream/results/grid_convergence_margin.json"
VECTORS = {"graded": json.loads(release.read_text())["design_parameters"]}
for name, (loading, *duty) in json.loads((ROOT / "candidates.json").read_text()).items():
    VECTORS[name] = [10 * loading] * 3 + duty
