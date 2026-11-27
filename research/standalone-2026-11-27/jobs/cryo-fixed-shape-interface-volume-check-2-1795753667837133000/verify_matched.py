"""Challenge the startup comparison using the same update for every model."""
import argparse
import json

import numpy as np

from analysis import ROOT, cut_graph, voxel_graph
from lattice import VECTORS
from verify import independent


def recompute():
    vector = VECTORS["uniform_strict"]
    temperature = -150.
    # At uniform initial temperature the physical conductive flux is zero.
    capacity = (2402.3 + 7.6 * temperature) * (1078 - .43 * temperature)
    source_only = 60 * 1000 * 1.47 * (323 - 4.6 * temperature) * vector[3] / capacity
    timestep = .0005
    cases = {
        "voxel_1.1": lambda: voxel_graph(1.1),
        "gap_1.1": lambda: voxel_graph(1.1, gap=True),
        "cut_1.1_q24": lambda: cut_graph(1.1, 24),
        "voxel_0.75": lambda: voxel_graph(.75),
    }
    rates = {name: independent(build(), vector, timestep, .5) for name, build in cases.items()}
    assert rates["voxel_1.1"] < 51.5 < min(rates[name] for name in cases if name != "voxel_1.1")
    assert 0 < source_only - rates["cut_1.1_q24"] < .001
    return {
        "dt_s": timestep,
        "diffusion_theta": .5,
        "initial_myocardial_source_only_c_per_min": source_only,
        "minima_c_per_min": rates,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="regenerate the matched-update record")
    args = parser.parse_args()
    record = recompute()
    path = ROOT / "matched-results.json"
    if args.write:
        path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    archived = json.loads(path.read_text())
    assert archived.keys() == record.keys()
    assert archived["minima_c_per_min"].keys() == record["minima_c_per_min"].keys()
    for key in ("dt_s", "diffusion_theta", "initial_myocardial_source_only_c_per_min"):
        np.testing.assert_allclose(record[key], archived[key], rtol=1e-12)
    for name, value in record["minima_c_per_min"].items():
        np.testing.assert_allclose(value, archived["minima_c_per_min"][name], atol=2e-5, rtol=0)
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
