"""Generate the matched cut-cell trajectories from the pinned schedules."""
import json
import numpy as np
from coupled import heat_history
from lattice import VECTORS
from thermal_graph import cut_graph

DT = .05
H = .75
Q = 24


def compute():
    graph = cut_graph(H, Q)
    histories, events = {}, {}
    for name, vector in VECTORS.items():
        histories[name], events[name] = heat_history(graph, vector, 'extrapolated', DT)
        print(name, events[name], flush=True)
    return graph, histories, events


if __name__ == '__main__':
    _, histories, events = compute()
    np.savez_compressed('thermal-cache.npz', **histories)
    with open('events.json', 'w') as handle:
        json.dump(events, handle, indent=2)
