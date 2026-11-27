"""Check whether filled chambers determine the stop on either thermal graph."""
import sys
from analysis import heat_history, voxel_graph, cut_graph, VECTORS

method, step = sys.argv[1:]
dt = float(step)
graph = voxel_graph(.75) if method == 'voxel' else cut_graph(.75,16 if dt==.1 else 24)
for mode in ('extrapolated', 'clamped'):
    for name, vector in VECTORS.items():
        history, event = heat_history(graph,vector,mode,dt)
        phase = graph[0]
        myocardium = history[-1,phase==2].min()
        chamber = history[-1,phase==3].min()
        assert abs(myocardium+35) < 1e-10
        assert chamber > -35
        print(method, dt, mode, name, f'chamber minimum {chamber:.6f} C', flush=True)
