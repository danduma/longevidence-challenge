"""Summarize full-horizon and 12 s paired gaps from archived stress traces."""
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


def scores(data, region, early=False):
    values = np.array(data[region])
    if early:
        values = values[np.array(data['times']) <= 12]
    return np.array([values[:,0].max(),20*logsumexp(values[:,1]/20)])


if __name__ == '__main__':
    result = json.loads(Path(__file__).with_name('results.json').read_text())
    fraction = result['retained_weight']
    for run in result['runs']:
        print(run['family'],run['cells'],run['q'],run['conditioning'],run['transfer'])
        for early in (False,True):
            by_region = {}
            for region in ('all','interior'):
                pair = {name:scores(trace,region,early) for name,trace in run['scores'].items()}
                by_region[region] = pair
                gap = pair['graded']-pair['uniform_strict']
                lead = np.minimum(pair['graded'],pair['uniform_strict'])-pair['uniform_released_policy']
                print(early,region,'gap',gap,'tuned lead',lead,'tensile peak',pair['graded'][0])
            print('retained exp fraction', [fraction*np.exp((by_region['interior'][name][1]-
                    by_region['all'][name][1])/20) for name in by_region['all']])
