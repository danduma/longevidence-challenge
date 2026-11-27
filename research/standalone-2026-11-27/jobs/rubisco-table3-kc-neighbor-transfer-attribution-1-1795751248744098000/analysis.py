"""Attribute a fixed-fold pooled KC rank contrast from cached identities."""
import json
from itertools import groupby
from pathlib import Path
from statistics import correlation, fmean, median

ROOT = Path(__file__).resolve().parent


def ranks(values):
    result = [0.] * len(values)
    ordered = sorted(range(len(values)), key=values.__getitem__)
    start = 0
    for _, group in groupby(ordered, key=values.__getitem__):
        tied = list(group)
        for index in tied:
            result[index] = start + (len(tied) + 1) / 2
        start += len(tied)
    return result


def rho(a, b):
    return correlation(ranks(a), ranks(b))


def predict(data, folds, k, include_ties=False):
    neighbors = [None] * len(data['rows'])
    for fold in folds:
        count = fold['chosen']['identity_knn'][0] if k is None else k
        for row in fold['test']:
            identity = data['identity'][row]
            ordered = sorted(fold['train'], key=lambda j: (-identity[j], j))
            cutoff = identity[ordered[count - 1]]
            neighbors[row] = ([j for j in ordered if identity[j] >= cutoff]
                              if include_ties else ordered[:count])
    y = [row['y'] for row in data['rows']]
    return [round(fmean(y[j] for j in selected), 10) for selected in neighbors], neighbors


def audit(data):
    y = [row['y'] for row in data['rows']]
    aggregate = {}
    for k in (1, 3, 5, 10):
        component, _ = predict(data, data['component'], k)
        direct, _ = predict(data, data['direct'], k)
        matched = [rho(y, predict(data, folds, k)[0]) for folds in data['thin']]
        aggregate[k] = dict(component=rho(y, component), direct=rho(y, direct),
                            matched_median=median(matched), matched_max=max(matched),
                            tie_inclusive=[rho(y, predict(data, data[arm], k, True)[0])
                                           for arm in ('component', 'direct')])
    c, ct = predict(data, data['component'], 1)
    d, dt = predict(data, data['direct'], 1)
    rc, rd = rho(y, c), rho(y, d)
    groups = {row['group'] for row in data['rows']}
    direct_only = {}
    for cf, df in zip(data['component'], data['direct']):
        extra = set(df['train']) - set(cf['train'])
        for i in cf['test']:
            direct_only[i] = extra
    traces = []
    for i, (a, b) in enumerate(zip(ct, dt)):
        if a != b:
            assert b[0] in direct_only[i]
            traces.append(dict(test=i, component=a[0], direct=b[0],
                               true_log_kc=y[i], component_log_kc=c[i], direct_log_kc=d[i],
                               component_identity=data['identity'][i][a[0]],
                               direct_identity=data['identity'][i][b[0]]))
    omissions = []
    hybrids = []
    for group in sorted(groups):
        drop = [i for i, row in enumerate(data['rows']) if row['group'] == group]
        keep = [i for i in range(len(y)) if i not in drop]
        keep_y, keep_c, keep_d = ([v[i] for i in keep] for v in (y, c, d))
        omissions.append(dict(group=group, rows=drop,
                              gap=rho(keep_y, keep_c) - rho(keep_y, keep_d)))
        if any(c[i] != d[i] for i in drop):
            swapped = [c[i] if i in drop else d[i] for i in range(len(y))]
            hybrids.append(dict(group=group, gain=rho(y, swapped) - rd))
    hybrids.sort(key=lambda item: (-item['gain'], item['group']))
    cumulative = []
    for count in (1, 3, 5):
        selected = {item['group'] for item in hybrids[:count]}
        swap = [c[i] if row['group'] in selected else d[i]
                for i, row in enumerate(data['rows'])]
        cumulative.append(rho(y, swap) - rd)
    return dict(fixed=aggregate, traces=traces, omissions=omissions, hybrids=hybrids,
                cumulative=cumulative, gap=rc-rd, direct_trimer=rho(y, data['trimer']),
                tuned={arm:rho(y, predict(data, data[arm], None, True)[0])
                       for arm in ('component', 'direct')})


if __name__ == '__main__':
    data = json.loads((ROOT / 'sources/cached-kc.json').read_text())
    (ROOT / 'results.json').write_text(json.dumps(audit(data), indent=2, allow_nan=False) + '\n')
