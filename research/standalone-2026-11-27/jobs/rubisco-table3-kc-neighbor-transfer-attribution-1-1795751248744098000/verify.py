"""Recompute the fixed-k and attribution results independently from the cache."""
import csv
import json
import math
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent
data = json.loads((ROOT / 'sources/cached-kc.json').read_text())
result = json.loads((ROOT / 'results.json').read_text())
rows = data['rows']
y = [row['y'] for row in rows]
indices = set(range(len(rows)))
with (ROOT / 'sources/pred_KC.csv').open() as handle:
    archive = list(csv.DictReader(handle))
assert len(archive) == len(rows) == 197
for row, archived in zip(rows, archive):
    assert row['species'] == archived['species']
    assert abs(row['y'] - float(archived['y_true'])) < 1e-12


def close(actual, expected):
    assert math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12), (actual, expected)


def centered_ranks(values):
    ordered = sorted(values)
    return [ordered.index(v) + (ordered.count(v) - len(values)) / 2 for v in values]


def spearman(a, b):
    x, z = centered_ranks(a), centered_ranks(b)
    return sum(p*q for p, q in zip(x, z)) / math.sqrt(sum(p*p for p in x) * sum(q*q for q in z))


def replay(folds, k, ties=False):
    predictions, neighbors = {}, {}
    for fold in folds:
        n = k if k is not None else fold['chosen']['identity_knn'][0]
        for i in fold['test']:
            left = set(fold['train'])
            selected = []
            for _ in range(n):
                j = max(left, key=lambda j: (data['identity'][i][j], -j))
                left.remove(j)
                selected.append(j)
            if ties:
                selected.extend(j for j in left if data['identity'][i][j] == data['identity'][i][selected[-1]])
            assert i not in predictions
            predictions[i] = round(math.fsum(y[j] for j in selected) / len(selected), 10)
            neighbors[i] = selected
    assert predictions.keys() == indices
    return [predictions[i] for i in range(len(y))], neighbors


genus = [r['species'].casefold().split()[0] for r in rows]
for fold, (component, direct) in enumerate(zip(data['component'], data['direct'])):
    test = component['test']
    base = {j for j in indices-set(test)
            if all(genus[j] != genus[i] and data['identity'][i][j] < .95
                   and rows[j]['group'] != rows[i]['group'] for i in test)}
    assert test == direct['test'] and base == set(direct['train'])
    assert set(component['train']) < base
    for draw in data['thin']:
        arm = draw[fold]
        assert arm['test'] == test
        assert len(arm['train']) == len(component['train'])
        assert set(arm['train']) <= base

for key, metrics in result['fixed'].items():
    k = int(key)
    for i, arm in enumerate(('component', 'direct')):
        close(spearman(y, replay(data[arm], k)[0]), metrics[arm])
        close(spearman(y, replay(data[arm], k, True)[0]), metrics['tie_inclusive'][i])
    scores = [spearman(y, replay(draw, k)[0]) for draw in data['thin']]
    close(median(scores), metrics['matched_median'])
    close(max(scores), metrics['matched_max'])
assert set(result['fixed']) == {'1', '3', '5', '10'}
for arm in ('component', 'direct'):
    predictions, _ = replay(data[arm], None, True)
    assert predictions == data['prior_knn'][arm]
    close(spearman(y, predictions), result['tuned'][arm])
close(spearman(y, data['trimer']), result['direct_trimer'])

c, cn = replay(data['component'], 1)
d, dn = replay(data['direct'], 1)
rd = spearman(y, d)
gap = spearman(y, c) - rd
close(gap, result['gap'])
transfers = {t['test']: t for t in result['traces']}
assert len(transfers) == len(result['traces'])
assert transfers.keys() == {i for i in indices if cn[i] != dn[i]}
for i, t in transfers.items():
    a, b = cn[i][0], dn[i][0]
    assert t['component'] == a and t['direct'] == b
    for key, value in [('true_log_kc', y[i]), ('component_log_kc', c[i]), ('direct_log_kc', d[i]),
                       ('component_identity', data['identity'][i][a]),
                       ('direct_identity', data['identity'][i][b])]:
        close(value, t[key])
    fold = next(f for f, arm in enumerate(data['component']) if i in arm['test'])
    assert b not in data['component'][fold]['train']

groups = {}
for i, row in enumerate(rows):
    groups.setdefault(row['group'], []).append(i)


def swap(selected):
    return spearman(y, [c[i] if row['group'] in selected else d[i]
                        for i, row in enumerate(rows)]) - rd


assert len(result['omissions']) == len(groups)
assert {o['group'] for o in result['omissions']} == groups.keys()
for record in result['omissions']:
    drop = groups[record['group']]
    assert drop == record['rows']
    a, b, z = ([v[i] for i in range(len(y)) if i not in drop] for v in (y, c, d))
    close(spearman(a, b) - spearman(a, z), record['gap'])
ranking = sorted((g for g, members in groups.items() if any(c[i] != d[i] for i in members)),
                 key=lambda g: (-swap({g}), g))
assert ranking == [h['group'] for h in result['hybrids']]
for record in result['hybrids']:
    close(swap({record['group']}), record['gain'])
for n, gain in zip((1, 3, 5), result['cumulative'], strict=True):
    close(swap(set(ranking[:n])), gain)
assert ranking[0] == min(result['omissions'], key=lambda o: o['gap'])['group']
assert all(o['gap'] > 0 for o in result['omissions'])
assert result['fixed']['1']['component'] > result['fixed']['1']['matched_max']
assert result['fixed']['3']['component'] > result['fixed']['3']['matched_max']
print(f'Verified {len(transfers)} changed donors, all group omissions, label swaps and fixed-k ranks.')
