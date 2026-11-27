"""Compute equal-loading window predictions from the publisher source workbook."""
import json
from itertools import combinations, count
from pathlib import Path
import re
from statistics import linear_regression, mean, median
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from workbook import worksheets

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'sources'
EDGES = tuple(range(-140, 1, 20))


def trace(sheet, time_col, temp_col):
    points = []
    for row in count(6):
        temp = sheet.get(f'{temp_col}{row}')
        if temp is None:
            break
        points.append((float(sheet[f'{time_col}{row}']), float(temp)))
    if not points or any(a[0] >= b[0] for a, b in zip(points, points[1:])):
        raise ValueError(f'Invalid time series in {temp_col}')
    return points


def transit(points, edges):
    crossings = []
    for edge in edges:
        for (t0, y0), (t1, y1) in zip(points, points[1:]):
            if y0 < edge <= y1:
                crossings.append(t0 + (t1 - t0) * (edge - y0) / (y1 - y0))
                break
        else:
            raise ValueError(f'No upward crossing at {edge}')
    return [b - a for a, b in zip(crossings, crossings[1:])]


def trends(profiles):
    points = [(j, value) for profile in profiles for j, value in enumerate(profile)]
    # A least-absolute-deviation optimum exists at an intersection of residual lines.
    candidates = []
    for (i, a), (j, b) in combinations(points, 2):
        if i != j:
            slope = (b - a) / (j - i)
            candidates.append((slope, a - slope * i))
    losses = [(sum(abs(value - (m * j + c)) for j, value in points), m, c)
              for m, c in candidates]
    best_loss = min(row[0] for row in losses)
    # Break equivalent fits toward the smallest slope and then intercept.
    lad = min((m, c) for loss, m, c in losses if loss <= best_loss + 1e-10)
    ols = linear_regression(*zip(*points))
    return [[m * j + c for j in range(len(profiles[0]))] for m, c in (lad, ols)]


def predictions(small, fixed):
    half, one = small
    pooled = list(map(mean, zip(half, one)))
    durations = half + one
    n = len(one)
    nearest_lad, nearest_ols = trends([one])
    pooled_lad, pooled_ols = trends(small)
    return {
        'nearest_linear_LAD': nearest_lad, 'nearest_linear_OLS': nearest_ols,
        'pooled_linear_LAD': pooled_lad, 'pooled_linear_OLS': pooled_ols,
        'half_liter': half, 'nearest_1L': one, 'pooled': pooled,
        'pooled_median_constant': [median(durations)] * n,
        'pooled_mean_constant': [mean(durations)] * n,
        'nearest_median_constant': [median(one)] * n,
        'nearest_mean_constant': [mean(one)] * n,
        'fixed_warm': [fixed] * n,
    }


def errors(actual, models):
    return {name: mean(abs(y - p) for y, p in zip(actual, pred, strict=True))
            for name, pred in models.items()}


def evaluate(actual, small, fixed):
    models = predictions(small, fixed)
    scores = errors(actual, models)
    # MAE is convex in a common blend weight, so its maximum is at an endpoint.
    worst_blend = max(scores['half_liter'], scores['nearest_1L'])
    best_constant = min(error for key, error in scores.items()
                        if key.endswith('_constant') or key == 'fixed_warm')
    return {'MAE_seconds': scores, 'predicted_seconds': models,
            'worst_blend_MAE_seconds': worst_blend,
            'worst_blend_margin_vs_best_constant_seconds': best_constant - worst_blend}


def exact_run(a, b):
    """Longest contiguous equal-temperature sequence at equal elapsed times."""
    longest = 0
    for offset in range(1 - len(b), len(a)):
        streak = 0
        for p, q in zip(a[max(0, offset):], b[max(0, -offset):]):
            streak = streak + 1 if abs(p[1] - q[1]) < 1e-9 else 0
            longest = max(longest, streak)
    return longest


def calculate():
    book = dict(worksheets(SOURCE / 'liter-data.xlsx'))
    sup, fig = book['Supplementary Information'], book['Fig.7']
    with ZipFile(SOURCE / 'liter-supp.docx') as z:
        xml = ET.fromstring(z.read('word/document.xml'))
    words = ''.join(xml.itertext())
    density = float(re.search(r'~(\d+) kg/m3 for M22 samples', words)[1])
    cp = float(re.search(r'~(\d+) J/kg.K for M22 samples', words)[1])
    sar = mean(float(sup[cell]) for cell in ('CX22', 'CY22', 'CZ22'))
    concentration = float(sup['DR5'].split('_')[1].split()[0])
    fixed_rate = sar * concentration * 1000 / (density * cp)
    small = [trace(sup, 'DQ', col) for col in ('DR', 'DS')]
    train = [transit(points, EDGES) for points in small]
    test = trace(sup, 'DQ', 'DT')
    target = transit(test, EDGES)
    probes = {fig[col+'5']: trace(fig, 'F', col) for col in 'GHI'}
    fixed = 20 / fixed_rate
    full = evaluate(target, train, fixed)
    return {
        'edges_C': EDGES, 'fixed_rate_C_per_min': fixed_rate * 60,
        'training_seconds': train, 'S22A_2L_seconds': target,
        'S22A_full': full,
        'S22A_warm': evaluate(target[2:], [row[2:] for row in train], fixed),
        'Fig7b_2L_probe_MAE_seconds': {
            name: errors(transit(points, EDGES), full['predicted_seconds'])
            for name, points in probes.items()
        },
        'max_contiguous_exact_temperature_match_S22A_vs_Fig7b': {
            col: {p: exact_run(trace(sup, 'DQ', col), trace(fig, 'F', p)) for p in 'GHI'}
            for col in ('DR', 'DS', 'DT')
        },
        'Fig7c_2L_triplet_means_C_per_min': [
            mean(float(fig[col+'7']) for col in cols) for cols in ('MPS', 'NQT', 'ORU')],
        'S22B_4p6_triplet_C_per_min': [float(sup[col+'7']) for col in ('EC', 'ED', 'EE')],
        'S22A_total_minus100_to_0_seconds': [sum(transit(points, EDGES[2:]))
                                             for points in [*small, test]],
    }


if __name__ == '__main__':
    (ROOT / 'results.json').write_text(json.dumps(calculate(), indent=2) + '\n')
