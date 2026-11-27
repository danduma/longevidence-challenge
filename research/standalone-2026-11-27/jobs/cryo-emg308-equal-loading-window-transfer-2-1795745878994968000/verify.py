"""Independently reconstruct all reported values using exact rational arithmetic."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import isclose
from pathlib import Path
from statistics import median, median_low
import re
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from workbook import worksheets

ROOT = Path(__file__).resolve().parent
S = ROOT / 'sources'


def require(condition, detail):
    if not condition:
        raise AssertionError(detail)


def check(saved, expected):
    if isinstance(expected, dict):
        require(saved.keys() == expected.keys(), 'Output keys differ')
        for key in expected:
            check(saved[key], expected[key])
    elif isinstance(expected, list):
        for got, want in zip(saved, expected, strict=True):
            check(got, want)
    else:
        require(isclose(saved, float(expected), rel_tol=1e-12, abs_tol=1e-12),
                f'Expected {expected}, found {saved}')


def average(values):
    values = list(values)
    return sum(values) / len(values)


def samples(sheet, timer, sensor):
    rows = sorted(int(m[1]) for key in sheet
                  if (m := re.fullmatch(sensor + r'(\d+)', key)) and int(m[1]) >= 6)
    require(rows == list(range(6, 6 + len(rows))), 'Trace contains missing rows')
    points = [(F(sheet[timer + str(i)]), F(sheet[sensor + str(i)])) for i in rows]
    require(all(b[0] - a[0] == 1 for a, b in zip(points, points[1:])), 'Unexpected cadence')
    return points


def durations(points, edges):
    times, temperatures = zip(*points)
    crossings = []
    for edge in edges:
        idx = next(i for i, temperature in enumerate(temperatures) if temperature >= edge)
        require(idx > 0, 'Crossing precedes observed trace')
        fraction_remaining = (temperatures[idx] - edge) / (temperatures[idx] - temperatures[idx-1])
        crossings.append(times[idx] - fraction_remaining * (times[idx] - times[idx-1]))
    return [right - left for left, right in zip(crossings, crossings[1:])]


def longest_shared(a, b):
    # Longest common contiguous substring on the 0.1 °C measurement grid.
    previous = [0] * (len(b) + 1)
    longest = 0
    for _, temp in a:
        current = [0]
        for j, (_, other) in enumerate(b):
            current.append(1 + previous[j] if round(temp * 10) == round(other * 10) else 0)
        longest = max(longest, max(current))
        previous = current
    return longest


def trend(profiles):
    points = [(i, y) for row in profiles for i, y in enumerate(row)]
    slopes = {(b - a) / (j - i) for (i, a), (j, b) in combinations(points, 2) if i != j}
    fits = [(m, median_low(y - m * i for i, y in points)) for m in slopes]
    losses = [(sum(abs(y - m * i - c) for i, y in points), m, c) for m, c in fits]
    best = min(row[0] for row in losses)
    lad = min((m, c) for loss, m, c in losses if loss - best <= F(1, 10**10))
    xbar, ybar = map(average, zip(*points))
    slope = sum((i - xbar) * (y - ybar) for i, y in points) / sum((i - xbar)**2 for i, _ in points)
    return [[m * i + c for i in range(len(profiles[0]))]
            for m, c in (lad, (slope, ybar - slope * xbar))]


def evaluate(actual, predictions):
    return {model: average(abs(p - y) for p, y in zip(values, actual, strict=True))
            for model, values in predictions.items()}


def reconstruct(a, b, actual, fixed):
    preds = {'half_liter': a, 'nearest_1L': b, 'pooled': [average(pair) for pair in zip(a, b)],
             'fixed_warm': [fixed] * len(a)}
    for label, profiles in (('pooled', [a, b]), ('nearest', [b])):
        samples = sum(profiles, [])
        for stat in (median, average):
            name = 'mean' if stat is average else 'median'
            preds[f'{label}_{name}_constant'] = [stat(samples)] * len(a)
        for loss, row in zip(('LAD', 'OLS'), trend(profiles)):
            preds[f'{label}_linear_{loss}'] = row

    # Enumerate residual zeroes and endpoints to verify the full blend bound.
    candidates = {F(0), F(1)}
    candidates.update((y - x) / (z - x) for x, z, y in zip(a, b, actual) if z != x)
    worst = max(average(abs(y - (x + w * (z - x))) for x, z, y in zip(a, b, actual))
                for w in candidates if 0 <= w <= 1)
    errors = evaluate(actual, preds)
    strongest_constant = min(err for name, err in errors.items()
                             if name.endswith('_constant') or name == 'fixed_warm')
    require(worst < strongest_constant, 'A blend lost to a tested constant')
    return {'predicted_seconds': preds, 'MAE_seconds': errors,
            'worst_blend_MAE_seconds': worst,
            'worst_blend_margin_vs_best_constant_seconds': strongest_constant - worst}


def main():
    workbook = S / 'liter-data.xlsx'
    require(sha256(workbook.read_bytes()).hexdigest() ==
            'c059e1beef528efd5a3aa04fc51d039d64e1c8583376a3596160e837069d2898', 'Workbook hash')
    sheets = dict(worksheets(workbook))
    s, f = sheets['Supplementary Information'], sheets['Fig.7']
    require([s[col+'5'] for col in ('DR', 'DS', 'DT')] ==
            [v + 'L_4.6 mgFe/mL' for v in ('0.5', '1', '2')], 'Volume/loading mapping')
    require(s['DQ2'] == 'Fig. S22A' and s['DT5'].split('_')[1] == f['K7'], 'Figure mapping')
    require(s['CX21'] == 'EMG308 in M22 without LM5' and s['CT22'] == '35', 'Assay mapping')
    supplement = S / 'liter-supp.docx'
    require(sha256(supplement.read_bytes()).hexdigest() ==
            '63336d3e3e002e1bfda535e3609aa2758e89a3c360176a79e28b324188ee4059', 'Supplement hash')
    with ZipFile(supplement) as archive:
        doc = ET.fromstring(archive.read('word/document.xml'))
    paragraph = next(''.join(p.itertext()) for p in doc.iter() if p.tag.endswith('}p')
                     and 'kg/m3 for M22 samples' in ''.join(p.itertext()))
    rho, cp = [F(x) for x in re.findall(r'~(\d+) (?:kg/m3|J/kg.K) for M22 samples', paragraph)]
    sar = sum(F(s[x+'22']) for x in ('CX', 'CY', 'CZ')) / 3
    seconds = 20 * rho * cp / (sar * F(s['DR5'].split('_')[1].split()[0]) * 1000)

    traces = {col: samples(s, 'DQ', col) for col in ('DR', 'DS', 'DT')}
    probes = {col: samples(f, 'F', col) for col in 'GHI'}
    edges = [-140, -120, -100, -80, -60, -40, -20, 0]
    a, b, actual = [durations(t, edges) for t in traces.values()]
    full = reconstruct(a, b, actual, seconds)
    warm = reconstruct(a[2:], b[2:], actual[2:], seconds)
    require(full['MAE_seconds']['nearest_1L'] == min(full['MAE_seconds'].values()), 'Full-span winner')
    require(warm['MAE_seconds']['nearest_linear_LAD'] == min(warm['MAE_seconds'].values()),
            'Warm-span winner')
    parallel = {f[col+'5']: evaluate(durations(t, edges), full['predicted_seconds'])
                for col, t in probes.items()}
    reuse = {col: {p: longest_shared(t, u) for p, u in probes.items()} for col, t in traces.items()}
    require(max(reuse['DT'].values()) == 3, 'Trace reuse pattern changed')
    triplets = [average(F(f[col+'7']) for col in cols) for cols in ('MPS', 'NQT', 'ORU')]
    supplemental = [F(s[k+'7']) for k in ('EC', 'ED', 'EE')]
    require([round(v, 5) for v in triplets] == [round(v, 5) for v in supplemental], 'Rate reuse pattern changed')
    check(json.loads((ROOT / 'results.json').read_text()), {
        'edges_C': edges,
        'fixed_rate_C_per_min': 1200 / seconds,
        'training_seconds': [a, b],
        'S22A_2L_seconds': actual,
        'S22A_full': full,
        'S22A_warm': warm,
        'Fig7b_2L_probe_MAE_seconds': parallel,
        'max_contiguous_exact_temperature_match_S22A_vs_Fig7b': reuse,
        'Fig7c_2L_triplet_means_C_per_min': triplets,
        'S22B_4p6_triplet_C_per_min': supplemental,
        'S22A_total_minus100_to_0_seconds': [sum(row[2:]) for row in (a, b, actual)],
    })
    print('Verified independent crossings, predictor ranking, and figure reuse against the workbook.')


if __name__ == '__main__':
    main()
