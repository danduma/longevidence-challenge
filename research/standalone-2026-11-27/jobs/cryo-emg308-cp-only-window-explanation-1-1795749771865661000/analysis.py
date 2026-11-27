"""Fit a common inverse power to smaller bags, with no 2 L data in fitting."""
import csv
from itertools import combinations, pairwise
import json
from pathlib import Path
import re
from statistics import mean, median
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from workbook import worksheets

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'sources'


def properties():
    with ZipFile(SOURCE / 'liter-supp.docx') as archive:
        tree = ET.fromstring(archive.read('word/document.xml'))
    rows = [
        [''.join(n.text or '' for n in cell.findall('.//{*}t'))
         for cell in row.findall('{*}tc')]
        for row in tree.findall('.//{*}tr')
    ]
    row, = [row for row in rows if row[:1] == ['M22'] and len(row) == 7]
    knots = sorted((int(temp), float(cp) * 1000)
                   for cp, temp in re.findall(r'(\d+\.\d+)\s*\((-?\d+)⁰C\)', row[2]))
    fixed_density = float(row[3].split('[')[0])
    with (SOURCE / 'solanki-table2.csv').open() as file:
        density = next(csv.DictReader(file))
    return knots, fixed_density, density


def cp_at(knots, t):
    for (x, y), (u, v) in pairwise(knots):
        if x <= t <= u:
            return y + (v - y) * (t - x) / (u - x)
    raise ValueError(f'Temperature outside specific heat support: {t}')


def integral(knots, density, lo, hi):
    # Two-point Gauss quadrature is exact for a cubic on each knot interval.
    cuts = [lo, *(x for x, _ in knots if lo < x < hi), hi]
    total = 0
    for a, b in pairwise(cuts):
        midpoint, radius = (a + b) / 2, (b - a) / 2
        for t in (midpoint - radius / 3**0.5, midpoint + radius / 3**0.5):
            total += radius * cp_at(knots, t) * density(t)
    return total / 1e6  # MJ/m³


def duration(sheet, col, edges):
    points = [(float(sheet[f'DQ{i}']), float(sheet[f'{col}{i}']))
              for i in range(6, 1000) if f'{col}{i}' in sheet]
    crossings = [
        next(t + (u - t) * (edge - y) / (z - y)
             for (t, y), (u, z) in pairwise(points) if y < edge <= z)
        for edge in edges
    ]
    return [b - a for a, b in pairwise(crossings)]


def score(small, target, energy):
    scale = sum(map(sum, small)) / (len(small) * sum(energy))
    xs = range(len(energy))
    trend_options = [
        [a + (b - a) * (x - i) / (j - i) for x in xs]
        for i, j in combinations(xs, 2)
        for a, b in [(small[1][i], small[1][j])]
    ]

    def mae(actual, predicted):
        return mean(abs(a - p) for a, p in zip(actual, predicted, strict=True))

    cp = [scale * e for e in energy]
    predictions = {
        'heat_capacity': cp,
        'pooled_mean_constant': [mean(small[0] + small[1])] * len(energy),
        'pooled_median_constant': [median(small[0] + small[1])] * len(energy),
        'nearest_median_constant': [median(small[1])] * len(energy),
        'nearest_mean_constant': [mean(small[1])] * len(energy),
        'nearest_1L': small[1],
        'equal_blend': [mean(pair) for pair in zip(*small)],
        'nearest_linear_LAD': min(trend_options, key=lambda v: mae(small[1], v)),
    }
    return {
        'energy_MJ_m3': energy, 'effective_power_MW_m3': 1 / scale,
        'seconds_0p5L_1L_2L': [*small, target],
        'predicted_seconds': predictions,
        'MAE_seconds': {k: mae(target, v) for k, v in predictions.items()},
        'residual_seconds_observed_minus_predicted': [a - p for a, p in zip(target, cp)],
        'last_to_first_duration_ratio': [row[-1] / row[0] for row in [*small, target, cp]],
        'relative_capacity_correction_last_to_first': (target[-1] / energy[-1]) / (target[0] / energy[0]),
    }


def calculate():
    knots, rho0, fit = properties()
    sheet = dict(worksheets(SOURCE / 'liter-data.xlsx'))['Supplementary Information']
    coefficients = [float(fit[key]) for key in ('b0_kg_m3', 'b1_kg_m3_K', 'b2_kg_m3_K2')]

    def variable_density(t):
        if not float(fit['min_C']) <= t <= float(fit['max_C']):
            raise ValueError(f'Temperature outside density support: {t}')
        a, b, c = coefficients
        return a + b*t + c*t*t

    results = {}
    for lo, name, rho in [
        (-140, 'full_constant_density', lambda _: rho0),
        (-100, 'paper_span_constant_density', lambda _: rho0),
        (-80, 'overlap_constant_density', lambda _: rho0),
        (-80, 'overlap_variable_density', variable_density),
    ]:
        edges = list(range(lo, 1, 20))
        half, one, two = [duration(sheet, col, edges) for col in ('DR', 'DS', 'DT')]
        energy = [integral(knots, rho, a, b) for a, b in pairwise(edges)]
        results[name] = {'edges_C': edges, **score([half, one], two, energy)}
    return results


if __name__ == '__main__':
    (ROOT / 'results.json').write_text(json.dumps(calculate(), indent=2) + '\n')
