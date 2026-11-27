"""Independently reconstruct the decisive comparisons from cached sources."""
import csv
from fractions import Fraction as F
from hashlib import sha256
from itertools import pairwise
import json
from math import isclose
from pathlib import Path
import re
from statistics import median
from zipfile import ZipFile

from workbook import worksheets

ROOT = Path(__file__).resolve().parent
SRC = ROOT / 'sources'


def close(expected, actual):
    assert isclose(float(expected), actual, rel_tol=1e-10, abs_tol=1e-10), (expected, actual)


def source_curve():
    with ZipFile(SRC / 'liter-supp.docx') as archive:
        document = archive.read('word/document.xml').decode()
    plain = ''.join(re.findall(r'<w:t(?: [^>]*)?>(.*?)</w:t>', document))
    assert 'Specific Heat: Cp (KJ/kg.K)' in plain
    block = plain.rsplit('Table S4: Thermo-physical properties', 1)[1].split('Table S5: Capabilities')[0].split('M220.3[25]')[1]
    values, density = block.split('[26]', 1)
    knots = sorted((F(t), F(cp) * 1000)
                   for cp, t in re.findall(r'([\d.]+)\s*\((-?\d+)⁰C\)', values))
    assert len(knots) == 7
    return knots, F(density.split('[27]')[0])


def crossings(sheet, column, edges):
    times = {}
    rows = sorted(int(key[2:]) for key in sheet if re.fullmatch(column + r'\d+', key)
                  and int(key[2:]) >= 6)
    for i, j in pairwise(rows):
        t, u = (F(sheet[f'DQ{r}']) for r in (i, j))
        assert t < u
        x, y = (F(sheet[f'{column}{r}']) for r in (i, j))
        for edge in edges:
            if edge not in times and x < edge <= y:
                times[edge] = (t * (y - edge) + u * (edge - x)) / (y - x)
    return [times[b] - times[a] for a, b in pairwise(edges)]


def simpson(knots, rho, start, end):
    result = F(0)
    for (left, cp_left), (right, cp_right) in pairwise(knots):
        a, b = max(start, left), min(end, right)
        if b <= a:
            continue
        def product(t):
            specific = (cp_left * (right - t) + cp_right * (t - left)) / (right - left)
            return specific * rho(t)
        result += (b - a) * (product(a) + 4 * product((a + b) / 2) + product(b)) / 6
    return result / 1000000


def mae(actual, prediction):
    return sum(abs(x - y) for x, y in zip(actual, prediction, strict=True)) / len(actual)


def check():
    input_hashes = {
        'liter-data.xlsx': 'c059e1beef528efd5a3aa04fc51d039d64e1c8583376a3596160e837069d2898',
        'liter-supp.docx': '63336d3e3e002e1bfda535e3609aa2758e89a3c360176a79e28b324188ee4059',
        'solanki-table2.csv': 'd157e4ff7a6f1044f9f235eb6aa38dcb0f3e266b390a32918bfd338108905d75',
    }
    for name, expected in input_hashes.items():
        assert sha256((SRC / name).read_bytes()).hexdigest() == expected, name
    knots, constant = source_curve()
    data = dict(worksheets(SRC / 'liter-data.xlsx'))['Supplementary Information']
    assert '4.6' in data['DR5'] and '4.6' in data['DS5'] and '4.6' in data['DT5']
    with (SRC / 'solanki-table2.csv').open() as infile:
        row = next(csv.DictReader(infile))
    b0, b1, b2 = (F(row[key]) for key in ('b0_kg_m3', 'b1_kg_m3_K', 'b2_kg_m3_K2'))
    low, high = F(row['min_C']), F(row['max_C'])

    def density(t):
        assert low <= t <= high
        return (b2*t + b1)*t + b0

    saved = json.loads((ROOT / 'results.json').read_text())
    for name, start, rho in (
        ('full_constant_density', -140, lambda _: constant),
        ('paper_span_constant_density', -100, lambda _: constant),
        ('overlap_constant_density', -80, lambda _: constant),
        ('overlap_variable_density', -80, density),
    ):
        record = saved[name]
        edges = list(range(start, 1, 20))
        assert record['edges_C'] == edges and knots[0][0] <= start and 0 <= knots[-1][0]
        half, one, two = [crossings(data, col, edges) for col in ('DR', 'DS', 'DT')]
        intervals = list(pairwise(edges))
        energies = [simpson(knots, rho, F(a), F(b)) for a, b in intervals]
        power = 2 * sum(energies) / (sum(half) + sum(one))
        cp_times = [energy / power for energy in energies]
        best_trend = None
        # For each slope candidate the median residual is the optimal intercept.
        for i, a in enumerate(one):
            for j, b in enumerate(one[i+1:], i+1):
                slope = (b-a)/(j-i)
                offset = median(y-slope*k for k, y in enumerate(one))
                trial = [offset+slope*k for k in range(len(one))]
                rank = mae(one, trial), slope, offset, trial
                best_trend = rank if best_trend is None else min(best_trend, rank)
        models = {
            'heat_capacity': cp_times, 'nearest_1L': one,
            'equal_blend': [(a+b)/2 for a,b in zip(half,one)],
            'nearest_linear_LAD': best_trend[-1],
            'pooled_mean_constant': [sum(half+one)/(2*len(one))]*len(one),
            'pooled_median_constant': [median(half+one)]*len(one),
            'nearest_mean_constant': [sum(one)/len(one)]*len(one),
            'nearest_median_constant': [median(one)]*len(one),
        }
        close(power, record['effective_power_MW_m3'])
        for a, b in zip(energies, record['energy_MJ_m3'], strict=True):
            close(a, b)
        for i, times in enumerate((half, one, two, cp_times)):
            close(times[-1]/times[0], record['last_to_first_duration_ratio'][i])
            if i < 3:
                for a, b in zip(times, record['seconds_0p5L_1L_2L'][i], strict=True):
                    close(a, b)
        for key, prediction in models.items():
            close(mae(two, prediction), record['MAE_seconds'][key])
            for a, b in zip(prediction, record['predicted_seconds'][key], strict=True):
                close(a, b)
        for a, b, c in zip(two, cp_times, record['residual_seconds_observed_minus_predicted'], strict=True):
            close(a-b, c)
        correction = (two[-1]/two[0]) / (energies[-1]/energies[0])
        close(correction, record['relative_capacity_correction_last_to_first'])
        assert mae(two, one) < mae(two, cp_times)
        if start == -140:
            assert mae(two, cp_times) < mae(two, models['nearest_median_constant'])
        if name == 'overlap_variable_density':
            assert correction > F('1.45') and mae(two, cp_times) > mae(two, models['nearest_median_constant'])
            assert two[0]-cp_times[0] < -3 and two[-1]-cp_times[-1] > 1
    print('Verified published knots, crossings, exact integrals, smaller-bag fits and held-out residuals.')


if __name__ == '__main__':
    check()
