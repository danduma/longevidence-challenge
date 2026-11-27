"""Signed matched contrasts and finite numerical sensitivity ranges."""

REFERENCE = 'uniform_released_policy'
PAIR = ('graded', 'uniform_strict')
REFINED = ('coarse', 'medium', 'fine')
METRICS = ('objective_mpa', 'peak_mpa')
OVERLAP = 'overlap_affine'


def contrasts(runs):
    output = {}
    for name, run in runs.items():
        schedules = run['schedules']
        output[name] = {}
        for mapping in schedules[REFERENCE]:
            output[name][mapping] = {}
            for metric in METRICS:
                value = {s: results[mapping][metric] for s, results in schedules.items()}
                output[name][mapping][metric] = {
                    'graded_minus_strict': value[PAIR[0]]-value[PAIR[1]],
                    'graded_minus_reference': value[PAIR[0]]-value[REFERENCE],
                    'strict_minus_reference': value[PAIR[1]]-value[REFERENCE],
                }
    return output


def sensitivity(runs, comparison):
    output = {}
    for metric in METRICS:
        key = 'graded_minus_strict'
        fixed = {grid: {mapping: values[metric][key]
                        for mapping, values in comparison[grid].items()}
                 for grid in REFINED}
        linear = {grid: values[OVERLAP] for grid, values in fixed.items()}
        fine = linear['fine']
        differences = {
            'affine_coarse_to_medium': abs(linear['medium']-linear['coarse']),
            'affine_medium_to_fine': abs(fine-linear['medium']),
            'affine_three_grid_range': max(linear.values())-min(linear.values()),
            'affine_quadrature_change': abs(fine-comparison['fine_quadrature'][OVERLAP][metric][key]),
            'fine_transfer_range': max(fixed['fine'].values())-min(fixed['fine'].values()),
        }
        shifts = {
            grid: {schedule: rows[OVERLAP][metric]-rows['idw8'][metric]
                   for schedule, rows in runs[grid]['schedules'].items()}
            for grid in REFINED
        }
        output[metric] = dict(affine_by_grid=linear, fine_abs_contrast=abs(fine),
                              sensitivity=differences, affine_minus_idw8_scores=shifts)
    return output


def record(runs, events):
    comparison = contrasts(runs)
    return dict(events=events, runs=runs, contrasts=comparison,
                sensitivity=sensitivity(runs,comparison))


if __name__ == '__main__':
    import json
    from pathlib import Path
    from analysis import SETTINGS
    runs = {name: json.loads(Path(f'{name}.json').read_text()) for name in SETTINGS}
    events = json.loads(Path('events.json').read_text())
    Path('results.json').write_text(json.dumps(record(runs,events),indent=2,allow_nan=False)+'\n')
