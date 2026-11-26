#!/usr/bin/env python3
"""Independent LP/eigenvalue verification of the scout's conditional Thg1 bound.

This is a model analysis of published single-addition kinetics, not a measured
fidelity, a general polymerase result, or an experimentally validated prediction.
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def main():
    root = Path(__file__).resolve().parents[1]
    source = root / "polymerase/thg1_table1_verified.csv"
    rows = list(csv.DictReader(source.open()))
    nts = "ACGU"
    complement = dict(zip(nts, "UGCA", strict=True))
    efficiency = np.array([[next(float(r["eff_uM_inv_min_inv"]) for r in rows
                                 if r["template"] == template and r["ntp"] == incoming)
                            for incoming in nts] for template in nts])
    results = []
    for censored in [10.0, 100.0, 10000.0]:
        e = efficiency.copy()
        e[nts.index("C"), nts.index("G")] = censored
        # Incoming-nucleotide indexing makes the correct efficiency diagonal.
        a = np.array([e[nts.index(complement[n])] for n in nts])
        b = a / np.diag(a)[:, None]
        np.fill_diagonal(b, 0.0)
        eigenvalues, eigenvectors = np.linalg.eig(b)
        ix = int(np.argmax(eigenvalues.real))
        rho = float(eigenvalues[ix].real)
        vector = np.abs(eigenvectors[:, ix].real)
        vector /= vector.sum()
        direct_probabilities = 1 / (1 + (b @ vector) / vector)
        eigen_optimum = 1 / (1 + rho)

        def feasible(probability):
            constraints = probability * a.copy()
            constraints[np.diag_indices(4)] -= np.diag(a)
            constraints /= np.max(np.abs(constraints), axis=1)[:, None]
            return linprog(np.zeros(4), A_ub=constraints, b_ub=np.zeros(4),
                           A_eq=np.ones((1, 4)), b_eq=[1.0], bounds=[(1e-14, None)] * 4,
                           method="highs", options={"primal_feasibility_tolerance": 1e-10,
                                                    "dual_feasibility_tolerance": 1e-10})

        low, high = 0.0, 1.0
        for _ in range(48):
            middle = (low + high) / 2
            if feasible(middle).success:
                low = middle
            else:
                high = middle
        assert abs(low - eigen_optimum) < 5e-8
        assert np.max(np.abs(direct_probabilities - eigen_optimum)) < 1e-10
        assert not feasible(0.99).success
        results.append({"censored_C_template_G_efficiency_assumption": censored,
                        "eigenvalue_minimax_probability": eigen_optimum,
                        "independent_LP_minimax_probability": low,
                        "LP_eigenvalue_difference": low - eigen_optimum,
                        "99_percent_feasible_under_model": False})
    # No assumption about the finite censored entry is needed for this ceiling.
    # rho(B) >= rho(B_principal) for a nonnegative matrix. Remove the row/column
    # indexed by the incoming G, whose correct efficiency is lower-censored.
    a = np.array([efficiency[nts.index(complement[n])] for n in nts])
    b = a / np.diag(a)[:, None]
    np.fill_diagonal(b, 0.0)
    keep = [nts.index(n) for n in "ACU"]
    rho_lower_bound = float(max(np.linalg.eigvals(b[np.ix_(keep, keep)]).real))
    cycle_product = (efficiency[nts.index("U"), nts.index("C")] / efficiency[nts.index("U"), nts.index("A")]
                     * efficiency[nts.index("G"), nts.index("A")] / efficiency[nts.index("G"), nts.index("C")])
    # Deterministic parameter rectangles, not confidence intervals. The source
    # reports +/- quantities; their statistical interpretation is not established
    # by the table, so do not relabel these envelopes as one/two-sigma intervals.
    row_data = {(r["template"], r["ntp"]): r for r in rows}
    rectangles = []
    for multiplier in [0, 1, 2]:
        lower = np.zeros((3, 3))
        for i, correct in enumerate("ACU"):
            template = complement[correct]
            r = row_data[template, correct]
            denominator = float(r["KD_uM"]) - multiplier * float(r["KD_reported_pm"])
            if denominator <= 0:
                raise ValueError("Unbounded correct efficiency in sensitivity envelope")
            correct_upper = (float(r["kpol_per_min"]) + multiplier * float(r["kpol_reported_pm"])) / denominator
            for j, incoming in enumerate("ACU"):
                if i == j:
                    continue
                r = row_data[template, incoming]
                mismatch_lower = max(0.0, float(r["kpol_per_min"]) - multiplier * float(r["kpol_reported_pm"])) / (float(r["KD_uM"]) + multiplier * float(r["KD_reported_pm"]))
                lower[i, j] = mismatch_lower / correct_upper
        rho_min = float(max(np.linalg.eigvals(lower).real))
        rectangles.append({"multiplier_of_reported_plus_minus_magnitudes": multiplier,
                           "minimum_offdiagonal_matrix": lower.tolist(),
                           "upper_bound_uniform_correct_probability": 1 / (1 + rho_min),
                           "interpretation": "Deterministic favorable-parameter envelope, not a statistical interval; efficiencies recomputed from kpol and KD."})
    report = {"source_doi": "10.1261/rna.078686.121", "table": "Table 1",
              "table_image_visually_checked": "sources/thg1-table1.jpg",
              "transcribed_data_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
              "independent_methods": "numpy eigensystem and scipy HiGHS LP feasibility bisection",
              "cases": results,
              "principal_submatrix_probability_upper_bound_for_any_censored_correct_G_rate": 1 / (1 + rho_lower_bound),
              "two_template_probability_upper_bound": 1 / (1 + cycle_product ** 0.5),
              "two_template_required_ratio_product_for_99_percent": 99 ** 2 * cycle_product,
              "favorable_parameter_rectangle_sensitivity": rectangles,
              "limitations": ["Dilute mutually exclusive substrate capture proportional to the published kpol/KD values is assumed.",
                              "The experiment measured single additions on preactivated tRNA substrates, not processive arbitrary-template synthesis.",
                              "Rounded point estimates are used; the bound is not a statistical confidence limit.",
                              "This does not apply automatically to Thg1-like proteins, alternative conditions, proofreading, or engineered variants."]}
    out = root / "polymerase/results"
    out.mkdir(exist_ok=True)
    (out / "thg1_independent_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
