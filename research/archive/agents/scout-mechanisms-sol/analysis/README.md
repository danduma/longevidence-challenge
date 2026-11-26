# Reproduce the exploratory analyses

Run from this job directory with Python 3. No third-party dependencies, network access, GPU, or wet laboratory work are required for these commands:

```sh
python3 analysis/cryo_source_audit.py
python3 analysis/quadruplet_audit.py
python3 analysis/thg1_competition.py
```

Outputs are deterministic except for explicit, fixed-seed heuristic/uncertainty sampling; all seeds are 20261126. JSON outputs contain numerical assumptions and warning fields. The cryopreservation script reads the downloaded original XLSX with Python's ZIP/XML libraries. The code audit exhaustively checks all 3,072 directed single-base substitutions in the 256-word space, filtering source words according to each code's sense set, and all ordered sense-codon pairs at frame offsets 1–3. The Thg1 analysis verifies its Perron eigenvalue through componentwise Collatz–Wielandt bounds and confirms equality of its optimized template probabilities.

These checks prove properties of the calculations and code maps under their stated assumptions. They do not validate living cells, polymerase processivity, tissue survival, or experimental challenge completion. The annealed quadruplet map is a heuristic example, not a certified optimum. The Thg1 parameter perturbation range is not a confidence interval.
