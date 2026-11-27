"""Print paired differences, keeping phase and stress metric explicit."""
import json
from pathlib import Path


def paired_rows(rows):
    for row in rows:
        grade = row["designs"]["graded"]
        uniform = row["designs"]["uniform_released_policy"]
        strict = row["designs"]["uniform_strict"]
        yield (row["spacing_mm"], row["phase"],
               uniform["objective_mpa"] - grade["objective_mpa"],
               uniform["peak_mpa"] - grade["peak_mpa"],
               strict["objective_mpa"] - grade["objective_mpa"],
               strict["peak_mpa"] - grade["peak_mpa"])


if __name__ == "__main__":
    path = Path(__file__).with_name("results.json")
    for spacing, phase, *differences in paired_rows(json.loads(path.read_text())):
        print(spacing, phase, *(f"{difference:+.7f}" for difference in differences))
