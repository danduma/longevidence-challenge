"""Build the coordinator record from the verified, fixed-schedule comparison."""
import json
from pathlib import Path

from summarize_results import paired_rows

ROOT = Path(__file__).resolve().parent
PIN = "1d498ee5c583e0112a762e60892d576ec2e27f57"
BASE = "https://github.com/Marc-Dvci/Vitrify/blob/" + PIN + "/"


def source(url, path, scope, role):
    assert (ROOT / path).stat().st_size > 0
    return {"url": url, "path": path, "read_scope": scope, "role": role}


def pinned(path, role, scope="Complete code"):
    return source(BASE + path, "sources/upstream/" + path,
                  f"{scope}, commit {PIN}, bytes checked 2026-11-27.", role)


def make():
    rows = json.loads((ROOT / "results.json").read_text())
    differences = list(paired_rows(rows))
    fine = [row for row in differences if row[0] == .75]
    coarse = [row for row in differences if row[0] == 1.1]

    def count_lower(group, index):
        return sum(row[index] < 0 for row in group)

    fine_policy = count_lower(fine, 2)
    coarse_policy = count_lower(coarse, 2)
    fine_peak = count_lower(fine, 3)
    coarse_peak = count_lower(coarse, 3)
    fine_strict = count_lower(fine, 4)
    coarse_strict = count_lower(coarse, 4)
    rates = [row["designs"]["uniform_strict"]["rate_70"]
             for row in rows if row["spacing_mm"] == 1.1]
    return {
        "kind": "contribution",
        "title": "Grid origin changes the strict rate verdict in Vitrify’s rat-heart phantom",
        "summary": (
            "A fixed half-cell shift can erase the strict uniform surround’s 1.1 mm rate margin. "
            "On the finest released grid, all three fixed schedules clear 51.5 °C/min in all five phases. "
            f"The release-policy uniform candidate has a lower sampled stress objective than the graded design "
            f"in {coarse_policy}/5 phases at 1.1 mm and {fine_policy}/5 at 0.75 mm. Peak-tension ordering reverses under some shifts."
        ),
        "prior_art": (
            "Checked the controlling reassessment and ledger, the earlier tuned-uniform comparison, pinned "
            "Vitrify geometry, convergence experiment, thermal and mechanical code, public issue/PR and margin "
            "history records, primary heart thermal and mechanical papers, the later artery mechanics paper, "
            "and specialist discussions at CryoNet and Biostasis. The released numerical precedent changes spacing "
            "at one origin per spacing. The previous coordinator result introduced the stronger tuned uniform "
            "baseline. The heart and artery studies used conventional mesh convergence, while the forum discussions "
            "address nonuniform heating and glass-transition stress. The surviving increment is a controlled "
            "within-spacing phase comparison of the same calibrated phantom and fixed protocols."
        ),
        "limitations": (
            "Five deterministic phases per spacing, one synthetic phantom, extrapolated SAR, and fixed schedules. "
            "Each spacing retains its own native shape calibration and loading centre. Mechanics are sampled "
            "Maxwell simulations with an unnormalized smooth maximum. The limiter occurs below glass transition. "
            "The results give neither an all-phase bound nor biological validation or an optimum over design classes."
        ),
        "sources": [
            pinned("src/vitrify/geometry.py", "Analytic phantom dimensions, material labels and native pillow fit."),
            pinned("results/grid_convergence_margin.json", "Original graded vector and native-grid thermal anchors.", "Complete JSON data"),
            pinned("experiments/grid_convergence.py", "Single-origin spacing sweep and fixed Case II comparator."),
            pinned("coupler/workflow.py", "Continuous radial loading and coupled objective conventions."),
            pinned("experiments/optimize_protocol.py", "Duty and loading bounds and rate policy."),
            pinned("src/vitrify/thermal.py", "Thermal stencil, rate and horizon."),
            pinned("src/vitrify/properties.py", "Material and extrapolated SAR laws."),
            pinned("tesseracts/mechanics/solver.py", "Q1 equilibrium, Maxwell history and smooth stress score."),
            pinned("LICENSE", "Apache-2.0 licensing of cached upstream files.", "Complete license"),
            source("https://api.github.com/repos/Marc-Dvci/Vitrify/commits?path=results/optimization_margin.json&per_page=20",
                   "sources/margin-history.json", "Complete first page of responses cached 2026-11-27.",
                   "Margin design history."),
            source("https://api.github.com/repos/Marc-Dvci/Vitrify/issues?state=all&per_page=100",
                   "sources/issues.json", "Complete first page of responses cached 2026-11-27.",
                   "Repository-specific prior-art check."),
            source("https://api.github.com/repos/Marc-Dvci/Vitrify/pulls?state=all&per_page=100",
                   "sources/pulls.json", "Complete first page of responses cached 2026-11-27.",
                   "Repository-specific prior-art check."),
            source("https://pmc.ncbi.nlm.nih.gov/articles/PMC8823202/?report=xml", "sources/thermal-excerpt.txt",
                   "2022 article, attributed excerpt and reading notes cached 2026-11-27.",
                   "Coldest-heart endpoint and rate interval interpretation."),
            source("https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0290063",
                   "sources/mechanics.xml", "Complete article XML cached 2026-11-27. Read thermal and numerical methods and relevant results.",
                   "Original heart stress and mesh convergence precedent."),
            source("https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0350078",
                   "sources/artery.xml", "Complete article XML cached 2026-11-27. Read numerical mesh and nanowarming sections.",
                   "Later fitted-mesh analysis and interface sensitivity precedent."),
            source("https://www.cryonet.org/cgi-bin/dsp.cgi?msg=18925", "sources/cryonet-excerpt.txt",
                   "Attributed excerpt and reading notes, 2002 message, cached 2026-11-27.",
                   "Specialist discussion of spatial heating."),
            source("https://www.biostasis.com/new-warming-breakthrough-for-cryopreserved-organs/",
                   "sources/biostasis-excerpt.txt", "Attributed excerpt and reading notes, 2017 article, cached 2026-11-27.",
                   "Specialist discussion of uniform heating around glass transition.")
        ],
        "artifacts": [
            "report.md", "analysis.py", "forward.py", "verify.py", "verify.sh", "results.json", "thermal.json",
            "protocol.json", "candidates.json", "requirements.txt", "attribution.md", "summarize_results.py",
            "sources/upstream-sha256.txt", "verification.log", "thermal-verification.log", "assemble_submission.py"
        ],
        "verification": [["sh", "verify.sh"]],
        "claim": (
            f"The strict witness passes 51.5 °C/min in 3/5 phases at 1.1 mm "
            f"({min(rates):.4f}–{max(rates):.4f} °C/min) and 5/5 at 0.75 mm. "
            "The selected limiter is always startup myocardium. Three 1.1 mm strict runs place it next to a wall voxel across an analytic 0.810–0.833 mm CPA gap. "
            f"At 1.1 and 0.75 mm the release-policy uniform score is lower than the released graded score "
            f"in {coarse_policy}/5 and {fine_policy}/5 phases for the smooth objective and "
            f"{coarse_peak}/5 and {fine_peak}/5 for peak tension. The strict uniform smooth score is lower in "
            f"{coarse_strict}/5 and {fine_strict}/5 phases, respectively. At 0.75 mm strict uniform has lower peak tension only at the native origin."
        ),
        "increment": (
            "Isolates sampling origin within each spacing by freezing the continuous shape and the graded loading "
            "field, while preserving the exact earlier uniform and graded loading and duty vectors. Reports paired "
            "mechanical ordering together with the rate decision and limiting voxel context."
        ),
        "decision_changed": (
            "Reject the 1.1 mm native strict pass as evidence of alignment-robust feasibility. Use the sampled "
            "finest-grid phase results and retain the release-policy uniform as the smooth-objective benchmark. Treat the peak-tension ordering as phase-sensitive."
        ),
        "verification_scope": (
            "Reruns all ten voxelizations, thirty finite-volume trajectories and thirty Maxwell stress trajectories "
            "from cached inputs. Asserts source native thermal anchors, the prior finest-grid stress anchor, "
            "pass/fail, limiter timing and adjacency, analytic gap, paired orderings, and numerical agreement of every cached score."
        ),
        "followups": [{
            "id": "cryo-heart-wall-contact-continuous-gap",
            "title": "Resolve the apparent myocardium–wall contacts in the rate stencil",
            "question": "Does resolving the CPA gap remove the low startup rates associated with the coarse wall-adjacent stencil?",
            "increment": (
                "Use the measured continuous axial gaps at the coarse contacts. Compare the same schedules and "
                "phantom with a conservative interface-resolved thermal discretization that preserves real contacts "
                "and resolves nonzero gaps. Use the native and failing phases as matched controls."
            ),
            "why_useful": (
                "Tests how much of the rate depression comes from numerically skipping a "
                "real thin gap. The analytic predicates and thermal coefficients permit the comparison now."
            ),
            "context": [BASE + "src/vitrify/geometry.py", BASE + "src/vitrify/thermal.py"],
            "depends_on": ["cryo-margin-grid-origin-sensitivity"],
            "watches": []
        }]
    }


if __name__ == "__main__":
    (ROOT / "submission.json").write_text(json.dumps(make(), indent=2, ensure_ascii=False) + "\n")
