"""Research contracts. No provider SDK or orchestrator imports."""
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urlparse


class ValidationError(ValueError):
    pass


def require_text(record, fields):
    if not isinstance(record, dict):
        raise ValidationError("Expected an object")
    for field in fields:
        if not isinstance(record.get(field), str) or not record[field].strip():
            raise ValidationError(f"Missing nonempty text: {field}")


def question_key(question):
    return hashlib.sha256(" ".join(question.casefold().split()).encode()).hexdigest()


def validate_task(task):
    require_text(task, ["id", "title", "question", "increment", "why_useful"])
    if not re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_.-]{0,95}", task["id"]):
        raise ValidationError("Task id must be a short filesystem-safe identifier")
    for field in ["context", "depends_on", "watches"]:
        if not isinstance(task.get(field), list) or any(not isinstance(v, str) for v in task[field]):
            raise ValidationError(f"{field} must be a list of strings")
    if not task["context"] or task["id"] in task["depends_on"]:
        raise ValidationError("Need source context and no self dependency")
    return task


def confined(root, relative):
    if not isinstance(relative, str) or Path(relative).is_absolute():
        raise ValidationError("Evidence paths must be relative")
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()) or path.is_symlink():
        raise ValidationError("Evidence path escapes job directory")
    return path


def evidence_manifest(root, paths):
    result = {}
    for relative in paths:
        path = confined(root, relative)
        if path.is_dir():
            children = [p for p in sorted(path.rglob("*")) if p.is_file()]
            if not children:
                raise ValidationError(f"Empty evidence directory: {relative}")
            child_paths = [p.relative_to(root.resolve()).as_posix() for p in children]
            result.update(evidence_manifest(root, child_paths))
            continue
        if not path.is_file():
            raise ValidationError(f"Missing evidence file: {relative}")
        result[path.relative_to(root.resolve()).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def validate_sources(job, sources):
    if not isinstance(sources, list) or not sources:
        raise ValidationError("Need source-backed evidence")
    for source in sources:
        require_text(source, ["url", "path", "read_scope", "role"])
        if urlparse(source["url"]).scheme not in ("https", "http"):
            raise ValidationError("Sources need original HTTP(S) URLs")
        path = confined(job, source["path"])
        if not path.is_file() or path.stat().st_size == 0:
            raise ValidationError("Source snapshot missing or empty")


def validate_gap(gap):
    if not isinstance(gap, dict):
        raise ValidationError("Gap must be an object")
    require_text(gap, ["question", "existing_evidence", "missing_data", "why_insufficient",
                       "supports", "refutes", "inconclusive", "minimum_design",
                       "experimental_unit", "controls", "precision_basis",
                       "work_possible_now", "reopen_trigger"])
    if not isinstance(gap.get("candidate_producers"), list):
        raise ValidationError("candidate_producers must be a list (empty if no supported match)")
    for producer in gap["candidate_producers"]:
        require_text(producer, ["name", "capability", "strategic_fit", "source_url", "confidence"])
        if not producer["source_url"].startswith("https://"):
            raise ValidationError("Producer fit requires a supporting source")
    return gap


def validate_result(job, result):
    if not isinstance(result, dict) or result.get("kind") not in ("contribution", "data_gap", "not_useful"):
        raise ValidationError("kind must be contribution, data_gap or not_useful")
    require_text(result, ["title", "summary", "prior_art", "limitations"])
    validate_sources(job, result.get("sources"))
    if not isinstance(result.get("artifacts"), list) or not result["artifacts"]:
        raise ValidationError("Need artifact files")
    evidence_manifest(job, result["artifacts"])
    for followup in result.get("followups", []):
        validate_task(followup)
    if result["kind"] == "contribution":
        require_text(result, ["claim", "increment", "decision_changed", "verification_scope"])
        if not isinstance(result.get("verification"), list) or not result["verification"]:
            raise ValidationError("A contribution needs fresh executable verification")
    if result["kind"] == "data_gap":
        validate_gap(result.get("gap"))
    return result


REVIEW_CHECKS = ("distinct", "source_supported", "reproducible", "useful", "scope_honest", "prior_art_checked")


def validate_review(review, result):
    if not isinstance(review, dict) or review.get("verdict") not in ("accept", "revise", "reject"):
        raise ValidationError("Invalid review verdict")
    require_text(review, ["reason", "strongest_objection"])
    if not isinstance(review.get("approved_followups"), list):
        raise ValidationError("Need explicit approved_followups list")
    offered = {task["id"] for task in result.get("followups", [])}
    if not set(review["approved_followups"]).issubset(offered):
        raise ValidationError("Reviewer approved an unknown followup")
    if review["verdict"] == "accept":
        if result["kind"] == "contribution":
            from investigations import validate_investigation
            validate_investigation(review.get("investigation"))
            if any(review.get("checks", {}).get(key) is not True for key in REVIEW_CHECKS):
                raise ValidationError("Acceptance requires every review check")
        elif result["kind"] == "data_gap" and review.get("data_gate_valid") is not True:
            raise ValidationError("Data gate not justified")
    return review


RESULT_CONTRACT = """
Return only a JSON object, also saved as submission.json in your job directory:
{
 "kind": "contribution | data_gap | not_useful", "title": "...", "summary": "...",
 "prior_art": "What you actually checked, closest precedent, and surviving increment",
 "limitations": "...", "sources": [{"url":"original https URL", "path":"sources/local-snapshot",
 "read_scope":"full text / code / excerpt / indexed snippet; retrieval date and version",
 "role":"what this source supports"}],
 "artifacts": ["report.md", "analysis.py", "verify.py", "results.json"],
 "verification": [["python3", "verify.py"]],
 "claim": "precise result", "increment": "difference from prior work",
 "decision_changed": "concrete decision or interpretation", "verification_scope": "what checks establish",
 "followups": [{"id":"stable-distinct-id", "title":"...", "question":"...",
 "increment":"...", "why_useful":"...", "context":["source URL or read-only context path"],
 "depends_on":[], "watches":[]}]
}
For data_gap include gap: {question, existing_evidence, missing_data, why_insufficient,
supports, refutes, inconclusive, minimum_design, experimental_unit, controls,
precision_basis, work_possible_now, reopen_trigger, candidate_producers:[{name,
capability, strategic_fit, source_url, confidence}]} (all other gap fields are strings).
Gap cards and not_useful outcomes DO NOT count as verified contributions. Their
verification may be empty; omit contribution-only claim/increment/decision fields.
Candidate fit is an inference, not a partnership or proof data exist. No invented power/sample size.
Every source needs a real nonempty snapshot in your job directory; distinguish downloaded
code/data, complete text, excerpts and inaccessible pages. Preserve attribution/licenses.
"""


REVIEW_CONTRACT = """
Return only JSON: {"verdict":"accept | revise | reject", "reason":"...",
"strongest_objection":"...", "checks":{"distinct":true,"source_supported":true,
"reproducible":true,"useful":true,"scope_honest":true,"prior_art_checked":true},
"data_gate_valid":false,"approved_followups":["offered task ids"],
"investigation":{"id":"stable-investigation-id","question":"the substantive research question",
"relation":"independent | validation | supersedes | excluded","rationale":"why this is a new investigation or belongs to an existing one"}}.
Accept a contribution only if ALL six checks are established. A successful script alone
is not scientific support. Check the actual code, outputs, cached primary sources,
strongest feasible comparator and the ledger for duplicates. Reject splitting sweeps,
old results relabelled as new, generic literature summaries or trivial arithmetic.
Do not demand positive results or human acceptance for a well-supported useful computational
finding. State limitations and do not equate agent review with human review or biological validation.
For an accepted computational contribution an investigation assessment is mandatory.
The distinct check concerns a useful new computational increment; the investigation
assessment determines whether it is also independent for the target count. Reuse
the existing investigation id for controls, numerical convergence, sensitivity checks,
solver repairs, tuning/scoring changes and corrections to that investigation. These may
be valuable accepted receipts but add NO new investigation count. Consult the ledger's
investigations list. A new independent id requires a substantive separate question and
useful answer that can stand apart from validating the earlier claim. Different task names
or a new metric alone do not establish independence. If new evidence defeats an earlier
conclusion, describe what is superseded; preserve its receipt instead of erasing history.
For a data gate, check computations/search alternatives were considered and exact evidence,
discriminating design, controls and reopening trigger are documented. Gap cards do not count.
Followups must be meaningfully distinct, feasible, source-grounded and worth doing;
do not approve quota-filling or recursive restatements. You may approve useful followups
even if the present contribution is rejected. Do not modify evidence or submission files.
"""
