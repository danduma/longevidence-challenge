"""Append-only counting assessments; computational receipts remain immutable."""
import datetime as dt
import json
import re

from protocol import ValidationError, require_text, question_key


def validate_investigation(value):
    require_text(value, ["id", "question", "relation", "rationale"])
    if not re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_.-]{0,95}", value["id"]):
        raise ValidationError("Invalid investigation id")
    if value["relation"] not in {"independent", "validation", "supersedes", "excluded"}:
        raise ValidationError("Invalid investigation relation")
    return value


LATEST = """SELECT r.* FROM investigation_reviews r JOIN
 (SELECT task_id,max(seq) AS seq FROM investigation_reviews GROUP BY task_id) latest
 ON latest.seq=r.seq"""


class InvestigationStore:
    def init_investigations(self):
        self.db.executescript("""
          CREATE TABLE IF NOT EXISTS investigations (
            id TEXT PRIMARY KEY, question TEXT NOT NULL, created TEXT NOT NULL);
          CREATE TABLE IF NOT EXISTS investigation_reviews (
            seq INTEGER PRIMARY KEY AUTOINCREMENT,
            task_id TEXT NOT NULL REFERENCES contributions(task_id),
            investigation_id TEXT NOT NULL REFERENCES investigations(id),
            relation TEXT NOT NULL, rationale TEXT NOT NULL, assessment TEXT NOT NULL,
            created TEXT NOT NULL);
        """)

    def _assess_investigation(self, task_id, assessment):
        validate_investigation(assessment)
        if not self.db.execute("SELECT 1 FROM contributions WHERE task_id=?", (task_id,)).fetchone():
            raise ValidationError("Investigation assessment requires an accepted receipt")
        existing = self.db.execute("SELECT question FROM investigations WHERE id=?", (assessment["id"],)).fetchone()
        if not existing:
            if any(question_key(r[0]) == question_key(assessment["question"])
                   for r in self.db.execute("SELECT question FROM investigations")):
                raise ValidationError("Reuse the existing id for this investigation question")
        if not existing and assessment["relation"] not in {"independent", "excluded"}:
            raise ValidationError("Validation must refer to an existing investigation")
        timestamp = dt.datetime.now(dt.timezone.utc).isoformat()
        self.db.execute("INSERT OR IGNORE INTO investigations VALUES (?,?,?)",
                        (assessment["id"], assessment["question"], timestamp))
        self.db.execute("""INSERT INTO investigation_reviews
          (task_id,investigation_id,relation,rationale,assessment,created) VALUES (?,?,?,?,?,?)""",
                        (task_id, assessment["id"], assessment["relation"], assessment["rationale"],
                         json.dumps(assessment, sort_keys=True), timestamp))
        self.event("investigation.assessed", task_id, investigation_id=assessment["id"],
                   relation=assessment["relation"], reason=assessment["rationale"], human_review=False)

    def assess_investigation(self, task_id, assessment):
        with self.db:
            self._assess_investigation(task_id, assessment)

    def investigation_counts(self):
        rows = list(self.db.execute(LATEST))
        active = {r["investigation_id"] for r in rows if r["relation"] != "excluded"}
        receipts = self.db.execute("SELECT count(*) FROM contributions").fetchone()[0]
        return dict(verified_contributions=len(active), verified_records=receipts,
                    unassessed_records=receipts-len(rows))

    def investigation_summaries(self):
        latest = list(self.db.execute(LATEST))
        created = dict(self.db.execute("SELECT task_id,created FROM contributions"))
        summaries = []
        for row in self.db.execute("SELECT * FROM investigations ORDER BY created,id"):
            members = [dict(task_id=r["task_id"], relation=r["relation"], rationale=r["rationale"])
                       for r in latest if r["investigation_id"] == row["id"]]
            active = [m for m in members if m["relation"] != "excluded"]
            summaries.append(dict(id=row["id"], question=row["question"], counted=bool(active),
                                  latest_receipt=max(active, key=lambda m: created[m["task_id"]])["task_id"] if active else None,
                                  members=members, human_review=False))
        return summaries
