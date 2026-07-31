"""Reconcile Python and TwoSampleMR single-lead Wald estimates."""

from __future__ import annotations

import csv
import json
from pathlib import Path


BASE = Path(r"E:/20241004_MCI/20260709")
PY = BASE / "mci_frozen_gca_lead_wald_20260725_142608" / "lead_wald_summary.csv"
R = BASE / "mci_frozen_gca_lead_wald_twosamplemr_20260725_142801" / "twosamplemr_lead_wald_summary.csv"
OUT = BASE / "mci_frozen_gca_lead_wald_twosamplemr_20260725_142801" / "reconciliation.json"


def rows(path, r_file=False):
    with path.open(newline="", encoding="utf-8") as fh:
        data = list(csv.DictReader(fh))
    if r_file:
        return {row["id.outcome"]: row for row in data}
    return {row["outcome"]: row for row in data}


p, r = rows(PY), rows(R, r_file=True)
if set(p) != set(r):
    raise SystemExit("Outcome sets differ")
fields = {"beta": "b", "se": "se", "pval": "pval", "OR": "OR", "OR_lower95": "OR_lci95", "OR_upper95": "OR_uci95"}
diffs = []
for outcome in p:
    rec = {"outcome": outcome}
    for pk, rk in fields.items():
        a, b = float(p[outcome][pk]), float(r[outcome][rk])
        rec[f"abs_diff_{pk}"] = abs(a - b)
    diffs.append(rec)
summary = {
    "python_summary": str(PY),
    "twosamplemr_summary": str(R),
    "outcome_count": len(diffs),
    "all_reconciled": all(max(v for k, v in d.items() if k != "outcome") < 1e-10 for d in diffs),
    "max_absolute_difference": max(v for d in diffs for k, v in d.items() if k != "outcome"),
    "lead_snp": "rs918928",
    "multi_snp_ivw_status": "rejected_due_to_strong_ld",
}
OUT.write_text(json.dumps({"summary": summary, "by_outcome": diffs}, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
