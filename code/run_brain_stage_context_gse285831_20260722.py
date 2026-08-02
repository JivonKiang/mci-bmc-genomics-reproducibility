from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_stage_discovery_gse63060_20260721 import fit_linear  # noqa: E402


ROOT = Path(os.environ.get("MCI_ROOT", Path(__file__).resolve().parents[2]))
BASE = ROOT / "20260709"
SC_DIR = ROOT / "20250816 revise" / "scRNA"
TAR_PATH = SC_DIR / "GSE285831_RAW.tar"
SERIES_PATH = SC_DIR / "GSE285831_series_matrix.txt.gz"
LOOKUP_PATH = BASE / "progression_concordance_20260721_1600" / "ensembl_symbol_lookup.json"
DISC_PATH = BASE / "stage_discovery_20260721_1510" / "stage_discovery_gene_effects.csv"
FREEZE_PATH = BASE / "stage_discovery_20260721_1510" / "stage_candidate_registry_monotonic_freeze_candidate.csv"
OUT = BASE / "brain_stage_context_gse285831_20260722"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_series_labels(path: Path) -> pd.DataFrame:
    rows: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession") or line.startswith("!Sample_title"):
                fields = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                rows[fields[0]] = [x.strip() for x in fields[1:]]
    sample_ids = rows["!Sample_geo_accession"]
    titles = rows["!Sample_title"]
    if len(sample_ids) != len(titles):
        raise ValueError("GSE285831 series metadata accession/title length mismatch")
    records = []
    for sid, title in zip(sample_ids, titles):
        upper = title.upper()
        if ", CTR," in upper:
            group = "HC"
        elif ", MCI," in upper:
            group = "MCI"
        elif ", AD," in upper:
            group = "AD"
        else:
            raise ValueError(f"Unrecognized GSE285831 title: {title}")
        records.append(
            {
                "sample_id": sid.strip('"'),
                "title": title.strip('"'),
                "group": group,
                "tissue": "frontal cortex (BA9)",
                "cell_type": "layer III pyramidal neurons",
            }
        )
    return pd.DataFrame(records).set_index("sample_id")


def read_counts(path: Path, sample_ids: list[str]) -> pd.DataFrame:
    cols: dict[str, pd.Series] = {}
    with tarfile.open(path, "r") as tar:
        members = [m for m in tar.getmembers() if m.name.endswith(".genes.results.gz")]
        if len(members) != len(sample_ids):
            raise ValueError(f"Expected {len(sample_ids)} count files, found {len(members)}")
        for member in members:
            m = re.match(r"(GSM\d+)_", Path(member.name).name)
            if not m:
                raise ValueError(f"Cannot parse sample id from {member.name}")
            sid = m.group(1)
            if sid not in sample_ids:
                raise ValueError(f"Count file {sid} is not in official sample metadata")
            raw = tar.extractfile(member)
            if raw is None:
                raise ValueError(f"Cannot extract {member.name}")
            with gzip.GzipFile(fileobj=raw) as gz:
                tab = pd.read_csv(gz, sep="\t", usecols=["gene_id", "expected_count"])
            tab["gene_id"] = tab["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
            tab["expected_count"] = pd.to_numeric(tab["expected_count"], errors="coerce").fillna(0.0)
            cols[sid] = tab.groupby("gene_id", sort=False)["expected_count"].sum()
    counts = pd.DataFrame(cols).fillna(0.0)
    return counts.loc[:, sample_ids]


def main() -> None:
    labels = parse_series_labels(SERIES_PATH)
    sample_ids = labels.index.tolist()
    counts = read_counts(TAR_PATH, sample_ids)
    labels = labels.loc[counts.columns]
    labels["stage_num"] = labels["group"].map({"HC": 0, "MCI": 1, "AD": 2})
    lib = counts.sum(axis=0).to_numpy(float)
    log_cpm = np.log2(counts.to_numpy(float) / lib[None, :] * 1e6 + 1.0)

    X_cat = np.column_stack(
        [
            np.ones(len(labels)),
            labels["group"].eq("MCI").to_numpy(float),
            labels["group"].eq("AD").to_numpy(float),
        ]
    )
    X_trend = np.column_stack([np.ones(len(labels)), labels["stage_num"].to_numpy(float)])
    contrasts_cat = {
        "MCI_vs_HC": np.array([0, 1, 0]),
        "AD_vs_MCI": np.array([0, -1, 1]),
        "AD_vs_HC": np.array([0, 0, 1]),
    }
    contrasts_trend = {"ordinal_HC_MCI_AD": np.array([0, 1])}
    effects = pd.concat(
        [
            fit_linear(log_cpm.T, X_cat, contrasts_cat, counts.index.astype(str).tolist()),
            fit_linear(log_cpm.T, X_trend, contrasts_trend, counts.index.astype(str).tolist()),
        ],
        ignore_index=True,
    )

    lookup = json.loads(LOOKUP_PATH.read_text(encoding="utf-8"))
    freeze = pd.read_csv(FREEZE_PATH)
    disc = pd.read_csv(DISC_PATH)
    trend = disc.loc[disc["contrast"].eq("ordinal_HC_MCI_AD"), ["gene", "estimate"]].rename(
        columns={"estimate": "discovery_trend_estimate"}
    )
    cand = freeze[["gene"]].merge(trend, on="gene", how="left")
    cand["ensembl_id"] = cand["gene"].map(lookup)
    cand = cand.dropna(subset=["ensembl_id"])
    cand = cand[cand["ensembl_id"].isin(counts.index)].copy()
    cand.to_csv(OUT / "frozen_stage_candidates_mapped_to_GSE285831.csv", index=False)
    if cand.empty:
        raise ValueError("No frozen candidates mapped to GSE285831")

    z_rows = []
    for ensembl_id in cand["ensembl_id"]:
        v = log_cpm[counts.index.get_loc(ensembl_id), :]
        sd = np.std(v, ddof=1)
        z_rows.append((v - np.mean(v)) / sd if sd > 0 else np.zeros_like(v))
    weights = cand["discovery_trend_estimate"].to_numpy(float)
    weights = weights / np.sum(np.abs(weights))
    labels["stage_program_score"] = weights @ np.vstack(z_rows)
    labels.reset_index().to_csv(OUT / "GSE285831_sample_stage_program_scores.csv", index=False)

    score_effects = pd.concat(
        [
            fit_linear(labels[["stage_program_score"]].to_numpy(float), X_cat, contrasts_cat, ["frozen_stage_program_score"]),
            fit_linear(labels[["stage_program_score"]].to_numpy(float), X_trend, contrasts_trend, ["frozen_stage_program_score"]),
        ],
        ignore_index=True,
    )
    score_effects.to_csv(OUT / "frozen_stage_program_brain_context_effects.csv", index=False)

    candidate_effects = effects.loc[effects["feature_id"].isin(cand["ensembl_id"])].merge(
        cand[["gene", "ensembl_id", "discovery_trend_estimate"]], left_on="feature_id", right_on="ensembl_id", how="left"
    )
    candidate_effects.to_csv(OUT / "frozen_stage_candidate_brain_context_effects.csv", index=False)
    labels.reset_index().to_csv(OUT / "GSE285831_sample_audit.csv", index=False)

    audit = {
        "dataset": "GSE285831",
        "source_tar": str(TAR_PATH),
        "source_tar_sha256": sha256(TAR_PATH),
        "series_matrix": str(SERIES_PATH),
        "n_samples": int(len(labels)),
        "sample_group_counts": labels["group"].value_counts().to_dict(),
        "n_gene_features": int(counts.shape[0]),
        "n_frozen_candidates": int(len(freeze)),
        "n_frozen_candidates_mapped": int(len(cand)),
        "tissue": "frontal cortex (BA9)",
        "cell_type": "layer III pyramidal neurons",
        "analysis_role": "cross-tissue stage-associated context; not single-cell and not progression validation",
        "model": "unadjusted group contrasts and ordinal trend; GEO metadata did not provide audited age/sex covariates",
        "old_single_cell_or_KO_claims_used": False,
    }
    (OUT / "brain_stage_context_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    pd.DataFrame([audit]).to_csv(OUT / "brain_stage_context_audit.csv", index=False)
    (OUT / "analysis_registry.md").write_text(
        "# GSE285831 brain stage-context registry\n\n"
        "- Official GEO series metadata define 20 LCM RNA-seq samples from frontal cortex BA9 layer III pyramidal neurons.\n"
        "- Group labels are parsed from official Sample_title fields: HC/CTR, MCI and AD.\n"
        "- The raw supplementary files are gene-level `genes.results` tables, not single-cell matrices.\n"
        "- The frozen blood stage candidates and ordinal weights are tested without re-selection.\n"
        "- No age/sex covariates were available in the audited series metadata; all tests are unadjusted context analyses.\n"
        "- This analysis cannot establish blood-to-brain causality, progression prediction or CellChat/KO effects.\n",
        encoding="utf-8",
    )
    print(json.dumps(audit, indent=2))
    print(score_effects.to_string(index=False))
    print("output_dir", OUT)


if __name__ == "__main__":
    main()
