"""Recompute orthogonal brain and TEMRA context for the frozen MCI panel.

The panel and its peak/trough component weights are read from the development
lock. No gene is selected from either context dataset.  Brain data are raw
GSE285831 LCM RNA-seq counts; TEMRA data are raw donor-level 10X matrices
aggregated to pseudobulk.  Outputs are context evidence, not validation or
causal claims.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_brain_stage_context_gse285831_20260722 import parse_series_labels, read_counts  # noqa: E402
from run_stage_discovery_gse63060_20260721 import fit_linear  # noqa: E402
from run_temra_stage_context_gse134578_20260722 import parse_temra_labels, read_temra_pseudobulk  # noqa: E402

ROOT = Path(os.environ.get("MCI_ROOT", Path(__file__).resolve().parents[2]))
BASE = ROOT / "20260709"
PANEL_PATH = BASE / "mci_development_optimized_panel_20260722_161919" / "development_locked_12_gene_panel.csv"
ANNOTATION_PATH = ROOT / "20250816 revise" / "GEO" / "Human.GRCh38.p13.annot.tsv.gz"
BRAIN_TAR = ROOT / "20250816 revise" / "scRNA" / "GSE285831_RAW.tar"
BRAIN_SERIES = ROOT / "20250816 revise" / "scRNA" / "GSE285831_series_matrix.txt.gz"
TEMRA_DIR = ROOT / "20250816 revise" / "scRNA" / "GSE134578_RAW"
TEMRA_SOFT = BASE / "data" / "GSE134578_family.soft.gz"
STAMP = "20260723_090020"
OUT = BASE / f"mci_current_panel_context_{STAMP}"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def resolve_panel(panel: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    exact: dict[str, list[str]] = {}
    aliases: dict[str, list[str]] = {}
    for _, row in annotation.fillna("").iterrows():
        symbol = str(row["Symbol"]).strip()
        if symbol:
            exact.setdefault(symbol.upper(), []).append(symbol)
        for token in str(row["Synonyms"]).split("|"):
            token = token.strip().upper()
            if token:
                aliases.setdefault(token, []).append(symbol)
    rows = []
    for _, row in panel.iterrows():
        gene = str(row["gene"])
        candidates = exact.get(gene.upper(), [])
        canonical = candidates[0] if candidates else ""
        status = "exact_symbol"
        note = ""
        if not canonical:
            uniq = sorted(set(aliases.get(gene.upper(), [])))
            if len(uniq) == 1:
                canonical = uniq[0]
                status = "unique_alias"
            else:
                status = "unresolved"
        hit = annotation[annotation["Symbol"].astype(str).str.upper().eq(canonical.upper())] if canonical else annotation.iloc[0:0]
        ensembl = str(hit["EnsemblGeneID"].iloc[0]) if len(hit) else ""
        if gene == "LAT1-3TM":
            canonical = ""
            ensembl = ""
            status = "unresolved_transcript_alias"
            note = "GPL10558 ILMN_138298 / RefSeq XR_001385.1; no safe Ensembl assignment"
        rows.append({
            "gene": gene,
            "trajectory_class": row["trajectory_class"],
            "canonical_symbol": canonical,
            "ensembl_id": ensembl,
            "mapping_status": status,
            "mapping_note": note,
        })
    return pd.DataFrame(rows)


def design(labels: pd.DataFrame) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    stage = labels["group"].astype(str)
    X_cat = np.column_stack([np.ones(len(labels)), stage.eq("MCI").to_numpy(float), stage.eq("AD").to_numpy(float)])
    X_trend = np.column_stack([np.ones(len(labels)), labels["stage_num"].to_numpy(float)])
    contrasts = {
        "MCI_vs_HC": np.array([0, 1, 0]),
        "AD_vs_MCI": np.array([0, -1, 1]),
        "AD_vs_HC": np.array([0, 0, 1]),
        "ordinal_HC_MCI_AD": np.array([0, 1]),
    }
    return (X_cat, X_trend), contrasts


def gene_effects(log_cpm: pd.DataFrame, labels: pd.DataFrame, feature_ids: list[str]) -> pd.DataFrame:
    labels = labels.copy()
    labels["stage_num"] = labels["group"].map({"HC": 0, "MCI": 1, "AD": 2})
    (X_cat, X_trend), contrasts = design(labels)
    cat = {k: v for k, v in contrasts.items() if k != "ordinal_HC_MCI_AD"}
    trend = {"ordinal_HC_MCI_AD": contrasts["ordinal_HC_MCI_AD"]}
    return pd.concat([
        fit_linear(log_cpm.to_numpy(float).T, X_cat, cat, feature_ids),
        fit_linear(log_cpm.to_numpy(float).T, X_trend, trend, feature_ids),
    ], ignore_index=True)


def score_from_matrix(values: pd.DataFrame, mapping: pd.DataFrame, id_col: str, mode: str) -> tuple[pd.Series, list[str]]:
    """Return fixed-12 and component-renormalized score using within-context z-scores."""
    rows = []
    ids = []
    for _, row in mapping.iterrows():
        identifier = str(row[id_col])
        if identifier and identifier in values.index:
            rows.append(row)
            ids.append(identifier)
    mapped = pd.DataFrame(rows)
    if mapped.empty:
        return pd.Series(dtype=float), []
    peak_n = int((mapping["trajectory_class"] == "MCI_peak").sum()) if mode == "fixed12" else int((mapped["trajectory_class"] == "MCI_peak").sum())
    trough_n = int((mapping["trajectory_class"] == "MCI_trough").sum()) if mode == "fixed12" else int((mapped["trajectory_class"] == "MCI_trough").sum())
    score = np.zeros(values.shape[1], dtype=float)
    for _, row in mapped.iterrows():
        v = values.loc[str(row[id_col])].to_numpy(float)
        sd = np.nanstd(v, ddof=1)
        z = (v - np.nanmean(v)) / sd if sd > 0 else np.zeros_like(v)
        w = 0.5 / peak_n if row["trajectory_class"] == "MCI_peak" else -0.5 / trough_n
        score += w * z
    return pd.Series(score, index=values.columns), mapped["gene"].tolist()


def run_brain(mapping: pd.DataFrame) -> dict:
    labels = parse_series_labels(BRAIN_SERIES)
    counts = read_counts(BRAIN_TAR, labels.index.tolist())
    labels = labels.loc[counts.columns].copy()
    labels["stage_num"] = labels["group"].map({"HC": 0, "MCI": 1, "AD": 2})
    lib = counts.sum(axis=0).to_numpy(float)
    log_cpm = pd.DataFrame(np.log2(counts.to_numpy(float) / lib[None, :] * 1e6 + 1), index=counts.index, columns=counts.columns)
    mapped = mapping[mapping["ensembl_id"].isin(counts.index)].copy()
    mapped.to_csv(OUT / "brain_current_panel_mapping_used.csv", index=False)
    effects = gene_effects(log_cpm.loc[mapped["ensembl_id"].tolist()], labels, mapped["ensembl_id"].tolist())
    effects = effects.merge(mapped[["gene", "trajectory_class", "ensembl_id"]], left_on="feature_id", right_on="ensembl_id", how="left")
    effects.to_csv(OUT / "brain_current_panel_gene_effects.csv", index=False)
    score_rows = []
    for mode in ["fixed12", "component_renorm"]:
        score, genes = score_from_matrix(log_cpm, mapping, "ensembl_id", mode)
        labels_mode = labels.copy()
        labels_mode["score"] = labels_mode.index.to_series().map(score)
        labels_mode = labels_mode.dropna(subset=["score"])
        score_df = pd.DataFrame([labels_mode["score"].to_numpy()], index=["frozen_score"], columns=labels_mode.index)
        eff = gene_effects(score_df, labels_mode, ["frozen_score"])
        eff["score_mode"] = mode
        eff["mapped_genes"] = ";".join(genes)
        score_rows.append(eff)
        labels_mode.reset_index().to_csv(OUT / f"brain_current_panel_sample_scores_{mode}.csv", index=False)
    pd.concat(score_rows, ignore_index=True).to_csv(OUT / "brain_current_panel_score_effects.csv", index=False)
    audit = {
        "dataset": "GSE285831",
        "n_samples": int(len(labels)),
        "sample_group_counts": labels["group"].value_counts().to_dict(),
        "n_gene_features": int(counts.shape[0]),
        "mapped_panel_genes": int(len(mapped)),
        "mapping_fraction": float(len(mapped) / len(mapping)),
        "tissue": "frontal cortex (BA9)",
        "cell_type": "layer III pyramidal neurons",
        "role": "brain stage-associated orthogonal context; unadjusted, not progression validation",
        "source_sha256": sha256(BRAIN_TAR),
    }
    return audit


def run_temra(mapping: pd.DataFrame) -> dict:
    labels = parse_temra_labels(TEMRA_SOFT)
    tar_paths = {p.name.split("_", 1)[0]: p for p in TEMRA_DIR.glob("GSM*_T*.tar.gz")}
    sample_rows = []
    vals_by_sample: dict[str, pd.Series] = {}
    symbols = set(mapping["canonical_symbol"].dropna().astype(str))
    for sid in labels.index:
        vals, n_cells, library_size = read_temra_pseudobulk(tar_paths[sid], symbols)
        vals_by_sample[sid] = vals
        sample_rows.append({"sample_id": sid, "n_cells": n_cells, "library_size": library_size})
    pseudobulk = pd.DataFrame(vals_by_sample).fillna(0.0).reindex(columns=labels.index).fillna(0.0)
    labels = labels.loc[pseudobulk.columns].copy()
    labels["stage_num"] = labels["group"].map({"HC": 0, "MCI": 1, "AD": 2})
    mapped = mapping[mapping["canonical_symbol"].isin(pseudobulk.index)].copy()
    mapped.to_csv(OUT / "temra_current_panel_mapping_used.csv", index=False)
    lib = pseudobulk.sum(axis=0).to_numpy(float)
    log_cpm = pd.DataFrame(np.log2(pseudobulk.to_numpy(float) / lib[None, :] * 1e6 + 1), index=pseudobulk.index, columns=pseudobulk.columns)
    effects = gene_effects(log_cpm.loc[mapped["canonical_symbol"].tolist()], labels, mapped["canonical_symbol"].tolist())
    effects = effects.merge(mapped[["gene", "trajectory_class", "canonical_symbol"]], left_on="feature_id", right_on="canonical_symbol", how="left")
    effects.to_csv(OUT / "temra_current_panel_gene_effects.csv", index=False)
    score_rows = []
    for mode in ["fixed12", "component_renorm"]:
        score, genes = score_from_matrix(log_cpm, mapping, "canonical_symbol", mode)
        labels_mode = labels.copy()
        labels_mode["score"] = labels_mode.index.to_series().map(score)
        labels_mode = labels_mode.dropna(subset=["score"])
        score_df = pd.DataFrame([labels_mode["score"].to_numpy()], index=["frozen_score"], columns=labels_mode.index)
        eff = gene_effects(score_df, labels_mode, ["frozen_score"])
        eff["score_mode"] = mode
        eff["mapped_genes"] = ";".join(genes)
        score_rows.append(eff)
        labels_mode.reset_index().to_csv(OUT / f"temra_current_panel_sample_scores_{mode}.csv", index=False)
    pd.concat(score_rows, ignore_index=True).to_csv(OUT / "temra_current_panel_score_effects.csv", index=False)
    pd.DataFrame(sample_rows).merge(labels.reset_index()[["sample_id", "group", "diagnosis", "title"]], on="sample_id").to_csv(OUT / "temra_current_panel_sample_audit.csv", index=False)
    audit = {
        "dataset": "GSE134578_TEMRA",
        "n_samples": int(len(labels)),
        "sample_group_counts": labels["group"].value_counts().to_dict(),
        "mapped_panel_genes": int(len(mapped)),
        "mapping_fraction": float(len(mapped) / len(mapping)),
        "median_cells_per_sample": float(np.median([x["n_cells"] for x in sample_rows])),
        "role": "donor-level TEMRA immune context; not bulk validation or progression validation",
        "normalization": "donor pseudobulk log2(CPM+1); cells are not independent",
    }
    return audit


def main() -> None:
    panel = pd.read_csv(PANEL_PATH)
    ann = pd.read_csv(ANNOTATION_PATH, sep="\t", compression="gzip", usecols=["Symbol", "Synonyms", "EnsemblGeneID"], dtype=str).fillna("")
    mapping = resolve_panel(panel, ann)
    mapping.to_csv(OUT / "current_panel_annotation_audit.csv", index=False)
    brain = run_brain(mapping)
    temra = run_temra(mapping)
    qa = {
        "generated": "2026-07-23 09:00:20 Asia/Shanghai",
        "panel": str(PANEL_PATH),
        "panel_size": int(len(panel)),
        "brain": brain,
        "temra": temra,
        "interpretation": "Context evidence only. No candidates or weights were selected from brain/TEMRA data; LAT1-3TM remains unresolved.",
    }
    (OUT / "context_QA.json").write_text(json.dumps(qa, indent=2, ensure_ascii=False), encoding="utf-8")
    lines = [
        f"# Current locked panel context audit ({qa['generated']})",
        "",
        "- The 12-gene panel and peak/trough roles were frozen before these analyses.",
        f"- GSE285831 mapped {brain['mapped_panel_genes']}/{len(panel)} panel genes; labels are HC/MCI/AD in frontal cortex BA9 layer III pyramidal neurons.",
        f"- GSE134578 mapped {temra['mapped_panel_genes']}/{len(panel)} panel genes; samples are donor-level CD8+ TEMRA pseudobulks.",
        "- Both datasets are unadjusted stage-context analyses. They cannot establish blood-to-brain transfer, progression prediction, causality, CellChat effects, or therapeutic efficacy.",
        "- Score files include fixed12 (missing genes retain the frozen denominator) and component_renorm sensitivity; gene-level files retain every mapped panel member and contrast.",
    ]
    (OUT / "context_QA_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(qa, indent=2, ensure_ascii=False))
    print("output_dir", OUT)


if __name__ == "__main__":
    main()
