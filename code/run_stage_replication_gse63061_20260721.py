from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from run_stage_discovery_gse63060_20260721 import (  # noqa: E402
    bh,
    fit_linear,
    read_bgx,
    read_series_matrix,
)


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "20260709" / "data"
DISCOVERY = ROOT / "20260709" / "stage_discovery_20260721_1510"
OUT = ROOT / "20260709" / "stage_replication_20260721_151320"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    expr_path = DATA / "GSE63061_series_matrix.txt.gz"
    design_path = DATA / "GSE63061_official_sample_labels.csv"
    ann_path = DATA / "GSE63061" / "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz"
    expr = read_series_matrix(expr_path)
    labels = pd.read_csv(design_path)
    labels["status"] = labels["status"].replace({"CTL": "HC"})
    labels["included"] = labels["included"].astype(str).str.lower()
    labels = labels.loc[labels["included"].eq("yes")].copy().set_index("sample")
    common = [s for s in expr.columns if s in labels.index]
    if len(common) != len(labels):
        raise ValueError(f"Expression/label mismatch: common={len(common)}, labels={len(labels)}")
    expr = expr.loc[:, common]
    labels = labels.loc[common]
    labels["stage_num"] = labels["status"].map({"HC": 0, "MCI": 1, "AD": 2})
    if labels["stage_num"].isna().any():
        raise ValueError(f"Unexpected labels: {labels.loc[labels['stage_num'].isna(), 'status'].unique()}")
    labels["male"] = labels["gender"].astype(str).str.lower().eq("male").astype(int)
    labels["age_c"] = labels["age"].astype(float) - labels["age"].astype(float).mean()
    X_cat = np.column_stack(
        [
            np.ones(len(labels)),
            labels["age_c"].to_numpy(float),
            labels["male"].to_numpy(float),
            labels["status"].eq("MCI").to_numpy(float),
            labels["status"].eq("AD").to_numpy(float),
        ]
    )
    X_trend = np.column_stack(
        [
            np.ones(len(labels)),
            labels["age_c"].to_numpy(float),
            labels["male"].to_numpy(float),
            labels["stage_num"].to_numpy(float),
        ]
    )
    contrasts_cat = {
        "MCI_vs_HC": np.array([0, 0, 0, 1, 0]),
        "AD_vs_MCI": np.array([0, 0, 0, -1, 1]),
        "AD_vs_HC": np.array([0, 0, 0, 0, 1]),
    }
    contrasts_trend = {"ordinal_HC_MCI_AD": np.array([0, 0, 0, 1])}
    Y_probe = expr.to_numpy(dtype=float).T
    probe_ids = expr.index.astype(str).tolist()
    probe_res = pd.concat(
        [
            fit_linear(Y_probe, X_cat, contrasts_cat, probe_ids),
            fit_linear(Y_probe, X_trend, contrasts_trend, probe_ids),
        ],
        ignore_index=True,
    )
    probe_res.to_csv(OUT / "stage_replication_probe_effects.csv", index=False)

    ann = read_bgx(ann_path)
    ann = ann.rename(columns={"Probe_Id": "feature_id", "Symbol": "symbol", "Entrez_Gene_ID": "entrez_id"})
    ann["feature_id"] = ann["feature_id"].astype(str)
    ann["symbol_clean"] = ann["symbol"].astype(str).str.strip()
    ann.loc[ann["symbol_clean"].isin(["", "---", "NA", "nan"]), "symbol_clean"] = np.nan
    ann.loc[ann["symbol_clean"].str.contains(r"[;,|]", regex=True, na=False), "symbol_clean"] = np.nan
    ann = ann.drop_duplicates("feature_id")
    probe_map = ann.set_index("feature_id")["symbol_clean"].reindex(probe_ids).dropna()
    gene_to_probes = probe_map.groupby(probe_map).groups
    genes = sorted(gene_to_probes)
    probe_pos = {probe: i for i, probe in enumerate(probe_ids)}
    gene_arrays = []
    for gene in genes:
        idx = [probe_pos[p] for p in gene_to_probes[gene] if p in probe_pos]
        gene_arrays.append(np.nanmedian(Y_probe[:, idx], axis=1))
    Y_gene = np.column_stack(gene_arrays)
    gene_res = pd.concat(
        [
            fit_linear(Y_gene, X_cat, contrasts_cat, genes),
            fit_linear(Y_gene, X_trend, contrasts_trend, genes),
        ],
        ignore_index=True,
    ).rename(columns={"feature_id": "gene"})
    gene_res["n_probes"] = gene_res["gene"].map({g: len(gene_to_probes[g]) for g in genes})
    gene_res.to_csv(OUT / "stage_replication_gene_effects.csv", index=False)

    freeze = pd.read_csv(DISCOVERY / "stage_candidate_registry_monotonic_freeze_candidate.csv")
    frozen_genes = freeze["gene"].astype(str).tolist()
    rep = gene_res.loc[gene_res["gene"].isin(frozen_genes)].copy()
    disc = pd.read_csv(DISCOVERY / "stage_discovery_gene_effects.csv")
    disc = disc.loc[disc["gene"].isin(frozen_genes)].copy()
    disc = disc.rename(columns={"estimate": "discovery_estimate", "se": "discovery_se", "p_value": "discovery_p_value", "fdr": "discovery_fdr"})
    rep = rep.rename(columns={"estimate": "replication_estimate", "se": "replication_se", "p_value": "replication_p_value", "fdr": "replication_fdr"})
    locked = disc[["gene", "contrast", "discovery_estimate", "discovery_se", "discovery_p_value", "discovery_fdr"]].merge(
        rep[["gene", "contrast", "replication_estimate", "replication_se", "replication_p_value", "replication_fdr"]],
        on=["gene", "contrast"],
        how="inner",
    )
    locked["same_direction"] = locked["discovery_estimate"] * locked["replication_estimate"] > 0
    locked["replication_nominal_p05"] = locked["replication_p_value"] < 0.05
    locked["replication_fdr05"] = locked["replication_fdr"] < 0.05
    locked.to_csv(OUT / "frozen_stage_candidate_replication_table.csv", index=False)
    primary_contrasts = ["MCI_vs_HC", "AD_vs_HC", "ordinal_HC_MCI_AD"]
    primary = locked.loc[locked["contrast"].isin(primary_contrasts)].copy()
    primary_wide = primary.pivot(index="gene", columns="contrast")
    replicated_nominal = (
        primary_wide["same_direction"].all(axis=1)
        & primary_wide["replication_nominal_p05"].all(axis=1)
    )
    replicated_fdr = (
        primary_wide["same_direction"].all(axis=1)
        & primary_wide["replication_fdr05"].all(axis=1)
    )
    summary = pd.DataFrame(
        {
            "gene": primary_wide.index,
            "all_three_same_direction_nominal_p05": replicated_nominal,
            "all_three_same_direction_fdr05": replicated_fdr,
        }
    )
    summary.to_csv(OUT / "stage_replication_summary_by_gene.csv", index=False)
    summary.loc[summary["all_three_same_direction_nominal_p05"]].to_csv(
        OUT / "stage_candidates_replicated_nominal_p05.csv", index=False
    )
    summary.loc[summary["all_three_same_direction_fdr05"]].to_csv(
        OUT / "stage_candidates_replicated_fdr05.csv", index=False
    )
    (OUT / "analysis_registry.md").write_text(
        "# GSE63061 stage replication analysis registry\n\n"
        "- Role: AddNeuroMed batch-level replication, not an independent external cohort.\n"
        "- Frozen input: stage_candidate_registry_monotonic_freeze_candidate.csv from GSE63060.\n"
        "- Model: expression ~ centered age + sex + categorical stage, HC reference; separate ordinal trend model.\n"
        "- Primary replication contrasts: MCI vs HC, AD vs HC and ordinal HC-MCI-AD trend.\n"
        "- Replication summaries report same-direction nominal P<0.05 and same-direction FDR<0.05 separately.\n"
        "- No GSE63061 result is fed back to change the GSE63060 candidate registry.\n",
        encoding="utf-8",
    )
    labels.reset_index().to_csv(OUT / "GSE63061_included_labels_used.csv", index=False)
    audit = {
        "dataset": "GSE63061",
        "series_matrix_sha256": sha256(expr_path),
        "design_sha256": sha256(design_path),
        "annotation_sha256": sha256(ann_path),
        "n_expression_probes": int(expr.shape[0]),
        "n_included_samples": int(expr.shape[1]),
        "stage_counts": labels["status"].value_counts().reindex(["HC", "MCI", "AD"]).fillna(0).astype(int).to_dict(),
        "n_gene_symbols_with_unambiguous_annotation": len(genes),
        "frozen_candidate_n": len(frozen_genes),
        "frozen_candidates_source": str(DISCOVERY / "stage_candidate_registry_monotonic_freeze_candidate.csv"),
        "model": "expression ~ age_centered + sex + stage; HC reference; separate ordinal trend model",
        "role": "AddNeuroMed batch-level replication; not external validation",
        "old_candidates_used": False,
    }
    (OUT / "stage_replication_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    pd.DataFrame([audit]).to_csv(OUT / "stage_replication_audit.csv", index=False)
    print(json.dumps(audit, indent=2))
    print("output_dir", OUT)


if __name__ == "__main__":
    main()
