from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_stage_discovery_gse63060_20260721 import fit_linear  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]
PUBLIC = ROOT / "20260709" / "external_public_data_20260721_151749"
DISCOVERY = ROOT / "20260709" / "stage_discovery_20260721_1510"
OUT = ROOT / "20260709" / "external_stage_validation_20260721_1700"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> None:
    count_path = PUBLIC / "GSE249477_count_clean_full.gz"
    label_path = PUBLIC / "GSE249477_sample_label_audit.csv"
    freeze_path = DISCOVERY / "stage_candidate_registry_monotonic_freeze_candidate.csv"
    disc_path = DISCOVERY / "stage_discovery_gene_effects.csv"
    if not count_path.exists():
        raise FileNotFoundError(count_path)

    labels = pd.read_csv(label_path)
    labels["dk_id"] = labels["title"].str.extract(r"\[(DK\d+_\d+)\s*\]", expand=False)
    if labels["dk_id"].isna().any():
        raise ValueError("Some GSE249477 official sample titles have no DK identifier")
    labels = labels.set_index("dk_id")

    header = pd.read_csv(count_path, compression="gzip", nrows=0)
    total_cols = [c for c in header.columns if c.endswith(" - Total counts")]
    base_cols = [c for c in ["Name", "Chromosome", "Region", "Identifier", "Feature ID"] if c in header.columns]
    if len(total_cols) != len(labels):
        raise ValueError(f"Count/label sample mismatch: count_total_columns={len(total_cols)}, labels={len(labels)}")
    raw = pd.read_csv(count_path, compression="gzip", usecols=base_cols + total_cols)
    raw["Identifier"] = raw["Identifier"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    raw = raw.loc[raw["Identifier"].ne("") & raw["Identifier"].ne("nan")].copy()
    raw[total_cols] = raw[total_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    counts = raw.groupby("Identifier", sort=False)[total_cols].sum()
    dk_ids = [re.match(r"^(DK\d+_\d+)", c).group(1) for c in total_cols]
    if len(set(dk_ids)) != len(dk_ids):
        raise ValueError("Duplicate DK identifiers in count matrix")
    counts.columns = dk_ids
    counts = counts.loc[:, [x for x in dk_ids if x in labels.index]]
    labels = labels.loc[counts.columns]
    lib = counts.sum(axis=0).to_numpy(float)
    log_cpm = np.log2(counts.to_numpy(float) / lib[None, :] * 1e6 + 1.0)

    labels["stage_num"] = labels["group"].map({"HC": 0, "MCI": 1, "AD": 2})
    labels["male"] = labels["sex"].astype(str).str.lower().eq("male").astype(int)
    labels["age_c"] = labels["age"].astype(float) - labels["age"].astype(float).mean()
    X_cat = np.column_stack(
        [
            np.ones(len(labels)),
            labels["age_c"].to_numpy(float),
            labels["male"].to_numpy(float),
            labels["group"].eq("MCI").to_numpy(float),
            labels["group"].eq("AD").to_numpy(float),
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
    gene_ids = counts.index.astype(str).tolist()
    effects = pd.concat(
        [
            fit_linear(log_cpm.T, X_cat, contrasts_cat, gene_ids),
            fit_linear(log_cpm.T, X_trend, contrasts_trend, gene_ids),
        ],
        ignore_index=True,
    ).rename(columns={"feature_id": "ensembl_id"})

    freeze = pd.read_csv(freeze_path)
    disc = pd.read_csv(disc_path)
    trend = disc.loc[disc["contrast"].eq("ordinal_HC_MCI_AD"), ["gene", "estimate"]].rename(columns={"estimate": "discovery_trend_estimate"})
    freeze = freeze.merge(trend, on="gene", how="left")
    cand = freeze.dropna(subset=["gene"]).copy()
    cand["ensembl_id"] = cand.get("ensembl_id", cand["gene"])
    # Use the locked symbol-to-Ensembl lookup generated for GSE282742.
    lookup_path = ROOT / "20260709" / "progression_concordance_20260721_1600" / "ensembl_symbol_lookup.json"
    lookup = json.loads(lookup_path.read_text(encoding="utf-8"))
    cand["ensembl_id"] = cand["gene"].map(lookup)
    cand = cand.dropna(subset=["ensembl_id"])
    cand = cand[cand["ensembl_id"].isin(counts.index)].copy()
    cand.to_csv(OUT / "frozen_stage_candidates_mapped_to_GSE249477.csv", index=False)
    if cand.empty:
        raise ValueError("No frozen candidates mapped to GSE249477 count matrix")
    score_rows = []
    for _, row in cand.iterrows():
        v = log_cpm[counts.index.get_loc(row["ensembl_id"]), :]
        sd = np.std(v, ddof=1)
        score_rows.append((v - np.mean(v)) / sd if sd > 0 else np.zeros_like(v))
    Z = np.vstack(score_rows)
    weights = cand["discovery_trend_estimate"].to_numpy(float)
    weights = weights / np.sum(np.abs(weights))
    labels["stage_program_score"] = weights @ Z
    labels.reset_index().to_csv(OUT / "GSE249477_sample_stage_program_scores.csv", index=False)

    score_effects = fit_linear(
        labels[["stage_program_score"]].to_numpy(float),
        X_cat,
        contrasts_cat,
        ["frozen_stage_program_score"],
    )
    score_effects = pd.concat(
        [
            score_effects,
            fit_linear(labels[["stage_program_score"]].to_numpy(float), X_trend, contrasts_trend, ["frozen_stage_program_score"]),
        ],
        ignore_index=True,
    )
    score_effects.to_csv(OUT / "frozen_stage_program_score_external_effects.csv", index=False)
    candidate_effects = effects.loc[effects["ensembl_id"].isin(cand["ensembl_id"])].merge(
        cand[["gene", "ensembl_id", "discovery_trend_estimate"]], on="ensembl_id", how="left"
    )
    candidate_effects.to_csv(OUT / "frozen_stage_candidate_external_effects.csv", index=False)
    audit = {
        "dataset": "GSE249477",
        "count_file": str(count_path),
        "n_rows_raw": int(len(raw)),
        "n_gene_rows_aggregated": int(counts.shape[0]),
        "n_samples": int(counts.shape[1]),
        "sample_group_counts": labels["group"].value_counts().to_dict(),
        "n_frozen_candidates": int(len(freeze)),
        "n_frozen_candidates_mapped": int(len(cand)),
        "score_definition": "gene-wise log2(CPM+1) z score; weights frozen from GSE63060 ordinal estimates",
        "role": "independent cross-sectional validation",
        "old_candidates_used": False,
    }
    (OUT / "external_stage_validation_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    pd.DataFrame([audit]).to_csv(OUT / "external_stage_validation_audit.csv", index=False)
    (OUT / "analysis_registry.md").write_text(
        "# GSE249477 external stage validation registry\n\n"
        "- Frozen input: GSE63060 monotonic candidates and ordinal effect weights.\n"
        "- GSE63061 and GSE282742 results are not used to change the frozen input.\n"
        "- Count source: official GEO supplementary total-count columns.\n"
        "- Normalization: log2(CPM+1); gene rows aggregated by Ensembl Identifier.\n"
        "- Main test: age/sex-adjusted HC, MCI and AD contrasts plus ordinal trend.\n"
        "- This is cross-sectional external validation, not longitudinal progression validation.\n",
        encoding="utf-8",
    )
    print(json.dumps(audit, indent=2))
    print(score_effects.to_string(index=False))
    print("output_dir", OUT)


if __name__ == "__main__":
    main()
