from __future__ import annotations

import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "20260709" / "data"
OUT = ROOT / "20260709" / "stage_discovery_20260721_1510"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_series_matrix(path: Path) -> pd.DataFrame:
    # GEO series matrices contain metadata lines beginning with '!'; pandas
    # can skip these while retaining the expression table.
    return pd.read_csv(path, sep="\t", comment="!", index_col=0, compression="gzip")


def read_bgx(path: Path) -> pd.DataFrame:
    header_line = None
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Probes]":
                header_line = fh.readline().rstrip("\n").split("\t")
                break
        if header_line is None:
            raise ValueError(f"[Probes] section not found in {path}")
        for line in fh:
            if not line or line.startswith("["):
                break
            values = line.split("\t")
            if len(values) < len(header_line):
                values.extend([""] * (len(header_line) - len(values)))
            rows.append(values[: len(header_line)])
    ann = pd.DataFrame(rows, columns=header_line)
    return ann


def bh(p: np.ndarray) -> np.ndarray:
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def fit_linear(Y: np.ndarray, X: np.ndarray, contrasts: dict[str, np.ndarray], feature_ids: list[str]) -> pd.DataFrame:
    """Vectorized OLS for samples x features, returning one row per feature/contrast."""
    n, p = X.shape
    rank = np.linalg.matrix_rank(X)
    if rank != p:
        raise ValueError(f"Design matrix is rank deficient: rank={rank}, p={p}")
    xtx_inv = np.linalg.inv(X.T @ X)
    beta = xtx_inv @ X.T @ Y
    resid = Y - X @ beta
    df_resid = n - rank
    sigma2 = np.sum(resid * resid, axis=0) / df_resid
    rows = []
    for name, contrast in contrasts.items():
        c = np.asarray(contrast, dtype=float)
        est = c @ beta
        se = np.sqrt(np.maximum(0, (c @ xtx_inv @ c) * sigma2))
        stat = np.divide(est, se, out=np.full_like(est, np.nan), where=se > 0)
        pval = 2 * t.sf(np.abs(stat), df_resid)
        rows.append(
            pd.DataFrame(
                {
                    "feature_id": feature_ids,
                    "contrast": name,
                    "estimate": est,
                    "se": se,
                    "t": stat,
                    "p_value": pval,
                    "fdr": bh(pval),
                    "n_samples": n,
                    "df_residual": df_resid,
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def clean_symbol(value: object) -> str | None:
    if pd.isna(value):
        return None
    s = str(value).strip()
    if not s or s in {"---", "NA", "nan"}:
        return None
    # Illumina annotations occasionally contain multiple symbols. Keep only
    # unambiguous symbols for the primary gene-level table.
    if re.search(r"[;,|]", s):
        return None
    return s


def main() -> None:
    expr_path = DATA / "GSE63060_series_matrix.txt.gz"
    design_path = DATA / "GSE63060_official_sample_labels.csv"
    ann_path = DATA / "GSE63063" / "GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx.gz"
    expr = read_series_matrix(expr_path)
    labels = pd.read_csv(design_path).rename(columns={"sample": "sample_id"})
    labels["status"] = labels["status"].replace({"CTL": "HC"})
    labels["included"] = labels["included"].astype(str).str.lower()
    labels = labels.loc[labels["included"].eq("yes")].copy()
    labels = labels.set_index("sample_id")
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

    # Intercept, age, sex, MCI and AD; HC is the reference level.
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
    probe_res.to_csv(OUT / "stage_discovery_probe_effects.csv", index=False)

    ann = read_bgx(ann_path)
    ann = ann.rename(columns={"Probe_Id": "feature_id", "Symbol": "symbol", "Entrez_Gene_ID": "entrez_id"})
    ann["feature_id"] = ann["feature_id"].astype(str)
    ann["symbol_clean"] = ann["symbol"].map(clean_symbol)
    ann = ann.drop_duplicates("feature_id")
    probe_map = ann.set_index("feature_id")["symbol_clean"].reindex(probe_ids)
    probe_map = probe_map.dropna()
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
    )
    gene_res = gene_res.rename(columns={"feature_id": "gene"})
    gene_res["n_probes"] = gene_res["gene"].map({g: len(gene_to_probes[g]) for g in genes})
    gene_res.to_csv(OUT / "stage_discovery_gene_effects.csv", index=False)

    # A transparent, non-hypothesis-driven candidate flag for downstream audit.
    # The thresholds are recorded here and must be frozen before opening any
    # validation cohort; they are not used to assert final biological hits.
    gene_res["candidate_primary"] = (
        gene_res["fdr"].lt(0.05)
        & gene_res["estimate"].abs().ge(0.20)
        & gene_res["n_probes"].ge(1)
    )
    gene_res.to_csv(OUT / "stage_discovery_gene_effects_with_candidate_flag.csv", index=False)

    # Collapse the four pre-registered contrasts into an auditable candidate
    # registry. The monotonic tier is only a discovery-stage freeze candidate;
    # it must still be tested without re-selection in GSE63061 and GSE249477.
    wide_est = gene_res.pivot(index="gene", columns="contrast", values="estimate")
    wide_fdr = gene_res.pivot(index="gene", columns="contrast", values="fdr")
    required = ["ordinal_HC_MCI_AD", "MCI_vs_HC", "AD_vs_MCI", "AD_vs_HC"]
    pass_flags = pd.DataFrame(index=wide_est.index)
    for contrast in required:
        pass_flags[f"{contrast}_pass"] = (
            wide_fdr[contrast].lt(0.05) & wide_est[contrast].abs().ge(0.20)
        )
    pass_flags["stage_monotonic_pass"] = (
        pass_flags["ordinal_HC_MCI_AD_pass"]
        & pass_flags["MCI_vs_HC_pass"]
        & pass_flags["AD_vs_HC_pass"]
        & (wide_est["MCI_vs_HC"] * wide_est["AD_vs_HC"] > 0)
    )
    pass_flags["early_transition_pass"] = pass_flags["MCI_vs_HC_pass"]
    pass_flags["late_transition_pass"] = pass_flags["AD_vs_MCI_pass"]
    pass_flags["any_stage_pass"] = pass_flags[[f"{c}_pass" for c in required]].any(axis=1)
    registry = pd.concat(
        [
            wide_est.add_prefix("estimate_"),
            wide_fdr.add_prefix("fdr_"),
            pass_flags,
        ],
        axis=1,
    ).reset_index(names="gene")
    registry["n_probes"] = registry["gene"].map({g: len(gene_to_probes[g]) for g in genes})
    registry.to_csv(OUT / "stage_candidate_registry_preliminary.csv", index=False)
    registry.loc[registry["stage_monotonic_pass"]].to_csv(
        OUT / "stage_candidate_registry_monotonic_freeze_candidate.csv", index=False
    )
    (OUT / "analysis_registry.md").write_text(
        "# GSE63060 stage discovery analysis registry\n\n"
        "- Discovery dataset: GSE63060, included=yes samples only.\n"
        "- Primary contrasts: ordinal HC-MCI-AD trend; MCI vs HC; AD vs MCI; AD vs HC.\n"
        "- Model: expression ~ centered age + sex + categorical stage, HC reference; ordinal trend uses numeric stage 0/1/2.\n"
        "- Gene aggregation: median of unambiguous GPL6947 probe symbols.\n"
        "- Exploratory candidate threshold: FDR < 0.05 and absolute estimate >= 0.20.\n"
        "- Monotonic freeze-candidate tier: ordinal, MCI-vs-HC and AD-vs-HC all pass, with the same HC-to-disease direction.\n"
        "- GSE63061 and GSE249477 results were not used to select this registry.\n"
        "- TOMM7, RPS24 and RPS27L were not pre-specified and were not used to restrict the search space.\n"
        "- This is a discovery-stage registry; validation and final candidate status require locked testing in independent data.\n",
        encoding="utf-8",
    )
    ann.loc[ann["feature_id"].isin(probe_ids)].to_csv(OUT / "GPL6947_probe_annotation_used.csv", index=False)
    labels.reset_index().to_csv(OUT / "GSE63060_included_labels_used.csv", index=False)

    counts = labels["status"].value_counts().reindex(["HC", "MCI", "AD"]).fillna(0).astype(int).to_dict()
    audit = {
        "dataset": "GSE63060",
        "series_matrix_sha256": sha256(expr_path),
        "design_sha256": sha256(design_path),
        "annotation_sha256": sha256(ann_path),
        "n_expression_probes": int(expr.shape[0]),
        "n_included_samples": int(expr.shape[1]),
        "stage_counts": counts,
        "n_gene_symbols_with_unambiguous_annotation": len(genes),
        "model": "expression ~ age_centered + sex + stage; HC reference; separate ordinal trend model",
        "gene_aggregation": "median of unambiguous GPL6947 probe symbols",
        "candidate_flag": "FDR < 0.05 and absolute gene-level estimate >= 0.20; exploratory until registry sign-off",
        "old_candidates_used": False,
    }
    (OUT / "stage_discovery_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    pd.DataFrame([audit]).to_csv(OUT / "stage_discovery_audit.csv", index=False)
    print(json.dumps(audit, indent=2))
    print("output_dir", OUT)


if __name__ == "__main__":
    main()
