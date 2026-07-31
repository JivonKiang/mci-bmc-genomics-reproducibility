from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
import concurrent.futures
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "20260709" / "data"
PUBLIC = ROOT / "20260709" / "external_public_data_20260721_151749"
DISCOVERY = ROOT / "20260709" / "stage_discovery_20260721_1510"
OUT = ROOT / "20260709" / "progression_concordance_20260721_1600"
OUT.mkdir(parents=True, exist_ok=True)


def bh(p: np.ndarray) -> np.ndarray:
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def fit_group_effect(Y: np.ndarray, age: np.ndarray, male: np.ndarray, group: np.ndarray, feature_ids: list[str]) -> pd.DataFrame:
    """Subject-level OLS: expression ~ centered age + sex + P-MCI indicator."""
    age_c = age - np.nanmean(age)
    X = np.column_stack([np.ones(len(age)), age_c, male, group.astype(float)])
    rank = np.linalg.matrix_rank(X)
    inv = np.linalg.inv(X.T @ X)
    beta = inv @ X.T @ Y
    resid = Y - X @ beta
    df_resid = len(age) - rank
    sigma2 = np.sum(resid * resid, axis=0) / df_resid
    c = np.array([0, 0, 0, 1.0])
    est = c @ beta
    se = np.sqrt(np.maximum(0, (c @ inv @ c) * sigma2))
    stat = np.divide(est, se, out=np.full_like(est, np.nan), where=se > 0)
    from scipy.stats import t

    p = 2 * t.sf(np.abs(stat), df_resid)
    return pd.DataFrame(
        {
            "feature_id": feature_ids,
            "estimate_P_MCI_vs_S_MCI": est,
            "se": se,
            "t": stat,
            "p_value": p,
            "fdr": bh(p),
            "n_subjects": len(age),
            "df_residual": df_resid,
        }
    )


def ensembl_lookup(symbol: str) -> str | None:
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/" + urllib.parse.quote(symbol)
    req = urllib.request.Request(url, headers={"User-Agent": "MCI-data-audit", "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as response:
            obj = json.loads(response.read().decode("utf-8"))
        return obj.get("id")
    except Exception:
        return None


def mygene_batch_lookup(symbols: list[str]) -> dict[str, str | None]:
    """Resolve a small locked symbol list without downloading a full annotation."""
    out: dict[str, str | None] = {s: None for s in symbols}
    try:
        payload = json.dumps(
            {
                "q": ",".join(symbols),
                "scopes": "symbol",
                "species": "human",
                "fields": "symbol,ensembl.gene",
                "size": 1000,
            }
        ).encode("utf-8")
        req = urllib.request.Request(
            "https://mygene.info/v3/query",
            data=payload,
            headers={"User-Agent": "MCI-data-audit", "Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=60) as response:
            payload = json.loads(response.read().decode("utf-8"))
        for hit in payload:
            symbol = hit.get("symbol")
            ens = hit.get("ensembl", {}).get("gene") if isinstance(hit.get("ensembl"), dict) else None
            if isinstance(ens, list):
                ens = ens[0] if ens else None
            if symbol in out and ens:
                out[symbol] = ens
    except Exception:
        pass
    return out


def mygene_alias_lookup(symbols: list[str]) -> dict[str, str | None]:
    out: dict[str, str | None] = {s: None for s in symbols}

    def one(symbol: str) -> tuple[str, str | None]:
        url = (
            "https://mygene.info/v3/query?q="
            + urllib.parse.quote(symbol)
            + "&species=human&fields=symbol,alias,ensembl.gene&size=10"
        )
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "MCI-data-audit"})
            with urllib.request.urlopen(req, timeout=30) as response:
                hits = json.loads(response.read().decode("utf-8")).get("hits", [])
            for hit in hits:
                aliases = hit.get("alias", [])
                aliases = [aliases] if isinstance(aliases, str) else aliases
                if hit.get("symbol", "").upper() != symbol.upper() and symbol.upper() not in {
                    str(a).upper() for a in aliases
                }:
                    continue
                ens = hit.get("ensembl", {})
                ens = ens.get("gene") if isinstance(ens, dict) else None
                if isinstance(ens, list):
                    ens = ens[0] if ens else None
                if ens:
                    return symbol, ens
        except Exception:
            pass
        return symbol, None

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        for symbol, ens in executor.map(one, symbols):
            out[symbol] = ens
    return out


def main() -> None:
    count_path = PUBLIC / "GSE282742_Expected_count.full.gz"
    audit_path = PUBLIC / "GSE282742_sample_subject_audit.csv"
    freeze_path = DISCOVERY / "stage_candidate_registry_monotonic_freeze_candidate.csv"
    disc_path = DISCOVERY / "stage_discovery_gene_effects.csv"
    if not count_path.exists():
        raise FileNotFoundError(count_path)

    sample_audit = pd.read_csv(audit_path)
    counts = pd.read_csv(count_path, sep="\t", compression="gzip", index_col=0)
    counts.index = counts.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    counts = counts.loc[~counts.index.duplicated(keep="first")]
    sample_audit = sample_audit.set_index("vgh_id")
    common_samples = [s for s in counts.columns if s in sample_audit.index]
    if len(common_samples) != len(counts.columns):
        raise ValueError(f"Count/audit sample mismatch: count={len(counts.columns)}, common={len(common_samples)}")
    counts = counts.loc[:, common_samples]
    sample_audit = sample_audit.loc[common_samples].copy()
    libsize = counts.sum(axis=0).to_numpy(float)
    log_cpm = np.log2(counts.to_numpy(float) / libsize[None, :] * 1e6 + 1.0)

    freeze = pd.read_csv(freeze_path)
    disc = pd.read_csv(disc_path)
    trend = disc.loc[disc["contrast"].eq("ordinal_HC_MCI_AD"), ["gene", "estimate"]].rename(columns={"estimate": "discovery_trend_estimate"})
    freeze = freeze.merge(trend, on="gene", how="left")
    cache_path = OUT / "ensembl_symbol_lookup.json"
    cache = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else {}
    missing_symbols = [s for s in freeze["gene"].astype(str) if s not in cache]
    if missing_symbols:
        cache.update(mygene_batch_lookup(missing_symbols))
    unresolved = [s for s in freeze["gene"].astype(str) if cache.get(s) is None]
    if unresolved:
        cache.update(mygene_alias_lookup(unresolved))
    cache_path.write_text(json.dumps(cache, indent=2), encoding="utf-8")
    freeze["ensembl_id"] = freeze["gene"].map(cache)
    freeze["ensembl_id"] = freeze["ensembl_id"].astype(object)
    freeze.to_csv(OUT / "frozen_stage_candidates_with_ensembl.csv", index=False)

    # Score is locked from discovery effect direction/weight. Counts are used
    # only after log-CPM conversion; gene-wise z scores make this a cross-platform
    # concordance test, not a clinical prediction score.
    count_ids = pd.Series(counts.index, index=counts.index)
    candidate = freeze.dropna(subset=["ensembl_id"]).copy()
    candidate = candidate[candidate["ensembl_id"].isin(count_ids.index)].copy()
    if candidate.empty:
        raise ValueError("No frozen candidate genes mapped to GSE282742 count matrix")
    gene_positions = {g: i for i, g in enumerate(counts.index)}
    candidate["row"] = candidate["ensembl_id"].map(gene_positions)
    candidate = candidate.dropna(subset=["row"]).copy()
    candidate["row"] = candidate["row"].astype(int)
    score_values = []
    for row in candidate["row"]:
        v = log_cpm[row, :]
        sd = np.nanstd(v, ddof=1)
        score_values.append((v - np.nanmean(v)) / sd if sd > 0 else np.zeros_like(v))
    Z = np.vstack(score_values)
    weights = candidate["discovery_trend_estimate"].to_numpy(float)
    weights = weights / np.sum(np.abs(weights))
    score = weights @ Z
    score_table = sample_audit.reset_index()[["vgh_id", "sample", "group", "subject_id", "age", "sex"]].copy()
    score_table["stage_program_score"] = score
    score_table.to_csv(OUT / "GSE282742_sample_stage_program_scores.csv", index=False)

    # Use the earliest available sample per subject for the cross-sectional P/S
    # comparison; repeated samples are retained in a separate transition table.
    base = score_table.sort_values(["subject_id", "age", "vgh_id"], na_position="last").groupby("subject_id", as_index=False).first()
    ps = base.loc[base["group"].isin(["P-MCI", "S-MCI"])].copy()
    ps["P_MCI_indicator"] = ps["group"].eq("P-MCI").astype(int)
    candidate_scores = []
    for j, row in candidate.reset_index(drop=True).iterrows():
        vals = Z[j, :]
        by_sample = pd.Series(vals, index=counts.columns)
        ps_vals = ps["vgh_id"].map(by_sample).to_numpy(float)
        candidate_scores.append(ps_vals)
    ps_matrix = np.column_stack(candidate_scores)
    ps_effects = fit_group_effect(
        ps_matrix,
        ps["age"].fillna(ps["age"].median()).to_numpy(float),
        ps["sex"].eq("M").astype(int).to_numpy(),
        ps["P_MCI_indicator"].to_numpy(),
        candidate["gene"].astype(str).tolist(),
    )
    ps_effects.to_csv(OUT / "frozen_candidate_progression_effects_subject_level.csv", index=False)

    # Independent secondary discovery: all expressed Ensembl features are
    # tested within GSE282742, but these results never feed back into the frozen
    # stage score. This is exploratory progression-specific discovery.
    sample_position = {s: i for i, s in enumerate(counts.columns)}
    ps_positions = [sample_position[s] for s in ps["vgh_id"]]
    Y_ps_all = log_cpm[:, ps_positions].T
    all_progression = fit_group_effect(
        Y_ps_all,
        ps["age"].fillna(ps["age"].median()).to_numpy(float),
        ps["sex"].eq("M").astype(int).to_numpy(),
        ps["P_MCI_indicator"].to_numpy(),
        counts.index.astype(str).tolist(),
    )
    all_progression.to_csv(OUT / "progression_discovery_all_ensembl_features_subject_level.csv", index=False)

    # Paired P-MCI -> AD transition summary for the 11 subjects represented in
    # both groups. Each subject contributes its earliest P-MCI and latest AD sample.
    transition_rows = []
    for subject, g in score_table.groupby("subject_id"):
        p = g.loc[g["group"].eq("P-MCI")].sort_values(["age", "vgh_id"], na_position="last")
        a = g.loc[g["group"].eq("AD")].sort_values(["age", "vgh_id"], na_position="last")
        if len(p) and len(a):
            p0, a1 = p.iloc[0], a.iloc[-1]
            transition_rows.append(
                {
                    "subject_id": subject,
                    "P_MCI_sample": p0["vgh_id"],
                    "AD_sample": a1["vgh_id"],
                    "P_MCI_age": p0["age"],
                    "AD_age": a1["age"],
                    "P_MCI_score": p0["stage_program_score"],
                    "AD_score": a1["stage_program_score"],
                    "delta_AD_minus_P_MCI": a1["stage_program_score"] - p0["stage_program_score"],
                }
            )
    transition = pd.DataFrame(transition_rows)
    transition.to_csv(OUT / "GSE282742_subject_level_P_MCI_to_AD_transitions.csv", index=False)
    if len(transition) >= 2:
        paired = ttest_rel(transition["AD_score"], transition["P_MCI_score"], nan_policy="omit")
        paired_summary = {
            "n_transition_subjects": int(len(transition)),
            "mean_delta_AD_minus_P_MCI": float(transition["delta_AD_minus_P_MCI"].mean()),
            "paired_t": float(paired.statistic),
            "paired_p": float(paired.pvalue),
        }
    else:
        paired_summary = {"n_transition_subjects": int(len(transition))}

    audit = {
        "dataset": "GSE282742",
        "expression_source": str(count_path),
        "expression_type": "GEO processed expected counts; raw human reads unavailable",
        "n_count_genes": int(counts.shape[0]),
        "n_samples": int(counts.shape[1]),
        "sample_group_counts": sample_audit["group"].value_counts().to_dict(),
        "n_unique_subjects": int(sample_audit["subject_id"].nunique()),
        "n_repeated_sample_rows": int(sample_audit["subject_id"].duplicated(keep=False).sum()),
        "n_frozen_candidates": int(len(freeze)),
        "n_frozen_candidates_mapped": int(len(candidate)),
        "n_subjects_P_or_S_baseline": int(len(ps)),
        "n_all_feature_progression_tests": int(len(all_progression)),
        "paired_transition": paired_summary,
        "score_definition": "gene-wise log-CPM z score within GSE282742; weights frozen from GSE63060 ordinal estimates; concordance test only",
        "old_candidates_used": False,
    }
    (OUT / "progression_concordance_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    pd.DataFrame([audit]).to_csv(OUT / "progression_concordance_audit.csv", index=False)
    (OUT / "analysis_registry.md").write_text(
        "# GSE282742 progression concordance analysis registry\n\n"
        "- Input program: GSE63060 monotonic freeze candidates and discovery ordinal effect weights.\n"
        "- No GSE282742 expression result was used to select genes, weights or cutoffs.\n"
        "- Expression source: GEO processed expected counts; human raw reads are not public.\n"
        "- Normalization: log2(CPM + 1), followed by gene-wise within-cohort z scoring.\n"
        "- Primary test: subject-level earliest-sample P-MCI vs S-MCI, adjusted for age and sex.\n"
        "- Repeated samples are not independent n; P-MCI-to-AD pairs are reported separately.\n"
        "- The score is a cross-platform progression-concordance score, not a clinical prediction model.\n",
        encoding="utf-8",
    )
    print(json.dumps(audit, indent=2))
    print("output_dir", OUT)


if __name__ == "__main__":
    main()
