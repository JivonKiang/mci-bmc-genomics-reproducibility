from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import os
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_stage_discovery_gse63060_20260721 import fit_linear  # noqa: E402


ROOT = Path(os.environ.get("MCI_ROOT", Path(__file__).resolve().parents[2]))
BASE = ROOT / "20260709"
RAW_DIR = ROOT / "20250816 revise" / "scRNA" / "GSE134578_RAW"
SOFT_PATH = BASE / "data" / "GSE134578_family.soft.gz"
LOOKUP_PATH = BASE / "progression_concordance_20260721_1600" / "ensembl_symbol_lookup.json"
DISC_PATH = BASE / "stage_discovery_20260721_1510" / "stage_discovery_gene_effects.csv"
FREEZE_PATH = BASE / "stage_discovery_20260721_1510" / "stage_candidate_registry_monotonic_freeze_candidate.csv"
OUT = BASE / "temra_stage_context_gse134578_20260722"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_temra_labels(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None and current.get("cell_type") == "CD8+ TEMRA cells":
                    records.append(current)
                current = {"sample_id": line.split("=", 1)[1].strip()}
            elif current is not None and line.startswith("!Sample_title"):
                current["title"] = line.split("=", 1)[1].strip().strip('"')
            elif current is not None and line.startswith("!Sample_characteristics_ch1"):
                value = line.split("=", 1)[1].strip().strip('"')
                if value.lower().startswith("diagnosis:"):
                    current["diagnosis"] = value.split(":", 1)[1].strip()
                elif value.lower().startswith("cell type:"):
                    current["cell_type"] = value.split(":", 1)[1].strip()
    if current is not None and current.get("cell_type") == "CD8+ TEMRA cells":
        records.append(current)
    labels = pd.DataFrame(records)
    labels["group"] = labels["diagnosis"].replace({"Healthy": "HC"})
    if set(labels["group"]) != {"HC", "MCI", "AD"}:
        raise ValueError(f"Unexpected TEMRA labels: {labels['group'].unique()}")
    return labels.set_index("sample_id").sort_index()


def read_temra_pseudobulk(tar_path: Path, gene_symbols: set[str]) -> tuple[pd.Series, int, int]:
    with tarfile.open(tar_path, "r:gz") as tar:
        genes = tar.extractfile("genes.tsv")
        barcodes = tar.extractfile("barcodes.tsv")
        matrix = tar.extractfile("matrix.mtx")
        if genes is None or barcodes is None or matrix is None:
            raise ValueError(f"Incomplete 10X tar: {tar_path.name}")
        gene_df = pd.read_csv(genes, sep="\t", header=None, usecols=[0, 1], names=["ensembl", "symbol"])
        cell_count = sum(1 for _ in barcodes)
        mat = mmread(io.BytesIO(matrix.read())).tocsr()
    if mat.shape[0] != len(gene_df) or mat.shape[1] != cell_count:
        raise ValueError(f"10X dimension mismatch in {tar_path.name}: matrix={mat.shape}, genes={len(gene_df)}, cells={cell_count}")
    library_size = np.asarray(mat.sum(axis=0)).ravel().astype(float)
    wanted = gene_df["symbol"].isin(gene_symbols)
    selected = gene_df.loc[wanted].copy()
    if selected.empty:
        return pd.Series(dtype=float), int(cell_count), int(library_size.sum())
    # Sum duplicate feature rows by gene symbol before sample-level normalization.
    values: dict[str, float] = {}
    for idx, symbol in zip(selected.index, selected["symbol"]):
        values[str(symbol)] = values.get(str(symbol), 0.0) + float(mat.getrow(int(idx)).sum())
    return pd.Series(values, dtype=float), int(cell_count), int(library_size.sum())


def main() -> None:
    labels = parse_temra_labels(SOFT_PATH)
    tar_paths = {p.name.split("_", 1)[0]: p for p in RAW_DIR.glob("GSM*_T*.tar.gz")}
    if set(labels.index) != set(tar_paths):
        raise ValueError(f"Official TEMRA samples and local tar files differ: labels={len(labels)}, files={len(tar_paths)}")

    freeze = pd.read_csv(FREEZE_PATH)
    disc = pd.read_csv(DISC_PATH)
    trend = disc.loc[disc["contrast"].eq("ordinal_HC_MCI_AD"), ["gene", "estimate"]].rename(
        columns={"estimate": "discovery_trend_estimate"}
    )
    cand = freeze[["gene"]].merge(trend, on="gene", how="left")
    gene_symbols = set(cand["gene"].dropna())
    sample_rows = []
    bulk_rows: dict[str, pd.Series] = {}
    for sid, tar_path in sorted(tar_paths.items()):
        vals, n_cells, total_counts = read_temra_pseudobulk(tar_path, gene_symbols)
        bulk_rows[sid] = vals
        sample_rows.append({"sample_id": sid, "n_cells": n_cells, "library_size": total_counts})
    pseudobulk = pd.DataFrame(bulk_rows).fillna(0.0).reindex(columns=labels.index).fillna(0.0)
    mapped = [g for g in gene_symbols if g in pseudobulk.index]
    cand = cand[cand["gene"].isin(mapped)].copy()
    cand.to_csv(OUT / "frozen_stage_candidates_mapped_to_TEMRA.csv", index=False)
    if cand.empty:
        raise ValueError("No frozen candidates mapped to TEMRA symbols")

    lib = pseudobulk.sum(axis=0).to_numpy(float)
    log_cpm = np.log2(pseudobulk.to_numpy(float) / lib[None, :] * 1e6 + 1.0)
    labels = labels.loc[pseudobulk.columns].copy()
    labels["stage_num"] = labels["group"].map({"HC": 0, "MCI": 1, "AD": 2})
    X_cat = np.column_stack(
        [np.ones(len(labels)), labels["group"].eq("MCI").to_numpy(float), labels["group"].eq("AD").to_numpy(float)]
    )
    X_trend = np.column_stack([np.ones(len(labels)), labels["stage_num"].to_numpy(float)])
    contrasts_cat = {
        "MCI_vs_HC": np.array([0, 1, 0]),
        "AD_vs_MCI": np.array([0, -1, 1]),
        "AD_vs_HC": np.array([0, 0, 1]),
    }
    contrasts_trend = {"ordinal_HC_MCI_AD": np.array([0, 1])}

    symbol_to_row = {symbol: i for i, symbol in enumerate(pseudobulk.index)}
    z_rows = []
    for gene in cand["gene"]:
        v = log_cpm[symbol_to_row[gene], :]
        sd = np.std(v, ddof=1)
        z_rows.append((v - np.mean(v)) / sd if sd > 0 else np.zeros_like(v))
    weights = cand["discovery_trend_estimate"].to_numpy(float)
    weights = weights / np.sum(np.abs(weights))
    labels["stage_program_score"] = weights @ np.vstack(z_rows)
    labels.reset_index().to_csv(OUT / "TEMRA_sample_stage_program_scores.csv", index=False)
    pd.DataFrame(sample_rows).merge(labels.reset_index()[["sample_id", "group", "diagnosis", "title"]], on="sample_id").to_csv(
        OUT / "TEMRA_sample_audit.csv", index=False
    )

    score_effects = pd.concat(
        [
            fit_linear(labels[["stage_program_score"]].to_numpy(float), X_cat, contrasts_cat, ["frozen_stage_program_score"]),
            fit_linear(labels[["stage_program_score"]].to_numpy(float), X_trend, contrasts_trend, ["frozen_stage_program_score"]),
        ],
        ignore_index=True,
    )
    score_effects.to_csv(OUT / "frozen_stage_program_TEMRA_context_effects.csv", index=False)

    audit = {
        "dataset": "GSE134578_TEMRA",
        "source_family_soft": str(SOFT_PATH),
        "n_samples": int(len(labels)),
        "sample_group_counts": labels["group"].value_counts().to_dict(),
        "n_mapped_frozen_candidates": int(len(cand)),
        "median_cells_per_sample": float(np.median([x["n_cells"] for x in sample_rows])),
        "total_local_tar_bytes": int(sum(p.stat().st_size for p in tar_paths.values())),
        "analysis_role": "donor-level TEMRA immune context; not blood bulk external validation and not progression validation",
        "normalization": "sample-level pseudobulk log2(CPM+1); cells are not treated as independent observations",
        "old_cell_level_claims_used": False,
    }
    (OUT / "temra_context_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    pd.DataFrame([audit]).to_csv(OUT / "temra_context_audit.csv", index=False)
    (OUT / "analysis_registry.md").write_text(
        "# GSE134578 TEMRA stage-context registry\n\n"
        "- Official GEO family metadata define 13 CD8+ TEMRA donor samples with HC/MCI/AD diagnosis fields.\n"
        "- Raw 10X matrices are aggregated to one pseudobulk profile per donor/sample.\n"
        "- Frozen GSE63060 stage candidates and ordinal weights are tested without re-selection.\n"
        "- No cell-level inferential test is used; the output is immune cellular-context evidence only.\n"
        "- This analysis does not validate blood bulk expression, MCI-to-AD progression, causality or CellChat effects.\n",
        encoding="utf-8",
    )
    print(json.dumps(audit, indent=2))
    print(score_effects.to_string(index=False))
    print("output_dir", OUT)


if __name__ == "__main__":
    main()
