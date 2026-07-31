"""Audit current-panel cis-eQTL availability through the eQTL Catalogue v3 API.

The v3 association routes do not expose a p-value filter.  This script therefore
pages each fixed dataset/gene query, filters returned rows locally at P<=5e-8,
and retains the raw returned rows plus a machine-readable QA manifest.  An
instrument is never promoted from an upstream error or an incomplete page.
"""
from __future__ import annotations

import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import pandas as pd

ROOT = Path(os.environ.get("MCI_PROJECT_ROOT", "."))
CLIENT = Path(os.environ.get("EQTL_CATALOGUE_REQUEST_CLIENT", "rest_request.py"))
MAPPING = ROOT / "mci_key_analyses_20260723_091500" / "locked_12_panel_local_annotation_audit.csv"
DATASETS = {
    "GTEx_brain_DLPFC": ("QTD000176", "GTEx", "brain (DLPFC)", 175),
    "GTEx_brain_cortex": ("QTD000171", "GTEx", "brain (cortex)", 205),
    "GTEx_brain_hippocampus": ("QTD000181", "GTEx", "brain (hippocampus)", 165),
    "BLUEPRINT_monocyte": ("QTD000021", "BLUEPRINT", "monocyte", 191),
    "BLUEPRINT_CD4_T": ("QTD000031", "BLUEPRINT", "CD4+ T cell", 167),
}
P_CUTOFF = 5e-8
PAGE_SIZE = 1000
MAX_PAGES = 25


def query_page(gene: str, ensembl: str, dataset_name: str, dataset_id: str, start: int) -> dict:
    payload = {
        "base_url": "https://www.ebi.ac.uk/eqtl/api",
        "path": f"v3/datasets/{dataset_id}/associations",
        "params": {"gene_id": ensembl, "start": start, "size": PAGE_SIZE},
        "max_items": PAGE_SIZE,
        "max_depth": 10,
        "timeout_sec": 90,
    }
    last_error = ""
    for attempt in range(3):
        try:
            run = subprocess.run(
                ["python", str(CLIENT)],
                input=json.dumps(payload),
                text=True,
                capture_output=True,
                timeout=120,
            )
            response = json.loads(run.stdout.strip() or "{}")
            if response.get("ok"):
                records = response.get("records", [])
                return {
                    "ok": True,
                    "records": records if isinstance(records, list) else [],
                    "available": response.get("record_count_available"),
                    "status_code": response.get("status_code"),
                    "error": "",
                }
            last_error = str(response.get("error", {}).get("message", "unknown API error"))
        except Exception as exc:  # network/client timeout; retained as a status
            last_error = str(exc)
        time.sleep(2 ** attempt)
    return {"ok": False, "records": [], "available": None, "status_code": None, "error": last_error}


def query_gene_dataset(job: tuple[str, str, str, str, str, int]) -> tuple[dict, list[dict]]:
    gene, ensembl, dataset_name, dataset_id, study, tissue, sample_size = job
    all_records: list[dict] = []
    starts: list[int] = []
    available = None
    status = "ok"
    error = ""
    for page in range(MAX_PAGES):
        start = page * PAGE_SIZE
        starts.append(start)
        result = query_page(gene, ensembl, dataset_name, dataset_id, start)
        if not result["ok"]:
            status = "client_or_upstream_error"
            error = result["error"]
            break
        page_records = result["records"]
        available = result.get("available") if result.get("available") is not None else available
        all_records.extend(page_records)
        if not page_records or len(page_records) < PAGE_SIZE:
            break
        if available is not None and start + len(page_records) >= int(available):
            break
    if status == "ok" and not all_records:
        status = "upstream_empty"
    filtered = []
    for row in all_records:
        if not isinstance(row, dict):
            continue
        try:
            pvalue = float(row.get("pvalue"))
        except (TypeError, ValueError):
            continue
        if pvalue <= P_CUTOFF:
            item = dict(row)
            item.update({"gene": gene, "ensembl_id": ensembl, "dataset": dataset_name, "dataset_id": dataset_id, "study": study, "tissue": tissue, "sample_size": sample_size})
            filtered.append(item)
    raw_rows = []
    for row in all_records:
        if isinstance(row, dict):
            item = dict(row)
            item.update({"gene": gene, "ensembl_id": ensembl, "dataset": dataset_name, "dataset_id": dataset_id, "study": study, "tissue": tissue, "sample_size": sample_size})
            raw_rows.append(item)
    manifest = {
        "gene": gene,
        "ensembl_id": ensembl,
        "dataset": dataset_name,
        "dataset_id": dataset_id,
        "study": study,
        "tissue": tissue,
        "sample_size": sample_size,
        "status": status,
        "pages_requested": len(starts),
        "records_returned": len(all_records),
        "record_count_available_last": available,
        "min_pvalue_returned": min((float(r["pvalue"]) for r in all_records if r.get("pvalue") is not None), default=None),
        "n_genomewide_rows": len(filtered),
        "error": error,
    }
    return manifest, raw_rows


def main() -> None:
    stamp = datetime.now(ZoneInfo("Asia/Shanghai")).strftime("%Y%m%d_%H%M%S")
    out = ROOT / f"mci_mr_gate_audit_v3_{stamp}"
    out.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv(MAPPING)
    genes = mapping[mapping["ensembl_id"].fillna("").astype(str).ne("")][["gene", "ensembl_id"]].drop_duplicates()
    jobs = []
    for _, row in genes.iterrows():
        for dataset_name, (dataset_id, study, tissue, sample_size) in DATASETS.items():
            jobs.append((str(row["gene"]), str(row["ensembl_id"]), dataset_name, dataset_id, study, tissue, sample_size))
    manifests: list[dict] = []
    raw_rows: list[dict] = []
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = [pool.submit(query_gene_dataset, job) for job in jobs]
        for future in as_completed(futures):
            manifest, rows = future.result()
            manifests.append(manifest)
            raw_rows.extend(rows)
    manifest_df = pd.DataFrame(manifests).sort_values(["gene", "dataset"])
    manifest_df.to_csv(out / "current_panel_mr_gate_v3_manifest.csv", index=False)
    with (out / "current_panel_mr_gate_v3_raw_records.jsonl").open("w", encoding="utf-8") as handle:
        for row in raw_rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")
    sig_df = pd.DataFrame([row for row in raw_rows if float(row.get("pvalue", 1.0)) <= P_CUTOFF])
    if not sig_df.empty:
        sig_df.to_csv(out / "current_panel_mr_gate_v3_genomewide_rows.csv", index=False)
    summary = {
        "generated": datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(),
        "genes_probed": int(genes.shape[0]),
        "datasets_probed": len(DATASETS),
        "queries": len(manifest_df),
        "successful_queries": int((manifest_df["status"] == "ok").sum()),
        "empty_queries": int((manifest_df["status"] == "upstream_empty").sum()),
        "error_queries": int((manifest_df["status"] == "client_or_upstream_error").sum()),
        "queries_with_genomewide_rows": int((manifest_df["n_genomewide_rows"] > 0).sum()),
        "genomewide_rows": int(len(sig_df)),
        "mr_run": False,
        "interpretation": "v3 has no server-side p-value filter; rows were locally filtered at P<=5e-8. MR remains gated until instruments are deduplicated, cis-window checked, LD-clumped and harmonised against a prespecified AD outcome.",
    }
    (out / "mr_gate_v3_QA.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (out / "mr_gate_v3_QA.md").write_text(
        f"# Current-panel MR gate audit via eQTL Catalogue v3 ({summary['generated']})\n\n"
        f"- Probed {summary['genes_probed']} mapped genes across {summary['datasets_probed']} fixed datasets ({summary['queries']} gene-dataset queries).\n"
        f"- Successful queries: {summary['successful_queries']}; empty: {summary['empty_queries']}; client/upstream errors: {summary['error_queries']}.\n"
        f"- Queries with locally filtered genome-wide rows: {summary['queries_with_genomewide_rows']}; rows: {summary['genomewide_rows']}.\n"
        "- v3 does not implement a p-value query filter. The raw JSONL is retained and filtering uses the returned `pvalue` field at P<=5e-8.\n"
        "- No outcome harmonisation or MR estimate was run in this audit. Any rows still require cis-window, duplicate-trait, LD-clumping and allele-harmonisation checks.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print("output_dir", out)


if __name__ == "__main__":
    main()
