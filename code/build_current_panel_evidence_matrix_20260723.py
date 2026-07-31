"""Build the auditable evidence matrix for the development-locked MCI panel."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("MCI_PROJECT_ROOT", "."))
STAMP = "20260723_091200"
OUT = ROOT / f"mci_evidence_matrix_{STAMP}"
OUT.mkdir(parents=True, exist_ok=True)

PANEL = ROOT / "mci_development_optimized_panel_20260722_161919" / "development_locked_12_gene_panel.csv"
DEV = ROOT / "mci_development_optimized_panel_20260722_161919" / "development_locked_12_contrasts.csv"
EXT = ROOT / "mci_key_analyses_20260723_091500" / "mapping_weight_sensitivity_contrasts.csv"
PROG_SCORE = ROOT / "mci_key_analyses_20260723_091500" / "GSE282742_locked_panel_progression_contrast.csv"
PROG_PAIR = ROOT / "mci_key_analyses_20260723_091500" / "GSE282742_paired_transition_summary.csv"
PROG_GENE = ROOT / "mci_key_analyses_20260723_091500" / "GSE282742_locked_panel_gene_level_progression_audit.csv"
PROG_MAP = ROOT / "mci_key_analyses_20260723_091500" / "GSE282742_panel_mapping_used.csv"
BRAIN = ROOT / "mci_current_panel_context_20260723_090020" / "brain_current_panel_gene_effects.csv"
BRAIN_MAP = ROOT / "mci_current_panel_context_20260723_090020" / "brain_current_panel_mapping_used.csv"
TEMRA = ROOT / "mci_current_panel_context_20260723_090020" / "temra_current_panel_gene_effects.csv"
TEMRA_MAP = ROOT / "mci_current_panel_context_20260723_090020" / "temra_current_panel_mapping_used.csv"
LONG = ROOT / "gse136243_current_panel_audit_20260723" / "GSE136243_current_panel_slope_summary.csv"
MR = ROOT / "audit_20260721_081827" / "MR_instrument_manifest.csv"

# Compact result from the UniProt REST `xref_pdb` query run on 2026-07-23.
# Empty PDB lists are retained as negative targetability evidence.
UNIPROT = {
    "GPI": ("P06744", "G6PI_HUMAN", "1IAT;1IRI;1JIQ;1JLH;1NUH;6XUH;6XUI;8BBH;8P2K;9FCW;9FHF;9FKC;9FKF"),
    "PHF15": ("Q9NQC1", "JADE2_HUMAN", ""),
    "VAT1": ("Q99536", "VAT1_HUMAN", "6K9Y;6LHR;6LII"),
    "EIF3B": ("P55884", "EIF3B_HUMAN", "2KRB;2NLW;5K1H;6YBT;6ZMW;6ZON;6ZP4;6ZVJ;7A09;7QP6;7QP7;8OZ0;8PJ1;8PJ2;8PJ3;8PJ4;8PJ5;8PJ6;8XXN;9BLN;9CPA"),
    "UCP2": ("P55851", "UCP2_HUMAN", ""),
    "C22ORF39": ("Q6P5X5", "CV039_HUMAN", ""),
    "USP38": ("Q8NB14", "UBP38_HUMAN", "4RXX"),
    "ZSWIM6": ("Q9HCJ5", "ZSWM6_HUMAN", ""),
    "PRDM10": ("Q9NQV6", "PRD10_HUMAN", "3IHX"),
    "RP2": ("O75695", "XRP2_HUMAN", "2BX6;3BH6;3BH7"),
    "GCA": ("P28676", "GRAN_HUMAN", "1F4O;1F4Q;1K94;1K95"),
}


def one_effect(table: pd.DataFrame, gene: str, contrast: str) -> dict[str, float]:
    d = table[(table["gene"].astype(str).str.upper() == gene.upper()) & table["contrast"].eq(contrast)]
    if d.empty:
        return {"estimate": np.nan, "p_value": np.nan, "fdr": np.nan}
    row = d.iloc[0]
    return {"estimate": float(row.get("estimate", row.get("estimate_P_MCI_vs_S_MCI", np.nan))), "p_value": float(row.get("p_value", np.nan)), "fdr": float(row.get("fdr", np.nan))}


def sign_ok(value: float, expected: int) -> bool:
    return bool(np.isfinite(value) and np.sign(value) == expected)


def main() -> None:
    panel = pd.read_csv(PANEL)
    dev = pd.read_csv(DEV)
    ext = pd.read_csv(EXT)
    prog_score = pd.read_csv(PROG_SCORE)
    prog_pair = pd.read_csv(PROG_PAIR)
    prog_gene_raw = pd.read_csv(PROG_GENE)
    prog_map = pd.read_csv(PROG_MAP)
    brain = pd.read_csv(BRAIN)
    brain_map = pd.read_csv(BRAIN_MAP)
    temra = pd.read_csv(TEMRA)
    temra_map = pd.read_csv(TEMRA_MAP)
    longitudinal = pd.read_csv(LONG).iloc[0]
    mr = pd.read_csv(MR)

    # GSE282742 gene table has Ensembl feature IDs; attach frozen symbols.
    prog_gene = prog_gene_raw.merge(prog_map[["gene", "trajectory_class", "ensembl_id"]], left_on="feature_id", right_on="ensembl_id", how="left")
    rows = []
    for _, p in panel.iterrows():
        gene = str(p["gene"])
        role = str(p["trajectory_class"])
        expected_mci_hc = 1 if role == "MCI_peak" else -1
        expected_ad_mci = -1 if role == "MCI_peak" else 1
        d60 = dev[(dev.dataset == "GSE63060") & (dev.contrast == "MCI_vs_HC")].iloc[0]
        d61 = dev[(dev.dataset == "GSE63061") & (dev.contrast == "MCI_vs_HC")].iloc[0]
        d60_adj = dev[(dev.dataset == "GSE63060") & (dev.contrast == "AD_vs_MCI")].iloc[0]
        d61_adj = dev[(dev.dataset == "GSE63061") & (dev.contrast == "AD_vs_MCI")].iloc[0]
        b_mci = one_effect(brain, gene, "MCI_vs_HC")
        b_ad = one_effect(brain, gene, "AD_vs_MCI")
        b_trend = one_effect(brain, gene, "ordinal_HC_MCI_AD")
        t_mci = one_effect(temra, gene, "MCI_vs_HC")
        t_ad = one_effect(temra, gene, "AD_vs_MCI")
        t_trend = one_effect(temra, gene, "ordinal_HC_MCI_AD")
        gp = prog_gene[prog_gene["gene"].astype(str).str.upper().eq(gene.upper())]
        gp_row = gp.iloc[0] if not gp.empty else pd.Series(dtype=float)
        # External and longitudinal results are score-level evidence, so repeat
        # their estimates on each row and explicitly label them as global.
        e = ext[(ext.dataset == "GSE249477") & ext.contrast.isin(["MCI_vs_HC", "AD_vs_MCI"]) & ext.scoring_mode.eq("component_renorm")]
        e_mci = e[e.contrast.eq("MCI_vs_HC")].iloc[0]
        e_ad = e[e.contrast.eq("AD_vs_MCI")].iloc[0]
        ps = prog_score.iloc[0]
        prot = UNIPROT.get(gene, ("", "", ""))
        mr_rows = mr[mr["gene"].astype(str).str.upper().eq(gene.upper())]
        brain_mci_ok = sign_ok(b_mci["estimate"], expected_mci_hc)
        brain_ad_ok = sign_ok(b_ad["estimate"], expected_ad_mci)
        temra_mci_ok = sign_ok(t_mci["estimate"], expected_mci_hc)
        temra_ad_ok = sign_ok(t_ad["estimate"], expected_ad_mci)
        context_count = int(sum([brain_mci_ok, brain_ad_ok, temra_mci_ok, temra_ad_ok]))
        brain_mapped = gene in set(brain_map["gene"])
        temra_mapped = gene in set(temra_map["gene"])
        structure_ids = prot[2]
        # The causal gate is intentionally strict: no current-panel gene is
        # allowed through without a current, gene-specific MR instrument row.
        causal_gate = bool(mr_rows.shape[0] > 0)
        context_gate = bool((brain_mapped or temra_mapped) and context_count >= 2)
        rows.append({
            "gene": gene,
            "trajectory_class": role,
            "expected_MCI_vs_HC_sign": expected_mci_hc,
            "expected_AD_vs_MCI_sign": expected_ad_mci,
            "development_GSE63060_MCI_vs_HC": float(d60[f"estimate"] if False else p["GSE63060_estimate_MCI_vs_HC"]),
            "development_GSE63061_MCI_vs_HC": float(p["GSE63061_estimate_MCI_vs_HC"]),
            "development_direction_concordant": bool(p["development_consistent"]),
            "brain_mapped": brain_mapped,
            "brain_MCI_vs_HC_estimate": b_mci["estimate"],
            "brain_MCI_vs_HC_p": b_mci["p_value"],
            "brain_AD_vs_MCI_estimate": b_ad["estimate"],
            "brain_AD_vs_MCI_p": b_ad["p_value"],
            "brain_ordinal_estimate": b_trend["estimate"],
            "brain_MCI_direction_ok": brain_mci_ok,
            "brain_AD_direction_ok": brain_ad_ok,
            "temra_mapped": temra_mapped,
            "temra_MCI_vs_HC_estimate": t_mci["estimate"],
            "temra_MCI_vs_HC_p": t_mci["p_value"],
            "temra_AD_vs_MCI_estimate": t_ad["estimate"],
            "temra_AD_vs_MCI_p": t_ad["p_value"],
            "temra_ordinal_estimate": t_trend["estimate"],
            "temra_MCI_direction_ok": temra_mci_ok,
            "temra_AD_direction_ok": temra_ad_ok,
            "orthogonal_context_direction_count": context_count,
            "GSE249477_score_MCI_vs_HC_estimate_global": float(e_mci["estimate"]),
            "GSE249477_score_MCI_vs_HC_p_global": float(e_mci["p_value_HC3"]),
            "GSE249477_score_AD_vs_MCI_estimate_global": float(e_ad["estimate"]),
            "GSE249477_score_AD_vs_MCI_p_global": float(e_ad["p_value_HC3"]),
            "GSE249477_transportability_interpretation": "direction_reversed_vs_development",
            "GSE282742_score_P_MCI_vs_S_MCI_estimate_global": float(ps["estimate"]),
            "GSE282742_score_P_MCI_vs_S_MCI_p_global": float(ps["p_value_HC3"]),
            "GSE282742_gene_P_MCI_vs_S_MCI_estimate": float(gp_row.get("estimate_P_MCI_vs_S_MCI", np.nan)) if not gp.empty else np.nan,
            "GSE282742_gene_P_MCI_vs_S_MCI_fdr": float(gp_row.get("fdr", np.nan)) if not gp.empty else np.nan,
            "GSE136243_C_vs_N_slope_difference_p_global": float(longitudinal["welch_p"]),
            "GSE136243_interpretation": "mixed normal-to-MCI/AD exRNA longitudinal context",
            "MR_current_panel_instrument_rows": int(len(mr_rows)),
            "MR_current_panel_status": "not_in_audited_manifest" if len(mr_rows) == 0 else "present_requires_reanalysis",
            "uniprot_accession": prot[0],
            "uniprot_id": prot[1],
            "PDB_ids_from_UniProt": structure_ids,
            "experimental_structure_available": bool(structure_ids),
            "context_gate_pass": context_gate,
            "causal_MR_gate_pass": causal_gate,
            "full_priority_gate_pass": bool(context_gate and causal_gate),
            "priority_status": "mechanistic_shortlist_pending_current_MR" if context_gate and structure_ids and not causal_gate else ("context_supported_but_MR_missing" if context_gate else "context_inconclusive"),
        })
    matrix = pd.DataFrame(rows)
    matrix.to_csv(OUT / "current_panel_evidence_matrix.csv", index=False)
    matrix.sort_values(["full_priority_gate_pass", "context_gate_pass", "orthogonal_context_direction_count", "experimental_structure_available"], ascending=[False, False, False, False]).to_csv(OUT / "current_panel_evidence_matrix_ranked.csv", index=False)
    structure = matrix[["gene", "uniprot_accession", "uniprot_id", "PDB_ids_from_UniProt", "experimental_structure_available"]]
    structure.to_csv(OUT / "current_panel_structure_targetability_audit.csv", index=False)
    summary = {
        "generated": "2026-07-23 09:12 Asia/Shanghai",
        "panel_size": int(len(matrix)),
        "brain_mapped": int(matrix["brain_mapped"].sum()),
        "temra_mapped": int(matrix["temra_mapped"].sum()),
        "context_gate_pass": int(matrix["context_gate_pass"].sum()),
        "causal_MR_gate_pass": int(matrix["causal_MR_gate_pass"].sum()),
        "full_priority_gate_pass": int(matrix["full_priority_gate_pass"].sum()),
        "structure_available": int(matrix["experimental_structure_available"].sum()),
        "interpretation": "No current-panel gene passes the full priority gate because the audited MR instrument manifest contains no current-panel gene. Mechanistic shortlist is therefore provisional and cannot justify docking/MD as a causal result.",
        "progression_score_context": {
            "estimate": float(ps["estimate"]),
            "p": float(ps["p_value_HC3"]),
            "n": int(ps["n"]),
            "paired_n": int(prog_pair.iloc[0]["n_pairs"]),
            "paired_p": float(prog_pair.iloc[0]["paired_p"]),
        },
    }
    (OUT / "evidence_matrix_QA.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    lines = [
        f"# Frozen MCI panel evidence matrix ({summary['generated']})",
        "",
        "## Gate definition",
        "- Context gate: mapped in brain or TEMRA and concordant with the frozen peak/trough direction in at least two of four orthogonal contrasts (brain/TEMRA MCI-vs-HC and AD-vs-MCI).",
        "- Causal MR gate: current gene must have a gene-specific instrument row in the audited MR manifest; old TOMM7/RPS24/RPS27L rows are not inherited.",
        "- Full priority gate = context gate AND causal MR gate. Structure availability is a feasibility flag, not proof of mechanism.",
        "",
        f"## Current result: {summary['context_gate_pass']} context-supported genes; {summary['causal_MR_gate_pass']} current-MR genes; {summary['full_priority_gate_pass']} full-gate genes.",
        "The absence of current-panel MR rows is a data/analysis gap, not a negative MR finding. It means the next executable task is a fresh cis-eQTL plus AD/clinical outcome MR audit for the frozen panel.",
        "",
        "## Boundary",
        "GSE249477 reverses the development score direction and therefore blocks a validated blood biomarker claim. GSE282742 and GSE136243 remain progression/longitudinal context, not complete independent MCI-to-AD mRNA validation. PDB availability only determines whether structural work is technically possible after a target clears the evidence gate.",
    ]
    (OUT / "evidence_matrix_QA_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("output_dir", OUT)


if __name__ == "__main__":
    main()
