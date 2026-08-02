# Code, data and figure mapping

Prepared: 3 August 2026

This mapping describes the public GitHub release for the BMC Genomics manuscript. The repository contains selected code and aggregate tables only; the complete audit package remains local to the submission workspace.

| Figure | Scientific role | Public aggregate source data | Public code location |
|---|---|---|---|
| Figure 1 | Study design, estimands and evidence boundary | `data/aggregate/evidence_unit_audit_20260730_114606.csv` | `code/redraw_figure1_evidence_architecture_20260730_1455.py` |
| Figure 2 | De novo 95-gene programme and locked 12-gene panel | `data/aggregate/source_figure2_95_gene_programme.csv`; `source_figure2_locked_12_gene_panel.csv` | `code/run_stage_discovery_gse63060_20260721.py`; `code/build_mci_centered_scores_20260722.R` |
| Figure 3 | Locked score, development reproducibility and GSE249477 test | `data/aggregate/source_figure3_discovery_only_compression.csv`; `source_figure3_locked_score_contrasts.csv` | `code/run_stage_replication_gse63061_20260721.py`; `code/run_external_stage_validation_gse249477_20260721.py` |
| Figure 4 | Exploratory progression and orthogonal context | `data/aggregate/source_figure4_progression_top100.csv`; `source_figure5_brain_context_effects.csv`; `source_figure5_temra_context_effects.csv`; `source_figure5_spi1_mr.csv` | `code/run_progression_concordance_gse282742_20260721.py`; `code/audit_gse282742_mapping_reconciliation_20260728_v2.R`; `code/run_current_panel_brain_temra_context_20260723.py` |
| Figure 5 | Single-cell immune context and perturbation gate | `data/aggregate/source_figure6_current_panel_virtual_perturbation_summary.csv`; `source_figure6_cell_proportion_MCI_vs_HC.csv`; `source_scRNA_current_panel_virtual_perturbation_summary.csv` | `code/run_mci_locked_panel_scTenifoldKnk_20260729.R` |
| Figure 6 | Integrated QC, communication and TF-activity context | `data/aggregate/source_figure6_cellchat_*`; `source_scRNA_cellchat_*`; `source_figure6_tf_regulon_audit.csv`; `source_scRNA_TF_regulon_audit.csv` | `code/generate_tables.R`; `code/test_locked_mci_panel_20260722.R` |

GSE285831 is a frontal-cortex brain expression context layer and GSE134578 is a donor-aware CSF single-cell context layer. Neither is an independent bulk-blood validation cohort. Their source references are recorded in `docs/REFERENCE_SOURCES_20260803.md`.

Sample-level, donor-level and cell-level tables are deliberately omitted from this public release.
