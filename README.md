# MCI blood transcriptomics reproducibility

Public code and aggregate result tables supporting the BMC Genomics submission prepared on 31 July 2026.

## Study scope

The study evaluates a de novo 95-gene MCI-centred blood expression programme and a locked 12-gene research score. It reports reproducibility within the AddNeuroMed development system, no observed reproduction in the independent GSE249477 blood RNA-seq cohort, and bounded progression, single-cell, MR and brain-context analyses.

The score is a research representation of a cross-sectional stage-associated state. This repository does not support a claim of a validated clinical biomarker, an MCI-to-AD prediction model, causal immune mediation or therapeutic target engagement.

## Contents

- `code/`: selected analysis, audit and figure-support scripts retained from the versioned project workflow.
- `data/aggregate/`: aggregate gene-, contrast-, pathway-, perturbation- and audit-level tables used for figures and evidence checks.
- `docs/`: code/data/figure mapping, data availability and public-release scope.

The scripts are provenance-preserving research code. Several scripts reference the original local project layout and may require path adaptation before independent execution. The repository is not presented as a one-command software package.

## Public data sources

The analyses use public GEO series: GSE63060, GSE63061, GSE249477, GSE282742, GSE136243, GSE150693, GSE285831 and GSE134578.

GEO is the primary source for the public datasets:

- https://www.ncbi.nlm.nih.gov/geo/
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285831
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134578

The repository does not redistribute controlled-access ADNI GO/2 data.

## Data protection boundary

This public release contains no sample-level, subject-level or cell-level records; no `sample_id`, `subject_id`, `cell_id` or donor-level metadata; no UMAP metadata; no RDA/RDS/H5/H5AD objects; and no controlled-access data. The aggregate tables retain only the summaries needed to inspect the reported analyses.

## Citation

See `CITATION.cff` and `docs/Data_and_Code_Availability_20260731.md`.

## Contact

Public profile: https://github.com/JivonKiang
