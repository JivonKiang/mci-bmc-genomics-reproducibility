# Public release scope

Prepared: 3 August 2026

## Included

- Selected Python, R and validation scripts from the BMC Genomics analysis workflow.
- Aggregate source tables for gene programmes, locked-score contrasts, pathway summaries, perturbation summaries, MR summaries and evidence-unit audits.
- Relative-path documentation linking code, aggregate tables and manuscript figures.
- Verified source references for the added GSE285831 brain-context and GSE134578 single-cell context layers.

## Excluded

The following classes were intentionally excluded from this public repository:

- `sample_id`, `subject_id`, `cell_id` and donor-level records.
- UMAP coordinates and cell-level expression or QC metadata.
- RDA, RDS, H5, H5AD, loom and other serialized analysis objects.
- Controlled-access ADNI GO/2 data.
- Manuscript DOCX files, supplementary ZIP archives and private correspondence.
- Temporary logs, Python caches and historical build artifacts not needed for the public code release.

The exclusion boundary is conservative. The public repository contains aggregate outputs only and should not be treated as a replacement for the original GEO records.
