# MCI blood transcriptomics reproducibility

Public research code and aggregate result tables supporting a BMC Genomics submission prepared on 3 August 2026.

[![Public release checks](https://github.com/JivonKiang/mci-bmc-genomics-reproducibility/actions/workflows/public-release-check.yml/badge.svg)](https://github.com/JivonKiang/mci-bmc-genomics-reproducibility/actions/workflows/public-release-check.yml)
[![License: MIT](https://img.shields.io/badge/Code%20license-MIT-blue.svg)](LICENSE)

## Release status

This is a provenance-preserving research release, not a clinical software package or a one-command reproduction of the complete local project. The public repository contains selected code and aggregate outputs; the original raw, controlled-access and sample-level inputs remain outside the repository.

## Study scope

The study evaluates a de novo 95-gene MCI-centred blood expression programme and a locked 12-gene research score. It reports reproducibility within the AddNeuroMed development system, no observed reproduction in the independent GSE249477 blood RNA-seq cohort, and bounded progression, single-cell, MR and brain-context analyses.

The score is a research representation of a cross-sectional stage-associated state. This repository does not support a claim of a validated clinical biomarker, an MCI-to-AD prediction model, causal immune mediation or therapeutic target engagement.

## Contents

- `code/`: selected Python and R analysis, audit and figure-support scripts retained from the versioned project workflow.
- `data/aggregate/`: aggregate gene-, contrast-, pathway-, perturbation- and audit-level tables used for figures and evidence checks.
- `docs/`: code/data/figure mapping, data availability, source references and public-release scope.
- `scripts/`: checks for the public release boundary and repository metadata.

The scripts retain provenance links to the original project. Several scripts reference the original local project layout and require path adaptation or excluded inputs before independent execution. See `code/README.md` and the mapping in `docs/Code_Data_Figure_Mapping_20260803.md`.

## Public-only validation

The repository can be checked without access to the excluded research inputs:

```text
python scripts/check_public_release.py
python -m compileall -q code
```

The GitHub Actions workflow runs these checks on pushes and pull requests. Passing these checks confirms repository structure and Python syntax only; it does not reproduce the scientific analyses.

## Reproducible environment

The Python scripts require Python 3.10 or newer and the packages listed in `requirements.txt`. Scripts that read local raw data require `MCI_ROOT` to point to the full local project; the public repository intentionally does not contain those raw or controlled-access files. R workflows require R 4.2 or newer and the packages listed in `R_environment.md`.

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

Please cite the repository when reusing its code or aggregate result tables. The preferred citation and the current citation status are documented in [`CITATION.md`](CITATION.md) and [`CITATION.cff`](CITATION.cff). The associated manuscript did not have a DOI at this release; the citation file should be updated when the final article is published.

Please do not describe this repository as a validated clinical biomarker or diagnostic package. See [`docs/Data_and_Code_Availability_20260803.md`](docs/Data_and_Code_Availability_20260803.md) and [`docs/REFERENCE_SOURCES_20260803.md`](docs/REFERENCE_SOURCES_20260803.md) for the release and source boundaries.

## License and reuse

The MIT License in [`LICENSE`](LICENSE) applies to original source code in this repository. The aggregate files under `data/aggregate/` are not automatically covered by the code license; see [`DATA_LICENSE.md`](DATA_LICENSE.md) and the terms of the original data sources before reusing or redistributing them.

## Contributing and data handling

See [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening an issue or pull request. Do not upload sample-level, subject-level, donor-level or controlled-access data, manuscripts, credentials or private correspondence to this repository. For sensitive reports, see [`SECURITY.md`](SECURITY.md).

## Contact

Public profile: https://github.com/JivonKiang
