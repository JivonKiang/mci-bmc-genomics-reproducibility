# Analysis code

The files in this directory are selected, versioned scripts from the MCI project audit and manuscript workflow. They cover discovery, development-system stability, the GSE249477 external test, progression mapping, single-cell context, MR sensitivity analysis, figure support and validation.

The scripts were developed against a local project workspace. Set `MCI_ROOT`, `MCI_PROJECT_ROOT`, `MN_SUBMISSION_PACKAGE` or `EQTL_CATALOGUE_REQUEST_CLIENT` as appropriate for the script before running it. The files are retained here as provenance-linked code rather than as a guaranteed standalone package. Use the relative source-data mapping in `../docs/Code_Data_Figure_Mapping_20260803.md` when adapting a script.

No credentials, controlled-access data or individual-level data are required by this public release. The two context helper modules used by `run_current_panel_brain_temra_context_20260723.py` are included in this release; raw inputs remain outside the repository and must be supplied through `MCI_ROOT`.
