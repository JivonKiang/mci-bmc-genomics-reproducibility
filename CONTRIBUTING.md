# Contributing

Contributions are welcome when they improve reproducibility, documentation,
testing, source provenance or public-release safety.

## Before opening a pull request

- Do not add sample-level, subject-level, donor-level or controlled-access
  data.
- Do not add manuscripts, supplementary archives, credentials, private
  correspondence or local absolute paths.
- Add or update the relevant source citation and explain the provenance of new
  aggregate outputs.
- Keep research claims bounded by the study scope and the evidence available
  in the repository.
- Run `python scripts/check_public_release.py` and
  `python -m compileall -q code` when applicable.

Pull requests should explain what changed, which public inputs are required,
and whether the change affects any figure or aggregate result mapping. The
maintainer may request a reproducibility note or a checksum for new public
inputs.

By submitting code or documentation, you confirm that you have the right to
contribute it and that it may be distributed under the repository's applicable
license. Do not submit copied material, institutional-only code or third-party
data without confirming its reuse terms.
