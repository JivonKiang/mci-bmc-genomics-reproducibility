"""Check the public repository boundary without requiring research inputs."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REQUIRED_FILES = (
    "README.md",
    "CITATION.cff",
    "CITATION.md",
    "DATA_LICENSE.md",
    "LICENSE",
    "CONTRIBUTING.md",
    "SECURITY.md",
    "docs/PUBLIC_RELEASE_SCOPE.md",
)
FORBIDDEN_SUFFIXES = {
    ".rda",
    ".rdata",
    ".rds",
    ".h5",
    ".h5ad",
    ".loom",
    ".docx",
    ".xlsx",
    ".zip",
}
FORBIDDEN_DATA_HEADERS = {
    "sample_id",
    "subject_id",
    "cell_id",
    "donor_id",
    "donor_metadata",
}


def repository_files() -> list[Path]:
    ignored_directories = {".git", ".venv", "__pycache__"}
    return [
        path
        for path in ROOT.rglob("*")
        if path.is_file()
        and not any(part in ignored_directories for part in path.parts)
        and path.suffix.lower() != ".pyc"
    ]


def main() -> int:
    errors: list[str] = []

    for relative in REQUIRED_FILES:
        if not (ROOT / relative).is_file():
            errors.append(f"missing required file: {relative}")

    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    if "license: MIT" not in citation:
        errors.append("CITATION.cff must declare license: MIT")
    if "repository-code:" not in citation:
        errors.append("CITATION.cff must declare repository-code")

    for path in repository_files():
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            errors.append(f"forbidden file type in public release: {path.relative_to(ROOT)}")

    aggregate_root = ROOT / "data" / "aggregate"
    for path in aggregate_root.glob("*.csv"):
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            header = next(csv.reader(handle), [])
        normalized = {column.strip().lower() for column in header}
        leaked = sorted(normalized & FORBIDDEN_DATA_HEADERS)
        if leaked:
            errors.append(f"forbidden aggregate header in {path.relative_to(ROOT)}: {', '.join(leaked)}")

    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1

    print(f"OK: public release checks passed for {len(repository_files())} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
