# Exact fiber elimination for the PCR3BP Levi-Civita convexity gates

Reproducibility sources and manuscript by César Augusto Grisa. Version **1.4.1**,
R07 editorial/reproducibility update, 27 October 2026. Version-specific Zenodo
identifier: https://doi.org/10.5281/zenodo.23512761. The Zenodo record shows
the registration and publication status of this identifier.

## Existing project and archival lineage

- Repository: https://github.com/cagrisa-creator/pcr3bp-fiber-elimination
- All-version Zenodo DOI: https://doi.org/10.5281/zenodo.23244861
- Earlier published version 1.4.0: https://doi.org/10.5281/zenodo.23244862

The version 1.4.0 DOI identifies the earlier archive. Version 1.4.1 uses the
existing record https://zenodo.org/records/23512761 in the same family; do not
create a separate record family or a duplicate draft.

## Current files

- `paper/main.tex` and `paper/main.pdf`: revised manuscript (14 pages).
- `paper/figures/`: exact assets for all three figures.
- `code/`, `audit/`, `results/`: revised implementations, exact checks,
  complete numerical samples, witness coordinates and scientific logs.
- `REPRODUCIBILITY.txt`: commands, dependencies, grids and result-to-file map.
- `CHANGELOG_AND_SCOPE.txt`: corrections and their numerical consequences.
- `frozen_inputs/`: historical code/data required for the comparison runs.
- `SHA256SUMS_R07.txt`: checksums for the listed current scientific files and metadata; obsolete Word/HAL variants are archived separately.

Only `paper/main.tex` and `paper/main.pdf` are current manuscript files.
The former Word/HAL mirrors and their deposit instructions are retained under
`historical/pre_R06/` for provenance only; they are obsolete and must not be
used as the current manuscript. No current Word or HAL mirror is supplied.
The obsolete output logs and their original environment/manifests are now
under `historical/pre_R07/`; `results/README_CURRENT.md` identifies the current
evidence. No old log has been rewritten to imply a fresh computation.

## Scientific status

The fiber identities and geometric reconstruction are exact. The numerical
reconnaissance consists of finite floating-point samples. No fixed parameter
pair or parameter region has an interval numerical convexity certificate.
An 80-digit witness check is a precision cross-check, not an interval enclosure.
The seven principal sections, numbered equations, central proofs and
bibliography are preserved. The R07 manuscript identifies Kim's LC result by
content and specifies the cited arXiv version for the nearby remark; it
changes no mathematical claim. Corrections include reflected Figure 2 domains,
the explicit geometric derivation, complete reproducibility data and precise
sampling qualifications. All 112 phase-map sign classifications are unchanged.

## Reproduce

Install `requirements.txt`, then follow `REPRODUCIBILITY.txt`.
The unchanged core exact audits, new 23-check reconstruction and packaged
diagnostic driver passed. All three figures were regenerated from this code.
Code/data files are identical to the audited snapshot except that one audit
summary replaces local machine paths with relative paths. Editorial letters,
referee reports and private continuity documents are excluded.

## License

MIT applies to the source code. The manuscript and figures remain copyright
the author, with all rights reserved pending publication; see `LICENSE`.
