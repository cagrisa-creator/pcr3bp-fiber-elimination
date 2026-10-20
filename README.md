# Exact fiber elimination for the PCR3BP Levi-Civita convexity gates

Reproducibility sources and manuscript by César Augusto Grisa. Release candidate **1.4.1**, prepared
20 October 2026. These sources update the existing GitHub repository. A version 1.4.1 Zenodo archive and its new DOI are still pending.

## Existing project and archival lineage

- Repository: https://github.com/cagrisa-creator/pcr3bp-fiber-elimination
- All-version Zenodo DOI: https://doi.org/10.5281/zenodo.23244861
- Earlier published version 1.4.0: https://doi.org/10.5281/zenodo.23244862

The earlier version DOI must not be used to identify this revised candidate.
A new version must be created within that same Zenodo record family.

## Current files

- `paper/main.tex` and `paper/main.pdf`: revised manuscript (14 pages).
- `paper/figures/`: exact assets for all three figures.
- `code/`, `audit/`, `results/`: revised implementations, exact checks,
  complete numerical samples, witness coordinates and scientific logs.
- `REPRODUCIBILITY.txt`: commands, dependencies, grids and result-to-file map.
- `CHANGELOG_AND_SCOPE.txt`: corrections and their numerical consequences.
- `frozen_inputs/`: historical code/data required for the comparison runs.
- `SHA256SUMS_R06.txt`: complete manifest for this public archive.

The current paper and its scientific data take precedence over legacy
Word/HAL/deposit instructions retained in older Git history or directories.
Those older files are not updated manuscript variants.

## Scientific status

The fiber identities and geometric reconstruction are exact. The numerical
reconnaissance consists of finite floating-point samples. No fixed parameter
pair or parameter region has an interval numerical convexity certificate.
An 80-digit witness check is a precision cross-check, not an interval enclosure.
The seven principal sections, numbered equations, central proofs and
bibliography are preserved. Corrections include reflected Figure 2 domains,
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
