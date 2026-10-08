# Exact fiber elimination for the PCR3BP Levi-Civita convexity gates

Reproducibility package and manuscript sources for

> **Exact fiber elimination for the Levi-Civita convexity gates of the planar
> restricted three-body problem, with a global convexity reconnaissance**
> C. A. Grisa (Universidade Federal Fluminense).

## Branch status

This branch is the **CMDA R04 author-review staging branch** prepared after the minor-revision report.

It now contains the R04 manuscript source, compiled PDF, and revised figure files that were uploaded through GitHub Desktop:

- `paper/main.tex`
- `paper/main.pdf`
- `paper/figures/phase_map.png`
- `paper/figures/dynamical_relevance.png`

It is still **not** the final public release and **not** yet the exact package to archive on Zenodo, because the Zenodo DOI and final deposit metadata have not yet been reserved/inserted.

## R04 scope

- Keeps the paper focused on exact fiber elimination and numerical reconnaissance.
- Removes the unsupported companion-manuscript dependency from the manuscript draft.
- Does not add material from outside the submitted-manuscript and referee-response scope.
- Addresses the referee's precision requests on dimension counting, the base/collar distinction, determinant lower bound versus true fiber minimum, numerical evidence versus interval certification, notation, bibliography, and identifiers.

## What is rigorous, and what is reconnaissance

- **Rigorous:** the exact algebraic fiber-elimination identities and the symbolic audit.
- **Rigorous sufficient bound:** the determinant gate is treated by an exact harmonic decomposition and a rigorous two-variable lower bound.
- **Numerical spot check:** the Earth-Moon values in the revised manuscript are pointwise numerical checks, not interval certificates.
- **Floating-point reconnaissance:** the global mass-energy map and dynamical-relevance observations are exploratory and are labeled as such.

## Files added for the R04 review workflow

- `docs/REVISION_R04_NOTES.md` - author-review staging notes.
- `R04_UPLOAD_SHA256SUMS.txt` - checksums for the generated R04 package files.
- `patches/R04_main_tex.patch` - staging notes for applying/checking the R04 source package.
- `scripts/APPLY_R04_FROM_ZIP.sh` - optional local helper script.

## Still pending before final release

1. Inspect the PR files on GitHub, especially the PDF, TeX source, figures, README, and CITATION metadata.
2. Re-run local QA/reproducibility checks if a full local scientific release is desired.
3. Create a Zenodo draft and reserve a DOI manually, without publishing yet.
4. Insert the reserved DOI into `paper/main.tex`, `CITATION.cff`, and any deposit metadata.
5. Rebuild `paper/main.pdf` after DOI insertion.
6. Commit the DOI-synchronized files.
7. Only then mark the PR ready, merge/tag/release, and publish the Zenodo record.

## Current DOI status

No Zenodo DOI is claimed in this branch yet. The manuscript contains a revision-stage note saying that the verified Zenodo DOI and archival release are still pending.

## Historical baseline

The historical public baseline remains v1.3.2. This R04 branch is a review-stage update and should not be cited as a final archived release until the DOI-synchronized release is tagged and deposited.
