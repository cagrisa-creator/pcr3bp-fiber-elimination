# CMDA R04 author-review branch

This branch stages the CMDA R04 author-review revision.

Status: **review-stage branch; not final Zenodo release; not yet ready for journal resubmission**.

## Scope

- The R04 revision removes the unsupported companion-manuscript dependency from the manuscript draft.
- It keeps the separate G459 manuscript out of this article.
- It addresses the referee's precision requests on dimension counting, the base/collar distinction, determinant lower bound versus true fiber minimum, numerical evidence versus interval certification, notation, bibliography, and identifiers.
- The real GitHub repository URL is used.
- A verified Zenodo DOI remains pending and must be inserted before final release/resubmission.

## Uploaded R04 manuscript files

The following R04 files have been uploaded to this branch through GitHub Desktop:

- `paper/main.tex`
- `paper/main.pdf`
- `paper/figures/phase_map.png`
- `paper/figures/dynamical_relevance.png`

The metadata/checklist files were then updated through the GitHub connector:

- `README.md`
- `CITATION.cff`
- `FILL_BEFORE_DEPOSIT.md`
- `R04_UPLOAD_SHA256SUMS.txt`

## Important remaining limitation

This is still not a public archival release. The manuscript currently says that the verified Zenodo DOI and archival release are pending. That is intentional.

## Next release workflow

1. Inspect the GitHub PR files.
2. Reserve a Zenodo DOI manually, without publishing yet.
3. Insert the reserved DOI into `paper/main.tex`, `CITATION.cff`, and deposit metadata.
4. Rebuild `paper/main.pdf` from the DOI-synchronized TeX source.
5. Commit the DOI-synchronized files.
6. Re-check the final package.
7. Mark the PR ready only after the final package is complete.
8. Merge/tag/release and publish Zenodo only after DOI and metadata are synchronized.
