# CMDA R04 author-review branch

This branch was created from `main` to stage the CMDA R04 author-review revision.

Status: **not final, not yet a Zenodo release, not yet ready for journal resubmission**.

## Scope

- The R04 revision removes the unsupported companion-manuscript dependency from the working draft.
- It keeps the separate G459 manuscript out of this article.
- It addresses the referee's precision requests on dimension counting, the base/collar distinction, determinant lower bound versus true fiber minimum, numerical evidence versus interval certification, notation, and bibliography.
- The real GitHub repository is already linked; a verified Zenodo DOI remains pending.

## Important limitation of this GitHub staging step

The available GitHub connector in this session can create and update UTF-8 text files directly. It does not provide a direct file-path upload operation for the R04 binary artifacts (`paper/main.pdf` and the revised PNG figures). Those binary artifacts remain in the local R04 package and must be uploaded through a full git client or another GitHub upload path before a public release or Zenodo archive is created.

Therefore, this branch is a staging marker and checklist, not the final synchronized R04 package.

## Files that still need synchronization before release

- `paper/main.tex` - R04 source.
- `paper/main.pdf` - compiled R04 manuscript.
- `paper/figures/phase_map.png` - R04 annotation-synchronized figure.
- `paper/figures/dynamical_relevance.png` - R04 annotation-synchronized figure.
- `README.md`, `CITATION.cff`, and any HAL/Zenodo metadata.
- Remove or regenerate stale DOCX files before any DOCX delivery.

## Next release workflow

1. Apply the R04 package with a normal git client.
2. Rebuild `paper/main.pdf` from `paper/main.tex`.
3. Confirm the figure files match the R04 package.
4. Re-run symbolic and numerical QA.
5. Create a final tag.
6. Reserve or create the Zenodo DOI from that exact tag.
7. Insert the DOI into the manuscript and metadata.
8. Rebuild and re-tag only after the DOI text is final.
