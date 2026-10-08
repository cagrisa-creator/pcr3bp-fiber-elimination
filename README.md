# Exact fiber elimination for the PCR3BP Levi-Civita convexity gates

Reproducibility package and manuscript sources for

> **Exact fiber elimination for the Levi-Civita convexity gates of the planar
> restricted three-body problem, with a global convexity reconnaissance**
> C. A. Grisa (Universidade Federal Fluminense).

## Branch status

This branch contains **CMDA R04 staging metadata** for the author-review draft prepared after the minor-revision report.

It is **not** the final public release and **not** yet the exact package to archive on Zenodo. The R04 manuscript/PDF/figure binaries remain in the private R04 package until they can be committed with a normal git client or another binary-safe GitHub upload path.

## Why this branch exists

The available connector in this session can safely create and update UTF-8 text files, but it does not provide a binary-safe upload path for the revised PDF and PNG figure files. To avoid a half-synchronized repository, this branch records the R04 state, checksums, and next steps without replacing the old PDF/figures in place.

## R04 scope

- Keeps the paper focused on exact fiber elimination and numerical reconnaissance.
- Removes the unsupported companion-manuscript dependency from the working draft.
- Does not import material from the separate G459 manuscript.
- Addresses the referee's precision requests on dimension counting, the base/collar distinction, determinant lower bound versus true fiber minimum, numerical evidence versus interval certification, notation, bibliography, and identifiers.

## Files added in this branch

- `docs/REVISION_R04_NOTES.md` - author-review staging notes.
- `R04_UPLOAD_SHA256SUMS.txt` - checksums for the R04 files that must be synchronized.
- `patches/R04_main_tex.patch` - staging instructions for applying the R04 source package.

## Files that must still be synchronized before release

- `paper/main.tex`
- `paper/main.pdf`
- `paper/figures/phase_map.png`
- `paper/figures/dynamical_relevance.png`
- `paper/figures/lobe_structure.png`
- `CITATION.cff`
- any HAL/Zenodo metadata
- DOCX files, if required, regenerated from the final approved LaTeX source

## Final release workflow

1. Apply the R04 package locally with a normal git client.
2. Rebuild `paper/main.pdf` from `paper/main.tex`.
3. Confirm file checksums against `R04_UPLOAD_SHA256SUMS.txt` or regenerate the manifest after any intentional change.
4. Re-run the symbolic and numerical QA.
5. Create the final tag.
6. Reserve/create the Zenodo DOI from that exact tag.
7. Insert the DOI into the manuscript and metadata.
8. Rebuild, re-check, and tag the final release.

## Historical baseline

The historical public baseline remains v1.3.2. It should not be cited as containing the R04 revision until the final synchronized files are committed and tagged.
