# Exact fiber elimination for the PCR3BP Levi-Civita convexity gates

Reproducibility package and manuscript sources for

> **Exact fiber elimination for the Levi-Civita convexity gates of the planar
> restricted three-body problem, with a global convexity reconnaissance**
> C. A. Grisa (Universidade Federal Fluminense).

## DOI

Zenodo DOI for this revised package:

- https://doi.org/10.5281/zenodo.23244862

GitHub repository:

- https://github.com/cagrisa-creator/pcr3bp-fiber-elimination

## Revision scope

- Keeps the paper focused on exact fiber elimination and numerical reconnaissance.
- Removes the unsupported companion-manuscript dependency from the manuscript draft.
- Does not add material from outside the submitted-manuscript and referee-response scope.
- Addresses the referee's precision requests on dimension counting, the base/collar distinction, determinant lower bound versus true fiber minimum, numerical evidence versus interval certification, notation, bibliography, and identifiers.

## What is rigorous, and what is reconnaissance

- **Rigorous:** the exact algebraic fiber-elimination identities and the symbolic audit.
- **Rigorous sufficient bound:** the determinant gate is treated by an exact harmonic decomposition and a rigorous two-variable lower bound.
- **Numerical spot check:** the Earth-Moon values in the revised manuscript are pointwise numerical checks, not interval certificates.
- **Floating-point reconnaissance:** the global mass-energy map and dynamical-relevance observations are exploratory and are labeled as such.

## Core files

- `paper/main.tex` - manuscript source.
- `paper/main.pdf` - compiled manuscript.
- `paper/figures/phase_map.png` and `paper/figures/dynamical_relevance.png` - revised figures.
- `code/fiber_elimination_audit.py` and `code/rotation_term_audit.py` - symbolic audit scripts.
- `CITATION.cff` - citation metadata.

## Notes

The exact algebraic identities and symbolic audits are the rigorous part of the package. The Earth-Moon numerical values, global mass-energy map, and dynamical-relevance observations are numerical checks/reconnaissance unless explicitly stated otherwise in the manuscript.
