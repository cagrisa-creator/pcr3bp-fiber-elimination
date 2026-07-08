# Exact fiber elimination for the PCR3BP Levi–Civita convexity gates

Reproducibility package for the paper

> **Exact fiber elimination for the Levi–Civita convexity gates of the planar
> restricted three-body problem, with a global convexity reconnaissance**
> C. A. Grisa (Universidade Federal Fluminense).

The paper source and compiled PDF are in [`paper/`](paper/).

## Erratum (version 2 / v1.3.0)

Version 1 misread the mass parameter of the defining function (μ is the
**distant**-primary fraction; the regularized primary has m_r = 1 − μ) and used
a worst-points-only scanning heuristic that produced false "convex" entries in
the phase map. Both are corrected in v1.3.0; see [`ERRATUM.md`](ERRATUM.md) for
the full account and `audit/audit_artigo1.py` for an independent, self-contained
verification (symbolic + numerical). The fiber-elimination lemmas and their
symbolic audit are unaffected.

## What this is

Certifying that a regularized energy surface of the planar circular restricted
three-body problem (PCR3BP) bounds a strictly convex domain reduces, after the
Levi–Civita regularization, to the positivity of two scalar **gates** `T` and `D`
over a four-dimensional domain: a two-dimensional base `z = (z1, z2)` times a
fiber angle `φ ∈ S¹`. This repository proves and reproduces an exact
**elimination of the fiber angle**:

- the tangential gate has an exact closed-form fiber minimum, `min_φ T = 𝒯(z)`;
- the curvature-determinant gate admits an exact five-term harmonic expansion
  with second-harmonic amplitude `128 F² s`, giving the rigorous, collar-sharp
  lower bound `D ≥ 𝒟(z) = A₀ − √(A₁²+B₁²) − 128 F² s`.

Both gates thereby reduce from four variables to two, **unconditionally**. Every
identity is confirmed by an exact symbolic audit whose eight residuals vanish
identically.

The package also contains a floating-point **global reconnaissance** of
Levi–Civita convexity across the mass–energy plane (whose non-convexity boundary
is consistent with Kim's 16/17 sufficient condition and extends well beyond it on the heavy side), a **dynamical-relevance observation** that
the non-convexity region avoids the collision-orbit axes and, on the lateral
lobes where it lives, is a single simple fold on a base that stays radially
star-shaped, and a **refinement at
the critical limit** showing that even the Moon component (m_r≈0.0122) of the Earth–Moon pair loses
LC-convexity on an ultrathin shell of width ≈ 9.5 × 10⁻⁶ above the critical
energy — the natural boundary of the LC route.

## Repository layout

```
.
├── paper/
│   ├── main.tex               # LaTeX source (amsart) — canonical
│   ├── main.pdf               # compiled paper (10 pp.) — canonical
│   ├── main.docx              # Word version, for venues requiring .docx
│   ├── main_docx_source.tex   # intermediate source used to generate main.docx
│   ├── DOCX_NOTES.md          # why main_docx_source.tex exists; how to regenerate main.docx
│   └── figures/
│       ├── phase_map.png     # global convexity reconnaissance (Fig. 1)
│       ├── dynamical_relevance.png  # off-axis lobes + on-axis D profile (Fig. 2)
│       └── lobe_structure.png       # simple-fold structure of the lobes (Fig. 3)
├── code/
│   ├── fiber_elimination_audit.py  # exact symbolic proof of the 4 lemmas
│   ├── convexity_map.py            # reduced-gate engine; regenerates the map data
│   ├── dynamical_relevance.py      # localizes the non-convexity lobes; makes Fig. 2
│   ├── lobe_structure.py           # fine structure of the lobes: simple fold + star-shape (Sec. 6.1, Fig. 3)
│   ├── rotation_term_audit.py      # exact closed form of the elliptic rotation term (Sec. 7, eq. rotation)
│   ├── critical_shell.py           # ultrathin non-convex shell at the critical limit (Sec. 5.5)
│   └── run_all.py                  # reproducibility driver (fast; checks all claims)
├── results/
│   ├── phase_map_data.csv    # (μ, Δc) grid with calT_min, calD_bound_min, D_true_min
│   ├── critical_shell_output.txt  # reference output of critical_shell.py (slow module)
│   ├── lobe_structure_output.txt  # reference output of lobe_structure.py
│   ├── expected_output.txt   # reference output of run_all.py
│   └── SHA256SUMS            # checksums of the tracked artifacts
├── docs/
│   ├── ROADMAP.md            # what is proven, what is open, and the routes past the ceiling
│   ├── ZENODO.md             # step-by-step: archiving this repo on Zenodo for a DOI
│   └── HAL.md                # step-by-step: submitting to HAL
├── hal_submission/           # self-contained bundle ready to upload to HAL
├── requirements.txt
├── CITATION.cff
├── LICENSE
└── README.md
```

## Quick start

```bash
pip install -r requirements.txt

cd code
python run_all.py            # ~1 min: runs the audit + checks the reproducible claims
```

Expected tail:

```
[1/5] Exact symbolic audit ...
      -> PASS
[2/5] Critical-energy sanity + spot reconnaissance ...
      -> PASS
[3/5] Dynamical-relevance observation (on-axis check) ...
      -> PASS
[4/5] Fine structure of the non-convexity lobes ...
      -> PASS
[5/5] Rotation-term closed form (elliptic coordinates) ...
      -> PASS
ALL CHECKS PASSED
```

The full reference output is in [`results/expected_output.txt`](results/expected_output.txt).

### Individual steps

```bash
cd code
python fiber_elimination_audit.py   # prints the 8 residuals (all 0) and PASS
python convexity_map.py             # ~2 min: regenerates results/phase_map_data.csv
python dynamical_relevance.py       # regenerates paper/figures/dynamical_relevance.png
python lobe_structure.py            # simple-fold + star-shape table; regenerates Fig. 3 data
python rotation_term_audit.py       # exact elliptic rotation term; prints PASS_ROTATION_TERM_AUDIT
python critical_shell.py            # ~20 min: shell width dc* and robustness checks (Sec. 5.5)
```

`critical_shell.py` is the slow module; its reference output is committed as
[`results/critical_shell_output.txt`](results/critical_shell_output.txt) so the
numbers quoted in the paper can be checked without rerunning it.

### Building the paper

The paper is provided in three formats, all already built in this repository:
`main.tex` (LaTeX source), `main.pdf` (compiled PDF), and `main.docx` (Word).
`main.tex`/`main.pdf` are canonical; `main.docx` is a courtesy version for
venues or collaborators that need `.docx` (see
[`paper/DOCX_NOTES.md`](paper/DOCX_NOTES.md) for how it is generated and how to
regenerate it after editing `main.tex`).

To rebuild the PDF from source:

```bash
cd paper
pdflatex main.tex && pdflatex main.tex
```

### Before you submit: fill in the placeholders

`main.tex` (and `main.docx`) contain a **Data and Code Availability** statement
with two placeholders that must be filled in before submission:

```
The computational code developed for this study is publicly available on
GitHub at https://github.com/seu-usuario/seu-repositorio. The exact version
used to generate the results presented here has been archived and is
available on Zenodo under DOI https://doi.org/10.5281/zenodo.XXXXXXX.
```

Replace the GitHub URL with this repository's real URL once pushed, and the
Zenodo DOI once archived — see [`docs/ZENODO.md`](docs/ZENODO.md). After
editing `main.tex`, re-derive `main_docx_source.tex` and regenerate `main.docx`
per [`paper/DOCX_NOTES.md`](paper/DOCX_NOTES.md), and recompile the PDF.

### Submitting to HAL

See [`docs/HAL.md`](docs/HAL.md) and the self-contained bundle in
[`hal_submission/`](hal_submission/).

## What is rigorous, and what is not

This distinction is made explicit throughout the paper and in
[`docs/ROADMAP.md`](docs/ROADMAP.md):

- **Rigorous (exact symbolic identities):** the four fiber-elimination lemmas and
  the reduction theorem (`code/fiber_elimination_audit.py`, all residuals zero).
- **Reproduction of a certified computation:** the reduced gates recover the
  interval-certified strict-convexity bounds near the Earth–Moon mass ratio
  (companion certificate manuscript).
- **Floating-point reconnaissance (not interval-certified):** the global
  convexity map and the dynamical-relevance observation. These are exploratory
  and are labeled as such.

The comparable-mass, near-critical regime is **provably out of reach** of the
Levi–Civita-convexity route (Kim 2017); the two natural routes past this ceiling
are discussed in the paper and the roadmap.

## Conventions

`μ` is the (large-primary) mass parameter; `c` is one half the Jacobi constant;
`s = z1² + z2²`, `a = z1² − z2²` in Levi–Civita base coordinates. The defining
function `F`, the pencil `L`, and the structure matrices `M₁, M₂` are as in the
paper.

## Requirements

Python 3.9+ with `sympy`, `numpy`, `scipy`, `matplotlib` (see
[`requirements.txt`](requirements.txt)), and a TeX distribution with `pdflatex`
to build the paper.

## Citing

See [`CITATION.cff`](CITATION.cff).

## License

Code: MIT (see [`LICENSE`](LICENSE)). The paper text and figures in `paper/` are
© the author, all rights reserved pending publication.


## Tested environments and timings (v1.3.1)

Verified on Python 3.12.3 / SymPy 1.14.0 / NumPy 2.4.4 / SciPy 1.17.1:
`fiber_elimination_audit.py` ~1 s (v1.3.1 uses exact cheap zero-tests before
`simplify`; identical verdicts); `run_all.py` ~3-4 min (the reconnaissance
spot-checks dominate); `convexity_map.py` full map ~4 min at 192x80;
`audit/audit_artigo1.py` ~3-4 min; `audit/v130_validation.py` ~3 min.
SymPy timings vary significantly across versions; if a script stalls, check
the SymPy version first.

### v1.3.2: fast import via verified expression cache

convexity_map.py now ships a symbolic-expression cache (code/_gates_cache.py):
import costs ~5 s (sympify + lambdify) instead of the full derivation, on every
tested Python. The cache is verified against the from-scratch derivation by
audit/verify_gates_cache.py (six exact zero residuals: PASS_GATES_CACHE).
Deleting _gates_cache.py falls back to the full derivation automatically.
