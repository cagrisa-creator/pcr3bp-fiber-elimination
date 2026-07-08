> **VERSION 2 (2026-07-08):** corrected mass convention (m_r = 1 - mu) and
> scanning heuristic; Kim cross-check claim withdrawn; shell reattributed to
> the Moon component; Earth-component macroscopic window added. Deposit as a
> new HAL version of the same record. See ERRATUM.md in the repository.

# HAL deposit metadata — copy-paste ready

## Title

Exact fiber elimination for the Levi–Civita convexity gates of the planar restricted three-body problem, with a global convexity reconnaissance

## Author

César Augusto Grisa
Affiliation: Instituto de Computação, Universidade Federal Fluminense (UFF), Niterói, RJ, Brazil
Email: cagrisa@id.uff.br
(add idHAL / ORCID during deposit if you have one)

## Document type

Pré-print / Document de travail (working paper) — not yet peer-reviewed at
time of writing. Update to "Article" with a link to the published version
once accepted by a journal, subject to that journal's self-archiving policy.

## Language

English

## Abstract (plain text, no LaTeX markup)

Computer-assisted strict-convexity certification of the regularized energy
surfaces of the planar circular restricted three-body problem (PCR3BP)
reduces to the simultaneous positivity of two scalar gates, T and D, over a
four-dimensional domain: a two-dimensional Levi–Civita base times a fiber
angle phi in S^1. We prove four exact algebraic identities that eliminate the
fiber angle in closed form. The tangential gate satisfies min_phi T = calT
with calT an explicit function of the base alone, attained at a computable
angle; the curvature-determinant gate admits an exact five-term harmonic
decomposition whose second-harmonic amplitude is 128 F^2 s for the structure
matrices at hand, yielding the rigorous lower bound D >= calD := A0 -
sqrt(A1^2+B1^2) - 128 F^2 s, which is sharp on the collar. Both gates thereby
reduce from four variables to two, unconditionally. All identities are
confirmed by an exact symbolic audit whose residuals vanish identically. As
an application, the reduced gates reproduce the certified convexity bounds
near the Earth–Moon mass ratio. Using the reduced gates we compute a global
reconnaissance of Levi–Civita convexity across the mass–energy plane; a global reconnaissance parametrized by the regularized-primary mass
fraction m_r = 1 - mu: near the critical energy the non-convex window spans
m_r from ~0.045 to beyond 0.999, consistent with Kim's 16/17 sufficient
condition for the heavy component and showing it is far from sharp. We record a numerical observation that the
non-convexity region avoids the collision-orbit axes and forms lateral lobes
on which the failure is a single simple fold -- one negative principal
curvature, the tangential gate staying positive -- while the base remains
radially star-shaped; and we delimit the (provably open) comparable-mass
regime. A refinement of the reconnaissance down to energy offsets of 1e-7
reveals that even the Moon component of the Earth-Moon system (m_r~0.0122)
loses LC-convexity on an
ultrathin shell of width about 9.5e-6 above the critical energy, where the
Lagrange point enters the surface; this locates the natural boundary of the
LC route at the critical limit.

## Keywords

restricted three-body problem; Levi-Civita regularization; strict convexity;
global surface of section; Birkhoff conjecture; computer-assisted proof;
Hofer-Wysocki-Zehnder; dynamical convexity

## Scientific domain (HAL classification)

Primary: Mathematics / Dynamical Systems (math.DS)
Secondary (if HAL's list allows more than one): Mathematical Physics
(math-ph); Computer Science / Symbolic Computation

MSC 2020 classification (from the paper): Primary 70F07; Secondary 37J06,
53A07, 68W30

## License

CC-BY (recommended by the CCSD; select during the file-upload step of the
deposit form — mandatory since February 2026). See `docs/HAL.md` for context.

## Files to attach

- Primary document: `main.pdf` (this folder)
- Optional supplementary file: `source_bundle.zip` (this folder) — contains
  `main.tex` and `figures/`, self-contained and verified to compile standalone
  with `pdflatex main.tex && pdflatex main.tex`.

## Related identifiers (fill in once available)

- GitHub repository: `https://github.com/seu-usuario/seu-repositorio`
- Zenodo DOI: `https://doi.org/10.5281/zenodo.XXXXXXX`
- HAL identifier: assigned by HAL after moderation (`hal-XXXXXXXX`)
