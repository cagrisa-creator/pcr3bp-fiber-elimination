# Roadmap: what is proven, what is open

This file states plainly what this repository establishes rigorously, what it
establishes as a reproduction of prior certified work, what it only reconnoiters
numerically, and what remains genuinely open. It is meant to be read alongside
Section 7 ("Limitations and outlook") of the paper.

## Proven (exact symbolic identities)

The four fiber-elimination lemmas and the reduction theorem (paper, Section 3):

- **Lemma 1.** `min_φ T(z,φ) = 𝒯(z)`, an exact closed-form minimum, attained at
  an explicit angle.
- **Lemma 2.** `𝒯` depends on `z` only through `(s,a)`; hence the quadrant
  reduction (Corollary) is a consequence, not an assumption.
- **Lemma 3.** `D(z,φ)` has an exact five-term harmonic expansion in `φ`, for
  *arbitrary* symmetric 2×2 data — a pure multilinear identity, proved in full
  generality before specializing to the PCR3BP structure matrices.
- **Lemma 4.** The second-harmonic amplitude of the specialized expansion is
  exactly `128 F² s`.
- **Theorem 1.** Combining the above: both convexity gates reduce from the
  four variables `(z, φ)` to the two base variables `z`, unconditionally.

All eight residuals of the symbolic audit (`code/fiber_elimination_audit.py`)
vanish identically — this is checked by exact computer algebra (`sympy`), not
by floating-point sampling.

## Reproduction of an existing certified computation

Section 4 of the paper applies the reduced gates on the interval-certified
Earth–Moon parameter rectangle from the companion certificate manuscript
(rigorous, `T ≥ 7.63×10⁻⁴`, `D ≥ 1.15×10⁻⁷` via a `129,238`-subcell
adaptive interval-arithmetic subdivision) and confirms the reduced surrogates
reproduce those bounds.

## Floating-point reconnaissance (not interval-certified)

- **Section 5 (global map).** A `float64` scan of `(μ, Δc)` locating where the
  determinant gate is negative. The corrected reconnaissance (v2) is consistent with Kim's 16/17 sufficient condition and shows it is far from sharp: the heavy-side non-convexity persists to m_r = 0.999 near the critical energy.
- **Section 6 (dynamical-relevance observation).** The non-convexity region
  forms two off-axis lobes; the two collision-orbit axes have strictly
  positive `D` in the sampled cases. Evidence, not proof: the Hryniewicz–
  Salomão hypothesis concerns *all* closed Reeb orbits, not just the two
  collision orbits.
- **Section 6.1 (fine structure of the lobes, new).** Inside the lobes, the
  tangential gate `𝒯` stays strictly positive with a uniform margin across the
  whole comparable-mass band (`μ ∈ [0.2, 0.8]`) — so the failure of convexity
  is always a *single simple fold* (`tr > 0`, `det < 0`, one negative
  principal curvature), never a double concavity — and the base stays
  radially star-shaped throughout (`∂F/∂r < 0` uniformly, including on rays
  through the lobes). Still floating-point reconnaissance, not a certificate,
  but a sharper picture of *how* the gate fails, not just *where*. See
  `code/lobe_structure.py` and Section "Fine structure of the lobes" of the
  paper.

## The critical shell (new, quantified)

Refining the reconnaissance down to energy offsets of 1e-7 shows that even the
MOON component (m_r ≈ 0.0122) of the Earth–Moon pair loses LC-convexity on an ultrathin shell (the Earth component, m_r ≈ 0.9878, has a macroscopic non-convex window, Δc ≲ 1.2e-2)
(c_L1, c_L1 + 9.5e-6] just above the critical energy, where L1 enters the
surface. The trace gate stays positive there; only the determinant fails.
Consequence: an LC-based certification of the full subcritical range stops at
c_L1 + 1e-5, and the remaining shell must be handed to the dynamical criterion
(fiberwise convexity / CZ indices ≥ 0). See `code/critical_shell.py` and
Section "Refinement at the critical limit" of the paper.

## What is provably out of reach for the LC-convexity route

By Kim's theorem (arXiv:1701.07258), the Levi–Civita embedding is non-convex
in the comparable-mass, near-critical corner of parameter space once the mass
ratio drops below `16/17`. No amount of tightening the gates changes this: the
surface itself is not convex there. This is a structural ceiling on the whole
LC-convexity approach, not a limitation of this paper's reduction.

## Two routes past the ceiling (both genuinely open)

### Route 1 — Conley–Zehnder indices (direct dynamical convexity)

Verifying the Hryniewicz–Salomão hypothesis directly, without going through
convexity of any embedding, via the Conley–Zehnder indices of the closed Reeb
orbits (shooting for the periodic retrograde orbit + integrating the
symplectic monodromy + extracting the index). This is the more powerful route
(it is the actual hypothesis of the theorem being invoked) but computationally
heavier, and correct only once the rotating-frame dynamics (Coriolis term) are
properly included in the shooting flow — an earlier attempt used an
inertial/two-center flow missing this term and did not converge. The
simple-fold structure of Section 6.1 narrows the target for this route
considerably: since the non-convexity is confined to thin off-axis lobes
(angular width a few degrees, radial extent within ~4% of the collar) and is
everywhere a single simple fold on a star-shaped base, only Reeb orbits that
actually enter those lobes can encounter the negative principal direction —
the search for problematic orbits does not need to cover the whole surface.
Not included in this repository; a future paper.

### Route 2 — Kim's convex elliptic (two-center) embedding, extended to a=1

Kim proves convexity of the doubly-covered elliptic embedding for the Euler
problem (`a=0`, no rotation; his Thm. 1.1) and, by a regular-perturbation
remark, for small rotation strength (his Rmk. 1.2). Extending this to the full
`a=1` PCR3BP is the open route. Facts established while exploring it:

- **Exact rotation term (proved, audited).** In elliptic coordinates the
  regularized rotation term has the closed form
  `(cosh²λ − cos²ν)(q₁p₂ − q₂p₁) = ½(p_λ sin 2ν + p_ν sinh 2λ)`,
  verified symbolically (residual exactly `0`) and numerically (`<10⁻¹³` over
  2×10⁵ points; `code/rotation_term_audit.py`). This term is **linear in the
  momenta**, so for any `a>0` the regularized Hamiltonian
  `Q_c^a = Q¹(λ,p_λ) + Q²(ν,p_ν) + a·½(p_λ sin 2ν + p_ν sinh 2λ)` is
  **non-mechanical**. Kim's proof reduces fiberwise convexity to the curvature
  of a planar curve, which is valid only in the mechanical case; that reduction
  therefore fails for `a>0`. This is sharper than "does not separate": it
  identifies exactly which structural property is lost and why. Now stated in
  the paper (Section 7, eq. rotation).
- **General 4D convexity tester (built, validated).** A direct projected-
  Hessian convexity test on `{Q_c^a=0}` was implemented and cross-validated
  against Salomão's mechanical criterion at `a=0` (100% sign agreement across
  `μ∈{0.5,0.4,0.3,0.2}`, `c` down to `c_J−10⁻¹⁰|c_J|`).
- **Branch-point analysis (proved, audited).** At the coordinate branch point
  `(λ,ν)=(0,π)`, `∂Q/∂x = a·p_ν`, `∂Q/∂y = −a·p_λ` (audited to `10⁻⁹`). At
  `a=0` the position-gradient vanishes for all momenta (the special critical
  structure Kim exploits); at `a≠0` it does not, so rotation turns the branch
  point into a *regular* point of `{Q=0}`. A preliminary numerical signal that
  small `μ` might lose convexity under rotation earlier than at `a=0` was found
  but **not confirmed**: it is symptomatic of an un-characterized physical
  domain for the compact component when `a≠0` (Kim's domain restriction was
  derived for the mechanical `a=0` case), not of a genuine interior effect.

Next concrete step: characterize the genuine "Earth" compact-component domain
in `(λ,ν,p_λ,p_ν)` for `a≠0` by continuation in `a` from the well-understood
`a=0` region, before any further convexity scan. Route-2 working files are in
the session deliverables (`ROTA_2_kim_euler/`), not in this repository; a future
paper.

## Bottom line

The comparable-mass, near-critical regime of the Birkhoff conjecture for the
PCR3BP remains open. This repository establishes, rigorously, a substantial
and unconditional simplification of the certification machinery everywhere
LC-convexity *does* hold (which — per the reconnaissance — is most of
parameter space, including the entire Earth–Moon regime), and delimits
precisely where the LC-convexity approach cannot succeed, with two concrete,
partially-scoped routes forward.
