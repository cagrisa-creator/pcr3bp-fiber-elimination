# ERRATUM — version 2 (repository v1.3.0), 2026-07-08

Version 1 of the paper (repository ≤ v1.2.0) contained two compounding errors
in the **floating-point reconnaissance** sections (Sections 5–6 of v1). The
exact fiber-elimination results (the four lemmas and their symbolic audit),
the reproduction of the Earth–Moon certificate values, the lobe fine-structure
observations at comparable masses, and the statement "LC-convex for
Δc ≳ 0.03" are **unaffected**.

## Error 1 — mass-convention mislabeling

In the defining function F of Eq. (2.1), the parameter μ is the mass fraction
of the **distant** primary (at q=(1,0)); the **regularized** primary, at the
origin, has mass fraction m_r = 1 − μ. This follows from the exact identity
F = |q|·(Ω(q) − c) with Ω the effective potential of the pair (mass 1−μ at the
origin, μ at (1,0), barycentre at (μ,0)); the identity is verified symbolically
in `audit/audit_artigo1.py` (Block A) and geometrically by the Hill-radius test
(Block B). Version 1 read μ as the regularized-primary fraction itself, and
consequently mislabeled every component in the reconnaissance narrative
(e.g. the run at μ = 0.987849, labeled "Earth component", regularizes the
**Moon**, m_r = 0.012151).

## Error 2 — scanning-heuristic false negatives (false "convex" entries)

The v1 scanner evaluated the exact fiber minimum of D only at the worst base
points of the loose amplitude bound 𝒟 = A0 − √(A1²+B1²) − √(A2²+B2²). When 𝒟
is looser in the interior than on the collar, all "worst points" are interior
(where the true D is positive) and a genuine collar lobe (where the bound is
sharp and D < 0) goes undetected. Exhaustive base×fiber scans show, e.g., at
μ = 0.05, Δc = 1e−4: true D_min = −2.507 at the collar (θ ≈ 5.7°), where the
v1 CSV reported +5.6. The corrected scanner evaluates D exactly on the union
of the worst points **and the entire collar band**, and was validated against
exhaustive scans at sample grid points (sign agreement everywhere; values
within grid resolution). Both `convexity_map.py` and `critical_shell.py`
shared the flaw; both are fixed.

## Corrected findings (main quantitative changes)

- Near the critical energy (Δc = 1e−4) the non-convex window in m_r spans
  ≈ (0.045, 0.999]: only components regularized at a light primary
  (m_r ≲ 0.045) are LC-convex there.
- Earth–Moon: the **Moon** component (m_r = 0.0122) is LC-convex for all
  scanned Δc ≥ 1e−4 and carries the ultrathin non-convex shell
  (Δc* ≈ 9.5×10⁻⁶) of Section 5.5 (v1 attributed the shell to the Earth
  component). The **Earth** component (m_r = 0.9878) is non-convex on the
  macroscopic window Δc ≲ 1.2×10⁻², recovering between Δc = 1.1×10⁻² and
  1.2×10⁻².
- The v1 claim that the empirical boundary "brackets Kim's 16/17 threshold"
  is **withdrawn**. It rested on the mislabeled variable — the arithmetic
  coincidence 1 − 16/17 = 1/17 placed the (mislabeled) light-side boundary
  near Kim's number — and on false-convex entries. The corrected map is
  consistent with Kim's sufficient condition (every heavy-side point with
  m_r < 16/17 is non-convex near critical) and shows the condition is far
  from sharp: non-convexity persists up to m_r = 0.999. The light-side
  boundary sits at m_r ≈ 0.045 at Δc = 1e−4 (not 1/17 ≈ 0.059).

## What changed in the repository (v1.2.0 → v1.3.0)

- `paper/main.tex` + `paper/main.pdf`: abstract, contribution (iv),
  new Remark 2.1 (mass convention), erratum note opening Section 5,
  Sections 5.2–5.5 rewritten, captions, outlook, data availability.
- `code/convexity_map.py`: convention docstring; collar band added to the
  exact-evaluation set; `m_r` column added to the CSV; map resolution raised
  to 192×80.
- `code/critical_shell.py`: same heuristic fix; component relabeled (Moon).
- `code/lobe_structure.py`: convention note (values unchanged).
- `code/run_all.py`: new regression checks, including one that fails on the
  v1 behavior (Earth component must be non-convex at Δc = 1e−4).
- `code/phase_map_figure.py` (new): regenerates Fig. 1 from the CSV.
- `results/phase_map_data.csv`: regenerated with the corrected scanner.
- `paper/figures/phase_map.png`: regenerated (m_r axis).
- `audit/audit_artigo1.py` (new): independent, self-contained verification of
  both errors (symbolic derivation, Hill-radius test, exhaustive-vs-heuristic
  comparison with N-convergence, calibration against v1's correct numbers).
  Run `python3 audit/audit_artigo1.py`; it must end with
  "AUDITORIA CONFIRMA T1 e T2".

## Note on the companion certificate paper

The interval certificate reproduced in Section 4 (cited as [GrisaCertificate])
uses the same defining function F. Its certified rectangle at μ ≈ 0.9878
therefore concerns the **Moon** component (m_r ≈ 0.0122). The certified
computation itself is unaffected; if that paper's text also labels the
component as the Earth's, the same relabeling applies there.


## Addendum — v1.3.1 (same day)

An independent external audit of v1.3.0 confirmed the mathematical correction
but found residual contamination in auxiliary files. v1.3.1 fixes: HAL
metadata abstract, README body phrases, stale docstrings in
critical_shell.py/convexity_map.py, regenerated results/expected_output.txt
and results/critical_shell_output.txt, regenerated main.docx from the
corrected source, a ~100x speedup of fiber_elimination_audit.py (exact cheap
zero-tests before simplify; identical verdicts), the Earth--Moon Earth
component (mu=0.012151) added to the phase-map grid (16x7), and a new
explicit validation artifact audit/v130_validation.py +
results/v130_validation_output.txt covering every dedicated boundary and
spot value quoted in the text.

## Addendum -- v1.3.2

Second external-audit round: docs/ROADMAP.md decontaminated (last file with
the withdrawn v1 narrative); code-availability version string synchronized;
symbolic-expression cache added so convexity_map.py imports in seconds on
every tested Python (the cache is verified against the full symbolic
derivation by audit/verify_gates_cache.py); FILL_BEFORE_DEPOSIT.md added
(GitHub/Zenodo/HAL identifiers are the depositor's inputs and remain
placeholders until deposit).
