PARAMETER REPRODUCTION AND AUDIT
==============================

Scope
-----
All numerical runs in this directory use finite float64 sampling, plus separate
80-decimal-digit checks of actual negative witnesses. No directed-rounding
interval arithmetic, parameter-box certificate, domain-wide positivity
certificate, or continuous transition theorem is supplied. A positive sample
does not prove positivity between samples. Even an 80-digit negative evaluation
is a high-precision numerical check, not an interval-certified witness.

Use immutable original and revised source directories separately. For a full
comparative reproduction, with numpy/scipy/sympy/mpmath installed:

  python audit_parameters.py --source ORIGINAL_FROZEN/code \
      --reference ORIGINAL_FROZEN/results/phase_map_data.csv \
      --corrected-source REVISED/code --output audit_results --suite all

The script can instead run one suite: phase, critical, spot, corrected,
threshold, or witnesses. The corrected suite requires the historical phase and
lunar JSON records produced by the phase and critical suites. The witnesses
suite should be run last. --deps DIRECTORY optionally adds an existing local
Python dependency directory to sys.path; it is not needed in a normal installed
environment. --export-only re-exports already computed lunar/spot JSON data.

No source file is modified. The historical algorithm's finite fiber sample,
worst-amplitude-bound candidates union outer three radial rows, first radial
root search, and boundary NaN omission are preserved in the historical runs.
The historical API and CSV name D_true_min mean a sampled estimate here, not
an exact continuous minimum over the fiber.

File map
--------
phase_map_rerun.csv:
  Original algorithm, 16 masses x 7 energy offsets, 192 first-quadrant base
  rays x 80 radial rows, 256 fiber angles, 80 worst-bound candidates union
  the outermost three radial rows. No sign differs from the frozen CSV.

phase_map_data_revised.csv:
  Same grid and selection, with the nonsingular F=0 formulas evaluated
  explicitly on every numerically located boundary point. This is the CSV
  corresponding to corrected_code/convexity_map.py and the revised phase map.

phase_revised_comparison.csv:
  All 112 rows, including T, determinant lower bound, sampled determinant
  minimum, historical rerun and frozen-CSV values, revised values, differences,
  and boundary NaN counts. The older phase_comparison.csv compares the frozen
  CSV with the unchanged algorithm on the current environment.

phase_details.json / phase_witnesses.csv:
  Historical sample locations and direct gate evaluations.
phase_revised_details.json / phase_revised_witnesses.csv:
  Revised sample locations and direct gate evaluations.
phase_boundary_complete_details.json:
  Additional diagnostic that also assembles the interior coefficients from
  cheap derivative blocks. This is an independent algebraically equivalent
  diagnostic, not the principal revised map. The corrected source keeps the
  audited original cached expressions in the interior.

lunar_details.json / lunar_track_witnesses.csv / lunar_bisection.json:
  Original-method lunar refinement, 600 rays x 400 radial rows, 1024 fiber
  angles, 120 worst-bound candidates union outer three rows, ten energy
  offsets and 30 floating-point energy bisection steps.
lunar_stability.json / lunar_stability_witnesses.csv:
  Radial sizes 60,150,400 and fiber sizes 512,1024 at delta_c=3e-6.
lunar_revised_details.json / lunar_revised_track_witnesses.csv /
lunar_revised_bisection.json:
  Actual corrected critical_shell.py compared with the historical track.
  This supplies a sampled sign transition, not a certified shell or limit.

spot_details.json / spot_witnesses.csv:
  Four corners of the Section-4 rectangle plus the preserved validation
  checks and an additional terrestrial delta_c=.01 check. The Section-4
  corner grid is newly declared explicitly: 400 rays x 160 radial rows,
  all 512 fiber samples for D. The reduced T is sampled over the base and
  collar D is evaluated from its F=0 formula. No dedicated file specifying
  the original Section-4 corner grid was recovered. The fresh grid supports
  the existing conservative rounded sampled values 2.29e-3 and 6.95e-4.

light_threshold_details.json / light_threshold_bisection_details.json /
light_threshold_witnesses.csv:
  Dedicated exhaustive finite base x fiber scans at delta_c=1e-4, 600 rays
  x 400 radial rows x 1024 fiber angles, with F=0 formulas at the boundary.
  The regularized mass fractions .04,.045,.05 are checked, followed by
  12 empirical mass-bisection steps in [.04,.05]. No continuous mass
  monotonicity is proved. The apparent transition depends on sampling.

negative_witnesses_80digits.json / negative_witnesses_80digits.csv:
  Direct evaluations of equations (2)-(5) using 80-digit arithmetic, not
  the determinant lower bound. A documented 1e-8 relative inward radius
  shift is used where a numerical boundary minimizer would otherwise have
  an ambiguous tiny F residual. Coordinates, angle, mu, c, F, actual D,
  actual T, reduced T, precision and shift are recorded.

Source corrections
------------------
convexity_map.py.patch and critical_shell.py.patch preserve the originals.
The boundary helper applies only to the last radial row:
  T=|grad F|^2,
  A0=(grad F)^T adj(-Hess F/2) grad F,
  A1=B1=A2=B2=0.
No negative interior F is clamped to a valid point. This removes NaNs caused
by tiny negative floating-point root residuals inside sqrt(F) on the boundary.
The actual boundary itself is still located numerically, without an enclosure.

critical_shell.py reuses the audited cached coefficient functions already used
by convexity_map.py, avoiding the expensive duplicate symbolic determinant
construction. Sampling grids and worst-point/collar selection remain unchanged.
Explanatory output labels were aligned with the manuscript: no claim that a
linear extrapolation proves a limiting value, no identification of a base ray
with a periodic orbit, and no presentation of sampled minima as exact minima.

Environment manifests record Python/dependency versions and original source
hashes. Summary JSON files record comparison magnitudes, boundary counts,
observed sign changes and explicit limitations. The final package should hash
the final source files and numerical outputs as a whole after assembly.
