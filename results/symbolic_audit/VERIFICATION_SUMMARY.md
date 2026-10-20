# Original audits and independent high-precision witness check

The three original scripts under `work/canonical` were executed unchanged. All exited with status zero:

|Script|Outcome|Runtime|
|---|---|---|
|`code/fiber_elimination_audit.py`|Eight residuals exactly zero; `PASS_C82_FIBER_ELIMINATION_AUDIT`|1.50745 s|
|`code/rotation_term_audit.py`|Symbolic residual exactly zero; 200,000 random checks, maximum absolute residual `4.263e-14`; `PASS_ROTATION_TERM_AUDIT`|4.40182 s|
|`audit/verify_gates_cache.py`|All six cached expressions have exactly zero residual against full symbolic reconstruction; `PASS_GATES_CACHE`|173.43035 s|

This is a total of fifteen exact zero-residual checks. The cache does correspond to the defining formulas; it is not merely consistent at a numerical sample. The rotation script's numerical pass threshold is `1e-9`; the measured error is much smaller. Its random seed and sample ranges are those of the unchanged original source.

`high_precision_witnesses.py` independently differentiates the defining formula `F` using mpmath at 80 decimal digits. It imports neither SymPy nor the expression cache, uses the four referee coordinates and the exported 4096-grid minimizing fiber angles, and constructs `L`, `T`, and `D` directly. At all four points, `F>0`, `calT>0`, `T>0`, `D<0`, and the matrix `C` has one negative and one positive eigenvalue. The maximum absolute differences from the stored double-precision direct-matrix calculation are:

|Quantity|Maximum absolute difference|
|---|---|
|F|`6.07310e-17`|
|calT|`3.52847e-15`|
|T|`3.52102e-15`|
|D|`2.32999e-14`|

For example, at `z=(.463842748,.0569527281)`, `phi=4.331961744989051`, the direct high-precision determinant is approximately `-.3117683021202375501`. Complete decimal values, inputs and per-point differences are in `high_precision_witnesses.json`. Runtime is approximately `.008 s` excluding imports.

High-precision numerical agreement is not an interval error enclosure and is not described as a validated numerical certificate. The exact symbolic checks certify the listed algebraic identities, not positivity over any parameter or base domain.

Use the three original sources, their `.log` files, `original_audits_summary.json`, and the independent high-precision source/result for the final reproducibility package. `verify_gates_cache_explicit.py` and `verify_cache_polynomial.py` were exploratory performance alternatives started while the original cache checker was still running. They are not needed, were abandoned after the original checker passed, and should remain scratch files rather than final verification tools.
