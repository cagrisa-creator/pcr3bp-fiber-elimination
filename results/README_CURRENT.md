# Current results and historical output records

The current numerical evidence consists of finite floating-point samples.
The exact symbolic checks verify identities. Neither positive samples nor
symbolic identities constitute an interval convexity certificate; no such
certificate is supplied for a parameter pair or parameter region.

Current evidence is recorded in:

- `phase_map_data.csv`: the revised sampled global map.
- `parameter_audit/`: spot checks, mass/energy and lunar refinements,
  before/after comparisons, and actual negative-gate numerical witnesses.
- `figure_audit/`: plotted samples, boundaries, ray profiles, reflection checks,
  and the data and comparisons underlying the revised figures.
- `symbolic_audit/`: exact-identity execution logs and related cross-checks.
- `package_driver.log`: the revised diagnostic-driver log, which uses base-axis
  terminology and "one negative direction" for sampled lobe points.

The root `REPRODUCIBILITY.txt` provides commands, grids, dependencies, the
result-to-file map, and the limitations of the computations. The additional
`v130_validation_output.txt` retains its explicitly versioned finite-grid
cross-check status; it is not a new interval certificate.

The files formerly named `critical_shell_output.txt`,
`lobe_structure_output.txt`, and `expected_output.txt` in this directory were
legacy outputs. Their exact bytes, together with their old `RUN_ENV.txt` and
old full-repository `SHA256SUMS`, are preserved in
`../historical/pre_R07/results/`. They are not the current expected output and
must not be used to restore withdrawn orbit, limiting-value, or singularity
claims. They were archived rather than edited to appear newly generated.

`ORIGINAL_FULL_TREE_MANIFEST_REFERENCE.txt` is an explicitly historical
manifest reference for the older full repository. It is identical to the
archived legacy manifest and is not the current archive manifest.

`../frozen_inputs/` contains immutable historical comparison inputs. The audit
commands depend on the preserved code and reference CSV there; the frozen logs
are historical, as explained in `../frozen_inputs/README_HISTORICAL.md`.
