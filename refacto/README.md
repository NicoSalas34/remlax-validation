# Safety nets of the remlax cleanup

- `golden.R`: 36 fits through R (engines, optimisers, structures, outputs,
  scans, Bayesian route), compared to the bit:
  `Rscript golden.R capture <dir> <remlax repo>` then
  `Rscript golden.R compare <reference> <dir>`.
- `golden_py.py`: Python engine values at fixed theta on every route (dense,
  low-memory, blocks, the three MME inversions), MME outputs and whole AI
  fits: `python3 golden_py.py capture <out.pkl> <remlax repo>` then
  `python3 golden_py.py compare <ref.pkl> <out.pkl>`. It pins itself to one
  CPU core: XLA on CPU is not deterministic across runs otherwise.

References:

- `reference_anglais/` (golden.R) and `reference_py_anglais.pkl`
  (golden_py.py), captured on remlax bc2af0c after the move to English names
  and messages. Every number equals those of the earlier references
  (before the cleanup, during it, and after the defect fixes of
  2026-10-30) once the old names are mapped to the new ones. The earlier
  references are in the history of this repository, up to commit b55f225.

- `reference_centrage/` (golden.R, remlax cfec55d) and
  `reference_py_centrage.pkl` (golden_py.py, remlax 73e81cf), captured
  after the likelihood moved to y centred on the span of X (remlax
  63b697a). Use these from now on. Against `reference_anglais`:
  - values at fixed theta move by 1e-9 to 1e-12; Hessians by at most
    3e-8 relative to their largest entry (the elementwise gaps near 2 come
    from entries close to zero);
  - the Bayesian draws differ by rounding (mean gap 0.007 for a posterior
    sd of 0.2), as any MCMC chain would;
  - `scan_exact` gains its `df` column (multi-trait scan, remlax fdd6948);
  - `ige_dense`: one component near the floor is now classed positive
    instead of at the bound, same logLik;
  - two optimisation paths change, both open issues: `ige_mme` stops
    6.8e-7 under the optimum on a face of the cone, without the dense
    polish (engine "mme" instead of "mme+dense"); `decalage` (spatial lag)
    ends at logLik -351.0 instead of -115.2. The reference had not
    converged either (max gradient 6086) and an older code reached
    -87.35: the lag optimiser is unstable on this case.
  golden_py at 73e81cf is identical to the bit to golden_py at 57eb468:
  the refactors of 2026-11-01 do not move any engine value.
