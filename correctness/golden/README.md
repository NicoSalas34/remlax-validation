# Golden references

Whole fits captured before a refactoring of remlax and compared after it. A
refactoring is accepted only when these numbers do not move.

- `golden.R`: 36 fits through R (engines, optimisers, structures, outputs,
  scans, Bayesian route), compared to the bit:
  `Rscript golden.R capture <dir> <remlax repo>` then
  `Rscript golden.R compare references/r <dir>`.
- `golden_py.py`: Python engine values at fixed theta on every route (dense,
  low-memory, blocks, the three MME inversions), MME outputs and whole AI
  fits: `python3 golden_py.py capture <out.pkl> <remlax repo>` then
  `python3 golden_py.py compare references/python.pkl <out.pkl>`. It pins itself to one
  CPU core: XLA on CPU is not deterministic across runs otherwise.

References:

- `references/r/` (golden.R, remlax cfec55d) and
  `references/python.pkl` (golden_py.py, remlax 73e81cf), captured
  after the likelihood moved to y centred on the span of X (remlax
  63b697a). Use these from now on. Against the previous references (`reference_anglais`,
remlax bc2af0c, kept in the git tag `archive/2026-11-17`):
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
