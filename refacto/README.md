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

- `reference_avant_nettoyage/` (golden.R) and `reference_py_bbc7340.pkl`
  (golden_py.py): captured before and during the cleanup of 2026-10. Every
  cleanup step compared identical to them.
- `reference_stabilisation/` and `reference_py_stabilisation.pkl`: the
  current references, after the defect fixes of 2026-10-30. Against the
  former ones they differ only where the fixes meant to: `n_iter` of the
  average information optimiser (accepted iterations instead of steps
  tried), and the constraint code B, without standard error, of variances at
  the raised floor of the MME engine (cases fa1_mme and ige_mme). Engine
  values, gradients and AI matrices are unchanged.
