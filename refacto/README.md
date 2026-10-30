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
