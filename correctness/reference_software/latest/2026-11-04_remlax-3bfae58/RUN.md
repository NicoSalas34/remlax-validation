# Validation campaign of 2026-11-04

| | |
|---|---|
| remlax | commit 3bfae58, version 0.2.0.9000 |
| remlax-validation | commit cc0e235 |
| Image | `ige_reml_v3.sif`: R 4.6.0, asreml 4.2.0.480, lme4 2.0.1, sommer 4.4.5, nlme 3.1.169, pbkrtest 0.5.5, JAX 0.11.1 |
| Hardware | every suite on CPU nodes of the cluster; CPU/GPU parity on a 1g.10gb slice of an A100 |
| Command | `cluster/validation_campaign/prepare.sh`, then `run_extra.sbatch` |

## Verdict

381 checks, 365 passed, 16 failed. These are the same verdicts as on
2026-10-30 with remlax 6dbada5. Each failure has a diagnosis in
[docs/validation-summary.md](../../../../docs/validation-summary.md): 12
concern ma2, sph and cir against asreml, 2 the 2D P-spline, 1 the
multimodal lvr kernel and 1 the random sweep.

The remlax suites ran in the same campaign. pytest gave 300 passed and 2
skipped. testthat gave 1029 expectations with no failure. The suites of the
removed RTMB engine (`10_tmb_parity`, `11_tmb_vs_dense`, `12_tmb_prec`) are
absent on purpose.

## Files

| File | Content |
|---|---|
| `logs/` | the logs of every suite, as written on the cluster |
| `checks.csv` | one row per printed verdict, extracted by `parse_logs.py` |
| `versions.json` | the version of each program, read from the logs |
| `pairs.csv` | every (remlax, reference) pair read from the checks |
| `fig_validation.png` | remlax against the references, and the log-likelihood gaps |
