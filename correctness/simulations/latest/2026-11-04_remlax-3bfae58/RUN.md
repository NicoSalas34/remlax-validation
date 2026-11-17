# Simulation checks of 2026-11-04

| | |
|---|---|
| remlax | commit 3bfae58, version 0.2.0.9000 |
| Hardware | CPU node of the cluster, 20 cores, 34.9 minutes |
| Command | `Rscript correctness/simulations/scripts/verif_simulations.R`, from `cluster/validation_campaign/run_extra.sbatch` |
| Criteria | log-likelihood within 1e-8 relative of the dense arbiter; within 1e-4 of the best route; every component within 3 standard errors of the truth |

## Verdict

256 fits: 240 passed, 15 failed, 1 refused as expected. The failures fall
on four cases. On `mt_g_fa2` (4 routes), one component lies 3.06 standard
errors from the truth while the log-likelihood equals the arbiter's. On
`boundary_ar1_rho_097` (6 routes) and `boundary_cor_097` (4 routes), a
correlation near 0.97 gives a component with |z| of 3 or more; four of these
fits also stop 4.5e-4 below the best route, and one differs from the arbiter
by 3.1e-8 relative. On `boundary_us8_near_rank2`, L-BFGS-B stops 2.3e-3
below the best route. `summary.txt` gives the reason of each failure.

## Files

| File | Content |
|---|---|
| `summary.txt` | every fit with its verdict, then each component against the truth |
| `fits.csv` | one row per fit and route |
| `components.csv` | one row per estimated component |
| `cases.csv` | the simulated cases |
