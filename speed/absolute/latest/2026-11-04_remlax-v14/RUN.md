# Absolute benchmark, remlax snapshot v14

| | |
|---|---|
| remlax | snapshot v14, frozen at the end of October 2026; its commit is not recorded in this repository |
| kmmr, asreml | every pass of the benchmark |
| Hardware | genotoul CPU nodes (8 cores) and one A100; a reserved meso node (8 cores) for asreml |
| Summarised | 2026-11-04, `python3 speed/absolute/scripts/summarise.py inputs .` |

## What the inputs hold

`summarise.py` keeps one remlax version per synthesis, v14 by default. The v14
pass redid the `us` cases and filled the missing routes. The other remlax
rows come from the v13 files. kmmr and asreml rows are read from every pass.
Three data sets are kept per point: g1, g2, g3 for one trait, g1, g4, g5 for
several traits.

| Files in `inputs/` | Content |
|---|---|
| `abs_cpu*`, `abs_gpu*` | genotoul passes: first pass, `_fixe`, `_rep`, `_v13`, `_v14` |
| `abs_meso_x13*` | the reserved meso node, with asreml |
| `*_chrono.csv`, `*_delais.csv` | step timings and timeouts (French names of the first passes) |
| `axes.csv` | the axis of each simulated case |

## Verdict

562 of 589 measurements succeeded. The tables and figures below are the
outputs.

| Output | Content |
|---|---|
| `measurements.csv` | one row per case, hardware and route: median time, failures, gap to the best log-likelihood |
| `tables.html` | the times side by side, per case |
| `figures/fig_time_n`, `fig_time_p`, `fig_time_C` | time against n, p and dim(C)/n |
| `figures/fig_rule` | the rule behind `engine = "auto"` |
| `figures/fig_cpu_gpu`, `fig_cpu_gpu_scatter` | CPU against A100 |
| `figures/fig_remlax_asreml` | remlax against asreml on meso |
| `figures/fig_decomposition_cpu`, `fig_decomposition_a100` | where the time of a fit goes |
