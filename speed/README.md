# Speed

How long does a fit take with remlax, against kmmr and asreml, on CPU and on
an A100? Only the latest tested versions of remlax are kept here.

| Folder | Question | Latest run |
|---|---|---|
| [absolute/](absolute/) | time against the number of observations n, of parameters p and dim(C)/n, on simulated cases | 2026-11-04, remlax snapshot v14 |
| [real_models/](real_models/) | time on seven DROPS and Krause models | 2026-10-31, remlax 6dbada5 (v15) |
| [engines/](engines/) | time of each remlax engine across model families | no run stored |
| [gpu_memory.md](gpu_memory.md) | GPU memory of a large direct and indirect effects model | note |

The automatic choice of engine in remlax (`engine = "auto"`) follows the rule
drawn in `absolute/latest/2026-11-04_remlax-v14/figures/fig_rule.png`.
