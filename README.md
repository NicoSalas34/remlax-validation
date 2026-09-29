# remlax-validation

Correctness and speed evidence for the R package
[remlax](https://github.com/NicoSalas34/remlax). This repository compares
remlax with asreml, lme4, sommer, nlme and pbkrtest on the same models and the
same data, and times it against asreml on CPU and GPU. It is kept apart from
the package because most of these comparisons need software that CRAN cannot
depend on (asreml is commercial) and a computing cluster.

The package ships the results: `vignette("validation", package = "remlax")`
reads a copy of `validation/results/checks_<date>.csv` and of the benchmark
table, and `remlax::rx_validate()` reruns the comparisons against lme4, nlme
and sommer on the user's machine.

## Contents

| path | what |
|---|---|
| `tests/R/` | the comparison suites, one line per check (`label`, measured values, verdict); `test_remlax_asreml*.R` need an asreml licence |
| `validation/` | log parsing (`parse_logs.py`), the table of checks (`results/checks_<date>.csv`), the figure, the regeneration of `docs/validation.md` |
| `benchmarks/` | speed at equal model: `bench_complexite.R` (size, model complexity, density), `summarise_benchc.py` (tables and figures), `verif_ecarts_asreml.R` (asreml continued where it stops short) |
| `docs/` | `validation-summary.md` (one page), `validation.md` (every check), `benchmarks.md` (speed), `vram.md` (GPU memory) |
| `paper/` | the draft of the software paper and its validation section |
| `R/charger_remlax.R` | loads the installed remlax package and exposes its namespace to the scripts |
| `moteur_remlax.py` | locates the Python engine shipped with the installed package, for the Python scripts |

## Running

Install remlax first (`remotes::install_github("NicoSalas34/remlax")`) and a
Python with jax, numpy and scipy (`remlax::rx_install_python()` or `RX_PY`).
Every script prints the remlax version it tested. From the root of this
repository:

```sh
Rscript tests/R/test_remlax_lme4.R       # lme4
Rscript tests/R/test_remlax_nlme.R       # nlme
Rscript tests/R/test_remlax_sommer.R     # sommer
Rscript tests/R/test_remlax_asreml3.R    # asreml (licence)
Rscript benchmarks/bench_complexite.R --sizes 500,2000 --software asreml,remlax --tag cpu4
python3 benchmarks/summarise_benchc.py --date <date>
```

The published campaign (2026-09-28/29) ran in the apptainer image
`ige_reml.sif` of the authors' cluster (R 4.6.0, asreml 4.2.0.480, lme4 2.0.1,
sommer 4.4.5, nlme 3.1.169, pbkrtest 0.5.5, jax 0.11.1) against remlax
revision `9a51c79`. The version of each program is recorded in
`validation/results/versions_<date>.json`.

## History

This repository was split from the remlax repository on 2026-09-29 with its
history (`git filter-repo`): the commits that built the validation are the
ones that were made there.