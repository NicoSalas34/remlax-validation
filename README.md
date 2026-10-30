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
| `campaign/` | the validation campaign on the cluster: `prepare.sh` copies remlax and this repository at their commits and submits `run_cpu.sbatch` (every suite, CPU) and `run_gpu.sbatch` (CPU/GPU parity) |
| `tests/R/` | the comparison suites, one line per check (`label`, measured values, verdict); `test_remlax_asreml*.R` need an asreml licence |
| `validation/` | log parsing (`parse_logs.py`), the table of checks (`results/checks_<date>.csv`), the figure, the regeneration of `docs/validation.md` |
| `benchmarks/` | speed at equal model: `bench_complexite.R` (size, model complexity, density), `summarise_benchc.py` (tables and figures), `verif_ecarts_asreml.R` (asreml continued where it stops short); `bench_engines.R` and `grille.R` time the dense engine against the mixed model equations engine with sparse inversion |
| `banc_moteurs/` | the engine benchmark of `docs/plan_banc_moteurs.md`: every remlax engine on CPU and GPU, with asreml, sommer and lme4, on the same simulated data and against a common dense REML arbiter |
| `banc_kmmr/` | remlax against kmmr and asreml-R on CPU and GPU; `R/` and `resultats_2026-10-27/` hold the seven DROPS and Krause models, `absolu/` the absolute benchmark on simulated cases (time against n, p and dim(C)/n), summarised by `absolu/synthese_abs.py` with figures in `absolu/sortie_v14/` |
| `verif_simulations/` | simulation checks: each model type is fitted by every applicable route, its REML log-likelihood is compared with a dense REML arbiter written in R, and its variance components with the simulated truth |
| `demo/` | simulated demonstrations: a field trial with a 2D P-spline surface and direct and indirect genetic effects, a one-step GWAS with LOCO, and studies of neighbourhood kernels and of the identifiability of the dilution |
| `refacto/` | the safety nets of remlax: `golden.R` captures and compares whole fits through R, `golden_py.py` the Python engine values at fixed theta, both pinned to one CPU core, against the references `reference_anglais/` and `reference_py_anglais.pkl` |
| `docs/` | `validation-summary.md` (one page), `validation.md` (every check), `validation-plan.md` (the claims and their evidence), `benchmarks.md` (speed), `vram.md` (GPU memory), `rapport_nuit_2026-10-27.md` (remlax against kmmr and asreml-R), `paysage_logiciels_reml.md` (survey of REML software), `plan_banc_moteurs.md` (protocol of the engine benchmark), `plan_branche_cran.md` (conditions of the CRAN release), `cran_policy_audit.md` (CRAN policy audit), `note_remlax.md` (the original note on the solver, translated from French) |
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
Rscript verif_simulations/verif_simulations.R --cores=10
Rscript demo/demo_remlax_spatial_cpu.R
Rscript banc_moteurs/R/executer.R --family=grm --n=2000 --rep=1 --out=results
python3 banc_kmmr/absolu/synthese_abs.py banc_kmmr/absolu/resultats banc_kmmr/absolu/output_v14
Rscript refacto/golden.R capture <dir> [<remlax sources>]
python3 refacto/golden_py.py capture <out.pkl> <remlax repository>
```

The `.sbatch` files in `banc_moteurs/` and `banc_kmmr/` launch the full
benchmarks on the authors' clusters.

The scripts write English file names and columns. The stored results predate
this choice and keep their original French names and columns (for instance
`banc_kmmr/absolu/resultats/`, `verif_simulations/resultats/ajustements.csv`,
the `sortie_*` folders of `demo/`). The scripts that read stored results map
the old names to the new ones on reading. `banc_kmmr/absolu/sortie_v14/` holds
the stored figures and tables of the absolute benchmark. The command above
regenerates them under English names in a new folder.

The published campaign (2026-09-28/29) ran in the apptainer image
`ige_reml.sif` of the authors' cluster (R 4.6.0, asreml 4.2.0.480, lme4 2.0.1,
sommer 4.4.5, nlme 3.1.169, pbkrtest 0.5.5, jax 0.11.1) against remlax
revision `9a51c79`. The version of each program is recorded in
`validation/results/versions_<date>.json`. The RTMB sparse engine was removed
from remlax after this campaign, together with its suites, and the mixed model
equations engine (`rx_fit(engine = "mme")`) now covers its scope.

## History

This repository was split from the remlax repository on 2026-09-29 with its
history (`git filter-repo`): the commits that built the validation are the
ones that were made there.