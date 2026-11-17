# remlax-validation

Evidence that the R package [remlax](https://github.com/NicoSalas34/remlax)
computes the right restricted likelihood, and measurements of how fast it does
so. remlax is compared with asreml, lme4, sommer, nlme, pbkrtest and kmmr on
the same models and the same data. This repository is kept apart from the
package because most comparisons need software that CRAN cannot depend on
(asreml is commercial) and a computing cluster.

## Latest results

Each benchmark keeps only its latest run, in `<benchmark>/latest/<date>_remlax-<version>/`.
Every run folder has a `RUN.md` that gives the remlax version, the hardware,
the command and the verdict. Older runs are in the git tag `archive/2026-11-17`
(see [docs/history.md](docs/history.md)).

| Question | Benchmark | remlax tested | Date | Verdict |
|---|---|---|---|---|
| Same likelihood as the reference programs? | [correctness/reference_software](correctness/reference_software/) | 3bfae58 (0.2.0.9000) | 2026-11-04 | 381 checks, 365 passed, 16 failed with a diagnosis |
| Recovers the simulated truth? | [correctness/simulations](correctness/simulations/) | 3bfae58 (0.2.0.9000) | 2026-11-04 | 256 fits: 240 passed, 15 failed, 1 expected refusal |
| Agrees on published data sets? | [correctness/published_data](correctness/published_data/) | 0.1.0 to 0.2.0.9000, per data set | 2026-11-01 | every counted check passed, see [docs/published-data.md](docs/published-data.md) |
| Unchanged by a refactoring? | [correctness/golden](correctness/golden/) | references captured at cfec55d (R) and 73e81cf (Python) | 2026-11-01 | reference fits for before/after comparisons |
| Time against n, p and dim(C)/n | [speed/absolute](speed/absolute/) | snapshot v14 | 2026-11-04 | 562 of 589 measurements; figures and tables |
| Time on seven real-data models | [speed/real_models](speed/real_models/) | 6dbada5 (v15) | 2026-10-31 | same optima as on 2026-10-27; tables |
| Time of each remlax engine | [speed/engines](speed/engines/) | none stored | | protocol and scripts only |

The summary for a reader is [docs/validation-summary.md](docs/validation-summary.md).
Every check of the latest campaign is in [docs/validation.md](docs/validation.md).

## Layout

| Folder | Content |
|---|---|
| `correctness/` | Does remlax compute the right thing? Comparisons with reference programs, simulations with a known truth, published data sets, golden references |
| `speed/` | How fast? Absolute benchmark on simulated cases, real-data models, remlax engines, GPU memory |
| `studies/` | Method studies built on remlax: co-GWAS power, demonstrations, genome scans, Neighbor GWAS, pathosystem |
| `cluster/` | SLURM jobs and task lists, one folder per benchmark |
| `docs/` | Summaries and reference notes: validation summary and plan, published data report, CRAN audit and release plan, REML software landscape |
| `paper/` | Draft of the software paper |
| `tools/` | Helpers shared by the scripts: `load_remlax.R`, `remlax_engine.py`, `complexity.py` |

Inside each benchmark, `scripts/` holds the code, `latest/` the latest run and
`README.md` the question, the method and how to rerun it.

## Running

Install remlax first (`remotes::install_github("NicoSalas34/remlax")`) and a
Python with jax, numpy and scipy (`remlax::rx_install_python()` or `RX_PY`).
Every script prints the remlax version it tested. From the root of this
repository:

```sh
Rscript correctness/reference_software/suites/test_remlax_lme4.R       # lme4
Rscript correctness/reference_software/suites/test_remlax_nlme.R       # nlme
Rscript correctness/reference_software/suites/test_remlax_sommer.R     # sommer
Rscript correctness/reference_software/suites/test_remlax_asreml3.R    # asreml (licence)
Rscript correctness/simulations/scripts/verif_simulations.R --cores=10
python3 correctness/published_data/make_report.py
python3 speed/absolute/scripts/summarise.py speed/absolute/latest/2026-11-04_remlax-v14/inputs <output folder>
Rscript speed/real_models/scripts/summarise.R --res speed/real_models/latest/2026-10-31_remlax-6dbada5/inputs \
        --asreml speed/real_models/latest/2026-10-31_remlax-6dbada5/inputs/banc_asreml.csv --versions v15
Rscript correctness/golden/golden.R capture <dir> [<remlax sources>]
```

The full campaigns run on the cluster from `cluster/` (see
[cluster/README.md](cluster/README.md)).

## Keeping the results clean

A new run goes to a new folder `latest/<date>_remlax-<version>/` with its
`RUN.md`. Once it is checked, the previous run folder is deleted in the same
commit. The git history keeps it, and a tag such as `archive/2026-11-17`
marks a state worth finding again.

The scripts write English file names and columns. Some stored inputs predate
this choice and keep French names and columns (`*_chrono.csv`,
`*_delais.csv`, `banc_*.csv`). The scripts that read them map the old names
to the new ones.

## History

This repository was split from the remlax repository on 2026-09-29 with its
history (`git filter-repo`). It was reorganised on 2026-11-17 around the
question each benchmark answers. [docs/history.md](docs/history.md) maps the
old folders to the new ones.
