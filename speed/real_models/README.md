# Real-data models

Time of seven models fitted on real data: four DROPS maize models
(`drops2_fa1`, `drops2_fa3`, `drops1_fa1`, `drops1_us`) and three Krause
models (`krause_gblup`, `krause_fa2`, `krause_us`). remlax runs on genotoul,
CPU and A100, against kmmr on the same nodes and asreml on meso.

| Folder | Content |
|---|---|
| `scripts/run_case.R` | fits one case with one tool and appends a row to the results |
| `scripts/summarise.R` | Markdown tables from the results of a run |
| `scripts/cas.R`, `sim.R`, `kmmr_cas.R`, `remlax_cas.R`, `asreml_cas.R` | the cases and the call of each tool |
| `latest/` | the latest run: `inputs/` (raw CSV) and `tables.md` |

`scripts/run_case.R` and `scripts/sim.R` also serve the absolute benchmark on
the cluster. The jobs are in [cluster/speed_real_models/](../../cluster/speed_real_models/).

## Regenerating the tables

```sh
RUN=speed/real_models/latest/2026-10-31_remlax-6dbada5
Rscript speed/real_models/scripts/summarise.R --res $RUN/inputs --asreml $RUN/inputs/banc_asreml.csv --versions v15 > $RUN/tables.md
```
