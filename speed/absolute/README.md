# Absolute benchmark

Time of a fit against the number of observations n, the number of variance
parameters p and the ratio dim(C)/n, on simulated GBLUP and multi-trait
cases. remlax is timed on its mixed model equations route and its V form
route, against kmmr. asreml runs on a reserved meso node with kmmr and remlax
as controls, because its licence exists on meso only.

| Folder | Content |
|---|---|
| `scripts/grid*.py` | the grids of simulated cases |
| `scripts/summarise.py` | reads the CSV files of a run and writes the tables and figures |
| `latest/` | the latest run: `inputs/` (raw CSV), `measurements.csv`, `tables.html`, `figures/` |

The cluster jobs are in [cluster/speed_absolute/](../../cluster/speed_absolute/).

## Regenerating the outputs

```sh
RUN=speed/absolute/latest/2026-11-04_remlax-v14
python3 speed/absolute/scripts/summarise.py $RUN/inputs $RUN
```

`summarise.py` needs the figure-style helpers in `~/.claude/skills/figure-style`.
