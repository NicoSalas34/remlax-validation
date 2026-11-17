# Published data sets

Each script of `scripts/` downloads a published data set, fits it with remlax
and with the programs of the source paper, and writes one table of checks.
`make_report.py` pools the tables into
[docs/published-data.md](../../docs/published-data.md) and a figure.

| Folder | Content |
|---|---|
| `scripts/` | one script per data set, `01_wheat_crossa2010.R` to `15_drops_multitrait.R` |
| `R/common.R` | shared tolerances, the check table and the folders |
| `latest/` | the latest run: one checks table and one log per data set, and the figure |
| `data_cache/` | downloaded data, not versioned |

## Running

```sh
for s in correctness/published_data/scripts/[0-9][0-9]_*.R; do Rscript "$s"; done
python3 correctness/published_data/make_report.py
```

The scripts write into the newest folder of `latest/`. For a new full run,
create `latest/<date>_remlax-<version>/` first, then delete the previous one.
