# Comparison with reference programs

remlax fits the same model on the same data as asreml, lme4, sommer, nlme and
pbkrtest. A check passes when the restricted log-likelihood at the optimum
agrees to 1e-6 absolute (1e-5 against asreml), the variance components to
1e-4 relative (1e-3 against asreml), the fixed effects to 1e-6 relative and
the BLUPs to 1e-4 absolute.

| Folder | Content |
|---|---|
| `suites/` | the R test suites, one line per check with its verdict; `test_remlax_asreml*.R` need an asreml licence |
| `report/` | `parse_logs.py` reads the logs into `checks.csv`; `make_validation_fig.py` draws the figure and writes `pairs.csv`; `make_validation_md.py` writes [docs/validation.md](../../docs/validation.md); `summary_stats.py` gives the largest gaps per reference |
| `latest/` | the latest campaign, with its logs, table of checks, pairs and figure |

## Running a campaign

The campaign runs on the cluster: `cluster/validation_campaign/prepare.sh`
copies remlax and this repository at their commits and submits the jobs. Then,
with `RUN=correctness/reference_software/latest/<date>_remlax-<commit>`:

```sh
rsync -a <host>:<work folder>/logs/ $RUN/logs/
python3 correctness/reference_software/report/parse_logs.py $RUN/logs=cluster --out $RUN/checks.csv --versions-out $RUN/versions.json
python3 correctness/reference_software/report/make_validation_fig.py --csv $RUN/checks.csv --out $RUN/fig_validation.png --pairs-out $RUN/pairs.csv
python3 correctness/reference_software/report/make_validation_md.py --csv $RUN/checks.csv --versions $RUN/versions.json --rev <commit> --date <date>
```

Write the `RUN.md` of the new run, then delete the previous run folder.
