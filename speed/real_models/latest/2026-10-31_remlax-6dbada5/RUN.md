# Real-data models, remlax 6dbada5

| | |
|---|---|
| remlax | commit 6dbada5, the released version, frozen as v15; run on 2026-10-30 |
| kmmr | the run of 2026-10-27 on the same genotoul nodes |
| asreml | meso, 8 cores |
| Hardware | genotoul CPU nodes (8 threads) and one A100 |
| Command | `cluster/speed_real_models/run.sbatch`, then `scripts/summarise.R --versions v15` |

## Verdict

Every remlax optimum equals the one of the run of 2026-10-27. Times vary on
the shared genotoul nodes. `tables.md` gives the times of each tool per case.

## What the inputs hold

| Files in `inputs/` | Content |
|---|---|
| `bench_cpu_v15*`, `bench_gpu_v15*` | the remlax rows of this run, with their step timings |
| `banc_cpu*`, `banc_gpu*` (`_v2`, `_v3`) | kmmr rows of 2026-10-27 (French names) |
| `banc_asreml.csv` | asreml on meso |
