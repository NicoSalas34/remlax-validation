# Cluster jobs

SLURM jobs and task lists, one folder per benchmark. The jobs ran on two
clusters: genotoul (`/work/user/nsalas`) and meso (`/scratch/users/salasn`).
They read the repository from a copy on the cluster, so paths inside them are
cluster paths.

| Folder | Benchmark | Entry point |
|---|---|---|
| `validation_campaign/` | [correctness/reference_software](../correctness/reference_software/) and [correctness/simulations](../correctness/simulations/) | `prepare.sh` copies remlax and this repository at their commits, then submits `run_cpu.sbatch` and `run_gpu.sbatch`; `run_extra.sbatch` runs the remlax testthat suite and the simulation checks |
| `speed_absolute/` | [speed/absolute](../speed/absolute/) | `prepare*.sbatch` simulate the cases, `run*.sbatch` time them, `fetch_results.sh <run>/inputs` brings the CSV back; `tasks/` holds the task lists |
| `speed_real_models/` | [speed/real_models](../speed/real_models/) | `run.sbatch` with `tasks.txt` |
| `speed_engines/` | [speed/engines](../speed/engines/) | `run_*.sbatch`; `rerun_missing.sh` resumes unfinished cases |

On genotoul, `$K/banc` is a copy of `speed/real_models/` and `$K/absolu`
holds the task lists and the simulated cases of the absolute benchmark.
