# History of the repository layout

On 2026-11-17 the repository was reorganised around the question each
benchmark answers. Each benchmark now keeps only its latest run. The state
before the reorganisation is the git tag `archive/2026-11-17`.

## Finding an old run

```sh
git ls-tree -r --name-only archive/2026-11-17 -- banc_kmmr/resultats_2026-10-27
git show archive/2026-11-17:docs/benchmarks.md
git restore --source=archive/2026-11-17 --worktree -- benchmarks/results    # copy back into the working tree
```

## Runs removed from the working tree

| Old path | Content | remlax |
|---|---|---|
| `benchmarks/`, `docs/benchmarks.md` | speed against asreml on size, complexity and density; the engine alone, CPU against a MIG slice | 9a51c79 and earlier (September 2026) |
| `banc_kmmr/resultats_2026-10-27/` (v1 to v10 rows), `VERSION_FIGEE_v*.txt`, `R/tableaux_html.py`, `docs/rapport_nuit_2026-10-27.md` | the night of benchmarks against kmmr, version by version | snapshots v1 to v10 |
| `banc_kmmr/absolu/resultats/abs_meso.csv`, `abs_meso2.csv`, `banc_cpu_cmp_*.csv` | discarded meso passes on overloaded nodes, comparisons of v10 and v13 | v10, v13 |
| `banc_kmmr/absolu/sortie/`, `sortie_v14/` | earlier outputs of the absolute benchmark, French names | v10 to v14 |
| `banc_kmmr/absolu/annuler_*.txt` | lists of jobs to cancel | |
| `validation/results/` (checks, pairs, versions, logs of 2026-09-01, 2026-09-28 and 2026-10-30) | earlier validation campaigns | up to 6dbada5 |
| `verif_simulations/resultats/`, `resultats_replicats/` | earlier simulation checks, French names | before 3bfae58 |
| `refacto/reference_anglais/`, `reference_py_anglais.pkl` | golden references before the centring of y | bc2af0c |

## Old and new folders

| Old | New |
|---|---|
| `tests/R/` | `correctness/reference_software/suites/` |
| `validation/parse_logs.py`, `make_validation_*.py`, `summary_stats.py` | `correctness/reference_software/report/` |
| `validation/results/` | `correctness/reference_software/latest/<run>/` |
| `verif_simulations/` | `correctness/simulations/` (`scripts/`, `latest/`) |
| `published_data/` | `correctness/published_data/` (`results/` became `latest/<run>/`) |
| `refacto/` | `correctness/golden/` (`reference_centrage/` became `references/r/`) |
| `banc_kmmr/absolu/` | `speed/absolute/` (scripts), `cluster/speed_absolute/` (jobs, task lists) |
| `banc_kmmr/absolu/synthese_abs.py`, `grille*.py` | `speed/absolute/scripts/summarise.py`, `grid*.py` |
| `banc_kmmr/R/` | `speed/real_models/scripts/` (`executer.R` became `run_case.R`, `synthese.R` became `summarise.R`) |
| `banc_kmmr/*.sbatch`, `taches.txt` | `cluster/speed_real_models/` |
| `banc_moteurs/` | `speed/engines/` (scripts), `cluster/speed_engines/` (jobs) |
| `docs/plan_banc_moteurs.md` | `speed/engines/PROTOCOL.md` |
| `docs/vram.md` | `speed/gpu_memory.md` |
| `campaign/` | `cluster/validation_campaign/` |
| `simulation/cogwas/`, `docs/rapport_qtl_additifs_2026-11-02.md`, `docs/figures_qtl_additifs/` | `studies/cogwas/` |
| `demo/` (`test_*.R`, `sortie_*/`) | `studies/demos/` (English script names, `output/`) |
| `validation/etude_scan_*` | `studies/genome_scan/study_scan_*` |
| `validation/neighbor_multitrait_sim.R` | `studies/neighbor_gwas/` |
| `R/charger_remlax.R`, `moteur_remlax.py` | `tools/load_remlax.R`, `tools/remlax_engine.py` |
| `docs/paysage_logiciels_reml.md`, `docs/plan_branche_cran.md` | `docs/reml_software_landscape.md`, `docs/cran_release_plan.md` |
