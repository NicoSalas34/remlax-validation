# Simulation checks of remlax

For each type of model that remlax can fit, this suite checks two things. The returned REML log-likelihood must be exact. The estimated variance components must recover the simulated values. It is not a speed benchmark.

## Running

From the root of the repository:

```sh
Rscript verif_simulations/verif_simulations.R --cores=10
```

The script restarts R with `OPENBLAS_NUM_THREADS=1`. Each JAX process is limited to two threads (`XLA_FLAGS` and `options(remlax.threads = 2)`). On 20 cores, the full suite takes about ten minutes.

Options:

- `--pkg=` path of the remlax sources, loaded by `pkgload::load_all`. Without this option, the script loads the installed package through `R/charger_remlax.R`.
- `--cases=` regular expression on case names (for example `--cases='^mt_'`). `--list` prints the catalogue.
- `--cores=` number of parallel processes (10 by default).
- `--reps=` number of replicates per case, on successive seeds (1 by default).
- `--engines=` regular expression on the routes (`lbfgsb`, `ai`, `mme`, `mme_sparse`, `mme_lbfgsb`).
- `--seed=` base seed. `--out=` folder of the results (`results/` by default).

## What each case does

The case builds a design and simulates y from the true V of the model. It then fits the model by every applicable route. The routes are `rx_fit(optimizer = "lbfgsb")` and `rx_fit(optimizer = "ai")`. The mixed model equations engine, `rx_fit(engine = "mme")`, is added when the model is within its scope, by AI-REML with automatic (`mme`) or sparse (`mme_sparse`) inversion. For a kinship given by its precision and for the ar1 × ar1 field, one more route fits the same engine by L-BFGS-B with sparse inversion (`mme_lbfgsb`). It replaces `rx_fit_sparse()`, the RTMB engine removed from remlax in October 2026. The backend is always the CPU.

The arbiter rebuilds V in R from the Sigma matrices returned by the fit (`fit$sigmas`, `fit$sigmas_res`) and from the level parameters of `fit$varcomp`. It uses the true geometry of the design: coordinates, field cells, incidences built by the case. It then evaluates REML in the V form, with the 2π constant, which is the convention of `fit$logLik`:

    -2 logL = (n - p) log 2π + log|V| + log|X'V⁻¹X| + y'Py

For the spatial lag (`rx_lag`), the arbiter applies the `reml_y` criterion: REML of y* = (I - M) y, Jacobian -2 log|det A| and mean term. The simulation goes through the same V function evaluated at the true parameters. The arbiter code is in `R/outils.R`, the cases in `R/cas.R` and `R/cas_limites.R`.

## Criteria

- arbiter: |logLik - arbiter| / |arbiter| < 1e-8.
- agreement of the routes: the logLik of each route is within 1e-4 of the best one for the case.
- convergence: `conv_decrement` or `conv_grad_rel` for every route. `scipy_success` is reported in the detail and decides nothing.
- recovery: |estimate - true| < 3 standard errors. The standard error comes from `fit$varcomp` (delta method on the Hessian). A route without a standard error takes the one of L-BFGS-B on the same model.

A fit is PASS if it meets the four criteria. EXPECTED_REFUSAL marks a route outside its scope that refuses the model with the expected message (AI-REML on a model with a lag).

## Reading the results

The output folder holds four files:

- `summary.txt`: per case, one line per route (status, logLik, arbiter, differences, max|z|, reasons for failure), then the components of L-BFGS-B.
- `fits.csv`: one line per case and route, with n, number of parameters, logLik, logLik_arbiter, rel_diff, diff_between_engines, convergence and its detail, duration.
- `components.csv`: one line per component, with true, estimate, se, the source of the standard error and z.
- `cases.csv`: family, description and duration of each case.

Component names follow `fit$varcomp$name`. An entry of Sigma is written `term[i,j]`, a level parameter `term!family` (with suffix 1, 2 when the family has several), and a residual section carries the name of its level in `dsum`.

With a single replicate and about 150 components, one or two values of |z| ≥ 3 can appear by chance. A `--reps=5` run on the case concerned settles it.

Two groups of components give high |z| without any computing fault. The variances of the `rx_spl2d` terms rest on 7 to 9 levels. Their distribution is very skewed, and the asymptotic standard error does not describe it. Over 24 replicates, the mean of each variance stays within two standard errors of the mean of the truth. The range of `lvr` enters a peaked correlation. The log-likelihood has a kink each time the range crosses a distance between two positions, and the finite-difference Hessian gives a meaningless standard error (0.01 for an actual difference of 0.1 to 1.4).

## Stored results

The folders `resultats/` (full suite) and `resultats_replicats/` (four replicates of the cases `mt_g_fa2`, `spline_spl2d`, `ige_mixture_two_species` and `metric_lvr`) were produced before the scripts were translated. They keep their French file names (`resume.txt`, `ajustements.csv`, `composantes.csv`, `cas.csv`), column names and values. The table below maps them to the current names.

| stored (French) | current (English) |
|---|---|
| `resume.txt`, `ajustements.csv`, `composantes.csv`, `cas.csv` | `summary.txt`, `fits.csv`, `components.csv`, `cases.csv` |
| `cas`, `moteur`, `modele`, `famille`, `replicat`, `graine`, `secondes` | `case`, `engine`, `model`, `family`, `replicate`, `seed`, `seconds` |
| `logLik_arbitre`, `ecart_rel`, `ecart_entre_moteurs` | `logLik_arbiter`, `rel_diff`, `diff_between_engines` |
| `detail_convergence`, `message_optimiseur`, `n_z_hors` | `convergence_detail`, `optimiser_message`, `n_z_out` |
| `statut`, `raisons` | `status`, `reasons` |
| `composante`, `vrai`, `estime` | `component`, `true`, `estimate` |
| status `PASSE`, `ECHEC`, `ERREUR`, `REFUS_ATTENDU` | `PASS`, `FAIL`, `ERROR`, `EXPECTED_REFUSAL` |
| se_source `propre`, `+bord` | `own`, `+bound` |
| model `principal`, `creux`, `jumeau` | `main`, `sparse`, `twin` |
| term names `champ`, `plante` | `field`, `plant` |
| case prefixes `parente_`, `niveaux_`, `residuelle_`, `metrique_`, `decalage_`, `limite_` | `kinship_`, `levels_`, `residual_`, `metric_`, `lag_`, `boundary_` |

The other case names were translated word for word (for example `ige_melange_deux_especes` became `ige_mixture_two_species`, `limite_us8_presque_rang2` became `boundary_us8_near_rank2`). The seed of a case depends on its position in the catalogue, not on its name, so the seeds are unchanged.
