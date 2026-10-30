remlax · night of 26 to 27 October 2026

# remlax against kmmr and asreml-R

Solver optimisation, full validation and a benchmark on seven quantitative genetics models, CPU and GPU, frozen version v10.

## Summary

On each of the seven cases, remlax reaches an optimum at least as high as kmmr and asreml-R. On four cases, it is strictly higher. The two one-stage DROPS models make kmmr break down, and it stops after two iterations. asreml caps their ar1 correlations at 0,98, pins specific variances at zero on drops2_fa3 and stops 3,8 units below remlax on krause_fa2.

On A100, remlax is the fastest on all Krause models: 22 times kmmr on krause_fa2, nearly 4 times on krause_us. On CPU, it beats kmmr on krause_us (316 s against 359 s) and drops2_fa3. It matches kmmr on the GBLUP and drops2_fa1, and it remains slower on krause_fa2 (995 s against 813 s), where it climbs 20 units higher. On the small models, most of its time is fixed: about 2,5 s to start the Python solver, once per R session.

×22 (krause_fa2 on A100: 44,4 s against 961 s for kmmr, higher optimum).

×3,7 (krause_us on A100: 22,1 s against 82,6 s).

357 → 29 s (drops1_fa1 on CPU, from v1 to v10, same nodes).

36 / 37 (cases of the simulation check; the last one is a known random failure).

## Benchmark results

kmmr and remlax run on the same genotoul nodes, AMD EPYC 7513 with 8 cores or A100 80 GB, in the same image. asreml-R runs on meso, where its licence is available. Its times are therefore only orders of magnitude. The kmmr time includes the inversion of K done in R before the call. The remlax time includes the whole `rx_fit()`, solver start included, at the second fit of the session (compiled programs read back from disk). ▼ marks an optimum lower than the best of the three.

| case | REML log-likelihood at the optimum |  |  | CPU, 8 cores (s) |  |  | A100 (s) |  |
|----|----|----|----|----|----|----|----|----|
|  | kmmr | asreml-R | remlax | kmmr | remlax | asreml | kmmr | remlax |
| DROPS two-stage, fa1 | -492.86 | -492.86 | -492.86 | 6.2 | 7.2 | 0.5 | 8.3 | 8.2 |
| DROPS two-stage, fa3 | ▼ -368.06 | ▼ -378.25 | -367.68 | 16.0 | 9.7 | 2.4 | 19.6 | 10.0 |
| DROPS one-stage, fa1 | ▼ -2154.65 | ▼ -1161.86 | -1161.21 | 3.9 | 29.1 | 15.3 | 5.0 | 42.7 |
| DROPS one-stage, us(5) | ▼ -2169.30 | ▼ -1146.72 | -1146.19 | 5.8 | 18.0 | 10.5 | 5.6 | 28.3 |
| Krause, GBLUP | 389.57 | 389.57 | 389.57 | 9.1 | 9.0 | 60.2 | 7.4 | 9.2 |
| Krause, fa2 × 5 treatments | ▼ 7109.92 | ▼ 7126.36 | 7130.18 | 813 | 995 | 7429 | 961 | 44.4 |
| Krause, us(5) × 5 treatments | 7138.93 | 7138.93 | 7138.93 | 359 | 316 | 7530 | 82.6 | 22.1 |

kmmr at a tolerance of 1e-6 on the log-likelihood. asreml-R with its default settings, with a 30 GB workspace for Krause. On drops1, asreml recovers the remlax optimum to within 0,02 when the correlations are fixed at 0,98 as it does.

### Iterations and cost of one evaluation

| case | kmmr it. | asreml it. | remlax AI steps tried | remlax s/eval CPU | remlax s/eval A100 | version |
|----|----|----|----|----|----|----|
| DROPS two-stage, fa1 | 23 | 13 | 24 | 0.071 | 0.030 | v10 |
| DROPS two-stage, fa3 | 38 | 30 | 138 | 0.070 | 0.033 | v10 |
| DROPS one-stage, fa1 | 2 | 24 | 76 | 0.435 | 0.607 | v10 |
| DROPS one-stage, us(5) | 2 | 23 | 30 | 0.416 | 0.608 | v10 |
| Krause, GBLUP | 11 | 7 | 4 | 0.676 | 0.056 | v10 |
| Krause, fa2 × 5 treatments | 14 | 15 | 56 | 28.122 | 1.252 | v10 |
| Krause, us(5) × 5 treatments | 12 | 9 | 9 | 30.670 | 1.242 | v10 |

The AI steps tried by remlax also count the rejected steps, evaluated on the likelihood alone. On the fa models, AI convergence stays linear near the end. This is the main limitation of remlax on CPU.

### Progress of remlax across the frozen versions

| case, hardware | v1 | v2 | v3 | v6 | v10 |
|----|----|----|----|----|----|
| DROPS two-stage, fa3, CPU | 186 s -367.68 | 12.9 s -367.68 | 10.4 s -367.68 | 10.4 s -367.68 | 9.7 s -367.68 |
| DROPS one-stage, fa1, CPU | 357 s -1161.21 | 25.4 s -1161.21 | 415 s -1161.21 | 34.5 s -1161.21 | 29.1 s -1161.21 |
| Krause, GBLUP, CPU | 11.7 s 389.57 | 14.4 s 389.57 | 11.1 s 389.57 | 10.0 s 389.57 | 9.0 s 389.57 |
| Krause, us(5) × 5 treatments, CPU |  | 1286 s 7138.93 |  | 444 s 7138.93 | 316 s 7138.93 |
| Krause, fa2 × 5 treatments, A100 | 480 s 7130.18 |  |  | 66.6 s 7130.18 | 44.4 s 7130.18 |
| Krause, us(5) × 5 treatments, A100 | 266 s 7138.93 |  |  | 28.0 s 7138.93 | 22.1 s 7138.93 |

v1 to v3: first optimisations (starting values, trials on the value alone, BLAS threads, full blocks). v6: AI step limits and a GPU route without transfers. v9: cache of the last evaluation and Heywood rule. v10: an accepted AI step reuses the factorisation of its trial. v3 found drops1_fa1 in 415 s on CPU and missed it on GPU. v10 finds it everywhere.

### Where remlax spends its time (v10, warm)

#### CPU, 8 cores

| case | start | preparation | 1st evaluation | optimisation | finish | outputs | solver total |
|----|----|----|----|----|----|----|----|
| DROPS two-stage, fa1 | 2.6 | 0.3 | 2.2 | 1.9 | 0.0 | 0.1 | 7.2 |
| DROPS two-stage, fa3 | 2.6 | 0.3 | 0.3 | 6.3 | 0.0 | 0.1 | 9.6 |
| DROPS one-stage, fa1 | 2.5 | 2.3 | 0.6 | 23.0 | 0.0 | 0.2 | 28.9 |
| DROPS one-stage, us(5) | 2.4 | 2.3 | 0.6 | 12.1 | 0.0 | 0.2 | 17.9 |
| Krause, GBLUP | 3.0 | 1.7 | 0.8 | 3.0 | 0.0 | 0.1 | 8.8 |
| Krause, fa2 × 5 treatments | 2.6 | 2.4 | 28.6 | 960 | 0.0 | 0.6 | 995 |
| Krause, us(5) × 5 treatments | 2.5 | 2.8 | 29.3 | 280 | 0.0 | 0.6 | 316 |

#### A100

| case | start | preparation | 1st evaluation | optimisation | finish | outputs | solver total |
|----|----|----|----|----|----|----|----|
| DROPS two-stage, fa1 | 2.8 | 0.3 | 0.5 | 0.8 | 0.0 | 0.1 | 8.1 |
| DROPS two-stage, fa3 | 2.8 | 0.3 | 0.5 | 2.7 | 0.0 | 0.1 | 10.0 |
| DROPS one-stage, fa1 | 2.7 | 2.9 | 0.8 | 32.2 | 0.0 | 0.3 | 42.6 |
| DROPS one-stage, us(5) | 2.8 | 3.0 | 0.8 | 17.7 | 0.0 | 0.3 | 28.1 |
| Krause, GBLUP | 2.8 | 1.5 | 0.5 | 0.4 | 0.0 | 0.1 | 9.0 |
| Krause, fa2 × 5 treatments | 2.7 | 2.8 | 1.8 | 33.1 | 0.0 | 0.1 | 44.3 |
| Krause, us(5) × 5 treatments | 2.8 | 2.8 | 1.8 | 10.7 | 0.0 | 0.1 | 22.0 |

## Changes in remlax

| change | defect fixed | measured effect |
|----|----|----|
| Step timer (`rx_fit()$chrono`) | no measurement of preparation, compilation or outputs | every step, from the launch by R to reading back the result |
| Starting values on the residuals of the fixed effects, by component and by section | the variance of raw y contains the mean differences between environments | krause_us (600 genotypes): 17 → 7 iterations |
| Exposure divided by the loaded (row, component) pairs | division by n × t: starting values t times too large in long format | drops2: 13 fewer damped iterations |
| AI step trials on the likelihood alone | each rejected trial cost a gradient and an AI matrix | drops2: 42 → 12 s |
| NumPy OpenBLAS on one thread during the MME | the NumPy and SciPy OpenBLAS competed for the cores | evaluation 1,46 → 0,25 s (m = 1428, 8 cores) |
| CUDA plugin not loaded on CPU | 2,0 s of backend initialisation at each launch | start 4,3 → 2,2 s |
| Persistent compilation cache | recompilation at each launch | GPU: first call 3,4 → 0,6 s |
| Vectorised component table, Python probe on demand | one JAX operation per entry of Σ; a Python process launched for nothing | 1,9 → 0,18 s; 1,5 s |
| Full blocks in the dense MME, deferred out of the sparse pattern | 178 million entries grouped by bincount at each evaluation | Krause preparation 18 → 3 s |
| On GPU, C and C⁻¹ stay on the card, blocks grouped by term | 2,8 GB round trip per evaluation; 406 operations to compile | Krause evaluation 6,9 → 1,25 s; drops2 cold 39 → 12 s |
| Last MME evaluation cached | finish and outputs refactorised C at the final point | krause_us on CPU: 444 → 384 s (v6 → v9) |
| An accepted AI step reuses the factorisation of its trial | the full evaluation refactorised C right after the trial on the value | about 10 s less Cholesky per accepted step on Krause on CPU |
| Limited AI: variance ÷10, level parameter at most 0,5 per step | jumps between basins of a multimodal likelihood | drops1: right optimum whatever the number of threads |
| A variance limited three steps in a row is free to reach the floor | Heywood case of an fa: ten more plateaus | krause_fa2: 101 → 56 steps tried |
| Leaving the floor after L-BFGS-B and after the AI | a variance at the floor has a zero gradient on the log scale | drops2 in V form: −2734,67 → −1753,64 |
| METIS and CHOLMOD, the latter forced into our ordering | SuperLU ran out of memory; CHOLMOD could return a factor in its own ordering | mvD7 in 26 GB on CPU; 7 errors of the simulation check → 0 |
| `engine = "auto"`, `optimizer = "auto"`, `theta_init` with NA | choice of engine left to the user; a preliminary call to fix a parameter | MME and AI when the model suits them |

Tried then dropped: a quasi-Newton correction (SR1) of the AI matrix, which increased the iterations everywhere, and a trust-region style ratio test, with no effect on the optima and up to 60 % slower (left as an option, `REMLAX_AI_RAPPORT`). A ÷1000 limit on the variances sped up the fa models but missed the optimum of drops1_us.

## Validation of the final version

| check | result |
|----|----|
| Python tests, CPU | 257 passed, 1 skipped |
| Python tests, GPU (L40S) | 253 passed, 1 skipped |
| testthat | all pass; 1 test skipped for lack of the IGE repository |
| `R CMD check` | Status: OK, with no warning or note |
| simulation check (37 cases, 4 routes, dense REML arbiter in R) | 36 cases passed; mt_g_fa2 fails on its seed 1 on the recovery criterion only (one component at 3,06 standard errors, same exact likelihood on all four routes) and passes on 4 seeds out of 5 |
| CPU/GPU parity (L40S, 7 benchmark bundles) | evaluations within 1e-13 or better, full fits within 2e-10 |

Four tests changed, each justified in its comment. On a data set with a flat and multimodal likelihood, the block engine and the trust-region finish now find a better point than the dense engine. The tests check that they do no worse. The trust-region finish no longer saves evaluations since the starting values were fixed, so its test now compares only the optimum. An analytical threshold of 1e-8 was too strict for GPU sums. It moves to 1e-7.

## Limitations and decisions to make

- The default remains `engine = "dense"` and `optimizer = "lbfgsb"`. Switching to `"auto"` would give the fastest route without tuning, but it would change the optimisation path of your scripts, including those of the IGE pipeline. I did not make this change without your consent.
- The new starting values change the path of every fit. The optima are the same on the simulation check and on the benchmark, but the fits of the IGE pipeline will have to be rerun before any comparison with the old ones.
- Multimodal likelihoods remain a risk. The AI step limits settled every observed case, without guaranteeing the global optimum. `n_restarts` remains the safeguard in case of doubt.
- On CPU, krause_fa2 remains slower than kmmr. The AI takes 19 iterations there against 14, plus reduced step trials early on. This is the next target.
- Starting the Python solver costs about 2,5 s once per R session. It accounts for most of the time on small models.
- asreml-R: the meso licence accepts only a few simultaneous runs. krause_us and krause_fa2 took two hours each there. The asreml fit of mvD7 crashes after the iteration ("double free or corruption").
- The folder `_archive_bancs_20261021` is excluded from the package build, but not from git. It remains to decide whether to keep it.

## Ready to push

Nothing is committed. In the remlax repository, 45 tracked files are modified, and the following ones are new and must be added:

    R/remlax_bayes.R  R/remlax_lag.R  R/remlax_mme.R  R/remlax_serveur.R
    inst/python/remlax/{ai,bayes,blocs,chrono,decoupage,econome,finition,lag,mme,mme_cli,mme_hybride,selinv,serveur}.py
    inst/python/remlax/selinv_c.c
    man/rx_bayes.Rd  man/rx_border.Rd  man/rx_fit_mme.Rd  man/rx_lag.Rd  man/rx_prior.Rd  man/rx_server_stop.Rd
    manual/ (sources; the LaTeX outputs are ignored)
    tests-python/test_{bayes,blocs,decalage,econome,finition,mme_hybride,plancher,poids}.py
    tests/testthat/test-decalage.R  vignettes/structural.Rmd

In remlax-validation, to add: `banc_kmmr/` (cases, kmmr, remlax and asreml drivers, synthesis), `verif_simulations/`, `docs/paysage_logiciels_reml.md` and this report (`docs/rapport_nuit_2026-10-27.md`). The image `container/ige_reml.def` of the IGE repository gains pymetis and cvxopt. The rebuilt image is on genotoul (`/work/user/nsalas/kmmr_banc/ige_reml_v2.sif`).

Frozen versions, sources and checksums: `/work/user/nsalas/kmmr_banc/` on genotoul (`rlib_remlax` for v1, then `rlib_remlax_v2` to `v10`, `VERSION_FIGEE_v*.txt`, `remlax_fige_v*_src.tgz`).
