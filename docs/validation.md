# Validation

Every number on this page was produced on 2026-10-30 by the test suites of this repository and of remlax, with remlax at git revision `bc2af0c`. The campaign ran on the user's SLURM cluster inside the apptainer image `ige_reml_v3.sif`, launched by `campaign/prepare.sh`: every suite on CPU, the CPU/GPU parity on a slice of an A100. The measurements were extracted from the logs (`validation/results/logs_2026-10-30/`) by `validation/parse_logs.py`, and this page is generated from the resulting table (`validation/results/checks_2026-10-30.csv`) by `validation/make_validation_md.py`. No value here was typed by hand.

**381 checks, 364 passed, 17 failed.**

Provenance of the rows:

- re-measured 2026-10-30, ssh:cluster cpu-dedicated and gpu A100 MIG, ige_reml_v3.sif, remlax bc2af0c: 381 checks

Software versions read from the log headers:

- R: 4.6.0 (ssh:cluster cpu-dedicated and gpu A100 MIG, ige_reml_v3.sif, remlax bc2af0c)
- nlme: 3.1.169 (ssh:cluster cpu-dedicated and gpu A100 MIG, ige_reml_v3.sif, remlax bc2af0c)
- lme4: 2.0.1 (ssh:cluster cpu-dedicated and gpu A100 MIG, ige_reml_v3.sif, remlax bc2af0c)
- sommer: 4.4.5 (ssh:cluster cpu-dedicated and gpu A100 MIG, ige_reml_v3.sif, remlax bc2af0c)
- asreml: 4.2.0.480 (ssh:cluster cpu-dedicated and gpu A100 MIG, ige_reml_v3.sif, remlax bc2af0c)
- pbkrtest: 0.5.5 (ssh:cluster cpu-dedicated and gpu A100 MIG, ige_reml_v3.sif, remlax bc2af0c)
- jax: 0.11.1 (ssh:cluster cpu-dedicated and gpu A100 MIG, ige_reml_v3.sif, remlax bc2af0c)

The check labels in the tables below are the strings the test suites printed, quoted verbatim. Suites run before their English translation printed French labels, and those labels are kept as they are. Translating them would give a more readable table and a trace that no longer matches the logs. The section headings and the commentary are in English.

Suites whose log was not available for this run, and therefore absent from the table: `10_tmb_parity.log`, `11_tmb_vs_dense.log`, `12_tmb_prec.log`. See the last section.

## Summary

| Suite | Checks | Compared against | Passed |
|---|---|---|---|
| Internal algebra | 25 | the solver against itself and against closed-form results | 25/25 |
| Unit tests (pytest) | 1 | closed forms and independent numpy references | 1/1 |
| Random sweep | 3 | properties over randomly drawn designs | 2/3 |
| Against lme4 and sommer | 31 | two open-source reference implementations | 31/31 |
| Against nlme | 61 | the reference implementation for correlation structures | 61/61 |
| Against lme4, random regressions, nesting and PEV | 30 | the reference implementation for crossed and nested random effects | 30/30 |
| Against sommer, genomic and multi-trait models | 19 | the reference implementation for genomic relationship matrices | 19/19 |
| Against asreml, structures and inference | 29 | the field reference | 27/29 |
| Against asreml, the rest of the catalogue | 22 | the field reference | 21/22 |
| Against asreml, the remaining structures | 111 | the field reference | 99/111 |
| Against asreml, BLUPs | 4 | the field reference | 4/4 |
| A full model of known truth | 11 | asreml, and the simulated parameters | 11/11 |
| CPU against GPU | 34 | the other backend, on the same input | 33/34 |

## Internal algebra

Gradients against central finite differences, the Hessian against a second differentiation, log-determinant and quadratic form against a direct computation, parameter counts, and the refusal of designs that would make V singular. These need no other software: they are the properties the implementation must have to be an implementation of REML at all.

**1. structures : nombre de parametres, positivite, aller-retour**

| Check | Measured | Verdict |
|---|---|---|
| iid(t=1,r=0) : p=1, symetrique, PSD | min eig 1.70e+00 | pass |
| iid(t=4,r=0) : p=1, symetrique, PSD | min eig 1.70e+00 | pass |
| diag(t=4,r=0) : p=4, symetrique, PSD | min eig 1.70e+00 | pass |
| us(t=2,r=0) : p=3, symetrique, PSD | min eig 1.70e+00 | pass |
| us(t=6,r=0) : p=21, symetrique, PSD | min eig 1.70e+00 | pass |
| fa(t=5,r=1) : p=10, symetrique, PSD | min eig 1.70e+00 | pass |
| fa(t=5,r=3) : p=17, symetrique, PSD | min eig 1.70e+00 | pass |
| fa(t=8,r=4) : p=34, symetrique, PSD | min eig 1.70e+00 | pass |
| us(t=2) represente n'importe quelle matrice PSD | ecart 1.78e-15 | pass |
| us(t=5) represente n'importe quelle matrice PSD | ecart 8.88e-16 | pass |
| us(t=9) represente n'importe quelle matrice PSD | ecart 1.78e-15 | pass |

**2. gradient autodiff contre differences finies centrees**

| Check | Measured | Verdict |
|---|---|---|
| gradient exact (t=1, K=False) | ecart relatif 3.91e-09 | pass |
| gradient exact (t=1, K=True) | ecart relatif 7.82e-10 | pass |
| gradient exact (t=3, K=False) | ecart relatif 1.88e-08 | pass |
| gradient exact (t=3, K=True) | ecart relatif 3.71e-08 | pass |

**3. assemblage de V : symetrie, positivite, coherence Kronecker**

| Check | Measured | Verdict |
|---|---|---|
| V symetrique | - | pass |
| V definie positive | min eig 4.000e-01 | pass |
| V == Z (Sigma (x) K) Z' + R  (voie independante) | ecart max 5.33e-15 | pass |

**4. REML analytique : plan a un facteur, equilibre**

| Check | Measured | Verdict |
|---|---|---|
| a=20 m=6 : sigma2 analytiques retrouvees | ecarts 7.2e-10 / 6.4e-10, decrement 3.2e-17 | pass |
| a=12 m=10 : sigma2 analytiques retrouvees | ecarts 4.0e-11 / 3.7e-11, decrement 1.2e-19 | pass |
| a=40 m=3 : sigma2 analytiques retrouvees | ecarts 4.1e-10 / 3.8e-09, decrement 5.5e-16 | pass |

**5. reproductibilite et invariance au peripherique**

| Check | Measured | Verdict |
|---|---|---|
| deux appels identiques donnent le meme optimum | ecart 0.00e+00 | pass |
| depart perturbe -> meme optimum | ecart 0.00e+00 | pass |

**6. cas degeneres : le solveur refuse plutot que de mentir**

| Check | Measured | Verdict |
|---|---|---|
| structure inconnue rejetee | - | pass |
| X singuliere : erreur explicite | ValueError | pass |

## Unit tests (pytest)

The pytest suite of remlax (tests-python): one test per feature, each against an independent reference (numpy Cholesky REML written out longhand, closed-form correlation matrices, finite differences). Collected here as a single line.

**pytest, tests/python**

| Check | Measured | Verdict |
|---|---|---|
| 274 tests collected by pytest | 273 passed, 1 skipped, 0 failed, 0 errors in 368.87 s | pass |

## Random sweep

Configurations are drawn at random (structure, dimensions, number of traits, seed) and fitted. Nothing is compared to a reference: the sweep asks only whether the solver ever produces something impossible - an exception, a non-finite log-likelihood, a covariance matrix that is not positive semi-definite, a run that does not converge.

**random configuration sweep, seed 0**

| Check | Measured | Verdict |
|---|---|---|
| 200 configurations: no exception, finite logLik, PSD Sigma, converged | clean 197/200; unstable optima 2 (best-of-restart gains 0.0044, 0.2745 in -2logL) | FAIL |

**random configuration sweep, seed 50**

| Check | Measured | Verdict |
|---|---|---|
| 200 configurations: no exception, finite logLik, PSD Sigma, converged | clean 197/200; unstable optima 3 (best-of-restart gains 0.0012, 0.0035, 0.0044 in -2logL) | pass |

**random configuration sweep, seed 100**

| Check | Measured | Verdict |
|---|---|---|
| 200 configurations: no exception, finite logLik, PSD Sigma, converged | clean 199/200; unstable optima 1 (best-of-restart gains 0.0364 in -2logL) | pass |

## Against lme4 and sommer

The models both packages can fit: a single random factor on three layouts, two crossed factors, a genomic relationship matrix, and a bivariate unstructured covariance. Agreement is required on the restricted log-likelihood first, then on the variance components, the fixed effects and the BLUPs.

**1. one random factor, against lme4**

| Check | Measured | Verdict |
|---|---|---|
| a=25 m=5 : identical -2logL | 451.397867339 vs 451.397867339 | pass |
| a=25 m=5 : variance components | s2g 1.361456/1.361456 | pass |
| a=25 m=5 : fixed effects | beta 1.7588 0.80093 | pass |
| a=25 m=5 : BLUPs | max diff 7.47e-09 | pass |
| a=10 m=12 : identical -2logL | 391.341158079 vs 391.341158079 | pass |
| a=10 m=12 : variance components | s2g 2.231905/2.231905 | pass |
| a=10 m=12 : fixed effects | beta 1.70765 0.74394 | pass |
| a=10 m=12 : BLUPs | max diff 5.01e-10 | pass |
| a=60 m=3 : identical -2logL | 653.887073882 vs 653.887073882 | pass |
| a=60 m=3 : variance components | s2g 2.502726/2.502726 | pass |
| a=60 m=3 : fixed effects | beta 2.30711 0.79389 | pass |
| a=60 m=3 : BLUPs | max diff 1.96e-09 | pass |

**2. two crossed factors, against lme4**

| Check | Measured | Verdict |
|---|---|---|
| crossed : identical -2logL | 380.944211187 vs 380.944211187 | pass |
| crossed : the three variances | 1.94207 0.20123 0.90019 | pass |

**3. kinship matrix, against sommer**

| Check | Measured | Verdict |
|---|---|---|
| kinship : variance components vs sommer | s2g 0.99215/0.99211 \| s2e 0.98662/0.98663 | pass |
| kinship : the GRM improves the fit (K vs I) | 515.7869 vs 516.7093 | pass |

**4. two traits, us structure, against sommer**

| Check | Measured | Verdict |
|---|---|---|
| us : 2x2 genetic covariance vs sommer | remlax 1.3169 0.8723 0.8723 1.0382 \| sommer 1.3169 0.8723 0.8723 1.0381 | pass |
| us : 2x2 residual covariance vs sommer | remlax 0.849 0.1019 0.1019 1.2783 \| sommer 0.849 0.1019 0.1019 1.2783 | pass |
| us : genetic correlation recovered | r 0.7460 vs 0.7460 | pass |

**5. nested structures: us must dominate diag, which dominates iid**

| Check | Measured | Verdict |
|---|---|---|
| nesting : -2logL decreasing iid >= diag >= us | iid=1324.4543 diag=1324.0042 us=1299.9284 | pass |

**6. WEIGHTED incidence (the case a formula cannot express)**

| Check | Measured | Verdict |
|---|---|---|
| weighted incidence : finite and defined fit | s2 1.5621 (true 1.21) \| s2e 0.6332 (true 0.64) \| decrement 1.9e-15 | pass |

**7. formula interface: same result as the explicit route**

| Check | Measured | Verdict |
|---|---|---|
| formula ~ g : -2logL identical to lme4 | diff 0.00e+00 | pass |
| formula ~ iid(g)+iid(b) : -2logL identical to lme4 | diff 5.68e-14 | pass |
| formula us(gid) == explicit route rx_term/rx_model | diff 0.00e+00 | pass |

**8. backend switch: the model does not change with the machine**

| Check | Measured | Verdict |
|---|---|---|
| backend auto and cpu give the same fit | cpu vs cpu, diff 0.00e+00 | pass |
| the chosen backend is reported in the result | cpu | pass |

**9. errors: usable messages instead of doubtful results**

| Check | Measured | Verdict |
|---|---|---|
| unknown term rejected | - | pass |
| missing column rejected | - | pass |
| multi-trait structure without trait rejected | - | pass |
| rank-deficient X rejected | - | pass |
| level missing from K rejected | - | pass |

## Against nlme

nlme::gls and nlme::lme on the same data: corAR1, corCAR1, corARMA (AR2 and ARMA(1,1)), corCompSymm, corExp, corGaus, corLin, corSpher, corSymm with varIdent (an unstructured covariance) and varIdent alone (one residual variance per group), then a random intercept with an AR1 residual. A grouped corStruct is written in remlax by declaring the group as the trait of the residual, residual = ~ id(g):ar1(t) with trait = g. Agreement on -2logL, the residual variance, the correlation parameters, the fixed effects and their standard errors.

**1. corAR1: AR1 within each group, against gls**

| Check | Measured | Verdict |
|---|---|---|
| logLik constant : nlme = remlax$logLik | nlme -182.486858 \| remlax -182.486858 \| logLik_asreml -74.052112 (+(n-p)/2 log 2pi = 108.4347) | pass |
| corAR1 : identical -2logL | 364.973716966 vs 364.973716966 | pass |
| corAR1 : residual variance | 1.627385 vs 1.627385 | pass |
| corAR1 : fixed effects | max rel diff 1.45e-09 | pass |
| corAR1 : standard errors of the fixed effects | max rel diff 1.67e-08 | pass |
| corAR1 : rho | 0.5338496 vs 0.5338496 | pass |

**2. corCAR1: continuous AR1 at irregular times, against gls**

| Check | Measured | Verdict |
|---|---|---|
| corCAR1 : identical -2logL | 232.839518984 vs 232.839518984 | pass |
| corCAR1 : residual variance | 1.152159 vs 1.152159 | pass |
| corCAR1 : fixed effects | max rel diff 5.27e-09 | pass |
| corCAR1 : standard errors of the fixed effects | max rel diff 3.40e-08 | pass |
| corCAR1 : phi (correlation per time unit) | 0.7287346 vs 0.7287346 (true 0.7) | pass |

**3. corARMA(p=2): AR2, against gls**

| Check | Measured | Verdict |
|---|---|---|
| ar2 : identical -2logL | 364.591914889 vs 364.591914889 | pass |
| ar2 : residual variance | 1.899831 vs 1.899831 | pass |
| ar2 : fixed effects | max rel diff 1.30e-09 | pass |
| ar2 : standard errors of the fixed effects | max rel diff 1.02e-08 | pass |
| ar2 : phi1 and phi2 | remlax 0.526414 0.145188 \| nlme 0.526414 0.145188 | pass |

**4. corARMA(p=1, q=1): ARMA(1,1), against gls**

| Check | Measured | Verdict |
|---|---|---|
| arma : identical -2logL | 348.507672823 vs 348.507672823 | pass |
| arma : residual variance | 1.461911 vs 1.461911 | pass |
| arma : fixed effects | max rel diff 1.43e-08 | pass |
| arma : standard errors of the fixed effects | max rel diff 2.45e-08 | pass |
| arma : phi and theta (remlax order: theta, phi) | remlax theta 0.389765 phi 0.259672 \| nlme phi 0.259672 theta 0.389765 | pass |

**5. corCompSymm: uniform correlation, against gls**

| Check | Measured | Verdict |
|---|---|---|
| cor : identical -2logL | 365.567018442 vs 365.567018442 | pass |
| cor : residual variance | 1.481197 vs 1.481197 | pass |
| cor : fixed effects | max rel diff 1.81e-10 | pass |
| cor : standard errors of the fixed effects | max rel diff 4.87e-09 | pass |
| cor : uniform rho | 0.2951065 vs 0.2951065 | pass |

**6. corExp: exponential decay in 1D, against gls**

| Check | Measured | Verdict |
|---|---|---|
| exp : identical -2logL | 229.403629489 vs 229.403629489 | pass |
| exp : residual variance | 0.633634 vs 0.633633 | pass |
| exp : fixed effects | max rel diff 6.23e-08 | pass |
| exp : standard errors of the fixed effects | max rel diff 9.81e-07 | pass |
| exp : rho = exp(-1/nlme range) | 0.6429754 vs 0.6429748 (range 2.2642, true 4) | pass |
| exp : the fit moved from its starting point | 8 iterations, rho 0.6430 | pass |

**7. corGaus: Gaussian decay in 1D, against gls**

| Check | Measured | Verdict |
|---|---|---|
| gau : identical -2logL | 235.903503011 vs 235.903503011 | pass |
| gau : residual variance | 1.106024 vs 1.106024 | pass |
| gau : fixed effects | max rel diff 5.05e-13 | pass |
| gau : standard errors of the fixed effects | max rel diff 5.76e-10 | pass |
| gau : rho = exp(-1/range^2) | 0.8895193 vs 0.8895193 (range 2.9226, true 3) | pass |

**8. corLin: truncated linear tent (lvr), against gls**

| Check | Measured | Verdict |
|---|---|---|
| lvr (gls range) : identical -2logL | 188.616692460 vs 188.616692460 | pass |
| lvr (gls range) : residual variance | 0.780188 vs 0.780188 | pass |
| lvr (gls range) : fixed effects | max rel diff 1.33e-14 | pass |
| lvr (gls range) : standard errors of the fixed effects | max rel diff 1.93e-14 | pass |
| lvr : range read = range imposed | 6.02241 vs 6.02241 (true 6) | pass |
| lvr : free remlax optimum no worse than gls | remlax 188.616692 (range 6.0224, 4 restarts) \| best gls over 2 starts 188.616692 (range 6.0224) | pass |

**9. corSpher: spherical in 2D, against gls**

| Check | Measured | Verdict |
|---|---|---|
| sph : identical -2logL | 297.607416797 vs 297.607416797 | pass |
| sph : residual variance | 0.982703 vs 0.982703 | pass |
| sph : fixed effects | max rel diff 4.67e-10 | pass |
| sph : standard errors of the fixed effects | max rel diff 5.43e-09 | pass |
| sph : range | 7.32780 vs 7.32780 (true 8) | pass |

**10. corSymm + varIdent: unstructured covariance (us), against gls**

| Check | Measured | Verdict |
|---|---|---|
| us : identical -2logL | 458.854629436 vs 458.854629445 | pass |
| us : the 10 components of Sigma | max abs diff / max var 5.79e-06 \| diag remlax 0.8898 1.5638 2.2088 0.953 \| nlme 0.8898 1.5638 2.2088 0.953 | pass |
| us : fixed effects | max rel diff 6.47e-07 | pass |

**11. varIdent alone: one residual variance per group (diag), against gls**

| Check | Measured | Verdict |
|---|---|---|
| diag : identical -2logL | 508.562150916 vs 508.562150916 | pass |
| diag : the five variances | remlax 0.206 0.9976 1.8157 4.1048 8.9001 \| nlme 0.206 0.9976 1.8157 4.1048 8.9001 | pass |
| diag : fixed effects | max rel diff 1.81e-07 | pass |

**12. lme: random intercept + within-group corAR1, against lme**

| Check | Measured | Verdict |
|---|---|---|
| lme AR1 : identical -2logL | 348.881810570 vs 348.881810570 | pass |
| lme AR1 : group variance | 1.545201 vs 1.545200 | pass |
| lme AR1 : residual variance | 0.888761 vs 0.888762 | pass |
| lme AR1 : rho | 0.3242477 vs 0.3242487 | pass |
| lme AR1 : fixed effects | max rel diff 3.65e-07 | pass |
| lme AR1 : standard errors of the fixed effects | max rel diff 2.90e-07 | pass |
| lme AR1 : BLUPs | max diff 8.76e-07 | pass |

## Against lme4, random regressions, nesting and PEV

Correlated and independent random slopes (a two-column us and diag term), a nested design, three crossed factors on an unbalanced layout, very unequal group sizes. Besides -2logL, components, fixed effects and their standard errors: the BLUPs, and the prediction error variance checked against the mixed-model-equations formula computed by hand (lme4's condVar follows a different convention, without the uncertainty on beta, and is checked against its own formula).

**1. correlated random slopes (1 + x | g): us with two columns**

| Check | Measured | Verdict |
|---|---|---|
| us slopes : identical -2logL | 708.105175507 vs 708.105175507 | pass |
| us slopes : fixed effects | max rel diff 2.00e-09 | pass |
| us slopes : standard errors of the fixed effects | max rel diff 6.00e-08 | pass |
| us slopes : the three components of Sigma | remlax 1.04933 0.31753 0.51161 \| lme4 1.04933 0.31753 0.51161 | pass |
| us slopes : residual variance | 0.649963 vs 0.649963 | pass |
| us slopes : BLUPs (intercepts and slopes) | max diff 3.11e-08 over 60 values | pass |
| us slopes : remlax PEV = MME formula (with beta) | max rel diff 1.02e-14 \| mean PEV 0.10677 | pass |
| us slopes : lme4 condVar = formula without beta | max rel diff 3.12e-08 \| mean condVar 0.08346 (the two conventions differ) | pass |

**2. independent random slopes (1 + x || g): diag with two columns**

| Check | Measured | Verdict |
|---|---|---|
| diag slopes : identical -2logL | 768.020691545 vs 768.020691545 | pass |
| diag slopes : fixed effects | max rel diff 6.07e-10 | pass |
| diag slopes : standard errors of the fixed effects | max rel diff 7.50e-09 | pass |
| diag slopes : the two variances | remlax 1.18772 0.49438 \| lme4 1.18772 0.49438 | pass |
| diag slopes : covariance zero by construction | 0.0e+00 | pass |

**3. nesting (1 | g/h): two iid terms**

| Check | Measured | Verdict |
|---|---|---|
| nested : identical -2logL | 764.300248907 vs 764.300248907 | pass |
| nested : fixed effects | max rel diff 2.58e-14 | pass |
| nested : standard errors of the fixed effects | max rel diff 4.65e-08 | pass |
| nested : variances of g and of h within g | g 1.53208/1.53208 \| h:g 0.26695/0.26695 \| res 1.05729/1.05729 | pass |
| nested : BLUPs of both terms | max diff 9.93e-09 / 3.52e-08 | pass |

**4. three crossed factors, unbalanced design**

| Check | Measured | Verdict |
|---|---|---|
| 3 crossed : identical -2logL | 1232.919208193 vs 1232.919208193 | pass |
| 3 crossed : fixed effects | max rel diff 1.01e-09 | pass |
| 3 crossed : standard errors of the fixed effects | max rel diff 7.85e-08 | pass |
| 3 crossed : the four variances | remlax 0.81044 0.69056 0.18664 1.02067 \| lme4 0.81044 0.69056 0.18664 1.02067 | pass |
| 3 crossed : BLUPs of the three terms | max diff 1.75e-08 | pass |
| 3 crossed : PEV of the three terms = MME formula | max rel diff 1.48e-14 | pass |
| 3 crossed : PEV of g decreases with replication | cor(PEV, n_rep) = -0.982 | pass |

**5. one factor, very unequal group sizes, two covariates**

| Check | Measured | Verdict |
|---|---|---|
| unbalanced : identical -2logL | 888.893080406 vs 888.893080406 | pass |
| unbalanced : fixed effects | max rel diff 2.66e-11 | pass |
| unbalanced : standard errors of the fixed effects | max rel diff 4.94e-10 | pass |
| unbalanced : variance components | s2g 1.603755/1.603755 \| s2e 1.214474/1.214474 | pass |
| unbalanced : BLUPs | max diff 6.98e-10 (sizes 1 to 12) | pass |

## Against sommer, genomic and multi-trait models

A GRM on one term with BLUPs and PEV; three traits with us genetic and us residual covariance and a GRM; the same with diag and the likelihood ratio between the two; additive plus dominance (two relationship matrices on the same factor); a multi-environment trial with one genetic variance per environment. sommer is tightened to tolParConvLL = 1e-12; its log-likelihood carries a data-dependent constant, so only differences of -2logL are compared.

**1. one GRM, univariate: components, fixed effects, BLUP, PEV**

| Check | Measured | Verdict |
|---|---|---|
| GRM : variance components | s2g 1.350769/1.350769 \| s2e 0.844907/0.844907 | pass |
| GRM : fixed effects | max rel diff 1.61e-13 | pass |
| GRM : standard errors of the fixed effects | max rel diff 1.33e-10 | pass |
| GRM : BLUPs | max diff 2.00e-11 | pass |
| GRM : PEV of the BLUPs | max rel diff 2.60e-10 \| mean PEV 0.17607 | pass |
| GRM : sommer logLik constant (measured) | sommer -42.2549 \| remlax logLik_asreml -135.2526 \| diff 92.9977 ((n-p)/2 log var(y) = 87.6572) | pass |

**2. three traits, genetic us with GRM and residual us**

| Check | Measured | Verdict |
|---|---|---|
| us3 + GRM : the 6 genetic components | max rel diff 7.74e-08 \| diag remlax 1.3404 1.0058 1.0875 \| sommer 1.3404 1.0058 1.0875 | pass |
| us3 + GRM : the 6 residual components | max rel diff 1.94e-08 | pass |
| us3 + GRM : the 3 genetic correlations | remlax 0.5164 0.1735 0.6232 \| sommer 0.5164 0.1735 0.6232 | pass |
| us3 + GRM : fixed effects (3 means) | max rel diff 1.91e-15 | pass |
| us3 + GRM : BLUPs of the 3 traits | max diff 6.65e-08 over 180 values | pass |

**3. three traits, genetic diag with GRM and residual diag; LRT us vs diag**

| Check | Measured | Verdict |
|---|---|---|
| diag3 + GRM : the 3 genetic variances | remlax 1.3399 0.9738 1.0889 \| sommer 1.3399 0.9738 1.0889 | pass |
| diag3 + GRM : the 3 residual variances | remlax 0.7849 1.3983 0.6668 \| sommer 0.7849 1.3983 0.6668 | pass |
| LRT us vs diag : same difference in -2logL | remlax 81.06687 \| sommer 81.06687 (6 df) | pass |

**4. additive + dominance: two kinship matrices on the same factor**

| Check | Measured | Verdict |
|---|---|---|
| A + D : the three components | remlax 0.93906 0.76455 0.78469 \| sommer 0.93906 0.76455 0.78469 | pass |
| A + D : additive and dominance BLUPs | max diff 6.47e-08 | pass |

**5. multi-environment: one genetic variance per environment with GRM**

| Check | Measured | Verdict |
|---|---|---|
| MET : the 4 genetic variances | remlax 0.5255 0.8048 1.1331 0.9474 \| sommer 0.5255 0.8048 1.1331 0.9474 | pass |
| MET : the 4 residual variances | remlax 0.7188 1.1889 1.5606 0.3662 \| sommer 0.7188 1.1889 1.5606 0.3662 | pass |
| MET : fixed effects | max rel diff 7.46e-11 | pass |

## Against asreml, structures and inference

The structures a plant breeder actually writes: a one-dimensional AR1, a separable AR1xAR1 field, two-dimensional splines, and a covariance shared between two terms through str(). Then the inference layer: heritability and its delta-method standard error, Wald tests, and the equivalences between saturated parameterisations.

**1. one-dimensional AR1**

| Check | Measured | Verdict |
|---|---|---|
| AR1 : rho | 0.6276741 vs 0.6276741 | pass |
| AR1 : variance | 1.1690526 vs 1.1690526 | pass |
| AR1 : residual | - | pass |
| AR1 : logLik (asreml convention) | -111.075133372 vs -111.075133372 | pass |

**2. separable AR1 x AR1 (spatial field)**

| Check | Measured | Verdict |
|---|---|---|
| AR1xAR1 : field variance | 1.2088279 vs 1.2088279 | pass |
| AR1xAR1 : row rho | 0.5940549 vs 0.5940549 | pass |
| AR1xAR1 : column rho | 0.3397100 vs 0.3397100 | pass |
| AR1xAR1 : residual | - | pass |
| AR1xAR1 : logLik | -154.818019348 vs -154.818019348 | pass |

**3. 2D spline (SAME basis matrices on both sides)**

| Check | Measured | Verdict |
|---|---|---|
| spline 2D : variance spl_fx | 1.063175e+01 vs 1.063180e+01 | pass |
| spline 2D : variance spl_fy | 1.552303e+01 vs 1.545232e+01 | FAIL |
| spline 2D : variance spl_fx_y | 6.737911e-11 vs 5.425611e-07 | pass |
| spline 2D : variance spl_x_fy | 3.775135e-11 vs 8.655429e-11 | pass |
| spline 2D : variance spl_fx_fy | 3.775135e-11 vs 3.454338e-11 | pass |
| spline 2D : residual | - | pass |
| spline 2D : logLik | -5.014079566 vs -5.014093202 | FAIL |

**4. str(): covariance SHARED between two terms**

| Check | Measured | Verdict |
|---|---|---|
| str : var term 1 | 1.1103369 vs 1.1103369 | pass |
| str : covariance | 0.5137383 vs 0.5137383 | pass |
| str : var term 2 | 0.6966671 vs 0.6966671 | pass |
| str : residual | - | pass |
| str : logLik | -133.657051701 vs -133.657051701 | pass |

**5. vpredict and Wald**

| Check | Measured | Verdict |
|---|---|---|
| vpredict : estimate of h2 | 0.6685333 vs 0.6685333 | pass |
| vpredict : STANDARD ERROR of h2 (delta method) | 0.0500972 vs 0.0500998 | pass |
| Wald : df of the trt term | 2 | pass |
| Wald : p-value of the last term (seq == cond) | 8.776e-09 vs 8.776e-09 | pass |

**6. SATURATED parametrisations: chol(t-1), ante(t-1), fa(t), rr(t)**

| Check | Measured | Verdict |
|---|---|---|
| chol(gid, rank=3)  == us | diff 4.55e-13 | pass |
| ante(gid, rank=3)  == us | diff 2.27e-13 | pass |
| fa(gid, rank=4)    == us | diff 5.68e-13 | pass |
| rr(gid, rank=4)    == us | diff 7.96e-13 | pass |

## Against asreml, the rest of the catalogue

The structures that were identified by matching asreml's output rather than from documentation - lvr, the anisotropic Matern, user-written structures, sectioned residuals - plus prediction and the Kenward-Roger denominator degrees of freedom against pbkrtest.

**1. lvr: truncated tent max(0, 1 - d/phi)**

| Check | Measured | Verdict |
|---|---|---|
| lvr : logLik = independent R profile | remlax -2.6786173 \| R -1.9143976 \| range 1.5055 vs 5.8896 | FAIL |
| lvr : asreml does not reach a better optimum | asreml -2.819575 <= remlax -2.678617 | pass |

**2. mtrn: anisotropic Matern, parameters FIXED on both sides**

| Check | Measured | Verdict |
|---|---|---|
| mtrn isotropic nu=0.8 : logLik = besselK reference | remlax 5.1067110 \| R 5.1067110 | pass |
| mtrn isotropic nu=0.8 : logLik = asreml | asreml 5.1067110 \| remlax 5.1067110 | pass |
| mtrn nu=1.5 delta=2 : logLik = besselK reference | remlax -6.9636274 \| R -6.9636274 | pass |
| mtrn nu=1.5 delta=2 : logLik = asreml | asreml -6.9636274 \| remlax -6.9636274 | pass |
| mtrn delta=2 alpha=.6 : logLik = besselK reference | remlax 4.9300437 \| R 4.9300437 | pass |
| mtrn nu=0.5 lambda=1 : logLik = besselK reference | remlax -1.0246717 \| R -1.0246717 | pass |
| mtrn nu=0.5 lambda=1 : logLik = asreml | asreml -1.0246717 \| remlax -1.0246717 | pass |

**3. mtrn: the range is indeed ESTIMATED when freed**

| Check | Measured | Verdict |
|---|---|---|
| mtrn : estimated phi = argmax of the R profile | phi 0.62086 vs 0.62086 (true 2.5) \| logLik 4.9738738 vs 4.9738738 | pass |

**5. dsum: one residual per section**

| Check | Measured | Verdict |
|---|---|---|
| dsum : two sections improve the likelihood | 1 section -146.2015 -> 2 sections -127.2089 (+18.99) | pass |
| dsum : the two residual variances are separated | 0.1388 / 1.2999 (true 0.25 and 1.96) | pass |
| dsum : logLik = asreml | asreml -39.909705 \| remlax -39.909705 | pass |

**5b. dsum: DIFFERENT structures by section**

| Check | Measured | Verdict |
|---|---|---|
| dsum : three sections, two AR1 and one iid | 3 sections \| 6 parameters \| rho on A+B | pass |

**6. predict: adjusted means and their standard errors**

| Check | Measured | Verdict |
|---|---|---|
| predict : values = asreml | max diff 1.31e-13 | pass |
| predict : standard errors = asreml | max diff 1.34e-10 | pass |
| predict : standard error of differences available | mean sed 0.18922 | pass |

**7. Kenward-Roger vs pbkrtest**

| Check | Measured | Verdict |
|---|---|---|
| K-R : denominator df = pbkrtest | remlax 21.76897 \| pbkrtest 21.76896 | pass |
| K-R : F statistic = pbkrtest | remlax 2.790148 \| pbkrtest 2.790148 | pass |
| K-R : adjusted standard errors = pbkrtest | max diff 1.31e-07 | pass |
| K-R : the adjustment INFLATES the standard errors | mean ratio 1.03753 | pass |

**8. predict with the random part (BLUP + prediction error)**

| Check | Measured | Verdict |
|---|---|---|
| predict : the BLUP enters the mean by block | standard deviation with BLUP 0.69886, without 0.00e+00 | pass |

## Against asreml, the remaining structures

Grouped residual series ar2, ar3, ma1, ma2, arma, sar, cor, corb and corg; the two-dimensional metric kernels iexp, igau, ieuc, aexp, agau, sph and cir on an irregular field; the non-saturated multi-trait structures fa(1), rr(1), corh, ante(1), chol(1) and diag on four traits with a GRM; the separable product id x ar1 x ar1 with shared correlations; and the multi-trait spatial residual us(trait):ar1(row):ar1(col).

**A1. ar2: grouped residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| ar2 : logLik (asreml convention) | -72.335196521 vs -72.335196521 | pass |
| ar2 : residual variance | 1.675628 vs 1.675628 | pass |
| ar2 : correlation parameters (as a set) | remlax 0.5325 0.19233 \| asreml 0.5325 0.19233 | pass |
| ar2 : remlax no worse than asreml | -72.335197 vs -72.335197 (diff -2.84e-14) | pass |
| ar2 : fixed effects (up to order) | remlax 0.6225 0.4419 \| asreml 0.6225 0.4419 | pass |

**A2. ar3: grouped residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| ar3 : logLik (asreml convention) | -75.186557795 vs -75.186557795 | pass |
| ar3 : residual variance | 1.658159 vs 1.658159 | pass |
| ar3 : correlation parameters (as a set) | remlax 0.4335 0.10803 0.21949 \| asreml 0.4335 0.10803 0.21949 | pass |
| ar3 : remlax no worse than asreml | -75.186558 vs -75.186558 (diff +1.42e-14) | pass |
| ar3 : fixed effects (up to order) | remlax 1.0652 0.4518 \| asreml 1.0652 0.4518 | pass |

**A3. ma1: grouped residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| ma1 : logLik (asreml convention) | -71.441165627 vs -71.441165627 | pass |
| ma1 : residual variance | 1.173080 vs 1.173080 | pass |
| ma1 : correlation parameters (as a set) | remlax -0.52845 \| asreml -0.52845 | pass |
| ma1 : remlax no worse than asreml | -71.441166 vs -71.441166 (diff +0.00e+00) | pass |
| ma1 : fixed effects (up to order) | remlax 1.0329 0.6044 \| asreml 1.0329 0.6044 | pass |

**A4. ma2: grouped residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| ma2 : logLik (asreml convention) | -72.524782986 vs -72.558053561 | FAIL |
| ma2 : residual variance | 1.400230 vs 1.369831 | FAIL |
| ma2 : correlation parameters (as a set) | remlax -0.58011 -0.44931 \| asreml -0.57924 -0.44714 | FAIL |
| ma2 : remlax no worse than asreml | -72.524783 vs -72.558054 (diff +3.33e-02) | pass |
| ma2 : fixed effects (up to order) | remlax 0.8983 0.4977 \| asreml 0.8994 0.4962 | FAIL |

**A5. arma: grouped residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| arma : logLik (asreml convention) | -73.124116717 vs -73.124116717 | pass |
| arma : residual variance | 1.376202 vs 1.376202 | pass |
| arma : correlation parameters (as a set) | remlax 0.28594 0.35587 \| asreml 0.35587 0.28594 | pass |
| arma : remlax no worse than asreml | -73.124117 vs -73.124117 (diff -8.53e-14) | pass |
| arma : fixed effects (up to order) | remlax 0.8109 0.3658 \| asreml 0.8109 0.3658 | pass |

**A6. sar: grouped residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| sar : logLik (asreml convention) | -70.706817237 vs -70.706817237 | pass |
| sar : residual variance | 1.369386 vs 1.369386 | pass |
| sar : correlation parameters (as a set) | remlax 0.63555 \| asreml 0.63555 | pass |
| sar : remlax no worse than asreml | -70.706817 vs -70.706817 (diff -8.53e-14) | pass |
| sar : fixed effects (up to order) | remlax 1.0847 0.3773 \| asreml 1.0847 0.3773 | pass |

**A7. cor: grouped residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| cor : logLik (asreml convention) | -96.140426816 vs -96.140426816 | pass |
| cor : residual variance | 1.515622 vs 1.515622 | pass |
| cor : correlation parameters (as a set) | remlax 0.21775 \| asreml 0.21775 | pass |
| cor : remlax no worse than asreml | -96.140427 vs -96.140427 (diff +1.42e-14) | pass |
| cor : fixed effects (up to order) | remlax 0.9673 0.6452 \| asreml 0.9673 0.6452 | pass |

**A8. corb: grouped residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| corb : logLik (asreml convention) | -69.307901067 vs -69.307901067 | pass |
| corb : residual variance | 1.092528 vs 1.092528 | pass |
| corb : correlation parameters (as a set) | remlax 0.41316 0.14577 \| asreml 0.41316 0.14577 | pass |
| corb : remlax no worse than asreml | -69.307901 vs -69.307901 (diff +1.42e-14) | pass |
| corb : fixed effects (up to order) | remlax 0.9803 0.4266 \| asreml 0.9803 0.4266 | pass |

**A9. corg: grouped residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| corg : logLik (asreml convention) | -97.912513106 vs -97.912513106 | pass |
| corg : residual variance | 0.905272 vs 0.905272 | pass |
| corg : correlation parameters (as a set) | remlax 0.38785 0.16896 0.22701 -0.07416 0.19222 0.43027 \| asreml 0.38785 0.16896 0.22701 -0.07416 0.19222 0.43027 | pass |
| corg : remlax no worse than asreml | -97.912513 vs -97.912513 (diff +7.11e-14) | pass |
| corg : fixed effects (up to order) | remlax 0.9558 0.4911 \| asreml 0.9558 0.4911 | pass |

**B1. iexp: 2D metric kernel, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| iexp : logLik (asreml convention) | -4.778255069 vs -4.778255069 | pass |
| iexp : variance | 1.476882 vs 1.476882 | pass |
| iexp : kernel parameters (as a set) | remlax 0.87186 \| asreml 0.87186 | pass |
| iexp : remlax no worse than asreml | -4.778255 vs -4.778255 (diff +5.33e-14) | pass |
| iexp : same logLik at the asreml parameters | -4.778255069 vs -4.778255069 (remlax evaluated at the asreml parameters) | pass |

**B2. igau: 2D metric kernel, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| igau : logLik (asreml convention) | -3.293847858 vs -3.293847858 | pass |
| igau : variance | 1.287447 vs 1.287447 | pass |
| igau : kernel parameters (as a set) | remlax 0.80114 \| asreml 0.80114 | pass |
| igau : remlax no worse than asreml | -3.293848 vs -3.293848 (diff +3.90e-12) | pass |
| igau : same logLik at the asreml parameters | -3.293847858 vs -3.293847858 (remlax evaluated at the asreml parameters) | pass |

**B3. ieuc: 2D metric kernel, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| ieuc : logLik (asreml convention) | -4.016578104 vs -4.016578104 | pass |
| ieuc : variance | 0.769838 vs 0.769838 | pass |
| ieuc : kernel parameters (as a set) | remlax 0.70585 \| asreml 0.70585 | pass |
| ieuc : remlax no worse than asreml | -4.016578 vs -4.016578 (diff +1.42e-14) | pass |
| ieuc : same logLik at the asreml parameters | -4.016578104 vs -4.016578104 (remlax evaluated at the asreml parameters) | pass |

**B4. aexp: 2D metric kernel, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| aexp : logLik (asreml convention) | -37.151545485 vs -37.151545489 | pass |
| aexp : variance | 1.424754 vs 1.424748 | pass |
| aexp : kernel parameters (as a set) | remlax 0.75991 0.65095 \| asreml 0.75991 0.65096 | pass |
| aexp : remlax no worse than asreml | -37.151545 vs -37.151545 (diff +4.57e-09) | pass |
| aexp : same logLik at the asreml parameters | -37.151545486 vs -37.151545489 (remlax evaluated at the asreml parameters) | pass |

**B5. agau: 2D metric kernel, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| agau : logLik (asreml convention) | 17.316795825 vs 17.316795825 | pass |
| agau : variance | 0.605183 vs 0.605183 | pass |
| agau : kernel parameters (as a set) | remlax 0.94404 0.24205 \| asreml 0.94404 0.24205 | pass |
| agau : remlax no worse than asreml | 17.316796 vs 17.316796 (diff +1.42e-10) | pass |
| agau : same logLik at the asreml parameters | 17.316795825 vs 17.316795825 (remlax evaluated at the asreml parameters) | pass |

**B6. sph: 2D metric kernel, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| sph : logLik (asreml convention) | 3.823231014 vs 3.753094938 | FAIL |
| sph : variance | 1.024529 vs 1.006357 | FAIL |
| sph : kernel parameters (as a set) | remlax 8.6514 \| asreml 8.75447 | FAIL |
| sph : remlax no worse than asreml | 3.823231 vs 3.753095 (diff +7.01e-02) | pass |
| sph : same logLik at the asreml parameters | 3.783476190 vs 3.753094938 (remlax evaluated at the asreml parameters) | FAIL |

**B7. cir: 2D metric kernel, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| cir : logLik (asreml convention) | 13.875897073 vs 10.847517451 | FAIL |
| cir : variance | 0.846601 vs 0.623517 | FAIL |
| cir : kernel parameters (as a set) | remlax 7.40768 \| asreml 4.98345 | FAIL |
| cir : remlax no worse than asreml | 13.875897 vs 10.847517 (diff +3.03e+00) | pass |
| cir : same logLik at the asreml parameters | 10.844976334 vs 10.847517451 (remlax evaluated at the asreml parameters) | FAIL |

**C1. us(trait) x GRM, 4 traits, diag residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| us : logLik (asreml convention) | -475.359451222 vs -475.359451222 | pass |
| us : remlax no worse than asreml | -475.359451 vs -475.359451 (diff +7.84e-12) | pass |
| us : genetic Sigma | max rel diff 4.23e-08 \| diag remlax 1.209 2.6058 0.711 1.3245 \| asreml 1.209 2.6058 0.711 1.3245 | pass |
| us : the 4 residual variances | remlax 0.798 0.8956 0.5819 0.9622 \| asreml 0.798 0.8956 0.5819 0.9622 | pass |

**C2. diag(trait) x GRM, 4 traits, diag residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| diag : logLik (asreml convention) | -493.193233887 vs -493.193233887 | pass |
| diag : remlax no worse than asreml | -493.193234 vs -493.193234 (diff -9.09e-13) | pass |
| diag : genetic Sigma | max rel diff 2.49e-08 \| diag remlax 1.2052 2.6244 0.7152 1.3477 \| asreml 1.2052 2.6244 0.7152 1.3477 | pass |
| diag : the 4 residual variances | remlax 0.7998 0.8921 0.5804 0.9565 \| asreml 0.7998 0.8921 0.5804 0.9565 | pass |

**C3. fa1(trait) x GRM, 4 traits, diag residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| fa1 : logLik (asreml convention) | -480.160993199 vs -480.160993199 | pass |
| fa1 : remlax no worse than asreml | -480.160993 vs -480.160993 (diff -5.50e-11) | pass |
| fa1 : the 4 residual variances | remlax 0.7967 0.8955 0.5807 0.9543 \| asreml 0.7967 0.8955 0.5807 0.9543 | pass |

**C4. rr1(trait) x GRM, 4 traits, diag residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| rr1 : logLik (asreml convention) | -564.439985070 vs -564.440028141 | pass |
| rr1 : remlax no worse than asreml | -564.439985 vs -564.440028 (diff +4.31e-05) | pass |
| rr1 : the 4 residual variances | remlax 1.9695 1.0314 1.3332 1.5949 \| asreml 1.9693 1.0321 1.3332 1.5944 | pass |

**C5. corh(trait) x GRM, 4 traits, diag residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| corh : logLik (asreml convention) | -492.756261237 vs -492.756261237 | pass |
| corh : remlax no worse than asreml | -492.756261 vs -492.756261 (diff +3.18e-12) | pass |
| corh : the 4 residual variances | remlax 0.7972 0.894 0.5795 0.9562 \| asreml 0.7972 0.894 0.5795 0.9562 | pass |

**C6. ante1(trait) x GRM, 4 traits, diag residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| ante1 : logLik (asreml convention) | -491.875918525 vs -491.875918525 | pass |
| ante1 : remlax no worse than asreml | -491.875919 vs -491.875919 (diff +3.58e-11) | pass |
| ante1 : the 4 residual variances | remlax 0.7999 0.8922 0.5821 0.9573 \| asreml 0.7999 0.8922 0.5821 0.9573 | pass |

**C7. chol1(trait) x GRM, 4 traits, diag residual, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| chol1 : logLik (asreml convention) | -492.112977719 vs -492.112977719 | pass |
| chol1 : remlax no worse than asreml | -492.112978 vs -492.112978 (diff -5.05e-11) | pass |
| chol1 : the 4 residual variances | remlax 0.7999 0.8922 0.581 0.9575 \| asreml 0.7999 0.8922 0.581 0.9575 | pass |

**D. sep: id(block) x ar1(col) x ar1(row), shared correlations, against asreml**

| Check | Measured | Verdict |
|---|---|---|
| sep : logLik (asreml convention) | -349.658684512 vs -349.658684512 | pass |
| sep : field variance | 1.742834 vs 1.742834 | pass |
| sep : column rho and row rho | remlax 0.30873 0.59191 \| asreml 0.30873 0.59191 | pass |
| sep : residual | 0.772377 vs 0.772377 | pass |

**E. spatial multi-trait residual us(trait):ar1(row):ar1(col), against asreml**

| Check | Measured | Verdict |
|---|---|---|
| us:ar1:ar1 : logLik (asreml convention) | -97.115201240 vs -97.115201240 | pass |
| us:ar1:ar1 : row and column rho | remlax 0.50089 0.36206 \| asreml 0.50089 0.36206 | pass |
| us:ar1:ar1 : 2x2 residual Sigma | remlax 1.1591 0.6273 0.6273 0.8517 \| asreml 1.1591 0.6273 0.6273 0.8517 | pass |
| us:ar1:ar1 : 2x2 genetic Sigma | remlax 0.5743 0.1596 0.1596 0.4802 \| asreml 0.5743 0.1596 0.1596 0.4802 | pass |

## Against asreml, BLUPs

BLUPs and their ordering on a one-factor model, level by level.

**1. fit on both sides, time measured**

| Check | Measured | Verdict |
|---|---|---|
| logLik in the asreml convention | -260.296848172 vs -260.296848182 | pass |

**3. BLUP: do both sides predict the same effects?**

| Check | Measured | Verdict |
|---|---|---|
| BLUP : maximum difference | max\|diff\| = 3.155e-06 | pass |
| BLUP : correlation | r = 1.000000000000 | pass |

**4. fixed effects and their variance**

| Check | Measured | Verdict |
|---|---|---|
| fixed effect (intercept) | 5.110533510 vs 5.110533510 | pass |

## A full model of known truth

The complete direct and indirect genetic effects model: a weighted neighbourhood incidence, a covariance shared between the direct and indirect terms, a separable field and a nugget - eight variance parameters - fitted on simulated data whose true values are known, by remlax and by asreml.

| Check | Measured | Verdict |
|---|---|---|
| the sign of cov(DGE, IGE intra) is recovered | -0.4221 vs -0.35 | pass |
| the fit is at an optimum | decrement 2.4e-24 | pass |
| remlax is no worse than asreml | diff +2.05e-02 in favour of remlax | pass |
| var_DGE : agreement with asreml | 1.09524 vs 1.09592 (SE 0.189) | pass |
| cov_DGE_IGE : agreement with asreml | -0.42211 vs -0.41939 (SE 0.180) | pass |
| var_IGEintra : agreement with asreml | 0.16268 vs 0.17627 (SE 0.354) | pass |
| var_IGEinter : agreement with asreml | 0.54369 vs 0.53748 (SE 0.341) | pass |
| var_field : agreement with asreml | 0.61009 vs 0.61022 (SE 0.164) | pass |
| rho_row : agreement with asreml | 0.53375 vs 0.53215 (SE 0.131) | pass |
| rho_column : agreement with asreml | 0.11931 vs 0.11850 (SE 0.111) | pass |
| nugget : agreement with asreml | 0.39824 vs 0.39694 (SE 0.156) | pass |

## CPU against GPU

Every design of the structure catalogue is serialised once, then fitted twice in the same process: on the CPU and on a slice of an A100. The model must not depend on the machine that fits it.

**structure catalogue, 34 serialised designs**

| Check | Measured | Verdict |
|---|---|---|
| design 'aexp': same -2logL on CPU and on an A100 MIG slice | -2logL 560.16430238 vs 560.16430238 \| rel gap 2.03e-16 \| max\|dtheta\| 2.78e-16 \| 1.1 s CPU vs 8.6 s MIG | pass |
| design 'ante1_mv': same -2logL on CPU and on an A100 MIG slice | -2logL 1284.22057669 vs 1284.22057669 \| rel gap 1.77e-16 \| max\|dtheta\| 2.67e-08 \| 1.5 s CPU vs 5.1 s MIG | pass |
| design 'ar1': same -2logL on CPU and on an A100 MIG slice | -2logL 585.15663362 vs 585.15663362 \| rel gap 0.00e+00 \| max\|dtheta\| 2.78e-16 \| 0.4 s CPU vs 2.3 s MIG | pass |
| design 'ar1ar1_res': same -2logL on CPU and on an A100 MIG slice | -2logL 513.75395410 vs 513.75395410 \| rel gap 0.00e+00 \| max\|dtheta\| 4.37e-16 \| 0.5 s CPU vs 2.7 s MIG | pass |
| design 'ar2': same -2logL on CPU and on an A100 MIG slice | -2logL 583.38706024 vs 583.38706024 \| rel gap 0.00e+00 \| max\|dtheta\| 6.66e-16 \| 0.8 s CPU vs 1.4 s MIG | pass |
| design 'ar3': same -2logL on CPU and on an A100 MIG slice | -2logL 581.02165146 vs 581.02165146 \| rel gap 0.00e+00 \| max\|dtheta\| 3.22e-15 \| 1.0 s CPU vs 1.5 s MIG | pass |
| design 'arma': same -2logL on CPU and on an A100 MIG slice | -2logL 583.22954846 vs 583.22954846 \| rel gap 0.00e+00 \| max\|dtheta\| 1.04e-04 \| 0.6 s CPU vs 0.8 s MIG | pass |
| design 'chol2_mv': same -2logL on CPU and on an A100 MIG slice | -2logL 1281.59614340 vs 1281.59614340 \| rel gap 1.77e-16 \| max\|dtheta\| 3.00e-14 \| 1.1 s CPU vs 1.2 s MIG | pass |
| design 'cir': same -2logL on CPU and on an A100 MIG slice | -2logL 591.99572228 vs 591.99572228 \| rel gap 0.00e+00 \| max\|dtheta\| 5.07e-09 \| 0.5 s CPU vs 0.8 s MIG | pass |
| design 'cor': same -2logL on CPU and on an A100 MIG slice | -2logL 585.55424610 vs 585.55424610 \| rel gap 1.94e-16 \| max\|dtheta\| 5.70e-15 \| 0.4 s CPU vs 0.6 s MIG | pass |
| design 'corb2': same -2logL on CPU and on an A100 MIG slice | -2logL 583.42821102 vs 583.42821102 \| rel gap 3.90e-16 \| max\|dtheta\| 3.22e-15 \| 0.4 s CPU vs 0.6 s MIG | pass |
| design 'corg': same -2logL on CPU and on an A100 MIG slice | -2logL 601.62640596 vs 601.62640596 \| rel gap 1.89e-16 \| max\|dtheta\| 8.48e-10 \| 0.8 s CPU vs 2.8 s MIG | pass |
| design 'corh_mv': same -2logL on CPU and on an A100 MIG slice | -2logL 1281.60150032 vs 1281.60150032 \| rel gap 1.77e-16 \| max\|dtheta\| 1.72e-15 \| 1.1 s CPU vs 1.1 s MIG | pass |
| design 'dsum': same -2logL on CPU and on an A100 MIG slice | -2logL 515.18686079 vs 515.18686079 \| rel gap 2.21e-16 \| max\|dtheta\| 1.18e-16 \| 0.3 s CPU vs 0.5 s MIG | pass |
| design 'dsum_ar1': same -2logL on CPU and on an A100 MIG slice | -2logL 52.15394381 vs 52.15394381 \| rel gap 6.81e-16 \| max\|dtheta\| 1.55e-15 \| 0.5 s CPU vs 0.9 s MIG | pass |
| design 'expo': same -2logL on CPU and on an A100 MIG slice | -2logL 585.55424610 vs 585.55424610 \| rel gap 0.00e+00 \| max\|dtheta\| 6.37e-15 \| 0.5 s CPU vs 1.2 s MIG | pass |
| design 'fa2_mv': same -2logL on CPU and on an A100 MIG slice | -2logL 1281.59614340 vs 1281.59614340 \| rel gap 0.00e+00 \| max\|dtheta\| 1.21e-08 \| 1.5 s CPU vs 3.7 s MIG | pass |
| design 'gau': same -2logL on CPU and on an A100 MIG slice | -2logL 585.55424610 vs 585.55424610 \| rel gap 0.00e+00 \| max\|dtheta\| 2.07e-14 \| 0.5 s CPU vs 0.6 s MIG | pass |
| design 'ieuc': same -2logL on CPU and on an A100 MIG slice | -2logL 590.00207600 vs 590.00207600 \| rel gap 0.00e+00 \| max\|dtheta\| 9.19e-07 \| 1.0 s CPU vs 0.9 s MIG | pass |
| design 'iexp': same -2logL on CPU and on an A100 MIG slice | -2logL 578.64147976 vs 578.64147976 \| rel gap 1.96e-16 \| max\|dtheta\| 5.55e-16 \| 0.4 s CPU vs 0.6 s MIG | pass |
| design 'iid': same -2logL on CPU and on an A100 MIG slice | -2logL 498.81465569 vs 498.81465569 \| rel gap 5.70e-16 \| max\|dtheta\| 7.49e-16 \| 0.3 s CPU vs 1.1 s MIG | pass |
| design 'lvr': same -2logL on CPU and on an A100 MIG slice | -2logL 585.55424610 vs 585.55424610 \| rel gap 3.88e-16 \| max\|dtheta\| 7.86e-12 \| 0.5 s CPU vs 0.7 s MIG | pass |
| design 'ma1': same -2logL on CPU and on an A100 MIG slice | -2logL 583.30849572 vs 583.30849572 \| rel gap 3.90e-16 \| max\|dtheta\| 1.13e-06 \| 0.6 s CPU vs 0.8 s MIG | pass |
| design 'mtrn_iso': same -2logL on CPU and on an A100 MIG slice | -2logL 589.40743274 vs 589.40743274 \| rel gap 1.93e-16 \| max\|dtheta\| 1.60e-12 \| 6.7 s CPU vs 288.1 s MIG | pass |
| design 'own_2par': same -2logL on CPU and on an A100 MIG slice | -2logL 607.25006661 vs 607.25006661 \| rel gap 0.00e+00 \| max\|dtheta\| 0.00e+00 \| 0.5 s CPU vs 0.8 s MIG | pass |
| design 'own_exp': same -2logL on CPU and on an A100 MIG slice | -2logL 585.55424610 vs 585.55424610 \| rel gap 0.00e+00 \| max\|dtheta\| 4.43e-07 \| 0.6 s CPU vs 0.8 s MIG | pass |
| design 'rr2_mv': same -2logL on CPU and on an A100 MIG slice | -2logL 1343.60530125 vs 1343.60530125 \| rel gap 1.69e-16 \| max\|dtheta\| 2.34e-08 \| 0.9 s CPU vs 3.2 s MIG | pass |
| design 'sar': same -2logL on CPU and on an A100 MIG slice | -2logL 585.04086512 vs 585.04086512 \| rel gap 1.94e-16 \| max\|dtheta\| 2.52e-08 \| 0.7 s CPU vs 0.9 s MIG | pass |
| design 'sph': same -2logL on CPU and on an A100 MIG slice | -2logL 587.42826916 vs 587.42826916 \| rel gap 1.94e-16 \| max\|dtheta\| 5.07e-09 \| 1.0 s CPU vs 1.0 s MIG | pass |
| design 'str_2': same -2logL on CPU and on an A100 MIG slice | -2logL 512.47455420 vs 512.47455420 \| rel gap 0.00e+00 \| max\|dtheta\| 2.64e-16 \| 0.7 s CPU vs 3.8 s MIG | pass |
| design 'us_mv': same -2logL on CPU and on an A100 MIG slice | -2logL 1281.59614340 vs 1281.59614340 \| rel gap 0.00e+00 \| max\|dtheta\| 9.30e-16 \| 1.2 s CPU vs 1.0 s MIG | pass |
| design 'us_res_mv': same -2logL on CPU and on an A100 MIG slice | -2logL 1395.95209251 vs 1395.95209251 \| rel gap 0.00e+00 \| max\|dtheta\| 8.88e-16 \| 0.8 s CPU vs 2.6 s MIG | pass |
| design 'vm_kin': same -2logL on CPU and on an A100 MIG slice | -2logL 494.64877180 vs 494.64877180 \| rel gap 5.75e-16 \| max\|dtheta\| 3.05e-09 \| 0.3 s CPU vs 0.5 s MIG | pass |
| design 'mtrn_aniso': same -2logL on CPU and on an A100 MIG slice | fit error: ValueError: -2logL or its gradient is not finite at the starting point (value 629.5077869074553, 2 non-finite gradient component(s)). The optimiser cannot start. Check the kernel or supply theta_init. | FAIL |

## The failed checks

A failed check stays failed on purpose: each one is a documented disagreement, not a reason to widen a tolerance.

**ma2, sph and cir against asreml (12 checks, as in 2026-09).** `tests/R/diag_asreml3_ecarts.R` settles them (its log sits with the others). A dense REML arbiter written in R finds the remlax value to the last digit, and asreml held at the parameters of remlax returns the remlax value. For sph and cir, asreml reports converge = FALSE (it oscillates between hills of a multimodal surface). For ma2, asreml declares convergence while its parameters still change by several per cent per iteration, and its two correlation parameters are not the MA coefficients of the standard form.

**2D P-spline against asreml (2 checks).** `rx_spl2d()` builds the PS-ANOVA of SpATS, five smooth components, and asreml receives the same five basis matrices. remlax reaches a restricted log-likelihood higher than asreml by 1.4e-5. asreml declares convergence and stays at its point when it is allowed to iterate further, so the two programs stop at different points of a flat direction. The variance of f(y) differs by 0.5 per cent, and the three components at the floor agree.

**lvr, truncated tent (1 check).** The likelihood is not unimodal: its support changes each time the range crosses an integer distance. remlax stops at a local optimum (range 1.51, log-likelihood -2.679) while a fine grid in R finds range 5.89 and -1.914. asreml stops at a worse point (-2.820). The remlax revision of 2026-09-28, which passed this check on the cluster then, gives the same local optimum on the authors' workstation: the outcome depends on the start, not on a change of the code.

**CPU/GPU parity, anisotropic Matern (1 check).** The fit of `mtrn_aniso` stops before its first iteration: two components of the gradient are not finite at the starting point. The anisotropic distance and the Matern function have an infinite derivative at zero distance, where the diagonal of the kernel sits. The check that refuses a non-finite start was added on 2026-09-28, after the parity of 2026-09-01 (34/34), which therefore did not see it. The 33 other designs agree to at most 6.8e-16 in relative value.

**Random sweep, seed 0 (1 check).** One configuration out of 200 ends with a Newton decrement of 0.011, above the threshold of 0.01. Two others reach a better optimum after a restart, with gains of 0.004 and 0.27 on -2logL. Seeds 50 and 100 have no unconverged configuration.

## What could not be measured here

**Structures with no external reference.** `mtrn` with all four parameters free, `own` (user-written kernels) and `sep` beyond three factors have no counterpart in asreml, lme4, sommer or nlme; they are checked against closed forms (besselK profile, the same kernel written as `exp`) and against the dense reference.

**What agreement means here.** Two REML implementations agree on the restricted log-likelihood at the optimum to 1e-6 absolute, on the variance components to 1e-4 relative (looser on purpose: on a flat surface two optimisers stop at slightly different theta for the same likelihood) and on fixed effects to 1e-6 relative. asreml stops on its own criterion and is compared at 1e-5 on the log-likelihood and 1e-3 on components; sommer is tightened to tolParConvLL = 1e-12 and compared at 1e-4. Iteration counts are never compared: average information, quasi-Newton and derivative-free optimisers do not count the same thing.

## Re-running all of it

```sh
# copies remlax and this repository at their commits to the cluster, submits both jobs
campaign/prepare.sh
# afterwards, from this repository
rsync -a cluster:<work folder>/logs/ validation/results/logs_<date>/
python3 validation/parse_logs.py validation/results/logs_<date>=<provenance> \
        --out validation/results/checks_<date>.csv
python3 validation/make_validation_md.py --csv validation/results/checks_<date>.csv \
        --versions validation/results/versions_<date>.json --rev <sha> --date <date>
```

`campaign/run_cpu.sbatch` lists every suite and the log it writes.

asreml requires a licence, and the `test_remlax_asreml*.R` suites check one out at run time. The other suites need no licence. The CPU/GPU parity check needs a machine with a GPU; everything else runs on a CPU.

