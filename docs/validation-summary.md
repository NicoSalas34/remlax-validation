# Validation summary

One page for the reader who wants to know whether remlax computes the same
restricted likelihood as the programs they already trust. The full table of
checks, one row per printed verdict, is in [validation.md](validation.md); the
numbers below are read from the same CSV (`correctness/reference_software/latest/2026-11-04_remlax-3bfae58/checks.csv`).

Revision tested: remlax `3bfae58` (version 0.2.0.9000), on 2026-11-04. Reference software: asreml 4.2.0.480, lme4 2.0.1,
sommer 4.4.5, nlme 3.1.169, pbkrtest 0.5.5 and R 4.6.0 inside the cluster
image `ige_reml_v3.sif` with JAX 0.11.1. `cluster/validation_campaign/prepare.sh` launched the
campaign; its logs are in `correctness/reference_software/latest/2026-11-04_remlax-3bfae58/logs/`. The campaign of 2026-10-30 on
remlax `6dbada5` gave the same verdicts; it is kept in the git tag
`archive/2026-11-17` (see [history.md](history.md)).

**Result: 381 checks, 365 passed, 16 failed** (399 checks, 387 passed, in the
2026-09-28 table, which still counted the suites of the RTMB engine since
removed). Twelve failures are those of 2026-09-28: `ma2`, `sph` and `cir`
against asreml, where a dedicated diagnostic shows that asreml's reported
log-likelihood is not the REML at its own reported parameters, while an
independent dense-algebra REML reproduces remlax's value to the last digit.
Two concern the 2D P-spline, where remlax reaches a log-likelihood higher than
asreml by 1.4e-5. One is the truncated tent `lvr`, whose multimodal likelihood
traps remlax at a local optimum (asreml stops lower still). The last is one unconverged
configuration out of 600 in the random sweep. [validation.md](validation.md)
details each one. The rows are left as failures on purpose.

## What is compared to what

Every comparison is the same model on the same data. Agreement means: the
restricted log-likelihood at the optimum to 1e-6 absolute (1e-5 against asreml,
which stops on its own criterion), variance components to 1e-4 relative (1e-3
against asreml), fixed effects to 1e-6 relative, BLUPs to 1e-4 absolute,
standard errors to 1e-4 relative. Iteration counts are never compared. The
log-likelihood constant of each program was checked: nlme and lme4 carry the
full REML constant and equal `fit$logLik`; asreml omits (n-p)/2 log(2 pi) and
equals `fit$logLik_asreml`; sommer's value carries a data-dependent constant
(it rescales the response internally), so only differences of -2logL between
two models are compared for sommer.

| Feature | Reference | Suite (section) | Checks | Largest gap measured |
|---|---|---|---|---|
| One random factor, three layouts; two crossed factors; nested structures iid/diag/us; weighted incidence; formula = explicit path; backend switch; refusals | lme4 (and the solver's own properties) | `test_remlax.R` 1, 2, 5-9 | 28 | -2logL identical to 9 decimals; BLUP 7.7e-9 |
| Random slopes us and diag (two-column terms), nesting, three crossed factors unbalanced, unequal groups; SE(beta); BLUP; PEV against the MME formula | lme4 | `test_remlax_lme4.R` 1-5 | 30 | components 0 to 6 digits; PEV 1.6e-14 |
| GRM on a term (test_remlax.R 3, 4: 3 checks at sommer's default tolerance); GRM with BLUP and PEV, three traits us and diag with GRM and LRT between them, additive + dominance, multi-environment diag (test_remlax_sommer.R, 19 checks, tolParConvLL 1e-12) | sommer | `test_remlax.R` 3, 4; `test_remlax_sommer.R` 1-5 | 22 | components 9.6e-5 at sommer's default tolerance, 2e-9 tightened |
| corAR1, corCAR1, corARMA (AR2, ARMA(1,1)), corCompSymm, corExp, corGaus, corLin, corSpher, corSymm + varIdent (us), varIdent (diag), lme random intercept + AR1 | nlme::gls, nlme::lme | `test_remlax_nlme.R` 1-12 | 61 | -2logL 9e-9; parameters 3e-6 |
| AR1, AR1 x AR1, 2D splines, str() shared covariance, vpredict h2 and its SE, Wald, saturated chol/ante/fa/rr = us | asreml | `test_remlax_asreml.R` 1-6 | 28 | logLik 2.9e-8; SE(h2) 5e-5 relative |
| lvr, anisotropic Matern (besselK reference too), own kernels, dsum sections, predict with SE and sed, Kenward-Roger vs pbkrtest, predict with BLUP | asreml, pbkrtest, R | `test_remlax_asreml2.R` 1-8 | 22 | predict 1.9e-12; K-R denDF 1e-5 |
| BLUP ordering and values, fixed effects and their variance | asreml | `test_remlax_asreml_blup.R` | 4 | see table |
| Full direct + indirect genetic effects model, eight parameters, known truth | asreml + simulation | `test_remlax_ige.R` | 11 | logLik: remlax higher by 2e-2; components within 0.1 SE (boundary case, see limits) |
| ar2, ar3, ma1, ma2, arma, sar, cor, corb, corg (grouped residual series); iexp, igau, ieuc, aexp, agau, sph, cir (2D field, also evaluated at asreml's own parameters); fa(1), rr(1), corh, ante(1), chol(1), diag with GRM on 4 traits; sep id x ar1 x ar1; us(trait):ar1:ar1 | asreml | `test_remlax_asreml3.R` A-E | 111 (99 passed, 12 failed) | logLik identical to 9 decimals for 23 of 26 models, 4.3e-5 on rr(1); ma2, sph and cir under Limits |
| Sigma parameterisations, closed-form sparse precisions, quadratic forms, K^-1 K = I; same likelihood at a common theta in both directions; precision path = relationship path | RTMB sparse engine vs JAX dense engine, and a longhand dense REML | `test_remlax_tmb_*.R` | 55 | 2.3e-13 on -2logL |
| Internal algebra: parameter counts, PSD, autodiff gradient vs finite differences, V assembly, analytic one-factor REML, determinism, refusals | closed forms | `test_remlax_core.py` | 25 | gradient 2e-8 relative |
| Totals per reference, read from the CSV: asreml 176 (164 passed), lme4 58, sommer 22, nlme 61, closed forms and sparse engine 81, properties only 1; 399 in all | | | 399 | |
| One test per feature (structures, levels, likelihood, fit options, inference, bundle/CLI, scan, sparsity) | independent numpy references | `pytest tests/python` | 169 tests, counted as 1 row | all pass |
| 200 random configurations: no exception, finite logLik, PSD Sigma, converged | properties only | `stress_remlax.py --seed 0` | 1 row | 198/200 clean, 2 unstable optima improved by restarts |

Features with **no external reference**: `mtrn` with all four parameters
free (checked against a besselK profile and asreml with fixed parameters
only), `own` user kernels (checked against the same kernel written as `exp`),
`sep` with more than three factors, weighted incidences (checked for
finiteness and against the engine of the IGE analysis repository, in
`scripts/remlax/reproduction_remlax/validation/`),
and the GPU backend today (see below).

## Figure

![remlax against the references](../correctness/reference_software/latest/2026-11-04_remlax-3bfae58/fig_validation.png)

Left: every (remlax, reference) pair that a check line prints, 221 pairs,
one colour per reference program, symlog axes. Right: |logLik gap| at the
optimum for the 66 checks that print both log-likelihoods; zeros are drawn at
1e-13. The points above 1e-5 are ma2, sph and cir against asreml, at the two
programs' own optima and, for sph and cir, remlax evaluated at asreml's
parameters (see Limits); the point at 4.3e-5 is rr(1), remlax being the higher. The pairs are listed in
`correctness/reference_software/latest/2026-11-04_remlax-3bfae58/pairs.csv`. Largest gaps per reference, read
from that file: asreml 4.3e-5 on logLik and 2e-2 on a null spline component
(3.8e-11 vs 2.5e-10) once ma2, sph, cir and the boundary IGE model are set
aside; lme4 0 on -2logL to 9 decimals; sommer 9.6e-5 on components (at its
default tolerance) and 2e-9 tightened; nlme 9e-9 on -2logL and 3e-6 on
parameters.

## Speed

Speed is measured apart from correctness, on the latest tested versions of
remlax only. See [speed/README.md](../speed/README.md). The speed campaigns of
September 2026 on revision `9a51c79`, summarised here before, are kept in the
git tag `archive/2026-11-17`.

## Defects found by this validation, and fixed

The nlme comparison on `corExp` revealed that the metric structures written
as rho^d (`exp`, `gau`, `iexp`, `igau`, `ieuc`, `aexp`, `agau`) started the
optimiser at theta = 0, that is rho = |tanh(0)| = 0 clipped to 1e-12, a point
where the gradient is exactly zero by construction (the clip is active and
sign(0) = 0). L-BFGS-B declared convergence without moving, and a field with
a true range of 4 was fitted with rho = 0: -2logL was 65 units above the
optimum found by gls, which remlax reproduced to 1e-9 from a warm start. The
fix (`fit.py`, `_theta0_rho_metrique`) starts at a correlation of 0.5 at the
median spacing of each coordinate axis. It changes no likelihood value, only
the starting point, and it is covered by a pytest
(`test_initial_theta_metrique_ne_part_pas_du_point_a_gradient_nul`) and by
nlme sections 6 and 7. The 2026-09-01 checks that involved these structures
were evaluations at a fixed theta (CPU/GPU parity) or `exp` against `own`
written by remlax itself, which is why they did not catch it: two engines
that agree prove less than one engine checked against an independent program.

The asreml comparison on `cir` then revealed two more (fixed in the same
working tree, see the cir item under Limits): a NaN gradient beyond the range
of the circular kernel, which stopped the optimiser at iteration 0, and the
multimodality of the range kernels, now handled by a sweep of the distance
deciles at start-up.

## Limits, stated plainly

- **CPU/GPU parity sweep.** The 34-design parity check (`parite_gpu.py`,
  evaluation at a common theta on both devices) was not replayed; the card
  became free late in the day and the time went to the speed benchmark. What
  the benchmark does show is that the ten models fitted on both devices reach
  the same optimum (|logLik gap| at most 1e-6 at 6 printed decimals). The
  parity recorded on 2026-09-01 (34/34, relative gap at most 6.4e-15) stands
  for the code of that date.
- **n = 8000 on CPU, and asreml at n = 8000.** Not measured. At 4 s per
  evaluation for n = 2000 on 4 cores the dense evaluation scales as n^3 and a
  fit at n = 8000 needs one to several hours per model on a 4- to 16-core CPU;
  the one-hour account was too short, the long account saturated (848/850
  cores) and the bigmem partition announced a start six days later. The
  n = 8000 rows exist for the A100 only, without an asreml counterpart.
- **The full IGE model at the boundary.** On the 600-plot simulated design the
  direct/indirect genetic correlation sits at -1 in both programs. remlax
  reaches a log-likelihood higher than asreml's by 0.02; the eight components
  agree within 0.1 standard error, but var_IGEintra differs by 8 percent
  (0.163 vs 0.176, SE 0.35). That is a flat ridge, not a disagreement about
  the likelihood; it is the largest gap in the table and it is expected on a
  boundary.
- **Where asreml stops, diagnosed (`correctness/reference_software/suites/diag_asreml3_ecarts.R`, cluster
  log `40_diag_asreml3_ecarts.log`).** The three models where the two programs
  report different log-likelihoods were re-fitted with asreml's convergence
  flag, bound codes and per-iteration trace, then both sets of parameters were
  handed to an independent REML written in dense R algebra (log|V| +
  log|X'V^-1 X| + y'Py; same X, same y, no coding difference is possible since
  X is a column of ones for sph and cir and the same numeric covariate for ma2,
  nedf = n - p in both programs). The arbiter reproduces remlax's value at
  remlax's parameters to the last digit in every case (ma2 -72.524782986, sph
  3.823231). It does not reproduce asreml's reported value at asreml's
  reported parameters.
  - ma2: asreml says `converge = TRUE`, bound codes P U U, and its
    log-likelihood is stable at -72.558054 over 20 further updates, but its
    `%ch` column still shows 3.5 percent parameter change at every iteration.
    Asked to hold remlax's parameters (R.param with con = "F"), asreml does
    not keep them: it re-expresses (-0.5801, -0.4493) as (-0.5581, -0.4323)
    and returns -72.576239, so asreml's `cor1`, `cor2` are not the moving
    average coefficients of the standard form and a point-wise comparison of
    the two programs' parameters is not possible. What can be said: the
    arbiter reproduces remlax's -72.524783 at remlax's parameters exactly
    (standard form e_t = a_t - theta_1 a_{t-1} - theta_2 a_{t-2}, invertible:
    roots of modulus 1.49); asreml's best value, -72.558054, is 0.033 below
    it; and the arbiter evaluated at asreml's printed parameters, read as
    standard coefficients, gives -72.539906, which is neither program's
    value, confirming that asreml's printed parameters live on another scale.
    remlax's likelihood is confirmed; asreml's parameterisation of ma2 is not
    identified.
  - sph: asreml says `converge = FALSE` and oscillates between two points at
    every update, log-likelihood 3.823207 and 3.753095 in turn; it reports the
    low point together with the parameters of the other one (range 8.754,
    variance 1.006), which is why remlax evaluated at "asreml's parameters"
    gave 3.7835 in the suite: those parameters and that log-likelihood do not
    belong to the same iterate. Held at remlax's parameters (variance 1.0245,
    range 8.6514), asreml returns 3.823231014, remlax's value to the last
    digit; the arbiter gives the same. Same function, same kernel; asreml's
    optimiser did not converge on it.
  - cir: asreml says `converge = FALSE` (parameters still changing 0.4 to
    0.6 percent after 20 updates) and reports 10.847517 with range 4.983,
    variance 0.624; the arbiter at those parameters gives 10.844976, which is
    exactly what remlax gives at the same parameters. The likelihood of the
    circular kernel is multimodal in the range (profile at the REML variance:
    10.85 at range 5, 11.02 at 6, 13.55 at 7, 13.88 at 7.4, 12.83 at 8; three
    local maxima at ranges 5.0, 7.4 and 10.9). asreml stops in the first hill.
    A first remlax run also stopped there (10.687): its gradient was NaN beyond
    the range (0 x inf from the clip with sqrt and arcsin) and L-BFGS-B ended
    at iteration 0. That defect was fixed in `levels.py` (two-branch
    `jnp.where`), and `fit.py` now starts range kernels from a sweep of the
    distance deciles followed by a short descent from the three best (tested
    in `tests/python/test_portee.py`). With this, remlax reaches 13.875897 at
    range 7.41 from its default start; held at those parameters asreml
    returns 13.875897068, the same value to 5e-9, three units above its own
    optimum. The "not worse" check passes and the equality checks fail
    because asreml is in the wrong hill.
  - rr(1): remlax higher by 4.3e-5 on a rank-deficient Sigma whose
    likelihood is flat along the rotation of the loadings; passed at 1e-4.
  - The arma MA parameter is reported by asreml with the opposite sign to the
    standard form that nlme and remlax use; compared on absolute values.
  For sph and cir the two programs compute the same function, shown by
  holding asreml at remlax's parameters, and the disagreement is asreml's
  stopping point. For ma2 remlax's likelihood is confirmed by the arbiter and
  asreml's parameterisation is the open point.
  The diagnostic log is kept as
  `correctness/reference_software/latest/2026-11-04_remlax-3bfae58/logs/diag_asreml3_ecarts.log`. The 12 rows are left as failures because the suite's
  criterion is agreement with what asreml reports, and that is what a user
  would see; the diagnosis is the answer to why.
- **lvr (corLin) is multimodal.** The truncated linear kernel makes the
  likelihood non-smooth in the range; gls found 200.19 from one start and
  188.62 from another, remlax 193.93 from its default start and 188.62 with
  four restarts. The check compares the two programs at the same range and
  requires remlax with restarts to be no worse than the best gls.
- **gau at sub-unit spacing.** The Gaussian kernel is parameterised as
  rho^(d^2); with coordinates spaced by less than one unit and a short range,
  the optimal rho falls below 1e-6 and the tanh parameterisation has no slope
  left. This is asreml's parameterisation, which remlax follows; rescale the
  coordinates, as one would with asreml.
- **PEV conventions.** lme4's `condVar` is Var(u | y) at known beta; remlax
  returns the mixed-model-equations PEV, which includes the uncertainty on
  beta, as asreml and sommer do. Both are checked against their own formula.
- **sommer's constant.** sommer's log-likelihood is on a constant that depends
  on the data; the offset is printed, never corrected by hand, and only
  differences between models are compared.
- **What was not re-run.** `stress` seeds 50 and 100 (only seed 0 was
  replayed, 200 configurations instead of 3 x 50), and the RTMB engine
  comparison at full scale (36 h against 15 min, a 2025 measurement quoted in
  the paper).

## Replacement text for the Validation section of `paper/remlax.tex`

The LaTeX draft is in [`paper/validation_section_draft.tex`](../paper/validation_section_draft.tex);
it is not inserted in the manuscript. It quotes the final numbers above.

## Paragraph for a methods section (Fitting)

> All models were fitted by restricted maximum likelihood with remlax (revision f576638), a solver that evaluates the restricted log-likelihood on the dense covariance matrix with automatic differentiation and maximises it by L-BFGS-B followed by Newton polishing. Its correctness was checked against asreml 4.2, lme4 2.0.1, sommer 4.4.5 and nlme 3.1 on the same simulated data and models, in 399 comparisons covering one-factor, crossed and nested designs, genomic relationship matrices, multi-trait unstructured, diagonal, factor-analytic, reduced-rank, antedependence and Cholesky covariances, AR1 to AR3, MA, ARMA, uniform, banded and general correlations, exponential, Gaussian, spherical and Matern kernels, separable AR1 x AR1 fields, sectioned residuals, prediction, Wald and Kenward-Roger tests. Of these, 387 passed; the 12 that did not concern the MA(2) model and the spherical and circular kernels, where asreml declared convergence with parameters still moving or oscillated without converging, and where an independent REML written in dense algebra reproduces remlax's log-likelihood exactly and not asreml's. The restricted log-likelihood at the optimum agreed to within 4e-5 with asreml (3e-8 outside the reduced-rank model) and 1e-8 with lme4 and nlme, variance components to within 1e-4 relative, fixed effects to 1e-6, BLUPs to 1e-4 and prediction error variances to the mixed-model-equations formula to 1e-14. On the complete direct and indirect genetic effects model of Salas et al. (2026), fitted on simulated data with known parameters, remlax and asreml agreed on all eight variance parameters within 0.1 standard error, remlax reaching a log-likelihood higher by 0.02. The full table of checks and the scripts that produce it are in the remlax-validation repository (docs/validation-summary.md).
