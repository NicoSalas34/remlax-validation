# remlax: a generic REML solver for CPU and GPU

**Project**: IGE Analysis 2024–2025 · **Written**: 31 August 2026 · **Extended and validated against asreml**: 1 September 2026 · **Structure catalogue and full IGE model**: 1 September 2026

remlax is a restricted likelihood solver **independent of the IGE model**. It is
driven from R with a syntax close to asreml and computed in JAX. It uses the GPU
when one is present, and the model does not change.

---

## 1. Scope

```
y = X β + Σ_k Z_k u_k + e,   u_k ~ N(0, Σ_k ⊗ K_k),   e ~ N(0, R)
```

A **term** is a set of `t` incidences (n × q) that share the same `q` levels,
with a covariance structure Σ on these `t` columns and a kinship K between
levels. This single object covers the following cases.

| case | t | K |
|---|---|---|
| simple random factor | 1 | I |
| genetic effect with a GRM | 1 | GRM |
| multi-trait `us` / `fa` | number of traits | I or GRM |
| IGE model of this project | DGE + within-species IGE + between-species IGE | GRM per species |

The last case forces the generality. Each of its columns has **its own**
incidence (a weighted neighbourhood), whereas a multi-trait model uses a trait
indicator.

## 2. Interface

Through a formula, in the spirit of asreml:

```r
source("R/remlax.R")

fit <- rx_reml(fixed    = y ~ 1 + traitement,
               random   = ~ vm(genotype, K = Kmat) + iid(bloc),
               residual = "units",
               data     = df,
               backend  = "auto")
```

Grammar of the random terms. The first argument is the grouping column.
`K =` and `rank =` are optional.

| syntax | Σ | parameters |
|---|---|---|
| `iid(f)` or `f` | σ² | 1 |
| `vm(f, K)` | σ², kinship K | 1 |
| `diag(f)` | diagonal between traits | t |
| `us(f)` | unstructured between traits | t(t+1)/2 |
| `fa(f, rank = k)` | ΛΛ' + diag(ψ) | nload + t |
| `mm(Z)` | **supplied** incidence | 1 |
| `ar1(f)` | σ² × AR1(ρ) between levels | 2 |
| `ar1(r, c)` | σ² × AR1(ρ_r) ⊗ AR1(ρ_c) | 3 |
| `str(~ a + b, struct = "us")` | **one** covariance for **several** terms | t(t+1)/2 |

### Catalogue of correlations between levels (appendix C of the ASReml-R 4.2 manual)

**Stationary 1D**: `cor` · `ar1` · `ar2` · `ar3` · `sar` · `ma1` · `ma2` ·
`arma` · `corb(order=)` · `corg`.
**Metric 1D**: `exp(coord=)` · `gau(coord=)` · `lvr`.
**Metric 2D**: `iexp` (city-block) · `igau` · `ieuc` (Euclidean) ·
`sph` · `cir` (with a range) · `aexp` · `agau` (anisotropic) · `mtrn` (Matérn).
**User-defined**: `own(f, expr=, n_par=)`.
**Separable**: `ar1(r, c)`.

The range models (`sph`, `cir`, `lvr`, `mtrn`) take a **distance** parameter.
It is therefore parameterised as `exp(θ)`, positive and unbounded, and not as a
correlation through `tanh`.

#### Two formulas missing from the manual

The help pages of `lvr` and `mtrn` give the signature but not the formula. Both
were **identified against asreml**, by comparing logLik curves at a fixed
parameter.

**`lvr`** is the truncated tent `C_ij = max(0, 1 − d_ij/φ)`. This is the linear
model with a sill from geostatistics. It is positive definite in 1D because it
is the autocorrelation of a box function. φ is a **range**, with constraint U in
asreml. Identification (q = 12, same data set, φ imposed on both sides):

| φ | 1 | 1,5 | 2 | 2,528 | 3 | 4 | 6 |
|---|---|---|---|---|---|---|---|
| asreml | −4,518854 | −3,301014 | −3,831719 | −3,419241 | −5,217616 | −4,338156 | −2,652383 |
| `max(0,1−d/φ)` | −4,518854 | −3,301014 | −3,831719 | −3,419241 | −5,217616 | −4,338156 | −2,652383 |

Six plausible families were ruled out before this one (`I + φ(max−d)`,
`min(i,j)`, spherical, `exp(−d/φ)`...). All of them missed by more than 0,6
point.

The curve is **not** unimodal. Its support changes each time φ crosses an
integer. On this data set asreml stops at φ = 2,53 (−3,419), while the global
optimum lies near φ = 47 (−2,723). This is a local optimum and not a formula
error. The regression test therefore compares against a **fine grid**. An
`optimize()` call would be trapped in the same way.

**`mtrn`** is the anisotropic Matérn of Haskard et al. (2007):

```
u =  dx·cos α + dy·sin α          h = ( |√δ·u|^λ + |v/√δ|^λ )^(1/λ)
v = −dx·sin α + dy·cos α          ρ = 2^(1−ν)/Γ(ν) · (h/φ)^ν · K_ν(h/φ)
```

δ acts **as a square root and in opposite directions on the two axes**. The
transformation preserves areas, so δ is not confounded with the range φ. Three
other conventions (δ on one axis only, δ without the square root) shift the
logLik by 0,5 to 4 points. Checked against asreml to 3×10⁻¹⁰ on (δ=2, α=0) and
(δ=2, α=0,6), and to 1,5×10⁻¹⁰ on ν ∈ {0,5 ; 0,8 ; 1,5} in the isotropic case.
λ = 1 (city-block) is checked exactly at ν = 0,5. **λ = 1 with ν ≠ 0,5 is not
positive definite** in two dimensions (eigenvalue −0,049 on a 6×5 grid). asreml
still returns a number. remlax refuses.

#### Any ν: K_ν in JAX

`jax.scipy.special` does not expose the modified Bessel function K_ν. A
callback to scipy would have no derivative, yet ν and φ are parameters to
estimate. `scripts/gpu/remlax/bessel.py` therefore computes directly

```
M(z,ν) = 2^(1−ν)/Γ(ν) · z^ν · K_ν(z) = ∫₀^∞ exp(ν log z + ν u − z cosh u) du  (symmetrised)
```

by **trapezoidal quadrature** on a fixed grid. The choice has a numerical
reason. The first version, a difference of hypergeometric series
(K_ν = π/2·(I₋ν − I_ν)/sin νπ), gave an absolute error of **3,5×10⁶ at ν = 3**.
That formula blows up at integer ν and already loses 8 digits at ν = 3,7. The
integral is analytic in ν and has no singularity in the order. Its trapezoidal
error decreases as exp(−π²/h) because the Poisson coefficients equal
K_{iξ}(z).

Checked against `scipy.special.kv` on ν ∈ [0,05 ; 8] (integers included) and
z ∈ [0 ; 80]: **maximum absolute error 1,5×10⁻¹⁴**, with derivatives in z and
in ν correct to 8×10⁻¹⁰.

#### `own`: a structure written by the user

asreml's `own(obj, fun)` calls an R function at each evaluation. This is
impossible here. R and JAX share no memory space, and above all the
**derivative** would be needed, which an opaque R function does not provide.
remlax therefore takes an **expression**. It crosses the boundary as a string
and is differentiated automatically:

```r
random = ~ own(col, expr = "exp(-lag*exp(p1))", n_par = 1)
```

Variables: `d` (Euclidean distance), `dx`, `dy`, `lag` (|i−j|), `I`, `J`,
`p1..pk` (**raw**, unconstrained parameters; the expression applies
`exp`/`tanh` itself if it needs positivity). `^` is translated to `**`. Without
this translation, an expression written in R produced
`xor does not accept dtype float64`. The result is normalised to a correlation
unless `normalise = FALSE`, since Σ already carries the scale.

**`ilv` is not provided.** asreml has one, but its formula could not be
recovered. None of the nine families tried reproduces its logLik curve. The
Euclidean tent, the natural candidate, is not positive definite in 2D (it
returns NaN). No structure is better than a wrong one under asreml's name.

### Structures of Σ

`iid`(=idv) · `diag`(=idh) · `us`(=corgh) · `fa(k)` · `rr(k)` · `chol(k)` ·
`ante(k)` · `corh`.

Parameter counts follow the manual. `chol(k)` and `ante(k)` count
(k+1)(ω−k/2). There is one documented exception. `rr(k)` counts
`n_loadings(ω,k)` instead of kω, because Γ is constrained to be trapezoidal.
Without this constraint Γ is defined only up to a rotation. The likelihood then
has a flat ridge of dimension k(k−1)/2 and the Hessian is singular by as much.
The fitted model is the same.

**Consistency check**: `chol(ω−1)`, `ante(ω−1)`, `fa(ω)` and `rr(ω)` are
**saturated** parameterisations of the same Σ. Checked: all of them reach the
`us` logLik to ≤ 3,4×10⁻¹³, with the same Σ to 7×10⁻¹⁴.

The **variance** lives in Σ. The structure between levels carries only a
**correlation**. This is the decomposition Σ_h = D C D of the manual. An asreml
`ar1v` is written `ar1(f)` here (Σ = σ²), and an `ar1h` is written
`ar1(f, struct = "diag")`.

Two choices avoid silent failures:

- **ar2 / ar3 are parameterised by partial correlations.** Imposing |φᵢ| < 1 is
  not enough. The stationarity region of an AR(2) is a triangle, and a φ outside
  it gives a "correlation" that is not positive definite. Levinson-Durbin maps
  (−1,1)^p exactly onto the admissible region. The reported φ remain those of
  asreml.
- **`cor` is bounded by θ > −1/(q−1)**, the exact positivity condition (at
  q = 10, θ = −0,3 gives an eigenvalue of −1,7).

One deliberate departure from the manual concerns `arma`. Its printed formula,
ρ₁ = (θ−φ)(1−θφ)/(1+θ²−2θφ), produces sequences that are **not** positive
definite for 204 draws out of 300 in the square |θ|<1, |φ|<1 (one of them gives
ρ₁ = −0,976 with ρ₂ = −0,944, which is impossible). remlax uses the standard
form of the ARMA(1,1). The fitted model is the same. Only the sign of the
reported MA parameter may differ.

### Residual: the same grammar as the random terms

```r
residual = ~ units                            # σ² I
residual = ~ us(trait):units                  # coupling between traits
residual = ~ ar1(row):ar1(col)                # separable residual field
residual = ~ diag(trait):ar1(row):ar1(col)    # both
```

Under the hood there is a single formula: `R[i,j] = Σ_trait[t_i,t_j] · C_unit[u_i,u_j]`.
The residual is a term like any other, with the identity as incidence. Nothing
had to be added to the solver. It only had to be exposed.

`str()` is the equivalent of asreml's `str(~ a + b, ~us(2):id(n))`. The grouped
terms share their levels and receive a common covariance. This correlates a
direct effect and a neighbour effect carried by the same genotypes, the case
that motivated the whole solver. Each term keeps **its own** incidence, `mm()`
included.

2D splines come from `rx_spl2d(x, y, nseg)`. This is a tensor product of
P-splines with a PS-ANOVA decomposition into three random blocks (`_x`, `_y`,
`_xy`) and a null part that **must** go into the fixed effects. The solver needs
no new machinery. A smooth surface is a random effect with a known incidence,
exactly as in sommer and asreml.

Multi-trait structures require `trait =` (data in **long format**, one row per
unit × trait). This format was preferred to a `cbind()` of responses because it
handles traits measured on different subsets of units without any special case.
Such subsets are the rule as soon as phenotyping is incomplete.

Some models cannot be written as a formula, for instance a **weighted**
incidence or heterogeneous columns. The explicit route uses the same machinery
without the syntactic sugar:

```r
mod <- rx_model(y, X,
                terms = list(rx_term("voisinage", list(W1, W2), K = Kb, struct = "us")),
                residual = rx_residual("us", trait = d$trait, unit = d$unit))
fit <- rx_fit(mod, backend = "auto")
```

## 3. CPU or GPU

`backend = "auto"` takes the GPU if `jax.devices()` sees one, otherwise the CPU,
and reports its choice. `"gpu"` **refuses explicitly** when no GPU is present
instead of falling back silently. The device is set once, through
`jax.default_device`. This is the only place in the code where the choice of
machine intervenes.

R and Python communicate through **files** (raw binary and a JSON manifest), not
through reticulate. On this cluster R and JAX live in two separate containers
(`ige_pipeline.sif`, `jax_gpu.sif`). Moreover, the file makes parity
verifiable, since both backends then read strictly the same input. `RX_PY`,
`IGE_JAX_CMD` or `IGE_JAX_SIF` say how to call Python.

## 4. Method, and what it shares with asreml

**Same objective**: the restricted likelihood. The optimum is therefore the
same, and the estimates coincide. This is checked below.

**Different path**: asreml uses average information (AI-REML) from the first
iteration. Here, L-BFGS-B on an **unconstrained** parameterisation makes the
approach, and a **regularised Newton polish** finishes the work. Σ is positive
by construction (the parameters are a Cholesky factor, never Σ itself), so no
constraint needs handling.

This mechanism proved decisive on two points:

- **Scaling of the objective.** `-2logL` is of order *n*, and so is its
  gradient. L-BFGS-B starts with an approximate Hessian equal to the identity.
  Its first step equals `-gradient`, hundreds of units in parameter space.
  Without division by *n*, measured on a one-factor model, the optimiser stops
  at `-2logL = 425,15` with |grad| = 74,6, against **414,3095** at the optimum.
  Eleven log-likelihood points are lost silently.
- **Polish.** From an L-BFGS-B run stopped early on purpose, the error on σ²
  falls from 2×10⁻² to 9×10⁻¹⁰ and the Newton decrement from 4×10⁻³ to
  9×10⁻¹⁸.

The reported stopping criterion is the **Newton decrement** `g' H⁻¹ g`. It is
the logLik rise still available, and is therefore comparable with the 3,84 of a
1 df LRT. `max|grad|` is not used, because it reaches several units at verified
optima.

## 5. Validation

### Inference: vpredict and Wald

```r
f <- rx_reml(y ~ trt, random = ~ gid, data = d,
             vpredict = c(h2 = "V1/(V1+V2)"), wald = TRUE)
```

The `Vi` are numbered in the order of `f$composantes_noms`. This order is
**fixed and documented**, because a silent renumbering would make every
expression wrong at the next run.

| | remlax | asreml |
|---|---|---|
| h² | 0,6685333 | 0,6685333 |
| SE(h²), delta method | 0,0500972 | 0,0500998 |
| Wald `trt`, 2 df | p = 8,776×10⁻⁹ | p = 8,776×10⁻⁹ |

**The remlax Wald test is CONDITIONAL** (type III: each term given all the
others). The asreml one is **sequential** by default (type I). Both coincide on
the last term of the model, the only case where the two notions meet. The test
therefore compares this term. Without `kenward_roger = TRUE` the denominator
degrees of freedom are not adjusted. The p-values are then those of a χ²/df,
valid asymptotically.

### Kenward-Roger, and predict

```r
f  <- rx_reml(y ~ trt + x, random = ~ iid(bloc), data = d, kenward_roger = TRUE)
pv <- rx_predict(f, classify = "trt", sed = TRUE)
```

Compared with `pbkrtest::KRmodcomp`, the reference implementation, on an
**unbalanced** block design (n = 27, 5 blocks of 3 to 7 plots):

| | remlax | pbkrtest |
|---|---|---|
| denominator df | 21,76897 | 21,76896 |
| F (`trt`, 2 df) | 2,790148 | 2,790148 |
| p | 0,0834186 | 0,0834186 |
| adjusted standard errors | to 1,3×10⁻⁷ | – |

Two points make the result wrong when missed. Both were found by comparison:

- **W is the inverse of the EXPECTED REML information** ½tr(P ∂V_i P ∂V_j), and
  not of the observed Hessian. The two coincide asymptotically but not at
  n = 27. The Hessian gave 21,83 df instead of 21,769.
- **A₁ and A₂ contract Θ against ∂Φ/∂θ = Φ P Φ**, and not against P. With P,
  they are of order n³ instead of O(1), and the degrees of freedom come out
  **negative** (−0,56 on this data set).

The **second-order term is omitted** on purpose. Measured, it is not invariant
under reparameterisation. On a block design, the same model gives
Φ_A = 0,4375 in variances against 0,3465 in log standard deviations (21 %
apart), while the first-order part coincides to 2×10⁻¹⁶. This term corrects the
second-order bias of the estimator *in the chosen parameterisation*. It makes
sense only in the natural components, where it vanishes as soon as V is linear
in them. pbkrtest and SAS make the same choice.

`predict` reproduces asreml: values to **1,7×10⁻¹⁰** and standard errors to
**1,2×10⁻⁸** on `classify = "trt"`, with a covariate (set to its mean) and an
averaged random block. With `include_random = TRUE`, the BLUP enters the
prediction and the error becomes a **prediction** error:

    L'ΦL − 2 sym(L'Φ X'V⁻¹ Z G M) + M'G M − M'G Z' P Z G M

All three terms matter. The second one is negative, because the correlation
between β̂ and û reduces the error of the total.

### CPU/GPU parity over the whole catalogue

`scripts/tests/export_bundles.R` serialises **34 designs** that cover every
structure, and `scripts/tests/parite_gpu.py` fits them on both backends. Both
read the *same* bundle, which is the point of the file format.

Result: identical −2logL, relative difference ≤ 1,1×10⁻¹⁴, most often
**exactly 0**. Replayed on an **A100 MIG slice** of the cluster (job 10326144,
image `ige_reml.sif`): the same 34 models, difference ≤ 4,1×10⁻¹⁵.

On a MIG slice (7 SMs out of 108), **XLA compilation dominates the computation
by far**. Anisotropic `mtrn` takes 252 s there against 2,8 s on CPU for a model
of 160 observations. The Bessel quadrature produces a kernel that is costly to
compile, not to run. For small models the CPU is the right choice. The GPU pays
off on the size of V, not on the number of models.

The sweep was useful mostly for the **three bugs that no suite had caught**.
A fourth one, found afterwards by rerunning the existing suites, deserves to be
known:

> **`$` does PARTIAL matching on R lists.** After a field `sigmas_res` was
> added to the result, the expression `r$sigmas`, not yet created at that point
> of `rx_read_result()`, started to match `sigmas_res`. The Σ of the terms were
> therefore appended to the list of residuals, and `fit$sigmas[[1]]` returned
> the **residual** instead of the first term. Seven failures across three suites
> came from this off-by-one shift. The logLik stayed correct, which made the
> diagnosis counter-intuitive. `rx_read_result()` now uses only `[[ ]]`, which
> matches exactly.


1. `own`: `^`, written in R, was the Python XOR
   (`xor does not accept dtype float64`). It is now translated to `**`.
2. `dsum(~ ar1(col) | site)` with several rows per column: R singular,
   −2logL = NaN **silently**. `validate()` checked duplicated units only under
   us/diag/fa. The check now covers every structure between units.
3. `ilv` returned NaN, which led to its removal.

**Fixed** parameters (asreml's `F` code): `fixed_theta = c(3, 5)` bounds these
parameters at their starting value. The diagnostic then sees them "at a bound"
and excludes them from the free subspace. This is the right behaviour for the
degrees of freedom of a later LRT. The `P` code (positive) is automatic, since
the parameters are log standard deviations, and `U` is the default.

### The full IGE model, on simulated data with a known truth

`scripts/tests/test_remlax_ige.R` reproduces the model of scripts 05/06 at a
small scale. It contains DGE and within-species IGE on a **weighted
neighbourhood** with a free `us(2)` covariance, between-species IGE, an
AR1×AR1 field and a nugget. There are eight variance parameters.

The design is 40×30 with 90/60 genotypes and species alternating by column. The
kernel is a power law d⁻¹ at order 6 or an exponential exp(−κ(d−δ₀)). asreml
fits the **same** model through
`str(~ geno + grp(nb), ~us(2):id(q)) + grp(xb) + ar1(row):ar1(col)`.

| component | remlax | asreml | difference | (asreml SE) |
|---|---|---|---|---|
| var_DGE | 1,09524 | 1,09592 | 6,7×10⁻⁴ | 0,189 |
| cov(DGE, within IGE) | −0,42211 | −0,41939 | 2,7×10⁻³ | 0,180 |
| var within IGE | 0,16268 | 0,17627 | 1,4×10⁻² | 0,354 |
| var between IGE | 0,54369 | 0,53748 | 6,2×10⁻³ | 0,341 |
| var field | 0,61009 | 0,61022 | 1,3×10⁻⁴ | 0,164 |
| ρ row | 0,53375 | 0,53215 | 1,6×10⁻³ | 0,131 |
| ρ column | 0,11931 | 0,11850 | 8,1×10⁻⁴ | 0,111 |
| nugget | 0,39824 | 0,39694 | 1,3×10⁻³ | 0,156 |

The eight components coincide to **less than 4 % of a standard error**, and
remlax stops 0,021 logLik point above asreml.

#### The size of the design decides what the test proves

The test used to run on 16×30 with **18/12 genotypes**, that is 96
observations for 8 variance parameters. At this size the `us(2)` block sits on
the boundary (r = −1,0000) **on both sides**. asreml returns `NA` as the
standard error of three field parameters, its way of saying that it sits on an
edge. The two solvers stop **0,305 logLik point** apart. Comparing their
components then proved nothing. They are two different points of the same flat
surface, and requiring them to coincide amounts to asking a maximum to equal a
non-maximum.

The **number of genotypes** separates DGE and within-species IGE, and the number
of plots does not, since both effects are carried by the same individuals.
Measured logLik difference between the two solvers:

| design | genotypes | n | logLik difference | components |
|---|---|---|---|---|
| 16×12 | 18 / 12 | 96 | 0,305 | not comparable |
| 40×30 | 18 / 12 | 600 | 0,015 | still on the boundary |
| 40×30 | 90 / 60 | 600 | **0,021** | **all 8 within < 4 % of an SE** |

The test therefore compares the components only when the two logLik values
agree to within 0,05, well below the 1,92 of a 1 df LRT. Otherwise it says so
explicitly instead of concluding to a disagreement.

**A trap uncovered by this test.** The AR1×AR1 field must be indexed on the
columns **actually occupied**, not on the full grid. Species alternate by
column, so asreml's lag 1 is lag 2 of the field on the full grid. Before the
fix: remlax ρ_col = −0,4856, asreml +0,2358, and (−0,4856)² = 0,2358
**exactly**, the square erasing the sign along the way. The same trap exists on
the real design.

### Against asreml, the reference of the field

`scripts/tests/test_remlax_asreml.R`: 20 checks, local enterprise licence.

| structure | agreement |
|---|---|
| `ar1(f)` | ρ, variance, residual and **logLik identical to 10⁻⁸** |
| separable `ar1(r, c)` | variance, ρ_row, ρ_column, residual, logLik to 10⁻⁸ |
| 2D spline, **same basis matrices** | variances and logLik to 3×10⁻⁸ |
| `str(~ a + b, us)` vs `str(~ gid + gid:w, ~us(2):id(40))` | the three elements of Σ and the **logLik to 1,9×10⁻¹²** |

**Log-likelihood convention**: asreml omits the constant
`(n−p)/2·log(2π)`, and remlax includes it (as lme4 does). The field
`logLik_asreml` makes the conversion, checked to 1,3×10⁻⁸. Without it the
difference is 219,6 on a data set with n = 240 and looks like a model
disagreement.

In one case **asreml refuses and remlax converges**: an AR1×AR1 field without
replication, where the field absorbs the whole residual. asreml stops on "1
singularities in the Average Information matrix". remlax returns an optimum
(decrement 1,8×10⁻²⁶) with the residual at 4×10⁻¹¹, and **flags** it through
`composantes_degenerees`. Both reach the same diagnosis, expressed differently.

### Against implementations written by others

`scripts/tests/test_remlax.R`: 31 checks.

| case | reference | result |
|---|---|---|
| one random factor (3 designs) | `lme4::lmer` | −2logL identical to **10⁻¹³**, BLUPs to 7×10⁻⁹ |
| two crossed factors | `lme4::lmer` | −2logL identical to 10⁻¹³ |
| genomic kinship | `sommer::mmer` | variances to 4×10⁻⁵ |
| 2 traits, genetic and residual `us` | `sommer::mmer` | G and R to the 4th decimal, genetic correlation **0,7460 vs 0,7460** |
| formula vs explicit route | itself | difference exactly 0 |

### Against algebra

`scripts/tests/test_remlax_core.py`: 22 checks. The analytical REML of a
balanced design is recovered to 10⁻⁹. The autodiff gradient matches finite
differences to 10⁻⁸. `V == Z (Σ ⊗ K) Z' + R`, computed by an independent route,
matches to **5×10⁻¹⁵**. `us` represents any PSD matrix to 10⁻¹⁵.

### Under a random sweep

`scripts/tests/stress_remlax.py` draws hundreds of configurations, including
awkward ones (zero variance, degenerate incidences, strong kinship, very
unbalanced designs). It checks **properties** instead of values: no exception,
finite logLik, negligible decrement, symmetric PSD Σ, and an optimum stable on
restart.

The sweep found three real defects, all fixed:

1. **15 crashes out of 120** in the spectral decomposition. The finite
   difference Hessian contains non-finite values when a step leaves the domain
   where V is positive definite. Fixed by switching to forward or backward
   differences, and by a diagnostic that declares itself unavailable instead of
   crashing.
2. **11 NaN logLik values**, all with a `us` residual. The generator produced
   duplicated (unit, trait) pairs, which makes R singular. The solver now
   **refuses** explicitly, since a user could make this error.
3. **`Singular matrix`** on `inv(V)`. Replaced by a solve with a fallback on the
   pseudo-inverse, with a warning that β and the BLUPs are then conditional.

After the fixes: **150 configurations, 0 exceptions, 0 non-finite logLik, 0 non
PSD Σ, 0 non-converged fits**. Two cases of **local optimum** remained, where a
restart from another point finds a better one. The Newton decrement cannot
detect them. It measures the rise available *locally*, so it is zero at the top
of a secondary hill. Hence the `n_restarts` option, which perturbs the retained
point and keeps the best one. It recovers 0,0085 and 0,017 of `-2logL` on the
two cases concerned.

## 6. Known limitations

- **V is formed explicitly** (dense n × n). At n = 16 000 this takes 2,1 GB in
  float64. This is the target regime, and nothing beyond it. Going larger would
  require a sparse factorisation or Henderson's equations. `dense_Z()` and
  `assemble_V()` are the only two places to change.
- **Restarts are off by default** (`n_restarts = 0`), because they multiply the
  cost. Turn them on for any fit meant for publication.
- ~~Checked on CPU only~~: **lifted on 01/09/2026**. JAX CUDA is installed on
  the workstation (RTX A1000 6 GB), and the suites run on GPU. Parity measured
  on the same model: **−2logL identical to the bit**, Σ to 2,8×10⁻¹⁵, BLUPs to
  3,5×10⁻¹⁴, β to 4,2×10⁻¹⁶.

  A trap was fixed along the way. `jax.devices()` **without an argument**
  returns only the default platform. As soon as a GPU is visible it lists no
  CPU, so filtering its output on `platform == "cpu"` returned an empty list and
  `--backend=cpu` fell back silently on the GPU. The correct call is
  `jax.devices("cpu")`.
- **Still out of scope**: non-Gaussian models (GLM). They would require an
  IRLS/Laplace layer above the solver, a project in its own right. Also `ilv`,
  whose asreml formula could not be recovered (see §2).
- **Kenward-Roger without the second-order term.** The full formula carries
  −R_ij/4 with R_ij = U′(∂²V/∂σᵢ∂σⱼ)U. This term is zero as soon as V is linear
  in the variance parameters. This holds for every structure here **except**
  the correlation parameters between levels (φ of an AR1, a range). pbkrtest and
  SAS make the same choice. The field `second_ordre_omis` says whether the model
  contains parameters for which the omission is not exact.
- **`predict()` covers the fixed effects and a single-trait random term.** A
  classify on a multi-trait term is flagged and ignored in the random part.
- **The range of the neighbourhood kernel is not estimated.** The chain
  κ → Z → V would require rebuilding the incidence at each iteration.
  `scripts/gpu/reml_uni_kernel.py` already does this, and the code would have to
  be ported here.
