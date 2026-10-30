# Benchmark plan: remlax engines, CPU and GPU, against asreml, sommer and lme4

Written on 2026-10-22 and revised the same day. This benchmark grounds the `engine = "auto"` rule of `rx_fit`, the inversion rule of the MME engine and the choice of device. It also compares remlax with asreml, sommer and lme4 on time, on the logLik reached and on the estimated parameters, whose true values are known by simulation.

## Questions

- When is the MME engine faster than the dense engine, and when is the dense engine required? The dense engine is required outside the MME scope (residual none, rr, rx_lag).
- In the MME engine, from which fill-in of the factor does the sparse inversion (SuperLU and supernodal selected inversion) become slower than the dense inversion of C?
- For each engine, from which size does the GPU beat the CPU, compilation and transfers included?
- Does remlax reach the same optimum as asreml, sommer and lme4, and in how much time?

## Routes compared

| Route | Engine | Inversion | Device |
|---|---|---|---|
| D-CPU | dense (V formed) | Cholesky of V | CPU |
| D-GPU | dense | Cholesky of V | GPU |
| M-dense-CPU | MME | dense C (LAPACK) | CPU |
| M-dense-GPU | MME | dense C (JAX) | GPU |
| M-sparse-CPU | MME | sparse C (SuperLU and selected inversion) | CPU |
| asreml | asreml-R 4.2 | sparse AI-REML | CPU |
| sommer | sommer 4.4 (`mmer`) | direct formulation in V | CPU |
| lme4 | lme4 2.0 (`lmer`, pedigreemm for the pedigree) | sparse Cholesky | CPU |

The remlax optimiser is AI on the variance scale (`optimizer = "ai"`) for every route, with the same bounds, the same starting point and the same stopping criterion. The diagnostic Hessian is switched off (`hessian = FALSE`).

## Model families

Each family covers one regime of the ratio between n (observations) and m (size of C) and of the fill-in of C. An external program runs only if it can fit exactly the same model.

| Family | Model | External |
|---|---|---|
| iid | iid genotype, iid block | asreml, sommer, lme4 |
| grm | genotype with a dense GRM, iid block | asreml, sommer, lme4 (Z L, K = L L') |
| us | 3 traits, genetic us, residual us per unit | asreml, sommer |
| ped | animal model, sparse A^-1, iid block | asreml, sommer, lme4 (pedigreemm) |
| ige | DGE and IGE us(2) on neighbourhood incidences, iid IEE | asreml (n ≤ 3000) |
| field | iid genotype, ar1 x ar1 field as a random term | asreml |
| resid | iid genotype, residual ar1(row):ar1(col) | asreml |
| exp | exponential correlation on 1D coordinates | asreml |
| fa | 5 traits, fa of rank 2, diag residual | asreml |

sommer and lme4 estimate neither an ar1 or exponential correlation nor an fa structure. They are absent from these families.

## Sizes and replicates

n takes the values 500, 1000, 1500, 2000, 3000, 4000, 6000 and 8000. n counts the observations, all traits together. n grows by adding genotypes, and the number of observations per genotype stays roughly fixed. m, the size of C, therefore grows with n. An effect of n must not be read as an effect of the number of replicates per genotype.

| Family | Levels of the main term | Observations per level and per trait | m at n = 500 | m at n = 8000 |
|---|---|---|---|---|
| iid | genotypes = n/5 | 5 | 126 | 1615 |
| grm | genotypes = n/4 | 4 | 152 | 2027 |
| us | genotypes = n/12 (3 traits) | 4 | 126 | 2001 |
| ped | animals = n | 1 | 522 | 8022 |
| ige | genotypes = n/5 | 5 | 701 | 11201 |
| field | genotypes = n/4 | 4 | 624 | 9967 |
| resid | genotypes = n/4 | 4 | 125 | 1963 |
| exp | positions = n/3 | 3 | 167 | 2667 |
| fa | genotypes = n/20 (5 traits) | 4 (2.5 at n = 500, floor of 40 genotypes) | 205 | 2005 |

Genotypes are drawn with replacement. The design is therefore unbalanced, and the number of observations per genotype fluctuates around its mean. `banc_moteurs/R/dimensions.R` rebuilds each case from its seed and writes `dimensions.csv`. The analysis draws from it a figure of the dimensions against n, puts the scaling rule of each family in the title of its panels, and draws a figure of time against m.

Each cell (family, n) is simulated five times. The seed depends on the family, on n and on the replicate, and every route and program of a replicate reads the same data. The 360 cases (9 families, 8 sizes, 5 replicates) form one SLURM array.

## Machine

Everything runs on node io-gpu-a100-01, MIG slice 1g.10gb (1/7 of the A100), 6 cores and 32 GB per task, container `ige_reml.sif`. The CPU routes and the external programs use the cores of the slice, with OMP, OpenBLAS and MKL set to 6 threads. CPU and GPU therefore share the same machine. Each fit runs in its own R process (callr), with a time limit of 1800 s. The time is measured inside that process, around the fitting call alone.

## Recorded measurements

One row per fit, in `fits_<family>_<n>_<rep>.csv`: family, n, replicate, seed, software, route, status (OK, ERROR, TIMEOUT, MEMORY), total time, `compile_s`, `eval_s`, numbers of iterations and evaluations, effective engine and inversion, declared convergence, native logLik and arbiter logLik. The arbiter logLik is the shared dense REML, evaluated at the estimates of each program, because the constants of the native logLiks differ from one program to the next. A second table, `parameters_<family>_<n>_<rep>.csv`, gives the true and the estimated value of each parameter for each fit. Campaigns run before the English pass wrote `ajustements_*.csv` and `parametres_*.csv` with French columns. The analysis scripts still read them.

## Analysis

Within each replicate, the reference is the best arbiter logLik among the successful fits. A fit is accurate within 1e-3 of this reference. It is kept for times and parameters within 0.1. Time ratios are computed within the replicate, then summarised by the median, minimum and maximum of the five replicates.

Figures produced by `banc_moteurs/R/analyser.R`:

- dimensions of the cases (genotypes, observations per genotype, m) against n;
- status of the fits by family and route;
- time against n, every route, as total time then as time of the remlax solver alone;
- time against m;
- GPU gain (CPU / GPU) for the dense and MME engines;
- MME against dense, and sparse against dense inversion;
- remlax against asreml;
- fastest route per cell, every route then remlax only;
- difference to the best logLik;
- parameter recovery (bias, relative RMSE);
- agreement of the external estimates with remlax;
- share of JAX compilation in the total time.

The `auto` rule is built on the map of the fastest remlax route, from predictors that can be computed before fitting (n, m, nnz of the pattern of C, symbolic fill-in, presence of a structured residual or of a parameterised dense K).

## Launch

```
W=/scratch/users/salasn/banc_moteurs_20261022
sbatch --array=0-359%7 $W/remlax-validation/banc_moteurs/banc_mig.sbatch
Rscript banc_moteurs/R/analyser.R --res=<local copy of $W/res>
```

A task whose CSV already exists stops at once. A resubmitted array therefore redoes only the missing cases.

## Known limits

- A MIG slice underestimates the GPU at large sizes. The full A100 is not used.
- asreml `expv` as a random term does not reach the REML optimum of the simulated model (ρ pushed towards 0, arbiter logLik 20 units below that of remlax at n = 500, even on a case with one observation per position). The gap is reported, not corrected.
- asreml fits the ige model only at n ≤ 3000, because the neighbourhood incidences go through dense columns of the data frame.
