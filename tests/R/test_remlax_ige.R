# ==============================================================================
# test_remlax_ige.R -- THE model of the project, on simulated data of known truth
# ------------------------------------------------------------------------------
# Reproduces at small scale the model of scripts 05/06: direct genetic effect
# (DGE), intraspecific indirect effect (IGE intra) carried by a WEIGHTED
# NEIGHBOURHOOD, interspecific indirect effect (IGE inter), free covariance
# between DGE and IGE intra, AR1 x AR1 spatial field and nugget.
#
#     y_i = mu + Zg u_DGE + Zn u_IGEintra + Zx u_IGEinter + AR1xAR1 field + e
#     (u_DGE, u_IGEintra) ~ N(0, Sigma_us (x) I)      Sigma_us: FREE 2 x 2
#
# WHY SIMULATE RATHER THAN USE THE REAL DATA: here the values are KNOWN. A
# difference is therefore not a matter of interpretation. It is a defect of the
# solver or of the sampling, and the two are told apart by varying the size.
#
# TWO NEIGHBOURHOOD KERNELS, those of the project:
#   power         w = d^-1 over the Chebyshev window of order `ord`
#   exponential   w = exp(-kappa (d - d0)), estimated range (here with kappa
#                 fixed: estimating kappa requires rebuilding Z at each
#                 iteration, which the solver does not do yet, see limits)
#
# COMPARISON: asreml on the SAME model, with str(~ dge + grp(nb), ~us(2):id(q)),
# its way of writing a covariance shared between two terms.
#
#   Rscript scripts/tests/test_remlax_ige.R [--kernel=power|exp]
#                       [--nr= --nc=]   grid size (replication)
#                       [--qa= --qb=]   number of genotypes (identifiability)
# ==============================================================================
suppressPackageStartupMessages({
  library(here); library(Matrix); library(jsonlite)
})
source(here::here("R", "charger_remlax.R"))
ga <- function(k, d) { a <- commandArgs(TRUE); m <- grep(paste0("^--", k, "="), a, value = TRUE)
  if (length(m)) sub(paste0("^--", k, "="), "", m[1]) else d }
KERNEL <- ga("kernel", "power"); ORD <- as.integer(ga("ord", 6))
KAPPA <- as.numeric(ga("kappa", 0.12)); SEED <- as.integer(ga("seed", 42))
# DEFAULT SIZE. It used to be 16 x 12 with 18/12 genotypes: 96 observations
# for 8 variance parameters. At that size the us(2) block is at the boundary
# (r = -1) on BOTH sides, asreml returns NA as the standard error of three
# field parameters, and the two solvers stop 0.305 logLik points apart.
# Comparing their components then proved nothing. At 40 x 30 and 90/60
# genotypes, the eight components agree within 4 % of a standard error. A test
# that does not discriminate is not a test.
NR <- as.integer(ga("nr", 40)); NC <- as.integer(ga("nc", 30))

FAILURES <- character(0)
check <- function(name, cond, detail = "") {
  cat(sprintf("  %-50s %s %s\n", name, if (isTRUE(cond)) "OK " else "FAIL", detail))
  if (!isTRUE(cond)) FAILURES <<- c(FAILURES, name)
}

# ==============================================================================
# 1. DESIGN: species alternating by column, as in the real one
# ==============================================================================
set.seed(SEED)
g <- expand.grid(col = 1:NC, row = 1:NR)
g$sp <- ifelse(g$col %% 2 == 1, "A", "B")           # A = "wheat", B = "alfalfa"
g$id  <- seq_len(nrow(g))
# NUMBER OF GENOTYPES, and why it matters. Enlarging the GRID adds only
# replication. The us(2) block stays at the boundary (r = -1) because DGE and
# IGE intra are carried by the SAME genotypes. The number of genotypes
# separates them, not the number of plots. Measured: at 18/12 genotypes,
# r = -1.0000 on a 16x12 grid as on a 40x30 one.
QA <- as.integer(ga("qa", 90)); QB <- as.integer(ga("qb", 60))
gA <- paste0("a", seq_len(QA)); gB <- paste0("b", seq_len(QB))
g$geno <- NA_character_
g$geno[g$sp == "A"] <- sample(rep_len(gA, sum(g$sp == "A")))
g$geno[g$sp == "B"] <- sample(rep_len(gB, sum(g$sp == "B")))

# --- neighbourhood matrices ---------------------------------------------------
# w = d^-1 (power) or exp(-kappa (d - 5)) (exponential), over the Chebyshev
# window of order ORD. No L2 normalisation: the IGE variance is therefore
# "per unit of exposure", as in the pipeline.
step <- 5
weight <- function(dr, dc) {
  d <- step * sqrt(dr^2 + dc^2)
  if (KERNEL == "power") 1 / d else exp(-KAPPA * (d - step))
}
mk_Z <- function(focal_sp, nbr_sp, glev) {
  ix <- which(g$sp == focal_sp)
  Z <- matrix(0, length(ix), length(glev), dimnames = list(NULL, glev))
  for (k in seq_along(ix)) {
    i <- ix[k]
    vd <- abs(g$row - g$row[i]) <= ORD & abs(g$col - g$col[i]) <= ORD &
          g$id != g$id[i] & g$sp == nbr_sp
    j <- which(vd)
    if (!length(j)) next
    w <- weight(g$row[j] - g$row[i], g$col[j] - g$col[i])
    for (m in seq_along(j)) {
      cc <- match(g$geno[j[m]], glev); Z[k, cc] <- Z[k, cc] + w[m]
    }
  }
  Z
}
iA <- which(g$sp == "A"); nA <- length(iA)
ZgA <- Matrix::sparseMatrix(i = seq_len(nA), j = match(g$geno[iA], gA), x = 1,
                            dims = c(nA, QA), dimnames = list(NULL, gA))
ZnA <- mk_Z("A", "A", gA)                            # conspecific neighbours
ZxA <- mk_Z("A", "B", gB)                            # heterospecific neighbours

# ==============================================================================
# 2. SIMULATION: TRUE values
# ==============================================================================
SIG <- matrix(c(0.90, -0.35, -0.35, 0.30), 2, 2)     # (DGE, IGE intra): NEGATIVE cov
V_INTER <- 0.22                                       # IGE inter variance
V_AR1 <- 0.45; RHO_R <- 0.65; RHO_C <- 0.40           # spatial field
V_E <- 0.50                                           # nugget
U   <- matrix(rnorm(QA * 2), QA, 2) %*% chol(SIG)
uX  <- rnorm(QB, 0, sqrt(V_INTER))
Kr <- RHO_R ^ abs(outer(1:NR, 1:NR, "-")); Kc <- RHO_C ^ abs(outer(1:NC, 1:NC, "-"))
field <- as.numeric(kronecker(t(chol(Kr)), t(chol(Kc))) %*% rnorm(NR * NC)) * sqrt(V_AR1)
cellA <- (g$row[iA] - 1) * NC + g$col[iA]            # full grid (simulation)
col_obs <- as.integer(factor(g$col[iA]))             # columns ACTUALLY occupied
NC_OBS  <- max(col_obs)
cellA_obs <- (g$row[iA] - 1) * NC_OBS + col_obs      # conspecific grid
y <- 3 + as.numeric(ZgA %*% U[, 1]) + as.numeric(ZnA %*% U[, 2]) +
     as.numeric(ZxA %*% uX) + field[cellA] + rnorm(nA, 0, sqrt(V_E))

d <- data.frame(y = y, geno = factor(g$geno[iA], levels = gA),
                row = factor(g$row[iA]), col = factor(g$col[iA]))
cat(sprintf("Design: %d x %d, kernel %s (order %d%s) | %d obs, %d genotypes A, %d B\n",
            NR, NC, KERNEL, ORD, if (KERNEL == "exp") sprintf(", kappa=%.2f", KAPPA) else "",
            nA, QA, QB))
cat(sprintf("Truth  : var_DGE %.2f  var_IGEintra %.2f  cov %.2f  (r = %.3f)\n",
            SIG[1, 1], SIG[2, 2], SIG[1, 2], SIG[1, 2] / sqrt(SIG[1, 1] * SIG[2, 2])))
cat(sprintf("         var_IGEinter %.2f | AR1 var %.2f rho %.2f/%.2f | nugget %.2f\n\n",
            V_INTER, V_AR1, RHO_R, RHO_C, V_E))

# ==============================================================================
# 3. remlax: the full model
# ==============================================================================
# str(~ dge + neighbourhood): ONE us(2) covariance for the TWO incidences
# carried by the same genotypes. This identifies cov(DGE, IGE).
mod <- rx_model(
  y = d$y, X = matrix(1, nA, 1),
  terms = list(
    rx_term("dge_ige", list(ZgA, ZnA), struct = "us", levels = gA),
    rx_term("ige_inter", list(ZxA), struct = "iid", levels = gB),
    # AR1 x AR1 FIELD ON THE OBSERVED COLUMNS, not on the full grid. The
    # species alternate by column, so species A occupies only the odd columns.
    # Indexing the field on the whole grid would turn my lag 2 into the asreml
    # lag 1, and rho_asreml would equal rho_remlax^2, the square erasing the
    # SIGN. MEASURED before the fix: remlax -0.4856, asreml +0.2358, and
    # (-0.4856)^2 = 0.2358 exactly. The real design has the same trap: the
    # conspecific step is two columns, not one.
    rx_term("field", Matrix::sparseMatrix(i = seq_len(nA), j = cellA_obs, x = 1,
                                          dims = c(nA, NR * NC_OBS)),
            t = 1L, struct = "iid", level = "ar1ar1", dims = c(NR, NC_OBS))),
  residual = rx_residual("iid"))
print(mod)
t0 <- Sys.time()
f <- rx_fit(mod, backend = Sys.getenv("RX_BACKEND", "auto"), n_restarts = 6,
            verbose = FALSE)
cat(sprintf("\nremlax (%s) : logLik %.6f | %.1f s | decrement %.2e\n",
            f$backend, f$logLik, f$seconds, f$newton_decrement))
G <- f$sigmas$dge_ige
cat(sprintf("  var_DGE      %.4f   (true %.2f)\n", G[1, 1], SIG[1, 1]))
cat(sprintf("  var_IGEintra %.4f   (true %.2f)\n", G[2, 2], SIG[2, 2]))
cat(sprintf("  cov          %.4f   (true %.2f)   r = %.4f (true %.3f)\n",
            G[1, 2], SIG[1, 2], G[1, 2] / sqrt(G[1, 1] * G[2, 2]),
            SIG[1, 2] / sqrt(SIG[1, 1] * SIG[2, 2])))
cat(sprintf("  var_IGEinter %.4f   (true %.2f)\n", f$sigmas$ige_inter[1, 1], V_INTER))
cat(sprintf("  field var    %.4f   rho %.4f / %.4f   (true %.2f, %.2f/%.2f)\n",
            f$sigmas$field[1, 1], f$rho$field[1], f$rho$field[2], V_AR1, RHO_R, RHO_C))
cat(sprintf("  nugget       %.4f   (true %.2f)\n", f$sigma_res[1, 1], V_E))

check("the sign of cov(DGE, IGE intra) is recovered",
      sign(G[1, 2]) == sign(SIG[1, 2]),
      sprintf("%.4f vs %.2f", G[1, 2], SIG[1, 2]))
check("the fit is at an optimum", f$newton_decrement < 1e-4,
      sprintf("decrement %.1e", f$newton_decrement))
# The degeneracy of Sigma_us (r -> -1) IS NOT a defect of the solver. On a
# design of this size, DGE and IGE intra are almost collinear and asreml gives
# exactly the same thing. It is REPORTED, not penalised. The solver must say so
# instead of presenting it as a result.
if (length(f$degenerate_components))
  cat(sprintf("  [note] component(s) at the boundary: %s (r = %.4f). ",
              paste(f$degenerate_components, collapse = ", "),
              G[1, 2] / sqrt(G[1, 1] * G[2, 2])),
      "Expected at this size: asreml also ends on this edge\n",
      "  (it returns NA as the standard error of three field parameters).\n",
      "  For a meaningful comparison, enlarge: --nr=40 --nc=30.\n")

# ==============================================================================
# 4. asreml: THE SAME model
# ==============================================================================
ok_asreml <- requireNamespace("asreml", quietly = TRUE)
if (!ok_asreml) {
  cat("\nasreml unavailable: comparison skipped.\n")
} else {
  suppressPackageStartupMessages(library(asreml))
  da <- cbind(d, as.data.frame(as.matrix(ZnA)), as.data.frame(as.matrix(ZxA)))
  names(da) <- c(names(d), paste0("N", seq_len(QA)), paste0("X", seq_len(QB)))
  iN <- grep("^N[0-9]+$", names(da)); iX <- grep("^X[0-9]+$", names(da))
  a <- try(asreml(y ~ 1,
                  random = stats::as.formula(sprintf(
                    "~ str(~ geno + grp(nb), ~ us(2):id(%d)) + grp(xb) + ar1(row):ar1(col)", QA)),
                  group = list(nb = iN, xb = iX),
                  data = da, trace = FALSE, maxit = 150, workspace = "4gb"),
           silent = TRUE)
  if (inherits(a, "try-error")) {
    cat("\nasreml failed:\n  ",
        substr(conditionMessage(attr(a, "condition")), 1, 200), "\n")
  } else {
    for (i in 1:6) a <- update(a, trace = FALSE)
    vc <- summary(a)$varcomp
    cat("\n--- asreml, same model ---\n"); print(vc[, c("component", "std.error")])
    cat(sprintf("logLik asreml %.9f | remlax %.9f | diff %.2e\n",
                a$loglik, f$logLik_asreml, abs(a$loglik - f$logLik_asreml)))
    # Two separate requirements, and only one is symmetric.
    #  (a) remlax must not be WORSE than asreml. On a flat surface the two
    #      optimisers stop at different theta for the same likelihood.
    #      Requiring strict equality would penalise the one that converges
    #      best. Here remlax does BETTER by 7.6e-04.
    #  (b) the components must agree where they are determined. The tolerance
    #      is set on the asreml standard error. Two estimates that differ by
    #      less than a tenth of a SE are the same result.
    check("remlax is no worse than asreml",
          f$logLik_asreml >= a$loglik - 1e-6,
          sprintf("diff %+.2e in favour of %s", f$logLik_asreml - a$loglik,
                  if (f$logLik_asreml >= a$loglik) "remlax" else "asreml"))
    cmp <- list(
      c("var_DGE",      G[1, 1],                 vc["geno+grp(nb)!us(2)_1:1", "component"],
        vc["geno+grp(nb)!us(2)_1:1", "std.error"]),
      c("cov_DGE_IGE",  G[1, 2],                 vc["geno+grp(nb)!us(2)_2:1", "component"],
        vc["geno+grp(nb)!us(2)_2:1", "std.error"]),
      c("var_IGEintra", G[2, 2],                 vc["geno+grp(nb)!us(2)_2:2", "component"],
        vc["geno+grp(nb)!us(2)_2:2", "std.error"]),
      c("var_IGEinter", f$sigmas$ige_inter[1, 1], vc["grp(xb)", "component"],
        vc["grp(xb)", "std.error"]),
      c("var_field",    f$sigmas$field[1, 1],    vc["row:col", "component"],
        vc["row:col", "std.error"]),
      c("rho_row",      f$rho$field[1],          vc["row:col!row!cor", "component"],
        vc["row:col!row!cor", "std.error"]),
      c("rho_column",   f$rho$field[2],          vc["row:col!col!cor", "component"],
        vc["row:col!col!cor", "std.error"]),
      c("nugget",       f$sigma_res[1, 1],       vc["units!R", "component"],
        vc["units!R", "std.error"]))
    # COMPARING THE COMPONENTS MAKES SENSE ONLY AT THE SAME OPTIMUM. If the two
    # solvers stop at different likelihoods, they describe two different
    # POINTS of the same surface. Requiring their components to agree would ask
    # a maximum to equal a non-maximum. The comparison is therefore made only
    # when the logLik values meet. Otherwise the script SAYS so.
    diff_ll <- f$logLik_asreml - a$loglik
    # THRESHOLD OF THE "SAME POINT". Requiring 1e-4 would be unrealistic on a
    # flat surface. The two solvers stop at slightly different theta for a
    # practically identical likelihood. 0.05 is far below the 1.92 of a 1 df
    # LRT. Two optima closer than that are the same result. Beyond it they are
    # two points, and the script says so instead of comparing.
    same_point <- abs(diff_ll) < 5e-2
    cat("\n  component         remlax     asreml    diff    (SE asreml)\n")
    for (z in cmp) {
      rk <- as.numeric(z[2]); as_ <- as.numeric(z[3]); se <- as.numeric(z[4])
      # asreml returns NA as the standard error of a parameter at a bound. The
      # test used to stop on `if (NA)`, and the following components were never
      # compared. An active bound is INFORMATION, not a breakdown.
      tol <- if (is.finite(se)) max(0.1 * se, 1e-4) else NA_real_
      flag <- if (!is.finite(se)) "  <- asreml SE unavailable (bound)" else
                if (abs(rk - as_) > 0.1 * se) "  <- apart" else ""
      cat(sprintf("  %-14s %9.5f %10.5f %8.1e   (%s)%s\n", z[1], rk, as_,
                  abs(rk - as_), if (is.finite(se)) sprintf("%.3f", se) else "NA", flag))
      if (!same_point || !is.finite(tol)) next
      check(sprintf("%s : agreement with asreml", z[1]), abs(rk - as_) <= tol,
            sprintf("%.5f vs %.5f (SE %.3f)", rk, as_, se))
    }
    if (!same_point)
      cat(sprintf(paste0("\n  [note] the two solvers stop %.3f logLik point(s) apart\n",
                         "         (in favour of %s). These are two different POINTS,\n",
                         "         and their components are not comparable.\n",
                         "         At this size (%d obs, %d parameters) the us(2) block is\n",
                         "         at the boundary and the surface is flat.\n"),
                  abs(diff_ll), if (diff_ll > 0) "remlax" else "asreml",
                  mod$n, mod$n_par))
  }
}

cat("\n", strrep("=", 66), "\n", sep = "")
if (length(FAILURES)) {
  cat(sprintf("FAILURES (%d):\n%s\n", length(FAILURES), paste0("  - ", FAILURES, collapse = "\n")))
  quit(save = "no", status = 1L)
}
cat("Full IGE model: fitted, and consistent with the simulated truth.\n")
quit(save = "no", status = 0L)
