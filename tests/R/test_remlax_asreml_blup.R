# =============================================================================
# BLUP, prediction errors and TIME: remlax against asreml
# =============================================================================
# WHY THIS FILE EXISTS. All the other validation is about the likelihood and
# the variance components. A breeder uses the BLUPs and their prediction error
# variance. A defect in this path would show in NO current test, exactly as the
# missing log|K| did not show in the estimates, since a constant does not touch
# the gradient. This is the most serious correctness gap of the repository.
#
# TIME IS MEASURED HERE AND NOWHERE ELSE. asreml exists only in the container,
# and this is the only speed comparison that is both available and expected by
# a reader. Both engines fit the SAME design, in the same process, one after
# the other.
#
# WRITTEN DEFENSIVELY. The exact names of the asreml outputs (coefficients,
# standard errors of the random effects) cannot be checked outside the
# container. This script PRINTS what it finds before comparing, and an
# unexpected structure produces a diagnostic, not a silent failure.
suppressPackageStartupMessages({library(Matrix); library(asreml)})
source("R/charger_remlax.R")
ok <- 0L; ko <- 0L
check <- function(name, cond, detail = "") {
  if (isTRUE(cond)) { ok <<- ok + 1L; cat(sprintf("  OK    %-44s %s\n", name, detail)) }
  else { ko <<- ko + 1L; cat(sprintf("  FAIL  %-44s %s\n", name, detail)) }
}

set.seed(20260901)
nr <- 30L; nc <- 24L; n <- nr * nc; q <- 60L
geno <- factor(rep_len(seq_len(q), n), levels = seq_len(q))
u <- rnorm(q, 0, 1.1)
row <- factor(rep(seq_len(nr), each = nc)); col <- factor(rep(seq_len(nc), nr))
y <- 5 + u[as.integer(geno)] + rnorm(n, 0, 0.8)
d <- data.frame(y = y, geno = geno, row = row, col = col)

cat("=== 1. fit on both sides, time measured ===\n")
t_a <- system.time(a <- asreml(y ~ 1, random = ~ geno, data = d, trace = FALSE))[["elapsed"]]
t_r <- system.time(f <- rx_reml(y ~ 1, random = ~ geno, data = d,
                                hessian = TRUE, blups = TRUE, verbose = FALSE))[["elapsed"]]
cat(sprintf("  asreml %7.2f s | remlax %7.2f s | ratio %.2f\n", t_a, t_r, t_r / t_a))
check("logLik in the asreml convention", abs(f$logLik_asreml - a$loglik) < 1e-5,
      sprintf("%.9f vs %.9f", f$logLik_asreml, a$loglik))

cat("\n=== 2. what each side exposes (diagnostic before comparison) ===\n")
cat("  remlax$blups :", paste(names(f$blups), collapse = ", "), "\n")
sa <- try(summary(a, coef = TRUE), silent = TRUE)
cat("  asreml summary(coef=TRUE) :",
    if (inherits(sa, "try-error")) "unavailable" else paste(names(sa), collapse = ", "), "\n")
cr <- if (!inherits(sa, "try-error") && !is.null(sa$coef.random)) sa$coef.random else NULL
if (is.null(cr)) cr <- try(a$coefficients$random, silent = TRUE)
cat("  columns of the asreml random effects:",
    if (inherits(cr, "try-error") || is.null(cr)) "not found"
    else paste(colnames(as.matrix(cr)), collapse = " | "), "\n")

cat("\n=== 3. BLUP: do both sides predict the same effects? ===\n")
bl_r <- as.numeric(f$blups[["geno"]])
cr_m <- if (!is.null(cr) && !inherits(cr, "try-error")) as.matrix(cr) else NULL
if (!is.null(cr_m)) {
  rows_g <- grep("^geno", rownames(cr_m))
  bl_a <- as.numeric(cr_m[rows_g, 1])
  # asreml orders its levels as the factor does; align by NAME, never by position
  names_a <- sub("^geno_?", "", rownames(cr_m)[rows_g])
  ord <- match(levels(geno), names_a)
  bl_a <- bl_a[ord]
  n_comm <- sum(!is.na(bl_a))
  cat(sprintf("  %d levels matched out of %d\n", n_comm, q))
  if (n_comm >= q %/% 2) {
    max_diff <- max(abs(bl_r - bl_a), na.rm = TRUE)
    corr <- cor(bl_r, bl_a, use = "complete.obs")
    check("BLUP : maximum difference", max_diff < 1e-4, sprintf("max|diff| = %.3e", max_diff))
    check("BLUP : correlation", corr > 1 - 1e-8, sprintf("r = %.12f", corr))
  } else {
    check("BLUP : level matching", FALSE,
          sprintf("only %d levels matched", n_comm))
  }
  if (ncol(cr_m) >= 2) {
    se_a <- as.numeric(cr_m[rows_g, 2])[ord]
    se_r <- if (!is.null(f$blup_se) && !is.null(f$blup_se[["geno"]]))
              as.numeric(f$blup_se[["geno"]]) else NULL
    if (!is.null(se_r)) {
      e <- max(abs(se_r - se_a), na.rm = TRUE)
      check("prediction error of the BLUPs", e < 1e-4, sprintf("max|diff| = %.3e", e))
    } else {
      cat("  NOTE: remlax exposes no BLUP standard error ($blup_se absent).\n")
      cat("         This is a GAP to fill, not a disagreement. Fields present:",
          paste(grep("blup|se", names(f), value = TRUE), collapse = ", "), "\n")
    }
  }
} else {
  check("retrieval of the asreml BLUPs", FALSE, "unexpected structure, see diagnostic")
}

cat("\n=== 4. fixed effects and their variance ===\n")
b_a <- as.numeric(a$coefficients$fixed)
b_r <- as.numeric(f$beta)
check("fixed effect (intercept)", abs(b_r[1] - b_a[length(b_a)]) < 1e-5,
      sprintf("%.9f vs %.9f", b_r[1], b_a[length(b_a)]))

cat(sprintf("\n=== %d checks, %d failure(s) | time asreml %.2f s, remlax %.2f s ===\n",
            ok + ko, ko, t_a, t_r))
if (ko > 0) quit(status = 1)
