# =============================================================================
# The decision grid: {sparse case, dense case} x {dense CPU, dense GPU, sparse CPU}
# =============================================================================
# The "sparse" engine is the remlax mixed model equations (MME) engine, with
# sparse inversion and the AI optimiser: rx_fit(engine = "mme",
# inversion = "sparse", optimizer = "ai"). It replaces the RTMB engine, removed
# from remlax in October 2026.
#
# SIX cells and not eight. The sparse factorisation runs on CPU, so the cell
# "sparse x GPU" DOES NOT EXIST. This absence makes the "sparse" and "GPU"
# choices exclusive. The decision is therefore a trade-off, and not an
# optimisation along two independent axes.
#
# The two cases are defined by STRUCTURE, never by size:
#   SPARSE: iid factor + ar1 field, whose inverses are diagonal and tridiagonal
#   DENSE:  genomic kinship, whose inverse is FULL
#
# The most instructive cell is the DENSE case forced on the sparse engine. It is
# forced here, through a full Kinv, to measure the size of the cost. This is the
# quantitative answer to "why not use sparse everywhere", and it calibrates the
# decision tool instead of leaving it on an uncalibrated flop ratio.
suppressPackageStartupMessages(library(Matrix))
`%||%` <- function(a, b) if (is.null(a)) b else a
source("R/charger_remlax.R")

a <- commandArgs(TRUE)
arg <- function(k, d) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1L] }
BK   <- arg("backend", "cpu")            # cpu | gpu, for the DENSE engine
CELL <- strsplit(arg("cells", "sparse-dense,sparse-sparse,dense-dense,dense-sparse"), ",")[[1]]
NS   <- as.integer(strsplit(arg("ns", "2000,4000"), ",")[[1]])
QS   <- as.integer(strsplit(arg("qs", "1000,2000"), ",")[[1]])
REPS <- as.integer(arg("reps", "3"))
OUT  <- arg("out", "bench_engines.csv")
TAG  <- arg("tag", "")
# --- sweep over traits, at constant n PER TRAIT ("traits-*" cells)
TS   <- as.integer(strsplit(arg("ts", "1,2,3,4,5,6,7,8"), ",")[[1]])
NU   <- as.integer(arg("n-unit", "1000"))      # observations per trait
QT   <- as.integer(arg("q-trait", "300"))      # levels of the genetic term
KIN  <- identical(arg("kinship", "no"), "yes")
# The kinship is built ONCE: it does not depend on t, and rebuilding it per
# cell would change the data between the points of the curve.
KMAT <- NULL
if (KIN) {
  set.seed(4242L)
  M <- matrix(rbinom(QT * 1000L, 2L, 0.3), QT, 1000L)
  fr <- colMeans(M) / 2; W <- sweep(M, 2, 2 * fr)
  KMAT <- (W %*% t(W)) / (2 * sum(fr * (1 - fr)))
  KMAT <- KMAT + diag(QT) * 1e-6 * mean(diag(KMAT))
}

# rx_fit exposes backend = c("auto", "gpu", "cpu"). It is set explicitly,
# never "auto", so that a GPU cell falling back to the CPU is visible.
fit_dense <- function(m, ...) rx_fit(m, backend = BK, ...)
# The MME engine with sparse inversion, on CPU.
fit_sparse <- function(m, ...) rx_fit(m, backend = "cpu", engine = "mme",
                                     inversion = "sparse", optimizer = "ai", ...)

# BIAS NOT TO MEASURE. Both engines go through the Python solver, OUTSIDE THE R
# PROCESS. Each call pays this invocation, and the first one the compilation.
# Comparing these wall times would give a false conclusion on the cost of the
# computation. The INTERNAL time reported by the solver itself (field
# `seconds`) is therefore recorded, with the wall time in a separate column.
# The gap between the two IS the invocation cost. It matters on its own, since
# a fit pays it once and a model comparison loop pays it at each call.
# WARM-UP OUTSIDE THE MEASUREMENT. The first call of a cell pays build costs
# that are not the measured computation, compilation included. The trap, as
# measured with the old RTMB engine: it showed 0.209 s at t = 1 then 0.012 s at
# t = 2, a factor of 17 against the expected slope. The median of two
# repetitions still contained the build. An untimed call therefore precedes
# the repetitions.
timed <- function(f, reps) {
  invisible(try(f(), silent = TRUE))
  ts <- numeric(reps)
  for (k in seq_len(reps)) { t0 <- proc.time()[["elapsed"]]; v <- f(); ts[k] <- proc.time()[["elapsed"]] - t0 }
  list(med = median(ts), min = min(ts), max = max(ts), val = v)
}

# --- the two design families ---------------------------------------------------
case_sparse <- function(n) {
  q1 <- n %/% 5L; nr <- 40L; nc <- n %/% (2L * nr)
  qf <- nr * nc
  f1 <- factor(sample.int(q1, n, replace = TRUE), levels = seq_len(q1))
  f2 <- factor(rep(seq_len(qf), length.out = n), levels = seq_len(qf))
  y <- 2 + rnorm(n, 0, 1.2); X <- matrix(1, n, 1)
  list(y = y, X = X,
       terms_d = list(rx_term("a", f1, struct = "iid"),
                      rx_term("s", f2, struct = "iid", level = "ar1")),
       terms_s = list(rx_term("a", f1, struct = "iid"),
                      rx_term("s", f2, struct = "iid", level = "ar1")),
       residual = rx_residual(), q = q1 + qf)
}

case_dense <- function(n, q) {
  # dense genomic kinship after VanRaden, on biallelic markers
  M <- matrix(rbinom(q * 500L, 2L, 0.3), q, 500L)
  p <- colMeans(M) / 2; W <- sweep(M, 2, 2 * p)
  K <- (W %*% t(W)) / (2 * sum(p * (1 - p)))
  K <- K + diag(q) * 1e-6 * mean(diag(K))
  lev <- factor(rep(seq_len(q), length.out = n), levels = seq_len(q))
  y <- 2 + rnorm(n, 0, 1.2); X <- matrix(1, n, 1)
  # Kinv is FULL: that is the point. It is supplied to force the sparse
  # factorisation of the MME engine on a dense kinship.
  Ki <- as(as(solve(K), "CsparseMatrix"), "generalMatrix")
  list(y = y, X = X,
       terms_d = list(rx_term("g", lev, K = K, struct = "iid")),
       terms_s = list(rx_term("g", lev, Kinv = Ki, struct = "iid")),
       residual = rx_residual(), q = q, K = K)
}

# --- MODEL COMPLEXITY AT CONSTANT n PER TRAIT ----------------------------------
# The design a breeder actually runs: the SAME plants, measured on more and more
# traits. It deliberately couples n and p, where the earlier sweep held n
# constant. Together, the two sweeps separate the effects.
#
# us(t) carries t(t+1)/2 parameters (the variances AND the covariances), that is
# 36 at t = 8. The residual is diag(t) so that p is dominated by the genetic
# term and stays readable. A us residual would double the count without
# changing the axis.
#
# The incidence is a LIST of t n x q matrices, one per column of Sigma, as
# rx_term expects. All observations of a unit share its level.
case_traits <- function(n_unit, t, q, K = NULL) {
  n <- n_unit * t
  lev <- rep(seq_len(q), length.out = n_unit)
  # one matrix per trait: row (i-1)*t + j carries trait j
  Zl <- lapply(seq_len(t), function(j) {
    M <- Matrix::sparseMatrix(i = seq_len(n_unit), j = lev, x = 1,
                              dims = c(n_unit, q))
    # duplicate each unit t times, but activate only the block of trait j
    idx <- (seq_len(n_unit) - 1L) * t + j
    Z <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                              dims = c(n, q))
    Z[idx, ] <- M
    Z
  })
  trait <- factor(rep(seq_len(t), times = n_unit))
  unit  <- rep(seq_len(n_unit), each = t)
  # genetic signal correlated between traits, so that us has something to
  # estimate. Without a real covariance, the optimum lies on the boundary and
  # the number of iterations measures a degeneracy instead of a difficulty.
  A <- matrix(rnorm(t * t), t, t) / sqrt(t)
  G <- A %*% t(A) + diag(t) * 0.5
  U <- matrix(rnorm(q * t), q, t) %*% chol(G)
  if (!is.null(K)) U <- chol(K) %*% U
  y <- 2 + U[cbind(lev[unit], as.integer(trait))] + rnorm(n, 0, 0.8)
  # t = 1 is the requested starting point, and model.matrix refuses a factor
  # with a single level. The overall mean is then enough.
  X <- if (t == 1L) matrix(1, n, 1) else model.matrix(~ trait)
  Ki <- if (is.null(K)) NULL else
        as(as(solve(K), "CsparseMatrix"), "generalMatrix")
  list(y = y, X = X, q = q, t = t, n = n,
       terms_d = list(rx_term("g", Zl, K = K, struct = "us")),
       terms_s = list(rx_term("g", Zl, Kinv = Ki, struct = "us")),
       # At t = 1 there is only one residual variance: `diag` has nothing to
       # diagonalise and the solver rightly refuses it.
       residual = if (t == 1L) rx_residual()
                  else rx_residual(struct = "diag", trait = trait, unit = unit))
}


rows <- list()
note <- function(...) cat(..., "\n", sep = "")

for (cl in CELL) {
  parts <- strsplit(cl, "-")[[1]]; case_ <- parts[1]; eng <- parts[2]
  grid <- if (case_ == "traits") lapply(TS, function(tt) list(n = NU * tt, q = QT, t = tt))
            else if (case_ == "sparse") lapply(NS, function(n) list(n = n, q = NA))
            else unlist(lapply(NS, function(n) lapply(QS, function(q) list(n = n, q = q))), recursive = FALSE)
  for (g in grid) {
    set.seed(1000L + g$n + (if (is.na(g$q)) 0L else g$q))
    # Same safeguard as the Python benchmark: below 2q units some levels are not
    # replicated, the system becomes nearly singular, and the extra cost would
    # read as an effect of complexity when the design is degenerate.
    if (identical(case_, "traits") && 2L * g$q > NU) {
      note(sprintf("[%s] t=%d SKIPPED: 2q=%d > n_unit=%d, insufficient replication",
                   cl, g$t, 2L * g$q, NU)); next
    }
    # q >= n: more effects than observations. The factor then carries
    # UNOBSERVED levels, which rx_term does not count, hence a K of size q
    # against a smaller term. The case is degenerate (the variance is not
    # identified). It is skipped and SAID, instead of interrupting the sweep.
    if (identical(case_, "traits") && g$q * 2L > NU) {
      note(sprintf("[%s] t=%d SKIPPED: q=%d > n_unit/2=%d, insufficient replication",
                   cl, g$t, g$q, NU %/% 2L)); next
    }
    if (!identical(case_, "traits") && !is.na(g$q) && g$q >= g$n) {
      note(sprintf("[%s] n=%d q=%d SKIPPED: q >= n, degenerate case", cl, g$n, g$q))
      next
    }
    d <- if (case_ == "traits") case_traits(NU, g$t, g$q, K = if (KIN) KMAT else NULL)
         else if (case_ == "sparse") case_sparse(g$n) else case_dense(g$n, g$q)
    tms <- if (eng == "dense") d$terms_d else d$terms_s
    m <- rx_model(d$y, d$X, terms = tms, residual = d$residual, name = cl)
    # IMPOSED THETA, IDENTICAL FOR BOTH ENGINES. Without it, each engine
    # evaluates at ITS initial theta and the logLik column cannot be compared
    # between engines. Trap checked: the values diverged by 1e-5 to 9e-3,
    # growing with t, which looked like a formulation defect. At a COMMON theta
    # the two engines agree to 1e-12. Zero on the transformed scale (variances
    # at 1, zero correlations) is the natural choice: well conditioned,
    # reproducible, and independent of the engine. The two parameter counts are
    # identical, as checked by the agreement suite.
    th0 <- rep(0, rx_n_theta(m))
    r <- tryCatch({
      # Both engines are called the same way, with the same outputs
      # (hessian = FALSE, blups = FALSE) and the same iteration cap. A call
      # asymmetry would read as an engine difference.
      fit_eng <- if (eng == "dense") fit_dense else fit_sparse
      ev <- timed(function() fit_eng(m, theta_init = th0, maxiter = 0L,
                                      polish = 0L, hessian = FALSE,
                                      blups = FALSE, verbose = FALSE), REPS)
      fi <- timed(function() fit_eng(m, hessian = FALSE, blups = FALSE,
                                      verbose = FALSE), 1L)
      list(ev = ev$val$seconds %||% NA_real_, ev_wall = ev$med,
           ev_min = ev$min, ev_max = ev$max, ll = ev$val$logLik,
           # THE LIKELIHOOD AT THE OPTIMUM, distinct from the one at the imposed
           # theta. Without it, comparing fit TIMES between engines makes no
           # sense: an engine that stops at a worse point is not faster, it is
           # worse. The ll column is the value at the common imposed theta and
           # does not answer this question.
           ll_fit = fi$val$logLik %||% NA_real_,
           fit = fi$val$seconds %||% NA_real_, fit_wall = fi$med,
           n_iter = fi$val$n_iter %||% NA_integer_,
           # A COUNT, not an inference. Dividing the fit time by the unit cost
           # gives a plausible and wrong number of evaluations. That error
           # made the rebuild of the cost model fail by a factor of two. The
           # solver counts them.
           n_eval = fi$val$n_eval %||% NA_integer_,
           compile_s = fi$val$compile_s %||% NA_real_, err = "")
    }, error = function(e) list(ev = NA_real_, ev_wall = NA_real_,
                                ev_min = NA_real_, ev_max = NA_real_,
                                ll = NA_real_, fit = NA_real_, fit_wall = NA_real_,
                                n_iter = NA_integer_, n_eval = NA_integer_,
                                compile_s = NA_real_, ll_fit = NA_real_,
                                err = substr(conditionMessage(e), 1, 200)))
    rows[[length(rows) + 1L]] <- data.frame(
      cell = cl, case = case_, engine = eng, backend = if (eng == "dense") BK else "cpu",
      n = g$n, t = if (is.null(g$t)) 1L else g$t, q_gen = g$q, q_total = d$q,
      p_us = if (is.null(g$t)) NA_integer_ else g$t * (g$t + 1L) / 2L, reps = REPS,
      eval_s = r$ev, eval_wall_s = r$ev_wall,
      eval_min_s = r$ev_min, eval_max_s = r$ev_max,
      fit_s = r$fit, fit_wall_s = r$fit_wall, n_iter = r$n_iter,
      n_eval = r$n_eval, compile_s = r$compile_s,
      logLik = r$ll, logLik_fit = r$ll_fit,
      rss_mb = as.numeric(gsub("[^0-9]", "", system("grep VmRSS /proc/self/status", intern = TRUE))) / 1024,
      tag = TAG, error = r$err, stringsAsFactors = FALSE)
    note(sprintf("[%s] t=%-2s n=%-6d q=%-5s eval %8.4f s (wall %8.3f) | fit %8.2f s | iter %-4s %s",
                 cl, ifelse(is.null(g$t), "-", g$t), g$n, ifelse(is.na(g$q), "-", g$q),
                 r$ev %||% NA_real_, r$ev_wall %||% NA_real_, r$fit %||% NA_real_,
                 r$n_iter, r$err))
  }
}
df <- do.call(rbind, rows)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.csv(df, OUT, row.names = FALSE)
note(sprintf("%d measurements written to %s", nrow(df), OUT))
