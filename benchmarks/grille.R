# =============================================================================
# ONE cell of the grid, in ONE process. Nothing else.
# =============================================================================
# WHY ONE PROCESS PER CELL. The earlier sweep measured every cell in one
# process, and its compilation column proved impossible to interpret: it fell
# from 79 s to 2.7 s on the last cells, a sign that an already compiled program
# was reused from one cell to the next. The order of execution thus changed the
# result. A fresh process per cell removes the cause instead of correcting the
# effect, at the price of one startup per measurement, which is precisely what
# this benchmark measures.
#
# THE THREE AXES, AND WHY THEY ARE SEPARABLE.
#   n  observations         -> cost of the dense factorisation, in n^3
#   q  random effects       -> cost of the sparse factorisation and assembly
#   p  estimated parameters -> NUMBER OF EVALUATIONS, not their unit cost
# p varies with the number of traits t under a us structure on the genetic
# term, the residual staying diag: p = t(t+1)/2 + t. Raising p through the
# RESIDUAL would confound it with the structure of R.
#
# THE FOUR CELLS. The first token is the CASE, the second the ENGINE.
#   sparse-sparse  ar1 field                     -> sparse engine (its ground)
#   sparse-dense   ar1 field                     -> dense engine
#   dense-dense    full genomic kinship          -> dense engine  (its ground)
#   dense-sparse   full genomic kinship          -> sparse engine THE WORST CASE
# The "sparse" engine is the remlax mixed model equations (MME) engine, with
# sparse inversion and the AI optimiser: rx_fit(engine = "mme",
# inversion = "sparse", optimizer = "ai"). It replaces the RTMB engine, removed
# from remlax in October 2026.
# The last cell forces a full kinship into the sparse factorisation, through
# the Kinv slot, to put a number on the cost of the wrong choice.
suppressPackageStartupMessages(library(Matrix))
a <- commandArgs(TRUE)
arg <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1L] }
CASE <- arg("case", "sparse"); ENG <- arg("engine", "dense"); BK <- arg("backend", "cpu")
N   <- as.integer(arg("n", "2000")); Q <- as.integer(arg("q", "500"))
TT  <- as.integer(arg("t", "2")); OUT <- arg("out", "grid.csv"); TAG <- arg("tag", "")
SEED <- as.integer(arg("seed", "20260901"))
MAXIT <- as.integer(arg("maxiter", "3000"))
REPO <- arg("repo", ".")
source(file.path(REPO, "R/charger_remlax.R"))
set.seed(SEED)

inc <- function(f, q) Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f),
                                           x = 1, dims = c(length(f), q))

# --- the design: t traits in long format, one incidence per trait
n_unit <- max(50L, N %/% TT); n <- n_unit * TT
lev <- if (CASE == "sparse") rep_len(seq_len(Q), n_unit) else sample.int(Q, n_unit, replace = TRUE)
trait <- rep(seq_len(TT), each = n_unit)
Zt <- inc(factor(rep(lev, TT), levels = seq_len(Q)), Q)
Zl <- lapply(seq_len(TT), function(k) { Z <- Zt; Z[trait != k, ] <- 0; Z })
y <- 2 + rnorm(n, 0, 1.2); X <- model.matrix(~ factor(trait))

K <- Ki <- NULL
if (CASE == "dense") {
  # dense genomic kinship after VanRaden, on biallelic markers
  M <- matrix(rbinom(Q * 500L, 2L, 0.3), Q, 500L)
  pf <- colMeans(M) / 2; W <- sweep(M, 2, 2 * pf)
  K <- (W %*% t(W)) / (2 * sum(pf * (1 - pf)))
  K <- K + diag(Q) * 1e-6 * mean(diag(K))
}

tm <- if (CASE == "sparse") {
  rx_term("g", Zl, struct = "us", level = "ar1")
} else if (ENG == "sparse") {
  # Kinv is FULL: that is the point of this cell.
  rx_term("g", Zl, struct = "us",
          Kinv = methods::as(methods::as(solve(K), "CsparseMatrix"), "generalMatrix"))
} else {
  rx_term("g", Zl, struct = "us", K = K)
}
m <- rx_model(y, X, terms = list(tm),
              residual = rx_residual(trait = factor(trait), struct = "diag"))

t0 <- proc.time()[["elapsed"]]
# SAME cap and same outputs on both sides, otherwise a truncation passes for a
# convergence. The MME engine runs on CPU (sparse factorisation).
r <- try(if (ENG == "sparse") rx_fit(m, backend = "cpu", engine = "mme", inversion = "sparse",
                                    optimizer = "ai", maxiter = MAXIT, hessian = FALSE,
                                    blups = FALSE, verbose = FALSE)
         else rx_fit(m, backend = BK, maxiter = MAXIT, hessian = FALSE,
                     blups = FALSE, verbose = FALSE),
         silent = TRUE)
wall <- proc.time()[["elapsed"]] - t0

row <- data.frame(
  cell = paste0(CASE, "-", ENG), case = CASE, engine = ENG,
  backend = if (ENG == "dense") BK else "cpu",
  n = n, q = Q, t = TT, tag = TAG,
  p = if (inherits(r, "try-error")) NA_integer_ else as.integer(r$n_par),
  # t(t+1)/2 for the us of the term, t for the diag residual, and ONE more
  # in the sparse case: the ar1 field carries its correlation.
  p_expected = as.integer(TT * (TT + 1) / 2 + TT + (CASE == "sparse")),
  n_iter = if (inherits(r, "try-error")) NA_integer_ else as.integer(r$n_iter %||% NA),
  # WALL against INTERNAL: both engines go through the Python solver, outside the
  # R process, and the wall time contains this invocation cost. Both are
  # recorded. Only the internal time compares the algebra; only the wall time
  # tells what the user experiences.
  fit_wall_s = wall,
  fit_internal_s = if (inherits(r, "try-error")) NA_real_
                   else as.numeric(r$seconds %||% wall),
  # Both engines report their own compilation and the cost of one evaluation.
  # The first call pays the compilation, so first - next = compilation. For
  # the MME engine, eval_s is the median of its complete evaluations.
  compile_s = if (inherits(r, "try-error")) NA_real_ else as.numeric(r$compile_s %||% NA),
  eval_s = if (inherits(r, "try-error")) NA_real_ else as.numeric(r$eval_s %||% NA),
  logLik = if (inherits(r, "try-error")) NA_real_ else as.numeric(r$logLik),
  maxiter = MAXIT,
  # A cell that reaches the cap is NOT a fit time.
  truncated = !inherits(r, "try-error") && !is.null(r$n_iter) && r$n_iter >= MAXIT,
  status = if (inherits(r, "try-error")) "fail" else "ok",
  error = if (inherits(r, "try-error"))
             substr(conditionMessage(attr(r, "condition")), 1, 160) else "",
  stringsAsFactors = FALSE)

write.table(row, OUT, sep = ",", row.names = FALSE, append = file.exists(OUT),
            col.names = !file.exists(OUT))
cat(sprintf("[%s-%s/%s] n=%-6d q=%-5d t=%-3d p=%-4s iter=%-5s wall %9.2f s | internal %9.2f s %s\n",
            CASE, ENG, if (ENG == "dense") BK else "cpu", n, Q, TT,
            row$p, row$n_iter, row$fit_wall_s, row$fit_internal_s, row$error))
