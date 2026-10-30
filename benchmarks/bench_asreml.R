# =============================================================================
# EXTERNAL TEST: the sparse / dense partition of ASReml, on the same designs
# =============================================================================
# WHY THIS FILE EXISTS
# --------------------
# Our central claim is that a DENSE kinship should not go through the sparse
# route, and that structure decides more than size. Measured on our two
# engines, it holds for OUR implementation. A reviewer will ask whether it
# generalises.
#
# ASReml can settle this in the best possible way, because the contrast is
# INTERNAL to the same program. Its manual 4.2, appendix A.1, describes the
# partition of the terms into a dense set and a sparse set. The terms of the
# fixed formula go into the dense set, those of the random and sparse formulae
# into the sparse set. The inverse coefficient matrix is FULLY formed for the
# dense set and only PARTIALLY for the sparse set. The option
#     asreml.options(dense = ~vm(...))
# moves a RANDOM term into the dense set.
#
# The SAME model is therefore compared, on the SAME data, in the same program,
# with the genomic term treated as sparse (default) then as dense. Same
# optimiser, same stopping rule, same code: only the solution route changes.
# None of our usual biases applies: no diverging paths between backends, no
# out-of-process invocation cost, no implementation difference.
#
# WHAT THE MANUAL ALREADY SAYS, TO BE CITED RATHER THAN REDISCOVERED
# ------------------------------------------------------------------
# Appendix A.1 itself recommends the dense treatment for kinship matrices: they
# are large, of the order of several thousands, and full, "and it can be more
# efficient to process such terms as dense". Our thesis is therefore stated by
# the reference program. This test PUTS A NUMBER ON IT.
#
# IT ALSO EXPLAINS WHY remlax DOES NOT RETURN THE PEV
# ---------------------------------------------------
# "the variance matrix of the BLUEs and BLUPs is only available for terms in the
# dense portion": the partial inverse does not give the prediction error
# variances of the sparse terms. This is the price of the partial inverse, and
# ASReml pays the same one, with the same output.
#
# WHAT THIS FILE DOES NOT MEASURE
# -------------------------------
# It does NOT compare the speed of remlax with ASReml. They are two different
# algorithms (average information against quasi-Newton), invoked differently,
# and the manual sets maxit = 13 by default where we use 3000, so the
# iterations are not comparable objects. Only the INTERNAL sparse/dense
# contrast is measured here, and it compares with OUR internal contrast, ratio
# to ratio.
# =============================================================================
suppressMessages({ library(asreml); library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  i <- match(paste0("--", name), args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
N_UNIT <- as.integer(opt("n-unit", "2000"))
QS     <- as.integer(strsplit(opt("qs", "200,500,1000"), ",")[[1]])
REPS   <- as.integer(opt("reps", "2"))
OUT    <- opt("out", "asreml_sparse_dense.csv")
TAG    <- opt("tag", "")
MAXIT  <- as.integer(opt("maxit", "30"))
# The default workspace is about 134 MB and saturates from q = 2000:
# "Insufficient workspace available when reordering matrices".
WS     <- opt("workspace", "16gb")

cat(sprintf("[asreml] version %s | n_unit = %d | q = %s | maxit = %d\n",
            as.character(packageVersion("asreml")), N_UNIT,
            paste(QS, collapse = ","), MAXIT))

# --- design: DENSE genomic kinship, built after VanRaden ----------------------
# The same construction as the genomic suite of bench.py, so that both series
# deal with the same statistical object.
design <- function(n_unit, q, n_marq = 1500, seed = 7) {
  set.seed(seed)
  p <- runif(n_marq, 0.05, 0.95)
  M <- matrix(rbinom(q * n_marq, 2, rep(p, each = q)), q, n_marq)
  Z <- scale(M, center = 2 * p, scale = FALSE)
  K <- tcrossprod(Z) / (2 * sum(p * (1 - p)))
  K <- K + diag(1e-6, q)                       # minimal regularisation
  gen <- factor(rep_len(seq_len(q), n_unit))
  u <- as.numeric(chol(K) %*% rnorm(q))
  y <- 1.5 + u[as.integer(gen)] + rnorm(n_unit, 0, 0.8)
  list(d = data.frame(y = y, gen = gen), K = K, q = q, n = n_unit)
}

# ASReml wants the INVERSE, with the INVERSE attribute set (manual 4.2 p.67).
inverse_set <- function(K) {
  Ki <- solve(K)
  dimnames(Ki) <- list(seq_len(nrow(K)), seq_len(nrow(K)))
  attr(Ki, "INVERSE") <- TRUE
  Ki
}

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(eval.parent(substitute(expr)), error = function(e) e)
  list(r = r, s = proc.time()[["elapsed"]] - t0)
}

rows <- list()
for (q in QS) {
  cat(sprintf("\n=== q = %d ===\n", q))
  dd <- design(N_UNIT, q)
  Ki <- inverse_set(dd$K)

  for (route in c("sparse", "dense")) {
    for (rep in seq_len(REPS)) {
      # The partition is set by asreml.options: by default a random term is
      # SPARSE; dense = ~vm(...) switches it.
      if (route == "dense") {
        asreml.options(dense = ~ vm(gen, Ki), trace = FALSE, maxit = MAXIT,
                       workspace = WS, pworkspace = WS)
      } else {
        asreml.options(dense = ~ NULL, trace = FALSE, maxit = MAXIT,
                       workspace = WS, pworkspace = WS)
      }
      z <- timed(asreml(fixed = y ~ 1, random = ~ vm(gen, Ki),
                        residual = ~ units, data = dd$d))
      if (inherits(z$r, "error")) {
        rows[[length(rows) + 1L]] <- data.frame(
          tag = TAG, q = q, n = dd$n, route = route, rep = rep, status = "fail",
          seconds = NA_real_, n_iter = NA_integer_, logLik = NA_real_,
          converged = NA, pev_available = NA, pev_se_median = NA_real_,
          error = substr(conditionMessage(z$r), 1, 180))
        cat(sprintf("  %-6s rep %d : FAIL %s\n", route, rep,
                    substr(conditionMessage(z$r), 1, 90)))
        next
      }
      a <- z$r
      # WHAT DOES THE OBJECT CARRY? Reading length(a$loglik) as an iteration
      # count returned 1 everywhere, so it is not the history. The fields are
      # printed ONCE instead of guessed, and the count is extracted robustly.
      if (!exists(".fields_printed")) {
        cat("  [fields of the asreml object] ",
            paste(names(a), collapse = " "), "\n")
        .fields_printed <<- TRUE
      }
      n_it <- suppressWarnings(as.integer(
        # Fields actually present (printed by this script): ifault converge
        # nedf nwv nsing noeff loglik sigma2 ... trace ... There is neither
        # `nit` nor `monitor`; the history is in `trace`, one column per
        # iteration.
        if (!is.null(a$trace)) ncol(as.matrix(a$trace))
        else if (length(a$loglik) > 1) length(a$loglik)
        else NA_integer_))
      # The prediction error variance is available ONLY for the dense part
      # (manual A.1). This is the second result of this test, and it is
      # checked instead of assumed.
      pev <- tryCatch({
        pv <- predict(a, classify = "vm(gen, Ki)", only = "vm(gen, Ki)",
                      sed = FALSE, trace = FALSE)
        # This contradicts my reading of appendix A.1, which says the BLUP
        # variance is available only for the dense part: predict() returned it
        # on BOTH sides. Either it recomputes it on demand whatever the
        # partition, or the term was not really switched. The median standard
        # error is therefore recorded, to see whether both routes return the
        # SAME thing or merely something.
        se <- pv$pvals$std.error
        if (is.null(se) || !any(is.finite(se))) FALSE else median(se, na.rm = TRUE)
      }, error = function(e) FALSE)
      rows[[length(rows) + 1L]] <- data.frame(
        tag = TAG, q = q, n = dd$n, route = route, rep = rep, status = "ok",
        seconds = z$s, n_iter = n_it, logLik = tail(a$loglik, 1),
        converged = isTRUE(a$converge), pev_available = !identical(pev, FALSE),
        pev_se_median = if (is.numeric(pev)) pev else NA_real_, error = "")
      cat(sprintf("  %-6s rep %d : %7.2f s | %3s iter | logLik %.6f | PEV %s\n",
                  route, rep, z$s, ifelse(is.na(n_it), "?", n_it),
                  tail(a$loglik, 1), pev))
    }
  }
}

res <- do.call(rbind, rows)
write.csv(res, OUT, row.names = FALSE)
cat(sprintf("\n[asreml] %d rows written to %s\n", nrow(res), OUT))

# --- the internal ratio, which is the object of the test ---------------------
ok <- res[res$status == "ok", ]
if (nrow(ok)) {
  cat("\n=== sparse / dense ratio, INTERNAL to ASReml ===\n")
  for (q in unique(ok$q)) {
    s <- ok[ok$q == q, ]
    tc <- median(s$seconds[s$route == "sparse"])
    td <- median(s$seconds[s$route == "dense"])
    lc <- unique(round(s$logLik[s$route == "sparse"], 6))
    ld <- unique(round(s$logLik[s$route == "dense"], 6))
    cat(sprintf("  q=%-5d sparse %7.2f s | dense %7.2f s | ratio %5.2f | ",
                q, tc, td, tc / td))
    cat(sprintf("same logLik: %s\n",
                isTRUE(all.equal(lc[1], ld[1], tolerance = 1e-6))))
  }
  cat("\nA ratio > 1 means the DENSE route is faster, as the manual states\n",
      "for a genomic kinship. The logLik values must be identical: the\n",
      "partition changes the solution route, not the model.\n", sep = "")
}
