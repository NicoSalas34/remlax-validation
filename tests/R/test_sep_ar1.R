# Does a three-factor separable product reach the solver?
# WHY THIS TEST. `ar1ar1` covers only TWO factors. A spatial field replicated
# independently by block is written id (x) ar1 (x) ar1. The identity factor
# carries no parameter, so the two correlations are SHARED between the
# replicates. The IGE_analysis engine estimates this model, and the R interface
# could not express it.
suppressMessages(source(file.path(Sys.getenv("RX_DIR", "."), "R", "charger_remlax.R")))
set.seed(4)
nb <- 3L; nr <- 4L; nc <- 5L; q <- nb * nr * nc
n <- 2L * q
# one observation per cell, duplicated, with the engine index:
# lin = tray*nr*nc + col*nr + row  -> the ROW varies fastest
cell <- rep(seq_len(q), 2L)
Z <- Matrix::sparseMatrix(i = seq_len(n), j = cell, x = 1, dims = c(n, q))
y <- rnorm(n); X <- matrix(1, n, 1)
tm <- rx_term("field", Z, struct = "iid", t = 1L, level = "sep",
              parts = list(list("id", nb), list("ar1", nc), list("ar1", nr)))
cat("n_par of the term:", tm$n_par, " (expected 3: one variance + two correlations)\n")
m <- rx_model(y, X, terms = list(field = tm),
              residual = rx_residual("iid", unit = seq_len(n)))
f <- rx_fit(m, hessian = FALSE, blups = FALSE, maxiter = 60L, verbose = FALSE)
cat("logLik :", format(f$logLik, digits = 8), "\n")
cat("theta  :", paste(round(f$theta, 4), collapse = " "), "\n")
cat("reported correlations:", paste(names(f$rho), signif(unlist(f$rho), 4),
                                    collapse = " | "), "\n")
# THE TWO CORRELATIONS IN PLAIN FORM: theta[2] is the COLUMN correlation and
# theta[3] the ROW correlation, in this order, because the engine index makes
# the row vary fastest. They are reported on the scale used, tanh(theta).
cat(sprintf("rho columns = %.4f | rho rows = %.4f\n",
            tanh(f$theta[2]), tanh(f$theta[3])))
# theta has 4 entries: 3 for the term (one variance + two correlations) and
# one for the residual. My first assertion miscounted.
ok <- length(f$theta) == 4L && tm$n_par == 3L && is.finite(f$logLik) &&
      # the R side flattens the keys into "<term>!<family>_<index>_<param>".
      # The factor index is in the name, which tells WHICH of the two
      # correlations belongs to the rows.
      length(grep("^field!ar1_", names(f$rho))) == 2L
cat(if (ok) "\nok: the three-factor product reaches the solver\n"
    else "\nFAIL\n")
if (!ok) quit(status = 1L)
