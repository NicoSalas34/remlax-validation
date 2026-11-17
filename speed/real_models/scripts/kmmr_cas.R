# ==============================================================================
# Fit of one case by kmmr. The formulas follow checks/checks.R of the package
# and section 5 of the paper (maxit = 50, default tolConvLL 0.001).
# ==============================================================================

inverse_kmmr <- function(K) {
  # kmmr expects a Matrix object (dsyMatrix), like kinshipDROPS
  Ki <- solve(Matrix::Matrix(K, sparse = FALSE))
  attr(Ki, "logDet") <- sum(log(eigen(K, only.values = TRUE, symmetric = TRUE)$values))
  Ki
}

# kmmr timings: kmmr prints (cat) the duration of its setup ("Setup") and of
# each iteration. This output is captured, printed again for the log, and the
# durations are read from it.
kmmr_trace <- function(expr) {
  out <- utils::capture.output(val <- force(expr), type = "output")
  cat(out, sep = "\n")
  hms <- function(x) { p <- as.numeric(strsplit(x, ":")[[1]]); sum(p * c(3600, 60, 1)[seq_along(p)]) }
  s <- gsub("\033\\[[0-9;]*m", "", out)
  setup <- regmatches(s, regexpr("^Setup +NA +[0-9:.]+", s))
  setup <- if (length(setup)) hms(sub(".* ", "", setup)) else NA_real_
  it <- grep("^[0-9]+ +-?[0-9.]+ +[0-9:.]+", s, value = TRUE)
  di <- vapply(strsplit(trimws(it), " +"), function(z) hms(z[3]), 0)
  list(val = val, setup = setup, iterations = di)
}

fit_kmmr <- function(cs, gpu = FALSE, maxit = 50L, tol = 0.001) {
  d <- cs$d
  tK <- Sys.time()
  Ki <- inverse_kmmr(cs$K)
  t_Ki <- as.numeric(difftime(Sys.time(), tK, units = "secs"))
  t_sv <- 0
  assign("Ki", Ki, envir = globalenv())        # kmmr formulas look up Ki by name
  t0 <- Sys.time()
  if (cs$type == "two_stage") {
    f_r <- stats::as.formula(sprintf("~ re(~ exfa(ManEnv, %d):lk(Gen, Ki), vars = ~ ManEnv)", cs$rank))
    tsv <- Sys.time()
    sv <- kmmr::kmmr(fixed = GY ~ -1 + ManEnv, random = f_r, residual = ~ units, weights = "weight",
                     data = d, getStartingValues = TRUE, useGPU = gpu, verbose = FALSE)
    t_sv <- as.numeric(difftime(Sys.time(), tsv, units = "secs"))
    # residual fixed at 1 (weights = inverses of the BLUE variances), starting
    # variances rescaled as in checks.R
    init <- sv$Init; fixed <- sv$Fixed
    ir <- length(init)
    inc <- mean(1 / d$weight) / init[ir]
    init[ir] <- 1; fixed[ir] <- TRUE
    nE <- nlevels(d$ManEnv)
    init[seq_len(nE)] <- init[seq_len(nE)] * inc
    t0 <- Sys.time()
    tr <- kmmr_trace(kmmr::kmmr(fixed = GY ~ -1 + ManEnv, random = f_r, residual = ~ units, weights = "weight",
                      init = init, isFixed = fixed, data = d, maxit = maxit, useGPU = gpu, tolConvLL = tol))
  } else if (cs$type == "one_stage") {
    ge <- if (cs$ge == "us") "corg(Env)" else "exfa(Env, 1)"
    f_r <- stats::as.formula(paste(
      "~ atEach('Env', ~ re(~ id(Rep), vars = ~ 1))",
      "+ atEach('Env', ~ re(~ id(Block), vars = ~ 1))",
      "+ atEach('Env', ~ re(~ ar1(Row):ar1(Col), vars = ~ 1))",
      sprintf("+ re(~ %s:lk(Gen, Ki), vars = ~ Env)", ge)))
    d$Row <- factor(d$Row); d$Col <- factor(d$Col)
    tr <- kmmr_trace(kmmr::kmmr(fixed = GY ~ 1 + Env, random = f_r,
                      residual = ~ re(~ id(Env):id(Gen):id(Rep), vars = ~ Env),
                      data = d, maxit = maxit, useGPU = gpu, tolConvLL = tol))
  } else if (cs$type == "sim" && cs$t == 1L) {
    tr <- kmmr_trace(kmmr::kmmr(fixed = GY ~ 1, random = ~ re(~ lk(Gen, Ki), vars = ~ 1),
                      residual = ~ units, data = d, maxit = maxit, useGPU = gpu, tolConvLL = tol))
  } else if (cs$type == "sim") {
    # r = 1 in all multi-environment cases: (Env, Gen) identifies the
    # observation. No other column, because the kmmr fa starting values
    # aggregate over all remaining columns.
    d <- d[, c("GY", "Env", "Gen")]
    f_r <- if (cs$ge == "diag") {
      # id(Env):lk(Gen, Ki) fails in kmmr (chol2inv on dCHMsimpl)
      ~ atEach("Env", ~ re(~ lk(Gen, Ki), vars = ~ 1))
    } else {
      ge <- if (cs$ge == "us") "corg(Env)" else sprintf("exfa(Env, %d)", cs$rank)
      stats::as.formula(sprintf("~ re(~ %s:lk(Gen, Ki), vars = ~ Env)", ge))
    }
    tr <- kmmr_trace(kmmr::kmmr(fixed = GY ~ -1 + Env, random = f_r,
                      residual = ~ re(~ id(Env):id(Gen), vars = ~ Env),
                      data = d, maxit = maxit, useGPU = gpu, tolConvLL = tol))
  } else if (cs$type == "gblup") {
    tr <- kmmr_trace(kmmr::kmmr(fixed = GY ~ 1 + Cycle, random = ~ re(~ lk(Gen, Ki), vars = ~ 1),
                      residual = ~ units, data = d, maxit = maxit, useGPU = gpu, tolConvLL = tol))
  } else {
    ge <- if (cs$ge == "us") "corg(Trt)" else "exfa(Trt, 2)"
    tr <- kmmr_trace(kmmr::kmmr(fixed = GY ~ -1 + TrtCycle,
                      random = stats::as.formula(sprintf("~ re(~ %s:lk(Gen, Ki), vars = ~ Trt)", ge)),
                      residual = ~ re(~ id(Trt):id(Gen), vars = ~ Trt),
                      data = d, maxit = maxit, useGPU = gpu, tolConvLL = tol))
  }
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fit <- tr$val
  it <- tr$iterations
  timings <- c(inverse_K = t_Ki, start_values = t_sv,
               preparation = tr$setup,
               first_iteration = if (length(it)) it[1] else NA_real_,
               later_iterations = if (length(it) > 1) sum(it[-1]) else 0,
               outside_solver = sec - sum(c(tr$setup, it), na.rm = TRUE))
  list(timings = timings, eval_s = if (length(it) > 1) stats::median(it[-1]) else NA_real_,
       logLik = utils::tail(fit$LL, 1), n_iter = length(fit$LL), converged = isTRUE(fit$converged),
       seconds = sec, n_par = sum(fit$isFixed == 0))
}
