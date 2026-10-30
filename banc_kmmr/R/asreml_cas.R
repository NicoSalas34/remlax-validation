# ==============================================================================
# Fit of one case by asreml-R 4.2, same model as remlax and kmmr.
#
#   two_stage : GY ~ -1 + ManEnv, fa(ManEnv, k):vm(Gen, Ki), weights,
#               residual fixed at 1 by asr_gaussian(dispersion = 1)
#               (asreml-R 4.2 manual, section 3.10, absolute weights)
#   one_stage : GY ~ 1 + Env. Per trial: at(Env):Rep, at(Env):Block, a RANDOM
#               ar1 x ar1 field (its variance, its two rho) and its own iid
#               residual, dsum(~ units | Env). G x E in
#               fa(Env, 1) or us(Env) x vm(Gen, Ki)
#   gblup     : GY ~ 1 + Cycle, vm(Gen, Ki), iid residual
#   multi_trt : GY ~ -1 + TrtCycle, fa(Trt, 2) or us(Trt) x vm(Gen, Ki),
#               dsum(~ units | Trt)
# No column is called `trait` (a reserved word of asreml).
# ==============================================================================

# dense inverse of K, in the format accepted by vm(): matrix with dimnames and
# the INVERSE attribute (manual, section 4.3.6). asreml computes log|K| itself.
inverse_asreml <- function(K) {
  Ki <- solve(K)
  dimnames(Ki) <- dimnames(K)
  attr(Ki, "INVERSE") <- TRUE
  Ki
}

# ar1 x ar1 field of one trial: row and column renumbered from 1 within the
# trial, complete levels 1..n (ar1 assumes contiguous, equally spaced levels,
# like the nr x nc grid of remlax). Outside the trial, level 1 is a
# placeholder, and at(Env, e) zeroes these rows.
trial_grid <- function(d, e) {
  ok <- d$Env == e
  r <- d$Row - min(d$Row[ok]) + 1L; cc <- d$Col - min(d$Col[ok]) + 1L
  nr <- max(r[ok]); nc <- max(cc[ok])
  r[!ok] <- 1L; cc[!ok] <- 1L
  list(R = factor(r, levels = seq_len(nr)), C = factor(cc, levels = seq_len(nc)))
}

# asreml trace: one line per iteration, wall-clock time to the second. It is
# captured, printed again for the log, and the time stamps are read from it.
asreml_trace <- function(expr) {
  out <- utils::capture.output(val <- force(expr), type = "output")
  cat(out, sep = "\n")
  # iteration line: number, logLik, ..., hh:mm:ss, sometimes "( 8 restrained)"
  re_wall <- "[0-9]{1,2}:[0-9]{2}:[0-9]{2}"
  it <- grep(paste0("^ *[0-9]+ +-?[0-9.]+ .* ", re_wall), out, value = TRUE)
  hms <- function(x) { p <- as.numeric(strsplit(x, ":")[[1]]); sum(p * c(3600, 60, 1)) }
  wall <- vapply(regmatches(it, regexpr(re_wall, it)), hms, 0, USE.NAMES = FALSE)
  list(val = val, wall = wall)
}

fit_asreml <- function(cs, maxit = 50L, workspace = NULL) {
  suppressPackageStartupMessages(library(asreml))
  # asreml really touches its workspace: 30gb suits the cluster nodes but kills
  # a 32 GB workstation. RX_ASREML_WORKSPACE overrides the default (e.g. "4gb").
  if (is.null(workspace)) workspace <- Sys.getenv("RX_ASREML_WORKSPACE",
    if (cs$type %in% c("gblup", "multi_trt", "sim")) "30gb" else "1gb")
  asreml.options(workspace = workspace, pworkspace = workspace, maxit = maxit, trace = TRUE)
  # threads: OMP_NUM_THREADS if set (the value that executer.R records), otherwise all
  threads <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", "")))
  if (!is.na(threads)) asreml.options(threads = threads)
  d <- cs$d
  # dsum requires data grouped by section. The row order does not change the
  # likelihood
  if (cs$type == "one_stage") d <- d[order(d$Env), ]
  if (cs$type == "multi_trt") d <- d[order(d$Trt), ]

  tK <- Sys.time()
  Ki <- inverse_asreml(cs$K)
  t_Ki <- as.numeric(difftime(Sys.time(), tK, units = "secs"))

  tp <- Sys.time()
  args <- list(data = d)
  if (cs$type == "two_stage") {
    args$fixed <- GY ~ -1 + ManEnv
    args$random <- stats::as.formula(sprintf("~ fa(ManEnv, %d):vm(Gen, Ki)", cs$rank))
    args$residual <- ~ units
    args$weights <- "weight"
    args$family <- asreml::asr_gaussian(dispersion = 1)
  } else if (cs$type == "one_stage") {
    fields <- character(0)
    for (e in levels(d$Env)) {
      g <- trial_grid(d, e)
      nR <- paste0("R_", e); nC <- paste0("C_", e)
      d[[nR]] <- g$R; d[[nC]] <- g$C
      fields <- c(fields, sprintf("at(Env, '%s'):ar1(%s):ar1(%s)", e, nR, nC))
    }
    ge <- if (cs$ge == "us") "us(Env)" else "fa(Env, 1)"
    args$data <- d
    args$fixed <- GY ~ 1 + Env
    args$random <- stats::as.formula(paste("~ at(Env):Rep + at(Env):Block +",
                                           paste(fields, collapse = " + "),
                                           sprintf("+ %s:vm(Gen, Ki)", ge)))
    args$residual <- ~ dsum(~ units | Env)
  } else if (cs$type == "sim" && cs$t == 1L) {
    args$fixed <- GY ~ 1
    args$random <- ~ vm(Gen, Ki)
    args$residual <- ~ units
  } else if (cs$type == "sim") {
    ge <- switch(cs$ge, us = "us(Env)", diag = "diag(Env)", sprintf("fa(Env, %d)", cs$rank))
    args$fixed <- GY ~ -1 + Env
    args$random <- stats::as.formula(sprintf("~ %s:vm(Gen, Ki)", ge))
    args$residual <- ~ dsum(~ units | Env)
  } else if (cs$type == "gblup") {
    args$fixed <- GY ~ 1 + Cycle
    args$random <- ~ vm(Gen, Ki)
    args$residual <- ~ units
  } else {
    ge <- if (cs$ge == "us") "us(Trt)" else "fa(Trt, 2)"
    args$fixed <- GY ~ -1 + TrtCycle
    args$random <- stats::as.formula(sprintf("~ %s:vm(Gen, Ki)", ge))
    args$residual <- ~ dsum(~ units | Trt)
  }
  t_prep <- as.numeric(difftime(Sys.time(), tp, units = "secs"))
  # the formulas are evaluated here, so Ki must be visible in this frame
  environment(args$random) <- environment()

  t0 <- Sys.time()
  tr <- asreml_trace(do.call(asreml::asreml, args))
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fit <- tr$val

  # wall-clock time stamps (to the second): start of the call, then end of each iteration
  m0 <- with(as.POSIXlt(t0), hour * 3600 + min * 60 + sec)
  wall <- tr$wall
  if (length(wall)) wall <- ifelse(wall < m0 - 1, wall + 86400, wall)   # past midnight
  n_it <- if (length(wall)) length(wall) else max(ncol(fit$trace) - 1L, NA_integer_)
  # trace span under 5 s: the second resolves nothing, steps set to NA
  if (length(wall) < 2 || wall[length(wall)] - wall[1] < 5) wall <- numeric(0)
  # cost of one iteration: median of the gaps between lines if they last at
  # least 5 s, otherwise the mean over the span 2..n (the second is too coarse)
  eval_s <- NA_real_
  if (length(wall) > 1) {
    di <- diff(wall)
    eval_s <- if (stats::median(di) >= 5) stats::median(di) else (wall[length(wall)] - wall[1]) / (length(wall) - 1)
  }
  timings <- c(inverse_K = t_Ki, preparation = t_prep, fit = sec,
               first_iteration = if (length(wall)) max(wall[1] - m0, 0) else NA_real_,
               later_iterations = if (length(wall) > 1) wall[length(wall)] - wall[1] else NA_real_,
               outside_iterations = if (length(wall)) max(sec - (wall[length(wall)] - m0), 0) else NA_real_)
  list(timings = timings, eval_s = eval_s, logLik = fit$loglik, n_iter = n_it,
       converged = isTRUE(fit$converge), seconds = sec,
       n_par = sum(fit$vparameters.con != "F"))
}
