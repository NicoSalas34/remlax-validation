# ==============================================================================
# Running a case: fits, arbiter, comparisons, verdicts.
# ==============================================================================

ARBITER_TOL <- 1e-8     # relative difference of logLik from the arbiter
ENGINE_TOL <- 1e-4      # absolute difference of logLik from the best route
Z_TOL <- 3              # |estimate - true| / se

.chr <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(x)[1]
.num <- function(x) if (is.null(x) || !length(x)) NA_real_ else suppressWarnings(as.numeric(x)[1])
.lgl <- function(x) if (is.null(x) || !length(x)) NA else as.logical(x)[1]

fit_spec <- function(spec, model, dir) {
  t0 <- Sys.time()
  out <- switch(spec$type,
    rx_fit = {
      fit <- rx_fit(model, backend = "cpu", verbose = FALSE, optimizer = spec$optimizer, dir = dir)
      # Verdict of the package: Newton decrement or projected relative gradient,
      # after polishing. scipy_success covers only the L-BFGS-B / AI phase
      # before polishing. It is reported in the detail and decides nothing.
      conv <- isTRUE(as.logical(fit$conv_decrement)) || isTRUE(as.logical(fit$conv_grad_rel))
      list(par = par_from_fit(fit), logLik = fit$logLik, conv = conv,
           detail = sprintf("success=%s decrement=%s grad_rel=%s hessian_ok=%s at_bound=%s iter=%s",
                            .lgl(fit$scipy_success), .lgl(fit$conv_decrement), .lgl(fit$conv_grad_rel),
                            .lgl(fit$conv_hessian_ok), .num(fit$n_at_bound), .num(fit$n_iter)),
           newton_decrement = .num(fit$newton_decrement), grad_rel = .num(fit$grad_rel),
           hessian_ok = .lgl(fit$conv_hessian_ok), n_at_bound = .num(fit$n_at_bound),
           message = .chr(fit$scipy_message))
    },
    # MME engine: AI-REML (automatic or sparse inversion) or L-BFGS-B with
    # sparse inversion, the route that replaces the former RTMB engine.
    mme = {
      fit <- rx_fit(model, backend = "cpu", verbose = FALSE, optimizer = spec$optimizer,
                    engine = "mme", inversion = spec$inversion, dir = dir)
      conv <- isTRUE(as.logical(fit$conv_decrement)) || isTRUE(as.logical(fit$conv_grad_rel))
      list(par = par_from_fit(fit), logLik = fit$logLik, conv = conv,
           detail = sprintf("engine=%s inversion=%s success=%s decrement=%s grad_rel=%s at_bound=%s iter=%s",
                            .chr(fit$engine), .chr(fit$inversion), .lgl(fit$scipy_success),
                            .lgl(fit$conv_decrement), .lgl(fit$conv_grad_rel), .num(fit$n_at_bound),
                            .num(fit$n_iter)),
           newton_decrement = .num(fit$newton_decrement), grad_rel = .num(fit$grad_rel),
           hessian_ok = .lgl(fit$conv_hessian_ok), n_at_bound = .num(fit$n_at_bound),
           message = .chr(fit$scipy_message))
    })
  out$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out
}

arbitrate <- function(entry, par) {
  m <- entry$model
  V <- entry$V_of(par)
  if (!is.null(entry$M_of)) reml_lag_dense(V, m$y, as.matrix(m$X), entry$M_of(par))$logLik
  else reml_dense(V, m$y, as.matrix(m$X))$logLik
}

run_case <- function(name, builder, seed, tmp) {
  t_case <- Sys.time()
  case <- tryCatch(builder(seed), error = function(e) e)
  if (inherits(case, "error"))
    return(list(fits = data.frame(case = name, engine = "construction", model = NA, n = NA, n_par = NA,
                                  logLik = NA, logLik_arbiter = NA, rel_diff = NA,
                                  diff_between_engines = NA, convergence = NA, convergence_detail = NA,
                                  optimiser_message = NA,
                                  n_z_out = NA, max_abs_z = NA, seconds = NA,
                                  status = "ERROR", reasons = paste("construction:", conditionMessage(case)),
                                  stringsAsFactors = FALSE),
                comps = NULL, note = NA, family = NA))
  res <- list()
  filt <- if (exists("ENGINE_FILTER")) ENGINE_FILTER else "."
  for (k in seq_along(case$engines)) {
    spec <- case$engines[[k]]
    if (!grepl(filt, spec$id)) next
    entry <- case$models[[spec$model]]
    lab <- if (spec$model == "main") spec$id else paste0(spec$id, "@", spec$model)
    dir <- file.path(tmp, paste0(name, "_g", seed, "_", k))
    a <- tryCatch(fit_spec(spec, entry$model, dir), error = function(e) e)
    unlink(dir, recursive = TRUE)
    if (inherits(a, "error")) {
      msg <- gsub("[\r\n]+", " ", conditionMessage(a))
      expected <- !is.null(spec$expected_error) && grepl(spec$expected_error, msg)
      res[[lab]] <- list(error = msg, expected = expected, spec = spec)
      next
    }
    a$arbiter <- tryCatch(arbitrate(entry, a$par), error = function(e) structure(NA_real_, msg = conditionMessage(e)))
    a$spec <- spec; a$model <- entry$model
    res[[lab]] <- a
  }
  ok <- vapply(res, function(r) is.null(r$error), TRUE)
  best <- if (any(ok)) max(vapply(res[ok], `[[`, 0, "logLik")) else NA_real_
  # Fallback standard errors (engines without SE): those of lbfgsb on the same model.
  se_ref <- function(model) {
    for (r in res) if (is.null(r$error) && r$spec$type == "rx_fit" && r$spec$model == model) return(r$par$se)
    NULL
  }
  rows_out <- list(); comps <- list()
  for (lab in names(res)) {
    r <- res[[lab]]
    if (!is.null(r$error)) {
      rows_out[[lab]] <- data.frame(case = name, engine = lab, model = r$spec$model, n = NA, n_par = NA,
                                  logLik = NA, logLik_arbiter = NA, rel_diff = NA,
                                  diff_between_engines = NA, convergence = NA, convergence_detail = NA,
                                  optimiser_message = NA,
                                  n_z_out = NA, max_abs_z = NA, seconds = NA,
                                  status = if (r$expected) "EXPECTED_REFUSAL" else "ERROR",
                                  reasons = substr(r$error, 1, 300), stringsAsFactors = FALSE)
      next
    }
    arb <- as.numeric(r$arbiter)
    rel <- (r$logLik - arb) / abs(arb)
    gap <- best - r$logLik
    # components
    est <- r$par$est; se <- r$par$se
    ren <- r$spec$rename
    if (!is.null(ren)) { i <- names(est) %in% names(ren); names(est)[i] <- ren[names(est)[i]]
                         names(se) <- names(est) }
    if (!is.null(case$derive)) {
      dv <- case$derive(r$par, r$spec$type)
      if (!is.null(dv)) { est <- c(est, stats::setNames(dv$est, dv$name)); se <- c(se, stats::setNames(dv$se, dv$name)) }
    }
    sr <- se_ref(r$spec$model)
    vr <- case$truth
    vr <- vr[names(vr) %in% names(est)]
    e_v <- est[names(vr)]; s_v <- se[names(vr)]
    src <- ifelse(is.finite(s_v), "own", NA_character_)
    if (!is.null(sr)) { j <- !is.finite(s_v) & names(vr) %in% names(sr)
      s_v[j] <- sr[names(vr)[j]]; src[j & is.finite(s_v)] <- "lbfgsb" }
    z <- (e_v - vr) / s_v
    # true component on the bound: |z| is reported and decides nothing
    at_bound <- names(vr) %in% (case$bound %||% character(0))
    src[at_bound] <- paste0(ifelse(is.na(src[at_bound]), "", src[at_bound]), "+bound")
    comps[[lab]] <- data.frame(case = name, engine = lab, component = names(vr), true = unname(vr),
                               estimate = unname(e_v), se = unname(s_v), se_source = src, z = unname(z),
                               stringsAsFactors = FALSE)
    absent <- setdiff(names(case$truth), names(est))
    reasons <- c()
    if (!is.finite(arb)) reasons <- c(reasons, paste("arbiter unavailable:", attr(r$arbiter, "msg") %||% "?"))
    else if (abs(rel) >= ARBITER_TOL) reasons <- c(reasons, sprintf("logLik != arbiter (rel diff %.2e)", rel))
    if (is.finite(gap) && gap >= ENGINE_TOL) reasons <- c(reasons, sprintf("optimum lower by %.2e", gap))
    if (!isTRUE(r$conv)) reasons <- c(reasons, "convergence not established")
    zf <- z[is.finite(z) & !at_bound]
    nz <- sum(abs(zf) >= Z_TOL)
    if (nz) reasons <- c(reasons, sprintf("%d component(s) with |z| >= %g", nz, Z_TOL))
    if (length(absent) && r$spec$type == "rx_fit" && r$spec$model == "main")
      reasons <- c(reasons, paste("true components missing:", paste(absent, collapse = ", ")))
    rows_out[[lab]] <- data.frame(case = name, engine = lab, model = r$spec$model, n = length(r$model$y),
                                n_par = r$model$n_par, logLik = r$logLik, logLik_arbiter = arb,
                                rel_diff = rel, diff_between_engines = gap, convergence = r$conv,
                                convergence_detail = r$detail, optimiser_message = r$message,
                                n_z_out = nz,
                                max_abs_z = if (length(zf)) max(abs(zf)) else NA_real_,
                                seconds = r$seconds,
                                status = if (length(reasons)) "FAIL" else "PASS",
                                reasons = paste(reasons, collapse = "; "), stringsAsFactors = FALSE)
  }
  list(fits = do.call(rbind, rows_out), comps = do.call(rbind, comps), note = case$note,
       family = case$family, seconds = as.numeric(difftime(Sys.time(), t_case, units = "secs")))
}
