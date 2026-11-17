# ==============================================================================
# Translation of a case into a remlax rx_model, and fit along one route.
#
# Routes: dense (V form, low memory), blocks (V form by blocks), mme
# (inversion chosen automatically), mme_hybrid (dense separator + blocks).
# The V routes use L-BFGS-B with the trust-region finish. The MME routes use
# AI-REML, like kmmr. The old French route names (blocs, mme_creux,
# mme_hybride) stay accepted, because stored task lists use them.
# ==============================================================================

indicators <- function(f, keep = NULL) {
  f <- factor(f)
  i <- seq_along(f); j <- as.integer(f)
  if (!is.null(keep)) { i <- i[keep]; j <- j[keep] }
  Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(length(f), nlevels(f)),
                       dimnames = list(NULL, levels(f)))
}

# indicators of a factor on the `ok` rows only (levels present)
sub_indicators <- function(f, ok) {
  g <- droplevels(factor(f)[ok])
  Matrix::sparseMatrix(i = which(ok), j = as.integer(g), x = 1, dims = c(length(f), nlevels(g)))
}

# G x E term: t components (one per environment), q genotypes, kinship K
term_ge <- function(name, env, gen, K, struct, rank = 0L) {
  env <- factor(env)
  Zl <- lapply(levels(env), function(e) indicators(gen, env == e))
  names(Zl) <- levels(env)
  remlax::rx_term(name, Zl, K = K, struct = struct, rank = rank)
}

remlax_model <- function(cs) {
  d <- cs$d
  if (cs$type == "two_stage") {
    X <- stats::model.matrix(~ -1 + ManEnv, d)
    tm <- term_ge("ge", d$ManEnv, d$Gen, cs$K, "fa", cs$rank)
    return(remlax::rx_model(d$GY, X, list(tm), remlax::rx_residual("iid"), weights = d$weight))
  }
  if (cs$type == "one_stage") {
    X <- stats::model.matrix(~ 1 + Env, d)
    terms <- list()
    secs <- list()
    for (e in levels(d$Env)) {
      ok <- d$Env == e
      terms[[length(terms) + 1L]] <- remlax::rx_term(paste0("rep_", e), sub_indicators(d$Rep, ok), t = 1L)
      terms[[length(terms) + 1L]] <- remlax::rx_term(paste0("blk_", e), sub_indicators(d$Block, ok), t = 1L)
      r <- d$Row[ok] - min(d$Row[ok]); cc <- d$Col[ok] - min(d$Col[ok])
      nr <- max(r) + 1L; nc <- max(cc) + 1L
      Zc <- Matrix::sparseMatrix(i = which(ok), j = r * nc + cc + 1L, x = 1, dims = c(nrow(d), nr * nc))
      terms[[length(terms) + 1L]] <- remlax::rx_term(paste0("field_", e), Zc, t = 1L,
                                                       level = "ar1ar1", dims = c(nr, nc))
      secs[[length(secs) + 1L]] <- remlax::rx_residual("iid", rows = which(ok), name = e)
    }
    terms[[length(terms) + 1L]] <- term_ge("ge", d$Env, d$Gen, cs$K,
                                           if (cs$ge == "us") "us" else "fa", if (cs$ge == "us") 0L else 1L)
    return(remlax::rx_model(d$GY, X, terms, remlax::rx_residual("iid", sections = secs)))
  }
  if (cs$type == "sim") {
    # indicators on ALL kinship levels, observed or not
    ind <- function(ok) Matrix::sparseMatrix(i = which(ok), j = as.integer(d$Gen)[ok], x = 1,
                                             dims = c(nrow(d), nlevels(d$Gen)))
    if (cs$t == 1L) {
      tm <- remlax::rx_term("g", ind(rep(TRUE, nrow(d))), K = cs$K, t = 1L)
      return(remlax::rx_model(d$GY, matrix(1, nrow(d), 1), list(tm), remlax::rx_residual("iid")))
    }
    Zl <- lapply(levels(d$Env), function(e) ind(d$Env == e)); names(Zl) <- levels(d$Env)
    st <- if (startsWith(cs$ge, "fa")) "fa" else cs$ge
    tm <- remlax::rx_term("ge", Zl, K = cs$K, struct = st, rank = cs$rank)
    secs <- lapply(levels(d$Env), function(e) remlax::rx_residual("iid", rows = which(d$Env == e), name = e))
    return(remlax::rx_model(d$GY, stats::model.matrix(~ -1 + Env, d), list(tm),
                            remlax::rx_residual("iid", sections = secs)))
  }
  if (cs$type == "gblup") {
    X <- stats::model.matrix(~ 1 + Cycle, d)
    tm <- remlax::rx_term("g", indicators(d$Gen), K = cs$K, t = 1L)
    return(remlax::rx_model(d$GY, X, list(tm), remlax::rx_residual("iid")))
  }
  X <- stats::model.matrix(~ -1 + TrtCycle, d)
  tm <- term_ge("ge", d$Trt, d$Gen, cs$K, if (cs$ge == "us") "us" else "fa", if (cs$ge == "us") 0L else 2L)
  secs <- lapply(levels(d$Trt), function(e) remlax::rx_residual("iid", rows = which(d$Trt == e), name = e))
  remlax::rx_model(d$GY, X, list(tm), remlax::rx_residual("iid", sections = secs))
}

fit_remlax <- function(cs, route, backend = "cpu", maxiter = 30000L) {
  t0 <- Sys.time()
  m <- remlax_model(cs)
  t_model <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  args <- switch(route,
    dense       = list(engine = "dense", memory = "low", optimizer = "lbfgsb"),
    # fully fixed routes of the absolute benchmark (engine x inversion x optimiser)
    V_ai        = list(engine = "dense", memory = "standard", optimizer = "ai"),
    V_lbfgsb    = list(engine = "dense", memory = "standard", optimizer = "lbfgsb"),
    mme_dense   = list(engine = "mme", inversion = "dense", optimizer = "ai"),
    blocks      = ,
    blocs       = list(engine = "blocks", optimizer = "lbfgsb"),
    mme         = list(engine = "mme", inversion = "auto", optimizer = "ai"),
    mme_sparse  = ,
    mme_creux   = list(engine = "mme", inversion = "sparse", optimizer = "ai"),
    mme_hybrid  = ,
    mme_hybride = list(engine = "mme", inversion = "hybrid", optimizer = "ai"),
    stop("unknown route: ", route))
  common <- list(model = m, backend = backend, hessian = FALSE, blups = FALSE, verbose = TRUE,
                 maxiter = maxiter)
  fixed <- NULL
  if (!is.null(m$weights)) {
    # residual fixed at 1 (weights = inverses of the BLUE variances). theta_init
    # holds NA everywhere except the residual (log sd 0). The NA entries take
    # the automatic start of the solver, without a preliminary call.
    th <- rep(NA_real_, m$n_par); th[m$n_par] <- 0
    common$theta_init <- th
    fixed <- m$n_par
  }
  t1 <- Sys.time()
  r <- do.call(remlax::rx_fit, c(common, list(fixed_theta = fixed), args))
  sec <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  list(logLik = r$logLik_asreml %||% r$logLik, logLik_remlax = r$logLik, n_iter = r$n_iter %||% NA,
       n_eval = r$n_eval %||% NA, converged = isTRUE(as.logical(r[["scipy_success"]])),
       seconds = sec, solver_seconds = r$seconds %||% NA, n_par = m$n_par,
       engine = r$engine %||% NA, inversion = r$inversion %||% NA,
       t_model = t_model, compile_s = r$compile_s %||% NA, eval_s = r$eval_s %||% NA,
       timings = r[["chrono"]])
}
