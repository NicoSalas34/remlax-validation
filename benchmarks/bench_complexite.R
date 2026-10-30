# =============================================================================
# bench_complexite.R: speed benchmark AT EQUAL MODEL on TWO axes, size (n) and
#                     model complexity (number of variance parameters,
#                     structures), for asreml, remlax CPU and remlax GPU
# =============================================================================
# WHAT IS MEASURED, per fit.
#   wall_s      wall clock of the complete call (R -> result)
#   startup_s   remlax: wall - compile - solver (Python launch, JAX import,
#               bundle write/read); asreml: NA
#   compile_s   remlax: XLA compilation, paid ONCE; asreml: NA
#   solver_s    remlax: optimisation only (L-BFGS-B + polishing); asreml: wall
#   n_iter      remlax: L-BFGS-B iterations; asreml: AI iterations
#   n_eval      remlax: evaluations of -2logL and gradient (line search and
#               polishing included); asreml: = n_iter
#   s_per_iter  solver_s / n_iter
#   s_per_eval  remlax: solver_s / n_eval; eval_s: one timed evaluation
#   logLik      asreml convention (fit$logLik_asreml), comparable to a$loglik
#   converged, n_par, n_obs, components
#
# CASES (n_par grows at fixed size):
#   iid     one random factor + covariate                                  2
#   grm     same with dense VanRaden kinship (q = n/4)                     2
#   ar1ar1  nr x nc field, iid genotype, ar1(row):ar1(col) residual        4
#   us3, us6, us9, us12   t traits, genetic us and residual us             t(t+1)
#   usK3, usK6            same with dense kinship on the genetic term      t(t+1)
#   ige     one-species reduction of the model of Salas et al. (2026): str(dge + neighbourhood) us(2) on GRM=I,
#           iid interspecific IGE, random ar1 x ar1 field, nugget          8
#
# TIME LIMIT PER FIT (--timeout-fit, s). remlax: `timeout` prefix on the Python
# command (RX_PY = "timeout <s> python3"); asreml: fork + kill. A stopped fit is
# recorded with status "timeout" and the limit in wall_s.
#
#   Rscript benchmarks/bench_complexite.R --sizes 500,2000,8000 \
#       --cases iid,grm,ar1ar1,us3,us6,us9,us12,usK3,usK6,ige \
#       --software asreml,remlax --tag cpu4 --timeout-fit 7200 --out benchmarks/results
#   RX_BACKEND=gpu Rscript benchmarks/bench_complexite.R --software remlax --tag a100
#
# Each measurement is printed as soon as it is taken (line "BENCH|..."), and the
# CSV is rewritten at each row, so a killed job leaves its data.
# =============================================================================
suppressPackageStartupMessages({ library(here); library(Matrix); library(jsonlite) })
source(here::here("R", "charger_remlax.R"))

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  i <- match(paste0("--", name), args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
SIZES <- as.integer(strsplit(opt("sizes", "500,2000"), ",")[[1]])
CASES <- strsplit(opt("cases", "iid,grm,ar1ar1,us3,us6,us9,us12,usK3,usK6,ige"), ",")[[1]]
SOFT  <- strsplit(opt("software", "asreml,remlax"), ",")[[1]]
TAG   <- opt("tag", "cpu")
OUT   <- opt("out", here::here("benchmarks", "results"))
TIMEOUT <- as.numeric(opt("timeout-fit", "7200"))
# Replicates per genotype: the DENSITY of the problem. q = n / REPS for iid, grm,
# ar1ar1 and ige; q = (n / t) / REPS for the multi-trait cases. REPS = 1 is
# estimable only with kinship (grm, usK) or neighbourhood (ige).
REPS <- as.numeric(opt("reps", "4"))
BK <- Sys.getenv("RX_BACKEND", "cpu")
MAXIT_REF <- 100L
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DATE <- format(Sys.Date(), "%Y-%m-%d")
OUT_FILE <- file.path(OUT, sprintf("benchc_%s_%s.csv", DATE, TAG))  # TAG must carry reps if != 4

# remlax: time limit through the timeout prefix on the Python command
if ("remlax" %in% SOFT && is.finite(TIMEOUT) && TIMEOUT > 0) {
  py <- rx_python_cmd()
  if (py[1] != "timeout") Sys.setenv(RX_PY = paste(c("timeout", as.integer(TIMEOUT), py), collapse = " "))
}
if ("asreml" %in% SOFT) {
  if (!requireNamespace("asreml", quietly = TRUE)) { cat("asreml not available\n"); SOFT <- setdiff(SOFT, "asreml") }
  else { suppressPackageStartupMessages(library(asreml)); asreml.options(trace = FALSE, workspace = "8gb", pworkspace = "2gb") }
}
sha <- Sys.getenv("RX_GIT_SHA", "")
if (!nzchar(sha)) sha <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(here::here()), "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE))[1], error = function(e) "?")
jaxv <- tryCatch(system2(rx_python_cmd()[length(rx_python_cmd())], c("-c", shQuote("import jax; print(jax.__version__)")), stdout = TRUE)[1], error = function(e) "?")
versions <- c(remlax = paste0(if (is.na(sha) || !nzchar(sha)) "?" else sha, "/jax", jaxv),
              asreml = if ("asreml" %in% SOFT) as.character(packageVersion("asreml")) else NA)
cat("versions:", paste(names(versions), versions, sep = "=", collapse = " | "), "| backend", BK,
    "| cores", Sys.getenv("SLURM_CPUS_PER_TASK", "?"), "| OMP", Sys.getenv("OMP_NUM_THREADS", "?"),
    "| timeout", TIMEOUT, "s\n")

COLS <- c("date", "tag", "backend", "case", "n", "reps", "q", "n_par", "software", "version", "status",
          "wall_s", "startup_s", "compile_s", "solver_s", "n_iter", "n_eval", "s_per_iter", "s_per_eval", "eval_s",
          "converged", "logLik", "components", "error")
RES <- list()
add_row <- function(...) {
  r <- list(...)
  for (k in setdiff(COLS, names(r))) r[[k]] <- NA
  r <- r[COLS]
  RES[[length(RES) + 1L]] <<- as.data.frame(r, stringsAsFactors = FALSE)
  cat("BENCH|", paste(sapply(r, function(v) if (is.null(v) || length(v) == 0) "" else as.character(v)), collapse = "|"), "\n", sep = "")
  write.csv(do.call(rbind, RES), OUT_FILE, row.names = FALSE)
}
comp_str <- function(v) paste(sprintf("%.6g", as.numeric(v)), collapse = ";")
timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; list(v = v, s = proc.time()[["elapsed"]] - t0) }
grm <- function(q, seed) {
  set.seed(seed); M <- matrix(rbinom(q * 2 * q, 2, 0.3), q, 2 * q)
  pfr <- colMeans(M) / 2; W <- sweep(M, 2, 2 * pfr)
  K <- tcrossprod(W) / (2 * sum(pfr * (1 - pfr))) + diag(1e-3, q)
  dimnames(K) <- list(paste0("g", 1:q), paste0("g", 1:q)); K
}
n_traits <- function(case_) as.integer(sub("^usK?", "", case_))

# ---- design generation ----------------------------------------------------------
gen_case <- function(case_, n) {
  set.seed(1000 + n)
  if (case_ %in% c("iid", "grm")) {
    q <- max(2L, round(n / REPS)); gid <- factor(rep(paste0("g", 1:q), length.out = n), levels = paste0("g", 1:q))
    K <- if (case_ == "grm") grm(q, n) else NULL
    u <- if (is.null(K)) rnorm(q, 0, 1.2) else as.numeric(t(chol(K)) %*% rnorm(q)) * 1.2
    x <- rnorm(n); y <- 5 + 0.4 * x + u[as.integer(gid)] + rnorm(n, 0, 1)
    list(d = data.frame(y = y, x = x, gid = gid), K = K, q = q)
  } else if (grepl("^usK?[0-9]+$", case_)) {
    tt <- n_traits(case_); n1 <- n %/% tt; q <- max(2L, round(n1 / REPS))
    gid <- factor(rep(paste0("g", 1:q), length.out = n1), levels = paste0("g", 1:q))
    K <- if (grepl("^usK", case_)) grm(q, n) else NULL
    G <- crossprod(matrix(rnorm(tt * tt), tt)) / tt + diag(0.5, tt)
    Rm <- crossprod(matrix(rnorm(tt * tt), tt)) / tt + diag(0.5, tt)
    Uw <- matrix(rnorm(q * tt), q, tt)
    if (!is.null(K)) Uw <- t(chol(K)) %*% Uw
    U <- Uw %*% chol(G)
    Y <- U[as.integer(gid), ] + matrix(rnorm(n1 * tt), n1, tt) %*% chol(Rm) + matrix(seq_len(tt), n1, tt, byrow = TRUE)
    colnames(Y) <- paste0("y", seq_len(tt))
    dw <- data.frame(gid = gid, Y)
    dl <- data.frame(y = c(Y), gid = rep(gid, tt), trait = factor(rep(colnames(Y), each = n1), levels = colnames(Y)),
                     unit = rep(seq_len(n1), tt))
    list(d = dw, dl = dl, t = tt, K = K, q = q)
  } else if (case_ == "ar1ar1") {
    nr <- round(sqrt(n / 1.25)); nc <- n %/% nr; ntot <- nr * nc; q <- max(2L, round(ntot / REPS))
    g <- expand.grid(col = seq_len(nc), row = seq_len(nr))
    Kr <- 0.6 ^ abs(outer(1:nr, 1:nr, "-")); Kc <- 0.4 ^ abs(outer(1:nc, 1:nc, "-"))
    e <- as.numeric(kronecker(t(chol(Kr)), t(chol(Kc))) %*% rnorm(ntot)) * 1.0
    gid <- factor(sample(rep(paste0("g", 1:q), length.out = ntot)), levels = paste0("g", 1:q))
    y <- 4 + rnorm(q, 0, 1.0)[as.integer(gid)] + e + rnorm(ntot, 0, 0.5)
    d <- data.frame(g, gid = gid, y = y); d$row <- factor(d$row); d$col <- factor(d$col)
    d <- d[order(d$row, d$col), ]
    list(d = d, dims = c(nr, nc), q = q)
  } else if (case_ == "ige") {
    # NR x NC grid, two species alternating by column: A (focal) on the odd
    # columns, B on the even ones. n plots of A. First-order neighbourhood
    # (Chebyshev window), weight 1/d, no normalisation.
    NC <- 2L * round(sqrt(n / 2 / 1.25)); NR <- ceiling(2 * n / NC)
    g <- expand.grid(col = seq_len(NC), row = seq_len(NR)); g$id <- seq_len(nrow(g))
    g$sp <- ifelse(g$col %% 2L == 1L, "A", "B")
    QA <- max(2L, round(sum(g$sp == "A") / REPS)); QB <- max(2L, round(sum(g$sp == "B") / REPS))
    gA <- paste0("a", seq_len(QA)); gB <- paste0("b", seq_len(QB))
    g$geno <- NA_character_
    g$geno[g$sp == "A"] <- sample(rep_len(gA, sum(g$sp == "A")))
    g$geno[g$sp == "B"] <- sample(rep_len(gB, sum(g$sp == "B")))
    iA <- which(g$sp == "A"); nA <- length(iA)
    idx <- matrix(NA_integer_, NR, NC); idx[cbind(g$row, g$col)] <- g$id
    mk_Z <- function(neighbour_sp, glev) {
      ii <- integer(0); jj <- integer(0); xx <- numeric(0)
      for (dr in -1:1) for (dc in -1:1) {
        if (dr == 0 && dc == 0) next
        r2 <- g$row[iA] + dr; c2 <- g$col[iA] + dc
        ok <- r2 >= 1 & r2 <= NR & c2 >= 1 & c2 <= NC
        j <- rep(NA_integer_, nA); j[ok] <- idx[cbind(r2[ok], c2[ok])]
        ok <- ok & !is.na(j) & g$sp[ifelse(is.na(j), 1L, j)] == neighbour_sp
        ii <- c(ii, which(ok)); jj <- c(jj, match(g$geno[j[ok]], glev)); xx <- c(xx, rep(1 / sqrt(dr^2 + dc^2), sum(ok)))
      }
      Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(nA, length(glev)), dimnames = list(NULL, glev))
    }
    ZgA <- Matrix::sparseMatrix(i = seq_len(nA), j = match(g$geno[iA], gA), x = 1, dims = c(nA, QA), dimnames = list(NULL, gA))
    ZnA <- mk_Z("A", gA); ZxA <- mk_Z("B", gB)
    SIG <- matrix(c(0.90, -0.35, -0.35, 0.30), 2, 2); V_INTER <- 0.22; V_AR1 <- 0.45; RHO_R <- 0.65; RHO_C <- 0.40; V_E <- 0.50
    U <- matrix(rnorm(QA * 2), QA, 2) %*% chol(SIG); uX <- rnorm(QB, 0, sqrt(V_INTER))
    col_obs <- as.integer(factor(g$col[iA])); NC_OBS <- max(col_obs)
    Kr <- RHO_R ^ abs(outer(1:NR, 1:NR, "-")); Kc <- RHO_C ^ abs(outer(1:NC_OBS, 1:NC_OBS, "-"))
    field <- as.numeric(kronecker(t(chol(Kr)), t(chol(Kc))) %*% rnorm(NR * NC_OBS)) * sqrt(V_AR1)
    cellA_obs <- (g$row[iA] - 1) * NC_OBS + col_obs
    y <- 3 + as.numeric(ZgA %*% U[, 1]) + as.numeric(ZnA %*% U[, 2]) + as.numeric(ZxA %*% uX) + field[cellA_obs] + rnorm(nA, 0, sqrt(V_E))
    d <- data.frame(y = y, geno = factor(g$geno[iA], levels = gA), row = factor(g$row[iA]), col = factor(col_obs))
    list(d = d, ZgA = ZgA, ZnA = ZnA, ZxA = ZxA, gA = gA, gB = gB, cell = cellA_obs, dims = c(NR, NC_OBS), QA = QA, QB = QB, q = QA)
  } else stop("unknown case ", case_)
}

# ---- one fit per software ----------------------------------------------------------
fit_remlax <- function(case_, D) {
  ctl <- list(backend = BK, verbose = FALSE, hessian = FALSE, blups = FALSE)
  if (case_ == "iid")  f <- do.call(rx_reml, c(list(y ~ x, random = ~ gid, data = D$d), ctl))
  if (case_ == "grm")  f <- do.call(rx_reml, c(list(y ~ x, random = ~ vm(gid, K = D$K), data = D$d), ctl))
  if (grepl("^us[0-9]+$", case_))
    f <- do.call(rx_reml, c(list(y ~ 0 + trait, random = ~ us(gid), residual = ~ us(trait):unit, data = D$dl, trait = "trait", unit = "unit"), ctl))
  if (grepl("^usK[0-9]+$", case_))
    f <- do.call(rx_reml, c(list(y ~ 0 + trait, random = ~ us(gid, K = D$K), residual = ~ us(trait):unit, data = D$dl, trait = "trait", unit = "unit"), ctl))
  if (case_ == "ar1ar1") f <- do.call(rx_reml, c(list(y ~ 1, random = ~ gid, residual = ~ ar1(row):ar1(col), data = D$d), ctl))
  if (case_ == "ige") {
    nA <- nrow(D$d)
    mod <- rx_model(y = D$d$y, X = matrix(1, nA, 1),
                    terms = list(rx_term("dge_ige", list(D$ZgA, D$ZnA), struct = "us", levels = D$gA),
                                 rx_term("ige_inter", list(D$ZxA), struct = "iid", levels = D$gB),
                                 rx_term("field", Matrix::sparseMatrix(i = seq_len(nA), j = D$cell, x = 1, dims = c(nA, prod(D$dims))),
                                         t = 1L, struct = "iid", level = "ar1ar1", dims = D$dims)),
                    residual = rx_residual("iid"))
    f <- do.call(rx_fit, c(list(mod), ctl))
  }
  comp <- c(unlist(lapply(f$sigmas, function(S) S[lower.tri(S, diag = TRUE)])), f$sigma_res[lower.tri(f$sigma_res, diag = TRUE)], unlist(f$rho))
  list(logLik = f$logLik_asreml, comp = comp, n_par = f$n_par, solver_s = f[["seconds"]], compile_s = f[["compile_s"]],
       n_eval = f[["n_eval"]], n_iter = f[["n_iter"]], eval_s = f[["eval_s"]],
       converged = isTRUE(as.logical(f$scipy_success)), version = versions[["remlax"]])
}
fit_asreml_raw <- function(case_, D) {
  d <- D$d
  if (case_ == "iid") a <- asreml(y ~ x, random = ~ gid, data = d, maxit = MAXIT_REF)
  if (case_ == "grm") { K <- D$K; a <- asreml(y ~ x, random = ~ vm(gid, K), data = d, maxit = MAXIT_REF) }
  if (grepl("^usK?[0-9]+$", case_)) {
    fml <- stats::as.formula(paste("cbind(", paste(grep("^y", names(d), value = TRUE), collapse = ","), ") ~ trait"))
    if (grepl("^usK", case_)) { K <- D$K; a <- asreml(fml, random = ~ us(trait):vm(gid, K), residual = ~ id(units):us(trait), data = d, maxit = MAXIT_REF) }
    else a <- asreml(fml, random = ~ us(trait):gid, residual = ~ id(units):us(trait), data = d, maxit = MAXIT_REF)
  }
  if (case_ == "ar1ar1") a <- asreml(y ~ 1, random = ~ gid, residual = ~ ar1(row):ar1(col), data = d, maxit = MAXIT_REF)
  if (case_ == "ige") {
    da <- cbind(d, as.data.frame(as.matrix(D$ZnA)), as.data.frame(as.matrix(D$ZxA)))
    names(da) <- c(names(d), paste0("N", seq_len(D$QA)), paste0("X", seq_len(D$QB)))
    # asreml evaluates `group` with eval(call$group) outside this function, so
    # the indices must be visible from the global environment.
    assign("iN", grep("^N[0-9]+$", names(da)), envir = globalenv())
    assign("iX", grep("^X[0-9]+$", names(da)), envir = globalenv())
    a <- asreml(y ~ 1, random = stats::as.formula(sprintf("~ str(~ geno + grp(nb), ~ us(2):id(%d)) + grp(xb) + ar1(row):ar1(col)", D$QA)),
                group = list(nb = iN, xb = iX), data = da, maxit = MAXIT_REF)
  }
  vc <- summary(a)$varcomp
  it <- tryCatch({ m <- a$monitor; if (is.matrix(m)) ncol(m) - 1L else if (is.matrix(a$trace)) nrow(a$trace) else NA_integer_ }, error = function(e) NA_integer_)
  list(logLik = a$loglik, comp = vc$component, n_par = nrow(vc), n_eval = it, n_iter = it,
       converged = isTRUE(a$converge), version = versions[["asreml"]])
}
fit_asreml <- function(case_, D) {
  # fork + time limit: setTimeLimit does not interrupt asreml
  if (!is.finite(TIMEOUT) || TIMEOUT <= 0) return(fit_asreml_raw(case_, D))
  p <- parallel::mcparallel(suppressWarnings(fit_asreml_raw(case_, D)))
  r <- parallel::mccollect(p, wait = FALSE, timeout = TIMEOUT)
  if (is.null(r)) { tools::pskill(p$pid, tools::SIGKILL); try(parallel::mccollect(p, wait = FALSE), silent = TRUE); stop("TIMEOUT") }
  r <- r[[1]]
  if (inherits(r, "try-error")) stop(paste(as.character(r), collapse = " "))
  if (inherits(r, "error")) stop(conditionMessage(r))
  r
}
FITS <- list(remlax = fit_remlax, asreml = fit_asreml)

# ---- loop -------------------------------------------------------------------------
for (n in SIZES) for (case_ in CASES) {
  D <- tryCatch(gen_case(case_, n), error = function(e) e)
  if (inherits(D, "error")) { add_row(date = DATE, tag = TAG, backend = BK, case = case_, n = n, software = "generation", status = "error", error = conditionMessage(D)); next }
  n_obs <- if (is.null(D$dl)) nrow(D$d) else nrow(D$dl)
  cat(sprintf("\n### %s n=%d (%d observations) ###\n", case_, n, n_obs)); flush.console()
  for (sw in SOFT) {
    # The asreml licence server is shared ("All licenses in use"): wait and
    # retry, the clock restarting from zero at each attempt.
    for (attempt in 1:30) {
      t0 <- proc.time()[["elapsed"]]
      r <- tryCatch(timed(suppressWarnings(FITS[[sw]](case_, D))), error = function(e) e)
      if (!(inherits(r, "error") && grepl("icense", conditionMessage(r)))) break
      cat("  [asreml licence busy, new attempt in 60 s]\n"); Sys.sleep(60)
    }
    bk <- if (sw == "remlax") BK else "cpu"
    if (inherits(r, "error")) {
      msg <- conditionMessage(r); elapsed <- proc.time()[["elapsed"]] - t0
      # remlax killed by `timeout`: the solver fails after ~TIMEOUT seconds
      to <- grepl("TIMEOUT", msg) || (is.finite(TIMEOUT) && TIMEOUT > 0 && elapsed >= 0.97 * TIMEOUT)
      add_row(date = DATE, tag = TAG, backend = bk, case = case_, n = n_obs, reps = REPS, q = D$q, software = sw, version = versions[[sw]],
            status = if (to) "timeout" else "error", wall_s = round(elapsed, 1),
            error = substr(gsub("[\r\n|]+", " ", msg), 1, 300))
      next
    }
    v <- r$v
    solver <- if (is.null(v$solver_s)) NA else v$solver_s
    compile <- if (is.null(v$compile_s)) NA else v$compile_s
    n_it <- if (length(v$n_iter) == 1L && is.finite(v$n_iter) && v$n_iter > 0) v$n_iter else NA
    n_ev <- if (length(v$n_eval) == 1L && is.finite(v$n_eval) && v$n_eval > 0) v$n_eval else NA
    base <- if (sw == "remlax") solver else r$s
    add_row(date = DATE, tag = TAG, backend = bk, case = case_, n = n_obs, reps = REPS, q = D$q, n_par = v$n_par, software = sw, version = v$version, status = "ok",
          wall_s = round(r$s, 3),
          startup_s = if (sw == "remlax") round(r$s - compile - solver, 3) else NA,
          compile_s = if (is.na(compile)) NA else round(compile, 3),
          solver_s = if (sw == "remlax") round(solver, 3) else round(r$s, 3),
          n_iter = n_it, n_eval = n_ev,
          s_per_iter = if (is.na(n_it)) NA else round(base / n_it, 4),
          s_per_eval = if (is.na(n_ev)) NA else round(base / n_ev, 4),
          eval_s = if (is.null(v$eval_s)) NA else round(v$eval_s, 4),
          converged = v$converged, logLik = sprintf("%.6f", v$logLik), components = comp_str(v$comp), error = "")
  }
}
cat("\nwritten:", OUT_FILE, "(", length(RES), "rows )\n")
