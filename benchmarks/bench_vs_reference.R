# =============================================================================
# bench_vs_reference.R: speed benchmark AT EQUAL MODEL: remlax against asreml,
#                       lme4 and sommer, same designs, same data
# =============================================================================
# WHAT IS MEASURED. For each case and each size: the total time of the fit call
# (wall clock), the number of iterations or evaluations each program reports,
# the time per evaluation derived from it, the logLik at the optimum (asreml
# convention for remlax: fit$logLik_asreml, the only one comparable to
# a$loglik; lme4 and sommer on their own constants) and the final variance
# components. The components check that the compared fits reached the same
# point.
#
# WHAT IS NOT COMPARABLE. The iterations: average information (asreml,
# sommer), quasi-Newton (remlax, L-BFGS-B) and derivative-free (lme4, bobyqa)
# do not count the same thing. The total time is comparable, and so is the
# optimum.
#
# CASES. iid      one random factor (q = n/4 genotypes, 4 replicates) + covariate
#        grm      the same with a VanRaden kinship matrix (dense q x q)
#        us3/us6  t traits on n/t units, genetic us and residual us
#        ar1ar1   nr x nc field, iid genotype, ar1(row):ar1(col) residual
#
#   Rscript benchmarks/bench_vs_reference.R --sizes 500,2000 --cases iid,grm,us3,us6,ar1ar1 \
#           --refs asreml,lme4,sommer --tag cpu4 --out benchmarks/results
#   RX_BACKEND=gpu Rscript benchmarks/bench_vs_reference.R --refs none --tag a100
#
# Each measurement is printed as soon as it is taken (line "BENCH|..."), so
# that a job killed by the walltime leaves its data in the log.
# =============================================================================
suppressPackageStartupMessages({ library(here); library(Matrix); library(jsonlite) })
source(here::here("R", "charger_remlax.R"))

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  i <- match(paste0("--", name), args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
SIZES <- as.integer(strsplit(opt("sizes", "500,2000"), ",")[[1]])
CASES <- strsplit(opt("cases", "iid,grm,us3,us6,ar1ar1"), ",")[[1]]
REFS  <- setdiff(strsplit(opt("refs", "asreml,lme4,sommer"), ",")[[1]], "none")
TAG   <- opt("tag", "cpu")
OUT   <- opt("out", here::here("benchmarks", "results"))
SOMMER_MAX_N <- as.integer(opt("sommer-max-n", "2000"))
BK <- Sys.getenv("RX_BACKEND", "cpu")
MAXIT_REF <- 100L
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DATE <- format(Sys.Date(), "%Y-%m-%d")
OUT_FILE <- file.path(OUT, sprintf("bench_%s_%s.csv", DATE, TAG))

for (p in REFS) if (!requireNamespace(p, quietly = TRUE)) { cat("package not available:", p, "\n"); REFS <- setdiff(REFS, p) }
if ("asreml" %in% REFS) { suppressPackageStartupMessages(library(asreml)); asreml.options(trace = FALSE, workspace = "6gb", pworkspace = "2gb") }
if ("lme4" %in% REFS) suppressPackageStartupMessages(library(lme4))
if ("sommer" %in% REFS) suppressPackageStartupMessages(library(sommer))
sha <- Sys.getenv("RX_GIT_SHA", "")
if (!nzchar(sha)) sha <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(here::here()), "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE))[1], error = function(e) "?")
versions <- c(remlax = paste0(if (is.na(sha) || !nzchar(sha)) "?" else sha, "/jax", tryCatch(system2(rx_python_cmd(), c("-c", shQuote("import jax; print(jax.__version__)")), stdout = TRUE)[1], error = function(e) "?")),
              sapply(REFS, function(p) as.character(packageVersion(p))))
cat("versions:", paste(names(versions), versions, sep = "=", collapse = " | "), "| backend", BK,
    "| cores", Sys.getenv("SLURM_CPUS_PER_TASK", "?"), "| OMP", Sys.getenv("OMP_NUM_THREADS", "?"), "\n")

COLS <- c("date", "tag", "backend", "case", "n", "n_par", "software", "version", "status",
          "wall_s", "solver_s", "compile_s", "n_eval", "s_per_eval", "converged", "logLik", "components", "error")
n_eval_ok <- function(v) length(v) == 1L && is.finite(v) && v > 0
RES <- list()
add_row <- function(...) {
  r <- list(...); r <- r[COLS[COLS %in% names(r)]]
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

# ---- design generation ----------------------------------------------------------
gen_case <- function(case_, n) {
  set.seed(1000 + n)
  if (case_ %in% c("iid", "grm")) {
    q <- n %/% 4; gid <- factor(rep(paste0("g", 1:q), each = 4), levels = paste0("g", 1:q))
    K <- if (case_ == "grm") grm(q, n) else NULL
    u <- if (is.null(K)) rnorm(q, 0, 1.2) else as.numeric(t(chol(K)) %*% rnorm(q)) * 1.2
    x <- rnorm(n); y <- 5 + 0.4 * x + u[as.integer(gid)] + rnorm(n, 0, 1)
    list(d = data.frame(y = y, x = x, gid = gid), K = K)
  } else if (case_ %in% c("us3", "us6")) {
    tt <- if (case_ == "us3") 3L else 6L; n1 <- n %/% tt; q <- n1 %/% 2
    gid <- factor(rep(paste0("g", 1:q), length.out = n1), levels = paste0("g", 1:q))
    G <- crossprod(matrix(rnorm(tt * tt), tt)) / tt + diag(0.5, tt)
    Rm <- crossprod(matrix(rnorm(tt * tt), tt)) / tt + diag(0.5, tt)
    U <- matrix(rnorm(q * tt), q, tt) %*% chol(G)
    Y <- U[as.integer(gid), ] + matrix(rnorm(n1 * tt), n1, tt) %*% chol(Rm) + matrix(seq_len(tt), n1, tt, byrow = TRUE)
    colnames(Y) <- paste0("y", seq_len(tt))
    dw <- data.frame(gid = gid, Y)
    dl <- data.frame(y = c(Y), gid = rep(gid, tt), trait = factor(rep(colnames(Y), each = n1), levels = colnames(Y)),
                     unit = rep(seq_len(n1), tt))
    list(d = dw, dl = dl, t = tt)
  } else if (case_ == "ar1ar1") {
    nr <- round(sqrt(n / 1.25)); nc <- n %/% nr; ntot <- nr * nc; q <- ntot %/% 4
    g <- expand.grid(col = seq_len(nc), row = seq_len(nr))
    Kr <- 0.6 ^ abs(outer(1:nr, 1:nr, "-")); Kc <- 0.4 ^ abs(outer(1:nc, 1:nc, "-"))
    e <- as.numeric(kronecker(t(chol(Kr)), t(chol(Kc))) %*% rnorm(ntot)) * 1.0
    gid <- factor(sample(rep(paste0("g", 1:q), length.out = ntot)), levels = paste0("g", 1:q))
    y <- 4 + rnorm(q, 0, 1.0)[as.integer(gid)] + e + rnorm(ntot, 0, 0.5)
    d <- data.frame(g, gid = gid, y = y); d$row <- factor(d$row); d$col <- factor(d$col)
    d <- d[order(d$row, d$col), ]
    list(d = d, dims = c(nr, nc))
  } else stop("unknown case ", case_)
}

# ---- one fit per software ----------------------------------------------------------
fit_remlax <- function(case_, D) {
  if (case_ == "iid")  f <- rx_reml(y ~ x, random = ~ gid, data = D$d, backend = BK, verbose = FALSE, hessian = FALSE, blups = FALSE)
  if (case_ == "grm")  f <- rx_reml(y ~ x, random = ~ vm(gid, K = D$K), data = D$d, backend = BK, verbose = FALSE, hessian = FALSE, blups = FALSE)
  if (case_ %in% c("us3", "us6")) {
    dl <- D$dl
    f <- rx_reml(y ~ 0 + trait, random = ~ us(gid), residual = ~ us(trait):unit, data = dl, trait = "trait", unit = "unit",
                 backend = BK, verbose = FALSE, hessian = FALSE, blups = FALSE)
  }
  if (case_ == "ar1ar1") f <- rx_reml(y ~ 1, random = ~ gid, residual = ~ ar1(row):ar1(col), data = D$d, backend = BK, verbose = FALSE, hessian = FALSE, blups = FALSE)
  comp <- c(unlist(lapply(f$sigmas, function(S) S[lower.tri(S, diag = TRUE)])), f$sigma_res[lower.tri(f$sigma_res, diag = TRUE)], unlist(f$rho))
  list(logLik = f$logLik_asreml, comp = comp, n_par = f$n_par, solver_s = f$seconds, compile_s = f$compile_s, n_eval = f$n_eval,
       converged = isTRUE(as.logical(f$scipy_success)), version = versions[["remlax"]])
}
fit_asreml <- function(case_, D) {
  d <- D$d
  if (case_ == "iid") a <- asreml(y ~ x, random = ~ gid, data = d, maxit = MAXIT_REF)
  if (case_ == "grm") { K <- D$K; a <- asreml(y ~ x, random = ~ vm(gid, K), data = d, maxit = MAXIT_REF) }
  if (case_ %in% c("us3", "us6")) {
    fml <- stats::as.formula(paste("cbind(", paste(grep("^y", names(d), value = TRUE), collapse = ","), ") ~ trait"))
    a <- asreml(fml, random = ~ us(trait):gid, residual = ~ id(units):us(trait), data = d, maxit = MAXIT_REF)
  }
  if (case_ == "ar1ar1") a <- asreml(y ~ 1, random = ~ gid, residual = ~ ar1(row):ar1(col), data = d, maxit = MAXIT_REF)
  vc <- summary(a)$varcomp
  it <- tryCatch({ m <- a$monitor; if (is.matrix(m)) ncol(m) - 1L else if (is.matrix(a$trace)) nrow(a$trace) else NA_integer_ }, error = function(e) NA_integer_)
  list(logLik = a$loglik, comp = vc$component, n_par = nrow(vc), n_eval = if (length(it) == 1L) it else NA_integer_,
       converged = isTRUE(a$converge), version = versions[["asreml"]])
}
fit_lme4 <- function(case_, D) {
  if (case_ != "iid") return(NULL)
  m <- lmer(y ~ x + (1 | gid), D$d, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  vc <- as.data.frame(VarCorr(m))
  list(logLik = as.numeric(logLik(m)), comp = vc$vcov, n_par = nrow(vc), n_eval = m@optinfo$feval,
       converged = length(m@optinfo$conv$lme4) == 0, version = versions[["lme4"]])
}
fit_sommer <- function(case_, D) {
  d <- D$d
  if (nrow(d) > SOMMER_MAX_N) return(NULL)
  if (case_ == "iid") so <- mmer(y ~ x, random = ~ vsr(gid), rcov = ~ units, data = d, verbose = FALSE, nIters = MAXIT_REF)
  else if (case_ == "grm") so <- mmer(y ~ x, random = ~ vsr(gid, Gu = D$K), rcov = ~ units, data = d, verbose = FALSE, nIters = MAXIT_REF)
  else if (case_ %in% c("us3", "us6")) {
    tt <- D$t
    fml <- stats::as.formula(paste("cbind(", paste(grep("^y", names(d), value = TRUE), collapse = ","), ") ~ 1"))
    so <- mmer(fml, random = ~ vsr(gid, Gtc = unsm(tt)), rcov = ~ vsr(units, Gtc = unsm(tt)), data = d, verbose = FALSE, nIters = MAXIT_REF)
  } else return(NULL)
  comp <- unlist(lapply(so$sigma, function(S) S[lower.tri(S, diag = TRUE)]))
  list(logLik = so$monitor[1, ncol(so$monitor)], comp = comp, n_par = length(comp), n_eval = ncol(so$monitor),
       converged = isTRUE(so$convergence), version = versions[["sommer"]])
}
FITS <- list(remlax = fit_remlax, asreml = fit_asreml, lme4 = fit_lme4, sommer = fit_sommer)

# ---- loop -------------------------------------------------------------------------
for (n in SIZES) for (case_ in CASES) {
  D <- gen_case(case_, n); n_obs <- if (is.null(D$dl)) nrow(D$d) else nrow(D$dl)
  cat(sprintf("\n### %s n=%d (%d observations) ###\n", case_, n, n_obs))
  for (sw in c("remlax", REFS)) {
    r <- tryCatch(timed(suppressWarnings(FITS[[sw]](case_, D))), error = function(e) e)
    if (inherits(r, "error")) { add_row(date = DATE, tag = TAG, backend = if (sw == "remlax") BK else "cpu", case = case_, n = n_obs,
                                      software = sw, status = "error", error = conditionMessage(r)); next }
    if (is.null(r$v)) { add_row(date = DATE, tag = TAG, backend = "cpu", case = case_, n = n_obs, software = sw, status = "skipped",
                              error = "model not expressible or size above cap"); next }
    v <- r$v
    add_row(date = DATE, tag = TAG, backend = if (sw == "remlax") BK else "cpu", case = case_, n = n_obs, n_par = v$n_par,
          software = sw, version = v$version, status = "ok", wall_s = round(r$s, 3),
          solver_s = if (is.null(v$solver_s)) NA else round(v$solver_s, 3),
          compile_s = if (is.null(v$compile_s)) NA else round(v$compile_s, 3),
          n_eval = if (n_eval_ok(v$n_eval)) v$n_eval else NA,
          s_per_eval = if (n_eval_ok(v$n_eval)) round((if (sw == "remlax" && !is.null(v$solver_s)) v$solver_s else r$s) / v$n_eval, 4) else NA,
          converged = v$converged, logLik = sprintf("%.6f", v$logLik), components = comp_str(v$comp), error = "")
  }
}
cat("\nwritten:", OUT_FILE, "(", length(RES), "rows )\n")
