# ==============================================================================
# Tools of the simulation checks: dense REML arbiter, correlation structures
# written in R, extraction of the components, comparisons.
#
# Nothing here calls the remlax engine. The arbiter rebuilds V from the Sigma
# matrices and level parameters RETURNED by a fit, with the TRUE geometry of the
# design (coordinates, cells, incidences built by the case). It then evaluates
# REML in the V form:
#     -2 logL = (n - p) log 2pi + log|V| + log|X'V^-1 X| + y'Py
# This is the convention of fit$logLik (2pi constant included). fit$logLik_asreml
# omits it.
# ==============================================================================

LOG2PI <- log(2 * pi)

# --- dense REML ---------------------------------------------------------------
reml_dense <- function(V, y, X) {
  n <- length(y); p <- ncol(X)
  L <- chol(V)
  ldV <- 2 * sum(log(diag(L)))
  sol <- function(B) backsolve(L, forwardsolve(t(L), B))
  ViX <- sol(X); Viy <- sol(y)
  A <- crossprod(X, ViX)
  La <- chol(A)
  ldA <- 2 * sum(log(diag(La)))
  beta <- backsolve(La, forwardsolve(t(La), crossprod(X, Viy)))
  r <- Viy - ViX %*% beta
  neg2 <- (n - p) * LOG2PI + ldV + ldA + sum(y * r)
  list(neg2 = neg2, logLik = -0.5 * neg2)
}

# reml_y criterion of the structural model (rx_lag): REML of y* = A y, Jacobian
# -2 log|det A|, and mean term -log|X_y'X_y| + log|X'X| with X_y = A^-1 X.
reml_lag_dense <- function(V, y, X, M) {
  n <- length(y)
  A <- diag(n) - M
  ys <- as.numeric(A %*% y)
  base <- reml_dense(V, ys, X)$neg2
  ld <- function(B) as.numeric(determinant(B, logarithm = TRUE)$modulus)
  Xy <- solve(A, X)
  neg2 <- base - 2 * ld(A) - ld(crossprod(Xy)) + ld(crossprod(X))
  list(neg2 = neg2, logLik = -0.5 * neg2)
}

# --- correlations between levels, parameters on the natural scale -----------
C_ar1 <- function(q, rho) rho^abs(outer(seq_len(q), seq_len(q), "-"))
C_ar1ar1 <- function(nr, nc, rr, rc) kronecker(C_ar1(nr, rr), C_ar1(nc, rc))

.dist_axes <- function(co) {
  co <- as.matrix(co)
  dx <- abs(outer(co[, 1], co[, 1], "-"))
  dy <- if (ncol(co) > 1) abs(outer(co[, 2], co[, 2], "-")) else 0 * dx
  list(dx = dx, dy = dy, d = sqrt(dx^2 + dy^2))
}

matern_R <- function(h, nu) {
  out <- matrix(1, nrow(h), ncol(h))
  nz <- h > 0
  out[nz] <- 2^(1 - nu) / gamma(nu) * h[nz]^nu * besselK(h[nz], nu)
  out
}

# Metric families, with the parameter meaning of levels.py: phi is a
# correlation raised to a distance (exp, gau, iexp, igau, ieuc, aexp, agau),
# a range (sph, cir, lvr) or a scale (mtrn).
C_metric <- function(kind, co, par, nu = 0.5) {
  D <- .dist_axes(co)
  switch(kind,
    exp  = par[1]^D$d,
    gau  = par[1]^(D$d^2),
    lvr  = pmax(1 - D$dx / par[1], 0),
    iexp = par[1]^(D$dx + D$dy),
    igau = par[1]^(D$dx^2 + D$dy^2),
    ieuc = par[1]^D$d,
    aexp = par[1]^D$dx * par[2]^D$dy,
    agau = par[1]^(D$dx^2) * par[2]^(D$dy^2),
    sph  = { t <- pmin(D$d / par[1], 1); ifelse(D$d < par[1], 1 - 1.5 * t + 0.5 * t^3, 0) },
    cir  = { t <- pmin(D$d / par[1], 1)
             ifelse(D$d < par[1], 1 - (2 / pi) * (t * sqrt(1 - t^2) + asin(t)), 0) },
    mtrn = matern_R(D$d / par[1], nu),
    stop("unknown metric family: ", kind))
}

# Uniform correlation and general correlation
C_cor <- function(q, r) (1 - r) * diag(q) + r * matrix(1, q, q)
# corg: unit lower triangular L with normalised rows, theta in the strict
# triangle row by row (parametrisation of levels.py, rewritten here for the
# delta method).
C_corg_theta <- function(th, q) {
  L <- diag(q); k <- 1L
  for (i in seq_len(q)) for (j in seq_len(i - 1L)) { L[i, j] <- th[k]; k <- k + 1L }
  L <- L / sqrt(rowSums(L^2))
  tcrossprod(L)
}

# --- assembly of V -------------------------------------------------------------
# Contribution of a term: Z (Sigma (x) K) Z', Z = [Z_1 ... Z_t] (trait slow).
V_term <- function(Zl, Sigma, K = NULL) {
  Z <- as.matrix(do.call(cbind, lapply(Zl, as.matrix)))
  q <- ncol(as.matrix(Zl[[1]]))
  if (is.null(K)) K <- diag(q)
  G <- kronecker(as.matrix(Sigma), as.matrix(K))
  Z %*% G %*% t(Z)
}

# Residual of a section: R[i,j] = Sigma[tr_i, tr_j] * C[i,j], where C is the
# correlation between the UNITS of observations i and j (n_s x n_s matrix).
R_section <- function(Sigma, tr, C) {
  Sigma <- as.matrix(Sigma)
  if (length(Sigma) == 1L) return(Sigma[1, 1] * C)
  Sigma[tr, tr, drop = FALSE] * C
}
C_same_unit <- function(unit) (outer(unit, unit, "==")) * 1

# --- fits: normalised extraction ------------------------------------------------
# Returns a list `par`: sig (Sigma per term), res (Sigma per section, in section
# order), lvl (level parameters on the natural scale, named as in varcomp), Psi
# (lag), est (estimated components named as in varcomp), se (matching standard
# errors).
par_from_fit <- function(fit) {
  vc <- fit$varcomp
  est <- stats::setNames(vc$value, vc$name)
  se  <- stats::setNames(vc$se, vc$name)
  list(sig = fit$sigmas, res = unname(fit$sigmas_res), res_names = names(fit$sigmas_res),
       lvl = est[grepl("!", names(est), fixed = TRUE)],
       Psi = if (!is.null(fit$lag)) fit$lag$Psi else NULL,
       est = est, se = se, theta = fit$theta, hessian = fit$hessian,
       par_floor = fit$par_floor, par_ceil = fit$par_ceil)
}

# Component names in the manner of component_names() (components.py).
sigma_names <- function(name, S) {
  S <- as.matrix(S); t <- nrow(S)
  if (t == 1L) return(stats::setNames(S[1, 1], name))
  v <- c(); nm <- c()
  for (i in seq_len(t)) for (j in seq_len(i)) { v <- c(v, S[i, j]); nm <- c(nm, sprintf("%s[%d,%d]", name, i, j)) }
  stats::setNames(v, nm)
}

# Diagonal only (diag structure: the covariances are zero by construction)
diag_names <- function(name, S) {
  v <- sigma_names(name, S)
  if (nrow(as.matrix(S)) == 1L) return(v)
  v[grepl("\\[(\\d+),\\1\\]$", names(v))]
}

residual_names <- function(model) {
  secs <- model$residual$sections
  if (identical(model$residual$struct, "none")) return(character(0))
  if (is.null(secs) || length(secs) <= 1L) return("residual")
  vapply(secs, function(s) s$name %||% "section", "")
}

par_from_mme <- function(mm, model) {
  rn <- residual_names(model)
  res <- unname(mm$sigma_res)
  est <- c(unlist(lapply(names(mm$sigmas), function(nm) sigma_names(nm, mm$sigmas[[nm]]))),
           unlist(unname(Map(function(nm, S) sigma_names(nm, S), rn, res))))
  # level parameters (rho of an ar1 or ar1 x ar1 term), named as for the other
  # engines: <term>!ar1, <term>!ar1ar11, <term>!ar1ar12
  lvl <- c()
  for (nm in names(mm$rho %||% list())) {
    v <- mm$rho[[nm]]; if (!length(v)) next
    lv <- model$terms[[nm]]$level
    lvl <- c(lvl, stats::setNames(v, if (length(v) == 1L) paste0(nm, "!", lv) else paste0(nm, "!", lv, seq_along(v))))
  }
  est <- c(est, lvl)
  list(sig = mm$sigmas, res = res, res_names = rn, lvl = lvl, Psi = NULL,
       est = est, se = stats::setNames(rep(NA_real_, length(est)), names(est)))
}

# --- numerical delta method on theta (derived components) ---------------------
cov_theta <- function(par) {
  H <- par$hessian; th <- par$theta
  if (is.null(H)) return(NULL)
  free <- !(th <= par$par_floor + 1e-7 | th >= par$par_ceil - 1e-7)
  Hf <- (H[free, free, drop = FALSE] + t(H[free, free, drop = FALSE])) / 2
  Vi <- tryCatch(solve(Hf), error = function(e) NULL)
  if (is.null(Vi)) return(NULL)
  Vt <- matrix(0, length(th), length(th)); Vt[free, free] <- 2 * Vi
  Vt
}
delta_se <- function(f, par, h = 1e-6) {
  Vt <- cov_theta(par); th <- par$theta
  f0 <- f(th)
  if (is.null(Vt)) return(list(val = f0, se = rep(NA_real_, length(f0))))
  J <- sapply(seq_along(th), function(k) { e <- rep(0, length(th)); e[k] <- h
    (f(th + e) - f(th - e)) / (2 * h) })
  J <- matrix(J, nrow = length(f0))
  list(val = f0, se = sqrt(pmax(diag(J %*% Vt %*% t(J)), 0)))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- simulation -----------------------------------------------------------------
# y = X beta + L z, L L' = V_true: a single simulation route, the same form of V
# as the arbiter, with the true parameters.
simulate_y <- function(V, X, beta) {
  L <- t(chol(V))
  as.numeric(X %*% beta + L %*% stats::rnorm(nrow(V)))
}

# GRM from markers simulated in families (half sibs), with blending to make it
# positive definite (a centred GRM has rank q - 1).
simulated_grm <- function(q, m = 800, n_fam = 15, blend = 0.05, prefix = "G") {
  p <- stats::runif(m, 0.1, 0.9)
  fam <- sample(rep(seq_len(n_fam), length.out = q))
  mothers <- t(sapply(seq_len(n_fam), function(f) stats::rbinom(m, 2, p)))
  M <- t(sapply(seq_len(q), function(i) {
    h <- mothers[fam[i], ] / 2
    stats::rbinom(m, 1, h) + stats::rbinom(m, 1, p)
  }))
  rownames(M) <- sprintf("%s%03d", prefix, seq_len(q))
  rx_grm(M, ploidy = 2, blend = blend, coding = "count")
}

# A^-1 by Henderson's rules, without inbreeding (unknown parents = 0).
ainv_henderson <- function(sire, dam) {
  n <- length(sire)
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  add <- function(i, j, x) { ii <<- c(ii, i); jj <<- c(jj, j); xx <<- c(xx, x) }
  for (k in seq_len(n)) {
    s <- sire[k]; d <- dam[k]
    ps <- c(s, d)[c(s, d) > 0]
    delta <- switch(length(ps) + 1L, 1, 4 / 3, 2)
    add(k, k, delta)
    for (a in ps) { add(k, a, -delta / 2); add(a, k, -delta / 2) }
    for (a in ps) for (b in ps) add(a, b, delta / 4)
  }
  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))
}

# n x q indicator incidence of a factor (levels in the order of the factor)
inc <- function(f, n = length(f), rows = seq_along(f)) {
  f <- factor(f)
  Matrix::sparseMatrix(i = rows, j = as.integer(f), x = 1, dims = c(n, nlevels(f)),
                       dimnames = list(NULL, levels(f)))
}
# Places an m x q matrix (rows `rows` of n) in an n x q matrix
place <- function(M, rows, n) {
  M <- methods::as(M, "CsparseMatrix")
  T <- methods::as(M, "TsparseMatrix")
  Matrix::sparseMatrix(i = rows[T@i + 1L], j = T@j + 1L, x = T@x, dims = c(n, ncol(M)),
                       dimnames = list(NULL, colnames(M)))
}
