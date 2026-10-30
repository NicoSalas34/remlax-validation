# ==============================================================================
# diag_asreml3_ecarts.R -- the three models where asreml and remlax do not
#                          return the same logLik (ma2, sph, cir): cause, not
#                          verdict
# ------------------------------------------------------------------------------
# Same data as test_remlax_asreml3.R (same seeds). For each model:
#   1. asreml: convergence message, bound codes, logLik trajectory, every
#      component; then 20 update() calls to see whether it oscillates;
#   2. an ARBITER written in dense R (log|V| + log|X'V^-1 X| + y'Py), evaluated
#      at the asreml parameters AND at the remlax ones, under each plausible
#      convention (sign of the MA theta; Euclidean or city-block distance).
#      The one that recovers the asreml value tells what asreml computes;
#   3. asreml refitted with its parameters FIXED at the remlax values (R.param,
#      con = "F"), to read the asreml logLik at the remlax point;
#   4. check of the constant: X (coding, columns, rank) on both sides, nedf of
#      asreml, n - p of remlax, log|X'V^-1 X| by hand.
#   Rscript tests/R/diag_asreml3_ecarts.R
# ==============================================================================
suppressPackageStartupMessages({ library(here); library(Matrix); library(jsonlite); library(asreml) })
source(here::here("R", "charger_remlax.R"))
BK <- Sys.getenv("RX_BACKEND", "cpu")
asreml.options(workspace = "2gb", pworkspace = "1gb")
cat("asreml", as.character(packageVersion("asreml")), "\n")

reml_R <- function(y, X, V) {
  # asreml convention: -1/2 (log|V| + log|X'V^-1X| + y'Py), without (n-p)/2 log 2pi
  cV <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(cV)) return(c(logLik = NA_real_, ldV = NA_real_, ldXVX = NA_real_, yPy = NA_real_))   # V not PD: the kernel is not valid on this metric
  Vi <- chol2inv(cV)
  XVX <- t(X) %*% Vi %*% X
  P <- Vi - Vi %*% X %*% solve(XVX, t(X) %*% Vi)
  ldV <- 2 * sum(log(diag(cV))); ldX <- determinant(XVX, log = TRUE)$modulus[1]; yPy <- as.numeric(t(y) %*% P %*% y)
  c(logLik = -0.5 * (ldV + ldX + yPy), ldV = ldV, ldXVX = ldX, yPy = yPy)
}
print_asreml <- function(a, tag) {
  cat("\n[", tag, "] asreml: converge =", a$converge, "| loglik =", sprintf("%.9f", a$loglik),
      "| nedf =", a$nedf, "| noeff =", paste(a$noeff, collapse = ","), "\n")
  vc <- summary(a)$varcomp; print(vc)
  cat("  vparameters :", paste(names(a$vparameters), sprintf("%.6f", a$vparameters), sep = "=", collapse = " | "), "\n")
  cat("  vparameters.con :", paste(a$vparameters.con, collapse = " "), "\n")
  if (!is.null(a$monitor)) { cat("  monitor (logLik per iteration):\n"); print(round(a$monitor[1:min(3, nrow(a$monitor)), , drop = FALSE], 6)) }
  cat("  fixed coef:", paste(round(as.numeric(coef(a)$fixed), 6), collapse = " "), "\n")
  invisible(vc)
}
follow_updates <- function(a, k = 20) {
  ll <- numeric(k)
  for (i in seq_len(k)) { a <- suppressWarnings(update(a, trace = FALSE)); ll[i] <- a$loglik }
  cat("  logLik after each update (", k, "):", paste(sprintf("%.6f", ll), collapse = " "), "\n")
  a
}
fix_R_param <- function(a, values) {
  # values: list of EXACT path ("/section/component") -> vector of values.
  # Each named component is fixed (con = "F") at these values. The paths are
  # those printed by str(a$R.param): "/g:tf/variance", "/g:tf/tf",
  # "/sph(cx, cy)/variance", "/sph(cx, cy)/sph(cx, cy)".
  rp <- a$R.param
  set_param <- function(x, path = "") {
    if (is.list(x) && !is.null(x$initial) && !is.null(x$con)) {
      if (path %in% names(values)) {
        v <- values[[path]]
        if (length(v) != length(x$initial)) stop("length ", length(v), " for ", path, " (expected ", length(x$initial), ")")
        x$initial[] <- v; x$con[] <- rep("F", length(x$con)); cat("    fixed", path, "=", paste(v, collapse = " "), "\n")
      }
      return(x)
    }
    if (is.list(x)) for (nm in names(x)) x[[nm]] <- set_param(x[[nm]], paste(path, nm, sep = "/"))
    x
  }
  set_param(rp)
}

# ------------------------------------------------------------------------------
# 1. ma2: grouped series (section A4 of the suite)
# ------------------------------------------------------------------------------
cat("\n=== ma2 ===\n")
make_series <- function(G, Tn, seed, noise) {
  set.seed(seed); n <- G * Tn
  g <- factor(rep(seq_len(G), each = Tn)); t <- rep(seq_len(Tn), G); x <- rnorm(n)
  e <- unlist(lapply(seq_len(G), function(i) noise(Tn)))
  data.frame(y = 1 + 0.5 * x + e, x = x, g = g, t = t, tf = factor(t))
}
d <- make_series(12L, 12L, 104, function(Tn) as.numeric(arima.sim(list(ma = c(0.5, 0.3)), Tn)))
d <- d[order(d$g, d$t), ]
fit <- rx_reml(y ~ x, residual = ~ id(g):ma2(tf), data = d, trait = "g", backend = BK, verbose = FALSE, n_restarts = 2)
cat("remlax : logLik_asreml", sprintf("%.9f", fit$logLik_asreml), "| s2", fit$sigma_res[1, 1], "| theta MA", paste(round(unlist(fit$rho), 6), collapse = " "),
    "| n_obs", fit$n_obs, "| p", ncol(fit$model$X), "| beta", paste(round(fit$beta, 6), collapse = " "), "| decrement", fit$newton_decrement, "\n")
a <- suppressWarnings(asreml(y ~ x, residual = ~ g:ma2(tf), data = d, trace = FALSE, maxit = 100))
vc <- print_asreml(a, "ma2, first call")
a <- follow_updates(a); vc <- print_asreml(a, "ma2, after 20 updates")
# R arbiter: MA(2) correlation under the two sign conventions
acf_ma2 <- function(th, Tn, sign_ = +1) {
  t1 <- sign_ * th[1]; t2 <- sign_ * th[2]; den <- 1 + t1^2 + t2^2
  r <- c(1, (t1 + t1 * t2) / den, t2 / den, rep(0, Tn - 3)); toeplitz(r)
}
X <- model.matrix(~ x, d); G <- nlevels(d$g); Tn <- nlevels(d$tf)
s2_as <- vc[grep("!R$", rownames(vc)), "component"]; th_as <- vc[grep("!theta|!cor|!lag|!ma", rownames(vc)), "component"]
cat("asreml read: s2", s2_as, " theta", paste(th_as, collapse = " "), "\n")
for (sign_ in c(+1, -1)) {
  V_as <- s2_as * kronecker(diag(G), acf_ma2(th_as, Tn, sign_))
  V_rx <- fit$sigma_res[1, 1] * kronecker(diag(G), acf_ma2(unlist(fit$rho)[1:2], Tn, sign_))
  r1 <- reml_R(d$y, X, V_as); r2 <- reml_R(d$y, X, V_rx)
  cat(sprintf("  R arbiter (sign %+d): at asreml params %.9f [asreml says %.9f] | at remlax params %.9f [remlax says %.9f]\n",
              sign_, r1["logLik"], a$loglik, r2["logLik"], fit$logLik_asreml))
  cat(sprintf("     pieces at asreml params: log|V| %.6f  log|X'V-1X| %.6f  y'Py %.6f\n", r1["ldV"], r1["ldXVX"], r1["yPy"]))
}
# asreml fixed at the remlax parameters
cat("  structure of R.param:\n"); str(a$R.param, max.level = 4, give.attr = FALSE)
rp <- tryCatch(fix_R_param(a, list("/g:tf/variance" = fit$sigma_res[1, 1], "/g:tf/tf" = as.numeric(unlist(fit$rho)[1:2]))),
               error = function(e) { cat("R.param :", conditionMessage(e), "\n"); NULL })
if (!is.null(rp)) {
  a2 <- tryCatch(suppressWarnings(asreml(y ~ x, residual = ~ g:ma2(tf), data = d, R.param = rp, trace = FALSE, maxit = 1)),
                 error = function(e) { cat("  asreml fixed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(a2)) { cat(sprintf("  asreml at the remlax params (fixed): %.9f  [remlax says %.9f; arbiter sign -1 says %.9f]\n", a2$loglik, fit$logLik_asreml,
                                  reml_R(d$y, X, fit$sigma_res[1, 1] * kronecker(diag(G), acf_ma2(unlist(fit$rho)[1:2], Tn, -1)))["logLik"])); print(summary(a2)$varcomp) }
}
# ------------------------------------------------------------------------------
# 2 and 3. sph and cir: irregular 2D field (sections B6, B7)
# ------------------------------------------------------------------------------
set.seed(200); nB <- 90
coordB <- data.frame(cx = round(runif(nB, 0, 20), 2), cy = round(runif(nB, 0, 20), 2))
coordB <- coordB[!duplicated(coordB), ]; nB <- nrow(coordB)
Dx <- abs(outer(coordB$cx, coordB$cx, "-")); Dy <- abs(outer(coordB$cy, coordB$cy, "-")); De <- sqrt(Dx^2 + Dy^2); Dcb <- Dx + Dy
sph_f <- function(D, r) { h <- pmin(D / r, 1); 1 - 1.5 * h + 0.5 * h^3 }
cir_f <- function(D, r) { h <- pmin(D / r, 1); 1 - (2 / pi) * (h * sqrt(1 - h^2) + asin(h)) }
kernels <- list(list(name = "sph", k = 6L, f = sph_f, C = sph_f(De, 8)),
                list(name = "cir", k = 7L, f = cir_f, C = cir_f(De, 8)))
for (s in kernels) {
  cat("\n===", s$name, "===\n")
  set.seed(300 + s$k)
  e <- as.numeric(t(chol(s$C + diag(1e-6, nB))) %*% rnorm(nB))
  d <- data.frame(y = 2 + e + rnorm(nB, 0, 0.1), coordB)
  fml <- stats::as.formula(sprintf("~ %s(cx, cy)", s$name))
  fit <- rx_reml(y ~ 1, residual = fml, data = d, backend = BK, verbose = FALSE, n_restarts = 3)
  rg_rx <- fit$rho[["residual!range"]]
  cat("remlax : logLik_asreml", sprintf("%.9f", fit$logLik_asreml), "| s2", fit$sigma_res[1, 1], "| range", rg_rx, "| n_iter", fit$n_iter, "| decrement", fit$newton_decrement, "\n")
  a <- suppressWarnings(asreml(y ~ 1, residual = fml, data = d, trace = FALSE, maxit = 100))
  vc <- print_asreml(a, paste(s$name, "first call"))
  a <- follow_updates(a); vc <- print_asreml(a, paste(s$name, "after 20 updates"))
  s2_as <- vc[grep("!R$", rownames(vc)), "component"][1]
  rg_as <- vc[grep("!phi|!range|!pow|!cor", rownames(vc)), "component"][1]
  cat("asreml read: s2", s2_as, " range", rg_as, "\n")
  X <- matrix(1, nB, 1)
  # unit order: remlax and the arbiter follow the row order of d (coordB). asreml sorts its
  # units by its own keys, but the likelihood does not depend on the order.
  for (dist in list(c("euclidean", "De"), c("city-block", "Dcb"))) {
    D <- get(dist[2])
    for (pars in list(c("asreml", s2_as, rg_as), c("remlax", fit$sigma_res[1, 1], rg_rx))) {
      V <- as.numeric(pars[2]) * s$f(D, as.numeric(pars[3]))
      r <- reml_R(d$y, X, V + diag(1e-10, nB))
      cat(sprintf("  R arbiter, distance %-11s at params %-6s : %.9f  (log|V| %.5f, log|X'V-1X| %.5f, y'Py %.5f)\n",
                  dist[1], pars[1], r["logLik"], r["ldV"], r["ldXVX"], r["yPy"]))
    }
  }
  cat(sprintf("  reminder: asreml says %.9f at its params; remlax says %.9f at its params\n", a$loglik, fit$logLik_asreml))
  # range on the asreml scale: what if its range were defined otherwise (phi = range/3, or sqrt)?
  for (fac in c(1, 3, 1 / 3, sqrt(3))) {
    V <- s2_as * s$f(De, rg_as * fac); r <- reml_R(d$y, X, V + diag(1e-10, nB))
    cat(sprintf("  Euclidean R arbiter, asreml range x %.3f : %.9f\n", fac, r["logLik"]))
  }
  # remlax fixed at the asreml params, and asreml fixed at the remlax params
  f2 <- rx_reml(y ~ 1, residual = fml, data = d, backend = BK, verbose = FALSE,
                theta_init = c(log(sqrt(s2_as)), log(rg_as)), maxiter = 0, polish = 0)
  cat(sprintf("  remlax at the asreml params: %.9f\n", f2$logLik_asreml))
  cat("  structure of R.param:\n"); str(a$R.param, max.level = 4, give.attr = FALSE)
  sec <- names(a$R.param)[1]
  rp <- tryCatch(fix_R_param(a, setNames(list(fit$sigma_res[1, 1], rg_rx), c(paste0("/", sec, "/variance"), paste0("/", sec, "/", sec)))),
                 error = function(e) { cat("R.param :", conditionMessage(e), "\n"); NULL })
  if (!is.null(rp)) {
    a2 <- tryCatch(suppressWarnings(asreml(y ~ 1, residual = fml, data = d, R.param = rp, trace = FALSE, maxit = 1)),
                   error = function(e) { cat("  asreml fixed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(a2)) { cat(sprintf("  asreml at the remlax params (fixed): %.9f\n", a2$loglik)); print(summary(a2)$varcomp) }
  }
}
cat("\nend of diagnostic\n")
