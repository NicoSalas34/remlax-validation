# ==============================================================================
# test_remlax_nlme.R -- correlation structures against nlme (gls and lme)
# ------------------------------------------------------------------------------
# nlme (Pinheiro & Bates) is the open reference for residual correlation
# structures: corAR1, corCAR1, corARMA, corCompSymm, corSymm, corExp, corGaus,
# corLin, corSpher, and varIdent for heteroscedasticity. Until now these
# structures had asreml (licensed) as their only reference, or none at all
# (ar2, arma, gau, cor). Here each one is compared with nlme on the SAME model
# and the SAME data.
#
# HOW A GROUPED corStruct IS WRITTEN IN remlax. nlme writes
# corAR1(form = ~ t | g) for "AR1 within each group g, no correlation between
# groups". In remlax the residual is R = Sigma_trait (x) C_unit. Declaring the
# GROUP as the "trait" (trait = "g") with the id(g) structure gives
# Sigma = s2 I_g, and C_unit carries the correlation between times or
# positions. residual = ~ id(g):ar1(tf) with trait = "g" is therefore exactly
# corAR1(form = ~ t | g).
#
# logLik CONSTANT. nlme returns the FULL restricted likelihood, with the term
# -(n-p)/2 log(2 pi), as lme4 does. It compares with fit$logLik of remlax, not
# with fit$logLik_asreml. The first check verifies this.
#
# TOLERANCES: those of the other suites (-2logL 1e-6 absolute, variances 1e-4
# relative, fixed effects 1e-6 relative), with nlme tightened by glsControl /
# lmeControl (tolerance 1e-10). At its default settings nlme stops a few 1e-6
# from its optimum, and the tolerance would then measure nlme, not remlax.
#
#   Rscript tests/R/test_remlax_nlme.R
# ==============================================================================
suppressPackageStartupMessages({
  library(here); library(Matrix); library(jsonlite); library(nlme)
})
source(here::here("R", "charger_remlax.R"))
BK <- Sys.getenv("RX_BACKEND", "cpu")
cat("nlme", as.character(packageVersion("nlme")), "| R", R.version$major, R.version$minor, "\n")

FAILURES <- character(0)
check <- function(name, cond, detail = "") {
  cat(sprintf("  %-56s %s %s\n", name, if (isTRUE(cond)) "OK " else "FAIL", detail))
  if (!isTRUE(cond)) FAILURES <<- c(FAILURES, name)
}
rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-8)
section <- function(title, expr) {
  cat("\n=== ", title, " ===\n", sep = "")
  tryCatch(expr, error = function(e) check(paste0(sub(" .*", "", title), " : execution"),
                                           FALSE, conditionMessage(e)))
}
CTL <- glsControl(tolerance = 1e-10, msTol = 1e-10, msMaxIter = 500, opt = "nlminb")
CTL_LME <- lmeControl(tolerance = 1e-10, msTol = 1e-10, msMaxIter = 500, opt = "nlminb",
                      niterEM = 50)
cor_par <- function(m) coef(m$modelStruct$corStruct, unconstrained = FALSE)
n2l <- function(m) -2 * as.numeric(logLik(m))
se_gls <- function(m) sqrt(diag(vcov(m)))

# Checks common to every section: -2logL, sigma2, beta, SE(beta).
common_checks <- function(tag, m, fit, s2_nlme = m$sigma^2) {
  check(paste0(tag, " : identical -2logL"), abs(-2 * fit$logLik - n2l(m)) < 1e-6,
        sprintf("%.9f vs %.9f", -2 * fit$logLik, n2l(m)))
  check(paste0(tag, " : residual variance"), rel(fit$sigma_res[1, 1], s2_nlme) < 1e-4,
        sprintf("%.6f vs %.6f", fit$sigma_res[1, 1], s2_nlme))
  check(paste0(tag, " : fixed effects"),
        max(rel(as.numeric(fit$beta), as.numeric(coef(m)))) < 1e-6,
        sprintf("max rel diff %.2e", max(rel(as.numeric(fit$beta), as.numeric(coef(m))))))
  check(paste0(tag, " : standard errors of the fixed effects"),
        max(rel(sqrt(diag(fit$vbeta)), se_gls(m))) < 1e-4,
        sprintf("max rel diff %.2e", max(rel(sqrt(diag(fit$vbeta)), se_gls(m)))))
}

# Grouped series: G groups of Tn times, one covariate, one fixed effect.
series <- function(G, Tn, seed, noise) {
  set.seed(seed)
  n <- G * Tn
  g <- factor(rep(seq_len(G), each = Tn)); t <- rep(seq_len(Tn), G); x <- rnorm(n)
  e <- unlist(lapply(seq_len(G), function(i) noise(Tn)))
  data.frame(y = 1 + 0.5 * x + e, x = x, g = g, t = t, tf = factor(t))
}

# ==============================================================================
section("1. corAR1: AR1 within each group, against gls", {
# ==============================================================================
  d <- series(12, 10, 1, function(Tn) as.numeric(arima.sim(list(ar = 0.6), Tn)))
  m <- gls(y ~ x, d, correlation = corAR1(form = ~ t | g), method = "REML", control = CTL)
  fit <- rx_reml(y ~ x, residual = ~ id(g):ar1(tf), data = d, trait = "g",
                 backend = BK, verbose = FALSE)
  n <- nrow(d); p <- 2
  check("logLik constant : nlme = remlax$logLik",
        abs(fit$logLik - as.numeric(logLik(m))) < 1e-6 &&
          abs(fit$logLik_asreml - fit$logLik - (n - p) / 2 * log(2 * pi)) < 1e-9,
        sprintf("nlme %.6f | remlax %.6f | logLik_asreml %.6f (+(n-p)/2 log 2pi = %.4f)",
                as.numeric(logLik(m)), fit$logLik, fit$logLik_asreml, (n - p) / 2 * log(2 * pi)))
  common_checks("corAR1", m, fit)
  check("corAR1 : rho", abs(fit$rho[["residual"]] - cor_par(m)) < 1e-4,
        sprintf("%.7f vs %.7f", fit$rho[["residual"]], cor_par(m)))
})

# ==============================================================================
section("2. corCAR1: continuous AR1 at irregular times, against gls", {
# ==============================================================================
  set.seed(2); G <- 10; Tn <- 12; n <- G * Tn
  g <- factor(rep(seq_len(G), each = Tn))
  # times differ from one group to the next (one unit = one position)
  pos <- as.numeric(unlist(lapply(seq_len(G), function(i) sort(runif(Tn, 0, 15)))))
  x <- rnorm(n)
  e <- unlist(lapply(seq_len(G), function(i) {
    D <- as.matrix(dist(pos[g == i])); as.numeric(t(chol(0.7 ^ D)) %*% rnorm(Tn)) }))
  d <- data.frame(y = 1 + 0.5 * x + e, x = x, g = g, pos = pos)
  m <- gls(y ~ x, d, correlation = corCAR1(form = ~ pos | g), method = "REML", control = CTL)
  fit <- rx_reml(y ~ x, residual = ~ id(g):exp(pos), data = d, trait = "g",
                 backend = BK, verbose = FALSE)
  common_checks("corCAR1", m, fit)
  check("corCAR1 : phi (correlation per time unit)",
        abs(fit$rho[["residual"]] - cor_par(m)) < 1e-4,
        sprintf("%.7f vs %.7f (true 0.7)", fit$rho[["residual"]], cor_par(m)))
})

# ==============================================================================
section("3. corARMA(p=2): AR2, against gls", {
# ==============================================================================
  d <- series(10, 12, 3, function(Tn) as.numeric(arima.sim(list(ar = c(0.5, 0.2)), Tn)))
  m <- gls(y ~ x, d, correlation = corARMA(form = ~ t | g, p = 2), method = "REML",
           control = CTL)
  fit <- rx_reml(y ~ x, residual = ~ id(g):ar2(tf), data = d, trait = "g",
                 backend = BK, verbose = FALSE)
  common_checks("ar2", m, fit)
  check("ar2 : phi1 and phi2", max(abs(fit$rho[["residual"]] - cor_par(m))) < 1e-4,
        sprintf("remlax %s | nlme %s", paste(round(fit$rho[["residual"]], 6), collapse = " "),
                paste(round(cor_par(m), 6), collapse = " ")))
})

# ==============================================================================
section("4. corARMA(p=1, q=1): ARMA(1,1), against gls", {
# ==============================================================================
  d <- series(10, 12, 4, function(Tn) as.numeric(arima.sim(list(ar = 0.5, ma = 0.3), Tn)))
  m <- gls(y ~ x, d, correlation = corARMA(form = ~ t | g, p = 1, q = 1), method = "REML",
           control = CTL)
  fit <- rx_reml(y ~ x, residual = ~ id(g):arma(tf), data = d, trait = "g",
                 backend = BK, verbose = FALSE)
  common_checks("arma", m, fit)
  # remlax orders its two parameters as (theta_MA, phi_AR), nlme as
  # (phi, theta). The two pairs are compared after reordering.
  r_rx <- fit$rho[["residual"]]; r_nl <- cor_par(m)
  check("arma : phi and theta (remlax order: theta, phi)",
        abs(r_rx[2] - r_nl[1]) < 1e-4 && abs(r_rx[1] - r_nl[2]) < 1e-4,
        sprintf("remlax theta %.6f phi %.6f | nlme phi %.6f theta %.6f",
                r_rx[1], r_rx[2], r_nl[1], r_nl[2]))
})

# ==============================================================================
section("5. corCompSymm: uniform correlation, against gls", {
# ==============================================================================
  d <- series(10, 12, 5, function(Tn) rnorm(1, 0, 0.8) + rnorm(Tn, 0, 1))
  m <- gls(y ~ x, d, correlation = corCompSymm(form = ~ 1 | g), method = "REML", control = CTL)
  fit <- rx_reml(y ~ x, residual = ~ id(g):cor(tf), data = d, trait = "g",
                 backend = BK, verbose = FALSE)
  common_checks("cor", m, fit)
  check("cor : uniform rho", abs(fit$rho[["residual"]] - cor_par(m)) < 1e-4,
        sprintf("%.7f vs %.7f", fit$rho[["residual"]], cor_par(m)))
})

# ==============================================================================
section("6. corExp: exponential decay in 1D, against gls", {
# ==============================================================================
  # WHY THIS CASE IS THE MOST IMPORTANT OF THE SUITE. On 2026-09-28 it showed
  # that the rho^d structures started at exactly rho = 0, a point where the
  # gradient is zero by construction. The optimiser declared convergence
  # without moving, and -2logL stayed 65 units away from the optimum found by
  # gls. The start is now rho = 0.5 at the nearest neighbour.
  set.seed(6); G <- 8; Tn <- 15; n <- G * Tn
  g <- factor(rep(seq_len(G), each = Tn))
  pos <- as.numeric(unlist(lapply(seq_len(G), function(i) sort(runif(Tn, 0, 30)))))
  x <- rnorm(n)
  e <- unlist(lapply(seq_len(G), function(i) {
    D <- as.matrix(dist(pos[g == i])); as.numeric(t(chol(exp(-D / 4))) %*% rnorm(Tn)) }))
  d <- data.frame(y = 1 + 0.5 * x + e, x = x, g = g, pos = pos)
  m <- gls(y ~ x, d, correlation = corExp(form = ~ pos | g), method = "REML", control = CTL)
  fit <- rx_reml(y ~ x, residual = ~ id(g):exp(pos), data = d, trait = "g",
                 backend = BK, verbose = FALSE)
  common_checks("exp", m, fit)
  check("exp : rho = exp(-1/nlme range)",
        abs(fit$rho[["residual"]] - exp(-1 / cor_par(m))) < 1e-4,
        sprintf("%.7f vs %.7f (range %.4f, true 4)", fit$rho[["residual"]],
                exp(-1 / cor_par(m)), cor_par(m)))
  check("exp : the fit moved from its starting point",
        fit$n_iter > 0 && fit$rho[["residual"]] > 0.05,
        sprintf("%d iterations, rho %.4f", fit$n_iter, fit$rho[["residual"]]))
})

# ==============================================================================
section("7. corGaus: Gaussian decay in 1D, against gls", {
# ==============================================================================
  # INTEGER positions with step >= 1. The Gaussian correlation is rho^(d^2)
  # with rho = exp(-1/range^2). With a step < 1 and a short range, the optimal
  # rho falls below 1e-6 and the tanh parametrisation has no slope left. This
  # is a limit of the asreml parametrisation that remlax adopts, not an
  # optimisation defect. The coordinates are rescaled, as with asreml.
  set.seed(7); G <- 8; Tn <- 15; n <- G * Tn
  g <- factor(rep(seq_len(G), each = Tn))
  pos <- as.numeric(unlist(lapply(seq_len(G), function(i) sort(sample(1:40, Tn)))))
  x <- rnorm(n)
  e <- unlist(lapply(seq_len(G), function(i) {
    D <- as.matrix(dist(pos[g == i]))
    as.numeric(t(chol(exp(-(D / 3) ^ 2) + diag(1e-8, Tn))) %*% rnorm(Tn)) }))
  d <- data.frame(y = 1 + 0.5 * x + e, x = x, g = g, pos = pos)
  m <- gls(y ~ x, d, correlation = corGaus(form = ~ pos | g), method = "REML", control = CTL)
  fit <- rx_reml(y ~ x, residual = ~ id(g):gau(pos), data = d, trait = "g",
                 backend = BK, verbose = FALSE)
  common_checks("gau", m, fit)
  check("gau : rho = exp(-1/range^2)",
        abs(fit$rho[["residual"]] - exp(-1 / cor_par(m) ^ 2)) < 1e-4,
        sprintf("%.7f vs %.7f (range %.4f, true 3)", fit$rho[["residual"]],
                exp(-1 / cor_par(m) ^ 2), cor_par(m)))
})

# ==============================================================================
section("8. corLin: truncated linear tent (lvr), against gls", {
# ==============================================================================
  set.seed(8); G <- 8; Tn <- 15; n <- G * Tn
  g <- factor(rep(seq_len(G), each = Tn))
  pos <- as.numeric(unlist(lapply(seq_len(G), function(i) sort(runif(Tn, 0, 30)))))
  x <- rnorm(n)
  e <- unlist(lapply(seq_len(G), function(i) {
    D <- as.matrix(dist(pos[g == i]))
    as.numeric(t(chol(pmax(1 - D / 6, 0) + diag(1e-8, Tn))) %*% rnorm(Tn)) }))
  d <- data.frame(y = 1 + 0.5 * x + e, x = x, g = g, pos = pos)
  # The truncated tent makes the likelihood NOT SMOOTH in the range (a kink at
  # each observed distance) and multimodal. gls finds 200.19 from value = 5 and
  # 188.62 from value = 2, remlax 193.93 from its own start. The test therefore
  # compares two things that do not depend on the mode reached. The first is
  # the likelihood of both programs AT THE SAME range (that of the best gls,
  # imposed on remlax through fixed_theta). The second is that the free remlax
  # optimum with restarts is no worse than the best gls.
  fits_nl <- lapply(c(2, 3, 4, 6, 8), function(v) tryCatch(
    gls(y ~ x, d, correlation = corLin(form = ~ pos | g, value = v), method = "REML",
        control = CTL), error = function(e) NULL))
  fits_nl <- Filter(Negate(is.null), fits_nl)
  m <- fits_nl[[which.min(sapply(fits_nl, n2l))]]
  fit <- rx_reml(y ~ x, residual = ~ id(g):lvr(pos), data = d, trait = "g",
                 backend = BK, verbose = FALSE, n_restarts = 4)
  fit_fixed <- rx_reml(y ~ x, residual = ~ id(g):lvr(pos), data = d, trait = "g",
                      backend = BK, verbose = FALSE,
                      theta_init = c(log(m$sigma), log(cor_par(m))), fixed_theta = 2)
  common_checks("lvr (gls range)", m, fit_fixed)
  check("lvr : range read = range imposed",
        rel(fit_fixed$rho[["residual!range"]], cor_par(m)) < 1e-6,
        sprintf("%.5f vs %.5f (true 6)", fit_fixed$rho[["residual!range"]], cor_par(m)))
  check("lvr : free remlax optimum no worse than gls",
        -2 * fit$logLik <= n2l(m) + 1e-6,
        sprintf("remlax %.6f (range %.4f, %d restarts) | best gls over %d starts %.6f (range %.4f)",
                -2 * fit$logLik, fit$rho[["residual!range"]], 4L, length(fits_nl), n2l(m), cor_par(m)))
})

# ==============================================================================
section("9. corSpher: spherical in 2D, against gls", {
# ==============================================================================
  set.seed(9); G <- 6; Tn <- 20; n <- G * Tn
  g <- factor(rep(seq_len(G), each = Tn))
  cx <- runif(n, 0, 20); cy <- runif(n, 0, 20)
  sph <- function(D, r) ifelse(D < r, 1 - 1.5 * D / r + 0.5 * (D / r) ^ 3, 0)
  e <- unlist(lapply(seq_len(G), function(i) {
    D <- as.matrix(dist(cbind(cx, cy)[g == i, ]))
    as.numeric(t(chol(sph(D, 8) + diag(1e-8, Tn))) %*% rnorm(Tn)) }))
  x <- rnorm(n)
  d <- data.frame(y = 1 + 0.5 * x + e, x = x, g = g, cx = cx, cy = cy)
  m <- gls(y ~ x, d, correlation = corSpher(form = ~ cx + cy | g, value = 6), method = "REML",
           control = CTL)
  fit <- rx_reml(y ~ x, residual = ~ id(g):sph(cx, cy), data = d, trait = "g",
                 backend = BK, verbose = FALSE)
  common_checks("sph", m, fit)
  check("sph : range", rel(fit$rho[["residual!range"]], cor_par(m)) < 1e-3,
        sprintf("%.5f vs %.5f (true 8)", fit$rho[["residual!range"]], cor_par(m)))
})

# ==============================================================================
section("10. corSymm + varIdent: unstructured covariance (us), against gls", {
# ==============================================================================
  set.seed(10); G <- 40; Tn <- 4; n <- G * Tn
  S <- matrix(c(1.0, 0.5, 0.3, 0.1,
                0.5, 1.5, 0.6, 0.2,
                0.3, 0.6, 2.0, 0.7,
                0.1, 0.2, 0.7, 0.8), 4, 4)
  g <- factor(rep(seq_len(G), each = Tn)); t <- rep(seq_len(Tn), G); x <- rnorm(n)
  e <- as.numeric(t(matrix(rnorm(n), G, Tn) %*% chol(S)))
  d <- data.frame(y = 1 + 0.5 * x + e, x = x, g = g, t = t, tf = factor(t))
  m <- gls(y ~ x, d, correlation = corSymm(form = ~ t | g), weights = varIdent(form = ~ 1 | tf),
           method = "REML", control = CTL)
  fit <- rx_reml(y ~ x, residual = ~ us(tf):g, data = d, trait = "tf",
                 backend = BK, verbose = FALSE)
  check("us : identical -2logL", abs(-2 * fit$logLik - n2l(m)) < 1e-6,
        sprintf("%.9f vs %.9f", -2 * fit$logLik, n2l(m)))
  # Sigma of gls: sigma^2 * D C D, D = varIdent weights (1 for the reference level).
  w <- c(1, coef(m$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))
  w <- w[order(c(levels(d$tf)[1], names(w)[-1]))]
  C <- corMatrix(m$modelStruct$corStruct)[[1]]
  S_nlme <- m$sigma ^ 2 * outer(w, w) * C
  # ABSOLUTE difference relative to the largest variance. At -2logL equal to
  # 1e-8, a covariance of 0.1 differs by 2.7e-3 in relative terms between two
  # optimisers stopped on a flat surface. The surface causes this, not the
  # likelihood.
  diff_us <- max(abs(as.numeric(fit$sigma_res) - as.numeric(S_nlme))) / max(abs(S_nlme))
  check("us : the 10 components of Sigma", diff_us < 1e-3,
        sprintf("max abs diff / max var %.2e | diag remlax %s | nlme %s", diff_us,
                paste(round(diag(fit$sigma_res), 4), collapse = " "),
                paste(round(diag(S_nlme), 4), collapse = " ")))
  check("us : fixed effects", max(rel(as.numeric(fit$beta), as.numeric(coef(m)))) < 1e-6,
        sprintf("max rel diff %.2e", max(rel(as.numeric(fit$beta), as.numeric(coef(m))))))
})

# ==============================================================================
section("11. varIdent alone: one residual variance per group (diag), against gls", {
# ==============================================================================
  set.seed(11); G <- 5; Tn <- 30; n <- G * Tn
  sds <- c(0.5, 1, 1.5, 2, 3)
  g <- factor(rep(seq_len(G), each = Tn)); x <- rnorm(n)
  d <- data.frame(y = 1 + 0.5 * x + rnorm(n, 0, rep(sds, each = Tn)), x = x, g = g,
                  u = factor(seq_len(n)))
  m <- gls(y ~ x, d, weights = varIdent(form = ~ 1 | g), method = "REML", control = CTL)
  fit <- rx_reml(y ~ x, residual = ~ diag(g):u, data = d, trait = "g",
                 backend = BK, verbose = FALSE)
  check("diag : identical -2logL", abs(-2 * fit$logLik - n2l(m)) < 1e-6,
        sprintf("%.9f vs %.9f", -2 * fit$logLik, n2l(m)))
  w <- c(1, coef(m$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))
  v_nlme <- m$sigma ^ 2 * w ^ 2
  check("diag : the five variances", max(rel(diag(fit$sigma_res), v_nlme)) < 1e-4,
        sprintf("remlax %s | nlme %s", paste(round(diag(fit$sigma_res), 4), collapse = " "),
                paste(round(v_nlme, 4), collapse = " ")))
  check("diag : fixed effects", max(rel(as.numeric(fit$beta), as.numeric(coef(m)))) < 1e-6,
        sprintf("max rel diff %.2e", max(rel(as.numeric(fit$beta), as.numeric(coef(m))))))
})

# ==============================================================================
section("12. lme: random intercept + within-group corAR1, against lme", {
# ==============================================================================
  d <- series(15, 8, 12, function(Tn) rnorm(1, 0, 1.2) + as.numeric(arima.sim(list(ar = 0.5), Tn)))
  m <- lme(y ~ x, random = ~ 1 | g, correlation = corAR1(form = ~ t | g), data = d,
           method = "REML", control = CTL_LME)
  fit <- rx_reml(y ~ x, random = ~ g, residual = ~ id(g):ar1(tf), data = d, trait = "g",
                 backend = BK, verbose = FALSE)
  s2g <- as.numeric(VarCorr(m)[1, 1]); s2e <- m$sigma ^ 2
  check("lme AR1 : identical -2logL", abs(-2 * fit$logLik - n2l(m)) < 1e-6,
        sprintf("%.9f vs %.9f", -2 * fit$logLik, n2l(m)))
  check("lme AR1 : group variance", rel(fit$sigmas$g[1, 1], s2g) < 1e-4,
        sprintf("%.6f vs %.6f", fit$sigmas$g[1, 1], s2g))
  check("lme AR1 : residual variance", rel(fit$sigma_res[1, 1], s2e) < 1e-4,
        sprintf("%.6f vs %.6f", fit$sigma_res[1, 1], s2e))
  check("lme AR1 : rho", abs(fit$rho[["residual"]] - cor_par(m)) < 1e-4,
        sprintf("%.7f vs %.7f", fit$rho[["residual"]], cor_par(m)))
  check("lme AR1 : fixed effects", max(rel(as.numeric(fit$beta), as.numeric(fixef(m)))) < 1e-6,
        sprintf("max rel diff %.2e", max(rel(as.numeric(fit$beta), as.numeric(fixef(m))))))
  check("lme AR1 : standard errors of the fixed effects",
        max(rel(sqrt(diag(fit$vbeta)), sqrt(diag(vcov(m))))) < 1e-4,
        sprintf("max rel diff %.2e", max(rel(sqrt(diag(fit$vbeta)), sqrt(diag(vcov(m)))))))
  b_rx <- as.numeric(fit$blups$g); b_nl <- ranef(m)[levels(d$g), 1]
  check("lme AR1 : BLUPs", max(abs(b_rx - b_nl)) < 1e-4,
        sprintf("max diff %.2e", max(abs(b_rx - b_nl))))
})

# ==============================================================================
cat("\n==============================================================\n")
if (length(FAILURES)) {
  cat("FAILURES (", length(FAILURES), "):\n", sep = ""); for (e in FAILURES) cat("  -", e, "\n")
  quit(status = 1)
}
cat("ALL CHECKS PASS: remlax returns the same restricted likelihood\n",
    "as nlme::gls and nlme::lme on the correlation structures tested.\n")
