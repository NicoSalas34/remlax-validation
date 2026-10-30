# ==============================================================================
# test_remlax_lme4.R -- crossed and nested terms, random slopes, BLUP and PEV
#                       against lme4
# ------------------------------------------------------------------------------
# test_remlax.R already covers one factor, two balanced crossed factors and the
# formula. This file covers the designs that lme4 fits and that the chapter
# uses under other names. They are a us covariance between two incidence
# columns of the same factor (random slopes, i.e. the two-column dge_ige term),
# the diag version of that term, a nesting, and three unbalanced crossed
# factors. On each design the BLUPs and the prediction error variance are
# checked as well.
#
# PEV: TWO CONVENTIONS. lme4 (ranef(condVar = TRUE)) returns Var(u | y) with
# beta KNOWN: sigma^2 (Z'Z + sigma^2 G^-1)^-1. remlax returns the PEV of the
# mixed model equations, diag(G - G Z' P Z G), which includes the uncertainty
# on beta, as asreml does. Both are computed here by hand from the lme4
# estimates, and each is compared with the software that returns it.
#
#   Rscript tests/R/test_remlax_lme4.R
# ==============================================================================
suppressPackageStartupMessages({
  library(here); library(Matrix); library(jsonlite); library(lme4)
})
source(here::here("R", "charger_remlax.R"))
BK <- Sys.getenv("RX_BACKEND", "cpu")
cat("lme4", as.character(packageVersion("lme4")), "| R", R.version$major, R.version$minor, "\n")

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
CTL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6),
                   check.conv.singular = "ignore")
n2l <- function(m) -2 * as.numeric(logLik(m))
inc <- function(f) Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
                                        dims = c(length(f), nlevels(f)))

# PEV by hand, from G (block diagonal by term), sigma2 and Z, X.
#   MME: sigma2 * [C^-1]_uu with C = [X'X X'Z ; Z'X Z'Z + sigma2 G^-1]
#   lme4: sigma2 * (Z'Z + sigma2 G^-1)^-1
pev_by_hand <- function(X, Z, Ginv, s2) {
  X <- as.matrix(X); Z <- as.matrix(Z)
  C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + s2 * Ginv))
  Ci <- solve(C); p <- ncol(X)
  list(mme = s2 * diag(Ci)[-(1:p)],
       cond = s2 * diag(solve(crossprod(Z) + s2 * Ginv)))
}
common_checks <- function(tag, m4, fit) {
  check(paste0(tag, " : identical -2logL"), abs(-2 * fit$logLik - n2l(m4)) < 1e-6,
        sprintf("%.9f vs %.9f", -2 * fit$logLik, n2l(m4)))
  check(paste0(tag, " : fixed effects"),
        max(rel(as.numeric(fit$beta), as.numeric(fixef(m4)))) < 1e-6,
        sprintf("max rel diff %.2e", max(rel(as.numeric(fit$beta), as.numeric(fixef(m4))))))
  check(paste0(tag, " : standard errors of the fixed effects"),
        max(rel(sqrt(diag(fit$vbeta)), sqrt(diag(as.matrix(vcov(m4)))))) < 1e-4,
        sprintf("max rel diff %.2e",
                max(rel(sqrt(diag(fit$vbeta)), sqrt(diag(as.matrix(vcov(m4))))))))
}

# ==============================================================================
section("1. correlated random slopes (1 + x | g): us with two columns", {
# ==============================================================================
  set.seed(41); q <- 30; mrep <- 8; n <- q * mrep
  g <- factor(rep(seq_len(q), each = mrep)); x <- rnorm(n)
  U <- matrix(rnorm(q * 2), q, 2) %*% chol(matrix(c(1, 0.3, 0.3, 0.5), 2))
  y <- 1 + 0.5 * x + U[g, 1] + U[g, 2] * x + rnorm(n, 0, 0.8)
  d <- data.frame(y = y, x = x, g = g)
  m4 <- lmer(y ~ x + (1 + x | g), d, REML = TRUE, control = CTL)
  G4 <- as.matrix(VarCorr(m4)$g)[1:2, 1:2]; s2 <- sigma(m4) ^ 2
  Z0 <- inc(g); Z1 <- Z0 * x
  mod <- rx_model(y, model.matrix(~ x, d),
                  list(rx_term("g", list(int = Z0, slope = Z1), struct = "us", levels = levels(g))))
  fit <- rx_fit(mod, backend = BK, verbose = FALSE, pev = TRUE)
  common_checks("us slopes", m4, fit)
  check("us slopes : the three components of Sigma",
        max(rel(as.numeric(fit$sigmas$g), as.numeric(G4))) < 1e-4,
        sprintf("remlax %s | lme4 %s", paste(round(as.numeric(fit$sigmas$g)[c(1, 2, 4)], 5), collapse = " "),
                paste(round(as.numeric(G4)[c(1, 2, 4)], 5), collapse = " ")))
  check("us slopes : residual variance", rel(fit$sigma_res[1, 1], s2) < 1e-4,
        sprintf("%.6f vs %.6f", fit$sigma_res[1, 1], s2))
  b4 <- as.matrix(ranef(m4)$g)[levels(g), ]
  check("us slopes : BLUPs (intercepts and slopes)",
        max(abs(as.numeric(fit$blups$g) - as.numeric(b4))) < 1e-4,
        sprintf("max diff %.2e over %d values", max(abs(as.numeric(fit$blups$g) - as.numeric(b4))), length(b4)))
  # PEV: G = Sigma (x) I_q, columns of Z ordered (column of Sigma, level)
  Z <- cbind(as.matrix(Z0), as.matrix(Z1))
  pv <- pev_by_hand(model.matrix(~ x, d), Z, kronecker(solve(fit$sigmas$g), diag(q)), fit$sigma_res[1, 1])
  pv4 <- attr(ranef(m4, condVar = TRUE)$g, "postVar")
  cond4 <- c(pv4[1, 1, ], pv4[2, 2, ])
  check("us slopes : remlax PEV = MME formula (with beta)",
        max(rel(as.numeric(fit$pev$g), pv$mme)) < 1e-4,
        sprintf("max rel diff %.2e | mean PEV %.5f", max(rel(as.numeric(fit$pev$g), pv$mme)), mean(pv$mme)))
  check("us slopes : lme4 condVar = formula without beta",
        max(rel(cond4, pv$cond)) < 1e-3,
        sprintf("max rel diff %.2e | mean condVar %.5f (the two conventions differ)",
                max(rel(cond4, pv$cond)), mean(cond4)))
})

# ==============================================================================
section("2. independent random slopes (1 + x || g): diag with two columns", {
# ==============================================================================
  set.seed(42); q <- 25; mrep <- 10; n <- q * mrep
  g <- factor(rep(seq_len(q), each = mrep)); x <- rnorm(n)
  y <- 2 + 0.3 * x + rnorm(q, 0, 1.1)[g] + rnorm(q, 0, 0.6)[g] * x + rnorm(n, 0, 0.9)
  d <- data.frame(y = y, x = x, g = g)
  m4 <- lmer(y ~ x + (1 + x || g), d, REML = TRUE, control = CTL)
  vc <- as.data.frame(VarCorr(m4))
  Z0 <- inc(g); Z1 <- Z0 * x
  mod <- rx_model(y, model.matrix(~ x, d),
                  list(rx_term("g", list(int = Z0, slope = Z1), struct = "diag", levels = levels(g))))
  fit <- rx_fit(mod, backend = BK, verbose = FALSE)
  common_checks("diag slopes", m4, fit)
  v4 <- c(vc$vcov[vc$grp == "g" & vc$var1 == "(Intercept)"], vc$vcov[vc$grp == "g.1" | (vc$grp == "g" & vc$var1 == "x")])
  v4 <- vc$vcov[vc$grp != "Residual"]
  check("diag slopes : the two variances", max(rel(diag(fit$sigmas$g), v4)) < 1e-4,
        sprintf("remlax %s | lme4 %s", paste(round(diag(fit$sigmas$g), 5), collapse = " "),
                paste(round(v4, 5), collapse = " ")))
  check("diag slopes : covariance zero by construction", fit$sigmas$g[1, 2] == 0,
        sprintf("%.1e", fit$sigmas$g[1, 2]))
})

# ==============================================================================
section("3. nesting (1 | g/h): two iid terms", {
# ==============================================================================
  set.seed(43); a <- 12; b <- 4; r <- 5; n <- a * b * r
  g <- factor(rep(seq_len(a), each = b * r)); h <- factor(rep(rep(seq_len(b), each = r), a))
  gh <- factor(paste(h, g, sep = ":"))     # h:g naming, as in ranef(lme4)
  y <- 3 + rnorm(a, 0, 1.2)[g] + rnorm(nlevels(gh), 0, 0.7)[gh] + rnorm(n, 0, 1.0)
  d <- data.frame(y = y, g = g, h = h, gh = gh)
  m4 <- lmer(y ~ 1 + (1 | g / h), d, REML = TRUE, control = CTL)
  vc <- as.data.frame(VarCorr(m4))
  fit <- rx_reml(y ~ 1, random = ~ g + gh, data = d, backend = BK, verbose = FALSE)
  common_checks("nested", m4, fit)
  check("nested : variances of g and of h within g",
        rel(fit$sigmas$g[1, 1], vc$vcov[vc$grp == "g"]) < 1e-4 &&
          rel(fit$sigmas$gh[1, 1], vc$vcov[vc$grp == "h:g"]) < 1e-4 &&
          rel(fit$sigma_res[1, 1], vc$vcov[vc$grp == "Residual"]) < 1e-4,
        sprintf("g %.5f/%.5f | h:g %.5f/%.5f | res %.5f/%.5f", fit$sigmas$g[1, 1],
                vc$vcov[vc$grp == "g"], fit$sigmas$gh[1, 1], vc$vcov[vc$grp == "h:g"],
                fit$sigma_res[1, 1], vc$vcov[vc$grp == "Residual"]))
  b_g <- ranef(m4)$g[levels(g), 1]; b_gh <- ranef(m4)$`h:g`[levels(gh), 1]
  check("nested : BLUPs of both terms",
        max(abs(as.numeric(fit$blups$g) - b_g)) < 1e-4 &&
          max(abs(as.numeric(fit$blups$gh) - b_gh)) < 1e-4,
        sprintf("max diff %.2e / %.2e", max(abs(as.numeric(fit$blups$g) - b_g)),
                max(abs(as.numeric(fit$blups$gh) - b_gh))))
})

# ==============================================================================
section("4. three crossed factors, unbalanced design", {
# ==============================================================================
  set.seed(44); a <- 20; b <- 6; cc <- 5
  full <- expand.grid(g = seq_len(a), bl = seq_len(b), an = seq_len(cc))
  full <- full[sample(nrow(full), 400), ]                     # 400 of the 600 cells
  g <- factor(full$g); bl <- factor(full$bl); an <- factor(full$an); n <- nrow(full)
  x <- rnorm(n)
  y <- 5 + 0.4 * x + rnorm(a, 0, 1.3)[g] + rnorm(b, 0, 0.8)[bl] + rnorm(cc, 0, 0.5)[an] + rnorm(n, 0, 1.0)
  d <- data.frame(y = y, x = x, g = g, bl = bl, an = an)
  m4 <- lmer(y ~ x + (1 | g) + (1 | bl) + (1 | an), d, REML = TRUE, control = CTL)
  vc <- as.data.frame(VarCorr(m4))
  fit <- rx_reml(y ~ x, random = ~ g + bl + an, data = d, backend = BK, verbose = FALSE, pev = TRUE)
  common_checks("3 crossed", m4, fit)
  v_rx <- c(fit$sigmas$g[1, 1], fit$sigmas$bl[1, 1], fit$sigmas$an[1, 1], fit$sigma_res[1, 1])
  v_4 <- vc$vcov[match(c("g", "bl", "an", "Residual"), vc$grp)]
  check("3 crossed : the four variances", max(rel(v_rx, v_4)) < 1e-4,
        sprintf("remlax %s | lme4 %s", paste(round(v_rx, 5), collapse = " "),
                paste(round(v_4, 5), collapse = " ")))
  diff_b <- max(sapply(c("g", "bl", "an"), function(f)
    max(abs(as.numeric(fit$blups[[f]]) - ranef(m4)[[f]][levels(d[[f]]), 1]))))
  check("3 crossed : BLUPs of the three terms", diff_b < 1e-4, sprintf("max diff %.2e", diff_b))
  # PEV of term g from the MME, with the three terms in Z
  Z <- cbind(as.matrix(inc(g)), as.matrix(inc(bl)), as.matrix(inc(an)))
  Ginv <- diag(c(rep(1 / v_rx[1], a), rep(1 / v_rx[2], b), rep(1 / v_rx[3], cc)))
  pv <- pev_by_hand(model.matrix(~ x, d), Z, Ginv, v_rx[4])
  pev_rx <- c(as.numeric(fit$pev$g), as.numeric(fit$pev$bl), as.numeric(fit$pev$an))
  check("3 crossed : PEV of the three terms = MME formula",
        max(rel(pev_rx, pv$mme)) < 1e-4,
        sprintf("max rel diff %.2e", max(rel(pev_rx, pv$mme))))
  check("3 crossed : PEV of g decreases with replication",
        cor(as.numeric(fit$pev$g), as.numeric(table(g))) < 0,
        sprintf("cor(PEV, n_rep) = %.3f", cor(as.numeric(fit$pev$g), as.numeric(table(g)))))
})

# ==============================================================================
section("5. one factor, very unequal group sizes, two covariates", {
# ==============================================================================
  set.seed(45); q <- 40
  size <- sample(1:12, q, replace = TRUE); n <- sum(size)
  g <- factor(rep(seq_len(q), size)); x1 <- rnorm(n); x2 <- runif(n)
  y <- 1 + 0.5 * x1 - 1.5 * x2 + rnorm(q, 0, 1.4)[g] + rnorm(n, 0, 1.0)
  d <- data.frame(y = y, x1 = x1, x2 = x2, g = g)
  m4 <- lmer(y ~ x1 + x2 + (1 | g), d, REML = TRUE, control = CTL)
  vc <- as.data.frame(VarCorr(m4))
  fit <- rx_reml(y ~ x1 + x2, random = ~ g, data = d, backend = BK, verbose = FALSE)
  common_checks("unbalanced", m4, fit)
  check("unbalanced : variance components",
        rel(fit$sigmas$g[1, 1], vc$vcov[vc$grp == "g"]) < 1e-4 &&
          rel(fit$sigma_res[1, 1], vc$vcov[vc$grp == "Residual"]) < 1e-4,
        sprintf("s2g %.6f/%.6f | s2e %.6f/%.6f", fit$sigmas$g[1, 1], vc$vcov[vc$grp == "g"],
                fit$sigma_res[1, 1], vc$vcov[vc$grp == "Residual"]))
  check("unbalanced : BLUPs",
        max(abs(as.numeric(fit$blups$g) - ranef(m4)$g[levels(g), 1])) < 1e-4,
        sprintf("max diff %.2e (sizes 1 to 12)", max(abs(as.numeric(fit$blups$g) - ranef(m4)$g[levels(g), 1]))))
})

# ==============================================================================
cat("\n==============================================================\n")
if (length(FAILURES)) {
  cat("FAILURES (", length(FAILURES), "):\n", sep = ""); for (e in FAILURES) cat("  -", e, "\n")
  quit(status = 1)
}
cat("ALL CHECKS PASS: same -2logL, components, fixed effects, BLUP\n",
    "and PEV (MME formula) as lme4 on every design tested.\n")
