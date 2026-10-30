# ==============================================================================
# test_remlax_sommer.R -- kinship matrices and multi-trait models against sommer
# ------------------------------------------------------------------------------
# sommer (Covarrubias-Pazaran) is the open reference for genomic models: a
# kinship matrix on one term, several traits with a us or diag covariance
# between traits, several kinship matrices on the same factor (additive +
# dominance), multi-environment. test_remlax.R covered only the univariate case
# with one GRM and the bivariate us case without GRM.
#
# CONVERGENCE OF sommer. At its default settings (tolParConvLL = 1e-3), sommer
# stops within 1e-4 of the components. The measured differences would then be
# those of its stopping rule. It is tightened here (tolParConvLL = 1e-12, large
# nIters), after which the components agree to 1e-6 relative.
#
# logLik CONSTANT. The likelihood returned by sommer (monitor[1, ]) is neither
# that of remlax nor that of asreml. It depends on an internal scaling of the
# response. Only DIFFERENCES in -2logL between two models fitted by the same
# software are comparable, and these are compared (section 3, us against diag).
# The constant offset is measured and printed, never corrected by hand.
#
#   Rscript tests/R/test_remlax_sommer.R
# ==============================================================================
suppressPackageStartupMessages({
  library(here); library(Matrix); library(jsonlite); library(sommer)
})
source(here::here("R", "charger_remlax.R"))
BK <- Sys.getenv("RX_BACKEND", "cpu")
cat("sommer", as.character(packageVersion("sommer")), "| R", R.version$major, R.version$minor, "\n")

FAILURES <- character(0)
check <- function(name, cond, detail = "") {
  cat(sprintf("  %-56s %s %s\n", name, if (isTRUE(cond)) "OK " else "FAIL", detail))
  if (!isTRUE(cond)) FAILURES <<- c(FAILURES, name)
}
rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-8)
relmax <- function(A, B) max(abs(as.numeric(A) - as.numeric(B))) / max(abs(as.numeric(B)))
section <- function(title, expr) {
  cat("\n=== ", title, " ===\n", sep = "")
  tryCatch(expr, error = function(e) check(paste0(sub(" .*", "", title), " : execution"),
                                           FALSE, conditionMessage(e)))
}
MMER <- function(...) suppressWarnings(mmer(..., verbose = FALSE, tolParInv = 1e-8,
                                            tolParConvLL = 1e-12, nIters = 300))
llik_so <- function(so) so$monitor[1, ncol(so$monitor)]
grm <- function(q, nm, seed) {
  set.seed(seed)
  M <- matrix(rbinom(q * nm, 2, 0.3), q, nm)
  K <- A.mat(M - 1) + diag(1e-4, q)
  dimnames(K) <- list(paste0("g", 1:q), paste0("g", 1:q))
  list(K = K, M = M)
}
long <- function(Y, gid) {
  n1 <- nrow(Y)
  data.frame(y = c(Y), gid = rep(gid, ncol(Y)),
             trait = factor(rep(colnames(Y), each = n1), levels = colnames(Y)),
             unit = rep(seq_len(n1), ncol(Y)))
}

# ==============================================================================
section("1. one GRM, univariate: components, fixed effects, BLUP, PEV", {
# ==============================================================================
  q <- 50; rep_ <- 4; n <- q * rep_
  K <- grm(q, 300, 61)$K
  set.seed(62)
  gid <- factor(rep(paste0("g", 1:q), each = rep_), levels = paste0("g", 1:q))
  x <- rnorm(n)
  u <- as.numeric(t(chol(K)) %*% rnorm(q)) * 1.2
  y <- 5 + 0.4 * x + u[as.integer(gid)] + rnorm(n, 0, 1.0)
  d <- data.frame(y = y, x = x, gid = gid)
  so <- MMER(y ~ x, random = ~ vsr(gid, Gu = K), rcov = ~ units, data = d)
  fit <- rx_reml(y ~ x, random = ~ vm(gid, K = K), data = d, backend = BK, verbose = FALSE, pev = TRUE)
  s_g <- so$sigma[[1]][1, 1]; s_e <- so$sigma[[2]][1, 1]
  check("GRM : variance components", rel(fit$sigmas$gid[1, 1], s_g) < 1e-4 && rel(fit$sigma_res[1, 1], s_e) < 1e-4,
        sprintf("s2g %.6f/%.6f | s2e %.6f/%.6f", fit$sigmas$gid[1, 1], s_g, fit$sigma_res[1, 1], s_e))
  check("GRM : fixed effects", max(rel(as.numeric(fit$beta), as.numeric(so$Beta$Estimate))) < 1e-6,
        sprintf("max rel diff %.2e", max(rel(as.numeric(fit$beta), as.numeric(so$Beta$Estimate)))))
  check("GRM : standard errors of the fixed effects",
        max(rel(sqrt(diag(fit$vbeta)), sqrt(diag(as.matrix(so$VarBeta))))) < 1e-4,
        sprintf("max rel diff %.2e", max(rel(sqrt(diag(fit$vbeta)), sqrt(diag(as.matrix(so$VarBeta)))))))
  b_so <- so$U[[1]][[1]][levels(gid)]
  check("GRM : BLUPs", max(abs(as.numeric(fit$blups$gid) - b_so)) < 1e-4,
        sprintf("max diff %.2e", max(abs(as.numeric(fit$blups$gid) - b_so))))
  pev_so <- diag(as.matrix(so$PevU[[1]][[1]]))[levels(gid)]
  check("GRM : PEV of the BLUPs", max(rel(as.numeric(fit$pev$gid), pev_so)) < 1e-4,
        sprintf("max rel diff %.2e | mean PEV %.5f", max(rel(as.numeric(fit$pev$gid), pev_so)), mean(pev_so)))
  # Constant: measured, not corrected. sommer returns a likelihood shifted by a
  # constant that depends on the data (scaling of the response).
  off <- llik_so(so) - fit$logLik_asreml
  check("GRM : sommer logLik constant (measured)", is.finite(off),
        sprintf("sommer %.4f | remlax logLik_asreml %.4f | diff %.4f ((n-p)/2 log var(y) = %.4f)",
                llik_so(so), fit$logLik_asreml, off, (n - 2) / 2 * log(var(y))))
})

# ==============================================================================
section("2. three traits, genetic us with GRM and residual us", {
# ==============================================================================
  q <- 60; rep_ <- 3; n1 <- q * rep_
  K <- grm(q, 400, 51)$K
  set.seed(52)
  G <- matrix(c(1.5, 0.8, 0.3, 0.8, 1.0, 0.4, 0.3, 0.4, 0.8), 3)
  Rm <- matrix(c(0.8, 0.2, 0.1, 0.2, 1.2, 0.3, 0.1, 0.3, 0.6), 3)
  U <- t(chol(K)) %*% matrix(rnorm(q * 3), q, 3) %*% chol(G)
  gid <- factor(rep(paste0("g", 1:q), each = rep_), levels = paste0("g", 1:q))
  Y <- U[as.integer(gid), ] + matrix(rnorm(n1 * 3), n1, 3) %*% chol(Rm) +
       matrix(c(3, 1, 2), n1, 3, byrow = TRUE)
  colnames(Y) <- c("y1", "y2", "y3")
  dw <- data.frame(gid = gid, Y)
  so <- MMER(cbind(y1, y2, y3) ~ 1, random = ~ vsr(gid, Gu = K, Gtc = unsm(3)),
             rcov = ~ vsr(units, Gtc = unsm(3)), data = dw)
  dl <- long(Y, gid)
  fit <- rx_reml(y ~ 0 + trait, random = ~ us(gid, K = K), residual = ~ us(trait):unit,
                 data = dl, trait = "trait", unit = "unit", backend = BK, verbose = FALSE)
  check("us3 + GRM : the 6 genetic components", relmax(fit$sigmas$gid, so$sigma[[1]]) < 1e-4,
        sprintf("max rel diff %.2e | diag remlax %s | sommer %s", relmax(fit$sigmas$gid, so$sigma[[1]]),
                paste(round(diag(fit$sigmas$gid), 4), collapse = " "), paste(round(diag(so$sigma[[1]]), 4), collapse = " ")))
  check("us3 + GRM : the 6 residual components", relmax(fit$sigma_res, so$sigma[[2]]) < 1e-4,
        sprintf("max rel diff %.2e", relmax(fit$sigma_res, so$sigma[[2]])))
  rg <- cov2cor(fit$sigmas$gid); rs <- cov2cor(so$sigma[[1]])
  check("us3 + GRM : the 3 genetic correlations", max(abs(rg[lower.tri(rg)] - rs[lower.tri(rs)])) < 1e-4,
        sprintf("remlax %s | sommer %s", paste(round(rg[lower.tri(rg)], 4), collapse = " "),
                paste(round(rs[lower.tri(rs)], 4), collapse = " ")))
  check("us3 + GRM : fixed effects (3 means)", max(rel(as.numeric(fit$beta), as.numeric(so$Beta$Estimate))) < 1e-6,
        sprintf("max rel diff %.2e", max(rel(as.numeric(fit$beta), as.numeric(so$Beta$Estimate)))))
  b_so <- sapply(c("y1", "y2", "y3"), function(tr) so$U[[1]][[tr]][levels(gid)])
  check("us3 + GRM : BLUPs of the 3 traits", max(abs(as.numeric(fit$blups$gid) - as.numeric(b_so))) < 1e-4,
        sprintf("max diff %.2e over %d values", max(abs(as.numeric(fit$blups$gid) - as.numeric(b_so))), length(b_so)))
  # The same design in diag (section 3): -2logL is kept for the difference.
  assign("LL_US", list(rx = -2 * fit$logLik, so = -2 * llik_so(so)), envir = globalenv())
  assign("DW", dw, envir = globalenv()); assign("DL", dl, envir = globalenv()); assign("K3", K, envir = globalenv())
})

# ==============================================================================
section("3. three traits, genetic diag with GRM and residual diag; LRT us vs diag", {
# ==============================================================================
  so <- MMER(cbind(y1, y2, y3) ~ 1, random = ~ vsr(gid, Gu = K3, Gtc = diag(3)),
             rcov = ~ vsr(units, Gtc = diag(3)), data = DW)
  fit <- rx_reml(y ~ 0 + trait, random = ~ diag(gid, K = K3), residual = ~ diag(trait):unit,
                 data = DL, trait = "trait", unit = "unit", backend = BK, verbose = FALSE)
  check("diag3 + GRM : the 3 genetic variances", max(rel(diag(fit$sigmas$gid), diag(so$sigma[[1]]))) < 1e-4,
        sprintf("remlax %s | sommer %s", paste(round(diag(fit$sigmas$gid), 4), collapse = " "),
                paste(round(diag(so$sigma[[1]]), 4), collapse = " ")))
  check("diag3 + GRM : the 3 residual variances", max(rel(diag(fit$sigma_res), diag(so$sigma[[2]]))) < 1e-4,
        sprintf("remlax %s | sommer %s", paste(round(diag(fit$sigma_res), 4), collapse = " "),
                paste(round(diag(so$sigma[[2]]), 4), collapse = " ")))
  d_rx <- -2 * fit$logLik - LL_US$rx; d_so <- -2 * llik_so(so) - LL_US$so
  check("LRT us vs diag : same difference in -2logL", abs(d_rx - d_so) < 1e-3,
        sprintf("remlax %.5f | sommer %.5f (6 df)", d_rx, d_so))
})

# ==============================================================================
section("4. additive + dominance: two kinship matrices on the same factor", {
# ==============================================================================
  q <- 80; rep_ <- 3; n <- q * rep_
  g0 <- grm(q, 500, 71); K <- g0$K
  D <- D.mat(g0$M - 1) + diag(1e-4, q); dimnames(D) <- dimnames(K)
  set.seed(72)
  gid <- factor(rep(paste0("g", 1:q), each = rep_), levels = paste0("g", 1:q))
  ua <- as.numeric(t(chol(K)) %*% rnorm(q)) * 1.0
  ud <- as.numeric(t(chol(D)) %*% rnorm(q)) * 0.7
  y <- 10 + ua[as.integer(gid)] + ud[as.integer(gid)] + rnorm(n, 0, 1.0)
  d <- data.frame(y = y, gid = gid, gid2 = gid)
  so <- MMER(y ~ 1, random = ~ vsr(gid, Gu = K) + vsr(gid2, Gu = D), rcov = ~ units, data = d)
  fit <- rx_reml(y ~ 1, random = ~ vm(gid, K = K) + vm(gid2, K = D), data = d, backend = BK, verbose = FALSE)
  v_rx <- c(fit$sigmas$gid[1, 1], fit$sigmas$gid2[1, 1], fit$sigma_res[1, 1])
  v_so <- c(so$sigma[[1]][1, 1], so$sigma[[2]][1, 1], so$sigma[[3]][1, 1])
  check("A + D : the three components", max(rel(v_rx, v_so)) < 1e-4,
        sprintf("remlax %s | sommer %s", paste(round(v_rx, 5), collapse = " "), paste(round(v_so, 5), collapse = " ")))
  b_rx <- c(as.numeric(fit$blups$gid), as.numeric(fit$blups$gid2))
  b_so <- c(so$U[[1]][[1]][levels(gid)], so$U[[2]][[1]][levels(gid)])
  check("A + D : additive and dominance BLUPs", max(abs(b_rx - b_so)) < 1e-4,
        sprintf("max diff %.2e", max(abs(b_rx - b_so))))
})

# ==============================================================================
section("5. multi-environment: one genetic variance per environment with GRM", {
# ==============================================================================
  q <- 50; ne <- 4
  K <- grm(q, 300, 81)$K
  set.seed(82)
  s_e <- c(0.6, 1.0, 1.4, 0.8)
  U <- t(chol(K)) %*% matrix(rnorm(q * ne), q, ne) %*% diag(sqrt(s_e))
  d <- expand.grid(gid = paste0("g", 1:q), env = paste0("E", 1:ne), rep = 1:2, stringsAsFactors = TRUE)
  d <- d[sample(nrow(d), 300), ]                       # incomplete design, 2 replicates
  d$gid <- factor(d$gid, levels = paste0("g", 1:q)); d$env <- factor(d$env)
  d$y <- 2 + as.numeric(d$env) + U[cbind(as.integer(d$gid), as.integer(d$env))] +
         rnorm(nrow(d), 0, c(0.8, 1.0, 1.2, 0.9)[as.integer(d$env)])
  d$unit <- seq_len(nrow(d))
  so <- MMER(y ~ env, random = ~ vsr(dsr(env), gid, Gu = K), rcov = ~ vsr(dsr(env), units), data = d)
  fit <- rx_reml(y ~ env, random = ~ diag(gid, K = K), residual = ~ diag(env):unit,
                 data = d, trait = "env", unit = "unit", backend = BK, verbose = FALSE)
  v_so_g <- sapply(1:ne, function(k) so$sigma[[k]][1, 1]); v_so_e <- sapply(1:ne, function(k) so$sigma[[ne + k]][1, 1])
  check("MET : the 4 genetic variances", relmax(diag(fit$sigmas$gid), v_so_g) < 1e-4,
        sprintf("remlax %s | sommer %s", paste(round(diag(fit$sigmas$gid), 4), collapse = " "),
                paste(round(v_so_g, 4), collapse = " ")))
  check("MET : the 4 residual variances", relmax(diag(fit$sigma_res), v_so_e) < 1e-4,
        sprintf("remlax %s | sommer %s", paste(round(diag(fit$sigma_res), 4), collapse = " "),
                paste(round(v_so_e, 4), collapse = " ")))
  check("MET : fixed effects", max(rel(as.numeric(fit$beta), as.numeric(so$Beta$Estimate))) < 1e-6,
        sprintf("max rel diff %.2e", max(rel(as.numeric(fit$beta), as.numeric(so$Beta$Estimate)))))
})

# ==============================================================================
cat("\n==============================================================\n")
if (length(FAILURES)) {
  cat("FAILURES (", length(FAILURES), "):\n", sep = ""); for (e in FAILURES) cat("  -", e, "\n")
  quit(status = 1)
}
cat("ALL CHECKS PASS: same components, fixed effects, BLUP and PEV\n",
    "as sommer on the genomic models tested.\n")
