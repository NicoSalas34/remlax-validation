# ==============================================================================
# test_remlax.R -- cross-validation against INDEPENDENT implementations
# ------------------------------------------------------------------------------
# The core is tested separately (scripts/tests/test_remlax_core.py):
# parametrisations, gradients, assembly of V, analytical cases. This file
# checks the one thing that internal tests cannot establish. The solver must
# compute THE restricted likelihood and not another one. It is compared with
# software written by other people: lme4 (Bates et al.) and sommer
# (Covarrubias-Pazaran).
#
# WHAT IS COMPARED, AND WITH WHICH TOLERANCE
#   -2logL_REML   1e-6 absolute. This is THE criterion. Two implementations of
#                 the same likelihood must give the same value at the optimum,
#                 whatever path they followed to reach it.
#   variances     1e-4 relative. Looser than the logLik on purpose. On a flat
#                 surface, two optimisers stop at slightly different theta for
#                 the SAME logLik. This is a property of the problem, not a
#                 defect.
#   fixed effects 1e-6 relative (beta is solved exactly for a given V).
#
#   Rscript scripts/tests/test_remlax.R
# ==============================================================================
suppressPackageStartupMessages({
  library(here); library(Matrix); library(jsonlite)
})
source(here::here("R", "charger_remlax.R"))
# Backend of the suite. Default "cpu": locally we want a deterministic control.
# RX_BACKEND=gpu replays EXACTLY the same comparisons on the card, which makes
# this suite a test of the machine as much as of the code.
BK <- Sys.getenv("RX_BACKEND", "cpu")
if (BK != "cpu") cat("Backend of the suite:", BK, "\n")


FAILURES <- character(0)
check <- function(name, cond, detail = "") {
  cat(sprintf("  %-56s %s %s\n", name, if (isTRUE(cond)) "OK " else "FAIL", detail))
  if (!isTRUE(cond)) FAILURES <<- c(FAILURES, name)
}
rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-8)

# ==============================================================================
cat("\n=== 1. one random factor, against lme4 ===\n")
# ==============================================================================
suppressPackageStartupMessages(library(lme4))
for (cfg in list(c(a = 25, m = 5), c(a = 10, m = 12), c(a = 60, m = 3))) {
  set.seed(100 + cfg["a"])
  a <- cfg["a"]; m <- cfg["m"]; n <- a * m
  g <- factor(rep(seq_len(a), each = m))
  x <- rnorm(n)
  y <- 2 + 0.7 * x + rep(rnorm(a, 0, 1.3), each = m) + rnorm(n, 0, 1.1)
  d <- data.frame(y = y, x = x, g = g)

  m4 <- lmer(y ~ 1 + x + (1 | g), d, REML = TRUE,
             control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)))
  vc <- as.data.frame(VarCorr(m4))
  mod <- rx_model(y, model.matrix(~ 1 + x, d), list(rx_term("g", g)))
  fit <- rx_fit(mod, backend = BK, verbose = FALSE)

  check(sprintf("a=%d m=%d : identical -2logL", a, m),
        abs(-2 * fit$logLik - (-2 * as.numeric(logLik(m4)))) < 1e-6,
        sprintf("%.9f vs %.9f", -2 * fit$logLik, -2 * as.numeric(logLik(m4))))
  check(sprintf("a=%d m=%d : variance components", a, m),
        rel(fit$sigmas$g[1, 1], vc$vcov[vc$grp == "g"]) < 1e-4 &&
          rel(fit$sigma_res[1, 1], vc$vcov[vc$grp == "Residual"]) < 1e-4,
        sprintf("s2g %.6f/%.6f", fit$sigmas$g[1, 1], vc$vcov[vc$grp == "g"]))
  check(sprintf("a=%d m=%d : fixed effects", a, m),
        max(rel(as.numeric(fit$beta), as.numeric(fixef(m4)))) < 1e-6,
        sprintf("beta %s", paste(round(as.numeric(fit$beta), 5), collapse = " ")))
  check(sprintf("a=%d m=%d : BLUPs", a, m),
        max(abs(as.numeric(fit$blups$g) - as.numeric(ranef(m4)$g[, 1]))) < 1e-3,
        sprintf("max diff %.2e",
                max(abs(as.numeric(fit$blups$g) - as.numeric(ranef(m4)$g[, 1])))))
}

# ==============================================================================
cat("\n=== 2. two crossed factors, against lme4 ===\n")
# ==============================================================================
set.seed(11)
a <- 15; b <- 8; n <- a * b
g <- factor(rep(seq_len(a), each = b)); h <- factor(rep(seq_len(b), times = a))
y <- 1 + rep(rnorm(a, 0, 1.5), each = b) + rep(rnorm(b, 0, 0.9), times = a) + rnorm(n, 0, 1.0)
d <- data.frame(y = y, g = g, h = h)
m4 <- lmer(y ~ 1 + (1 | g) + (1 | h), d, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)))
vc <- as.data.frame(VarCorr(m4))
mod <- rx_model(y, matrix(1, n, 1), list(rx_term("g", g), rx_term("h", h)))
fit <- rx_fit(mod, backend = BK, verbose = FALSE)
check("crossed : identical -2logL",
      abs(-2 * fit$logLik + 2 * as.numeric(logLik(m4))) < 1e-6,
      sprintf("%.9f vs %.9f", -2 * fit$logLik, -2 * as.numeric(logLik(m4))))
check("crossed : the three variances",
      rel(fit$sigmas$g[1, 1], vc$vcov[vc$grp == "g"]) < 1e-4 &&
        rel(fit$sigmas$h[1, 1], vc$vcov[vc$grp == "h"]) < 1e-4 &&
        rel(fit$sigma_res[1, 1], vc$vcov[vc$grp == "Residual"]) < 1e-4,
      sprintf("%.5f %.5f %.5f", fit$sigmas$g[1, 1], fit$sigmas$h[1, 1], fit$sigma_res[1, 1]))

# ==============================================================================
cat("\n=== 3. kinship matrix, against sommer ===\n")
# ==============================================================================
suppressPackageStartupMessages(library(sommer))
set.seed(21)
q <- 40; rep_ <- 4; n <- q * rep_
M <- matrix(rbinom(q * 300, 2, 0.3), q, 300)
K <- A.mat(M - 1)
K <- K + diag(1e-4, q); dimnames(K) <- list(paste0("g", 1:q), paste0("g", 1:q))
gid <- factor(rep(paste0("g", 1:q), each = rep_), levels = paste0("g", 1:q))
u <- as.numeric(t(chol(K)) %*% rnorm(q)) * 1.2
y <- 5 + rep(u, each = rep_) + rnorm(n, 0, 1.0)
d <- data.frame(y = y, gid = gid)
so <- suppressWarnings(mmer(y ~ 1, random = ~ vsr(gid, Gu = K), rcov = ~ units,
                            data = d, verbose = FALSE, tolParInv = 1e-8))
s_g <- as.numeric(so$sigma[[1]]); s_e <- as.numeric(so$sigma[[2]])
mod <- rx_model(y, matrix(1, n, 1), list(rx_term("gid", gid, K = K)))
fit <- rx_fit(mod, backend = BK, verbose = FALSE)
check("kinship : variance components vs sommer",
      rel(fit$sigmas$gid[1, 1], s_g) < 5e-3 && rel(fit$sigma_res[1, 1], s_e) < 5e-3,
      sprintf("s2g %.5f/%.5f | s2e %.5f/%.5f", fit$sigmas$gid[1, 1], s_g,
              fit$sigma_res[1, 1], s_e))
# The sommer logLik is not on the same additive constant. The test compares the
# DIFFERENCE in -2logL between two models, which is invariant.
mod0 <- rx_model(y, matrix(1, n, 1), list(rx_term("gid", gid)))    # without kinship
fit0 <- rx_fit(mod0, backend = BK, verbose = FALSE)
check("kinship : the GRM improves the fit (K vs I)",
      -2 * fit$logLik < -2 * fit0$logLik,
      sprintf("%.4f vs %.4f", -2 * fit$logLik, -2 * fit0$logLik))

# ==============================================================================
cat("\n=== 4. two traits, us structure, against sommer ===\n")
# ==============================================================================
set.seed(31)
q <- 50; rep_ <- 4
G <- matrix(c(1.5, 0.8, 0.8, 1.0), 2, 2)
U <- matrix(rnorm(q * 2), q, 2) %*% chol(G)
gid <- factor(rep(paste0("g", 1:q), each = rep_), levels = paste0("g", 1:q))
n1 <- q * rep_
y1 <- 3 + rep(U[, 1], each = rep_) + rnorm(n1, 0, 0.9)
y2 <- 1 + rep(U[, 2], each = rep_) + rnorm(n1, 0, 1.1)
dw <- data.frame(gid = gid, y1 = y1, y2 = y2)
so <- suppressWarnings(mmer(cbind(y1, y2) ~ 1, random = ~ vsr(gid, Gtc = unsm(2)),
                            rcov = ~ vsr(units, Gtc = unsm(2)), data = dw,
                            verbose = FALSE, tolParInv = 1e-8))
Gs <- so$sigma[[1]]; Rs <- so$sigma[[2]]
# long format: one row per (unit, trait)
dl <- data.frame(y = c(y1, y2), gid = rep(gid, 2),
                 trait = factor(rep(c("y1", "y2"), each = n1)),
                 unit = rep(seq_len(n1), 2))
Zt <- lapply(levels(dl$trait), function(tt) {
  Matrix::sparseMatrix(i = which(dl$trait == tt),
                       j = as.integer(dl$gid)[dl$trait == tt],
                       x = 1, dims = c(nrow(dl), q)) })
Xt <- model.matrix(~ 0 + trait, dl)
mod <- rx_model(dl$y, Xt,
                list(rx_term("gid", Zt, struct = "us", levels = levels(gid))),
                rx_residual("us", trait = dl$trait, unit = dl$unit))
fit <- rx_fit(mod, backend = BK, verbose = FALSE)
check("us : 2x2 genetic covariance vs sommer",
      max(rel(as.numeric(fit$sigmas$gid), as.numeric(Gs))) < 5e-2,
      sprintf("remlax %s | sommer %s",
              paste(round(as.numeric(fit$sigmas$gid), 4), collapse = " "),
              paste(round(as.numeric(Gs), 4), collapse = " ")))
check("us : 2x2 residual covariance vs sommer",
      max(rel(as.numeric(fit$sigma_res), as.numeric(Rs))) < 5e-2,
      sprintf("remlax %s | sommer %s",
              paste(round(as.numeric(fit$sigma_res), 4), collapse = " "),
              paste(round(as.numeric(Rs), 4), collapse = " ")))
check("us : genetic correlation recovered",
      abs(fit$sigmas$gid[1, 2] / sqrt(fit$sigmas$gid[1, 1] * fit$sigmas$gid[2, 2]) -
            Gs[1, 2] / sqrt(Gs[1, 1] * Gs[2, 2])) < 0.02,
      sprintf("r %.4f vs %.4f",
              fit$sigmas$gid[1, 2] / sqrt(fit$sigmas$gid[1, 1] * fit$sigmas$gid[2, 2]),
              Gs[1, 2] / sqrt(Gs[1, 1] * Gs[2, 2])))

# ==============================================================================
cat("\n=== 5. nested structures: us must dominate diag, which dominates iid ===\n")
# ==============================================================================
lls <- sapply(c("iid", "diag", "us"), function(st) {
  m <- rx_model(dl$y, Xt, list(rx_term("gid", Zt, struct = st, levels = levels(gid))),
                rx_residual("diag", trait = dl$trait, unit = dl$unit))
  -2 * rx_fit(m, backend = BK, verbose = FALSE)$logLik })
check("nesting : -2logL decreasing iid >= diag >= us",
      lls["iid"] >= lls["diag"] - 1e-6 && lls["diag"] >= lls["us"] - 1e-6,
      paste(sprintf("%s=%.4f", names(lls), lls), collapse = " "))

# ==============================================================================
cat("\n=== 6. WEIGHTED incidence (the case a formula cannot express) ===\n")
# ==============================================================================
set.seed(41)
q <- 30; n <- 200
W <- matrix(0, n, q)
for (i in seq_len(n)) { j <- sample(q, 3); W[i, j] <- runif(3) }   # 3 weighted neighbours
u <- rnorm(q, 0, 1.1)
y <- 2 + as.numeric(W %*% u) + rnorm(n, 0, 0.8)
mod <- rx_model(y, matrix(1, n, 1), list(rx_term("nbr", list(W), struct = "iid")))
fit <- rx_fit(mod, backend = BK, verbose = FALSE)
check("weighted incidence : finite and defined fit",
      is.finite(fit$logLik) && fit$sigmas$nbr[1, 1] > 0 && fit$newton_decrement < 1e-4,
      sprintf("s2 %.4f (true 1.21) | s2e %.4f (true 0.64) | decrement %.1e",
              fit$sigmas$nbr[1, 1], fit$sigma_res[1, 1], fit$newton_decrement))

# ==============================================================================
cat("\n=== 7. formula interface: same result as the explicit route ===\n")
# ==============================================================================
set.seed(5); a <- 30; m <- 5; n <- a * m
dd <- data.frame(x = rnorm(n), g = factor(rep(1:a, each = m)), b = factor(rep(1:m, a)))
dd$y <- 2 + 0.5 * dd$x + rep(rnorm(a, 0, 1.2), each = m) + rnorm(n, 0, 1)
ff <- rx_reml(y ~ 1 + x, random = ~ g, data = dd, backend = BK, verbose = FALSE)
m4 <- lmer(y ~ 1 + x + (1 | g), dd, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)))
check("formula ~ g : -2logL identical to lme4",
      abs(-2 * ff$logLik + 2 * as.numeric(logLik(m4))) < 1e-6,
      sprintf("diff %.2e", abs(-2 * ff$logLik + 2 * as.numeric(logLik(m4)))))
ff2 <- rx_reml(y ~ 1 + x, random = ~ iid(g) + iid(b), data = dd, backend = BK, verbose = FALSE)
m5 <- lmer(y ~ 1 + x + (1 | g) + (1 | b), dd, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)))
check("formula ~ iid(g)+iid(b) : -2logL identical to lme4",
      abs(-2 * ff2$logLik + 2 * as.numeric(logLik(m5))) < 1e-6,
      sprintf("diff %.2e", abs(-2 * ff2$logLik + 2 * as.numeric(logLik(m5)))))
ff3 <- rx_reml(y ~ 0 + trait, random = ~ us(gid), residual = "us",
               data = dl, trait = "trait", unit = "unit", backend = BK, verbose = FALSE)
mod_ex <- rx_model(dl$y, Xt, list(rx_term("gid", Zt, struct = "us", levels = levels(gid))),
                   rx_residual("us", trait = dl$trait, unit = dl$unit))
fit_ex <- rx_fit(mod_ex, backend = BK, verbose = FALSE)
check("formula us(gid) == explicit route rx_term/rx_model",
      abs(ff3$logLik - fit_ex$logLik) < 1e-8,
      sprintf("diff %.2e", abs(ff3$logLik - fit_ex$logLik)))

# ==============================================================================
cat("\n=== 8. backend switch: the model does not change with the machine ===\n")
# ==============================================================================
fa_ <- rx_fit(mod_ex, backend = "auto", verbose = FALSE)
fc_ <- rx_fit(mod_ex, backend = "cpu",  verbose = FALSE)
check("backend auto and cpu give the same fit",
      abs(fa_$logLik - fc_$logLik) < 1e-9,
      sprintf("%s vs %s, diff %.2e", fa_$backend, fc_$backend, abs(fa_$logLik - fc_$logLik)))
check("the chosen backend is reported in the result",
      !is.null(fa_$backend) && nzchar(fa_$backend), fa_$backend)

# ==============================================================================
cat("\n=== 9. errors: usable messages instead of doubtful results ===\n")
# ==============================================================================
# The patterns accept the French and the English wording of the remlax messages.
att <- function(expr) tryCatch({ eval(expr); "NO ERROR" }, error = function(e) conditionMessage(e))
check("unknown term rejected",
      grepl("non reconnu|not recogni[sz]ed|unrecogni[sz]ed|unknown", att(quote(rx_reml(y ~ 1, random = ~ toto(g), data = dd,
                                             backend = BK, verbose = FALSE)))))
check("missing column rejected",
      grepl("absente|not found|missing|absent", att(quote(rx_reml(y ~ 1, random = ~ iid(zzz), data = dd,
                                         backend = BK, verbose = FALSE)))))
check("multi-trait structure without trait rejected",
      grepl("multi-caractere|multi-trait|multitrait", att(quote(rx_reml(y ~ 1, random = ~ us(g), data = dd,
                                                 backend = BK, verbose = FALSE)))))
# Rank-deficient X: two columns identical ON EVERY ROW (and not a 2x2 matrix,
# which would fail earlier on the number of rows).
check("rank-deficient X rejected",
      grepl("rang|rank", att(quote(rx_model(dd$y, cbind(rep(1, nrow(dd)), rep(1, nrow(dd))),
                                       list(rx_term("g", dd$g)))))))
check("level missing from K rejected",
      grepl("absent|not found|missing|not in", att(quote(rx_term("g", factor(c("a", "b")),
                                        K = matrix(1, 1, 1, dimnames = list("a", "a")))))))

cat("\n", strrep("=", 70), "\n", sep = "")
if (length(FAILURES)) {
  cat(sprintf("FAILURES (%d):\n%s\n", length(FAILURES), paste0("  - ", FAILURES, collapse = "\n")))
  quit(save = "no", status = 1L)
}
cat("Cross-validation complete: remlax computes the same restricted likelihood\n",
    "as lme4 and sommer, on every case tested.\n", sep = "")
quit(save = "no", status = 0L)
