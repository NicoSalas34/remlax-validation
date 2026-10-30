# ==============================================================================
# test_remlax_asreml2.R -- validation of the structures and inference ADDED on
# 01/09/2026: lvr, ilv, anisotropic mtrn, own, dsum, predict, Kenward-Roger.
#
# THREE KINDS OF CONTROL, and all three are needed:
#   asreml     for what it can do (lvr, mtrn, dsum, predict)
#   pbkrtest   for Kenward-Roger, whose reference implementation it is
#   plain R    an independent reimplementation of the likelihood profile,
#              which tests the WHOLE CHAIN (R formula -> package -> solver) and
#              not only the final number
#
# The formulas of lvr and of the mtrn anisotropy are NOT in the manual. They
# were identified by comparing logLik curves with asreml (see
# docs/note_remlax.md). These tests pin them down.
#
#   Rscript scripts/tests/test_remlax_asreml2.R
# ==============================================================================
suppressPackageStartupMessages({ library(Matrix) })
ROOT <- Sys.getenv("IGE_ROOT", Sys.getenv("IGE_RACINE", "."))
source(file.path(ROOT, "R", "charger_remlax.R"))
BACKEND <- Sys.getenv("RX_BACKEND", "auto")

.n_ok <- 0L; .n_ko <- 0L
ok <- function(lbl, cond, det = "") {
  if (isTRUE(cond)) { .n_ok <<- .n_ok + 1L; cat(sprintf("  %-52s OK   %s\n", lbl, det)) }
  else { .n_ko <<- .n_ko + 1L; cat(sprintf("  %-52s FAIL  %s\n", lbl, det)) }
}
a_asreml <- requireNamespace("asreml", quietly = TRUE) &&
  !inherits(try(suppressMessages(library(asreml)), silent = TRUE), "try-error")
if (a_asreml) asreml.options(trace = FALSE)
a_pbkr <- requireNamespace("pbkrtest", quietly = TRUE) && requireNamespace("lme4", quietly = TRUE)

# --- REML likelihood of a GIVEN correlation, sigma2 profiled ----------------
# Control fully independent of the solver. It tells whether the structure built
# by remlax is the one we believe it is.
ll_profile <- function(C, y, X) {
  n <- length(y); p <- qr(X)$rank
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12) return(NA_real_)
  Ci <- chol2inv(chol(C)); A <- t(X) %*% Ci %*% X
  b <- solve(A, t(X) %*% Ci %*% y); r <- y - X %*% b
  s2 <- as.numeric(t(r) %*% Ci %*% r / (n - p))
  as.numeric(-0.5 * (sum(log(ev)) + (n - p) * log(s2) + log(det(A)) + (n - p)))
}

set.seed(20260901)

cat("\n=== 1. lvr: truncated tent max(0, 1 - d/phi) ===\n")
q <- 14
d1 <- data.frame(pos = 1:q)
H <- abs(outer(d1$pos, d1$pos, "-"))
d1$y <- as.numeric(t(chol(pmax(0, 1 - H / 6) + diag(0.25, q))) %*% rnorm(q)) + 5
f1 <- rx_reml(y ~ 1, random = NULL, residual = ~ lvr(pos), data = d1,
              backend = BACKEND, verbose = FALSE, n_restarts = 4)
range_ <- f1$rho[["residual!range"]]
# FINE GRID, not `optimize`. The likelihood of a truncated tent is not
# unimodal. Its support changes each time the range crosses an integer, which
# creates as many local optima. A golden-section search stays trapped there,
# and this is what happens to asreml on this same data set.
.ll_lvr <- function(p) ll_profile(matrix(pmax(0, 1 - H / p), q, q), matrix(d1$y),
                                 matrix(1, q, 1))
gr_p <- seq(1.02, 60, by = 0.002)
gr_v <- vapply(gr_p, .ll_lvr, 0)
i0 <- which.max(replace(gr_v, !is.finite(gr_v), -Inf))
refine <- optimize(.ll_lvr, gr_p[c(max(i0 - 1, 1), min(i0 + 1, length(gr_p)))],
                maximum = TRUE, tol = 1e-11)
gr <- if (refine$objective > gr_v[i0]) refine else list(maximum = gr_p[i0], objective = gr_v[i0])
ok("lvr : logLik = independent R profile",
   abs(f1$logLik_asreml - gr$objective) < 1e-5,
   sprintf("remlax %.7f | R %.7f | range %.4f vs %.4f",
           f1$logLik_asreml, gr$objective, range_ %||% NA, gr$maximum))
if (a_asreml) {
  d1$posf <- factor(d1$pos)
  fa1 <- try(asreml(y ~ 1, residual = ~ lvr(posf), data = d1, maxit = 60), silent = TRUE)
  ok("lvr : asreml does not reach a better optimum",
     !inherits(fa1, "try-error") && fa1$loglik <= f1$logLik_asreml + 1e-5,
     if (inherits(fa1, "try-error")) "asreml failed" else
       sprintf("asreml %.6f <= remlax %.6f", fa1$loglik, f1$logLik_asreml))
}

cat("\n=== 2. mtrn: anisotropic Matern, parameters FIXED on both sides ===\n")
nr <- 7; nc <- 6
g2 <- expand.grid(x = 1:nr, y = 1:nc)
dx <- outer(g2$x, g2$x, "-"); dy <- outer(g2$y, g2$y, "-")
matern_ref <- function(phi, nu, delta, alpha, lambda) {
  u <-  dx * cos(alpha) + dy * sin(alpha)
  v <- -dx * sin(alpha) + dy * cos(alpha)
  sd_ <- sqrt(delta)
  h <- (abs(sd_ * u)^lambda + abs(v / sd_)^lambda)^(1 / lambda)
  z <- h / phi
  matrix(ifelse(z < 1e-12, 1, (2^(1 - nu) / gamma(nu)) * z^nu * besselK(pmax(z, 1e-12), nu)),
         nrow(g2), nrow(g2))
}
# Response drawn from a TRUE Matern field (range 2.5, nu 1.5) plus a nugget.
# On pure noise the range is not identifiable. The likelihood is flat and test 3
# used to pass without checking anything.
g2$resp <- as.numeric(t(chol(matern_ref(2.5, 1.5, 1, 0, 2) + diag(0.15, nrow(g2)))) %*%
                      rnorm(nrow(g2))) + 3

for (case_ in list(list(lbl = "isotropic nu=0.8", phi = 3, nu = 0.8, de = 1, al = 0, la = 2),
                 list(lbl = "nu=1.5 delta=2",  phi = 3, nu = 1.5, de = 2, al = 0, la = 2),
                 list(lbl = "delta=2 alpha=.6", phi = 3, nu = 0.8, de = 2, al = 0.6, la = 2),
                 list(lbl = "nu=0.5 lambda=1",  phi = 3, nu = 0.5, de = 1, al = 0, la = 1))) {
  # The formula is built as TEXT. substitute() would return a `call`, which the
  # parser refuses (rightly: a formula carries its environment).
  fml <- stats::as.formula(sprintf(
    "~ mtrn(x, y, phi = '%s F', nu = '%s F', delta = '%s F', alpha = '%s F', lambda = %s)",
    case_$phi, case_$nu, case_$de, case_$al, case_$la))
  fm <- rx_reml(resp ~ 1, residual = fml, data = g2, backend = BACKEND, verbose = FALSE)
  ref <- ll_profile(matern_ref(case_$phi, case_$nu, case_$de, case_$al, case_$la),
                   matrix(g2$resp), matrix(1, nrow(g2), 1))
  ok(sprintf("mtrn %s : logLik = besselK reference", case_$lbl),
     is.finite(ref) && abs(fm$logLik_asreml - ref) < 1e-6,
     sprintf("remlax %.7f | R %.7f", fm$logLik_asreml, ref))
  if (a_asreml) {
    fam <- try(asreml(resp ~ 1, residual = fml, data = g2, maxit = 30), silent = TRUE)
    if (!inherits(fam, "try-error") && is.finite(ref))
      ok(sprintf("mtrn %s : logLik = asreml", case_$lbl),
         abs(fam$loglik - fm$logLik_asreml) < 1e-5,
         sprintf("asreml %.7f | remlax %.7f", fam$loglik, fm$logLik_asreml))
  }
}

cat("\n=== 3. mtrn: the range is indeed ESTIMATED when freed ===\n")
fm2 <- rx_reml(resp ~ 1, residual = ~ mtrn(x, y, phi = 2, nu = "1.5 F"),
               data = g2, backend = BACKEND, verbose = FALSE, n_restarts = 4)
ph <- fm2$rho[["residual!phi"]]   # mtrn: prefixed keys (phi/nu/delta/alpha)
# Fine grid, for the same reason as for lvr. Nothing guarantees that the
# profile of a Matern is unimodal in phi, and `optimize` stopped on the edge of
# the interval while reporting an optimum.
.ll_m <- function(p) ll_profile(matern_ref(p, 1.5, 1, 0, 2), matrix(g2$resp),
                               matrix(1, nrow(g2), 1))
g_p <- exp(seq(log(0.05), log(200), length.out = 3000))
g_v <- vapply(g_p, .ll_m, 0)
j0 <- which.max(replace(g_v, !is.finite(g_v), -Inf))
refine <- optimize(.ll_m, g_p[c(max(j0 - 1, 1), min(j0 + 1, length(g_p)))],
                maximum = TRUE, tol = 1e-11)
prof <- if (refine$objective > g_v[j0]) refine else list(maximum = g_p[j0], objective = g_v[j0])
ok("mtrn : estimated phi = argmax of the R profile",
   is.finite(prof$objective) && abs(fm2$logLik_asreml - prof$objective) < 1e-5 &&
     abs(log(ph %||% NA) - log(prof$maximum)) < 0.02,
   sprintf("phi %.5f vs %.5f (true 2.5) | logLik %.7f vs %.7f",
           ph %||% NA, prof$maximum, fm2$logLik_asreml, prof$objective))

cat("\n=== 4. own: a user-written structure reproduces exp() ===\n")
q4 <- 16; d4 <- data.frame(pos = factor(1:q4), posn = 1:q4)
H4 <- abs(outer(1:q4, 1:q4, "-"))
d4$y <- as.numeric(t(chol(0.7^H4 + diag(0.3, q4))) %*% rnorm(q4)) + 2
fe <- rx_reml(y ~ 1, residual = ~ exp(posn), data = d4, backend = BACKEND,
              verbose = FALSE, n_restarts = 3)
fo <- rx_reml(y ~ 1, residual = ~ own(pos, expr = "exp(-lag*exp(p1))", n_par = 1),
              data = d4, backend = BACKEND, verbose = FALSE, n_restarts = 3)
ok("own : phi^d written as exp(-d*exp(p)) gives the same logLik",
   abs(fe$logLik - fo$logLik) < 1e-6,
   sprintf("exp() %.8f | own() %.8f", fe$logLik, fo$logLik))

cat("\n=== 5. dsum: one residual per section ===\n")
n5 <- 96
d5 <- data.frame(site = factor(rep(c("A", "B"), each = n5 / 2)),
                 gid = factor(rep(sprintf("g%02d", 1:24), each = 4)),
                 col = factor(rep(1:8, length.out = n5)))
u5 <- rnorm(24, 0, 1.1)
d5$y <- 4 + u5[as.integer(d5$gid)] +
  rnorm(n5, 0, ifelse(d5$site == "A", 0.5, 1.4))          # CONTRASTING variances
f5a <- rx_reml(y ~ 1, random = ~ iid(gid), residual = ~ units, data = d5,
               backend = BACKEND, verbose = FALSE)
f5b <- rx_reml(y ~ 1, random = ~ iid(gid), residual = ~ dsum(~ units | site),
               data = d5, backend = BACKEND, verbose = FALSE)
vr <- vapply(f5b$sigmas_res, function(S) S[1, 1], 0)
ok("dsum : two sections improve the likelihood",
   f5b$logLik > f5a$logLik + 1,
   sprintf("1 section %.4f -> 2 sections %.4f (%+.2f)", f5a$logLik, f5b$logLik,
           f5b$logLik - f5a$logLik))
ok("dsum : the two residual variances are separated",
   length(vr) == 2L && max(vr) / min(vr) > 3,
   sprintf("%s (true 0.25 and 1.96)", paste(sprintf("%.4f", vr), collapse = " / ")))
if (a_asreml) {
  fa5 <- try(asreml(y ~ 1, random = ~ gid, residual = ~ dsum(~ units | site),
                    data = d5, maxit = 60), silent = TRUE)
  if (!inherits(fa5, "try-error"))
    ok("dsum : logLik = asreml", abs(fa5$loglik - f5b$logLik_asreml) < 1e-4,
       sprintf("asreml %.6f | remlax %.6f", fa5$loglik, f5b$logLik_asreml))
}

cat("\n=== 5b. dsum: DIFFERENT structures by section ===\n")
d5b <- expand.grid(col = 1:8, site = c("A", "B", "C"))
d5b$col <- factor(d5b$col); d5b$site <- factor(d5b$site)
d5b$gid <- factor(rep(1:6, length.out = nrow(d5b)))
d5b$y <- rnorm(nrow(d5b)) + rnorm(6)[as.integer(d5b$gid)]
f5c <- rx_reml(y ~ 1, random = ~ iid(gid),
               residual = ~ dsum(~ ar1(col) + ar1(col) + units | site,
                                 levels = list("A", "B", "C")),
               data = d5b, backend = BACKEND, verbose = FALSE)
ok("dsum : three sections, two AR1 and one iid",
   length(f5c$sigmas_res) == 3L && sum(grepl("^[AB]$", names(f5c$rho))) == 2L &&
     f5c$n_par == 6L,
   sprintf("%d sections | %d parameters | rho on %s", length(f5c$sigmas_res),
           f5c$n_par, paste(setdiff(names(f5c$rho), ""), collapse = "+")))
# An AR1 section whose units are duplicated makes R singular. The solver must
# REFUSE it, not return NaN. The parity sweep unearthed this defect.
bad <- try(rx_reml(y ~ 1, random = ~ iid(gid),
                       residual = ~ dsum(~ ar1(col) + units | site,
                                         levels = list(c("A", "B"), "C")),
                       data = d5b, backend = BACKEND, verbose = FALSE), silent = TRUE)
ok("dsum : duplicated units in a section -> explicit refusal",
   inherits(bad, "try-error"), "error raised instead of -2logL = NaN")

cat("\n=== 6. predict: adjusted means and their standard errors ===\n")
n6 <- 90
d6 <- data.frame(trt = factor(rep(c("A", "B", "C"), each = n6 / 3)),
                 block = factor(rep(1:6, length.out = n6)),
                 x = rnorm(n6))
ub <- rnorm(6, 0, 0.9)
d6$y <- 5 + c(A = 0, B = 1.2, C = -0.4)[as.character(d6$trt)] + 0.5 * d6$x +
  ub[as.integer(d6$block)] + rnorm(n6, 0, 0.7)
f6 <- rx_reml(y ~ trt + x, random = ~ iid(block), data = d6, backend = BACKEND,
              verbose = FALSE, wald = TRUE)
pv <- rx_predict(f6, classify = "trt", sed = TRUE)
if (a_asreml) {
  fa6 <- try(asreml(y ~ trt + x, random = ~ block, data = d6, maxit = 60), silent = TRUE)
  if (!inherits(fa6, "try-error")) {
    pa <- predict(fa6, classify = "trt")$pvals
    ok("predict : values = asreml",
       max(abs(pv$predicted.value - pa$predicted.value)) < 1e-4,
       sprintf("max diff %.2e", max(abs(pv$predicted.value - pa$predicted.value))))
    ok("predict : standard errors = asreml",
       max(abs(pv$std.error - pa$std.error)) < 1e-4,
       sprintf("max diff %.2e", max(abs(pv$std.error - pa$std.error))))
  }
}
ok("predict : covariate at its mean, predictions estimable",
   all(pv$estimable) && nrow(pv) == 3L,
   sprintf("%s", paste(sprintf("%s=%.4f(%.4f)", pv$trt, pv$predicted.value,
                               pv$std.error), collapse = " ")))
ok("predict : standard error of differences available",
   is.finite(attr(pv, "sed.mean")),
   sprintf("mean sed %.5f", attr(pv, "sed.mean")))

cat("\n=== 7. Kenward-Roger vs pbkrtest ===\n")
if (a_pbkr) {
  set.seed(11); nb <- 5
  d7 <- do.call(rbind, lapply(1:nb, function(b) {
    k <- sample(3:7, 1)                                  # deliberate IMBALANCE
    data.frame(block = b, trt = sample(c("A", "B", "C"), k, replace = TRUE))
  }))
  d7$block <- factor(d7$block); d7$trt <- factor(d7$trt); d7$x <- rnorm(nrow(d7))
  u7 <- rnorm(nb, 0, 1.1)
  d7$y <- 2 + c(A = 0, B = 0.7, C = 1.3)[as.character(d7$trt)] + 0.4 * d7$x +
    u7[as.integer(d7$block)] + rnorm(nrow(d7), 0, 0.9)
  m1 <- lme4::lmer(y ~ trt + x + (1 | block), data = d7, REML = TRUE)
  m0 <- lme4::lmer(y ~ x + (1 | block), data = d7, REML = TRUE)
  kr <- pbkrtest::KRmodcomp(m1, m0)$stats
  se_ref <- sqrt(diag(as.matrix(pbkrtest::vcovAdj(m1))))
  f7 <- rx_reml(y ~ trt + x, random = ~ iid(block), data = d7, backend = BACKEND,
                verbose = FALSE, kenward_roger = TRUE)
  tw <- f7$kenward_roger$tests
  i <- which(tw$term == "trt")
  ok("K-R : denominator df = pbkrtest",
     abs(tw$denDF[i] - kr$ddf) < 1e-3,
     sprintf("remlax %.5f | pbkrtest %.5f", tw$denDF[i], kr$ddf))
  ok("K-R : F statistic = pbkrtest", abs(tw$F[i] - kr$Fstat) < 1e-5,
     sprintf("remlax %.6f | pbkrtest %.6f", tw$F[i], kr$Fstat))
  ok("K-R : adjusted standard errors = pbkrtest",
     max(abs(f7$kenward_roger$se_beta - se_ref)) < 1e-5,
     sprintf("max diff %.2e", max(abs(f7$kenward_roger$se_beta - se_ref))))
  ok("K-R : the adjustment INFLATES the standard errors",
     all(f7$kenward_roger$se_beta >= f7$kenward_roger$se_beta_raw - 1e-12),
     sprintf("mean ratio %.5f",
             mean(f7$kenward_roger$se_beta / f7$kenward_roger$se_beta_raw)))
} else cat("  (pbkrtest/lme4 missing: section skipped)\n")

cat("\n=== 8. predict with the random part (BLUP + prediction error) ===\n")
f8 <- rx_reml(y ~ 1 + x, random = ~ iid(block), data = d6, backend = BACKEND, verbose = FALSE)
pb <- rx_predict(f8, classify = "block", include_random = TRUE, verbose = FALSE)
pf <- rx_predict(f8, classify = "block", include_random = FALSE, verbose = FALSE)
ok("predict : the BLUP enters the mean by block",
   sd(pb$predicted.value) > 1e-6 && sd(pf$predicted.value) < 1e-9,
   sprintf("standard deviation with BLUP %.5f, without %.2e",
           sd(pb$predicted.value), sd(pf$predicted.value)))
ok("predict : the PREDICTION error stays finite and positive",
   all(is.finite(pb$std.error)) && all(pb$std.error > 0),
   sprintf("SE %s", paste(sprintf("%.4f", pb$std.error), collapse = " ")))

cat(sprintf("\n%s\n%d test(s) OK, %d failure(s)\n", strrep("=", 70), .n_ok, .n_ko))
if (.n_ko > 0L) quit(status = 1)
