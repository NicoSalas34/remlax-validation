# ==============================================================================
# test_remlax_asreml.R -- validation against ASREML, the reference of the field
# ------------------------------------------------------------------------------
# The other suites compare remlax with lme4 and sommer. This one compares it
# with asreml, on the structures that asreml can fit and the two others cannot:
# AR1, separable AR1, 2D splines, and str() (covariance SHARED between several
# terms).
#
# LOG-LIKELIHOOD CONVENTION. asreml omits the constant (n-p)/2 * log(2*pi),
# remlax includes it (as lme4 does). The comparison therefore uses
# `logLik_asreml`, a field provided for this purpose. Without this precaution
# the difference is several hundred and looks like a model disagreement: it is
# 219.6 on a data set with n=240, p=1.
#
# TOLERANCES. logLik to 1e-5 (the two optimisers stop differently).
# Components to 1e-3 relative, EXCEPT those at the floor on both sides, where a
# relative difference has no meaning (comparing 3e-11 with 2e-10 gives 0.85).
#
#   Rscript scripts/tests/test_remlax_asreml.R
# ==============================================================================
suppressPackageStartupMessages({
  library(here); library(Matrix); library(jsonlite); library(splines); library(asreml)
})
source(here::here("R", "charger_remlax.R"))

FAILURES <- character(0)
check <- function(name, cond, detail = "") {
  cat(sprintf("  %-52s %s %s\n", name, if (isTRUE(cond)) "OK " else "FAIL", detail))
  if (!isTRUE(cond)) FAILURES <<- c(FAILURES, name)
}
relc <- function(a, b, floor_ = 1e-6) {
  if (max(abs(a), abs(b)) < floor_) return(0)   # two numerical zeros
  abs(a - b) / max(abs(b), floor_)
}
stab <- function(a, k = 5) { for (i in seq_len(k)) a <- update(a, trace = FALSE); a }

# ==============================================================================
cat("\n=== 1. one-dimensional AR1 ===\n")
# ==============================================================================
set.seed(1); q <- 40; rp <- 6; n <- q * rp
K <- 0.75 ^ abs(outer(1:q, 1:q, "-"))
u <- as.numeric(t(chol(K)) %*% rnorm(q)) * 1.4
lev <- rep(1:q, each = rp)
d <- data.frame(y = 2 + u[lev] + rnorm(n, 0, 0.8), lev = factor(lev, levels = 1:q))
d <- d[order(d$lev), ]
a <- stab(asreml(y ~ 1, random = ~ ar1v(lev), data = d, trace = FALSE,
                 maxit = 100, workspace = "2gb"))
vc <- summary(a)$varcomp
f <- rx_reml(y ~ 1, random = ~ ar1(lev), data = d, backend = "auto",
             n_restarts = 5, verbose = FALSE)
check("AR1 : rho",  relc(f$rho[[1]][1], vc["lev!lev!cor", "component"]) < 1e-3,
      sprintf("%.7f vs %.7f", f$rho[[1]][1], vc["lev!lev!cor", "component"]))
check("AR1 : variance", relc(f$sigmas[[1]][1, 1], vc["lev!lev!var", "component"]) < 1e-3,
      sprintf("%.7f vs %.7f", f$sigmas[[1]][1, 1], vc["lev!lev!var", "component"]))
check("AR1 : residual", relc(f$sigma_res[1, 1], vc["units!R", "component"]) < 1e-3)
check("AR1 : logLik (asreml convention)", abs(f$logLik_asreml - a$loglik) < 1e-5,
      sprintf("%.9f vs %.9f", f$logLik_asreml, a$loglik))

# ==============================================================================
cat("\n=== 2. separable AR1 x AR1 (spatial field) ===\n")
# ==============================================================================
set.seed(4); nr <- 12; nc <- 10; K2 <- 2
g <- expand.grid(col = 1:nc, row = 1:nr)
Kr <- 0.7 ^ abs(outer(1:nr, 1:nr, "-")); Kc <- 0.45 ^ abs(outer(1:nc, 1:nc, "-"))
u <- as.numeric(kronecker(t(chol(Kr)), t(chol(Kc))) %*% rnorm(nr * nc)) * 1.2
d <- do.call(rbind, replicate(K2, g, simplify = FALSE))
d$y <- 3 + u[(d$row - 1) * nc + d$col] + rnorm(nrow(d), 0, 0.9)
d$row <- factor(d$row); d$col <- factor(d$col); d <- d[order(d$row, d$col), ]
a <- stab(asreml(y ~ 1, random = ~ ar1(row):ar1(col), data = d, trace = FALSE,
                 maxit = 100, workspace = "2gb"))
vc <- summary(a)$varcomp
f <- rx_reml(y ~ 1, random = ~ ar1(row, col), data = d, backend = "auto",
             n_restarts = 5, verbose = FALSE)
check("AR1xAR1 : field variance",
      relc(f$sigmas[[1]][1, 1], vc["row:col", "component"]) < 1e-3,
      sprintf("%.7f vs %.7f", f$sigmas[[1]][1, 1], vc["row:col", "component"]))
check("AR1xAR1 : row rho",
      relc(f$rho[[1]][1], vc["row:col!row!cor", "component"]) < 1e-3,
      sprintf("%.7f vs %.7f", f$rho[[1]][1], vc["row:col!row!cor", "component"]))
check("AR1xAR1 : column rho",
      relc(f$rho[[1]][2], vc["row:col!col!cor", "component"]) < 1e-3,
      sprintf("%.7f vs %.7f", f$rho[[1]][2], vc["row:col!col!cor", "component"]))
check("AR1xAR1 : residual", relc(f$sigma_res[1, 1], vc["units!R", "component"]) < 1e-3)
check("AR1xAR1 : logLik", abs(f$logLik_asreml - a$loglik) < 1e-5,
      sprintf("%.9f vs %.9f", f$logLik_asreml, a$loglik))

# ==============================================================================
cat("\n=== 3. 2D spline (SAME basis matrices on both sides) ===\n")
# ==============================================================================
set.seed(6); nr <- 20; nc <- 18
g <- expand.grid(col = 1:nc, row = 1:nr)
surf <- 2 * sin(g$row / 4) + 1.5 * cos(g$col / 3) + 0.02 * g$row * g$col
y <- 5 + surf + rnorm(nrow(g), 0, 0.6)
sp <- rx_spl2d(g$row, g$col, nseg = c(6, 6))
Zs <- lapply(sp$terms, function(t) as.matrix(t$Zl[[1]]))
dd <- data.frame(y = y, sp$X, Zs[[1]], Zs[[2]], Zs[[3]])
colnames(dd) <- c("y", paste0("L", seq_len(ncol(sp$X))),
                  paste0("A", seq_len(ncol(Zs[[1]]))), paste0("B", seq_len(ncol(Zs[[2]]))),
                  paste0("C", seq_len(ncol(Zs[[3]]))))
a <- stab(asreml(stats::as.formula(paste("y ~", paste(paste0("L", seq_len(ncol(sp$X))), collapse = " + "))),
                 random = ~ grp(gA) + grp(gB) + grp(gC),
                 group = list(gA = grep("^A", names(dd)), gB = grep("^B", names(dd)),
                              gC = grep("^C", names(dd))),
                 data = dd, trace = FALSE, maxit = 100, workspace = "2gb"))
vc <- summary(a)$varcomp
f <- rx_fit(rx_model(y, cbind(1, sp$X), sp$terms), backend = "auto",
            n_restarts = 3, verbose = FALSE)
for (i in seq_along(sp$terms)) {
  nm <- sp$terms[[i]]$name; ref <- vc[c("grp(gA)", "grp(gB)", "grp(gC)")[i], "component"]
  check(sprintf("spline 2D : variance %s", nm),
        relc(f$sigmas[[nm]][1, 1], ref) < 1e-3,
        sprintf("%.6e vs %.6e", f$sigmas[[nm]][1, 1], ref))
}
check("spline 2D : residual", relc(f$sigma_res[1, 1], vc["units!R", "component"]) < 1e-3)
check("spline 2D : logLik", abs(f$logLik_asreml - a$loglik) < 1e-5,
      sprintf("%.9f vs %.9f", f$logLik_asreml, a$loglik))
check("spline 2D : zero component flagged as degenerate",
      "spl_xy" %in% f$degenerate_components,
      paste(f$degenerate_components, collapse = ", "))

# ==============================================================================
cat("\n=== 4. str(): covariance SHARED between two terms ===\n")
# ==============================================================================
set.seed(12); q <- 40; rp <- 6; n <- q * rp
S <- matrix(c(1.4, 0.7, 0.7, 0.9), 2, 2)
U <- matrix(rnorm(q * 2), q, 2) %*% chol(S)
gid <- factor(rep(paste0("g", 1:q), each = rp), levels = paste0("g", 1:q))
w <- rnorm(n)
d <- data.frame(y = 3 + U[as.integer(gid), 1] + w * U[as.integer(gid), 2] + rnorm(n, 0, 0.8),
                gid = gid, w = w)
Zw <- Matrix::sparseMatrix(i = seq_len(n), j = as.integer(gid), x = w, dims = c(n, q),
                           dimnames = list(NULL, levels(gid)))
a <- stab(asreml(y ~ 1 + w, random = ~ str(~ gid + gid:w, ~ us(2):id(40)),
                 data = d, trace = FALSE, maxit = 100, workspace = "2gb"))
vc <- summary(a)$varcomp
f <- rx_reml(y ~ 1 + w, random = ~ str(~ gid + mm(Zw, name = "slope"), struct = "us"),
             data = d, backend = "auto", n_restarts = 3, verbose = FALSE)
G <- f$sigmas[[1]]
check("str : var term 1", relc(G[1, 1], vc[1, "component"]) < 1e-3,
      sprintf("%.7f vs %.7f", G[1, 1], vc[1, "component"]))
check("str : covariance", relc(G[1, 2], vc[2, "component"]) < 1e-3,
      sprintf("%.7f vs %.7f", G[1, 2], vc[2, "component"]))
check("str : var term 2", relc(G[2, 2], vc[3, "component"]) < 1e-3,
      sprintf("%.7f vs %.7f", G[2, 2], vc[3, "component"]))
check("str : residual", relc(f$sigma_res[1, 1], vc["units!R", "component"]) < 1e-3)
check("str : logLik", abs(f$logLik_asreml - a$loglik) < 1e-5,
      sprintf("%.9f vs %.9f", f$logLik_asreml, a$loglik))

# ==============================================================================
cat("\n=== 5. vpredict and Wald ===\n")
# ==============================================================================
set.seed(31); q <- 60; rp <- 5; n <- q * rp
dv <- data.frame(gid = factor(rep(paste0("g", 1:q), each = rp)),
                 trt = factor(rep(1:3, length.out = n)))
dv$y <- 4 + rep(rnorm(q, 0, 1.3), each = rp) + rnorm(n, 0, 0.9) + as.numeric(dv$trt) * 0.5
a <- stab(asreml(y ~ trt, random = ~ gid, data = dv, trace = FALSE, maxit = 100), 4)
h_a <- asreml::vpredict(a, h2 ~ V1 / (V1 + V2))
w_a <- asreml::wald(a)
f <- rx_reml(y ~ trt, random = ~ gid, data = dv, backend = "auto", n_restarts = 3,
             verbose = FALSE, vpredict = c(h2 = "V1/(V1+V2)"), wald = TRUE)
pv <- f$vpredict$predictions
check("vpredict : estimate of h2",
      abs(pv$value[1] - h_a$Estimate) < 1e-5,
      sprintf("%.7f vs %.7f", pv$value[1], h_a$Estimate))
check("vpredict : STANDARD ERROR of h2 (delta method)",
      abs(pv$se[1] - h_a$SE) < 1e-5,
      sprintf("%.7f vs %.7f", pv$se[1], h_a$SE))
# asreml gives a SEQUENTIAL Wald test (type I) by default, remlax a CONDITIONAL
# one (type III). They agree on the LAST term, the only case where the two
# notions meet, so this term is compared.
tw <- f$wald$tests
i_trt <- which(tw$term == "trt")
p_asr <- w_a[rownames(w_a) == "trt", "Pr(Chisq)"]
check("Wald : df of the trt term", tw$df[i_trt] == 2L, sprintf("%d", tw$df[i_trt]))
check("Wald : p-value of the last term (seq == cond)",
      abs(tw$p[i_trt] - p_asr) < 1e-9,
      sprintf("%.4g vs %.4g", tw$p[i_trt], p_asr))

# ==============================================================================
cat("\n=== 6. SATURATED parametrisations: chol(t-1), ante(t-1), fa(t), rr(t) ===\n")
# ==============================================================================
# They all describe the SAME unstructured Sigma. A logLik divergence would
# signal a parametrisation error, not a model choice.
set.seed(77); qq <- 45; TT <- 4; rr_ <- 3
gg <- factor(rep(paste0("g", 1:qq), each = rr_)); nn <- qq * rr_
GG <- crossprod(matrix(rnorm(TT * TT), TT, TT)) + diag(TT)
UU <- matrix(rnorm(qq * TT), qq, TT) %*% chol(GG)
YY <- sapply(1:TT, function(k) 2 + k * 0.3 + rep(UU[, k], each = rr_) + rnorm(nn, 0, 0.9))
dsat <- data.frame(y = as.vector(YY), gid = rep(gg, TT),
                   trait = factor(rep(paste0("t", 1:TT), each = nn)), unit = rep(1:nn, TT))
fsat <- function(r) rx_reml(y ~ 0 + trait, random = stats::as.formula(paste("~", r)),
                            residual = ~ diag(trait):units, data = dsat, trait = "trait",
                            unit = "unit", backend = "auto", n_restarts = 4, verbose = FALSE)
fu <- fsat("us(gid)")
for (r in c("chol(gid, rank=3)", "ante(gid, rank=3)", "fa(gid, rank=4)", "rr(gid, rank=4)"))
  check(sprintf("%-18s == us", r), abs(fsat(r)$logLik - fu$logLik) < 1e-6,
        sprintf("diff %.2e", abs(fsat(r)$logLik - fu$logLik)))

cat("\n", strrep("=", 68), "\n", sep = "")
if (length(FAILURES)) {
  cat(sprintf("FAILURES (%d):\n%s\n", length(FAILURES), paste0("  - ", FAILURES, collapse = "\n")))
  quit(save = "no", status = 1L)
}
cat("remlax reproduces asreml on AR1, AR1xAR1, 2D splines, str(), vpredict and Wald.\n")
quit(save = "no", status = 0L)
