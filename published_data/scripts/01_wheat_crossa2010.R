# ==============================================================================
# 01_wheat_crossa2010.R -- GBLUP, pedigree BLUP and multi-environment GBLUP on
#                          the CIMMYT wheat data of Crossa et al. (2010)
# ------------------------------------------------------------------------------
# DATA. BGLR::wheat: 599 inbred lines, 1279 DArT markers coded 0/1, grain yield
# in four mega-environments (columns "1", "2", "4", "5" of wheat.Y, already
# standardised within environment), and the pedigree coancestry matrix wheat.A
# (Crossa et al. 2010, Genetics 186:713-724, doi:10.1534/genetics.110.118521).
#
# RELATIONSHIP MATRIX. The lines are inbred, so a DArT presence is a
# homozygous dose: the markers enter VanRaden's formula as 0/2. G is built
# ONCE by vanraden() and the same matrix is passed to every program, so the
# comparison tests the solvers and not the builders of G. wheat.A is positive
# definite as shipped (smallest eigenvalue about 1e-3) and is used unchanged.
#
# MODELS
#   1. GBLUP in each environment            y = 1mu + g + e, g ~ N(0, s2g G)
#   2. pedigree BLUP in each environment    same with A
#   3. A + G in one model (environment 1)   two kernels on the same lines
#   4. multi-environment, long format       us(env) x G, diag residual by env
#   5. same with fa(env, 1) x G
#
# REFERENCES. asreml for every model (components, logLik without the 2 pi
# constant, BLUPs); sommer::mmes, rrBLUP::mixed.solve and gaston::lmm.aireml
# for the single-kernel GBLUP. Their convergence rules differ: sommer is run
# with a tight tolerance on the change of log-likelihood.
#
#   systemd-run --user --scope -q -p MemoryMax=7G \
#     Rscript published_data/scripts/01_wheat_crossa2010.R
# ==============================================================================
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(.self), "..", "R", "common.R"))
suppressPackageStartupMessages({
  library(asreml); library(sommer); library(rrBLUP); library(gaston)
})
asreml.options(workspace = "4gb", trace = FALSE, maxit = 100)

log_file <- file.path(pd_root(), "results", "wheat_crossa2010.log")
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
sink(log_file, split = TRUE)

start_checks("wheat_crossa2010", "Crossa et al. 2010, doi:10.1534/genetics.110.118521")
data(wheat, package = "BGLR")
ids <- rownames(wheat.A)
stopifnot(nrow(wheat.X) == length(ids), nrow(wheat.Y) == length(ids))
G <- vanraden(2 * wheat.X)
dimnames(G) <- list(ids, ids)
A <- wheat.A
cat(sprintf("G: %d x %d, mean diagonal %.3f; A: smallest eigenvalue %.2e\n",
            nrow(G), ncol(G), mean(diag(G)),
            min(eigen(A, only.values = TRUE, symmetric = TRUE)$values)))
envs <- colnames(wheat.Y)

# asreml random coefficients of one term, named by line, in the order `ids`.
asr_blup <- function(a, pattern, ids) {
  cr <- coef(a)$random
  rn <- rownames(cr)
  sel <- grepl(pattern, rn, fixed = TRUE)
  v <- cr[sel, 1]
  names(v) <- sub(".*_", "", rn[sel])
  v[ids]
}
asr_vc <- function(a) summary(a)$varcomp

timing <- list()

# ------------------------------------------------------------------------------
# 1 and 2. single-environment GBLUP and pedigree BLUP
# ------------------------------------------------------------------------------
for (kern in c("G", "A")) {
  K <- if (kern == "G") G else A
  Ka <- as_asreml_G(K)
  for (e in envs) {
    mname <- sprintf("%sBLUP env %s", if (kern == "G") "G" else "A", e)
    d <- data.frame(gid = factor(ids, levels = ids), y = wheat.Y[, e])
    t0 <- proc.time()[3]
    f <- rx_cpu(y ~ 1, random = ~ vm(gid, K = K), data = d)
    t_rx <- proc.time()[3] - t0
    t0 <- proc.time()[3]
    a <- stab(asreml(y ~ 1, random = ~ vm(gid, Ka), data = d))
    t_as <- proc.time()[3] - t0
    vc <- asr_vc(a)
    sg <- f$sigmas$gid[1, 1]; se <- f$sigma_res[1, 1]
    add_check(mname, "sigma2_g", "asreml", sg, vc[1, "component"])
    add_check(mname, "sigma2_e", "asreml", se, vc["units!R", "component"])
    add_check(mname, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
              tol = 1e-5, metric = "abs")
    add_blup_checks(mname, "asreml", f$blups$gid[ids, 1], asr_blup(a, "vm(gid, Ka)_", ids))
    timing[[mname]] <- c(remlax = t_rx, asreml = t_as)

    if (kern == "G") {
      # sommer, direct inversion with the raw relationship matrix
      s <- mmes(y ~ 1, random = ~ vsm(ism(gid), Gu = K), rcov = ~ units, data = d,
                nIters = 200, tolParConvLL = 1e-9, verbose = FALSE, dateWarning = FALSE)
      th <- unlist(lapply(s$theta, function(m) m[1, 1]))
      add_check(mname, "sigma2_g", "sommer", sg, th[1])
      add_check(mname, "sigma2_e", "sommer", se, th[length(th)])
      u_s <- s$uList[[1]]
      u_s <- setNames(as.numeric(u_s[, 1]), rownames(u_s))
      add_blup_checks(mname, "sommer", f$blups$gid[ids, 1], u_s[ids])

      # rrBLUP, REML by default
      r <- mixed.solve(d$y, K = K, method = "REML")
      add_check(mname, "sigma2_g", "rrBLUP", sg, r$Vu)
      add_check(mname, "sigma2_e", "rrBLUP", se, r$Ve)
      add_blup_checks(mname, "rrBLUP", f$blups$gid[ids, 1], r$u[ids])

      # gaston, AI-REML
      gz <- lmm.aireml(d$y, X = matrix(1, nrow(d), 1), K = K, verbose = FALSE,
                       eps = 1e-10, max_iter = 200L)
      add_check(mname, "sigma2_g", "gaston", sg, gz$tau)
      add_check(mname, "sigma2_e", "gaston", se, gz$sigma2)
      add_blup_checks(mname, "gaston", f$blups$gid[ids, 1],
                      setNames(as.numeric(gz$BLUP_omega), ids))
    }
  }
}

# ------------------------------------------------------------------------------
# 3. A + G in one model, environment 1
# ------------------------------------------------------------------------------
# A and G describe the same lines; the pedigree part captures what the
# markers miss (Crossa et al. 2010 fitted the same combination in a Bayesian
# setting). The two kernels are partly collinear, which makes this the
# hardest single-environment surface of the script.
d <- data.frame(gid = factor(ids, levels = ids), y = wheat.Y[, envs[1]])
d$ped <- d$gid
Ga <- as_asreml_G(G); Aa <- as_asreml_G(A)
mname <- sprintf("A+G env %s", envs[1])
t0 <- proc.time()[3]
f <- rx_cpu(y ~ 1, random = ~ vm(gid, K = G) + vm(ped, K = A), data = d, n_restarts = 3)
t_rx <- proc.time()[3] - t0
t0 <- proc.time()[3]
a <- stab(asreml(y ~ 1, random = ~ vm(gid, Ga) + vm(ped, Aa), data = d))
t_as <- proc.time()[3] - t0
vc <- asr_vc(a)
add_check(mname, "sigma2_G", "asreml", f$sigmas$gid[1, 1], vc["vm(gid, Ga)", "component"])
add_check(mname, "sigma2_A", "asreml", f$sigmas$ped[1, 1], vc["vm(ped, Aa)", "component"])
add_check(mname, "sigma2_e", "asreml", f$sigma_res[1, 1], vc["units!R", "component"])
add_check(mname, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
          tol = 1e-5, metric = "abs")
add_blup_checks(paste(mname, "G part"), "asreml", f$blups$gid[ids, 1],
                asr_blup(a, "vm(gid, Ga)_", ids))
add_blup_checks(paste(mname, "A part"), "asreml", f$blups$ped[ids, 1],
                asr_blup(a, "vm(ped, Aa)_", ids))
timing[[mname]] <- c(remlax = t_rx, asreml = t_as)

s <- mmes(y ~ 1, random = ~ vsm(ism(gid), Gu = G) + vsm(ism(ped), Gu = A), rcov = ~ units,
          data = d, nIters = 200, tolParConvLL = 1e-9, verbose = FALSE, dateWarning = FALSE)
th <- unlist(lapply(s$theta, function(m) m[1, 1]))
add_check(mname, "sigma2_G", "sommer", f$sigmas$gid[1, 1], th[1])
add_check(mname, "sigma2_A", "sommer", f$sigmas$ped[1, 1], th[2])
add_check(mname, "sigma2_e", "sommer", f$sigma_res[1, 1], th[3])

# ------------------------------------------------------------------------------
# 4 and 5. multi-environment GBLUP in long format
# ------------------------------------------------------------------------------
# One record per line and environment. The genetic covariance between
# environments is us (10 parameters) or fa(1) (4 loadings + 4 specific
# variances); the residual is a separate variance per environment.
dl <- data.frame(gid = factor(rep(ids, length(envs)), levels = ids),
                 env = factor(rep(paste0("E", envs), each = length(ids)),
                              levels = paste0("E", envs)),
                 y = as.vector(wheat.Y))
dl <- dl[order(dl$env, dl$gid), ]
env_lab <- levels(dl$env)

# Genetic covariance of an asreml multi-environment fit, rebuilt from varcomp.
# Labels: "env:vm(gid, Ga)!env_E2:E1" (us), "fa(env, 1):vm(gid, Ga)!E2!fa1" and
# "...!E2!var" (fa), "env_E2!R" (residual section).
us_label <- function(i, j) sprintf("env:vm(gid, Ga)!env_%s:%s", env_lab[i], env_lab[j])
asr_us <- function(vc, nt) {
  S <- matrix(NA_real_, nt, nt)
  for (i in seq_len(nt)) for (j in seq_len(i))
    S[i, j] <- S[j, i] <- vc[us_label(i, j), "component"]
  S
}
asr_fa <- function(vc, nt) {
  lam <- vc[sprintf("fa(env, 1):vm(gid, Ga)!%s!fa1", env_lab), "component"]
  psi <- vc[sprintf("fa(env, 1):vm(gid, Ga)!%s!var", env_lab), "component"]
  tcrossprod(lam) + diag(psi)
}
asr_res <- function(vc) vc[sprintf("env_%s!R", env_lab), "component"]
# BLUP matrix (lines x environments) of a multi-environment asreml fit.
# Coefficient names: "env_E1:vm(gid, Ga)_775" (us), "fa(env, 1)_E1:vm(gid, Ga)_775" (fa).
asr_blup_me <- function(a, prefix) {
  cr <- coef(a)$random
  rn <- rownames(cr)
  sapply(env_lab, function(e) {
    sel <- startsWith(rn, sprintf("%s%s:vm(gid, Ga)_", prefix, e))
    v <- cr[sel, 1]; names(v) <- sub(".*_", "", rn[sel]); v[ids]
  })
}

nt <- length(env_lab)
# Independent REML log-likelihood (asreml convention, without the 2 pi
# constant) of the multi-environment us model, by dense algebra in R.
X_me <- model.matrix(~ env, dl)
Zg_me <- model.matrix(~ 0 + gid, dl)
E_me <- model.matrix(~ 0 + env, dl)
dense_reml_us <- function(S, r) {
  V <- matrix(0, nrow(dl), nrow(dl))
  for (i in seq_len(nt)) for (j in seq_len(nt))
    V <- V + S[i, j] * (Zg_me * E_me[, i]) %*% G %*% t(Zg_me * E_me[, j])
  diag(V) <- diag(V) + as.vector(E_me %*% r)
  L <- chol(V)
  Vi <- chol2inv(L)
  XtViX <- crossprod(X_me, Vi %*% X_me)
  Py <- Vi %*% dl$y - Vi %*% X_me %*% solve(XtViX, crossprod(X_me, Vi %*% dl$y))
  -0.5 * (2 * sum(log(diag(L))) + as.numeric(determinant(XtViX)$modulus) + sum(dl$y * Py))
}
us_call <- quote(asreml(y ~ env, random = ~ us(env):vm(gid, Ga),
                        residual = ~ dsum(~ units | env), data = dl))
fa_call <- quote(asreml(y ~ env, random = ~ fa(env, 1):vm(gid, Ga),
                        residual = ~ dsum(~ units | env), data = dl))
for (st in c("us", "fa1")) {
  mname <- sprintf("multi-env %s x G", st)
  t0 <- proc.time()[3]
  f <- if (st == "us")
    rx_cpu(y ~ env, random = ~ us(gid, K = G), residual = ~ diag(env):units,
           data = dl, trait = "env", unit = "gid")
  else
    rx_cpu(y ~ env, random = ~ fa(gid, K = G, rank = 1), residual = ~ diag(env):units,
           data = dl, trait = "env", unit = "gid")
  t_rx <- proc.time()[3] - t0
  t0 <- proc.time()[3]
  a <- eval(if (st == "us") us_call else fa_call)
  a <- stab(a)
  t_as <- proc.time()[3] - t0
  cat(sprintf("%s: asreml from its default start: logLik %.6f, converged %s\n",
              mname, a$loglik, a$converge))
  if (st == "us" && !isTRUE(a$converge)) {
    # asreml keeps the us matrix positive definite by shrinking its AI updates
    # ("US updates modified ... to remain positive definite") and stops
    # unconverged at maxit. The remlax optimum sits on the boundary of the
    # cone (smallest eigenvalue printed below), which asreml cannot reach. The
    # counted checks are therefore:
    #   - the remlax log-likelihood against an independent dense evaluation of
    #     the REML formula at the remlax estimates (asreml convention);
    #   - asreml finds no higher likelihood than remlax, from its default
    #     start or started at the remlax estimates. asreml refuses a singular
    #     start ("Variance structure is not positive definite"), so the
    #     estimates are given with a 1e-7 ridge on the diagonal of Sigma.
    # asreml evaluated with every parameter fixed at the remlax estimates
    # (ridge 1e-9) and the components of its restarted run are context rows:
    # they describe where asreml goes, not a property of remlax.
    cat(sprintf("smallest eigenvalue of the remlax genetic us matrix: %.3e\n",
                min(eigen(f$sigmas$gid, symmetric = TRUE, only.values = TRUE)$values)))
    ll_dense <- dense_reml_us(f$sigmas$gid, diag(f$sigma_res))
    add_check(mname, "logLik vs dense REML at remlax estimates", "dense R",
              f$logLik_asreml, ll_dense, tol = 1e-5, metric = "abs",
              note = "independent evaluation of -0.5(log|V| + log|X'V^-1X| + y'Py)")
    add_check(mname, "logLik asreml (default start) minus remlax", "asreml",
              a$loglik - f$logLik_asreml, NA, tol = 1e-5, metric = "value",
              note = sprintf("asreml converged = %s; pass if asreml is not higher", a$converge))
    sv0 <- eval(as.call(c(as.list(us_call), start.values = TRUE)))$vparameters.table
    rownames(sv0) <- sv0$Component
    sv_at <- function(ridge, constraint) {
      sv <- sv0
      Se <- f$sigmas$gid + diag(ridge, nt)
      for (i in seq_len(nt)) for (j in seq_len(i)) sv[us_label(i, j), "Value"] <- Se[i, j]
      sv[sprintf("env_%s!R", env_lab), "Value"] <- diag(f$sigma_res)
      sv$Constraint <- constraint
      sv
    }
    sv <- sv_at(1e-7, "P")
    a2 <- stab(eval(as.call(c(as.list(us_call), G.param = quote(sv), R.param = quote(sv)))))
    cat(sprintf("%s: asreml started at the remlax estimates (+1e-7 I): logLik %.6f, converged %s\n",
                mname, a2$loglik, a2$converge))
    add_check(mname, "logLik asreml (remlax start) minus remlax", "asreml",
              a2$loglik - f$logLik_asreml, NA, tol = 1e-5, metric = "value",
              note = sprintf("asreml converged = %s; pass if asreml is not higher", a2$converge))
    sv <- sv_at(1e-9, "F")
    a_fix <- eval(as.call(c(as.list(us_call), G.param = quote(sv), R.param = quote(sv),
                            maxit = 1)))
    add_check(mname, "logLik asreml, all fixed at remlax estimates", "asreml",
              f$logLik_asreml, a_fix$loglik, metric = "abs", kind = "context",
              note = sprintf("Sigma + 1e-9 I; dense REML at the same point gives %.8f",
                             dense_reml_us(f$sigmas$gid + diag(1e-9, nt), diag(f$sigma_res))))
    add_blup_checks(paste(mname, "(asreml fixed at remlax estimates)"), "asreml",
                    as.vector(f$blups$gid[ids, , drop = FALSE]),
                    as.vector(asr_blup_me(a_fix, "env_")))
    vc <- asr_vc(a2)
    S_as <- asr_us(vc, nt)
    S_rx <- f$sigmas$gid
    for (i in seq_len(nt)) for (j in seq_len(i))
      add_check(paste(mname, "(asreml from remlax start)"),
                sprintf("G %s %s-%s", if (i == j) "var" else "cov", env_lab[j], env_lab[i]),
                "asreml", S_rx[i, j], S_as[i, j], kind = "context",
                note = "asreml drifts away from the boundary optimum to a lower likelihood")
    timing[[sprintf("multi-env %s x G", st)]] <- c(remlax = t_rx, asreml = t_as)
    next
  }
  vc <- asr_vc(a)
  S_as <- if (st == "us") asr_us(vc, nt) else asr_fa(vc, nt)
  S_rx <- f$sigmas$gid
  for (i in seq_len(nt)) for (j in seq_len(i)) {
    # a covariance can sit near zero: its difference is scaled by the
    # geometric mean of the two variances (a correlation scale)
    if (i == j) add_check(mname, sprintf("G var %s", env_lab[i]), "asreml", S_rx[i, j], S_as[i, j])
    else add_check(mname, sprintf("G cov %s-%s (/ sqrt(var var))", env_lab[j], env_lab[i]),
                   "asreml", abs(S_rx[i, j] - S_as[i, j]) / sqrt(S_as[i, i] * S_as[j, j]), NA,
                   tol = 1e-3, metric = "value")
  }
  R_rx <- diag(f$sigma_res)
  R_as <- asr_res(vc)
  for (i in seq_len(nt))
    add_check(mname, sprintf("residual var %s", env_lab[i]), "asreml", R_rx[i], R_as[i])
  add_check(mname, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
            tol = 1e-5, metric = "abs")
  B_as <- asr_blup_me(a, if (st == "us") "env_" else "fa(env, 1)_")
  B_rx <- f$blups$gid[ids, , drop = FALSE]
  add_blup_checks(mname, "asreml", as.vector(B_rx), as.vector(B_as))
  timing[[sprintf("multi-env %s x G", st)]] <- c(remlax = t_rx, asreml = t_as)
}

cat("\nRuntimes (s, elapsed, remlax CPU vs asreml):\n")
print(round(do.call(rbind, timing), 2))
tab <- write_checks()
sink()
