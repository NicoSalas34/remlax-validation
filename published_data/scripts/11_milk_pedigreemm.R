# ==============================================================================
# 11_milk_pedigreemm.R -- repeated lactations: animal model with permanent
# environment, on the Holstein data of pedigreemm
# ------------------------------------------------------------------------------
# Data: milk and pedCowsR shipped with pedigreemm (Vazquez et al. 2010, J.
# Anim. Sci. 88:497-504, doi:10.2527/jas.2009-1952). 3397 lactations (first to
# fifth) of 1359 Holstein cows, daughters of 38 sires in 57 herds, downloaded
# by the authors from the USDA site; the pedigree has 6547 animals. Milk and
# protein yields at 305 days are divided by 1000 so that the variances are of
# order one. The help page says kg, but the values look like pounds; the unit
# does not matter for the comparison.
#
# Models, with the additive relationship of the full pedigree (A^-1 from
# pedigreemm::getAInv, passed to every program):
#   animal         y ~ lact + log(dim) + a + herd
#   repeatability  y ~ lact + log(dim) + a + pe + herd
#                  pe: permanent environment of the cow, iid, on the same id
#   bivariate      milk and protein: us(trait) on a and pe, iid herd per
#                  trait, us(trait) residual
#
# References, same model and data, REML: asreml (ginv), pedigreemm (lme4 REML
# on the Cholesky factor of A) for the univariate models.
#
#   Rscript published_data/scripts/11_milk_pedigreemm.R
# ==============================================================================
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(.here, "..", "R", "common.R"))
suppressPackageStartupMessages({ library(asreml); library(pedigreemm) })

start_checks("milk_pedigreemm", "Vazquez et al. 2010, doi:10.2527/jas.2009-1952")
cat(sprintf("  %-12s %s\n", "pedigreemm", as.character(utils::packageVersion("pedigreemm"))))

data(milk)
labs <- pedCowsR@label
Ai <- getAInv(pedCowsR)
dimnames(Ai) <- list(labs, labs)
d <- milk
d$id <- factor(as.character(d$id))
stopifnot(all(levels(d$id) %in% labs))
d$pe <- d$id
d$lact <- factor(d$lact)
d$ldim <- log(d$dim)
d$milk_t <- d$milk / 1000; d$prot_t <- d$prot / 1000
cat(sprintf("records %d, cows %d, herds %d, pedigree %d\n", nrow(d), nlevels(d$id),
            nlevels(d$herd), length(labs)))

ginv <- function(M) {
  T <- methods::as(methods::as(Matrix::forceSymmetric(M, "L"), "generalMatrix"), "TsparseMatrix")
  k <- T@i >= T@j
  g <- cbind(Row = T@i[k] + 1L, Column = T@j[k] + 1L, Value = T@x[k])
  g <- g[order(g[, "Row"], g[, "Column"]), ]
  attr(g, "rowNames") <- rownames(M); attr(g, "INVERSE") <- TRUE
  g
}
gA <- ginv(Ai)
vcr <- function(vc, row) {
  if (!row %in% rownames(vc)) stop("no varcomp row ", row, " in: ", paste(rownames(vc), collapse = ", "))
  vc[row, "component"]
}
asr_blup <- function(a, prefix) {
  cr <- summary(a, coef = TRUE)$coef.random
  k <- startsWith(rownames(cr), prefix)
  setNames(cr[k, "solution"], sub(prefix, "", rownames(cr)[k], fixed = TRUE))
}
timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; attr(v, "secs") <- proc.time()[["elapsed"]] - t0; v }
TIMES <- list()

# ------------------------------------------------------------------------------
# Univariate: animal and repeatability models, milk and protein
# ------------------------------------------------------------------------------
for (tr in c("milk_t", "prot_t")) for (rep_model in c(FALSE, TRUE)) {
  m <- paste(tr, if (rep_model) "repeatability" else "animal")
  fx <- stats::as.formula(paste(tr, "~ lact + ldim"))
  rx_r <- if (rep_model) ~ vm(id, Kinv = Ai) + pe + herd else ~ vm(id, Kinv = Ai) + herd
  as_r <- if (rep_model) ~ vm(id, gA) + pe + herd else ~ vm(id, gA) + herd
  pm_f <- stats::as.formula(paste(tr, "~ lact + ldim + (1 | id) + (1 | herd)",
                                  if (rep_model) "+ (1 | pe)" else ""))
  cat(sprintf("\n=== %s ===\n", m))
  f <- timed(rx_cpu(fx, random = rx_r, data = d, engine = "mme"))
  print(f)
  a <- timed(stab(asreml(fx, random = as_r, data = d, trace = FALSE, maxit = 100,
                         workspace = "2gb")))
  vc <- summary(a)$varcomp; print(vc)
  p <- timed(pedigreemm(pm_f, data = d, pedigree = list(id = pedCowsR), REML = TRUE))
  vp <- as.data.frame(lme4::VarCorr(p))
  vpg <- function(g) vp$vcov[vp$grp == g]
  TIMES[[m]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"),
                  pedigreemm = attr(p, "secs"))
  comps <- c(A = "id", herd = "herd", if (rep_model) c(PE = "pe"))
  as_nm <- c(A = "vm(id, gA)", herd = "herd", PE = "pe")
  for (cn in names(comps)) {
    v <- f$sigmas[[comps[[cn]]]][1, 1]
    add_check(m, paste0("sigma2 ", cn), "asreml", v, vcr(vc, as_nm[[cn]]))
    add_check(m, paste0("sigma2 ", cn), "pedigreemm", v, vpg(comps[[cn]]))
  }
  add_check(m, "sigma2 residual", "asreml", f$sigma_res[1, 1], vcr(vc, "units!R"))
  add_check(m, "sigma2 residual", "pedigreemm", f$sigma_res[1, 1], vpg("Residual"))
  add_check(m, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
            tol = 1e-5, metric = "abs")
  # lme4 reports the REML log-likelihood with its 2 pi constant.
  add_check(m, "logLik (with constant)", "pedigreemm", f$logLik,
            as.numeric(stats::logLik(p)), tol = 1e-4, metric = "abs")
  u_rx <- setNames(f$blups$id[, 1], rownames(f$blups$id))
  u_as <- asr_blup(a, "vm(id, gA)_")[names(u_rx)]
  add_blup_checks(paste(m, "A"), "asreml", u_rx, u_as)
  bx <- coef(summary(p))
  bx_rx <- as.numeric(f$beta)
  add_check(m, "fixed effects max |diff| / se", "pedigreemm",
            max(abs(bx_rx - bx[, "Estimate"]) / bx[, "Std. Error"]), NA,
            tol = 1e-2, metric = "value")
}

# ------------------------------------------------------------------------------
# Bivariate milk + protein, repeatability model
# ------------------------------------------------------------------------------
mname <- "bivariate milk+prot repeatability"
cat(sprintf("\n=== %s ===\n", mname))
d$rec <- factor(seq_len(nrow(d)))
dl <- rbind(data.frame(d[, c("id", "pe", "herd", "lact", "ldim", "rec")], car = "milk", y = d$milk_t),
            data.frame(d[, c("id", "pe", "herd", "lact", "ldim", "rec")], car = "prot", y = d$prot_t))
dl$car <- factor(dl$car, levels = c("milk", "prot"))
dl <- dl[order(dl$rec, dl$car), ]
f <- timed(rx_cpu(y ~ car + car:lact + car:ldim,
                  random = ~ us(id, Kinv = Ai) + us(pe) + diag(herd),
                  residual = ~ us(car):units, data = dl, trait = "car", unit = "rec",
                  engine = "mme"))
print(f)
biv_call <- quote(asreml(cbind(milk_t, prot_t) ~ trait + trait:lact + trait:ldim,
                         random = ~ us(trait):vm(id, gA) + us(trait):pe + diag(trait):herd,
                         residual = ~ units:us(trait), data = d, trace = FALSE,
                         maxit = 100, workspace = "4gb"))
a <- timed(stab(eval(biv_call)))
vc <- summary(a)$varcomp; print(vc)
TIMES[[mname]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"))
tl <- c("milk_t", "prot_t")
get_us <- function(term) {
  S <- matrix(NA_real_, 2, 2)
  for (i in 1:2) for (j in 1:i) {
    k <- match(sprintf("%s!trait_%s:%s", term, tl[i], tl[j]), rownames(vc))
    S[i, j] <- S[j, i] <- vc[k[1], "component"]
  }
  S
}
# The additive and permanent environment terms of a cow are nearly confounded
# in the repeatability model: the likelihood is flat along their sum. The
# arbiter is the REML log-likelihood computed densely in R (asreml
# convention) at the estimates of each program.
for (blk in list(list("A", f$sigmas$id, get_us("trait:vm(id, gA)")),
                 list("PE", f$sigmas$pe, get_us("trait:pe")),
                 list("residual", f$sigma_res, get_us("units:trait")))) {
  Srx <- blk[[2]]; Sas <- blk[[3]]
  add_check(mname, paste(blk[[1]], "var milk"), "asreml", Srx[1, 1], Sas[1, 1])
  add_check(mname, paste(blk[[1]], "var prot"), "asreml", Srx[2, 2], Sas[2, 2])
  add_check(mname, paste(blk[[1]], "cov (/ sqrt(var var))"), "asreml",
            abs(Srx[1, 2] - Sas[1, 2]) / sqrt(Sas[1, 1] * Sas[2, 2]), NA,
            tol = 1e-3, metric = "value")
}
A_full <- as.matrix(getA(pedCowsR)); dimnames(A_full) <- list(labs, labs)
cows <- levels(d$id); Ac <- A_full[cows, cows]; rm(A_full)
Xl <- stats::model.matrix(~ car + car:lact + car:ldim, dl); yl <- dl$y; nl <- nrow(dl)
tix <- as.integer(dl$car)
Zc <- Matrix::sparseMatrix(i = seq_len(nl), j = match(as.character(dl$id), cows), x = 1,
                           dims = c(nl, length(cows)))
Zh <- Matrix::sparseMatrix(i = seq_len(nl), j = as.integer(dl$herd), x = 1)
Zr <- Matrix::sparseMatrix(i = seq_len(nl), j = as.integer(dl$rec), x = 1)
dense_reml_biv <- function(SA, SP, SH, SR) {
  o <- function(S) outer(tix, tix, function(i, j) S[cbind(i, j)])
  V <- o(SA) * as.matrix(Zc %*% Ac %*% Matrix::t(Zc)) + o(SP) * as.matrix(Matrix::tcrossprod(Zc)) +
       o(SH) * as.matrix(Matrix::tcrossprod(Zh)) + o(SR) * as.matrix(Matrix::tcrossprod(Zr))
  R <- chol(V); Vi <- chol2inv(R); XtViX <- crossprod(Xl, Vi %*% Xl)
  P <- Vi - Vi %*% Xl %*% solve(XtViX, crossprod(Xl, Vi))
  -0.5 * (2 * sum(log(diag(R))) + as.numeric(determinant(XtViX)$modulus) + drop(t(yl) %*% P %*% yl))
}
ll_rx <- dense_reml_biv(f$sigmas$id, f$sigmas$pe, f$sigmas$herd, f$sigma_res)
ll_as <- dense_reml_biv(get_us("trait:vm(id, gA)"), get_us("trait:pe"),
                        diag(c(vcr(vc, "trait:herd!trait_milk_t"), vcr(vc, "trait:herd!trait_prot_t"))),
                        get_us("units:trait"))
add_check(mname, "logLik reported vs dense R at the remlax estimates", "dense R",
          f$logLik_asreml, ll_rx, tol = 1e-6, metric = "abs")
add_check(mname, "logLik asreml vs dense R at the asreml estimates", "dense R",
          a$loglik, ll_as, tol = 1e-6, metric = "abs", kind = "context",
          note = "checks the arbiter on the asreml side")
add_check(mname, "dense R: remlax estimates minus asreml estimates", "dense R",
          ll_rx - ll_as, 0, tol = -1e-6, metric = "ge",
          note = "pass if the remlax point is not lower than the asreml point")
add_check(mname, "herd var milk", "asreml", f$sigmas$herd[1, 1],
          vcr(vc, sprintf("trait:herd!trait_%s", tl[1])))
add_check(mname, "herd var prot", "asreml", f$sigmas$herd[2, 2],
          vcr(vc, sprintf("trait:herd!trait_%s", tl[2])))
add_check(mname, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
          tol = 1e-5, metric = "abs")

cat("\nRuntimes (s, elapsed):\n")
print(TIMES)
write_checks()
