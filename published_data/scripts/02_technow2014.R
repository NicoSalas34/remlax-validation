# ==============================================================================
# 02_technow2014.R -- hybrid GBLUP with two GCA kernels on the maize single
#                     crosses of Technow et al. (2014)
# ------------------------------------------------------------------------------
# DATA. enhancer::DT_technow (the data sets formerly shipped with sommer):
# 1254 single-cross hybrids between 123 dent and 86 flint inbred lines, with
# adjusted means for grain yield (GY) and grain moisture (GM), and 35478 SNPs
# for each heterotic group (Md_technow, Mf_technow, coded 0/1). Technow et al.
# 2014, Genetics 197:1343-1355, doi:10.1534/genetics.114.165860.
#
# RELATIONSHIP MATRICES. The parents are inbred, so a coded allele is a
# homozygous dose: the markers enter VanRaden's formula as 0/2, one matrix per
# heterotic group, built once and passed unchanged to every program.
#
# MODELS
#   1. GY and GM, each:   y = 1mu + gca_dent + gca_flint + e
#                         gca_dent ~ N(0, s2d Gd), gca_flint ~ N(0, s2f Gf)
#   2. GY + GM bivariate: us(trait) x Gd + us(trait) x Gf, us residual
#
# With one record per hybrid, the specific combining ability is confounded
# with the residual and is not fitted.
#
# REFERENCES. asreml for every model, sommer::mmes for the univariate ones.
#
#   systemd-run --user --scope -q -p MemoryMax=7G \
#     Rscript published_data/scripts/02_technow2014.R
# ==============================================================================
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(.self), "..", "R", "common.R"))
suppressPackageStartupMessages({
  library(asreml); library(sommer); library(enhancer)
})
asreml.options(workspace = "4gb", trace = FALSE, maxit = 100)

log_file <- file.path(pd_root(), "results", "maize_technow2014.log")
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
sink(log_file, split = TRUE)

start_checks("maize_technow2014", "Technow et al. 2014, doi:10.1534/genetics.114.165860")
data(DT_technow, package = "enhancer")
to_dose <- function(M) {
  X <- matrix(as.integer(M), nrow(M), ncol(M), dimnames = dimnames(M))
  2 * X
}
Gd <- vanraden(to_dose(Md_technow))
Gf <- vanraden(to_dose(Mf_technow))
cat(sprintf("Gd: %d x %d, Gf: %d x %d\n", nrow(Gd), ncol(Gd), nrow(Gf), ncol(Gf)))

d <- data.frame(hy = factor(DT_technow$hy),
                dent = factor(as.character(DT_technow$dent), levels = rownames(Gd)),
                flint = factor(as.character(DT_technow$flint), levels = rownames(Gf)),
                GY = DT_technow$GY, GM = DT_technow$GM)
stopifnot(!anyNA(d$dent), !anyNA(d$flint))
cat(sprintf("%d hybrids, %d dent and %d flint parents used\n", nrow(d),
            length(unique(d$dent)), length(unique(d$flint))))
Gda <- as_asreml_G(Gd); Gfa <- as_asreml_G(Gf)
dent_ids <- rownames(Gd); flint_ids <- rownames(Gf)

asr_blup <- function(a, pattern, lev) {
  cr <- coef(a)$random
  rn <- rownames(cr)
  sel <- grepl(pattern, rn, fixed = TRUE)
  v <- cr[sel, 1]; names(v) <- sub(".*_", "", rn[sel]); v[lev]
}
timing <- list()

# ------------------------------------------------------------------------------
# 1. univariate, two GCA kernels
# ------------------------------------------------------------------------------
for (tr in c("GY", "GM")) {
  mname <- sprintf("2 GCA kernels %s", tr)
  dd <- d; dd$y <- dd[[tr]]
  t0 <- proc.time()[3]
  f <- rx_cpu(y ~ 1, random = ~ vm(dent, K = Gd) + vm(flint, K = Gf), data = dd,
              n_restarts = 3)
  t_rx <- proc.time()[3] - t0
  t0 <- proc.time()[3]
  a <- stab(asreml(y ~ 1, random = ~ vm(dent, Gda) + vm(flint, Gfa), data = dd))
  t_as <- proc.time()[3] - t0
  vc <- summary(a)$varcomp
  add_check(mname, "sigma2_GCA dent", "asreml", f$sigmas$dent[1, 1], vc["vm(dent, Gda)", "component"])
  add_check(mname, "sigma2_GCA flint", "asreml", f$sigmas$flint[1, 1], vc["vm(flint, Gfa)", "component"])
  add_check(mname, "sigma2_e", "asreml", f$sigma_res[1, 1], vc["units!R", "component"])
  add_check(mname, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
            tol = 1e-5, metric = "abs")
  add_blup_checks(paste(mname, "dent"), "asreml", f$blups$dent[dent_ids, 1],
                  asr_blup(a, "vm(dent, Gda)_", dent_ids))
  add_blup_checks(paste(mname, "flint"), "asreml", f$blups$flint[flint_ids, 1],
                  asr_blup(a, "vm(flint, Gfa)_", flint_ids))
  timing[[mname]] <- c(remlax = t_rx, asreml = t_as)

  s <- mmes(y ~ 1, random = ~ vsm(ism(dent), Gu = Gd) + vsm(ism(flint), Gu = Gf),
            rcov = ~ units, data = dd, nIters = 200, tolParConvLL = 1e-9,
            verbose = FALSE, dateWarning = FALSE)
  th <- unlist(lapply(s$theta, function(m) m[1, 1]))
  add_check(mname, "sigma2_GCA dent", "sommer", f$sigmas$dent[1, 1], th[1])
  add_check(mname, "sigma2_GCA flint", "sommer", f$sigmas$flint[1, 1], th[2])
  add_check(mname, "sigma2_e", "sommer", f$sigma_res[1, 1], th[3])
  ud <- s$uList[[1]]; uf <- s$uList[[2]]
  add_blup_checks(paste(mname, "dent"), "sommer", f$blups$dent[dent_ids, 1],
                  setNames(as.numeric(ud[, 1]), rownames(ud))[dent_ids])
  add_blup_checks(paste(mname, "flint"), "sommer", f$blups$flint[flint_ids, 1],
                  setNames(as.numeric(uf[, 1]), rownames(uf))[flint_ids])
}

# ------------------------------------------------------------------------------
# 2. bivariate GY + GM, us genetic covariances in each group, us residual
# ------------------------------------------------------------------------------
dl <- rbind(data.frame(d[, c("hy", "dent", "flint")], trait = "GY", y = d$GY),
            data.frame(d[, c("hy", "dent", "flint")], trait = "GM", y = d$GM))
dl$trait <- factor(dl$trait, levels = c("GY", "GM"))
dl <- dl[order(dl$hy, dl$trait), ]
mname <- "bivariate GY+GM"
t0 <- proc.time()[3]
f <- rx_cpu(y ~ trait, random = ~ us(dent, K = Gd) + us(flint, K = Gf),
            residual = ~ us(trait):units, data = dl, trait = "trait", unit = "hy",
            n_restarts = 3)
t_rx <- proc.time()[3] - t0
t0 <- proc.time()[3]
# asreml side in WIDE format, its native multivariate form. In long format
# a column named `trait` collides with the reserved asreml keyword and gives a
# wrong model without any error (measured here: logLik -1552.34 for the us
# model and -1568.03 for the diag model, whose value must equal the sum of the
# two univariate fits, -3143.57).
a <- stab(asreml(cbind(GY, GM) ~ trait,
                 random = ~ us(trait):vm(dent, Gda) + us(trait):vm(flint, Gfa),
                 residual = ~ units:us(trait), data = d))
t_as <- proc.time()[3] - t0
vc <- summary(a)$varcomp
print(vc)
tl <- levels(dl$trait)
get_us <- function(term) {
  S <- matrix(NA_real_, 2, 2)
  for (i in 1:2) for (j in 1:i) {
    k <- match(sprintf("%s!trait_%s:%s", term, tl[i], tl[j]), rownames(vc))
    S[i, j] <- S[j, i] <- vc[k[1], "component"]
  }
  S
}
S_d <- get_us("trait:vm(dent, Gda)")
S_f <- get_us("trait:vm(flint, Gfa)")
S_r <- get_us("units:trait")
for (blk in list(list("GCA dent", f$sigmas$dent, S_d), list("GCA flint", f$sigmas$flint, S_f),
                 list("residual", f$sigma_res, S_r))) {
  Srx <- blk[[2]]; Sas <- blk[[3]]
  add_check(mname, paste(blk[[1]], "var GY"), "asreml", Srx[1, 1], Sas[1, 1])
  add_check(mname, paste(blk[[1]], "var GM"), "asreml", Srx[2, 2], Sas[2, 2])
  add_check(mname, paste(blk[[1]], "cov (/ sqrt(var var))"), "asreml",
            abs(Srx[1, 2] - Sas[1, 2]) / sqrt(Sas[1, 1] * Sas[2, 2]), NA,
            tol = 1e-3, metric = "value")
}
add_check(mname, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
          tol = 1e-5, metric = "abs")
cr <- coef(a)$random
for (grp in c("dent", "flint")) {
  lev <- if (grp == "dent") dent_ids else flint_ids
  term <- sprintf("vm(%s, %s)", grp, if (grp == "dent") "Gda" else "Gfa")
  B_as <- sapply(tl, function(t) {
    sel <- startsWith(rownames(cr), sprintf("trait_%s:%s_", t, term))
    v <- cr[sel, 1]; names(v) <- sub(".*_", "", rownames(cr)[sel]); v[lev]
  })
  add_blup_checks(paste(mname, grp), "asreml", as.vector(f$blups[[grp]][lev, ]),
                  as.vector(B_as))
}
timing[[mname]] <- c(remlax = t_rx, asreml = t_as)

cat("\nRuntimes (s, elapsed, remlax CPU vs asreml):\n")
print(round(do.call(rbind, timing), 2))
tab <- write_checks()
sink()
