# ==============================================================================
# 12_tilapia_joshi2018.R -- additive, dominance and maternal variance in Nile
# tilapia, against the ASReml values printed by the authors
# ------------------------------------------------------------------------------
# Data: Joshi, Woolliams, Meuwissen & Gjoen (2018, Heredity 120:452-462,
# doi:10.1038/s41437-017-0046-x, CC BY 4.0), supplementary files 5
# (phenotypes, 2524 fish of 155 full-sib families) and 6 (pedigree, 4050
# animals). GIFT strain, factorial mating with reciprocal crosses.
#
# Model of the paper (Statistical analysis section), fitted with ASReml-4:
#   ADM  y = Xb + Z1 a + Z2 d + Z3 m + e
#   a ~ N(0, A s2A), d ~ N(0, D s2D) on every animal of the pedigree,
#   m ~ N(0, I s2M), e ~ N(0, I s2E).
# A is the numerator relationship matrix, D the matrix of coefficients of
# fraternity, its inverse computed with nadiv. Sub-models AD, AM and A drop
# m, d or both.
#
# Two details of the published fits differ from the text, and were found by
# refitting: the paper's values are reproduced only with them.
#   - the maternal effect has one level per dam (18 dams), not per full-sib
#     family (155 families) as the text says;
#   - the day enters as a numeric covariate, with batch as a factor: rc +
#     batch + Day1 for body weight, length and depth; rc + batch + Day2 for
#     body thickness (measured on the filleting day); rc + batch + Day2 +
#     filleter for fillet weight and yield. The text gives day 7 (or 10) degrees of
#     freedom; with day as an 8-level factor the variances come out about
#     3 % lower and the log-likelihoods do not match.
#
# Published values compared, at the rounding of the paper:
#   Table 2      log-likelihood of every model for the six traits. The printed
#                values equal the REML log-likelihood (asreml convention)
#                plus an exact multiple of 100 (200 for BD, 12200 for BWH):
#                the check is on the remainder modulo 100, which also makes
#                the differences between models exact.
#   Table S3.1   variance components (additive, dominance, maternal,
#                residual) of every model; "0 [0]" marks a component bound to
#                zero, which remlax must also put at its floor.
# The A*D*M* models (pedigree re-based at generation 20) and the sire-dam
# models are not refitted: the files do not give the generation of each fish.
#
# Software reference, same model and data: asreml-R, with the same inverses.
#
#   Rscript published_data/scripts/12_tilapia_joshi2018.R
# ==============================================================================
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(.here, "..", "R", "common.R"))
suppressPackageStartupMessages({ library(asreml); library(nadiv) })

start_checks("tilapia_joshi2018", "Joshi et al. 2018, doi:10.1038/s41437-017-0046-x")
cat(sprintf("  %-12s %s\n", "nadiv", as.character(utils::packageVersion("nadiv"))))
PAPER <- "Joshi et al. 2018"

base <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41437-017-0046-x/MediaObjects/"
files <- c(pheno = "41437_2017_46_MOESM5_ESM.csv", ped = "41437_2017_46_MOESM6_ESM.csv")
for (k in names(files)) {
  dest <- pd_cache("tilapia_joshi2018", files[[k]])
  if (!file.exists(dest)) utils::download.file(paste0(base, files[[k]]), dest, mode = "wb", quiet = TRUE)
}
ph  <- utils::read.csv(pd_cache("tilapia_joshi2018", files[["pheno"]]), skip = 1, check.names = FALSE)
ped <- utils::read.csv(pd_cache("tilapia_joshi2018", files[["ped"]]), skip = 1)
names(ph) <- c("ID", "Sire", "Dam", "batch", "BWH", "BL", "BD", "BT", "FW", "Filleter",
               "FY", "Day1", "Day2", "fam", "rc")
ped <- ped[, c("ID", "Sire", "Dam")]
stopifnot(all(ph$ID %in% ped$ID))

pp  <- prepPed(ped)
ids <- as.character(pp$ID)
Ai  <- makeAinv(pp)$Ainv
Dl  <- suppressWarnings(makeD(pp, invertD = TRUE))
Di  <- Dl$Dinv
A   <- as.matrix(makeA(pp)); D <- as.matrix(Dl$D)
dimnames(Ai) <- dimnames(Di) <- dimnames(A) <- dimnames(D) <- list(ids, ids)
cat(sprintf("density of D^-1: %.3f (the pedigree is deep: most pairs share a coefficient of fraternity)\n",
            Matrix::nnzero(Di) / length(ids)^2))
cat(sprintf("fish %d, families %d, pedigree %d\n", nrow(ph), length(unique(ph$fam)), length(ids)))

ginv <- function(M) {
  T <- methods::as(methods::as(Matrix::forceSymmetric(M, "L"), "generalMatrix"), "TsparseMatrix")
  k <- T@i >= T@j
  g <- cbind(Row = T@i[k] + 1L, Column = T@j[k] + 1L, Value = T@x[k])
  g <- g[order(g[, "Row"], g[, "Column"]), ]
  attr(g, "rowNames") <- rownames(M); attr(g, "INVERSE") <- TRUE
  g
}
gA <- ginv(Ai); gD <- ginv(Di)

d <- ph
d$id   <- factor(as.character(d$ID), levels = ids)
d$idD  <- d$id
d$fam  <- factor(d$fam); d$dam <- factor(d$Dam)
d$rc   <- factor(d$rc); d$batch <- factor(d$batch)
d$Filleter <- factor(d$Filleter)   # Day1 and Day2 stay numeric (see the header)

# Printed values. Components: Table S3.1; NA = not in the model. Log-likelihoods: Table 2.
pub_vc <- list(
  BD  = list(A = c(0.26, NA, NA, 0.32), AD = c(0.15, 0.15, NA, 0.27),
             ADM = c(0, 0.14, 0.05, 0.35), AM = c(0.09, NA, 0.05, 0.40)),
  BL  = list(A = c(0.97, NA, NA, 2.33), AD = c(0.74, 0.33, NA, 2.20),
             ADM = c(0.30, 0.36, 0.12, 2.39), AM = c(0.55, NA, 0.12, 2.52)),
  BT  = list(A = c(1.97, NA, NA, 7.44), AD = c(0.96, 1.51, NA, 6.86),
             ADM = c(0, 1.53, 0.29, 7.30), AM = c(1.05, NA, 0.28, 7.87)),
  BWH = list(A = c(3216, NA, NA, 4312), AD = c(2085, 1688, NA, 3673),
             ADM = c(0, 1548, 637, 4739), AM = c(1008, NA, 635, 5342)),
  FW  = list(A = c(503, NA, NA, 626), AD = c(425, 126, NA, 575),
             ADM = c(0, 152, 119, 755), AM = c(104, NA, 116, 812)),
  FY  = list(A = c(2.28, NA, NA, 7.34), AD = c(2.28, 0, NA, 7.34),
             ADM = c(2.11, 0, 0.05, 7.42), AM = c(2.11, NA, 0.05, 7.42)))
# Unit of the last printed digit of the components.
unit_vc <- c(BD = 0.01, BL = 0.01, BT = 0.01, BWH = 1, FW = 1, FY = 0.01)
pub_ll <- list(
  BD  = c(A = -81.786, AD = -79.494, AM = -78.220, ADM = -75.660),
  BL  = c(A = -90.60,  AD = -90.34,  AM = -89.38,  ADM = -89.05),
  BT  = c(A = -56.94,  AD = -56.14,  AM = -55.55,  ADM = -54.57),
  BWH = c(A = -28.34,  AD = -26.66,  AM = -24.15,  ADM = -22.19),
  FW  = c(A = -21.95,  AD = -21.46,  AM = -15.59,  ADM = -14.80),
  FY  = c(A = -68.37,  AD = -68.37,  AM = -68.33,  ADM = -68.33))
unit_ll <- c(BD = 0.001, BL = 0.01, BT = 0.01, BWH = 0.01, FW = 0.01, FY = 0.01)

# D^-1 is almost full here, so the sparse mixed model equations gain nothing.
# remlax uses its dense engine on the 2524 records, with A and D themselves.
rand_rx <- list(A = ~ vm(id, K = A), AD = ~ vm(id, K = A) + vm(idD, K = D),
                AM = ~ vm(id, K = A) + dam,
                ADM = ~ vm(id, K = A) + vm(idD, K = D) + dam)
rand_as <- list(A = ~ vm(id, gA), AD = ~ vm(id, gA) + vm(idD, gD), AM = ~ vm(id, gA) + dam,
                ADM = ~ vm(id, gA) + vm(idD, gD) + dam)
comp_rx <- c(A = "id", D = "idD", M = "dam")
comp_as <- c(A = "vm(id, gA)", D = "vm(idD, gD)", M = "dam")
# Independent arbiter: the REML log-likelihood (asreml convention) computed
# densely in R at the remlax estimates.
dense_reml <- function(dt, fx, f, m) {
  X <- stats::model.matrix(fx, dt); y <- dt[[all.vars(fx)[1]]]
  ii <- as.character(dt$ID); V <- diag(f$sigma_res[1, 1], nrow(dt))
  V <- V + f$sigmas$id[1, 1] * A[ii, ii]
  if (grepl("D", m)) V <- V + f$sigmas$idD[1, 1] * D[ii, ii]
  if (grepl("M", m)) { Zm <- stats::model.matrix(~ 0 + dam, dt); V <- V + f$sigmas$dam[1, 1] * tcrossprod(Zm) }
  R <- chol(V); Vi <- chol2inv(R); XtViX <- crossprod(X, Vi %*% X)
  P <- Vi - Vi %*% X %*% solve(XtViX, crossprod(X, Vi))
  -0.5 * (2 * sum(log(diag(R))) + as.numeric(determinant(XtViX)$modulus) + drop(t(y) %*% P %*% y))
}
vc_or0 <- function(vc, row) if (row %in% rownames(vc)) vc[row, "component"] else NA_real_
asr_blup <- function(a, prefix) {
  cr <- summary(a, coef = TRUE)$coef.random
  k <- startsWith(rownames(cr), prefix)
  setNames(cr[k, "solution"], sub(prefix, "", rownames(cr)[k], fixed = TRUE))
}
TIMES <- list()

for (tr in names(pub_vc)) {
  fx <- stats::as.formula(paste(tr, switch(tr, FW = , FY = "~ rc + batch + Day2 + Filleter",
                                           BT = "~ rc + batch + Day2", "~ rc + batch + Day1")))
  dt <- d[!is.na(d[[tr]]), ]
  for (m in c("A", "AD", "AM", "ADM")) {
    lab <- paste(tr, m)
    cat(sprintf("\n=== %s ===\n", lab))
    t0 <- proc.time()[["elapsed"]]
    f <- rx_cpu(fx, random = rand_rx[[m]], data = dt)
    t_rx <- proc.time()[["elapsed"]] - t0
    print(f)
    t0 <- proc.time()[["elapsed"]]
    a <- stab(asreml(fx, random = rand_as[[m]], data = dt, trace = FALSE, maxit = 100,
                     workspace = "2gb"))
    t_as <- proc.time()[["elapsed"]] - t0
    vc <- summary(a)$varcomp; print(vc)
    TIMES[[lab]] <- c(remlax = t_rx, asreml = t_as)
    pv <- setNames(pub_vc[[tr]][[m]], c("A", "D", "M", "residual"))
    for (cn in names(pv)) {
      if (is.na(pv[[cn]])) next
      v <- if (cn == "residual") f$sigma_res[1, 1] else f$sigmas[[comp_rx[[cn]]]][1, 1]
      ref <- if (cn == "residual") vc["units!R", "component"] else vc_or0(vc, comp_as[[cn]])
      # A component at its bound on both sides is compared on the scale of
      # the phenotypic variance: asreml stops at a small positive value
      # (0.0012 for a phenotypic variance of 7000), remlax at its floor.
      vp <- sum(vc[, "component"])
      if (max(v, ref) < 1e-5 * vp)
        add_check(lab, paste0("sigma2 ", cn, " (at bound, / phenotypic var)"), "asreml",
                  abs(v - ref) / vp, NA, tol = 1e-6, metric = "value")
      else add_check(lab, paste0("sigma2 ", cn), "asreml", v, ref)
      # At the rounding of the paper: half a unit of the last printed digit,
      # plus 0.2 % for the stopping rule of ASReml-4.
      add_check(lab, paste0("sigma2 ", cn), "paper (Table S3.1)", v, pv[[cn]],
                tol = 0.5 * unit_vc[[tr]] + 2e-3 * abs(pv[[cn]]), metric = "abs",
                kind = "published", note = PAPER)
    }
    # asreml's log-likelihood with these supplied inverses sits 0.7e-5 to
    # 1.4e-5 away from the dense REML value at the same estimates; remlax
    # matches the dense value to 1e-9. The dense row is the counted check.
    add_check(lab, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
              tol = 1e-5, metric = "abs", kind = "context",
              note = "asreml loglik with the nadiv ginv is about 1e-5 off the dense REML value")
    add_check(lab, "logLik (asreml convention)", "dense R", f$logLik_asreml,
              dense_reml(dt, fx, f, m), tol = 1e-6, metric = "abs")
    off <- f$logLik_asreml - pub_ll[[tr]][[m]]
    add_check(lab, "logLik remainder modulo 100", "paper (Table 2)", abs(off - 100 * round(off / 100)),
              NA, tol = 0.5 * unit_ll[[tr]] + 2e-3, metric = "value", kind = "published",
              note = sprintf("%s; printed value = remlax + %g", PAPER, -100 * round(off / 100)))
    u_rx <- setNames(f$blups$id[, 1], rownames(f$blups$id))
    if (stats::sd(u_rx) > 1e-6)
      add_blup_checks(paste(lab, "A"), "asreml", u_rx, asr_blup(a, "vm(id, gA)_")[names(u_rx)])
  }
}

cat("\nRuntimes (s, elapsed):\n")
print(round(do.call(rbind, TIMES), 2))
write_checks()
