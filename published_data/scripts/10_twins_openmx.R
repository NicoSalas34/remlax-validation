# ==============================================================================
# 10_twins_openmx.R -- ACE and ADE twin models on the Australian twin data
# ------------------------------------------------------------------------------
# Data: twinData shipped with OpenMx (Neale et al. 2016, Psychometrika
# 81:535-549, doi:10.1007/s11336-014-9435-8), Australian twin pairs of the
# NH&MRC twin register with weight, height and body mass index. Female pairs,
# monozygotic (MZFF) and dizygotic (DZFF), body mass index (bmi1, bmi2).
#
# The classical twin models written as linear mixed models:
#   ACE   y = mu + a + c + e    a ~ N(0, s2A K_A), c ~ N(0, s2C I_pair)
#   ADE   y = mu + a + d + e    d ~ N(0, s2D K_D)
# K_A and K_D are defined on genotypes, not on individuals: the two twins of a
# monozygotic pair share one genotype level; a dizygotic pair has two levels,
# related by 1/2 in K_A and 1/4 in K_D. Both matrices are then positive
# definite and block diagonal, and so are their inverses. The covariance
# between the twins of a pair is s2A + s2C (MZ) and s2A/2 + s2C (DZ) under
# ACE, s2A + s2D and s2A/2 + s2D/4 under ADE, as in the classical model.
#
# References, same model and data, REML: asreml (block diagonal K), gaston
# lmm.aireml (individual-level matrices Z K Z'), sommer::mmes. OpenMx fits the
# same models by maximum likelihood in the wide (pair) format; its estimates
# differ from REML by the degrees of freedom of the mean and are reported as
# context, never counted.
#
#   Rscript published_data/scripts/10_twins_openmx.R
# ==============================================================================
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(.here, "..", "R", "common.R"))
suppressPackageStartupMessages({ library(asreml); library(sommer); library(gaston); library(OpenMx) })

start_checks("twins_openmx", "OpenMx twinData, Neale et al. 2016, doi:10.1007/s11336-014-9435-8")
cat(sprintf("  %-12s %s\n", "OpenMx", as.character(utils::packageVersion("OpenMx"))))
MLNOTE <- "OpenMx maximum likelihood on the wide pair format; REML differs by the df of the mean"

data(twinData)
tw <- twinData[twinData$zygosity %in% c("MZFF", "DZFF"), ]
tw$pair <- paste0("p", seq_len(nrow(tw)))
long <- rbind(data.frame(pair = tw$pair, zyg = tw$zygosity, twin = 1L, y = tw$bmi1),
              data.frame(pair = tw$pair, zyg = tw$zygosity, twin = 2L, y = tw$bmi2))
long <- long[!is.na(long$y), ]
long <- long[order(long$pair, long$twin), ]
mz <- long$zyg == "MZFF"
long$geno <- ifelse(mz, paste0(long$pair, "_mz"), paste0(long$pair, "_", long$twin))
cat(sprintf("pairs MZ %d, DZ %d; individuals with BMI %d\n", sum(tw$zygosity == "MZFF"),
            sum(tw$zygosity == "DZFF"), nrow(long)))

# Genotype levels and block diagonal K_A, K_D (and their inverses).
glev <- unique(long$geno)
pair_of <- sub("_(mz|1|2)$", "", glev)
build_K <- function(r) {
  K <- Matrix::Diagonal(length(glev))
  K <- methods::as(K, "generalMatrix")
  dz <- split(seq_along(glev), pair_of)
  dz <- dz[lengths(dz) == 2L]
  i <- vapply(dz, `[`, 1L, 1L); j <- vapply(dz, `[`, 1L, 2L)
  K <- K + Matrix::sparseMatrix(i = c(i, j), j = c(j, i), x = r, dims = dim(K))
  dimnames(K) <- list(glev, glev)
  K
}
KA <- build_K(0.5); KD <- build_K(0.25)
KAi <- solve(KA); KDi <- solve(KD)
dimnames(KAi) <- dimnames(KA); dimnames(KDi) <- dimnames(KD)
long$geno <- factor(long$geno, levels = glev)
long$genoD <- long$geno
long$pair <- factor(long$pair)
d <- long
cat(sprintf("genotype levels %d (MZ pairs count once)\n", length(glev)))

# Individual-level matrices for gaston: Z K Z'.
Zg <- Matrix::sparseMatrix(i = seq_len(nrow(d)), j = as.integer(d$geno), x = 1,
                           dims = c(nrow(d), length(glev)))
Zp <- Matrix::sparseMatrix(i = seq_len(nrow(d)), j = as.integer(d$pair), x = 1)
KA_ind <- as.matrix(Zg %*% KA %*% Matrix::t(Zg))
KD_ind <- as.matrix(Zg %*% KD %*% Matrix::t(Zg))
KC_ind <- as.matrix(Matrix::tcrossprod(Zp))

KA_as <- as_asreml_G(as.matrix(KA)); KD_as <- as_asreml_G(as.matrix(KD))
vcr <- function(vc, row) {
  if (!row %in% rownames(vc)) stop("no varcomp row ", row, " in: ", paste(rownames(vc), collapse = ", "))
  vc[row, "component"]
}
timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; attr(v, "secs") <- proc.time()[["elapsed"]] - t0; v }
TIMES <- list()

# ------------------------------------------------------------------------------
# OpenMx, maximum likelihood, wide format (context)
# ------------------------------------------------------------------------------
openmx_fit <- function(type = c("ACE", "ADE")) {
  type <- match.arg(type)
  sv <- c("bmi1", "bmi2")
  mzd <- tw[tw$zygosity == "MZFF", sv]; dzd <- tw[tw$zygosity == "DZFF", sv]
  rd <- if (type == "ACE") 1 else 0.25
  top <- mxModel("top",
    mxMatrix("Full", 1, 1, TRUE, 0.6, name = "a"),
    mxMatrix("Full", 1, 1, TRUE, 0.3, name = "x"),   # c under ACE, d under ADE
    mxMatrix("Full", 1, 1, TRUE, 0.6, name = "e"),
    mxMatrix("Full", 1, 2, TRUE, 21, labels = c("mu", "mu"), name = "M"),
    mxAlgebra(a %*% t(a), name = "A"), mxAlgebra(x %*% t(x), name = "X"),
    mxAlgebra(e %*% t(e), name = "E"),
    mxAlgebra(rbind(cbind(A + X + E, A + X), cbind(A + X, A + X + E)), name = "Smz"),
    mxAlgebra(rbind(cbind(A + X + E, 0.5 %x% A + rd %x% X),
                    cbind(0.5 %x% A + rd %x% X, A + X + E)), name = "Sdz"),
    mxMatrix("Full", 1, 1, FALSE, rd, name = "rd"))
  mzm <- mxModel("MZ", mxData(mzd, type = "raw"),
                 mxExpectationNormal("top.Smz", "top.M", dimnames = sv), mxFitFunctionML())
  dzm <- mxModel("DZ", mxData(dzd, type = "raw"),
                 mxExpectationNormal("top.Sdz", "top.M", dimnames = sv), mxFitFunctionML())
  m <- mxModel(type, top, mzm, dzm, mxFitFunctionMultigroup(c("MZ", "DZ")))
  r <- mxRun(m, silent = TRUE)
  c(A = mxEval(top.A, r)[1], X = mxEval(top.X, r)[1], E = mxEval(top.E, r)[1],
    minus2LL = r$output$minimum, status = r$output$status$code)
}

# ------------------------------------------------------------------------------
# ACE and ADE
# ------------------------------------------------------------------------------
models <- list(
  ACE = list(rx = ~ vm(geno, Kinv = KAi) + pair, as = ~ vm(geno, KA_as) + pair,
             so = ~ vsm(ism(geno), Gu = as.matrix(KA)) + pair,
             gaston = list(KA_ind, KC_ind),
             comp = c(A = "geno", C = "pair"),
             as_names = c(A = "vm(geno, KA_as)", C = "pair"),
             so_names = c(A = "vsm(ism(geno), Gu = as.matrix(KA))",
                          C = "sommer::vsm( sommer::ism( pair ) )")),
  ADE = list(rx = ~ vm(geno, Kinv = KAi) + vm(genoD, Kinv = KDi),
             as = ~ vm(geno, KA_as) + vm(genoD, KD_as),
             so = ~ vsm(ism(geno), Gu = as.matrix(KA)) + vsm(ism(genoD), Gu = as.matrix(KD)),
             gaston = list(KA_ind, KD_ind),
             comp = c(A = "geno", D = "genoD"),
             as_names = c(A = "vm(geno, KA_as)", D = "vm(genoD, KD_as)"),
             so_names = c(A = "vsm(ism(geno), Gu = as.matrix(KA))",
                          D = "vsm(ism(genoD), Gu = as.matrix(KD))")))

for (m in names(models)) {
  s <- models[[m]]
  cat(sprintf("\n=== %s: %s ===\n", m, deparse(s$rx)))
  f <- timed(rx_cpu(y ~ 1, random = s$rx, data = d, engine = "mme"))
  print(f)
  a <- timed(stab(asreml(y ~ 1, random = s$as, data = d, trace = FALSE, maxit = 100,
                         workspace = "2gb")))
  vc <- summary(a)$varcomp; print(vc)
  g <- timed(lmm.aireml(d$y, X = matrix(1, nrow(d)), K = s$gaston, verbose = FALSE,
                        eps = 1e-10, max_iter = 200, min_tau = 1e-10))
  so <- timed(mmes(y ~ 1, random = s$so, rcov = ~ units, data = d, verbose = FALSE,
                   dateWarning = FALSE))
  TIMES[[m]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"),
                  gaston = attr(g, "secs"), sommer = attr(so, "secs"))
  cn <- names(s$comp)
  for (k in seq_along(cn)) {
    v <- f$sigmas[[s$comp[[k]]]][1, 1]
    add_check(m, paste0("sigma2 ", cn[k]), "asreml", v, vcr(vc, s$as_names[[cn[k]]]))
    add_check(m, paste0("sigma2 ", cn[k]), "gaston", v, g$tau[k])
    j <- match(s$so_names[[cn[k]]], names(so$theta))
    add_check(m, paste0("sigma2 ", cn[k]), "sommer", v,
              if (is.na(j)) NA else as.numeric(so$theta[[j]]))
  }
  add_check(m, "sigma2 E", "asreml", f$sigma_res[1, 1], vcr(vc, "units!R"))
  add_check(m, "sigma2 E", "gaston", f$sigma_res[1, 1], g$sigma2)
  add_check(m, "sigma2 E", "sommer", f$sigma_res[1, 1],
            as.numeric(so$theta[[match("units", names(so$theta))]]))
  add_check(m, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
            tol = 1e-5, metric = "abs")
  # gaston's logL leaves out the 2 pi constant, as asreml does.
  add_check(m, "logLik (asreml convention)", "gaston", f$logLik_asreml, g$logL,
            tol = 1e-5, metric = "abs")

  om <- openmx_fit(m)
  cat(sprintf("OpenMx %s (ML): A %.5f  %s %.5f  E %.5f  -2LL %.4f  status %d\n", m,
              om[["A"]], cn[2], om[["X"]], om[["E"]], om[["minus2LL"]], om[["status"]]))
  add_check(m, paste0("sigma2 ", cn[1]), "OpenMx (ML)", f$sigmas[[s$comp[[1]]]][1, 1], om[["A"]],
            kind = "context", note = MLNOTE)
  add_check(m, paste0("sigma2 ", cn[2]), "OpenMx (ML)", f$sigmas[[s$comp[[2]]]][1, 1], om[["X"]],
            kind = "context", note = MLNOTE)
  add_check(m, "sigma2 E", "OpenMx (ML)", f$sigma_res[1, 1], om[["E"]],
            kind = "context", note = MLNOTE)
}

cat("\nRuntimes (s, elapsed):\n")
print(do.call(rbind, TIMES))
write_checks()
