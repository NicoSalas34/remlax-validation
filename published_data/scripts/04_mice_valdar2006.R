# ==============================================================================
# 04_mice_valdar2006.R -- genomic BLUP on the heterogeneous-stock mice
# ------------------------------------------------------------------------------
# Data: the outbred heterogeneous-stock mice of Valdar et al. (2006, Nat.
# Genet. 38:879-887, doi:10.1038/ng1840), as shipped in BGLR (mice.X: 1814
# mice x 10346 SNPs coded 0/1/2; mice.A: pedigree relationship matrix;
# mice.pheno: phenotypes, sex, litter and cage). The same population served for
# the first empirical test of genomic selection (Legarra et al. 2008, Genetics
# 180:611-618, doi:10.1534/genetics.108.088575) and for the BGLR examples.
#
# Models, identical in every program (the genomic matrix G is built once, by
# VanRaden's first method, and passed to all of them):
#   gblup_cage   y ~ GENDER + Litter + cage (iid) + G
#   ped_gen      y ~ GENDER + Litter + cage + A (pedigree) + G (two kernels)
#   gblup_only   y ~ GENDER + G, the single-kernel model of gaston and rrBLUP
# for body mass index (Obesity.BMI) and body length (Obesity.BodyLength).
#
# References: asreml (all models), sommer::mmes (gblup_cage, ped_gen),
# gaston::lmm.aireml and rrBLUP::mixed.solve (gblup_only).
#
# Published values: Legarra et al. (2008, Materials and methods) quote the
# pedigree heritabilities of Valdar et al. (2006b): 0.21 for body mass index
# and 0.13 for body length. They come from a pedigree model on phenotypes
# pre-corrected for other fixed effects (body weight, season, month, day).
# They are reported as context, beside the genomic heritability of
# gblup_cage, and never counted.
#
#   Rscript published_data/scripts/04_mice_valdar2006.R
# ==============================================================================
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(.here, "..", "R", "common.R"))
suppressPackageStartupMessages({
  library(asreml); library(sommer); library(BGLR); library(gaston); library(rrBLUP)
})

start_checks("mice_valdar2006", "Valdar et al. 2006, doi:10.1038/ng1840")
LEGARRA <- "pedigree h2 of Valdar et al. 2006b quoted by Legarra et al. 2008, doi:10.1534/genetics.108.088575"

data(mice)
stopifnot(identical(rownames(mice.X), as.character(mice.pheno$SUBJECT.NAME)),
          identical(rownames(mice.A), rownames(mice.X)))
G <- vanraden(mice.X)
A <- mice.A + diag(1e-4, nrow(mice.A))
dimnames(A) <- dimnames(mice.A)
ids <- rownames(G)
d0 <- data.frame(id = factor(ids, levels = ids), idA = factor(ids, levels = ids),
                 GENDER = factor(mice.pheno$GENDER), Litter = factor(mice.pheno$Litter),
                 cage = factor(mice.pheno$cage),
                 BMI = mice.pheno$Obesity.BMI, BodyLength = mice.pheno$Obesity.BodyLength)
cat(sprintf("mice %d, SNPs %d, cages %d\n", nrow(d0), ncol(mice.X), nlevels(d0$cage)))
G_as <- as_asreml_G(G); A_as <- as_asreml_G(A)
WS <- "4gb"

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
so_cage <- "sommer::vsm( sommer::ism( cage ) )"
TIMES <- list()
# The remlax solver stays alive between fits (a JAX server holding about 2 GB
# here). It is stopped after each fit so that asreml, with its 4 GB workspace,
# fits beside it under the 7 GB cap of the run.
free_rx <- function() invisible(try(rx_server_stop(), silent = TRUE))
H2_CONTEXT <- c(BMI = 0.21, BodyLength = 0.13)

for (tr in c("BMI", "BodyLength")) {
  d <- d0[!is.na(d0[[tr]]), ]
  d$y <- d[[tr]]
  cat(sprintf("\n================ trait %s: %d records ================\n", tr, nrow(d)))

  # --- gblup_cage -------------------------------------------------------------
  m <- paste0(tr, " gblup_cage")
  cat("\n===", m, "===\n")
  f <- timed(rx_cpu(y ~ GENDER + Litter, random = ~ iid(cage) + vm(id, K = G), data = d,
                    vpredict = c(h2 = "V2/(V1+V2+V3)")))
  print(f); free_rx()
  a <- timed(stab(asreml(y ~ GENDER + Litter, random = ~ cage + vm(id, G_as), data = d,
                         trace = FALSE, maxit = 100, workspace = WS)))
  vc <- summary(a)$varcomp; print(vc)
  so <- timed(mmes(y ~ GENDER + Litter, random = ~ cage + vsm(ism(id), Gu = G),
                   rcov = ~ units, data = d, verbose = FALSE, dateWarning = FALSE, nIters = 200,
                   tolParConvLL = 1e-8, tolParConvNorm = 1e-8))
  TIMES[[m]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"), sommer = attr(so, "secs"))
  rx <- c(cage = f$sigmas$cage[1, 1], G = f$sigmas$id[1, 1], residual = f$sigma_res[1, 1])
  asr <- c(cage = vcr(vc, "cage"), G = vcr(vc, "vm(id, G_as)"), residual = vcr(vc, "units!R"))
  som <- c(cage = so$theta[[so_cage]][1, 1], G = so$theta[["vsm(ism(id), Gu = G)"]][1, 1],
           residual = so$theta[["units"]][1, 1])
  for (cn in names(rx)) {
    add_check(m, paste0("sigma2 ", cn), "asreml", rx[[cn]], asr[[cn]])
    add_check(m, paste0("sigma2 ", cn), "sommer", rx[[cn]], som[[cn]])
  }
  add_check(m, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik, tol = 1e-5, metric = "abs")
  u_rx <- setNames(f$blups$id[, 1], rownames(f$blups$id))
  add_blup_checks(m, "asreml", u_rx, asr_blup(a, "vm(id, G_as)_")[names(u_rx)])
  u_so <- so$uList[["vsm(ism(id), Gu = G)"]]
  add_blup_checks(m, "sommer", u_rx, setNames(u_so[, 1], rownames(u_so))[names(u_rx)])
  add_check(m, "h2 (genomic, cage in the denominator)", "Valdar 2006b",
            f$vpredict$predictions$value[1], H2_CONTEXT[[tr]], kind = "context", note = LEGARRA)

  # --- ped_gen: pedigree and genomic kernels together --------------------------
  m <- paste0(tr, " ped_gen")
  cat("\n===", m, "===\n")
  f2 <- timed(rx_cpu(y ~ GENDER + Litter, random = ~ iid(cage) + vm(idA, K = A) + vm(id, K = G),
                     data = d, n_restarts = 3))
  print(f2); free_rx()
  a2 <- timed(stab(asreml(y ~ GENDER + Litter, random = ~ cage + vm(idA, A_as) + vm(id, G_as),
                          data = d, trace = FALSE, maxit = 100, workspace = WS)))
  vc2 <- summary(a2)$varcomp; print(vc2)
  so2 <- timed(mmes(y ~ GENDER + Litter,
                    random = ~ cage + vsm(ism(idA), Gu = A) + vsm(ism(id), Gu = G),
                    rcov = ~ units, data = d, verbose = FALSE, dateWarning = FALSE, nIters = 200,
                   tolParConvLL = 1e-8, tolParConvNorm = 1e-8))
  TIMES[[m]] <- c(remlax = attr(f2, "secs"), asreml = attr(a2, "secs"), sommer = attr(so2, "secs"))
  rx <- c(cage = f2$sigmas$cage[1, 1], A = f2$sigmas$idA[1, 1], G = f2$sigmas$id[1, 1],
          residual = f2$sigma_res[1, 1])
  asr <- c(cage = vcr(vc2, "cage"), A = vcr(vc2, "vm(idA, A_as)"), G = vcr(vc2, "vm(id, G_as)"),
           residual = vcr(vc2, "units!R"))
  som <- c(cage = so2$theta[[so_cage]][1, 1], A = so2$theta[["vsm(ism(idA), Gu = A)"]][1, 1],
           G = so2$theta[["vsm(ism(id), Gu = G)"]][1, 1], residual = so2$theta[["units"]][1, 1])
  for (cn in names(rx)) {
    add_check(m, paste0("sigma2 ", cn), "asreml", rx[[cn]], asr[[cn]])
    add_check(m, paste0("sigma2 ", cn), "sommer", rx[[cn]], som[[cn]])
  }
  add_check(m, "logLik (asreml convention)", "asreml", f2$logLik_asreml, a2$loglik, tol = 1e-5, metric = "abs")
  u_rx <- setNames(f2$blups$id[, 1], rownames(f2$blups$id))
  add_blup_checks(m, "asreml", u_rx, asr_blup(a2, "vm(id, G_as)_")[names(u_rx)])

  # --- gblup_only: one kernel, the model of gaston and rrBLUP ------------------
  m <- paste0(tr, " gblup_only")
  cat("\n===", m, "===\n")
  f3 <- timed(rx_cpu(y ~ GENDER, random = ~ vm(id, K = G), data = d))
  print(f3); free_rx()
  a3 <- timed(stab(asreml(y ~ GENDER, random = ~ vm(id, G_as), data = d,
                          trace = FALSE, maxit = 100, workspace = WS)))
  vc3 <- summary(a3)$varcomp; print(vc3)
  Gd <- G[as.character(d$id), as.character(d$id)]
  X <- model.matrix(~ GENDER, d)
  ga <- timed(lmm.aireml(d$y, X = X, K = Gd, verbose = FALSE, eps = 1e-10, max_iter = 200L))
  rr <- timed(mixed.solve(d$y, X = X, K = Gd, method = "REML"))
  TIMES[[m]] <- c(remlax = attr(f3, "secs"), asreml = attr(a3, "secs"),
                  gaston = attr(ga, "secs"), rrBLUP = attr(rr, "secs"))
  rx <- c(G = f3$sigmas$id[1, 1], residual = f3$sigma_res[1, 1])
  add_check(m, "sigma2 G", "asreml", rx[["G"]], vcr(vc3, "vm(id, G_as)"))
  add_check(m, "sigma2 residual", "asreml", rx[["residual"]], vcr(vc3, "units!R"))
  add_check(m, "sigma2 G", "gaston", rx[["G"]], ga$tau)
  add_check(m, "sigma2 residual", "gaston", rx[["residual"]], ga$sigma2)
  add_check(m, "sigma2 G", "rrBLUP", rx[["G"]], rr$Vu)
  add_check(m, "sigma2 residual", "rrBLUP", rx[["residual"]], rr$Ve)
  add_check(m, "logLik (asreml convention)", "asreml", f3$logLik_asreml, a3$loglik, tol = 1e-5, metric = "abs")
  u_rx <- setNames(f3$blups$id[, 1], rownames(f3$blups$id))
  add_blup_checks(m, "asreml", u_rx, asr_blup(a3, "vm(id, G_as)_")[names(u_rx)])
  add_blup_checks(m, "gaston", u_rx, setNames(ga$BLUP_omega, as.character(d$id))[names(u_rx)])
  add_blup_checks(m, "rrBLUP", u_rx, setNames(as.numeric(rr$u), as.character(d$id))[names(u_rx)])
}

cat("\nElapsed seconds (remlax includes process start-up):\n")
for (m in names(TIMES)) cat(sprintf("  %-26s %s\n", m,
  paste(sprintf("%s %.1f", names(TIMES[[m]]), TIMES[[m]]), collapse = " | ")))
write_checks()
