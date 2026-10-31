# ==============================================================================
# 07_drops_gwas.R -- GWAS on the DROPS maize panel: remlax scan against
#                    statgenGWAS, gaston and asreml
# ------------------------------------------------------------------------------
# DATA. The DROPS panel (Millet et al. 2016, Plant Physiol 172:749-764,
# doi:10.1104/pp.16.00621; Millet et al. 2019, Sci Data 6:190011,
# doi:10.1038/sdata.2019.11): 246 maize hybrids (dent lines x common flint
# tester), 41,722 SNPs of the 50K Illumina and 600K Affymetrix arrays, one
# genotypic mean (BLUE) per hybrid and trial. It ships with statgenGWAS
# (dropsPheno, dropsMarkers, dropsMap) and its vignette runs single-trait GWAS
# on the trial Mur13W, the one used here.
#
# MODEL. One observation per hybrid, y = 1 mu + g + e, g ~ N(0, s2_g K),
# e ~ N(0, s2_e I). This is the EMMAX / P3D model, the only one where an
# external reference for the scan exists (see validation/etude_scan_vs_rrblup.R
# and docs/scan of remlax). K is the VanRaden kinship computed by statgenGWAS,
# plus a ridge of 1e-4 on the diagonal: the centred VanRaden matrix has rank
# q - 1 and asreml refuses it as a G-structure. The SAME matrix is passed to
# every program, so the comparison tests the solvers and not the builders.
#
# WHAT EACH PROGRAM COMPUTES, AND HOW THE CONVENTIONS ARE RECONCILED.
#   remlax  rx_scan: beta_j = num_j/den_j, var(beta_j) = 1/den_j and
#           chi2_j = num_j^2/den_j, at the scale of the null model (fixed V).
#   statgenGWAS (fastGLS) and rrBLUP: same GLS beta, but the test is an
#           F(1, n - 2) with the residual scale re-estimated at each SNP on
#           the model that contains it. The exact relation (vignette scan,
#           section 5) is
#               s2_j  = (y'Py - num_j^2/den_j) / (n - 2)
#               F_j   = chi2_j / s2_j
#           where y'Py is computed once with the same V. statgenGWAS reports
#           effectSe at the scale of the null model, sqrt(1/den_j) with its V
#           (measured: ratio 1 to the remlax se_dir), and only its p-value uses
#           s2_j.
#   gaston  association.test(test = "score"): chi2 at fixed V, the remlax
#           convention, with gaston's own AI-REML components. The statistic is
#           computed by gg_GWAS_lmm_score_f; its output differs by about 2e-6
#           relative from the same statistic recomputed in double precision
#           from gaston's own P and Py. The counted check is against that
#           double-precision recomputation; the native output is context.
#           association.test(test = "wald"): one REML refit per SNP, the
#           counterpart of rx_scan_exact().
#
# TWO LEVELS OF COMPARISON.
#   (A) components imposed: the scan of remlax is evaluated at the components
#       of the other program (theta = log(sd), maxiter 0). Any gap is then an
#       algebra gap. Tolerance 1e-6: rx_scan writes its table with 10
#       significant digits, and K has a condition number of about 3e5.
#   (B) each program with its own REML: the component gap with asreml and
#       gaston is a counted check (1e-3, the campaign tolerance); the induced
#       p-value gap is reported as context. statgenGWAS EMMA stops short of
#       the REML optimum (its delta comes from a grid refined by a 1-D search),
#       so its component gap is context and the counted check is that the
#       REML log-likelihood of remlax, evaluated by remlax at both points,
#       is not lower than at the statgenGWAS components. The top SNP and the set of SNPs above the Bonferroni threshold
#       are counted: they must not depend on the program.
#
#   Rscript published_data/scripts/07_drops_gwas.R
# ==============================================================================
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(f[1]), "..", "R", "common.R"))
})
suppressPackageStartupMessages({
  library(statgenGWAS); library(gaston); library(asreml)
})
options(gaston.verbose = FALSE)
RcppParallel::setThreadOptions(numThreads = 4)

start_checks("maize_drops_gwas",
             "Millet et al. 2016 Plant Physiol; Millet et al. 2019 Sci Data")

TRIAL  <- "Mur13W"
TRAITS <- c("anthesis", "grain.yield")
RIDGE  <- 1e-4
MAF    <- 0.05

# ---- data --------------------------------------------------------------------
data(dropsMarkers); data(dropsMap); data(dropsPheno)
M <- as.matrix(dropsMarkers[, -1])
rownames(M) <- dropsMarkers$Ind
storage.mode(M) <- "double"
map0 <- data.frame(chr = dropsMap$Chromosome, pos = dropsMap$Position,
                   row.names = dropsMap$SNP.names)
map0 <- map0[colnames(M), ]
ph <- dropsPheno[dropsPheno$Experiment == TRIAL, ]
ph <- data.frame(genotype = as.character(ph$Variety_ID), ph[, TRAITS])
ph <- ph[match(rownames(M), ph$genotype), ]
stopifnot(identical(ph$genotype, rownames(M)))
q <- nrow(M)

# SNPs kept: MAF >= 0.05 on these 246 hybrids, as statgenGWAS computes it
# (allele frequency = mean dose / 2). Every program then sees the same columns.
frq <- colMeans(M) / 2
M <- M[, pmin(frq, 1 - frq) >= MAF]
map0 <- map0[colnames(M), ]
cat(sprintf("trial %s: %d hybrids, %d SNPs after MAF >= %.2f\n", TRIAL, q, ncol(M), MAF))

kin_of <- function(Mx) {
  K <- statgenGWAS::kinship(Mx, method = "vanRaden")
  K <- K + diag(RIDGE, nrow(K))
  dimnames(K) <- list(rownames(Mx), rownames(Mx))
  K
}
K <- kin_of(M)
cat(sprintf("K: condition number %.3g\n", kappa(K, exact = TRUE)))

map_rx <- data.frame(snp = colnames(M), chr = as.character(map0$chr), pos = map0$pos)
gd <- createGData(geno = M, map = map0, pheno = list(ph = ph), kin = K)

# rx_model for a kinship K and optional extra fixed columns: used by the scan
# with imposed components, by rx_scan_exact() and by rx_scan_loco().
Zid <- Matrix::Diagonal(q)
make_build <- function(y) function(K, X_sup = NULL) {
  X <- cbind(intercept = rep(1, length(y)), X_sup)
  rx_model(y, X, list(rx_term("gid", Zid, K = K, t = 1L, levels = rownames(K))),
           rx_residual(), name = "y")
}

# y'Py for V = vg K + ve I and X = 1 (dimension q, computed once).
yPy_of <- function(y, K, vg, ve) {
  Vi <- solve(vg * K + diag(ve, length(y)))
  one <- rep(1, length(y))
  Vi1 <- Vi %*% one
  Py <- Vi %*% y - Vi1 * as.numeric(crossprod(Vi1, y) / crossprod(one, Vi1))
  sum(y * Py)
}
# statgenGWAS quantities rebuilt from the raw scan outputs.
to_statgen <- function(r, yPy, n) {
  s2 <- (yPy - r$num_dir^2 / r$den_dir_dir) / (n - 2)
  Fst <- r$chi2_dir / s2
  data.frame(snp = r$snp, effect = r$beta_dir, se = sqrt(s2 / r$den_dir_dir),
             mlp = -stats::pf(Fst, 1, n - 2, lower.tail = FALSE, log.p = TRUE) / log(10))
}
mlp_chi2 <- function(chi2) -stats::pchisq(chi2, 1, lower.tail = FALSE, log.p = TRUE) / log(10)
lambda_of <- function(mlp) stats::median(stats::qchisq(10^-mlp, 1, lower.tail = FALSE)) /
  stats::qchisq(0.5, 1)
scan_rx <- function(fit, model, Mx, mp)
  rx_scan(fit, model, markers = Mx, map = mp, incidences = c(dir = "gid"),
          tests = "dir", maf = NULL, coding = "012", backend = "cpu", verbose = FALSE)

# gaston needs a bed.matrix.
bm <- as.bed.matrix(x = M, fam = data.frame(famid = rownames(M), id = rownames(M),
                                           father = 0, mother = 0, sex = 0, pheno = NA),
                    bim = data.frame(chr = map0$chr, id = colnames(M), dist = 0,
                                     pos = map0$pos, A1 = "A", A2 = "B"))
eigK <- eigen(K, symmetric = TRUE)

timing <- list()
plots <- list()

for (tr in TRAITS) {
  y <- ph[[tr]]
  n <- length(y)
  mdl <- paste0(TRIAL, " ", tr)
  cat(sprintf("\n==== %s ====\n", mdl))
  build <- make_build(y)
  m0 <- build(K)
  d <- data.frame(gid = factor(rownames(M), levels = rownames(K)), y = y)

  # ---- null model: remlax, asreml, gaston, statgenGWAS ------------------------
  t_fit <- system.time(f0 <- rx_fit(m0, backend = "cpu", verbose = FALSE))[["elapsed"]]
  vg_rx <- f0$sigmas[[1]][1, 1]; ve_rx <- f0$sigma_res[1, 1]
  cat(sprintf("remlax: s2_g %.10g s2_e %.10g logLik_asreml %.10f (%.1f s)\n",
              vg_rx, ve_rx, f0$logLik_asreml, t_fit))

  a <- stab(asreml(y ~ 1, random = ~ vm(gid, as_asreml_G(K)), data = d,
                   trace = FALSE, workspace = "2gb", maxit = 100))
  vc <- summary(a)$varcomp
  add_check(mdl, "null s2_g", "asreml", vg_rx, vc[1, "component"])
  add_check(mdl, "null s2_e", "asreml", ve_rx, vc["units!R", "component"])
  add_check(mdl, "null logLik", "asreml", f0$logLik_asreml, a$loglik, tol = 1e-5, metric = "abs")

  t_ga0 <- system.time(
    ga <- lmm.aireml(y, X = matrix(1, n), K = K, eps = 1e-10, max_iter = 200, verbose = FALSE)
  )[["elapsed"]]
  add_check(mdl, "null s2_g", "gaston", vg_rx, ga$tau)
  add_check(mdl, "null s2_e", "gaston", ve_rx, ga$sigma2)

  t_sg <- system.time(
    sg <- runSingleTraitGwas(gd, trials = "ph", traits = tr, kin = K,
                             GLSMethod = "single", remlAlgo = "EMMA",
                             useMAF = TRUE, MAF = 0, MAC = 0,
                             thrType = "bonf", alpha = 0.05)
  )[["elapsed"]]
  vcs <- sg$GWASInfo$varComp$ph[[tr]]
  vg_sg <- unname(vcs[1]); ve_sg <- unname(vcs[2])
  cat(sprintf("statgenGWAS EMMA: s2_g %.10g s2_e %.10g (%.1f s)\n", vg_sg, ve_sg, t_sg))
  f_at_sg <- rx_fit(m0, backend = "cpu", verbose = FALSE, maxiter = 0, polish = 0,
                    theta_init = 0.5 * log(c(vg_sg, ve_sg)), hessian = FALSE)
  add_check(mdl, "null s2_g", "statgenGWAS", vg_rx, vg_sg, kind = "context",
            note = "statgenGWAS EMMA does not reach the REML optimum, see logLik row")
  add_check(mdl, "null s2_e", "statgenGWAS", ve_rx, ve_sg, kind = "context",
            note = "statgenGWAS EMMA does not reach the REML optimum, see logLik row")
  add_check(mdl, "REML logLik gain over statgen comps", "statgenGWAS",
            f0$logLik - f_at_sg$logLik, NA, tol = 0, metric = "ge",
            note = sprintf("remlax REML logLik %.8f at its optimum, %.8f at the statgenGWAS components",
                           f0$logLik, f_at_sg$logLik))
  gw <- as.data.frame(sg$GWAResult$ph)
  gw <- gw[match(colnames(M), gw$snp), ]
  stopifnot(identical(gw$snp, colnames(M)))
  gw_mlp <- -log10(gw$pValue)

  # ---- (A) statgenGWAS components imposed on the remlax scan -------------------
  f_sg <- list(theta = 0.5 * log(c(vg_sg, ve_sg)))
  r_sg <- scan_rx(f_sg, m0, M, map_rx)
  rb <- to_statgen(r_sg, yPy_of(y, K, vg_sg, ve_sg), n)
  sdb <- stats::sd(gw$effect)
  add_check(mdl, "effect, statgen V (max|d|/sd)", "statgenGWAS",
            max(abs(rb$effect - gw$effect)) / sdb, NA, tol = 1e-6, metric = "value")
  add_check(mdl, "effectSe, statgen V (max rel)", "statgenGWAS",
            max(abs(r_sg$se_dir - gw$effectSe) / gw$effectSe), NA, tol = 1e-6, metric = "value",
            note = "statgenGWAS reports the SE at the null-model scale, as remlax")
  add_check(mdl, "-log10 p, statgen V (max|d|)", "statgenGWAS",
            max(abs(rb$mlp - gw_mlp)), NA, tol = 1e-6, metric = "value",
            note = "remlax chi2 converted to the F(1, n-2) of statgenGWAS")

  # ---- (A) gaston components imposed: score test at fixed V --------------------
  t_gs <- system.time(
    gs <- association.test(bm, Y = y, X = matrix(1, n), method = "lmm",
                           response = "quantitative", test = "score", K = K,
                           eigenK = eigK, eps = 1e-10, max_iter = 200)
  )[["elapsed"]]
  # The same score in double precision, from gaston's own P and Py.
  gaP <- lmm.aireml(y, X = matrix(1, n), K = K, eps = 1e-10, max_iter = 200,
                    verbose = FALSE, get.P = TRUE)
  Mc <- sweep(M, 2, colMeans(M))
  T_dbl <- as.numeric(crossprod(Mc, gaP$Py))^2 / colSums(Mc * (gaP$P %*% Mc))
  rm(Mc)
  # The null that association.test fits internally, on trans.X(X).
  gaI <- lmm.aireml(y, X = gaston:::trans.X(matrix(1, n), mean.y = mean(y)), K = K, eps = 1e-10,
                   max_iter = 200, verbose = FALSE)
  f_ga <- list(theta = 0.5 * log(c(gaP$tau, gaP$sigma2)))
  r_ga <- scan_rx(f_ga, m0, M, map_rx)
  add_check(mdl, "score chi2, gaston V (max rel)", "gaston",
            max(abs(r_ga$chi2_dir - T_dbl) / pmax(T_dbl, 1)), NA, tol = 1e-6,
            metric = "value", note = "score recomputed in double from gaston's P and Py")
  add_check(mdl, "score chi2 native output (max rel)", "gaston",
            max(abs(r_ga$chi2_dir - gs$score) / pmax(gs$score, 1)), NA, kind = "context",
            metric = "value",
            note = sprintf("gaston native output vs its own double recomputation: %.2e",
                           max(abs(gs$score - T_dbl) / pmax(T_dbl, 1))))

  # ---- (B) each program with its own REML -------------------------------------
  t_scan <- system.time(r_rx <- scan_rx(f0, m0, M, map_rx))[["elapsed"]]
  rx_F <- to_statgen(r_rx, yPy_of(y, K, vg_rx, ve_rx), n)
  rx_mlp <- mlp_chi2(r_rx$chi2_dir)
  add_check(mdl, "-log10 p own REML (max|d|)", "statgenGWAS",
            max(abs(rx_F$mlp - gw_mlp)), NA, kind = "context", metric = "value",
            note = "remlax REML vs statgenGWAS EMMA components, F convention on both sides")
  add_check(mdl, "-log10 p own REML (max|d|)", "gaston",
            max(abs(rx_mlp + log10(gs$p))), NA, kind = "context", metric = "value",
            note = "remlax REML vs gaston AI-REML components, chi2 at fixed V on both sides")
  thr <- -log10(0.05 / ncol(M))
  top_rx <- colnames(M)[which.max(rx_F$mlp)]
  top_sg <- colnames(M)[which.max(gw_mlp)]
  # statgenGWAS EMMA stops short of the REML optimum, so a SNP near the
  # threshold can change side: the top SNP and the set are context there and
  # counted against gaston, whose AI-REML reaches the optimum of remlax.
  add_check(mdl, "top SNP identical", "statgenGWAS", as.numeric(top_rx == top_sg), 1,
            kind = "context", metric = "abs",
            note = paste("remlax", top_rx, "| statgenGWAS", top_sg))
  top_gs <- colnames(M)[which.max(T_dbl)]
  add_check(mdl, "top SNP identical", "gaston", as.numeric(colnames(M)[which.max(rx_mlp)] == top_gs),
            1, tol = 0, metric = "abs",
            note = paste("remlax", colnames(M)[which.max(rx_mlp)], "| gaston", top_gs))
  set_rx <- colnames(M)[rx_F$mlp > thr]; set_sg <- colnames(M)[gw_mlp > thr]
  flip <- union(setdiff(set_rx, set_sg), setdiff(set_sg, set_rx))
  fj <- match(flip, colnames(M))
  add_check(mdl, "Bonferroni set symmetric diff", "statgenGWAS", length(flip), 0,
            kind = "context", metric = "abs",
            note = paste0(sprintf("%d SNPs above %.2f in remlax (F), %d in statgenGWAS",
                                  length(set_rx), thr, length(set_sg)),
                          if (length(flip)) paste0("; differing: ", paste(sprintf(
                            "%s remlax %.3f statgen %.3f", flip, rx_F$mlp[fj], gw_mlp[fj]),
                            collapse = ", ")) else ""))
  # gaston's own null (lmm.aireml on X, double precision). association.test
  # refits the null on trans.X(X), and on flowering time in the rice panel
  # that fit ends at another point with logL = Inf; its native output is
  # therefore context.
  set_chi <- colnames(M)[rx_mlp > thr]; set_gs <- colnames(M)[mlp_chi2(T_dbl) > thr]
  set_nat <- colnames(M)[-log10(gs$p) > thr]
  add_check(mdl, "Bonferroni set symmetric diff", "gaston",
            length(union(setdiff(set_chi, set_gs), setdiff(set_gs, set_chi))), 0,
            tol = 0, metric = "abs",
            note = sprintf("%d SNPs in remlax (chi2 fixed V), %d with gaston's null (double), %d in association.test native output",
                           length(set_chi), length(set_gs), length(set_nat)))
  add_check(mdl, "Bonferroni set, native output", "gaston",
            length(union(setdiff(set_chi, set_nat), setdiff(set_nat, set_chi))), 0,
            kind = "context", metric = "abs",
            note = sprintf("association.test null refit on trans.X(X): tau %.6g sigma2 %.6g logL %s",
                           gaI$tau, gaI$sigma2, format(gaI$logL)))
  add_check(mdl, "lambda GC (F convention)", "statgenGWAS", lambda_of(rx_F$mlp),
            lambda_of(gw_mlp), kind = "context")
  add_check(mdl, "lambda GC (chi2 fixed V)", "remlax", attr(r_rx, "meta")$lambda_gc$dir,
            NA, kind = "context", note = "lambda reported by rx_scan")

  # ---- exact refit of the top SNPs against gaston's per-SNP Wald -------------
  t_ex <- system.time(
    ex <- rx_scan_exact(r_rx, "dir", build = build, K_of = function(Mx) K,
                        markers = M, map = map_rx, incidences = c(dir = "gid"),
                        top = 10L, loco = FALSE, fit_args = list(backend = "cpu"),
                        verbose = FALSE)
  )[["elapsed"]]
  jj <- match(ex$snp, colnames(M))
  gw_w <- association.test(bm[, jj], Y = y, X = matrix(1, n), method = "lmm",
                           response = "quantitative", test = "wald", K = K,
                           eigenK = eigK, tol = 1e-10)
  add_check(mdl, "top-10 refit beta (max rel)", "gaston",
            max(abs(ex$beta_dir - gw_w$beta) / abs(gw_w$beta)), NA, tol = 1e-3,
            metric = "value", note = "rx_scan_exact vs gaston Wald, one REML per SNP")
  add_check(mdl, "top-10 refit -log10 p (max|d|)", "gaston",
            max(abs(-log10(ex$p_exact) + log10(gw_w$p))), NA, tol = 1e-3, metric = "value",
            note = "rx_scan_exact vs gaston Wald, one REML per SNP")
  add_check(mdl, "top-10 refit all converged", "remlax", mean(ex$converged), NA,
            tol = 1, metric = "ge")

  timing[[tr]] <- data.frame(trait = tr, remlax_null_fit_s = t_fit, remlax_scan_s = t_scan,
                             remlax_exact_top10_s = t_ex, statgenGWAS_s = t_sg,
                             gaston_aireml_s = t_ga0, gaston_score_s = t_gs)
  plots[[tr]] <- list(chr = map_rx$chr, pos = map_rx$pos, rx = rx_F$mlp, sg = gw_mlp, thr = thr)
}

# ---- LOCO on anthesis: chromosome-specific kinships --------------------------
tr <- "anthesis"; y <- ph[[tr]]; n <- length(y); build <- make_build(y)
mdl <- paste0(TRIAL, " ", tr, " LOCO")
cat(sprintf("\n==== %s ====\n", mdl))
chrs <- unique(map_rx$chr)
Kloco <- lapply(stats::setNames(chrs, chrs), function(cc) kin_of(M[, map_rx$chr != cc]))
t_sgl <- system.time(
  sgl <- runSingleTraitGwas(gd, trials = "ph", traits = tr, kin = unname(Kloco),
                            GLSMethod = "multi", remlAlgo = "EMMA",
                            useMAF = TRUE, MAF = 0, MAC = 0, thrType = "bonf")
)[["elapsed"]]
gwl <- as.data.frame(sgl$GWAResult$ph)
gwl <- gwl[match(colnames(M), gwl$snp), ]
vcl <- sgl$GWASInfo$varComp$ph[[tr]]
cat("statgenGWAS varComp per chromosome:\n"); print(str(vcl))

# (A) statgenGWAS chromosome-specific components imposed on remlax.
rb_l <- list()
for (k in seq_along(chrs)) {
  cc <- chrs[k]; inside <- map_rx$chr == cc
  vk <- unlist(vcl[[k]])
  r_k <- scan_rx(list(theta = 0.5 * log(vk[1:2])), build(Kloco[[cc]]),
                 M[, inside, drop = FALSE], map_rx[inside, ])
  rb_l[[cc]] <- to_statgen(r_k, yPy_of(y, Kloco[[cc]], vk[1], vk[2]), n)
}
rb_l <- do.call(rbind, rb_l)
rb_l <- rb_l[match(colnames(M), rb_l$snp), ]
add_check(mdl, "effect, statgen V (max|d|/sd)", "statgenGWAS",
          max(abs(rb_l$effect - gwl$effect)) / stats::sd(gwl$effect), NA, tol = 1e-6,
          metric = "value")
add_check(mdl, "-log10 p, statgen V (max|d|)", "statgenGWAS",
          max(abs(rb_l$mlp + log10(gwl$pValue))), NA, tol = 1e-6, metric = "value")

# (B) rx_scan_loco with its own REML per chromosome.
t_loco <- system.time(
  rl <- rx_scan_loco(build, K_of = kin_of, markers = M, map = map_rx,
                     incidences = c(dir = "gid"), tests = "dir", maf = NULL, coding = "012",
                     fit_args = list(backend = "cpu"), verbose = FALSE)
)[["elapsed"]]
rl <- rl[match(colnames(M), rl$snp), ]
loco <- attr(rl, "loco")
gap_vc <- 0; gain_min <- Inf; mlp_l <- numeric(nrow(rl))
for (k in seq_along(chrs)) {
  cc <- chrs[k]; inside <- map_rx$chr == cc
  th <- loco$theta[[match(cc, loco$chr)]]
  vk <- exp(2 * th)
  vs <- unlist(vcl[[k]])[1:2]
  gap_vc <- max(gap_vc, relc(vk[1], vs[1]), relc(vk[2], vs[2]))
  mk <- build(Kloco[[cc]])
  l_rx <- rx_fit(mk, backend = "cpu", verbose = FALSE, maxiter = 0, polish = 0,
                 theta_init = th, hessian = FALSE)$logLik
  l_sg <- rx_fit(mk, backend = "cpu", verbose = FALSE, maxiter = 0, polish = 0,
                 theta_init = 0.5 * log(vs), hessian = FALSE)$logLik
  gain_min <- min(gain_min, l_rx - l_sg)
  mlp_l[inside] <- to_statgen(rl[inside, ], yPy_of(y, Kloco[[cc]], vk[1], vk[2]), n)$mlp
}
add_check(mdl, "per-chromosome components (max rel)", "statgenGWAS", gap_vc, NA,
          kind = "context", metric = "value", note = "rx_scan_loco REML vs statgenGWAS EMMA")
add_check(mdl, "min REML logLik gain over statgen", "statgenGWAS", gain_min, NA,
          tol = 0, metric = "ge",
          note = "smallest over chromosomes of logLik(remlax theta) - logLik(statgenGWAS components)")
add_check(mdl, "all chromosome nulls converged", "remlax", mean(loco$converged), NA,
          tol = 1, metric = "ge")
add_check(mdl, "-log10 p own REML (max|d|)", "statgenGWAS",
          max(abs(mlp_l + log10(gwl$pValue))), NA, kind = "context", metric = "value")
thr <- -log10(0.05 / ncol(M))
s1 <- colnames(M)[mlp_l > thr]; s2 <- colnames(M)[-log10(gwl$pValue) > thr]
add_check(mdl, "Bonferroni set symmetric diff", "statgenGWAS",
          length(union(setdiff(s1, s2), setdiff(s2, s1))), 0, kind = "context", metric = "abs",
          note = sprintf("%d SNPs in remlax LOCO, %d in statgenGWAS multi", length(s1), length(s2)))
timing[["loco"]] <- data.frame(trait = "anthesis LOCO", remlax_null_fit_s = NA,
                               remlax_scan_s = t_loco, remlax_exact_top10_s = NA,
                               statgenGWAS_s = t_sgl, gaston_aireml_s = NA, gaston_score_s = NA)
plots[["anthesis LOCO"]] <- list(chr = map_rx$chr, pos = map_rx$pos, rx = mlp_l,
                                 sg = -log10(gwl$pValue), thr = thr)

# ---- outputs -------------------------------------------------------------------
tim <- do.call(rbind, timing)
cat("\nTimings (elapsed seconds):\n"); print(tim, row.names = FALSE)
utils::write.csv(tim, file.path(pd_root(), "results", "maize_drops_gwas_timing.csv"),
                 row.names = FALSE)

png(file.path(pd_root(), "results", "maize_drops_gwas_manhattan.png"),
    width = 2400, height = 600 * length(plots), res = 200)
op <- par(mfrow = c(length(plots), 1), mar = c(3.5, 4.2, 2, 1))
for (nm in names(plots)) {
  pl <- plots[[nm]]
  o <- order(as.numeric(pl$chr), pl$pos)
  x <- seq_along(o)
  col <- ifelse(as.numeric(pl$chr[o]) %% 2 == 0, "grey55", "grey20")
  plot(x, pl$sg[o], pch = 16, cex = 0.5, col = col, xaxt = "n", xlab = "chromosome",
       ylab = expression(-log[10](p)), main = paste(TRIAL, nm),
       ylim = c(0, max(pl$sg, pl$rx, pl$thr) * 1.05))
  points(x, pl$rx[o], pch = 1, cex = 0.7, col = "#D55E00")
  abline(h = pl$thr, lty = 2)
  mids <- tapply(x, pl$chr[o], mean)
  axis(1, at = mids, labels = names(mids), tick = FALSE)
  legend("topleft", c("statgenGWAS", "remlax (F convention)"), pch = c(16, 1),
         col = c("grey20", "#D55E00"), bty = "n", cex = 0.8)
}
par(op); dev.off()

write_checks()
