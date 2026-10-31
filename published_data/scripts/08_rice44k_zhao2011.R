# ==============================================================================
# 08_rice44k_zhao2011.R -- GWAS on the rice diversity panel (44K chip):
#                          remlax scan against statgenGWAS, gaston and asreml
# ------------------------------------------------------------------------------
# DATA. Zhao et al. (2011) Genome-wide association mapping reveals a rich
# genetic architecture of complex traits in Oryza sativa. Nat Commun 2:467,
# doi:10.1038/ncomms1467. 413 inbred accessions from 82 countries, 36,901 SNPs
# (call rate > 70 %, MAF > 0.01), 34 traits. Public files from
# www.ricediversity.org/44kgwas, downloaded into data_cache/rice44k/:
#   RiceDiversity.44K.MSU6.Genotypes_PLINK.zip  (sativas413.ped/.map/.fam)
#   RiceDiversity_44K_Phenotypes_34traits_PLINK.txt
#
# MODEL. The model of the paper (Methods, "Genome-wide association"): EMMA
# mixed model, an IBS relatedness matrix as random effect, and "for all
# samples, SNPs and the top four PCs were used as fixed effects". Here:
#   y = X b + g + e,  X = [1, PC1..PC4],  g ~ N(0, s2_g K),  e ~ N(0, s2_e I)
# with K the IBS kinship of statgenGWAS plus a ridge of 1e-4 (the same matrix
# for every program). Missing genotypes (< 4.5 % in the paper) are imputed by
# the most frequent genotype of the SNP, so that the doses stay integers for
# gaston. SNPs with MAF < 0.05 among the phenotyped accessions are dropped.
#
# CONVENTIONS. The same as 07_drops_gwas.R: statgenGWAS re-estimates the
# residual scale per SNP and tests with an F(1, n - c - 1), c = ncol(X); the
# remlax chi2 at fixed V is converted with
#     s2_j = (y'Py - num_j^2/den_j) / (n - c - 1),   F_j = chi2_j / s2_j.
# statgenGWAS reports effectSe at the null-model scale, as remlax. gaston's
# score test is a chi2 at fixed V; the counted check is against the same
# statistic recomputed in double precision from gaston's P and Py, and
# gaston's native output (about 2e-6 away from it on DROPS) is context.
# statgenGWAS EMMA stops short of the REML optimum (07_drops_gwas.R): its
# component gap is context, and the counted check is that remlax's REML
# log-likelihood at its optimum is not lower than at the EMMA components.
#
# PUBLISHED ASSOCIATIONS (Zhao et al. 2011, Results, verified in the PMC full
# text PMC3195253): "the same known genes showed clear signal for the same
# phenotypes", with GS3 for grain length and Waxy for amylose content; SD1 for
# plant height, detected by "both the naive and the mixed model"; "the most
# significant signal was observed under very long days in Aberdeen around
# HD1" for flowering time. The paper gives no coordinates. The gene positions
# below are MSU Rice Genome Annotation loci (LOC_Os06g04200, LOC_Os03g29380,
# LOC_Os01g66100, LOC_Os06g16370), rounded to 10 kb, and are an input of this
# script, not a value of the paper. Whether the lead SNP falls within 1 Mb of
# the gene is reported as context: the model reproduces the paper's, but the
# software, the kinship builder and the imputation differ.
#
#   Rscript published_data/scripts/08_rice44k_zhao2011.R
# ==============================================================================
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(f[1]), "..", "R", "common.R"))
})
suppressPackageStartupMessages({
  library(statgenGWAS); library(gaston); library(asreml); library(data.table)
})
options(gaston.verbose = FALSE)
RcppParallel::setThreadOptions(numThreads = 4)

start_checks("rice44k_zhao2011", "Zhao et al. 2011 Nat Commun 2:467")

RIDGE <- 1e-4
MAF   <- 0.05
WIN   <- 1e6
TRAITS <- list(
  "Amylose content"            = list(gene = "Waxy", chr = "6", pos = 1.77e6),
  "Seed length"                = list(gene = "GS3",  chr = "3", pos = 16.73e6),
  "Plant height"               = list(gene = "SD1",  chr = "1", pos = 38.38e6),
  "Flowering time at Aberdeen" = list(gene = "HD1",  chr = "6", pos = 9.34e6))

# ---- data ----------------------------------------------------------------------
gdir <- pd_cache("rice44k", "geno", "RiceDiversity_44K_Genotypes_PLINK")
pfile <- pd_cache("rice44k", "RiceDiversity_44K_Phenotypes_34traits_PLINK.txt")
if (!file.exists(file.path(gdir, "sativas413.ped")))
  stop("download and unzip the ricediversity.org files into ", dirname(gdir))

ped <- fread(file.path(gdir, "sativas413.ped"), header = FALSE, sep = " ",
             colClasses = "character")
mp <- fread(file.path(gdir, "sativas413.map"), header = FALSE,
            col.names = c("chr", "snp", "cM", "pos"))
ids <- ped[[2]]
A <- as.matrix(ped[, -(1:6)])
rm(ped)
p0 <- nrow(mp)
stopifnot(ncol(A) == 2 * p0)
A1 <- A[, seq(1, 2 * p0, 2)]; A2 <- A[, seq(2, 2 * p0, 2)]
rm(A)
# Dose of the first allele observed at each SNP; "0" is a missing call.
M <- matrix(NA_real_, length(ids), p0, dimnames = list(ids, mp$snp))
for (j in seq_len(p0)) {
  a1 <- A1[, j]; a2 <- A2[, j]
  ok <- a1 != "0" & a2 != "0"
  if (!any(ok)) next
  ref <- a1[ok][1]
  M[ok, j] <- (a1[ok] == ref) + (a2[ok] == ref)
}
rm(A1, A2)
miss <- mean(is.na(M))
# Imputation by the most frequent genotype.
for (j in which(colSums(is.na(M)) > 0)) {
  tb <- tabulate(M[, j] + 1, nbins = 3)
  M[is.na(M[, j]), j] <- which.max(tb) - 1
}
cat(sprintf("%d accessions, %d SNPs, %.2f %% missing calls imputed\n",
            nrow(M), ncol(M), 100 * miss))

ph <- fread(cmd = paste("tr '\\r' '\\n' <", shQuote(pfile)), sep = "\t",
            na.strings = c("NA", ""), check.names = FALSE)
ph <- as.data.frame(ph)
ph$NSFTVID <- as.character(ph$NSFTVID)
ph <- ph[match(ids, ph$NSFTVID), ]
stopifnot(identical(ph$NSFTVID, ids))

# Kinship on the whole panel (IBS, as in the paper) and the top four PCs of the
# standardised doses (EIGENSOFT standardises them).
frq_all <- colMeans(M) / 2
Mall <- M[, pmin(frq_all, 1 - frq_all) >= 0.01]
Kall <- statgenGWAS::kinship(Mall, method = "IBS")
dimnames(Kall) <- list(ids, ids)
Ws <- scale(Mall)
pcs <- prcomp(Ws, center = FALSE, rank. = 4)$x
colnames(pcs) <- paste0("PC", 1:4)
rm(Ws)
map_all <- data.frame(snp = mp$snp, chr = as.character(mp$chr), pos = mp$pos,
                      row.names = mp$snp)

yPy_of <- function(y, X, K, vg, ve) {
  Vi <- solve(vg * K + diag(ve, length(y)))
  ViX <- Vi %*% X
  Py <- Vi %*% y - ViX %*% solve(crossprod(X, ViX), crossprod(ViX, y))
  sum(y * Py)
}
to_statgen <- function(r, yPy, n, c) {
  v2 <- n - c - 1
  s2 <- (yPy - r$num_dir^2 / r$den_dir_dir) / v2
  Fst <- r$chi2_dir / s2
  data.frame(snp = r$snp, effect = r$beta_dir, se = sqrt(s2 / r$den_dir_dir),
             mlp = -stats::pf(Fst, 1, v2, lower.tail = FALSE, log.p = TRUE) / log(10))
}
mlp_chi2 <- function(chi2) -stats::pchisq(chi2, 1, lower.tail = FALSE, log.p = TRUE) / log(10)
lambda_of <- function(mlp) stats::median(stats::qchisq(10^-mlp, 1, lower.tail = FALSE)) /
  stats::qchisq(0.5, 1)

timing <- list(); plots <- list()

for (tr in names(TRAITS)) {
  info <- TRAITS[[tr]]
  keep <- is.finite(ph[[tr]])
  y <- ph[[tr]][keep]; n <- length(y)
  idt <- ids[keep]
  Mt <- M[keep, ]
  ft <- colMeans(Mt) / 2
  Mt <- Mt[, pmin(ft, 1 - ft) >= MAF]
  mapt <- map_all[colnames(Mt), ]
  K <- Kall[keep, keep] + diag(RIDGE, n)
  X <- cbind(intercept = 1, pcs[keep, ])
  c_x <- ncol(X)
  mdl <- tr
  cat(sprintf("\n==== %s: %d accessions, %d SNPs ====\n", tr, n, ncol(Mt)))

  Zid <- Matrix::Diagonal(n)
  build <- function(K, X_sup = NULL)
    rx_model(y, cbind(X, X_sup), list(rx_term("gid", Zid, K = K, t = 1L, levels = idt)),
             rx_residual(), name = "y")
  m0 <- build(K)
  scan_rx <- function(fit)
    rx_scan(fit, m0, markers = Mt, map = mapt, incidences = c(dir = "gid"),
            tests = "dir", maf = NULL, coding = "012", backend = "cpu", verbose = FALSE)

  # ---- null model -------------------------------------------------------------
  t_fit <- system.time(f0 <- rx_fit(m0, backend = "cpu", verbose = FALSE))[["elapsed"]]
  vg_rx <- f0$sigmas[[1]][1, 1]; ve_rx <- f0$sigma_res[1, 1]
  cat(sprintf("remlax: s2_g %.10g s2_e %.10g logLik_asreml %.10f (%.1f s)\n",
              vg_rx, ve_rx, f0$logLik_asreml, t_fit))
  d <- data.frame(gid = factor(idt, levels = idt), y = y, pcs[keep, ])
  a <- stab(asreml(y ~ PC1 + PC2 + PC3 + PC4, random = ~ vm(gid, as_asreml_G(K)),
                   data = d, trace = FALSE, workspace = "2gb", maxit = 100))
  vc <- summary(a)$varcomp
  add_check(mdl, "null s2_g", "asreml", vg_rx, vc[1, "component"])
  add_check(mdl, "null s2_e", "asreml", ve_rx, vc["units!R", "component"])
  add_check(mdl, "null logLik", "asreml", f0$logLik_asreml, a$loglik, tol = 1e-5, metric = "abs")
  add_check(mdl, "heritability s2_g/(s2_g+s2_e)", "asreml", vg_rx / (vg_rx + ve_rx),
            vc[1, "component"] / (vc[1, "component"] + vc["units!R", "component"]))

  t_ga0 <- system.time(
    ga <- lmm.aireml(y, X = X, K = K, eps = 1e-10, max_iter = 200, verbose = FALSE)
  )[["elapsed"]]
  add_check(mdl, "null s2_g", "gaston", vg_rx, ga$tau)
  add_check(mdl, "null s2_e", "gaston", ve_rx, ga$sigma2)

  phs <- data.frame(genotype = idt, y = y)
  cvs <- data.frame(pcs[keep, ], row.names = idt)
  gd <- createGData(geno = Mt, map = data.frame(chr = mapt$chr, pos = mapt$pos,
                                                row.names = mapt$snp),
                    pheno = list(ph = phs), kin = K, covar = cvs)
  t_sg <- system.time(
    sg <- runSingleTraitGwas(gd, trials = "ph", traits = "y", covar = paste0("PC", 1:4),
                             kin = K, GLSMethod = "single", remlAlgo = "EMMA",
                             useMAF = TRUE, MAF = 0, MAC = 0,
                             thrType = "bonf", alpha = 0.05)
  )[["elapsed"]]
  vcs <- sg$GWASInfo$varComp$ph$y
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
  gw <- gw[match(colnames(Mt), gw$snp), ]
  stopifnot(identical(gw$snp, colnames(Mt)))
  gw_mlp <- -log10(gw$pValue)

  # ---- (A) components imposed ------------------------------------------------
  r_sg <- scan_rx(list(theta = 0.5 * log(c(vg_sg, ve_sg))))
  rb <- to_statgen(r_sg, yPy_of(y, X, K, vg_sg, ve_sg), n, c_x)
  add_check(mdl, "effect, statgen V (max|d|/sd)", "statgenGWAS",
            max(abs(rb$effect - gw$effect)) / stats::sd(gw$effect), NA, tol = 1e-6,
            metric = "value")
  add_check(mdl, "effectSe, statgen V (max rel)", "statgenGWAS",
            max(abs(r_sg$se_dir - gw$effectSe) / gw$effectSe), NA, tol = 1e-6, metric = "value",
            note = "statgenGWAS reports the SE at the null-model scale, as remlax")
  add_check(mdl, "-log10 p, statgen V (max|d|)", "statgenGWAS",
            max(abs(rb$mlp - gw_mlp)), NA, tol = 1e-6, metric = "value",
            note = "remlax chi2 converted to the F(1, n-c-1) of statgenGWAS")

  bm <- as.bed.matrix(x = Mt, fam = data.frame(famid = idt, id = idt, father = 0, mother = 0,
                                              sex = 0, pheno = NA),
                      bim = data.frame(chr = as.integer(mapt$chr), id = mapt$snp, dist = 0,
                                       pos = mapt$pos, A1 = "A", A2 = "B"))
  eigK <- eigen(K, symmetric = TRUE)
  t_gs <- system.time(
    gs <- association.test(bm, Y = y, X = X, method = "lmm", response = "quantitative",
                           test = "score", K = K, eigenK = eigK, eps = 1e-10, max_iter = 200)
  )[["elapsed"]]
  gaP <- lmm.aireml(y, X = X, K = K, eps = 1e-10, max_iter = 200, verbose = FALSE,
                    get.P = TRUE)
  Mc <- sweep(Mt, 2, colMeans(Mt))
  T_dbl <- as.numeric(crossprod(Mc, gaP$Py))^2 / colSums(Mc * (gaP$P %*% Mc))
  rm(Mc)
  # The null that association.test fits internally, on trans.X(X).
  gaI <- lmm.aireml(y, X = gaston:::trans.X(X, mean.y = mean(y)), K = K, eps = 1e-10,
                   max_iter = 200, verbose = FALSE)
  r_ga <- scan_rx(list(theta = 0.5 * log(c(gaP$tau, gaP$sigma2))))
  add_check(mdl, "score chi2, gaston V (max rel)", "gaston",
            max(abs(r_ga$chi2_dir - T_dbl) / pmax(T_dbl, 1)), NA, tol = 1e-6,
            metric = "value", note = "score recomputed in double from gaston's P and Py")
  add_check(mdl, "score chi2 native output (max rel)", "gaston",
            max(abs(r_ga$chi2_dir - gs$score) / pmax(gs$score, 1)), NA, kind = "context",
            metric = "value",
            note = sprintf("gaston native output vs its own double recomputation: %.2e",
                           max(abs(gs$score - T_dbl) / pmax(T_dbl, 1))))

  # ---- (B) own REML -----------------------------------------------------------
  t_scan <- system.time(r_rx <- scan_rx(f0))[["elapsed"]]
  rx_F <- to_statgen(r_rx, yPy_of(y, X, K, vg_rx, ve_rx), n, c_x)
  rx_mlp <- mlp_chi2(r_rx$chi2_dir)
  add_check(mdl, "-log10 p own REML (max|d|)", "statgenGWAS",
            max(abs(rx_F$mlp - gw_mlp)), NA, kind = "context", metric = "value",
            note = "remlax REML vs statgenGWAS EMMA components, F convention on both sides")
  add_check(mdl, "-log10 p own REML (max|d|)", "gaston",
            max(abs(rx_mlp + log10(gs$p))), NA, kind = "context", metric = "value",
            note = "remlax REML vs gaston AI-REML components, chi2 at fixed V on both sides")
  thr <- -log10(0.05 / ncol(Mt))
  top_rx <- colnames(Mt)[which.max(rx_F$mlp)]
  top_sg <- colnames(Mt)[which.max(gw_mlp)]
  # statgenGWAS EMMA stops short of the REML optimum, so a SNP near the
  # threshold can change side: the top SNP and the set are context there and
  # counted against gaston, whose AI-REML reaches the optimum of remlax.
  add_check(mdl, "top SNP identical", "statgenGWAS", as.numeric(top_rx == top_sg), 1,
            kind = "context", metric = "abs",
            note = paste("remlax", top_rx, "| statgenGWAS", top_sg))
  top_gs <- colnames(Mt)[which.max(T_dbl)]
  add_check(mdl, "top SNP identical", "gaston", as.numeric(colnames(Mt)[which.max(rx_mlp)] == top_gs),
            1, tol = 0, metric = "abs",
            note = paste("remlax", colnames(Mt)[which.max(rx_mlp)], "| gaston", top_gs))
  set_rx <- colnames(Mt)[rx_F$mlp > thr]; set_sg <- colnames(Mt)[gw_mlp > thr]
  flip <- union(setdiff(set_rx, set_sg), setdiff(set_sg, set_rx))
  fj <- match(flip, colnames(Mt))
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
  set_chi <- colnames(Mt)[rx_mlp > thr]; set_gs <- colnames(Mt)[mlp_chi2(T_dbl) > thr]
  set_nat <- colnames(Mt)[-log10(gs$p) > thr]
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

  # ---- published association: the lead SNP near the known gene -------------
  j <- which.max(rx_F$mlp)
  near <- mapt$chr == info$chr & abs(mapt$pos - info$pos) <= WIN
  best_near <- max(rx_F$mlp[near])
  dist_lead <- if (mapt$chr[j] == info$chr) abs(mapt$pos[j] - info$pos) else Inf
  add_check(mdl, sprintf("lead SNP within 1 Mb of %s", info$gene), "Zhao 2011",
            as.numeric(dist_lead <= WIN), 1, kind = "context", metric = "abs",
            note = sprintf("lead %s chr%s:%d -log10p %.2f; best within 1 Mb of %s (chr%s ~%.2f Mb) %.2f; threshold %.2f",
                           top_rx, mapt$chr[j], as.integer(mapt$pos[j]), rx_F$mlp[j],
                           info$gene, info$chr, info$pos / 1e6, best_near, thr))

  timing[[tr]] <- data.frame(trait = tr, n = n, p = ncol(Mt), remlax_null_fit_s = t_fit,
                             remlax_scan_s = t_scan, statgenGWAS_s = t_sg,
                             gaston_aireml_s = t_ga0, gaston_score_s = t_gs)
  plots[[tr]] <- list(chr = mapt$chr, pos = mapt$pos, rx = rx_F$mlp, sg = gw_mlp, thr = thr,
                      gene = info)
}

# ---- outputs -------------------------------------------------------------------
tim <- do.call(rbind, timing)
cat("\nTimings (elapsed seconds):\n"); print(tim, row.names = FALSE)
utils::write.csv(tim, file.path(pd_root(), "results", "rice44k_zhao2011_timing.csv"),
                 row.names = FALSE)

png(file.path(pd_root(), "results", "rice44k_zhao2011_manhattan.png"),
    width = 2400, height = 600 * length(plots), res = 200)
op <- par(mfrow = c(length(plots), 1), mar = c(3.5, 4.2, 2, 1))
for (nm in names(plots)) {
  pl <- plots[[nm]]
  o <- order(as.numeric(pl$chr), pl$pos)
  x <- seq_along(o)
  col <- ifelse(as.numeric(pl$chr[o]) %% 2 == 0, "grey55", "grey20")
  plot(x, pl$sg[o], pch = 16, cex = 0.5, col = col, xaxt = "n", xlab = "chromosome",
       ylab = expression(-log[10](p)),
       main = sprintf("%s (known gene %s, chr %s)", nm, pl$gene$gene, pl$gene$chr),
       ylim = c(0, max(pl$sg, pl$rx, pl$thr) * 1.05))
  points(x, pl$rx[o], pch = 1, cex = 0.7, col = "#D55E00")
  abline(h = pl$thr, lty = 2)
  g <- which(pl$chr[o] == pl$gene$chr)
  abline(v = x[g][which.min(abs(pl$pos[o][g] - pl$gene$pos))], col = "#0072B2", lty = 3)
  mids <- tapply(x, pl$chr[o], mean)
  axis(1, at = mids, labels = names(mids), tick = FALSE)
  legend("topleft", c("statgenGWAS", "remlax (F convention)", "known gene"),
         pch = c(16, 1, NA), lty = c(NA, NA, 3), col = c("grey20", "#D55E00", "#0072B2"),
         bty = "n", cex = 0.8)
}
par(op); dev.off()

write_checks()
