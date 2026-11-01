# ==============================================================================
# run_cascade.R -- simulate one mixture, then recover its parameters step by step
# ------------------------------------------------------------------------------
#   Rscript simulation/cogwas/run_cascade.R --rep 1 [--backend gpu] [--out dir]
#
# Simulation: sim_genomes + sim_design + sim_pheno_full (sim_cogwas.R), default
# parameters of sim_par_default(), wheat phenotyped.
# Recovery, all with remlax:
#   1. mixed model: DGE and IGE intra (us 2x2, kinship of the wheat chip), IGE
#      inter (kinship of the alfalfa family doses), bin, bin:row, bin:col, residual
#      -> variances and r_g against the realised values;
#   2. GWAS at fixed V (rx_scan): dir, ind, sim, sim|dir+ind on the wheat chip;
#      ind on the alfalfa family chip -> each simulated QTL found or not;
#   3. exact refit of the top SNPs (rx_scan_exact, batched on the device);
#   4. co-GWAS (rx_cogwas) within and between species, test pair|dir+ind.
# A causal locus is "found" when a chip SNP within +/- win (5 cM) of it is below
# the Bonferroni threshold of the test; significant SNPs outside every causal
# window of that test are counted as false positives (by peak: one per
# chromosome segment of width win).
# ==============================================================================
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = NULL) { i <- match(paste0("--", k), a)
  if (is.na(i) || i == length(a)) d else a[i + 1L] }
REP <- as.integer(arg("rep", "1")); BACKEND <- arg("backend", "cpu")
OUT <- arg("out", "cascade_out"); WIN <- as.numeric(arg("win", "5e6"))
R2 <- as.numeric(arg("r2", "0.8")); CHIP <- as.numeric(arg("chip_frac", "0.5"))
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
suppressPackageStartupMessages({ library(Matrix); library(remlax) })
source(file.path(HERE, "sim_cogwas.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
T0 <- Sys.time()
say <- function(...) cat(sprintf(...), "\n", sep = "")

# ---------------------------------------------------------------- simulation
par <- sim_par_default()
gen <- sim_genomes(n_causal_wheat = 1500L, n_causal_alf = 700L, chip_frac = CHIP, seed = REP)
des <- sim_design(gen, rank_intra = 2L, rank_inter = 1L, seed = REP)
ph <- sim_pheno_full(gen, des, par, seed = REP)
tr <- ph$truth
cw <- des$cells[des$wheat, ]; n <- length(ph$y)
say("simulation: %d wheat plants, rho_shared %.3f, realised r_g %.3f", n, tr$rho_shared, tr$realised[["r_g"]])

# ---------------------------------------------------------------- analysis data
lines <- rownames(gen$wheat$G); fams <- rownames(gen$alfalfa$G_fam)
Zd <- sparseMatrix(seq_len(n), match(cw$line, lines), x = 1, dims = c(n, length(lines)))
alf <- des$cells[des$alfalfa, ]
Zfa <- sparseMatrix(seq_len(nrow(alf)), match(alf$family, fams), x = 1, dims = c(nrow(alf), length(fams)))
Zw <- des$W_intra %*% Zd; Zb <- des$W_inter %*% Zfa
maf_ok <- function(M) { f <- colMeans(M) / 2; pmin(f, 1 - f) >= 0.05 & apply(M, 2, stats::sd) > 0 }
Mw <- gen$wheat$G[, gen$wheat$chip]; mw <- gen$wheat$map[gen$wheat$chip, ]
k <- maf_ok(Mw); Mw <- sweep(Mw[, k], 2, colMeans(Mw[, k]), "-"); mw <- mw[k, ]
Ma <- gen$alfalfa$G_fam[, gen$alfalfa$chip]; ma <- gen$alfalfa$map[gen$alfalfa$chip, ]
k <- maf_ok(Ma); Ma <- sweep(Ma[, k], 2, colMeans(Ma[, k]), "-"); ma <- ma[k, ]
kin <- function(Mc) { K <- tcrossprod(Mc); K <- K / mean(diag(K)); 0.98 * K + 0.02 * diag(nrow(K)) }
Kw <- kin(Mw); Ka <- kin(Ma)
dimnames(Kw) <- list(lines, lines); dimnames(Ka) <- list(fams, fams)
fac <- function(x) sparse.model.matrix(~ x - 1, data.frame(x = factor(x)))
build <- function(K, X_sup = NULL) {
  rx_model(ph$y, cbind(matrix(1, n, 1), X_sup),
           terms = list(rx_term("gen", Z = cbind(Zd, Zw), K = K, t = 2L, struct = "us"),
                        rx_term("ige_inter", Z = Zb, K = Ka, t = 1L),
                        rx_term("bin", Z = fac(cw$bin), t = 1L),
                        rx_term("bin_row", Z = fac(paste(cw$bin, cw$row)), t = 1L),
                        rx_term("bin_col", Z = fac(paste(cw$bin, cw$col)), t = 1L)),
           residual = rx_residual())
}

# ---------------------------------------------------------------- 1. mixed model
mod <- build(Kw)
t1 <- Sys.time()
fit <- rx_fit(mod, backend = BACKEND, hessian = TRUE, blups = FALSE, verbose = FALSE)
t_fit <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
S <- fit$sigmas$gen; se <- fit$se_theta
v_inter_true <- tr$realised[["var_inter_fam"]]
nonadd <- tr$share_y[c("sim", "pair_intra", "pair_inter")] * tr$realised[["var_y"]]
vc <- data.frame(
  parameter = c("var_DGE", "var_IGE_intra", "r_g", "var_IGE_inter (family)", "var_bin", "var_row",
                "var_col", "var_residual (+ sim, pairs, within-family IGE inter)"),
  simulated = c(tr$realised[["var_dge"]], tr$realised[["var_ige"]], tr$realised[["r_g"]], v_inter_true,
                tr$realised[["var_bin"]], tr$realised[["var_row"]], tr$realised[["var_col"]],
                tr$realised[["var_res"]] + sum(nonadd)),
  estimated = c(S[1, 1], S[2, 2], S[1, 2] / sqrt(S[1, 1] * S[2, 2]), fit$sigmas$ige_inter[1, 1],
                fit$sigmas$bin[1, 1], fit$sigmas$bin_row[1, 1], fit$sigmas$bin_col[1, 1],
                as.numeric(fit$sigma_res)[1]))
say("\n1. mixed model (%.0f s, converged %s)", t_fit,
    isTRUE(as.logical(fit$conv_grad_rel)) || isTRUE(as.logical(fit$conv_decrement)))
print(transform(vc, simulated = round(simulated, 3), estimated = round(estimated, 3)), row.names = FALSE)

# ---------------------------------------------------------------- 2. GWAS
t1 <- Sys.time()
TESTS <- c("dir", "ind", "sim", "sim|dir+ind")
sc <- rx_scan(fit, mod, Mw, map = data.frame(snp = colnames(Mw), chr = mw$chr, pos = mw$pos),
              incidences = c(dir = "gen:1", ind = "gen:2"), tests = TESTS, maf = NULL,
              backend = BACKEND, verbose = FALSE)
sci <- rx_scan(fit, mod, list(ind = Ma), map = data.frame(snp = colnames(Ma), chr = ma$chr, pos = ma$pos),
               incidences = c(ind = "ige_inter"), tests = "ind", maf = NULL, backend = BACKEND,
               verbose = FALSE)
t_scan <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
in_win <- function(map, site_map, s) map$chr == site_map$chr[s] & abs(map$pos - site_map$pos[s]) <= WIN
find <- function(scan, col, map, site_map, sites, label, beta) {
  thr <- 0.05 / sum(is.finite(scan[[col]]))
  do.call(rbind, lapply(seq_along(sites), function(i) {
    w <- which(in_win(map, site_map, sites[i]) & is.finite(scan[[col]]))
    j <- if (length(w)) w[which.min(scan[[col]][w])] else NA
    data.frame(effect = label, locus = sprintf("%d:%.0f", site_map$chr[sites[i]], site_map$pos[sites[i]]),
               simulated_beta = beta[i], test = col, best_p = if (is.na(j)) NA else scan[[col]][j],
               found = !is.na(j) && scan[[col]][j] < thr, threshold = thr)
  }))
}
false_peaks <- function(scan, col, map, site_map, sites) {
  thr <- 0.05 / sum(is.finite(scan[[col]]))
  sig <- which(scan[[col]] < thr)
  if (!length(sig)) return(0L)
  causal <- Reduce(`|`, lapply(sites, function(s) in_win(map, site_map, s)), rep(FALSE, nrow(map)))
  out <- sig[!causal[sig]]
  if (!length(out)) return(0L)
  length(unique(paste(map$chr[out], floor(map$pos[out] / WIN))))
}
wm <- gen$wheat$map; am <- gen$alfalfa$map
st <- tr$sets
gw <- rbind(find(sc, "p_dir", mw, wm, st$major_dir, "major direct", tr$beta$major_dir),
            find(sc, "p_ind", mw, wm, st$major_ige, "major indirect intra", tr$beta$major_ige),
            find(sc, "p_sim", mw, wm, st$sim, "similarity", tr$beta$sim),
            find(sc, "p_sim|dir+ind", mw, wm, st$sim, "similarity | dir+ind", tr$beta$sim),
            find(sci, "p_ind", ma, am, st$major_inter, "major indirect inter", tr$beta$major_inter))
fp <- c(dir = false_peaks(sc, "p_dir", mw, wm, st$major_dir),
        ind = false_peaks(sc, "p_ind", mw, wm, st$major_ige),
        sim = false_peaks(sc, "p_sim|dir+ind", mw, wm, st$sim),
        inter = false_peaks(sci, "p_ind", ma, am, st$major_inter))
lam <- attr(sc, "meta")$lambda_gc
say("\n2. GWAS (%.0f s) | lambda dir %.2f ind %.2f sim %.2f sim|dir+ind %.2f inter %.2f",
    t_scan, lam[["dir"]], lam[["ind"]], lam[["sim"]], lam[["sim|dir+ind"]], attr(sci, "meta")$lambda_gc[["ind"]])
print(transform(gw, simulated_beta = signif(simulated_beta, 3), best_p = signif(best_p, 3),
                threshold = signif(threshold, 3)), row.names = FALSE)
say("   false peaks (significant outside causal windows): %s",
    paste(sprintf("%s %d", names(fp), fp), collapse = ", "))

# ---------------------------------------------------------------- 3. exact refit
t1 <- Sys.time()
map_w <- data.frame(snp = colnames(Mw), chr = mw$chr, pos = mw$pos)
ex <- lapply(c("dir", "ind"), function(tt)
  rx_scan_exact(sc, tt, build, function(W) Kw, Mw, map_w, c(dir = "gen:1", ind = "gen:2"),
                top = 5L, loco = FALSE, batch = 0L, backend = BACKEND, verbose = FALSE))
t_ex <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
say("\n3. exact refit of the 5 top SNPs per test (%.0f s)", t_ex)
for (i in 1:2) print(ex[[i]][, c("snp", "chr", "pos", "p_fixed_V", "p_exact", "converged")], row.names = FALSE)

# ---------------------------------------------------------------- 4. co-GWAS
prune <- function(M, map) { k <- rx_ld_prune(M, map, r2 = R2, window = 500L); list(M = M[, k], map = map[k, ]) }
PW <- prune(Mw, mw); PA <- prune(Ma, ma)
t1 <- Sys.time()
cg <- list(
  inter = rx_cogwas(fit, mod, PW$M, PA$M, incidences = c(dir = "gen:1", ind = "ige_inter"),
                    tests = "pair|dir+ind", top = 5000L, backend = BACKEND, verbose = FALSE),
  intra = rx_cogwas(fit, mod, PW$M, PW$M, incidences = c(dir = "gen:1", ind = "gen:2"),
                    tests = "pair|dir+ind", top = 5000L, backend = BACKEND, verbose = FALSE))
t_cg <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
pair_row <- function(kind, a_site, b_site, b_map, PB, beta) {
  top <- cg[[kind]]$top; thr <- 0.05 / cg[[kind]]$meta$n_pairs
  hit <- in_win(PW$map, wm, a_site)[top$focal] & in_win(PB$map, b_map, b_site)[top$neighbour]
  data.frame(effect = paste("pair", kind), focal = sprintf("%d:%.0f", wm$chr[a_site], wm$pos[a_site]),
             neighbour = sprintf("%d:%.0f", b_map$chr[b_site], b_map$pos[b_site]), simulated_beta = beta,
             n_pairs = cg[[kind]]$meta$n_pairs, lambda = cg[[kind]]$meta$lambda[["pair|dir+ind"]],
             best_p = if (any(hit)) min(top$p[hit]) else NA, found = any(hit & top$p < thr),
             sig_outside = sum(!hit & top$p < thr))
}
pr <- rbind(pair_row("intra", st$pair_intra[1], st$pair_intra[2], wm, PW, tr$beta$pair_intra),
            pair_row("inter", st$pair_inter[["a"]], st$pair_inter[["b"]], am, PA, tr$beta$pair_inter))
sim_diag <- do.call(rbind, lapply(seq_along(st$sim), function(i) {
  pr_ <- pair_row("intra", st$sim[i], st$sim[i], wm, PW, tr$beta$sim[i]); pr_$effect <- "similarity (diagonal pair)"; pr_ }))
say("\n4. co-GWAS (%.0f s), %d x %d wheat SNP and %d alfalfa SNP after pruning", t_cg,
    ncol(PW$M), ncol(PW$M), ncol(PA$M))
print(transform(rbind(pr, sim_diag), simulated_beta = signif(simulated_beta, 3), best_p = signif(best_p, 3),
                lambda = round(lambda, 3)), row.names = FALSE)

# ---------------------------------------------------------------- outputs
utils::write.csv(vc, file.path(OUT, sprintf("rep%03d_variances.csv", REP)), row.names = FALSE)
utils::write.csv(gw, file.path(OUT, sprintf("rep%03d_gwas.csv", REP)), row.names = FALSE)
utils::write.csv(rbind(pr, sim_diag), file.path(OUT, sprintf("rep%03d_cogwas.csv", REP)), row.names = FALSE)
utils::write.csv(data.frame(rep = REP, t(fp), lambda_dir = lam[["dir"]], lambda_ind = lam[["ind"]],
                            lambda_sim = lam[["sim|dir+ind"]], t_fit = t_fit, t_scan = t_scan, t_exact = t_ex,
                            t_cogwas = t_cg, t_total = as.numeric(difftime(Sys.time(), T0, units = "secs"))),
                 file.path(OUT, sprintf("rep%03d_summary.csv", REP)), row.names = FALSE)
saveRDS(list(truth = tr, par = par), file.path(OUT, sprintf("rep%03d_truth.rds", REP)))
say("\nwritten to %s (%.0f s in total)", OUT, as.numeric(difftime(Sys.time(), T0, units = "secs")))
