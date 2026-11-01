# ==============================================================================
# run_cascade.R -- simulate one mixture, then recover its parameters step by step
# ------------------------------------------------------------------------------
#   Rscript simulation/cogwas/run_cascade.R --rep 1 --scenario S2 [--backend gpu] [--out dir]
#
# Simulation: sim_genomes + sim_design + sim_pheno_full with sim_scenario()
# (S0 additive only, S1 + major similarity and pair QTL, S2 + polygenic
# interaction background). Wheat phenotyped.
# Recovery with remlax, for two null models:
#   no_kernel  intercept; DGE and IGE intra (us 2x2, kinship of the wheat chip);
#              IGE inter (kinship of the alfalfa family doses); bin, bin:row,
#              bin:col iid; residual
#   kernels    the same plus three interaction kernels: between species
#              (Z K_w Z') * (Zb K_a Zb'), within (Z K_w Z') * (Zw K_w Zw'), and
#              similarity S S' / p with S = (Z M) * (Zw M) (Sato's K2)
# For each: 1. variances against the realised values (with standard errors);
# 2. GWAS at fixed V (dir, ind, sim|dir+ind on the wheat chip, ind on the alfalfa
# family chip); 4. co-GWAS within and between species (pair|dir+ind). Step 3,
# the exact refit of the top SNPs, runs for the no_kernel model only.
# A causal locus is found when a chip SNP within +/- win (5 cM) is below the
# Bonferroni threshold; false peaks are significant SNPs outside every causal
# window of the test, counted once per chromosome segment of width win.
# ==============================================================================
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = NULL) { i <- match(paste0("--", k), a)
  if (is.na(i) || i == length(a)) d else a[i + 1L] }
REP <- as.integer(arg("rep", "1")); SCEN <- arg("scenario", "S2"); BACKEND <- arg("backend", "cpu")
OUT <- arg("out", "cascade_out"); WIN <- as.numeric(arg("win", "5e6"))
R2 <- as.numeric(arg("r2", "0.8")); CHIP <- as.numeric(arg("chip_frac", "0.5"))
VARIANTS <- strsplit(arg("variants", "no_kernel,kernels"), ",")[[1]]
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
suppressPackageStartupMessages({ library(Matrix); library(remlax) })
source(file.path(HERE, "sim_cogwas.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
T0 <- Sys.time()
say <- function(...) cat(sprintf(...), "\n", sep = "")
tag <- sprintf("%s_rep%03d", SCEN, REP)
seed <- REP + 1000L * match(SCEN, c("S0", "S1", "S2"))

# ---------------------------------------------------------------- simulation
par <- sim_scenario(SCEN)
gen <- sim_genomes(n_causal_wheat = 1500L, n_causal_alf = 700L, chip_frac = CHIP, seed = seed)
des <- sim_design(gen, rank_intra = 2L, rank_inter = 1L, seed = seed)
ph <- sim_pheno_full(gen, des, par, seed = seed)
tr <- ph$truth; vp <- tr$var_parts
cw <- des$cells[des$wheat, ]; n <- length(ph$y)
say("%s: %d wheat plants, rho_shared %.3f, realised r_g %.3f", tag, n, tr$rho_shared, tr$realised[["r_g"]])

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
norm_k <- function(K) { K <- K / mean(diag(K)); diag(K) <- diag(K) + 1e-6; K }
Ssim <- as.matrix(Zd %*% Mw) * as.matrix(Zw %*% Mw)
KERN <- list(k_inter = rx_cogwas_kernel(Zd, Kw, Zb, Ka), k_intra = rx_cogwas_kernel(Zd, Kw, Zw, Kw),
             k_sim = norm_k(tcrossprod(Ssim) / ncol(Ssim)))
rm(Ssim)
fac <- function(x) sparse.model.matrix(~ x - 1, data.frame(x = factor(x)))
build_for <- function(variant) function(K, X_sup = NULL) {
  terms <- list(rx_term("gen", Z = cbind(Zd, Zw), K = K, t = 2L, struct = "us"),
                rx_term("ige_inter", Z = Zb, K = Ka, t = 1L),
                rx_term("bin", Z = fac(cw$bin), t = 1L),
                rx_term("bin_row", Z = fac(paste(cw$bin, cw$row)), t = 1L),
                rx_term("bin_col", Z = fac(paste(cw$bin, cw$col)), t = 1L))
  if (variant == "kernels")
    terms <- c(terms, lapply(names(KERN), function(nm) rx_term(nm, Z = Diagonal(n), K = KERN[[nm]], t = 1L)))
  rx_model(ph$y, cbind(matrix(1, n, 1), X_sup), terms = terms, residual = rx_residual())
}
# within-family part of the alfalfa effect on wheat: not in the family-level term
fam_of <- gen$alfalfa$fam_of_plant
dev_a <- tr$g$inter - tr$g$inter_fam[fam_of]
v_within <- stats::var(as.numeric(des$W_inter %*% dev_a[des$cells$plant[des$alfalfa]]))
wm <- gen$wheat$map; am <- gen$alfalfa$map; st <- tr$sets
in_win <- function(map, site_map, s) map$chr == site_map$chr[s] & abs(map$pos - site_map$pos[s]) <= WIN
find <- function(scan, col, map, site_map, sites, label) {
  if (!length(sites)) return(NULL)
  thr <- 0.05 / sum(is.finite(scan[[col]]))
  do.call(rbind, lapply(seq_along(sites), function(i) {
    w <- which(in_win(map, site_map, sites[i]) & is.finite(scan[[col]]))
    j <- if (length(w)) w[which.min(scan[[col]][w])] else NA
    data.frame(effect = label, test = col, best_p = if (is.na(j)) NA else scan[[col]][j],
               found = !is.na(j) && scan[[col]][j] < thr)
  }))
}
false_peaks <- function(scan, col, map, site_map, sites) {
  thr <- 0.05 / sum(is.finite(scan[[col]])); sig <- which(scan[[col]] < thr)
  if (!length(sig)) return(0L)
  causal <- Reduce(`|`, lapply(sites, function(s) in_win(map, site_map, s)), rep(FALSE, nrow(map)))
  out <- sig[!causal[sig]]
  if (!length(out)) 0L else length(unique(paste(map$chr[out], floor(map$pos[out] / WIN))))
}
prune <- function(M, map) { k <- rx_ld_prune(M, map, r2 = R2, window = 500L); list(M = M[, k], map = map[k, ]) }
PW <- prune(Mw, mw); PA <- prune(Ma, ma)
map_w <- data.frame(snp = colnames(Mw), chr = mw$chr, pos = mw$pos)
map_a <- data.frame(snp = colnames(Ma), chr = ma$chr, pos = ma$pos)

# ---------------------------------------------------------------- one null model
analyse <- function(variant) {
  build <- build_for(variant); mod <- build(Kw)
  t1 <- Sys.time()
  fit <- rx_fit(mod, backend = BACKEND, hessian = TRUE, blups = FALSE, verbose = FALSE)
  t_fit <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  conv <- isTRUE(as.logical(fit$conv_grad_rel)) || isTRUE(as.logical(fit$conv_decrement))
  S <- fit$sigmas$gen
  kern_on <- variant == "kernels"
  inter_parts <- c("sim", "pair_intra", "pair_inter", "poly_sim", "poly_pair_intra", "poly_pair_inter")
  in_resid <- if (kern_on) c("sim", "pair_intra", "pair_inter") else inter_parts
  truth <- c(var_DGE = tr$realised[["var_dge"]], var_IGE_intra = tr$realised[["var_ige"]],
             r_g = tr$realised[["r_g"]], var_IGE_inter = tr$realised[["var_inter_fam"]],
             var_bin = tr$realised[["var_bin"]], var_row = tr$realised[["var_row"]],
             var_col = tr$realised[["var_col"]],
             var_residual = tr$realised[["var_res"]] + sum(vp[in_resid]) + v_within)
  est <- c(S[1, 1], S[2, 2], S[1, 2] / sqrt(S[1, 1] * S[2, 2]), fit$sigmas$ige_inter[1, 1],
           fit$sigmas$bin[1, 1], fit$sigmas$bin_row[1, 1], fit$sigmas$bin_col[1, 1], as.numeric(fit$sigma_res)[1])
  if (kern_on) {
    truth <- c(truth, var_k_inter = vp[["poly_pair_inter"]], var_k_intra = vp[["poly_pair_intra"]],
               var_k_sim = vp[["poly_sim"]])
    est <- c(est, fit$sigmas$k_inter[1, 1], fit$sigmas$k_intra[1, 1], fit$sigmas$k_sim[1, 1])
  }
  vc <- data.frame(tag = tag, scenario = SCEN, rep = REP, variant = variant, parameter = names(truth),
                   simulated = unname(truth), estimated = est, converged = conv)
  vcomp <- cbind(tag = tag, scenario = SCEN, rep = REP, variant = variant, as.data.frame(fit$varcomp))
  say("\n[%s] mixed model %.0f s, converged %s", variant, t_fit, conv)
  print(transform(vc[, c("parameter", "simulated", "estimated")], simulated = round(simulated, 3),
                  estimated = round(estimated, 3)), row.names = FALSE)
  # GWAS
  t1 <- Sys.time()
  sc <- rx_scan(fit, mod, Mw, map = map_w, incidences = c(dir = "gen:1", ind = "gen:2"),
                tests = c("dir", "ind", "sim|dir+ind"), maf = NULL, backend = BACKEND, verbose = FALSE)
  sci <- rx_scan(fit, mod, list(ind = Ma), map = map_a, incidences = c(ind = "ige_inter"),
                 tests = "ind", maf = NULL, backend = BACKEND, verbose = FALSE)
  t_scan <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  gw <- rbind(find(sc, "p_dir", mw, wm, st$major_dir, "major direct"),
              find(sc, "p_ind", mw, wm, st$major_ige, "major indirect intra"),
              find(sc, "p_sim|dir+ind", mw, wm, st$sim, "similarity"),
              find(sci, "p_ind", ma, am, st$major_inter, "major indirect inter"))
  gw <- cbind(tag = tag, scenario = SCEN, rep = REP, variant = variant, gw)
  lam <- attr(sc, "meta")$lambda_gc
  gsum <- data.frame(tag = tag, scenario = SCEN, rep = REP, variant = variant,
                     test = c("dir", "ind", "sim|dir+ind", "inter ind"),
                     lambda = c(lam[["dir"]], lam[["ind"]], lam[["sim|dir+ind"]], attr(sci, "meta")$lambda_gc[["ind"]]),
                     false_peaks = c(false_peaks(sc, "p_dir", mw, wm, st$major_dir),
                                     false_peaks(sc, "p_ind", mw, wm, st$major_ige),
                                     false_peaks(sc, "p_sim|dir+ind", mw, wm, st$sim),
                                     false_peaks(sci, "p_ind", ma, am, st$major_inter)))
  say("[%s] GWAS %.0f s | lambda %s | found %d of %d", variant, t_scan,
      paste(sprintf("%.2f", gsum$lambda), collapse = " "), sum(gw$found), nrow(gw))
  # exact refit
  ex <- NULL
  if (variant == "no_kernel") {
    t1 <- Sys.time()
    ex <- do.call(rbind, lapply(c("dir", "ind"), function(tt) {
      e <- rx_scan_exact(sc, tt, build, function(W) Kw, Mw, map_w, c(dir = "gen:1", ind = "gen:2"),
                         top = 5L, loco = FALSE, batch = 0L, backend = BACKEND, verbose = FALSE)
      data.frame(tag = tag, scenario = SCEN, rep = REP, test = tt, snp = e$snp, p_fixed_V = e$p_fixed_V,
                 p_exact = e$p_exact, converged = e$converged)
    }))
    say("[%s] exact refit %.0f s", variant, as.numeric(difftime(Sys.time(), t1, units = "secs")))
  }
  # co-GWAS
  t1 <- Sys.time()
  cg <- list(inter = rx_cogwas(fit, mod, PW$M, PA$M, incidences = c(dir = "gen:1", ind = "ige_inter"),
                               tests = "pair|dir+ind", top = 5000L, backend = BACKEND, verbose = FALSE),
             intra = rx_cogwas(fit, mod, PW$M, PW$M, incidences = c(dir = "gen:1", ind = "gen:2"),
                               tests = "pair|dir+ind", top = 5000L, backend = BACKEND, verbose = FALSE))
  t_cg <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  pair_row <- function(kind, a_site, b_site, b_map, PB, label) {
    top <- cg[[kind]]$top; thr <- 0.05 / cg[[kind]]$meta$n_pairs
    hit <- if (length(a_site)) in_win(PW$map, wm, a_site)[top$focal] & in_win(PB$map, b_map, b_site)[top$neighbour]
           else rep(FALSE, nrow(top))
    data.frame(tag = tag, scenario = SCEN, rep = REP, variant = variant, effect = label,
               simulated = length(a_site) > 0, n_pairs = cg[[kind]]$meta$n_pairs,
               lambda = cg[[kind]]$meta$lambda[["pair|dir+ind"]],
               best_p = if (any(hit)) min(top$p[hit]) else NA, found = any(hit & top$p < thr),
               n_sig = sum(top$p < thr), sig_outside = sum(!hit & top$p < thr))
  }
  has_pi <- par$var_pair_intra > 0; has_pe <- par$var_pair_inter > 0
  pr <- rbind(pair_row("intra", if (has_pi) st$pair_intra[1] else integer(0), st$pair_intra[2], wm, PW, "pair intra"),
              pair_row("inter", if (has_pe) st$pair_inter[["a"]] else integer(0), st$pair_inter[["b"]], am, PA, "pair inter"))
  say("[%s] co-GWAS %.0f s | lambda intra %.3f inter %.3f | found intra %s inter %s | sig outside %d / %d",
      variant, t_cg, pr$lambda[1], pr$lambda[2], pr$found[1], pr$found[2], pr$sig_outside[1], pr$sig_outside[2])
  list(vc = vc, vcomp = vcomp, gw = gw, gsum = gsum, ex = ex, pr = pr,
       time = data.frame(tag = tag, scenario = SCEN, rep = REP, variant = variant, t_fit = t_fit,
                         t_scan = t_scan, t_cogwas = t_cg))
}

res <- lapply(VARIANTS, analyse)
bind <- function(k) do.call(rbind, lapply(res, `[[`, k))
for (k in c("vc", "vcomp", "gw", "gsum", "ex", "pr", "time")) {
  x <- bind(k)
  if (!is.null(x) && nrow(x)) utils::write.csv(x, file.path(OUT, sprintf("%s_%s.csv", tag, k)), row.names = FALSE)
}
utils::write.csv(data.frame(tag = tag, scenario = SCEN, rep = REP, t(vp), t(tr$realised),
                            rho_shared = tr$rho_shared, v_within_inter = v_within),
                 file.path(OUT, sprintf("%s_truth.csv", tag)), row.names = FALSE)
say("\n%s written to %s (%.0f s)", tag, OUT, as.numeric(difftime(Sys.time(), T0, units = "secs")))
