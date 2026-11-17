# ==============================================================================
# run_cascade.R -- simulate one mixture, then recover its parameters step by step
# ------------------------------------------------------------------------------
#   Rscript studies/cogwas/run_cascade.R --rep 1 --scenario S2 [--backend gpu] [--out dir]
#
# Simulation: sim_genomes + sim_design + sim_pheno_full with sim_scenario(),
# sum framework, order-1 neighbourhood (2 wheat and 6 alfalfa neighbours at
# the centre of a bin), parameters of Salas et al. 2026. Wheat phenotyped.
# Recovery with remlax, for two null models:
#   additive     intercept; DGE and IGE intra (us 2x2, kinship of the wheat
#                chip); IGE inter (kinship of the alfalfa family doses); bin,
#                bin:row, bin:col iid; residual
#   interaction  the same plus the DGE x IGE interaction within species,
#                (Z K_w Z') o (Zw K_w Zw'), and between species,
#                (Z K_w Z') o (Zb K_a Zb'). Both are the plant-level form of a
#                genotype-pair effect gamma(f, g) with covariance K x K.
# For each: 1. variances against the realised values, BLUP accuracies (lines,
# families, genotype pairs); 2. GWAS at fixed V (dir, ind, sim|dir+ind on the
# wheat chip, ind on the alfalfa family chip); 3. exact refit of the top SNPs
# (interaction model only); 4. co-GWAS within and between species
# (pair|dir+ind) on the LD-pruned chips.
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
VARIANTS <- strsplit(arg("variants", "additive,interaction"), ",")[[1]]
# design and effect sizes (defaults: Salas et al. 2026 values, real panel sizes)
N_LINES <- as.integer(arg("n_lines", "181")); N_FAM <- as.integer(arg("n_fam", "106"))
VAR_GXG <- as.numeric(arg("var_gxg", "0.05")); PVE_MAJOR <- as.numeric(arg("pve_major", "0.10"))
EXACT <- arg("exact", "yes") == "yes"
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
suppressPackageStartupMessages({ library(Matrix); library(remlax) })
source(file.path(HERE, "sim_cogwas.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
T0 <- Sys.time()
say <- function(...) cat(sprintf(...), "\n", sep = "")
default_cell <- N_LINES == 181L && N_FAM == 106L && VAR_GXG == 0.05 && PVE_MAJOR == 0.10
cell <- if (default_cell) "" else sprintf("L%03d_F%03d_g%.2f_q%.2f_", N_LINES, N_FAM, VAR_GXG, PVE_MAJOR)
tag <- sprintf("%s_%srep%03d", SCEN, cell, REP)
seed <- REP + 1000L * match(SCEN, c("S0", "S1", "S2")) + 7L * N_LINES
save_csv <- function(x, k) utils::write.csv(cbind(x, n_lines = N_LINES, n_fam = N_FAM, var_gxg = VAR_GXG,
                                                    pve_major = PVE_MAJOR),
                                              file.path(OUT, sprintf("%s_%s.csv", tag, k)), row.names = FALSE)

# ---------------------------------------------------------------- simulation
par <- sim_scenario(SCEN)
par$pve_major <- PVE_MAJOR
if (par$var_gxg_intra > 0) { par$var_gxg_intra <- VAR_GXG; par$var_gxg_inter <- VAR_GXG }
# about 2,120 alfalfa plants whatever the number of families
gen <- sim_genomes(n_lines = N_LINES, n_fam = N_FAM, plants_per_fam = as.integer(ceiling(2120 / N_FAM)),
                   n_causal_wheat = 1500L, n_causal_alf = 700L, chip_frac = CHIP, seed = seed)
des <- sim_design(gen, rank_intra = 1L, rank_inter = 1L, scale = "sum", seed = seed)
ph <- sim_pheno_full(gen, des, par, seed = seed)
tr <- ph$truth; vp <- tr$var_parts
cw <- des$cells[des$wheat, ]; n <- length(ph$y)
say("%s: %d wheat plants, rho_shared %.3f, realised r_g %.3f", tag, n, tr$rho_shared, tr$realised[["r_g"]])
print(round(vp / sum(vp), 3))

# ---------------------------------------------------------------- analysis data
lines <- rownames(gen$wheat$G); fams <- rownames(gen$alfalfa$G_fam)
Zd <- sparseMatrix(seq_len(n), match(cw$line, lines), x = 1, dims = c(n, length(lines)))
alf <- des$cells[des$alfalfa, ]
Zfa <- sparseMatrix(seq_len(nrow(alf)), match(alf$family, fams), x = 1, dims = c(nrow(alf), length(fams)))
Zp <- sparseMatrix(seq_len(nrow(alf)), alf$plant, x = 1, dims = c(nrow(alf), nrow(gen$alfalfa$G_plant)))
Zw <- des$W_intra %*% Zd; Zb <- des$W_inter %*% Zfa
maf_ok <- function(M) { f <- colMeans(M) / 2; pmin(f, 1 - f) >= 0.05 & apply(M, 2, stats::sd) > 0 }
Mw <- gen$wheat$G[, gen$wheat$chip]; mw <- gen$wheat$map[gen$wheat$chip, ]
k <- maf_ok(Mw); Mw <- sweep(Mw[, k], 2, colMeans(Mw[, k]), "-"); mw <- mw[k, ]
Ma <- gen$alfalfa$G_fam[, gen$alfalfa$chip]; ma <- gen$alfalfa$map[gen$alfalfa$chip, ]
k <- maf_ok(Ma); Ma <- sweep(Ma[, k], 2, colMeans(Ma[, k]), "-"); ma <- ma[k, ]
kin <- function(Mc) { K <- tcrossprod(Mc); K <- K / mean(diag(K)); 0.98 * K + 0.02 * diag(nrow(K)) }
Kw <- kin(Mw); Ka <- kin(Ma)
dimnames(Kw) <- list(lines, lines); dimnames(Ka) <- list(fams, fams)
# not normalised: the variance stays that of gamma(f, g) per neighbour
KERN <- list(k_intra = rx_cogwas_kernel(Zd, Kw, Zw, Kw, normalize = FALSE),
             k_inter = rx_cogwas_kernel(Zd, Kw, Zb, Ka, normalize = FALSE))
fac <- function(x) sparse.model.matrix(~ x - 1, data.frame(x = factor(x)))
Zsp <- list(bin = fac(cw$bin), bin_row = fac(paste(cw$bin, cw$row)), bin_col = fac(paste(cw$bin, cw$col)))
build_for <- function(variant) function(K, X_sup = NULL) {
  terms <- c(list(rx_term("gen", Z = cbind(Zd, Zw), K = K, t = 2L, struct = "us"),
                  rx_term("ige_inter", Z = Zb, K = Ka, t = 1L)),
             lapply(names(Zsp), function(nm) rx_term(nm, Z = Zsp[[nm]], t = 1L)))
  if (variant == "interaction")
    terms <- c(terms, lapply(names(KERN), function(nm) rx_term(nm, Z = Diagonal(n), K = KERN[[nm]], t = 1L)))
  rx_model(ph$y, cbind(matrix(1, n, 1), X_sup), terms = terms, residual = rx_residual())
}
# parts of the simulated signal that no term of a model can carry
fam_of <- gen$alfalfa$fam_of_plant
dev_a <- tr$g$inter - tr$g$inter_fam[fam_of]
v_within <- stats::var(as.numeric(des$W_inter %*% dev_a[alf$plant]))
u_gxg_inter <- rowSums(as.matrix(Zd %*% tr$gxg$inter) * as.matrix(des$W_inter %*% Zp))
u_gxg_inter_fam <- rowSums(as.matrix(Zd %*% tr$gxg$inter_fam) * as.matrix(Zb))
v_gxg_within <- stats::var(u_gxg_inter - u_gxg_inter_fam)
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

# BLUPs of the genotype-pair interaction effects, gamma_hat = s2 K_f (Zf' diag(alpha) Zn) K_n,
# with alpha = V^-1 (y - X beta) = e_hat / s2_e.
gxg_blups <- function(fit, mod) {
  B <- fit$blups
  fitted <- drop(mod$X %*% fit$beta) + as.numeric(Zd %*% B$gen[, 1] + Zw %*% B$gen[, 2]) +
    as.numeric(Zb %*% B$ige_inter[, 1]) +
    Reduce(`+`, lapply(names(Zsp), function(nm) as.numeric(Zsp[[nm]] %*% B[[nm]][, 1])))
  if (!is.null(B$k_intra)) fitted <- fitted + B$k_intra[, 1] + B$k_inter[, 1]
  alpha <- (ph$y - fitted) / as.numeric(fit$sigma_res)[1]
  chk <- stats::cor(B$k_intra[, 1], fit$sigmas$k_intra[1, 1] * drop(KERN$k_intra %*% alpha))
  say("  check of alpha: cor(BLUP k_intra, s2 K alpha) = %.6f", chk)
  Da <- Diagonal(x = alpha)
  list(intra = fit$sigmas$k_intra[1, 1] * Kw %*% as.matrix(t(Zd) %*% Da %*% Zw) %*% Kw,
       inter = fit$sigmas$k_inter[1, 1] * Kw %*% as.matrix(t(Zd) %*% Da %*% Zb) %*% Ka)
}

# ---------------------------------------------------------------- one null model
analyse <- function(variant) {
  build <- build_for(variant); mod <- build(Kw)
  t1 <- Sys.time()
  fit <- rx_fit(mod, backend = BACKEND, hessian = TRUE, blups = TRUE, verbose = FALSE)
  t_fit <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  conv <- isTRUE(as.logical(fit$conv_grad_rel)) || isTRUE(as.logical(fit$conv_decrement))
  S <- fit$sigmas$gen
  inter_on <- variant == "interaction"
  r_ <- tr$realised
  v_res <- r_[["var_res"]] + v_within + if (inter_on) v_gxg_within else vp[["gxg_intra"]] + vp[["gxg_inter"]]
  truth <- c(var_DGE = r_[["var_dge"]], var_IGE_intra = r_[["var_ige"]], r_g = r_[["r_g"]],
             var_IGE_inter = r_[["var_inter_fam"]], var_bin = r_[["var_bin"]], var_row = r_[["var_row"]],
             var_col = r_[["var_col"]], var_residual = v_res)
  est <- c(S[1, 1], S[2, 2], S[1, 2] / sqrt(S[1, 1] * S[2, 2]), fit$sigmas$ige_inter[1, 1],
           fit$sigmas$bin[1, 1], fit$sigmas$bin_row[1, 1], fit$sigmas$bin_col[1, 1], as.numeric(fit$sigma_res)[1])
  if (inter_on) {
    truth <- c(truth, var_GxG_intra = r_[["var_gxg_intra"]], var_GxG_inter = r_[["var_gxg_inter_fam"]])
    est <- c(est, fit$sigmas$k_intra[1, 1], fit$sigmas$k_inter[1, 1])
  }
  vc <- data.frame(tag = tag, scenario = SCEN, rep = REP, variant = variant, parameter = names(truth),
                   simulated = unname(truth), estimated = est, converged = conv)
  vcomp <- cbind(tag = tag, scenario = SCEN, rep = REP, variant = variant, as.data.frame(fit$varcomp))
  say("\n[%s] mixed model %.0f s, converged %s", variant, t_fit, conv)
  print(transform(vc[, c("parameter", "simulated", "estimated")], simulated = round(simulated, 4),
                  estimated = round(estimated, 4)), row.names = FALSE)
  # BLUPs against the simulated values
  B <- fit$blups
  bl <- data.frame(tag = tag, variant = variant, line = lines, dge = tr$g$dge, ige = tr$g$ige,
                   dge_hat = B$gen[, 1], ige_hat = B$gen[, 2])
  bf <- data.frame(tag = tag, variant = variant, family = fams, inter = as.numeric(tr$g$inter_fam[fams]),
                   inter_hat = B$ige_inter[, 1])
  acc <- data.frame(effect = c("DGE", "IGE intra", "IGE inter"),
                    accuracy = c(stats::cor(bl$dge, bl$dge_hat), stats::cor(bl$ige, bl$ige_hat),
                                 stats::cor(bf$inter, bf$inter_hat)))
  gx <- NULL
  if (inter_on && par$var_gxg_intra > 0) {
    gh <- gxg_blups(fit, mod)
    gx <- rbind(data.frame(tag = tag, effect = "GxG intra", truth = as.vector(tr$gxg$intra), blup = as.vector(gh$intra)),
                data.frame(tag = tag, effect = "GxG inter", truth = as.vector(tr$gxg$inter_fam), blup = as.vector(gh$inter)))
    acc <- rbind(acc, data.frame(effect = c("GxG intra", "GxG inter"),
                                 accuracy = c(stats::cor(as.vector(tr$gxg$intra), as.vector(gh$intra)),
                                              stats::cor(as.vector(tr$gxg$inter_fam), as.vector(gh$inter)))))
  }
  acc <- cbind(tag = tag, scenario = SCEN, rep = REP, variant = variant, acc)
  say("[%s] BLUP accuracy: %s", variant, paste(sprintf("%s %.2f", acc$effect, acc$accuracy), collapse = ", "))
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
  scan_out <- rbind(data.frame(variant = variant, species = "wheat", map_w, p_dir = sc$p_dir, p_ind = sc$p_ind,
                               p_sim = sc[["p_sim|dir+ind"]]),
                    data.frame(variant = variant, species = "alfalfa", map_a, p_dir = NA, p_ind = sci$p_ind, p_sim = NA))
  # exact refit
  ex <- NULL
  if (inter_on && EXACT) {
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
  # causal pairs of each co-GWAS: similarity loci are the diagonal pairs (s, s)
  causal_pairs <- list(
    intra = c(if (par$var_gxg_intra > 0) list(list(a = st$pair_intra[1], b = st$pair_intra[2], label = "pair intra")),
              lapply(st$sim, function(s) list(a = s, b = s, label = "similarity (diagonal pair)"))),
    inter = if (par$var_gxg_inter > 0) list(list(a = st$pair_inter[["a"]], b = st$pair_inter[["b"]], label = "pair inter")))
  nb_map <- list(intra = list(map = PW$map, site = wm), inter = list(map = PA$map, site = am))
  pr <- NULL; cg_out <- NULL
  for (kind in c("intra", "inter")) {
    top <- cg[[kind]]$top; thr <- 0.05 / cg[[kind]]$meta$n_pairs
    hits <- lapply(causal_pairs[[kind]], function(cp)
      in_win(PW$map, wm, cp$a)[top$focal] & in_win(nb_map[[kind]]$map, nb_map[[kind]]$site, cp$b)[top$neighbour])
    any_hit <- Reduce(`|`, hits, rep(FALSE, nrow(top)))
    lam_k <- cg[[kind]]$meta$lambda[["pair|dir+ind"]]
    rows <- lapply(seq_along(hits), function(i) data.frame(
      effect = causal_pairs[[kind]][[i]]$label, best_p = if (any(hits[[i]])) min(top$p[hits[[i]]]) else NA,
      found = any(hits[[i]] & top$p < thr)))
    rows <- c(rows, list(data.frame(effect = paste("outside", kind), best_p = NA, found = NA)))
    pr <- rbind(pr, cbind(tag = tag, scenario = SCEN, rep = REP, variant = variant, kind = kind,
                          do.call(rbind, rows), n_pairs = cg[[kind]]$meta$n_pairs, lambda = lam_k,
                          n_sig = sum(top$p < thr), sig_outside = sum(!any_hit & top$p < thr)))
    cg_out <- rbind(cg_out, data.frame(variant = variant, kind = kind,
                                       focal_chr = PW$map$chr[top$focal], focal_pos = PW$map$pos[top$focal],
                                       neigh_chr = nb_map[[kind]]$map$chr[top$neighbour],
                                       neigh_pos = nb_map[[kind]]$map$pos[top$neighbour],
                                       p = top$p, causal = any_hit, threshold = thr))
  }
  say("[%s] co-GWAS %.0f s | lambda intra %.3f inter %.3f | %s", variant, t_cg,
      pr$lambda[pr$kind == "intra"][1], pr$lambda[pr$kind == "inter"][1],
      paste(sprintf("%s %s", pr$effect, pr$found), collapse = ", "))
  list(vc = vc, vcomp = vcomp, acc = acc, blup_lines = bl, blup_fam = bf, gxg = gx, gw = gw, gsum = gsum,
       scan = scan_out, ex = ex, pr = pr, cg = cg_out,
       time = data.frame(tag = tag, scenario = SCEN, rep = REP, variant = variant, t_fit = t_fit,
                         t_scan = t_scan, t_cogwas = t_cg))
}

res <- lapply(VARIANTS, analyse)
for (k in names(res[[1]])) {
  x <- do.call(rbind, lapply(res, `[[`, k))
  if (!is.null(x) && nrow(x)) save_csv(x, k)
}
save_csv(data.frame(tag = tag, scenario = SCEN, rep = REP, t(vp), t(tr$realised),
                    rho_shared = tr$rho_shared, v_within_inter = v_within, v_gxg_within = v_gxg_within), "truth")
# simulated causal loci, for the figures
loci <- function(s, map, label) if (length(s)) data.frame(effect = label, chr = map$chr[s], pos = map$pos[s]) else NULL
save_csv(rbind(loci(st$major_dir, wm, "major direct"), loci(st$major_ige, wm, "major indirect intra"),
               loci(st$sim, wm, "similarity"), loci(st$major_inter, am, "major indirect inter"),
               loci(st$pair_intra, wm, c("pair intra focal", "pair intra neighbour")),
               loci(st$pair_inter[["a"]], wm, "pair inter focal"), loci(st$pair_inter[["b"]], am, "pair inter neighbour")),
         "loci")
say("\n%s written to %s (%.0f s)", tag, OUT, as.numeric(difftime(Sys.time(), T0, units = "secs")))
