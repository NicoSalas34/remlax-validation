# ==============================================================================
# run_bench.R -- one replicate of the co-GWAS benchmark on simulated mixtures
# ------------------------------------------------------------------------------
#   Rscript simulation/cogwas/run_bench.R --rep 1 --pair_inter 0.02 --pair_intra 0.02 \
#           [--kernel yes] [--backend gpu] [--chip_frac 0.5] [--r2 0.8] [--out results.csv]
#
# Analysis model (what a user would fit): intercept; DGE and IGE intra (us 2x2,
# VanRaden kinship of the wheat chip), IGE inter (kinship of the alfalfa family
# doses); bin, bin:row and bin:column iid; optionally the interaction kernels
# (rx_cogwas_kernel) inter and intra; iid residual. Then rx_cogwas between
# species (wheat chip x alfalfa family doses) and within (wheat x wheat), test
# pair|dir+ind, after MAF >= 0.05 and LD pruning.
#
# Detection of a simulated pair: a "window pair" is a pair of chip SNPs within
# `win` (pseudo-bp, 5e6 = 5 cM) of the causal focal and neighbour sites. Reported:
# the smallest p among window pairs in the top table, whether it is below the
# Bonferroni threshold, and the p of the best-tag pair (chip SNPs of highest r^2
# with each causal site), computed by a direct call on those two SNPs.
# ==============================================================================
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = NULL) { i <- match(paste0("--", k), a)
  if (is.na(i) || i == length(a)) d else a[i + 1L] }
REP <- as.integer(arg("rep", "1"))
P_INTER <- as.numeric(arg("pair_inter", "0")); P_INTRA <- as.numeric(arg("pair_intra", "0"))
KERNEL <- tolower(arg("kernel", "yes")) %in% c("yes", "oui", "true")
BACKEND <- arg("backend", "cpu"); CHIP <- as.numeric(arg("chip_frac", "0.5"))
R2 <- as.numeric(arg("r2", "0.8")); TOP <- as.integer(arg("top", "5000"))
WIN <- as.numeric(arg("win", "5e6")); OUT <- arg("out", "bench_cogwas.csv")
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))

suppressPackageStartupMessages({ library(Matrix); library(remlax) })
source(file.path(HERE, "sim_cogwas.R"))
t0 <- Sys.time()

gen <- sim_genomes(chip_frac = CHIP, seed = REP)
des <- sim_design(gen, rank_intra = 2L, rank_inter = 1L, seed = REP)
ph  <- sim_pheno(gen, des, pair_intra = P_INTRA, pair_inter = P_INTER, seed = REP)
cw <- des$cells[des$wheat, ]
n <- length(ph$y)

# ---------------------------------------------------------------- analysis data
lines <- rownames(gen$wheat$G); fams <- rownames(gen$alfalfa$G_fam)
Zl <- sparseMatrix(seq_len(n), match(cw$line, lines), x = 1, dims = c(n, length(lines)))
alf <- des$cells[des$alfalfa, ]
Zf_alf <- sparseMatrix(seq_len(nrow(alf)), match(alf$family, fams), x = 1,
                       dims = c(nrow(alf), length(fams)))
Zd <- Zl; Zw <- des$W_intra %*% Zl; Zb <- des$W_inter %*% Zf_alf
prep <- function(M, map) {
  f <- colMeans(M) / 2; g <- pmin(f, 1 - f) >= 0.05 & apply(M, 2, stats::sd) > 0
  M <- M[, g, drop = FALSE]; map <- map[g, , drop = FALSE]
  k <- rx_ld_prune(M, map, r2 = R2, window = 500L)
  list(M = sweep(M[, k, drop = FALSE], 2, colMeans(M[, k, drop = FALSE]), "-"), map = map[k, ],
       idx = which(g)[k])
}
Mw <- gen$wheat$G[, gen$wheat$chip]; Ma <- gen$alfalfa$G_fam[, gen$alfalfa$chip]
DW <- prep(Mw, gen$wheat$map[gen$wheat$chip, ]); DA <- prep(Ma, gen$alfalfa$map[gen$alfalfa$chip, ])
kin <- function(Mc) { K <- tcrossprod(Mc); K <- K / mean(diag(K)); 0.98 * K + 0.02 * diag(nrow(K)) }
Kw <- kin(DW$M); Ka <- kin(DA$M)

# ---------------------------------------------------------------- null model
fac <- function(x) sparse.model.matrix(~ x - 1, data.frame(x = factor(x)))
terms <- list(rx_term("gen", Z = cbind(Zd, Zw), K = Kw, t = 2L, struct = "us"),
              rx_term("ige_inter", Z = Zb, K = Ka, t = 1L),
              rx_term("bin", Z = fac(cw$bin), t = 1L),
              rx_term("bin_row", Z = fac(paste(cw$bin, cw$row)), t = 1L),
              rx_term("bin_col", Z = fac(paste(cw$bin, cw$col)), t = 1L))
if (KERNEL) terms <- c(terms, list(
  rx_term("pair_inter", Z = Diagonal(n), K = rx_cogwas_kernel(Zd, Kw, Zb, Ka), t = 1L),
  rx_term("pair_intra", Z = Diagonal(n), K = rx_cogwas_kernel(Zd, Kw, Zw, Kw), t = 1L)))
mod <- rx_model(ph$y, matrix(1, n, 1), terms = terms, residual = rx_residual())
t1 <- Sys.time()
fit <- rx_fit(mod, backend = BACKEND, hessian = FALSE, blups = FALSE, verbose = FALSE)
t_fit <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
conv <- isTRUE(as.logical(fit$conv_grad_rel)) || isTRUE(as.logical(fit$conv_decrement))

# ---------------------------------------------------------------- co-GWAS
causal_tags <- function(G, causal_site, D) {
  # chip SNP (after pruning) of highest r^2 with the causal site, and the window
  x <- G[, causal_site]; r2 <- suppressWarnings(stats::cor(D$M, x)^2); r2[!is.finite(r2)] <- 0
  list(best = which.max(r2), r2 = max(r2))
}
win_of <- function(map_all, D, site) {
  ch <- map_all$chr[site]; ps <- map_all$pos[site]
  which(D$map$chr == ch & abs(D$map$pos - ps) <= WIN)
}
one <- function(kind) {
  if (kind == "inter") {
    Mb <- DA; inc <- c(dir = "gen:1", ind = "ige_inter")
    Gb <- gen$alfalfa$G_fam; map_b <- gen$alfalfa$map
  } else {
    Mb <- DW; inc <- c(dir = "gen:1", ind = "gen:2")
    Gb <- gen$wheat$G; map_b <- gen$wheat$map
  }
  t1 <- Sys.time()
  cg <- rx_cogwas(fit, mod, DW$M, Mb$M, incidences = inc, tests = "pair|dir+ind", top = TOP,
                  backend = BACKEND, verbose = FALSE)
  tt <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  top <- cg$top; thr <- 0.05 / cg$meta$n_pairs
  out <- data.frame(kind = kind, n_pairs = cg$meta$n_pairs, lambda = cg$meta$lambda[["pair|dir+ind"]],
                    p_min = min(top$p), n_sig = sum(top$p < thr), seconds = tt,
                    simulated = !is.null(ph$truth[[kind]]))
  tr <- ph$truth[[kind]]
  if (!is.null(tr)) {
    wa <- win_of(gen$wheat$map, DW, tr[["a"]]); wb <- win_of(map_b, Mb, tr[["b"]])
    hit <- top$focal %in% wa & top$neighbour %in% wb
    out$p_window <- if (any(hit)) min(top$p[hit]) else NA
    out$sig_window <- any(hit & top$p < thr)
    out$false_sig <- sum(!hit & top$p < thr)
    ta <- causal_tags(gen$wheat$G, tr[["a"]], DW); tb <- causal_tags(Gb, tr[["b"]], Mb)
    cb <- rx_cogwas(fit, mod, DW$M[, ta$best, drop = FALSE], Mb$M[, tb$best, drop = FALSE],
                    incidences = inc, tests = "pair|dir+ind", top = 1, backend = BACKEND,
                    verbose = FALSE)
    out$r2_tag_focal <- ta$r2; out$r2_tag_neighbour <- tb$r2; out$p_best_tag <- cb$top$p[1]
  } else {
    out$p_window <- NA; out$sig_window <- NA; out$false_sig <- out$n_sig
    out$r2_tag_focal <- NA; out$r2_tag_neighbour <- NA; out$p_best_tag <- NA
  }
  out
}
res <- rbind(one("inter"), one("intra"))
res <- cbind(rep = REP, pair_inter = P_INTER, pair_intra = P_INTRA, kernel = KERNEL, chip_frac = CHIP,
             r2 = R2, snp_wheat = ncol(DW$M), snp_alf = ncol(DA$M), converged = conv,
             seconds_fit = t_fit, res,
             seconds_total = as.numeric(difftime(Sys.time(), t0, units = "secs")))
utils::write.table(res, OUT, sep = ",", row.names = FALSE, append = file.exists(OUT),
                   col.names = !file.exists(OUT))
print(res[, c("rep", "kind", "n_pairs", "lambda", "p_min", "n_sig", "p_window", "sig_window",
              "false_sig", "p_best_tag", "r2_tag_focal", "r2_tag_neighbour", "seconds")], row.names = FALSE)
