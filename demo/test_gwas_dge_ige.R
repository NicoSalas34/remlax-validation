# ==============================================================================
# One-step GWAS with remlax on a simple DGE + IGE data set: rx_scan, LOCO
# (rx_scan_loco) and exact refit of the top hits (rx_scan_exact)
#
# 300 genotypes, 4 replicates, field of 80 rows x 15 columns, blocks of 5.
# Markers: 5 chromosomes x 500 SNPs, linkage disequilibrium along the
# chromosome (latent AR(1) haplotypes) and two subpopulations. Each plot
# receives the IGE of its 4 row and column neighbours (L2-normalised incidence).
# Truth: correlated polygenic DGE and IGE (r = -0.3), plus two causal SNPs, one
# with a DIRECT effect on chromosome 2 (4 % of the phenotypic variance) and one
# with an INDIRECT effect (dose of the neighbours) on chromosome 4 (3 %).
#
#   Rscript demo/test_gwas_dge_ige.R [output_folder]
#
# remlax is the installed package, loaded by R/charger_remlax.R (path relative
# to this script), which also exposes the internal functions.
# ==============================================================================

here_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(here_dir, "..", "R", "charger_remlax.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
out_dir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "output_gwas_dge_ige"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(11)

# ------------------------------------------------------------------------------
# Markers
# ------------------------------------------------------------------------------
q <- 300L; n_chr <- 5L; m_chr <- 500L
gen_names <- sprintf("G%03d", seq_len(q))
pop <- rep(1:2, length.out = q)
haplotypes <- function(m) {                   # 2q latent AR(1) haplotypes, rho = 0.9
  z <- matrix(0, 2 * q, m); z[, 1] <- rnorm(2 * q)
  for (j in 2:m) z[, j] <- 0.9 * z[, j - 1] + sqrt(1 - 0.81) * rnorm(2 * q)
  z
}
M <- do.call(cbind, lapply(seq_len(n_chr), function(cc) {
  p0 <- runif(m_chr, 0.1, 0.5)
  pp <- sapply(p0, function(p) rbeta(2, p * 19, (1 - p) * 19))       # Fst ~ 0.05
  z <- haplotypes(m_chr)
  a <- pnorm(z) < pp[rep(pop, 2), ]                                  # allele per haplotype
  a[seq_len(q), ] + a[q + seq_len(q), ]
}))
dimnames(M) <- list(gen_names, sprintf("snp%d_%03d", rep(seq_len(n_chr), each = m_chr), rep(seq_len(m_chr), n_chr)))
snp_map <- data.frame(snp = colnames(M), chr = rep(seq_len(n_chr), each = m_chr), pos = rep(seq_len(m_chr), n_chr))
f <- colMeans(M) / 2
keep <- pmin(f, 1 - f) >= 0.05
M <- M[, keep]; snp_map <- snp_map[keep, ]
Mc <- sweep(M, 2, colMeans(M))                                       # centred doses
cat(sprintf("%d SNPs after MAF >= 0.05\n", ncol(M)))

# Kinship: cross-product of the centred doses, mean diagonal scaled to 1,
# blended with 2 % identity. The same function serves the full null and the LOCO nulls.
K_of <- function(W) {
  G <- tcrossprod(W); G <- G / mean(diag(G))
  G <- 0.98 * G + 0.02 * diag(nrow(G)); dimnames(G) <- list(gen_names, gen_names); G
}
K <- K_of(Mc)

# ------------------------------------------------------------------------------
# Design, neighbourhood and phenotypes
# ------------------------------------------------------------------------------
d <- expand.grid(col = 1:15, row = 1:80)
d$rep <- factor(sprintf("R%d", (d$row - 1L) %/% 20L + 1L))
d$block <- factor(sprintf("%s_r%02d_b%d", d$rep, d$row, (d$col - 1L) %/% 5L + 1L))
d$plot <- sprintf("P%04d", seq_len(nrow(d)))
d$gid <- factor(NA, levels = gen_names)
for (r in levels(d$rep)) d$gid[d$rep == r] <- sample(rep(gen_names, length.out = sum(d$rep == r)))
n <- nrow(d)
nb <- rx_neighbourhood(cbind(d$row, d$col), group = rep("A", n), level = d$gid, id = d$plot,
                       rank = 1, window = "euclidean", kernel = "none", normalise = TRUE)
Zn <- nb$level[["A<-A"]][, gen_names]; dimnames(Zn) <- list(NULL, gen_names)
Zd <- Matrix::sparseMatrix(i = seq_len(n), j = as.integer(d$gid), x = 1, dims = c(n, q),
                           dimnames = list(NULL, gen_names))
k_ige <- rx_exposure(Zd, Zn, K = K)$k

# Polygenic: DGE 0.5; IGE 0.15 of phenotypic variance; r = -0.3
s2i <- 0.15 / k_ige
Sig <- matrix(c(0.5, rep(-0.3 * sqrt(0.5 * s2i), 2), s2i), 2)
U <- t(chol(K)) %*% matrix(rnorm(2 * q), q, 2) %*% chol(Sig)
# Causal SNPs: intermediate frequency, middle of chromosomes 2 and 4
pick <- function(cc) { i <- which(snp_map$chr == cc & abs(colMeans(M) / 2 - 0.4) < 0.1)
  i[which.min(abs(snp_map$pos[i] - 250))] }
j_dir <- pick(2); j_ind <- pick(4)
x_dir <- as.numeric(Zd %*% Mc[, j_dir]); x_ind <- as.numeric(Zn %*% Mc[, j_ind])
b_dir <- sqrt(0.08 / var(x_dir)); b_ind <- sqrt(0.06 / var(x_ind))   # 4 % and 3 % of Vp ~ 2
block_eff <- setNames(rnorm(nlevels(d$block), 0, sqrt(0.2)), levels(d$block))
d$y <- 10 + as.numeric(Zd %*% U[, 1]) + as.numeric(Zn %*% U[, 2]) + b_dir * x_dir + b_ind * x_ind +
  block_eff[as.character(d$block)] + rnorm(n, 0, 1)
cat(sprintf("Causal SNP direct %s (chr 2), indirect %s (chr 4)\n", snp_map$snp[j_dir], snp_map$snp[j_ind]))

build_model <- function(K, X_sup = NULL) {
  X <- cbind(stats::model.matrix(~ rep, d), X_sup)
  rx_model(d$y, X, list(rx_term("dge_ige", list(Zd, Zn), K = K, struct = "us"),
                        rx_term("block", d$block)), rx_residual())
}
incid <- c(dir = "dge_ige:1", ind = "dge_ige:2")
tests <- c("dir", "ind", "dir+ind", "ind|dir", "dir|ind")
fa <- list(backend = "cpu", hessian = FALSE)

# ------------------------------------------------------------------------------
# Scans
# ------------------------------------------------------------------------------
t0 <- Sys.time()
fit0 <- do.call(rx_fit, c(list(build_model(K), verbose = FALSE), fa))
cat(sprintf("full null: logLik %.3f, %.1f s\n", fit0$logLik, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
print(summary(fit0)$varcomp)
t0 <- Sys.time()
r_std <- rx_scan(fit0, build_model(K), Mc, map = snp_map, incidences = incid, tests = tests,
                 maf = NULL, verbose = FALSE)
cat(sprintf("standard scan: %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
t0 <- Sys.time()
r_loco <- rx_scan_loco(build_model, K_of, Mc, snp_map, incidences = incid, tests = tests,
                       fit_args = fa, theta_init = fit0$theta)
cat(sprintf("LOCO scan: %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
ex_dir <- rx_scan_exact(r_loco, "dir", build_model, K_of, Mc, snp_map, incid, top = 8, fit_args = fa,
                        theta_init = fit0$theta)
ex_ind <- rx_scan_exact(r_loco, "ind", build_model, K_of, Mc, snp_map, incid, top = 8, fit_args = fa,
                        theta_init = fit0$theta)

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
threshold <- 0.05 / nrow(r_std)
rank_stats <- function(res, snp, test) { p <- res[[paste0("p_", test)]]
  c(mlog10p = -log10(p[res$snp == snp]), rank = rank(p)[res$snp == snp]) }
tab <- do.call(rbind, lapply(list(c(snp_map$snp[j_dir], "dir"), c(snp_map$snp[j_ind], "ind")), function(z) {
  ex <- if (z[2] == "dir") ex_dir else ex_ind
  data.frame(snp_causal = z[1], test = z[2],
             mlog10p_standard = rank_stats(r_std, z[1], z[2])[1], rank_standard = rank_stats(r_std, z[1], z[2])[2],
             mlog10p_loco = rank_stats(r_loco, z[1], z[2])[1], rank_loco = rank_stats(r_loco, z[1], z[2])[2],
             mlog10p_exact = if (z[1] %in% ex$snp) -log10(ex$p_exact[ex$snp == z[1]]) else NA)
}))
cat(sprintf("\nBonferroni threshold: -log10 p = %.2f\n", -log10(threshold)))
print(tab, digits = 3, row.names = FALSE)
cat("\nInflation lambda: standard\n"); print(round(.rx_scan_lambda(r_std, tests), 3))
cat("LOCO\n"); print(round(attr(r_loco, "lambda"), 3))
cat("\nRefitted top hits (dir)\n"); print(ex_dir[, c("snp", "chr", "pos", "p_fixed_V", "p_exact", "beta_dir", "se_dir", "converged")], digits = 3)
cat("\nRefitted top hits (ind)\n"); print(ex_ind[, c("snp", "chr", "pos", "p_fixed_V", "p_exact", "beta_ind", "se_ind", "converged")], digits = 3)
cat(sprintf("\nTrue effects: beta_dir = %.3f, beta_ind = %.3f\n", b_dir, b_ind))
write.csv(r_std, file.path(out_dir, "scan_standard.csv"), row.names = FALSE)
write.csv(r_loco, file.path(out_dir, "scan_loco.csv"), row.names = FALSE)
write.csv(rbind(cbind(test = "dir", ex_dir[, c("snp", "chr", "pos", "p_fixed_V", "p_exact")]),
                cbind(test = "ind", ex_ind[, c("snp", "chr", "pos", "p_fixed_V", "p_exact")])),
          file.path(out_dir, "exact_top_hits.csv"), row.names = FALSE)

# Manhattan
theme_set(theme_minimal(base_size = 10) +
          theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10)))
manhattan <- function(test, title) {
  df <- rbind(data.frame(method = "standard", r_std[, c("snp", paste0("p_", test))]),
              data.frame(method = "LOCO", r_loco[match(r_std$snp, r_loco$snp), c("snp", paste0("p_", test))]))
  names(df)[3] <- "p"
  df$x <- match(df$snp, snp_map$snp); df$chr <- snp_map$chr[df$x]
  df$method <- factor(df$method, levels = c("standard", "LOCO"))
  causal <- snp_map$snp[if (test == "dir") j_dir else j_ind]
  ex <- if (test == "dir") ex_dir else ex_ind
  exd <- data.frame(x = match(ex$snp, snp_map$snp), p = ex$p_exact, method = factor("LOCO", levels = c("standard", "LOCO")))
  ggplot(df, aes(x, -log10(p), colour = factor(chr %% 2))) + geom_point(size = 0.6) +
    geom_hline(yintercept = -log10(threshold), linetype = 2, colour = "grey40") +
    geom_vline(xintercept = match(causal, snp_map$snp), colour = "#d6604d", alpha = 0.4) +
    geom_point(data = exd, aes(x, -log10(p)), inherit.aes = FALSE, shape = 1, size = 2.2, colour = "black") +
    facet_wrap(~ method, nrow = 1) + scale_colour_manual(values = c("#2166ac", "#67a9cf"), guide = "none") +
    labs(title = title, x = "SNP (chromosomes 1 to 5)", y = "-log10 p")
}
fig <- manhattan("dir", "a  Direct test (dir); red: causal SNP; circles: exact refit of the 8 top hits") /
  manhattan("ind", "b  Indirect test (ind, dose of the neighbours)")
ggsave(file.path(out_dir, "fig_gwas_dge_ige.png"), fig, width = 12, height = 7, dpi = 150, bg = "white")
cat("\nOutputs in", normalizePath(out_dir), "\n")
