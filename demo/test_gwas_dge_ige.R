# ==============================================================================
# GWAS en une etape avec remlax sur un jeu DGE + IGE simple : rx_scan, LOCO
# (rx_scan_loco) et reajustement exact des tetes (rx_scan_exact)
#
# 300 genotypes, 4 repetitions, champ de 80 rangs x 15 colonnes, blocs de 5.
# Marqueurs : 5 chromosomes x 500 SNP, desequilibre de liaison le long du
# chromosome (haplotypes latents AR(1)) et deux sous-populations. Chaque parcelle
# recoit l'IGE de ses 4 voisines de rang et de colonne (incidence normalisee L2).
# Verite : DGE et IGE polygeniques correles (r = -0,3), plus deux SNP causaux,
# un a effet DIRECT sur le chromosome 2 (4 % de la variance phenotypique) et un a
# effet INDIRECT (dose des voisines) sur le chromosome 4 (3 %).
#
#   Rscript demo/test_gwas_dge_ige.R [dossier_sortie]
#
# remlax est le paquet installe, charge par R/charger_remlax.R (chemin relatif a
# ce script), qui expose aussi les fonctions internes.
# ==============================================================================

ici <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(ici, "..", "R", "charger_remlax.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
sortie <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "sortie_gwas_dge_ige"
dir.create(sortie, showWarnings = FALSE, recursive = TRUE)
set.seed(11)

# ------------------------------------------------------------------------------
# Marqueurs
# ------------------------------------------------------------------------------
q <- 300L; n_chr <- 5L; m_chr <- 500L
noms_gen <- sprintf("G%03d", seq_len(q))
pop <- rep(1:2, length.out = q)
haplos <- function(m) {                   # 2q haplotypes latents AR(1), rho = 0,9
  z <- matrix(0, 2 * q, m); z[, 1] <- rnorm(2 * q)
  for (j in 2:m) z[, j] <- 0.9 * z[, j - 1] + sqrt(1 - 0.81) * rnorm(2 * q)
  z
}
M <- do.call(cbind, lapply(seq_len(n_chr), function(cc) {
  p0 <- runif(m_chr, 0.1, 0.5)
  pp <- sapply(p0, function(p) rbeta(2, p * 19, (1 - p) * 19))       # Fst ~ 0,05
  z <- haplos(m_chr)
  a <- pnorm(z) < pp[rep(pop, 2), ]                                  # allele par haplotype
  a[seq_len(q), ] + a[q + seq_len(q), ]
}))
dimnames(M) <- list(noms_gen, sprintf("snp%d_%03d", rep(seq_len(n_chr), each = m_chr), rep(seq_len(m_chr), n_chr)))
carte <- data.frame(snp = colnames(M), chr = rep(seq_len(n_chr), each = m_chr), pos = rep(seq_len(m_chr), n_chr))
f <- colMeans(M) / 2
garder <- pmin(f, 1 - f) >= 0.05
M <- M[, garder]; carte <- carte[garder, ]
Mc <- sweep(M, 2, colMeans(M))                                       # doses centrees
cat(sprintf("%d SNP apres MAF >= 0,05\n", ncol(M)))

# Parente : produit croise des doses centrees, diagonale moyenne ramenee a 1,
# melangee a 2 % d'identite. La meme fonction sert au nul complet et aux nuls LOCO.
K_de <- function(W) {
  G <- tcrossprod(W); G <- G / mean(diag(G))
  G <- 0.98 * G + 0.02 * diag(nrow(G)); dimnames(G) <- list(noms_gen, noms_gen); G
}
K <- K_de(Mc)

# ------------------------------------------------------------------------------
# Dispositif, voisinage et phenotypes
# ------------------------------------------------------------------------------
d <- expand.grid(col = 1:15, rang = 1:80)
d$rep <- factor(sprintf("R%d", (d$rang - 1L) %/% 20L + 1L))
d$bloc <- factor(sprintf("%s_r%02d_b%d", d$rep, d$rang, (d$col - 1L) %/% 5L + 1L))
d$parcelle <- sprintf("P%04d", seq_len(nrow(d)))
d$gid <- factor(NA, levels = noms_gen)
for (r in levels(d$rep)) d$gid[d$rep == r] <- sample(rep(noms_gen, length.out = sum(d$rep == r)))
n <- nrow(d)
nb <- rx_neighbourhood(cbind(d$rang, d$col), group = rep("A", n), level = d$gid, id = d$parcelle,
                       rank = 1, window = "euclidean", kernel = "none", normalise = TRUE)
Zn <- nb$level[["A<-A"]][, noms_gen]; dimnames(Zn) <- list(NULL, noms_gen)
Zd <- Matrix::sparseMatrix(i = seq_len(n), j = as.integer(d$gid), x = 1, dims = c(n, q),
                           dimnames = list(NULL, noms_gen))
k_ige <- rx_exposure(Zd, Zn, K = K)$k

# Polygenique : DGE 0,5 ; IGE 0,15 de variance phenotypique ; r = -0,3
s2i <- 0.15 / k_ige
Sig <- matrix(c(0.5, rep(-0.3 * sqrt(0.5 * s2i), 2), s2i), 2)
U <- t(chol(K)) %*% matrix(rnorm(2 * q), q, 2) %*% chol(Sig)
# SNP causaux : frequence moyenne, milieu des chromosomes 2 et 4
choisir <- function(cc) { i <- which(carte$chr == cc & abs(colMeans(M) / 2 - 0.4) < 0.1)
  i[which.min(abs(carte$pos[i] - 250))] }
j_dir <- choisir(2); j_ind <- choisir(4)
x_dir <- as.numeric(Zd %*% Mc[, j_dir]); x_ind <- as.numeric(Zn %*% Mc[, j_ind])
b_dir <- sqrt(0.08 / var(x_dir)); b_ind <- sqrt(0.06 / var(x_ind))   # 4 % et 3 % de Vp ~ 2
bloc_eff <- setNames(rnorm(nlevels(d$bloc), 0, sqrt(0.2)), levels(d$bloc))
d$y <- 10 + as.numeric(Zd %*% U[, 1]) + as.numeric(Zn %*% U[, 2]) + b_dir * x_dir + b_ind * x_ind +
  bloc_eff[as.character(d$bloc)] + rnorm(n, 0, 1)
cat(sprintf("SNP causal direct %s (chr 2), indirect %s (chr 4)\n", carte$snp[j_dir], carte$snp[j_ind]))

construire <- function(K, X_sup = NULL) {
  X <- cbind(stats::model.matrix(~ rep, d), X_sup)
  rx_model(d$y, X, list(rx_term("dge_ige", list(Zd, Zn), K = K, struct = "us"),
                        rx_term("bloc", d$bloc)), rx_residual())
}
incid <- c(dir = "dge_ige:1", ind = "dge_ige:2")
tests <- c("dir", "ind", "dir+ind", "ind|dir", "dir|ind")
fa <- list(backend = "cpu", hessian = FALSE)

# ------------------------------------------------------------------------------
# Scans
# ------------------------------------------------------------------------------
t0 <- Sys.time()
fit0 <- do.call(rx_fit, c(list(construire(K), verbose = FALSE), fa))
cat(sprintf("nul complet : logLik %.3f, %.1f s\n", fit0$logLik, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
print(summary(fit0)$varcomp)
t0 <- Sys.time()
r_std <- rx_scan(fit0, construire(K), Mc, carte = carte, incidences = incid, tests = tests,
                 maf = NULL, verbose = FALSE)
cat(sprintf("scan standard : %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
t0 <- Sys.time()
r_loco <- rx_scan_loco(construire, K_de, Mc, carte, incidences = incid, tests = tests,
                       fit_args = fa, theta_init = fit0$theta)
cat(sprintf("scan LOCO : %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
ex_dir <- rx_scan_exact(r_loco, "dir", construire, K_de, Mc, carte, incid, top = 8, fit_args = fa,
                        theta_init = fit0$theta)
ex_ind <- rx_scan_exact(r_loco, "ind", construire, K_de, Mc, carte, incid, top = 8, fit_args = fa,
                        theta_init = fit0$theta)

# ------------------------------------------------------------------------------
# Bilan
# ------------------------------------------------------------------------------
seuil <- 0.05 / nrow(r_std)
bil <- function(res, snp, test) { p <- res[[paste0("p_", test)]]
  c(mlog10p = -log10(p[res$snp == snp]), rang = rank(p)[res$snp == snp]) }
tab <- do.call(rbind, lapply(list(c(carte$snp[j_dir], "dir"), c(carte$snp[j_ind], "ind")), function(z) {
  ex <- if (z[2] == "dir") ex_dir else ex_ind
  data.frame(snp_causal = z[1], test = z[2],
             mlog10p_standard = bil(r_std, z[1], z[2])[1], rang_standard = bil(r_std, z[1], z[2])[2],
             mlog10p_loco = bil(r_loco, z[1], z[2])[1], rang_loco = bil(r_loco, z[1], z[2])[2],
             mlog10p_exact = if (z[1] %in% ex$snp) -log10(ex$p_exact[ex$snp == z[1]]) else NA)
}))
cat(sprintf("\nSeuil de Bonferroni : -log10 p = %.2f\n", -log10(seuil)))
print(tab, digits = 3, row.names = FALSE)
cat("\nInflation lambda : standard\n"); print(round(.rx_scan_lambda(r_std, tests), 3))
cat("LOCO\n"); print(round(attr(r_loco, "lambda"), 3))
cat("\nTetes reajustees (dir)\n"); print(ex_dir[, c("snp", "chr", "pos", "p_fige", "p_exact", "beta_dir", "se_dir", "converge")], digits = 3)
cat("\nTetes reajustees (ind)\n"); print(ex_ind[, c("snp", "chr", "pos", "p_fige", "p_exact", "beta_ind", "se_ind", "converge")], digits = 3)
cat(sprintf("\nEffets vrais : beta_dir = %.3f, beta_ind = %.3f\n", b_dir, b_ind))
write.csv(r_std, file.path(sortie, "scan_standard.csv"), row.names = FALSE)
write.csv(r_loco, file.path(sortie, "scan_loco.csv"), row.names = FALSE)
write.csv(rbind(cbind(test = "dir", ex_dir[, c("snp", "chr", "pos", "p_fige", "p_exact")]),
                cbind(test = "ind", ex_ind[, c("snp", "chr", "pos", "p_fige", "p_exact")])),
          file.path(sortie, "tetes_exactes.csv"), row.names = FALSE)

# Manhattan
theme_set(theme_minimal(base_size = 10) +
          theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10)))
man <- function(test, titre) {
  df <- rbind(data.frame(methode = "standard", r_std[, c("snp", paste0("p_", test))]),
              data.frame(methode = "LOCO", r_loco[match(r_std$snp, r_loco$snp), c("snp", paste0("p_", test))]))
  names(df)[3] <- "p"
  df$x <- match(df$snp, carte$snp); df$chr <- carte$chr[df$x]
  df$methode <- factor(df$methode, levels = c("standard", "LOCO"))
  causal <- carte$snp[if (test == "dir") j_dir else j_ind]
  ex <- if (test == "dir") ex_dir else ex_ind
  exd <- data.frame(x = match(ex$snp, carte$snp), p = ex$p_exact, methode = factor("LOCO", levels = c("standard", "LOCO")))
  ggplot(df, aes(x, -log10(p), colour = factor(chr %% 2))) + geom_point(size = 0.6) +
    geom_hline(yintercept = -log10(seuil), linetype = 2, colour = "grey40") +
    geom_vline(xintercept = match(causal, carte$snp), colour = "#d6604d", alpha = 0.4) +
    geom_point(data = exd, aes(x, -log10(p)), inherit.aes = FALSE, shape = 1, size = 2.2, colour = "black") +
    facet_wrap(~ methode, nrow = 1) + scale_colour_manual(values = c("#2166ac", "#67a9cf"), guide = "none") +
    labs(title = titre, x = "SNP (chromosomes 1 a 5)", y = "-log10 p")
}
fig <- man("dir", "a  Test direct (dir) ; rouge : SNP causal ; cercles : reajustement exact des 8 tetes") /
  man("ind", "b  Test indirect (ind, dose des voisines)")
ggsave(file.path(sortie, "fig_gwas_dge_ige.png"), fig, width = 12, height = 7, dpi = 150, bg = "white")
cat("\nSorties dans", normalizePath(sortie), "\n")
