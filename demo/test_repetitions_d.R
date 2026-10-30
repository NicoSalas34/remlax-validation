# ==============================================================================
# Identifiabilite de la dilution d selon le nombre de repetitions
#
# Dispositif de la demo agrandi : 150 genotypes, chaque repetition occupe 10
# rangs x 15 colonnes (un genotype par repetition), avec 3, 6 ou 12 repetitions.
# Sous-blocs disjoints heterogenes (bandes de 1 a 5 rangs sur 5 colonnes, le
# meilleur cas de test_identifiabilite_d.R), 10 % de parcelles vides, noyau
# puissance r = 3, lambda = 1,5, d = 0,5, IGE 0,4 de variance phenotypique,
# r(DGE, IGE) = -0,5. La tendance spatiale de la demo est etiree sur la longueur
# du champ. Profil REML du modele B (DGE + IGE independants, bloc, spl2d) sur
# lambda in {1 ; 1,5 ; 2 ; 2,5} et d in {0 ; 0,25 ; ... ; 1,25}.
#
#   Rscript test_repetitions_d.R liste
#   Rscript test_repetitions_d.R tache <i> <dossier>   (6 ajustements : un lambda, toute la grille en d)
#   Rscript test_repetitions_d.R bilan <dossier>
# ==============================================================================

a0 <- commandArgs(trailingOnly = TRUE)
mode <- if (length(a0)) a0[1] else "liste"
vraiR <- list(r = 3L, lambda = 1.5, d = 0.5, dge = 1.0, ige = 0.4, r_dge_ige = -0.5,
              bloc = 0.3, spatial = 1.2, residuelle = 0.8)
grille_l <- c(1, 1.5, 2, 2.5)
grille_d <- seq(0, 1.25, 0.25)
taches <- expand.grid(lambda = grille_l, graine = 101:104, n_rep = c(3L, 6L, 12L))
taches <- taches[, c("n_rep", "graine", "lambda")]
if (mode == "liste") { cat(nrow(taches), "taches\n"); quit(save = "no") }

if (mode == "tache") {
  suppressPackageStartupMessages(library(remlax))
  i <- as.integer(a0[2]); dossier <- a0[3]
  dir.create(dossier, showWarnings = FALSE, recursive = TRUE)
  t <- taches[i, ]
  n_gen <- 150L; n_col <- 15L; n_rang <- 10L * t$n_rep
  noms_gen <- sprintf("G%03d", seq_len(n_gen))

  # Parente : marqueurs de la demo (graine 5), deux sous-populations
  set.seed(5L)
  p0 <- runif(1500L, 0.05, 0.5); pop <- rep(1:2, length.out = n_gen)
  P <- sapply(p0, function(p) rbeta(2, p * 9, (1 - p) * 9))
  M <- matrix(rbinom(n_gen * 1500L, 2L, P[pop, ]), n_gen, 1500L); rownames(M) <- noms_gen
  K <- rx_grm(M, blend = 0.02, coding = "count")

  # Dispositif et donnees
  set.seed(t$graine)
  d <- expand.grid(col = seq_len(n_col), rang = seq_len(n_rang))
  d$rep <- factor(sprintf("R%02d", (d$rang - 1L) %/% 10L + 1L))
  d$bloc <- factor(sprintf("%s_r%03d_b%d", d$rep, d$rang, (d$col - 1L) %/% 5L + 1L))
  d$parcelle <- sprintf("P%04d", seq_len(nrow(d)))
  d$gid <- factor(NA, levels = noms_gen)
  for (r in levels(d$rep)) d$gid[d$rep == r] <- sample(noms_gen)
  hauteurs <- rep(1:5, length.out = ceiling(n_rang / 3))
  bandes <- rep(seq_along(hauteurs), times = hauteurs)[seq_len(n_rang)]
  d$sb <- factor(paste((d$col - 1L) %/% 5L, bandes[d$rang]))
  d <- d[sort(sample(nrow(d), round(0.9 * nrow(d)))), ]; rownames(d) <- NULL
  d$bloc <- droplevels(d$bloc)
  # Tendance de la demo, etiree sur la longueur du champ
  tb <- function(x, y) { y <- 1 + (y - 1) * 29 / (n_rang - 1)
    0.9 * (y - 15.5) / 14.5 - 1.6 * exp(-((x - 11)^2 / 10 + (y - 9)^2 / 30)) + 0.35 * sin(2 * pi * x / 10) }
  br <- tb(d$col, d$rang); d$spatial <- (br - mean(br)) / sd(br) * sqrt(vraiR$spatial)

  voisinage <- function(lambda, dd) {
    nb <- rx_neighbourhood(coord = cbind(d$rang, d$col), group = rep("A", nrow(d)), block = d$sb,
                           level = d$gid, id = d$parcelle, rank = vraiR$r, reach = as.numeric(lambda),
                           dilution = as.numeric(dd), window = "euclidean", kernel = "power")
    Z <- nb$level[["A<-A"]]
    Pm <- Matrix::sparseMatrix(i = match(colnames(Z), noms_gen), j = seq_len(ncol(Z)), x = 1,
                               dims = c(n_gen, ncol(Z)))
    Z <- Z %*% Matrix::t(Pm); dimnames(Z) <- list(NULL, noms_gen)
    methods::as(Z, "CsparseMatrix")
  }
  Zd <- Matrix::sparseMatrix(i = seq_len(nrow(d)), j = as.integer(d$gid), x = 1,
                             dims = c(nrow(d), n_gen), dimnames = list(NULL, noms_gen))
  Zv <- voisinage(vraiR$lambda, vraiR$d)
  s2i <- vraiR$ige / rx_exposure(Zd, Zv, K = K)$k
  Sig <- matrix(c(vraiR$dge, rep(vraiR$r_dge_ige * sqrt(vraiR$dge * s2i), 2), s2i), 2)
  U <- t(chol(K)) %*% matrix(rnorm(2 * n_gen), n_gen, 2) %*% chol(Sig)
  b <- setNames(rnorm(nlevels(d$bloc), 0, sqrt(vraiR$bloc)), levels(d$bloc))
  d$y <- 6 + d$spatial + U[as.integer(d$gid), 1] + as.numeric(Zv %*% U[, 2]) +
    b[as.character(d$bloc)] + rnorm(nrow(d), 0, sqrt(vraiR$residuelle))
  nseg_r <- round(10 * t$n_rep / 3)

  sorties <- lapply(grille_d, function(dd) {
    Z <- voisinage(t$lambda, dd)
    random <- stats::as.formula(sprintf(
      "~ vm(gid, K = K) + mm(Z, K = K, name = 'ige') + bloc + spl2d(col, rang, nseg = c(7, %d))", nseg_r),
      env = environment())
    t0 <- Sys.time()
    f <- tryCatch(rx_reml(y ~ rep, random = random, data = d, backend = "cpu", verbose = FALSE,
                          hessian = FALSE, blups = FALSE, pev = FALSE), error = function(e) conditionMessage(e))
    if (is.character(f)) return(data.frame(t, d = dd, n = nrow(d), logLik = NA, conv = NA, secondes = NA,
                                           erreur = substr(f, 1, 150)))
    data.frame(t, d = dd, n = nrow(d), logLik = f$logLik,
               conv = isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement)),
               secondes = as.numeric(difftime(Sys.time(), t0, units = "secs")), erreur = "")
  })
  res <- do.call(rbind, sorties)
  print(res)
  utils::write.csv(res, file.path(dossier, sprintf("tache_%03d.csv", i)), row.names = FALSE)
  quit(save = "no")
}

if (mode == "bilan") {
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  dossier <- a0[2]
  res <- do.call(rbind, lapply(list.files(dossier, "^tache_.*\\.csv$", full.names = TRUE), utils::read.csv))
  cat(sprintf("%d ajustements, %d echec(s), %d non converge(s)\n", nrow(res), sum(is.na(res$logLik)),
              sum(!res$conv, na.rm = TRUE)))
  bilan <- do.call(rbind, lapply(split(res, list(res$n_rep, res$graine), drop = TRUE), function(q) {
    pd <- stats::aggregate(logLik ~ d, q, max); pd <- pd[order(pd$d), ]
    j <- which.max(pd$logLik); jj <- min(max(j, 2L), nrow(pd) - 1L)
    cf <- stats::coef(stats::lm(pd$logLik[(jj - 1):(jj + 1)] ~ poly(pd$d[(jj - 1):(jj + 1)], 2, raw = TRUE)))
    d_hat <- if (cf[3] < 0) min(max(-cf[2] / (2 * cf[3]), 0), 1.25) else pd$d[j]
    fine <- seq(0, 1.25, length.out = 501)
    yf <- stats::spline(pd$d, pd$logLik - max(pd$logLik), xout = fine)$y
    dedans <- fine[yf >= -stats::qchisq(0.95, 1) / 2]
    data.frame(n_rep = q$n_rep[1], n = q$n[1], graine = q$graine[1], d_chapeau = d_hat,
               lambda_chapeau = q$lambda[which.max(q$logLik)], ic_bas = min(dedans), ic_haut = max(dedans),
               et_d = if (cf[3] < 0) 1 / sqrt(-2 * cf[3]) else NA_real_,
               chute_au_vrai = pd$logLik[pd$d == vraiR$d] - max(pd$logLik), chute_a_0 = pd$logLik[1] - max(pd$logLik))
  }))
  bilan <- bilan[order(bilan$n_rep, bilan$graine), ]
  print(bilan, digits = 3, row.names = FALSE)
  cat("\nMoyennes par nombre de repetitions\n")
  print(stats::aggregate(cbind(n, d_chapeau, largeur_ic = ic_haut - ic_bas, et_d, chute_a_0) ~ n_rep,
                         bilan, mean), digits = 3)
  utils::write.csv(bilan, file.path(dossier, "bilan_repetitions.csv"), row.names = FALSE)
  pr <- do.call(rbind, lapply(split(res, list(res$n_rep, res$graine), drop = TRUE), function(q) {
    pd <- stats::aggregate(logLik ~ d, q, max); pd$dll <- pd$logLik - max(pd$logLik)
    cbind(pd, n_rep = q$n_rep[1], graine = q$graine[1]) }))
  theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10)))
  p1 <- ggplot(pr, aes(d, dll, group = graine)) + geom_hline(yintercept = -stats::qchisq(0.95, 1) / 2,
                                                             linetype = 3, colour = "grey50") +
    geom_vline(xintercept = vraiR$d, linetype = 2, colour = "grey50") +
    geom_line(colour = "#1b7837") + geom_point(colour = "#1b7837", size = 1) +
    facet_wrap(~ sprintf("%d repetitions", n_rep), nrow = 1) + coord_cartesian(ylim = c(-8, 0)) +
    labs(title = "a  Profil de vraisemblance en d (maximum sur lambda), une courbe par graine",
         x = "d", y = "delta logLik") + theme(panel.grid.major.y = element_line(colour = "grey90"))
  p2 <- ggplot(bilan, aes(factor(graine), d_chapeau)) +
    geom_hline(yintercept = vraiR$d, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(ymin = ic_bas, ymax = ic_haut), colour = "#1b7837", size = 0.3) +
    facet_wrap(~ sprintf("%d repetitions (n = %d)", n_rep, n), nrow = 1) +
    labs(title = "b  d estime et intervalle a 95 % (tirets : vrai d = 0,5)", x = "graine", y = "d") +
    theme(panel.grid.major.y = element_line(colour = "grey90"))
  ggsave(file.path(dossier, "fig_repetitions_d.png"), p1 / p2 +
           plot_annotation(title = "Dilution d et nombre de repetitions (sous-blocs heterogenes, noyau puissance)"),
         width = 12, height = 7.5, dpi = 150, bg = "white")
  cat("Figure :", normalizePath(file.path(dossier, "fig_repetitions_d.png")), "\n")
}
