# ==============================================================================
# Voisinage sans rayon a estimer : toutes les voisines dans une fenetre large,
# poids exponentiels, dilution
#
# Modele ajuste : w_ij = exp(-lambda (delta_ij - 1)) / n_i^d pour toutes les
# voisines a distance <= R = 6 (dans le meme sous-bloc), n_i = nombre de ces
# voisines. Seuls lambda et d restent a estimer ; la portee passe par lambda.
#
# Donnees : exactement les jeux de test_identifiabilite_d.R (memes graines, memes
# cas), simules avec un rayon r = 3, d = 0,5, noyau puissance (lambda 1,5) ou
# exponentiel (lambda 0,8). Le modele sans rayon est donc mal specifie sur la
# coupure a 3 et sur le n_i de la dilution. Les logLik se comparent directement a
# celles du modele a r = 3 (profils.csv du test d'identifiabilite).
#
#   Rscript test_sans_rayon.R liste
#   Rscript test_sans_rayon.R tache <i> <dossier>   (un cas, une graine, un lambda : 6 ajustements en d)
#   Rscript test_sans_rayon.R bilan <dossier> [profils_r3.csv]
# ==============================================================================

a0 <- commandArgs(trailingOnly = TRUE)
mode <- if (length(a0)) a0[1] else "liste"
R_FENETRE <- 6L
vraiD <- list(r = 3L, d = 0.5, lambda = c(puissance = 1.5, exponentiel = 0.8),
              dge = 1.0, ige = 0.4, r_dge_ige = -0.5, bloc = 0.3, residuelle = 0.8)
noyaux <- c(puissance = "power", exponentiel = "exponential")
grille_l <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3)
grille_d <- seq(0, 1.25, 0.25)
taches <- expand.grid(lambda = grille_l, graine = 101:104, noyau = names(noyaux),
                      sous_blocs = c("continu", "5x5", "heterogene"), stringsAsFactors = FALSE)
taches <- taches[, c("sous_blocs", "noyau", "graine", "lambda")]
if (mode == "liste") { cat(nrow(taches), "taches\n"); quit(save = "no") }

if (mode == "tache") {
  suppressPackageStartupMessages(library(remlax))
  i <- as.integer(a0[2]); dossier <- a0[3]
  dir.create(dossier, showWarnings = FALSE, recursive = TRUE)
  l <- readLines("demo_remlax_spatial_cpu.R")
  pre <- l[1:(grep("^# 3. Ajustements", l) - 2)]
  pre <- sub("^args <- commandArgs.*", "args <- character(0)", pre)
  pre <- pre[!grepl("library\\((ggplot2|patchwork)\\)", pre)]
  sortie <- tempdir()
  eval(parse(text = pre))

  bandes <- rep(seq_len(10), times = rep(1:5, 2))
  sous_bloc <- function(d2, type) switch(type,
    continu = factor(rep(1L, nrow(d2))),
    "5x5" = factor(paste((d2$col - 1L) %/% 5L, (d2$rang - 1L) %/% 5L)),
    heterogene = factor(paste((d2$col - 1L) %/% 5L, bandes[d2$rang])))
  voisinage <- function(d2, r, lambda, dd, noyau) {
    nb <- rx_neighbourhood(coord = cbind(d2$rang, d2$col), group = rep("A", nrow(d2)), block = d2$sb,
                           level = d2$gid, id = d2$parcelle, rank = as.integer(r),
                           reach = as.numeric(lambda), dilution = as.numeric(dd), window = "euclidean",
                           kernel = noyaux[[noyau]], normalise = FALSE)
    Z <- nb$level[["A<-A"]]
    P <- Matrix::sparseMatrix(i = match(colnames(Z), noms_gen), j = seq_len(ncol(Z)), x = 1,
                              dims = c(n_gen, ncol(Z)))
    Z <- Z %*% Matrix::t(P)
    dimnames(Z) <- list(NULL, noms_gen)
    methods::as(Z, "CsparseMatrix")
  }
  # Simulation identique a test_identifiabilite_d.R (meme ordre des tirages)
  simuler <- function(graine, sous_blocs, noyau) {
    set.seed(graine)
    d2 <- d[sort(sample(nrow(d), round(0.9 * nrow(d)))), c("col", "rang", "rep", "bloc", "parcelle", "gid")]
    rownames(d2) <- NULL
    d2$bloc <- droplevels(d2$bloc)
    d2$sb <- sous_bloc(d2, sous_blocs)
    Zv <- voisinage(d2, vraiD$r, vraiD$lambda[[noyau]], vraiD$d, noyau)
    Zd <- Matrix::sparseMatrix(i = seq_len(nrow(d2)), j = as.integer(d2$gid), x = 1,
                               dims = c(nrow(d2), n_gen), dimnames = list(NULL, noms_gen))
    k <- rx_exposure(Zd, Zv, K = K)$k
    s2i <- vraiD$ige / k
    Sig <- matrix(c(vraiD$dge, rep(vraiD$r_dge_ige * sqrt(vraiD$dge * s2i), 2), s2i), 2)
    U <- t(chol(K)) %*% matrix(rnorm(2 * n_gen), n_gen, 2) %*% chol(Sig)
    b <- setNames(rnorm(nlevels(d2$bloc), 0, sqrt(vraiD$bloc)), levels(d2$bloc))
    d2$y <- 6 + 0.2 * (d2$rep == "R2") - 0.1 * (d2$rep == "R3") + tendance_vraie(d2$col, d2$rang) +
      U[as.integer(d2$gid), 1] + as.numeric(Zv %*% U[, 2]) + b[as.character(d2$bloc)] +
      rnorm(nrow(d2), 0, sqrt(vraiD$residuelle))
    d2
  }

  t <- taches[i, ]
  d2 <- simuler(t$graine, t$sous_blocs, t$noyau)
  res <- do.call(rbind, lapply(grille_d, function(dd) {
    Z <- voisinage(d2, R_FENETRE, t$lambda, dd, "exponentiel")
    random <- stats::as.formula(
      "~ vm(gid, K = K) + mm(Z, K = K, name = 'ige') + bloc + spl2d(col, rang, nseg = c(7, 10))",
      env = environment())
    f <- tryCatch(rx_reml(y ~ rep, random = random, data = d2, backend = "cpu", verbose = FALSE,
                          hessian = FALSE, blups = FALSE, pev = FALSE), error = function(e) conditionMessage(e))
    if (is.character(f)) return(data.frame(t, d = dd, logLik = NA, conv = NA, erreur = substr(f, 1, 150)))
    data.frame(t, d = dd, logLik = f$logLik,
               conv = isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement)), erreur = "")
  }))
  print(res)
  utils::write.csv(res, file.path(dossier, sprintf("tache_%03d.csv", i)), row.names = FALSE)
  quit(save = "no")
}

if (mode == "bilan") {
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  dossier <- a0[2]
  f_r3 <- if (length(a0) >= 3) a0[3] else NA
  res <- do.call(rbind, lapply(list.files(dossier, "^tache_.*\\.csv$", full.names = TRUE), utils::read.csv))
  cat(sprintf("%d ajustements, %d echec(s), %d non converge(s)\n", nrow(res), sum(is.na(res$logLik)),
              sum(!res$conv, na.rm = TRUE)))
  nom_cas <- function(sb, nk) paste(c(continu = "champ continu", "5x5" = "sous-blocs 5 x 5",
                                      heterogene = "sous-blocs heterogenes")[sb], nk, sep = ", ")
  res$cas <- nom_cas(res$sous_blocs, res$noyau)
  r3 <- if (!is.na(f_r3) && file.exists(f_r3)) utils::read.csv(f_r3) else NULL
  parabole <- function(x, y) {
    j <- which.max(y); jj <- min(max(j, 2L), length(y) - 1L)
    cf <- stats::coef(stats::lm(y[(jj - 1):(jj + 1)] ~ poly(x[(jj - 1):(jj + 1)], 2, raw = TRUE)))
    list(sommet = if (cf[3] < 0) min(max(-cf[2] / (2 * cf[3]), min(x)), max(x)) else x[j],
         et = if (cf[3] < 0) 1 / sqrt(-2 * cf[3]) else NA_real_)
  }
  bilan <- do.call(rbind, lapply(split(res, list(res$cas, res$graine), drop = TRUE), function(q) {
    pd <- stats::aggregate(logLik ~ d, q, max); pd <- pd[order(pd$d), ]
    pl <- stats::aggregate(logLik ~ lambda, q, max); pl <- pl[order(pl$lambda), ]
    sd_ <- parabole(pd$d, pd$logLik); sl <- parabole(pl$lambda, pl$logLik)
    fine <- seq(0, 1.25, length.out = 501)
    yf <- stats::spline(pd$d, pd$logLik - max(pd$logLik), xout = fine)$y
    dedans <- fine[yf >= -stats::qchisq(0.95, 1) / 2]
    fl <- seq(0, 3, length.out = 601)
    yl <- stats::spline(pl$lambda, pl$logLik - max(pl$logLik), xout = fl)$y
    dl <- fl[yl >= -stats::qchisq(0.95, 1) / 2]
    ll_r3 <- if (is.null(r3)) NA_real_ else {
      m <- r3[r3$sous_blocs == q$sous_blocs[1] & r3$noyau == q$noyau[1] & r3$graine == q$graine[1], ]
      if (nrow(m)) max(m$logLik, na.rm = TRUE) else NA_real_ }
    data.frame(cas = q$cas[1], noyau = q$noyau[1], graine = q$graine[1],
               d_chapeau = sd_$sommet, d_bas = min(dedans), d_haut = max(dedans), et_d = sd_$et,
               lambda_chapeau = sl$sommet, l_bas = min(dl), l_haut = max(dl),
               chute_d_vrai = pd$logLik[pd$d == vraiD$d] - max(pd$logLik),
               logLik_max = max(q$logLik), logLik_r3 = ll_r3, delta_vs_r3 = max(q$logLik) - ll_r3)
  }))
  bilan <- bilan[order(bilan$cas, bilan$graine), ]
  print(bilan, digits = 3, row.names = FALSE)
  cat("\nMoyennes par cas\n")
  print(stats::aggregate(cbind(d_chapeau, largeur_d = d_haut - d_bas, et_d, lambda_chapeau,
                               largeur_l = l_haut - l_bas, delta_vs_r3) ~ cas, bilan, mean), digits = 3)
  utils::write.csv(bilan, file.path(dossier, "bilan_sans_rayon.csv"), row.names = FALSE)

  theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10)))
  q1 <- res[res$graine == 101, ]
  q1 <- do.call(rbind, lapply(split(q1, q1$cas), function(q) { q$dll <- pmax(q$logLik - max(q$logLik), -12); q }))
  p_h <- ggplot(q1, aes(lambda, d)) + geom_tile(aes(fill = dll)) +
    geom_contour(aes(z = dll), breaks = -stats::qchisq(0.95, 2) / 2, colour = "black", linewidth = 0.5) +
    geom_hline(yintercept = vraiD$d, linetype = 2, colour = "white") +
    facet_wrap(~ cas, nrow = 1) + scale_fill_viridis_c(name = "delta\nlogLik") + coord_cartesian(expand = FALSE) +
    labs(title = "a  Profil (lambda, d) du modele sans rayon, graine 101 ; contour = region a 95 %, tirets : vrai d",
         x = "lambda (exponentiel, fenetre R = 6)", y = "d")
  pdl <- do.call(rbind, lapply(split(res, list(res$cas, res$graine), drop = TRUE), function(q) {
    pd <- stats::aggregate(logLik ~ d, q, max); pd$dll <- pd$logLik - max(pd$logLik)
    cbind(pd, cas = q$cas[1], graine = q$graine[1]) }))
  p_d <- ggplot(pdl, aes(d, dll, group = graine)) +
    geom_hline(yintercept = -stats::qchisq(0.95, 1) / 2, linetype = 3, colour = "grey50") +
    geom_vline(xintercept = vraiD$d, linetype = 2, colour = "grey50") +
    geom_line(colour = "#1b7837") + geom_point(colour = "#1b7837", size = 0.8) +
    facet_wrap(~ cas, nrow = 1) + coord_cartesian(ylim = c(-6, 0)) +
    labs(title = "b  Profil en d (maximum sur lambda), une courbe par graine", x = "d", y = "delta logLik") +
    theme(panel.grid.major.y = element_line(colour = "grey90"))
  p_c <- ggplot(bilan, aes(factor(graine), delta_vs_r3)) + geom_hline(yintercept = 0, colour = "grey50") +
    geom_col(fill = "#2166ac", width = 0.6) + facet_wrap(~ cas, nrow = 1) +
    labs(title = "c  logLik maximale : modele sans rayon moins modele a r = 3 (profil lambda, d du test precedent)",
         x = "graine", y = "delta logLik") + theme(panel.grid.major.y = element_line(colour = "grey90"))
  ggsave(file.path(dossier, "fig_sans_rayon.png"), p_h / p_d / p_c +
           plot_annotation(title = "Voisinage sans rayon : noyau exponentiel sur fenetre R = 6, lambda et d estimes"),
         width = 16, height = 11, dpi = 150, bg = "white")
  cat("Figure :", normalizePath(file.path(dossier, "fig_sans_rayon.png")), "\n")
}
