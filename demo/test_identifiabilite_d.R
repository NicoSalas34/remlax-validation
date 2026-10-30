# ==============================================================================
# Identifiabilite de la dilution d : sous-blocs disjoints et noyau exponentiel
#
# Reprend le dispositif de demo_remlax_spatial_cpu.R (memes genotypes, meme
# parente, meme tendance spatiale, 10 % de parcelles vides) et croise deux
# facteurs :
#   - champ continu, ou sous-blocs disjoints de 5 x 5 parcelles (une parcelle ne
#     voit jamais une parcelle d'un autre sous-bloc ; le nombre de voisines n_i
#     varie alors beaucoup plus) ;
#   - noyau puissance w = delta^-lambda (lambda vrai 1,5), ou exponentiel
#     w = exp(-lambda (delta - 1)) (lambda vrai 0,8, meme poids a delta = 3).
# Le rayon est fixe a sa vraie valeur r = 3 et d vaut 0,5. Pour chaque cas et
# chaque graine, la vraisemblance REML du modele B (DGE + IGE independants,
# bloc, spl2d) est profilee sur une grille (lambda, d). Le profil en d (maximum
# sur lambda a d fixe) donne d estime, son intervalle a 95 % (chute de 1,92) et
# une erreur-type par la courbure.
#
#   Rscript test_identifiabilite_d.R [n_graines] [coeurs] [dossier_sortie]
# ==============================================================================

suppressPackageStartupMessages({
  library(remlax)
  library(ggplot2)
  library(patchwork)
})
a0 <- commandArgs(trailingOnly = TRUE)
n_graines <- if (length(a0) >= 1) as.integer(a0[1]) else 4L
coeurs <- if (length(a0) >= 2) as.integer(a0[2]) else 12L
sortie <- if (length(a0) >= 3) a0[3] else "sortie_identifiabilite_d"
dir.create(sortie, showWarnings = FALSE, recursive = TRUE)

# Dispositif, parente et tendance de la demo (sections 1 et 2, sans ajustement)
l <- readLines("demo_remlax_spatial_cpu.R")
pre <- l[1:(grep("^# 3. Ajustements", l) - 2)]
pre <- sub("^args <- commandArgs.*", "args <- character(0)", pre)
sortie_demo <- sortie
eval(parse(text = pre))
sortie <- sortie_demo

vraiD <- list(r = 3L, d = 0.5, lambda = c(puissance = 1.5, exponentiel = 0.8),
             dge = 1.0, ige = 0.4, r_dge_ige = -0.5, bloc = 0.3, residuelle = 0.8)
grilles <- list(puissance = seq(0.5, 3, 0.25), exponentiel = seq(0, 2, 0.2))
grille_d <- seq(0, 1.25, 0.25)
noyaux <- c(puissance = "power", exponentiel = "exponential")

# Sous-blocs disjoints : "continu" (aucun), "5x5" (carres de 5 x 5 parcelles),
# "heterogene" (bandes de 1, 2, 3, 4 puis 5 rangs sur 5 colonnes, deux fois).
bandes <- rep(seq_len(10), times = rep(1:5, 2))           # 30 rangs
sous_bloc <- function(d2, type) switch(type,
  continu = factor(rep(1L, nrow(d2))),
  "5x5" = factor(paste((d2$col - 1L) %/% 5L, (d2$rang - 1L) %/% 5L)),
  heterogene = factor(paste((d2$col - 1L) %/% 5L, bandes[d2$rang])))

# Jeu simule pour (graine, sous_blocs, noyau) : parcelles vides, valeurs
# genetiques, blocs, residus. Deterministe, donc refait dans chaque processus.
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

ajuster <- function(d2, Z) {
  random <- stats::as.formula("~ vm(gid, K = K) + mm(Z, K = K, name = 'ige') + bloc + spl2d(col, rang, nseg = c(7, 10))",
                              env = environment())
  rx_reml(y ~ rep, random = random, data = d2, backend = "cpu", verbose = FALSE,
          hessian = FALSE, blups = FALSE, pev = FALSE)
}

# Taches : cas x graine x lambda x d. Les jeux sont mis en cache par processus.
cas <- expand.grid(sous_blocs = c("continu", "5x5", "heterogene"), noyau = names(noyaux),
                   stringsAsFactors = FALSE)
taches <- do.call(rbind, lapply(seq_len(nrow(cas)), function(i) do.call(rbind, lapply(seq_len(n_graines), function(g)
  expand.grid(sous_blocs = cas$sous_blocs[i], noyau = cas$noyau[i], graine = 100L + g,
              lambda = grilles[[cas$noyau[i]]], d = grille_d, stringsAsFactors = FALSE)))))
cache <- new.env()
une_tache <- function(i) {
  t <- taches[i, ]
  cle <- paste(t$graine, t$sous_blocs, t$noyau)
  if (is.null(cache[[cle]])) cache[[cle]] <- simuler(t$graine, t$sous_blocs, t$noyau)
  d2 <- cache[[cle]]
  f <- tryCatch(ajuster(d2, voisinage(d2, vraiD$r, t$lambda, t$d, t$noyau)), error = function(e) conditionMessage(e))
  if (is.character(f)) return(cbind(t, logLik = NA_real_, conv = NA, erreur = substr(f, 1, 150)))
  cbind(t, logLik = f$logLik, conv = isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement)),
        erreur = "")
}
# Essai a blanc : une tache par cas, puis arret.
if (n_graines == 0L) {
  taches <- do.call(rbind, lapply(seq_len(nrow(cas)), function(i) data.frame(
    sous_blocs = cas$sous_blocs[i], noyau = cas$noyau[i], graine = 101L,
    lambda = vraiD$lambda[[cas$noyau[i]]], d = vraiD$d, stringsAsFactors = FALSE)))
  print(do.call(rbind, lapply(seq_len(nrow(taches)), une_tache)))
  for (ty in cas$sous_blocs[!duplicated(cas$sous_blocs)]) {
    d2 <- simuler(101L, ty, "puissance")
    nb <- rx_neighbourhood(cbind(d2$rang, d2$col), group = rep("A", nrow(d2)), block = d2$sb,
                           id = d2$parcelle, rank = 3L, window = "euclidean", output = "unit")
    n_i <- nb$n_neighbours[[1]]
    cat(sprintf("%-11s n_i de %d a %d, moyenne %.1f, sd(log n_i) = %.3f, %d sans voisine\n", ty,
                min(n_i), max(n_i), mean(n_i), sd(log(n_i[n_i > 0])), sum(n_i == 0)))
  }
  quit(save = "no")
}
Sys.setenv(XLA_FLAGS = "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2")
options(remlax.threads = 2L)
cat(sprintf("%d ajustements sur %d processus\n", nrow(taches), coeurs))
t0 <- Sys.time()
# Les taches d'un meme jeu vont au meme processus (prepartition par blocs contigus).
morceaux <- split(seq_len(nrow(taches)), cut(seq_len(nrow(taches)), coeurs, labels = FALSE))
res <- do.call(rbind, parallel::mclapply(morceaux, function(ii) do.call(rbind, lapply(ii, une_tache)),
                                          mc.cores = coeurs, mc.preschedule = FALSE))
cat(sprintf("fait en %.1f min ; %d echec(s), %d non convergé(s)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(is.na(res$logLik)),
            sum(!res$conv, na.rm = TRUE)))
if (any(is.na(res$logLik))) print(unique(res$erreur[is.na(res$logLik)]))
res$cas <- paste(c(continu = "champ continu", "5x5" = "sous-blocs 5 x 5",
                   heterogene = "sous-blocs heterogenes")[res$sous_blocs], res$noyau, sep = ", ")
write.csv(res, file.path(sortie, "profils.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Profil en d : maximum sur lambda a d fixe ; estimation, intervalle, courbure
# ------------------------------------------------------------------------------
parabole <- function(x, y) {           # sommet et courbure par les 3 points autour du max
  j <- which.max(y); j <- min(max(j, 2L), length(y) - 1L)
  cf <- stats::coef(stats::lm(y[(j - 1):(j + 1)] ~ poly(x[(j - 1):(j + 1)], 2, raw = TRUE)))
  list(sommet = if (cf[3] < 0) -cf[2] / (2 * cf[3]) else x[which.max(y)], courbure = 2 * cf[3])
}
intervalle <- function(x, y, chute = stats::qchisq(0.95, 1) / 2) {
  y <- y - max(y); fine <- seq(min(x), max(x), length.out = 501)
  yf <- stats::spline(x, y, xout = fine)$y
  dedans <- fine[yf >= -chute]
  c(bas = min(dedans), haut = max(dedans))
}
bilan <- do.call(rbind, lapply(split(res, list(res$cas, res$graine), drop = TRUE), function(q) {
  pd <- stats::aggregate(logLik ~ d, q, max); pl <- stats::aggregate(logLik ~ lambda, q, max)
  sd_ <- parabole(pd$d, pd$logLik); sl <- parabole(pl$lambda, pl$logLik)
  ic <- intervalle(pd$d, pd$logLik)
  data.frame(cas = q$cas[1], noyau = q$noyau[1], graine = q$graine[1],
             d_chapeau = min(max(sd_$sommet, 0), 1.25), lambda_chapeau = sl$sommet,
             ic_bas = ic[["bas"]], ic_haut = ic[["haut"]], largeur_ic = ic[["haut"]] - ic[["bas"]],
             et_d = if (sd_$courbure < 0) 1 / sqrt(-sd_$courbure) else NA_real_,
             chute_au_vrai = pd$logLik[pd$d == vraiD$d] - max(pd$logLik))
}))
bilan <- bilan[order(bilan$cas, bilan$graine), ]
write.csv(bilan, file.path(sortie, "bilan_d.csv"), row.names = FALSE)
cat("\nEstimation de d par cas et par graine (vrai d = 0,5)\n")
print(bilan, digits = 3, row.names = FALSE)
cat("\nMoyennes par cas\n")
print(stats::aggregate(cbind(d_chapeau, largeur_ic, et_d, chute_au_vrai) ~ cas, bilan, mean), digits = 3)

# ------------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 10) +
          theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10)))
q1 <- res[res$graine == min(res$graine), ]
q1 <- do.call(rbind, lapply(split(q1, q1$cas), function(q) { q$dll <- pmax(q$logLik - max(q$logLik), -12); q }))
vv <- unique(data.frame(cas = q1$cas, lambda = vraiD$lambda[q1$noyau], d = vraiD$d))
p_h <- ggplot(q1, aes(lambda, d)) + geom_raster(aes(fill = dll), interpolate = TRUE) +
  geom_contour(aes(z = dll), breaks = -stats::qchisq(0.95, 2) / 2, colour = "black", linewidth = 0.5) +
  geom_point(data = vv, shape = 4, size = 3) +
  facet_wrap(~ cas, nrow = 1, scales = "free_x") +
  scale_fill_viridis_c(name = "delta\nlogLik") + coord_cartesian(expand = FALSE) +
  labs(title = sprintf("a  Profil (lambda, d), graine %d ; croix = vrai, contour = region a 95 %%", min(res$graine)),
       x = "lambda", y = "d")
p_d <- ggplot(bilan, aes(factor(graine), d_chapeau, colour = noyau)) +
  geom_hline(yintercept = vraiD$d, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = ic_bas, ymax = ic_haut), size = 0.3) +
  facet_wrap(~ cas, nrow = 1) +
  scale_colour_manual(values = c(puissance = "#2166ac", exponentiel = "#d6604d"), guide = "none") +
  labs(title = "b  d estime et intervalle a 95 % du profil, par graine (tirets : vrai d = 0,5)",
       x = "graine", y = "d") + theme(panel.grid.major.y = element_line(colour = "grey90"))
fig <- p_h / p_d + plot_annotation(title = "Identifiabilite de la dilution : sous-blocs disjoints et forme du noyau")
ggsave(file.path(sortie, "fig_identifiabilite_d.png"), fig, width = 15, height = 8.5, dpi = 150, bg = "white")
cat("\nSorties dans", normalizePath(sortie), "\n")
