# ==============================================================================
# Noyau de voisinage libre : un anneau par distance, IGE correle au DGE
#
# Meme dispositif et memes simulations que test_identifiabilite_d.R (r = 3,
# d = 0,5, noyau puissance lambda = 1,5 ou exponentiel lambda = 0,8, champ
# continu ou sous-blocs heterogenes, r(DGE, IGE) = -0,5). Le modele ne connait
# plus la forme du noyau. Les voisines sont rangees en 6 anneaux de distance
# (1, rac(2), 2, rac(5), rac(8), 3), avec une incidence Z_k par anneau, divisee
# par n_i^d.
#
# Un seul terme couvre le DGE et les 6 anneaux : une structure fa sur 7
# "caracteres" (DGE, anneau 1, ..., anneau 6), de parente K. Avec un facteur f
# par genotype :
#   DGE       = lambda_D f + e_D      (variance specifique psi_D libre)
#   anneau k  = l_k f                 (variance specifique fixee a zero)
# Les l_k donnent le noyau a une echelle pres, et lambda_D porte la covariance
# DGE-IGE. Trois ajustements par point :
#   libre  fa1, toutes les variances specifiques libres (controle)
#   fa1    fa1, variances specifiques des anneaux fixees a zero (le modele)
#   fa2    fa2, idem : la forme du noyau peut varier entre genotypes (controle,
#          la verite est de rang 1)
# La dilution d est profilee sur une grille.
#
#   Rscript test_noyau_anneaux.R tache <i> <dossier_resultats>   (une tache)
#   Rscript test_noyau_anneaux.R bayes <i> <dossier_resultats>   (loi a posteriori,
#     cas x graine i parmi 16, a d = 0,5, depuis le fa1 contraint)
#   Rscript test_noyau_anneaux.R liste                            (nombre de taches)
#   Rscript test_noyau_anneaux.R bilan <dossier_resultats>        (bilan et figure)
# ==============================================================================

a0 <- commandArgs(trailingOnly = TRUE)
mode <- if (length(a0)) a0[1] else "liste"
`%||%` <- function(a, b) if (is.null(a)) b else a

vraiA <- list(r = 3L, d = 0.5, lambda = c(puissance = 1.5, exponentiel = 0.8),
              dge = 1.0, ige = 0.4, r_dge_ige = -0.5, bloc = 0.3, residuelle = 0.8)
noyaux <- c(puissance = "power", exponentiel = "exponential")
forme <- function(delta, noyau, lambda) if (noyau == "puissance") delta^(-lambda) else exp(-lambda * (delta - 1))
d2_anneaux <- c(1, 2, 4, 5, 8, 9)                 # delta^2 des 6 anneaux a r = 3
grille_d <- seq(0, 1.25, 0.25)
taches <- expand.grid(d = grille_d, graine = 101:104, noyau = names(noyaux),
                      sous_blocs = c("continu", "heterogene"), stringsAsFactors = FALSE)
taches <- taches[, c("sous_blocs", "noyau", "graine", "d")]

taches_b <- unique(taches[, c("sous_blocs", "noyau", "graine")]); taches_b$d <- vraiA$d
if (mode == "liste") { cat(nrow(taches), "taches REML,", nrow(taches_b), "taches bayesiennes\n"); quit(save = "no") }

# ------------------------------------------------------------------------------
# Mode tache
# ------------------------------------------------------------------------------
if (mode %in% c("tache", "bayes", "variances")) {
  i <- as.integer(a0[2]); dossier <- a0[3]
  if (mode == "variances") { i <- 1L; dossier <- a0[2] }
  dir.create(dossier, showWarnings = FALSE, recursive = TRUE)
  suppressPackageStartupMessages(library(remlax))
  # Dispositif, parente et tendance de la demo (sans ses graphiques)
  l <- readLines("demo_remlax_spatial_cpu.R")
  pre <- l[1:(grep("^# 3. Ajustements", l) - 2)]
  pre <- sub("^args <- commandArgs.*", "args <- character(0)", pre)
  pre <- pre[!grepl("library\\((ggplot2|patchwork)\\)", pre)]
  sortie <- tempdir()
  eval(parse(text = pre))

  bandes <- rep(seq_len(10), times = rep(1:5, 2))
  sous_bloc <- function(d2, type) switch(type,
    continu = factor(rep(1L, nrow(d2))),
    heterogene = factor(paste((d2$col - 1L) %/% 5L, bandes[d2$rang])))

  anneaux <- function(d2, dd) {
    n <- nrow(d2)
    q2 <- outer(d2$rang, d2$rang, "-")^2 + outer(d2$col, d2$col, "-")^2
    vois <- outer(as.integer(d2$sb), as.integer(d2$sb), "==") & q2 > 0 & q2 <= vraiA$r^2
    n_i <- rowSums(vois); div <- ifelse(n_i > 0, n_i^dd, 1)
    g <- as.integer(d2$gid)
    lapply(stats::setNames(d2_anneaux, sprintf("a%d", seq_along(d2_anneaux))), function(a) {
      ij <- which(vois & q2 == a, arr.ind = TRUE)
      methods::as(Matrix::sparseMatrix(i = ij[, 1], j = g[ij[, 2]], x = 1 / div[ij[, 1]],
                                       dims = c(n, n_gen), dimnames = list(NULL, noms_gen)), "CsparseMatrix")
    })
  }

  simuler <- function(graine, sous_blocs, noyau) {
    set.seed(graine)
    d2 <- d[sort(sample(nrow(d), round(0.9 * nrow(d)))), c("col", "rang", "rep", "bloc", "parcelle", "gid")]
    rownames(d2) <- NULL
    d2$bloc <- droplevels(d2$bloc)
    d2$sb <- sous_bloc(d2, sous_blocs)
    nb <- rx_neighbourhood(coord = cbind(d2$rang, d2$col), group = rep("A", nrow(d2)), block = d2$sb,
                           level = d2$gid, id = d2$parcelle, rank = vraiA$r,
                           reach = vraiA$lambda[[noyau]], dilution = vraiA$d, window = "euclidean",
                           kernel = noyaux[[noyau]], normalise = FALSE)
    Z <- nb$level[["A<-A"]]
    P <- Matrix::sparseMatrix(i = match(colnames(Z), noms_gen), j = seq_len(ncol(Z)), x = 1,
                              dims = c(n_gen, ncol(Z)))
    Zv <- Z %*% Matrix::t(P)
    # Controle : la somme des anneaux ponderes par la forme vraie redonne Zv
    Za <- anneaux(d2, vraiA$d)
    Zs <- Reduce(`+`, Map(function(Zk, a) Zk * forme(sqrt(a), noyau, vraiA$lambda[[noyau]]), Za, d2_anneaux))
    stopifnot(max(abs(as.matrix(Zs - Zv))) < 1e-10)
    Zd <- Matrix::sparseMatrix(i = seq_len(nrow(d2)), j = as.integer(d2$gid), x = 1,
                               dims = c(nrow(d2), n_gen), dimnames = list(NULL, noms_gen))
    k <- rx_exposure(Zd, Zv, K = K)$k
    s2i <- vraiA$ige / k
    Sig <- matrix(c(vraiA$dge, rep(vraiA$r_dge_ige * sqrt(vraiA$dge * s2i), 2), s2i), 2)
    U <- t(chol(K)) %*% matrix(rnorm(2 * n_gen), n_gen, 2) %*% chol(Sig)
    b <- setNames(rnorm(nlevels(d2$bloc), 0, sqrt(vraiA$bloc)), levels(d2$bloc))
    d2$y <- 6 + 0.2 * (d2$rep == "R2") - 0.1 * (d2$rep == "R3") + tendance_vraie(d2$col, d2$rang) +
      U[as.integer(d2$gid), 1] + as.numeric(Zv %*% U[, 2]) + b[as.character(d2$bloc)] +
      rnorm(nrow(d2), 0, sqrt(vraiA$residuelle))
    d2
  }

  # theta : le terme fa vient en premier (loadings trapezoidaux par ligne, puis
  # log psi des 7 caracteres), puis bloc, les 5 variances de spl2d, la residuelle.
  # Le controle ci-dessous verifie cette disposition sur la Sigma rendue.
  n_load <- function(t, r) sum(pmin(seq_len(t), r))
  ajuster <- function(d2, Za, rang, fixer_anneaux, theta0 = NULL, bayes = NULL) {
    for (k in names(Za)) assign(paste0("Z_", k), Za[[k]])
    # incidence DGE sur les 150 genotypes : un genotype sans parcelle semee
    # garde sa colonne, comme dans les incidences des anneaux
    Z_dge <- Matrix::sparseMatrix(i = seq_len(nrow(d2)), j = as.integer(d2$gid), x = 1,
                                  dims = c(nrow(d2), n_gen), dimnames = list(NULL, noms_gen))
    ige <- paste0("str(~ mm(Z_dge, name = 'dge') + ", paste(sprintf("mm(Z_%s, name = '%s')", names(Za), names(Za)), collapse = " + "),
                  sprintf(", struct = 'fa', rank = %d, K = K, name = 'ige')", rang))
    random <- stats::as.formula(paste("~", ige, "+ bloc + spl2d(col, rang, nseg = c(7, 10))"),
                                env = environment())
    nl <- n_load(7L, rang)
    fixe <- if (fixer_anneaux) nl + 2:7 else NULL
    t0 <- Sys.time()
    f <- rx_reml(y ~ rep, random = random, data = d2, backend = "cpu", verbose = FALSE,
                 hessian = FALSE, blups = FALSE, pev = FALSE, theta_init = theta0, fixed_theta = fixe,
                 bayes = bayes)
    th <- f$theta
    stopifnot(length(th) == nl + 7L + 7L)
    L <- matrix(0, 7, rang); kk <- 0
    for (a in 1:7) for (b in seq_len(min(a, rang))) { kk <- kk + 1; L[a, b] <- th[kk] }
    psi2 <- exp(2 * th[nl + 1:7])
    S <- as.matrix(f$sigmas$ige)
    stopifnot(max(abs(S - (L %*% t(L) + diag(psi2)))) < 1e-6 * max(abs(S)))
    list(f = f, L = L, psi2 = psi2, theta = th, nl = nl,
         secondes = as.numeric(difftime(Sys.time(), t0, units = "secs")),
         conv = isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement)))
  }

  # Variances sur l'echelle phenotypique, depuis les tirages bayesiens et le REML
  # a d = 0,5 : DGE = lambda_D^2 + psi_D ; IGE x k = moyenne de diag(M K M'),
  # M = somme_k l_k Z_k ; bloc et residuelle = exp(2 theta).
  if (mode == "variances") {
    fs <- list.files(dossier, "^bayes_[0-9]+\\.rds$", full.names = TRUE)
    out <- list()
    for (f in fs) {
      ib <- as.integer(sub("^bayes_([0-9]+)\\.rds$", "\\1", basename(f)))
      t <- taches_b[ib, ]
      d2 <- simuler(t$graine, t$sous_blocs, t$noyau)
      Za <- anneaux(d2, t$d)
      ige_k <- function(l) { M <- Reduce(`+`, Map(`*`, Za, l)); mean(Matrix::rowSums((M %*% K) * M)) }
      th <- readRDS(f)$theta
      sel <- unique(round(seq(1, nrow(th), length.out = 300)))
      nl <- 7L; p <- ncol(th)
      dr <- data.frame(dge = th[sel, 1]^2 + exp(2 * th[sel, nl + 1]),
                       ige_k = vapply(sel, function(s) ige_k(th[s, 2:7]), 0),
                       bloc = exp(2 * th[sel, nl + 8]), residuelle = exp(2 * th[sel, p]))
      q <- function(x) stats::quantile(x, c(0.05, 0.5, 0.95))
      # REML au meme point (tache fa1 a d = 0,5)
      it <- which(taches$sous_blocs == t$sous_blocs & taches$noyau == t$noyau &
                  taches$graine == t$graine & taches$d == t$d)
      fr <- file.path(dossier, sprintf("tache_%03d.csv", it))
      reml <- if (file.exists(fr)) { r <- utils::read.csv(fr)
        c(dge = r$lambda_D^2 + r$psi_D, ige_k = ige_k(unlist(r[sprintf("l%d", 1:6)]))) } else c(dge = NA, ige_k = NA)
      for (v in names(dr)) out[[length(out) + 1]] <- data.frame(
        t[, 1:3], composante = v, q05 = q(dr[[v]])[1], mediane = q(dr[[v]])[2], q95 = q(dr[[v]])[3],
        reml = if (v %in% names(reml)) reml[[v]] else NA_real_,
        vrai = c(dge = vraiA$dge, ige_k = vraiA$ige, bloc = vraiA$bloc, residuelle = vraiA$residuelle)[[v]])
    }
    res_v <- do.call(rbind, out); rownames(res_v) <- NULL
    res_v$couvre <- res_v$q05 <= res_v$vrai & res_v$vrai <= res_v$q95
    print(res_v, digits = 3)
    cat("\nCouverture des intervalles a 90 % par composante :\n")
    print(tapply(res_v$couvre, res_v$composante, function(x) sprintf("%d/%d", sum(x), length(x))))
    utils::write.csv(res_v, file.path(dossier, "variances_bayes.csv"), row.names = FALSE)
    quit(save = "no")
  }

  if (mode == "bayes") {
    stopifnot(rx_bayes_available())
    t <- taches_b[i, ]
    cat(sprintf("bayes %d : %s, %s, graine %d, d = %g\n", i, t$sous_blocs, t$noyau, t$graine, t$d))
    d2 <- simuler(t$graine, t$sous_blocs, t$noyau)
    Za <- anneaux(d2, t$d)
    libre <- ajuster(d2, Za, 1L, FALSE)
    th1 <- libre$theta; th1[libre$nl + 2:7] <- -12
    fa1 <- ajuster(d2, Za, 1L, TRUE, th1)
    bb <- ajuster(d2, Za, 1L, TRUE, fa1$theta,
                  bayes = rx_bayes(start = "map", chains = 8L, superchains = 4L, warmup = 400L,
                                   draws = 400L, seed = i, effects = 0L))
    th <- bb$f$bayes$theta                       # tirages x p
    L <- th[, 1:7]                               # DGE puis anneaux 1 a 6 (rang 1)
    sg <- sign(L[, 2]); sg[sg == 0] <- 1; L <- L * sg
    psiD <- exp(2 * th[, bb$nl + 1])
    w <- L[, 2:7] / L[, 2]
    q <- function(x, p) unname(stats::quantile(x, p))
    res_b <- data.frame(t, anneau = 1:6, delta = sqrt(d2_anneaux),
                        w_med = apply(w, 2, q, 0.5), w_q05 = apply(w, 2, q, 0.05), w_q95 = apply(w, 2, q, 0.95),
                        l_med = apply(L[, 2:7], 2, q, 0.5), l_q05 = apply(L[, 2:7], 2, q, 0.05),
                        l_q95 = apply(L[, 2:7], 2, q, 0.95), p_positif = colMeans(L[, 2:7] > 0),
                        w_vrai = forme(sqrt(d2_anneaux), t$noyau, vraiA$lambda[[t$noyau]]),
                        r_med = q(L[, 1] / sqrt(L[, 1]^2 + psiD), 0.5),
                        r_q05 = q(L[, 1] / sqrt(L[, 1]^2 + psiD), 0.05),
                        r_q95 = q(L[, 1] / sqrt(L[, 1]^2 + psiD), 0.95),
                        divergences = bb$f$bayes$diagnostics$divergences %||% NA,
                        rhat_max = bb$f$bayes$diagnostics$rhat_emboite_max %||% NA,
                        ess_min = bb$f$bayes$diagnostics$ess_min %||% NA,
                        l_reml = fa1$L[2:7, 1] * sign(fa1$L[2, 1]), secondes = bb$secondes)
    print(res_b, digits = 3)
    str(bb$f$bayes$diagnostics)
    utils::write.csv(res_b, file.path(dossier, sprintf("bayes_%02d.csv", i)), row.names = FALSE)
    saveRDS(list(theta = th, diagnostics = bb$f$bayes$diagnostics, summary = bb$f$bayes$summary),
            file.path(dossier, sprintf("bayes_%02d.rds", i)))
    quit(save = "no")
  }

  t <- taches[i, ]
  cat(sprintf("tache %d : %s, %s, graine %d, d = %g\n", i, t$sous_blocs, t$noyau, t$graine, t$d))
  d2 <- simuler(t$graine, t$sous_blocs, t$noyau)
  Za <- anneaux(d2, t$d)
  libre <- ajuster(d2, Za, 1L, FALSE)
  # fa1 contraint : depart a l'optimum libre, log psi des anneaux au plancher
  th1 <- libre$theta; th1[libre$nl + 2:7] <- -12
  fa1 <- ajuster(d2, Za, 1L, TRUE, th1)
  # fa2 contraint : premier facteur = celui de fa1, second facteur petit
  L2 <- cbind(fa1$L[, 1], c(0, rep(0.05, 6)))
  th2 <- c(unlist(lapply(1:7, function(a) L2[a, seq_len(min(a, 2L))])),
           fa1$theta[fa1$nl + 1:7], fa1$theta[(fa1$nl + 8):length(fa1$theta)])
  fa2 <- ajuster(d2, Za, 2L, TRUE, th2)
  # signe : chargement de l'anneau 1 positif
  l <- fa1$L[, 1]; if (l[2] < 0) l <- -l
  ligne <- data.frame(t, logLik_libre = libre$f$logLik, logLik_fa1 = fa1$f$logLik, logLik_fa2 = fa2$f$logLik,
                      conv_libre = libre$conv, conv_fa1 = fa1$conv, conv_fa2 = fa2$conv,
                      lambda_D = l[1], psi_D = fa1$psi2[1], r_dge_ige = l[1] / sqrt(l[1]^2 + fa1$psi2[1]),
                      t(stats::setNames(l[2:7], sprintf("l%d", 1:6))),
                      psi_anneaux_libre_max = max(libre$psi2[2:7]),
                      secondes = libre$secondes + fa1$secondes + fa2$secondes)
  print(ligne)
  utils::write.csv(ligne, file.path(dossier, sprintf("tache_%03d.csv", i)), row.names = FALSE)
  quit(save = "no")
}

# ------------------------------------------------------------------------------
# Mode bilan
# ------------------------------------------------------------------------------
if (mode == "bilan") {
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  dossier <- a0[2]
  fs <- list.files(dossier, "^tache_.*\\.csv$", full.names = TRUE)
  res <- do.call(rbind, lapply(fs, utils::read.csv))
  cat(sprintf("%d taches sur %d\n", nrow(res), nrow(taches)))
  res$cas <- paste(c(continu = "champ continu", heterogene = "sous-blocs heterogenes")[res$sous_blocs],
                   res$noyau, sep = ", ")
  bilan <- do.call(rbind, lapply(split(res, list(res$cas, res$graine), drop = TRUE), function(q) {
    q <- q[order(q$d), ]; j <- which.max(q$logLik_fa1); jj <- min(max(j, 2L), nrow(q) - 1L)
    cf <- stats::coef(stats::lm(q$logLik_fa1[(jj - 1):(jj + 1)] ~ poly(q$d[(jj - 1):(jj + 1)], 2, raw = TRUE)))
    d_hat <- if (cf[3] < 0) min(max(-cf[2] / (2 * cf[3]), 0), 1.25) else q$d[j]
    fine <- seq(0, 1.25, length.out = 501)
    yf <- stats::spline(q$d, q$logLik_fa1 - max(q$logLik_fa1), xout = fine)$y
    dedans <- fine[yf >= -stats::qchisq(0.95, 1) / 2]
    data.frame(cas = q$cas[1], graine = q$graine[1], d_chapeau = d_hat, ic_bas = min(dedans),
               ic_haut = max(dedans), chute_au_vrai = q$logLik_fa1[q$d == vraiA$d] - max(q$logLik_fa1),
               r_dge_ige = q$r_dge_ige[j], lrt_rang2 = 2 * (q$logLik_fa2[j] - q$logLik_fa1[j]),
               lrt_psi_anneaux = 2 * (q$logLik_libre[j] - q$logLik_fa1[j]),
               n_non_conv = sum(!q$conv_fa1))
  }))
  bilan <- bilan[order(bilan$cas, bilan$graine), ]
  cat("\nNoyau libre par anneaux (vrai d = 0,5, vrai r(DGE, IGE) = -0,5, vrai rang 1)\n")
  print(bilan, digits = 3, row.names = FALSE)
  cat("\nMoyennes par cas\n")
  print(stats::aggregate(cbind(d_chapeau, largeur_ic = ic_haut - ic_bas, r_dge_ige, lrt_rang2) ~ cas,
                         bilan, mean), digits = 3)
  utils::write.csv(bilan, file.path(dossier, "bilan_anneaux.csv"), row.names = FALSE)

  fo <- do.call(rbind, lapply(split(res, list(res$cas, res$graine), drop = TRUE), function(q) {
    q <- q[which.max(q$logLik_fa1), ]; l <- unlist(q[sprintf("l%d", 1:6)])
    data.frame(cas = q$cas, noyau = q$noyau, graine = q$graine, delta = sqrt(d2_anneaux), w = l / l[1])
  }))
  moy <- stats::aggregate(w ~ cas + delta, fo, mean)
  vf <- do.call(rbind, lapply(names(noyaux), function(nk) data.frame(
    noyau = nk, delta = seq(1, 3, 0.02), w = forme(seq(1, 3, 0.02), nk, vraiA$lambda[[nk]]))))
  vf <- merge(vf, unique(fo[, c("cas", "noyau")]))
  theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10)))
  p_f <- ggplot(fo, aes(delta, w)) + geom_hline(yintercept = 0, colour = "grey70") +
    geom_line(data = vf, colour = "black", linewidth = 0.7) +
    geom_line(aes(group = graine), colour = "#1b7837", alpha = 0.35) +
    geom_point(colour = "#1b7837", size = 1.2, alpha = 0.6) +
    geom_line(data = moy, colour = "#1b7837", linewidth = 1, linetype = 2) +
    facet_wrap(~ cas, nrow = 1) +
    labs(title = "a  Noyau libre au meilleur d : chargements relatifs a l'anneau 1",
         subtitle = "trait fin : une graine ; tirets : moyenne des graines ; noir : forme vraie",
         x = "distance (pas de grille)", y = "poids relatif") +
    theme(panel.grid.major.y = element_line(colour = "grey90"))
  p_d <- ggplot(bilan, aes(factor(graine), d_chapeau)) +
    geom_hline(yintercept = vraiA$d, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(ymin = ic_bas, ymax = ic_haut), colour = "#1b7837", size = 0.3) +
    facet_wrap(~ cas, nrow = 1) +
    labs(title = "b  d estime avec le noyau libre, intervalle a 95 % du profil (tirets : vrai d = 0,5)",
         x = "graine", y = "d") + theme(panel.grid.major.y = element_line(colour = "grey90"))
  fig <- p_f / p_d + plot_annotation(title = "Noyau de voisinage libre par anneaux, IGE correle au DGE (fa, anneaux sans variance specifique)")
  ggsave(file.path(dossier, "fig_noyau_anneaux.png"), fig, width = 15, height = 8.5, dpi = 150, bg = "white")
  cat("\nFigure :", normalizePath(file.path(dossier, "fig_noyau_anneaux.png")), "\n")
}
