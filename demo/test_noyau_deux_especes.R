# ==============================================================================
# Estimation du noyau de voisinage (lambda, d) et effet d'un noyau faux sur les
# IGE, dans le cadre du chapitre 2 de la these (Salas et al., Genetics)
#
# Simulation a deux especes A et B, comme dans l'article : 100 genotypes par
# espece, DGE, IGE intra et IGE inter correles (r(DGE, IGE intra) = r(DGE, IGE
# inter) = -0,6, r(IGE intra, IGE inter) = 0,6), genotypes non apparentes
# (M = I), H2_DGE = 0,3, pas d'IEE. Le dispositif reprend la geometrie du
# dispositif reel : bacs de 20 colonnes x 16 rangs, especes alternees par
# colonne, pas de grille 1 (5 cm), 10 % de plantes absentes. Les bacs et les
# plantes absentes font varier le nombre de voisins, sans quoi d ne serait pas
# identifiable (un champ sans bord donne a toutes les plantes le meme n_i).
#
# Noyau vrai : fenetre euclidienne de rayon r = 3, poids delta^-lambda avec
# lambda = 1,5, dilution n_i^-d avec d = 0,5, n_i compte par espece emettrice.
# Les variances IGE sont fixees sur l'echelle phenotypique : k x V(IGE intra) =
# 0,2 et k x V(IGE inter) = 0,3 (V(DGE) = 1), soit l'ordre de grandeur de
# l'article (rapport 0,05 a 0,1 par voisin, 2 voisins intra et 6 inter).
#
# Par jeu de donnees (n_rep in {5, 10, 20}, 30 graines), le modele de l'article
# (us 3 x 3 par espece, residuelle par espece, bac fixe) est ajuste :
#   - "r3"   : grille lambda x d a r = 3 (le vrai rayon) ;
#   - "R6"   : la meme grille dans une fenetre large R = 6 ;
#   - "r"    : r in {1, 2, 4, 5} aux vrais lambda et d ;
#   - "naif" : rang 1 de Chebyshev (8 voisins), poids egaux, sans dilution,
#              le voisinage des simulations de l'article.
# Chaque ajustement donne sa logLik (profil) et, contre la verite simulee, la
# precision des BLUP (DGE, IGE intra, IGE inter, TBV intra) et l'erreur
# log10(V estimee x k / V vraie x k) de chaque composante.
#
#   Rscript test_noyau_deux_especes.R liste
#   Rscript test_noyau_deux_especes.R tache <i> <dossier> [coeurs]
#   Rscript test_noyau_deux_especes.R bilan <dossier>
# ==============================================================================

a0 <- commandArgs(trailingOnly = TRUE)
mode <- if (length(a0)) a0[1] else "liste"
vrai <- list(r = 3, lambda = 1.5, d = 0.5, h2 = 0.3, dge = 1, ige_intra_k = 0.2, ige_inter_k = 0.3,
             r_dge_intra = -0.6, r_dge_inter = -0.6, r_intra_inter = 0.6, absentes = 0.10)
grille_l <- c(0.5, 1, 1.5, 2, 2.5, 3)
grille_d <- c(0, 0.25, 0.5, 0.75, 1)
r_faux <- c(1, 2, 4, 5)
R_LARGE <- 6
taches <- expand.grid(graine = 1:30, n_rep = c(5L, 10L, 20L))[, c("n_rep", "graine")]
if (mode == "liste") { cat(nrow(taches), "taches\n"); quit(save = "no") }

# Liste des ajustements d'une tache
plan_fits <- function() {
  g <- expand.grid(lambda = grille_l, d = grille_d)
  rbind(data.frame(grille = "r3", r = vrai$r, g, fenetre = "euclidean"),
        data.frame(grille = "R6", r = R_LARGE, g, fenetre = "euclidean"),
        data.frame(grille = "r", r = r_faux, lambda = vrai$lambda, d = vrai$d, fenetre = "euclidean"),
        data.frame(grille = "naif", r = 1, lambda = 0, d = 0, fenetre = "chebyshev"))
}

if (mode == "tache") {
  suppressPackageStartupMessages({ library(remlax); library(Matrix) })
  i <- as.integer(a0[2]); dossier <- a0[3]
  coeurs <- if (length(a0) >= 4) as.integer(a0[4]) else 1L
  dir.create(dossier, showWarnings = FALSE, recursive = TRUE)
  t <- taches[i, ]
  n_gen <- 100L
  gA <- sprintf("A%03d", seq_len(n_gen)); gB <- sprintf("B%03d", seq_len(n_gen))

  # --- dispositif : bacs 20 x 16, colonnes impaires A, paires B -----------------
  set.seed(1000L * t$n_rep + t$graine)
  n_bac <- ceiling(n_gen * t$n_rep / 160)
  d <- do.call(rbind, lapply(seq_len(n_bac), function(b)
    data.frame(bac = sprintf("bac%02d", b), col = rep(1:20, 16), rang = rep(1:16, each = 20))))
  d$espece <- ifelse(d$col %% 2 == 1, "A", "B")
  d$gen <- NA_character_
  for (sp in c("A", "B")) {
    w <- which(d$espece == sp)
    d$gen[w] <- sample(rep_len(if (sp == "A") gA else gB, length(w)))
  }
  d <- d[sort(sample(nrow(d), round((1 - vrai$absentes) * nrow(d)))), ]
  d <- d[order(d$espece != "A"), ]                   # lignes A puis B
  rownames(d) <- NULL
  d$id <- sprintf("P%05d", seq_len(nrow(d)))
  d$espece <- factor(d$espece); d$bac <- factor(d$bac)
  d$esp_bac <- factor(paste(d$espece, d$bac))
  n <- nrow(d); lA <- d$espece == "A"; lB <- !lA

  # Incidence n x 100 d'une paire recevant<-emettant, colonnes = genotypes emetteurs
  incid <- function(nb, paire, niv) {
    Z <- nb$level[[paire]]
    Zt <- methods::as(Z, "TsparseMatrix")
    Matrix::sparseMatrix(i = match(rownames(Z), d$id)[Zt@i + 1L], j = match(colnames(Z), niv)[Zt@j + 1L],
                         x = Zt@x, dims = c(n, length(niv)), dimnames = list(NULL, niv))
  }
  voisinage <- function(r, lambda, dd, fenetre = "euclidean") {
    nb <- rx_neighbourhood(coord = cbind(d$rang, d$col), group = as.character(d$espece), block = d$bac,
                           level = d$gen, id = d$id, rank = r, reach = unname(lambda),
                           dilution = unname(dd), window = fenetre, kernel = "power", output = "level")
    list(AA = incid(nb, "A<-A", gA), BA = incid(nb, "B<-A", gA),   # emis par A
         BB = incid(nb, "B<-B", gB), AB = incid(nb, "A<-B", gB))   # emis par B
  }
  kexp <- function(Z, lignes) mean(Matrix::rowSums(Z^2)[lignes])  # k = moyenne diag(ZZ'), M = I
  ZdA <- Matrix::sparseMatrix(i = which(lA), j = match(d$gen[lA], gA), x = 1, dims = c(n, n_gen),
                              dimnames = list(NULL, gA))
  ZdB <- Matrix::sparseMatrix(i = which(lB), j = match(d$gen[lB], gB), x = 1, dims = c(n, n_gen),
                              dimnames = list(NULL, gB))

  # --- verite -----------------------------------------------------------------
  Vv <- voisinage(vrai$r, vrai$lambda, vrai$d)
  k_vrai <- c(AA = kexp(Vv$AA, lA), BA = kexp(Vv$BA, lB), BB = kexp(Vv$BB, lB), AB = kexp(Vv$AB, lA))
  tirer <- function(k_intra, k_inter) {
    v <- c(vrai$dge, vrai$ige_intra_k / k_intra, vrai$ige_inter_k / k_inter)
    R <- matrix(c(1, vrai$r_dge_intra, vrai$r_dge_inter, vrai$r_dge_intra, 1, vrai$r_intra_inter,
                  vrai$r_dge_inter, vrai$r_intra_inter, 1), 3)
    S <- diag(sqrt(v)) %*% R %*% diag(sqrt(v))
    list(U = matrix(rnorm(3 * n_gen), n_gen, 3) %*% chol(S), S = S)
  }
  tA <- tirer(k_vrai[["AA"]], k_vrai[["BA"]]); tB <- tirer(k_vrai[["BB"]], k_vrai[["AB"]])
  UA <- tA$U; UB <- tB$U; dimnames(UA) <- dimnames(UB) <- NULL
  v_dee <- vrai$dge * (1 / vrai$h2 - 1)
  g <- as.numeric(ZdA %*% UA[, 1] + Vv$AA %*% UA[, 2] + Vv$BA %*% UA[, 3] +
                  ZdB %*% UB[, 1] + Vv$BB %*% UB[, 2] + Vv$AB %*% UB[, 3])
  d$y <- ifelse(lA, 10, 5) + g + rnorm(n, 0, sqrt(v_dee))
  # TBV intra vraie : DGE + (somme moyenne des poids intra) x IGE intra
  s_vrai <- c(A = mean(Matrix::rowSums(Vv$AA)[lA]), B = mean(Matrix::rowSums(Vv$BB)[lB]))
  vraies <- list(A = cbind(UA, UA[, 1] + s_vrai[["A"]] * UA[, 2]),
                 B = cbind(UB, UB[, 1] + s_vrai[["B"]] * UB[, 2]))
  vk <- list(A = c(vrai$dge, vrai$ige_intra_k, vrai$ige_inter_k),
             B = c(vrai$dge, vrai$ige_intra_k, vrai$ige_inter_k))

  # --- un ajustement ----------------------------------------------------------
  ajuster <- function(p) {
    t0 <- Sys.time()
    V <- if (p$grille == "r3" && p$lambda == vrai$lambda && p$d == vrai$d) Vv
         else voisinage(p$r, p$lambda, p$d, p$fenetre)
    ZaD <- ZdA; ZaI <- V$AA; ZaX <- V$BA; ZbD <- ZdB; ZbI <- V$BB; ZbX <- V$AB
    random <- ~ str(~ mm(ZaD) + mm(ZaI) + mm(ZaX), struct = "us", name = "gA") +
                str(~ mm(ZbD) + mm(ZbI) + mm(ZbX), struct = "us", name = "gB")
    environment(random) <- environment()
    f <- tryCatch(rx_reml(y ~ esp_bac, random = random, residual = ~ dsum(~ units | espece), data = d,
                          backend = "cpu", verbose = FALSE, hessian = FALSE, pev = FALSE),
                  error = function(e) conditionMessage(e))
    base <- data.frame(p, n = n, logLik = NA_real_, conv = NA, secondes = NA_real_, erreur = "")
    if (is.character(f)) { base$erreur <- substr(gsub("\n", " ", f), 1, 200); return(base) }
    base$logLik <- f$logLik
    base$conv <- isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement))
    v <- summary(f)$varcomp
    k_f <- list(A = c(1, kexp(V$AA, lA), kexp(V$BA, lB)), B = c(1, kexp(V$BB, lB), kexp(V$AB, lA)))
    s_f <- c(A = mean(Matrix::rowSums(V$AA)[lA]), B = mean(Matrix::rowSums(V$BB)[lB]))
    eff <- c("dge", "intra", "inter")
    for (sp in c("A", "B")) {
      nm <- paste0("g", sp)
      S <- matrix(NA_real_, 3, 3)
      for (a in 1:3) for (b in 1:a) {
        lab <- sprintf("%s[%d,%d]", nm, a, b)
        if (lab %in% rownames(v)) S[a, b] <- S[b, a] <- v[lab, "component"]
      }
      Bl <- f$blups[[nm]]
      Bl <- Bl[if (sp == "A") gA else gB, , drop = FALSE]
      Bl <- cbind(Bl, Bl[, 1] + s_f[[sp]] * Bl[, 2])
      for (a in 1:3) {
        base[[sprintf("vk_%s_%s", eff[a], sp)]] <- S[a, a] * k_f[[sp]][a]
        base[[sprintf("err_%s_%s", eff[a], sp)]] <- log10(S[a, a] * k_f[[sp]][a] / vk[[sp]][a])
      }
      base[[sprintf("rho_dge_intra_%s", sp)]] <- S[2, 1] / sqrt(S[1, 1] * S[2, 2])
      base[[sprintf("rho_dge_inter_%s", sp)]] <- S[3, 1] / sqrt(S[1, 1] * S[3, 3])
      base[[sprintf("rho_intra_inter_%s", sp)]] <- S[3, 2] / sqrt(S[2, 2] * S[3, 3])
      for (a in 1:4) base[[sprintf("prec_%s_%s", c(eff, "tbv")[a], sp)]] <- stats::cor(Bl[, a], vraies[[sp]][, a])
    }
    base$secondes <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    base
  }

  pl <- plan_fits()
  if (identical(Sys.getenv("ESSAI"), "1"))           # essai rapide : verite, naif, un rayon faux
    pl <- pl[(pl$grille == "r3" & pl$lambda == vrai$lambda & pl$d == vrai$d) | pl$grille == "naif" |
             (pl$grille == "r" & pl$r == 1), ]
  options(remlax.threads = 1L)
  cat(sprintf("tache %d : n_rep %d, graine %d, %d bacs, n = %d, %d ajustements sur %d coeurs\n",
              i, t$n_rep, t$graine, n_bac, n, nrow(pl), coeurs))
  cat(sprintf("k vrais : AA %.3f BA %.3f BB %.3f AB %.3f\n", k_vrai[1], k_vrai[2], k_vrai[3], k_vrai[4]))
  t_tache <- Sys.time()
  res <- parallel::mclapply(seq_len(nrow(pl)), function(j) ajuster(pl[j, ]), mc.cores = coeurs,
                            mc.preschedule = FALSE)
  ok <- vapply(res, is.data.frame, TRUE)
  if (!all(ok)) {
    cat("ajustements perdus :", which(!ok), "\n")
    res[!ok] <- lapply(which(!ok), function(j) data.frame(pl[j, ], n = n, logLik = NA_real_, conv = NA,
                                                         secondes = NA_real_, erreur = "processus perdu"))
  }
  nm_all <- unique(unlist(lapply(res, names)))
  res <- do.call(rbind, lapply(res, function(x) { x[setdiff(nm_all, names(x))] <- NA; x[nm_all] }))
  res <- cbind(n_rep = t$n_rep, graine = t$graine, res)
  utils::write.csv(res, file.path(dossier, sprintf("tache_%03d.csv", i)), row.names = FALSE)
  cat(sprintf("fini en %.1f min ; %d echec(s), %d non converge(s)\n",
              as.numeric(difftime(Sys.time(), t_tache, units = "mins")), sum(is.na(res$logLik)),
              sum(!res$conv, na.rm = TRUE)))
  print(res[res$grille != "r3" | (res$lambda == vrai$lambda), c("grille", "r", "lambda", "d", "logLik",
          "conv", "prec_intra_A", "prec_inter_A", "err_intra_A", "err_inter_A", "secondes")], digits = 3)
  quit(save = "no")
}

if (mode == "bilan") {
  dossier <- a0[2]
  source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))),
                   "test_noyau_deux_especes_bilan.R"))
}
