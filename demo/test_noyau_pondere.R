# ==============================================================================
# Dilution par la masse du noyau : essai sur un jeu de donnees
#
# Meme simulation que test_noyau_deux_especes.R (chapitre 2 de la these) : deux
# especes, 100 genotypes chacune, DGE, IGE intra et IGE inter correles, bacs de
# 20 colonnes x 16 rangs, especes alternees par colonne, 10 % d'absentes.
#
# La dilution de rx_neighbourhood divise chaque ligne par n_i^d, n_i etant le
# NOMBRE de voisins dans la fenetre de rayon r : d depend alors de r. Ici la
# taille du groupe est la MASSE du noyau, m_i = somme_j w_ij, et chaque ligne
# est divisee par m_i^d (d = 0 : somme ; d = 1 : moyenne ponderee). Avec un
# noyau exponentiel, la masse converge et la fenetre ne sert qu'au calcul.
#
# Deux verites sur le meme dispositif (n_rep = 10, graine 1) :
#   A : noyau exponentiel w = exp(-lambda (delta - 1)), lambda = 1 par pas de
#       grille, dilution par la masse d = 0,5, fenetre de calcul R = 8 ;
#   B : l'ancienne verite, puissance delta^-1,5 tronquee a r = 3, dilution par
#       le nombre d = 0,5.
# Ajustements pour chaque verite :
#   - "pond"    : grille lambda x d, exponentiel, masse, R = 8 ;
#   - "fenetre" : lambda = 1, d = 0,5, masse, R in {4, 6, 12} ;
#   - "nombre"  : exponentiel lambda = 1, dilution par le NOMBRE, R = 8, d en grille ;
#   - "puis"    : puissance r = 3, lambda = 1,5, d = 0,5, nombre (verite B) ;
#   - "chap3"   : puissance r = 3, lambda = 1, d = 0 ;
#   - "naif"    : rang 1 de Chebyshev, poids egaux.
#
# Une tache = une verite et une graine (30 graines x 2 verites = 60 taches) ;
# la graine fixe le dispositif et les valeurs genetiques.
#
#   Rscript test_noyau_pondere.R tache <i> <dossier> [coeurs]
# ==============================================================================

a0 <- commandArgs(trailingOnly = TRUE)
i_tache <- as.integer(a0[2]); dossier <- a0[3]
coeurs <- if (length(a0) >= 4) as.integer(a0[4]) else 1L
nv <- c("A", "B")[(i_tache - 1L) %/% 30L + 1L]
graine <- (i_tache - 1L) %% 30L + 1L
dir.create(dossier, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(remlax); library(Matrix) })

vrai <- list(h2 = 0.3, dge = 1, ige_intra_k = 0.2, ige_inter_k = 0.3,
             r_dge_intra = -0.6, r_dge_inter = -0.6, r_intra_inter = 0.6, absentes = 0.10)
n_rep <- 10L; n_gen <- 100L
gA <- sprintf("A%03d", seq_len(n_gen)); gB <- sprintf("B%03d", seq_len(n_gen))

# --- dispositif (identique a test_noyau_deux_especes.R) -----------------------
set.seed(1000L * n_rep + graine)
n_bac <- ceiling(n_gen * n_rep / 160)
d <- do.call(rbind, lapply(seq_len(n_bac), function(b)
  data.frame(bac = sprintf("bac%02d", b), col = rep(1:20, 16), rang = rep(1:16, each = 20))))
d$espece <- ifelse(d$col %% 2 == 1, "A", "B")
d$gen <- NA_character_
for (sp in c("A", "B")) {
  w <- which(d$espece == sp)
  d$gen[w] <- sample(rep_len(if (sp == "A") gA else gB, length(w)))
}
d <- d[sort(sample(nrow(d), round((1 - vrai$absentes) * nrow(d)))), ]
d <- d[order(d$espece != "A"), ]
rownames(d) <- NULL
d$id <- sprintf("P%05d", seq_len(nrow(d)))
d$espece <- factor(d$espece); d$bac <- factor(d$bac)
d$esp_bac <- factor(paste(d$espece, d$bac))
n <- nrow(d); lA <- d$espece == "A"; lB <- !lA

incid <- function(nb, paire, niv) {
  Z <- nb$level[[paire]]
  Zt <- methods::as(Z, "TsparseMatrix")
  Matrix::sparseMatrix(i = match(rownames(Z), d$id)[Zt@i + 1L], j = match(colnames(Z), niv)[Zt@j + 1L],
                       x = Zt@x, dims = c(n, length(niv)), dimnames = list(NULL, niv))
}
# Voisinage. dilution = "nombre" : n_i^-d de rx_neighbourhood ; "masse" : m_i^-d,
# m_i = somme des poids de la ligne (rx_neighbourhood appele sans dilution).
voisinage <- function(noyau, r, lambda, dd, dilution, fenetre = "euclidean") {
  nb <- rx_neighbourhood(coord = cbind(d$rang, d$col), group = as.character(d$espece), block = d$bac,
                         level = d$gen, id = d$id, rank = r, reach = lambda,
                         dilution = if (dilution == "nombre") dd else 0, window = fenetre,
                         kernel = noyau, output = "level")
  V <- list(AA = incid(nb, "A<-A", gA), BA = incid(nb, "B<-A", gA),
            BB = incid(nb, "B<-B", gB), AB = incid(nb, "A<-B", gB))
  if (dilution == "masse" && dd > 0) V <- lapply(V, function(Z) {
    m <- Matrix::rowSums(Z); m[m == 0] <- 1
    Zd <- methods::as(Matrix::Diagonal(x = m^(-dd)) %*% Z, "CsparseMatrix")
    dimnames(Zd) <- dimnames(Z); Zd
  })
  V
}
kexp <- function(Z, lignes) mean(Matrix::rowSums(Z^2)[lignes])
ZdA <- Matrix::sparseMatrix(i = which(lA), j = match(d$gen[lA], gA), x = 1, dims = c(n, n_gen),
                            dimnames = list(NULL, gA))
ZdB <- Matrix::sparseMatrix(i = which(lB), j = match(d$gen[lB], gB), x = 1, dims = c(n, n_gen),
                            dimnames = list(NULL, gB))

verites <- list(
  A = list(noyau = "exponential", r = 8, lambda = 1, d = 0.5, dilution = "masse"),
  B = list(noyau = "power", r = 3, lambda = 1.5, d = 0.5, dilution = "nombre"))

# Tirage des valeurs genetiques et du phenotype sous une verite, graine propre
# a chaque jeu.
simuler <- function(vt) {
  Vv <- voisinage(vt$noyau, vt$r, vt$lambda, vt$d, vt$dilution)
  kv <- c(AA = kexp(Vv$AA, lA), BA = kexp(Vv$BA, lB), BB = kexp(Vv$BB, lB), AB = kexp(Vv$AB, lA))
  set.seed(5000L + graine)
  R <- matrix(c(1, vrai$r_dge_intra, vrai$r_dge_inter, vrai$r_dge_intra, 1, vrai$r_intra_inter,
                vrai$r_dge_inter, vrai$r_intra_inter, 1), 3)
  tirer <- function(ki, kx) {
    v <- c(vrai$dge, vrai$ige_intra_k / ki, vrai$ige_inter_k / kx)
    matrix(rnorm(3 * n_gen), n_gen, 3) %*% chol(diag(sqrt(v)) %*% R %*% diag(sqrt(v)))
  }
  UA <- tirer(kv[["AA"]], kv[["BA"]]); UB <- tirer(kv[["BB"]], kv[["AB"]])
  g <- as.numeric(ZdA %*% UA[, 1] + Vv$AA %*% UA[, 2] + Vv$BA %*% UA[, 3] +
                  ZdB %*% UB[, 1] + Vv$BB %*% UB[, 2] + Vv$AB %*% UB[, 3])
  y <- ifelse(lA, 10, 5) + g + rnorm(n, 0, sqrt(vrai$dge * (1 / vrai$h2 - 1)))
  s <- c(A = mean(Matrix::rowSums(Vv$AA)[lA]), B = mean(Matrix::rowSums(Vv$BB)[lB]))
  list(y = y, k = kv, s = s,
       vraies = list(A = cbind(UA, UA[, 1] + s[["A"]] * UA[, 2]), B = cbind(UB, UB[, 1] + s[["B"]] * UB[, 2])))
}

plan <- function() {
  g <- expand.grid(lambda = c(0.5, 0.75, 1, 1.5, 2), d = c(0, 0.25, 0.5, 0.75, 1))
  rbind(data.frame(grille = "pond", noyau = "exponential", r = 8, g, dilution = "masse", fenetre = "euclidean"),
        data.frame(grille = "fenetre", noyau = "exponential", r = c(4, 6, 12), lambda = 1, d = 0.5,
                   dilution = "masse", fenetre = "euclidean"),
        data.frame(grille = "nombre", noyau = "exponential", r = 8, lambda = 1, d = c(0, 0.25, 0.5, 0.75, 1),
                   dilution = "nombre", fenetre = "euclidean"),
        data.frame(grille = "puis", noyau = "power", r = 3, lambda = 1.5, d = 0.5, dilution = "nombre",
                   fenetre = "euclidean"),
        data.frame(grille = "chap3", noyau = "power", r = 3, lambda = 1, d = 0, dilution = "nombre",
                   fenetre = "euclidean"),
        data.frame(grille = "naif", noyau = "none", r = 1, lambda = 0, d = 0, dilution = "nombre",
                   fenetre = "chebyshev"))
}

ajuster <- function(p, sim) {
  t0 <- Sys.time()
  dd <- d; dd$y <- sim$y
  V <- voisinage(p$noyau, p$r, p$lambda, p$d, p$dilution, p$fenetre)
  ZaD <- ZdA; ZaI <- V$AA; ZaX <- V$BA; ZbD <- ZdB; ZbI <- V$BB; ZbX <- V$AB
  random <- ~ str(~ mm(ZaD) + mm(ZaI) + mm(ZaX), struct = "us", name = "gA") +
              str(~ mm(ZbD) + mm(ZbI) + mm(ZbX), struct = "us", name = "gB")
  environment(random) <- environment()
  f <- tryCatch(rx_reml(y ~ esp_bac, random = random, residual = ~ dsum(~ units | espece), data = dd,
                        backend = "cpu", verbose = FALSE, hessian = FALSE, pev = FALSE),
                error = function(e) conditionMessage(e))
  base <- data.frame(p, logLik = NA_real_, conv = NA, secondes = NA_real_, erreur = "")
  if (is.character(f)) { base$erreur <- substr(gsub("\n", " ", f), 1, 200); return(base) }
  base$logLik <- f$logLik
  base$conv <- isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement))
  v <- summary(f)$varcomp
  k_f <- list(A = c(1, kexp(V$AA, lA), kexp(V$BA, lB)), B = c(1, kexp(V$BB, lB), kexp(V$AB, lA)))
  s_f <- c(A = mean(Matrix::rowSums(V$AA)[lA]), B = mean(Matrix::rowSums(V$BB)[lB]))
  vk <- c(vrai$dge, vrai$ige_intra_k, vrai$ige_inter_k)
  eff <- c("dge", "intra", "inter")
  for (sp in c("A", "B")) {
    nm <- paste0("g", sp)
    S <- matrix(NA_real_, 3, 3)
    for (a in 1:3) for (b in 1:a) {
      lab <- sprintf("%s[%d,%d]", nm, a, b)
      if (lab %in% rownames(v)) S[a, b] <- S[b, a] <- v[lab, "component"]
    }
    Bl <- f$blups[[nm]][if (sp == "A") gA else gB, , drop = FALSE]
    Bl <- cbind(Bl, Bl[, 1] + s_f[[sp]] * Bl[, 2])
    for (a in 1:3) base[[sprintf("err_%s_%s", eff[a], sp)]] <- log10(S[a, a] * k_f[[sp]][a] / vk[a])
    for (a in 1:4) base[[sprintf("prec_%s_%s", c(eff, "tbv")[a], sp)]] <- stats::cor(Bl[, a], sim$vraies[[sp]][, a])
  }
  base$secondes <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  base
}

options(remlax.threads = 1L)
pl <- plan()
sim <- simuler(verites[[nv]])
cat(sprintf("tache %d : verite %s, graine %d, n = %d, %d bacs, %d ajustements sur %d coeurs\n",
            i_tache, nv, graine, n, n_bac, nrow(pl), coeurs))
cat(sprintf("k AA %.3f BA %.3f BB %.3f AB %.3f ; s A %.3f B %.3f\n",
            sim$k[1], sim$k[2], sim$k[3], sim$k[4], sim$s[1], sim$s[2]))
t1 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(pl)), function(j) ajuster(pl[j, ], sim), mc.cores = coeurs,
                          mc.preschedule = FALSE)
ok <- vapply(res, is.data.frame, TRUE)
res[!ok] <- lapply(which(!ok), function(j) data.frame(pl[j, ], logLik = NA_real_, conv = NA,
                                                     secondes = NA_real_, erreur = "processus perdu"))
nm_all <- unique(unlist(lapply(res, names)))
res <- do.call(rbind, lapply(res, function(x) { x[setdiff(nm_all, names(x))] <- NA; x[nm_all] }))
res <- cbind(verite = nv, graine = graine, res)
utils::write.csv(res, file.path(dossier, sprintf("verite_%s_graine_%02d.csv", nv, graine)), row.names = FALSE)
cat(sprintf("fini en %.1f min ; %d echec(s), %d non converge(s)\n",
            as.numeric(difftime(Sys.time(), t1, units = "mins")), sum(is.na(res$logLik)),
            sum(!res$conv, na.rm = TRUE)))
