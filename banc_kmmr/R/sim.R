# ==============================================================================
# Cas simules du banc absolu : temps contre n, contre p, et regime dim(C) >> n.
#
# Nom d'un cas : sim_q<q>_t<t>_r<r>_s<struct>_x<x>_o<o>
#   q       genotypes phenotypes
#   t       environnements (caracteres) ; t = 1 donne un GBLUP
#   r       repetitions par genotype et par environnement observe
#   struct  structure G x E : id (t = 1), diag, fa1, fa2, fa3, us
#   x       genotypes NON phenotypes par genotype phenotype (ensemble de
#           prediction) : la parente porte q (1 + x) genotypes, C aussi
#   o       pourcentage des couples (genotype, environnement) observes
#           (essais clairsemes) ; 100 = plan complet
#
# Modele commun aux trois outils :
#   y = mu_env + g + e,  g ~ N(0, Sigma_g (x) K),  e ~ N(0, sigma2_env) par env.
# Parente : VanRaden sur 3000 marqueurs tires dans 10 sous-populations,
# K = 0,98 G / mean(diag G) + 0,02 I (definie positive meme si q > marqueurs).
#
# Les jeux sont ecrits une fois par `preparer_sim()` dans <dossier>/sim/, pour
# que les trois outils et les deux grappes lisent les memes nombres.
# ==============================================================================

lire_nom_sim <- function(cas) {
  v <- function(cle) sub(paste0(".*_", cle, "([^_]+).*"), "\\1", cas)
  list(q = as.integer(v("q")), t = as.integer(v("t")), r = as.integer(v("r")),
       s = v("s"), x = as.numeric(v("x")), o = as.numeric(v("o")))
}

# parente de q_tot genotypes, graine fixee par q_tot : tous les cas de meme
# q_tot partagent la meme matrice
simuler_K <- function(q_tot, m = 3000L, n_pop = 10L) {
  set.seed(1000L + q_tot)
  p0 <- stats::runif(m, 0.05, 0.5)
  pop <- sample.int(n_pop, q_tot, replace = TRUE)
  # frequences par sous-population : derive de type Balding-Nichols, Fst 0,1
  Fst <- 0.1
  P <- sapply(p0, function(p) stats::rbeta(n_pop, p * (1 - Fst) / Fst, (1 - p) * (1 - Fst) / Fst))
  P <- pmin(pmax(P, 0.01), 0.99)
  M <- matrix(stats::rbinom(q_tot * m, 2L, P[pop, ]), q_tot, m)
  f <- colMeans(M) / 2
  W <- sweep(M, 2L, 2 * f)
  G <- tcrossprod(W) / (2 * sum(f * (1 - f)))
  K <- 0.98 * G / mean(diag(G)) + 0.02 * diag(q_tot)
  noms <- sprintf("G%05d", seq_len(q_tot))
  dimnames(K) <- list(noms, noms)
  K
}

sigma_g <- function(t) {
  if (t == 1L) return(matrix(1))
  set.seed(77L + t)
  L <- cbind(stats::runif(t, 0.6, 0.95), stats::rnorm(t, 0, 0.35))
  S <- tcrossprod(L)
  S <- S + diag(pmax(1 - diag(S), 0.05))
  D <- diag(1 / sqrt(diag(S)))
  D %*% S %*% D                                  # variances genetiques 1
}

simuler_jeu <- function(cas, K) {
  a <- lire_nom_sim(cas)
  q <- a$q; t <- a$t
  q_tot <- nrow(K)
  set.seed(sum(utf8ToInt(cas)))
  noms <- rownames(K)
  obs <- noms[seq_len(q)]                       # les q premiers sont phenotypes
  # valeurs genetiques : chol(K)' Z chol(Sigma)
  Lk <- chol(K[obs, obs])
  Sg <- sigma_g(t)
  # sigma_g() refixe la graine quand t > 1 : sans cette ligne, les jeux _gN
  # multicaracteres etaient des copies de g1. Les noms sans suffixe gardent
  # l'ancien tirage pour que les mesures g1 restent valides. Les fichiers _g2 et
  # _g3 multicaracteres ecrits avant ce correctif sont ranges dans sim/copies_g1.
  if (grepl("_g[0-9]+$", cas)) set.seed(sum(utf8ToInt(cas)) + 1L)
  U <-crossprod(Lk, matrix(stats::rnorm(q * t), q, t)) %*% chol(Sg)
  # couples observes : chaque genotype dans max(1, round(o t / 100)) env
  k_env <- max(1L, round(a$o * t / 100))
  cellules <- do.call(rbind, lapply(seq_len(q), function(i)
    cbind(i, if (k_env >= t) seq_len(t) else sort(sample.int(t, k_env)))))
  cellules <- cellules[rep(seq_len(nrow(cellules)), each = a$r), , drop = FALSE]
  sig_e <- seq(1.5, 2.5, length.out = t)
  env <- factor(sprintf("E%02d", cellules[, 2]), levels = sprintf("E%02d", seq_len(t)))
  mu <- seq(5, 10, length.out = t)
  y <- mu[cellules[, 2]] + U[cellules] + stats::rnorm(nrow(cellules), 0, sqrt(sig_e[cellules[, 2]]))
  d <- data.frame(GY = y, Env = env, Gen = factor(obs[cellules[, 1]], levels = noms))
  d <- d[order(d$Env, d$Gen), ]
  d$Unit <- factor(sprintf("U%07d", seq_len(nrow(d))))
  rownames(d) <- NULL
  d
}

fichier_K <- function(dossier, q_tot) file.path(dossier, "sim", sprintf("K_q%d.rds", q_tot))
fichier_d <- function(dossier, cas) file.path(dossier, "sim", paste0(cas, ".rds"))

preparer_sim <- function(cas_liste, dossier) {
  dir.create(file.path(dossier, "sim"), showWarnings = FALSE, recursive = TRUE)
  for (cas in cas_liste) {
    a <- lire_nom_sim(cas)
    q_tot <- as.integer(round(a$q * (1 + a$x)))
    fk <- fichier_K(dossier, q_tot)
    if (!file.exists(fk)) {
      t0 <- Sys.time()
      saveRDS(simuler_K(q_tot), fk, compress = FALSE)
      cat(sprintf("K q_tot = %d : %.1f s\n", q_tot, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    fd <- fichier_d(dossier, cas)
    if (!file.exists(fd)) {
      d <- simuler_jeu(cas, readRDS(fk))
      saveRDS(d, fd)
      cat(sprintf("%s : n = %d, q_tot = %d\n", cas, nrow(d), q_tot))
    }
  }
}

charger_sim <- function(cas, dossier) {
  a <- lire_nom_sim(cas)
  q_tot <- as.integer(round(a$q * (1 + a$x)))
  d <- readRDS(fichier_d(dossier, cas))
  K <- readRDS(fichier_K(dossier, q_tot))
  rang <- if (startsWith(a$s, "fa")) as.integer(sub("fa", "", a$s)) else 0L
  list(cas = cas, d = d, K = K, type = "sim", t = a$t, ge = a$s, rang = rang)
}
