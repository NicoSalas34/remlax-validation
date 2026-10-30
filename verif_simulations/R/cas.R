# ==============================================================================
# Catalogue des cas simules. Chaque constructeur prend une graine et rend :
#   nom, famille, note
#   modeles  liste nommee ; chaque entree porte `model` (rx_model), `V_de`
#            (fonction par -> V, ecrite en R avec la geometrie vraie) et, pour un
#            decalage, `M_de` (par -> M)
#   moteurs  liste de specifications (id, type, optimizer, modele, erreur_attendue)
#   vrai     composantes vraies, nommees comme varcomp (fit$varcomp$nom)
#   derive   NULL, ou fonction (par, type) -> data.frame(nom, est, se) pour des
#            quantites derivees (correlations de corg)
# La simulation passe par simuler_y(V_vrai) : la V de l'arbitre, evaluee aux
# parametres vrais.
# ==============================================================================

m_fit <- function(opt, modele = "principal", erreur_attendue = NULL)
  list(id = opt, type = "rx_fit", optimizer = opt, modele = modele, erreur_attendue = erreur_attendue)
# moteur MME : rx_fit(engine = "mme"), inversion automatique ou creuse forcee
m_mme <- function(modele = "principal", inversion = "auto")
  list(id = if (inversion == "auto") "mme" else paste0("mme_", inversion), type = "mme",
       optimizer = "ai", inversion = inversion, modele = modele)
# moteur MME par L-BFGS-B, inversion creuse. Cette voie remplace l'ancien moteur
# RTMB (rx_fit_sparse), retire de remlax en octobre 2026, sur son perimetre :
# parente donnee par sa precision (Kinv) et champ ar1 separable.
m_mme_lbfgsb <- function(modele = "principal")
  list(id = "mme_lbfgsb", type = "mme", optimizer = "lbfgsb", inversion = "creuse", modele = modele)
M_STD <- list(m_fit("lbfgsb"), m_fit("ai"), m_mme(), m_mme(inversion = "creuse"))
M_STD_MME <- M_STD

cas_def <- function(nom, famille, modeles, moteurs, vrai, derive = NULL, note = "")
  list(nom = nom, famille = famille, modeles = modeles, moteurs = moteurs, vrai = vrai,
       derive = derive, note = note)

# Sigma vraie depuis ses composantes nommees
S_de <- function(vals, nom, t) {
  if (t == 1L) return(matrix(vals[[nom]], 1, 1))
  S <- matrix(0, t, t)
  for (i in seq_len(t)) for (j in seq_len(i)) S[i, j] <- S[j, i] <- vals[[sprintf("%s[%d,%d]", nom, i, j)]]
  S
}

# ------------------------------------------------------------------------------
# A. Termes simples, parente
# ------------------------------------------------------------------------------
cas_iid1 <- function(graine) {
  set.seed(graine)
  q <- 60; nk <- sample(6:10, q, replace = TRUE)
  g <- factor(rep(sprintf("g%02d", seq_len(q)), nk)); n <- length(g)
  x <- stats::rnorm(n); X <- cbind(1, x)
  Zg <- list(inc(g))
  V_de <- function(par) V_terme(Zg, par$sig$g) + R_section(par$res[[1]], NULL, diag(n))
  vrai <- c(g = 0.6, residuelle = 1.0)
  y <- simuler_y(V_de(list(sig = list(g = 0.6), res = list(1.0))), X, c(10, 0.5))
  m <- rx_model(y, X, list(rx_term("g", g)), rx_residual())
  cas_def("iid_un_terme", "termes simples", list(principal = list(model = m, V_de = V_de)),
          M_STD_MME, vrai, note = "un facteur aleatoire, effectifs desequilibres")
}

cas_iid2 <- function(graine) {
  set.seed(graine)
  d <- expand.grid(g = sprintf("g%02d", 1:50), env = sprintf("e%02d", 1:12))
  d <- d[sort(sample(nrow(d), round(0.92 * nrow(d)))), ]
  n <- nrow(d); X <- matrix(1, n, 1)
  Zg <- list(inc(d$g)); Ze <- list(inc(d$env))
  V_de <- function(par) V_terme(Zg, par$sig$g) + V_terme(Ze, par$sig$env) + par$res[[1]][1, 1] * diag(n)
  vrai <- c(g = 0.5, env = 0.8, residuelle = 1.0)
  y <- simuler_y(V_de(list(sig = list(g = 0.5, env = 0.8), res = list(matrix(1)))), X, 5)
  m <- rx_model(y, X, list(rx_term("g", d$g), rx_term("env", d$env)), rx_residual())
  cas_def("iid_deux_termes_croises", "termes simples", list(principal = list(model = m, V_de = V_de)),
          M_STD_MME, vrai, note = "genotype x milieu croises, 8 % de cellules manquantes")
}

cas_grm <- function(graine) {
  set.seed(graine)
  q <- 200; K <- grm_simulee(q, prefixe = "G")
  g <- factor(rep(rownames(K), each = 3), levels = rownames(K))
  keep <- sort(sample(length(g), 570)); g <- droplevels(g[keep])
  n <- length(g); X <- matrix(1, n, 1)
  Kg <- K[levels(g), levels(g)]
  Zg <- list(inc(g))
  V_de <- function(par) V_terme(Zg, par$sig$g, Kg) + par$res[[1]][1, 1] * diag(n)
  vrai <- c(g = 0.5, residuelle = 0.5)
  y <- simuler_y(V_de(list(sig = list(g = 0.5), res = list(matrix(0.5)))), X, 2)
  m <- rx_model(y, X, list(rx_term("g", g, K = K)), rx_residual())
  cas_def("parente_grm_dense", "parente", list(principal = list(model = m, V_de = V_de)),
          M_STD_MME, vrai, note = "GRM VanRaden (rx_grm, blending 0,05), 200 genotypes x 3 repetitions")
}

cas_pedigree <- function(graine) {
  set.seed(graine)
  # 100 fondateurs, 50 couples x 4 descendants, puis 100 couples entre familles
  # differentes x 2 : aucun apparentement entre parents, donc A^-1 de Henderson
  # sans consanguinite est exacte.
  nf <- 100; sire <- rep(0L, nf); dam <- rep(0L, nf); fam <- rep(NA_integer_, nf)
  for (k in 1:50) { sire <- c(sire, rep(2L * k - 1L, 4)); dam <- c(dam, rep(2L * k, 4)); fam <- c(fam, rep(k, 4)) }
  g2 <- (nf + 1):(nf + 200)
  f2 <- fam[g2]
  for (k in 1:100) {
    a <- sample(g2[f2 == ((k - 1) %% 50) + 1], 1)
    b <- sample(g2[f2 == (k %% 50) + 1], 1)
    sire <- c(sire, rep(a, 2)); dam <- c(dam, rep(b, 2))
  }
  N <- length(sire); ids <- sprintf("I%03d", seq_len(N))
  Ainv <- ainv_henderson(sire, dam); dimnames(Ainv) <- list(ids, ids)
  A <- as.matrix(solve(Ainv)); A <- (A + t(A)) / 2; dimnames(A) <- list(ids, ids)
  a <- factor(ids, levels = ids); n <- N; X <- matrix(1, n, 1)
  Za <- list(inc(a))
  V_de <- function(par) V_terme(Za, par$sig$a, A) + par$res[[1]][1, 1] * diag(n)
  vrai <- c(a = 0.5, residuelle = 0.5)
  y <- simuler_y(V_de(list(sig = list(a = 0.5), res = list(matrix(0.5)))), X, 3)
  md <- rx_model(y, X, list(rx_term("a", a, K = A)), rx_residual())
  mc <- rx_model(y, X, list(rx_term("a", a, Kinv = Ainv)), rx_residual())
  cas_def("parente_pedigree_kinv", "parente",
          list(principal = list(model = md, V_de = V_de), creux = list(model = mc, V_de = V_de)),
          list(m_fit("lbfgsb"), m_fit("ai"), m_mme("creux"), m_mme("creux", "creuse"), m_mme_lbfgsb("creux")), vrai,
          note = paste("modele animal, 500 individus sur 3 generations ; A^-1 par Henderson",
                       "(fonction R, asreml::ainverse non utilise) ; Kinv seul : moteur MME",
                       "(AI et L-BFGS-B) ; rx_fit sur le meme modele avec K = A dense"))
}

# ------------------------------------------------------------------------------
# B. Multi-caracteres : structures de Sigma
# ------------------------------------------------------------------------------
# Format long, t caracteres, n_g genotypes x 3 plantes, ~5 % de lignes absentes.
mt_dispositif <- function(t, n_g = 120, reps = 3, frac_manq = 0.05) {
  tn <- paste0("T", seq_len(t))
  pl <- data.frame(plante = sprintf("p%04d", seq_len(n_g * reps)),
                   g = rep(sprintf("g%03d", seq_len(n_g)), each = reps))
  d <- merge(pl, data.frame(trait = tn), by = NULL)
  d <- d[order(d$trait, d$plante), ]
  d <- d[sort(sample(nrow(d), round((1 - frac_manq) * nrow(d)))), ]
  d$trait <- factor(d$trait, levels = tn); d$g <- factor(d$g); d$plante <- factor(d$plante)
  n <- nrow(d)
  # niveaux complets pour chaque caractere, meme si un genotype y manque
  Zl <- lapply(tn, function(k) { r <- which(d$trait == k)
    Matrix::sparseMatrix(i = r, j = as.integer(d$g[r]), x = 1, dims = c(n, nlevels(d$g)),
                         dimnames = list(NULL, levels(d$g))) })
  names(Zl) <- tn
  list(d = d, n = n, t = t, tn = tn, Zl = Zl, X = stats::model.matrix(~ trait - 1, d),
       tr = as.integer(d$trait), Cu = C_meme_unite(as.integer(d$plante)))
}

cas_mt <- function(graine, nom, t, struct, rank, Sg, res_struct, Sr, res_rank = 0L,
                   moteurs = M_STD, note = "", n_g = 120) {
  set.seed(graine)
  D <- mt_dispositif(t, n_g = n_g)
  V_de <- function(par) V_terme(D$Zl, par$sig$g) + R_section(par$res[[1]], D$tr, D$Cu)
  y <- simuler_y(V_de(list(sig = list(g = Sg), res = list(Sr))), D$X, seq_len(t))
  m <- rx_model(y, D$X, list(rx_term("g", D$Zl, struct = struct, rank = rank)),
                rx_residual(res_struct, trait = D$d$trait, unit = D$d$plante, rank = res_rank))
  vrai <- c((if (struct == "diag") noms_diag else noms_sigma)("g", Sg),
            (if (res_struct == "diag") noms_diag else noms_sigma)("residuelle", Sr))
  cas_def(nom, "multi-caracteres", list(principal = list(model = m, V_de = V_de)), moteurs, vrai,
          note = note)
}

fa_sigma <- function(L, psi) L %*% t(L) + diag(psi)
cas_mt_diag <- function(graine) cas_mt(graine, "mt_g_diag_r_diag", 3, "diag", 0L, diag(c(0.5, 0.8, 0.3)),
                                       "diag", diag(c(1, 0.7, 0.5)), moteurs = M_STD_MME,
                                       note = "3 caracteres, 120 genotypes x 3 plantes, 5 % manquants")
cas_mt_us <- function(graine) {
  Sg <- matrix(c(0.6, 0.2, -0.15, 0.2, 0.5, 0.1, -0.15, 0.1, 0.4), 3)
  Sr <- matrix(c(1, 0.3, 0.2, 0.3, 0.8, -0.1, 0.2, -0.1, 0.6), 3)
  cas_mt(graine, "mt_g_us_r_us", 3, "us", 0L, Sg, "us", Sr, moteurs = M_STD_MME,
         note = "Sigma genetique us, residuelle us entre caracteres d'une meme plante")
}
cas_mt_fa1 <- function(graine) cas_mt(graine, "mt_g_fa1", 4, "fa", 1L,
                                      fa_sigma(matrix(c(0.7, 0.5, -0.4, 0.6)), c(0.2, 0.3, 0.25, 0.15)),
                                      "diag", diag(c(0.8, 0.6, 0.7, 0.5)), moteurs = M_STD_MME,
                                      note = "fa rang 1, 4 caracteres")
cas_mt_fa2 <- function(graine) cas_mt(graine, "mt_g_fa2", 5, "fa", 2L,
                                      fa_sigma(matrix(c(0.7, 0.5, -0.3, 0.6, 0.4, 0, 0.5, 0.6, -0.4, 0.3), 5),
                                               c(0.15, 0.2, 0.2, 0.15, 0.25)),
                                      "diag", diag(c(0.8, 0.6, 0.7, 0.5, 0.6)), n_g = 100, moteurs = M_STD_MME,
                                      note = "fa rang 2, 5 caracteres (rang 2 non identifie a 4 caracteres)")
cas_mt_rr2 <- function(graine) {
  G <- matrix(c(0.7, 0.5, -0.3, 0.6, 0, 0.4, 0.5, -0.3), 4)
  cas_mt(graine, "mt_g_rr2", 4, "rr", 2L, G %*% t(G), "diag", diag(c(0.8, 0.6, 0.7, 0.5)),
         note = "rang reduit 2 sans variance specifique : Sigma genetique singuliere")
}
cas_mt_chol1 <- function(graine) {
  L <- diag(4); L[2, 1] <- 0.5; L[3, 2] <- -0.4; L[4, 3] <- 0.3
  cas_mt(graine, "mt_g_chol1", 4, "chol", 1L, L %*% diag(c(0.6, 0.4, 0.5, 0.3)) %*% t(L),
         "diag", diag(c(0.8, 0.6, 0.7, 0.5)), note = "Cholesky en bande d'ordre 1")
}
cas_mt_ante1 <- function(graine) {
  Ut <- diag(4); Ut[2, 1] <- 0.5; Ut[3, 2] <- -0.4; Ut[4, 3] <- 0.3
  Ai <- solve(Ut)
  cas_mt(graine, "mt_g_ante1", 4, "ante", 1L, Ai %*% diag(1 / c(1.2, 1.5, 1.3, 1.6)^2) %*% t(Ai),
         "diag", diag(c(0.8, 0.6, 0.7, 0.5)), note = "antedependance d'ordre 1")
}
cas_mt_corh <- function(graine) {
  sd <- c(0.8, 0.6, 0.7)
  cas_mt(graine, "mt_g_corh", 3, "corh", 0L, diag(sd) %*% C_cor(3, 0.4) %*% diag(sd),
         "diag", diag(c(0.8, 0.6, 0.7)), note = "variances heterogenes, correlation uniforme 0,4")
}
cas_mt_resfa1 <- function(graine) cas_mt(graine, "mt_residuelle_fa1", 4, "diag", 0L, diag(c(0.4, 0.5, 0.3, 0.6)),
                                         "fa", fa_sigma(matrix(c(0.6, 0.5, 0.4, -0.3)), c(0.5, 0.4, 0.6, 0.5)),
                                         res_rank = 1L, note = "residuelle fa rang 1 entre caracteres d'une plante")

# ------------------------------------------------------------------------------
# C. Residuelles structurees
# ------------------------------------------------------------------------------
parse_res <- function(...) getFromNamespace(".rx_parse_residual", "remlax")(...)

# Champ 20 x 30, un genotype par cellule (100 genotypes x 6), residuelle
# ar1(ligne) x ar1(colonne). `manq` : fraction de cellules retirees.
cas_res_ar1ar1 <- function(graine, manq = 0, nom = "residuelle_ar1xar1") {
  set.seed(graine)
  d <- expand.grid(col = 1:20, row = 1:30)
  d$g <- factor(sample(rep(sprintf("g%03d", 1:100), 6)))
  if (manq > 0) d <- d[sort(sample(nrow(d), round((1 - manq) * nrow(d)))), ]
  d$row <- factor(d$row); d$col <- factor(d$col)
  n <- nrow(d); X <- matrix(1, n, 1); Zg <- list(inc(d$g))
  rr <- as.integer(as.character(d$row)); cc <- as.integer(as.character(d$col))
  C_obs <- function(r1, r2) r1^abs(outer(rr, rr, "-")) * r2^abs(outer(cc, cc, "-"))
  V_de <- function(par) V_terme(Zg, par$sig$g) +
    R_section(par$res[[1]], NULL, C_obs(par$niv[["residuelle!ar1ar11"]], par$niv[["residuelle!ar1ar12"]]))
  vrai <- c(g = 0.4, residuelle = 1, "residuelle!ar1ar11" = 0.6, "residuelle!ar1ar12" = 0.4)
  y <- simuler_y(V_de(list(sig = list(g = 0.4), res = list(1),
                           niv = vrai[c("residuelle!ar1ar11", "residuelle!ar1ar12")])), X, 1)
  d$y <- y
  m <- rx_model(y, X, list(rx_term("g", d$g)), parse_res(~ ar1(row):ar1(col), d))
  cas_def(nom, "residuelles", list(principal = list(model = m, V_de = V_de)), M_STD, vrai,
          note = if (manq > 0) sprintf("%.0f %% de cellules absentes du champ 30 x 20", 100 * manq)
                 else "champ complet 30 x 20, residual = ~ ar1(row):ar1(col)")
}
cas_res_ar1ar1_manq <- function(graine) cas_res_ar1ar1(graine, manq = 0.05, nom = "residuelle_ar1xar1_manquants")

cas_dsum_ar1_iid <- function(graine) {
  set.seed(graine)
  a <- expand.grid(col = 1:15, row = 1:20); a$site <- "A"
  # site B dans l'etendue des lignes et colonnes de A : les facteurs row et col
  # sont crees sur toutes les donnees, et B n'y ajoute aucun niveau.
  b <- data.frame(col = rep(1:13, each = 20), row = rep(1:20, 13), site = "B")
  d <- rbind(a, b)
  d$g <- factor(sample(rep(sprintf("g%03d", 1:112), 5)))
  d$site <- factor(d$site); n <- nrow(d)
  X <- stats::model.matrix(~ site, d); Zg <- list(inc(d$g))
  ra <- which(d$site == "A"); rb <- which(d$site == "B")
  rr <- d$row[ra]; cc <- d$col[ra]
  V_de <- function(par) {
    V <- V_terme(Zg, par$sig$g)
    V[ra, ra] <- V[ra, ra] + R_section(par$res[[1]], NULL,
      par$niv[["A!ar1ar11"]]^abs(outer(rr, rr, "-")) * par$niv[["A!ar1ar12"]]^abs(outer(cc, cc, "-")))
    V[rb, rb] <- V[rb, rb] + par$res[[2]][1, 1] * diag(length(rb))
    V
  }
  vrai <- c(g = 0.5, A = 1.0, "A!ar1ar11" = 0.5, "A!ar1ar12" = 0.3, B = 0.6)
  y <- simuler_y(V_de(list(sig = list(g = 0.5), res = list(matrix(1), matrix(0.6)),
                           niv = vrai[c("A!ar1ar11", "A!ar1ar12")])), X, c(2, 1))
  d$row <- factor(d$row); d$col <- factor(d$col)
  m <- rx_model(y, X, list(rx_term("g", d$g)),
                parse_res(~ dsum(~ ar1(row):ar1(col) + units | site, levels = list("A", "B")), d))
  cas_def("dsum_ar1xar1_et_iid", "residuelles", list(principal = list(model = m, V_de = V_de)), M_STD, vrai,
          note = "deux sites, dsum(~ ar1(row):ar1(col) + units | site) : structures differentes")
}

cas_dsum_iid_iid <- function(graine) {
  set.seed(graine)
  d <- data.frame(site = factor(rep(c("A", "B"), c(300, 260))))
  d$g <- factor(sample(rep(sprintf("g%03d", 1:80), 7)))
  n <- nrow(d); X <- stats::model.matrix(~ site, d); Zg <- list(inc(d$g))
  ra <- which(d$site == "A"); rb <- which(d$site == "B")
  V_de <- function(par) {
    V <- V_terme(Zg, par$sig$g)
    V[ra, ra] <- V[ra, ra] + par$res[[1]][1, 1] * diag(length(ra))
    V[rb, rb] <- V[rb, rb] + par$res[[2]][1, 1] * diag(length(rb))
    V
  }
  vrai <- c(g = 0.5, A = 1.0, B = 0.4)
  y <- simuler_y(V_de(list(sig = list(g = 0.5), res = list(matrix(1), matrix(0.4)))), X, c(2, 1))
  m <- rx_model(y, X, list(rx_term("g", d$g)), parse_res(~ dsum(~ units | site), d))
  cas_def("dsum_iid_heterogene", "residuelles", list(principal = list(model = m, V_de = V_de)),
          M_STD_MME, vrai, note = "residuelle iid propre a chaque site, dsum(~ units | site)")
}

# Residuelle "none" : un terme de plante [I sur T1, I sur T2] us porte la
# residuelle. Jumeau : meme modele avec residuelle us(trait):units ; les deux
# doivent atteindre la meme logLik.
cas_res_none <- function(graine) {
  set.seed(graine)
  t <- 2; np <- 300; tn <- c("T1", "T2")
  d <- data.frame(plante = factor(rep(sprintf("p%03d", 1:np), t)), trait = factor(rep(tn, each = np)))
  d$g <- factor(rep(sprintf("g%03d", rep(1:100, each = 3)), t))
  n <- nrow(d); X <- stats::model.matrix(~ trait - 1, d)
  Zg <- lapply(tn, function(k) { r <- which(d$trait == k)
    Matrix::sparseMatrix(i = r, j = as.integer(d$g[r]), x = 1, dims = c(n, nlevels(d$g))) })
  Zp <- lapply(tn, function(k) { r <- which(d$trait == k)
    Matrix::sparseMatrix(i = r, j = as.integer(d$plante[r]), x = 1, dims = c(n, np)) })
  Sg <- matrix(c(0.5, 0.15, 0.15, 0.4), 2); Sp <- matrix(c(1, 0.35, 0.35, 0.8), 2)
  V_de <- function(par) V_terme(Zg, par$sig$g) + V_terme(Zp, par$sig$plante)
  V_jum <- function(par) V_terme(Zg, par$sig$g) +
    R_section(par$res[[1]], as.integer(d$trait), C_meme_unite(as.integer(d$plante)))
  y <- simuler_y(V_de(list(sig = list(g = Sg, plante = Sp))), X, c(1, 2))
  m <- rx_model(y, X, list(rx_term("g", Zg, struct = "us"), rx_term("plante", Zp, struct = "us")),
                rx_residual("none"))
  ren <- c("residuelle[1,1]" = "plante[1,1]", "residuelle[2,1]" = "plante[2,1]",
           "residuelle[2,2]" = "plante[2,2]")
  mj <- rx_model(y, X, list(rx_term("g", Zg, struct = "us")),
                 rx_residual("us", trait = d$trait, unit = d$plante))
  cas_def("residuelle_none_terme_plante", "residuelles",
          list(principal = list(model = m, V_de = V_de), jumeau = list(model = mj, V_de = V_jum)),
          list(m_fit("lbfgsb"), m_fit("ai"), m_mme(),
               c(m_fit("lbfgsb", "jumeau"), list(renommer = ren)), c(m_mme("jumeau"), list(renommer = ren))),
          c(noms_sigma("g", Sg), noms_sigma("plante", Sp)),
          note = "residual none, la plante us porte la residuelle ; jumeau a residuelle us(trait):units")
}

# ------------------------------------------------------------------------------
# D. Structures entre niveaux
# ------------------------------------------------------------------------------
cas_lvl_ar1 <- function(graine) {
  set.seed(graine)
  q <- 100; plot <- factor(rep(sprintf("x%03d", 1:q), each = 4)); n <- length(plot)
  X <- matrix(1, n, 1); Zp <- list(inc(plot))
  V_de <- function(par) V_terme(Zp, par$sig$plot, C_ar1(q, par$niv[["plot!ar1"]])) +
    par$res[[1]][1, 1] * diag(n)
  vrai <- c(plot = 1, "plot!ar1" = 0.7, residuelle = 0.5)
  y <- simuler_y(V_de(list(sig = list(plot = 1), res = list(matrix(0.5)), niv = vrai["plot!ar1"])), X, 0)
  m <- rx_model(y, X, list(rx_term("plot", plot, level = "ar1")), rx_residual())
  cas_def("niveaux_ar1", "structures entre niveaux", list(principal = list(model = m, V_de = V_de)),
          M_STD_MME, vrai, note = "transect de 100 parcelles x 4 sous-echantillons, terme ar1")
}

cas_lvl_ar1ar1 <- function(graine) {
  set.seed(graine)
  nr <- 20; nc <- 25
  d <- expand.grid(col = 1:nc, row = 1:nr)            # cellule = (row-1)*nc + col
  d$g <- factor(sample(rep(sprintf("g%03d", 1:100), 5)))
  n <- nrow(d); X <- matrix(1, n, 1)
  cell <- factor((d$row - 1) * nc + d$col, levels = seq_len(nr * nc))
  Zs <- list(inc(cell)); Zg <- list(inc(d$g))
  V_de <- function(par) V_terme(Zg, par$sig$g) +
    V_terme(Zs, par$sig$champ, C_ar1ar1(nr, nc, par$niv[["champ!ar1ar11"]], par$niv[["champ!ar1ar12"]])) +
    par$res[[1]][1, 1] * diag(n)
  vrai <- c(g = 0.4, champ = 0.8, "champ!ar1ar11" = 0.6, "champ!ar1ar12" = 0.4, residuelle = 0.3)
  y <- simuler_y(V_de(list(sig = list(g = 0.4, champ = 0.8), res = list(matrix(0.3)),
                           niv = vrai[c("champ!ar1ar11", "champ!ar1ar12")])), X, 0)
  m <- rx_model(y, X, list(rx_term("g", d$g), rx_term("champ", cell, level = "ar1ar1", dims = c(nr, nc))),
                rx_residual())
  cas_def("niveaux_ar1xar1_champ", "structures entre niveaux",
          list(principal = list(model = m, V_de = V_de)), c(M_STD_MME, list(m_mme_lbfgsb())), vrai,
          note = "champ aleatoire ar1 x ar1 sur 20 x 25 cellules, avec pepite et genotype")
}

# Familles metriques : 150 positions x 3 observations, terme `s` metrique.
cas_metrique <- function(graine, kind, par_vrai, dim2 = TRUE, nu = 0.5) {
  set.seed(graine)
  q <- 150
  co <- if (dim2) cbind(stats::runif(q, 0, 15), stats::runif(q, 0, 15)) else cbind(sort(stats::runif(q, 0, 50)))
  pos <- factor(rep(sprintf("s%03d", 1:q), each = 3)); n <- length(pos)
  X <- matrix(1, n, 1); Zs <- list(inc(pos))
  np <- length(par_vrai)
  noms_niv <- if (np == 1L) paste0("s!", kind) else paste0("s!", kind, seq_len(np))
  V_de <- function(par) V_terme(Zs, par$sig$s, C_metrique(kind, co, unname(par$niv[noms_niv]), nu = nu)) +
    par$res[[1]][1, 1] * diag(n)
  vrai <- c(s = 1, stats::setNames(par_vrai, noms_niv), residuelle = 0.4)
  y <- simuler_y(V_de(list(sig = list(s = 1), res = list(matrix(0.4)), niv = vrai[noms_niv])), X, 0)
  opts <- if (kind == "mtrn") getFromNamespace(".rx_mtrn_opts", "remlax")(list(nu = sprintf("%g F", nu))) else NULL
  m <- rx_model(y, X, list(rx_term("s", pos, level = kind, coord = co, opts = opts)), rx_residual())
  cas_def(paste0("metrique_", kind), "structures entre niveaux",
          list(principal = list(model = m, V_de = V_de)), M_STD, vrai,
          note = sprintf("%s, %s, 150 positions x 3 obs%s", kind, if (dim2) "2D" else "1D",
                         if (kind == "mtrn") sprintf(", nu fixe a %g", nu) else ""))
}
CAS_METRIQUES <- list(
  exp  = list(par = 0.85, d2 = FALSE), gau = list(par = 0.97, d2 = FALSE), lvr = list(par = 8, d2 = FALSE),
  iexp = list(par = 0.7, d2 = TRUE), igau = list(par = 0.9, d2 = TRUE), ieuc = list(par = 0.7, d2 = TRUE),
  sph  = list(par = 5, d2 = TRUE), cir = list(par = 5, d2 = TRUE),
  aexp = list(par = c(0.6, 0.8), d2 = TRUE), agau = list(par = c(0.85, 0.93), d2 = TRUE),
  mtrn = list(par = 1.5, d2 = TRUE, nu = 1.5))

# Genotype x milieu : t = 60 genotypes (Sigma iid partagee), q = 4 milieux
# correles entre eux par une structure de niveaux (cor ou corg).
ge_dispositif <- function() {
  d <- expand.grid(rep = 1:2, env = sprintf("e%d", 1:4), g = sprintf("g%02d", 1:60))
  d$env <- factor(d$env); d$g <- factor(d$g)
  n <- nrow(d)
  Zl <- lapply(levels(d$g), function(k) { r <- which(d$g == k)
    Matrix::sparseMatrix(i = r, j = as.integer(d$env[r]), x = 1, dims = c(n, 4)) })
  list(d = d, n = n, Zl = Zl, X = stats::model.matrix(~ env, d))
}
cas_lvl_cor <- function(graine) {
  set.seed(graine)
  D <- ge_dispositif()
  V_de <- function(par) V_terme(D$Zl, par$sig$ge[1, 1] * diag(60), C_cor(4, par$niv[["ge!cor"]])) +
    par$res[[1]][1, 1] * diag(D$n)
  vrai <- c("ge[1,1]" = 0.8, "ge!cor" = 0.5, residuelle = 0.5)
  y <- simuler_y(V_de(list(sig = list(ge = matrix(0.8)), res = list(matrix(0.5)), niv = vrai["ge!cor"])),
                 D$X, c(1, 0.5, -0.5, 0.2))
  m <- rx_model(y, D$X, list(suppressMessages(rx_term("ge", D$Zl, struct = "iid", level = "cor"))), rx_residual())
  cas_def("niveaux_cor_ge", "structures entre niveaux", list(principal = list(model = m, V_de = V_de)),
          M_STD, vrai, note = "60 genotypes x 4 milieux x 2 rep., correlation uniforme entre milieux")
}
cas_lvl_corg <- function(graine) {
  set.seed(graine)
  D <- ge_dispositif()
  Cv <- matrix(c(1, .6, .3, .1, .6, 1, .5, .2, .3, .5, 1, .4, .1, .2, .4, 1), 4)
  Cde <- function(th) C_corg_theta(th, 4)
  V_de <- function(par) V_terme(D$Zl, par$sig$ge[1, 1] * diag(60), Cde(par$theta[2:7])) +
    par$res[[1]][1, 1] * diag(D$n)
  # simulation : V vraie directement (corg n'a pas d'echelle naturelle dans varcomp)
  Vv <- V_terme(D$Zl, 0.8 * diag(60), Cv) + 0.5 * diag(D$n)
  y <- simuler_y(Vv, D$X, c(1, 0.5, -0.5, 0.2))
  m <- rx_model(y, D$X, list(suppressMessages(rx_term("ge", D$Zl, struct = "iid", level = "corg", order = 4L))),
                rx_residual())
  ij <- which(lower.tri(Cv), arr.ind = TRUE)
  vrai <- c("ge[1,1]" = 0.8, residuelle = 0.5,
            stats::setNames(Cv[ij], sprintf("C[%d,%d]", ij[, 1], ij[, 2])))
  derive <- function(par, type) {
    if (type != "rx_fit") return(NULL)
    f <- function(th) Cde(th[2:7])[ij]
    r <- delta_se(f, par)
    data.frame(nom = sprintf("C[%d,%d]", ij[, 1], ij[, 2]), est = r$val, se = r$se)
  }
  cas_def("niveaux_corg_ge", "structures entre niveaux", list(principal = list(model = m, V_de = V_de)),
          M_STD, vrai, derive = derive,
          note = "correlation generale entre 4 milieux ; correlations par delta-methode sur theta")
}

# ------------------------------------------------------------------------------
# E. Splines 2D (PS-ANOVA)
# ------------------------------------------------------------------------------
cas_spl2d <- function(graine) {
  set.seed(graine)
  d <- expand.grid(col = 1:20, row = 1:30)
  d$g <- factor(sample(rep(sprintf("g%03d", 1:100), 6)))
  n <- nrow(d)
  sp <- rx_spl2d(d$col, d$row, nseg = c(6, 8))
  X <- cbind(1, sp$X)
  Zg <- list(inc(d$g))
  noms <- vapply(sp$terms, `[[`, "", "name")
  Zsp <- lapply(sp$terms, function(tm) tm$Zl)
  names(Zsp) <- noms
  ech <- vapply(Zsp, function(z) mean(Matrix::rowSums(z[[1]]^2)), 0)
  part <- c(0.3, 0.3, 0.15, 0.15, 0.3)
  s2 <- stats::setNames(part / ech, noms)
  V_de <- function(par) {
    V <- V_terme(Zg, par$sig$g) + par$res[[1]][1, 1] * diag(n)
    for (k in noms) V <- V + V_terme(Zsp[[k]], par$sig[[k]])
    V
  }
  sig_v <- c(list(g = 0.4), lapply(as.list(s2), function(v) matrix(v)))
  y <- simuler_y(V_de(list(sig = sig_v, res = list(matrix(0.5)))), X, c(5, rep(0.3, ncol(sp$X))))
  m <- rx_model(y, X, c(list(rx_term("g", d$g)), sp$terms), rx_residual())
  cas_def("spline_spl2d", "splines", list(principal = list(model = m, V_de = V_de)), M_STD_MME,
          c(g = 0.4, s2, residuelle = 0.5),
          note = "rx_spl2d(nseg = c(6, 8)) : cinq termes PS-ANOVA, partie nulle en effets fixes")
}

# ------------------------------------------------------------------------------
# F. Decalage spatial (rx_lag)
# ------------------------------------------------------------------------------
cas_lag <- function(graine) {
  set.seed(graine)
  nr <- 15; nc <- 20; N <- nr * nc; q <- 2; tn <- c("T1", "T2")
  co <- expand.grid(c = seq_len(nc), r = seq_len(nr))
  W <- Matrix::Matrix((as.matrix(stats::dist(co, method = "maximum")) == 1) * 1, sparse = TRUE)
  dimnames(W) <- NULL
  g <- sample(rep(seq_len(60), length.out = N))
  trait <- factor(rep(tn, each = N), levels = tn); unit <- rep(seq_len(N), q)
  n <- N * q; X <- stats::model.matrix(~ trait - 1)
  Zg1 <- Matrix::sparseMatrix(i = seq_len(N), j = g, x = 1, dims = c(N, 60))
  zero <- Matrix::Matrix(0, N, 60, sparse = TRUE)
  ZD <- list(T1 = rbind(Zg1, zero), T2 = rbind(zero, Zg1))
  Wd <- as.matrix(W)
  M_de <- function(par) { Psi <- par$Psi; M <- matrix(0, n, n)
    for (a in 1:q) for (b in 1:q) M[(a - 1) * N + seq_len(N), (b - 1) * N + seq_len(N)] <- Psi[a, b] * Wd
    M }
  Cu <- C_meme_unite(unit); tri <- as.integer(trait)
  V_de <- function(par) V_terme(ZD, par$sig$gen) + R_section(par$res[[1]], tri, Cu)
  Sg <- matrix(c(1, 0.3, 0.3, 0.8), 2); Sr <- matrix(c(1, -0.2, -0.2, 0.9), 2)
  Psi <- matrix(c(0.06, 0.05, -0.04, 0.03), 2, dimnames = list(tn, tn))
  pv <- list(sig = list(gen = Sg), res = list(Sr), Psi = Psi)
  eta <- simuler_y(V_de(pv), X, c(3, 6))
  y <- as.numeric(solve(diag(n) - M_de(pv), eta))
  L <- rx_lag(W, trait, unit)
  m <- rx_model(y, X, list(rx_term("gen", ZD, struct = "us")),
                rx_residual("us", trait = trait, unit = unit), lag = L)
  vrai <- c(noms_sigma("gen", Sg), noms_sigma("residuelle", Sr),
            "Psi[T1<-T1]" = 0.06, "Psi[T1<-T2]" = -0.04, "Psi[T2<-T1]" = 0.05, "Psi[T2<-T2]" = 0.03)
  cas_def("decalage_rx_lag", "decalage spatial",
          list(principal = list(model = m, V_de = V_de, M_de = M_de)),
          list(m_fit("lbfgsb"), m_fit("ai", erreur_attendue = "lag"),
               modifyList(m_mme(), list(optimizer = "lbfgsb"))), vrai,
          note = "y = (I - Psi (x) W)^-1 (Xb + Zu + e), W reine 15 x 20, 2 caracteres, critere reml_y")
}

# ------------------------------------------------------------------------------
# G. Modeles DGE / IGE
# ------------------------------------------------------------------------------
# Melange de deux especes sur 2 bacs de 30 rangs x 20 colonnes, colonnes
# alternees (A impaires, B paires) ; 150 genotypes par espece x 4 plantes.
cas_ige_melange <- function(graine, nom = "ige_melange_deux_especes", part_iee = 0.05, Rg = NULL) {
  set.seed(graine)
  bac <- expand.grid(col = 1:20, row = 1:30)
  d <- rbind(transform(bac, bloc = "b1"), transform(bac, bloc = "b2"))
  d$sp <- ifelse(d$col %% 2 == 1, "A", "B")
  d$id <- sprintf("pl%04d", seq_len(nrow(d)))
  d$geno <- NA_character_
  for (s in c("A", "B")) { r <- which(d$sp == s); d$geno[r] <- sample(rep(sprintf("%s%03d", s, 1:150), 4)) }
  n <- nrow(d); ia <- which(d$sp == "A"); ib <- which(d$sp == "B")
  K <- list(A = grm_simulee(150, prefixe = "A"), B = grm_simulee(150, prefixe = "B"))
  nb <- rx_neighbourhood(coord = cbind(d$row, d$col), group = d$sp, block = d$bloc, level = d$geno,
                         id = d$id, rank = list("A<-A" = 2, "A<-B" = 1, "B<-A" = 1, "B<-B" = 2),
                         kernel = "none", normalise = TRUE, output = "both")
  lev <- list(A = sort(unique(d$geno[ia])), B = sort(unique(d$geno[ib])))
  rows <- list(A = ia, B = ib); autre <- c(A = "B", B = "A")
  Zgen <- lapply(c(A = "A", B = "B"), function(s) {
    r <- rows[[s]]; o <- autre[[s]]
    DGE <- Matrix::sparseMatrix(i = r, j = match(d$geno[r], lev[[s]]), x = 1, dims = c(n, 150),
                                dimnames = list(NULL, lev[[s]]))
    list(DGE = DGE,
         IGE_intra = placer(nb$level[[paste0(s, "<-", s)]][, lev[[s]]], r, n),
         IGE_inter = placer(nb$level[[paste0(o, "<-", s)]][, lev[[s]]], rows[[o]], n))
  })
  # IEE : un terme iid au niveau plante, toutes paires confondues
  Zi <- Matrix::Matrix(0, n, n, sparse = TRUE)
  for (a in c("A", "B")) for (b in c("A", "B"))
    Zi[rows[[a]], rows[[b]]] <- nb$unit[[paste0(a, "<-", b)]]
  colnames(Zi) <- d$id
  # Mise a l'echelle d'exposition : part de variance phenotypique visee par
  # composante, divisee par sa fonctionnelle d'exposition (rx_exposure).
  Sig <- list()
  for (s in c("A", "B")) {
    r <- rows[[s]]; o <- autre[[s]]
    ex_d <- rx_exposure(Zgen[[s]]$DGE[r, ], Zgen[[s]]$IGE_intra[r, ], K = K[[s]])
    ex_x <- rx_exposure(NULL, Zgen[[s]]$IGE_inter[rows[[o]], ], K = K[[s]])
    v <- c(0.30 / ex_d$d, 0.07 / ex_d$k, 0.05 / ex_x$k)
    R <- if (is.null(Rg)) matrix(c(1, -0.4, -0.2, -0.4, 1, 0.3, -0.2, 0.3, 1), 3) else Rg
    Sig[[paste0("g", s)]] <- diag(sqrt(v)) %*% R %*% diag(sqrt(v))
  }
  k_iee <- mean(Matrix::rowSums(Zi^2))
  s2_iee <- part_iee / k_iee
  Ks <- list(A = K$A[lev$A, lev$A], B = K$B[lev$B, lev$B])
  ZiZ <- as.matrix(Matrix::tcrossprod(Zi))
  tri <- as.integer(factor(d$sp))
  V_de <- function(par) V_terme(Zgen$A, par$sig$gA, Ks$A) + V_terme(Zgen$B, par$sig$gB, Ks$B) +
    par$sig$iee[1, 1] * ZiZ + R_section(par$res[[1]], tri, diag(n))
  Sr <- diag(c(0.5, 0.55))
  pv <- list(sig = list(gA = Sig$gA, gB = Sig$gB, iee = matrix(s2_iee)), res = list(Sr))
  X <- stats::model.matrix(~ sp + bloc, d)
  y <- simuler_y(V_de(pv), X, c(10, 4, 0.5))
  m <- rx_model(y, X, list(rx_term("gA", Zgen$A, K = K$A, struct = "us"),
                           rx_term("gB", Zgen$B, K = K$B, struct = "us"),
                           rx_term("iee", list(Zi), struct = "iid")),
                rx_residual("diag", trait = factor(d$sp)))
  vrai <- c(noms_sigma("gA", Sig$gA), noms_sigma("gB", Sig$gB), iee = s2_iee, noms_diag("residuelle", Sr))
  cas_def(nom, "DGE / IGE",
          list(principal = list(model = m, V_de = V_de)), M_STD_MME, vrai,
          note = paste("2 bacs 30 x 20, colonnes alternees ; DGE + IGE intra + IGE inter us (3 x 3) par",
                       "espece avec GRM ; IEE iid ; residuelle diag par espece ; rangs 2 (intra) et 1",
                       "(inter), normalisation L2 ; h2 0,30, IGE 7 et 5 %, IEE 5 %, r(DGE, IGE) -0,4"))
}

cas_ige_mono <- function(graine) {
  set.seed(graine)
  d <- expand.grid(col = 1:20, row = 1:30)
  d$id <- sprintf("pl%03d", seq_len(nrow(d)))
  d$geno <- sample(rep(sprintf("g%03d", 1:120), 5))
  n <- nrow(d)
  nb <- rx_neighbourhood(coord = cbind(d$row, d$col), group = rep("A", n), level = d$geno, id = d$id,
                         rank = 1, kernel = "none", normalise = TRUE)
  lev <- sort(unique(d$geno))
  Zd <- Matrix::sparseMatrix(i = seq_len(n), j = match(d$geno, lev), x = 1, dims = c(n, 120),
                             dimnames = list(NULL, lev))
  Zn <- nb$level[["A<-A"]][, lev]; dimnames(Zn) <- list(NULL, lev)
  k <- rx_exposure(Zd, Zn)$k
  s2 <- c(dge = 0.3, ige = 0.08 / k, residuelle = 0.6)
  X <- matrix(1, n, 1)
  V_de <- function(par) V_terme(list(Zd), par$sig$dge) + V_terme(list(Zn), par$sig$ige) +
    par$res[[1]][1, 1] * diag(n)
  y <- simuler_y(V_de(list(sig = list(dge = s2[["dge"]], ige = s2[["ige"]]), res = list(matrix(s2[["residuelle"]])))),
                 X, 10)
  m <- rx_model(y, X, list(rx_term("dge", list(Zd)), rx_term("ige", list(Zn))), rx_residual())
  cas_def("ige_univarie_une_espece", "DGE / IGE", list(principal = list(model = m, V_de = V_de)),
          M_STD_MME, s2, note = "monoculture 30 x 20, 120 genotypes x 5, DGE iid + IGE iid (rang 1, L2)")
}

# ------------------------------------------------------------------------------
CATALOGUE <- c(
  list(iid_un_terme = cas_iid1, iid_deux_termes_croises = cas_iid2,
       parente_grm_dense = cas_grm, parente_pedigree_kinv = cas_pedigree,
       mt_g_diag_r_diag = cas_mt_diag, mt_g_us_r_us = cas_mt_us, mt_g_fa1 = cas_mt_fa1,
       mt_g_fa2 = cas_mt_fa2, mt_g_rr2 = cas_mt_rr2, mt_g_chol1 = cas_mt_chol1,
       mt_g_ante1 = cas_mt_ante1, mt_g_corh = cas_mt_corh, mt_residuelle_fa1 = cas_mt_resfa1,
       residuelle_ar1xar1 = cas_res_ar1ar1, residuelle_ar1xar1_manquants = cas_res_ar1ar1_manq,
       dsum_ar1xar1_et_iid = cas_dsum_ar1_iid, dsum_iid_heterogene = cas_dsum_iid_iid,
       residuelle_none_terme_plante = cas_res_none,
       niveaux_ar1 = cas_lvl_ar1, niveaux_ar1xar1_champ = cas_lvl_ar1ar1,
       niveaux_cor_ge = cas_lvl_cor, niveaux_corg_ge = cas_lvl_corg),
  stats::setNames(lapply(names(CAS_METRIQUES), function(k) {
    s <- CAS_METRIQUES[[k]]
    function(graine) cas_metrique(graine, k, s$par, dim2 = s$d2, nu = s$nu %||% 0.5)
  }), paste0("metrique_", names(CAS_METRIQUES))),
  list(spline_spl2d = cas_spl2d, decalage_rx_lag = cas_lag,
       ige_melange_deux_especes = cas_ige_melange, ige_univarie_une_espece = cas_ige_mono))
