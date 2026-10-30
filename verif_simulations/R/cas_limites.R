# ==============================================================================
# Cas limites : chaque structure placee a la frontiere de son espace de
# parametres, ou pres d'elle. Rang incomplet d'une us, variance specifique nulle
# d'un fa, variance nulle, correlation ar1 ou cor proche de 1, IGE et IEE nuls,
# rang surestime (rr, fa), echelle de y extreme, terme a 4 niveaux.
#
# Les criteres restent ceux de la suite. Une composante dont la valeur vraie est
# sur la borne (variance nulle) est listee dans `bord` : son erreur-type par le
# Hessien n'a pas de sens, et son |z| est rapporte sans trancher.
#
# Toutes les voies MME sont essayees : inversion automatique, dense, creuse et
# hybride.
# ==============================================================================

M_LIM <- list(m_fit("lbfgsb"), m_fit("ai"), m_mme(), m_mme(inversion = "dense"),
              m_mme(inversion = "creuse"), m_mme(inversion = "hybrid"))

limite <- function(cas, bord = character(0), note = NULL) {
  cas$famille <- "cas limites"
  cas$bord <- bord
  if (!is.null(note)) cas$note <- note
  cas
}

# --- us de rang incomplet ---------------------------------------------------------
lim_us_rang1 <- function(graine) {
  l <- c(0.7, 0.5, -0.4, 0.6)
  Sr <- matrix(0.25, 4, 4); diag(Sr) <- 1
  limite(cas_mt(graine, "limite_us_rang1_t4", 4, "us", 0L, l %*% t(l), "us", Sr, moteurs = M_LIM),
         note = "Sigma genetique us de rang 1, residuelle us : optimum au bord du cone")
}

lim_us8 <- function(graine) {
  set.seed(graine + 7L)
  L <- matrix(stats::rnorm(16, sd = 0.5), 8)
  Sr <- matrix(0.3, 8, 8); diag(Sr) <- 1
  limite(cas_mt(graine, "limite_us8_presque_rang2", 8, "us", 0L, L %*% t(L) + diag(0.01, 8), "us", Sr,
                moteurs = M_LIM, n_g = 100),
         note = "us(8) genetique presque de rang 2 (LL' + 0,01 I), comme le banc absolu ; 100 genotypes")
}

lim_deux_us <- function(graine) {
  set.seed(graine)
  t <- 3; D <- mt_dispositif(t, n_g = 150)
  fam <- factor(sprintf("f%02d", (as.integer(D$d$g) - 1L) %/% 5L + 1L))
  Zf <- lapply(D$tn, function(k) { r <- which(D$d$trait == k)
    Matrix::sparseMatrix(i = r, j = as.integer(fam[r]), x = 1, dims = c(D$n, nlevels(fam)),
                         dimnames = list(NULL, levels(fam))) })
  names(Zf) <- D$tn
  Lg <- matrix(c(0.6, 0.4, -0.3, 0, 0.4, 0.3), 3); Sg <- Lg %*% t(Lg)
  lf <- c(0.5, -0.4, 0.3); Sf <- lf %*% t(lf)
  Sr <- matrix(c(1, 0.3, 0.2, 0.3, 0.8, -0.1, 0.2, -0.1, 0.6), 3)
  V_de <- function(par) V_terme(D$Zl, par$sig$g) + V_terme(Zf, par$sig$fam) +
    R_section(par$res[[1]], D$tr, D$Cu)
  y <- simuler_y(V_de(list(sig = list(g = Sg, fam = Sf), res = list(Sr))), D$X, seq_len(t))
  m <- rx_model(y, D$X, list(rx_term("g", D$Zl, struct = "us"), rx_term("fam", Zf, struct = "us")),
                rx_residual("us", trait = D$d$trait, unit = D$d$plante))
  vrai <- c(noms_sigma("g", Sg), noms_sigma("fam", Sf), noms_sigma("residuelle", Sr))
  limite(cas_def("limite_deux_us_au_bord", "cas limites", list(principal = list(model = m, V_de = V_de)),
                 M_LIM, vrai),
         note = "deux termes us au bord a la fois : genotype us(3) de rang 2, famille us(3) de rang 1")
}

lim_us_cor_forte <- function(graine) {
  R <- matrix(0.97, 3, 3); diag(R) <- 1; sd <- c(0.8, 0.6, 0.7)
  limite(cas_mt(graine, "limite_us_correlation_097", 3, "us", 0L, diag(sd) %*% R %*% diag(sd),
                "us", matrix(c(1, 0.3, 0.2, 0.3, 0.8, -0.1, 0.2, -0.1, 0.6), 3), moteurs = M_LIM),
         note = "us(3) de rang plein mais correlations genetiques de 0,97 : interieur mal conditionne")
}

# --- fa, rr : variance specifique nulle, rang surestime ----------------------------
lim_fa1_psi0 <- function(graine)
  limite(cas_mt(graine, "limite_fa1_psi_nul", 4, "fa", 1L,
                fa_sigma(matrix(c(0.7, 0.5, -0.4, 0.6)), c(0.2, 0, 0.25, 0.15)),
                "diag", diag(c(0.8, 0.6, 0.7, 0.5)), moteurs = M_LIM),
         note = "fa rang 1, variance specifique du caractere 2 nulle (cas de Heywood)")

lim_fa2_rang1 <- function(graine)
  limite(cas_mt(graine, "limite_fa2_vrai_rang1", 5, "fa", 2L,
                fa_sigma(matrix(c(0.7, 0.5, -0.3, 0.6, 0.4)), c(0.15, 0.2, 0.2, 0.15, 0.25)),
                "diag", diag(c(0.8, 0.6, 0.7, 0.5, 0.6)), moteurs = M_LIM, n_g = 100),
         note = "fa rang 2 ajuste, verite de rang 1 : second facteur nul, direction plate")

lim_rr2_rang1 <- function(graine) {
  l <- c(0.7, 0.5, -0.3, 0.6)
  limite(cas_mt(graine, "limite_rr2_vrai_rang1", 4, "rr", 2L, l %*% t(l),
                "diag", diag(c(0.8, 0.6, 0.7, 0.5)), moteurs = M_LIM),
         note = "rr rang 2 ajuste, verite de rang 1 sans variance specifique")
}

# --- variances nulles ----------------------------------------------------------------
lim_diag0 <- function(graine)
  limite(cas_mt(graine, "limite_diag_variance_nulle", 3, "diag", 0L, diag(c(0.5, 0, 0.3)),
                "diag", diag(c(1, 0.7, 0.5)), moteurs = M_LIM),
         bord = "g[2,2]", note = "diag genetique, variance du caractere 2 nulle")

lim_iid0 <- function(graine) {
  set.seed(graine)
  q <- 60; nk <- sample(6:10, q, replace = TRUE)
  g <- factor(rep(sprintf("g%02d", seq_len(q)), nk)); n <- length(g)
  X <- cbind(1, stats::rnorm(n)); Zg <- list(inc(g))
  V_de <- function(par) V_terme(Zg, par$sig$g) + R_section(par$res[[1]], NULL, diag(n))
  y <- simuler_y(V_de(list(sig = list(g = 0), res = list(1.0))), X, c(10, 0.5))
  m <- rx_model(y, X, list(rx_term("g", g)), rx_residual())
  limite(cas_def("limite_iid_variance_nulle", "cas limites", list(principal = list(model = m, V_de = V_de)),
                 M_LIM, c(g = 0, residuelle = 1)),
         bord = "g", note = "un facteur aleatoire de variance vraie nulle")
}

# --- correlations de niveaux proches de la borne -----------------------------------
lim_ar1 <- function(graine, rho, nom) {
  set.seed(graine)
  q <- 100; plot <- factor(rep(sprintf("x%03d", 1:q), each = 4)); n <- length(plot)
  X <- matrix(1, n, 1); Zp <- list(inc(plot))
  V_de <- function(par) V_terme(Zp, par$sig$plot, C_ar1(q, par$niv[["plot!ar1"]])) +
    par$res[[1]][1, 1] * diag(n)
  vrai <- c(plot = 1, "plot!ar1" = rho, residuelle = 0.5)
  y <- simuler_y(V_de(list(sig = list(plot = 1), res = list(matrix(0.5)), niv = vrai["plot!ar1"])), X, 0)
  m <- rx_model(y, X, list(rx_term("plot", plot, level = "ar1")), rx_residual())
  limite(cas_def(nom, "cas limites", list(principal = list(model = m, V_de = V_de)), M_LIM, vrai),
         note = sprintf("terme ar1 sur un transect de 100 parcelles, rho = %g", rho))
}
lim_ar1_fort <- function(graine) lim_ar1(graine, 0.97, "limite_ar1_rho_097")
lim_ar1_neg <- function(graine) lim_ar1(graine, -0.6, "limite_ar1_rho_negatif")

lim_res_ar1ar1 <- function(graine) {
  set.seed(graine)
  d <- expand.grid(col = 1:20, row = 1:30)
  d$g <- factor(sample(rep(sprintf("g%03d", 1:100), 6)))
  d$row <- factor(d$row); d$col <- factor(d$col)
  n <- nrow(d); X <- matrix(1, n, 1); Zg <- list(inc(d$g))
  rr <- as.integer(as.character(d$row)); cc <- as.integer(as.character(d$col))
  C_obs <- function(r1, r2) r1^abs(outer(rr, rr, "-")) * r2^abs(outer(cc, cc, "-"))
  V_de <- function(par) V_terme(Zg, par$sig$g) +
    R_section(par$res[[1]], NULL, C_obs(par$niv[["residuelle!ar1ar11"]], par$niv[["residuelle!ar1ar12"]]))
  vrai <- c(g = 0.4, residuelle = 1, "residuelle!ar1ar11" = 0.95, "residuelle!ar1ar12" = 0.9)
  y <- simuler_y(V_de(list(sig = list(g = 0.4), res = list(1),
                           niv = vrai[c("residuelle!ar1ar11", "residuelle!ar1ar12")])), X, 1)
  d$y <- y
  m <- rx_model(y, X, list(rx_term("g", d$g)), parse_res(~ ar1(row):ar1(col), d))
  limite(cas_def("limite_residuelle_ar1xar1_forte", "cas limites",
                 list(principal = list(model = m, V_de = V_de)), M_LIM, vrai),
         note = "residuelle ar1(row):ar1(col) avec rho = 0,95 et 0,90 sur un champ 30 x 20")
}

lim_cor_forte <- function(graine) {
  set.seed(graine)
  D <- ge_dispositif()
  V_de <- function(par) V_terme(D$Zl, par$sig$ge[1, 1] * diag(60), C_cor(4, par$niv[["ge!cor"]])) +
    par$res[[1]][1, 1] * diag(D$n)
  vrai <- c("ge[1,1]" = 0.8, "ge!cor" = 0.97, residuelle = 0.5)
  y <- simuler_y(V_de(list(sig = list(ge = matrix(0.8)), res = list(matrix(0.5)), niv = vrai["ge!cor"])),
                 D$X, c(1, 0.5, -0.5, 0.2))
  m <- rx_model(y, D$X, list(suppressMessages(rx_term("ge", D$Zl, struct = "iid", level = "cor"))), rx_residual())
  limite(cas_def("limite_cor_097", "cas limites", list(principal = list(model = m, V_de = V_de)), M_LIM, vrai),
         note = "correlation uniforme de 0,97 entre 4 milieux")
}

# --- IGE ---------------------------------------------------------------------------
lim_ige0 <- function(graine) {
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
  s2 <- c(dge = 0.3, ige = 0, residuelle = 0.6)
  X <- matrix(1, n, 1)
  V_de <- function(par) V_terme(list(Zd), par$sig$dge) + V_terme(list(Zn), par$sig$ige) +
    par$res[[1]][1, 1] * diag(n)
  y <- simuler_y(V_de(list(sig = list(dge = 0.3, ige = 0), res = list(matrix(0.6)))), X, 10)
  m <- rx_model(y, X, list(rx_term("dge", list(Zd)), rx_term("ige", list(Zn))), rx_residual())
  limite(cas_def("limite_ige_nul", "cas limites", list(principal = list(model = m, V_de = V_de)), M_LIM, s2),
         bord = "ige", note = "monoculture 30 x 20, DGE iid + IGE iid de variance vraie nulle")
}

lim_melange <- function(graine) {
  a <- c(1, 0); b <- c(-0.6, 0.8); cc <- c(0.6, 0.8)
  R <- rbind(a, b, cc) %*% t(rbind(a, b, cc))
  limite(cas_ige_melange(graine, nom = "limite_melange_iee_nul_us_rang2", part_iee = 0, Rg = R),
         bord = "iee",
         note = paste("melange deux especes : IEE de variance nulle, et (DGE, IGE intra, IGE inter)",
                      "us(3) de rang 2 pour chaque espece"))
}

# --- echelle de y et peu de niveaux ------------------------------------------------
lim_echelle <- function(graine, f, nom, n_env = 12) {
  set.seed(graine)
  d <- expand.grid(g = sprintf("g%02d", 1:50), env = sprintf("e%02d", seq_len(n_env)))
  d <- d[sort(sample(nrow(d), round(0.92 * nrow(d)))), ]
  n <- nrow(d); X <- matrix(1, n, 1)
  Zg <- list(inc(d$g)); Ze <- list(inc(d$env))
  V_de <- function(par) V_terme(Zg, par$sig$g) + V_terme(Ze, par$sig$env) + par$res[[1]][1, 1] * diag(n)
  vrai <- c(g = 0.5, env = 0.8, residuelle = 1.0) * f^2
  y <- simuler_y(V_de(list(sig = list(g = vrai[["g"]], env = vrai[["env"]]), res = list(matrix(vrai[["residuelle"]])))),
                 X, 5 * f)
  m <- rx_model(y, X, list(rx_term("g", d$g), rx_term("env", d$env)), rx_residual())
  limite(cas_def(nom, "cas limites", list(principal = list(model = m, V_de = V_de)), M_LIM, vrai),
         note = sprintf("genotype x milieu croises, %d milieux, y a l'echelle %g", n_env, f))
}
lim_echelle_grande <- function(graine) lim_echelle(graine, 1e3, "limite_echelle_1e3")
lim_echelle_petite <- function(graine) lim_echelle(graine, 1e-3, "limite_echelle_1em3")
lim_peu_niveaux <- function(graine) lim_echelle(graine, 1, "limite_terme_4_niveaux", n_env = 4)

CATALOGUE_LIMITES <- list(
  limite_us_rang1_t4 = lim_us_rang1, limite_us8_presque_rang2 = lim_us8,
  limite_deux_us_au_bord = lim_deux_us, limite_us_correlation_097 = lim_us_cor_forte,
  limite_fa1_psi_nul = lim_fa1_psi0, limite_fa2_vrai_rang1 = lim_fa2_rang1,
  limite_rr2_vrai_rang1 = lim_rr2_rang1, limite_diag_variance_nulle = lim_diag0,
  limite_iid_variance_nulle = lim_iid0, limite_ar1_rho_097 = lim_ar1_fort,
  limite_ar1_rho_negatif = lim_ar1_neg, limite_residuelle_ar1xar1_forte = lim_res_ar1ar1,
  limite_cor_097 = lim_cor_forte, limite_ige_nul = lim_ige0,
  limite_melange_iee_nul_us_rang2 = lim_melange, limite_echelle_1e3 = lim_echelle_grande,
  limite_echelle_1em3 = lim_echelle_petite, limite_terme_4_niveaux = lim_peu_niveaux)
CATALOGUE <- c(CATALOGUE, CATALOGUE_LIMITES)
