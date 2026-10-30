# ==============================================================================
# Familles du banc "nombre de parametres" : modeles multi-caracteres dont le
# nombre de parametres de variance croit avec le nombre t de caracteres, a
# nombre d'observations n FIXE. Meme contrat que familles.R.
#
#   mtus<t>   genotype us(t) + residuelle us(t) par unite     t(t+1) parametres
#   mtige<t>  DGE et IGE dans un bloc us(2t) sur les incidences
#             directe et de voisinage (a la maniere de mvD7),
#             + residuelle us(t) par plante                   t(2t+1) + t(t+1)/2
#   mtigee<t> mtige + IEE : effet environnemental indirect au
#             niveau de la PLANTE (incidence de voisinage
#             plante -> plante, us(t) entre caracteres). Autant
#             de niveaux que de plantes : m de l'ordre de n et
#             C presque pleine, comme dans mvD7               mtige + t(t+1)/2
#
# Sources apres familles.R (qui definit FAMILLES, graine_de, inc_niv, ...).
# ==============================================================================

T_PARAMS <- c(1L, 2L, 3L, 4L, 6L, 8L)
FAMILLES <- c(FAMILLES, paste0("mtus", T_PARAMS), paste0("mtige", T_PARAMS), paste0("mtigee", T_PARAMS))

# matrice de covariance : variances v, correlation uniforme r (definie positive
# pour r > -1/(t-1))
cov_unif <- function(v, r) { t <- length(v); R <- matrix(r, t, t); diag(R) <- 1; D <- diag(sqrt(v), t); D %*% R %*% D }

# V = Z (G (x) I_q) Z' + R par unite, R[trait_i, trait_j] pour deux lignes de la meme unite
.V_mt <- function(Zall, G, q, R, tr, un, t) {
  Vg <- as.matrix(Zall %*% Matrix::kronecker(methods::as(G, "CsparseMatrix"), Matrix::Diagonal(q)) %*% Matrix::t(Zall))
  Rm <- Matrix::sparseMatrix(i = seq_along(tr), j = tr, x = 1, dims = c(length(tr), t))
  Zu <- Matrix::sparseMatrix(i = seq_along(un), j = un, x = 1)
  Vg + as.matrix(Matrix::tcrossprod(Zu)) * as.matrix(Rm %*% R %*% Matrix::t(Rm))
}

fam_mtus <- function(n, graine, t) {
  set.seed(graine)
  nu <- n %/% t; q <- max(30L, nu %/% 4L)
  tn <- sprintf("y%d", seq_len(t))
  Sg <- cov_unif(seq(0.5, 1, length.out = t), 0.4)
  Sr <- cov_unif(rep(1, t), 0.2)
  gen <- factor(sample(sprintf("g%04d", 1:q), nu, TRUE), levels = sprintf("g%04d", 1:q))
  dl <- data.frame(unit = rep(seq_len(nu), t), trait = factor(rep(tn, each = nu), levels = tn), gen = rep(gen, t))
  X <- if (t > 1) stats::model.matrix(~ trait - 1, dl) else matrix(1, nrow(dl), 1)
  Zl <- lapply(tn, function(a) inc_niv(dl$gen) * (dl$trait == a))
  vrai <- c(vec_sigma("g", Sg), if (t > 1) vec_sigma("residuelle", Sr) else c(residuelle = Sr[1, 1]))
  Zall <- do.call(cbind, Zl); tr <- as.integer(dl$trait)
  V_de <- function(e) {
    R <- if (t > 1) mat_sigma(e, "residuelle", t) else matrix(e[["residuelle"]], 1, 1)
    .V_mt(Zall, mat_sigma(e, "g", t), q, R, tr, dl$unit, t)
  }
  dl$y <- simuler_y(V_de(vrai), X, seq_len(ncol(X)))
  res <- if (t > 1) rx_residual("us", trait = dl$trait, unit = dl$unit) else rx_residual()
  m <- rx_model(dl$y, X, list(if (t > 1) rx_term("g", Zl, struct = "us") else rx_term("g", Zl[[1]], t = 1L)),
                residual = res)
  w <- data.frame(gen = gen); for (a in tn) w[[a]] <- dl$y[dl$trait == a]
  list(modeles = list(dense = m, mme = m), vrai = vrai, V_de = V_de, y = dl$y, X = X,
       externes = list(asreml = function() {
         if (t == 1) {
           a <- asreml::asreml(y1 ~ 1, random = ~ gen, data = w, trace = FALSE)
           v <- asreml_vc(a)
           return(list(est = c(g = pick(v, "^gen$"), residuelle = pick(v, "^units!R$")),
                       logLik = a$loglik, conv = isTRUE(a$converge)))
         }
         f <- stats::as.formula(sprintf("cbind(%s) ~ trait", paste(tn, collapse = ", ")))
         a <- asreml::asreml(f, random = ~ us(trait):gen, residual = ~ id(units):us(trait),
                             data = w, trace = FALSE, maxit = 100)
         v <- asreml_vc(a); e <- c()
         for (i in seq_len(t)) for (j in seq_len(i)) {
           e[sprintf("g[%d,%d]", i, j)] <- pick(v, sprintf("^trait:gen!trait_%s:%s$", tn[i], tn[j]))
           e[sprintf("residuelle[%d,%d]", i, j)] <- pick(v, sprintf("^units:trait!trait_%s:%s$", tn[i], tn[j]))
         }
         list(est = e, logLik = a$loglik, conv = isTRUE(a$converge)) }))
}

fam_mtige <- function(n, graine, t, iee = FALSE) {
  set.seed(graine)
  np_ <- n %/% t
  cote <- 20L; nbac <- ceiling(np_ / cote^2)
  g <- expand.grid(col = seq_len(cote), row = seq_len(cote), bac = seq_len(nbac))[seq_len(np_), ]
  q <- max(30L, np_ %/% 5L)
  niv <- sprintf("g%04d", 1:q)
  gen <- factor(sample(niv, np_, TRUE), levels = niv)
  cle <- paste(g$bac, g$row, g$col)
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (dr in -1:1) for (dc in -1:1) {
    if (dr == 0 && dc == 0) next
    j <- match(paste(g$bac, g$row + dr, g$col + dc), cle); ok <- !is.na(j)
    ii <- c(ii, which(ok)); jj <- c(jj, j[ok]); xx <- c(xx, rep(1 / sqrt(dr^2 + dc^2), sum(ok)))
  }
  Wp <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(np_, np_))
  Zd_p <- inc_niv(gen); Zn_p <- Wp %*% Zd_p
  tn <- sprintf("y%d", seq_len(t))
  dl <- data.frame(unit = rep(seq_len(np_), t), trait = factor(rep(tn, each = np_), levels = tn))
  P <- Matrix::sparseMatrix(i = seq_len(nrow(dl)), j = dl$unit, x = 1, dims = c(nrow(dl), np_))
  Zd <- lapply(tn, function(a) (P %*% Zd_p) * (dl$trait == a))
  Zn <- lapply(tn, function(a) (P %*% Zn_p) * (dl$trait == a))
  # G sur [DGE_1..DGE_t, IGE_1..IGE_t] : correlations DGE 0,5, IGE 0,5,
  # DGE-IGE -0,2 sur le meme caractere, -0,1 entre caracteres
  Rg <- matrix(0, 2 * t, 2 * t)
  Rg[1:t, 1:t] <- 0.5; Rg[t + 1:t, t + 1:t] <- 0.5
  Rg[1:t, t + 1:t] <- -0.1; Rg[t + 1:t, 1:t] <- -0.1
  for (k in seq_len(t)) { Rg[k, t + k] <- Rg[t + k, k] <- -0.2 }
  diag(Rg) <- 1
  sdg <- c(rep(sqrt(0.5), t), rep(sqrt(0.1), t))
  G <- diag(sdg, 2 * t) %*% Rg %*% diag(sdg, 2 * t)
  if (min(eigen(G, only.values = TRUE)$values) <= 1e-8) stop("mtige : G simulee non definie positive")
  Sr <- cov_unif(rep(1, t), 0.3)
  Se <- cov_unif(rep(0.2, t), 0.3)                       # IEE : variance 0,2, correlation 0,3
  vrai <- c(vec_sigma("g", G), if (iee) vec_sigma("iee", Se),
            if (t > 1) vec_sigma("residuelle", Sr) else c(residuelle = 1))
  # IEE : incidence plante -> plantes voisines, une par caractere
  Ze <- if (iee) lapply(tn, function(a) (P %*% Wp) * (dl$trait == a)) else NULL
  Zeall <- if (iee) do.call(cbind, Ze) else NULL
  X <- if (t > 1) stats::model.matrix(~ trait - 1, dl) else matrix(1, nrow(dl), 1)
  Zall <- do.call(cbind, c(Zd, Zn)); tr <- as.integer(dl$trait)
  V_de <- function(e) {
    R <- if (t > 1) mat_sigma(e, "residuelle", t) else matrix(e[["residuelle"]], 1, 1)
    V <- .V_mt(Zall, mat_sigma(e, "g", 2 * t), q, R, tr, dl$unit, t)
    if (iee) V <- V + as.matrix(Zeall %*% Matrix::kronecker(methods::as(mat_sigma(e, "iee", t), "CsparseMatrix"),
                                                           Matrix::Diagonal(np_)) %*% Matrix::t(Zeall))
    V
  }
  y <- simuler_y(V_de(vrai), X, seq_len(ncol(X)))
  res <- if (t > 1) rx_residual("us", trait = dl$trait, unit = dl$unit) else rx_residual()
  termes <- list(rx_term("g", c(Zd, Zn), struct = "us"))
  if (iee) termes[[2]] <- if (t > 1) rx_term("iee", Ze, struct = "us") else rx_term("iee", Ze[[1]], t = 1L)
  m <- rx_model(y, X, termes, residual = res)
  ext <- list()
  ext$asreml <- function() {
    # asreml : incidences en colonnes denses (une par niveau et par bloc), donnees
    # triees par unite puis caractere pour la residuelle id(unit):us(car).
    # LA COLONNE DES CARACTERES NE S'APPELLE PAS `trait` : c'est le facteur
    # reserve qu'asreml cree lui-meme en format large (cbind). En format long,
    # la collision faussait son modele sans erreur (residuelle seule : Sigma
    # 0,818/0,330/1,548 au lieu de la covariance empirique 0,929/0,399/1,631 ;
    # arrets sous l'optimum et variances negatives sur mtige et mtigee).
    o <- order(dl$unit, dl$trait)
    blocs <- c(Zd, Zn)
    dd <- data.frame(y = y[o], car = dl$trait[o], unit = factor(dl$unit[o]))
    M <- do.call(cbind, lapply(blocs, function(B) as.matrix(B[o, , drop = FALSE])))
    if (iee) M <- cbind(M, do.call(cbind, lapply(Ze, function(B) as.matrix(B[o, , drop = FALSE]))))
    colnames(M) <- sprintf("c%05d", seq_len(ncol(M)))
    dd <- cbind(dd, M)
    grp <- lapply(seq_along(blocs), function(b) 3L + (b - 1L) * q + seq_len(q))
    names(grp) <- sprintf("b%02d", seq_along(blocs))
    if (iee) {
      ge <- lapply(seq_len(t), function(b) 3L + 2L * t * q + (b - 1L) * np_ + seq_len(np_))
      names(ge) <- sprintf("e%02d", seq_len(t))
      grp <- c(grp, ge)
    }
    list2env(list(grp = grp, dd = dd), globalenv())
    fe <- if (!iee) "" else if (t > 1)
      sprintf(" + str(~ %s, ~ us(%d):id(%d))", paste(sprintf("grp(e%02d)", seq_len(t)), collapse = " + "), t, np_)
    else " + grp(e01)"
    fr <- stats::as.formula(sprintf("~ str(~ %s, ~ us(%d):id(%d))%s",
                                    paste(sprintf("grp(b%02d)", seq_len(2L * t)), collapse = " + "), 2L * t, q, fe))
    a <- if (t > 1)
      asreml::asreml(y ~ car, random = fr, residual = ~ id(unit):us(car), group = grp, data = dd,
                     trace = FALSE, maxit = 100, workspace = "4gb")
    else asreml::asreml(y ~ 1, random = fr, group = grp, data = dd, trace = FALSE, maxit = 100, workspace = "4gb")
    v <- asreml_vc(a)
    us_ <- v[grep(sprintf("^grp\\(b01\\).*!us\\(%d\\)_", 2L * t), names(v))]
    e <- c()
    for (i in seq_len(2L * t)) for (j in seq_len(i)) e[sprintf("g[%d,%d]", i, j)] <- pick(us_, sprintf("_%d:%d$", i, j))
    if (iee) {
      ve <- v[grep("^grp\\(e01\\)", names(v))]
      if (t > 1) { for (i in seq_len(t)) for (j in seq_len(i)) e[sprintf("iee[%d,%d]", i, j)] <- pick(ve, sprintf("_%d:%d$", i, j)) }
      else e["iee"] <- pick(v, "^grp\\(e01\\)$")
    }
    if (t > 1) {
      for (i in seq_len(t)) for (j in seq_len(i))
        e[sprintf("residuelle[%d,%d]", i, j)] <- pick(v, sprintf("^unit:car!car_%s:%s$", tn[i], tn[j]))
    } else e["residuelle"] <- pick(v, "^units!R$")
    list(est = e, logLik = a$loglik, conv = isTRUE(a$converge))
  }
  list(modeles = list(dense = m, mme = m), vrai = vrai, V_de = V_de, y = y, X = X, externes = ext)
}

for (t_ in T_PARAMS) local({
  tt <- t_
  CONSTRUCTEURS[[paste0("mtus", tt)]] <<- function(n, graine) fam_mtus(n, graine, tt)
  CONSTRUCTEURS[[paste0("mtige", tt)]] <<- function(n, graine) fam_mtige(n, graine, tt)
  CONSTRUCTEURS[[paste0("mtigee", tt)]] <<- function(n, graine) fam_mtige(n, graine, tt, iee = TRUE)
})

# nombre de parametres de variance d'une famille du banc "parametres"
n_params_famille <- function(fam) {
  t <- as.integer(sub("^mt(us|igee|ige)", "", fam))
  if (startsWith(fam, "mtus")) return(t * (t + 1L))
  base <- (2L * t) * (2L * t + 1L) / 2L + t * (t + 1L) / 2L
  if (startsWith(fam, "mtigee")) base + t * (t + 1L) / 2L else base
}
