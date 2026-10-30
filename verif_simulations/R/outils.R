# ==============================================================================
# Outils de la verification par simulation : arbitre REML dense, structures de
# correlation ecrites en R, extraction des composantes, comparaisons.
#
# Rien ici n'appelle le moteur de remlax. L'arbitre recompose V a partir des
# matrices Sigma et des parametres de niveaux RENDUS par un ajustement, avec la
# geometrie VRAIE du dispositif (coordonnees, cellules, incidences construites
# par le cas), puis evalue la REML par la forme en V :
#     -2 logL = (n - p) log 2pi + log|V| + log|X'V^-1 X| + y'Py
# C'est la convention de fit$logLik (constante 2pi incluse) ; fit$logLik_asreml
# l'omet.
# ==============================================================================

LOG2PI <- log(2 * pi)

# --- REML dense ---------------------------------------------------------------
reml_dense <- function(V, y, X) {
  n <- length(y); p <- ncol(X)
  L <- chol(V)
  ldV <- 2 * sum(log(diag(L)))
  sol <- function(B) backsolve(L, forwardsolve(t(L), B))
  ViX <- sol(X); Viy <- sol(y)
  A <- crossprod(X, ViX)
  La <- chol(A)
  ldA <- 2 * sum(log(diag(La)))
  beta <- backsolve(La, forwardsolve(t(La), crossprod(X, Viy)))
  r <- Viy - ViX %*% beta
  neg2 <- (n - p) * LOG2PI + ldV + ldA + sum(y * r)
  list(neg2 = neg2, logLik = -0.5 * neg2)
}

# Critere reml_y du modele structurel (rx_lag) : REML de y* = A y, jacobien
# -2 log|det A|, et terme de moyenne -log|X_y'X_y| + log|X'X| avec X_y = A^-1 X.
reml_lag_dense <- function(V, y, X, M) {
  n <- length(y)
  A <- diag(n) - M
  ys <- as.numeric(A %*% y)
  base <- reml_dense(V, ys, X)$neg2
  ld <- function(B) as.numeric(determinant(B, logarithm = TRUE)$modulus)
  Xy <- solve(A, X)
  neg2 <- base - 2 * ld(A) - ld(crossprod(Xy)) + ld(crossprod(X))
  list(neg2 = neg2, logLik = -0.5 * neg2)
}

# --- correlations entre niveaux, parametres sur l'echelle naturelle -----------
C_ar1 <- function(q, rho) rho^abs(outer(seq_len(q), seq_len(q), "-"))
C_ar1ar1 <- function(nr, nc, rr, rc) kronecker(C_ar1(nr, rr), C_ar1(nc, rc))

.dist_axes <- function(co) {
  co <- as.matrix(co)
  dx <- abs(outer(co[, 1], co[, 1], "-"))
  dy <- if (ncol(co) > 1) abs(outer(co[, 2], co[, 2], "-")) else 0 * dx
  list(dx = dx, dy = dy, d = sqrt(dx^2 + dy^2))
}

matern_R <- function(h, nu) {
  out <- matrix(1, nrow(h), ncol(h))
  nz <- h > 0
  out[nz] <- 2^(1 - nu) / gamma(nu) * h[nz]^nu * besselK(h[nz], nu)
  out
}

# Familles metriques, avec le sens des parametres de levels.py : phi est une
# correlation elevee a une distance (exp, gau, iexp, igau, ieuc, aexp, agau),
# une portee (sph, cir, lvr) ou une echelle (mtrn).
C_metrique <- function(kind, co, par, nu = 0.5) {
  D <- .dist_axes(co)
  switch(kind,
    exp  = par[1]^D$d,
    gau  = par[1]^(D$d^2),
    lvr  = pmax(1 - D$dx / par[1], 0),
    iexp = par[1]^(D$dx + D$dy),
    igau = par[1]^(D$dx^2 + D$dy^2),
    ieuc = par[1]^D$d,
    aexp = par[1]^D$dx * par[2]^D$dy,
    agau = par[1]^(D$dx^2) * par[2]^(D$dy^2),
    sph  = { t <- pmin(D$d / par[1], 1); ifelse(D$d < par[1], 1 - 1.5 * t + 0.5 * t^3, 0) },
    cir  = { t <- pmin(D$d / par[1], 1)
             ifelse(D$d < par[1], 1 - (2 / pi) * (t * sqrt(1 - t^2) + asin(t)), 0) },
    mtrn = matern_R(D$d / par[1], nu),
    stop("famille metrique inconnue : ", kind))
}

# Correlation uniforme et correlation generale
C_cor <- function(q, r) (1 - r) * diag(q) + r * matrix(1, q, q)
# corg : L unitriangulaire a lignes normees, theta en triangle strict ligne par
# ligne (parametrisation de levels.py, reecrite ici pour la delta-methode).
C_corg_theta <- function(th, q) {
  L <- diag(q); k <- 1L
  for (i in seq_len(q)) for (j in seq_len(i - 1L)) { L[i, j] <- th[k]; k <- k + 1L }
  L <- L / sqrt(rowSums(L^2))
  tcrossprod(L)
}

# --- assemblage de V ------------------------------------------------------------
# Contribution d'un terme : Z (Sigma (x) K) Z', Z = [Z_1 ... Z_t] (caractere lent).
V_terme <- function(Zl, Sigma, K = NULL) {
  Z <- as.matrix(do.call(cbind, lapply(Zl, as.matrix)))
  q <- ncol(as.matrix(Zl[[1]]))
  if (is.null(K)) K <- diag(q)
  G <- kronecker(as.matrix(Sigma), as.matrix(K))
  Z %*% G %*% t(Z)
}

# Residuelle d'une section : R[i,j] = Sigma[tr_i, tr_j] * C[i,j], ou C est la
# correlation entre les UNITES des observations i et j (matrice n_s x n_s).
R_section <- function(Sigma, tr, C) {
  Sigma <- as.matrix(Sigma)
  if (length(Sigma) == 1L) return(Sigma[1, 1] * C)
  Sigma[tr, tr, drop = FALSE] * C
}
C_meme_unite <- function(unit) (outer(unit, unit, "==")) * 1

# --- ajustements : extraction normalisee ----------------------------------------
# Rend une liste `par` : sig (Sigma par terme), res (Sigma par section, dans
# l'ordre des sections), niv (parametres de niveaux sur l'echelle naturelle,
# nommes comme dans varcomp), Psi (decalage), est (composantes estimees nommees
# comme varcomp), se (erreurs-types correspondantes).
par_de_fit <- function(fit) {
  vc <- fit$varcomp
  est <- stats::setNames(vc$valeur, vc$nom)
  se  <- stats::setNames(vc$se, vc$nom)
  list(sig = fit$sigmas, res = unname(fit$sigmas_res), res_noms = names(fit$sigmas_res),
       niv = est[grepl("!", names(est), fixed = TRUE)],
       Psi = if (!is.null(fit$lag)) fit$lag$Psi else NULL,
       est = est, se = se, theta = fit$theta, hessian = fit$hessian,
       par_floor = fit$par_floor, par_ceil = fit$par_ceil)
}

# Noms de composantes a la maniere de component_names() (components.py).
noms_sigma <- function(nom, S) {
  S <- as.matrix(S); t <- nrow(S)
  if (t == 1L) return(stats::setNames(S[1, 1], nom))
  v <- c(); nm <- c()
  for (i in seq_len(t)) for (j in seq_len(i)) { v <- c(v, S[i, j]); nm <- c(nm, sprintf("%s[%d,%d]", nom, i, j)) }
  stats::setNames(v, nm)
}

# Diagonale seule (structure diag : les covariances sont nulles par construction)
noms_diag <- function(nom, S) {
  v <- noms_sigma(nom, S)
  if (nrow(as.matrix(S)) == 1L) return(v)
  v[grepl("\\[(\\d+),\\1\\]$", names(v))]
}

noms_residuelle <- function(model) {
  secs <- model$residual$sections
  if (identical(model$residual$struct, "none")) return(character(0))
  if (is.null(secs) || length(secs) <= 1L) return("residuelle")
  vapply(secs, function(s) s$name %||% "section", "")
}

par_de_mme <- function(mm, model) {
  rn <- noms_residuelle(model)
  res <- unname(mm$sigma_res)
  est <- c(unlist(lapply(names(mm$sigmas), function(nm) noms_sigma(nm, mm$sigmas[[nm]]))),
           unlist(unname(Map(function(nm, S) noms_sigma(nm, S), rn, res))))
  # parametres de niveaux (rho d'un terme ar1 ou ar1 x ar1), nommes comme pour
  # les autres moteurs : <terme>!ar1, <terme>!ar1ar11, <terme>!ar1ar12
  niv <- c()
  for (nm in names(mm$rho %||% list())) {
    v <- mm$rho[[nm]]; if (!length(v)) next
    lv <- model$terms[[nm]]$level
    niv <- c(niv, stats::setNames(v, if (length(v) == 1L) paste0(nm, "!", lv) else paste0(nm, "!", lv, seq_along(v))))
  }
  est <- c(est, niv)
  list(sig = mm$sigmas, res = res, res_noms = rn, niv = niv, Psi = NULL,
       est = est, se = stats::setNames(rep(NA_real_, length(est)), names(est)))
}

# --- delta-methode numerique sur theta (composantes derivees) -----------------
cov_theta <- function(par) {
  H <- par$hessian; th <- par$theta
  if (is.null(H)) return(NULL)
  libre <- !(th <= par$par_floor + 1e-7 | th >= par$par_ceil - 1e-7)
  Hf <- (H[libre, libre, drop = FALSE] + t(H[libre, libre, drop = FALSE])) / 2
  Vi <- tryCatch(solve(Hf), error = function(e) NULL)
  if (is.null(Vi)) return(NULL)
  Vt <- matrix(0, length(th), length(th)); Vt[libre, libre] <- 2 * Vi
  Vt
}
delta_se <- function(f, par, h = 1e-6) {
  Vt <- cov_theta(par); th <- par$theta
  f0 <- f(th)
  if (is.null(Vt)) return(list(val = f0, se = rep(NA_real_, length(f0))))
  J <- sapply(seq_along(th), function(k) { e <- rep(0, length(th)); e[k] <- h
    (f(th + e) - f(th - e)) / (2 * h) })
  J <- matrix(J, nrow = length(f0))
  list(val = f0, se = sqrt(pmax(diag(J %*% Vt %*% t(J)), 0)))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- simulation -----------------------------------------------------------------
# y = X beta + L z, L L' = V_vrai : une seule voie de simulation, la meme forme
# de V que l'arbitre, avec les parametres vrais.
simuler_y <- function(V, X, beta) {
  L <- t(chol(V))
  as.numeric(X %*% beta + L %*% stats::rnorm(nrow(V)))
}

# GRM a partir de marqueurs simules en familles (demi-freres), avec blending
# pour la rendre definie positive (une GRM centree est de rang q - 1).
grm_simulee <- function(q, m = 800, n_fam = 15, blend = 0.05, prefixe = "G") {
  p <- stats::runif(m, 0.1, 0.9)
  fam <- sample(rep(seq_len(n_fam), length.out = q))
  meres <- t(sapply(seq_len(n_fam), function(f) stats::rbinom(m, 2, p)))
  M <- t(sapply(seq_len(q), function(i) {
    h <- meres[fam[i], ] / 2
    stats::rbinom(m, 1, h) + stats::rbinom(m, 1, p)
  }))
  rownames(M) <- sprintf("%s%03d", prefixe, seq_len(q))
  rx_grm(M, ploidy = 2, blend = blend, coding = "count")
}

# A^-1 par les regles de Henderson, sans consanguinite (parents inconnus = 0).
ainv_henderson <- function(sire, dam) {
  n <- length(sire)
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  aj <- function(i, j, x) { ii <<- c(ii, i); jj <<- c(jj, j); xx <<- c(xx, x) }
  for (k in seq_len(n)) {
    s <- sire[k]; d <- dam[k]
    ps <- c(s, d)[c(s, d) > 0]
    delta <- switch(length(ps) + 1L, 1, 4 / 3, 2)
    aj(k, k, delta)
    for (a in ps) { aj(k, a, -delta / 2); aj(a, k, -delta / 2) }
    for (a in ps) for (b in ps) aj(a, b, delta / 4)
  }
  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))
}

# Incidence indicatrice n x q d'un facteur (niveaux dans l'ordre du facteur)
inc <- function(f, n = length(f), rows = seq_along(f)) {
  f <- factor(f)
  Matrix::sparseMatrix(i = rows, j = as.integer(f), x = 1, dims = c(n, nlevels(f)),
                       dimnames = list(NULL, levels(f)))
}
# Place une matrice m x q (lignes `rows` de n) dans une matrice n x q
placer <- function(M, rows, n) {
  M <- methods::as(M, "CsparseMatrix")
  T <- methods::as(M, "TsparseMatrix")
  Matrix::sparseMatrix(i = rows[T@i + 1L], j = T@j + 1L, x = T@x, dims = c(n, ncol(M)),
                       dimnames = list(NULL, colnames(M)))
}
