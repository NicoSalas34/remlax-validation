# ==============================================================================
# Families of the "number of parameters" benchmark: multi-trait models whose
# number of variance parameters grows with the number t of traits, at a FIXED
# number of observations n. Same contract as familles.R.
#
#   mtus<t>   genotype us(t) + residual us(t) per unit        t(t+1) parameters
#   mtige<t>  DGE and IGE in one us(2t) block on the direct
#             and neighbourhood incidences (as in mvD7),
#             + residual us(t) per plant                      t(2t+1) + t(t+1)/2
#   mtigee<t> mtige + IEE: indirect environmental effect at
#             PLANT level (neighbourhood incidence plant ->
#             plant, us(t) between traits). As many levels as
#             plants: m of the order of n and C nearly full,
#             as in mvD7                                      mtige + t(t+1)/2
#
# Sourced after familles.R (which defines FAMILIES, seed_of, inc_levels, ...).
# ==============================================================================

T_PARAMS <- c(1L, 2L, 3L, 4L, 6L, 8L)
FAMILIES <- c(FAMILIES, paste0("mtus", T_PARAMS), paste0("mtige", T_PARAMS), paste0("mtigee", T_PARAMS))

# covariance matrix: variances v, uniform correlation r (positive definite
# for r > -1/(t-1))
cov_unif <- function(v, r) { t <- length(v); R <- matrix(r, t, t); diag(R) <- 1; D <- diag(sqrt(v), t); D %*% R %*% D }

# V = Z (G (x) I_q) Z' + R per unit, R[trait_i, trait_j] for two rows of the same unit
.V_mt <- function(Zall, G, q, R, tr, un, t) {
  Vg <- as.matrix(Zall %*% Matrix::kronecker(methods::as(G, "CsparseMatrix"), Matrix::Diagonal(q)) %*% Matrix::t(Zall))
  Rm <- Matrix::sparseMatrix(i = seq_along(tr), j = tr, x = 1, dims = c(length(tr), t))
  Zu <- Matrix::sparseMatrix(i = seq_along(un), j = un, x = 1)
  Vg + as.matrix(Matrix::tcrossprod(Zu)) * as.matrix(Rm %*% R %*% Matrix::t(Rm))
}

fam_mtus <- function(n, seed, t) {
  set.seed(seed)
  nu <- n %/% t; q <- max(30L, nu %/% 4L)
  tn <- sprintf("y%d", seq_len(t))
  Sg <- cov_unif(seq(0.5, 1, length.out = t), 0.4)
  Sr <- cov_unif(rep(1, t), 0.2)
  gen <- factor(sample(sprintf("g%04d", 1:q), nu, TRUE), levels = sprintf("g%04d", 1:q))
  dl <- data.frame(unit = rep(seq_len(nu), t), trait = factor(rep(tn, each = nu), levels = tn), gen = rep(gen, t))
  X <- if (t > 1) stats::model.matrix(~ trait - 1, dl) else matrix(1, nrow(dl), 1)
  Zl <- lapply(tn, function(a) inc_levels(dl$gen) * (dl$trait == a))
  truth <- c(vec_sigma("g", Sg), if (t > 1) vec_sigma("residual", Sr) else c(residual = Sr[1, 1]))
  Zall <- do.call(cbind, Zl); tr <- as.integer(dl$trait)
  V_of <- function(e) {
    R <- if (t > 1) mat_sigma(e, "residual", t) else matrix(e[["residual"]], 1, 1)
    .V_mt(Zall, mat_sigma(e, "g", t), q, R, tr, dl$unit, t)
  }
  dl$y <- simulate_y(V_of(truth), X, seq_len(ncol(X)))
  res <- if (t > 1) rx_residual("us", trait = dl$trait, unit = dl$unit) else rx_residual()
  m <- rx_model(dl$y, X, list(if (t > 1) rx_term("g", Zl, struct = "us") else rx_term("g", Zl[[1]], t = 1L)),
                residual = res)
  w <- data.frame(gen = gen); for (a in tn) w[[a]] <- dl$y[dl$trait == a]
  list(models = list(dense = m, mme = m), truth = truth, V_of = V_of, y = dl$y, X = X,
       external = list(asreml = function() {
         if (t == 1) {
           a <- asreml::asreml(y1 ~ 1, random = ~ gen, data = w, trace = FALSE)
           v <- asreml_vc(a)
           return(list(est = c(g = pick(v, "^gen$"), residual = pick(v, "^units!R$")),
                       logLik = a$loglik, conv = isTRUE(a$converge)))
         }
         f <- stats::as.formula(sprintf("cbind(%s) ~ trait", paste(tn, collapse = ", ")))
         a <- asreml::asreml(f, random = ~ us(trait):gen, residual = ~ id(units):us(trait),
                             data = w, trace = FALSE, maxit = 100)
         v <- asreml_vc(a); e <- c()
         for (i in seq_len(t)) for (j in seq_len(i)) {
           e[sprintf("g[%d,%d]", i, j)] <- pick(v, sprintf("^trait:gen!trait_%s:%s$", tn[i], tn[j]))
           e[sprintf("residual[%d,%d]", i, j)] <- pick(v, sprintf("^units:trait!trait_%s:%s$", tn[i], tn[j]))
         }
         list(est = e, logLik = a$loglik, conv = isTRUE(a$converge)) }))
}

fam_mtige <- function(n, seed, t, iee = FALSE) {
  set.seed(seed)
  np_ <- n %/% t
  side <- 20L; n_trays <- ceiling(np_ / side^2)
  g <- expand.grid(col = seq_len(side), row = seq_len(side), tray = seq_len(n_trays))[seq_len(np_), ]
  q <- max(30L, np_ %/% 5L)
  lev <- sprintf("g%04d", 1:q)
  gen <- factor(sample(lev, np_, TRUE), levels = lev)
  key <- paste(g$tray, g$row, g$col)
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (dr in -1:1) for (dc in -1:1) {
    if (dr == 0 && dc == 0) next
    j <- match(paste(g$tray, g$row + dr, g$col + dc), key); ok <- !is.na(j)
    ii <- c(ii, which(ok)); jj <- c(jj, j[ok]); xx <- c(xx, rep(1 / sqrt(dr^2 + dc^2), sum(ok)))
  }
  Wp <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(np_, np_))
  Zd_p <- inc_levels(gen); Zn_p <- Wp %*% Zd_p
  tn <- sprintf("y%d", seq_len(t))
  dl <- data.frame(unit = rep(seq_len(np_), t), trait = factor(rep(tn, each = np_), levels = tn))
  P <- Matrix::sparseMatrix(i = seq_len(nrow(dl)), j = dl$unit, x = 1, dims = c(nrow(dl), np_))
  Zd <- lapply(tn, function(a) (P %*% Zd_p) * (dl$trait == a))
  Zn <- lapply(tn, function(a) (P %*% Zn_p) * (dl$trait == a))
  # G on [DGE_1..DGE_t, IGE_1..IGE_t]: correlations DGE 0.5, IGE 0.5,
  # DGE-IGE -0.2 on the same trait, -0.1 between traits
  Rg <- matrix(0, 2 * t, 2 * t)
  Rg[1:t, 1:t] <- 0.5; Rg[t + 1:t, t + 1:t] <- 0.5
  Rg[1:t, t + 1:t] <- -0.1; Rg[t + 1:t, 1:t] <- -0.1
  for (k in seq_len(t)) { Rg[k, t + k] <- Rg[t + k, k] <- -0.2 }
  diag(Rg) <- 1
  sdg <- c(rep(sqrt(0.5), t), rep(sqrt(0.1), t))
  G <- diag(sdg, 2 * t) %*% Rg %*% diag(sdg, 2 * t)
  if (min(eigen(G, only.values = TRUE)$values) <= 1e-8) stop("mtige: simulated G is not positive definite")
  Sr <- cov_unif(rep(1, t), 0.3)
  Se <- cov_unif(rep(0.2, t), 0.3)                       # IEE: variance 0.2, correlation 0.3
  truth <- c(vec_sigma("g", G), if (iee) vec_sigma("iee", Se),
            if (t > 1) vec_sigma("residual", Sr) else c(residual = 1))
  # IEE: incidence plant -> neighbouring plants, one per trait
  Ze <- if (iee) lapply(tn, function(a) (P %*% Wp) * (dl$trait == a)) else NULL
  Zeall <- if (iee) do.call(cbind, Ze) else NULL
  X <- if (t > 1) stats::model.matrix(~ trait - 1, dl) else matrix(1, nrow(dl), 1)
  Zall <- do.call(cbind, c(Zd, Zn)); tr <- as.integer(dl$trait)
  V_of <- function(e) {
    R <- if (t > 1) mat_sigma(e, "residual", t) else matrix(e[["residual"]], 1, 1)
    V <- .V_mt(Zall, mat_sigma(e, "g", 2 * t), q, R, tr, dl$unit, t)
    if (iee) V <- V + as.matrix(Zeall %*% Matrix::kronecker(methods::as(mat_sigma(e, "iee", t), "CsparseMatrix"),
                                                           Matrix::Diagonal(np_)) %*% Matrix::t(Zeall))
    V
  }
  y <- simulate_y(V_of(truth), X, seq_len(ncol(X)))
  res <- if (t > 1) rx_residual("us", trait = dl$trait, unit = dl$unit) else rx_residual()
  terms_ <- list(rx_term("g", c(Zd, Zn), struct = "us"))
  if (iee) terms_[[2]] <- if (t > 1) rx_term("iee", Ze, struct = "us") else rx_term("iee", Ze[[1]], t = 1L)
  m <- rx_model(y, X, terms_, residual = res)
  ext <- list()
  ext$asreml <- function() {
    # asreml: incidences as dense columns (one per level and per block), data
    # sorted by unit then trait for the residual id(unit):us(car).
    # THE TRAIT COLUMN IS NOT CALLED `trait`. asreml reserves that factor and
    # creates it itself in wide format (cbind). In long format the clash
    # corrupted its model without any error (residual only: Sigma
    # 0.818/0.330/1.548 instead of the empirical covariance 0.929/0.399/1.631;
    # stops below the optimum and negative variances on mtige and mtigee).
    o <- order(dl$unit, dl$trait)
    blocks <- c(Zd, Zn)
    dd <- data.frame(y = y[o], car = dl$trait[o], unit = factor(dl$unit[o]))
    M <- do.call(cbind, lapply(blocks, function(B) as.matrix(B[o, , drop = FALSE])))
    if (iee) M <- cbind(M, do.call(cbind, lapply(Ze, function(B) as.matrix(B[o, , drop = FALSE]))))
    colnames(M) <- sprintf("c%05d", seq_len(ncol(M)))
    dd <- cbind(dd, M)
    grp <- lapply(seq_along(blocks), function(b) 3L + (b - 1L) * q + seq_len(q))
    names(grp) <- sprintf("b%02d", seq_along(blocks))
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
        e[sprintf("residual[%d,%d]", i, j)] <- pick(v, sprintf("^unit:car!car_%s:%s$", tn[i], tn[j]))
    } else e["residual"] <- pick(v, "^units!R$")
    list(est = e, logLik = a$loglik, conv = isTRUE(a$converge))
  }
  list(models = list(dense = m, mme = m), truth = truth, V_of = V_of, y = y, X = X, external = ext)
}

for (t_ in T_PARAMS) local({
  tt <- t_
  BUILDERS[[paste0("mtus", tt)]] <<- function(n, seed) fam_mtus(n, seed, tt)
  BUILDERS[[paste0("mtige", tt)]] <<- function(n, seed) fam_mtige(n, seed, tt)
  BUILDERS[[paste0("mtigee", tt)]] <<- function(n, seed) fam_mtige(n, seed, tt, iee = TRUE)
})

# number of variance parameters of a family of the "parameters" benchmark
n_params_family <- function(fam) {
  t <- as.integer(sub("^mt(us|igee|ige)", "", fam))
  if (startsWith(fam, "mtus")) return(t * (t + 1L))
  base <- (2L * t) * (2L * t + 1L) / 2L + t * (t + 1L) / 2L
  if (startsWith(fam, "mtigee")) base + t * (t + 1L) / 2L else base
}
