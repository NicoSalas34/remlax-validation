# ==============================================================================
# Catalogue of the simulated cases. Each builder takes a seed and returns:
#   name, family, note
#   models   named list; each entry holds `model` (rx_model), `V_of` (function
#            par -> V, written in R with the true geometry) and, for a lag,
#            `M_of` (par -> M)
#   engines  list of specifications (id, type, optimizer, model, expected_error)
#   truth    true components, named as in varcomp (fit$varcomp$name)
#   derive   NULL, or a function (par, type) -> data.frame(name, est, se) for
#            derived quantities (correlations of corg)
# The simulation goes through simulate_y(V_true): the arbiter's V, evaluated at
# the true parameters.
# ==============================================================================

m_fit <- function(opt, model = "main", expected_error = NULL)
  list(id = opt, type = "rx_fit", optimizer = opt, model = model, expected_error = expected_error)
# MME engine: rx_fit(engine = "mme"), automatic or forced sparse inversion
m_mme <- function(model = "main", inversion = "auto")
  list(id = if (inversion == "auto") "mme" else paste0("mme_", inversion), type = "mme",
       optimizer = "ai", inversion = inversion, model = model)
# MME engine by L-BFGS-B, sparse inversion. This route replaces the former RTMB
# engine (rx_fit_sparse), removed from remlax in October 2026, on its scope:
# kinship given by its precision (Kinv) and separable ar1 field.
m_mme_lbfgsb <- function(model = "main")
  list(id = "mme_lbfgsb", type = "mme", optimizer = "lbfgsb", inversion = "sparse", model = model)
M_STD <- list(m_fit("lbfgsb"), m_fit("ai"), m_mme(), m_mme(inversion = "sparse"))
M_STD_MME <- M_STD

case_def <- function(name, family, models, engines, truth, derive = NULL, note = "")
  list(name = name, family = family, models = models, engines = engines, truth = truth,
       derive = derive, note = note)

# True Sigma from its named components
S_from <- function(vals, name, t) {
  if (t == 1L) return(matrix(vals[[name]], 1, 1))
  S <- matrix(0, t, t)
  for (i in seq_len(t)) for (j in seq_len(i)) S[i, j] <- S[j, i] <- vals[[sprintf("%s[%d,%d]", name, i, j)]]
  S
}

# ------------------------------------------------------------------------------
# A. Simple terms, kinship
# ------------------------------------------------------------------------------
case_iid1 <- function(seed) {
  set.seed(seed)
  q <- 60; nk <- sample(6:10, q, replace = TRUE)
  g <- factor(rep(sprintf("g%02d", seq_len(q)), nk)); n <- length(g)
  x <- stats::rnorm(n); X <- cbind(1, x)
  Zg <- list(inc(g))
  V_of <- function(par) V_term(Zg, par$sig$g) + R_section(par$res[[1]], NULL, diag(n))
  truth <- c(g = 0.6, residual = 1.0)
  y <- simulate_y(V_of(list(sig = list(g = 0.6), res = list(1.0))), X, c(10, 0.5))
  m <- rx_model(y, X, list(rx_term("g", g)), rx_residual())
  case_def("iid_one_term", "simple terms", list(main = list(model = m, V_of = V_of)),
          M_STD_MME, truth, note = "one random factor, unbalanced group sizes")
}

case_iid2 <- function(seed) {
  set.seed(seed)
  d <- expand.grid(g = sprintf("g%02d", 1:50), env = sprintf("e%02d", 1:12))
  d <- d[sort(sample(nrow(d), round(0.92 * nrow(d)))), ]
  n <- nrow(d); X <- matrix(1, n, 1)
  Zg <- list(inc(d$g)); Ze <- list(inc(d$env))
  V_of <- function(par) V_term(Zg, par$sig$g) + V_term(Ze, par$sig$env) + par$res[[1]][1, 1] * diag(n)
  truth <- c(g = 0.5, env = 0.8, residual = 1.0)
  y <- simulate_y(V_of(list(sig = list(g = 0.5, env = 0.8), res = list(matrix(1)))), X, 5)
  m <- rx_model(y, X, list(rx_term("g", d$g), rx_term("env", d$env)), rx_residual())
  case_def("iid_two_crossed_terms", "simple terms", list(main = list(model = m, V_of = V_of)),
          M_STD_MME, truth, note = "genotype x environment crossed, 8 % of cells missing")
}

case_grm <- function(seed) {
  set.seed(seed)
  q <- 200; K <- simulated_grm(q, prefix = "G")
  g <- factor(rep(rownames(K), each = 3), levels = rownames(K))
  keep <- sort(sample(length(g), 570)); g <- droplevels(g[keep])
  n <- length(g); X <- matrix(1, n, 1)
  Kg <- K[levels(g), levels(g)]
  Zg <- list(inc(g))
  V_of <- function(par) V_term(Zg, par$sig$g, Kg) + par$res[[1]][1, 1] * diag(n)
  truth <- c(g = 0.5, residual = 0.5)
  y <- simulate_y(V_of(list(sig = list(g = 0.5), res = list(matrix(0.5)))), X, 2)
  m <- rx_model(y, X, list(rx_term("g", g, K = K)), rx_residual())
  case_def("kinship_grm_dense", "kinship", list(main = list(model = m, V_of = V_of)),
          M_STD_MME, truth, note = "VanRaden GRM (rx_grm, blending 0.05), 200 genotypes x 3 replicates")
}

case_pedigree <- function(seed) {
  set.seed(seed)
  # 100 founders, 50 pairs x 4 offspring, then 100 pairs between different
  # families x 2: no relationship between parents, so Henderson's A^-1 without
  # inbreeding is exact.
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
  V_of <- function(par) V_term(Za, par$sig$a, A) + par$res[[1]][1, 1] * diag(n)
  truth <- c(a = 0.5, residual = 0.5)
  y <- simulate_y(V_of(list(sig = list(a = 0.5), res = list(matrix(0.5)))), X, 3)
  md <- rx_model(y, X, list(rx_term("a", a, K = A)), rx_residual())
  mc <- rx_model(y, X, list(rx_term("a", a, Kinv = Ainv)), rx_residual())
  case_def("kinship_pedigree_kinv", "kinship",
          list(main = list(model = md, V_of = V_of), sparse = list(model = mc, V_of = V_of)),
          list(m_fit("lbfgsb"), m_fit("ai"), m_mme("sparse"), m_mme("sparse", "sparse"), m_mme_lbfgsb("sparse")), truth,
          note = paste("animal model, 500 individuals over 3 generations; A^-1 by Henderson",
                       "(R function, asreml::ainverse not used); Kinv alone: MME engine",
                       "(AI and L-BFGS-B); rx_fit on the same model with dense K = A"))
}

# ------------------------------------------------------------------------------
# B. Multi-trait: structures of Sigma
# ------------------------------------------------------------------------------
# Long format, t traits, n_g genotypes x 3 plants, ~5 % of rows missing.
mt_design <- function(t, n_g = 120, reps = 3, frac_missing = 0.05) {
  tn <- paste0("T", seq_len(t))
  pl <- data.frame(plant = sprintf("p%04d", seq_len(n_g * reps)),
                   g = rep(sprintf("g%03d", seq_len(n_g)), each = reps))
  d <- merge(pl, data.frame(trait = tn), by = NULL)
  d <- d[order(d$trait, d$plant), ]
  d <- d[sort(sample(nrow(d), round((1 - frac_missing) * nrow(d)))), ]
  d$trait <- factor(d$trait, levels = tn); d$g <- factor(d$g); d$plant <- factor(d$plant)
  n <- nrow(d)
  # complete levels for each trait, even if a genotype is missing from it
  Zl <- lapply(tn, function(k) { r <- which(d$trait == k)
    Matrix::sparseMatrix(i = r, j = as.integer(d$g[r]), x = 1, dims = c(n, nlevels(d$g)),
                         dimnames = list(NULL, levels(d$g))) })
  names(Zl) <- tn
  list(d = d, n = n, t = t, tn = tn, Zl = Zl, X = stats::model.matrix(~ trait - 1, d),
       tr = as.integer(d$trait), Cu = C_same_unit(as.integer(d$plant)))
}

case_mt <- function(seed, name, t, struct, rank, Sg, res_struct, Sr, res_rank = 0L,
                   engines = M_STD, note = "", n_g = 120) {
  set.seed(seed)
  D <- mt_design(t, n_g = n_g)
  V_of <- function(par) V_term(D$Zl, par$sig$g) + R_section(par$res[[1]], D$tr, D$Cu)
  y <- simulate_y(V_of(list(sig = list(g = Sg), res = list(Sr))), D$X, seq_len(t))
  m <- rx_model(y, D$X, list(rx_term("g", D$Zl, struct = struct, rank = rank)),
                rx_residual(res_struct, trait = D$d$trait, unit = D$d$plant, rank = res_rank))
  truth <- c((if (struct == "diag") diag_names else sigma_names)("g", Sg),
            (if (res_struct == "diag") diag_names else sigma_names)("residual", Sr))
  case_def(name, "multi-trait", list(main = list(model = m, V_of = V_of)), engines, truth,
          note = note)
}

fa_sigma <- function(L, psi) L %*% t(L) + diag(psi)
case_mt_diag <- function(seed) case_mt(seed, "mt_g_diag_r_diag", 3, "diag", 0L, diag(c(0.5, 0.8, 0.3)),
                                       "diag", diag(c(1, 0.7, 0.5)), engines = M_STD_MME,
                                       note = "3 traits, 120 genotypes x 3 plants, 5 % missing")
case_mt_us <- function(seed) {
  Sg <- matrix(c(0.6, 0.2, -0.15, 0.2, 0.5, 0.1, -0.15, 0.1, 0.4), 3)
  Sr <- matrix(c(1, 0.3, 0.2, 0.3, 0.8, -0.1, 0.2, -0.1, 0.6), 3)
  case_mt(seed, "mt_g_us_r_us", 3, "us", 0L, Sg, "us", Sr, engines = M_STD_MME,
         note = "genetic Sigma us, residual us between traits of the same plant")
}
case_mt_fa1 <- function(seed) case_mt(seed, "mt_g_fa1", 4, "fa", 1L,
                                      fa_sigma(matrix(c(0.7, 0.5, -0.4, 0.6)), c(0.2, 0.3, 0.25, 0.15)),
                                      "diag", diag(c(0.8, 0.6, 0.7, 0.5)), engines = M_STD_MME,
                                      note = "fa rank 1, 4 traits")
case_mt_fa2 <- function(seed) case_mt(seed, "mt_g_fa2", 5, "fa", 2L,
                                      fa_sigma(matrix(c(0.7, 0.5, -0.3, 0.6, 0.4, 0, 0.5, 0.6, -0.4, 0.3), 5),
                                               c(0.15, 0.2, 0.2, 0.15, 0.25)),
                                      "diag", diag(c(0.8, 0.6, 0.7, 0.5, 0.6)), n_g = 100, engines = M_STD_MME,
                                      note = "fa rank 2, 5 traits (rank 2 is not identified with 4 traits)")
case_mt_rr2 <- function(seed) {
  G <- matrix(c(0.7, 0.5, -0.3, 0.6, 0, 0.4, 0.5, -0.3), 4)
  case_mt(seed, "mt_g_rr2", 4, "rr", 2L, G %*% t(G), "diag", diag(c(0.8, 0.6, 0.7, 0.5)),
         note = "reduced rank 2 without specific variance: singular genetic Sigma")
}
case_mt_chol1 <- function(seed) {
  L <- diag(4); L[2, 1] <- 0.5; L[3, 2] <- -0.4; L[4, 3] <- 0.3
  case_mt(seed, "mt_g_chol1", 4, "chol", 1L, L %*% diag(c(0.6, 0.4, 0.5, 0.3)) %*% t(L),
         "diag", diag(c(0.8, 0.6, 0.7, 0.5)), note = "banded Cholesky of order 1")
}
case_mt_ante1 <- function(seed) {
  Ut <- diag(4); Ut[2, 1] <- 0.5; Ut[3, 2] <- -0.4; Ut[4, 3] <- 0.3
  Ai <- solve(Ut)
  case_mt(seed, "mt_g_ante1", 4, "ante", 1L, Ai %*% diag(1 / c(1.2, 1.5, 1.3, 1.6)^2) %*% t(Ai),
         "diag", diag(c(0.8, 0.6, 0.7, 0.5)), note = "antedependence of order 1")
}
case_mt_corh <- function(seed) {
  sd <- c(0.8, 0.6, 0.7)
  case_mt(seed, "mt_g_corh", 3, "corh", 0L, diag(sd) %*% C_cor(3, 0.4) %*% diag(sd),
         "diag", diag(c(0.8, 0.6, 0.7)), note = "heterogeneous variances, uniform correlation 0.4")
}
case_mt_resfa1 <- function(seed) case_mt(seed, "mt_residual_fa1", 4, "diag", 0L, diag(c(0.4, 0.5, 0.3, 0.6)),
                                         "fa", fa_sigma(matrix(c(0.6, 0.5, 0.4, -0.3)), c(0.5, 0.4, 0.6, 0.5)),
                                         res_rank = 1L, note = "residual fa rank 1 between traits of a plant")

# ------------------------------------------------------------------------------
# C. Structured residuals
# ------------------------------------------------------------------------------
parse_res <- function(...) getFromNamespace(".rx_parse_residual", "remlax")(...)

# Field 20 x 30, one genotype per cell (100 genotypes x 6), residual
# ar1(row) x ar1(column). `missing_frac`: fraction of cells removed.
case_res_ar1ar1 <- function(seed, missing_frac = 0, name = "residual_ar1xar1") {
  set.seed(seed)
  d <- expand.grid(col = 1:20, row = 1:30)
  d$g <- factor(sample(rep(sprintf("g%03d", 1:100), 6)))
  if (missing_frac > 0) d <- d[sort(sample(nrow(d), round((1 - missing_frac) * nrow(d)))), ]
  d$row <- factor(d$row); d$col <- factor(d$col)
  n <- nrow(d); X <- matrix(1, n, 1); Zg <- list(inc(d$g))
  rr <- as.integer(as.character(d$row)); cc <- as.integer(as.character(d$col))
  C_obs <- function(r1, r2) r1^abs(outer(rr, rr, "-")) * r2^abs(outer(cc, cc, "-"))
  V_of <- function(par) V_term(Zg, par$sig$g) +
    R_section(par$res[[1]], NULL, C_obs(par$lvl[["residual!ar1ar11"]], par$lvl[["residual!ar1ar12"]]))
  truth <- c(g = 0.4, residual = 1, "residual!ar1ar11" = 0.6, "residual!ar1ar12" = 0.4)
  y <- simulate_y(V_of(list(sig = list(g = 0.4), res = list(1),
                           lvl = truth[c("residual!ar1ar11", "residual!ar1ar12")])), X, 1)
  d$y <- y
  m <- rx_model(y, X, list(rx_term("g", d$g)), parse_res(~ ar1(row):ar1(col), d))
  case_def(name, "residuals", list(main = list(model = m, V_of = V_of)), M_STD, truth,
          note = if (missing_frac > 0) sprintf("%.0f %% of cells missing from the 30 x 20 field", 100 * missing_frac)
                 else "complete 30 x 20 field, residual = ~ ar1(row):ar1(col)")
}
case_res_ar1ar1_miss <- function(seed) case_res_ar1ar1(seed, missing_frac = 0.05, name = "residual_ar1xar1_missing")

case_dsum_ar1_iid <- function(seed) {
  set.seed(seed)
  a <- expand.grid(col = 1:15, row = 1:20); a$site <- "A"
  # site B within the row and column range of A: the row and col factors are
  # built on all the data, and B adds no level to them.
  b <- data.frame(col = rep(1:13, each = 20), row = rep(1:20, 13), site = "B")
  d <- rbind(a, b)
  d$g <- factor(sample(rep(sprintf("g%03d", 1:112), 5)))
  d$site <- factor(d$site); n <- nrow(d)
  X <- stats::model.matrix(~ site, d); Zg <- list(inc(d$g))
  ra <- which(d$site == "A"); rb <- which(d$site == "B")
  rr <- d$row[ra]; cc <- d$col[ra]
  V_of <- function(par) {
    V <- V_term(Zg, par$sig$g)
    V[ra, ra] <- V[ra, ra] + R_section(par$res[[1]], NULL,
      par$lvl[["A!ar1ar11"]]^abs(outer(rr, rr, "-")) * par$lvl[["A!ar1ar12"]]^abs(outer(cc, cc, "-")))
    V[rb, rb] <- V[rb, rb] + par$res[[2]][1, 1] * diag(length(rb))
    V
  }
  truth <- c(g = 0.5, A = 1.0, "A!ar1ar11" = 0.5, "A!ar1ar12" = 0.3, B = 0.6)
  y <- simulate_y(V_of(list(sig = list(g = 0.5), res = list(matrix(1), matrix(0.6)),
                           lvl = truth[c("A!ar1ar11", "A!ar1ar12")])), X, c(2, 1))
  d$row <- factor(d$row); d$col <- factor(d$col)
  m <- rx_model(y, X, list(rx_term("g", d$g)),
                parse_res(~ dsum(~ ar1(row):ar1(col) + units | site, levels = list("A", "B")), d))
  case_def("dsum_ar1xar1_and_iid", "residuals", list(main = list(model = m, V_of = V_of)), M_STD, truth,
          note = "two sites, dsum(~ ar1(row):ar1(col) + units | site): different structures")
}

case_dsum_iid_iid <- function(seed) {
  set.seed(seed)
  d <- data.frame(site = factor(rep(c("A", "B"), c(300, 260))))
  d$g <- factor(sample(rep(sprintf("g%03d", 1:80), 7)))
  n <- nrow(d); X <- stats::model.matrix(~ site, d); Zg <- list(inc(d$g))
  ra <- which(d$site == "A"); rb <- which(d$site == "B")
  V_of <- function(par) {
    V <- V_term(Zg, par$sig$g)
    V[ra, ra] <- V[ra, ra] + par$res[[1]][1, 1] * diag(length(ra))
    V[rb, rb] <- V[rb, rb] + par$res[[2]][1, 1] * diag(length(rb))
    V
  }
  truth <- c(g = 0.5, A = 1.0, B = 0.4)
  y <- simulate_y(V_of(list(sig = list(g = 0.5), res = list(matrix(1), matrix(0.4)))), X, c(2, 1))
  m <- rx_model(y, X, list(rx_term("g", d$g)), parse_res(~ dsum(~ units | site), d))
  case_def("dsum_iid_heterogeneous", "residuals", list(main = list(model = m, V_of = V_of)),
          M_STD_MME, truth, note = "iid residual specific to each site, dsum(~ units | site)")
}

# Residual "none": a plant term [I on T1, I on T2] us carries the residual.
# Twin: same model with residual us(trait):units. Both must reach the same
# logLik.
case_res_none <- function(seed) {
  set.seed(seed)
  t <- 2; np <- 300; tn <- c("T1", "T2")
  d <- data.frame(plant = factor(rep(sprintf("p%03d", 1:np), t)), trait = factor(rep(tn, each = np)))
  d$g <- factor(rep(sprintf("g%03d", rep(1:100, each = 3)), t))
  n <- nrow(d); X <- stats::model.matrix(~ trait - 1, d)
  Zg <- lapply(tn, function(k) { r <- which(d$trait == k)
    Matrix::sparseMatrix(i = r, j = as.integer(d$g[r]), x = 1, dims = c(n, nlevels(d$g))) })
  Zp <- lapply(tn, function(k) { r <- which(d$trait == k)
    Matrix::sparseMatrix(i = r, j = as.integer(d$plant[r]), x = 1, dims = c(n, np)) })
  Sg <- matrix(c(0.5, 0.15, 0.15, 0.4), 2); Sp <- matrix(c(1, 0.35, 0.35, 0.8), 2)
  V_of <- function(par) V_term(Zg, par$sig$g) + V_term(Zp, par$sig$plant)
  V_twin <- function(par) V_term(Zg, par$sig$g) +
    R_section(par$res[[1]], as.integer(d$trait), C_same_unit(as.integer(d$plant)))
  y <- simulate_y(V_of(list(sig = list(g = Sg, plant = Sp))), X, c(1, 2))
  m <- rx_model(y, X, list(rx_term("g", Zg, struct = "us"), rx_term("plant", Zp, struct = "us")),
                rx_residual("none"))
  ren <- c("residual[1,1]" = "plant[1,1]", "residual[2,1]" = "plant[2,1]",
           "residual[2,2]" = "plant[2,2]")
  mj <- rx_model(y, X, list(rx_term("g", Zg, struct = "us")),
                 rx_residual("us", trait = d$trait, unit = d$plant))
  case_def("residual_none_plant_term", "residuals",
          list(main = list(model = m, V_of = V_of), twin = list(model = mj, V_of = V_twin)),
          list(m_fit("lbfgsb"), m_fit("ai"), m_mme(),
               c(m_fit("lbfgsb", "twin"), list(rename = ren)), c(m_mme("twin"), list(rename = ren))),
          c(sigma_names("g", Sg), sigma_names("plant", Sp)),
          note = "residual none, the us plant term carries the residual; twin with residual us(trait):units")
}

# ------------------------------------------------------------------------------
# D. Structures between levels
# ------------------------------------------------------------------------------
case_lvl_ar1 <- function(seed) {
  set.seed(seed)
  q <- 100; plot <- factor(rep(sprintf("x%03d", 1:q), each = 4)); n <- length(plot)
  X <- matrix(1, n, 1); Zp <- list(inc(plot))
  V_of <- function(par) V_term(Zp, par$sig$plot, C_ar1(q, par$lvl[["plot!ar1"]])) +
    par$res[[1]][1, 1] * diag(n)
  truth <- c(plot = 1, "plot!ar1" = 0.7, residual = 0.5)
  y <- simulate_y(V_of(list(sig = list(plot = 1), res = list(matrix(0.5)), lvl = truth["plot!ar1"])), X, 0)
  m <- rx_model(y, X, list(rx_term("plot", plot, level = "ar1")), rx_residual())
  case_def("levels_ar1", "structures between levels", list(main = list(model = m, V_of = V_of)),
          M_STD_MME, truth, note = "transect of 100 plots x 4 subsamples, ar1 term")
}

case_lvl_ar1ar1 <- function(seed) {
  set.seed(seed)
  nr <- 20; nc <- 25
  d <- expand.grid(col = 1:nc, row = 1:nr)            # cell = (row-1)*nc + col
  d$g <- factor(sample(rep(sprintf("g%03d", 1:100), 5)))
  n <- nrow(d); X <- matrix(1, n, 1)
  cell <- factor((d$row - 1) * nc + d$col, levels = seq_len(nr * nc))
  Zs <- list(inc(cell)); Zg <- list(inc(d$g))
  V_of <- function(par) V_term(Zg, par$sig$g) +
    V_term(Zs, par$sig$field, C_ar1ar1(nr, nc, par$lvl[["field!ar1ar11"]], par$lvl[["field!ar1ar12"]])) +
    par$res[[1]][1, 1] * diag(n)
  truth <- c(g = 0.4, field = 0.8, "field!ar1ar11" = 0.6, "field!ar1ar12" = 0.4, residual = 0.3)
  y <- simulate_y(V_of(list(sig = list(g = 0.4, field = 0.8), res = list(matrix(0.3)),
                           lvl = truth[c("field!ar1ar11", "field!ar1ar12")])), X, 0)
  m <- rx_model(y, X, list(rx_term("g", d$g), rx_term("field", cell, level = "ar1ar1", dims = c(nr, nc))),
                rx_residual())
  case_def("levels_ar1xar1_field", "structures between levels",
          list(main = list(model = m, V_of = V_of)), c(M_STD_MME, list(m_mme_lbfgsb())), truth,
          note = "ar1 x ar1 random field on 20 x 25 cells, with nugget and genotype")
}

# Metric families: 150 positions x 3 observations, metric term `s`.
case_metric <- function(seed, kind, par_true, dim2 = TRUE, nu = 0.5) {
  set.seed(seed)
  q <- 150
  co <- if (dim2) cbind(stats::runif(q, 0, 15), stats::runif(q, 0, 15)) else cbind(sort(stats::runif(q, 0, 50)))
  pos <- factor(rep(sprintf("s%03d", 1:q), each = 3)); n <- length(pos)
  X <- matrix(1, n, 1); Zs <- list(inc(pos))
  np <- length(par_true)
  lvl_names <- if (np == 1L) paste0("s!", kind) else paste0("s!", kind, seq_len(np))
  V_of <- function(par) V_term(Zs, par$sig$s, C_metric(kind, co, unname(par$lvl[lvl_names]), nu = nu)) +
    par$res[[1]][1, 1] * diag(n)
  truth <- c(s = 1, stats::setNames(par_true, lvl_names), residual = 0.4)
  y <- simulate_y(V_of(list(sig = list(s = 1), res = list(matrix(0.4)), lvl = truth[lvl_names])), X, 0)
  opts <- if (kind == "mtrn") getFromNamespace(".rx_mtrn_opts", "remlax")(list(nu = sprintf("%g F", nu))) else NULL
  m <- rx_model(y, X, list(rx_term("s", pos, level = kind, coord = co, opts = opts)), rx_residual())
  case_def(paste0("metric_", kind), "structures between levels",
          list(main = list(model = m, V_of = V_of)), M_STD, truth,
          note = sprintf("%s, %s, 150 positions x 3 obs%s", kind, if (dim2) "2D" else "1D",
                         if (kind == "mtrn") sprintf(", nu fixed at %g", nu) else ""))
}
METRIC_CASES <- list(
  exp  = list(par = 0.85, d2 = FALSE), gau = list(par = 0.97, d2 = FALSE), lvr = list(par = 8, d2 = FALSE),
  iexp = list(par = 0.7, d2 = TRUE), igau = list(par = 0.9, d2 = TRUE), ieuc = list(par = 0.7, d2 = TRUE),
  sph  = list(par = 5, d2 = TRUE), cir = list(par = 5, d2 = TRUE),
  aexp = list(par = c(0.6, 0.8), d2 = TRUE), agau = list(par = c(0.85, 0.93), d2 = TRUE),
  mtrn = list(par = 1.5, d2 = TRUE, nu = 1.5))

# Genotype x environment: t = 60 genotypes (shared iid Sigma), q = 4
# environments correlated by a level structure (cor or corg).
ge_design <- function() {
  d <- expand.grid(rep = 1:2, env = sprintf("e%d", 1:4), g = sprintf("g%02d", 1:60))
  d$env <- factor(d$env); d$g <- factor(d$g)
  n <- nrow(d)
  Zl <- lapply(levels(d$g), function(k) { r <- which(d$g == k)
    Matrix::sparseMatrix(i = r, j = as.integer(d$env[r]), x = 1, dims = c(n, 4)) })
  list(d = d, n = n, Zl = Zl, X = stats::model.matrix(~ env, d))
}
case_lvl_cor <- function(seed) {
  set.seed(seed)
  D <- ge_design()
  V_of <- function(par) V_term(D$Zl, par$sig$ge[1, 1] * diag(60), C_cor(4, par$lvl[["ge!cor"]])) +
    par$res[[1]][1, 1] * diag(D$n)
  truth <- c("ge[1,1]" = 0.8, "ge!cor" = 0.5, residual = 0.5)
  y <- simulate_y(V_of(list(sig = list(ge = matrix(0.8)), res = list(matrix(0.5)), lvl = truth["ge!cor"])),
                 D$X, c(1, 0.5, -0.5, 0.2))
  m <- rx_model(y, D$X, list(suppressMessages(rx_term("ge", D$Zl, struct = "iid", level = "cor"))), rx_residual())
  case_def("levels_cor_ge", "structures between levels", list(main = list(model = m, V_of = V_of)),
          M_STD, truth, note = "60 genotypes x 4 environments x 2 reps, uniform correlation between environments")
}
case_lvl_corg <- function(seed) {
  set.seed(seed)
  D <- ge_design()
  Cv <- matrix(c(1, .6, .3, .1, .6, 1, .5, .2, .3, .5, 1, .4, .1, .2, .4, 1), 4)
  C_of <- function(th) C_corg_theta(th, 4)
  V_of <- function(par) V_term(D$Zl, par$sig$ge[1, 1] * diag(60), C_of(par$theta[2:7])) +
    par$res[[1]][1, 1] * diag(D$n)
  # simulation: true V directly (corg has no natural scale in varcomp)
  Vv <- V_term(D$Zl, 0.8 * diag(60), Cv) + 0.5 * diag(D$n)
  y <- simulate_y(Vv, D$X, c(1, 0.5, -0.5, 0.2))
  m <- rx_model(y, D$X, list(suppressMessages(rx_term("ge", D$Zl, struct = "iid", level = "corg", order = 4L))),
                rx_residual())
  ij <- which(lower.tri(Cv), arr.ind = TRUE)
  truth <- c("ge[1,1]" = 0.8, residual = 0.5,
            stats::setNames(Cv[ij], sprintf("C[%d,%d]", ij[, 1], ij[, 2])))
  derive <- function(par, type) {
    if (type != "rx_fit") return(NULL)
    f <- function(th) C_of(th[2:7])[ij]
    r <- delta_se(f, par)
    data.frame(name = sprintf("C[%d,%d]", ij[, 1], ij[, 2]), est = r$val, se = r$se)
  }
  case_def("levels_corg_ge", "structures between levels", list(main = list(model = m, V_of = V_of)),
          M_STD, truth, derive = derive,
          note = "general correlation between 4 environments; correlations by the delta method on theta")
}

# ------------------------------------------------------------------------------
# E. 2D splines (PS-ANOVA)
# ------------------------------------------------------------------------------
case_spl2d <- function(seed) {
  set.seed(seed)
  d <- expand.grid(col = 1:20, row = 1:30)
  d$g <- factor(sample(rep(sprintf("g%03d", 1:100), 6)))
  n <- nrow(d)
  sp <- rx_spl2d(d$col, d$row, nseg = c(6, 8))
  X <- cbind(1, sp$X)
  Zg <- list(inc(d$g))
  term_names <- vapply(sp$terms, `[[`, "", "name")
  Zsp <- lapply(sp$terms, function(tm) tm$Zl)
  names(Zsp) <- term_names
  scale <- vapply(Zsp, function(z) mean(Matrix::rowSums(z[[1]]^2)), 0)
  share <- c(0.3, 0.3, 0.15, 0.15, 0.3)
  s2 <- stats::setNames(share / scale, term_names)
  V_of <- function(par) {
    V <- V_term(Zg, par$sig$g) + par$res[[1]][1, 1] * diag(n)
    for (k in term_names) V <- V + V_term(Zsp[[k]], par$sig[[k]])
    V
  }
  sig_v <- c(list(g = 0.4), lapply(as.list(s2), function(v) matrix(v)))
  y <- simulate_y(V_of(list(sig = sig_v, res = list(matrix(0.5)))), X, c(5, rep(0.3, ncol(sp$X))))
  m <- rx_model(y, X, c(list(rx_term("g", d$g)), sp$terms), rx_residual())
  case_def("spline_spl2d", "splines", list(main = list(model = m, V_of = V_of)), M_STD_MME,
          c(g = 0.4, s2, residual = 0.5),
          note = "rx_spl2d(nseg = c(6, 8)): five PS-ANOVA terms, null part as fixed effects")
}

# ------------------------------------------------------------------------------
# F. Spatial lag (rx_lag)
# ------------------------------------------------------------------------------
case_lag <- function(seed) {
  set.seed(seed)
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
  M_of <- function(par) { Psi <- par$Psi; M <- matrix(0, n, n)
    for (a in 1:q) for (b in 1:q) M[(a - 1) * N + seq_len(N), (b - 1) * N + seq_len(N)] <- Psi[a, b] * Wd
    M }
  Cu <- C_same_unit(unit); tri <- as.integer(trait)
  V_of <- function(par) V_term(ZD, par$sig$gen) + R_section(par$res[[1]], tri, Cu)
  Sg <- matrix(c(1, 0.3, 0.3, 0.8), 2); Sr <- matrix(c(1, -0.2, -0.2, 0.9), 2)
  Psi <- matrix(c(0.06, 0.05, -0.04, 0.03), 2, dimnames = list(tn, tn))
  pv <- list(sig = list(gen = Sg), res = list(Sr), Psi = Psi)
  eta <- simulate_y(V_of(pv), X, c(3, 6))
  y <- as.numeric(solve(diag(n) - M_of(pv), eta))
  L <- rx_lag(W, trait, unit)
  m <- rx_model(y, X, list(rx_term("gen", ZD, struct = "us")),
                rx_residual("us", trait = trait, unit = unit), lag = L)
  truth <- c(sigma_names("gen", Sg), sigma_names("residual", Sr),
            "Psi[T1<-T1]" = 0.06, "Psi[T1<-T2]" = -0.04, "Psi[T2<-T1]" = 0.05, "Psi[T2<-T2]" = 0.03)
  case_def("lag_rx_lag", "spatial lag",
          list(main = list(model = m, V_of = V_of, M_of = M_of)),
          list(m_fit("lbfgsb"), m_fit("ai", expected_error = "lag"),
               modifyList(m_mme(), list(optimizer = "lbfgsb"))), truth,
          note = "y = (I - Psi (x) W)^-1 (Xb + Zu + e), queen W 15 x 20, 2 traits, reml_y criterion")
}

# ------------------------------------------------------------------------------
# G. DGE / IGE models
# ------------------------------------------------------------------------------
# Mixture of two species on 2 trays of 30 rows x 20 columns, alternate columns
# (A odd, B even); 150 genotypes per species x 4 plants.
case_ige_mixture <- function(seed, name = "ige_mixture_two_species", share_iee = 0.05, Rg = NULL) {
  set.seed(seed)
  tray <- expand.grid(col = 1:20, row = 1:30)
  d <- rbind(transform(tray, block = "b1"), transform(tray, block = "b2"))
  d$sp <- ifelse(d$col %% 2 == 1, "A", "B")
  d$id <- sprintf("pl%04d", seq_len(nrow(d)))
  d$geno <- NA_character_
  for (s in c("A", "B")) { r <- which(d$sp == s); d$geno[r] <- sample(rep(sprintf("%s%03d", s, 1:150), 4)) }
  n <- nrow(d); ia <- which(d$sp == "A"); ib <- which(d$sp == "B")
  K <- list(A = simulated_grm(150, prefix = "A"), B = simulated_grm(150, prefix = "B"))
  nb <- rx_neighbourhood(coord = cbind(d$row, d$col), group = d$sp, block = d$block, level = d$geno,
                         id = d$id, rank = list("A<-A" = 2, "A<-B" = 1, "B<-A" = 1, "B<-B" = 2),
                         kernel = "none", normalise = TRUE, output = "both")
  lev <- list(A = sort(unique(d$geno[ia])), B = sort(unique(d$geno[ib])))
  rows <- list(A = ia, B = ib); other <- c(A = "B", B = "A")
  Zgen <- lapply(c(A = "A", B = "B"), function(s) {
    r <- rows[[s]]; o <- other[[s]]
    DGE <- Matrix::sparseMatrix(i = r, j = match(d$geno[r], lev[[s]]), x = 1, dims = c(n, 150),
                                dimnames = list(NULL, lev[[s]]))
    list(DGE = DGE,
         IGE_intra = place(nb$level[[paste0(s, "<-", s)]][, lev[[s]]], r, n),
         IGE_inter = place(nb$level[[paste0(o, "<-", s)]][, lev[[s]]], rows[[o]], n))
  })
  # IEE: one iid term at plant level, all pairs pooled
  Zi <- Matrix::Matrix(0, n, n, sparse = TRUE)
  for (a in c("A", "B")) for (b in c("A", "B"))
    Zi[rows[[a]], rows[[b]]] <- nb$unit[[paste0(a, "<-", b)]]
  colnames(Zi) <- d$id
  # Exposure scaling: share of phenotypic variance targeted per component,
  # divided by its exposure functional (rx_exposure).
  Sig <- list()
  for (s in c("A", "B")) {
    r <- rows[[s]]; o <- other[[s]]
    ex_d <- rx_exposure(Zgen[[s]]$DGE[r, ], Zgen[[s]]$IGE_intra[r, ], K = K[[s]])
    ex_x <- rx_exposure(NULL, Zgen[[s]]$IGE_inter[rows[[o]], ], K = K[[s]])
    v <- c(0.30 / ex_d$d, 0.07 / ex_d$k, 0.05 / ex_x$k)
    R <- if (is.null(Rg)) matrix(c(1, -0.4, -0.2, -0.4, 1, 0.3, -0.2, 0.3, 1), 3) else Rg
    Sig[[paste0("g", s)]] <- diag(sqrt(v)) %*% R %*% diag(sqrt(v))
  }
  k_iee <- mean(Matrix::rowSums(Zi^2))
  s2_iee <- share_iee / k_iee
  Ks <- list(A = K$A[lev$A, lev$A], B = K$B[lev$B, lev$B])
  ZiZ <- as.matrix(Matrix::tcrossprod(Zi))
  tri <- as.integer(factor(d$sp))
  V_of <- function(par) V_term(Zgen$A, par$sig$gA, Ks$A) + V_term(Zgen$B, par$sig$gB, Ks$B) +
    par$sig$iee[1, 1] * ZiZ + R_section(par$res[[1]], tri, diag(n))
  Sr <- diag(c(0.5, 0.55))
  pv <- list(sig = list(gA = Sig$gA, gB = Sig$gB, iee = matrix(s2_iee)), res = list(Sr))
  X <- stats::model.matrix(~ sp + block, d)
  y <- simulate_y(V_of(pv), X, c(10, 4, 0.5))
  m <- rx_model(y, X, list(rx_term("gA", Zgen$A, K = K$A, struct = "us"),
                           rx_term("gB", Zgen$B, K = K$B, struct = "us"),
                           rx_term("iee", list(Zi), struct = "iid")),
                rx_residual("diag", trait = factor(d$sp)))
  truth <- c(sigma_names("gA", Sig$gA), sigma_names("gB", Sig$gB), iee = s2_iee, diag_names("residual", Sr))
  case_def(name, "DGE / IGE",
          list(main = list(model = m, V_of = V_of)), M_STD_MME, truth,
          note = paste("2 trays 30 x 20, alternate columns; DGE + IGE intra + IGE inter us (3 x 3) per",
                       "species with GRM; IEE iid; diag residual per species; ranks 2 (intra) and 1",
                       "(inter), L2 normalisation; h2 0.30, IGE 7 and 5 %, IEE 5 %, r(DGE, IGE) -0.4"))
}

case_ige_mono <- function(seed) {
  set.seed(seed)
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
  s2 <- c(dge = 0.3, ige = 0.08 / k, residual = 0.6)
  X <- matrix(1, n, 1)
  V_of <- function(par) V_term(list(Zd), par$sig$dge) + V_term(list(Zn), par$sig$ige) +
    par$res[[1]][1, 1] * diag(n)
  y <- simulate_y(V_of(list(sig = list(dge = s2[["dge"]], ige = s2[["ige"]]), res = list(matrix(s2[["residual"]])))),
                 X, 10)
  m <- rx_model(y, X, list(rx_term("dge", list(Zd)), rx_term("ige", list(Zn))), rx_residual())
  case_def("ige_univariate_one_species", "DGE / IGE", list(main = list(model = m, V_of = V_of)),
          M_STD_MME, s2, note = "monoculture 30 x 20, 120 genotypes x 5, DGE iid + IGE iid (rank 1, L2)")
}

# ------------------------------------------------------------------------------
CATALOGUE <- c(
  list(iid_one_term = case_iid1, iid_two_crossed_terms = case_iid2,
       kinship_grm_dense = case_grm, kinship_pedigree_kinv = case_pedigree,
       mt_g_diag_r_diag = case_mt_diag, mt_g_us_r_us = case_mt_us, mt_g_fa1 = case_mt_fa1,
       mt_g_fa2 = case_mt_fa2, mt_g_rr2 = case_mt_rr2, mt_g_chol1 = case_mt_chol1,
       mt_g_ante1 = case_mt_ante1, mt_g_corh = case_mt_corh, mt_residual_fa1 = case_mt_resfa1,
       residual_ar1xar1 = case_res_ar1ar1, residual_ar1xar1_missing = case_res_ar1ar1_miss,
       dsum_ar1xar1_and_iid = case_dsum_ar1_iid, dsum_iid_heterogeneous = case_dsum_iid_iid,
       residual_none_plant_term = case_res_none,
       levels_ar1 = case_lvl_ar1, levels_ar1xar1_field = case_lvl_ar1ar1,
       levels_cor_ge = case_lvl_cor, levels_corg_ge = case_lvl_corg),
  stats::setNames(lapply(names(METRIC_CASES), function(k) {
    s <- METRIC_CASES[[k]]
    function(seed) case_metric(seed, k, s$par, dim2 = s$d2, nu = s$nu %||% 0.5)
  }), paste0("metric_", names(METRIC_CASES))),
  list(spline_spl2d = case_spl2d, lag_rx_lag = case_lag,
       ige_mixture_two_species = case_ige_mixture, ige_univariate_one_species = case_ige_mono))
