# ==============================================================================
# Families of the engine benchmark. Each builder takes (n, seed) and returns:
#   models    rx_model for remlax: `dense` (dense routes) and `mme` (MME routes,
#             often the same object; pedigree: K for dense, Kinv for MME)
#   truth     true parameters, under NAMES SHARED by all software
#   V_of(est) V built from a vector of named parameters (arbiter)
#   y, X      response and fixed effects, for the arbiter
#   external  named list of functions without arguments -> list(est, logLik, conv)
#             for asreml, sommer and lme4 when they apply
# Shared names: simple component "gen", matrix "g[i,j]" (i >= j),
# residual "residual" or "residual[i,j]", level parameters
# "<term>!<name>".
# ==============================================================================

SIZES <- c(500L, 1000L, 1500L, 2000L, 3000L, 4000L, 6000L, 8000L)
FAMILIES <- c("iid", "grm", "us", "ped", "ige", "field", "resid", "exp", "fa")
seed_of <- function(fam, n, rep) as.integer(20261022L + 100000L * match(fam, FAMILIES) + 10L * n + rep)

# Incidence that KEEPS every level of the factor. inc() of the suite
# refactors the factor and drops levels never drawn, which misaligns Z and K.
inc_levels <- function(f) Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
                                            dims = c(length(f), nlevels(f)), dimnames = list(NULL, levels(f)))

vec_sigma <- function(name, S) {
  S <- as.matrix(S); t <- nrow(S)
  if (t == 1L) return(stats::setNames(S[1, 1], name))
  out <- c()
  for (i in seq_len(t)) for (j in seq_len(i)) out[sprintf("%s[%d,%d]", name, i, j)] <- S[i, j]
  out
}
mat_sigma <- function(est, name, t) {
  if (t == 1L) return(matrix(est[[name]], 1, 1))
  S <- matrix(0, t, t)
  for (i in seq_len(t)) for (j in seq_len(i)) S[i, j] <- S[j, i] <- est[[sprintf("%s[%d,%d]", name, i, j)]]
  S
}

# asreml variance components: table (row name -> component)
asreml_vc <- function(a) { vc <- summary(a)$varcomp; stats::setNames(vc$component, rownames(vc)) }
pick <- function(v, pattern) {
  i <- grep(pattern, names(v), perl = TRUE)
  if (length(i) != 1L) stop("component missing or ambiguous: ", pattern, " among ", paste(names(v), collapse = ", "))
  unname(v[i])
}
sommer_vc <- function(m) { vc <- summary(m)$varcomp; stats::setNames(vc$VarComp, rownames(vc)) }

# lme4 with a kinship K = L L': Z u = (Z L) v, v iid (the marginal model is the same)
lmer_kinship <- function(formula, data, factor_name, K) {
  lf <- lme4::lFormula(formula, data = data, REML = TRUE)
  idx <- which(names(lf$reTrms$cnms) == factor_name)
  start <- c(0, cumsum(vapply(lf$reTrms$Ztlist, nrow, 0)))
  rows <- (start[idx] + 1):start[idx + 1]
  lev <- rownames(lf$reTrms$Zt)[rows]          # lme4 drops absent levels
  L <- t(chol(K[lev, lev]))
  Zt <- lf$reTrms$Zt
  Zt[rows, ] <- methods::as(Matrix::t(L) %*% Zt[rows, ], "CsparseMatrix")
  lf$reTrms$Zt <- methods::as(Zt, "CsparseMatrix")
  dev <- do.call(lme4::mkLmerDevfun, lf)
  opt <- lme4::optimizeLmer(dev)
  lme4::mkMerMod(environment(dev), opt, lf$reTrms, fr = lf$fr)
}
lme4_vc <- function(f) {
  vc <- as.data.frame(lme4::VarCorr(f))
  stats::setNames(vc$vcov, ifelse(vc$grp == "Residual", "residual", vc$grp))
}

# --------------------------------------------------------------------------- iid
fam_iid <- function(n, seed) {
  set.seed(seed)
  q <- max(20L, n %/% 5L)
  d <- data.frame(gen = factor(sample(sprintf("g%04d", 1:q), n, TRUE), levels = sprintf("g%04d", 1:q)),
                  blk = factor(sample(sprintf("b%02d", 1:25), n, TRUE), levels = sprintf("b%02d", 1:25)),
                  x = stats::rnorm(n))
  X <- cbind(1, d$x); Zg <- inc(d$gen); Zb <- inc(d$blk)
  truth <- c(gen = 0.5, blk = 0.3, residual = 1.0)
  V_of <- function(e) as.matrix(e[["gen"]] * Matrix::tcrossprod(Zg) + e[["blk"]] * Matrix::tcrossprod(Zb)) +
    e[["residual"]] * diag(n)
  d$y <- simulate_y(V_of(truth), X, c(2, 0.5))
  m <- rx_model(d$y, X, list(rx_term("gen", d$gen), rx_term("blk", d$blk)))
  list(models = list(dense = m, mme = m), truth = truth, V_of = V_of, y = d$y, X = X,
       external = list(
         asreml = function() { a <- asreml::asreml(y ~ x, random = ~ gen + blk, data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(gen = pick(v, "^gen$"), blk = pick(v, "^blk$"), residual = pick(v, "^units!R$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) },
         sommer = function() { s <- sommer::mmer(y ~ x, random = ~ gen + blk, rcov = ~ units, data = d, verbose = FALSE)
           v <- sommer_vc(s)
           list(est = c(gen = pick(v, "^gen"), blk = pick(v, "^blk"), residual = pick(v, "^units")),
                logLik = as.numeric(utils::tail(s$monitor[1, ], 1)), conv = isTRUE(s$convergence)) },
         lme4 = function() { f <- lme4::lmer(y ~ x + (1 | gen) + (1 | blk), data = d, REML = TRUE)
           v <- lme4_vc(f)
           list(est = c(gen = v[["gen"]], blk = v[["blk"]], residual = v[["residual"]]),
                logLik = as.numeric(stats::logLik(f)), conv = length(f@optinfo$conv$lme4$messages) == 0) }))
}

# --------------------------------------------------------------------------- grm
fam_grm <- function(n, seed) {
  set.seed(seed)
  q <- max(50L, n %/% 4L)
  K <- simulated_grm(q, 800, 15, 0.05, "g")
  lev <- rownames(K)
  d <- data.frame(gen = factor(sample(lev, n, TRUE), levels = lev),
                  blk = factor(sample(sprintf("b%02d", 1:25), n, TRUE), levels = sprintf("b%02d", 1:25)),
                  x = stats::rnorm(n))
  X <- cbind(1, d$x); Zg <- inc_levels(d$gen); Zb <- inc(d$blk)
  truth <- c(gen = 0.5, blk = 0.3, residual = 1.0)
  V_of <- function(e) e[["gen"]] * as.matrix(Zg %*% K %*% Matrix::t(Zg)) +
    as.matrix(e[["blk"]] * Matrix::tcrossprod(Zb)) + e[["residual"]] * diag(n)
  d$y <- simulate_y(V_of(truth), X, c(2, 0.5))
  Zl <- list(Zg); colnames(Zl[[1]]) <- lev
  m <- rx_model(d$y, X, list(rx_term("gen", Zl, K = K), rx_term("blk", d$blk)))
  list(models = list(dense = m, mme = m), truth = truth, V_of = V_of, y = d$y, X = X,
       external = list(
         asreml = function() { a <- asreml::asreml(y ~ x, random = ~ vm(gen, K) + blk, data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(gen = pick(v, "^vm\\(gen, K\\)$"), blk = pick(v, "^blk$"), residual = pick(v, "^units!R$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) },
         sommer = function() { s <- sommer::mmer(y ~ x, random = ~ sommer::vsr(gen, Gu = K) + blk, rcov = ~ units,
                                                 data = d, verbose = FALSE)
           v <- sommer_vc(s)
           list(est = c(gen = pick(v, "^u:gen"), blk = pick(v, "^blk"), residual = pick(v, "^units")),
                logLik = as.numeric(utils::tail(s$monitor[1, ], 1)), conv = isTRUE(s$convergence)) },
         lme4 = function() { f <- lmer_kinship(y ~ x + (1 | gen) + (1 | blk), d, "gen", K)
           v <- lme4_vc(f)
           list(est = c(gen = v[["gen"]], blk = v[["blk"]], residual = v[["residual"]]),
                logLik = as.numeric(stats::logLik(f)), conv = TRUE) }))
}

# ---------------------------------------------------------------------------- us
fam_us <- function(n, seed) {
  set.seed(seed)
  nu <- n %/% 3L; q <- max(30L, nu %/% 4L)
  Sg <- matrix(c(1, .5, .2, .5, .8, .3, .2, .3, .6), 3)
  Sr <- matrix(c(1, .3, .1, .3, .9, .2, .1, .2, .7), 3)
  gen <- factor(sample(sprintf("g%04d", 1:q), nu, TRUE), levels = sprintf("g%04d", 1:q))
  dl <- data.frame(unit = rep(seq_len(nu), 3), trait = factor(rep(c("y1", "y2", "y3"), each = nu)),
                   gen = rep(gen, 3))
  X <- stats::model.matrix(~ trait - 1, dl)
  Zl <- lapply(c("y1", "y2", "y3"), function(a) inc_levels(dl$gen) * (dl$trait == a))
  truth <- c(vec_sigma("g", Sg), vec_sigma("residual", Sr))
  Zall <- do.call(cbind, Zl); tr <- as.integer(dl$trait); Zu <- inc(factor(dl$unit))
  V_of <- function(e) {
    G <- mat_sigma(e, "g", 3); R <- mat_sigma(e, "residual", 3)
    # V = Z (G (x) I) Z' + R per unit: R[tr_i, tr_j] for two rows of the same unit
    Vg <- as.matrix(Zall %*% Matrix::kronecker(methods::as(G, "CsparseMatrix"), Matrix::Diagonal(q)) %*% Matrix::t(Zall))
    Rm <- Matrix::sparseMatrix(i = seq_along(tr), j = tr, x = 1, dims = c(length(tr), 3))
    Vg + as.matrix(Matrix::tcrossprod(Zu)) * as.matrix(Rm %*% R %*% Matrix::t(Rm))
  }
  dl$y <- simulate_y(V_of(truth), X, c(1, 2, 3))
  w <- data.frame(gen = gen, y1 = dl$y[dl$trait == "y1"], y2 = dl$y[dl$trait == "y2"], y3 = dl$y[dl$trait == "y3"])
  m <- rx_model(dl$y, X, list(rx_term("g", Zl, struct = "us")),
                residual = rx_residual("us", trait = dl$trait, unit = dl$unit))
  list(models = list(dense = m, mme = m), truth = truth, V_of = V_of, y = dl$y, X = X,
       external = list(
         asreml = function() { a <- asreml::asreml(cbind(y1, y2, y3) ~ trait, random = ~ us(trait):gen,
                                                   residual = ~ id(units):us(trait), data = w, trace = FALSE)
           v <- asreml_vc(a); e <- c()
           tn <- c("y1", "y2", "y3")
           for (i in 1:3) for (j in 1:i) {
             e[sprintf("g[%d,%d]", i, j)] <- pick(v, sprintf("^trait:gen!trait_%s:%s$", tn[i], tn[j]))
             e[sprintf("residual[%d,%d]", i, j)] <- pick(v, sprintf("^units:trait!trait_%s:%s$", tn[i], tn[j]))
           }
           list(est = e, logLik = a$loglik, conv = isTRUE(a$converge)) },
         sommer = function() { s <- sommer::mmer(cbind(y1, y2, y3) ~ 1, random = ~ sommer::vsr(gen, Gtc = sommer::unsm(3)),
                                                 rcov = ~ sommer::vsr(units, Gtc = sommer::unsm(3)), data = w, verbose = FALSE)
           e <- c(vec_sigma("g", s$sigma[[1]]), vec_sigma("residual", s$sigma[[2]]))
           list(est = e, logLik = as.numeric(utils::tail(s$monitor[1, ], 1)), conv = isTRUE(s$convergence)) }))
}

# --------------------------------------------------------------------------- ped
fam_ped <- function(n, seed) {
  set.seed(seed)
  f0 <- max(50L, n %/% 10L); ngen <- 4L
  per_gen <- ceiling((n - f0) / ngen)
  ped <- data.frame(id = seq_len(f0), sire = 0L, dam = 0L)
  for (g in seq_len(ngen)) {
    prev <- ped$id; k <- min(per_gen, n - nrow(ped)); if (k <= 0) break
    ids <- max(ped$id) + seq_len(k)
    ped <- rbind(ped, data.frame(id = ids, sire = sample(prev, k, TRUE), dam = sample(prev, k, TRUE)))
  }
  ped <- ped[seq_len(n), ]
  ai <- asreml::ainverse(ped)
  lab <- as.character(ped$id)
  Ainv <- Matrix::sparseMatrix(i = ai[, 1], j = ai[, 2], x = ai[, 3], dims = c(n, n), symmetric = TRUE)
  Ainv <- methods::as(methods::as(Ainv, "CsparseMatrix"), "generalMatrix"); dimnames(Ainv) <- list(lab, lab)
  ldA <- attr(ai, "logdet")
  A <- as.matrix(solve(Matrix::forceSymmetric(Ainv))); dimnames(A) <- list(lab, lab)
  d <- data.frame(animal = factor(lab, levels = lab),
                  blk = factor(sample(sprintf("b%02d", 1:20), n, TRUE), levels = sprintf("b%02d", 1:20)),
                  x = stats::rnorm(n))
  X <- cbind(1, d$x); Zb <- inc(d$blk)
  truth <- c(animal = 0.4, blk = 0.2, residual = 0.6)
  V_of <- function(e) e[["animal"]] * A + e[["blk"]] * as.matrix(Matrix::tcrossprod(Zb)) + e[["residual"]] * diag(n)
  d$y <- simulate_y(V_of(truth), X, c(3, 0.5))
  Za <- list(Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1, dims = c(n, n), dimnames = list(NULL, lab)))
  m_dense <- rx_model(d$y, X, list(rx_term("animal", Za, K = A), rx_term("blk", d$blk)))
  m_mme <- rx_model(d$y, X, list(rx_term("animal", Za, Kinv = Ainv, Kinv_logdet = ldA), rx_term("blk", d$blk)))
  list(models = list(dense = m_dense, mme = m_mme), truth = truth, V_of = V_of, y = d$y, X = X,
       external = list(
         asreml = function() { a <- asreml::asreml(y ~ x, random = ~ vm(animal, ai) + blk, data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(animal = pick(v, "^vm\\(animal, ai\\)$"), blk = pick(v, "^blk$"), residual = pick(v, "^units!R$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) },
         sommer = function() { s <- sommer::mmer(y ~ x, random = ~ sommer::vsr(animal, Gu = A) + blk, rcov = ~ units,
                                                 data = d, verbose = FALSE)
           v <- sommer_vc(s)
           list(est = c(animal = pick(v, "^u:animal"), blk = pick(v, "^blk"), residual = pick(v, "^units")),
                logLik = as.numeric(utils::tail(s$monitor[1, ], 1)), conv = isTRUE(s$convergence)) },
         lme4 = function() {
           suppressPackageStartupMessages(library(lme4))   # pedigreemm calls lmer without a prefix
           pd <- pedigreemm::pedigree(sire = ifelse(ped$sire == 0, NA, ped$sire), dam = ifelse(ped$dam == 0, NA, ped$dam), label = ped$id)
           f <- pedigreemm::pedigreemm(y ~ x + (1 | animal) + (1 | blk), pedigree = list(animal = pd), data = d, REML = TRUE,
             control = lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
           v <- lme4_vc(f)
           list(est = c(animal = v[["animal"]], blk = v[["blk"]], residual = v[["residual"]]),
                logLik = as.numeric(stats::logLik(f)), conv = TRUE) }))
}

# --------------------------------------------------------------------------- ige
# One species in 20 x 20 trays; first-order neighbourhood (8 neighbours, weight 1/d).
# DGE and IGE us (2 x 2) on iid genotypes, IEE iid at plant level.
fam_ige <- function(n, seed) {
  set.seed(seed)
  side <- 20L; n_trays <- ceiling(n / side^2)
  g <- expand.grid(col = seq_len(side), row = seq_len(side), tray = seq_len(n_trays))[seq_len(n), ]
  q <- max(30L, n %/% 5L)
  lev <- sprintf("g%04d", 1:q)
  g$gen <- factor(sample(lev, n, TRUE), levels = lev)
  key <- paste(g$tray, g$row, g$col)
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (dr in -1:1) for (dc in -1:1) {
    if (dr == 0 && dc == 0) next
    j <- match(paste(g$tray, g$row + dr, g$col + dc), key); ok <- !is.na(j)
    ii <- c(ii, which(ok)); jj <- c(jj, j[ok]); xx <- c(xx, rep(1 / sqrt(dr^2 + dc^2), sum(ok)))
  }
  Wp <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))       # plant x neighbouring plant
  Zd <- inc_levels(g$gen)
  Zn <- Wp %*% Zd; colnames(Zn) <- lev                                      # neighbouring genotypes
  X <- matrix(1, n, 1)
  truth <- c("g[1,1]" = 0.5, "g[2,1]" = -0.1, "g[2,2]" = 0.1, iee = 0.1, residual = 1.0)
  Zdn <- cbind(Zd, Zn)
  V_of <- function(e) {
    G <- mat_sigma(e, "g", 2)
    as.matrix(Zdn %*% Matrix::kronecker(methods::as(G, "CsparseMatrix"), Matrix::Diagonal(q)) %*% Matrix::t(Zdn) +
                e[["iee"]] * Matrix::tcrossprod(Wp)) + e[["residual"]] * diag(n)
  }
  y <- simulate_y(V_of(truth), X, 5)
  m <- rx_model(y, X, list(rx_term("g", list(Zd, Zn), struct = "us"), rx_term("iee", Wp, t = 1L)))
  ext <- list()
  if (n <= 3000L) ext$asreml <- function() {
    dd <- data.frame(y = y, as.matrix(Zd), as.matrix(Zn), as.matrix(Wp), check.names = FALSE)
    names(dd) <- c("y", paste0("d_", lev), paste0("n_", lev), sprintf("w_%05d", seq_len(n)))
    grp <- list(gd = 1 + seq_len(q), gn = 1 + q + seq_len(q), gw = 1 + 2 * q + seq_len(n))
    list2env(list(grp = grp, dd = dd, q = q), globalenv())   # asreml evaluates group and str() in the global environment
    fr <- stats::as.formula(sprintf("~ str(~ grp(gd) + grp(gn), ~ us(2):id(%d)) + grp(gw)", q))
    a <- asreml::asreml(y ~ 1, random = fr, group = grp, data = dd, trace = FALSE, workspace = "4gb")
    v <- asreml_vc(a)
    us_ <- v[grep("!us\\(2\\)_", names(v))]
    list(est = c("g[1,1]" = pick(us_, "1:1$"), "g[2,1]" = pick(us_, "2:1$"), "g[2,2]" = pick(us_, "2:2$"),
                 iee = pick(v, "^grp\\(gw\\)$"), residual = pick(v, "^units!R$")),
         logLik = a$loglik, conv = isTRUE(a$converge))
  }
  list(models = list(dense = m, mme = m), truth = truth, V_of = V_of, y = y, X = X, external = ext)
}

# --------------------------------------------------------------- ar1 x ar1 grids
grid_dims <- function(n) { nr <- max(5L, round(sqrt(n / 1.25))); nc <- ceiling(n / nr); c(nr, nc) }

fam_field <- function(n, seed) {
  set.seed(seed)
  dm <- grid_dims(n); nr <- dm[1]; nc <- dm[2]; n <- nr * nc
  d <- expand.grid(col = seq_len(nc), row = seq_len(nr))           # cell = (row-1)*nc + col
  q <- max(20L, n %/% 4L)
  d$gen <- factor(sample(sprintf("g%04d", 1:q), n, TRUE), levels = sprintf("g%04d", 1:q))
  d$row <- factor(d$row, levels = seq_len(nr)); d$col <- factor(d$col, levels = seq_len(nc))
  X <- matrix(1, n, 1); Zg <- inc(d$gen)
  truth <- c(gen = 0.4, field = 0.8, "field!ar1ar11" = 0.6, "field!ar1ar12" = 0.4, residual = 0.3)
  V_of <- function(e) e[["gen"]] * as.matrix(Matrix::tcrossprod(Zg)) + e[["field"]] * C_ar1ar1(nr, nc, e[["field!ar1ar11"]], e[["field!ar1ar12"]]) +
    e[["residual"]] * diag(n)
  d$y <- simulate_y(V_of(truth), X, 4)
  cell <- factor(seq_len(n), levels = seq_len(n))
  m <- rx_model(d$y, X, list(rx_term("gen", d$gen), rx_term("field", cell, level = "ar1ar1", dims = c(nr, nc))))
  list(models = list(dense = m, mme = m), truth = truth, V_of = V_of, y = d$y, X = X,
       external = list(
         asreml = function() { a <- asreml::asreml(y ~ 1, random = ~ gen + ar1v(row):ar1(col), data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(gen = pick(v, "^gen$"), field = pick(v, "^row:col!row!var$"),
                        "field!ar1ar11" = pick(v, "^row:col!row!cor$"), "field!ar1ar12" = pick(v, "^row:col!col!cor$"),
                        residual = pick(v, "^units!R$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) }))
}

fam_resid <- function(n, seed) {
  set.seed(seed)
  dm <- grid_dims(n); nr <- dm[1]; nc <- dm[2]; n <- nr * nc
  d <- expand.grid(col = seq_len(nc), row = seq_len(nr))
  q <- max(20L, n %/% 4L)
  d$gen <- factor(sample(sprintf("g%04d", 1:q), n, TRUE), levels = sprintf("g%04d", 1:q))
  d$row <- factor(d$row, levels = seq_len(nr)); d$col <- factor(d$col, levels = seq_len(nc))
  d <- d[order(d$row, d$col), ]
  X <- matrix(1, n, 1); Zg <- inc(d$gen)
  truth <- c(gen = 0.4, residual = 1.0, "residual!ar1ar11" = 0.6, "residual!ar1ar12" = 0.4)
  V_of <- function(e) e[["gen"]] * as.matrix(Matrix::tcrossprod(Zg)) + e[["residual"]] * C_ar1ar1(nr, nc, e[["residual!ar1ar11"]], e[["residual!ar1ar12"]])
  d$y <- simulate_y(V_of(truth), X, 4)
  res <- remlax:::.rx_parse_residual(~ ar1(row):ar1(col), d)
  m <- rx_model(d$y, X, list(rx_term("gen", d$gen)), residual = res)
  list(models = list(dense = m, mme = m), truth = truth, V_of = V_of, y = d$y, X = X,
       external = list(
         asreml = function() { a <- asreml::asreml(y ~ 1, random = ~ gen, residual = ~ ar1(row):ar1(col), data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(gen = pick(v, "^gen$"), residual = pick(v, "^row:col!R$"),
                        "residual!ar1ar11" = pick(v, "^row:col!row!cor$"), "residual!ar1ar12" = pick(v, "^row:col!col!cor$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) }))
}

# --------------------------------------------------------------------------- exp
# Transect: q positions (3 observations per position), correlation phi^d.
fam_exp <- function(n, seed) {
  set.seed(seed)
  q <- n %/% 3L; n <- 3L * q
  x <- sort(stats::runif(q, 0, q / 4))
  pos <- factor(rep(sprintf("s%05d", 1:q), each = 3), levels = sprintf("s%05d", 1:q))
  d <- data.frame(pos = pos, xpos = rep(x, each = 3))
  X <- matrix(1, n, 1); Zs <- inc(pos)
  truth <- c(s = 1.0, "s!exp" = 0.7, residual = 0.5)
  D <- abs(outer(x, x, "-"))
  ZZ <- function(C) as.matrix(Zs %*% C %*% Matrix::t(Zs))
  V_of <- function(e) e[["s"]] * ZZ(e[["s!exp"]]^D) + e[["residual"]] * diag(n)
  d$y <- simulate_y(V_of(truth), X, 1)
  m <- rx_model(d$y, X, list(rx_term("s", pos, level = "exp", coord = x)))
  list(models = list(dense = m, mme = m), truth = truth, V_of = V_of, y = d$y, X = X,
       external = list(
         asreml = function() { a <- asreml::asreml(y ~ 1, random = ~ expv(xpos), data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(s = pick(v, "^xpos!xpos!var$"), "s!exp" = pick(v, "^xpos!xpos!pow$"),
                        residual = pick(v, "^units!R$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) }))
}

# ---------------------------------------------------------------------------- fa
fam_fa <- function(n, seed) {
  set.seed(seed)
  t <- 5L; nu <- n %/% t; q <- max(40L, nu %/% 4L)
  Lam <- matrix(c(0.7, 0.5, -0.3, 0.6, 0.4, 0, 0.5, 0.6, -0.4, 0.3), 5)
  Sg <- Lam %*% t(Lam) + diag(c(0.15, 0.2, 0.2, 0.15, 0.25))
  sr <- c(0.8, 0.6, 0.7, 0.5, 0.6)
  tn <- sprintf("y%d", 1:t)
  gen <- factor(sample(sprintf("g%04d", 1:q), nu, TRUE), levels = sprintf("g%04d", 1:q))
  dl <- data.frame(unit = rep(seq_len(nu), t), trait = factor(rep(tn, each = nu), levels = tn), gen = rep(gen, t))
  X <- stats::model.matrix(~ trait - 1, dl)
  Zl <- lapply(tn, function(a) inc_levels(dl$gen) * (dl$trait == a))
  truth <- c(vec_sigma("g", Sg), stats::setNames(sr, sprintf("residual[%d,%d]", 1:t, 1:t)))
  Zall <- do.call(cbind, Zl); tr <- as.integer(dl$trait)
  V_of <- function(e) {
    G <- mat_sigma(e, "g", t)
    as.matrix(Zall %*% Matrix::kronecker(methods::as(G, "CsparseMatrix"), Matrix::Diagonal(q)) %*% Matrix::t(Zall)) +
      diag(vapply(tr, function(a) e[[sprintf("residual[%d,%d]", a, a)]], 0))
  }
  dl$y <- simulate_y(V_of(truth), X, seq_len(t))
  w <- data.frame(gen = gen); for (a in tn) w[[a]] <- dl$y[dl$trait == a]
  m <- rx_model(dl$y, X, list(rx_term("g", Zl, struct = "fa", rank = 2L)),
                residual = rx_residual("diag", trait = dl$trait, unit = dl$unit))
  list(models = list(dense = m, mme = m), truth = truth, V_of = V_of, y = dl$y, X = X,
       external = list(
         asreml = function() { a <- asreml::asreml(cbind(y1, y2, y3, y4, y5) ~ trait, random = ~ fa(trait, 2):gen,
                                                   residual = ~ id(units):diag(trait), data = w, trace = FALSE, maxit = 100)
           v <- asreml_vc(a)
           psi <- vapply(tn, function(a_) pick(v, sprintf("^fa\\(trait, 2\\):gen!%s!var$", a_)), 0)
           L <- cbind(vapply(tn, function(a_) pick(v, sprintf("^fa\\(trait, 2\\):gen!%s!fa1$", a_)), 0),
                      vapply(tn, function(a_) pick(v, sprintf("^fa\\(trait, 2\\):gen!%s!fa2$", a_)), 0))
           G <- L %*% t(L) + diag(psi)
           r <- vapply(tn, function(a_) pick(v, sprintf("^units:trait!trait_%s$", a_)), 0)
           list(est = c(vec_sigma("g", G), stats::setNames(r, sprintf("residual[%d,%d]", 1:t, 1:t))),
                logLik = a$loglik, conv = isTRUE(a$converge)) }))
}

BUILDERS <- list(iid = fam_iid, grm = fam_grm, us = fam_us, ped = fam_ped, ige = fam_ige,
                      field = fam_field, resid = fam_resid, exp = fam_exp, fa = fam_fa)

# families of the "number of parameters" benchmark (mtus<t>, mtige<t>), in the
# same logical file: loaded here so that child processes which source this
# file know them too
local({
  # path of this file: the source() frame that carries `ofile`, wherever it is
  # in the stack (source() called inside a function, as in a child process)
  f <- ""
  for (i in rev(seq_len(sys.nframe()))) {
    o <- sys.frame(i)$ofile
    if (!is.null(o)) { f <- file.path(dirname(normalizePath(o)), "familles_params.R"); break }
  }
  if (nzchar(f) && file.exists(f)) source(f, local = FALSE)
})
