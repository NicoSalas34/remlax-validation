# ==============================================================================
# Familles du banc des moteurs. Chaque constructeur prend (n, graine) et rend :
#   modeles   rx_model pour remlax : `dense` (voies denses) et `mme` (voies MME,
#             souvent le meme objet ; pedigree : K pour le dense, Kinv pour le MME)
#   vrai      parametres vrais, NOMS COMMUNS a tous les logiciels
#   V_de(est) V construite a partir d'un vecteur de parametres nommes (arbitre)
#   y, X      reponse et effets fixes, pour l'arbitre
#   externes  liste nommee de fonctions sans argument -> list(est, logLik, conv)
#             pour asreml, sommer et lme4 quand ils s'appliquent
# Noms communs : composante simple "gen", matrice "g[i,j]" (i >= j),
# residuelle "residuelle" ou "residuelle[i,j]", parametres de niveaux
# "<terme>!<nom>".
# ==============================================================================

TAILLES <- c(500L, 1000L, 1500L, 2000L, 3000L, 4000L, 6000L, 8000L)
FAMILLES <- c("iid", "grm", "us", "ped", "ige", "champ", "resid", "exp", "fa")
graine_de <- function(fam, n, rep) as.integer(20261022L + 100000L * match(fam, FAMILLES) + 10L * n + rep)

# Incidence qui GARDE tous les niveaux du facteur (inc() de la suite les
# refactorise et perd les niveaux non tires, ce qui desaligne Z et K).
inc_niv <- function(f) Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
                                            dims = c(length(f), nlevels(f)), dimnames = list(NULL, levels(f)))

vec_sigma <- function(nom, S) {
  S <- as.matrix(S); t <- nrow(S)
  if (t == 1L) return(stats::setNames(S[1, 1], nom))
  out <- c()
  for (i in seq_len(t)) for (j in seq_len(i)) out[sprintf("%s[%d,%d]", nom, i, j)] <- S[i, j]
  out
}
mat_sigma <- function(est, nom, t) {
  if (t == 1L) return(matrix(est[[nom]], 1, 1))
  S <- matrix(0, t, t)
  for (i in seq_len(t)) for (j in seq_len(i)) S[i, j] <- S[j, i] <- est[[sprintf("%s[%d,%d]", nom, i, j)]]
  S
}

# Composantes de variance d'asreml : tableau (nom de ligne -> composante)
asreml_vc <- function(a) { vc <- summary(a)$varcomp; stats::setNames(vc$component, rownames(vc)) }
pick <- function(v, motif) {
  i <- grep(motif, names(v), perl = TRUE)
  if (length(i) != 1L) stop("composante introuvable ou ambigue : ", motif, " parmi ", paste(names(v), collapse = ", "))
  unname(v[i])
}
sommer_vc <- function(m) { vc <- summary(m)$varcomp; stats::setNames(vc$VarComp, rownames(vc)) }

# lme4 avec une parente K = L L' : Z u = (Z L) v, v iid (le modele marginal est le meme)
lmer_parente <- function(formule, data, facteur, K) {
  lf <- lme4::lFormula(formule, data = data, REML = TRUE)
  idx <- which(names(lf$reTrms$cnms) == facteur)
  deb <- c(0, cumsum(vapply(lf$reTrms$Ztlist, nrow, 0)))
  rows <- (deb[idx] + 1):deb[idx + 1]
  niv <- rownames(lf$reTrms$Zt)[rows]          # lme4 ecarte les niveaux absents
  L <- t(chol(K[niv, niv]))
  Zt <- lf$reTrms$Zt
  Zt[rows, ] <- methods::as(Matrix::t(L) %*% Zt[rows, ], "CsparseMatrix")
  lf$reTrms$Zt <- methods::as(Zt, "CsparseMatrix")
  dev <- do.call(lme4::mkLmerDevfun, lf)
  opt <- lme4::optimizeLmer(dev)
  lme4::mkMerMod(environment(dev), opt, lf$reTrms, fr = lf$fr)
}
lme4_vc <- function(f) {
  vc <- as.data.frame(lme4::VarCorr(f))
  stats::setNames(vc$vcov, ifelse(vc$grp == "Residual", "residuelle", vc$grp))
}

# --------------------------------------------------------------------------- iid
fam_iid <- function(n, graine) {
  set.seed(graine)
  q <- max(20L, n %/% 5L)
  d <- data.frame(gen = factor(sample(sprintf("g%04d", 1:q), n, TRUE), levels = sprintf("g%04d", 1:q)),
                  blk = factor(sample(sprintf("b%02d", 1:25), n, TRUE), levels = sprintf("b%02d", 1:25)),
                  x = stats::rnorm(n))
  X <- cbind(1, d$x); Zg <- inc(d$gen); Zb <- inc(d$blk)
  vrai <- c(gen = 0.5, blk = 0.3, residuelle = 1.0)
  V_de <- function(e) as.matrix(e[["gen"]] * Matrix::tcrossprod(Zg) + e[["blk"]] * Matrix::tcrossprod(Zb)) +
    e[["residuelle"]] * diag(n)
  d$y <- simuler_y(V_de(vrai), X, c(2, 0.5))
  m <- rx_model(d$y, X, list(rx_term("gen", d$gen), rx_term("blk", d$blk)))
  list(modeles = list(dense = m, mme = m), vrai = vrai, V_de = V_de, y = d$y, X = X,
       externes = list(
         asreml = function() { a <- asreml::asreml(y ~ x, random = ~ gen + blk, data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(gen = pick(v, "^gen$"), blk = pick(v, "^blk$"), residuelle = pick(v, "^units!R$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) },
         sommer = function() { s <- sommer::mmer(y ~ x, random = ~ gen + blk, rcov = ~ units, data = d, verbose = FALSE)
           v <- sommer_vc(s)
           list(est = c(gen = pick(v, "^gen"), blk = pick(v, "^blk"), residuelle = pick(v, "^units")),
                logLik = as.numeric(utils::tail(s$monitor[1, ], 1)), conv = isTRUE(s$convergence)) },
         lme4 = function() { f <- lme4::lmer(y ~ x + (1 | gen) + (1 | blk), data = d, REML = TRUE)
           v <- lme4_vc(f)
           list(est = c(gen = v[["gen"]], blk = v[["blk"]], residuelle = v[["residuelle"]]),
                logLik = as.numeric(stats::logLik(f)), conv = length(f@optinfo$conv$lme4$messages) == 0) }))
}

# --------------------------------------------------------------------------- grm
fam_grm <- function(n, graine) {
  set.seed(graine)
  q <- max(50L, n %/% 4L)
  K <- grm_simulee(q, prefixe = "g")
  niv <- rownames(K)
  d <- data.frame(gen = factor(sample(niv, n, TRUE), levels = niv),
                  blk = factor(sample(sprintf("b%02d", 1:25), n, TRUE), levels = sprintf("b%02d", 1:25)),
                  x = stats::rnorm(n))
  X <- cbind(1, d$x); Zg <- inc_niv(d$gen); Zb <- inc(d$blk)
  vrai <- c(gen = 0.5, blk = 0.3, residuelle = 1.0)
  V_de <- function(e) e[["gen"]] * as.matrix(Zg %*% K %*% Matrix::t(Zg)) +
    as.matrix(e[["blk"]] * Matrix::tcrossprod(Zb)) + e[["residuelle"]] * diag(n)
  d$y <- simuler_y(V_de(vrai), X, c(2, 0.5))
  Zl <- list(Zg); colnames(Zl[[1]]) <- niv
  m <- rx_model(d$y, X, list(rx_term("gen", Zl, K = K), rx_term("blk", d$blk)))
  list(modeles = list(dense = m, mme = m), vrai = vrai, V_de = V_de, y = d$y, X = X,
       externes = list(
         asreml = function() { a <- asreml::asreml(y ~ x, random = ~ vm(gen, K) + blk, data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(gen = pick(v, "^vm\\(gen, K\\)$"), blk = pick(v, "^blk$"), residuelle = pick(v, "^units!R$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) },
         sommer = function() { s <- sommer::mmer(y ~ x, random = ~ sommer::vsr(gen, Gu = K) + blk, rcov = ~ units,
                                                 data = d, verbose = FALSE)
           v <- sommer_vc(s)
           list(est = c(gen = pick(v, "^u:gen"), blk = pick(v, "^blk"), residuelle = pick(v, "^units")),
                logLik = as.numeric(utils::tail(s$monitor[1, ], 1)), conv = isTRUE(s$convergence)) },
         lme4 = function() { f <- lmer_parente(y ~ x + (1 | gen) + (1 | blk), d, "gen", K)
           v <- lme4_vc(f)
           list(est = c(gen = v[["gen"]], blk = v[["blk"]], residuelle = v[["residuelle"]]),
                logLik = as.numeric(stats::logLik(f)), conv = TRUE) }))
}

# ---------------------------------------------------------------------------- us
fam_us <- function(n, graine) {
  set.seed(graine)
  nu <- n %/% 3L; q <- max(30L, nu %/% 4L)
  Sg <- matrix(c(1, .5, .2, .5, .8, .3, .2, .3, .6), 3)
  Sr <- matrix(c(1, .3, .1, .3, .9, .2, .1, .2, .7), 3)
  gen <- factor(sample(sprintf("g%04d", 1:q), nu, TRUE), levels = sprintf("g%04d", 1:q))
  dl <- data.frame(unit = rep(seq_len(nu), 3), trait = factor(rep(c("y1", "y2", "y3"), each = nu)),
                   gen = rep(gen, 3))
  X <- stats::model.matrix(~ trait - 1, dl)
  Zl <- lapply(c("y1", "y2", "y3"), function(a) inc_niv(dl$gen) * (dl$trait == a))
  vrai <- c(vec_sigma("g", Sg), vec_sigma("residuelle", Sr))
  Zall <- do.call(cbind, Zl); tr <- as.integer(dl$trait); Zu <- inc(factor(dl$unit))
  V_de <- function(e) {
    G <- mat_sigma(e, "g", 3); R <- mat_sigma(e, "residuelle", 3)
    # V = Z (G (x) I) Z' + R par unite : R[tr_i, tr_j] pour deux lignes de la meme unite
    Vg <- as.matrix(Zall %*% Matrix::kronecker(methods::as(G, "CsparseMatrix"), Matrix::Diagonal(q)) %*% Matrix::t(Zall))
    Rm <- Matrix::sparseMatrix(i = seq_along(tr), j = tr, x = 1, dims = c(length(tr), 3))
    Vg + as.matrix(Matrix::tcrossprod(Zu)) * as.matrix(Rm %*% R %*% Matrix::t(Rm))
  }
  dl$y <- simuler_y(V_de(vrai), X, c(1, 2, 3))
  w <- data.frame(gen = gen, y1 = dl$y[dl$trait == "y1"], y2 = dl$y[dl$trait == "y2"], y3 = dl$y[dl$trait == "y3"])
  m <- rx_model(dl$y, X, list(rx_term("g", Zl, struct = "us")),
                residual = rx_residual("us", trait = dl$trait, unit = dl$unit))
  list(modeles = list(dense = m, mme = m), vrai = vrai, V_de = V_de, y = dl$y, X = X,
       externes = list(
         asreml = function() { a <- asreml::asreml(cbind(y1, y2, y3) ~ trait, random = ~ us(trait):gen,
                                                   residual = ~ id(units):us(trait), data = w, trace = FALSE)
           v <- asreml_vc(a); e <- c()
           tn <- c("y1", "y2", "y3")
           for (i in 1:3) for (j in 1:i) {
             e[sprintf("g[%d,%d]", i, j)] <- pick(v, sprintf("^trait:gen!trait_%s:%s$", tn[i], tn[j]))
             e[sprintf("residuelle[%d,%d]", i, j)] <- pick(v, sprintf("^units:trait!trait_%s:%s$", tn[i], tn[j]))
           }
           list(est = e, logLik = a$loglik, conv = isTRUE(a$converge)) },
         sommer = function() { s <- sommer::mmer(cbind(y1, y2, y3) ~ 1, random = ~ sommer::vsr(gen, Gtc = sommer::unsm(3)),
                                                 rcov = ~ sommer::vsr(units, Gtc = sommer::unsm(3)), data = w, verbose = FALSE)
           e <- c(vec_sigma("g", s$sigma[[1]]), vec_sigma("residuelle", s$sigma[[2]]))
           list(est = e, logLik = as.numeric(utils::tail(s$monitor[1, ], 1)), conv = isTRUE(s$convergence)) }))
}

# --------------------------------------------------------------------------- ped
fam_ped <- function(n, graine) {
  set.seed(graine)
  f0 <- max(50L, n %/% 10L); ngen <- 4L
  par_gen <- ceiling((n - f0) / ngen)
  ped <- data.frame(id = seq_len(f0), sire = 0L, dam = 0L)
  for (g in seq_len(ngen)) {
    prev <- ped$id; k <- min(par_gen, n - nrow(ped)); if (k <= 0) break
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
  vrai <- c(animal = 0.4, blk = 0.2, residuelle = 0.6)
  V_de <- function(e) e[["animal"]] * A + e[["blk"]] * as.matrix(Matrix::tcrossprod(Zb)) + e[["residuelle"]] * diag(n)
  d$y <- simuler_y(V_de(vrai), X, c(3, 0.5))
  Za <- list(Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1, dims = c(n, n), dimnames = list(NULL, lab)))
  m_dense <- rx_model(d$y, X, list(rx_term("animal", Za, K = A), rx_term("blk", d$blk)))
  m_mme <- rx_model(d$y, X, list(rx_term("animal", Za, Kinv = Ainv, Kinv_logdet = ldA), rx_term("blk", d$blk)))
  list(modeles = list(dense = m_dense, mme = m_mme), vrai = vrai, V_de = V_de, y = d$y, X = X,
       externes = list(
         asreml = function() { a <- asreml::asreml(y ~ x, random = ~ vm(animal, ai) + blk, data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(animal = pick(v, "^vm\\(animal, ai\\)$"), blk = pick(v, "^blk$"), residuelle = pick(v, "^units!R$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) },
         sommer = function() { s <- sommer::mmer(y ~ x, random = ~ sommer::vsr(animal, Gu = A) + blk, rcov = ~ units,
                                                 data = d, verbose = FALSE)
           v <- sommer_vc(s)
           list(est = c(animal = pick(v, "^u:animal"), blk = pick(v, "^blk"), residuelle = pick(v, "^units")),
                logLik = as.numeric(utils::tail(s$monitor[1, ], 1)), conv = isTRUE(s$convergence)) },
         lme4 = function() {
           suppressPackageStartupMessages(library(lme4))   # pedigreemm appelle lmer sans prefixe
           pd <- pedigreemm::pedigree(sire = ifelse(ped$sire == 0, NA, ped$sire), dam = ifelse(ped$dam == 0, NA, ped$dam), label = ped$id)
           f <- pedigreemm::pedigreemm(y ~ x + (1 | animal) + (1 | blk), pedigree = list(animal = pd), data = d, REML = TRUE,
             control = lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
           v <- lme4_vc(f)
           list(est = c(animal = v[["animal"]], blk = v[["blk"]], residuelle = v[["residuelle"]]),
                logLik = as.numeric(stats::logLik(f)), conv = TRUE) }))
}

# --------------------------------------------------------------------------- ige
# Une espece en bacs de 20 x 20 ; voisinage d'ordre 1 (8 voisins, poids 1/d).
# DGE et IGE us (2 x 2) sur genotypes iid, IEE iid au niveau plante.
fam_ige <- function(n, graine) {
  set.seed(graine)
  cote <- 20L; nbac <- ceiling(n / cote^2)
  g <- expand.grid(col = seq_len(cote), row = seq_len(cote), bac = seq_len(nbac))[seq_len(n), ]
  q <- max(30L, n %/% 5L)
  niv <- sprintf("g%04d", 1:q)
  g$gen <- factor(sample(niv, n, TRUE), levels = niv)
  cle <- paste(g$bac, g$row, g$col)
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (dr in -1:1) for (dc in -1:1) {
    if (dr == 0 && dc == 0) next
    j <- match(paste(g$bac, g$row + dr, g$col + dc), cle); ok <- !is.na(j)
    ii <- c(ii, which(ok)); jj <- c(jj, j[ok]); xx <- c(xx, rep(1 / sqrt(dr^2 + dc^2), sum(ok)))
  }
  Wp <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))       # plante x plante voisine
  Zd <- inc_niv(g$gen)
  Zn <- Wp %*% Zd; colnames(Zn) <- niv                                      # genotypes voisins
  X <- matrix(1, n, 1)
  vrai <- c("g[1,1]" = 0.5, "g[2,1]" = -0.1, "g[2,2]" = 0.1, iee = 0.1, residuelle = 1.0)
  Zdn <- cbind(Zd, Zn)
  V_de <- function(e) {
    G <- mat_sigma(e, "g", 2)
    as.matrix(Zdn %*% Matrix::kronecker(methods::as(G, "CsparseMatrix"), Matrix::Diagonal(q)) %*% Matrix::t(Zdn) +
                e[["iee"]] * Matrix::tcrossprod(Wp)) + e[["residuelle"]] * diag(n)
  }
  y <- simuler_y(V_de(vrai), X, 5)
  m <- rx_model(y, X, list(rx_term("g", list(Zd, Zn), struct = "us"), rx_term("iee", Wp, t = 1L)))
  ext <- list()
  if (n <= 3000L) ext$asreml <- function() {
    dd <- data.frame(y = y, as.matrix(Zd), as.matrix(Zn), as.matrix(Wp), check.names = FALSE)
    names(dd) <- c("y", paste0("d_", niv), paste0("n_", niv), sprintf("w_%05d", seq_len(n)))
    grp <- list(gd = 1 + seq_len(q), gn = 1 + q + seq_len(q), gw = 1 + 2 * q + seq_len(n))
    list2env(list(grp = grp, dd = dd, q = q), globalenv())   # asreml evalue group et str() dans l'environnement global
    fr <- stats::as.formula(sprintf("~ str(~ grp(gd) + grp(gn), ~ us(2):id(%d)) + grp(gw)", q))
    a <- asreml::asreml(y ~ 1, random = fr, group = grp, data = dd, trace = FALSE, workspace = "4gb")
    v <- asreml_vc(a)
    us_ <- v[grep("!us\\(2\\)_", names(v))]
    list(est = c("g[1,1]" = pick(us_, "1:1$"), "g[2,1]" = pick(us_, "2:1$"), "g[2,2]" = pick(us_, "2:2$"),
                 iee = pick(v, "^grp\\(gw\\)$"), residuelle = pick(v, "^units!R$")),
         logLik = a$loglik, conv = isTRUE(a$converge))
  }
  list(modeles = list(dense = m, mme = m), vrai = vrai, V_de = V_de, y = y, X = X, externes = ext)
}

# ------------------------------------------------------------- grilles ar1 x ar1
grille <- function(n) { nr <- max(5L, round(sqrt(n / 1.25))); nc <- ceiling(n / nr); c(nr, nc) }

fam_champ <- function(n, graine) {
  set.seed(graine)
  dm <- grille(n); nr <- dm[1]; nc <- dm[2]; n <- nr * nc
  d <- expand.grid(col = seq_len(nc), row = seq_len(nr))           # cellule = (row-1)*nc + col
  q <- max(20L, n %/% 4L)
  d$gen <- factor(sample(sprintf("g%04d", 1:q), n, TRUE), levels = sprintf("g%04d", 1:q))
  d$row <- factor(d$row, levels = seq_len(nr)); d$col <- factor(d$col, levels = seq_len(nc))
  X <- matrix(1, n, 1); Zg <- inc(d$gen)
  vrai <- c(gen = 0.4, champ = 0.8, "champ!ar1ar11" = 0.6, "champ!ar1ar12" = 0.4, residuelle = 0.3)
  V_de <- function(e) e[["gen"]] * as.matrix(Matrix::tcrossprod(Zg)) + e[["champ"]] * C_ar1ar1(nr, nc, e[["champ!ar1ar11"]], e[["champ!ar1ar12"]]) +
    e[["residuelle"]] * diag(n)
  d$y <- simuler_y(V_de(vrai), X, 4)
  cell <- factor(seq_len(n), levels = seq_len(n))
  m <- rx_model(d$y, X, list(rx_term("gen", d$gen), rx_term("champ", cell, level = "ar1ar1", dims = c(nr, nc))))
  list(modeles = list(dense = m, mme = m), vrai = vrai, V_de = V_de, y = d$y, X = X,
       externes = list(
         asreml = function() { a <- asreml::asreml(y ~ 1, random = ~ gen + ar1v(row):ar1(col), data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(gen = pick(v, "^gen$"), champ = pick(v, "^row:col!row!var$"),
                        "champ!ar1ar11" = pick(v, "^row:col!row!cor$"), "champ!ar1ar12" = pick(v, "^row:col!col!cor$"),
                        residuelle = pick(v, "^units!R$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) }))
}

fam_resid <- function(n, graine) {
  set.seed(graine)
  dm <- grille(n); nr <- dm[1]; nc <- dm[2]; n <- nr * nc
  d <- expand.grid(col = seq_len(nc), row = seq_len(nr))
  q <- max(20L, n %/% 4L)
  d$gen <- factor(sample(sprintf("g%04d", 1:q), n, TRUE), levels = sprintf("g%04d", 1:q))
  d$row <- factor(d$row, levels = seq_len(nr)); d$col <- factor(d$col, levels = seq_len(nc))
  d <- d[order(d$row, d$col), ]
  X <- matrix(1, n, 1); Zg <- inc(d$gen)
  vrai <- c(gen = 0.4, residuelle = 1.0, "residuelle!ar1ar11" = 0.6, "residuelle!ar1ar12" = 0.4)
  V_de <- function(e) e[["gen"]] * as.matrix(Matrix::tcrossprod(Zg)) + e[["residuelle"]] * C_ar1ar1(nr, nc, e[["residuelle!ar1ar11"]], e[["residuelle!ar1ar12"]])
  d$y <- simuler_y(V_de(vrai), X, 4)
  res <- remlax:::.rx_parse_residual(~ ar1(row):ar1(col), d)
  m <- rx_model(d$y, X, list(rx_term("gen", d$gen)), residual = res)
  list(modeles = list(dense = m, mme = m), vrai = vrai, V_de = V_de, y = d$y, X = X,
       externes = list(
         asreml = function() { a <- asreml::asreml(y ~ 1, random = ~ gen, residual = ~ ar1(row):ar1(col), data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(gen = pick(v, "^gen$"), residuelle = pick(v, "^row:col!R$"),
                        "residuelle!ar1ar11" = pick(v, "^row:col!row!cor$"), "residuelle!ar1ar12" = pick(v, "^row:col!col!cor$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) }))
}

# --------------------------------------------------------------------------- exp
# Transect : q positions (3 observations par position), correlation phi^d.
fam_exp <- function(n, graine) {
  set.seed(graine)
  q <- n %/% 3L; n <- 3L * q
  x <- sort(stats::runif(q, 0, q / 4))
  pos <- factor(rep(sprintf("s%05d", 1:q), each = 3), levels = sprintf("s%05d", 1:q))
  d <- data.frame(pos = pos, xpos = rep(x, each = 3))
  X <- matrix(1, n, 1); Zs <- inc(pos)
  vrai <- c(s = 1.0, "s!exp" = 0.7, residuelle = 0.5)
  D <- abs(outer(x, x, "-"))
  ZZ <- function(C) as.matrix(Zs %*% C %*% Matrix::t(Zs))
  V_de <- function(e) e[["s"]] * ZZ(e[["s!exp"]]^D) + e[["residuelle"]] * diag(n)
  d$y <- simuler_y(V_de(vrai), X, 1)
  m <- rx_model(d$y, X, list(rx_term("s", pos, level = "exp", coord = x)))
  list(modeles = list(dense = m, mme = m), vrai = vrai, V_de = V_de, y = d$y, X = X,
       externes = list(
         asreml = function() { a <- asreml::asreml(y ~ 1, random = ~ expv(xpos), data = d, trace = FALSE)
           v <- asreml_vc(a)
           list(est = c(s = pick(v, "^xpos!xpos!var$"), "s!exp" = pick(v, "^xpos!xpos!pow$"),
                        residuelle = pick(v, "^units!R$")),
                logLik = a$loglik, conv = isTRUE(a$converge)) }))
}

# ---------------------------------------------------------------------------- fa
fam_fa <- function(n, graine) {
  set.seed(graine)
  t <- 5L; nu <- n %/% t; q <- max(40L, nu %/% 4L)
  Lam <- matrix(c(0.7, 0.5, -0.3, 0.6, 0.4, 0, 0.5, 0.6, -0.4, 0.3), 5)
  Sg <- Lam %*% t(Lam) + diag(c(0.15, 0.2, 0.2, 0.15, 0.25))
  sr <- c(0.8, 0.6, 0.7, 0.5, 0.6)
  tn <- sprintf("y%d", 1:t)
  gen <- factor(sample(sprintf("g%04d", 1:q), nu, TRUE), levels = sprintf("g%04d", 1:q))
  dl <- data.frame(unit = rep(seq_len(nu), t), trait = factor(rep(tn, each = nu), levels = tn), gen = rep(gen, t))
  X <- stats::model.matrix(~ trait - 1, dl)
  Zl <- lapply(tn, function(a) inc_niv(dl$gen) * (dl$trait == a))
  vrai <- c(vec_sigma("g", Sg), stats::setNames(sr, sprintf("residuelle[%d,%d]", 1:t, 1:t)))
  Zall <- do.call(cbind, Zl); tr <- as.integer(dl$trait)
  V_de <- function(e) {
    G <- mat_sigma(e, "g", t)
    as.matrix(Zall %*% Matrix::kronecker(methods::as(G, "CsparseMatrix"), Matrix::Diagonal(q)) %*% Matrix::t(Zall)) +
      diag(vapply(tr, function(a) e[[sprintf("residuelle[%d,%d]", a, a)]], 0))
  }
  dl$y <- simuler_y(V_de(vrai), X, seq_len(t))
  w <- data.frame(gen = gen); for (a in tn) w[[a]] <- dl$y[dl$trait == a]
  m <- rx_model(dl$y, X, list(rx_term("g", Zl, struct = "fa", rank = 2L)),
                residual = rx_residual("diag", trait = dl$trait, unit = dl$unit))
  list(modeles = list(dense = m, mme = m), vrai = vrai, V_de = V_de, y = dl$y, X = X,
       externes = list(
         asreml = function() { a <- asreml::asreml(cbind(y1, y2, y3, y4, y5) ~ trait, random = ~ fa(trait, 2):gen,
                                                   residual = ~ id(units):diag(trait), data = w, trace = FALSE, maxit = 100)
           v <- asreml_vc(a)
           psi <- vapply(tn, function(a_) pick(v, sprintf("^fa\\(trait, 2\\):gen!%s!var$", a_)), 0)
           L <- cbind(vapply(tn, function(a_) pick(v, sprintf("^fa\\(trait, 2\\):gen!%s!fa1$", a_)), 0),
                      vapply(tn, function(a_) pick(v, sprintf("^fa\\(trait, 2\\):gen!%s!fa2$", a_)), 0))
           G <- L %*% t(L) + diag(psi)
           r <- vapply(tn, function(a_) pick(v, sprintf("^units:trait!trait_%s$", a_)), 0)
           list(est = c(vec_sigma("g", G), stats::setNames(r, sprintf("residuelle[%d,%d]", 1:t, 1:t))),
                logLik = a$loglik, conv = isTRUE(a$converge)) }))
}

CONSTRUCTEURS <- list(iid = fam_iid, grm = fam_grm, us = fam_us, ped = fam_ped, ige = fam_ige,
                      champ = fam_champ, resid = fam_resid, exp = fam_exp, fa = fam_fa)

# familles du banc "nombre de parametres" (mtus<t>, mtige<t>), dans le meme
# fichier logique : chargees ici pour que les processus fils qui sourcent ce
# fichier les connaissent aussi
local({
  # chemin de ce fichier : le cadre de source() qui porte `ofile`, ou qu'il soit
  # dans la pile (source() appele dans une fonction, comme dans un processus fils)
  f <- ""
  for (i in rev(seq_len(sys.nframe()))) {
    o <- sys.frame(i)$ofile
    if (!is.null(o)) { f <- file.path(dirname(normalizePath(o)), "familles_params.R"); break }
  }
  if (nzchar(f) && file.exists(f)) source(f, local = FALSE)
})
