# ==============================================================================
# Filet de securite du nettoyage de remlax : resultats de reference
#
# Ajuste un panel de modeles qui parcourt les chemins du solveur (moteurs dense,
# econome, MME dense et creux, blocs ; optimiseurs L-BFGS-B et AI ; structures
# us, fa, ar1 x ar1, spl2d, termes explicites, decalage spatial ; sorties PEV,
# Wald, Kenward-Roger, vpredict, prediction ; scan, LOCO, reajustement exact ;
# voie bayesienne ; evaluation a theta fixe ; aides de voisinage) et enregistre
# chaque objet rendu. Le nettoyage ne doit changer AUCUN nombre : on compare
# ensuite deux captures champ par champ.
#
# Les champs qui varient d'une execution a l'autre du MEME code (temps, chemins)
# sont reperes en capturant deux fois l'etat de reference, puis ignores.
#
#   Rscript golden.R capture <dossier> [chemin_du_paquet]
#   Rscript golden.R compare <dossier_reference> <dossier_nouveau> [<dossier_reference_bis>]
#
# CPU seulement, un fil de calcul : la vraisemblance depend au sixieme chiffre
# du nombre de fils BLAS, et le GPU n'est pas deterministe.
# ==============================================================================

a0 <- commandArgs(trailingOnly = TRUE)
mode <- a0[1]

# ------------------------------------------------------------------------------
# Capture
# ------------------------------------------------------------------------------
if (mode == "capture") {
  dossier <- a0[2]
  paquet <- if (length(a0) >= 3) a0[3] else "/home/salas/Documents/Git/remlax"
  dir.create(dossier, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
             REMLAX_THREADS = "1", JAX_PLATFORMS = "cpu", REMLAX_JAX_CACHE = "0")
  suppressPackageStartupMessages({ pkgload::load_all(paquet, quiet = TRUE); library(Matrix) })
  options(remlax.threads = 1L)
  `%||%` <- function(a, b) if (is.null(a)) b else a

  journal <- character()
  cas <- function(nom, expr) {
    t0 <- Sys.time()
    sortie <- character()
    val <- tryCatch({
      sortie <- utils::capture.output(v <- force(expr), type = "output")
      v
    }, error = function(e) structure(list(message = conditionMessage(e)), class = "golden_erreur"))
    affichage <- if (inherits(val, "rx_fit"))
      tryCatch(c(utils::capture.output(print(val)), utils::capture.output(print(summary(val)))),
               error = function(e) paste("affichage :", conditionMessage(e)))
    else NULL
    saveRDS(list(valeur = val, affichage = affichage), file.path(dossier, paste0(nom, ".rds")))
    etat <- if (inherits(val, "golden_erreur")) paste("ERREUR :", val$message) else "ok"
    ligne <- sprintf("%-28s %6.1f s  %s", nom, as.numeric(difftime(Sys.time(), t0, units = "secs")), etat)
    journal <<- c(journal, ligne); message(ligne)
    invisible(val)
  }
  fa <- function(...) c(list(backend = "cpu", verbose = FALSE), list(...))

  # --- donnees -----------------------------------------------------------------
  set.seed(2026)
  d1 <- expand.grid(gid = factor(sprintf("G%02d", 1:60)), bloc = factor(1:4))
  d1$x <- rnorm(nrow(d1))
  d1$y <- 12 + 0.4 * d1$x + rnorm(60, 0, 1.2)[d1$gid] + rnorm(4, 0, 0.6)[d1$bloc] + rnorm(nrow(d1))

  set.seed(11)
  q <- 80; p <- 300
  M <- matrix(sample(c(-1, 1), q * p, replace = TRUE), q, p,
              dimnames = list(sprintf("L%02d", 1:q), sprintf("snp%03d", 1:p)))
  carte <- data.frame(snp = colnames(M), chr = rep(1:3, each = 100), pos = rep(1:100, 3))
  Wc <- scale(M, scale = FALSE); K <- tcrossprod(Wc) / p + diag(0.01, q)
  dimnames(K) <- list(rownames(M), rownames(M))
  u <- as.numeric(t(chol(K)) %*% rnorm(q))
  d2 <- data.frame(gid = factor(rep(rownames(M), each = 3), levels = rownames(M)))
  d2$y <- 3 + 0.8 * M[as.character(d2$gid), 40] + u[as.integer(d2$gid)] + rnorm(nrow(d2))

  set.seed(5)
  t3 <- 3; q3 <- 50
  G3 <- matrix(c(1, .5, .2, .5, 1, .4, .2, .4, 1), 3)
  u3 <- matrix(rnorm(q3 * t3), q3) %*% chol(G3)
  d3 <- expand.grid(rep = 1:2, gid = factor(sprintf("H%02d", 1:q3)), trait = factor(paste0("T", 1:t3)))
  d3$unite <- factor(paste(d3$gid, d3$rep))
  d3$y <- as.numeric(d3$trait) * 2 + u3[cbind(as.integer(d3$gid), as.integer(d3$trait))] + rnorm(nrow(d3))
  K3 <- diag(q3); dimnames(K3) <- list(levels(d3$gid), levels(d3$gid))

  set.seed(8)
  ne <- 5; q4 <- 40
  L4 <- matrix(rnorm(ne), ne, 1)
  u4 <- rnorm(q4) %o% as.numeric(L4) + matrix(rnorm(q4 * ne, 0, .3), q4)
  d4 <- expand.grid(gid = factor(1:q4), env = factor(paste0("E", 1:ne)), rep = 1:2)
  d4$y <- as.numeric(d4$env) + u4[cbind(as.integer(d4$gid), as.integer(d4$env))] + rnorm(nrow(d4))

  set.seed(7); g5 <- expand.grid(row = 1:16, col = 1:10)
  g5$gid <- factor(sample(rep(1:40, 4)))
  g5$trt <- factor(rep(c("a", "b"), length.out = nrow(g5)))
  g5$y <- 5 + 0.03 * (g5$row - 8)^2 - 0.4 * cos(g5$col / 3) + 0.3 * (g5$trt == "b") +
    rnorm(40)[g5$gid] + rnorm(nrow(g5), 0, 0.5)

  # dispositif a deux especes en bacs (termes explicites, moteur par blocs)
  set.seed(21)
  ng <- 30
  dd <- do.call(rbind, lapply(1:3, function(b) data.frame(bac = sprintf("b%d", b),
                                                           col = rep(1:12, 10), rang = rep(1:10, each = 12))))
  dd$esp <- ifelse(dd$col %% 2 == 1, "A", "B")
  dd$gen <- NA_character_
  for (sp in c("A", "B")) { w <- which(dd$esp == sp); dd$gen[w] <- sample(rep_len(paste0(sp, 1:ng), length(w))) }
  dd <- dd[sort(sample(nrow(dd), round(0.92 * nrow(dd)))), ]
  dd <- dd[order(dd$esp != "A"), ]; rownames(dd) <- NULL
  dd$id <- sprintf("P%04d", seq_len(nrow(dd)))
  nb <- rx_neighbourhood(coord = cbind(dd$rang, dd$col), group = dd$esp, block = dd$bac, level = dd$gen,
                         id = dd$id, rank = 2, reach = 1, dilution = 0.5, window = "euclidean",
                         kernel = "exponential", output = "level")
  inc <- function(paire, niv) {
    Z <- nb$level[[paire]]; Zt <- methods::as(Z, "TsparseMatrix")
    Matrix::sparseMatrix(i = match(rownames(Z), dd$id)[Zt@i + 1L], j = match(colnames(Z), niv)[Zt@j + 1L],
                         x = Zt@x, dims = c(nrow(dd), length(niv)), dimnames = list(NULL, niv))
  }
  gA <- paste0("A", 1:ng); gB <- paste0("B", 1:ng)
  ZaD <- Matrix::sparseMatrix(i = which(dd$esp == "A"), j = match(dd$gen[dd$esp == "A"], gA), x = 1,
                              dims = c(nrow(dd), ng), dimnames = list(NULL, gA))
  ZbD <- Matrix::sparseMatrix(i = which(dd$esp == "B"), j = match(dd$gen[dd$esp == "B"], gB), x = 1,
                              dims = c(nrow(dd), ng), dimnames = list(NULL, gB))
  ZaI <- inc("A<-A", gA); ZaX <- inc("B<-A", gA); ZbI <- inc("B<-B", gB); ZbX <- inc("A<-B", gB)
  ua <- matrix(rnorm(ng * 3), ng) %*% chol(matrix(c(1, -.3, -.3, -.3, .5, .2, -.3, .2, .5), 3))
  ub <- matrix(rnorm(ng * 3), ng) %*% chol(matrix(c(1, -.3, -.3, -.3, .5, .2, -.3, .2, .5), 3))
  dd$y <- ifelse(dd$esp == "A", 10, 5) +
    as.numeric(ZaD %*% ua[, 1] + ZaI %*% ua[, 2] + ZaX %*% ua[, 3] + ZbD %*% ub[, 1] + ZbI %*% ub[, 2] +
                 ZbX %*% ub[, 3]) + rnorm(3, 0, .5)[factor(dd$bac)] + rnorm(nrow(dd))
  dd$esp <- factor(dd$esp); dd$bac <- factor(dd$bac); dd$esp_bac <- factor(paste(dd$esp, dd$bac))
  r_ige <- ~ str(~ mm(ZaD) + mm(ZaI) + mm(ZaX), struct = "us", name = "gA") +
    str(~ mm(ZbD) + mm(ZbI) + mm(ZbX), struct = "us", name = "gB")
  r_ige_bloc <- ~ str(~ mm(ZaD) + mm(ZaI) + mm(ZaX), struct = "us", name = "gA") +
    str(~ mm(ZbD) + mm(ZbI) + mm(ZbX), struct = "us", name = "gB") + iid(esp_bac)

  # decalage spatial (jeu de test-decalage.R)
  reine <- function(nr, nc) {
    co <- expand.grid(c = seq_len(nc), r = seq_len(nr))
    D <- unname(as.matrix(stats::dist(co, method = "maximum")))
    Matrix::Matrix((D == 1) * 1, sparse = TRUE)
  }
  set.seed(1)
  nr <- 6; nc <- 5; ngl <- 10; ql <- 2; N <- nr * nc
  Wl <- reine(nr, nc); Psi <- matrix(c(.15, .05, 0, .1), 2)
  gl <- sample(rep(seq_len(ngl), length.out = N))
  al <- matrix(rnorm(ngl * ql), ngl) %*% chol(matrix(c(1, .3, .3, 1), 2))
  el <- matrix(rnorm(N * ql), N) %*% chol(matrix(c(1, -.2, -.2, 1), 2))
  eta <- as.numeric(al[gl, ] + el) + rep(1:ql * 3, each = N)
  yl <- as.numeric(solve(diag(N * ql) - kronecker(Psi, as.matrix(Wl)), eta))
  trl <- factor(rep(c("T1", "T2"), each = N)); unl <- rep(seq_len(N), ql)
  Xl <- stats::model.matrix(~ trl - 1)
  Zg <- Matrix::sparseMatrix(i = seq_len(N), j = gl, x = 1, dims = c(N, ngl))
  z0 <- Matrix::Matrix(0, N, ngl, sparse = TRUE)
  ZDl <- list(T1 = rbind(Zg, z0), T2 = rbind(z0, Zg))

  # --- cas -----------------------------------------------------------------------
  cas("simple_lbfgsb", do.call("rx_reml", c(list(y ~ x, random = ~ gid + iid(bloc), data = d1,
       vpredict = c(h2 = "V1/(V1+V2+V3)"), wald = TRUE, kenward_roger = TRUE, pev = TRUE), fa())))
  f_simple <- readRDS(file.path(dossier, "simple_lbfgsb.rds"))$valeur
  cas("simple_ai", do.call("rx_reml", c(list(y ~ x, random = ~ gid + iid(bloc), data = d1, pev = TRUE,
       vpredict = c(h2 = "V1/(V1+V2+V3)")), fa(optimizer = "ai"))))
  cas("simple_mme_ai", do.call("rx_reml", c(list(y ~ x, random = ~ gid + iid(bloc), data = d1, pev = TRUE),
       fa(engine = "mme", optimizer = "ai"))))
  cas("simple_auto", do.call("rx_reml", c(list(y ~ x, random = ~ gid + iid(bloc), data = d1),
       fa(engine = "auto", optimizer = "auto"))))
  cas("simple_restarts", do.call("rx_reml", c(list(y ~ x, random = ~ gid + iid(bloc), data = d1),
       fa(n_restarts = 2L))))
  options(remlax.server = FALSE)
  cas("simple_sans_serveur", do.call("rx_reml", c(list(y ~ x, random = ~ gid + iid(bloc), data = d1), fa())))
  options(remlax.server = NULL)
  cas("predict_gid", rx_predict(f_simple, classify = "gid"))
  cas("predict_x", rx_predict(f_simple, classify = "x"))
  cas("evaluation_theta", { m <- rx_model(d1$y, cbind(1, d1$x), list(rx_term("gid", d1$gid), rx_term("bloc", d1$bloc)))
       do.call("rx_fit", c(list(m, theta_init = f_simple$theta, maxiter = 0L, polish = 0L), fa())) })
  cas("theta_fixe", { m <- rx_model(d1$y, cbind(1, d1$x), list(rx_term("gid", d1$gid), rx_term("bloc", d1$bloc)))
       th0 <- f_simple$theta; th0[2] <- log(0.5) / 2
       do.call("rx_fit", c(list(m, theta_init = th0, fixed_theta = 2L), fa())) })

  cas("gblup_dense", do.call("rx_reml", c(list(y ~ 1, random = ~ vm(gid, K = K), data = d2, pev = TRUE,
       vpredict = c(h2 = "V1/(V1+V2)")), fa())))
  cas("gblup_mme_ai", do.call("rx_reml", c(list(y ~ 1, random = ~ vm(gid, K = K), data = d2, pev = TRUE),
       fa(engine = "mme", optimizer = "ai"))))
  cas("gblup_mme_creux", do.call("rx_reml", c(list(y ~ 1, random = ~ vm(gid, K = K), data = d2),
       fa(engine = "mme", optimizer = "ai", inversion = "sparse"))))
  cas("gblup_econome", do.call("rx_reml", c(list(y ~ 1, random = ~ vm(gid, K = K), data = d2, pev = TRUE),
       fa(memory = "low"))))
  cas("gblup_econome_ai", do.call("rx_reml", c(list(y ~ 1, random = ~ vm(gid, K = K), data = d2),
       fa(memory = "low", optimizer = "ai"))))

  cas("mt_us_lbfgsb", do.call("rx_reml", c(list(y ~ trait, random = ~ us(gid, K = K3), residual = ~ us(trait):units,
       data = d3, trait = "trait", unit = "unite", pev = TRUE,
       vpredict = c(rg12 = "V2/sqrt(V1*V3)")), fa())))
  cas("mt_us_ai", do.call("rx_reml", c(list(y ~ trait, random = ~ us(gid, K = K3), residual = ~ us(trait):units,
       data = d3, trait = "trait", unit = "unite"), fa(optimizer = "ai"))))
  cas("mt_us_mme", do.call("rx_reml", c(list(y ~ trait, random = ~ us(gid, K = K3), residual = ~ us(trait):units,
       data = d3, trait = "trait", unit = "unite"), fa(engine = "mme", optimizer = "ai"))))
  cas("fa1_ai", do.call("rx_reml", c(list(y ~ env, random = ~ fa(gid, rank = 1), residual = ~ diag(env):units,
       data = d4, trait = "env"), fa(optimizer = "ai"))))
  cas("fa1_mme", do.call("rx_reml", c(list(y ~ env, random = ~ fa(gid, rank = 1), residual = ~ diag(env):units,
       data = d4, trait = "env"), fa(engine = "mme", optimizer = "ai"))))
  cas("diag_lbfgsb", do.call("rx_reml", c(list(y ~ env, random = ~ diag(gid), residual = ~ diag(env):units,
       data = d4, trait = "env"), fa())))

  cas("ar1_champ", do.call("rx_reml", c(list(y ~ trt, random = ~ gid + ar1(row, col), residual = ~ units,
       data = g5), fa())))
  cas("spl2d", do.call("rx_reml", c(list(y ~ trt, random = ~ gid + spl2d(col, row, nseg = c(5, 8)), data = g5,
       pev = TRUE), fa())))
  f_spl <- readRDS(file.path(dossier, "spl2d.rds"))$valeur
  cas("spl2d_sorties", list(dim = rx_dimensions(f_spl), h2 = rx_heritability(f_spl, "gid"),
       tendance = rx_spatial_trend(f_spl, grid = c(20L, 20L))))

  cas("ige_dense", do.call("rx_reml", c(list(y ~ esp_bac, random = r_ige, residual = ~ dsum(~ units | esp),
       data = dd), fa())))
  cas("ige_ai", do.call("rx_reml", c(list(y ~ esp_bac, random = r_ige, residual = ~ dsum(~ units | esp),
       data = dd), fa(optimizer = "ai"))))
  cas("ige_mme", do.call("rx_reml", c(list(y ~ esp_bac, random = r_ige, residual = ~ dsum(~ units | esp),
       data = dd), fa(engine = "mme", optimizer = "ai"))))
  cas("ige_blocs", do.call("rx_reml", c(list(y ~ esp, random = r_ige_bloc, residual = ~ dsum(~ units | esp),
       data = dd), fa(engine = "blocks"))))
  cas("ige_blocs_ref_dense", do.call("rx_reml", c(list(y ~ esp, random = r_ige_bloc, residual = ~ dsum(~ units | esp),
       data = dd), fa())))
  cas("voisinage", list(nb = nb, exposition = rx_exposure(ZaD, ZaI), grm = rx_grm((M + 1)[1:20, 1:50], coding = "count")))

  cas("decalage", { m <- rx_model(yl, Xl, list(rx_term("gen", ZDl, struct = "us")),
       rx_residual("us", trait = trl, unit = unl), lag = rx_lag(Wl, trl, unl))
       do.call("rx_fit", c(list(m), fa())) })
  cas("residuelle_none", { m <- rx_model(yl, Xl, list(rx_term("gen", ZDl, struct = "us"),
       rx_term("plante", list(T1 = rbind(Diagonal(N), Matrix(0, N, N)), T2 = rbind(Matrix(0, N, N), Diagonal(N))),
               struct = "us")), rx_residual("none"))
       do.call("rx_fit", c(list(m), fa())) })

  cas("scan", { fit <- do.call("rx_reml", c(list(y ~ 1, random = ~ vm(gid, K = K), data = d2), fa()))
       m <- rx_model(d2$y, matrix(1, nrow(d2), 1), list(rx_term("gid", d2$gid, K = K)))
       f2 <- do.call("rx_fit", c(list(m), fa()))
       rx_scan(f2, m, M, carte = carte, incidences = c(dir = "gid"), tests = "dir", maf = NULL,
               backend = "cpu", verbose = FALSE) })
  construire <- function(Kk, X_sup = NULL)
    rx_model(d2$y, cbind(rep(1, nrow(d2)), X_sup), list(rx_term("gid", d2$gid, K = Kk)), rx_residual())
  K_de <- function(W) { W <- scale(W, scale = FALSE); Kk <- tcrossprod(W) / ncol(W) + diag(0.01, nrow(W))
    dimnames(Kk) <- list(rownames(W), rownames(W)); Kk }
  cas("scan_loco", rx_scan_loco(construire, K_de, M, carte, incidences = c(dir = "gid"), tests = "dir",
       fit_args = list(backend = "cpu"), verbose = FALSE, backend = "cpu"))
  lo <- readRDS(file.path(dossier, "scan_loco.rds"))$valeur
  cas("scan_exact", rx_scan_exact(lo, "dir", construire, K_de, M, carte, c(dir = "gid"), top = 3,
       fit_args = list(backend = "cpu"), verbose = FALSE))

  if (isTRUE(rx_bayes_available()))
    cas("bayes", do.call("rx_reml", c(list(y ~ 1, random = ~ gid, data = d1[1:160, ], vpredict = c(h2 = "V1/(V1+V2)"),
         bayes = rx_bayes(chains = 4, superchains = 2, warmup = 100, draws = 100, seed = 3)), fa())))

  rx_server_stop()
  writeLines(journal, file.path(dossier, "journal.txt"))
  quit(save = "no")
}

# ------------------------------------------------------------------------------
# Comparaison
# ------------------------------------------------------------------------------
if (mode == "compare") {
  ref <- a0[2]; nouv <- a0[3]; bis <- if (length(a0) >= 4) a0[4] else NULL
  # Aplatit un objet en liste nommee de feuilles : chemin -> valeur
  aplatir <- function(x, chemin = "") {
    if (is.environment(x) || is.function(x)) return(list())
    if (is(x, "Matrix")) x <- as.matrix(x)
    if (is.list(x) && !is.data.frame(x)) {
      noms <- names(x) %||% rep("", length(x))
      noms[noms == ""] <- paste0("[", which(noms == ""), "]")
      return(do.call(c, c(list(list()), lapply(seq_along(x), function(i)
        aplatir(x[[i]], paste0(chemin, "$", noms[i]))))))
    }
    if (is.data.frame(x)) return(do.call(c, c(list(list()), lapply(names(x), function(n)
      aplatir(x[[n]], paste0(chemin, "$", n))))))
    stats::setNames(list(x), chemin)
  }
  `%||%` <- function(a, b) if (is.null(a)) b else a
  ecart <- function(a, b) {
    if (is.null(a) || is.null(b)) return(if (is.null(a) && is.null(b)) 0 else Inf)
    if (is.numeric(a) && is.numeric(b)) {
      if (length(a) != length(b) || !identical(dim(a), dim(b))) return(Inf)
      na <- is.na(a) | is.na(b)
      if (any(is.na(a) != is.na(b))) return(Inf)
      if (all(na)) return(0)
      a <- a[!na]; b <- b[!na]
      inf <- is.infinite(a) | is.infinite(b)
      if (any(inf) && !identical(a[inf], b[inf])) return(Inf)
      a <- a[!inf]; b <- b[!inf]
      if (!length(a)) return(0)
      return(max(abs(a - b) / pmax(abs(a), abs(b), 1e-300)))
    }
    if (identical(a, b)) 0 else Inf
  }
  lire <- function(f) {                # les durees imprimees changent a chaque execution
    x <- readRDS(f)
    if (!is.null(x$affichage)) x$affichage <- gsub("[0-9]+([.,][0-9]+)? ?s\\b", "<duree>", x$affichage)
    x
  }
  volatils <- character()
  if (!is.null(bis)) {                 # champs qui changent entre deux captures du meme code
    for (f in list.files(ref, "\\.rds$")) {
      A <- aplatir(lire(file.path(ref, f))); B <- aplatir(lire(file.path(bis, f)))
      for (k in union(names(A), names(B))) if (ecart(A[[k]], B[[k]]) > 0)
        volatils <- c(volatils, paste(f, k))
    }
    writeLines(volatils, file.path(nouv, "volatils.txt"))
  } else if (file.exists(file.path(ref, "volatils.txt"))) volatils <- readLines(file.path(ref, "volatils.txt"))
  motif_temps <- "secondes|chrono|compile_s|eval_s|\\$dir$|temps|t_optim|elapsed|duree|fils_bornes"
  res <- list()
  for (f in list.files(ref, "\\.rds$")) {
    if (!file.exists(file.path(nouv, f))) { res[[f]] <- data.frame(cas = f, champ = "(absent)", ecart = Inf); next }
    A <- aplatir(lire(file.path(ref, f))); B <- aplatir(lire(file.path(nouv, f)))
    for (k in union(names(A), names(B))) {
      if (paste(f, k) %in% volatils || grepl(motif_temps, k)) next
      e <- ecart(A[[k]], B[[k]])
      if (e > 0) res[[paste(f, k)]] <- data.frame(cas = f, champ = k, ecart = e)
    }
  }
  tab <- if (length(res)) do.call(rbind, res) else data.frame(cas = character(), champ = character(), ecart = numeric())
  rownames(tab) <- NULL
  cat(sprintf("%d fichiers compares, %d champs volatils ignores\n", length(list.files(ref, "\\.rds$")),
              length(volatils)))
  if (!nrow(tab)) cat("AUCUN ECART : resultats identiques au bit pres.\n") else {
    tab <- tab[order(-tab$ecart), ]
    cat(sprintf("%d champs differents ; ecart relatif max %.3g\n", nrow(tab), max(tab$ecart)))
    print(utils::head(tab, 40), row.names = FALSE)
  }
  utils::write.csv(tab, file.path(nouv, "ecarts.csv"), row.names = FALSE)
}
