# ==============================================================================
# Traduction d'un cas en rx_model remlax, et ajustement selon une voie.
#
# Voies : dense (forme V, memoire econome), blocs (forme V par blocs), mme
# (inversion choisie automatiquement), mme_hybride (separateur dense + blocs).
# Les voies V utilisent L-BFGS-B avec la finition a region de confiance ; les
# voies MME l'AI-REML, comme kmmr.
# ==============================================================================

indicatrices <- function(f, filtre = NULL) {
  f <- factor(f)
  i <- seq_along(f); j <- as.integer(f)
  if (!is.null(filtre)) { i <- i[filtre]; j <- j[filtre] }
  Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(length(f), nlevels(f)),
                       dimnames = list(NULL, levels(f)))
}

# indicatrices d'un facteur sur les seules lignes `ok` (niveaux presents)
sous_indicatrices <- function(f, ok) {
  g <- droplevels(factor(f)[ok])
  Matrix::sparseMatrix(i = which(ok), j = as.integer(g), x = 1, dims = c(length(f), nlevels(g)))
}

# terme G x E : t composantes (une par environnement), q genotypes, parente K
terme_ge <- function(nom, env, gen, K, struct, rang = 0L) {
  env <- factor(env)
  Zl <- lapply(levels(env), function(e) indicatrices(gen, env == e))
  names(Zl) <- levels(env)
  remlax::rx_term(nom, Zl, K = K, struct = struct, rank = rang)
}

modele_remlax <- function(cs) {
  d <- cs$d
  if (cs$type == "deux_etapes") {
    X <- stats::model.matrix(~ -1 + ManEnv, d)
    tm <- terme_ge("ge", d$ManEnv, d$Gen, cs$K, "fa", cs$rang)
    return(remlax::rx_model(d$GY, X, list(tm), remlax::rx_residual("iid"), weights = d$weight))
  }
  if (cs$type == "une_etape") {
    X <- stats::model.matrix(~ 1 + Env, d)
    termes <- list()
    secs <- list()
    for (e in levels(d$Env)) {
      ok <- d$Env == e
      termes[[length(termes) + 1L]] <- remlax::rx_term(paste0("rep_", e), sous_indicatrices(d$Rep, ok), t = 1L)
      termes[[length(termes) + 1L]] <- remlax::rx_term(paste0("blk_", e), sous_indicatrices(d$Block, ok), t = 1L)
      r <- d$Row[ok] - min(d$Row[ok]); cc <- d$Col[ok] - min(d$Col[ok])
      nr <- max(r) + 1L; nc <- max(cc) + 1L
      Zc <- Matrix::sparseMatrix(i = which(ok), j = r * nc + cc + 1L, x = 1, dims = c(nrow(d), nr * nc))
      termes[[length(termes) + 1L]] <- remlax::rx_term(paste0("champ_", e), Zc, t = 1L,
                                                       level = "ar1ar1", dims = c(nr, nc))
      secs[[length(secs) + 1L]] <- remlax::rx_residual("iid", rows = which(ok), name = e)
    }
    termes[[length(termes) + 1L]] <- terme_ge("ge", d$Env, d$Gen, cs$K,
                                              if (cs$ge == "us") "us" else "fa", if (cs$ge == "us") 0L else 1L)
    return(remlax::rx_model(d$GY, X, termes, remlax::rx_residual("iid", sections = secs)))
  }
  if (cs$type == "sim") {
    # indicatrices sur TOUS les niveaux de la parente, observes ou non
    ind <- function(ok) Matrix::sparseMatrix(i = which(ok), j = as.integer(d$Gen)[ok], x = 1,
                                             dims = c(nrow(d), nlevels(d$Gen)))
    if (cs$t == 1L) {
      tm <- remlax::rx_term("g", ind(rep(TRUE, nrow(d))), K = cs$K, t = 1L)
      return(remlax::rx_model(d$GY, matrix(1, nrow(d), 1), list(tm), remlax::rx_residual("iid")))
    }
    Zl <- lapply(levels(d$Env), function(e) ind(d$Env == e)); names(Zl) <- levels(d$Env)
    st <- if (startsWith(cs$ge, "fa")) "fa" else cs$ge
    tm <- remlax::rx_term("ge", Zl, K = cs$K, struct = st, rank = cs$rang)
    secs <- lapply(levels(d$Env), function(e) remlax::rx_residual("iid", rows = which(d$Env == e), name = e))
    return(remlax::rx_model(d$GY, stats::model.matrix(~ -1 + Env, d), list(tm),
                            remlax::rx_residual("iid", sections = secs)))
  }
  if (cs$type == "gblup") {
    X <- stats::model.matrix(~ 1 + Cycle, d)
    tm <- remlax::rx_term("g", indicatrices(d$Gen), K = cs$K, t = 1L)
    return(remlax::rx_model(d$GY, X, list(tm), remlax::rx_residual("iid")))
  }
  X <- stats::model.matrix(~ -1 + TrtCycle, d)
  tm <- terme_ge("ge", d$Trt, d$Gen, cs$K, if (cs$ge == "us") "us" else "fa", if (cs$ge == "us") 0L else 2L)
  secs <- lapply(levels(d$Trt), function(e) remlax::rx_residual("iid", rows = which(d$Trt == e), name = e))
  remlax::rx_model(d$GY, X, list(tm), remlax::rx_residual("iid", sections = secs))
}

ajuster_remlax <- function(cs, voie, backend = "cpu", maxiter = 30000L) {
  t0 <- Sys.time()
  m <- modele_remlax(cs)
  t_modele <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  args <- switch(voie,
    dense       = list(engine = "dense", memory = "low", optimizer = "lbfgsb"),
    # voies entierement fixees du banc absolu (moteur x inversion x optimiseur)
    V_ai        = list(engine = "dense", memory = "standard", optimizer = "ai"),
    V_lbfgsb    = list(engine = "dense", memory = "standard", optimizer = "lbfgsb"),
    mme_dense   = list(engine = "mme", inversion = "dense", optimizer = "ai"),
    blocs       = list(engine = "blocks", optimizer = "lbfgsb"),
    mme         = list(engine = "mme", inversion = "auto", optimizer = "ai"),
    mme_creux   = list(engine = "mme", inversion = "sparse", optimizer = "ai"),
    mme_hybride = list(engine = "mme", inversion = "hybrid", optimizer = "ai"),
    stop("voie inconnue : ", voie))
  commun <- list(model = m, backend = backend, hessian = FALSE, blups = FALSE, verbose = TRUE,
                 maxiter = maxiter)
  fixe <- NULL
  if (!is.null(m$weights)) {
    # residuelle fixee a 1 (poids = inverses des variances des BLUE) : theta_init
    # porte NA partout sauf la residuelle (log sd 0) ; les NA prennent le depart
    # automatique du solveur, sans appel preliminaire.
    th <- rep(NA_real_, m$n_par); th[m$n_par] <- 0
    commun$theta_init <- th
    fixe <- m$n_par
  }
  t1 <- Sys.time()
  r <- do.call(remlax::rx_fit, c(commun, list(fixed_theta = fixe), args))
  sec <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  list(logLik = r$logLik_asreml %||% r$logLik, logLik_remlax = r$logLik, n_iter = r$n_iter %||% NA,
       n_eval = r$n_eval %||% NA, converge = isTRUE(as.logical(r[["scipy_success"]])),
       secondes = sec, secondes_solveur = r$secondes %||% NA, n_par = m$n_par,
       moteur = r$moteur %||% NA, inversion = r$inversion %||% NA,
       t_modele = t_modele, compile_s = r$compile_s %||% NA, eval_s = r$eval_s %||% NA,
       chrono = r[["chrono"]])
}
