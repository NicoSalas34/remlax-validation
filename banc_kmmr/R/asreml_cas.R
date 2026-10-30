# ==============================================================================
# Ajustement d'un cas par asreml-R 4.2, meme modele que remlax et kmmr.
#
#   deux_etapes : GY ~ -1 + ManEnv, fa(ManEnv, k):vm(Gen, Ki), poids,
#                 residuelle fixee a 1 par asr_gaussian(dispersion = 1)
#                 (manuel asreml-R 4.2, section 3.10, poids absolus)
#   une_etape   : GY ~ 1 + Env ; par essai : at(Env):Rep, at(Env):Block, un
#                 champ ar1 x ar1 ALEATOIRE (sa variance, ses deux rho) et une
#                 residuelle iid propre, dsum(~ units | Env) ; G x E en
#                 fa(Env, 1) ou us(Env) x vm(Gen, Ki)
#   gblup       : GY ~ 1 + Cycle, vm(Gen, Ki), residuelle iid
#   multi_trt   : GY ~ -1 + TrtCycle, fa(Trt, 2) ou us(Trt) x vm(Gen, Ki),
#                 dsum(~ units | Trt)
# Aucune colonne ne s'appelle `trait` (mot reserve d'asreml).
# ==============================================================================

# inverse dense de K, au format accepte par vm() : matrice avec dimnames et
# attribut INVERSE (manuel, section 4.3.6). asreml calcule lui-meme log|K|.
inverse_asreml <- function(K) {
  Ki <- solve(K)
  dimnames(Ki) <- dimnames(K)
  attr(Ki, "INVERSE") <- TRUE
  Ki
}

# Champ ar1 x ar1 d'un essai : rang et colonne renumerotes a partir de 1 dans
# l'essai, niveaux 1..n complets (ar1 suppose des niveaux contigus et
# equidistants, comme la grille nr x nc de remlax). Hors de l'essai, le niveau
# 1 sert de bouche-trou : at(Env, e) annule ces lignes.
grille_essai <- function(d, e) {
  ok <- d$Env == e
  r <- d$Row - min(d$Row[ok]) + 1L; cc <- d$Col - min(d$Col[ok]) + 1L
  nr <- max(r[ok]); nc <- max(cc[ok])
  r[!ok] <- 1L; cc[!ok] <- 1L
  list(R = factor(r, levels = seq_len(nr)), C = factor(cc, levels = seq_len(nc)))
}

# Trace d'asreml : une ligne par iteration, heure murale a la seconde. On la
# capture, on la reimprime pour le journal, et on en tire les instants.
asreml_trace <- function(expr) {
  sortie <- utils::capture.output(val <- force(expr), type = "output")
  cat(sortie, sep = "\n")
  # ligne d'iteration : numero, logLik, ..., hh:mm:ss, parfois "( 8 restrained)"
  re_mur <- "[0-9]{1,2}:[0-9]{2}:[0-9]{2}"
  it <- grep(paste0("^ *[0-9]+ +-?[0-9.]+ .* ", re_mur), sortie, value = TRUE)
  hms <- function(x) { p <- as.numeric(strsplit(x, ":")[[1]]); sum(p * c(3600, 60, 1)) }
  mur <- vapply(regmatches(it, regexpr(re_mur, it)), hms, 0, USE.NAMES = FALSE)
  list(val = val, mur = mur)
}

ajuster_asreml <- function(cs, maxit = 50L, workspace = NULL) {
  suppressPackageStartupMessages(library(asreml))
  if (is.null(workspace)) workspace <- if (cs$type %in% c("gblup", "multi_trt", "sim")) "30gb" else "1gb"
  asreml.options(workspace = workspace, pworkspace = workspace, maxit = maxit, trace = TRUE)
  # fils : OMP_NUM_THREADS s'il est pose (celui que note executer.R), sinon tous
  fils <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", "")))
  if (!is.na(fils)) asreml.options(threads = fils)
  d <- cs$d
  # dsum exige des donnees groupees par section ; l'ordre des lignes ne change
  # pas la vraisemblance
  if (cs$type == "une_etape") d <- d[order(d$Env), ]
  if (cs$type == "multi_trt") d <- d[order(d$Trt), ]

  tK <- Sys.time()
  Ki <- inverse_asreml(cs$K)
  t_Ki <- as.numeric(difftime(Sys.time(), tK, units = "secs"))

  tp <- Sys.time()
  args <- list(data = d)
  if (cs$type == "deux_etapes") {
    args$fixed <- GY ~ -1 + ManEnv
    args$random <- stats::as.formula(sprintf("~ fa(ManEnv, %d):vm(Gen, Ki)", cs$rang))
    args$residual <- ~ units
    args$weights <- "weight"
    args$family <- asreml::asr_gaussian(dispersion = 1)
  } else if (cs$type == "une_etape") {
    champs <- character(0)
    for (e in levels(d$Env)) {
      g <- grille_essai(d, e)
      nR <- paste0("R_", e); nC <- paste0("C_", e)
      d[[nR]] <- g$R; d[[nC]] <- g$C
      champs <- c(champs, sprintf("at(Env, '%s'):ar1(%s):ar1(%s)", e, nR, nC))
    }
    ge <- if (cs$ge == "us") "us(Env)" else "fa(Env, 1)"
    args$data <- d
    args$fixed <- GY ~ 1 + Env
    args$random <- stats::as.formula(paste("~ at(Env):Rep + at(Env):Block +",
                                           paste(champs, collapse = " + "),
                                           sprintf("+ %s:vm(Gen, Ki)", ge)))
    args$residual <- ~ dsum(~ units | Env)
  } else if (cs$type == "sim" && cs$t == 1L) {
    args$fixed <- GY ~ 1
    args$random <- ~ vm(Gen, Ki)
    args$residual <- ~ units
  } else if (cs$type == "sim") {
    ge <- switch(cs$ge, us = "us(Env)", diag = "diag(Env)", sprintf("fa(Env, %d)", cs$rang))
    args$fixed <- GY ~ -1 + Env
    args$random <- stats::as.formula(sprintf("~ %s:vm(Gen, Ki)", ge))
    args$residual <- ~ dsum(~ units | Env)
  } else if (cs$type == "gblup") {
    args$fixed <- GY ~ 1 + Cycle
    args$random <- ~ vm(Gen, Ki)
    args$residual <- ~ units
  } else {
    ge <- if (cs$ge == "us") "us(Trt)" else "fa(Trt, 2)"
    args$fixed <- GY ~ -1 + TrtCycle
    args$random <- stats::as.formula(sprintf("~ %s:vm(Gen, Ki)", ge))
    args$residual <- ~ dsum(~ units | Trt)
  }
  t_prep <- as.numeric(difftime(Sys.time(), tp, units = "secs"))
  # les formules sont evaluees ici : Ki doit etre visible dans ce cadre
  environment(args$random) <- environment()

  t0 <- Sys.time()
  tr <- asreml_trace(do.call(asreml::asreml, args))
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fit <- tr$val

  # instants muraux (a la seconde) : debut de l'appel, puis fin de chaque iteration
  m0 <- with(as.POSIXlt(t0), hour * 3600 + min * 60 + sec)
  mur <- tr$mur
  if (length(mur)) mur <- ifelse(mur < m0 - 1, mur + 86400, mur)   # passage de minuit
  n_it <- if (length(mur)) length(mur) else max(ncol(fit$trace) - 1L, NA_integer_)
  # empan de la trace sous 5 s : la seconde ne resout plus rien, etapes en NA
  if (length(mur) < 2 || mur[length(mur)] - mur[1] < 5) mur <- numeric(0)
  # cout d'une iteration : mediane des ecarts entre lignes si elles durent au
  # moins 5 s, sinon moyenne sur l'empan 2..n (la seconde est trop grossiere)
  eval_s <- NA_real_
  if (length(mur) > 1) {
    di <- diff(mur)
    eval_s <- if (stats::median(di) >= 5) stats::median(di) else (mur[length(mur)] - mur[1]) / (length(mur) - 1)
  }
  chrono <- c(inverse_K = t_Ki, preparation = t_prep, ajustement = sec,
              premiere_iteration = if (length(mur)) max(mur[1] - m0, 0) else NA_real_,
              iterations_suivantes = if (length(mur) > 1) mur[length(mur)] - mur[1] else NA_real_,
              hors_iterations = if (length(mur)) max(sec - (mur[length(mur)] - m0), 0) else NA_real_)
  list(chrono = chrono, eval_s = eval_s, logLik = fit$loglik, n_iter = n_it,
       converge = isTRUE(fit$converge), secondes = sec,
       n_par = sum(fit$vparameters.con != "F"))
}
