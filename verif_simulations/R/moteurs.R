# ==============================================================================
# Execution d'un cas : ajustements, arbitre, comparaisons, verdicts.
# ==============================================================================

SEUIL_ARBITRE <- 1e-8     # ecart relatif de logLik a l'arbitre
SEUIL_MOTEURS <- 1e-4     # ecart absolu de logLik a la meilleure des voies
SEUIL_Z <- 3              # |estime - vrai| / se

.chr <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(x)[1]
.num <- function(x) if (is.null(x) || !length(x)) NA_real_ else suppressWarnings(as.numeric(x)[1])
.lgl <- function(x) if (is.null(x) || !length(x)) NA else as.logical(x)[1]

ajuster <- function(spec, model, dir) {
  t0 <- Sys.time()
  out <- switch(spec$type,
    rx_fit = {
      fit <- rx_fit(model, backend = "cpu", verbose = FALSE, optimizer = spec$optimizer, dir = dir)
      # Verdict du paquet : decrement de Newton ou gradient relatif projete,
      # apres polissage. scipy_success ne porte que sur la phase L-BFGS-B / AI
      # avant polissage ; il est rapporte dans le detail, sans trancher.
      conv <- isTRUE(as.logical(fit$conv_decrement)) || isTRUE(as.logical(fit$conv_grad_rel))
      list(par = par_de_fit(fit), logLik = fit$logLik, conv = conv,
           detail = sprintf("succes=%s decrement=%s grad_rel=%s hessien_ok=%s bornes=%s iter=%s",
                            .lgl(fit$scipy_success), .lgl(fit$conv_decrement), .lgl(fit$conv_grad_rel),
                            .lgl(fit$conv_hessien_ok), .num(fit$n_at_bound), .num(fit$n_iter)),
           newton_decrement = .num(fit$newton_decrement), grad_rel = .num(fit$grad_rel),
           hessien_ok = .lgl(fit$conv_hessien_ok), n_at_bound = .num(fit$n_at_bound),
           message = .chr(fit$scipy_message))
    },
    # moteur MME : AI-REML (inversion automatique ou creuse) ou L-BFGS-B en
    # inversion creuse, la voie qui remplace l'ancien moteur RTMB.
    mme = {
      fit <- rx_fit(model, backend = "cpu", verbose = FALSE, optimizer = spec$optimizer,
                    engine = "mme", inversion = spec$inversion, dir = dir)
      conv <- isTRUE(as.logical(fit$conv_decrement)) || isTRUE(as.logical(fit$conv_grad_rel))
      list(par = par_de_fit(fit), logLik = fit$logLik, conv = conv,
           detail = sprintf("moteur=%s inversion=%s succes=%s decrement=%s grad_rel=%s bornes=%s iter=%s",
                            .chr(fit$moteur), .chr(fit$inversion), .lgl(fit$scipy_success),
                            .lgl(fit$conv_decrement), .lgl(fit$conv_grad_rel), .num(fit$n_at_bound),
                            .num(fit$n_iter)),
           newton_decrement = .num(fit$newton_decrement), grad_rel = .num(fit$grad_rel),
           hessien_ok = .lgl(fit$conv_hessien_ok), n_at_bound = .num(fit$n_at_bound),
           message = .chr(fit$scipy_message))
    })
  out$secondes <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out
}

arbitrer <- function(entree, par) {
  m <- entree$model
  V <- entree$V_de(par)
  if (!is.null(entree$M_de)) reml_lag_dense(V, m$y, as.matrix(m$X), entree$M_de(par))$logLik
  else reml_dense(V, m$y, as.matrix(m$X))$logLik
}

executer_cas <- function(nom, constructeur, graine, tmp) {
  t_cas <- Sys.time()
  cas <- tryCatch(constructeur(graine), error = function(e) e)
  if (inherits(cas, "error"))
    return(list(fits = data.frame(cas = nom, moteur = "construction", modele = NA, n = NA, n_par = NA,
                                  logLik = NA, logLik_arbitre = NA, ecart_rel = NA,
                                  ecart_entre_moteurs = NA, convergence = NA, detail_convergence = NA,
                                  message_optimiseur = NA,
                                  n_z_hors = NA, max_abs_z = NA, secondes = NA,
                                  statut = "ERREUR", raisons = paste("construction :", conditionMessage(cas)),
                                  stringsAsFactors = FALSE),
                comps = NULL, note = NA, famille = NA))
  res <- list()
  filtre <- if (exists("MOTEURS_FILTRE")) MOTEURS_FILTRE else "."
  for (k in seq_along(cas$moteurs)) {
    spec <- cas$moteurs[[k]]
    if (!grepl(filtre, spec$id)) next
    entree <- cas$modeles[[spec$modele]]
    lab <- if (spec$modele == "principal") spec$id else paste0(spec$id, "@", spec$modele)
    dir <- file.path(tmp, paste0(nom, "_g", graine, "_", k))
    a <- tryCatch(ajuster(spec, entree$model, dir), error = function(e) e)
    unlink(dir, recursive = TRUE)
    if (inherits(a, "error")) {
      msg <- gsub("[\r\n]+", " ", conditionMessage(a))
      attendu <- !is.null(spec$erreur_attendue) && grepl(spec$erreur_attendue, msg)
      res[[lab]] <- list(erreur = msg, attendu = attendu, spec = spec)
      next
    }
    a$arbitre <- tryCatch(arbitrer(entree, a$par), error = function(e) structure(NA_real_, msg = conditionMessage(e)))
    a$spec <- spec; a$model <- entree$model
    res[[lab]] <- a
  }
  ok <- vapply(res, function(r) is.null(r$erreur), TRUE)
  best <- if (any(ok)) max(vapply(res[ok], `[[`, 0, "logLik")) else NA_real_
  # Erreurs-types de repli (moteurs sans SE) : celles de lbfgsb sur le meme modele.
  se_ref <- function(modele) {
    for (r in res) if (is.null(r$erreur) && r$spec$type == "rx_fit" && r$spec$modele == modele) return(r$par$se)
    NULL
  }
  lignes <- list(); comps <- list()
  for (lab in names(res)) {
    r <- res[[lab]]
    if (!is.null(r$erreur)) {
      lignes[[lab]] <- data.frame(cas = nom, moteur = lab, modele = r$spec$modele, n = NA, n_par = NA,
                                  logLik = NA, logLik_arbitre = NA, ecart_rel = NA,
                                  ecart_entre_moteurs = NA, convergence = NA, detail_convergence = NA,
                                  message_optimiseur = NA,
                                  n_z_hors = NA, max_abs_z = NA, secondes = NA,
                                  statut = if (r$attendu) "REFUS_ATTENDU" else "ERREUR",
                                  raisons = substr(r$erreur, 1, 300), stringsAsFactors = FALSE)
      next
    }
    arb <- as.numeric(r$arbitre)
    er <- (r$logLik - arb) / abs(arb)
    de <- best - r$logLik
    # composantes
    est <- r$par$est; se <- r$par$se
    ren <- r$spec$renommer
    if (!is.null(ren)) { i <- names(est) %in% names(ren); names(est)[i] <- ren[names(est)[i]]
                         names(se) <- names(est) }
    if (!is.null(cas$derive)) {
      dv <- cas$derive(r$par, r$spec$type)
      if (!is.null(dv)) { est <- c(est, stats::setNames(dv$est, dv$nom)); se <- c(se, stats::setNames(dv$se, dv$nom)) }
    }
    sr <- se_ref(r$spec$modele)
    vr <- cas$vrai
    vr <- vr[names(vr) %in% names(est)]
    e_v <- est[names(vr)]; s_v <- se[names(vr)]
    src <- ifelse(is.finite(s_v), "propre", NA_character_)
    if (!is.null(sr)) { j <- !is.finite(s_v) & names(vr) %in% names(sr)
      s_v[j] <- sr[names(vr)[j]]; src[j & is.finite(s_v)] <- "lbfgsb" }
    z <- (e_v - vr) / s_v
    # composante vraie sur la borne : |z| rapporte, sans trancher
    sur_bord <- names(vr) %in% (cas$bord %||% character(0))
    src[sur_bord] <- paste0(ifelse(is.na(src[sur_bord]), "", src[sur_bord]), "+bord")
    comps[[lab]] <- data.frame(cas = nom, moteur = lab, composante = names(vr), vrai = unname(vr),
                               estime = unname(e_v), se = unname(s_v), se_source = src, z = unname(z),
                               stringsAsFactors = FALSE)
    manque <- setdiff(names(cas$vrai), names(est))
    raisons <- c()
    if (!is.finite(arb)) raisons <- c(raisons, paste("arbitre indisponible :", attr(r$arbitre, "msg") %||% "?"))
    else if (abs(er) >= SEUIL_ARBITRE) raisons <- c(raisons, sprintf("logLik != arbitre (ecart rel %.2e)", er))
    if (is.finite(de) && de >= SEUIL_MOTEURS) raisons <- c(raisons, sprintf("optimum inferieur de %.2e", de))
    if (!isTRUE(r$conv)) raisons <- c(raisons, "convergence non etablie")
    zf <- z[is.finite(z) & !sur_bord]
    nz <- sum(abs(zf) >= SEUIL_Z)
    if (nz) raisons <- c(raisons, sprintf("%d composante(s) a |z| >= %g", nz, SEUIL_Z))
    if (length(manque) && r$spec$type == "rx_fit" && r$spec$modele == "principal")
      raisons <- c(raisons, paste("composantes vraies absentes :", paste(manque, collapse = ", ")))
    lignes[[lab]] <- data.frame(cas = nom, moteur = lab, modele = r$spec$modele, n = length(r$model$y),
                                n_par = r$model$n_par, logLik = r$logLik, logLik_arbitre = arb,
                                ecart_rel = er, ecart_entre_moteurs = de, convergence = r$conv,
                                detail_convergence = r$detail, message_optimiseur = r$message,
                                n_z_hors = nz,
                                max_abs_z = if (length(zf)) max(abs(zf)) else NA_real_,
                                secondes = r$secondes,
                                statut = if (length(raisons)) "ECHEC" else "PASSE",
                                raisons = paste(raisons, collapse = " ; "), stringsAsFactors = FALSE)
  }
  list(fits = do.call(rbind, lignes), comps = do.call(rbind, comps), note = cas$note,
       famille = cas$famille, secondes = as.numeric(difftime(Sys.time(), t_cas, units = "secs")))
}
