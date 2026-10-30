#!/usr/bin/env Rscript
# ==============================================================================
# Synthese du banc remlax / kmmr / asreml-R : tableaux Markdown.
#
#   Rscript R/synthese.R --res <dossier des CSV genotoul> --asreml <banc_asreml.csv>
#
# Temps compares :
#   kmmr   : appel kmmr() + inversion de K faite dans R avant l'appel
#            (+ valeurs de depart pour les cas en deux etapes) ;
#   remlax : appel rx_fit() entier, demarrage du solveur Python compris ;
#            "a chaud" = second ajustement, programmes compiles relus sur disque ;
#   asreml : appel asreml() seul (inversion de K comprise), noeuds de meso.
# Pour remlax, on prend pour chaque tache la version la plus recente disponible
# parmi les fichiers donnes dans l'ordre de --versions.
# ==============================================================================
a <- commandArgs(TRUE)
arg <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1L] }
res <- arg("res", ".")
versions <- strsplit(arg("versions", "v5,v3"), ",")[[1]]
f_asr <- arg("asreml", NULL)

lire <- function(f) if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
chrono_tot <- function(ch, cas_, outil_, voie_, etiq) {
  if (is.null(ch)) return(NA_real_)
  x <- ch[ch$cas == cas_ & ch$outil == outil_ & ch$voie == voie_ & ch$etiquette == etiq, ]
  x
}

cas_ordre <- c("drops2_fa1", "drops2_fa3", "drops1_fa1", "drops1_us",
               "krause_gblup", "krause_fa2", "krause_us")

tableau <- function(cible) {
  # kmmr : le fichier v3 (meme noeuds, memes conditions)
  km <- NULL; kch <- NULL
  for (v in c("v3", "v2", "")) {
    f <- file.path(res, sprintf("banc_%s%s.csv", cible, if (nzchar(v)) paste0("_", v) else ""))
    d <- lire(f)
    if (!is.null(d)) {
      d <- d[d$outil == "kmmr" & d$voie == "ai_tol0.000001", ]
      ch <- lire(sub("\\.csv$", "_chrono.csv", f))
      d <- d[!(d$cas %in% km$cas), ]
      km <- rbind(km, d); kch <- rbind(kch, ch)
    }
  }
  rx <- NULL; rch <- NULL
  for (v in versions) {
    f <- file.path(res, sprintf("banc_%s_%s.csv", cible, v))
    d <- lire(f)
    if (is.null(d)) next
    d <- d[d$outil == "remlax" & d$voie == "mme" & d$etiquette == "chaud", ]
    d <- d[!(d$cas %in% rx$cas), ]
    if (nrow(d)) d$version <- v
    rx <- rbind(rx, d); rch <- rbind(rch, lire(sub("\\.csv$", "_chrono.csv", f)))
  }
  lignes <- c(sprintf("| cas | kmmr logLik | kmmr it. | kmmr s | remlax logLik | remlax it. | remlax s | remlax s/eval | version |"),
              "|---|---|---|---|---|---|---|---|---|")
  for (cs in cas_ordre) {
    k <- km[km$cas == cs, ][1, ]
    r <- rx[rx$cas == cs, ][1, ]
    kt <- NA
    if (!is.na(k$secondes)) {
      x <- kch[kch$cas == cs & kch$outil == "kmmr" & kch$voie == "ai_tol0.000001", ]
      extra <- sum(x$secondes[x$etape %in% c("inverse_K", "valeurs_depart")])
      kt <- k$secondes + extra
    }
    lignes <- c(lignes, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s |", cs,
      if (is.na(k$logLik)) "" else sprintf("%.4f", k$logLik), if (is.na(k$n_iter)) "" else k$n_iter,
      if (is.na(kt)) "" else sprintf("%.1f", kt),
      if (is.null(r) || is.na(r$logLik)) "" else sprintf("%.4f", r$logLik),
      if (is.null(r) || is.na(r$n_iter)) "" else r$n_iter,
      if (is.null(r) || is.na(r$secondes)) "" else sprintf("%.1f", r$secondes),
      if (is.null(r) || is.na(r$eval_s)) "" else sprintf("%.3f", r$eval_s),
      if (is.null(r) || is.null(r$version)) "" else r$version))
  }
  lignes
}

for (cible in c("cpu", "gpu")) {
  cat(sprintf("\n### %s\n\n", toupper(cible)))
  cat(tableau(cible), sep = "\n")
}

if (!is.null(f_asr) && file.exists(f_asr)) {
  d <- lire(f_asr)
  d <- d[!is.na(d$logLik), ]
  cat("\n### asreml-R (CPU, meso, 8 coeurs)\n\n| cas | logLik | it. | s | converge |\n|---|---|---|---|---|\n")
  for (cs in cas_ordre) {
    x <- d[d$cas == cs, ]
    if (!nrow(x)) { cat(sprintf("| %s | | | | |\n", cs)); next }
    x <- x[nrow(x), ]
    cat(sprintf("| %s | %.4f | %s | %.1f | %s |\n", cs, x$logLik, x$n_iter, x$secondes, x$converge))
  }
}

# phases remlax (a chaud), CPU
for (cible in c("cpu", "gpu")) {
  for (v in versions) {
    f <- file.path(res, sprintf("banc_%s_%s_chrono.csv", cible, v))
    ch <- lire(f)
    if (is.null(ch)) next
    ch <- ch[ch$outil == "remlax" & ch$voie == "mme" & ch$etiquette == "chaud", ]
    if (!nrow(ch)) next
    et <- c("demarrage_serveur", "preparation", "compilation", "optimisation", "finition",
            "effets", "inference", "export_R")
    cat(sprintf("\n### Phases remlax a chaud (%s, %s), secondes\n\n| cas | %s |\n|%s\n", toupper(cible), v,
                paste(et, collapse = " | "), paste(rep("---|", length(et) + 1), collapse = "")))
    for (cs in intersect(cas_ordre, unique(ch$cas))) {
      x <- ch[ch$cas == cs, ]
      cat(sprintf("| %s | %s |\n", cs, paste(sprintf("%.2f", vapply(et, function(e) sum(x$secondes[x$etape == e]), 0)), collapse = " | ")))
    }
  }
}
