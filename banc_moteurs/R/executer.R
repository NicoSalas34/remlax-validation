# ==============================================================================
# Un cas du banc : famille, n, repetition. Toutes les voies remlax et les
# logiciels externes sur les MEMES donnees ; arbitre REML dense commun.
#
#   Rscript R/executer.R --famille=grm --n=2000 --rep=1 --sortie=resultats [--pkg=<remlax>]
#            [--voies=D-CPU,D-GPU,M-dense-CPU,M-dense-GPU,M-creux-CPU] [--externes=asreml,sommer,lme4]
#            [--delai=1800]
#
# Chaque ajustement tourne dans un processus R a part (callr), sous delai :
# un logiciel qui ne converge pas ou sature la memoire ne bloque pas le cas.
# Sorties (une ligne par ajustement) : <sortie>/ajustements_<famille>_<n>_<rep>.csv
# et <sortie>/parametres_<famille>_<n>_<rep>.csv.
# ==============================================================================
args <- commandArgs(TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", k, "="), "", a[1]) else d }
fam <- opt("famille", "iid"); n <- as.integer(opt("n", "500")); rep <- as.integer(opt("rep", "1"))
sortie <- opt("sortie", "resultats"); pkg <- opt("pkg", "")
voies <- strsplit(opt("voies", "D-CPU,D-GPU,M-dense-CPU,M-dense-GPU,M-creux-CPU"), ",")[[1]]
externes <- strsplit(opt("externes", "asreml,sommer,lme4"), ",")[[1]]
delai <- as.numeric(opt("delai", "1800"))
ici <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
outils <- normalizePath(file.path(ici, "..", "verif_simulations", "R", "outils.R"))
dir.create(sortie, recursive = TRUE, showWarnings = FALSE)

charger <- function(pkg) {
  if (nzchar(pkg)) suppressMessages(pkgload::load_all(pkg, quiet = TRUE)) else suppressMessages(library(remlax))
}
charger(pkg)
source(outils); source(file.path(ici, "R", "familles.R"))
graine <- graine_de(fam, n, rep)
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

# REPRISE. Chaque ajustement est ecrit des qu'il finit (fichiers .partiel), et
# un temoin en_cours_<cas> nomme celui qui tourne. Si la tache meurt (memoire
# depassee : SLURM tue tout le groupe de processus), la relance garde ce qui
# est fait, note MEMOIRE pour l'ajustement du temoin et fait le reste. Les
# echecs de licence ou d'espace de travail asreml sont refaits.
cle <- sprintf("%s_%d_%d", fam, n, rep)
f_aj <- file.path(sortie, paste0("ajustements_", cle, ".csv"))
f_pa <- file.path(sortie, paste0("parametres_", cle, ".csv"))
f_aj_p <- paste0(f_aj, ".partiel"); f_pa_p <- paste0(f_pa, ".partiel")
f_cours <- file.path(sortie, paste0("en_cours_", cle))
lire_si <- function(a, b) if (file.exists(a)) utils::read.csv(a, stringsAsFactors = FALSE) else
  if (file.exists(b)) utils::read.csv(b, stringsAsFactors = FALSE) else NULL
lignes0 <- lire_si(f_aj, f_aj_p); params0 <- lire_si(f_pa, f_pa_p)
if (!is.null(lignes0)) {
  refaire <- lignes0$statut == "ERREUR" & grepl("License error|licenses in use|workspace", lignes0$message)
  lignes0 <- lignes0[!refaire, , drop = FALSE]
  if (!is.null(params0)) params0 <- params0[paste(params0$logiciel, params0$voie) %in% paste(lignes0$logiciel, lignes0$voie), , drop = FALSE]
}
deja <- if (is.null(lignes0)) character(0) else lignes0$voie
mort <- if (file.exists(f_cours)) readLines(f_cours, warn = FALSE)[1] else NA_character_

t0 <- Sys.time()
cas <- CONSTRUCTEURS[[fam]](n, graine)
t_gen <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
n_obs <- length(cas$y)
cat(sprintf("== %s n=%d (%d obs) rep=%d graine=%d | generation %.1f s\n", fam, n, n_obs, rep, graine, t_gen))

VOIES <- list("D-CPU" = list(m = "dense", backend = "cpu", engine = "dense", inversion = "auto"),
              "D-CPU-LBFGS" = list(m = "dense", backend = "cpu", engine = "dense", inversion = "auto", optimizer = "lbfgsb"),
              "D-GPU-LBFGS" = list(m = "dense", backend = "gpu", engine = "dense", inversion = "auto", optimizer = "lbfgsb"),
              "D-GPU" = list(m = "dense", backend = "gpu", engine = "dense", inversion = "auto"),
              "M-dense-CPU" = list(m = "mme", backend = "cpu", engine = "mme", inversion = "dense"),
              "M-dense-GPU" = list(m = "mme", backend = "gpu", engine = "mme", inversion = "dense"),
              "M-creux-CPU" = list(m = "mme", backend = "cpu", engine = "mme", inversion = "creuse"))

# "--voies=none" : aucune voie remlax (une tache par voie dans banc_iee_long)
voies <- intersect(voies, names(VOIES))

# un ajustement remlax dans un processus a part ; le temps est mesure DANS ce
# processus, autour du seul appel a rx_fit
lancer_remlax <- function(v) {
  model <- cas$modeles[[v$m]]
  callr::r(function(model, v, pkg, outils) {
    if (nzchar(pkg)) suppressMessages(pkgload::load_all(pkg, quiet = TRUE)) else suppressMessages(library(remlax))
    source(outils)
    t0 <- Sys.time()
    f <- rx_fit(model, backend = v$backend, engine = v$engine, inversion = v$inversion,
                optimizer = if (is.null(v$optimizer)) "ai" else v$optimizer,
                hessian = FALSE, blups = FALSE, verbose = FALSE)
    s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    p <- par_de_fit(f)
    list(est = p$est, logLik = f$logLik, conv = isTRUE(as.logical(f$scipy_success)), grad_rel = f$grad_rel %||% NA_real_,
         secondes = s, solveur_s = f$secondes, compile_s = f$compile_s, eval_s = f$eval_s,
         n_iter = f$n_iter, n_eval = f$n_eval, moteur = f$moteur %||% v$engine, inversion = f$inversion %||% "",
         n_par = f$n_par,
         message = f$scipy_message)
  }, args = list(model, v, pkg, outils), timeout = delai)
}
lancer_externe <- function(nom) {
  callr::r(function(cas_fun, n, graine, pkg, outils, familles, nom) {
    if (nzchar(pkg)) suppressMessages(pkgload::load_all(pkg, quiet = TRUE)) else suppressMessages(library(remlax))
    source(outils); source(familles)
    if (nom == "asreml") { suppressMessages(library(asreml)); asreml::asreml.options(trace = FALSE, workspace = "4gb", pworkspace = "1gb") }
    cas <- CONSTRUCTEURS[[cas_fun]](n, graine)
    t0 <- Sys.time()
    r <- cas$externes[[nom]]()
    r$secondes <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    r
  }, args = list(fam, n, graine, pkg, outils, file.path(ici, "R", "familles.R"), nom), timeout = delai)
}

lignes <- if (is.null(lignes0)) list() else list(lignes0)
params <- if (is.null(params0) || !nrow(params0)) list() else list(params0)
ecrire <- function(final = FALSE) {
  a <- if (final) f_aj else f_aj_p; b <- if (final) f_pa else f_pa_p
  if (length(lignes)) utils::write.csv(do.call(rbind, lignes), a, row.names = FALSE)
  if (length(params)) utils::write.csv(do.call(rbind, params), b, row.names = FALSE)
  if (final) unlink(c(f_aj_p, f_pa_p, f_cours))
}
noter <- function(logiciel, voie, r, statut, erreur = "") {
  la <- NA_real_
  if (!is.null(r$est) && all(names(cas$vrai) %in% names(r$est)))
    la <- tryCatch(reml_dense(cas$V_de(r$est[names(cas$vrai)]), cas$y, cas$X)$logLik, error = function(e) NA_real_)
  lignes[[length(lignes) + 1L]] <<- data.frame(
    famille = fam, n = n, n_obs = n_obs, rep = rep, graine = graine, logiciel = logiciel, voie = voie,
    statut = statut, secondes = r$secondes %||% NA_real_, solveur_s = r$solveur_s %||% NA_real_,
    compile_s = r$compile_s %||% NA_real_, eval_s = r$eval_s %||% NA_real_,
    n_iter = r$n_iter %||% NA_integer_, n_eval = r$n_eval %||% NA_integer_,
    n_par = r$n_par %||% NA_integer_,
    moteur = r$moteur %||% "", inversion = r$inversion %||% "", conv = r$conv %||% NA, grad_rel = r$grad_rel %||% NA_real_,
    logLik_natif = r$logLik %||% NA_real_, logLik_arbitre = la,
    message = substr(gsub("[\r\n,]+", " ", paste(r$message %||% "", erreur)), 1, 300), stringsAsFactors = FALSE)
  if (!is.null(r$est)) {
    nm <- intersect(names(cas$vrai), names(r$est))
    params[[length(params) + 1L]] <<- data.frame(famille = fam, n = n, rep = rep, logiciel = logiciel, voie = voie,
                                                 parametre = nm, vrai = unname(cas$vrai[nm]), estime = unname(r$est[nm]),
                                                 stringsAsFactors = FALSE)
  }
  ecrire()
}
# statut d'un echec : delai, memoire (GPU ou processus tue), ou erreur
statut_de <- function(e) {
  m <- conditionMessage(e)
  if (grepl("timeout|timed out", m, ignore.case = TRUE)) "DELAI"
  else if (grepl("RESOURCE_EXHAUSTED|out of memory|ran out of memory|bfc_allocator|cannot allocate|killed", m, ignore.case = TRUE)) "MEMOIRE"
  else "ERREUR"
}
ajuster <- function(logiciel, voie, lancer) {
  if (voie %in% deja) { cat(sprintf("   %-26s deja fait\n", paste(logiciel, voie))); return(invisible()) }
  if (identical(voie, mort)) {
    noter(logiciel, voie, list(), "MEMOIRE", "processus tue par SLURM (memoire de la tache depassee)")
  } else {
    writeLines(voie, f_cours)
    r <- tryCatch(lancer(), error = function(e) e)
    if (inherits(r, "error")) noter(logiciel, voie, list(), statut_de(r), conditionMessage(r))
    else noter(logiciel, voie, r, "OK")
    unlink(f_cours)
  }
  cat(sprintf("   %-26s %s\n", paste(logiciel, voie), utils::tail(lignes, 1)[[1]]$statut))
}
# asreml : licence en ligne a jetons comptes ; on attend qu'un jeton se libere
lancer_asreml <- function() {
  for (k in 1:40) {
    r <- tryCatch(lancer_externe("asreml"), error = function(e) e)
    if (!inherits(r, "error") || !grepl("License error|licenses in use", conditionMessage(r))) break
    Sys.sleep(30)
  }
  if (inherits(r, "error")) stop(r)
  r
}

for (vn in voies) ajuster("remlax", vn, function() lancer_remlax(VOIES[[vn]]))
for (nom in intersect(externes, names(cas$externes)))
  ajuster(nom, nom, if (nom == "asreml") lancer_asreml else function() lancer_externe(nom))
ecrire(final = TRUE)
cat(sprintf("== fin %s en %.1f s\n", cle, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
