#!/usr/bin/env Rscript
# ==============================================================================
# Verification par simulation de remlax : logLik contre un arbitre REML dense
# ecrit en R, accord des voies d'ajustement (L-BFGS-B, AI-REML, equations du
# modele mixte en inversion automatique ou creuse) et recuperation des
# composantes simulees.
#
#   Rscript verif_simulations/verif_simulations.R [--pkg=/chemin/remlax] [--cas=motif]
#           [--coeurs=10] [--graine=20261021] [--sortie=dossier] [--rep=1]
#
# --pkg     sources de remlax chargees par pkgload::load_all ; sans --pkg, le paquet
#           installe est charge par R/charger_remlax.R a la racine du depot
# --cas     expression reguliere sur les noms de cas (defaut : tous) ; --liste les affiche
# --coeurs  processus paralleles (chacun limite a 2 fils JAX/BLAS)
# --rep     replicats par cas (graines successives)
# --moteurs expression reguliere sur les voies (lbfgsb, ai, mme, mme_creuse,
#           mme_lbfgsb ; defaut : toutes)
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
opt <- function(nom, defaut) {
  a <- grep(paste0("^--", nom, "="), args, value = TRUE)
  if (length(a)) sub(paste0("^--", nom, "="), "", a[1]) else defaut
}
ici <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})

# Les fils BLAS de R se fixent au demarrage : on relance le script avec un fil
# par processus si l'appelant ne l'a pas fait.
if (!nzchar(Sys.getenv("OPENBLAS_NUM_THREADS"))) {
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "2", MKL_NUM_THREADS = "1")
  st <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(ici, "verif_simulations.R")), shQuote(args)))
  quit(save = "no", status = st)
}

pkg     <- opt("pkg", "")
motif   <- opt("cas", ".")
coeurs  <- as.integer(opt("coeurs", "10"))
graine0 <- as.integer(opt("graine", "20261021"))
nrep    <- as.integer(opt("rep", "1"))
MOTEURS_FILTRE <- opt("moteurs", ".")
sortie  <- opt("sortie", file.path(ici, "resultats"))

Sys.setenv(XLA_FLAGS = "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2",
           JAX_PLATFORMS = "cpu")
if (nzchar(pkg)) {
  suppressMessages(pkgload::load_all(pkg, quiet = TRUE, export_all = TRUE))
} else {
  source(file.path(ici, "..", "R", "charger_remlax.R"))
  pkg <- sprintf("remlax %s installe dans %s", utils::packageVersion("remlax"),
                 dirname(system.file(package = "remlax")))
}
options(remlax.threads = 2L, warn = 1)
for (f in c("outils.R", "cas.R", "cas_limites.R", "moteurs.R")) source(file.path(ici, "R", f))

if ("--liste" %in% args) { cat(names(CATALOGUE), sep = "\n"); quit(save = "no") }
noms <- grep(motif, names(CATALOGUE), value = TRUE)
if (!length(noms)) stop("aucun cas ne correspond a --cas=", motif)
# les cas lourds d'abord, pour equilibrer les processus
lourds <- c("ige_melange_deux_especes", "mt_g_fa2", "decalage_rx_lag", "mt_residuelle_fa1", "mt_g_rr2",
            "mt_g_ante1", "mt_g_chol1", "mt_g_fa1", "spline_spl2d", "niveaux_corg_ge",
            "limite_melange_iee_nul_us_rang2", "limite_us8_presque_rang2", "limite_fa2_vrai_rang1",
            "limite_deux_us_au_bord")
noms <- c(intersect(lourds, noms), setdiff(noms, lourds))
taches <- expand.grid(rep = seq_len(nrep), cas = noms, stringsAsFactors = FALSE)

dir.create(sortie, recursive = TRUE, showWarnings = FALSE)
tmp <- tempfile("verif_"); dir.create(tmp)
py <- rx_python_check(quiet = TRUE)
cat(sprintf("remlax charge depuis %s\nPython : %s (jax %s)\n%d cas x %d replicat(s), %d processus\n",
            pkg, py$versions[["executable"]], py$versions[["jax"]], length(noms), nrep, coeurs))

t0 <- Sys.time()
sorties <- parallel::mclapply(seq_len(nrow(taches)), function(i) {
  nom <- taches$cas[i]; r <- taches$rep[i]
  graine <- graine0 + 1000L * r + match(nom, names(CATALOGUE))
  o <- tryCatch(executer_cas(nom, CATALOGUE[[nom]], graine, tmp), error = function(e)
    list(fits = data.frame(cas = nom, moteur = "execution", statut = "ERREUR",
                           raisons = conditionMessage(e), stringsAsFactors = FALSE)))
  o$fits$replicat <- r; o$fits$graine <- graine
  if (!is.null(o$comps)) { o$comps$replicat <- r; o$comps$graine <- graine }
  message(sprintf("[%s] %-32s rep %d : %s", format(Sys.time(), "%H:%M:%S"), nom, r,
                  paste(sprintf("%s=%s", o$fits$moteur, o$fits$statut), collapse = " ")))
  o
}, mc.cores = coeurs, mc.preschedule = FALSE)
unlink(tmp, recursive = TRUE)

plat <- function(l, k) { x <- lapply(l, function(o) if (is.list(o)) o[[k]] else NULL)
  x <- x[!vapply(x, is.null, TRUE)]
  if (!length(x)) return(NULL)
  cols <- unique(unlist(lapply(x, names)))
  do.call(rbind, lapply(x, function(d) { for (c in setdiff(cols, names(d))) d[[c]] <- NA; d[cols] })) }
fits <- plat(sorties, "fits"); comps <- plat(sorties, "comps")
infos <- do.call(rbind, lapply(seq_along(sorties), function(i) {
  o <- sorties[[i]]
  data.frame(cas = taches$cas[i], replicat = taches$rep[i], famille = .chr(o$famille),
             note = .chr(o$note), secondes = .num(o$secondes), stringsAsFactors = FALSE)
}))
utils::write.csv(fits, file.path(sortie, "ajustements.csv"), row.names = FALSE)
if (!is.null(comps)) utils::write.csv(comps, file.path(sortie, "composantes.csv"), row.names = FALSE)
utils::write.csv(infos, file.path(sortie, "cas.csv"), row.names = FALSE)

# --- resume texte -----------------------------------------------------------------
L <- c(sprintf("Verification par simulation de remlax (%s)", format(Sys.time(), "%Y-%m-%d %H:%M")),
       sprintf("Paquet : %s", pkg),
       sprintf("Criteres : |logLik - arbitre| / |arbitre| < %g ; ecart a la meilleure voie < %g ; |z| < %g",
               SEUIL_ARBITRE, SEUIL_MOTEURS, SEUIL_Z),
       sprintf("Duree totale : %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))), "")
tab <- table(fits$statut)
L <- c(L, paste("Ajustements :", paste(sprintf("%s %d", names(tab), tab), collapse = ", ")), "")
for (nm in unique(fits$cas)) {
  f <- fits[fits$cas == nm, ]
  inf <- infos[infos$cas == nm, ][1, ]
  L <- c(L, sprintf("== %s  [%s]  %s", nm, inf$famille, inf$note))
  for (i in seq_len(nrow(f))) {
    L <- c(L, sprintf("   %-18s rep %d  %-13s logLik %s  arbitre %s  ecart_rel %s  ecart_voies %s  max|z| %s  %s",
                      f$moteur[i], f$replicat[i], f$statut[i],
                      formatC(f$logLik[i], format = "f", digits = 6), formatC(f$logLik_arbitre[i], format = "f", digits = 6),
                      formatC(f$ecart_rel[i], format = "e", digits = 1),
                      formatC(f$ecart_entre_moteurs[i], format = "e", digits = 1),
                      formatC(f$max_abs_z[i], format = "f", digits = 2), f$raisons[i]))
  }
  if (!is.null(comps)) {
    cp <- comps[comps$cas == nm & comps$moteur == "lbfgsb", ]
    if (nrow(cp)) L <- c(L, sprintf("      %-22s vrai %9.4f  estime %9.4f  se %8.4f  z %6.2f",
                                    cp$composante, cp$vrai, cp$estime, cp$se, cp$z))
  }
}
writeLines(L, file.path(sortie, "resume.txt"))
cat(L, sep = "\n")
