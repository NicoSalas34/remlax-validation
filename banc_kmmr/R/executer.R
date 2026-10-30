#!/usr/bin/env Rscript
# ==============================================================================
# Un ajustement du banc remlax contre kmmr, ajoute en une ligne a un CSV.
#
#   Rscript R/executer.R --cas drops2_fa1 --outil kmmr --gpu 1 --kmmr <src kmmr> --sortie res.csv
#   Rscript R/executer.R --cas drops2_fa1 --outil remlax --voie mme_hybride --gpu 1 ...
#   Rscript R/executer.R --cas drops2_fa1 --outil asreml --kmmr <src kmmr> --sortie res.csv
# ==============================================================================
a <- commandArgs(TRUE)
arg <- function(k, d = NULL) { i <- match(paste0("--", k), a)
  if (is.na(i) || i == length(a)) d else a[i + 1L] }
ici <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(ici, "cas.R")); source(file.path(ici, "kmmr_cas.R")); source(file.path(ici, "remlax_cas.R"))
source(file.path(ici, "asreml_cas.R")); source(file.path(ici, "sim.R"))
`%||%` <- function(x, y) if (is.null(x)) y else x
suppressPackageStartupMessages(library(Matrix))

cas <- arg("cas"); outil <- arg("outil"); voie <- arg("voie", "")
gpu <- arg("gpu", "0") == "1"
cs <- charger_cas(cas, arg("kmmr"), arg("sim"))
cat(sprintf("== %s | %s %s | GPU %s | n = %d | %d genotypes\n", cas, outil, voie, gpu, nrow(cs$d), nrow(cs$K)))
res <- tryCatch(
  if (outil == "kmmr") ajuster_kmmr(cs, gpu = gpu, tol = as.numeric(arg("tol", "0.001"))) else if (outil == "asreml") ajuster_asreml(cs) else ajuster_remlax(cs, voie, backend = if (gpu) "gpu" else "cpu"),
  error = function(e) list(erreur = conditionMessage(e)))
ligne <- data.frame(etiquette = arg("etiquette", ""), cas = cas, outil = outil, voie = if (nzchar(voie)) voie else if (outil == "asreml") "ai" else paste0("ai_tol", arg("tol", "0.001")), gpu = gpu,
                    noeud = Sys.info()[["nodename"]], fils = Sys.getenv("OMP_NUM_THREADS", NA),
                    n = nrow(cs$d), q_K = nrow(cs$K), logLik = res$logLik %||% NA, n_iter = res$n_iter %||% NA,
                    converge = res$converge %||% NA, secondes = res$secondes %||% NA,
                    n_par = res$n_par %||% NA, moteur = res$moteur %||% NA,
                    eval_s = res$eval_s %||% NA, compile_s = res$compile_s %||% NA,
                    t_modele = res$t_modele %||% NA,
                    inversion = res$inversion %||% NA, erreur = res$erreur %||% "",
                    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
print(ligne)
f <- arg("sortie")
utils::write.table(ligne, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
# CHRONO PAR ETAPES, en format long a cote du CSV : une ligne par etape
ch <- res$chrono
if (length(ch)) {
  fc <- sub("\\.csv$", "_chrono.csv", f)
  lc <- data.frame(etiquette = arg("etiquette", ""), cas = cas, outil = outil, voie = ligne$voie, gpu = gpu, noeud = ligne$noeud,
                   etape = names(ch), secondes = as.numeric(ch), date = ligne$date)
  print(lc[, c("etape", "secondes")], row.names = FALSE)
  utils::write.table(lc, fc, sep = ",", row.names = FALSE, col.names = !file.exists(fc), append = file.exists(fc))
}
