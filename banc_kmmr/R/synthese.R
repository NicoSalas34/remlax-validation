#!/usr/bin/env Rscript
# ==============================================================================
# Synthesis of the remlax / kmmr / asreml-R benchmark: Markdown tables.
#
#   Rscript R/synthese.R --res <folder of the genotoul CSV files> --asreml <banc_asreml.csv>
#
# Times compared:
#   kmmr   : kmmr() call + inversion of K done in R before the call
#            (+ starting values for the two-stage cases);
#   remlax : whole rx_fit() call, start of the Python solver included.
#            "warm" = second fit, compiled programs read back from disk;
#   asreml : asreml() call alone (inversion of K included), meso nodes.
# For remlax, each task takes the most recent version available among the
# files given in the order of --versions.
#
# The files are bench_<target>[_<version>].csv with bench_..._timings.csv.
# The stored files of 2026-10-27 use the old French names (banc_*.csv,
# banc_*_chrono.csv) and French columns and values: they are read too, and
# mapped to English right after reading.
# ==============================================================================
a <- commandArgs(TRUE)
arg <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1L] }
res <- arg("res", ".")
versions <- strsplit(arg("versions", "v5,v3"), ",")[[1]]
f_asr <- arg("asreml", NULL)

# old French columns and values -> English
FR_COLS <- c(etiquette = "label", cas = "case", outil = "tool", voie = "route", noeud = "node",
             fils = "threads", converge = "converged", secondes = "seconds", moteur = "engine",
             t_modele = "t_model", erreur = "error", etape = "step", delai = "timeout")
FR_ROUTES <- c(blocs = "blocks", mme_creux = "mme_sparse", mme_hybride = "mme_hybrid")
# timing steps: kmmr and asreml steps of this benchmark, and remlax steps
# (old French names of the CSV files already written -> current names)
FR_STEPS <- c(valeurs_depart = "start_values", premiere_iteration = "first_iteration",
              iterations_suivantes = "later_iterations", hors_solveur = "outside_solver",
              ajustement = "fit", hors_iterations = "outside_iterations",
              demarrage_serveur = "server_startup", solveur_total = "solver_total", lecture_R = "read_R",
              origine = "origin", demarrage = "startup", peripherique = "device", lecture = "read",
              finition = "finish", effets = "effects", hessien = "hessian", ecriture = "write",
              transmission = "transfer")
recode <- function(x, map) ifelse(x %in% names(map), map[x], x)
to_english <- function(d) {
  names(d) <- recode(names(d), FR_COLS)
  if (!is.null(d$label)) {        # cold / warm / single, with an optional r2_ prefix
    d$label <- sub("(^|_)froid$", "\\1cold", d$label)
    d$label <- sub("(^|_)chaud$", "\\1warm", d$label)
    d$label <- sub("(^|_)unique$", "\\1single", d$label)
  }
  if (!is.null(d$route)) d$route <- recode(d$route, FR_ROUTES)
  if (!is.null(d$step)) d$step <- recode(d$step, FR_STEPS)
  d
}
read_res <- function(f) if (!is.null(f) && file.exists(f)) to_english(utils::read.csv(f, stringsAsFactors = FALSE)) else NULL
# results file of a target and version: English name first, then the old French one
results_file <- function(target, v) {
  stem <- sprintf("%s%s.csv", target, if (nzchar(v)) paste0("_", v) else "")
  f <- file.path(res, paste0(c("bench_", "banc_"), stem))
  if (file.exists(f[1])) f[1] else f[2]
}
# timings file next to a results file
timings_file <- function(f) {
  g <- c(sub("\\.csv$", "_timings.csv", f), sub("\\.csv$", "_chrono.csv", f))
  if (file.exists(g[1])) g[1] else g[2]
}

case_order <- c("drops2_fa1", "drops2_fa3", "drops1_fa1", "drops1_us",
                "krause_gblup", "krause_fa2", "krause_us")

table_md <- function(target) {
  # kmmr: the v3 file (same nodes, same conditions)
  km <- NULL; kch <- NULL
  for (v in c("v3", "v2", "")) {
    f <- results_file(target, v)
    d <- read_res(f)
    if (!is.null(d)) {
      d <- d[d$tool == "kmmr" & d$route == "ai_tol0.000001", ]
      ch <- read_res(timings_file(f))
      d <- d[!(d$case %in% km$case), ]
      km <- rbind(km, d); kch <- rbind(kch, ch)
    }
  }
  rx <- NULL; rch <- NULL
  for (v in versions) {
    f <- results_file(target, v)
    d <- read_res(f)
    if (is.null(d)) next
    d <- d[d$tool == "remlax" & d$route == "mme" & d$label == "warm", ]
    d <- d[!(d$case %in% rx$case), ]
    if (nrow(d)) d$version <- v
    rx <- rbind(rx, d); rch <- rbind(rch, read_res(timings_file(f)))
  }
  lines <- c(sprintf("| case | kmmr logLik | kmmr it. | kmmr s | remlax logLik | remlax it. | remlax s | remlax s/eval | version |"),
             "|---|---|---|---|---|---|---|---|---|")
  for (cs in case_order) {
    k <- km[km$case == cs, ][1, ]
    r <- rx[rx$case == cs, ][1, ]
    kt <- NA
    if (!is.na(k$seconds)) {
      x <- kch[kch$case == cs & kch$tool == "kmmr" & kch$route == "ai_tol0.000001", ]
      extra <- sum(x$seconds[x$step %in% c("inverse_K", "start_values")])
      kt <- k$seconds + extra
    }
    lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s |", cs,
      if (is.na(k$logLik)) "" else sprintf("%.4f", k$logLik), if (is.na(k$n_iter)) "" else k$n_iter,
      if (is.na(kt)) "" else sprintf("%.1f", kt),
      if (is.null(r) || is.na(r$logLik)) "" else sprintf("%.4f", r$logLik),
      if (is.null(r) || is.na(r$n_iter)) "" else r$n_iter,
      if (is.null(r) || is.na(r$seconds)) "" else sprintf("%.1f", r$seconds),
      if (is.null(r) || is.na(r$eval_s)) "" else sprintf("%.3f", r$eval_s),
      if (is.null(r) || is.null(r$version)) "" else r$version))
  }
  lines
}

for (target in c("cpu", "gpu")) {
  cat(sprintf("\n### %s\n\n", toupper(target)))
  cat(table_md(target), sep = "\n")
}

if (!is.null(f_asr) && file.exists(f_asr)) {
  d <- read_res(f_asr)
  d <- d[!is.na(d$logLik), ]
  cat("\n### asreml-R (CPU, meso, 8 cores)\n\n| case | logLik | it. | s | converged |\n|---|---|---|---|---|\n")
  for (cs in case_order) {
    x <- d[d$case == cs, ]
    if (!nrow(x)) { cat(sprintf("| %s | | | | |\n", cs)); next }
    x <- x[nrow(x), ]
    cat(sprintf("| %s | %.4f | %s | %.1f | %s |\n", cs, x$logLik, x$n_iter, x$seconds, x$converged))
  }
}

# remlax phases (warm)
for (target in c("cpu", "gpu")) {
  for (v in versions) {
    ch <- read_res(timings_file(results_file(target, v)))
    if (is.null(ch)) next
    ch <- ch[ch$tool == "remlax" & ch$route == "mme" & ch$label == "warm", ]
    if (!nrow(ch)) next
    steps <- c("server_startup", "preparation", "compilation", "optimisation", "finish",
            "effects", "inference", "export_R")
    cat(sprintf("\n### Warm remlax phases (%s, %s), seconds\n\n| case | %s |\n|%s\n", toupper(target), v,
                paste(steps, collapse = " | "), paste(rep("---|", length(steps) + 1), collapse = "")))
    for (cs in intersect(case_order, unique(ch$case))) {
      x <- ch[ch$case == cs, ]
      cat(sprintf("| %s | %s |\n", cs, paste(sprintf("%.2f", vapply(steps, function(e) sum(x$seconds[x$step == e]), 0)), collapse = " | ")))
    }
  }
}
