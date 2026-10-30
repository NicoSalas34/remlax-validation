#!/usr/bin/env Rscript
# ==============================================================================
# One fit of the remlax against kmmr benchmark, appended as one line to a CSV.
#
#   Rscript R/executer.R --case drops2_fa1 --tool kmmr --gpu 1 --kmmr <kmmr src> --out res.csv
#   Rscript R/executer.R --case drops2_fa1 --tool remlax --route mme_hybrid --gpu 1 ...
#   Rscript R/executer.R --case drops2_fa1 --tool asreml --kmmr <kmmr src> --out res.csv
#
# Timings by step go in long format to <out>_timings.csv, next to <out>.
# ==============================================================================
a <- commandArgs(TRUE)
arg <- function(k, d = NULL) { i <- match(paste0("--", k), a)
  if (is.na(i) || i == length(a)) d else a[i + 1L] }
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "cas.R")); source(file.path(here, "kmmr_cas.R")); source(file.path(here, "remlax_cas.R"))
source(file.path(here, "asreml_cas.R")); source(file.path(here, "sim.R"))
`%||%` <- function(x, y) if (is.null(x)) y else x
suppressPackageStartupMessages(library(Matrix))

case <- arg("case"); tool <- arg("tool"); route <- arg("route", "")
# old French route names of the stored task lists -> English names
FR_ROUTES <- c(blocs = "blocks", mme_creux = "mme_sparse", mme_hybride = "mme_hybrid")
if (route %in% names(FR_ROUTES)) route <- FR_ROUTES[[route]]
gpu <- arg("gpu", "0") == "1"
cs <- load_case(case, arg("kmmr"), arg("sim"))
cat(sprintf("== %s | %s %s | GPU %s | n = %d | %d genotypes\n", case, tool, route, gpu, nrow(cs$d), nrow(cs$K)))
res <- tryCatch(
  if (tool == "kmmr") fit_kmmr(cs, gpu = gpu, tol = as.numeric(arg("tol", "0.001"))) else if (tool == "asreml") fit_asreml(cs) else fit_remlax(cs, route, backend = if (gpu) "gpu" else "cpu"),
  error = function(e) list(error = conditionMessage(e)))
row <- data.frame(label = arg("label", ""), case = case, tool = tool, route = if (nzchar(route)) route else if (tool == "asreml") "ai" else paste0("ai_tol", arg("tol", "0.001")), gpu = gpu,
                  node = Sys.info()[["nodename"]], threads = Sys.getenv("OMP_NUM_THREADS", NA),
                  n = nrow(cs$d), q_K = nrow(cs$K), logLik = res$logLik %||% NA, n_iter = res$n_iter %||% NA,
                  converged = res$converged %||% NA, seconds = res$seconds %||% NA,
                  n_par = res$n_par %||% NA, engine = res$engine %||% NA,
                  eval_s = res$eval_s %||% NA, compile_s = res$compile_s %||% NA,
                  t_model = res$t_model %||% NA,
                  inversion = res$inversion %||% NA, error = res$error %||% "",
                  date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
print(row)
f <- arg("out")
utils::write.table(row, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
# TIMINGS BY STEP, in long format next to the CSV: one line per step
ch <- res$timings
if (length(ch)) {
  fc <- sub("\\.csv$", "_timings.csv", f)
  lc <- data.frame(label = arg("label", ""), case = case, tool = tool, route = row$route, gpu = gpu, node = row$node,
                   step = names(ch), seconds = as.numeric(ch), date = row$date)
  print(lc[, c("step", "seconds")], row.names = FALSE)
  utils::write.table(lc, fc, sep = ",", row.names = FALSE, col.names = !file.exists(fc), append = file.exists(fc))
}
