# =============================================================================
# verif_ecarts_asreml.R: the benchmark cases where remlax returns a logLik
#                        higher than asreml by more than 1e-3. Does asreml stop
#                        early, or is the function different?
# -----------------------------------------------------------------------------
# For each case (same data as bench_complexite.R, same seeds), the script fits
# asreml with maxit = 100 (the benchmark call), then repeats update() until the
# logLik moves by less than 1e-6 or 40 updates. It then fits remlax (backend
# RX_BACKEND) and prints the asreml logLik path and the final gap to remlax.
# If asreml reaches remlax by iterating further, the benchmark gap comes from
# its stopping rule (relative change of logLik < 0.002 x iteration and
# components < 1 %) and the likelihoods agree.
#   Rscript benchmarks/verif_ecarts_asreml.R [--cases us9:500,us12:500,us12:2000,ige:500,usK6:2000:40]
# =============================================================================
suppressPackageStartupMessages({ library(here); library(Matrix); library(jsonlite); library(asreml) })
asreml.options(trace = FALSE, workspace = "8gb", pworkspace = "2gb")
src <- readLines(here::here("benchmarks", "bench_complexite.R"))
# reuse the benchmark definitions (generators, fits) without its loop
i0 <- grep("^# ---- (loop|boucle)", src)
args_saved <- commandArgs(trailingOnly = TRUE)
commandArgs <- function(trailingOnly = TRUE) c("--software", "asreml,remlax", "--timeout-fit", "0")
eval(parse(text = src[seq_len(i0 - 1)]))
opt2 <- function(name, default) { i <- match(paste0("--", name), args_saved); if (is.na(i) || i == length(args_saved)) default else args_saved[i + 1L] }
# The third field is the number of replicates per genotype. The size axis of
# the benchmark set q = n1 %/% 2 for the multi-trait cases (us9 at n = 500:
# n1 = 55, q = 27). round(n1 / reps) in the current generator reproduces it with
# reps = n1 / q.
CASES_CHECK <- strsplit(opt2("cases", "us9:500:2.037037,us12:500:2.05,us12:2000:2,ige:500:4,usK6:2000:40"), ",")[[1]]
OUT <- opt2("out", here::here("benchmarks", "results"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
res <- list()
for (spec in CASES_CHECK) tryCatch({
  p <- strsplit(spec, ":")[[1]]; case_ <- p[1]; n <- as.integer(p[2]); REPS <<- if (length(p) > 2) as.numeric(p[3]) else 4
  D <- gen_case(case_, n)
  cat(sprintf("\n=== %s n=%d reps=%g q=%d ===\n", case_, n, REPS, D$q))
  # remlax
  t0 <- proc.time()[["elapsed"]]; f <- fit_remlax(case_, D); t_rx <- proc.time()[["elapsed"]] - t0
  cat(sprintf("remlax: logLik %.6f | %d iterations | %d evaluations | %.1f s\n", f$logLik, f$n_iter, f$n_eval, t_rx))
  # asreml: benchmark call, then updates
  d <- D$d
  obj <- NULL
  if (grepl("^us[0-9]+$", case_)) {
    fml <- stats::as.formula(paste("cbind(", paste(grep("^y", names(d), value = TRUE), collapse = ","), ") ~ trait"))
    obj <- asreml(fml, random = ~ us(trait):gid, residual = ~ id(units):us(trait), data = d, maxit = MAXIT_REF)
  } else if (grepl("^usK[0-9]+$", case_)) {
    fml <- stats::as.formula(paste("cbind(", paste(grep("^y", names(d), value = TRUE), collapse = ","), ") ~ trait"))
    K <- D$K; obj <- asreml(fml, random = ~ us(trait):vm(gid, K), residual = ~ id(units):us(trait), data = d, maxit = MAXIT_REF)
  } else if (case_ == "ige") {
    da <- cbind(d, as.data.frame(as.matrix(D$ZnA)), as.data.frame(as.matrix(D$ZxA)))
    names(da) <- c(names(d), paste0("N", seq_len(D$QA)), paste0("X", seq_len(D$QB)))
    assign("iN", grep("^N[0-9]+$", names(da)), envir = globalenv()); assign("iX", grep("^X[0-9]+$", names(da)), envir = globalenv())
    obj <- asreml(y ~ 1, random = stats::as.formula(sprintf("~ str(~ geno + grp(nb), ~ us(2):id(%d)) + grp(xb) + ar1(row):ar1(col)", D$QA)),
                  group = list(nb = iN, xb = iX), data = da, maxit = MAXIT_REF)
  } else stop("case not handled here: ", case_)
  traj <- obj$loglik; conv <- obj$converge
  cat(sprintf("asreml, benchmark call (maxit %d): logLik %.6f | converged %s | gap to remlax %+.4e\n", MAXIT_REF, obj$loglik, conv, obj$loglik - f$logLik))
  for (k in 1:40) {
    prev <- obj$loglik
    obj <- suppressWarnings(update(obj, trace = FALSE, maxit = MAXIT_REF))
    traj <- c(traj, obj$loglik)
    if (abs(obj$loglik - prev) < 1e-6) break
  }
  cat(sprintf("asreml, after %d update(s): logLik %.6f | converged %s | gap to remlax %+.4e\n", length(traj) - 1, obj$loglik, obj$converge, obj$loglik - f$logLik))
  cat("  path:", paste(sprintf("%.4f", traj), collapse = " "), "\n")
  res[[spec]] <- list(case = case_, n = n, reps = REPS, remlax_logLik = f$logLik, remlax_backend = BK,
                      asreml_logLik_bench = traj[1], asreml_converge_bench = conv,
                      asreml_logLik_updates = obj$loglik, n_updates = length(traj) - 1, asreml_converge_final = obj$converge,
                      gap_bench = traj[1] - f$logLik, gap_final = obj$loglik - f$logLik, trajectory = traj)
  write_json(res, file.path(OUT, "verif_ecarts_asreml.json"), auto_unbox = TRUE, digits = 10, pretty = TRUE)
}, error = function(e) {
  cat("FAIL", spec, ":", conditionMessage(e), "\n")
  res[[spec]] <<- list(spec = spec, error = conditionMessage(e))
  write_json(res, file.path(OUT, "verif_ecarts_asreml.json"), auto_unbox = TRUE, digits = 10, pretty = TRUE)
})
cat("\nwritten:", file.path(OUT, "verif_ecarts_asreml.json"), "\n")