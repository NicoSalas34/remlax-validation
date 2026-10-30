# ==============================================================================
# One benchmark case: family, n, replicate. Every remlax route and the external
# software on the SAME data; one shared dense REML arbiter.
#
#   Rscript R/executer.R --family=grm --n=2000 --rep=1 --out=results [--pkg=<remlax>]
#            [--routes=D-CPU,D-GPU,M-dense-CPU,M-dense-GPU,M-sparse-CPU] [--external=asreml,sommer,lme4]
#            [--timeout=1800]
#
# Each fit runs in its own R process (callr), under a time limit. Software that
# does not converge or fills the memory does not block the case.
# Outputs (one row per fit): <out>/fits_<family>_<n>_<rep>.csv
# and <out>/parameters_<family>_<n>_<rep>.csv.
# ==============================================================================
args <- commandArgs(TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", k, "="), "", a[1]) else d }
fam <- opt("family", "iid"); n <- as.integer(opt("n", "500")); rep <- as.integer(opt("rep", "1"))
out_dir <- opt("out", "results"); pkg <- opt("pkg", "")
routes <- strsplit(opt("routes", "D-CPU,D-GPU,M-dense-CPU,M-dense-GPU,M-sparse-CPU"), ",")[[1]]
external <- strsplit(opt("external", "asreml,sommer,lme4"), ",")[[1]]
timeout <- as.numeric(opt("timeout", "1800"))
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
tools_file <- normalizePath(file.path(here, "..", "verif_simulations", "R", "outils.R"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_pkg <- function(pkg) {
  if (nzchar(pkg)) suppressMessages(pkgload::load_all(pkg, quiet = TRUE)) else suppressMessages(library(remlax))
}
load_pkg(pkg)
source(tools_file); source(file.path(here, "R", "familles.R"))
seed <- seed_of(fam, n, rep)
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

# RESUME. Each fit is written as soon as it ends (.partial files), and a
# marker running_<case> names the fit in progress. If the task dies (memory
# exceeded: SLURM kills the whole process group), the rerun keeps what is
# done, records MEMORY for the fit named by the marker and does the rest.
# asreml licence or workspace failures are redone.
key <- sprintf("%s_%d_%d", fam, n, rep)
f_fits <- file.path(out_dir, paste0("fits_", key, ".csv"))
f_pars <- file.path(out_dir, paste0("parameters_", key, ".csv"))
f_fits_p <- paste0(f_fits, ".partial"); f_pars_p <- paste0(f_pars, ".partial")
f_running <- file.path(out_dir, paste0("running_", key))
read_if <- function(a, b) if (file.exists(a)) utils::read.csv(a, stringsAsFactors = FALSE) else
  if (file.exists(b)) utils::read.csv(b, stringsAsFactors = FALSE) else NULL
rows0 <- read_if(f_fits, f_fits_p); params0 <- read_if(f_pars, f_pars_p)
if (!is.null(rows0)) {
  redo <- rows0$status == "ERROR" & grepl("License error|licenses in use|workspace", rows0$message)
  rows0 <- rows0[!redo, , drop = FALSE]
  if (!is.null(params0)) params0 <- params0[paste(params0$software, params0$route) %in% paste(rows0$software, rows0$route), , drop = FALSE]
}
done <- if (is.null(rows0)) character(0) else rows0$route
killed <- if (file.exists(f_running)) readLines(f_running, warn = FALSE)[1] else NA_character_

t0 <- Sys.time()
sim <- BUILDERS[[fam]](n, seed)
t_gen <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
n_obs <- length(sim$y)
cat(sprintf("== %s n=%d (%d obs) rep=%d seed=%d | simulation %.1f s\n", fam, n, n_obs, rep, seed, t_gen))

ROUTES <- list("D-CPU" = list(m = "dense", backend = "cpu", engine = "dense", inversion = "auto"),
               "D-CPU-LBFGS" = list(m = "dense", backend = "cpu", engine = "dense", inversion = "auto", optimizer = "lbfgsb"),
               "D-GPU-LBFGS" = list(m = "dense", backend = "gpu", engine = "dense", inversion = "auto", optimizer = "lbfgsb"),
               "D-GPU" = list(m = "dense", backend = "gpu", engine = "dense", inversion = "auto"),
               "M-dense-CPU" = list(m = "mme", backend = "cpu", engine = "mme", inversion = "dense"),
               "M-dense-GPU" = list(m = "mme", backend = "gpu", engine = "mme", inversion = "dense"),
               "M-sparse-CPU" = list(m = "mme", backend = "cpu", engine = "mme", inversion = "sparse"))

# "--routes=none": no remlax route (one task per route in banc_iee_long)
routes <- intersect(routes, names(ROUTES))

# one remlax fit in its own process; the time is measured INSIDE that
# process, around the rx_fit call alone
run_remlax <- function(v) {
  model <- sim$models[[v$m]]
  callr::r(function(model, v, pkg, tools_file) {
    if (nzchar(pkg)) suppressMessages(pkgload::load_all(pkg, quiet = TRUE)) else suppressMessages(library(remlax))
    source(tools_file)
    t0 <- Sys.time()
    f <- rx_fit(model, backend = v$backend, engine = v$engine, inversion = v$inversion,
                optimizer = if (is.null(v$optimizer)) "ai" else v$optimizer,
                hessian = FALSE, blups = FALSE, verbose = FALSE)
    s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    p <- par_from_fit(f)
    list(est = p$est, logLik = f$logLik, conv = isTRUE(as.logical(f$scipy_success)), grad_rel = f$grad_rel %||% NA_real_,
         seconds = s, solver_s = f$seconds, compile_s = f$compile_s, eval_s = f$eval_s,
         n_iter = f$n_iter, n_eval = f$n_eval, engine = f$engine %||% v$engine, inversion = f$inversion %||% "",
         n_par = f$n_par,
         message = f$scipy_message)
  }, args = list(model, v, pkg, tools_file), timeout = timeout)
}
run_external <- function(name) {
  callr::r(function(fam, n, seed, pkg, tools_file, families_file, name) {
    if (nzchar(pkg)) suppressMessages(pkgload::load_all(pkg, quiet = TRUE)) else suppressMessages(library(remlax))
    source(tools_file); source(families_file)
    if (name == "asreml") { suppressMessages(library(asreml)); asreml::asreml.options(trace = FALSE, workspace = "4gb", pworkspace = "1gb") }
    sim <- BUILDERS[[fam]](n, seed)
    t0 <- Sys.time()
    r <- sim$external[[name]]()
    r$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    r
  }, args = list(fam, n, seed, pkg, tools_file, file.path(here, "R", "familles.R"), name), timeout = timeout)
}

rows <- if (is.null(rows0)) list() else list(rows0)
params <- if (is.null(params0) || !nrow(params0)) list() else list(params0)
write_out <- function(final = FALSE) {
  a <- if (final) f_fits else f_fits_p; b <- if (final) f_pars else f_pars_p
  if (length(rows)) utils::write.csv(do.call(rbind, rows), a, row.names = FALSE)
  if (length(params)) utils::write.csv(do.call(rbind, params), b, row.names = FALSE)
  if (final) unlink(c(f_fits_p, f_pars_p, f_running))
}
record <- function(software, route, r, status, error = "") {
  la <- NA_real_
  if (!is.null(r$est) && all(names(sim$truth) %in% names(r$est)))
    la <- tryCatch(reml_dense(sim$V_of(r$est[names(sim$truth)]), sim$y, sim$X)$logLik, error = function(e) NA_real_)
  rows[[length(rows) + 1L]] <<- data.frame(
    family = fam, n = n, n_obs = n_obs, rep = rep, seed = seed, software = software, route = route,
    status = status, seconds = r$seconds %||% NA_real_, solver_s = r$solver_s %||% NA_real_,
    compile_s = r$compile_s %||% NA_real_, eval_s = r$eval_s %||% NA_real_,
    n_iter = r$n_iter %||% NA_integer_, n_eval = r$n_eval %||% NA_integer_,
    n_par = r$n_par %||% NA_integer_,
    engine = r$engine %||% "", inversion = r$inversion %||% "", conv = r$conv %||% NA, grad_rel = r$grad_rel %||% NA_real_,
    logLik_native = r$logLik %||% NA_real_, logLik_arbiter = la,
    message = substr(gsub("[\r\n,]+", " ", paste(r$message %||% "", error)), 1, 300), stringsAsFactors = FALSE)
  if (!is.null(r$est)) {
    nm <- intersect(names(sim$truth), names(r$est))
    params[[length(params) + 1L]] <<- data.frame(family = fam, n = n, rep = rep, software = software, route = route,
                                                 parameter = nm, true = unname(sim$truth[nm]), estimate = unname(r$est[nm]),
                                                 stringsAsFactors = FALSE)
  }
  write_out()
}
# status of a failure: timeout, memory (GPU or killed process), or error
status_of <- function(e) {
  m <- conditionMessage(e)
  if (grepl("timeout|timed out", m, ignore.case = TRUE)) "TIMEOUT"
  else if (grepl("RESOURCE_EXHAUSTED|out of memory|ran out of memory|bfc_allocator|cannot allocate|killed", m, ignore.case = TRUE)) "MEMORY"
  else "ERROR"
}
fit_route <- function(software, route, launch) {
  if (route %in% done) { cat(sprintf("   %-26s already done\n", paste(software, route))); return(invisible()) }
  if (identical(route, killed)) {
    record(software, route, list(), "MEMORY", "process killed by SLURM (task memory exceeded)")
  } else {
    writeLines(route, f_running)
    r <- tryCatch(launch(), error = function(e) e)
    if (inherits(r, "error")) record(software, route, list(), status_of(r), conditionMessage(r))
    else record(software, route, r, "OK")
    unlink(f_running)
  }
  cat(sprintf("   %-26s %s\n", paste(software, route), utils::tail(rows, 1)[[1]]$status))
}
# asreml: online licence with counted tokens; wait until a token is free
run_asreml <- function() {
  for (k in 1:40) {
    r <- tryCatch(run_external("asreml"), error = function(e) e)
    if (!inherits(r, "error") || !grepl("License error|licenses in use", conditionMessage(r))) break
    Sys.sleep(30)
  }
  if (inherits(r, "error")) stop(r)
  r
}

for (vn in routes) fit_route("remlax", vn, function() run_remlax(ROUTES[[vn]]))
for (name in intersect(external, names(sim$external)))
  fit_route(name, name, if (name == "asreml") run_asreml else function() run_external(name))
write_out(final = TRUE)
cat(sprintf("== end %s in %.1f s\n", key, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
