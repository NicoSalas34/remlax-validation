#!/usr/bin/env Rscript
# ==============================================================================
# Simulation checks of remlax: logLik against a dense REML arbiter written in
# R, agreement of the fitting routes (L-BFGS-B, AI-REML, mixed model equations
# with automatic or sparse inversion) and recovery of the simulated components.
#
#   Rscript verif_simulations/verif_simulations.R [--pkg=/path/remlax] [--cases=pattern]
#           [--cores=10] [--seed=20261021] [--out=folder] [--reps=1]
#
# --pkg     remlax sources loaded by pkgload::load_all; without --pkg, the
#           installed package is loaded by R/charger_remlax.R at the repository root
# --cases   regular expression on case names (default: all); --list prints them
# --cores   parallel processes (each limited to 2 JAX/BLAS threads)
# --reps    replicates per case (successive seeds)
# --engines regular expression on the routes (lbfgsb, ai, mme, mme_sparse,
#           mme_lbfgsb; default: all)
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  a <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(a)) sub(paste0("^--", name, "="), "", a[1]) else default
}
here <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})

# R fixes its BLAS threads at start-up. The script restarts itself with one
# thread per process if the caller did not set it.
if (!nzchar(Sys.getenv("OPENBLAS_NUM_THREADS"))) {
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "2", MKL_NUM_THREADS = "1")
  st <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(here, "verif_simulations.R")), shQuote(args)))
  quit(save = "no", status = st)
}

pkg     <- opt("pkg", "")
pattern <- opt("cases", ".")
cores   <- as.integer(opt("cores", "10"))
seed0   <- as.integer(opt("seed", "20261021"))
nrep    <- as.integer(opt("reps", "1"))
ENGINE_FILTER <- opt("engines", ".")
out_dir <- opt("out", file.path(here, "results"))

Sys.setenv(XLA_FLAGS = "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2",
           JAX_PLATFORMS = "cpu")
if (nzchar(pkg)) {
  suppressMessages(pkgload::load_all(pkg, quiet = TRUE, export_all = TRUE))
} else {
  source(file.path(here, "..", "R", "charger_remlax.R"))
  pkg <- sprintf("remlax %s installed in %s", utils::packageVersion("remlax"),
                 dirname(system.file(package = "remlax")))
}
options(remlax.threads = 2L, warn = 1)
for (f in c("outils.R", "cas.R", "cas_limites.R", "moteurs.R")) source(file.path(here, "R", f))

if ("--list" %in% args) { cat(names(CATALOGUE), sep = "\n"); quit(save = "no") }
case_names <- grep(pattern, names(CATALOGUE), value = TRUE)
if (!length(case_names)) stop("no case matches --cases=", pattern)
# heavy cases first, to balance the processes
heavy <- c("ige_mixture_two_species", "mt_g_fa2", "lag_rx_lag", "mt_residual_fa1", "mt_g_rr2",
           "mt_g_ante1", "mt_g_chol1", "mt_g_fa1", "spline_spl2d", "levels_corg_ge",
           "boundary_mixture_iee_zero_us_rank2", "boundary_us8_near_rank2", "boundary_fa2_true_rank1",
           "boundary_two_us_at_bound")
case_names <- c(intersect(heavy, case_names), setdiff(case_names, heavy))
tasks <- expand.grid(rep = seq_len(nrep), case = case_names, stringsAsFactors = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tmp <- tempfile("verif_"); dir.create(tmp)
py <- rx_python_check(quiet = TRUE)
cat(sprintf("remlax loaded from %s\nPython: %s (jax %s)\n%d cases x %d replicate(s), %d processes\n",
            pkg, py$versions[["executable"]], py$versions[["jax"]], length(case_names), nrep, cores))

t0 <- Sys.time()
outputs <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
  name <- tasks$case[i]; r <- tasks$rep[i]
  seed <- seed0 + 1000L * r + match(name, names(CATALOGUE))
  o <- tryCatch(run_case(name, CATALOGUE[[name]], seed, tmp), error = function(e)
    list(fits = data.frame(case = name, engine = "execution", status = "ERROR",
                           reasons = conditionMessage(e), stringsAsFactors = FALSE)))
  o$fits$replicate <- r; o$fits$seed <- seed
  if (!is.null(o$comps)) { o$comps$replicate <- r; o$comps$seed <- seed }
  message(sprintf("[%s] %-32s rep %d: %s", format(Sys.time(), "%H:%M:%S"), name, r,
                  paste(sprintf("%s=%s", o$fits$engine, o$fits$status), collapse = " ")))
  o
}, mc.cores = cores, mc.preschedule = FALSE)
unlink(tmp, recursive = TRUE)

flatten <- function(l, k) { x <- lapply(l, function(o) if (is.list(o)) o[[k]] else NULL)
  x <- x[!vapply(x, is.null, TRUE)]
  if (!length(x)) return(NULL)
  cols <- unique(unlist(lapply(x, names)))
  do.call(rbind, lapply(x, function(d) { for (c in setdiff(cols, names(d))) d[[c]] <- NA; d[cols] })) }
fits <- flatten(outputs, "fits"); comps <- flatten(outputs, "comps")
infos <- do.call(rbind, lapply(seq_along(outputs), function(i) {
  o <- outputs[[i]]
  data.frame(case = tasks$case[i], replicate = tasks$rep[i], family = .chr(o$family),
             note = .chr(o$note), seconds = .num(o$seconds), stringsAsFactors = FALSE)
}))
utils::write.csv(fits, file.path(out_dir, "fits.csv"), row.names = FALSE)
if (!is.null(comps)) utils::write.csv(comps, file.path(out_dir, "components.csv"), row.names = FALSE)
utils::write.csv(infos, file.path(out_dir, "cases.csv"), row.names = FALSE)

# --- text summary ------------------------------------------------------------------
L <- c(sprintf("Simulation checks of remlax (%s)", format(Sys.time(), "%Y-%m-%d %H:%M")),
       sprintf("Package: %s", pkg),
       sprintf("Criteria: |logLik - arbiter| / |arbiter| < %g; difference from the best route < %g; |z| < %g",
               ARBITER_TOL, ENGINE_TOL, Z_TOL),
       sprintf("Total time: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))), "")
tab <- table(fits$status)
L <- c(L, paste("Fits:", paste(sprintf("%s %d", names(tab), tab), collapse = ", ")), "")
for (nm in unique(fits$case)) {
  f <- fits[fits$case == nm, ]
  inf <- infos[infos$case == nm, ][1, ]
  L <- c(L, sprintf("== %s  [%s]  %s", nm, inf$family, inf$note))
  for (i in seq_len(nrow(f))) {
    L <- c(L, sprintf("   %-18s rep %d  %-16s logLik %s  arbiter %s  rel_diff %s  route_diff %s  max|z| %s  %s",
                      f$engine[i], f$replicate[i], f$status[i],
                      formatC(f$logLik[i], format = "f", digits = 6), formatC(f$logLik_arbiter[i], format = "f", digits = 6),
                      formatC(f$rel_diff[i], format = "e", digits = 1),
                      formatC(f$diff_between_engines[i], format = "e", digits = 1),
                      formatC(f$max_abs_z[i], format = "f", digits = 2), f$reasons[i]))
  }
  if (!is.null(comps)) {
    cp <- comps[comps$case == nm & comps$engine == "lbfgsb", ]
    if (nrow(cp)) L <- c(L, sprintf("      %-22s true %9.4f  estimate %9.4f  se %8.4f  z %6.2f",
                                    cp$component, cp$true, cp$estimate, cp$se, cp$z))
  }
}
writeLines(L, file.path(out_dir, "summary.txt"))
cat(L, sep = "\n")
