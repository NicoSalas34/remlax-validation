# ==============================================================================
# Dimensions of each benchmark case, rebuilt from the seed (the data are
# deterministic): number of genotypes (or of levels of the main genetic term),
# observations per level and size m of the mixed model equations. n grows by
# adding genotypes, at a roughly constant number of observations per genotype.
# The analysis must show it, so that an effect of n is not read as an effect
# of replication.
#
#   Rscript R/dimensions.R --out=<folder of the benchmark CSVs> [--pkg=<remlax>] [--cores=6]
#
# Output: <out>/dimensions.csv, one row per (family, n, replicate).
# ==============================================================================
args <- commandArgs(TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", k, "="), "", a[1]) else d }
out_dir <- opt("out", "results"); pkg <- opt("pkg", ""); cores <- as.integer(opt("cores", "6"))
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
if (nzchar(pkg)) suppressMessages(pkgload::load_all(pkg, quiet = TRUE)) else suppressMessages(library(remlax))
source(file.path(here, "..", "verif_simulations", "R", "outils.R")); source(file.path(here, "R", "familles.R"))

# main genetic term: the first term of the model. For exp, its levels are
# positions on the transect, not genotypes
LEVEL <- c(iid = "genotypes", grm = "genotypes", us = "genotypes", ped = "animals", ige = "genotypes",
           field = "genotypes", resid = "genotypes", exp = "positions", fa = "genotypes")
SIZES <- c(500, 1000, 1500, 2000, 3000, 4000, 6000, 8000)
case_grid <- expand.grid(rep = 1:5, n = SIZES, family = names(BUILDERS), stringsAsFactors = FALSE)

one <- function(i) {
  f <- case_grid$family[i]; n <- case_grid$n[i]; r <- case_grid$rep[i]
  sim <- BUILDERS[[f]](n, seed_of(f, n, r))
  m <- sim$models$mme
  Zs <- lapply(m$terms, function(t) if (is.list(t$Z)) t$Z else list(t$Z))
  z1 <- Zs[[1]][[1]]
  nt <- if (f %in% c("us", "fa")) length(Zs[[1]]) else 1L       # traits
  nlev <- ncol(z1)
  obs_lev <- Matrix::colSums(z1 != 0)                            # observations of each level (first trait)
  data.frame(family = f, n = n, rep = r, n_obs = length(sim$y), n_traits = nt,
             level = LEVEL[[f]], n_levels = nlev, levels_observed = sum(obs_lev > 0),
             obs_per_level = length(sim$y) / (nlev * nt), obs_per_level_min = min(obs_lev),
             obs_per_level_max = max(obs_lev), p = ncol(sim$X),
             m = ncol(sim$X) + sum(vapply(unlist(Zs, recursive = FALSE), ncol, 0L)),
             stringsAsFactors = FALSE)
}
rows <- parallel::mclapply(seq_len(nrow(case_grid)), function(i) tryCatch(one(i), error = function(e) {
  message("failure ", paste(case_grid[i, ], collapse = " "), ": ", conditionMessage(e)); NULL }), mc.cores = cores)
d <- do.call(rbind, rows)
utils::write.csv(d, file.path(out_dir, "dimensions.csv"), row.names = FALSE)
print(stats::aggregate(cbind(n_levels, obs_per_level, m) ~ family + n, d, stats::median))
