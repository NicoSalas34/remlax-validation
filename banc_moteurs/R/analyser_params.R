# ==============================================================================
# Analysis of the "number of parameters" benchmark (banc_params.sbatch): fitting
# time against the number of variance parameters, at n = 4000, for each remlax
# route and asreml.
#
#   Rscript R/analyser_params.R --res=<res_params folder> [--out=<folder>/figures]
#
# Results written before the English pass (ajustements_*.csv, French columns and
# values) are read as well and translated on reading.
# ==============================================================================
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", k, "="), "", a[1]) else d }
res <- opt("res", "res_params"); out_dir <- opt("out", file.path(res, "figures"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Old French column names and values -> English
FR_EN_COLS <- c(famille = "family", graine = "seed", logiciel = "software", voie = "route", statut = "status",
                secondes = "seconds", solveur_s = "solver_s", moteur = "engine", logLik_natif = "logLik_native",
                logLik_arbitre = "logLik_arbiter")
FR_EN_VALUES <- list(status = c(ERREUR = "ERROR", DELAI = "TIMEOUT", MEMOIRE = "MEMORY"),
                     route = c("M-creux-CPU" = "M-sparse-CPU"))
to_english <- function(d) {
  names(d) <- ifelse(names(d) %in% names(FR_EN_COLS), FR_EN_COLS[names(d)], names(d))
  for (k in intersect(names(FR_EN_VALUES), names(d))) {
    v <- FR_EN_VALUES[[k]]; i <- d[[k]] %in% names(v); d[[k]][i] <- unname(v[d[[k]][i]])
  }
  d
}

# English files first; an old French file is read only when no English file has its name
f_en <- list.files(res, "^fits_.*\\.csv$", full.names = TRUE, recursive = TRUE)
f_fr <- list.files(res, "^ajustements_.*\\.csv$", full.names = TRUE, recursive = TRUE)
f_fr <- f_fr[!sub("/ajustements_", "/fits_", f_fr) %in% f_en]
f <- c(f_en, f_fr)
l <- lapply(f, function(x) to_english(utils::read.csv(x, stringsAsFactors = FALSE)))
cols <- unique(unlist(lapply(l, names)))
fits <- do.call(rbind, lapply(l, function(d) { for (k in setdiff(cols, names(d))) d[[k]] <- NA; d[cols] }))

npar <- function(fam) {
  t <- as.integer(sub("^mt(us|igee|ige)", "", fam))
  base <- (2L * t) * (2L * t + 1L) / 2L + t * (t + 1L) / 2L
  ifelse(startsWith(fam, "mtus"), t * (t + 1L),
         ifelse(startsWith(fam, "mtigee"), base + t * (t + 1L) / 2L, base))
}
fits$family_type <- ifelse(startsWith(fits$family, "mtus"), "mtus: genotype us + residual us",
                    ifelse(startsWith(fits$family, "mtigee"), "mtigee: DGE/IGE us(2t) + plant IEE us(t) + residual us",
                           "mtige: DGE/IGE us(2t) + residual us"))
fits$p <- npar(fits$family)
fits$t <- as.integer(sub("^mt(us|igee|ige)", "", fits$family))
fits$cell <- paste(fits$family, fits$n, fits$rep)
fits$success <- fits$status == "OK" & is.finite(fits$logLik_arbiter)
ref <- tapply(ifelse(fits$success, fits$logLik_arbiter, NA), fits$cell, max, na.rm = TRUE)
ref[!is.finite(ref)] <- NA
fits$diff <- ref[fits$cell] - fits$logLik_arbiter
fits$kept <- fits$success & fits$diff <= 0.1
ROUTES <- c("D-CPU", "D-GPU", "D-CPU-LBFGS", "D-GPU-LBFGS", "M-dense-CPU", "M-dense-GPU", "M-sparse-CPU", "asreml")
COL <- c("D-CPU" = "#1f78b4", "D-GPU" = "#a6cee3", "D-CPU-LBFGS" = "#6a3d9a", "D-GPU-LBFGS" = "#cab2d6",
         "M-dense-CPU" = "#33a02c", "M-dense-GPU" = "#b2df8a", "M-sparse-CPU" = "#e31a1c", asreml = "#ff7f00")
fits$route <- factor(fits$route, ROUTES)
theme_set(theme_bw(base_size = 10) + theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95")))
save_fig <- function(p, name, w = 11, h = 5.5) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
}

tj <- fits[fits$kept, ]
# "number of parameters" benchmark: n = 1000 t (1000 units whatever t); the
# "size" benchmark (mtigee3, variable n) is drawn separately
tj$bench_par <- tj$n == 1000L * tj$t
s <- stats::aggregate(seconds ~ family_type + p + t + route, tj[tj$bench_par, ],
                      function(x) c(med = stats::median(x), lo = min(x), hi = max(x)))
s <- cbind(s[1:4], as.data.frame(s$seconds))
p1 <- ggplot(s, aes(p, med, colour = route, fill = route)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point() +
  scale_x_log10() + scale_y_log10() + facet_wrap(~ family_type) +
  scale_colour_manual(values = COL, drop = FALSE) + scale_fill_manual(values = COL, drop = FALSE) +
  labs(x = "number of variance parameters", y = "total time (s)", colour = NULL, fill = NULL,
       title = "Fitting time against the number of parameters, 1000 units, n = 1000 t (MIG slice)",
       subtitle = "median of 3 replicates, min-max band; fits within 0.1 of the best logLik")
save_fig(p1, "figP1_time_vs_parameters")
utils::write.csv(s, file.path(out_dir, "time_parameters.csv"), row.names = FALSE)

# iterations and evaluations: what drives the cost of each optimiser
it <- stats::aggregate(cbind(n_iter, n_eval) ~ family_type + p + route,
                       fits[fits$kept & fits$route != "asreml" & fits$n == 1000L * fits$t, ], stats::median)
p2 <- ggplot(it, aes(p, n_eval, colour = route)) + geom_line() + geom_point() + scale_x_log10() + scale_y_log10() +
  facet_wrap(~ family_type) + scale_colour_manual(values = COL, drop = FALSE) +
  labs(x = "number of variance parameters", y = "likelihood evaluations", colour = NULL,
       title = "Evaluations per fit: L-BFGS-B against AI")
save_fig(p2, "figP2_evaluations")

if (any(tj$family == "mtigee3")) {
  st <- stats::aggregate(seconds ~ family + n + route, tj[tj$family == "mtigee3", ], function(x) c(med = stats::median(x), lo = min(x), hi = max(x)))
  st <- cbind(st[1:3], as.data.frame(st$seconds))
  p3 <- ggplot(st, aes(n, med, colour = route, fill = route)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point() +
    scale_x_log10() + scale_y_log10() + facet_wrap(~ family) +
    scale_colour_manual(values = COL, drop = FALSE) + scale_fill_manual(values = COL, drop = FALSE) +
    labs(x = "n (observations)", y = "total time (s)", colour = NULL, fill = NULL,
         title = "Fitting time against n (MIG slice)",
         subtitle = "median of the replicates, min-max band; fits within 0.1 of the best logLik")
  save_fig(p3, "figP3_time_vs_n", 8, 5.5)
  utils::write.csv(st, file.path(out_dir, "time_size.csv"), row.names = FALSE)
}

status_tab <- as.data.frame(table(family = fits$family, route = fits$route,
                                  status = ifelse(fits$kept, "kept", ifelse(fits$success, "below optimum", fits$status))))
status_tab <- status_tab[status_tab$Freq > 0, ]
utils::write.csv(status_tab, file.path(out_dir, "status_parameters.csv"), row.names = FALSE)
cat(sprintf("%d fits, %d kept; figures in %s\n", nrow(fits), sum(fits$kept), out_dir))
