# ==============================================================================
# Analysis of the engine benchmark: gathers the CSVs of executer.R, then draws
# the comparison figures (time, CPU against GPU, MME against dense, inversion,
# accuracy of the logLik, parameter recovery, agreement between software).
#
#   Rscript R/analyser.R --res=<folder of the CSVs> [--out=<folder>/figures]
#
# A cell = (family, n, replicate). The same seed serves every route of a cell.
# Time ratios and logLik differences are computed WITHIN the cell, then the
# five replicates are summarised (median, min, max).
# The reference logLik of a cell is the best ARBITER logLik (shared dense REML,
# evaluated at the estimates of each software) among the successful fits. The
# native logLiks do not share the same constant.
#
# Results written before the English pass (ajustements_*.csv, parametres_*.csv,
# French columns and values) are read as well and translated on reading.
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
args <- commandArgs(TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", k, "="), "", a[1]) else d }
res <- opt("res", "results"); out_dir <- opt("out", file.path(res, "figures"))
fdim <- opt("dimensions", file.path(res, "dimensions.csv"))      # written by R/dimensions.R
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Old French column names and values -> English
FR_EN_COLS <- c(famille = "family", graine = "seed", logiciel = "software", voie = "route", statut = "status",
                secondes = "seconds", solveur_s = "solver_s", moteur = "engine", logLik_natif = "logLik_native",
                logLik_arbitre = "logLik_arbiter", parametre = "parameter", vrai = "true", estime = "estimate",
                n_caracteres = "n_traits", niveau = "level", n_niveaux = "n_levels", niveaux_observes = "levels_observed",
                obs_par_niveau = "obs_per_level", obs_par_niveau_min = "obs_per_level_min",
                obs_par_niveau_max = "obs_per_level_max")
FR_EN_VALUES <- list(status = c(ERREUR = "ERROR", DELAI = "TIMEOUT", MEMOIRE = "MEMORY"),
                     route = c("M-creux-CPU" = "M-sparse-CPU"), family = c(champ = "field"),
                     level = c(animaux = "animals"))
to_english <- function(d) {
  names(d) <- ifelse(names(d) %in% names(FR_EN_COLS), FR_EN_COLS[names(d)], names(d))
  for (k in intersect(names(FR_EN_VALUES), names(d))) {
    v <- FR_EN_VALUES[[k]]; i <- d[[k]] %in% names(v); d[[k]][i] <- unname(v[d[[k]][i]])
  }
  if ("parameter" %in% names(d)) d$parameter <- sub("^residuelle", "residual", sub("^champ", "field", d$parameter))
  d
}

# English files first; an old French file is read only when its case has no English file
read_all <- function(prefix_en, prefix_fr) {
  f_en <- list.files(res, paste0("^", prefix_en, "_.*\\.csv$"), full.names = TRUE)
  f_fr <- list.files(res, paste0("^", prefix_fr, "_.*\\.csv$"), full.names = TRUE)
  key_fr <- sub("^champ_", "field_", sub(paste0("^", prefix_fr, "_"), "", basename(f_fr)))
  f_fr <- f_fr[!key_fr %in% sub(paste0("^", prefix_en, "_"), "", basename(f_en))]
  f <- c(f_en, f_fr)
  if (!length(f)) stop("no ", prefix_en, "_*.csv (or ", prefix_fr, "_*.csv) file in ", res)
  l <- lapply(f, function(x) to_english(utils::read.csv(x, stringsAsFactors = FALSE)))
  # columns added over the versions of the benchmark (n_par, grad_rel): fill them in
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { for (k in setdiff(cols, names(d))) d[[k]] <- NA; d[cols] }))
}
fits <- read_all("fits", "ajustements")
pa <- read_all("parameters", "parametres")

FAMILIES <- c("iid", "grm", "us", "ped", "ige", "field", "resid", "exp", "fa")
ROUTES <- c("D-CPU", "D-GPU", "M-dense-CPU", "M-dense-GPU", "M-sparse-CPU", "asreml", "sommer", "lme4")
COL <- c("D-CPU" = "#1f78b4", "D-GPU" = "#a6cee3", "M-dense-CPU" = "#33a02c", "M-dense-GPU" = "#b2df8a",
         "M-sparse-CPU" = "#e31a1c", asreml = "#ff7f00", sommer = "#6a3d9a", lme4 = "#b15928")
fits$family <- factor(fits$family, FAMILIES); pa$family <- factor(pa$family, FAMILIES)
fits$route <- factor(fits$route, ROUTES); pa$route <- factor(pa$route, ROUTES)
fits$cell <- paste(fits$family, fits$n, fits$rep)
pa$cell <- paste(pa$family, pa$n, pa$rep)

# DIMENSIONS. n grows by adding genotypes (number of observations per genotype
# roughly fixed), not by adding replicates: m, the size of C, grows with n.
# Each panel carries the rule of its family, so that an effect of n is not read
# as an effect of the number of replicates per genotype.
dims <- if (file.exists(fdim)) to_english(utils::read.csv(fdim, stringsAsFactors = FALSE)) else NULL
if (!is.null(dims)) {
  dims$family <- factor(dims$family, FAMILIES)
  rule <- vapply(FAMILIES, function(f) {
    d <- dims[dims$family == f, ]
    if (!nrow(d)) return(f)
    k <- stats::median(d$n_obs / d$n_levels)
    k <- if (abs(k - round(k)) < 0.15) round(k) else round(k, 1)
    nt <- d$n_traits[1]
    sprintf("%s: %s = %s%s", f, d$level[1], if (k == 1) "n" else paste0("n/", k),
            if (nt > 1) sprintf(" (%d traits)", nt) else "")
  }, "")
  lab_fam <- as_labeller(rule)
  fits <- merge(fits, dims[, c("family", "n", "rep", "n_levels", "m")], by = c("family", "n", "rep"), all.x = TRUE)
} else {
  message("dimensions.csv missing: run R/dimensions.R; panels without a scaling rule")
  lab_fam <- label_value
}

# GPU memory failures written as ERROR by the first tasks of the array
fits$status[fits$status == "ERROR" & grepl("bfc_allocator|RESOURCE_EXHAUSTED|out of memory", fits$message)] <- "MEMORY"
# a fit counts if it ended and its arbiter logLik is known
fits$success <- fits$status == "OK" & is.finite(fits$logLik_arbiter)
ref <- tapply(ifelse(fits$success, fits$logLik_arbiter, NA), fits$cell, max, na.rm = TRUE)
ref[!is.finite(ref)] <- NA
fits$diff_logLik <- ref[fits$cell] - fits$logLik_arbiter          # >= 0, in logLik units
# tolerance: below 1e-3 the difference is numerical; above, the fit is below the optimum
TOL <- 1e-3
fits$accurate <- fits$success & fits$diff_logLik <= TOL
# time and parameters: a fit within 0.1 logLik unit is kept (asreml sometimes
# stops a few 1e-3 short of the optimum); beyond that it is discarded
TOL_USE <- 0.1
fits$kept <- fits$success & fits$diff_logLik <= TOL_USE

summarise_by <- function(d, val, by = c("family", "n", "route")) {
  s <- stats::aggregate(d[[val]], d[by], function(x) c(med = stats::median(x), lo = min(x), hi = max(x), k = length(x)))
  cbind(s[by], as.data.frame(s$x))
}
theme_bench <- theme_bw(base_size = 10) + theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"))
save_fig <- function(p, name, w = 11, h = 8) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
}

# -------------------------------------------------------------- status summary
status_tab <- as.data.frame(table(family = fits$family, route = fits$route, status = ifelse(fits$accurate, "accurate",
                            ifelse(fits$success, "below optimum", fits$status))))
status_tab <- status_tab[status_tab$Freq > 0, ]
utils::write.csv(status_tab, file.path(out_dir, "summary_status.csv"), row.names = FALSE)
p <- ggplot(status_tab, aes(route, Freq, fill = status)) + geom_col() + facet_wrap(~ family, labeller = lab_fam, nrow = 3) +
  scale_fill_manual(values = c(accurate = "#4daf4a", "below optimum" = "#ff7f00", ERROR = "#e41a1c", TIMEOUT = "#984ea3", MEMORY = "#377eb8")) +
  labs(x = NULL, y = "fits", fill = NULL, title = "Status of the fits by family and route",
       subtitle = sprintf("accurate: arbiter logLik within %g of the best of the cell", TOL)) +
  theme_bench + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p, "fig0_status")

# ------------------------------------------------------- dimensions against n
if (!is.null(dims)) {
  dl <- rbind(data.frame(dims[c("family", "n")], measure = "genotypes (or levels) of the main term", v = dims$n_levels),
              data.frame(dims[c("family", "n")], measure = "observations per level and per trait", v = dims$obs_per_level),
              data.frame(dims[c("family", "n")], measure = "m (size of C)", v = dims$m))
  dl <- stats::aggregate(v ~ family + n + measure, dl, stats::median)
  p <- ggplot(dl, aes(n, v, colour = family)) + geom_line() + geom_point(size = 1) +
    scale_x_log10() + scale_y_log10() + facet_wrap(~ measure, scales = "free_y", nrow = 1) +
    labs(x = "n (observations)", y = NULL, colour = NULL, title = "Dimensions of the simulated cases",
         subtitle = paste("n grows by adding genotypes; the number of observations per genotype stays roughly fixed.",
                          "Median of 5 replicates.")) +
    theme_bench
  save_fig(p, "figD_dimensions", 12, 4.5)
  dt <- stats::aggregate(cbind(n_levels, obs_per_level, obs_per_level_min, obs_per_level_max, m) ~ family + level + n,
                         dims, stats::median)
  utils::write.csv(dt[order(dt$family, dt$n), ], file.path(out_dir, "dimensions_summary.csv"), row.names = FALSE)
}

# ------------------------------------------------------------- time against n
tj <- fits[fits$kept, ]
st <- summarise_by(tj, "seconds")
if (!is.null(dims)) st <- merge(st, stats::aggregate(cbind(n_levels, m) ~ family + n, fits, stats::median),
                                by = c("family", "n"), all.x = TRUE)
p <- ggplot(st, aes(n, med, colour = route, fill = route)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point(size = 1) +
  scale_x_log10() + scale_y_log10() + facet_wrap(~ family, labeller = lab_fam, nrow = 3, scales = "free_y") +
  scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
  labs(x = "n (observations)", y = "total time (s)", colour = NULL, fill = NULL,
       title = "Fitting time", subtitle = sprintf("median of 5 replicates, min-max band; fits within %g of the best logLik", TOL_USE)) +
  theme_bench
save_fig(p, "fig1_time")
utils::write.csv(st, file.path(out_dir, "time_summary.csv"), row.names = FALSE)
# remlax: time of the Python solver alone (without process start-up, JAX import
# and package export), beside the total time of the external software
ts <- tj; rx <- ts$route %in% ROUTES[1:5]
ts$seconds[rx] <- ts$solver_s[rx]
ts <- ts[is.finite(ts$seconds), ]
p <- ggplot(summarise_by(ts, "seconds"), aes(n, med, colour = route, fill = route)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point(size = 1) +
  scale_x_log10() + scale_y_log10() + facet_wrap(~ family, labeller = lab_fam, nrow = 3, scales = "free_y") +
  scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
  labs(x = "n (observations)", y = "time (s)", colour = NULL, fill = NULL,
       title = "Time of the remlax solver alone",
       subtitle = "remlax: Python solver, compilation included, without start-up or export; external: total time") +
  theme_bench
save_fig(p, "fig1b_solver_time")
# same time against m: the cost of the MME engine follows m, that of the dense engine follows n
if (!is.null(dims)) {
  # m varies slightly from one replicate to the next (unobserved levels): median m per (family, n)
  sm <- merge(summarise_by(tj, "seconds"), stats::aggregate(m ~ family + n, dims, stats::median), by = c("family", "n"))
  p <- ggplot(sm, aes(m, med, colour = route, fill = route)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point(size = 1) +
    scale_x_log10() + scale_y_log10() + facet_wrap(~ family, labeller = lab_fam, nrow = 3, scales = "free") +
    scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
    labs(x = "m (size of C: fixed effects + random levels)", y = "total time (s)", colour = NULL, fill = NULL,
         title = "Fitting time against the size of the mixed model equations") + theme_bench
  save_fig(p, "fig1c_time_vs_m")
}

# ratio of two routes within each cell, then summary over replicates
ratio <- function(num, den) {
  a <- tj[tj$route == num, c("family", "n", "rep", "cell", "seconds")]
  b <- tj[tj$route == den, c("cell", "seconds")]
  m <- merge(a, b, by = "cell", suffixes = c("_num", "_den"))
  if (!nrow(m)) return(NULL)
  m$r <- m$seconds_num / m$seconds_den
  s <- stats::aggregate(r ~ family + n, m, function(x) c(med = stats::median(x), lo = min(x), hi = max(x)))
  cbind(s[c("family", "n")], as.data.frame(s$r), comparison = sprintf("%s / %s", num, den))
}
plot_ratios <- function(r, title, sub, name) {
  r <- do.call(rbind, r)
  if (is.null(r) || !nrow(r)) return(invisible())
  p <- ggplot(r, aes(n, med, colour = comparison, fill = comparison)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point(size = 1) +
    scale_x_log10() + scale_y_log10() + facet_wrap(~ family, labeller = lab_fam, nrow = 3) +
    labs(x = "n (observations)", y = "time ratio", colour = NULL, fill = NULL, title = title, subtitle = sub) +
    theme_bench
  save_fig(p, name)
  utils::write.csv(r, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
}
plot_ratios(list(ratio("D-CPU", "D-GPU"), ratio("M-dense-CPU", "M-dense-GPU")),
            "GPU gain", "CPU time / GPU time, same node; above 1 the GPU wins", "fig2_gpu_vs_cpu")
plot_ratios(list(ratio("D-CPU", "M-dense-CPU"), ratio("D-CPU", "M-sparse-CPU"), ratio("D-GPU", "M-dense-GPU")),
            "MME engine against dense engine", "dense time / MME time; above 1 the MME wins", "fig3_mme_vs_dense")
plot_ratios(list(ratio("M-dense-CPU", "M-sparse-CPU")),
            "Inversion of C in the MME engine", "dense inversion time / sparse inversion time (CPU)", "fig4_inversion")
plot_ratios(list(ratio("asreml", "D-CPU"), ratio("asreml", "M-sparse-CPU"), ratio("asreml", "M-dense-GPU"),
                 ratio("asreml", "D-GPU")),
            "remlax against asreml", "asreml time / remlax time; above 1 remlax wins", "fig5_remlax_vs_asreml")

# ------------------------------------------------------ fastest route per cell
fastest <- do.call(rbind, lapply(split(tj, list(tj$family, tj$n), drop = TRUE), function(d) {
  s <- tapply(d$seconds, droplevels(d$route), stats::median)
  data.frame(family = d$family[1], n = d$n[1], route = names(s)[which.min(s)], seconds = min(s))
}))
fastest$route <- factor(fastest$route, ROUTES)
p <- ggplot(fastest, aes(factor(n), family, fill = route)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.3g", seconds)), size = 2.6) + scale_fill_manual(values = COL, drop = FALSE) +
  labs(x = "n (observations)", y = NULL, fill = NULL, title = "Fastest route per cell",
       subtitle = "median of the kept replicates; label = time (s)") + theme_bench
save_fig(p, "fig6_fastest_route", 10, 5)
fastest_rx <- do.call(rbind, lapply(split(tj[tj$route %in% ROUTES[1:5], ], list(tj$family[tj$route %in% ROUTES[1:5]], tj$n[tj$route %in% ROUTES[1:5]]), drop = TRUE), function(d) {
  s <- tapply(d$seconds, droplevels(d$route), stats::median)
  data.frame(family = d$family[1], n = d$n[1], route = names(s)[which.min(s)], seconds = min(s))
}))
fastest_rx$route <- factor(fastest_rx$route, ROUTES)
p <- ggplot(fastest_rx, aes(factor(n), family, fill = route)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.3g", seconds)), size = 2.6) + scale_fill_manual(values = COL[1:5], drop = FALSE) +
  labs(x = "n (observations)", y = NULL, fill = NULL, title = "Fastest remlax route per cell",
       subtitle = "basis of the engine = \"auto\" rule") + theme_bench
save_fig(p, "fig7_fastest_remlax_route", 10, 5)

# ------------------------------------------------------------ logLik accuracy
ok <- fits[fits$success, ]
p <- ggplot(ok, aes(factor(n), pmax(diff_logLik, 1e-10), colour = route)) +
  geom_hline(yintercept = TOL, linetype = 2) +
  geom_point(position = position_dodge(width = 0.7), size = 0.9, alpha = 0.8) +
  scale_y_log10() + facet_wrap(~ family, labeller = lab_fam, nrow = 3) + scale_colour_manual(values = COL) +
  labs(x = "n (observations)", y = "difference to the best logLik (arbiter)", colour = NULL,
       title = "Accuracy of the optimum", subtitle = sprintf("differences floored at 1e-10; dashes: tolerance %g", TOL)) +
  theme_bench
save_fig(p, "fig8_logLik_diff")

# ------------------------------------------------------ parameter recovery
# remlax: the routes give the same optimum; keep D-CPU (or the first accurate one)
pj <- merge(pa, fits[, c("cell", "route", "kept")], by = c("cell", "route"))
pj <- pj[pj$kept, ]
pj$software <- ifelse(pj$route %in% ROUTES[1:5], "remlax", as.character(pj$route))
pr <- pj[!duplicated(pj[c("cell", "software", "parameter")]), ]
pr$error <- pr$estimate - pr$true
rec <- stats::aggregate(error ~ family + n + software + parameter + true, pr,
                        function(x) c(bias = mean(x), rmse = sqrt(mean(x^2)), k = length(x)))
rec <- cbind(rec[c("family", "n", "software", "parameter", "true")], as.data.frame(rec$error))
utils::write.csv(rec, file.path(out_dir, "parameter_recovery.csv"), row.names = FALSE)
# one value per family: median relative RMSE over the variance parameters
rec$rrmse <- rec$rmse / pmax(abs(rec$true), 0.05)
rf <- stats::aggregate(rrmse ~ family + n + software, rec, stats::median)
COL_S <- c(remlax = "#1f78b4", asreml = "#ff7f00", sommer = "#6a3d9a", lme4 = "#b15928")
p <- ggplot(rf, aes(n, rrmse, colour = software)) + geom_line() + geom_point(size = 1) +
  scale_x_log10() + scale_y_log10() + facet_wrap(~ family, labeller = lab_fam, nrow = 3) + scale_colour_manual(values = COL_S) +
  labs(x = "n (observations)", y = "relative RMSE (median over the parameters)", colour = NULL,
       title = "Recovery of the simulated parameters",
       subtitle = "RMSE over 5 replicates divided by max(|true|, 0.05); overlapping curves indicate the same estimator") +
  theme_bench
save_fig(p, "fig9_recovery")
p <- ggplot(pr, aes(factor(n), error, fill = software)) + geom_hline(yintercept = 0) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.25) + facet_wrap(~ family + parameter, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = COL_S) +
  labs(x = "n", y = "estimate - true", fill = NULL, title = "Estimation error by parameter") +
  theme_bench + theme(strip.text = element_text(size = 6), axis.text = element_text(size = 5))
save_fig(p, "fig10_error_by_parameter", 16, 22)

# --------------------------------------------------- agreement between software
rx <- pr[pr$software == "remlax", c("cell", "parameter", "estimate")]
ac <- merge(pr[pr$software != "remlax", ], rx, by = c("cell", "parameter"), suffixes = c("", "_remlax"))
if (nrow(ac)) {
  ac$diff <- abs(ac$estimate - ac$estimate_remlax) / pmax(abs(ac$estimate_remlax), 0.05)
  ad <- stats::aggregate(diff ~ family + n + rep + software, ac, max)
  p <- ggplot(ad, aes(factor(n), pmax(diff, 1e-8), colour = software)) +
    geom_point(position = position_dodge(width = 0.6), size = 1) + scale_y_log10() +
    facet_wrap(~ family, labeller = lab_fam, nrow = 3) + scale_colour_manual(values = COL_S) +
    labs(x = "n (observations)", y = "maximum relative difference to remlax", colour = NULL,
         title = "Agreement of the estimates with remlax",
         subtitle = "per cell, max over the parameters of |estimate - remlax| / max(|remlax|, 0.05)") + theme_bench
  save_fig(p, "fig11_software_agreement")
  utils::write.csv(ad, file.path(out_dir, "software_agreement.csv"), row.names = FALSE)
}

# --------------------------------------------- compilation against computation (remlax)
rc <- tj[tj$route %in% ROUTES[1:5] & is.finite(tj$compile_s), ]
if (nrow(rc)) {
  rc$compile_share <- rc$compile_s / rc$seconds
  s <- summarise_by(rc, "compile_share")
  p <- ggplot(s, aes(n, med, colour = route)) + geom_line() + geom_point(size = 1) + scale_x_log10() +
    facet_wrap(~ family, labeller = lab_fam, nrow = 3) + scale_colour_manual(values = COL) +
    labs(x = "n (observations)", y = "share of compilation in the total time", colour = NULL,
         title = "Weight of JAX compilation") + theme_bench
  save_fig(p, "fig12_compilation")
}
cat(sprintf("%d fits, %d accurate, %d kept; figures in %s\n", nrow(fits), sum(fits$accurate), sum(fits$kept), out_dir))
