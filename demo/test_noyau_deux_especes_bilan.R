# ==============================================================================
# Summary of test_noyau_deux_especes.R (sourced by its "summary" mode, or alone:
#   Rscript test_noyau_deux_especes_bilan.R <folder>)
#
# Estimation of lambda and d: profile of the logLik on the grid, in d (maximum
# over lambda) and in lambda (maximum over d), interpolated by a spline; 95 %
# interval = drop of qchisq(0.95, 1) / 2. Grid "r3" (true radius) and "R6"
# (wide window). Accuracy and variance error: against the simulated truth,
# species A and B stacked (same parameters).
#
# Reads the task files task_*.csv, or the tache_*.csv files of the earlier
# French version, whose columns and values are mapped to English on reading.
# ==============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
if (!exists("dir") || !is.character(dir)) dir <- commandArgs(trailingOnly = TRUE)[1]
if (!exists("truth")) truth <- list(r = 3, lambda = 1.5, d = 0.5)

fr_cols <- c(graine = "seed", grille = "grid", fenetre = "window", secondes = "seconds", erreur = "error")
to_english <- function(x) {
  nm <- names(x)
  nm <- ifelse(nm %in% names(fr_cols), fr_cols[nm], nm)
  names(x) <- sub("^prec_", "acc_", nm)
  if ("grid" %in% names(x)) x$grid[x$grid == "naif"] <- "naive"
  x
}
fs <- list.files(dir, "^task_.*\\.csv$", full.names = TRUE)
if (!length(fs)) fs <- list.files(dir, "^tache_.*\\.csv$", full.names = TRUE)
res <- do.call(rbind, lapply(fs, function(f) {
  x <- to_english(utils::read.csv(f, stringsAsFactors = FALSE)); x$error <- as.character(x$error); x }))
out_dir <- file.path(dir, "summary"); dir.create(out_dir, showWarnings = FALSE)
log_con <- file(file.path(out_dir, "log.txt"), "w")
say <- function(...) { txt <- sprintf(...); cat(txt, "\n"); cat(txt, "\n", file = log_con) }
say("%d data sets, %d fits, %d failure(s), %d not converged",
    length(fs), nrow(res), sum(is.na(res$logLik)), sum(res$conv %in% FALSE))
say("Median time per fit (s) by n_rep: %s",
    paste(sprintf("%d: %.0f", sort(unique(res$n_rep)),
                  tapply(res$seconds, res$n_rep, stats::median, na.rm = TRUE)), collapse = ", "))

# --- long format: one row per (fit, species, effect) --------------------------
eff <- c(dge = "DGE", intra = "IGE intra", inter = "IGE inter", tbv = "TBV intra")
long <- do.call(rbind, lapply(c("A", "B"), function(sp) do.call(rbind, lapply(names(eff), function(e) {
  data.frame(res[, c("n_rep", "seed", "grid", "r", "lambda", "d", "logLik", "conv")], species = sp,
             effect = eff[[e]], acc = res[[sprintf("acc_%s_%s", e, sp)]],
             err = if (e == "tbv") NA_real_ else res[[sprintf("err_%s_%s", e, sp)]])
}))))
long$effect <- factor(long$effect, levels = eff)

# --- 1. estimation of lambda and d by profile ---------------------------------
profile1 <- function(x, ll, bound) {
  o <- order(x); x <- x[o]; ll <- ll[o] - max(ll)
  fine <- seq(min(x), max(x), length.out = 801)
  yf <- stats::spline(x, ll, xout = fine, method = "natural")$y
  inside <- fine[yf >= max(yf) - stats::qchisq(0.95, 1) / 2]
  c(est = fine[which.max(yf)], low = min(inside), high = max(inside),
    at_edge = as.numeric(fine[which.max(yf)] %in% range(x)))
}
est <- do.call(rbind, lapply(split(res[res$grid %in% c("r3", "R6") & !is.na(res$logLik), ],
                                   ~ n_rep + seed + grid, drop = TRUE), function(q) {
  pd <- stats::aggregate(logLik ~ d, q, max); pl <- stats::aggregate(logLik ~ lambda, q, max)
  ed <- profile1(pd$d, pd$logLik); el <- profile1(pl$lambda, pl$logLik)
  j <- which.max(q$logLik)
  data.frame(n_rep = q$n_rep[1], seed = q$seed[1], grid = q$grid[1],
             d_hat = ed[["est"]], d_low = ed[["low"]], d_high = ed[["high"]], d_at_edge = ed[["at_edge"]],
             l_hat = el[["est"]], l_low = el[["low"]], l_high = el[["high"]], l_at_edge = el[["at_edge"]],
             l_grid = q$lambda[j], d_grid = q$d[j], logLik_max = q$logLik[j])
}))
# logLik at the true kernel and at the maximum of R6, to compare the two windows
trueLL <- res[res$grid == "r3" & res$lambda == truth$lambda & res$d == truth$d, c("n_rep", "seed", "logLik")]
names(trueLL)[3] <- "logLik_true"
est <- merge(est, trueLL, by = c("n_rep", "seed"))
est$gain_over_true <- est$logLik_max - est$logLik_true
utils::write.csv(est, file.path(out_dir, "estimates_lambda_d.csv"), row.names = FALSE)

synth <- do.call(rbind, lapply(split(est, ~ grid + n_rep, drop = TRUE), function(q) data.frame(
  grid = q$grid[1], n_rep = q$n_rep[1], datasets = nrow(q),
  d_mean = mean(q$d_hat), d_bias = mean(q$d_hat) - truth$d,
  d_rmse = sqrt(mean((q$d_hat - truth$d)^2)), d_ci_width = mean(q$d_high - q$d_low),
  d_coverage = mean(q$d_low <= truth$d & q$d_high >= truth$d), d0_excluded = mean(q$d_low > 0),
  d_at_edge = mean(q$d_at_edge == 1),
  l_mean = mean(q$l_hat), l_bias = mean(q$l_hat) - truth$lambda,
  l_rmse = sqrt(mean((q$l_hat - truth$lambda)^2)), l_ci_width = mean(q$l_high - q$l_low),
  l_coverage = mean(q$l_low <= truth$lambda & q$l_high >= truth$lambda), l_at_edge = mean(q$l_at_edge == 1),
  gain_logLik = stats::median(q$gain_over_true))))
synth <- synth[order(synth$grid != "r3", synth$n_rep), ]
say("\n1. Estimation of lambda (true %.2f) and d (true %.2f) by profile", truth$lambda, truth$d)
captured <- utils::capture.output(print(synth, digits = 3, row.names = FALSE))
say("%s", paste(captured, collapse = "\n"))
utils::write.csv(synth, file.path(out_dir, "synthesis_lambda_d.csv"), row.names = FALSE)

# --- 2. accuracy of the IGE by fitted kernel -----------------------------------
# Estimated kernel: highest grid point of r3 (and of R6)
sel_est <- merge(res[res$grid %in% c("r3", "R6"), c("n_rep", "seed", "grid", "lambda", "d")],
                 est[, c("n_rep", "seed", "grid", "l_grid", "d_grid")])
sel_est <- sel_est[sel_est$lambda == sel_est$l_grid & sel_est$d == sel_est$d_grid, ]
key <- function(x) paste(x$n_rep, x$seed, x$grid, x$r, x$lambda, x$d)
L_est <- long[key(long) %in% key(merge(sel_est, unique(res[, c("n_rep", "seed", "grid", "r", "lambda", "d")]))), ]
L_est$model <- ifelse(L_est$grid == "r3", "r = 3, lambda and d estimated", "R = 6, lambda and d estimated")
L_true <- long[long$grid == "r3" & long$lambda == truth$lambda & long$d == truth$d, ]
L_true$model <- "true kernel"
L_naive <- long[long$grid == "naive", ]; L_naive$model <- "naive (rank 1, equal weights)"
L_unif <- long[long$grid == "r3" & long$lambda == 1 & long$d == 0, ]
L_unif$model <- "r = 3, lambda = 1, d = 0 (chapter 3)"
cmp <- rbind(L_true, L_est, L_unif, L_naive)
cmp$model <- factor(cmp$model, levels = c("true kernel", "r = 3, lambda and d estimated",
                                          "R = 6, lambda and d estimated", "r = 3, lambda = 1, d = 0 (chapter 3)",
                                          "naive (rank 1, equal weights)"))
tab_cmp <- stats::aggregate(cbind(acc, err) ~ model + effect + n_rep, cmp, stats::median, na.action = na.pass)
say("\n2. Median accuracy and median log10 error by fitted kernel")
say("%s", paste(utils::capture.output(print(tab_cmp[order(tab_cmp$n_rep, tab_cmp$effect, tab_cmp$model), ],
                                            digits = 3, row.names = FALSE)), collapse = "\n"))
utils::write.csv(tab_cmp, file.path(out_dir, "accuracy_by_model.csv"), row.names = FALSE)

# --- 3. misspecification curves ------------------------------------------------
L_l <- long[long$grid == "r3" & long$d == truth$d, ]; L_l$axis <- "lambda (r = 3, d = 0.5)"; L_l$x <- L_l$lambda
L_d <- long[long$grid == "r3" & long$lambda == truth$lambda, ]; L_d$axis <- "d (r = 3, lambda = 1.5)"; L_d$x <- L_d$d
L_r <- long[(long$grid %in% c("r", "R6") & long$lambda == truth$lambda & long$d == truth$d) |
            (long$grid == "r3" & long$lambda == truth$lambda & long$d == truth$d), ]
L_r$axis <- "r (lambda = 1.5, d = 0.5)"; L_r$x <- L_r$r
mis <- rbind(L_l, L_d, L_r)
mis$axis <- factor(mis$axis, levels = c("lambda (r = 3, d = 0.5)", "d (r = 3, lambda = 1.5)", "r (lambda = 1.5, d = 0.5)"))
q3 <- function(v) c(med = stats::median(v, na.rm = TRUE), q1 = unname(stats::quantile(v, 0.25, na.rm = TRUE)),
                    q3 = unname(stats::quantile(v, 0.75, na.rm = TRUE)))
agg <- function(var) {
  m <- data.frame(mis[, c("axis", "effect", "n_rep")], xv = mis$x, val = mis[[var]])
  a <- stats::aggregate(val ~ axis + xv + effect + n_rep, m, q3)
  data.frame(a[, c("axis", "effect", "n_rep")], x = a$xv, med = a$val[, "med"], q1 = a$val[, "q1"],
             q3 = a$val[, "q3"])
}
mp <- agg("acc"); me <- agg("err"); me <- me[me$effect != "TBV intra", ]
utils::write.csv(merge(mp, me, by = c("axis", "x", "effect", "n_rep"), all.x = TRUE, suffixes = c("_acc", "_err")),
                 file.path(out_dir, "misspecification.csv"), row.names = FALSE)
vx <- data.frame(axis = levels(mis$axis), x = c(truth$lambda, truth$d, truth$r))
vx$axis <- factor(vx$axis, levels = levels(mis$axis))

theme_set(theme_minimal(base_size = 10) +
          theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                plot.title = element_text(face = "bold", size = 10), strip.text = element_text(size = 9)))
col_rep <- c("5" = "#9ecae1", "10" = "#4292c6", "20" = "#08306b")

fig_mis <- function(m, ylab, title, zero = NULL) {
  p <- ggplot(m, aes(x, med, colour = factor(n_rep), fill = factor(n_rep))) +
    geom_vline(data = vx, aes(xintercept = x), linetype = 2, colour = "grey55") +
    geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.15, colour = NA) +
    geom_line() + geom_point(size = 1.2) +
    facet_grid(effect ~ axis, scales = "free") +
    scale_colour_manual(values = col_rep, name = "replicates") +
    scale_fill_manual(values = col_rep, name = "replicates") +
    labs(title = title, x = NULL, y = ylab)
  if (!is.null(zero)) p <- p + geom_hline(yintercept = zero, colour = "grey40", linewidth = 0.3)
  p
}
ggsave(file.path(out_dir, "fig1_accuracy_wrong_kernel.png"),
       fig_mis(mp, "correlation BLUP / true value (median and interquartile range)",
               "Accuracy of the genetic values when one kernel parameter is wrong (dashes: true value)"),
       width = 10, height = 8.5, dpi = 150, bg = "white")
ggsave(file.path(out_dir, "fig2_variance_error_wrong_kernel.png"),
       fig_mis(me, "log10(estimated V x k / true V x k)",
               "Error on the variances, phenotypic scale, when one kernel parameter is wrong", zero = 0),
       width = 10, height = 6.8, dpi = 150, bg = "white")

# --- 4. estimates of lambda and d ----------------------------------------------
est$window <- factor(ifelse(est$grid == "r3", "r = 3 (true radius)", "R = 6 (wide window)"),
                     levels = c("r = 3 (true radius)", "R = 6 (wide window)"))
pe <- function(v, tv, name) ggplot(est, aes(factor(n_rep), .data[[v]], fill = window)) +
  geom_hline(yintercept = tv, linetype = 2, colour = "grey45") +
  geom_boxplot(outlier.size = 0.6, width = 0.6, position = position_dodge(0.7)) +
  scale_fill_manual(values = c("#41ab5d", "#fd8d3c"), name = NULL) +
  labs(x = "replicates per genotype", y = name)
ggsave(file.path(out_dir, "fig3_estimates_lambda_d.png"),
       (pe("l_hat", truth$lambda, "estimated lambda") + ggtitle("a  lambda (dashes: true 1.5)")) +
       (pe("d_hat", truth$d, "estimated d") + ggtitle("b  d (dashes: true 0.5)")) +
       plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 10, height = 4.5, dpi = 150, bg = "white")

# --- 5. lambda x d maps: mean relative logLik, IGE accuracy and error ----------
r3 <- res[res$grid == "r3" & !is.na(res$logLik), ]
r3$dll <- r3$logLik - stats::ave(r3$logLik, r3$n_rep, r3$seed, FUN = max)
r3$acc_ige <- (r3$acc_intra_A + r3$acc_intra_B + r3$acc_inter_A + r3$acc_inter_B) / 4
r3$err_ige <- (r3$err_intra_A + r3$err_intra_B + r3$err_inter_A + r3$err_inter_B) / 4
map_df <- stats::aggregate(cbind(dll, acc_ige, err_ige) ~ lambda + d + n_rep, r3, stats::median)
ridge <- data.frame(lambda = seq(0.5, 3, 0.05)); ridge$d <- truth$d + (truth$lambda - ridge$lambda) / 2
ridge <- ridge[ridge$d >= 0 & ridge$d <= 1, ]
pc <- function(v, name, pal, lim = NULL) ggplot(map_df, aes(lambda, d, fill = .data[[v]])) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", .data[[v]])), size = 2.4) +
  geom_line(data = ridge, aes(lambda, d), inherit.aes = FALSE, linetype = 3) +
  annotate("point", x = truth$lambda, y = truth$d, shape = 4, size = 3) +
  facet_wrap(~ factor(sprintf("%d replicates", n_rep), levels = sprintf("%d replicates", sort(unique(n_rep)))), nrow = 1) +
  (if (pal == "RdBu") scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = name)
   else scale_fill_distiller(palette = pal, direction = 1, limits = lim, name = name)) +
  labs(x = "lambda", y = "d") + theme(panel.grid = element_blank())
ggsave(file.path(out_dir, "fig4_maps_lambda_d.png"),
       (pc("dll", "delta logLik", "Greens", c(-15, 0)) + ggtitle("a  logLik relative to the maximum of the data set (median)")) /
       (pc("acc_ige", "accuracy", "Blues") + ggtitle("b  mean accuracy of intra and inter IGE")) /
       (pc("err_ige", "log10 error", "RdBu") + ggtitle("c  mean log10 error of the IGE x k variances")) +
       plot_annotation(caption = "cross: true kernel; dotted: d + lambda / 2 = constant (disc approximation)"),
       width = 11, height = 10.5, dpi = 150, bg = "white")
close(log_con)
cat("Outputs in", normalizePath(out_dir), "\n")
