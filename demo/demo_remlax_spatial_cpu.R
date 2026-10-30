# ==============================================================================
# remlax demonstration on CPU: a simulated field trial with a spatial trend
# smoothed by 2D P-splines, direct genetic values (DGE) and indirect genetic
# effects (IGE) between neighbouring plots, with genomic kinship.
#
# Design: 150 genotypes, 3 replicates, incomplete blocks of 5 plots
# (resolvable alpha-type design), field of 30 rows x 15 columns.
# Each plot receives the IGE of the genotypes of its 4 immediate neighbours
# (same row or same column), as in the competition between adjacent plots of a
# variety trial. The neighbourhood incidence is normalised by row (L2 norm),
# and the IGE variance is set by its share of phenotypic variance, through the
# exposure k of rx_exposure().
#
# Every model carries rep (fixed), block and the surface spl2d(col, row).
# They differ by their genetic part:
#   A : DGE only, vm(gid, K)
#   B : independent DGE + IGE, vm(gid, K) + mm(Zn, K)
#   C : correlated DGE and IGE, str(~ gid + mm(Zn), struct = "us", K)
#
# Second part (section 6): the same trial, with 10 % empty plots and an
# unknown neighbourhood. The weights are delta^-lambda within a radius of r
# grid steps, divided by n_i^d (dilution, n_i = number of neighbours). The
# second part recovers r, lambda and d by profiling the REML likelihood (coarse
# grid over r, lambda, d, then fine grid over lambda, d at the best r), then by
# direct maximisation of this profile in lambda and d. It finally fits model C
# with the estimated neighbourhood.
#
# Run from a folder without renv (remlax is in the user's R library):
#   Rscript demo_remlax_spatial_cpu.R [seed] [output_folder] [part2 TRUE/FALSE]
# Outputs: four PNG figures and three CSV tables in the output folder.
# ==============================================================================

suppressPackageStartupMessages({
  library(remlax)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
seed <- if (length(args) >= 1) as.integer(args[1]) else 5L
out_dir <- if (length(args) >= 2) args[2] else "output_demo_ige"
part2 <- if (length(args) >= 3) as.logical(args[3]) else TRUE
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

# ------------------------------------------------------------------------------
# 1. True values
# ------------------------------------------------------------------------------
# `ige` is the share of phenotypic variance carried by the IGE (sigma2_IGE x k),
# `r_dge_ige` the genetic correlation between direct and indirect effects. It
# is negative: a vigorous genotype depresses its neighbours.
truth <- list(dge = 1.0, ige = 0.4, r_dge_ige = -0.5, block = 0.3, spatial = 1.2, residual = 0.8)
n_gen <- 150L; n_rep <- 3L; k_block <- 5L
n_col <- 15L; n_row <- 30L                      # each replicate: 10 rows
rows_per_rep <- n_row / n_rep

# Markers: two subpopulations (Balding-Nichols drift, Fst = 0.1), then VanRaden
# genomic kinship.
n_snp <- 1500L
p0 <- runif(n_snp, 0.05, 0.5)
pop <- rep(1:2, length.out = n_gen)
P <- sapply(p0, function(p) rbeta(2, p * 9, (1 - p) * 9))
M <- matrix(rbinom(n_gen * n_snp, 2L, P[pop, ]), n_gen, n_snp)
gen_names <- sprintf("G%03d", seq_len(n_gen))
rownames(M) <- gen_names
K <- rx_grm(M, blend = 0.02, coding = "count")

# True spatial trend: a fertility gradient along the rows, a localised wet
# depression and a ripple between columns. It is centred and scaled on the
# plots so that its variance equals truth$spatial. The same transformation
# applies on the fine grid.
raw_trend <- function(x, y)
  0.9 * (y - 15.5) / 14.5 -
  1.6 * exp(-((x - 11)^2 / 10 + (y - 9)^2 / 30)) +
  0.35 * sin(2 * pi * x / 10)

# ------------------------------------------------------------------------------
# 2. Design and neighbourhood
# ------------------------------------------------------------------------------
d <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))
d$rep <- factor(sprintf("R%d", (d$row - 1L) %/% rows_per_rep + 1L))
d$block <- factor(sprintf("%s_r%02d_b%d", d$rep, d$row, (d$col - 1L) %/% k_block + 1L))
d$plot <- sprintf("P%03d", seq_len(nrow(d)))
# Each genotype once per replicate, drawn at random within the replicate. A
# real alpha design also optimises the encounters between genotypes. This draw
# keeps the structure (complete replicates, incomplete blocks) without that
# optimisation.
d$gid <- factor(NA, levels = gen_names)
for (r in levels(d$rep)) d$gid[d$rep == r] <- sample(gen_names)

# Neighbours: the plots one grid step away in Euclidean distance, hence the 4
# row and column neighbours (3 or 2 on the edges). The replicates touch in the
# field, so neighbours from another replicate count too.
nb <- rx_neighbourhood(coord = cbind(d$row, d$col), group = rep("A", nrow(d)), level = d$gid,
                       id = d$plot, rank = 1, window = "euclidean", kernel = "none",
                       normalise = TRUE)
Zn <- nb$level[["A<-A"]][, gen_names]           # plots x neighbouring genotypes
Zd <- Matrix::sparseMatrix(i = seq_len(nrow(d)), j = as.integer(d$gid), x = 1,
                           dims = c(nrow(d), n_gen), dimnames = list(NULL, gen_names))
# Exposure: the phenotypic variance brought by a term equals sigma2 x k.
expo <- rx_exposure(Zd, Zn, K = K)
k_ige <- expo$k

# True genetic values: (DGE, IGE) ~ N(0, Sigma (x) K)
s2_ige <- truth$ige / k_ige
Sigma_true <- matrix(c(truth$dge, rep(truth$r_dge_ige * sqrt(truth$dge * s2_ige), 2), s2_ige), 2)
U <- t(chol(K)) %*% matrix(rnorm(2 * n_gen), n_gen, 2) %*% chol(Sigma_true)
dge_true <- setNames(U[, 1], gen_names)
ige_true <- setNames(U[, 2], gen_names)
# Total value of a genotype in pure stand: its direct effect, plus the effect
# it exerts on every plot that sees it. s = total weight emitted per plot, on
# average (about 2 with 4 neighbours normalised in L2).
s_emitted <- sum(Zn) / nrow(d)
tbv_true <- dge_true + s_emitted * ige_true

br <- raw_trend(d$col, d$row)
m_br <- mean(br); s_br <- sd(br)
true_trend <- function(x, y) (raw_trend(x, y) - m_br) / s_br * sqrt(truth$spatial)

d$spatial_true <- true_trend(d$col, d$row)
d$g_true <- dge_true[as.character(d$gid)]
d$ige_received <- as.numeric(Zn %*% ige_true)
b_true <- setNames(rnorm(nlevels(d$block), 0, sqrt(truth$block)), levels(d$block))
d$block_true <- b_true[as.character(d$block)]
d$e_true <- rnorm(nrow(d), 0, sqrt(truth$residual))
d$y <- 6 + 0.2 * (d$rep == "R2") - 0.1 * (d$rep == "R3") +
  d$spatial_true + d$g_true + d$ige_received + d$block_true + d$e_true

cat(sprintf("Design: %d plots, %d genotypes, %d replicates, %d blocks of %d\n",
            nrow(d), n_gen, n_rep, nlevels(d$block), k_block))
cat(sprintf("Neighbourhood: %s neighbours per plot; exposure k = %.3f; emitted weight s = %.2f\n",
            paste(names(table(nb$n_neighbours[[1]])), collapse = "/"), k_ige, s_emitted))

# ------------------------------------------------------------------------------
# 3. Fits
# ------------------------------------------------------------------------------
fit_model <- function(genetic) {
  t0 <- Sys.time()
  random <- stats::as.formula(paste("~", genetic, "+ block + spl2d(col, row, nseg = c(7, 10))"))
  f <- rx_reml(y ~ rep, random = random, data = d, backend = "cpu", verbose = FALSE)
  f$duration <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  f
}
fits <- list(
  "A  DGE only"              = fit_model("vm(gid, K = K)"),
  "B  DGE + IGE"             = fit_model("vm(gid, K = K) + mm(Zn, K = K, name = 'ige')"),
  "C  DGE + IGE correlated"  = fit_model("str(~ gid + mm(Zn, name = 'ige'), struct = 'us', K = K, name = 'dge_ige')"))
model_names <- names(fits)
colours <- setNames(c("#c2a5cf", "#999999", "#1b7837"), model_names)
fC <- fits[[3]]
print(summary(fC))
# The three criteria behind the Convergence line of summary()
for (m in model_names) with(fits[[m]], cat(sprintf(
  "%-24s criteria: decrement %s, gradient %s, Hessian %s (%s negative eigenvalue(s), %s flat direction(s)); %d variance(s) at the bound\n",
  m, conv_decrement, conv_grad_rel, conv_hessian_ok, n_neg_eig, n_null_dir %||% NA, n_at_bound)))

# ------------------------------------------------------------------------------
# 4. What the inference recovers
# ------------------------------------------------------------------------------
# BLUP per genotype: DGE and IGE (zero for A, which has no IGE).
blups_of <- function(f) {
  b <- f$blups
  if (!is.null(b$dge_ige)) return(list(dge = b$dge_ige[gen_names, 1], ige = b$dge_ige[gen_names, 2]))
  list(dge = setNames(as.numeric(b$gid), rownames(b$gid))[gen_names],
       ige = if (is.null(b$ige)) setNames(rep(0, n_gen), gen_names)
             else setNames(as.numeric(b$ige), rownames(b$ige))[gen_names])
}
g_est <- lapply(fits, blups_of)
top <- function(v, p = 0.1) names(sort(v, decreasing = TRUE))[seq_len(round(p * length(v)))]

# Genetic components on the phenotypic scale: sigma2_DGE, sigma2_IGE x k, and
# the correlation r(DGE, IGE).
genetic_components <- function(f) {
  v <- summary(f)$varcomp
  val <- function(nm) if (nm %in% rownames(v)) v[nm, "component"] else NA_real_
  se_of <- function(nm) if (nm %in% rownames(v)) v[nm, "std.error"] else NA_real_
  if ("dge_ige[1,1]" %in% rownames(v)) {
    s11 <- val("dge_ige[1,1]"); s21 <- val("dge_ige[2,1]"); s22 <- val("dge_ige[2,2]")
    return(data.frame(component = c("DGE", "IGE x k"), estimate = c(s11, s22 * k_ige),
                      se = c(se_of("dge_ige[1,1]"), se_of("dge_ige[2,2]") * k_ige),
                      r = s21 / sqrt(s11 * s22)))
  }
  data.frame(component = c("DGE", "IGE x k"), estimate = c(val("gid"), val("ige") * k_ige),
             se = c(se_of("gid"), se_of("ige") * k_ige), r = NA_real_)
}

# Environment estimated by C: spline surface + replicate effects, centred on
# the plots. The replicates absorb part of the gradient along the rows, hence
# their addition before the comparison with the true trend.
cf <- summary(fC)$coef.fixed
beta_rep <- c(R1 = 0, R2 = cf["repR2", "solution"], R3 = cf["repR3", "solution"])
trend_C <- function(nx, ny) {
  tr <- rx_spatial_trend(fC, grid = c(nx, ny))
  tr$trend <- tr$trend + beta_rep[pmin(n_rep, floor((tr$row - 0.5) / rows_per_rep) + 1L)]
  tr
}
tr_p <- trend_C(n_col, n_row)
mean_C <- mean(tr_p$trend)
d$spatial_est <- tr_p$trend[match(paste(d$col, d$row), paste(tr_p$col, tr_p$row))] - mean_C
for (i in seq_along(fits)) d[[paste0("res", LETTERS[i])]] <- residuals(fits[[i]])

tab <- do.call(rbind, lapply(model_names, function(m) {
  f <- fits[[m]]; g <- g_est[[m]]; s <- summary(f)
  tbv <- g$dge + s_emitted * g$ige
  data.frame(model = m, logLik = f$logLik, AIC = s$aic, converged = isTRUE(s$converged),
             seconds = round(f$duration, 1),
             r_dge = cor(g$dge, dge_true),
             r_ige = if (all(g$ige == 0)) NA_real_ else cor(g$ige, ige_true),
             r_tbv = cor(tbv, tbv_true),
             selection_10pc = length(intersect(top(tbv), top(tbv_true))) / length(top(tbv_true)))
}))
write.csv(tab, file.path(out_dir, "model_comparison.csv"), row.names = FALSE)

comp <- do.call(rbind, lapply(model_names, function(m) {
  v <- summary(fits[[m]])$varcomp
  data.frame(model = m, component = rownames(v), estimate = v$component, se = v$std.error)
}))
write.csv(comp, file.path(out_dir, "variance_components.csv"), row.names = FALSE)
gen <- do.call(rbind, lapply(model_names, function(m) cbind(model = m, genetic_components(fits[[m]]))))

cat("\nModel comparison\n"); print(tab, digits = 4, row.names = FALSE)
cat("\nGenetic components, phenotypic scale (true: DGE", truth$dge, ", IGE x k", truth$ige,
    ", r", truth$r_dge_ige, ")\n"); print(gen, digits = 3, row.names = FALSE)
cat(sprintf("\nCorrelation between estimated (C) and true trend on the plots: %.3f\n",
            cor(d$spatial_est, d$spatial_true)))

# ------------------------------------------------------------------------------
# 5. Figures
# ------------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 10) +
          theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10),
                legend.key.height = unit(0.5, "cm")))
lim_sym <- function(...) { m <- max(abs(c(...))); c(-m, m) }
map_plot <- function(df, var, title, lim = NULL) {
  ggplot(df, aes(col, row)) + geom_tile(aes(fill = .data[[var]])) +
    scale_fill_distiller(palette = "RdBu", direction = 1, limits = lim, name = NULL) +
    coord_equal(expand = FALSE) + labs(title = title, x = "column", y = "row")
}
surface <- function(df, z, title, lim) {
  ggplot(df, aes(col, row)) + geom_raster(aes(fill = .data[[z]])) +
    geom_contour(aes(z = .data[[z]]), colour = "grey20", linewidth = 0.2, bins = 8) +
    scale_fill_distiller(palette = "RdBu", direction = 1, limits = lim, name = NULL) +
    coord_equal(expand = FALSE) + labs(title = title, x = "column", y = "row")
}

# Outlines of blocks and replicates
bb <- aggregate(cbind(col, row) ~ block, d, range)
blocks_df <- data.frame(xmin = bb$col[, 1] - 0.5, xmax = bb$col[, 2] + 0.5,
                        ymin = bb$row[, 1] - 0.5, ymax = bb$row[, 2] + 0.5)
reps <- data.frame(xmin = 0.5, xmax = n_col + 0.5, ymin = (seq_len(n_rep) - 1) * rows_per_rep + 0.5,
                   ymax = seq_len(n_rep) * rows_per_rep + 0.5, rep = levels(d$rep))
grid_df <- expand.grid(col = seq(1, n_col, length.out = 90), row = seq(1, n_row, length.out = 180))
grid_df$truth <- true_trend(grid_df$col, grid_df$row)

# Figure 1: simulated data. The layout also shows one plot and its neighbours.
p0_ <- d[d$col == 8 & d$row == 15, ]
nbrs <- d[abs(d$col - p0_$col) + abs(d$row - p0_$row) == 1, ]
p_layout <- ggplot() +
  geom_rect(data = reps, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = rep), alpha = 0.35) +
  geom_rect(data = blocks_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "grey30", linewidth = 0.3) +
  geom_tile(data = nbrs, aes(col, row), fill = "#d6604d", alpha = 0.8) +
  geom_tile(data = p0_, aes(col, row), fill = "black") +
  scale_fill_manual(values = c(R1 = "#8fbcd4", R2 = "#e3b778", R3 = "#b5d6a7"), name = NULL) +
  coord_equal(expand = FALSE) +
  labs(title = "a  Layout: 3 replicates, blocks of 5;\n   one plot (black) and its 4 neighbours (red)",
       x = "column", y = "row")
p_trend <- surface(grid_df, "truth", "b  True spatial trend", lim_sym(grid_df$truth))
p_g <- map_plot(d, "g_true", "c  True DGE of the plot", lim_sym(d$g_true))
p_i <- map_plot(d, "ige_received", "d  IGE received from the neighbours", lim_sym(d$ige_received))
p_y <- ggplot(d, aes(col, row, fill = y)) + geom_tile() +
  scale_fill_viridis_c(name = NULL) + coord_equal(expand = FALSE) +
  labs(title = "e  Observed phenotype", x = "column", y = "row")
sources <- c("DGE", "IGE received", "block", "spatial", "residual")
parts <- data.frame(source = factor(sources, levels = sources),
                    variance = c(var(d$g_true), var(d$ige_received), var(d$block_true), var(d$spatial_true),
                                 var(d$e_true)))
p_parts <- ggplot(parts, aes(source, variance, fill = source)) + geom_col(width = 0.6) +
  scale_fill_manual(values = c("#1b7837", "#d6604d", "#999999", "#2166ac", "#bababa"), guide = "none") +
  labs(title = "f  Realised variances on the plots", x = NULL, y = "variance") +
  theme(panel.grid.major.y = element_line(colour = "grey90"))
fig1 <- (p_layout | p_trend | p_g) / (p_i | p_y | p_parts) +
  plot_annotation(title = sprintf("Simulated data: %d plots, %d genotypes, r(DGE, IGE) = %.1f (seed %d)",
                                  nrow(d), n_gen, truth$r_dge_ige, seed))
ggsave(file.path(out_dir, "fig1_simulated_data.png"), fig1, width = 12, height = 10, dpi = 150, bg = "white")

# Figure 2: the surface recovered by C
tr_fine <- trend_C(90, 180)
tr_fine$trend <- tr_fine$trend - mean_C
L <- lim_sym(grid_df$truth, tr_fine$trend)
p_v <- surface(grid_df, "truth", "a  True trend", L)
p_e <- surface(tr_fine, "trend", "b  Trend estimated by C (spl2d + rep)", L)
p_sc <- ggplot(d, aes(spatial_true, spatial_est)) +
  geom_abline(colour = "grey60", linetype = 2) + geom_point(size = 0.9, alpha = 0.6, colour = "#2166ac") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3,
           label = sprintf("r = %.3f", cor(d$spatial_true, d$spatial_est))) +
  coord_equal() + labs(title = "c  Per plot", x = "true trend", y = "estimated trend")
Lr <- lim_sym(d$resA, d$resC)
p_rA <- map_plot(d, "resA", "d  Residuals of A (no IGE)", Lr)
p_rC <- map_plot(d, "resC", "e  Residuals of C (correlated IGE)", Lr)
dens <- data.frame(model = factor(rep(model_names, each = nrow(d)), levels = model_names),
                   res = c(d$resA, d$resB, d$resC))
p_dens <- ggplot(dens, aes(res, colour = model)) + geom_density(linewidth = 0.7) +
  scale_colour_manual(values = colours, name = NULL) +
  labs(title = "f  Distribution of residuals", x = "residual", y = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.85))
fig2 <- (p_v | p_e | p_sc) / (p_rA | p_rC | p_dens) +
  plot_annotation(title = "Spatial trend: true against estimated")
ggsave(file.path(out_dir, "fig2_surface.png"), fig2, width = 12, height = 10, dpi = 150, bg = "white")

# Figure 3: the genetic part
scatter_plot <- function(tv, field, mods, title) {
  df <- do.call(rbind, lapply(mods, function(m)
    data.frame(truth = tv, blup = g_est[[m]][[field]], model = m)))
  df$model <- factor(df$model, levels = model_names)
  ggplot(df, aes(truth, blup, colour = model)) +
    geom_abline(colour = "grey60", linetype = 2) + geom_point(size = 1, alpha = 0.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
    scale_colour_manual(values = colours, name = NULL) +
    labs(title = title, x = "true value", y = "BLUP") +
    theme(legend.position = "inside", legend.position.inside = c(0.25, 0.88))
}
p_bd <- scatter_plot(dge_true, "dge", model_names[c(1, 3)], "a  DGE: BLUP against true")
p_bi <- scatter_plot(ige_true, "ige", model_names[2:3], "b  IGE: BLUP against true")
gen$component <- factor(gen$component, levels = c("DGE", "IGE x k"))
gen$model <- factor(gen$model, levels = model_names)
vv <- data.frame(component = factor(c("DGE", "IGE x k"), levels = levels(gen$component)),
                 truth = c(truth$dge, truth$ige))
rC <- gen$r[gen$model == model_names[3]][1]
p_vc <- ggplot(gen[!is.na(gen$estimate), ], aes(component, estimate, colour = model)) +
  geom_point(data = vv, aes(component, truth), inherit.aes = FALSE, shape = 95, size = 14) +
  geom_pointrange(aes(ymin = pmax(estimate - 2 * se, 0), ymax = estimate + 2 * se),
                  position = position_dodge(width = 0.6), size = 0.3) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.5, size = 3,
           label = sprintf("r(DGE, IGE) estimated by C: %.2f (true %.2f)", rC, truth$r_dge_ige)) +
  scale_colour_manual(values = colours, name = NULL) +
  labs(title = "c  Genetic variances, phenotypic scale", subtitle = "black bar = true value",
       x = NULL, y = "variance (+/- 2 s.e.)") +
  theme(legend.position = "none", panel.grid.major.y = element_line(colour = "grey90"))
ind <- rbind(data.frame(model = tab$model, criterion = "r(DGE)", value = tab$r_dge),
             data.frame(model = tab$model, criterion = "r(IGE)", value = tab$r_ige),
             data.frame(model = tab$model, criterion = "r(total value)", value = tab$r_tbv),
             data.frame(model = tab$model, criterion = "15 best recovered", value = tab$selection_10pc))
ind$model <- factor(ind$model, levels = model_names)
ind$criterion <- factor(ind$criterion, levels = unique(ind$criterion))
p_ind <- ggplot(ind[!is.na(ind$value), ], aes(model, value, fill = model)) + geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.3, size = 2.8) +
  facet_wrap(~ criterion, nrow = 1) + scale_fill_manual(values = colours, guide = "none") +
  scale_y_continuous(limits = c(0, 1.08)) + scale_x_discrete(labels = c("A", "B", "C")) +
  labs(title = "d  Accuracy and selection on the total value DGE + s IGE", x = NULL, y = NULL)
p_aic <- ggplot(tab, aes(factor(model, levels = model_names), AIC - min(AIC), fill = model)) +
  geom_col(width = 0.65) + scale_fill_manual(values = colours, guide = "none") +
  scale_x_discrete(labels = c("A", "B", "C")) +
  labs(title = "e  AIC relative to the best model", x = NULL, y = "delta AIC") +
  theme(panel.grid.major.y = element_line(colour = "grey90"))
fig3 <- (p_bd | p_bi | p_vc) / (p_ind + p_aic + plot_layout(widths = c(3, 1))) +
  plot_annotation(title = "Genetic part: A = DGE only, B = DGE + IGE, C = correlated DGE and IGE")
ggsave(file.path(out_dir, "fig3_genetics.png"), fig3, width = 14, height = 9, dpi = 150, bg = "white")

# ==============================================================================
# 6. Second part: unknown neighbourhood (r, lambda, d)
# ==============================================================================
# Row-wise L2 normalisation would cancel the dilution: it divides each row by
# its norm, after the division by n_i^d. The weights therefore stay raw here,
# and d is identifiable only because n_i varies: field edges and empty plots.
# Each point of the profile is a complete REML fit of model B (independent DGE
# and IGE), without Hessian or BLUP for speed.
if (part2) {
  cat("\n=== Part 2: neighbourhood to estimate ===\n")
  truth2 <- list(r = 3L, lambda = 1.5, d = 0.5, dge = 1.0, ige = 0.4, r_dge_ige = -0.5,
                 block = 0.3, residual = 0.8)
  set.seed(seed + 1000L)
  d2 <- d[sort(sample(nrow(d), round(0.9 * nrow(d)))), c("col", "row", "rep", "block", "plot", "gid")]
  rownames(d2) <- NULL
  d2$block <- droplevels(d2$block)

  # Incidence plots x neighbouring genotypes for (r, lambda, d). A genotype
  # absent from the neighbourhood keeps a zero column.
  neighbourhood_Z <- function(r, lambda, dd, normalise = FALSE) {
    # a named scalar would be read as a list named by pair
    r <- as.integer(unname(r)); lambda <- as.numeric(unname(lambda)); dd <- as.numeric(unname(dd))
    nb <- rx_neighbourhood(coord = cbind(d2$row, d2$col), group = rep("A", nrow(d2)),
                           level = d2$gid, id = d2$plot, rank = r, reach = lambda,
                           dilution = dd, window = "euclidean", kernel = "power",
                           normalise = normalise)
    Z <- nb$level[["A<-A"]]
    P <- Matrix::sparseMatrix(i = match(colnames(Z), gen_names), j = seq_len(ncol(Z)), x = 1,
                              dims = c(n_gen, ncol(Z)))
    Z <- Z %*% Matrix::t(P)
    dimnames(Z) <- list(NULL, gen_names)
    methods::as(Z, "CsparseMatrix")
  }
  Zd2 <- Matrix::sparseMatrix(i = seq_len(nrow(d2)), j = as.integer(d2$gid), x = 1,
                              dims = c(nrow(d2), n_gen), dimnames = list(NULL, gen_names))
  Zv <- neighbourhood_Z(truth2$r, truth2$lambda, truth2$d)
  k2 <- rx_exposure(Zd2, Zv, K = K)$k
  s2i <- truth2$ige / k2
  Sig2 <- matrix(c(truth2$dge, rep(truth2$r_dge_ige * sqrt(truth2$dge * s2i), 2), s2i), 2)
  U2 <- t(chol(K)) %*% matrix(rnorm(2 * n_gen), n_gen, 2) %*% chol(Sig2)
  dge2 <- setNames(U2[, 1], gen_names); ige2 <- setNames(U2[, 2], gen_names)
  s2_emitted <- sum(Zv) / nrow(d2)
  tbv2 <- dge2 + s2_emitted * ige2
  b2 <- setNames(rnorm(nlevels(d2$block), 0, sqrt(truth2$block)), levels(d2$block))
  d2$ige_received <- as.numeric(Zv %*% ige2)
  d2$y <- 6 + 0.2 * (d2$rep == "R2") - 0.1 * (d2$rep == "R3") + true_trend(d2$col, d2$row) +
    dge2[as.character(d2$gid)] + d2$ige_received + b2[as.character(d2$block)] +
    rnorm(nrow(d2), 0, sqrt(truth2$residual))
  cat(sprintf("%d sown plots; true neighbourhood r = %d, lambda = %g, d = %g; k = %.3f\n",
              nrow(d2), truth2$r, truth2$lambda, truth2$d, k2))

  fit_model2 <- function(Z, correlated = FALSE, light = TRUE) {
    gen <- if (correlated) "str(~ gid + mm(Z, name = 'ige'), struct = 'us', K = K, name = 'dge_ige')"
           else "vm(gid, K = K) + mm(Z, K = K, name = 'ige')"
    random <- stats::as.formula(paste("~", gen, "+ block + spl2d(col, row, nseg = c(7, 10))"),
                                env = environment())
    t0 <- Sys.time()
    f <- if (light) rx_reml(y ~ rep, random = random, data = d2, backend = "cpu", verbose = FALSE,
                            hessian = FALSE, blups = FALSE, pev = FALSE)
         else rx_reml(y ~ rep, random = random, data = d2, backend = "cpu", verbose = FALSE)
    f$duration <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    f
  }
  profile_point <- function(r, lambda, dd) {
    f <- tryCatch(fit_model2(neighbourhood_Z(r, lambda, dd)), error = function(e) conditionMessage(e))
    if (is.character(f)) return(data.frame(r = r, lambda = lambda, d = dd, logLik = NA, conv = NA,
                                           error = substr(gsub("[\r\n]+", " ", f), 1, 200)))
    data.frame(r = r, lambda = lambda, d = dd, logLik = f$logLik,
               conv = isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement)),
               error = "")
  }
  # Parallel fits: one remlax server per process, two threads each.
  cores <- max(1L, min(8L, parallel::detectCores() %/% 2L))
  Sys.setenv(XLA_FLAGS = "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2")
  options(remlax.threads = 2L)
  in_parallel <- function(g) do.call(rbind, parallel::mclapply(seq_len(nrow(g)), function(i)
    profile_point(g$r[i], g$lambda[i], g$d[i]), mc.cores = cores))

  # Step 1: coarse grid. At r = 1, the 4 neighbours are at the same distance
  # and lambda changes nothing.
  t0 <- Sys.time()
  g1 <- rbind(data.frame(r = 1L, lambda = 0, d = c(0, 0.5, 1)),
              expand.grid(r = 2:4, lambda = c(0, 1, 2, 3), d = c(0, 0.5, 1)))
  p1 <- in_parallel(g1)
  r_hat <- p1$r[which.max(p1$logLik)]
  cat(sprintf("Coarse grid: %d fits in %.0f s; best r = %d\n", nrow(p1),
              as.numeric(difftime(Sys.time(), t0, units = "secs")), r_hat))

  # Step 2: fine grid over (lambda, d) at the best r.
  t0 <- Sys.time()
  g2 <- expand.grid(r = r_hat, lambda = if (r_hat == 1L) 0 else seq(0, 3, 0.25), d = seq(0, 1.25, 0.25))
  p2 <- in_parallel(g2)
  i2 <- which.max(p2$logLik)
  cat(sprintf("Fine grid: %d fits in %.0f s; best point lambda = %g, d = %g\n",
              nrow(p2), as.numeric(difftime(Sys.time(), t0, units = "secs")), p2$lambda[i2], p2$d[i2]))

  write.csv(rbind(cbind(step = "coarse", p1), cbind(step = "fine", p2)),
            file.path(out_dir, "neighbourhood_profile.csv"), row.names = FALSE)
  if (any(!is.na(c(p1$error, p2$error)) & nzchar(c(p1$error, p2$error))))
    print(unique(c(p1$error, p2$error)))

  # Step 3: direct maximisation of the profile in (lambda, d), from the best
  # point of the fine grid. The gradient comes from forward differences, whose
  # three fits run in parallel.
  t0 <- Sys.time()
  memo <- new.env()
  eval_fg <- function(par) {
    key <- paste(sprintf("%.6f", par), collapse = "_")
    if (!is.null(memo[[key]])) return(memo[[key]])
    h <- 0.02
    pts <- rbind(par, par + c(h, 0), par + c(0, h))
    res <- parallel::mclapply(1:3, function(i) profile_point(r_hat, pts[i, 1], pts[i, 2]),
                              mc.cores = min(3L, cores))
    ll <- vapply(res, function(x) if (is.data.frame(x)) x$logLik else NA_real_, 0)
    if (any(!is.finite(ll))) {
      msg <- vapply(res, function(x) if (is.data.frame(x)) x$error else as.character(x), "")
      stop("direct maximisation: fit failed at (", paste(round(par, 3), collapse = ", "),
           "): ", paste(unique(msg), collapse = " | "))
    }
    memo[[key]] <- list(f = -ll[1], g = -(ll[2:3] - ll[1]) / h)
    memo[[key]]
  }
  if (r_hat == 1L) {
    o <- stats::optimize(function(x) -profile_point(1L, 0, x)$logLik, c(0, 1.5))
    est <- c(lambda = 0, d = o$minimum); ll_est <- -o$objective; n_eval <- NA
  } else {
    o <- stats::optim(c(p2$lambda[i2], p2$d[i2]), function(p) eval_fg(p)$f, function(p) eval_fg(p)$g,
                      method = "L-BFGS-B", lower = c(0, 0), upper = c(4, 1.5),
                      control = list(maxit = 25, factr = 1e9))
    est <- c(lambda = o$par[1], d = o$par[2]); ll_est <- -o$value; n_eval <- length(ls(memo)) * 3L
  }
  cat(sprintf("Direct maximisation: lambda = %.2f, d = %.2f, logLik %.3f (%s fits, %.0f s)\n",
              est[["lambda"]], est[["d"]], ll_est, n_eval,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  options(remlax.threads = NULL)

  # Step 4: full model C under three neighbourhoods. The naive one is that of
  # the first part (4 neighbours, L2). The oracle knows the true neighbourhood.
  Zs <- list("naive: r = 1, L2" = neighbourhood_Z(1L, 0, 0, normalise = TRUE),
             "estimated" = neighbourhood_Z(r_hat, est[["lambda"]], est[["d"]]),
             "true (oracle)" = Zv)
  final_fits <- lapply(Zs, function(Z) fit_model2(Z, correlated = TRUE, light = FALSE))
  summary_tab <- do.call(rbind, lapply(names(Zs), function(m) {
    f <- final_fits[[m]]; Z <- Zs[[m]]; v <- summary(f)$varcomp
    k_m <- rx_exposure(Zd2, Z, K = K)$k; s_m <- sum(Z) / nrow(d2)
    b <- f$blups$dge_ige[gen_names, ]
    s11 <- v["dge_ige[1,1]", "component"]; s21 <- v["dge_ige[2,1]", "component"]
    s22 <- v["dge_ige[2,2]", "component"]
    data.frame(neighbourhood = m, logLik = f$logLik, dge = s11, ige_k = s22 * k_m,
               r_dge_ige = s21 / sqrt(s11 * s22), acc_dge = cor(b[, 1], dge2),
               acc_ige = cor(b[, 2], ige2), acc_tbv = cor(b[, 1] + s_m * b[, 2], tbv2),
               converged = isTRUE(summary(f)$converged), seconds = round(f$duration, 1))
  }))
  cat("\nModel C under three neighbourhoods (true: DGE 1, IGE x k 0.4, r -0.5)\n")
  print(summary_tab, digits = 3, row.names = FALSE)

  # Relation between lambda and d at fixed r. The variance of the IGE received
  # by plot i is proportional to ||w_i||^2 = sum_j delta_ij^(-2 lambda) / n_i^(2 d).
  # The common scale goes into sigma2_IGE: the data only see how ||w_i||
  # varies with n_i. Write sum_j delta_ij^(-2 lambda) = n_i M_i(lambda), with
  # M_i the mean of delta^(-2 lambda) over the neighbours of i, and
  #   log ||w_i|| = (1/2 - d) log n_i + (1/2) log M_i(lambda).
  # If b(lambda) is the slope of log M_i against log n_i, the exponent seen by
  # the data is e = 1/2 - d + b(lambda)/2, and the curve that keeps it constant
  # passes through the true point: d(lambda) = d_true + (b(lambda) - b(lambda_true)) / 2.
  # Continuous approximation: if n_i varies as the area of a disc of radius
  # rho_i, n_i ~ rho_i^2 and sum delta^(-2 lambda) ~ rho_i^(2 - 2 lambda), hence
  # b = -lambda and the curve d + lambda / 2 = constant. On the design, edges
  # and empty plots cut the neighbourhood like a half-plane, so b is computed.
  slope_b <- function(r, lambda) {
    nb <- rx_neighbourhood(cbind(d2$row, d2$col), group = rep("A", nrow(d2)), id = d2$plot,
                           rank = r, reach = lambda, window = "euclidean", kernel = "power",
                           output = "unit")
    n_i <- nb$n_neighbours[[1]]; s2 <- Matrix::rowSums(nb$unit[[1]]^2); ok <- n_i > 0
    if (lambda == 0) return(0)
    unname(stats::coef(stats::lm(log(s2[ok] / n_i[ok]) ~ log(n_i[ok])))[2])
  }
  lam_g <- seq(0, 3, 0.05)
  b_g <- vapply(lam_g, function(l) slope_b(r_hat, l), 0)
  b_v <- slope_b(r_hat, truth2$lambda); b_e <- slope_b(r_hat, est[["lambda"]])
  curves <- rbind(
    data.frame(lambda = lam_g, d = truth2$d + (b_g - b_v) / 2, what = "constant exponent, design"),
    data.frame(lambda = lam_g, d = truth2$d - (lam_g - truth2$lambda) / 2, what = "d + lambda/2 constant (disc)"))
  curves <- curves[curves$d >= 0 & curves$d <= 1.25, ]
  # Ridge of the profile: best d for each lambda, by a parabola through the
  # three grid points around the maximum (at the edge of the grid, the edge point).
  ridge <- do.call(rbind, lapply(split(p2, p2$lambda), function(q) {
    q <- q[order(q$d), ]; j <- which.max(q$logLik)
    dd <- q$d[j]
    if (j > 1 && j < nrow(q)) {
      y <- q$logLik[(j - 1):(j + 1)]; h <- diff(q$d)[1]
      dd <- q$d[j] + h * (y[1] - y[3]) / (2 * (y[1] - 2 * y[2] + y[3]))
    }
    data.frame(lambda = q$lambda[1], d = dd)
  }))
  e_v <- 0.5 - truth2$d + b_v / 2; e_e <- 0.5 - est[["d"]] + b_e / 2
  cat(sprintf(paste0("\nRelation lambda-d at r = %d: b(lambda) = %.3f at the truth, %.3f at the estimate;\n",
                     "  exponent e = 1/2 - d + b/2: %.3f at the truth, %.3f at the estimate;\n",
                     "  d on the curve at the estimated lambda: %.3f (design), %.3f (disc); estimated d %.3f\n"),
              r_hat, b_v, b_e, e_v, e_e, truth2$d + (b_e - b_v) / 2,
              truth2$d - (est[["lambda"]] - truth2$lambda) / 2, est[["d"]]))

  # Figure 4
  p1$dll <- pmax(p1$logLik - max(p1$logLik), -40)
  p_g <- ggplot(p1, aes(lambda, dll, colour = factor(d))) + geom_line(data = p1[p1$r > 1, ]) +
    geom_point(size = 1.5) +
    facet_wrap(~ paste("r =", r), nrow = 1) +
    scale_colour_manual(values = c("#d6604d", "#4393c3", "#1b7837"), name = "d") +
    labs(title = "a  Coarse profile: logLik relative to the maximum (truncated at -40)",
         x = "lambda", y = "delta logLik") +
    theme(panel.grid.major.y = element_line(colour = "grey90"))
  p2$dll <- pmax(p2$logLik - max(p2$logLik, ll_est), -15)
  threshold <- -stats::qchisq(0.95, 2) / 2
  pts <- data.frame(lambda = c(truth2$lambda, est[["lambda"]]), d = c(truth2$d, est[["d"]]),
                    what = c(sprintf("true (r = %d)", truth2$r), "estimated"))
  p_c <- ggplot(p2, aes(lambda, d)) + geom_raster(aes(fill = dll), interpolate = TRUE) +
    geom_contour(aes(z = dll), breaks = threshold, colour = "black", linewidth = 0.5) +
    geom_line(data = curves, aes(linetype = what), colour = "white", linewidth = 0.8) +
    geom_point(data = ridge, colour = "white", size = 1.2) +
    geom_point(data = pts, aes(shape = what), size = 3, colour = "black", fill = "white") +
    scale_shape_manual(values = c(21, 4), name = NULL) +
    scale_linetype_manual(values = c("solid", "22"), name = NULL,
                          guide = guide_legend(override.aes = list(colour = "black"))) +
    scale_fill_viridis_c(name = "delta\nlogLik") + coord_cartesian(expand = FALSE) +
    labs(title = sprintf("b  Fine profile at r = %d; black contour: 95 %% confidence region", r_hat),
         subtitle = "white curves: constant exponent through the truth; white points: ridge of the profile",
         x = "lambda", y = "d") + theme(legend.position = "bottom", legend.box = "vertical")
  kernel_shape <- function(r, lambda, dd, what) {
    o <- expand.grid(dx = -r:r, dy = -r:r); o$delta <- sqrt(o$dx^2 + o$dy^2)
    o <- o[o$delta > 0 & o$delta <= r, ]
    # weight relative to the nearest neighbour: n_i^-d is a scale common to the
    # inner plots, absorbed by the IGE variance
    data.frame(delta = sort(unique(o$delta)), w = sort(unique(o$delta))^(-lambda), what = what)
  }
  ker <- rbind(kernel_shape(truth2$r, truth2$lambda, truth2$d, "true"),
               kernel_shape(r_hat, est[["lambda"]], est[["d"]], "estimated"))
  p_k <- ggplot(ker, aes(delta, w, colour = what)) + geom_line() + geom_point(size = 2) +
    scale_colour_manual(values = c(true = "black", estimated = "#1b7837"), name = NULL) +
    labs(title = "c  Kernel shape: weight relative to the nearest neighbour",
         x = "distance (grid steps)", y = "relative weight") + theme(panel.grid.major.y = element_line(colour = "grey90"))
  nu <- do.call(rbind, lapply(names(Zs), function(m)
    data.frame(truth = ige2, blup = final_fits[[m]]$blups$dge_ige[gen_names, 2], neighbourhood = m)))
  nu$neighbourhood <- factor(nu$neighbourhood, levels = names(Zs))
  p_n <- ggplot(nu, aes(truth, blup)) + geom_point(size = 0.8, alpha = 0.6, colour = "#d6604d") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.5) +
    facet_wrap(~ neighbourhood, nrow = 1, scales = "free_y") +
    labs(title = "d  IGE: BLUP against true by fitted neighbourhood", x = "true IGE", y = "BLUP")
  ib <- rbind(data.frame(neighbourhood = summary_tab$neighbourhood, criterion = "r(IGE)", value = summary_tab$acc_ige),
              data.frame(neighbourhood = summary_tab$neighbourhood, criterion = "r(total value)", value = summary_tab$acc_tbv),
              data.frame(neighbourhood = summary_tab$neighbourhood, criterion = "IGE x k (true 0.4)", value = summary_tab$ige_k))
  ib$neighbourhood <- factor(ib$neighbourhood, levels = names(Zs))
  ib$criterion <- factor(ib$criterion, levels = unique(ib$criterion))
  p_b <- ggplot(ib, aes(neighbourhood, value, fill = neighbourhood)) + geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.3, size = 2.8) +
    facet_wrap(~ criterion, nrow = 1) +
    scale_fill_manual(values = c("#c2a5cf", "#1b7837", "#999999"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    scale_x_discrete(labels = c("naive", "estimated", "oracle")) +
    labs(title = "e  Model C: accuracy and IGE variance", x = NULL, y = NULL)
  fig4 <- p_g / (p_c | p_k) / (p_n | p_b) +
    plot_annotation(title = sprintf(
      "Unknown neighbourhood: true r = %d, lambda = %g, d = %g; estimated r = %d, lambda = %.2f, d = %.2f",
      truth2$r, truth2$lambda, truth2$d, r_hat, est[["lambda"]], est[["d"]]))
  ggsave(file.path(out_dir, "fig4_estimated_neighbourhood.png"), fig4, width = 14, height = 13, dpi = 150,
         bg = "white")
}

cat("\nFigures and tables written to", normalizePath(out_dir), "\n")
