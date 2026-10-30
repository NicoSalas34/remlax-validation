# ==============================================================================
# Identifiability of the dilution d: disjoint sub-blocks and exponential kernel
#
# Reuses the design of demo_remlax_spatial_cpu.R (same genotypes, same kinship,
# same spatial trend, 10 % empty plots) and crosses two factors:
#   - continuous field, or disjoint sub-blocks of 5 x 5 plots (a plot never
#     sees a plot of another sub-block, so the number of neighbours n_i varies
#     much more);
#   - power kernel w = delta^-lambda (true lambda 1.5), or exponential kernel
#     w = exp(-lambda (delta - 1)) (true lambda 0.8, same weight at delta = 3).
# The radius is fixed at its true value r = 3 and d is 0.5. For each case and
# each seed, the REML likelihood of model B (independent DGE + IGE, block,
# spl2d) is profiled on a (lambda, d) grid. The profile in d (maximum over
# lambda at fixed d) gives the estimated d, its 95 % interval (drop of 1.92)
# and a standard error from the curvature.
#
#   Rscript test_identifiabilite_d.R [n_seeds] [cores] [output_folder]
# ==============================================================================

suppressPackageStartupMessages({
  library(remlax)
  library(ggplot2)
  library(patchwork)
})
a0 <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(a0) >= 1) as.integer(a0[1]) else 4L
cores <- if (length(a0) >= 2) as.integer(a0[2]) else 12L
out_dir <- if (length(a0) >= 3) a0[3] else "output_identifiability_d"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Design, kinship and trend of the demo (sections 1 and 2, without fitting)
l <- readLines("demo_remlax_spatial_cpu.R")
pre <- l[1:(grep("^# 3. Fits", l) - 2)]
pre <- sub("^args <- commandArgs.*", "args <- character(0)", pre)
out_dir_own <- out_dir
eval(parse(text = pre))
out_dir <- out_dir_own

truthD <- list(r = 3L, d = 0.5, lambda = c(power = 1.5, exponential = 0.8),
               dge = 1.0, ige = 0.4, r_dge_ige = -0.5, block = 0.3, residual = 0.8)
grids <- list(power = seq(0.5, 3, 0.25), exponential = seq(0, 2, 0.2))
grid_d <- seq(0, 1.25, 0.25)
kernels <- c(power = "power", exponential = "exponential")

# Disjoint sub-blocks: "continuous" (none), "5x5" (squares of 5 x 5 plots),
# "heterogeneous" (bands of 1, 2, 3, 4 then 5 rows over 5 columns, twice).
bands <- rep(seq_len(10), times = rep(1:5, 2))            # 30 rows
sub_block <- function(d2, type) switch(type,
  continuous = factor(rep(1L, nrow(d2))),
  "5x5" = factor(paste((d2$col - 1L) %/% 5L, (d2$row - 1L) %/% 5L)),
  heterogeneous = factor(paste((d2$col - 1L) %/% 5L, bands[d2$row])))

# Simulated data set for (seed, sub_blocks, kernel): empty plots, genetic
# values, blocks, residuals. Deterministic, hence rebuilt in each process.
simulate <- function(seed, sub_blocks, kernel) {
  set.seed(seed)
  d2 <- d[sort(sample(nrow(d), round(0.9 * nrow(d)))), c("col", "row", "rep", "block", "plot", "gid")]
  rownames(d2) <- NULL
  d2$block <- droplevels(d2$block)
  d2$sb <- sub_block(d2, sub_blocks)
  Zv <- neighbourhood_Z(d2, truthD$r, truthD$lambda[[kernel]], truthD$d, kernel)
  Zd <- Matrix::sparseMatrix(i = seq_len(nrow(d2)), j = as.integer(d2$gid), x = 1,
                             dims = c(nrow(d2), n_gen), dimnames = list(NULL, gen_names))
  k <- rx_exposure(Zd, Zv, K = K)$k
  s2i <- truthD$ige / k
  Sig <- matrix(c(truthD$dge, rep(truthD$r_dge_ige * sqrt(truthD$dge * s2i), 2), s2i), 2)
  U <- t(chol(K)) %*% matrix(rnorm(2 * n_gen), n_gen, 2) %*% chol(Sig)
  b <- setNames(rnorm(nlevels(d2$block), 0, sqrt(truthD$block)), levels(d2$block))
  d2$y <- 6 + 0.2 * (d2$rep == "R2") - 0.1 * (d2$rep == "R3") + true_trend(d2$col, d2$row) +
    U[as.integer(d2$gid), 1] + as.numeric(Zv %*% U[, 2]) + b[as.character(d2$block)] +
    rnorm(nrow(d2), 0, sqrt(truthD$residual))
  d2
}

neighbourhood_Z <- function(d2, r, lambda, dd, kernel) {
  nb <- rx_neighbourhood(coord = cbind(d2$row, d2$col), group = rep("A", nrow(d2)), block = d2$sb,
                         level = d2$gid, id = d2$plot, rank = as.integer(r),
                         reach = as.numeric(lambda), dilution = as.numeric(dd), window = "euclidean",
                         kernel = kernels[[kernel]], normalise = FALSE)
  Z <- nb$level[["A<-A"]]
  P <- Matrix::sparseMatrix(i = match(colnames(Z), gen_names), j = seq_len(ncol(Z)), x = 1,
                            dims = c(n_gen, ncol(Z)))
  Z <- Z %*% Matrix::t(P)
  dimnames(Z) <- list(NULL, gen_names)
  methods::as(Z, "CsparseMatrix")
}

fit_model <- function(d2, Z) {
  random <- stats::as.formula("~ vm(gid, K = K) + mm(Z, K = K, name = 'ige') + block + spl2d(col, row, nseg = c(7, 10))",
                              env = environment())
  rx_reml(y ~ rep, random = random, data = d2, backend = "cpu", verbose = FALSE,
          hessian = FALSE, blups = FALSE, pev = FALSE)
}

# Tasks: case x seed x lambda x d. The data sets are cached per process.
cases <- expand.grid(sub_blocks = c("continuous", "5x5", "heterogeneous"), kernel = names(kernels),
                     stringsAsFactors = FALSE)
tasks <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) do.call(rbind, lapply(seq_len(n_seeds), function(g)
  expand.grid(sub_blocks = cases$sub_blocks[i], kernel = cases$kernel[i], seed = 100L + g,
              lambda = grids[[cases$kernel[i]]], d = grid_d, stringsAsFactors = FALSE)))))
cache <- new.env()
one_task <- function(i) {
  t <- tasks[i, ]
  key <- paste(t$seed, t$sub_blocks, t$kernel)
  if (is.null(cache[[key]])) cache[[key]] <- simulate(t$seed, t$sub_blocks, t$kernel)
  d2 <- cache[[key]]
  f <- tryCatch(fit_model(d2, neighbourhood_Z(d2, truthD$r, t$lambda, t$d, t$kernel)), error = function(e) conditionMessage(e))
  if (is.character(f)) return(cbind(t, logLik = NA_real_, conv = NA, error = substr(f, 1, 150)))
  cbind(t, logLik = f$logLik, conv = isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement)),
        error = "")
}
# Dry run: one task per case, then stop.
if (n_seeds == 0L) {
  tasks <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) data.frame(
    sub_blocks = cases$sub_blocks[i], kernel = cases$kernel[i], seed = 101L,
    lambda = truthD$lambda[[cases$kernel[i]]], d = truthD$d, stringsAsFactors = FALSE)))
  print(do.call(rbind, lapply(seq_len(nrow(tasks)), one_task)))
  for (ty in cases$sub_blocks[!duplicated(cases$sub_blocks)]) {
    d2 <- simulate(101L, ty, "power")
    nb <- rx_neighbourhood(cbind(d2$row, d2$col), group = rep("A", nrow(d2)), block = d2$sb,
                           id = d2$plot, rank = 3L, window = "euclidean", output = "unit")
    n_i <- nb$n_neighbours[[1]]
    cat(sprintf("%-13s n_i from %d to %d, mean %.1f, sd(log n_i) = %.3f, %d without neighbour\n", ty,
                min(n_i), max(n_i), mean(n_i), sd(log(n_i[n_i > 0])), sum(n_i == 0)))
  }
  quit(save = "no")
}
Sys.setenv(XLA_FLAGS = "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2")
options(remlax.threads = 2L)
cat(sprintf("%d fits over %d processes\n", nrow(tasks), cores))
t0 <- Sys.time()
# The tasks of one data set go to the same process (prepartition by contiguous blocks).
chunks <- split(seq_len(nrow(tasks)), cut(seq_len(nrow(tasks)), cores, labels = FALSE))
res <- do.call(rbind, parallel::mclapply(chunks, function(ii) do.call(rbind, lapply(ii, one_task)),
                                          mc.cores = cores, mc.preschedule = FALSE))
cat(sprintf("done in %.1f min; %d failure(s), %d not converged\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(is.na(res$logLik)),
            sum(!res$conv, na.rm = TRUE)))
if (any(is.na(res$logLik))) print(unique(res$error[is.na(res$logLik)]))
res$case <- paste(c(continuous = "continuous field", "5x5" = "5 x 5 sub-blocks",
                    heterogeneous = "heterogeneous sub-blocks")[res$sub_blocks], res$kernel, sep = ", ")
write.csv(res, file.path(out_dir, "profiles.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Profile in d: maximum over lambda at fixed d; estimate, interval, curvature
# ------------------------------------------------------------------------------
parabola <- function(x, y) {           # vertex and curvature from the 3 points around the max
  j <- which.max(y); j <- min(max(j, 2L), length(y) - 1L)
  cf <- stats::coef(stats::lm(y[(j - 1):(j + 1)] ~ poly(x[(j - 1):(j + 1)], 2, raw = TRUE)))
  list(vertex = if (cf[3] < 0) -cf[2] / (2 * cf[3]) else x[which.max(y)], curvature = 2 * cf[3])
}
interval <- function(x, y, drop = stats::qchisq(0.95, 1) / 2) {
  y <- y - max(y); fine <- seq(min(x), max(x), length.out = 501)
  yf <- stats::spline(x, y, xout = fine)$y
  inside <- fine[yf >= -drop]
  c(low = min(inside), high = max(inside))
}
summary_d <- do.call(rbind, lapply(split(res, list(res$case, res$seed), drop = TRUE), function(q) {
  pd <- stats::aggregate(logLik ~ d, q, max); pl <- stats::aggregate(logLik ~ lambda, q, max)
  sd_ <- parabola(pd$d, pd$logLik); sl <- parabola(pl$lambda, pl$logLik)
  ci <- interval(pd$d, pd$logLik)
  data.frame(case = q$case[1], kernel = q$kernel[1], seed = q$seed[1],
             d_hat = min(max(sd_$vertex, 0), 1.25), lambda_hat = sl$vertex,
             ci_low = ci[["low"]], ci_high = ci[["high"]], ci_width = ci[["high"]] - ci[["low"]],
             sd_d = if (sd_$curvature < 0) 1 / sqrt(-sd_$curvature) else NA_real_,
             drop_at_true = pd$logLik[pd$d == truthD$d] - max(pd$logLik))
}))
summary_d <- summary_d[order(summary_d$case, summary_d$seed), ]
write.csv(summary_d, file.path(out_dir, "summary_d.csv"), row.names = FALSE)
cat("\nEstimate of d by case and seed (true d = 0.5)\n")
print(summary_d, digits = 3, row.names = FALSE)
cat("\nMeans by case\n")
print(stats::aggregate(cbind(d_hat, ci_width, sd_d, drop_at_true) ~ case, summary_d, mean), digits = 3)

# ------------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 10) +
          theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10)))
q1 <- res[res$seed == min(res$seed), ]
q1 <- do.call(rbind, lapply(split(q1, q1$case), function(q) { q$dll <- pmax(q$logLik - max(q$logLik), -12); q }))
vv <- unique(data.frame(case = q1$case, lambda = truthD$lambda[q1$kernel], d = truthD$d))
p_h <- ggplot(q1, aes(lambda, d)) + geom_raster(aes(fill = dll), interpolate = TRUE) +
  geom_contour(aes(z = dll), breaks = -stats::qchisq(0.95, 2) / 2, colour = "black", linewidth = 0.5) +
  geom_point(data = vv, shape = 4, size = 3) +
  facet_wrap(~ case, nrow = 1, scales = "free_x") +
  scale_fill_viridis_c(name = "delta\nlogLik") + coord_cartesian(expand = FALSE) +
  labs(title = sprintf("a  Profile (lambda, d), seed %d; cross = truth, contour = 95 %% region", min(res$seed)),
       x = "lambda", y = "d")
p_d <- ggplot(summary_d, aes(factor(seed), d_hat, colour = kernel)) +
  geom_hline(yintercept = truthD$d, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high), size = 0.3) +
  facet_wrap(~ case, nrow = 1) +
  scale_colour_manual(values = c(power = "#2166ac", exponential = "#d6604d"), guide = "none") +
  labs(title = "b  Estimated d and 95 % profile interval, by seed (dashes: true d = 0.5)",
       x = "seed", y = "d") + theme(panel.grid.major.y = element_line(colour = "grey90"))
fig <- p_h / p_d + plot_annotation(title = "Identifiability of the dilution: disjoint sub-blocks and kernel shape")
ggsave(file.path(out_dir, "fig_identifiability_d.png"), fig, width = 15, height = 8.5, dpi = 150, bg = "white")
cat("\nOutputs in", normalizePath(out_dir), "\n")
