# ==============================================================================
# Identifiability of the dilution d by number of replicates
#
# Enlarged design of the demo: 150 genotypes, each replicate fills 10 rows x 15
# columns (one genotype per replicate), with 3, 6 or 12 replicates.
# Heterogeneous disjoint sub-blocks (bands of 1 to 5 rows over 5 columns, the
# best case of test_identifiabilite_d.R), 10 % empty plots, power kernel r = 3,
# lambda = 1.5, d = 0.5, IGE 0.4 of phenotypic variance, r(DGE, IGE) = -0.5.
# The spatial trend of the demo is stretched over the length of the field.
# REML profile of model B (independent DGE + IGE, block, spl2d) over
# lambda in {1; 1.5; 2; 2.5} and d in {0; 0.25; ...; 1.25}.
#
#   Rscript test_repetitions_d.R list
#   Rscript test_repetitions_d.R task <i> <folder>   (6 fits: one lambda, the whole grid in d)
#   Rscript test_repetitions_d.R summary <folder>
# ==============================================================================

a0 <- commandArgs(trailingOnly = TRUE)
mode <- if (length(a0)) a0[1] else "list"
truthR <- list(r = 3L, lambda = 1.5, d = 0.5, dge = 1.0, ige = 0.4, r_dge_ige = -0.5,
               block = 0.3, spatial = 1.2, residual = 0.8)
grid_l <- c(1, 1.5, 2, 2.5)
grid_d <- seq(0, 1.25, 0.25)
tasks <- expand.grid(lambda = grid_l, seed = 101:104, n_rep = c(3L, 6L, 12L))
tasks <- tasks[, c("n_rep", "seed", "lambda")]
# Task files of the earlier French version (tache_*.csv) carry French column
# names; they are mapped to English on reading.
fr_cols <- c(graine = "seed", secondes = "seconds", erreur = "error")
to_english <- function(x) {
  names(x) <- ifelse(names(x) %in% names(fr_cols), fr_cols[names(x)], names(x))
  x
}
if (mode == "list") { cat(nrow(tasks), "tasks\n"); quit(save = "no") }

if (mode == "task") {
  suppressPackageStartupMessages(library(remlax))
  i <- as.integer(a0[2]); dir <- a0[3]
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  t <- tasks[i, ]
  n_gen <- 150L; n_col <- 15L; n_row <- 10L * t$n_rep
  gen_names <- sprintf("G%03d", seq_len(n_gen))

  # Kinship: markers of the demo (seed 5), two subpopulations
  set.seed(5L)
  p0 <- runif(1500L, 0.05, 0.5); pop <- rep(1:2, length.out = n_gen)
  P <- sapply(p0, function(p) rbeta(2, p * 9, (1 - p) * 9))
  M <- matrix(rbinom(n_gen * 1500L, 2L, P[pop, ]), n_gen, 1500L); rownames(M) <- gen_names
  K <- rx_grm(M, blend = 0.02, coding = "count")

  # Design and data
  set.seed(t$seed)
  d <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))
  d$rep <- factor(sprintf("R%02d", (d$row - 1L) %/% 10L + 1L))
  d$block <- factor(sprintf("%s_r%03d_b%d", d$rep, d$row, (d$col - 1L) %/% 5L + 1L))
  d$plot <- sprintf("P%04d", seq_len(nrow(d)))
  d$gid <- factor(NA, levels = gen_names)
  for (r in levels(d$rep)) d$gid[d$rep == r] <- sample(gen_names)
  heights <- rep(1:5, length.out = ceiling(n_row / 3))
  bands <- rep(seq_along(heights), times = heights)[seq_len(n_row)]
  d$sb <- factor(paste((d$col - 1L) %/% 5L, bands[d$row]))
  d <- d[sort(sample(nrow(d), round(0.9 * nrow(d)))), ]; rownames(d) <- NULL
  d$block <- droplevels(d$block)
  # Trend of the demo, stretched over the length of the field
  tb <- function(x, y) { y <- 1 + (y - 1) * 29 / (n_row - 1)
    0.9 * (y - 15.5) / 14.5 - 1.6 * exp(-((x - 11)^2 / 10 + (y - 9)^2 / 30)) + 0.35 * sin(2 * pi * x / 10) }
  br <- tb(d$col, d$row); d$spatial <- (br - mean(br)) / sd(br) * sqrt(truthR$spatial)

  neighbourhood_Z <- function(lambda, dd) {
    nb <- rx_neighbourhood(coord = cbind(d$row, d$col), group = rep("A", nrow(d)), block = d$sb,
                           level = d$gid, id = d$plot, rank = truthR$r, reach = as.numeric(lambda),
                           dilution = as.numeric(dd), window = "euclidean", kernel = "power")
    Z <- nb$level[["A<-A"]]
    Pm <- Matrix::sparseMatrix(i = match(colnames(Z), gen_names), j = seq_len(ncol(Z)), x = 1,
                               dims = c(n_gen, ncol(Z)))
    Z <- Z %*% Matrix::t(Pm); dimnames(Z) <- list(NULL, gen_names)
    methods::as(Z, "CsparseMatrix")
  }
  Zd <- Matrix::sparseMatrix(i = seq_len(nrow(d)), j = as.integer(d$gid), x = 1,
                             dims = c(nrow(d), n_gen), dimnames = list(NULL, gen_names))
  Zv <- neighbourhood_Z(truthR$lambda, truthR$d)
  s2i <- truthR$ige / rx_exposure(Zd, Zv, K = K)$k
  Sig <- matrix(c(truthR$dge, rep(truthR$r_dge_ige * sqrt(truthR$dge * s2i), 2), s2i), 2)
  U <- t(chol(K)) %*% matrix(rnorm(2 * n_gen), n_gen, 2) %*% chol(Sig)
  b <- setNames(rnorm(nlevels(d$block), 0, sqrt(truthR$block)), levels(d$block))
  d$y <- 6 + d$spatial + U[as.integer(d$gid), 1] + as.numeric(Zv %*% U[, 2]) +
    b[as.character(d$block)] + rnorm(nrow(d), 0, sqrt(truthR$residual))
  nseg_r <- round(10 * t$n_rep / 3)

  outputs <- lapply(grid_d, function(dd) {
    Z <- neighbourhood_Z(t$lambda, dd)
    random <- stats::as.formula(sprintf(
      "~ vm(gid, K = K) + mm(Z, K = K, name = 'ige') + block + spl2d(col, row, nseg = c(7, %d))", nseg_r),
      env = environment())
    t0 <- Sys.time()
    f <- tryCatch(rx_reml(y ~ rep, random = random, data = d, backend = "cpu", verbose = FALSE,
                          hessian = FALSE, blups = FALSE, pev = FALSE), error = function(e) conditionMessage(e))
    if (is.character(f)) return(data.frame(t, d = dd, n = nrow(d), logLik = NA, conv = NA, seconds = NA,
                                           error = substr(f, 1, 150)))
    data.frame(t, d = dd, n = nrow(d), logLik = f$logLik,
               conv = isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement)),
               seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")), error = "")
  })
  res <- do.call(rbind, outputs)
  print(res)
  utils::write.csv(res, file.path(dir, sprintf("task_%03d.csv", i)), row.names = FALSE)
  quit(save = "no")
}

if (mode == "summary") {
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  dir <- a0[2]
  fs <- list.files(dir, "^task_.*\\.csv$", full.names = TRUE)
  if (!length(fs)) fs <- list.files(dir, "^tache_.*\\.csv$", full.names = TRUE)
  res <- do.call(rbind, lapply(fs, function(f) to_english(utils::read.csv(f))))
  cat(sprintf("%d fits, %d failure(s), %d not converged\n", nrow(res), sum(is.na(res$logLik)),
              sum(!res$conv, na.rm = TRUE)))
  summary_rep <- do.call(rbind, lapply(split(res, list(res$n_rep, res$seed), drop = TRUE), function(q) {
    pd <- stats::aggregate(logLik ~ d, q, max); pd <- pd[order(pd$d), ]
    j <- which.max(pd$logLik); jj <- min(max(j, 2L), nrow(pd) - 1L)
    cf <- stats::coef(stats::lm(pd$logLik[(jj - 1):(jj + 1)] ~ poly(pd$d[(jj - 1):(jj + 1)], 2, raw = TRUE)))
    d_hat <- if (cf[3] < 0) min(max(-cf[2] / (2 * cf[3]), 0), 1.25) else pd$d[j]
    fine <- seq(0, 1.25, length.out = 501)
    yf <- stats::spline(pd$d, pd$logLik - max(pd$logLik), xout = fine)$y
    inside <- fine[yf >= -stats::qchisq(0.95, 1) / 2]
    data.frame(n_rep = q$n_rep[1], n = q$n[1], seed = q$seed[1], d_hat = d_hat,
               lambda_hat = q$lambda[which.max(q$logLik)], ci_low = min(inside), ci_high = max(inside),
               sd_d = if (cf[3] < 0) 1 / sqrt(-2 * cf[3]) else NA_real_,
               drop_at_true = pd$logLik[pd$d == truthR$d] - max(pd$logLik), drop_at_0 = pd$logLik[1] - max(pd$logLik))
  }))
  summary_rep <- summary_rep[order(summary_rep$n_rep, summary_rep$seed), ]
  print(summary_rep, digits = 3, row.names = FALSE)
  cat("\nMeans by number of replicates\n")
  print(stats::aggregate(cbind(n, d_hat, ci_width = ci_high - ci_low, sd_d, drop_at_0) ~ n_rep,
                         summary_rep, mean), digits = 3)
  utils::write.csv(summary_rep, file.path(dir, "summary_replicates.csv"), row.names = FALSE)
  pr <- do.call(rbind, lapply(split(res, list(res$n_rep, res$seed), drop = TRUE), function(q) {
    pd <- stats::aggregate(logLik ~ d, q, max); pd$dll <- pd$logLik - max(pd$logLik)
    cbind(pd, n_rep = q$n_rep[1], seed = q$seed[1]) }))
  theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10)))
  p1 <- ggplot(pr, aes(d, dll, group = seed)) + geom_hline(yintercept = -stats::qchisq(0.95, 1) / 2,
                                                             linetype = 3, colour = "grey50") +
    geom_vline(xintercept = truthR$d, linetype = 2, colour = "grey50") +
    geom_line(colour = "#1b7837") + geom_point(colour = "#1b7837", size = 1) +
    facet_wrap(~ sprintf("%d replicates", n_rep), nrow = 1) + coord_cartesian(ylim = c(-8, 0)) +
    labs(title = "a  Likelihood profile in d (maximum over lambda), one curve per seed",
         x = "d", y = "delta logLik") + theme(panel.grid.major.y = element_line(colour = "grey90"))
  p2 <- ggplot(summary_rep, aes(factor(seed), d_hat)) +
    geom_hline(yintercept = truthR$d, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high), colour = "#1b7837", size = 0.3) +
    facet_wrap(~ sprintf("%d replicates (n = %d)", n_rep, n), nrow = 1) +
    labs(title = "b  Estimated d and 95 % interval (dashes: true d = 0.5)", x = "seed", y = "d") +
    theme(panel.grid.major.y = element_line(colour = "grey90"))
  ggsave(file.path(dir, "fig_replicates_d.png"), p1 / p2 +
           plot_annotation(title = "Dilution d and number of replicates (heterogeneous sub-blocks, power kernel)"),
         width = 12, height = 7.5, dpi = 150, bg = "white")
  cat("Figure:", normalizePath(file.path(dir, "fig_replicates_d.png")), "\n")
}
