# ==============================================================================
# Neighbourhood without a radius to estimate: every neighbour within a wide
# window, exponential weights, dilution
#
# Fitted model: w_ij = exp(-lambda (delta_ij - 1)) / n_i^d for every neighbour
# at distance <= R = 6 (in the same sub-block), n_i = number of these
# neighbours. Only lambda and d remain to estimate. The reach goes through
# lambda.
#
# Data: exactly the data sets of test_identifiabilite_d.R (same seeds, same
# cases), simulated with a radius r = 3, d = 0.5, power kernel (lambda 1.5) or
# exponential kernel (lambda 0.8). The model without a radius is therefore
# misspecified on the cut at 3 and on the n_i of the dilution. The logLik values
# compare directly with those of the model at r = 3 (profiles.csv of the
# identifiability test).
#
#   Rscript test_sans_rayon.R list
#   Rscript test_sans_rayon.R task <i> <folder>   (one case, one seed, one lambda: 6 fits in d)
#   Rscript test_sans_rayon.R summary <folder> [profiles_r3.csv]
# The r = 3 profiles may also be the profils.csv of the earlier French version
# (demo/sortie_identifiabilite_d/profils.csv); French columns and values are
# mapped to English on reading, as are the French task files tache_*.csv.
# ==============================================================================

a0 <- commandArgs(trailingOnly = TRUE)
mode <- if (length(a0)) a0[1] else "list"
R_WINDOW <- 6L
truthD <- list(r = 3L, d = 0.5, lambda = c(power = 1.5, exponential = 0.8),
               dge = 1.0, ige = 0.4, r_dge_ige = -0.5, block = 0.3, residual = 0.8)
kernels <- c(power = "power", exponential = "exponential")
grid_l <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3)
grid_d <- seq(0, 1.25, 0.25)
tasks <- expand.grid(lambda = grid_l, seed = 101:104, kernel = names(kernels),
                     sub_blocks = c("continuous", "5x5", "heterogeneous"), stringsAsFactors = FALSE)
tasks <- tasks[, c("sub_blocks", "kernel", "seed", "lambda")]
if (mode == "list") { cat(nrow(tasks), "tasks\n"); quit(save = "no") }

# Old French columns and values, mapped to English on reading.
fr_cols <- c(sous_blocs = "sub_blocks", noyau = "kernel", graine = "seed", erreur = "error", cas = "case")
fr_vals <- c(continu = "continuous", heterogene = "heterogeneous", puissance = "power",
             exponentiel = "exponential")
to_english <- function(x) {
  names(x) <- ifelse(names(x) %in% names(fr_cols), fr_cols[names(x)], names(x))
  for (v in intersect(c("sub_blocks", "kernel"), names(x)))
    x[[v]] <- ifelse(x[[v]] %in% names(fr_vals), fr_vals[x[[v]]], x[[v]])
  x
}

if (mode == "task") {
  suppressPackageStartupMessages(library(remlax))
  i <- as.integer(a0[2]); dir <- a0[3]
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  l <- readLines("demo_remlax_spatial_cpu.R")
  pre <- l[1:(grep("^# 3. Fits", l) - 2)]
  pre <- sub("^args <- commandArgs.*", "args <- character(0)", pre)
  pre <- pre[!grepl("library\\((ggplot2|patchwork)\\)", pre)]
  out_dir <- tempdir()
  eval(parse(text = pre))

  bands <- rep(seq_len(10), times = rep(1:5, 2))
  sub_block <- function(d2, type) switch(type,
    continuous = factor(rep(1L, nrow(d2))),
    "5x5" = factor(paste((d2$col - 1L) %/% 5L, (d2$row - 1L) %/% 5L)),
    heterogeneous = factor(paste((d2$col - 1L) %/% 5L, bands[d2$row])))
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
  # Simulation identical to test_identifiabilite_d.R (same order of draws)
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

  t <- tasks[i, ]
  d2 <- simulate(t$seed, t$sub_blocks, t$kernel)
  res <- do.call(rbind, lapply(grid_d, function(dd) {
    Z <- neighbourhood_Z(d2, R_WINDOW, t$lambda, dd, "exponential")
    random <- stats::as.formula(
      "~ vm(gid, K = K) + mm(Z, K = K, name = 'ige') + block + spl2d(col, row, nseg = c(7, 10))",
      env = environment())
    f <- tryCatch(rx_reml(y ~ rep, random = random, data = d2, backend = "cpu", verbose = FALSE,
                          hessian = FALSE, blups = FALSE, pev = FALSE), error = function(e) conditionMessage(e))
    if (is.character(f)) return(data.frame(t, d = dd, logLik = NA, conv = NA, error = substr(f, 1, 150)))
    data.frame(t, d = dd, logLik = f$logLik,
               conv = isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement)), error = "")
  }))
  print(res)
  utils::write.csv(res, file.path(dir, sprintf("task_%03d.csv", i)), row.names = FALSE)
  quit(save = "no")
}

if (mode == "summary") {
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  dir <- a0[2]
  f_r3 <- if (length(a0) >= 3) a0[3] else NA
  fs <- list.files(dir, "^task_.*\\.csv$", full.names = TRUE)
  if (!length(fs)) fs <- list.files(dir, "^tache_.*\\.csv$", full.names = TRUE)
  res <- do.call(rbind, lapply(fs, function(f) to_english(utils::read.csv(f))))
  cat(sprintf("%d fits, %d failure(s), %d not converged\n", nrow(res), sum(is.na(res$logLik)),
              sum(!res$conv, na.rm = TRUE)))
  case_name <- function(sb, nk) paste(c(continuous = "continuous field", "5x5" = "5 x 5 sub-blocks",
                                        heterogeneous = "heterogeneous sub-blocks")[sb], nk, sep = ", ")
  res$case <- case_name(res$sub_blocks, res$kernel)
  r3 <- if (!is.na(f_r3) && file.exists(f_r3)) to_english(utils::read.csv(f_r3)) else NULL
  parabola <- function(x, y) {
    j <- which.max(y); jj <- min(max(j, 2L), length(y) - 1L)
    cf <- stats::coef(stats::lm(y[(jj - 1):(jj + 1)] ~ poly(x[(jj - 1):(jj + 1)], 2, raw = TRUE)))
    list(vertex = if (cf[3] < 0) min(max(-cf[2] / (2 * cf[3]), min(x)), max(x)) else x[j],
         se = if (cf[3] < 0) 1 / sqrt(-2 * cf[3]) else NA_real_)
  }
  summary_nr <- do.call(rbind, lapply(split(res, list(res$case, res$seed), drop = TRUE), function(q) {
    pd <- stats::aggregate(logLik ~ d, q, max); pd <- pd[order(pd$d), ]
    pl <- stats::aggregate(logLik ~ lambda, q, max); pl <- pl[order(pl$lambda), ]
    sd_ <- parabola(pd$d, pd$logLik); sl <- parabola(pl$lambda, pl$logLik)
    fine <- seq(0, 1.25, length.out = 501)
    yf <- stats::spline(pd$d, pd$logLik - max(pd$logLik), xout = fine)$y
    inside <- fine[yf >= -stats::qchisq(0.95, 1) / 2]
    fl <- seq(0, 3, length.out = 601)
    yl <- stats::spline(pl$lambda, pl$logLik - max(pl$logLik), xout = fl)$y
    dl <- fl[yl >= -stats::qchisq(0.95, 1) / 2]
    ll_r3 <- if (is.null(r3)) NA_real_ else {
      m <- r3[r3$sub_blocks == q$sub_blocks[1] & r3$kernel == q$kernel[1] & r3$seed == q$seed[1], ]
      if (nrow(m)) max(m$logLik, na.rm = TRUE) else NA_real_ }
    data.frame(case = q$case[1], kernel = q$kernel[1], seed = q$seed[1],
               d_hat = sd_$vertex, d_low = min(inside), d_high = max(inside), sd_d = sd_$se,
               lambda_hat = sl$vertex, l_low = min(dl), l_high = max(dl),
               drop_d_true = pd$logLik[pd$d == truthD$d] - max(pd$logLik),
               logLik_max = max(q$logLik), logLik_r3 = ll_r3, delta_vs_r3 = max(q$logLik) - ll_r3)
  }))
  summary_nr <- summary_nr[order(summary_nr$case, summary_nr$seed), ]
  print(summary_nr, digits = 3, row.names = FALSE)
  cat("\nMeans by case\n")
  print(stats::aggregate(cbind(d_hat, width_d = d_high - d_low, sd_d, lambda_hat,
                               width_l = l_high - l_low, delta_vs_r3) ~ case, summary_nr, mean), digits = 3)
  utils::write.csv(summary_nr, file.path(dir, "summary_no_radius.csv"), row.names = FALSE)

  theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10)))
  q1 <- res[res$seed == 101, ]
  q1 <- do.call(rbind, lapply(split(q1, q1$case), function(q) { q$dll <- pmax(q$logLik - max(q$logLik), -12); q }))
  p_h <- ggplot(q1, aes(lambda, d)) + geom_tile(aes(fill = dll)) +
    geom_contour(aes(z = dll), breaks = -stats::qchisq(0.95, 2) / 2, colour = "black", linewidth = 0.5) +
    geom_hline(yintercept = truthD$d, linetype = 2, colour = "white") +
    facet_wrap(~ case, nrow = 1) + scale_fill_viridis_c(name = "delta\nlogLik") + coord_cartesian(expand = FALSE) +
    labs(title = "a  Profile (lambda, d) of the model without a radius, seed 101; contour = 95 % region, dashes: true d",
         x = "lambda (exponential, window R = 6)", y = "d")
  pdl <- do.call(rbind, lapply(split(res, list(res$case, res$seed), drop = TRUE), function(q) {
    pd <- stats::aggregate(logLik ~ d, q, max); pd$dll <- pd$logLik - max(pd$logLik)
    cbind(pd, case = q$case[1], seed = q$seed[1]) }))
  p_d <- ggplot(pdl, aes(d, dll, group = seed)) +
    geom_hline(yintercept = -stats::qchisq(0.95, 1) / 2, linetype = 3, colour = "grey50") +
    geom_vline(xintercept = truthD$d, linetype = 2, colour = "grey50") +
    geom_line(colour = "#1b7837") + geom_point(colour = "#1b7837", size = 0.8) +
    facet_wrap(~ case, nrow = 1) + coord_cartesian(ylim = c(-6, 0)) +
    labs(title = "b  Profile in d (maximum over lambda), one curve per seed", x = "d", y = "delta logLik") +
    theme(panel.grid.major.y = element_line(colour = "grey90"))
  p_c <- ggplot(summary_nr, aes(factor(seed), delta_vs_r3)) + geom_hline(yintercept = 0, colour = "grey50") +
    geom_col(fill = "#2166ac", width = 0.6) + facet_wrap(~ case, nrow = 1) +
    labs(title = "c  Maximum logLik: model without a radius minus model at r = 3 (lambda, d profile of the previous test)",
         x = "seed", y = "delta logLik") + theme(panel.grid.major.y = element_line(colour = "grey90"))
  ggsave(file.path(dir, "fig_no_radius.png"), p_h / p_d / p_c +
           plot_annotation(title = "Neighbourhood without a radius: exponential kernel on a window R = 6, lambda and d estimated"),
         width = 16, height = 11, dpi = 150, bg = "white")
  cat("Figure:", normalizePath(file.path(dir, "fig_no_radius.png")), "\n")
}
