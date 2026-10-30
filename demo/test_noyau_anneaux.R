# ==============================================================================
# Free neighbourhood kernel: one ring per distance, IGE correlated with DGE
#
# Same design and same simulations as test_identifiabilite_d.R (r = 3, d = 0.5,
# power kernel lambda = 1.5 or exponential kernel lambda = 0.8, continuous
# field or heterogeneous sub-blocks, r(DGE, IGE) = -0.5). The model no longer
# knows the kernel shape. The neighbours are sorted into 6 distance rings
# (1, sqrt(2), 2, sqrt(5), sqrt(8), 3), with one incidence Z_k per ring, divided
# by n_i^d.
#
# A single term covers the DGE and the 6 rings: an fa structure on 7 "traits"
# (DGE, ring 1, ..., ring 6), with kinship K. With one factor f per genotype:
#   DGE       = lambda_D f + e_D      (specific variance psi_D free)
#   ring k    = l_k f                 (specific variance fixed at zero)
# The l_k give the kernel up to a scale, and lambda_D carries the DGE-IGE
# covariance. Three fits per point:
#   free   fa1, all specific variances free (control)
#   fa1    fa1, specific variances of the rings fixed at zero (the model)
#   fa2    fa2, same: the kernel shape may vary between genotypes (control,
#          the truth has rank 1)
# The dilution d is profiled on a grid.
#
#   Rscript test_noyau_anneaux.R task <i> <results_folder>      (one task)
#   Rscript test_noyau_anneaux.R bayes <i> <results_folder>     (posterior,
#     case x seed i among 16, at d = 0.5, from the constrained fa1)
#   Rscript test_noyau_anneaux.R variances <results_folder>     (variances from the posteriors)
#   Rscript test_noyau_anneaux.R list                           (number of tasks)
#   Rscript test_noyau_anneaux.R summary <results_folder>       (summary and figure)
# ==============================================================================

a0 <- commandArgs(trailingOnly = TRUE)
mode <- if (length(a0)) a0[1] else "list"
`%||%` <- function(a, b) if (is.null(a)) b else a

truthA <- list(r = 3L, d = 0.5, lambda = c(power = 1.5, exponential = 0.8),
               dge = 1.0, ige = 0.4, r_dge_ige = -0.5, block = 0.3, residual = 0.8)
kernels <- c(power = "power", exponential = "exponential")
shape_w <- function(delta, kernel, lambda) if (kernel == "power") delta^(-lambda) else exp(-lambda * (delta - 1))
d2_rings <- c(1, 2, 4, 5, 8, 9)                   # delta^2 of the 6 rings at r = 3
grid_d <- seq(0, 1.25, 0.25)
tasks <- expand.grid(d = grid_d, seed = 101:104, kernel = names(kernels),
                     sub_blocks = c("continuous", "heterogeneous"), stringsAsFactors = FALSE)
tasks <- tasks[, c("sub_blocks", "kernel", "seed", "d")]

tasks_b <- unique(tasks[, c("sub_blocks", "kernel", "seed")]); tasks_b$d <- truthA$d
if (mode == "list") { cat(nrow(tasks), "REML tasks,", nrow(tasks_b), "Bayesian tasks\n"); quit(save = "no") }

# Results written by an earlier, French version of this script: task files
# named tache_*.csv, French columns and values. Mapped to English on reading.
fr_cols <- c(sous_blocs = "sub_blocks", noyau = "kernel", graine = "seed", logLik_libre = "logLik_free",
             conv_libre = "conv_free", psi_anneaux_libre_max = "psi_rings_free_max", secondes = "seconds")
fr_vals <- c(continu = "continuous", heterogene = "heterogeneous", puissance = "power",
             exponentiel = "exponential")
to_english <- function(df) {
  names(df) <- ifelse(names(df) %in% names(fr_cols), fr_cols[names(df)], names(df))
  for (v in intersect(c("sub_blocks", "kernel"), names(df)))
    df[[v]] <- ifelse(df[[v]] %in% names(fr_vals), fr_vals[df[[v]]], df[[v]])
  df
}
task_file <- function(dir, i) {
  f <- file.path(dir, sprintf("task_%03d.csv", i))
  if (file.exists(f)) f else file.path(dir, sprintf("tache_%03d.csv", i))
}

# ------------------------------------------------------------------------------
# Task mode
# ------------------------------------------------------------------------------
if (mode %in% c("task", "bayes", "variances")) {
  i <- as.integer(a0[2]); dir <- a0[3]
  if (mode == "variances") { i <- 1L; dir <- a0[2] }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  suppressPackageStartupMessages(library(remlax))
  # Design, kinship and trend of the demo (without its plots)
  l <- readLines("demo_remlax_spatial_cpu.R")
  pre <- l[1:(grep("^# 3. Fits", l) - 2)]
  pre <- sub("^args <- commandArgs.*", "args <- character(0)", pre)
  pre <- pre[!grepl("library\\((ggplot2|patchwork)\\)", pre)]
  out_dir <- tempdir()
  eval(parse(text = pre))

  bands <- rep(seq_len(10), times = rep(1:5, 2))
  sub_block <- function(d2, type) switch(type,
    continuous = factor(rep(1L, nrow(d2))),
    heterogeneous = factor(paste((d2$col - 1L) %/% 5L, bands[d2$row])))

  rings <- function(d2, dd) {
    n <- nrow(d2)
    q2 <- outer(d2$row, d2$row, "-")^2 + outer(d2$col, d2$col, "-")^2
    nbr <- outer(as.integer(d2$sb), as.integer(d2$sb), "==") & q2 > 0 & q2 <= truthA$r^2
    n_i <- rowSums(nbr); div <- ifelse(n_i > 0, n_i^dd, 1)
    g <- as.integer(d2$gid)
    lapply(stats::setNames(d2_rings, sprintf("a%d", seq_along(d2_rings))), function(a) {
      ij <- which(nbr & q2 == a, arr.ind = TRUE)
      methods::as(Matrix::sparseMatrix(i = ij[, 1], j = g[ij[, 2]], x = 1 / div[ij[, 1]],
                                       dims = c(n, n_gen), dimnames = list(NULL, gen_names)), "CsparseMatrix")
    })
  }

  simulate <- function(seed, sub_blocks, kernel) {
    set.seed(seed)
    d2 <- d[sort(sample(nrow(d), round(0.9 * nrow(d)))), c("col", "row", "rep", "block", "plot", "gid")]
    rownames(d2) <- NULL
    d2$block <- droplevels(d2$block)
    d2$sb <- sub_block(d2, sub_blocks)
    nb <- rx_neighbourhood(coord = cbind(d2$row, d2$col), group = rep("A", nrow(d2)), block = d2$sb,
                           level = d2$gid, id = d2$plot, rank = truthA$r,
                           reach = truthA$lambda[[kernel]], dilution = truthA$d, window = "euclidean",
                           kernel = kernels[[kernel]], normalise = FALSE)
    Z <- nb$level[["A<-A"]]
    P <- Matrix::sparseMatrix(i = match(colnames(Z), gen_names), j = seq_len(ncol(Z)), x = 1,
                              dims = c(n_gen, ncol(Z)))
    Zv <- Z %*% Matrix::t(P)
    # Check: the sum of the rings weighted by the true shape gives back Zv
    Za <- rings(d2, truthA$d)
    Zs <- Reduce(`+`, Map(function(Zk, a) Zk * shape_w(sqrt(a), kernel, truthA$lambda[[kernel]]), Za, d2_rings))
    stopifnot(max(abs(as.matrix(Zs - Zv))) < 1e-10)
    Zd <- Matrix::sparseMatrix(i = seq_len(nrow(d2)), j = as.integer(d2$gid), x = 1,
                               dims = c(nrow(d2), n_gen), dimnames = list(NULL, gen_names))
    k <- rx_exposure(Zd, Zv, K = K)$k
    s2i <- truthA$ige / k
    Sig <- matrix(c(truthA$dge, rep(truthA$r_dge_ige * sqrt(truthA$dge * s2i), 2), s2i), 2)
    U <- t(chol(K)) %*% matrix(rnorm(2 * n_gen), n_gen, 2) %*% chol(Sig)
    b <- setNames(rnorm(nlevels(d2$block), 0, sqrt(truthA$block)), levels(d2$block))
    d2$y <- 6 + 0.2 * (d2$rep == "R2") - 0.1 * (d2$rep == "R3") + true_trend(d2$col, d2$row) +
      U[as.integer(d2$gid), 1] + as.numeric(Zv %*% U[, 2]) + b[as.character(d2$block)] +
      rnorm(nrow(d2), 0, sqrt(truthA$residual))
    d2
  }

  # theta: the fa term comes first (trapezoidal loadings by row, then log psi
  # of the 7 traits), then block, the 5 variances of spl2d, the residual.
  # The check below verifies this layout on the returned Sigma.
  n_load <- function(t, r) sum(pmin(seq_len(t), r))
  fit_model <- function(d2, Za, fa_rank, fix_rings, theta0 = NULL, bayes = NULL) {
    for (k in names(Za)) assign(paste0("Z_", k), Za[[k]])
    # DGE incidence on the 150 genotypes: a genotype without a sown plot keeps
    # its column, as in the ring incidences
    Z_dge <- Matrix::sparseMatrix(i = seq_len(nrow(d2)), j = as.integer(d2$gid), x = 1,
                                  dims = c(nrow(d2), n_gen), dimnames = list(NULL, gen_names))
    ige <- paste0("str(~ mm(Z_dge, name = 'dge') + ", paste(sprintf("mm(Z_%s, name = '%s')", names(Za), names(Za)), collapse = " + "),
                  sprintf(", struct = 'fa', rank = %d, K = K, name = 'ige')", fa_rank))
    random <- stats::as.formula(paste("~", ige, "+ block + spl2d(col, row, nseg = c(7, 10))"),
                                env = environment())
    nl <- n_load(7L, fa_rank)
    fixed <- if (fix_rings) nl + 2:7 else NULL
    t0 <- Sys.time()
    f <- rx_reml(y ~ rep, random = random, data = d2, backend = "cpu", verbose = FALSE,
                 hessian = FALSE, blups = FALSE, pev = FALSE, theta_init = theta0, fixed_theta = fixed,
                 bayes = bayes)
    th <- f$theta
    stopifnot(length(th) == nl + 7L + 7L)
    L <- matrix(0, 7, fa_rank); kk <- 0
    for (a in 1:7) for (b in seq_len(min(a, fa_rank))) { kk <- kk + 1; L[a, b] <- th[kk] }
    psi2 <- exp(2 * th[nl + 1:7])
    S <- as.matrix(f$sigmas$ige)
    stopifnot(max(abs(S - (L %*% t(L) + diag(psi2)))) < 1e-6 * max(abs(S)))
    list(f = f, L = L, psi2 = psi2, theta = th, nl = nl,
         seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
         conv = isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement)))
  }

  # Variances on the phenotypic scale, from the Bayesian draws and REML at
  # d = 0.5: DGE = lambda_D^2 + psi_D; IGE x k = mean of diag(M K M'),
  # M = sum_k l_k Z_k; block and residual = exp(2 theta).
  if (mode == "variances") {
    fs <- list.files(dir, "^bayes_[0-9]+\\.rds$", full.names = TRUE)
    out <- list()
    for (f in fs) {
      ib <- as.integer(sub("^bayes_([0-9]+)\\.rds$", "\\1", basename(f)))
      t <- tasks_b[ib, ]
      d2 <- simulate(t$seed, t$sub_blocks, t$kernel)
      Za <- rings(d2, t$d)
      ige_k <- function(l) { M <- Reduce(`+`, Map(`*`, Za, l)); mean(Matrix::rowSums((M %*% K) * M)) }
      th <- readRDS(f)$theta
      sel <- unique(round(seq(1, nrow(th), length.out = 300)))
      nl <- 7L; p <- ncol(th)
      dr <- data.frame(dge = th[sel, 1]^2 + exp(2 * th[sel, nl + 1]),
                       ige_k = vapply(sel, function(s) ige_k(th[s, 2:7]), 0),
                       block = exp(2 * th[sel, nl + 8]), residual = exp(2 * th[sel, p]))
      q <- function(x) stats::quantile(x, c(0.05, 0.5, 0.95))
      # REML at the same point (fa1 task at d = 0.5)
      it <- which(tasks$sub_blocks == t$sub_blocks & tasks$kernel == t$kernel &
                  tasks$seed == t$seed & tasks$d == t$d)
      fr <- task_file(dir, it)
      reml <- if (file.exists(fr)) { r <- utils::read.csv(fr)
        c(dge = r$lambda_D^2 + r$psi_D, ige_k = ige_k(unlist(r[sprintf("l%d", 1:6)]))) } else c(dge = NA, ige_k = NA)
      for (v in names(dr)) out[[length(out) + 1]] <- data.frame(
        t[, 1:3], component = v, q05 = q(dr[[v]])[1], median = q(dr[[v]])[2], q95 = q(dr[[v]])[3],
        reml = if (v %in% names(reml)) reml[[v]] else NA_real_,
        true = c(dge = truthA$dge, ige_k = truthA$ige, block = truthA$block, residual = truthA$residual)[[v]])
    }
    res_v <- do.call(rbind, out); rownames(res_v) <- NULL
    res_v$covers <- res_v$q05 <= res_v$true & res_v$true <= res_v$q95
    print(res_v, digits = 3)
    cat("\nCoverage of the 90 % intervals by component:\n")
    print(tapply(res_v$covers, res_v$component, function(x) sprintf("%d/%d", sum(x), length(x))))
    utils::write.csv(res_v, file.path(dir, "bayes_variances.csv"), row.names = FALSE)
    quit(save = "no")
  }

  if (mode == "bayes") {
    stopifnot(rx_bayes_available())
    t <- tasks_b[i, ]
    cat(sprintf("bayes %d: %s, %s, seed %d, d = %g\n", i, t$sub_blocks, t$kernel, t$seed, t$d))
    d2 <- simulate(t$seed, t$sub_blocks, t$kernel)
    Za <- rings(d2, t$d)
    free <- fit_model(d2, Za, 1L, FALSE)
    th1 <- free$theta; th1[free$nl + 2:7] <- -12
    fa1 <- fit_model(d2, Za, 1L, TRUE, th1)
    bb <- fit_model(d2, Za, 1L, TRUE, fa1$theta,
                    bayes = rx_bayes(start = "map", chains = 8L, superchains = 4L, warmup = 400L,
                                     draws = 400L, seed = i, effects = 0L))
    th <- bb$f$bayes$theta                       # draws x p
    L <- th[, 1:7]                               # DGE then rings 1 to 6 (rank 1)
    sg <- sign(L[, 2]); sg[sg == 0] <- 1; L <- L * sg
    psiD <- exp(2 * th[, bb$nl + 1])
    w <- L[, 2:7] / L[, 2]
    q <- function(x, p) unname(stats::quantile(x, p))
    res_b <- data.frame(t, ring = 1:6, delta = sqrt(d2_rings),
                        w_med = apply(w, 2, q, 0.5), w_q05 = apply(w, 2, q, 0.05), w_q95 = apply(w, 2, q, 0.95),
                        l_med = apply(L[, 2:7], 2, q, 0.5), l_q05 = apply(L[, 2:7], 2, q, 0.05),
                        l_q95 = apply(L[, 2:7], 2, q, 0.95), p_positive = colMeans(L[, 2:7] > 0),
                        w_true = shape_w(sqrt(d2_rings), t$kernel, truthA$lambda[[t$kernel]]),
                        r_med = q(L[, 1] / sqrt(L[, 1]^2 + psiD), 0.5),
                        r_q05 = q(L[, 1] / sqrt(L[, 1]^2 + psiD), 0.05),
                        r_q95 = q(L[, 1] / sqrt(L[, 1]^2 + psiD), 0.95),
                        divergences = bb$f$bayes$diagnostics$divergences %||% NA,
                        rhat_max = bb$f$bayes$diagnostics$nested_rhat_max %||% NA,
                        ess_min = bb$f$bayes$diagnostics$ess_min %||% NA,
                        l_reml = fa1$L[2:7, 1] * sign(fa1$L[2, 1]), seconds = bb$seconds)
    print(res_b, digits = 3)
    str(bb$f$bayes$diagnostics)
    utils::write.csv(res_b, file.path(dir, sprintf("bayes_%02d.csv", i)), row.names = FALSE)
    saveRDS(list(theta = th, diagnostics = bb$f$bayes$diagnostics, summary = bb$f$bayes$summary),
            file.path(dir, sprintf("bayes_%02d.rds", i)))
    quit(save = "no")
  }

  t <- tasks[i, ]
  cat(sprintf("task %d: %s, %s, seed %d, d = %g\n", i, t$sub_blocks, t$kernel, t$seed, t$d))
  d2 <- simulate(t$seed, t$sub_blocks, t$kernel)
  Za <- rings(d2, t$d)
  free <- fit_model(d2, Za, 1L, FALSE)
  # constrained fa1: start at the free optimum, log psi of the rings at the floor
  th1 <- free$theta; th1[free$nl + 2:7] <- -12
  fa1 <- fit_model(d2, Za, 1L, TRUE, th1)
  # constrained fa2: first factor = that of fa1, small second factor
  L2 <- cbind(fa1$L[, 1], c(0, rep(0.05, 6)))
  th2 <- c(unlist(lapply(1:7, function(a) L2[a, seq_len(min(a, 2L))])),
           fa1$theta[fa1$nl + 1:7], fa1$theta[(fa1$nl + 8):length(fa1$theta)])
  fa2 <- fit_model(d2, Za, 2L, TRUE, th2)
  # sign: loading of ring 1 positive
  l <- fa1$L[, 1]; if (l[2] < 0) l <- -l
  row_out <- data.frame(t, logLik_free = free$f$logLik, logLik_fa1 = fa1$f$logLik, logLik_fa2 = fa2$f$logLik,
                        conv_free = free$conv, conv_fa1 = fa1$conv, conv_fa2 = fa2$conv,
                        lambda_D = l[1], psi_D = fa1$psi2[1], r_dge_ige = l[1] / sqrt(l[1]^2 + fa1$psi2[1]),
                        t(stats::setNames(l[2:7], sprintf("l%d", 1:6))),
                        psi_rings_free_max = max(free$psi2[2:7]),
                        seconds = free$seconds + fa1$seconds + fa2$seconds)
  print(row_out)
  utils::write.csv(row_out, file.path(dir, sprintf("task_%03d.csv", i)), row.names = FALSE)
  quit(save = "no")
}

# ------------------------------------------------------------------------------
# Summary mode
# ------------------------------------------------------------------------------
if (mode == "summary") {
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  dir <- a0[2]
  fs <- list.files(dir, "^task_.*\\.csv$", full.names = TRUE)
  if (!length(fs)) fs <- list.files(dir, "^tache_.*\\.csv$", full.names = TRUE)
  res <- do.call(rbind, lapply(fs, function(f) to_english(utils::read.csv(f))))
  cat(sprintf("%d tasks out of %d\n", nrow(res), nrow(tasks)))
  res$case <- paste(c(continuous = "continuous field", heterogeneous = "heterogeneous sub-blocks")[res$sub_blocks],
                    res$kernel, sep = ", ")
  summary_rings <- do.call(rbind, lapply(split(res, list(res$case, res$seed), drop = TRUE), function(q) {
    q <- q[order(q$d), ]; j <- which.max(q$logLik_fa1); jj <- min(max(j, 2L), nrow(q) - 1L)
    cf <- stats::coef(stats::lm(q$logLik_fa1[(jj - 1):(jj + 1)] ~ poly(q$d[(jj - 1):(jj + 1)], 2, raw = TRUE)))
    d_hat <- if (cf[3] < 0) min(max(-cf[2] / (2 * cf[3]), 0), 1.25) else q$d[j]
    fine <- seq(0, 1.25, length.out = 501)
    yf <- stats::spline(q$d, q$logLik_fa1 - max(q$logLik_fa1), xout = fine)$y
    inside <- fine[yf >= -stats::qchisq(0.95, 1) / 2]
    data.frame(case = q$case[1], seed = q$seed[1], d_hat = d_hat, ci_low = min(inside),
               ci_high = max(inside), drop_at_true = q$logLik_fa1[q$d == truthA$d] - max(q$logLik_fa1),
               r_dge_ige = q$r_dge_ige[j], lrt_rank2 = 2 * (q$logLik_fa2[j] - q$logLik_fa1[j]),
               lrt_psi_rings = 2 * (q$logLik_free[j] - q$logLik_fa1[j]),
               n_not_conv = sum(!q$conv_fa1))
  }))
  summary_rings <- summary_rings[order(summary_rings$case, summary_rings$seed), ]
  cat("\nFree kernel by rings (true d = 0.5, true r(DGE, IGE) = -0.5, true rank 1)\n")
  print(summary_rings, digits = 3, row.names = FALSE)
  cat("\nMeans by case\n")
  print(stats::aggregate(cbind(d_hat, ci_width = ci_high - ci_low, r_dge_ige, lrt_rank2) ~ case,
                         summary_rings, mean), digits = 3)
  utils::write.csv(summary_rings, file.path(dir, "summary_rings.csv"), row.names = FALSE)

  fo <- do.call(rbind, lapply(split(res, list(res$case, res$seed), drop = TRUE), function(q) {
    q <- q[which.max(q$logLik_fa1), ]; l <- unlist(q[sprintf("l%d", 1:6)])
    data.frame(case = q$case, kernel = q$kernel, seed = q$seed, delta = sqrt(d2_rings), w = l / l[1])
  }))
  avg <- stats::aggregate(w ~ case + delta, fo, mean)
  vf <- do.call(rbind, lapply(names(kernels), function(nk) data.frame(
    kernel = nk, delta = seq(1, 3, 0.02), w = shape_w(seq(1, 3, 0.02), nk, truthA$lambda[[nk]]))))
  vf <- merge(vf, unique(fo[, c("case", "kernel")]))
  theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10)))
  p_f <- ggplot(fo, aes(delta, w)) + geom_hline(yintercept = 0, colour = "grey70") +
    geom_line(data = vf, colour = "black", linewidth = 0.7) +
    geom_line(aes(group = seed), colour = "#1b7837", alpha = 0.35) +
    geom_point(colour = "#1b7837", size = 1.2, alpha = 0.6) +
    geom_line(data = avg, colour = "#1b7837", linewidth = 1, linetype = 2) +
    facet_wrap(~ case, nrow = 1) +
    labs(title = "a  Free kernel at the best d: loadings relative to ring 1",
         subtitle = "thin line: one seed; dashes: mean over seeds; black: true shape",
         x = "distance (grid steps)", y = "relative weight") +
    theme(panel.grid.major.y = element_line(colour = "grey90"))
  p_d <- ggplot(summary_rings, aes(factor(seed), d_hat)) +
    geom_hline(yintercept = truthA$d, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high), colour = "#1b7837", size = 0.3) +
    facet_wrap(~ case, nrow = 1) +
    labs(title = "b  d estimated with the free kernel, 95 % profile interval (dashes: true d = 0.5)",
         x = "seed", y = "d") + theme(panel.grid.major.y = element_line(colour = "grey90"))
  fig <- p_f / p_d + plot_annotation(title = "Free neighbourhood kernel by rings, IGE correlated with DGE (fa, rings without specific variance)")
  ggsave(file.path(dir, "fig_ring_kernel.png"), fig, width = 15, height = 8.5, dpi = 150, bg = "white")
  cat("\nFigure:", normalizePath(file.path(dir, "fig_ring_kernel.png")), "\n")
}
