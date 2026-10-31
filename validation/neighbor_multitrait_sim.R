# ==============================================================================
# neighbor_multitrait_sim.R -- calibration and power of the multi-trait
# Neighbor GWAS (rx_neighbor_gwas) by simulation
# ------------------------------------------------------------------------------
# Plants on a grid in two groups, markers coded -1/+1. Two traits:
#   y = mu_t + X1 b1_t + X2 b2_t + e
# with polygenic self effects on all markers (X1 = doses) and neighbour
# effects on all markers (X2 = Sato's neighbour covariates), correlated
# between traits (Sigma_1, Sigma_2), plus correlated residuals. One causal
# SNP per replicate carries a neighbour (similarity) effect on trait 1 only.
#
# Per replicate: rx_neighbor_gwas with struct "us", tests
#   dir, sim|dir, sim[y1]|dir, sim[y2]|dir, sim@common|dir, sim@het|dir.
# Calibration: the null SNPs (all but the causal one, and outside +/- 5
# markers around it): proportion of p < 0.05 and genomic control lambda per
# test. Power: rank and p-value of the causal SNP per test.
#
#   Rscript validation/neighbor_multitrait_sim.R [n_rep] [out.csv]
# ==============================================================================
suppressPackageStartupMessages({ library(remlax); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 1) as.integer(args[1]) else 20L
out <- if (length(args) >= 2) args[2] else "neighbor_multitrait_sim.csv"
backend <- Sys.getenv("PD_BACKEND", "cpu")

one_rep <- function(r) {
  set.seed(1000 + r)
  nr <- 20; nc <- 15; n <- nr * nc; q <- 400
  coord <- as.matrix(expand.grid(x = seq_len(nc), y = seq_len(nr)))
  grp <- ifelse(coord[, 2] <= nr / 2, 1L, 2L)
  M <- matrix(sample(c(-1, 1), n * q, TRUE), n, q,
              dimnames = list(paste0("p", 1:n), paste0("m", 1:q)))
  W <- rx_neighbor_weights(coord, 1.5, grp)          # queen neighbourhood
  X2 <- M * as.matrix(W %*% M)
  S1 <- matrix(c(1, 0.5, 0.5, 1), 2); S2 <- matrix(c(1, -0.4, -0.4, 1), 2)
  b1 <- matrix(rnorm(2 * q), q) %*% chol(S1) * sqrt(0.3 / q)
  b2 <- matrix(rnorm(2 * q), q) %*% chol(S2) * sqrt(0.6 / q)
  E <- matrix(rnorm(2 * n), n) %*% chol(matrix(c(1, 0.3, 0.3, 1), 2))
  causal <- sample(20:(q - 20), 1)
  Y <- M %*% b1 + X2 %*% b2 + E
  Y[, 1] <- Y[, 1] + 0.35 * X2[, causal] / sd(X2[, causal])
  d <- data.frame(plant = rownames(M), y1 = Y[, 1], y2 = Y[, 2], grp = factor(grp))
  tests <- c("dir", "sim|dir", "sim[y1]|dir", "sim[y2]|dir", "sim@common|dir", "sim@het|dir")
  t0 <- proc.time()[3]
  nb <- rx_neighbor_gwas(d, c("y1", "y2"), "plant", M, W, fixed = ~ grp, tests = tests,
                         fit_args = list(backend = backend), verbose = FALSE)
  secs <- proc.time()[3] - t0
  null <- setdiff(seq_len(q), (causal - 5):(causal + 5))
  do.call(rbind, lapply(tests, function(tt) {
    p <- nb$scan[[paste0("p_", tt)]]; df <- nb$scan[[paste0("df_", tt)]][1]
    chi <- stats::qchisq(p[null], df, lower.tail = FALSE)
    data.frame(rep = r, test = tt, df = df, type1_05 = mean(p[null] < 0.05),
               lambda = stats::median(chi) / stats::qchisq(0.5, df),
               p_causal = p[causal], rank_causal = rank(p)[causal],
               cor_self = nb$cor_self[1, 2], cor_nei = nb$cor_neighbour[1, 2],
               seconds = secs, check.names = FALSE)
  }))
}

res <- do.call(rbind, lapply(seq_len(n_rep), function(r) {
  x <- one_rep(r)
  cat(sprintf("rep %d: %s\n", r, paste(sprintf("%s lambda %.2f", x$test, x$lambda), collapse = " | ")))
  x
}))
utils::write.csv(res, out, row.names = FALSE)
agg <- aggregate(cbind(type1_05, lambda, p_causal_lt_1e3 = p_causal < 1e-3) ~ test, res, mean)
print(agg)
cat(sprintf("simulated cor self 0.5 / estimated %.2f; cor neighbour -0.4 / estimated %.2f\n",
            mean(res$cor_self), mean(res$cor_nei)))
