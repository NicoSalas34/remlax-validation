# ==============================================================================
# Simulated cases of the absolute benchmark: time against n, against p, and the
# dim(C) >> n regime.
#
# Case name: sim_q<q>_t<t>_r<r>_s<struct>_x<x>_o<o>
#   q       phenotyped genotypes
#   t       environments (traits). t = 1 gives a GBLUP
#   r       replicates per genotype and per observed environment
#   struct  G x E structure: id (t = 1), diag, fa1, fa2, fa3, us
#   x       NON-phenotyped genotypes per phenotyped genotype (prediction set).
#           The kinship holds q (1 + x) genotypes, and so does C
#   o       percentage of observed (genotype, environment) pairs (sparse
#           trials). 100 = complete design
#
# Model shared by the three tools:
#   y = mu_env + g + e,  g ~ N(0, Sigma_g (x) K),  e ~ N(0, sigma2_env) per env.
# Kinship: VanRaden on 3000 markers drawn in 10 subpopulations,
# K = 0.98 G / mean(diag G) + 0.02 I (positive definite even if q > markers).
#
# `prepare_sim()` writes the data sets once in <dir>/sim/, so that the three
# tools and the two clusters read the same numbers.
# ==============================================================================

parse_sim_name <- function(case) {
  v <- function(key) sub(paste0(".*_", key, "([^_]+).*"), "\\1", case)
  list(q = as.integer(v("q")), t = as.integer(v("t")), r = as.integer(v("r")),
       s = v("s"), x = as.numeric(v("x")), o = as.numeric(v("o")))
}

# kinship of q_tot genotypes, seed fixed by q_tot: all cases with the same
# q_tot share the same matrix
simulate_K <- function(q_tot, m = 3000L, n_pop = 10L) {
  set.seed(1000L + q_tot)
  p0 <- stats::runif(m, 0.05, 0.5)
  pop <- sample.int(n_pop, q_tot, replace = TRUE)
  # frequencies per subpopulation: Balding-Nichols type drift, Fst 0.1
  Fst <- 0.1
  P <- sapply(p0, function(p) stats::rbeta(n_pop, p * (1 - Fst) / Fst, (1 - p) * (1 - Fst) / Fst))
  P <- pmin(pmax(P, 0.01), 0.99)
  M <- matrix(stats::rbinom(q_tot * m, 2L, P[pop, ]), q_tot, m)
  f <- colMeans(M) / 2
  W <- sweep(M, 2L, 2 * f)
  G <- tcrossprod(W) / (2 * sum(f * (1 - f)))
  K <- 0.98 * G / mean(diag(G)) + 0.02 * diag(q_tot)
  ids <- sprintf("G%05d", seq_len(q_tot))
  dimnames(K) <- list(ids, ids)
  K
}

sigma_g <- function(t) {
  if (t == 1L) return(matrix(1))
  set.seed(77L + t)
  L <- cbind(stats::runif(t, 0.6, 0.95), stats::rnorm(t, 0, 0.35))
  S <- tcrossprod(L)
  S <- S + diag(pmax(1 - diag(S), 0.05))
  D <- diag(1 / sqrt(diag(S)))
  D %*% S %*% D                                  # genetic variances 1
}

simulate_data <- function(case, K) {
  a <- parse_sim_name(case)
  q <- a$q; t <- a$t
  q_tot <- nrow(K)
  set.seed(sum(utf8ToInt(case)))
  ids <- rownames(K)
  obs <- ids[seq_len(q)]                        # the first q are phenotyped
  # genetic values: chol(K)' Z chol(Sigma)
  Lk <- chol(K[obs, obs])
  Sg <- sigma_g(t)
  # sigma_g() resets the seed when t > 1. Without this line, the multi-trait
  # _gN data sets were copies of g1. Names without a suffix keep the old draw,
  # so that the g1 measurements remain valid. The multi-trait _g2 and _g3 files
  # written before this fix are stored in sim/copies_g1.
  if (grepl("_g[0-9]+$", case)) set.seed(sum(utf8ToInt(case)) + 1L)
  U <-crossprod(Lk, matrix(stats::rnorm(q * t), q, t)) %*% chol(Sg)
  # observed pairs: each genotype in max(1, round(o t / 100)) environments
  k_env <- max(1L, round(a$o * t / 100))
  cells <- do.call(rbind, lapply(seq_len(q), function(i)
    cbind(i, if (k_env >= t) seq_len(t) else sort(sample.int(t, k_env)))))
  cells <- cells[rep(seq_len(nrow(cells)), each = a$r), , drop = FALSE]
  sig_e <- seq(1.5, 2.5, length.out = t)
  env <- factor(sprintf("E%02d", cells[, 2]), levels = sprintf("E%02d", seq_len(t)))
  mu <- seq(5, 10, length.out = t)
  y <- mu[cells[, 2]] + U[cells] + stats::rnorm(nrow(cells), 0, sqrt(sig_e[cells[, 2]]))
  d <- data.frame(GY = y, Env = env, Gen = factor(obs[cells[, 1]], levels = ids))
  d <- d[order(d$Env, d$Gen), ]
  d$Unit <- factor(sprintf("U%07d", seq_len(nrow(d))))
  rownames(d) <- NULL
  d
}

K_file <- function(dir, q_tot) file.path(dir, "sim", sprintf("K_q%d.rds", q_tot))
data_file <- function(dir, case) file.path(dir, "sim", paste0(case, ".rds"))

prepare_sim <- function(case_list, dir) {
  dir.create(file.path(dir, "sim"), showWarnings = FALSE, recursive = TRUE)
  for (case in case_list) {
    a <- parse_sim_name(case)
    q_tot <- as.integer(round(a$q * (1 + a$x)))
    fk <- K_file(dir, q_tot)
    if (!file.exists(fk)) {
      t0 <- Sys.time()
      saveRDS(simulate_K(q_tot), fk, compress = FALSE)
      cat(sprintf("K q_tot = %d: %.1f s\n", q_tot, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    fd <- data_file(dir, case)
    if (!file.exists(fd)) {
      d <- simulate_data(case, readRDS(fk))
      saveRDS(d, fd)
      cat(sprintf("%s: n = %d, q_tot = %d\n", case, nrow(d), q_tot))
    }
  }
}

load_sim <- function(case, dir) {
  a <- parse_sim_name(case)
  q_tot <- as.integer(round(a$q * (1 + a$x)))
  d <- readRDS(data_file(dir, case))
  K <- readRDS(K_file(dir, q_tot))
  rank <- if (startsWith(a$s, "fa")) as.integer(sub("fa", "", a$s)) else 0L
  list(case = case, d = d, K = K, type = "sim", t = a$t, ge = a$s, rank = rank)
}
