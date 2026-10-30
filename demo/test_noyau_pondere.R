# ==============================================================================
# Dilution by the kernel mass: test on one data set
#
# Same simulation as test_noyau_deux_especes.R (chapter 2 of the thesis): two
# species, 100 genotypes each, correlated DGE, intra IGE and inter IGE, trays of
# 20 columns x 16 rows, species alternating by column, 10 % missing plants.
#
# The dilution of rx_neighbourhood divides each row by n_i^d, where n_i is the
# NUMBER of neighbours in the window of radius r, so d depends on r. Here the
# group size is the MASS of the kernel, m_i = sum_j w_ij, and each row is
# divided by m_i^d (d = 0: sum; d = 1: weighted mean). With an exponential
# kernel, the mass converges and the window only serves the computation.
#
# Two truths on the same design (n_rep = 10, seed 1):
#   A : exponential kernel w = exp(-lambda (delta - 1)), lambda = 1 per grid
#       step, dilution by the mass d = 0.5, computation window R = 8;
#   B : the former truth, power delta^-1.5 truncated at r = 3, dilution by the
#       count d = 0.5.
# Fits for each truth:
#   - "weighted" : lambda x d grid, exponential, mass, R = 8;
#   - "window"   : lambda = 1, d = 0.5, mass, R in {4, 6, 12};
#   - "count"    : exponential lambda = 1, dilution by the COUNT, R = 8, d on a grid;
#   - "power"    : power r = 3, lambda = 1.5, d = 0.5, count (truth B);
#   - "ch3"      : power r = 3, lambda = 1, d = 0;
#   - "naive"    : Chebyshev rank 1, equal weights.
#
# One task = one truth and one seed (30 seeds x 2 truths = 60 tasks). The seed
# sets the design and the genetic values.
#
#   Rscript test_noyau_pondere.R task <i> <folder> [cores]
# ==============================================================================

a0 <- commandArgs(trailingOnly = TRUE)
i_task <- as.integer(a0[2]); dir <- a0[3]
cores <- if (length(a0) >= 4) as.integer(a0[4]) else 1L
nv <- c("A", "B")[(i_task - 1L) %/% 30L + 1L]
seed <- (i_task - 1L) %% 30L + 1L
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(remlax); library(Matrix) })

truth <- list(h2 = 0.3, dge = 1, ige_intra_k = 0.2, ige_inter_k = 0.3,
              r_dge_intra = -0.6, r_dge_inter = -0.6, r_intra_inter = 0.6, missing = 0.10)
n_rep <- 10L; n_gen <- 100L
gA <- sprintf("A%03d", seq_len(n_gen)); gB <- sprintf("B%03d", seq_len(n_gen))

# --- design (identical to test_noyau_deux_especes.R) ------------------------
set.seed(1000L * n_rep + seed)
n_tray <- ceiling(n_gen * n_rep / 160)
d <- do.call(rbind, lapply(seq_len(n_tray), function(b)
  data.frame(tray = sprintf("tray%02d", b), col = rep(1:20, 16), row = rep(1:16, each = 20))))
d$species <- ifelse(d$col %% 2 == 1, "A", "B")
d$gen <- NA_character_
for (sp in c("A", "B")) {
  w <- which(d$species == sp)
  d$gen[w] <- sample(rep_len(if (sp == "A") gA else gB, length(w)))
}
d <- d[sort(sample(nrow(d), round((1 - truth$missing) * nrow(d)))), ]
d <- d[order(d$species != "A"), ]
rownames(d) <- NULL
d$id <- sprintf("P%05d", seq_len(nrow(d)))
d$species <- factor(d$species); d$tray <- factor(d$tray)
d$sp_tray <- factor(paste(d$species, d$tray))
n <- nrow(d); lA <- d$species == "A"; lB <- !lA

incid <- function(nb, pair, lev) {
  Z <- nb$level[[pair]]
  Zt <- methods::as(Z, "TsparseMatrix")
  Matrix::sparseMatrix(i = match(rownames(Z), d$id)[Zt@i + 1L], j = match(colnames(Z), lev)[Zt@j + 1L],
                       x = Zt@x, dims = c(n, length(lev)), dimnames = list(NULL, lev))
}
# Neighbourhood. dilution = "count": n_i^-d of rx_neighbourhood; "mass": m_i^-d,
# m_i = sum of the weights of the row (rx_neighbourhood called without dilution).
neighbourhood_Z <- function(kernel, r, lambda, dd, dilution, window = "euclidean") {
  nb <- rx_neighbourhood(coord = cbind(d$row, d$col), group = as.character(d$species), block = d$tray,
                         level = d$gen, id = d$id, rank = r, reach = lambda,
                         dilution = if (dilution == "count") dd else 0, window = window,
                         kernel = kernel, output = "level")
  V <- list(AA = incid(nb, "A<-A", gA), BA = incid(nb, "B<-A", gA),
            BB = incid(nb, "B<-B", gB), AB = incid(nb, "A<-B", gB))
  if (dilution == "mass" && dd > 0) V <- lapply(V, function(Z) {
    m <- Matrix::rowSums(Z); m[m == 0] <- 1
    Zd <- methods::as(Matrix::Diagonal(x = m^(-dd)) %*% Z, "CsparseMatrix")
    dimnames(Zd) <- dimnames(Z); Zd
  })
  V
}
kexp <- function(Z, rows) mean(Matrix::rowSums(Z^2)[rows])
ZdA <- Matrix::sparseMatrix(i = which(lA), j = match(d$gen[lA], gA), x = 1, dims = c(n, n_gen),
                            dimnames = list(NULL, gA))
ZdB <- Matrix::sparseMatrix(i = which(lB), j = match(d$gen[lB], gB), x = 1, dims = c(n, n_gen),
                            dimnames = list(NULL, gB))

truths <- list(
  A = list(kernel = "exponential", r = 8, lambda = 1, d = 0.5, dilution = "mass"),
  B = list(kernel = "power", r = 3, lambda = 1.5, d = 0.5, dilution = "count"))

# Draw of the genetic values and of the phenotype under one truth, with a seed
# specific to each data set.
simulate <- function(vt) {
  Vv <- neighbourhood_Z(vt$kernel, vt$r, vt$lambda, vt$d, vt$dilution)
  kv <- c(AA = kexp(Vv$AA, lA), BA = kexp(Vv$BA, lB), BB = kexp(Vv$BB, lB), AB = kexp(Vv$AB, lA))
  set.seed(5000L + seed)
  R <- matrix(c(1, truth$r_dge_intra, truth$r_dge_inter, truth$r_dge_intra, 1, truth$r_intra_inter,
                truth$r_dge_inter, truth$r_intra_inter, 1), 3)
  draw <- function(ki, kx) {
    v <- c(truth$dge, truth$ige_intra_k / ki, truth$ige_inter_k / kx)
    matrix(rnorm(3 * n_gen), n_gen, 3) %*% chol(diag(sqrt(v)) %*% R %*% diag(sqrt(v)))
  }
  UA <- draw(kv[["AA"]], kv[["BA"]]); UB <- draw(kv[["BB"]], kv[["AB"]])
  g <- as.numeric(ZdA %*% UA[, 1] + Vv$AA %*% UA[, 2] + Vv$BA %*% UA[, 3] +
                  ZdB %*% UB[, 1] + Vv$BB %*% UB[, 2] + Vv$AB %*% UB[, 3])
  y <- ifelse(lA, 10, 5) + g + rnorm(n, 0, sqrt(truth$dge * (1 / truth$h2 - 1)))
  s <- c(A = mean(Matrix::rowSums(Vv$AA)[lA]), B = mean(Matrix::rowSums(Vv$BB)[lB]))
  list(y = y, k = kv, s = s,
       true_vals = list(A = cbind(UA, UA[, 1] + s[["A"]] * UA[, 2]), B = cbind(UB, UB[, 1] + s[["B"]] * UB[, 2])))
}

plan <- function() {
  g <- expand.grid(lambda = c(0.5, 0.75, 1, 1.5, 2), d = c(0, 0.25, 0.5, 0.75, 1))
  rbind(data.frame(grid = "weighted", kernel = "exponential", r = 8, g, dilution = "mass", window = "euclidean"),
        data.frame(grid = "window", kernel = "exponential", r = c(4, 6, 12), lambda = 1, d = 0.5,
                   dilution = "mass", window = "euclidean"),
        data.frame(grid = "count", kernel = "exponential", r = 8, lambda = 1, d = c(0, 0.25, 0.5, 0.75, 1),
                   dilution = "count", window = "euclidean"),
        data.frame(grid = "power", kernel = "power", r = 3, lambda = 1.5, d = 0.5, dilution = "count",
                   window = "euclidean"),
        data.frame(grid = "ch3", kernel = "power", r = 3, lambda = 1, d = 0, dilution = "count",
                   window = "euclidean"),
        data.frame(grid = "naive", kernel = "none", r = 1, lambda = 0, d = 0, dilution = "count",
                   window = "chebyshev"))
}

fit_model <- function(p, sim) {
  t0 <- Sys.time()
  dd <- d; dd$y <- sim$y
  V <- neighbourhood_Z(p$kernel, p$r, p$lambda, p$d, p$dilution, p$window)
  ZaD <- ZdA; ZaI <- V$AA; ZaX <- V$BA; ZbD <- ZdB; ZbI <- V$BB; ZbX <- V$AB
  random <- ~ str(~ mm(ZaD) + mm(ZaI) + mm(ZaX), struct = "us", name = "gA") +
              str(~ mm(ZbD) + mm(ZbI) + mm(ZbX), struct = "us", name = "gB")
  environment(random) <- environment()
  f <- tryCatch(rx_reml(y ~ sp_tray, random = random, residual = ~ dsum(~ units | species), data = dd,
                        backend = "cpu", verbose = FALSE, hessian = FALSE, pev = FALSE),
                error = function(e) conditionMessage(e))
  base <- data.frame(p, logLik = NA_real_, conv = NA, seconds = NA_real_, error = "")
  if (is.character(f)) { base$error <- substr(gsub("\n", " ", f), 1, 200); return(base) }
  base$logLik <- f$logLik
  base$conv <- isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement))
  v <- summary(f)$varcomp
  k_f <- list(A = c(1, kexp(V$AA, lA), kexp(V$BA, lB)), B = c(1, kexp(V$BB, lB), kexp(V$AB, lA)))
  s_f <- c(A = mean(Matrix::rowSums(V$AA)[lA]), B = mean(Matrix::rowSums(V$BB)[lB]))
  vk <- c(truth$dge, truth$ige_intra_k, truth$ige_inter_k)
  eff <- c("dge", "intra", "inter")
  for (sp in c("A", "B")) {
    nm <- paste0("g", sp)
    S <- matrix(NA_real_, 3, 3)
    for (a in 1:3) for (b in 1:a) {
      lab <- sprintf("%s[%d,%d]", nm, a, b)
      if (lab %in% rownames(v)) S[a, b] <- S[b, a] <- v[lab, "component"]
    }
    Bl <- f$blups[[nm]][if (sp == "A") gA else gB, , drop = FALSE]
    Bl <- cbind(Bl, Bl[, 1] + s_f[[sp]] * Bl[, 2])
    for (a in 1:3) base[[sprintf("err_%s_%s", eff[a], sp)]] <- log10(S[a, a] * k_f[[sp]][a] / vk[a])
    for (a in 1:4) base[[sprintf("acc_%s_%s", c(eff, "tbv")[a], sp)]] <- stats::cor(Bl[, a], sim$true_vals[[sp]][, a])
  }
  base$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  base
}

options(remlax.threads = 1L)
pl <- plan()
sim <- simulate(truths[[nv]])
cat(sprintf("task %d: truth %s, seed %d, n = %d, %d trays, %d fits on %d cores\n",
            i_task, nv, seed, n, n_tray, nrow(pl), cores))
cat(sprintf("k AA %.3f BA %.3f BB %.3f AB %.3f; s A %.3f B %.3f\n",
            sim$k[1], sim$k[2], sim$k[3], sim$k[4], sim$s[1], sim$s[2]))
t1 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(pl)), function(j) fit_model(pl[j, ], sim), mc.cores = cores,
                          mc.preschedule = FALSE)
ok <- vapply(res, is.data.frame, TRUE)
res[!ok] <- lapply(which(!ok), function(j) data.frame(pl[j, ], logLik = NA_real_, conv = NA,
                                                     seconds = NA_real_, error = "process lost"))
nm_all <- unique(unlist(lapply(res, names)))
res <- do.call(rbind, lapply(res, function(x) { x[setdiff(nm_all, names(x))] <- NA; x[nm_all] }))
res <- cbind(truth = nv, seed = seed, res)
utils::write.csv(res, file.path(dir, sprintf("truth_%s_seed_%02d.csv", nv, seed)), row.names = FALSE)
cat(sprintf("finished in %.1f min; %d failure(s), %d not converged\n",
            as.numeric(difftime(Sys.time(), t1, units = "mins")), sum(is.na(res$logLik)),
            sum(!res$conv, na.rm = TRUE)))
