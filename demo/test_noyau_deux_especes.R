# ==============================================================================
# Estimation of the neighbourhood kernel (lambda, d) and effect of a wrong kernel
# on the IGE, in the setting of chapter 2 of the thesis (Salas et al., Genetics)
#
# Two-species simulation A and B, as in the article: 100 genotypes per species,
# correlated DGE, intra IGE and inter IGE (r(DGE, intra IGE) = r(DGE, inter
# IGE) = -0.6, r(intra IGE, inter IGE) = 0.6), unrelated genotypes (M = I),
# H2_DGE = 0.3, no IEE. The design follows the geometry of the real design:
# trays of 20 columns x 16 rows, species alternating by column, grid step 1
# (5 cm), 10 % missing plants. Trays and missing plants make the number of
# neighbours vary. Without them d would not be identifiable (a field without
# edges gives every plant the same n_i).
#
# True kernel: Euclidean window of radius r = 3, weights delta^-lambda with
# lambda = 1.5, dilution n_i^-d with d = 0.5, n_i counted per emitting species.
# The IGE variances are set on the phenotypic scale: k x V(intra IGE) = 0.2 and
# k x V(inter IGE) = 0.3 (V(DGE) = 1), the order of magnitude of the article
# (ratio 0.05 to 0.1 per neighbour, 2 intra and 6 inter neighbours).
#
# For each data set (n_rep in {5, 10, 20}, 30 seeds), the model of the article
# (us 3 x 3 per species, residual per species, fixed tray) is fitted:
#   - "r3"    : lambda x d grid at r = 3 (the true radius);
#   - "R6"    : the same grid in a wide window R = 6;
#   - "r"     : r in {1, 2, 4, 5} at the true lambda and d;
#   - "naive" : Chebyshev rank 1 (8 neighbours), equal weights, no dilution,
#               the neighbourhood of the simulations of the article.
# Each fit gives its logLik (profile) and, against the simulated truth, the
# accuracy of the BLUPs (DGE, intra IGE, inter IGE, intra TBV) and the error
# log10(estimated V x k / true V x k) of each component.
#
#   Rscript test_noyau_deux_especes.R list
#   Rscript test_noyau_deux_especes.R task <i> <folder> [cores]
#   Rscript test_noyau_deux_especes.R summary <folder>
# QUICK=1 restricts a task to the truth, the naive kernel and one wrong radius.
# ==============================================================================

a0 <- commandArgs(trailingOnly = TRUE)
mode <- if (length(a0)) a0[1] else "list"
truth <- list(r = 3, lambda = 1.5, d = 0.5, h2 = 0.3, dge = 1, ige_intra_k = 0.2, ige_inter_k = 0.3,
              r_dge_intra = -0.6, r_dge_inter = -0.6, r_intra_inter = 0.6, missing = 0.10)
grid_l <- c(0.5, 1, 1.5, 2, 2.5, 3)
grid_d <- c(0, 0.25, 0.5, 0.75, 1)
r_wrong <- c(1, 2, 4, 5)
R_WIDE <- 6
tasks <- expand.grid(seed = 1:30, n_rep = c(5L, 10L, 20L))[, c("n_rep", "seed")]
if (mode == "list") { cat(nrow(tasks), "tasks\n"); quit(save = "no") }

# List of the fits of one task
plan_fits <- function() {
  g <- expand.grid(lambda = grid_l, d = grid_d)
  rbind(data.frame(grid = "r3", r = truth$r, g, window = "euclidean"),
        data.frame(grid = "R6", r = R_WIDE, g, window = "euclidean"),
        data.frame(grid = "r", r = r_wrong, lambda = truth$lambda, d = truth$d, window = "euclidean"),
        data.frame(grid = "naive", r = 1, lambda = 0, d = 0, window = "chebyshev"))
}

if (mode == "task") {
  suppressPackageStartupMessages({ library(remlax); library(Matrix) })
  i <- as.integer(a0[2]); dir <- a0[3]
  cores <- if (length(a0) >= 4) as.integer(a0[4]) else 1L
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  t <- tasks[i, ]
  n_gen <- 100L
  gA <- sprintf("A%03d", seq_len(n_gen)); gB <- sprintf("B%03d", seq_len(n_gen))

  # --- design: trays 20 x 16, odd columns A, even columns B -------------------
  set.seed(1000L * t$n_rep + t$seed)
  n_tray <- ceiling(n_gen * t$n_rep / 160)
  d <- do.call(rbind, lapply(seq_len(n_tray), function(b)
    data.frame(tray = sprintf("tray%02d", b), col = rep(1:20, 16), row = rep(1:16, each = 20))))
  d$species <- ifelse(d$col %% 2 == 1, "A", "B")
  d$gen <- NA_character_
  for (sp in c("A", "B")) {
    w <- which(d$species == sp)
    d$gen[w] <- sample(rep_len(if (sp == "A") gA else gB, length(w)))
  }
  d <- d[sort(sample(nrow(d), round((1 - truth$missing) * nrow(d)))), ]
  d <- d[order(d$species != "A"), ]                   # rows A then B
  rownames(d) <- NULL
  d$id <- sprintf("P%05d", seq_len(nrow(d)))
  d$species <- factor(d$species); d$tray <- factor(d$tray)
  d$sp_tray <- factor(paste(d$species, d$tray))
  n <- nrow(d); lA <- d$species == "A"; lB <- !lA

  # Incidence n x 100 of a receiving<-emitting pair, columns = emitting genotypes
  incid <- function(nb, pair, lev) {
    Z <- nb$level[[pair]]
    Zt <- methods::as(Z, "TsparseMatrix")
    Matrix::sparseMatrix(i = match(rownames(Z), d$id)[Zt@i + 1L], j = match(colnames(Z), lev)[Zt@j + 1L],
                         x = Zt@x, dims = c(n, length(lev)), dimnames = list(NULL, lev))
  }
  neighbourhood_Z <- function(r, lambda, dd, window = "euclidean") {
    nb <- rx_neighbourhood(coord = cbind(d$row, d$col), group = as.character(d$species), block = d$tray,
                           level = d$gen, id = d$id, rank = r, reach = unname(lambda),
                           dilution = unname(dd), window = window, kernel = "power", output = "level")
    list(AA = incid(nb, "A<-A", gA), BA = incid(nb, "B<-A", gA),   # emitted by A
         BB = incid(nb, "B<-B", gB), AB = incid(nb, "A<-B", gB))   # emitted by B
  }
  kexp <- function(Z, rows) mean(Matrix::rowSums(Z^2)[rows])  # k = mean diag(ZZ'), M = I
  ZdA <- Matrix::sparseMatrix(i = which(lA), j = match(d$gen[lA], gA), x = 1, dims = c(n, n_gen),
                              dimnames = list(NULL, gA))
  ZdB <- Matrix::sparseMatrix(i = which(lB), j = match(d$gen[lB], gB), x = 1, dims = c(n, n_gen),
                              dimnames = list(NULL, gB))

  # --- truth ------------------------------------------------------------------
  Vv <- neighbourhood_Z(truth$r, truth$lambda, truth$d)
  k_true <- c(AA = kexp(Vv$AA, lA), BA = kexp(Vv$BA, lB), BB = kexp(Vv$BB, lB), AB = kexp(Vv$AB, lA))
  draw <- function(k_intra, k_inter) {
    v <- c(truth$dge, truth$ige_intra_k / k_intra, truth$ige_inter_k / k_inter)
    R <- matrix(c(1, truth$r_dge_intra, truth$r_dge_inter, truth$r_dge_intra, 1, truth$r_intra_inter,
                  truth$r_dge_inter, truth$r_intra_inter, 1), 3)
    S <- diag(sqrt(v)) %*% R %*% diag(sqrt(v))
    list(U = matrix(rnorm(3 * n_gen), n_gen, 3) %*% chol(S), S = S)
  }
  tA <- draw(k_true[["AA"]], k_true[["BA"]]); tB <- draw(k_true[["BB"]], k_true[["AB"]])
  UA <- tA$U; UB <- tB$U; dimnames(UA) <- dimnames(UB) <- NULL
  v_dee <- truth$dge * (1 / truth$h2 - 1)
  g <- as.numeric(ZdA %*% UA[, 1] + Vv$AA %*% UA[, 2] + Vv$BA %*% UA[, 3] +
                  ZdB %*% UB[, 1] + Vv$BB %*% UB[, 2] + Vv$AB %*% UB[, 3])
  d$y <- ifelse(lA, 10, 5) + g + rnorm(n, 0, sqrt(v_dee))
  # true intra TBV: DGE + (mean sum of intra weights) x intra IGE
  s_true <- c(A = mean(Matrix::rowSums(Vv$AA)[lA]), B = mean(Matrix::rowSums(Vv$BB)[lB]))
  true_vals <- list(A = cbind(UA, UA[, 1] + s_true[["A"]] * UA[, 2]),
                    B = cbind(UB, UB[, 1] + s_true[["B"]] * UB[, 2]))
  vk <- list(A = c(truth$dge, truth$ige_intra_k, truth$ige_inter_k),
             B = c(truth$dge, truth$ige_intra_k, truth$ige_inter_k))

  # --- one fit ----------------------------------------------------------------
  fit_model <- function(p) {
    t0 <- Sys.time()
    V <- if (p$grid == "r3" && p$lambda == truth$lambda && p$d == truth$d) Vv
         else neighbourhood_Z(p$r, p$lambda, p$d, p$window)
    ZaD <- ZdA; ZaI <- V$AA; ZaX <- V$BA; ZbD <- ZdB; ZbI <- V$BB; ZbX <- V$AB
    random <- ~ str(~ mm(ZaD) + mm(ZaI) + mm(ZaX), struct = "us", name = "gA") +
                str(~ mm(ZbD) + mm(ZbI) + mm(ZbX), struct = "us", name = "gB")
    environment(random) <- environment()
    f <- tryCatch(rx_reml(y ~ sp_tray, random = random, residual = ~ dsum(~ units | species), data = d,
                          backend = "cpu", verbose = FALSE, hessian = FALSE, pev = FALSE),
                  error = function(e) conditionMessage(e))
    base <- data.frame(p, n = n, logLik = NA_real_, conv = NA, seconds = NA_real_, error = "")
    if (is.character(f)) { base$error <- substr(gsub("\n", " ", f), 1, 200); return(base) }
    base$logLik <- f$logLik
    base$conv <- isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement))
    v <- summary(f)$varcomp
    k_f <- list(A = c(1, kexp(V$AA, lA), kexp(V$BA, lB)), B = c(1, kexp(V$BB, lB), kexp(V$AB, lA)))
    s_f <- c(A = mean(Matrix::rowSums(V$AA)[lA]), B = mean(Matrix::rowSums(V$BB)[lB]))
    eff <- c("dge", "intra", "inter")
    for (sp in c("A", "B")) {
      nm <- paste0("g", sp)
      S <- matrix(NA_real_, 3, 3)
      for (a in 1:3) for (b in 1:a) {
        lab <- sprintf("%s[%d,%d]", nm, a, b)
        if (lab %in% rownames(v)) S[a, b] <- S[b, a] <- v[lab, "component"]
      }
      Bl <- f$blups[[nm]]
      Bl <- Bl[if (sp == "A") gA else gB, , drop = FALSE]
      Bl <- cbind(Bl, Bl[, 1] + s_f[[sp]] * Bl[, 2])
      for (a in 1:3) {
        base[[sprintf("vk_%s_%s", eff[a], sp)]] <- S[a, a] * k_f[[sp]][a]
        base[[sprintf("err_%s_%s", eff[a], sp)]] <- log10(S[a, a] * k_f[[sp]][a] / vk[[sp]][a])
      }
      base[[sprintf("rho_dge_intra_%s", sp)]] <- S[2, 1] / sqrt(S[1, 1] * S[2, 2])
      base[[sprintf("rho_dge_inter_%s", sp)]] <- S[3, 1] / sqrt(S[1, 1] * S[3, 3])
      base[[sprintf("rho_intra_inter_%s", sp)]] <- S[3, 2] / sqrt(S[2, 2] * S[3, 3])
      for (a in 1:4) base[[sprintf("acc_%s_%s", c(eff, "tbv")[a], sp)]] <- stats::cor(Bl[, a], true_vals[[sp]][, a])
    }
    base$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    base
  }

  pl <- plan_fits()
  if (identical(Sys.getenv("QUICK"), "1"))           # quick test: truth, naive, one wrong radius
    pl <- pl[(pl$grid == "r3" & pl$lambda == truth$lambda & pl$d == truth$d) | pl$grid == "naive" |
             (pl$grid == "r" & pl$r == 1), ]
  options(remlax.threads = 1L)
  cat(sprintf("task %d: n_rep %d, seed %d, %d trays, n = %d, %d fits on %d cores\n",
              i, t$n_rep, t$seed, n_tray, n, nrow(pl), cores))
  cat(sprintf("true k: AA %.3f BA %.3f BB %.3f AB %.3f\n", k_true[1], k_true[2], k_true[3], k_true[4]))
  t_task <- Sys.time()
  res <- parallel::mclapply(seq_len(nrow(pl)), function(j) fit_model(pl[j, ]), mc.cores = cores,
                            mc.preschedule = FALSE)
  ok <- vapply(res, is.data.frame, TRUE)
  if (!all(ok)) {
    cat("lost fits:", which(!ok), "\n")
    res[!ok] <- lapply(which(!ok), function(j) data.frame(pl[j, ], n = n, logLik = NA_real_, conv = NA,
                                                         seconds = NA_real_, error = "process lost"))
  }
  nm_all <- unique(unlist(lapply(res, names)))
  res <- do.call(rbind, lapply(res, function(x) { x[setdiff(nm_all, names(x))] <- NA; x[nm_all] }))
  res <- cbind(n_rep = t$n_rep, seed = t$seed, res)
  utils::write.csv(res, file.path(dir, sprintf("task_%03d.csv", i)), row.names = FALSE)
  cat(sprintf("finished in %.1f min; %d failure(s), %d not converged\n",
              as.numeric(difftime(Sys.time(), t_task, units = "mins")), sum(is.na(res$logLik)),
              sum(!res$conv, na.rm = TRUE)))
  print(res[res$grid != "r3" | (res$lambda == truth$lambda), c("grid", "r", "lambda", "d", "logLik",
          "conv", "acc_intra_A", "acc_inter_A", "err_intra_A", "err_inter_A", "seconds")], digits = 3)
  quit(save = "no")
}

if (mode == "summary") {
  dir <- a0[2]
  source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))),
                   "test_noyau_deux_especes_bilan.R"))
}
