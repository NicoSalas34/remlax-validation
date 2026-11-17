# ==============================================================================
# diag_qtl_power.R -- detection of the additive major QTL: one-step against
# two-step GWAS, on the same simulated wheat x alfalfa mixture
# ------------------------------------------------------------------------------
#   Rscript diag_qtl_power.R --rep 1 --out qtl_diag [--pve_major 0.10] [--r_g -0.6]
#                            [--n_lines 181 --n_fam 106 --var_gxg 0.05] [--backend cpu]
#
# Simulation and null model as in run_cascade.R (scenario S2, sum framework,
# interaction null model). Each target (DGE, IGE wheat, IGE alfalfa) is tested
# by seven methods:
#   one-step   p3d    rx_scan at the variances of the null (fixed V)
#              loco   rx_scan_loco / per-chromosome nulls, K without the chromosome
#              exact  refit with the SNP as a fixed effect, top 10 SNPs of p3d
#   two-step   blup_k   BLUPs of the null, line-level model with K
#              blup_ols BLUPs of the null, regression without K
#              blue_k   BLUEs (genotype effects fixed, spatial and GxG random), with K
#              blue_ols BLUEs, regression without K
# Every method also tests the causal SNPs themselves (not on the chip), which
# separates the power at the QTL from the loss through linkage disequilibrium.
# Found: Bonferroni 5 % over the chip SNPs; also after genomic control.
# ==============================================================================
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = NULL) { i <- match(paste0("--", k), a)
  if (is.na(i) || i == length(a)) d else a[i + 1L] }
REP <- as.integer(arg("rep", "1")); OUT <- arg("out", "qtl_diag"); BACKEND <- arg("backend", "cpu")
N_LINES <- as.integer(arg("n_lines", "181")); N_FAM <- as.integer(arg("n_fam", "106"))
VAR_GXG <- as.numeric(arg("var_gxg", "0.05")); PVE_MAJOR <- as.numeric(arg("pve_major", "0.10"))
R_G <- as.numeric(arg("r_g", "-0.6")); WIN <- 5e6; CHIP <- 0.5
DO_LOCO <- arg("loco", "yes") == "yes"; DO_EXACT <- arg("exact", "yes") == "yes"
DO_TWO <- arg("two_step", "yes") == "yes"; EXACT_TOP <- as.integer(arg("exact_top", "10"))
# design size: 12 bins of 160 wheat cells (1,920 plants); more bins = more plants per genotype
BINS <- as.integer(arg("bins", "12")); N_MAJOR <- as.integer(arg("n_major", "3"))
# null model: "interaction" (DGE x IGE kernels, as run_cascade.R) or "additive" (interaction left in the residual)
NULL_MODEL <- arg("null", "interaction")
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
suppressPackageStartupMessages({ library(Matrix); library(remlax) })
source(file.path(HERE, "sim_cogwas.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
T0 <- Sys.time()
say <- function(...) cat(sprintf(...), "\n", sep = "")
size <- if (BINS == 12L && N_MAJOR == 3L) "" else sprintf("B%03d_M%d_", BINS, N_MAJOR)
if (NULL_MODEL != "interaction") size <- paste0(size, "add_")
tag <- sprintf("L%03d_q%.2f_r%+.1f_g%.2f_%srep%03d", N_LINES, PVE_MAJOR, R_G, VAR_GXG, size, REP)
seed <- REP + 50000L + 7L * N_LINES + round(1000 * PVE_MAJOR) + round(100 * R_G) +
  13L * (BINS - 12L) + 17L * (N_MAJOR - 3L)
save_csv <- function(x, k) utils::write.csv(cbind(tag = tag, rep = REP, n_lines = N_LINES, pve_major = PVE_MAJOR,
                                                  r_g_target = R_G, var_gxg = VAR_GXG, bins = BINS,
                                                  n_major = N_MAJOR, x),
                                            file.path(OUT, sprintf("%s_%s.csv", tag, k)), row.names = FALSE)

# ---------------------------------------------------------------- simulation (as run_cascade.R)
par <- sim_scenario("S2")
par$pve_major <- PVE_MAJOR; par$r_g <- R_G
par$var_gxg_intra <- VAR_GXG; par$var_gxg_inter <- VAR_GXG
par$n_major_dir <- N_MAJOR; par$n_major_ige <- N_MAJOR
# enough alfalfa plants for the 160 alfalfa cells of each bin (2,120 by default, as run_cascade.R)
n_alf <- max(2120, ceiling(1.1 * 160 * BINS))
gen <- sim_genomes(n_lines = N_LINES, n_fam = N_FAM, plants_per_fam = as.integer(ceiling(n_alf / N_FAM)),
                   n_causal_wheat = 1500L, n_causal_alf = 700L, chip_frac = CHIP, seed = seed)
des <- sim_design(gen, n_bins = BINS, rank_intra = 1L, rank_inter = 1L, scale = "sum", seed = seed)
ph <- sim_pheno_full(gen, des, par, seed = seed)
tr <- ph$truth; vp <- tr$var_parts; st <- tr$sets
cw <- des$cells[des$wheat, ]; n <- length(ph$y)
say("%s: %d wheat plants, realised r_g %.3f, phenotypic variance %.2f", tag, n, tr$realised[["r_g"]], stats::var(ph$y))

lines <- rownames(gen$wheat$G); fams <- rownames(gen$alfalfa$G_fam)
Zd <- sparseMatrix(seq_len(n), match(cw$line, lines), x = 1, dims = c(n, length(lines)))
alf <- des$cells[des$alfalfa, ]
Zfa <- sparseMatrix(seq_len(nrow(alf)), match(alf$family, fams), x = 1, dims = c(nrow(alf), length(fams)))
Zw <- des$W_intra %*% Zd; Zb <- des$W_inter %*% Zfa
maf_ok <- function(M) { f <- colMeans(M) / 2; pmin(f, 1 - f) >= 0.05 & apply(M, 2, stats::sd) > 0 }
ctr <- function(M) sweep(M, 2, colMeans(M), "-")
Mw <- gen$wheat$G[, gen$wheat$chip]; mw <- gen$wheat$map[gen$wheat$chip, ]
k <- maf_ok(Mw); Mw <- ctr(Mw[, k]); mw <- mw[k, ]
Ma <- gen$alfalfa$G_fam[, gen$alfalfa$chip]; ma <- gen$alfalfa$map[gen$alfalfa$chip, ]
k <- maf_ok(Ma); Ma <- ctr(Ma[, k]); ma <- ma[k, ]
kin <- function(Mc) { K <- tcrossprod(Mc); K <- K / mean(diag(K)); 0.98 * K + 0.02 * diag(nrow(K)) }
Kw <- kin(Mw); Ka <- kin(Ma)
dimnames(Kw) <- list(lines, lines); dimnames(Ka) <- list(fams, fams)
t_k <- Sys.time()
KERN <- if (NULL_MODEL == "interaction")
  list(k_intra = rx_cogwas_kernel(Zd, Kw, Zw, Kw, normalize = FALSE),
       k_inter = rx_cogwas_kernel(Zd, Kw, Zb, Ka, normalize = FALSE)) else list()
fac <- function(x) sparse.model.matrix(~ x - 1, data.frame(x = factor(x)))
Zsp <- list(bin = fac(cw$bin), bin_row = fac(paste(cw$bin, cw$row)), bin_col = fac(paste(cw$bin, cw$col)))
other_terms <- function() c(lapply(names(Zsp), function(nm) rx_term(nm, Z = Zsp[[nm]], t = 1L)),
                            lapply(names(KERN), function(nm) rx_term(nm, Z = Diagonal(n), K = KERN[[nm]], t = 1L)))
build2 <- function(K, KA = Ka, X_sup = NULL)
  rx_model(ph$y, cbind(matrix(1, n, 1), X_sup),
           terms = c(list(rx_term("gen", Z = cbind(Zd, Zw), K = K, t = 2L, struct = "us"),
                          rx_term("ige_inter", Z = Zb, K = KA, t = 1L)), other_terms()),
           residual = rx_residual())
build <- function(K, X_sup = NULL) build2(K, Ka, X_sup)
FIT_ARGS <- list(backend = BACKEND, hessian = FALSE, blups = TRUE, verbose = FALSE)
fit_of <- function(m, theta_init = NULL) do.call(rx_fit, c(list(m, theta_init = theta_init), FIT_ARGS))
conv_of <- function(f) isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement))

# markers: chip, and the causal major SNPs (oracle), same centring
mk <- function(G, s, map, pre) { if (!length(s)) return(NULL)
  M <- ctr(G[, s, drop = FALSE]); colnames(M) <- sprintf("%s%d", pre, seq_along(s))
  list(M = M, map = data.frame(snp = colnames(M), chr = map$chr[s], pos = map$pos[s])) }
OR <- list(DGE = mk(gen$wheat$G, st$major_dir, gen$wheat$map, "causal_dge"),
           IGE_wheat = mk(gen$wheat$G, st$major_ige, gen$wheat$map, "causal_igew"),
           IGE_alfalfa = mk(gen$alfalfa$G_fam, st$major_inter, gen$alfalfa$map, "causal_igea"))
map_w <- data.frame(snp = colnames(Mw), chr = mw$chr, pos = mw$pos)
map_a <- data.frame(snp = colnames(Ma), chr = ma$chr, pos = ma$pos)
TG <- list(DGE = list(M = Mw, map = map_w, K = Kw, inc = "dir", test = "dir", sp = "wheat"),
           IGE_wheat = list(M = Mw, map = map_w, K = Kw, inc = "ind", test = "ind", sp = "wheat"),
           IGE_alfalfa = list(M = Ma, map = map_a, K = Ka, inc = "ind", test = "ind", sp = "alfalfa"))
# realised share of each major in its genotypic component (line or family level)
g_true <- list(DGE = tr$g$dge, IGE_wheat = tr$g$ige, IGE_alfalfa = as.numeric(tr$g$inter_fam[fams]))
beta_true <- list(DGE = tr$beta$major_dir, IGE_wheat = tr$beta$major_ige, IGE_alfalfa = tr$beta$major_inter)
share <- lapply(names(OR), function(t) apply(OR[[t]]$M, 2, stats::var) * beta_true[[t]]^2 / stats::var(g_true[[t]]))
names(share) <- names(OR)
r2_tag <- function(t, j) { w <- which(TG[[t]]$map$chr == OR[[t]]$map$chr[j] & abs(TG[[t]]$map$pos - OR[[t]]$map$pos[j]) <= WIN)
  if (!length(w)) NA else max(stats::cor(OR[[t]]$M[, j], TG[[t]]$M[, w])^2) }

# ---------------------------------------------------------------- scoring
chi2 <- function(p) stats::qchisq(p, 1, lower.tail = FALSE)
QTL <- list(); CAL <- list()
score <- function(method, t, p_chip, p_causal, chip_snp = TG[[t]]$map$snp) {
  map <- TG[[t]]$map[match(chip_snp, TG[[t]]$map$snp), ]
  ok <- is.finite(p_chip); thr <- 0.05 / sum(ok)
  lam <- stats::median(chi2(p_chip[ok])) / stats::qchisq(0.5, 1)
  om <- OR[[t]]$map
  in_w <- function(j) ok & map$chr == om$chr[j] & abs(map$pos - om$pos[j]) <= WIN
  rows <- lapply(seq_len(nrow(om)), function(j) {
    w <- which(in_w(j)); b <- if (length(w)) min(p_chip[w]) else NA
    data.frame(method = method, target = t, qtl = j, share = share[[t]][j], r2_tag = r2_tag(t, j),
               p_tag = b, found_tag = isTRUE(b < thr), found_tag_gc = isTRUE(chi2(b) / max(lam, 1) > chi2(thr)),
               p_causal = if (length(p_causal) >= j) p_causal[j] else NA,
               found_causal = isTRUE(p_causal[j] < thr))
  })
  causal <- Reduce(`|`, lapply(seq_len(nrow(om)), in_w), rep(FALSE, nrow(map)))
  sig <- which(ok & p_chip < thr & !causal)
  fp <- if (length(sig)) length(unique(paste(map$chr[sig], floor(map$pos[sig] / WIN)))) else 0L
  QTL[[length(QTL) + 1L]] <<- do.call(rbind, rows)
  CAL[[length(CAL) + 1L]] <<- data.frame(method = method, target = t, lambda = lam, false_peaks = fp, n_chip = sum(ok))
}
# per-SNP regression without K (two-step "ols")
p_ols <- function(y, M) { r <- drop(stats::cor(y, M)); df <- length(y) - 2
  2 * stats::pt(-abs(r * sqrt(df / (1 - r^2))), df) }
# line- or family-level mixed model y = 1 + u(K) + e, then rx_scan
two_step_k <- function(y, t) {
  m <- rx_model(y, matrix(1, length(y), 1), terms = list(rx_term("g", Z = Diagonal(length(y)), K = TG[[t]]$K, t = 1L)),
                residual = rx_residual())
  f <- fit_of(m)
  sc <- function(M, map) rx_scan(f, m, M, map = map, incidences = c(dir = "g"), tests = "dir", maf = NULL,
                                 backend = BACKEND, verbose = FALSE)$p_dir
  list(chip = sc(TG[[t]]$M, TG[[t]]$map), causal = sc(OR[[t]]$M, OR[[t]]$map), h2 = f$sigmas$g[1, 1] /
         (f$sigmas$g[1, 1] + as.numeric(f$sigma_res)[1]), conv = conv_of(f))
}
H2 <- list()
TIME <- list(); tick <- function(k, t1) TIME[[k]] <<- as.numeric(difftime(Sys.time(), t1, units = "secs"))

# ---------------------------------------------------------------- null model and one-step
tick("kernels", t_k)
t1 <- Sys.time(); mod <- build(Kw); fit <- fit_of(mod); tick("null", t1)
say("null model %.0f s, converged %s", TIME$null, conv_of(fit))
B <- fit$blups
blup <- list(DGE = B$gen[, 1], IGE_wheat = B$gen[, 2], IGE_alfalfa = B$ige_inter[, 1])
INC <- list(wheat = c(dir = "gen:1", ind = "gen:2"), alfalfa = c(ind = "ige_inter"))
scan1 <- function(f, m, t, M, map) {
  mk_ <- if (TG[[t]]$sp == "wheat") M else list(ind = M)
  rx_scan(f, m, mk_, map = map, incidences = INC[[TG[[t]]$sp]], tests = TG[[t]]$test, maf = NULL,
          backend = BACKEND, verbose = FALSE)[[paste0("p_", TG[[t]]$test)]]
}
t1 <- Sys.time()
P3D <- lapply(names(TG), function(t) list(chip = scan1(fit, mod, t, TG[[t]]$M, TG[[t]]$map),
                                          causal = scan1(fit, mod, t, OR[[t]]$M, OR[[t]]$map)))
names(P3D) <- names(TG); tick("p3d", t1)
for (t in names(TG)) score("p3d", t, P3D[[t]]$chip, P3D[[t]]$causal)

# LOCO: one null per chromosome, K of the species without that chromosome
if (DO_LOCO) {
  t1 <- Sys.time()
  loco_species <- function(sp) {
    Mc <- if (sp == "wheat") Mw else Ma; mc <- if (sp == "wheat") map_w else map_a
    ts <- names(TG)[vapply(TG, function(x) x$sp == sp, TRUE)]
    out <- setNames(lapply(ts, function(t) list(chip = rep(NA_real_, ncol(Mc)), causal = rep(NA_real_, ncol(OR[[t]]$M)))), ts)
    for (cc in unique(as.character(mc$chr))) {
      inside <- as.character(mc$chr) == cc
      Kc <- kin(Mc[, !inside, drop = FALSE])
      m_c <- if (sp == "wheat") build2(Kc, Ka) else build2(Kw, Kc)
      f_c <- fit_of(m_c, theta_init = fit$theta)
      for (t in ts) {
        out[[t]]$chip[inside] <- scan1(f_c, m_c, t, Mc[, inside, drop = FALSE], mc[inside, ])
        oc <- as.character(OR[[t]]$map$chr) == cc
        if (any(oc)) out[[t]]$causal[oc] <- scan1(f_c, m_c, t, OR[[t]]$M[, oc, drop = FALSE], OR[[t]]$map[oc, ])
      }
    }
    out
  }
  LOCO <- c(loco_species("wheat"), loco_species("alfalfa")); tick("loco", t1)
  for (t in names(TG)) score("loco", t, LOCO[[t]]$chip, LOCO[[t]]$causal)
  say("LOCO %.0f s", TIME$loco)
}

# exact refit (fixed V replaced by a refit with the SNP fixed), top 10 SNPs of p3d and the causal SNPs
if (DO_EXACT) {
  t1 <- Sys.time()
  for (t in names(TG)) {
    res_of <- function(M, map, p) { r <- map; r[[paste0("p_", TG[[t]]$test)]] <- p; r }
    ex <- function(M, map, p, top) rx_scan_exact(res_of(M, map, p), TG[[t]]$test, build, function(W) Kw, M, map,
                                                 INC[[TG[[t]]$sp]], top = top, loco = FALSE, batch = 0L,
                                                 theta_init = fit$theta,
                                                 backend = BACKEND, verbose = FALSE)
    e_chip <- tryCatch(ex(TG[[t]]$M, TG[[t]]$map, P3D[[t]]$chip, EXACT_TOP), error = function(e) NULL)
    e_cau <- tryCatch(ex(OR[[t]]$M, OR[[t]]$map, P3D[[t]]$causal, ncol(OR[[t]]$M)), error = function(e) NULL)
    p_chip <- P3D[[t]]$chip
    if (!is.null(e_chip)) p_chip[match(e_chip$snp, TG[[t]]$map$snp)] <- e_chip$p_exact
    p_cau <- if (is.null(e_cau)) rep(NA, ncol(OR[[t]]$M)) else e_cau$p_exact[match(OR[[t]]$map$snp, e_cau$snp)]
    if (is.null(e_chip) || is.null(e_cau)) say("  exact refit failed for %s", t)
    score("exact", t, p_chip, p_cau)
  }
  tick("exact", t1); say("exact refit %.0f s", TIME$exact)
}

# ---------------------------------------------------------------- two-step
if (DO_TWO) {
t1 <- Sys.time()
for (t in names(TG)) {
  y <- blup[[t]]
  r <- tryCatch(two_step_k(y, t), error = function(e) { say("  blup_k failed for %s: %s", t, conditionMessage(e)); NULL })
  if (!is.null(r)) { score("blup_k", t, r$chip, r$causal); H2[[length(H2) + 1L]] <- data.frame(method = "blup_k", target = t, h2_stage2 = r$h2, conv = r$conv) }
  score("blup_ols", t, p_ols(y, TG[[t]]$M), p_ols(y, OR[[t]]$M))
}
tick("blup", t1)
# stage 1 with the genotype effects fixed: X = [1, Zd(-1), Zw, Zb], spatial and GxG random
t1 <- Sys.time()
X1 <- as.matrix(cbind(1, Zd[, -1], Zw, Zb))
idx <- list(DGE = 1 + seq_len(ncol(Zd) - 1), IGE_wheat = ncol(Zd) + seq_len(ncol(Zw)),
            IGE_alfalfa = ncol(Zd) + ncol(Zw) + seq_len(ncol(Zb)))
rk <- qr(X1)$rank
say("stage 1: %d fixed columns, rank %d", ncol(X1), rk)
m1 <- rx_model(ph$y, X1, terms = other_terms(), residual = rx_residual())
f1 <- tryCatch(fit_of(m1), error = function(e) { say("  stage 1 failed: %s", conditionMessage(e)); NULL })
if (!is.null(f1)) {
  be <- as.numeric(f1$beta)
  blue <- list(DGE = c(0, be[idx$DGE]), IGE_wheat = be[idx$IGE_wheat], IGE_alfalfa = be[idx$IGE_alfalfa])
  for (t in names(TG)) {
    y <- blue[[t]]
    r <- tryCatch(two_step_k(y, t), error = function(e) { say("  blue_k failed for %s: %s", t, conditionMessage(e)); NULL })
    if (!is.null(r)) { score("blue_k", t, r$chip, r$causal); H2[[length(H2) + 1L]] <- data.frame(method = "blue_k", target = t, h2_stage2 = r$h2, conv = r$conv) }
    score("blue_ols", t, p_ols(y, TG[[t]]$M), p_ols(y, OR[[t]]$M))
  }
  ACC <- data.frame(target = names(TG),
                    blup = vapply(names(TG), function(t) stats::cor(blup[[t]], g_true[[t]]), 0),
                    blue = vapply(names(TG), function(t) stats::cor(blue[[t]], g_true[[t]]), 0))
  print(ACC, row.names = FALSE)
  save_csv(ACC, "acc")
}
tick("blue", t1)
}

qtl <- do.call(rbind, QTL); cal <- do.call(rbind, CAL)
save_csv(qtl, "qtl"); save_csv(cal, "cal")
if (length(H2)) save_csv(do.call(rbind, H2), "h2")
save_csv(data.frame(step = names(TIME), seconds = unlist(TIME)), "time")
save_csv(data.frame(null_converged = conv_of(fit), realised_r_g = tr$realised[["r_g"]], var_y = stats::var(ph$y),
                    var_dge = tr$realised[["var_dge"]], var_ige = tr$realised[["var_ige"]],
                    var_inter_fam = tr$realised[["var_inter_fam"]]), "truth")
agg <- aggregate(cbind(found_tag, found_causal) ~ method + target, qtl, sum)
print(merge(agg, cal[, c("method", "target", "lambda")]), row.names = FALSE)
say("%s written to %s (%.0f s)", tag, OUT, as.numeric(difftime(Sys.time(), T0, units = "secs")))
