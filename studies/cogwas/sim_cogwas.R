# ==============================================================================
# sim_cogwas.R -- simulated wheat x alfalfa mixtures for co-GWAS validation
# ------------------------------------------------------------------------------
# Kept outside remlax on purpose: remlax is a solver, simulation lives here.
#
# sim_genomes()  AlphaSimR genomes
#   wheat   : n_lines fixed lines, MaCS coalescent with the "WHEAT" history,
#             21 chromosomes; doses 0/2.
#   alfalfa : 2 x n_fam autotetraploid parents ("GENERIC" history, ploidy 4,
#             8 chromosomes); each family is the cross of two parents, and
#             every plant has its own genotype (dose 0..4). The chip only sees
#             the family's expected dose, (mother + father) / 4 on 0..2, as in
#             the real data (genotyped mothers, F2 families).
#   Causal sites are drawn among the segregating sites; the chip is drawn among
#   the others, so a QTL is seen only through markers in LD with it.
# sim_design()   12 bins of 16 rows x 20 columns, wheat in odd columns, alfalfa
#   in even ones; lines and families randomised; neighbours by Chebyshev
#   distance within a bin, summed (default, Bijma) or averaged.
# sim_pheno()    wheat phenotypes:
#   y = DGE + IGE_intra + IGE_inter + pair_intra + pair_inter + bin + row + col + e
#   with polygenic additive effects, pair QTL (focal dose at a x mean neighbour
#   dose at b, wheat or alfalfa neighbours), each part scaled to a share of the
#   phenotypic variance.
# ==============================================================================

suppressPackageStartupMessages(library(AlphaSimR))

#' AlphaSimR genomes of the wheat lines and the alfalfa families.
sim_genomes <- function(n_lines = 181L, n_fam = 106L, plants_per_fam = 20L,
                        seg_wheat = 2000L, seg_alf = 3000L, chip_frac = 0.5,
                        n_causal_wheat = 400L, n_causal_alf = 400L, seed = 1L) {
  set.seed(seed)
  # wheat: fixed lines
  fw <- runMacs(nInd = n_lines, nChr = 21L, segSites = seg_wheat, inbred = TRUE, species = "WHEAT")
  SPw <- SimParam$new(fw)
  lines <- newPop(fw, simParam = SPw)
  Gw <- pullSegSiteGeno(lines, simParam = SPw)                 # lines x sites, 0/2
  rownames(Gw) <- sprintf("W%03d", seq_len(n_lines))
  mapw <- .sim_map(SPw, "W")
  # alfalfa: tetraploid parents, one cross per family
  fa <- runMacs(nInd = 2L * n_fam, nChr = 8L, segSites = seg_alf, ploidy = 4L, species = "GENERIC")
  SPa <- SimParam$new(fa)
  par <- newPop(fa, simParam = SPa)
  plan <- cbind(seq(1L, 2L * n_fam, 2L), seq(2L, 2L * n_fam, 2L))
  kids <- makeCross(par, plan, nProgeny = plants_per_fam, simParam = SPa)
  Ga_plant <- pullSegSiteGeno(kids, simParam = SPa)            # plants x sites, 0..4
  fam_of_plant <- rep(sprintf("F%03d", seq_len(n_fam)), each = plants_per_fam)
  Gp <- pullSegSiteGeno(par, simParam = SPa)
  Ga_fam <- (Gp[plan[, 1], ] + Gp[plan[, 2], ]) / 4             # expected dose, 0..2
  rownames(Ga_fam) <- sprintf("F%03d", seq_len(n_fam))
  mapa <- .sim_map(SPa, "A")
  # causal sites, then the chip among the others (polymorphic only)
  pick <- function(G, n) sample(which(apply(G, 2, stats::var) > 0), n)
  cw <- pick(Gw, n_causal_wheat); ca <- pick(Ga_plant, n_causal_alf)
  chip <- function(G, causal, frac) {
    cand <- setdiff(which(apply(G, 2, stats::var) > 0), causal)
    sort(sample(cand, round(frac * length(cand))))
  }
  list(wheat = list(G = Gw, map = mapw, causal = cw, chip = chip(Gw, cw, chip_frac)),
       alfalfa = list(G_plant = Ga_plant, G_fam = Ga_fam, fam_of_plant = fam_of_plant,
                      map = mapa, causal = ca, chip = chip(Ga_fam, ca, chip_frac)))
}

# Site names, chromosome and pseudo-physical position (genetic position x 1e8).
.sim_map <- function(SP, prefix) {
  gm <- SP$genMap
  do.call(rbind, lapply(seq_along(gm), function(k)
    data.frame(snp = sprintf("%s%02d_%05d", prefix, k, seq_along(gm[[k]])), chr = k,
               pos = round(as.numeric(gm[[k]]) * 1e8), stringsAsFactors = FALSE)))
}

#' Field design and neighbour weights.
sim_design <- function(gen, n_bins = 12L, n_row = 16L, n_col = 20L, rank_intra = 1L,
                       rank_inter = 1L, scale = c("sum", "mean"), seed = 1L) {
  scale <- match.arg(scale)
  set.seed(seed)
  cells <- expand.grid(row = seq_len(n_row), col = seq_len(n_col), bin = seq_len(n_bins))
  cells$species <- ifelse(cells$col %% 2L == 1L, "wheat", "alfalfa")
  wi <- which(cells$species == "wheat"); ai <- which(cells$species == "alfalfa")
  lines <- rownames(gen$wheat$G); n_alf_plants <- nrow(gen$alfalfa$G_plant)
  cells$line <- NA_character_; cells$line[wi] <- sample(rep_len(sample(lines), length(wi)))
  cells$plant <- NA_integer_
  cells$plant[ai] <- sample(n_alf_plants, length(ai), replace = length(ai) > n_alf_plants)
  cells$family <- NA_character_; cells$family[ai] <- gen$alfalfa$fam_of_plant[cells$plant[ai]]
  W <- function(from, to, rank) .sim_weights(cells, from, to, rank, scale)
  list(cells = cells, wheat = wi, alfalfa = ai,
       W_intra = W(wi, wi, rank_intra),       # wheat neighbours of each wheat plant
       W_inter = W(wi, ai, rank_inter))       # alfalfa neighbours of each wheat plant
}

# Neighbours of kind `to` within Chebyshev distance `rank` (same bin), for each
# plant of `from`: their sum (Bijma) or their mean. Rows of plants without
# neighbour stay zero.
.sim_weights <- function(cells, from, to, rank, scale = "sum") {
  i <- integer(0); j <- integer(0)
  for (b in unique(cells$bin)) {
    f <- from[cells$bin[from] == b]; t <- to[cells$bin[to] == b]
    d <- pmax(abs(outer(cells$row[f], cells$row[t], "-")), abs(outer(cells$col[f], cells$col[t], "-")))
    hit <- which(d > 0 & d <= rank, arr.ind = TRUE)
    i <- c(i, match(f[hit[, 1]], from)); j <- c(j, match(t[hit[, 2]], to))
  }
  W <- Matrix::sparseMatrix(i, j, x = 1, dims = c(length(from), length(to)))
  if (scale == "sum") return(W)
  rs <- Matrix::rowSums(W); rs[rs == 0] <- 1
  Matrix::Diagonal(x = 1 / rs) %*% W
}

#' Wheat phenotypes. Shares of the phenotypic variance: dge, ige_intra, ige_inter
#' (polygenic, on the causal sites), pair_intra and pair_inter (one pair QTL each,
#' 0 for none), spatial (bin + row + column), the rest is residual.
sim_pheno <- function(gen, des, dge = 0.3, ige_intra = 0.05, ige_inter = 0.05,
                      pair_intra = 0, pair_inter = 0, spatial = 0.1, seed = 1L) {
  set.seed(seed)
  c0 <- des$cells
  Gw <- gen$wheat$G; Ga <- gen$alfalfa$G_plant
  cw <- gen$wheat$causal; ca <- gen$alfalfa$causal
  centre <- function(x) x - mean(x)
  line_of <- c0$line[des$wheat]; plant_of <- c0$plant[des$alfalfa]
  # polygenic values of lines and alfalfa plants
  bw_d <- rnorm(length(cw)); bw_i <- rnorm(length(cw)); ba <- rnorm(length(ca))
  g_dge <- drop(scale(Gw[, cw], scale = FALSE) %*% bw_d)
  g_igei <- drop(scale(Gw[, cw], scale = FALSE) %*% bw_i)
  g_igea <- drop(scale(Ga[, ca], scale = FALSE) %*% ba)
  parts <- list(
    dge = g_dge[line_of],
    ige_intra = as.numeric(des$W_intra %*% g_igei[line_of]),
    ige_inter = as.numeric(des$W_inter %*% g_igea[plant_of]))
  # pair QTL: focal causal site a x neighbour causal site b
  truth <- list()
  if (pair_intra > 0) {
    a <- sample(cw, 1); b <- sample(setdiff(cw, a), 1)
    xa <- centre(Gw[line_of, a]); xb <- centre(Gw[line_of, b])
    parts$pair_intra <- xa * as.numeric(des$W_intra %*% xb)
    truth$intra <- c(a = unname(a), b = unname(b))
  }
  if (pair_inter > 0) {
    a <- sample(cw, 1); b <- sample(ca, 1)
    xa <- centre(Gw[line_of, a]); zb <- centre(Ga[plant_of, b])
    parts$pair_inter <- xa * as.numeric(des$W_inter %*% zb)
    truth$inter <- c(a = unname(a), b = unname(b))
  }
  cw_ <- des$cells[des$wheat, ]
  parts$spatial <- rnorm(max(cw_$bin))[cw_$bin] +
    rnorm(max(cw_$bin) * 16)[(cw_$bin - 1) * 16 + cw_$row] +
    rnorm(max(cw_$bin) * 20)[(cw_$bin - 1) * 20 + cw_$col]
  share <- c(dge = dge, ige_intra = ige_intra, ige_inter = ige_inter, pair_intra = pair_intra,
             pair_inter = pair_inter, spatial = spatial)
  share <- share[names(parts)]
  if (sum(share) >= 1) stop("sim_pheno: the shares add up to 1 or more.")
  y <- rnorm(length(des$wheat), sd = sqrt(1 - sum(share)))
  for (k in names(parts)) {
    v <- stats::var(parts[[k]])
    if (share[[k]] > 0 && v > 0) y <- y + parts[[k]] * sqrt(share[[k]] / v)
  }
  list(y = y, truth = truth, share = share)
}

# ==============================================================================
# sim_pheno_full() -- wheat phenotypes in the sum framework (Bijma)
# ------------------------------------------------------------------------------
# Each wheat plant receives the SUM of its neighbours' effects (W not
# normalised), with the parameters of Salas et al. 2026
# (simulate_multitrait_neighbors.R): var_dge = 1, var(IGE) = alpha x var_dge
# per neighbour, r(DGE, IGE intra) = -0.6, and the environmental variance set
# by the individual DGE heritability H2 = var_dge / (var_dge + var_env).
# Additive genetic values (what the mixed model estimates):
#   g_dge   (lines)          polygenic + major direct QTL
#   g_ige   (lines)          polygenic + major indirect QTL, summed over wheat neighbours
#   g_inter (alfalfa plants) polygenic + major QTL, summed over alfalfa neighbours
# Partial pleiotropy: the polygenic causal SNPs are split into direct only,
# indirect only and shared. On shared SNPs the indirect effect is
# rho z_d + sqrt(1 - rho^2) z_i; rho is set by bisection so that the REALISED
# correlation between g_dge and g_ige (majors and LD included) equals r_g.
# DGE x IGE interaction: one effect gamma(f, g) per ordered pair (focal line
# f, neighbour genotype g), summed over the neighbours. Its covariance over
# pairs is K x K (Kronecker); on the plants it is (Z K Z') o (Zw K Zw').
# gamma = sum over explicit pairs of causal loci (a, b) of theta_ab x_a(f) x_b(g),
# with, at the allelic level, major QTL: similarity (a = b) and pair (a != b)
# within wheat, pair between wheat and alfalfa. var(gamma) is per neighbour,
# in units of var_dge. All causal loci are drawn among gen$*$causal, hence off
# the chip.
# ==============================================================================
sim_par_default <- function() list(
  var_dge = 1, alpha_intra = 0.10, alpha_inter = 0.10, r_g = -0.6, h2 = 0.30,
  env_share = c(bin = 0.10, row = 0.05, col = 0.05),    # of var_env, the rest is residual
  n_poly = 1000L, frac_shared = 0.7, frac_dir_only = 0.15, n_poly_inter = 500L,
  n_major_dir = 3L, n_major_ige = 3L, n_major_inter = 2L, pve_major = 0.10,
  var_gxg_intra = 0.05, var_gxg_inter = 0.05,           # per neighbour, x var_dge
  n_gxg_intra = 2000L, n_gxg_inter = 2000L, poly_gxg = TRUE,
  n_sim = 2L, pve_sim = 0.10, pve_pair_intra = 0.25, pve_pair_inter = 0.30)

#' Scenarios: S0 additive only; S1 adds the allelic interaction QTL (similarity
#' and pairs) without polygenic interaction; S2 adds the polygenic DGE x IGE
#' interaction background.
sim_scenario <- function(name) {
  p <- sim_par_default()
  if (!name %in% c("S0", "S1", "S2")) stop("unknown scenario ", name)
  if (name == "S0") { p$var_gxg_intra <- 0; p$var_gxg_inter <- 0 }
  if (name == "S1") p$poly_gxg <- FALSE
  p
}

# Genotype-pair interaction matrix (focal genotypes x neighbour genotypes):
# majors as outer products scaled to pve x v, polygenic part to the rest.
.sim_gxg <- function(XF, XN, majors, pves, pool_f, pool_n, n_poly, v, poly) {
  if (v <= 0) return(list(G = matrix(0, nrow(XF), nrow(XN)), beta = numeric(0)))
  scale_mat <- function(G, w) G * sqrt(w / stats::var(as.vector(G)))
  G <- matrix(0, nrow(XF), nrow(XN)); beta <- numeric(0)
  for (k in seq_along(majors)) {
    ab <- majors[[k]]; o <- outer(XF[, ab[1]], XN[, ab[2]])
    s <- sample(c(-1, 1), 1) * sqrt(pves[k] * v / stats::var(as.vector(o)))
    G <- G + s * o; beta <- c(beta, s)
  }
  if (poly && sum(pves) < 1) {
    a <- sample(pool_f, n_poly, TRUE); b <- sample(pool_n, n_poly, TRUE)
    G <- G + scale_mat(XF[, a, drop = FALSE] %*% (stats::rnorm(n_poly) * t(XN[, b, drop = FALSE])),
                       v * (1 - sum(pves)))
  }
  list(G = G, beta = beta)
}

sim_pheno_full <- function(gen, des, par = sim_par_default(), seed = 1L) {
  set.seed(seed)
  cells <- des$cells; cw <- cells[des$wheat, ]; n <- nrow(cw)
  line_of <- match(cw$line, rownames(gen$wheat$G))
  plant_of <- cells$plant[des$alfalfa]
  Gw <- gen$wheat$G; Ga <- gen$alfalfa$G_plant
  Xw <- sweep(Gw, 2, colMeans(Gw), "-"); Xa <- sweep(Ga, 2, colMeans(Ga), "-")
  var_i <- par$alpha_intra * par$var_dge; var_a <- par$alpha_inter * par$var_dge
  pool <- sample(gen$wheat$causal)
  take <- function(k) { if (k == 0) return(integer(0)); s <- pool[seq_len(k)]; pool <<- pool[-seq_len(k)]; unname(s) }
  maj_d <- take(par$n_major_dir); maj_i <- take(par$n_major_ige)
  sim_s <- take(if (par$var_gxg_intra > 0) par$n_sim else 0L)
  pair_intra <- take(2L); pair_inter_a <- take(1L)
  if (length(pool) < par$n_poly) stop("sim_pheno_full: not enough causal wheat sites.")
  poly <- pool[seq_len(par$n_poly)]
  n_sh <- round(par$frac_shared * par$n_poly); n_d <- round(par$frac_dir_only * par$n_poly)
  sh <- poly[seq_len(n_sh)]; d_only <- poly[n_sh + seq_len(n_d)]
  i_only <- poly[-seq_len(n_sh + n_d)]
  # majors: effect so that var(X_s beta) = pve x target variance, random sign
  major <- function(X, s, v) sample(c(-1, 1), length(s), TRUE) * sqrt(par$pve_major * v / apply(X[, s, drop = FALSE], 2, stats::var))
  b_maj_d <- major(Xw, maj_d, par$var_dge); b_maj_i <- major(Xw, maj_i, var_i)
  g_maj_d <- drop(Xw[, maj_d, drop = FALSE] %*% b_maj_d)
  g_maj_i <- drop(Xw[, maj_i, drop = FALSE] %*% b_maj_i)
  scale_to <- function(g, v) g * sqrt(v / stats::var(g))
  # The reachable range of the realised r_g depends on the draws (LD between
  # sets, majors). Out of reach: new polygenic draws, up to 50 attempts.
  for (attempt in 1:50) {
    zd_sh <- stats::rnorm(n_sh); zi_sh <- stats::rnorm(n_sh)
    bd <- c(zd_sh, stats::rnorm(length(d_only)))
    g_d <- scale_to(drop(Xw[, c(sh, d_only)] %*% bd), par$var_dge * (1 - par$n_major_dir * par$pve_major)) + g_maj_d
    zi_only <- stats::rnorm(length(i_only))
    g_i_of <- function(rho) {
      bi <- c(rho * zd_sh + sqrt(1 - rho^2) * zi_sh, zi_only)
      scale_to(drop(Xw[, c(sh, i_only)] %*% bi), var_i * (1 - par$n_major_ige * par$pve_major)) + g_maj_i
    }
    r_of <- function(rho) stats::cor(g_d, g_i_of(rho))
    lo <- -0.999; hi <- 0.999; r_lo <- r_of(lo); r_hi <- r_of(hi)
    if (par$r_g >= min(r_lo, r_hi) && par$r_g <= max(r_lo, r_hi)) break
    if (attempt == 50)
      stop(sprintf("sim_pheno_full: r_g = %.2f not reachable with this pleiotropy (last range %.2f to %.2f).",
                   par$r_g, min(r_lo, r_hi), max(r_lo, r_hi)))
  }
  for (it in 1:60) {                       # r_of is monotone in rho
    mid <- (lo + hi) / 2
    if ((r_of(mid) - par$r_g) * (r_lo - par$r_g) > 0) { lo <- mid; r_lo <- r_of(lo) } else hi <- mid
  }
  rho <- (lo + hi) / 2; g_i <- g_i_of(rho)
  # LD between majors and polygenic SNPs moves the total off target: rescale
  # (correlation unchanged), majors' effects with it
  c_d <- sqrt(par$var_dge / stats::var(g_d)); c_i <- sqrt(var_i / stats::var(g_i))
  g_d <- g_d * c_d; b_maj_d <- b_maj_d * c_d; g_i <- g_i * c_i; b_maj_i <- b_maj_i * c_i
  # IGE inter: alfalfa plants
  pool_a <- sample(gen$alfalfa$causal)
  maj_a <- unname(pool_a[seq_len(par$n_major_inter)])
  poly_a <- unname(pool_a[par$n_major_inter + seq_len(par$n_poly_inter)])
  pair_inter_b <- unname(pool_a[par$n_major_inter + par$n_poly_inter + 1L])
  b_maj_a <- major(Xa, maj_a, var_a)
  g_a <- scale_to(drop(Xa[, poly_a] %*% stats::rnorm(length(poly_a))), var_a * (1 - par$n_major_inter * par$pve_major)) +
    drop(Xa[, maj_a, drop = FALSE] %*% b_maj_a)
  c_a <- sqrt(var_a / stats::var(g_a)); g_a <- g_a * c_a; b_maj_a <- b_maj_a * c_a
  # DGE x IGE interaction, genotype-pair matrices
  gi <- .sim_gxg(Xw, Xw, c(lapply(sim_s, function(s) c(s, s)), list(pair_intra)),
                 c(rep(par$pve_sim, length(sim_s)), par$pve_pair_intra),
                 gen$wheat$causal, gen$wheat$causal, par$n_gxg_intra, par$var_gxg_intra, par$poly_gxg)
  ge <- .sim_gxg(Xw, Xa, list(c(pair_inter_a, pair_inter_b)), par$pve_pair_inter,
                 gen$wheat$causal, gen$alfalfa$causal, par$n_gxg_inter, par$var_gxg_inter, par$poly_gxg)
  fam <- gen$alfalfa$fam_of_plant
  G_inter_fam <- t(apply(ge$G, 1, function(r) tapply(r, fam, mean)))       # lines x families
  # phenotype contributions on wheat plants (sums over neighbours)
  nl <- nrow(Gw); Zd <- Matrix::sparseMatrix(seq_len(n), line_of, x = 1, dims = c(n, nl))
  Zp <- Matrix::sparseMatrix(seq_along(plant_of), plant_of, x = 1, dims = c(length(plant_of), nrow(Ga)))
  Zw <- des$W_intra %*% Zd; Zn <- des$W_inter %*% Zp
  parts <- list(dge = g_d[line_of], ige = as.numeric(Zw %*% g_i), inter = as.numeric(Zn %*% g_a),
                gxg_intra = Matrix::rowSums(as.matrix(Zd %*% gi$G) * as.matrix(Zw)),
                gxg_inter = Matrix::rowSums(as.matrix(Zd %*% ge$G) * as.matrix(Zn)))
  var_env <- par$var_dge * (1 - par$h2) / par$h2
  es <- par$env_share; nb <- max(cw$bin)
  e_bin <- stats::rnorm(nb, sd = sqrt(es[["bin"]] * var_env))
  e_row <- stats::rnorm(nb * 16, sd = sqrt(es[["row"]] * var_env))
  e_col <- stats::rnorm(nb * 20, sd = sqrt(es[["col"]] * var_env))
  parts$spatial <- e_bin[cw$bin] + e_row[(cw$bin - 1) * 16 + cw$row] + e_col[(cw$bin - 1) * 20 + cw$col]
  parts$residual <- stats::rnorm(n, sd = sqrt((1 - sum(es)) * var_env))
  y <- Reduce(`+`, parts)
  g_a_fam <- tapply(g_a, fam, mean)
  vv <- function(x) if (length(x) > 1) stats::var(x) else 0
  truth <- list(
    par = par, rho_shared = rho, attempts = attempt,
    sets = list(shared = sh, dir_only = d_only, ind_only = i_only, major_dir = maj_d, major_ige = maj_i,
                sim = sim_s, pair_intra = pair_intra, pair_inter = c(a = pair_inter_a, b = pair_inter_b),
                major_inter = maj_a),
    beta = list(major_dir = b_maj_d, major_ige = b_maj_i, major_inter = b_maj_a,
                gxg_intra = gi$beta, gxg_inter = ge$beta),
    g = list(dge = g_d, ige = g_i, inter = g_a, inter_fam = g_a_fam),
    gxg = list(intra = gi$G, inter = ge$G, inter_fam = G_inter_fam),
    realised = c(var_dge = vv(g_d), var_ige = vv(g_i), r_g = stats::cor(g_d, g_i),
                 var_inter_plant = vv(g_a), var_inter_fam = vv(g_a_fam),
                 var_gxg_intra = vv(as.vector(gi$G)), var_gxg_inter_plant = vv(as.vector(ge$G)),
                 var_gxg_inter_fam = vv(as.vector(G_inter_fam)),
                 var_bin = vv(e_bin), var_row = vv(e_row), var_col = vv(e_col),
                 var_res = vv(parts$residual), var_y = vv(y)),
    var_parts = vapply(parts, vv, 0))
  list(y = y, truth = truth)
}
