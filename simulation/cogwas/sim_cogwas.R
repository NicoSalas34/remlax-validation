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
#   in even ones; lines and families randomised; neighbour weights by
#   Chebyshev distance within a bin, averaged over the neighbours of each kind.
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
                       rank_inter = 1L, seed = 1L) {
  set.seed(seed)
  cells <- expand.grid(row = seq_len(n_row), col = seq_len(n_col), bin = seq_len(n_bins))
  cells$species <- ifelse(cells$col %% 2L == 1L, "wheat", "alfalfa")
  wi <- which(cells$species == "wheat"); ai <- which(cells$species == "alfalfa")
  lines <- rownames(gen$wheat$G); n_alf_plants <- nrow(gen$alfalfa$G_plant)
  cells$line <- NA_character_; cells$line[wi] <- sample(rep_len(sample(lines), length(wi)))
  cells$plant <- NA_integer_
  cells$plant[ai] <- sample(n_alf_plants, length(ai), replace = length(ai) > n_alf_plants)
  cells$family <- NA_character_; cells$family[ai] <- gen$alfalfa$fam_of_plant[cells$plant[ai]]
  W <- function(from, to, rank) .sim_weights(cells, from, to, rank)
  list(cells = cells, wheat = wi, alfalfa = ai,
       W_intra = W(wi, wi, rank_intra),       # wheat neighbours of each wheat plant
       W_inter = W(wi, ai, rank_inter))       # alfalfa neighbours of each wheat plant
}

# Mean over the neighbours of kind `to` within Chebyshev distance `rank` (same
# bin), for each plant of `from`; rows of plants without neighbour stay zero.
.sim_weights <- function(cells, from, to, rank) {
  i <- integer(0); j <- integer(0)
  for (b in unique(cells$bin)) {
    f <- from[cells$bin[from] == b]; t <- to[cells$bin[to] == b]
    d <- pmax(abs(outer(cells$row[f], cells$row[t], "-")), abs(outer(cells$col[f], cells$col[t], "-")))
    hit <- which(d > 0 & d <= rank, arr.ind = TRUE)
    i <- c(i, match(f[hit[, 1]], from)); j <- c(j, match(t[hit[, 2]], to))
  }
  W <- Matrix::sparseMatrix(i, j, x = 1, dims = c(length(from), length(to)))
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
# sim_pheno_full() -- wheat phenotypes with every kind of effect
# ------------------------------------------------------------------------------
# Genotype-level variances (what the mixed model estimates):
#   g_dge   (lines)          polygenic + major direct QTL
#   g_ige   (lines)          polygenic + major indirect QTL, acting on wheat neighbours
#   g_inter (alfalfa plants) polygenic + major QTL, acting on wheat neighbours
# Partial pleiotropy: the polygenic causal SNPs are split into direct only,
# indirect only and shared. On shared SNPs the indirect effect is
# rho z_d + sqrt(1 - rho^2) z_i; rho is set by bisection so that the REALISED
# correlation between g_dge and g_ige (majors included, LD included) equals r_g.
# Non-additive effects, in phenotype units (variance of the term):
#   sim        focal dose at s x mean wheat-neighbour dose at the same s (Sato)
#   pair_intra focal dose at a x mean wheat-neighbour dose at b (a != b)
#   pair_inter focal dose at a x mean alfalfa-neighbour dose at b (plant genotype)
# Spatial (bin, bin:row, bin:col) and residual, also in phenotype units.
# All causal loci are drawn among gen$*$causal, hence off the chip.
# ==============================================================================
sim_par_default <- function() list(
  n_poly = 1000L, frac_shared = 0.7, frac_dir_only = 0.15,
  var_dge = 0.30, var_ige = 0.30, r_g = -0.4, var_inter = 0.30, n_poly_inter = 500L,
  n_major_dir = 3L, n_major_ige = 3L, n_major_inter = 2L, pve_major = 0.10,
  n_sim = 2L, var_sim = 0.02, var_pair_intra = 0.03, var_pair_inter = 0.03,
  n_poly_sim = 500L, var_poly_sim = 0.02, n_poly_pair_intra = 2000L, var_poly_pair_intra = 0.03,
  n_poly_pair_inter = 2000L, var_poly_pair_inter = 0.03,
  var_bin = 0.05, var_row = 0.03, var_col = 0.03, var_res = 0.50)

#' Scenarios of the campaign: S0 additive only (polygenic and major additive
#' QTL), S1 adds the major similarity and pair QTL, S2 adds the polygenic
#' interaction background (explicit pairs).
sim_scenario <- function(name) {
  p <- sim_par_default()
  if (name %in% c("S0", "S1")) { p$var_poly_sim <- 0; p$var_poly_pair_intra <- 0; p$var_poly_pair_inter <- 0 }
  if (name == "S0") { p$n_sim <- 0L; p$var_pair_intra <- 0; p$var_pair_inter <- 0 }
  if (!name %in% c("S0", "S1", "S2")) stop("unknown scenario ", name)
  p
}

sim_pheno_full <- function(gen, des, par = sim_par_default(), seed = 1L) {
  set.seed(seed)
  cells <- des$cells; cw <- cells[des$wheat, ]
  line_of <- match(cw$line, rownames(gen$wheat$G))
  plant_of <- cells$plant[des$alfalfa]
  Gw <- gen$wheat$G; Ga <- gen$alfalfa$G_plant
  Xw <- sweep(Gw, 2, colMeans(Gw), "-"); Xa <- sweep(Ga, 2, colMeans(Ga), "-")
  pool <- sample(gen$wheat$causal)
  take <- function(k) { if (k == 0) return(integer(0)); s <- pool[seq_len(k)]; pool <<- pool[-seq_len(k)]; unname(s) }
  maj_d <- take(par$n_major_dir); maj_i <- take(par$n_major_ige); sim_s <- take(par$n_sim)
  pair_intra <- take(2L); pair_inter_a <- take(1L)
  if (length(pool) < par$n_poly) stop("sim_pheno_full: not enough causal wheat sites.")
  poly <- pool[seq_len(par$n_poly)]
  n_sh <- round(par$frac_shared * par$n_poly); n_d <- round(par$frac_dir_only * par$n_poly)
  sh <- poly[seq_len(n_sh)]; d_only <- poly[n_sh + seq_len(n_d)]
  i_only <- poly[-seq_len(n_sh + n_d)]
  # majors: effect so that var(X_s beta) = pve x target variance, random sign
  major <- function(X, s, v) sample(c(-1, 1), length(s), TRUE) * sqrt(par$pve_major * v / apply(X[, s, drop = FALSE], 2, stats::var))
  b_maj_d <- major(Xw, maj_d, par$var_dge); b_maj_i <- major(Xw, maj_i, par$var_ige)
  g_maj_d <- drop(Xw[, maj_d, drop = FALSE] %*% b_maj_d)
  g_maj_i <- drop(Xw[, maj_i, drop = FALSE] %*% b_maj_i)
  scale_to <- function(g, v) g * sqrt(v / stats::var(g))
  # The reachable range of the realised r_g depends on the draws (LD between
  # sets, majors). Out of reach: new polygenic draws, up to 20 attempts.
  for (attempt in 1:20) {
    zd_sh <- rnorm(n_sh); zi_sh <- rnorm(n_sh)
    bd <- c(zd_sh, rnorm(length(d_only))); names(bd) <- c(sh, d_only)
    g_poly_d <- scale_to(drop(Xw[, c(sh, d_only)] %*% bd), par$var_dge * (1 - par$n_major_dir * par$pve_major))
    g_d <- g_poly_d + g_maj_d
    zi_only <- rnorm(length(i_only))
    g_i_of <- function(rho) {
      bi <- c(rho * zd_sh + sqrt(1 - rho^2) * zi_sh, zi_only)
      g_poly_i <- scale_to(drop(Xw[, c(sh, i_only)] %*% bi), par$var_ige * (1 - par$n_major_ige * par$pve_major))
      list(g = g_poly_i + g_maj_i, b = bi)
    }
    r_of <- function(rho) stats::cor(g_d, g_i_of(rho)$g)
    lo <- -0.999; hi <- 0.999; r_lo <- r_of(lo); r_hi <- r_of(hi)
    if (par$r_g >= min(r_lo, r_hi) && par$r_g <= max(r_lo, r_hi)) break
    if (attempt == 20)
      stop(sprintf("sim_pheno_full: r_g = %.2f not reachable with this pleiotropy (last range %.2f to %.2f).",
                   par$r_g, min(r_lo, r_hi), max(r_lo, r_hi)))
  }
  for (it in 1:60) {                       # r_of is monotone in rho
    mid <- (lo + hi) / 2
    if ((r_of(mid) - par$r_g) * (r_lo - par$r_g) > 0) { lo <- mid; r_lo <- r_of(lo) } else hi <- mid
  }
  rho <- (lo + hi) / 2; gi <- g_i_of(rho); g_i <- gi$g
  # IGE inter: alfalfa plants
  pool_a <- sample(gen$alfalfa$causal)
  maj_a <- unname(pool_a[seq_len(par$n_major_inter)])
  poly_a <- unname(pool_a[par$n_major_inter + seq_len(par$n_poly_inter)])
  pair_inter_b <- unname(pool_a[par$n_major_inter + par$n_poly_inter + 1L])
  b_maj_a <- major(Xa, maj_a, par$var_inter)
  g_a <- scale_to(drop(Xa[, poly_a] %*% rnorm(length(poly_a))), par$var_inter * (1 - par$n_major_inter * par$pve_major)) +
    drop(Xa[, maj_a, drop = FALSE] %*% b_maj_a)
  # phenotype contributions on wheat plants
  Wi <- des$W_intra; Wa <- des$W_inter
  parts <- list(dge = g_d[line_of], ige = as.numeric(Wi %*% g_i[line_of]),
                inter = as.numeric(Wa %*% g_a[plant_of]))
  inter_term <- function(x_focal, x_neigh, W, v) {
    if (v <= 0) return(list(u = 0, beta = 0))
    u <- x_focal * as.numeric(W %*% x_neigh); s <- sample(c(-1, 1), 1)
    list(u = s * u * sqrt(v / stats::var(u)), beta = s * sqrt(v / stats::var(u)))
  }
  sims <- lapply(sim_s, function(s) inter_term(Xw[line_of, s], Xw[line_of, s], Wi, par$var_sim))
  parts$sim <- if (length(sims)) Reduce(`+`, lapply(sims, `[[`, "u")) else 0
  pi_ <- inter_term(Xw[line_of, pair_intra[1]], Xw[line_of, pair_intra[2]], Wi, par$var_pair_intra)
  pe_ <- inter_term(Xw[line_of, pair_inter_a], Xa[plant_of, pair_inter_b], Wa, par$var_pair_inter)
  parts$pair_intra <- pi_$u; parts$pair_inter <- pe_$u
  # polygenic interaction background: explicit pairs of causal loci, small
  # normal effects, the sum scaled to its variance
  poly_pairs <- function(F, N, v) {          # F, N: plants x pairs (focal, neighbour mean)
    if (v <= 0) return(0)
    scale_to(drop((F * N) %*% rnorm(ncol(F))), v)
  }
  cwz <- gen$wheat$causal; caz <- gen$alfalfa$causal
  ps <- if (par$var_poly_sim > 0) sample(cwz, par$n_poly_sim) else integer(0)
  parts$poly_sim <- poly_pairs(Xw[line_of, ps, drop = FALSE], as.matrix(Wi %*% Xw[line_of, ps, drop = FALSE]),
                               par$var_poly_sim)
  pia <- sample(cwz, par$n_poly_pair_intra, TRUE); pib <- sample(cwz, par$n_poly_pair_intra, TRUE)
  parts$poly_pair_intra <- poly_pairs(Xw[line_of, pia, drop = FALSE], as.matrix(Wi %*% Xw[line_of, pib, drop = FALSE]),
                                      par$var_poly_pair_intra)
  pea <- sample(cwz, par$n_poly_pair_inter, TRUE); peb <- sample(caz, par$n_poly_pair_inter, TRUE)
  parts$poly_pair_inter <- poly_pairs(Xw[line_of, pea, drop = FALSE], as.matrix(Wa %*% Xa[plant_of, peb, drop = FALSE]),
                                      par$var_poly_pair_inter)
  nb <- max(cw$bin)
  e_bin <- rnorm(nb, sd = sqrt(par$var_bin)); e_row <- rnorm(nb * 16, sd = sqrt(par$var_row))
  e_col <- rnorm(nb * 20, sd = sqrt(par$var_col))
  parts$spatial <- e_bin[cw$bin] + e_row[(cw$bin - 1) * 16 + cw$row] + e_col[(cw$bin - 1) * 20 + cw$col]
  parts$residual <- rnorm(nrow(cw), sd = sqrt(par$var_res))
  y <- Reduce(`+`, parts)
  if (length(y) == 1L) stop("sim_pheno_full: empty phenotype.")
  fam <- gen$alfalfa$fam_of_plant
  g_a_fam <- tapply(g_a, fam, mean)
  truth <- list(
    par = par, rho_shared = rho, attempts = attempt,
    sets = list(shared = sh, dir_only = d_only, ind_only = i_only, major_dir = maj_d, major_ige = maj_i,
                sim = sim_s, pair_intra = pair_intra, pair_inter = c(a = pair_inter_a, b = pair_inter_b),
                major_inter = maj_a),
    beta = list(major_dir = b_maj_d, major_ige = b_maj_i, major_inter = b_maj_a,
                sim = vapply(sims, `[[`, 0, "beta"), pair_intra = pi_$beta, pair_inter = pe_$beta),
    g = list(dge = g_d, ige = g_i, inter = g_a, inter_fam = g_a_fam),
    realised = c(var_dge = stats::var(g_d), var_ige = stats::var(g_i), r_g = stats::cor(g_d, g_i),
                 var_inter_plant = stats::var(g_a), var_inter_fam = stats::var(g_a_fam),
                 var_bin = stats::var(e_bin), var_row = stats::var(e_row), var_col = stats::var(e_col),
                 var_res = stats::var(parts$residual), var_y = stats::var(y)),
    share_y = vapply(parts, function(x) if (length(x) > 1) stats::var(x) else 0, 0) / stats::var(y),
    var_parts = vapply(parts, function(x) if (length(x) > 1) stats::var(x) else 0, 0))
  list(y = y, truth = truth)
}
