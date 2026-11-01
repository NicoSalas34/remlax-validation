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
