# ==============================================================================
# Cases of the remlax against kmmr benchmark: same data, same models.
#
# The data sets of the kmmr paper (Melsen et al. 2026, arXiv:2610.16671) are
# not all public. The cases use the data shipped with the kmmr package (reduced
# DROPS) and the public data set of Krause et al. (2018) found in its dev/
# folder.
#
#   drops2_fa1, drops2_fa3 : two stages, 28 management x environment
#                            combinations, 50 hybrids, weights, residual fixed
#                            at 1, fa(k) x K
#   drops1_fa1, drops1_us  : one stage, 5 trials, 245 genotypes. Per trial:
#                            replicate, block, ar1(row) x ar1(column), own
#                            residual. G x E in fa1 or us(5) x K
#   krause_gblup           : one treatment, all cycles, ~3 768 genotypes, dense K
#   krause_fa2, krause_us  : 5 treatments x 3 768 genotypes, 18 840 observations,
#                            G x treatment fa2 or us(5) x K, fixed treatment x cycle
# ==============================================================================

case_names <- c("drops2_fa1", "drops2_fa3", "drops1_fa1", "drops1_us",
                "krause_gblup", "krause_fa2", "krause_us")

load_case <- function(case, kmmr_dir, sim_dir = NULL) {
  if (startsWith(case, "sim_")) return(load_sim(case, sim_dir))
  if (startsWith(case, "drops2")) {
    e <- new.env()
    load(file.path(kmmr_dir, "data", "dataDROPS.rda"), e)
    load(file.path(kmmr_dir, "data", "kinshipDROPS.rda"), e)
    d <- droplevels(e$dataDROPS[!is.na(e$dataDROPS$GY), ])
    K <- as.matrix(e$kinshipDROPS)
    d$Gen <- factor(d$Gen, levels = rownames(K))
    return(list(case = case, d = d, K = K, type = "two_stage",
                env = "ManEnv", rank = as.integer(sub("drops2_fa", "", case))))
  }
  if (startsWith(case, "drops1")) {
    d <- readRDS(file.path(kmmr_dir, "checks", "dataPlotsDROPS.rds"))
    K <- as.matrix(readRDS(file.path(kmmr_dir, "checks", "kinshipPlotsDROPS.rds")))
    d <- droplevels(d[!is.na(d$GY), ])
    d$Gen <- factor(as.character(d$Gen), levels = rownames(K))
    d$Rep <- factor(paste(d$Env, d$Rep, sep = ":"))
    d$Block <- factor(paste(d$Rep, d$Block, sep = ":"))
    d$Row <- as.integer(as.character(d$Row)); d$Col <- as.integer(as.character(d$Col))
    return(list(case = case, d = d, K = K, type = "one_stage",
                ge = if (case == "drops1_us") "us" else "fa1"))
  }
  # Krause et al. 2018
  b <- utils::read.csv(file.path(kmmr_dir, "dev", "Krause_et_al_2018_Yield_BLUEs.csv"))
  K <- readRDS(file.path(kmmr_dir, "dev", "Krause_K.rds"))
  b <- b[as.character(b$GID) %in% rownames(K) & !is.na(b$Grain_Yield_BLUE), ]
  names(b)[names(b) == "Grain_Yield_BLUE"] <- "GY"
  b$Trt <- factor(b$Managed_Treatment)
  b$Cycle <- factor(b$Breeding.Cycle)
  b$TrtCycle <- factor(paste(b$Trt, b$Cycle, sep = ":"))
  if (case == "krause_gblup") {
    b <- b[b$Trt == "Optimal Flat", ]
    g <- unique(as.character(b$GID))
    K <- K[g, g]
    b$Gen <- factor(as.character(b$GID), levels = g)
    return(list(case = case, d = droplevels(b), K = K, type = "gblup"))
  }
  g <- unique(as.character(b$GID))
  K <- K[g, g]
  b$Gen <- factor(as.character(b$GID), levels = g)
  list(case = case, d = droplevels(b), K = K, type = "multi_trt",
       ge = if (case == "krause_us") "us" else "fa2")
}
