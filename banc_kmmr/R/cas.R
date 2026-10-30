# ==============================================================================
# Cas du banc remlax contre kmmr : memes donnees, memes modeles.
#
# Les jeux de l'article de kmmr (Melsen et al. 2026, arXiv:2610.16671) ne sont
# pas tous publics. Les cas reprennent les donnees embarquees par le paquet
# kmmr (DROPS reduit) et le jeu public de Krause et al. (2018) present dans son
# dossier dev/.
#
#   drops2_fa1, drops2_fa3 : deux etapes, 28 combinaisons gestion x environnement,
#                            50 hybrides, poids, residuelle fixee a 1, fa(k) x K
#   drops1_fa1, drops1_us  : une etape, 5 essais, 245 genotypes ; par essai :
#                            repetition, bloc, ar1(rang) x ar1(colonne), residuelle
#                            propre ; G x E en fa1 ou us(5) x K
#   krause_gblup           : un traitement, tous les cycles, ~3 768 genotypes, K dense
#   krause_fa2, krause_us  : 5 traitements x 3 768 genotypes, 18 840 observations,
#                            G x traitement fa2 ou us(5) x K, traitement x cycle fixe
# ==============================================================================

cas_noms <- c("drops2_fa1", "drops2_fa3", "drops1_fa1", "drops1_us",
              "krause_gblup", "krause_fa2", "krause_us")

charger_cas <- function(cas, dossier_kmmr, dossier_sim = NULL) {
  if (startsWith(cas, "sim_")) return(charger_sim(cas, dossier_sim))
  if (startsWith(cas, "drops2")) {
    e <- new.env()
    load(file.path(dossier_kmmr, "data", "dataDROPS.rda"), e)
    load(file.path(dossier_kmmr, "data", "kinshipDROPS.rda"), e)
    d <- droplevels(e$dataDROPS[!is.na(e$dataDROPS$GY), ])
    K <- as.matrix(e$kinshipDROPS)
    d$Gen <- factor(d$Gen, levels = rownames(K))
    return(list(cas = cas, d = d, K = K, type = "deux_etapes",
                env = "ManEnv", rang = as.integer(sub("drops2_fa", "", cas))))
  }
  if (startsWith(cas, "drops1")) {
    d <- readRDS(file.path(dossier_kmmr, "checks", "dataPlotsDROPS.rds"))
    K <- as.matrix(readRDS(file.path(dossier_kmmr, "checks", "kinshipPlotsDROPS.rds")))
    d <- droplevels(d[!is.na(d$GY), ])
    d$Gen <- factor(as.character(d$Gen), levels = rownames(K))
    d$Rep <- factor(paste(d$Env, d$Rep, sep = ":"))
    d$Block <- factor(paste(d$Rep, d$Block, sep = ":"))
    d$Row <- as.integer(as.character(d$Row)); d$Col <- as.integer(as.character(d$Col))
    return(list(cas = cas, d = d, K = K, type = "une_etape",
                ge = if (cas == "drops1_us") "us" else "fa1"))
  }
  # Krause et al. 2018
  b <- utils::read.csv(file.path(dossier_kmmr, "dev", "Krause_et_al_2018_Yield_BLUEs.csv"))
  K <- readRDS(file.path(dossier_kmmr, "dev", "Krause_K.rds"))
  b <- b[as.character(b$GID) %in% rownames(K) & !is.na(b$Grain_Yield_BLUE), ]
  names(b)[names(b) == "Grain_Yield_BLUE"] <- "GY"
  b$Trt <- factor(b$Managed_Treatment)
  b$Cycle <- factor(b$Breeding.Cycle)
  b$TrtCycle <- factor(paste(b$Trt, b$Cycle, sep = ":"))
  if (cas == "krause_gblup") {
    b <- b[b$Trt == "Optimal Flat", ]
    g <- unique(as.character(b$GID))
    K <- K[g, g]
    b$Gen <- factor(as.character(b$GID), levels = g)
    return(list(cas = cas, d = droplevels(b), K = K, type = "gblup"))
  }
  g <- unique(as.character(b$GID))
  K <- K[g, g]
  b$Gen <- factor(as.character(b$GID), levels = g)
  list(cas = cas, d = droplevels(b), K = K, type = "multi_trt",
       ge = if (cas == "krause_us") "us" else "fa2")
}
