# ==============================================================================
# Dimensions de chaque cas du banc, reconstruites a partir de la graine (les
# donnees sont deterministes) : nombre de genotypes (ou de niveaux du terme
# genetique principal), observations par niveau et taille m des equations du
# modele mixte. n croit en ajoutant des genotypes, a nombre d'observations par
# genotype a peu pres constant : l'analyse doit le montrer pour qu'un effet de
# n ne soit pas lu comme un effet de repetition.
#
#   Rscript R/dimensions.R --sortie=<dossier des CSV du banc> [--pkg=<remlax>] [--coeurs=6]
#
# Sortie : <sortie>/dimensions.csv, une ligne par (famille, n, repetition).
# ==============================================================================
args <- commandArgs(TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", k, "="), "", a[1]) else d }
sortie <- opt("sortie", "resultats"); pkg <- opt("pkg", ""); coeurs <- as.integer(opt("coeurs", "6"))
ici <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
if (nzchar(pkg)) suppressMessages(pkgload::load_all(pkg, quiet = TRUE)) else suppressMessages(library(remlax))
source(file.path(ici, "..", "verif_simulations", "R", "outils.R")); source(file.path(ici, "R", "familles.R"))

# terme genetique principal : le premier terme du modele ; pour exp, ce sont
# des positions sur le transect et non des genotypes
NIVEAU <- c(iid = "genotypes", grm = "genotypes", us = "genotypes", ped = "animaux", ige = "genotypes",
            champ = "genotypes", resid = "genotypes", exp = "positions", fa = "genotypes")
TAILLES <- c(500, 1000, 1500, 2000, 3000, 4000, 6000, 8000)
grille_cas <- expand.grid(rep = 1:5, n = TAILLES, famille = names(CONSTRUCTEURS), stringsAsFactors = FALSE)

une <- function(i) {
  f <- grille_cas$famille[i]; n <- grille_cas$n[i]; r <- grille_cas$rep[i]
  cas <- CONSTRUCTEURS[[f]](n, graine_de(f, n, r))
  m <- cas$modeles$mme
  Zs <- lapply(m$terms, function(t) if (is.list(t$Z)) t$Z else list(t$Z))
  z1 <- Zs[[1]][[1]]
  nt <- if (f %in% c("us", "fa")) length(Zs[[1]]) else 1L       # caracteres
  niv <- ncol(z1)
  obs_niv <- Matrix::colSums(z1 != 0)                            # observations de chaque niveau (1er caractere)
  data.frame(famille = f, n = n, rep = r, n_obs = length(cas$y), n_caracteres = nt,
             niveau = NIVEAU[[f]], n_niveaux = niv, niveaux_observes = sum(obs_niv > 0),
             obs_par_niveau = length(cas$y) / (niv * nt), obs_par_niveau_min = min(obs_niv),
             obs_par_niveau_max = max(obs_niv), p = ncol(cas$X),
             m = ncol(cas$X) + sum(vapply(unlist(Zs, recursive = FALSE), ncol, 0L)),
             stringsAsFactors = FALSE)
}
lignes <- parallel::mclapply(seq_len(nrow(grille_cas)), function(i) tryCatch(une(i), error = function(e) {
  message("echec ", paste(grille_cas[i, ], collapse = " "), " : ", conditionMessage(e)); NULL }), mc.cores = coeurs)
d <- do.call(rbind, lignes)
utils::write.csv(d, file.path(sortie, "dimensions.csv"), row.names = FALSE)
print(stats::aggregate(cbind(n_niveaux, obs_par_niveau, m) ~ famille + n, d, stats::median))
