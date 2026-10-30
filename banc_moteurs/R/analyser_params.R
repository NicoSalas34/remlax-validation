# ==============================================================================
# Analyse du banc "nombre de parametres" (banc_params.sbatch) : temps
# d'ajustement en fonction du nombre de parametres de variance, a n = 4000, pour
# chaque voie remlax et asreml.
#
#   Rscript R/analyser_params.R --res=<dossier res_params> [--sortie=<dossier>/figures]
# ==============================================================================
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", k, "="), "", a[1]) else d }
res <- opt("res", "res_params"); sortie <- opt("sortie", file.path(res, "figures"))
dir.create(sortie, recursive = TRUE, showWarnings = FALSE)
f <- list.files(res, "^ajustements_.*\\.csv$", full.names = TRUE, recursive = TRUE)
l <- lapply(f, utils::read.csv, stringsAsFactors = FALSE)
cols <- unique(unlist(lapply(l, names)))
aj <- do.call(rbind, lapply(l, function(d) { for (k in setdiff(cols, names(d))) d[[k]] <- NA; d[cols] }))

npar <- function(fam) {
  t <- as.integer(sub("^mt(us|igee|ige)", "", fam))
  base <- (2L * t) * (2L * t + 1L) / 2L + t * (t + 1L) / 2L
  ifelse(startsWith(fam, "mtus"), t * (t + 1L),
         ifelse(startsWith(fam, "mtigee"), base + t * (t + 1L) / 2L, base))
}
aj$famille_type <- ifelse(startsWith(aj$famille, "mtus"), "mtus : genotype us + residuelle us",
                   ifelse(startsWith(aj$famille, "mtigee"), "mtigee : DGE/IGE us(2t) + IEE plante us(t) + residuelle us",
                          "mtige : DGE/IGE us(2t) + residuelle us"))
aj$p <- npar(aj$famille)
aj$t <- as.integer(sub("^mt(us|igee|ige)", "", aj$famille))
aj$cellule <- paste(aj$famille, aj$n, aj$rep)
aj$reussi <- aj$statut == "OK" & is.finite(aj$logLik_arbitre)
ref <- tapply(ifelse(aj$reussi, aj$logLik_arbitre, NA), aj$cellule, max, na.rm = TRUE)
ref[!is.finite(ref)] <- NA
aj$ecart <- ref[aj$cellule] - aj$logLik_arbitre
aj$retenu <- aj$reussi & aj$ecart <= 0.1
VOIES <- c("D-CPU", "D-GPU", "D-CPU-LBFGS", "D-GPU-LBFGS", "M-dense-CPU", "M-dense-GPU", "M-creux-CPU", "asreml")
COUL <- c("D-CPU" = "#1f78b4", "D-GPU" = "#a6cee3", "D-CPU-LBFGS" = "#6a3d9a", "D-GPU-LBFGS" = "#cab2d6",
          "M-dense-CPU" = "#33a02c", "M-dense-GPU" = "#b2df8a", "M-creux-CPU" = "#e31a1c", asreml = "#ff7f00")
aj$voie <- factor(aj$voie, VOIES)
theme_set(theme_bw(base_size = 10) + theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95")))
sauver <- function(p, nom, l = 11, h = 5.5) {
  ggsave(file.path(sortie, paste0(nom, ".pdf")), p, width = l, height = h)
  ggsave(file.path(sortie, paste0(nom, ".png")), p, width = l, height = h, dpi = 150)
}

tj <- aj[aj$retenu, ]
# banc "nombre de parametres" : n = 1000 t (1000 unites quel que soit t) ; le
# banc "taille" (mtigee3, n variable) est trace a part
tj$banc_par <- tj$n == 1000L * tj$t
s <- stats::aggregate(secondes ~ famille_type + p + t + voie, tj[tj$banc_par, ],
                      function(x) c(med = stats::median(x), lo = min(x), hi = max(x)))
s <- cbind(s[1:4], as.data.frame(s$secondes))
p1 <- ggplot(s, aes(p, med, colour = voie, fill = voie)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point() +
  scale_x_log10() + scale_y_log10() + facet_wrap(~ famille_type) +
  scale_colour_manual(values = COUL, drop = FALSE) + scale_fill_manual(values = COUL, drop = FALSE) +
  labs(x = "nombre de parametres de variance", y = "temps total (s)", colour = NULL, fill = NULL,
       title = "Temps d'ajustement contre le nombre de parametres, 1000 unites, n = 1000 t (tranche MIG)",
       subtitle = "mediane de 3 repetitions, bande min-max ; ajustements a moins de 0,1 de la meilleure logLik")
sauver(p1, "figP1_temps_contre_parametres")
utils::write.csv(s, file.path(sortie, "temps_parametres.csv"), row.names = FALSE)

# iterations et evaluations : ce qui fait le cout par optimiseur
it <- stats::aggregate(cbind(n_iter, n_eval) ~ famille_type + p + voie,
                       aj[aj$retenu & aj$voie != "asreml" & aj$n == 1000L * aj$t, ], stats::median)
p2 <- ggplot(it, aes(p, n_eval, colour = voie)) + geom_line() + geom_point() + scale_x_log10() + scale_y_log10() +
  facet_wrap(~ famille_type) + scale_colour_manual(values = COUL, drop = FALSE) +
  labs(x = "nombre de parametres de variance", y = "evaluations de la vraisemblance", colour = NULL,
       title = "Evaluations par ajustement : L-BFGS-B contre AI")
sauver(p2, "figP2_evaluations")

if (any(tj$famille == "mtigee3")) {
  st <- stats::aggregate(secondes ~ famille + n + voie, tj[tj$famille == "mtigee3", ], function(x) c(med = stats::median(x), lo = min(x), hi = max(x)))
  st <- cbind(st[1:3], as.data.frame(st$secondes))
  p3 <- ggplot(st, aes(n, med, colour = voie, fill = voie)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point() +
    scale_x_log10() + scale_y_log10() + facet_wrap(~ famille) +
    scale_colour_manual(values = COUL, drop = FALSE) + scale_fill_manual(values = COUL, drop = FALSE) +
    labs(x = "n (observations)", y = "temps total (s)", colour = NULL, fill = NULL,
         title = "Temps d'ajustement contre n (tranche MIG)",
         subtitle = "mediane des repetitions, bande min-max ; ajustements a moins de 0,1 de la meilleure logLik")
  sauver(p3, "figP3_temps_contre_n", 8, 5.5)
  utils::write.csv(st, file.path(sortie, "temps_taille.csv"), row.names = FALSE)
}

etat <- as.data.frame(table(famille = aj$famille, voie = aj$voie,
                            statut = ifelse(aj$retenu, "retenu", ifelse(aj$reussi, "sous l'optimum", aj$statut))))
etat <- etat[etat$Freq > 0, ]
utils::write.csv(etat, file.path(sortie, "etats_parametres.csv"), row.names = FALSE)
cat(sprintf("%d ajustements, %d retenus ; figures dans %s\n", nrow(aj), sum(aj$retenu), sortie))
