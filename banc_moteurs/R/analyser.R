# ==============================================================================
# Analyse du banc des moteurs : rassemble les CSV de executer.R, puis trace les
# figures de comparaison (temps, CPU contre GPU, MME contre dense, inversion,
# justesse de la logLik, recuperation des parametres, accord entre logiciels).
#
#   Rscript R/analyser.R --res=<dossier des CSV> [--sortie=<dossier>/figures]
#
# Une cellule = (famille, n, repetition). La meme graine sert a toutes les voies
# d'une cellule : les rapports de temps et les ecarts de logLik se calculent
# DANS la cellule, puis on resume les cinq repetitions (mediane, min, max).
# La logLik de reference d'une cellule est la meilleure logLik ARBITRE (REML
# dense commun, evalue aux estimations de chaque logiciel) parmi les
# ajustements reussis : les logLik natives n'ont pas la meme constante.
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
args <- commandArgs(TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", k, "="), "", a[1]) else d }
res <- opt("res", "resultats"); sortie <- opt("sortie", file.path(res, "figures"))
fdim <- opt("dimensions", file.path(res, "dimensions.csv"))      # produit par R/dimensions.R
dir.create(sortie, recursive = TRUE, showWarnings = FALSE)

lire <- function(motif) {
  f <- list.files(res, motif, full.names = TRUE)
  if (!length(f)) stop("aucun fichier ", motif, " dans ", res)
  l <- lapply(f, utils::read.csv, stringsAsFactors = FALSE)
  # colonnes ajoutees au fil des versions du banc (n_par, grad_rel) : on complete
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { for (k in setdiff(cols, names(d))) d[[k]] <- NA; d[cols] }))
}
aj <- lire("^ajustements_.*\\.csv$")
pa <- lire("^parametres_.*\\.csv$")

FAMILLES <- c("iid", "grm", "us", "ped", "ige", "champ", "resid", "exp", "fa")
VOIES <- c("D-CPU", "D-GPU", "M-dense-CPU", "M-dense-GPU", "M-creux-CPU", "asreml", "sommer", "lme4")
COUL <- c("D-CPU" = "#1f78b4", "D-GPU" = "#a6cee3", "M-dense-CPU" = "#33a02c", "M-dense-GPU" = "#b2df8a",
          "M-creux-CPU" = "#e31a1c", asreml = "#ff7f00", sommer = "#6a3d9a", lme4 = "#b15928")
aj$famille <- factor(aj$famille, FAMILLES); pa$famille <- factor(pa$famille, FAMILLES)
aj$voie <- factor(aj$voie, VOIES); pa$voie <- factor(pa$voie, VOIES)
aj$cellule <- paste(aj$famille, aj$n, aj$rep)
pa$cellule <- paste(pa$famille, pa$n, pa$rep)

# DIMENSIONS. n croit en ajoutant des genotypes (nombre d'observations par
# genotype a peu pres fixe), pas en ajoutant des repetitions : m, la taille de
# C, croit avec n. Chaque panneau porte la regle de sa famille, pour qu'un effet
# de n ne soit pas lu comme un effet du nombre de repetitions par genotype.
dims <- if (file.exists(fdim)) utils::read.csv(fdim, stringsAsFactors = FALSE) else NULL
if (!is.null(dims)) {
  dims$famille <- factor(dims$famille, FAMILLES)
  regle <- vapply(FAMILLES, function(f) {
    d <- dims[dims$famille == f, ]
    if (!nrow(d)) return(f)
    k <- stats::median(d$n_obs / d$n_niveaux)
    k <- if (abs(k - round(k)) < 0.15) round(k) else round(k, 1)
    nt <- d$n_caracteres[1]
    sprintf("%s : %s = %s%s", f, d$niveau[1], if (k == 1) "n" else paste0("n/", k),
            if (nt > 1) sprintf(" (%d caracteres)", nt) else "")
  }, "")
  lab_fam <- as_labeller(regle)
  aj <- merge(aj, dims[, c("famille", "n", "rep", "n_niveaux", "m")], by = c("famille", "n", "rep"), all.x = TRUE)
} else {
  message("dimensions.csv absent : lancer R/dimensions.R ; panneaux sans regle d'echelle")
  lab_fam <- label_value
}

# echecs de memoire GPU ecrits ERREUR par les premieres taches du tableau
aj$statut[aj$statut == "ERREUR" & grepl("bfc_allocator|RESOURCE_EXHAUSTED|out of memory", aj$message)] <- "MEMOIRE"
# un ajustement compte s'il a abouti et si sa logLik arbitre est connue
aj$reussi <- aj$statut == "OK" & is.finite(aj$logLik_arbitre)
ref <- tapply(ifelse(aj$reussi, aj$logLik_arbitre, NA), aj$cellule, max, na.rm = TRUE)
ref[!is.finite(ref)] <- NA
aj$ecart_logLik <- ref[aj$cellule] - aj$logLik_arbitre          # >= 0, en unites de logLik
# tolerance : sous 1e-3 l'ecart est numerique ; au-dela l'ajustement est sous l'optimum
TOL <- 1e-3
aj$juste <- aj$reussi & aj$ecart_logLik <= TOL
# temps et parametres : un ajustement a moins de 0,1 unite de logLik est retenu
# (asreml s'arrete parfois a quelques 1e-3 de l'optimum) ; au-dela il est ecarte
TOL_USAGE <- 0.1
aj$retenu <- aj$reussi & aj$ecart_logLik <= TOL_USAGE

resumer <- function(d, val, par = c("famille", "n", "voie")) {
  s <- stats::aggregate(d[[val]], d[par], function(x) c(med = stats::median(x), lo = min(x), hi = max(x), k = length(x)))
  cbind(s[par], as.data.frame(s$x))
}
theme_banc <- theme_bw(base_size = 10) + theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"))
enregistrer <- function(p, nom, l = 11, h = 8) {
  ggsave(file.path(sortie, paste0(nom, ".pdf")), p, width = l, height = h)
  ggsave(file.path(sortie, paste0(nom, ".png")), p, width = l, height = h, dpi = 150)
}

# ---------------------------------------------------------------- bilan d'etat
etat <- as.data.frame(table(famille = aj$famille, voie = aj$voie, statut = ifelse(aj$juste, "juste",
                      ifelse(aj$reussi, "sous l'optimum", aj$statut))))
etat <- etat[etat$Freq > 0, ]
utils::write.csv(etat, file.path(sortie, "bilan_etats.csv"), row.names = FALSE)
p <- ggplot(etat, aes(voie, Freq, fill = statut)) + geom_col() + facet_wrap(~ famille, labeller = lab_fam, nrow = 3) +
  scale_fill_manual(values = c(juste = "#4daf4a", "sous l'optimum" = "#ff7f00", ERREUR = "#e41a1c", DELAI = "#984ea3", MEMOIRE = "#377eb8")) +
  labs(x = NULL, y = "ajustements", fill = NULL, title = "Etat des ajustements par famille et par voie",
       subtitle = sprintf("juste : logLik arbitre a moins de %g de la meilleure de la cellule", TOL)) +
  theme_banc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
enregistrer(p, "fig0_etats")

# ------------------------------------------------------- dimensions contre n
if (!is.null(dims)) {
  dl <- rbind(data.frame(dims[c("famille", "n")], mesure = "genotypes (ou niveaux) du terme principal", v = dims$n_niveaux),
              data.frame(dims[c("famille", "n")], mesure = "observations par niveau et par caractere", v = dims$obs_par_niveau),
              data.frame(dims[c("famille", "n")], mesure = "m (taille de C)", v = dims$m))
  dl <- stats::aggregate(v ~ famille + n + mesure, dl, stats::median)
  p <- ggplot(dl, aes(n, v, colour = famille)) + geom_line() + geom_point(size = 1) +
    scale_x_log10() + scale_y_log10() + facet_wrap(~ mesure, scales = "free_y", nrow = 1) +
    labs(x = "n (observations)", y = NULL, colour = NULL, title = "Dimensions des cas simules",
         subtitle = paste("n croit en ajoutant des genotypes ; le nombre d'observations par genotype reste a peu pres fixe.",
                          "Mediane des 5 repetitions.")) +
    theme_banc
  enregistrer(p, "figD_dimensions", 12, 4.5)
  dt <- stats::aggregate(cbind(n_niveaux, obs_par_niveau, obs_par_niveau_min, obs_par_niveau_max, m) ~ famille + niveau + n,
                         dims, stats::median)
  utils::write.csv(dt[order(dt$famille, dt$n), ], file.path(sortie, "dimensions_resume.csv"), row.names = FALSE)
}

# ------------------------------------------------------------- temps contre n
tj <- aj[aj$retenu, ]
st <- resumer(tj, "secondes")
if (!is.null(dims)) st <- merge(st, stats::aggregate(cbind(n_niveaux, m) ~ famille + n, aj, stats::median),
                                by = c("famille", "n"), all.x = TRUE)
p <- ggplot(st, aes(n, med, colour = voie, fill = voie)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point(size = 1) +
  scale_x_log10() + scale_y_log10() + facet_wrap(~ famille, labeller = lab_fam, nrow = 3, scales = "free_y") +
  scale_colour_manual(values = COUL) + scale_fill_manual(values = COUL) +
  labs(x = "n (observations)", y = "temps total (s)", colour = NULL, fill = NULL,
       title = "Temps d'ajustement", subtitle = sprintf("mediane de 5 repetitions, bande min-max ; ajustements a moins de %g de la meilleure logLik", TOL_USAGE)) +
  theme_banc
enregistrer(p, "fig1_temps")
utils::write.csv(st, file.path(sortie, "temps_resume.csv"), row.names = FALSE)
# remlax : temps du seul solveur Python (hors lancement du processus, import de
# JAX et export du paquet), a cote du temps total des logiciels externes
ts <- tj; rx <- ts$voie %in% VOIES[1:5]
ts$secondes[rx] <- ts$solveur_s[rx]
ts <- ts[is.finite(ts$secondes), ]
p <- ggplot(resumer(ts, "secondes"), aes(n, med, colour = voie, fill = voie)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point(size = 1) +
  scale_x_log10() + scale_y_log10() + facet_wrap(~ famille, labeller = lab_fam, nrow = 3, scales = "free_y") +
  scale_colour_manual(values = COUL) + scale_fill_manual(values = COUL) +
  labs(x = "n (observations)", y = "temps (s)", colour = NULL, fill = NULL,
       title = "Temps du solveur remlax seul",
       subtitle = "remlax : solveur Python, compilation comprise, sans lancement ni export ; externes : temps total") +
  theme_banc
enregistrer(p, "fig1b_temps_solveur")
# meme temps contre m : le cout du moteur MME suit m, celui du dense suit n
if (!is.null(dims)) {
  # m varie un peu d'une repetition a l'autre (niveaux non observes) : m median par (famille, n)
  sm <- merge(resumer(tj, "secondes"), stats::aggregate(m ~ famille + n, dims, stats::median), by = c("famille", "n"))
  p <- ggplot(sm, aes(m, med, colour = voie, fill = voie)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point(size = 1) +
    scale_x_log10() + scale_y_log10() + facet_wrap(~ famille, labeller = lab_fam, nrow = 3, scales = "free") +
    scale_colour_manual(values = COUL) + scale_fill_manual(values = COUL) +
    labs(x = "m (taille de C : effets fixes + niveaux aleatoires)", y = "temps total (s)", colour = NULL, fill = NULL,
         title = "Temps d'ajustement contre la taille des equations du modele mixte") + theme_banc
  enregistrer(p, "fig1c_temps_contre_m")
}

# rapport de deux voies dans chaque cellule, puis resume des repetitions
rapport <- function(num, den) {
  a <- tj[tj$voie == num, c("famille", "n", "rep", "cellule", "secondes")]
  b <- tj[tj$voie == den, c("cellule", "secondes")]
  m <- merge(a, b, by = "cellule", suffixes = c("_num", "_den"))
  if (!nrow(m)) return(NULL)
  m$r <- m$secondes_num / m$secondes_den
  s <- stats::aggregate(r ~ famille + n, m, function(x) c(med = stats::median(x), lo = min(x), hi = max(x)))
  cbind(s[c("famille", "n")], as.data.frame(s$r), comparaison = sprintf("%s / %s", num, den))
}
trace_rapports <- function(r, titre, sous, nom) {
  r <- do.call(rbind, r)
  if (is.null(r) || !nrow(r)) return(invisible())
  p <- ggplot(r, aes(n, med, colour = comparaison, fill = comparaison)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point(size = 1) +
    scale_x_log10() + scale_y_log10() + facet_wrap(~ famille, labeller = lab_fam, nrow = 3) +
    labs(x = "n (observations)", y = "rapport des temps", colour = NULL, fill = NULL, title = titre, subtitle = sous) +
    theme_banc
  enregistrer(p, nom)
  utils::write.csv(r, file.path(sortie, paste0(nom, ".csv")), row.names = FALSE)
}
trace_rapports(list(rapport("D-CPU", "D-GPU"), rapport("M-dense-CPU", "M-dense-GPU")),
               "Gain du GPU", "temps CPU / temps GPU, meme noeud ; au-dessus de 1 le GPU l'emporte", "fig2_gpu_contre_cpu")
trace_rapports(list(rapport("D-CPU", "M-dense-CPU"), rapport("D-CPU", "M-creux-CPU"), rapport("D-GPU", "M-dense-GPU")),
               "Moteur MME contre moteur dense", "temps dense / temps MME ; au-dessus de 1 le MME l'emporte", "fig3_mme_contre_dense")
trace_rapports(list(rapport("M-dense-CPU", "M-creux-CPU")),
               "Inversion de C dans le moteur MME", "temps inversion dense / temps inversion creuse (CPU)", "fig4_inversion")
trace_rapports(list(rapport("asreml", "D-CPU"), rapport("asreml", "M-creux-CPU"), rapport("asreml", "M-dense-GPU"),
                    rapport("asreml", "D-GPU")),
               "remlax contre asreml", "temps asreml / temps remlax ; au-dessus de 1 remlax l'emporte", "fig5_remlax_contre_asreml")

# ------------------------------------------------- meilleure voie par cellule
mv <- do.call(rbind, lapply(split(tj, list(tj$famille, tj$n), drop = TRUE), function(d) {
  s <- tapply(d$secondes, droplevels(d$voie), stats::median)
  data.frame(famille = d$famille[1], n = d$n[1], voie = names(s)[which.min(s)], secondes = min(s))
}))
mv$voie <- factor(mv$voie, VOIES)
p <- ggplot(mv, aes(factor(n), famille, fill = voie)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.3g", secondes)), size = 2.6) + scale_fill_manual(values = COUL, drop = FALSE) +
  labs(x = "n (observations)", y = NULL, fill = NULL, title = "Voie la plus rapide par cellule",
       subtitle = "mediane des repetitions retenues ; etiquette = temps (s)") + theme_banc
enregistrer(p, "fig6_meilleure_voie", 10, 5)
mvr <- do.call(rbind, lapply(split(tj[tj$voie %in% VOIES[1:5], ], list(tj$famille[tj$voie %in% VOIES[1:5]], tj$n[tj$voie %in% VOIES[1:5]]), drop = TRUE), function(d) {
  s <- tapply(d$secondes, droplevels(d$voie), stats::median)
  data.frame(famille = d$famille[1], n = d$n[1], voie = names(s)[which.min(s)], secondes = min(s))
}))
mvr$voie <- factor(mvr$voie, VOIES)
p <- ggplot(mvr, aes(factor(n), famille, fill = voie)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.3g", secondes)), size = 2.6) + scale_fill_manual(values = COUL[1:5], drop = FALSE) +
  labs(x = "n (observations)", y = NULL, fill = NULL, title = "Voie remlax la plus rapide par cellule",
       subtitle = "base de la regle engine = \"auto\"") + theme_banc
enregistrer(p, "fig7_meilleure_voie_remlax", 10, 5)

# ----------------------------------------------------------- justesse logLik
ok <- aj[aj$reussi, ]
p <- ggplot(ok, aes(factor(n), pmax(ecart_logLik, 1e-10), colour = voie)) +
  geom_hline(yintercept = TOL, linetype = 2) +
  geom_point(position = position_dodge(width = 0.7), size = 0.9, alpha = 0.8) +
  scale_y_log10() + facet_wrap(~ famille, labeller = lab_fam, nrow = 3) + scale_colour_manual(values = COUL) +
  labs(x = "n (observations)", y = "ecart a la meilleure logLik (arbitre)", colour = NULL,
       title = "Justesse de l'optimum", subtitle = sprintf("ecarts plafonnes en bas a 1e-10 ; tirets : tolerance %g", TOL)) +
  theme_banc
enregistrer(p, "fig8_ecart_logLik")

# --------------------------------------------- recuperation des parametres
# remlax : les voies donnent le meme optimum, on garde D-CPU (ou la premiere juste)
pj <- merge(pa, aj[, c("cellule", "voie", "retenu")], by = c("cellule", "voie"))
pj <- pj[pj$retenu, ]
pj$logiciel <- ifelse(pj$voie %in% VOIES[1:5], "remlax", as.character(pj$voie))
pr <- pj[!duplicated(pj[c("cellule", "logiciel", "parametre")]), ]
pr$erreur <- pr$estime - pr$vrai
rec <- stats::aggregate(erreur ~ famille + n + logiciel + parametre + vrai, pr,
                        function(x) c(biais = mean(x), rmse = sqrt(mean(x^2)), k = length(x)))
rec <- cbind(rec[c("famille", "n", "logiciel", "parametre", "vrai")], as.data.frame(rec$erreur))
utils::write.csv(rec, file.path(sortie, "recuperation_parametres.csv"), row.names = FALSE)
# une valeur par famille : RMSE relative mediane sur les parametres de variance
rec$rrmse <- rec$rmse / pmax(abs(rec$vrai), 0.05)
rf <- stats::aggregate(rrmse ~ famille + n + logiciel, rec, stats::median)
COUL_L <- c(remlax = "#1f78b4", asreml = "#ff7f00", sommer = "#6a3d9a", lme4 = "#b15928")
p <- ggplot(rf, aes(n, rrmse, colour = logiciel)) + geom_line() + geom_point(size = 1) +
  scale_x_log10() + scale_y_log10() + facet_wrap(~ famille, labeller = lab_fam, nrow = 3) + scale_colour_manual(values = COUL_L) +
  labs(x = "n (observations)", y = "RMSE relative (mediane sur les parametres)", colour = NULL,
       title = "Recuperation des parametres simules",
       subtitle = "RMSE sur 5 repetitions divisee par max(|vrai|, 0.05) ; les courbes superposees signalent le meme estimateur") +
  theme_banc
enregistrer(p, "fig9_recuperation")
p <- ggplot(pr, aes(factor(n), erreur, fill = logiciel)) + geom_hline(yintercept = 0) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.25) + facet_wrap(~ famille + parametre, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = COUL_L) +
  labs(x = "n", y = "estime - vrai", fill = NULL, title = "Erreur d'estimation par parametre") +
  theme_banc + theme(strip.text = element_text(size = 6), axis.text = element_text(size = 5))
enregistrer(p, "fig10_erreur_par_parametre", 16, 22)

# -------------------------------------------------- accord entre logiciels
rx <- pr[pr$logiciel == "remlax", c("cellule", "parametre", "estime")]
ac <- merge(pr[pr$logiciel != "remlax", ], rx, by = c("cellule", "parametre"), suffixes = c("", "_remlax"))
if (nrow(ac)) {
  ac$diff <- abs(ac$estime - ac$estime_remlax) / pmax(abs(ac$estime_remlax), 0.05)
  ad <- stats::aggregate(diff ~ famille + n + rep + logiciel, ac, max)
  p <- ggplot(ad, aes(factor(n), pmax(diff, 1e-8), colour = logiciel)) +
    geom_point(position = position_dodge(width = 0.6), size = 1) + scale_y_log10() +
    facet_wrap(~ famille, labeller = lab_fam, nrow = 3) + scale_colour_manual(values = COUL_L) +
    labs(x = "n (observations)", y = "ecart relatif maximal a remlax", colour = NULL,
         title = "Accord des estimations avec remlax",
         subtitle = "par cellule, max sur les parametres de |estime - remlax| / max(|remlax|, 0.05)") + theme_banc
  enregistrer(p, "fig11_accord_logiciels")
  utils::write.csv(ad, file.path(sortie, "accord_logiciels.csv"), row.names = FALSE)
}

# ------------------------------------------- compilation contre calcul (remlax)
rc <- tj[tj$voie %in% VOIES[1:5] & is.finite(tj$compile_s), ]
if (nrow(rc)) {
  rc$part_compilation <- rc$compile_s / rc$secondes
  s <- resumer(rc, "part_compilation")
  p <- ggplot(s, aes(n, med, colour = voie)) + geom_line() + geom_point(size = 1) + scale_x_log10() +
    facet_wrap(~ famille, labeller = lab_fam, nrow = 3) + scale_colour_manual(values = COUL) +
    labs(x = "n (observations)", y = "part de la compilation dans le temps total", colour = NULL,
         title = "Poids de la compilation JAX") + theme_banc
  enregistrer(p, "fig12_compilation")
}
cat(sprintf("%d ajustements, %d justes, %d retenus ; figures dans %s\n", nrow(aj), sum(aj$juste), sum(aj$retenu), sortie))
