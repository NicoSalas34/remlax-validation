# ==============================================================================
# Bilan de test_noyau_deux_especes.R (source par son mode "bilan", ou seul :
#   Rscript test_noyau_deux_especes_bilan.R <dossier>)
#
# Estimation de lambda et d : profil de la logLik sur la grille, en d (maximum
# sur lambda) et en lambda (maximum sur d), interpole par spline ; intervalle a
# 95 % = chute de qchisq(0,95, 1) / 2. Grille "r3" (vrai rayon) et "R6"
# (fenetre large). Precision et erreur de variance : contre la verite simulee,
# especes A et B empilees (meme parametres).
# ==============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
if (!exists("dossier")) dossier <- commandArgs(trailingOnly = TRUE)[1]
if (!exists("vrai")) vrai <- list(r = 3, lambda = 1.5, d = 0.5)
fs <- list.files(dossier, "^tache_.*\\.csv$", full.names = TRUE)
res <- do.call(rbind, lapply(fs, function(f) {
  x <- utils::read.csv(f, stringsAsFactors = FALSE); x$erreur <- as.character(x$erreur); x }))
sortie <- file.path(dossier, "bilan"); dir.create(sortie, showWarnings = FALSE)
journal <- file(file.path(sortie, "journal.txt"), "w")
dire <- function(...) { txt <- sprintf(...); cat(txt, "\n"); cat(txt, "\n", file = journal) }
dire("%d jeux de donnees, %d ajustements, %d echec(s), %d non converge(s)",
     length(fs), nrow(res), sum(is.na(res$logLik)), sum(res$conv %in% FALSE))
dire("Temps median par ajustement (s) selon n_rep : %s",
     paste(sprintf("%d : %.0f", sort(unique(res$n_rep)),
                   tapply(res$secondes, res$n_rep, stats::median, na.rm = TRUE)), collapse = ", "))

# --- format long : une ligne par (ajustement, espece, effet) ------------------
eff <- c(dge = "DGE", intra = "IGE intra", inter = "IGE inter", tbv = "TBV intra")
long <- do.call(rbind, lapply(c("A", "B"), function(sp) do.call(rbind, lapply(names(eff), function(e) {
  data.frame(res[, c("n_rep", "graine", "grille", "r", "lambda", "d", "logLik", "conv")], espece = sp,
             effet = eff[[e]], prec = res[[sprintf("prec_%s_%s", e, sp)]],
             err = if (e == "tbv") NA_real_ else res[[sprintf("err_%s_%s", e, sp)]])
}))))
long$effet <- factor(long$effet, levels = eff)

# --- 1. estimation de lambda et d par profil ---------------------------------
profil1 <- function(x, ll, borne) {
  o <- order(x); x <- x[o]; ll <- ll[o] - max(ll)
  fin <- seq(min(x), max(x), length.out = 801)
  yf <- stats::spline(x, ll, xout = fin, method = "natural")$y
  dedans <- fin[yf >= max(yf) - stats::qchisq(0.95, 1) / 2]
  c(est = fin[which.max(yf)], bas = min(dedans), haut = max(dedans),
    au_bord = as.numeric(fin[which.max(yf)] %in% range(x)))
}
est <- do.call(rbind, lapply(split(res[res$grille %in% c("r3", "R6") & !is.na(res$logLik), ],
                                   ~ n_rep + graine + grille, drop = TRUE), function(q) {
  pd <- stats::aggregate(logLik ~ d, q, max); pl <- stats::aggregate(logLik ~ lambda, q, max)
  ed <- profil1(pd$d, pd$logLik); el <- profil1(pl$lambda, pl$logLik)
  j <- which.max(q$logLik)
  data.frame(n_rep = q$n_rep[1], graine = q$graine[1], grille = q$grille[1],
             d_chapeau = ed[["est"]], d_bas = ed[["bas"]], d_haut = ed[["haut"]], d_au_bord = ed[["au_bord"]],
             l_chapeau = el[["est"]], l_bas = el[["bas"]], l_haut = el[["haut"]], l_au_bord = el[["au_bord"]],
             l_grille = q$lambda[j], d_grille = q$d[j], logLik_max = q$logLik[j])
}))
# logLik au vrai noyau et au maximum de R6, pour comparer les deux fenetres
vraiLL <- res[res$grille == "r3" & res$lambda == vrai$lambda & res$d == vrai$d, c("n_rep", "graine", "logLik")]
names(vraiLL)[3] <- "logLik_vrai"
est <- merge(est, vraiLL, by = c("n_rep", "graine"))
est$gain_sur_vrai <- est$logLik_max - est$logLik_vrai
utils::write.csv(est, file.path(sortie, "estimations_lambda_d.csv"), row.names = FALSE)

synth <- do.call(rbind, lapply(split(est, ~ grille + n_rep, drop = TRUE), function(q) data.frame(
  grille = q$grille[1], n_rep = q$n_rep[1], jeux = nrow(q),
  d_moyen = mean(q$d_chapeau), d_biais = mean(q$d_chapeau) - vrai$d,
  d_rmse = sqrt(mean((q$d_chapeau - vrai$d)^2)), d_largeur_ic = mean(q$d_haut - q$d_bas),
  d_couverture = mean(q$d_bas <= vrai$d & q$d_haut >= vrai$d), d0_exclu = mean(q$d_bas > 0),
  d_au_bord = mean(q$d_au_bord == 1),
  l_moyen = mean(q$l_chapeau), l_biais = mean(q$l_chapeau) - vrai$lambda,
  l_rmse = sqrt(mean((q$l_chapeau - vrai$lambda)^2)), l_largeur_ic = mean(q$l_haut - q$l_bas),
  l_couverture = mean(q$l_bas <= vrai$lambda & q$l_haut >= vrai$lambda), l_au_bord = mean(q$l_au_bord == 1),
  gain_logLik = stats::median(q$gain_sur_vrai))))
synth <- synth[order(synth$grille != "r3", synth$n_rep), ]
dire("\n1. Estimation de lambda (vrai %.2f) et d (vrai %.2f) par profil", vrai$lambda, vrai$d)
capture <- utils::capture.output(print(synth, digits = 3, row.names = FALSE))
dire("%s", paste(capture, collapse = "\n"))
utils::write.csv(synth, file.path(sortie, "synthese_lambda_d.csv"), row.names = FALSE)

# --- 2. precision des IGE selon le noyau ajuste --------------------------------
# Noyau estime : point de grille le plus haut de r3 (et de R6)
sel_est <- merge(res[res$grille %in% c("r3", "R6"), c("n_rep", "graine", "grille", "lambda", "d")],
                 est[, c("n_rep", "graine", "grille", "l_grille", "d_grille")])
sel_est <- sel_est[sel_est$lambda == sel_est$l_grille & sel_est$d == sel_est$d_grille, ]
cle <- function(x) paste(x$n_rep, x$graine, x$grille, x$r, x$lambda, x$d)
L_est <- long[cle(long) %in% cle(merge(sel_est, unique(res[, c("n_rep", "graine", "grille", "r", "lambda", "d")]))), ]
L_est$modele <- ifelse(L_est$grille == "r3", "r = 3, lambda et d estimes", "R = 6, lambda et d estimes")
L_vrai <- long[long$grille == "r3" & long$lambda == vrai$lambda & long$d == vrai$d, ]
L_vrai$modele <- "noyau vrai"
L_naif <- long[long$grille == "naif", ]; L_naif$modele <- "naif (rang 1, poids egaux)"
L_unif <- long[long$grille == "r3" & long$lambda == 1 & long$d == 0, ]
L_unif$modele <- "r = 3, lambda = 1, d = 0 (chapitre 3)"
cmp <- rbind(L_vrai, L_est, L_unif, L_naif)
cmp$modele <- factor(cmp$modele, levels = c("noyau vrai", "r = 3, lambda et d estimes",
                                            "R = 6, lambda et d estimes", "r = 3, lambda = 1, d = 0 (chapitre 3)",
                                            "naif (rang 1, poids egaux)"))
tab_cmp <- stats::aggregate(cbind(prec, err) ~ modele + effet + n_rep, cmp, stats::median, na.action = na.pass)
dire("\n2. Precision mediane et erreur log10 mediane selon le noyau ajuste")
dire("%s", paste(utils::capture.output(print(tab_cmp[order(tab_cmp$n_rep, tab_cmp$effet, tab_cmp$modele), ],
                                             digits = 3, row.names = FALSE)), collapse = "\n"))
utils::write.csv(tab_cmp, file.path(sortie, "precision_par_modele.csv"), row.names = FALSE)

# --- 3. courbes de mauvaise specification -------------------------------------
L_l <- long[long$grille == "r3" & long$d == vrai$d, ]; L_l$axe <- "lambda (r = 3, d = 0,5)"; L_l$x <- L_l$lambda
L_d <- long[long$grille == "r3" & long$lambda == vrai$lambda, ]; L_d$axe <- "d (r = 3, lambda = 1,5)"; L_d$x <- L_d$d
L_r <- long[(long$grille %in% c("r", "R6") & long$lambda == vrai$lambda & long$d == vrai$d) |
            (long$grille == "r3" & long$lambda == vrai$lambda & long$d == vrai$d), ]
L_r$axe <- "r (lambda = 1,5, d = 0,5)"; L_r$x <- L_r$r
mis <- rbind(L_l, L_d, L_r)
mis$axe <- factor(mis$axe, levels = c("lambda (r = 3, d = 0,5)", "d (r = 3, lambda = 1,5)", "r (lambda = 1,5, d = 0,5)"))
q3 <- function(v) c(med = stats::median(v, na.rm = TRUE), q1 = unname(stats::quantile(v, 0.25, na.rm = TRUE)),
                    q3 = unname(stats::quantile(v, 0.75, na.rm = TRUE)))
agr <- function(var) {
  m <- data.frame(mis[, c("axe", "effet", "n_rep")], xv = mis$x, val = mis[[var]])
  a <- stats::aggregate(val ~ axe + xv + effet + n_rep, m, q3)
  data.frame(a[, c("axe", "effet", "n_rep")], x = a$xv, med = a$val[, "med"], q1 = a$val[, "q1"],
             q3 = a$val[, "q3"])
}
mp <- agr("prec"); me <- agr("err"); me <- me[me$effet != "TBV intra", ]
utils::write.csv(merge(mp, me, by = c("axe", "x", "effet", "n_rep"), all.x = TRUE, suffixes = c("_prec", "_err")),
                 file.path(sortie, "mauvaise_specification.csv"), row.names = FALSE)
vx <- data.frame(axe = levels(mis$axe), x = c(vrai$lambda, vrai$d, vrai$r))
vx$axe <- factor(vx$axe, levels = levels(mis$axe))

theme_set(theme_minimal(base_size = 10) +
          theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                plot.title = element_text(face = "bold", size = 10), strip.text = element_text(size = 9)))
col_rep <- c("5" = "#9ecae1", "10" = "#4292c6", "20" = "#08306b")

fig_mis <- function(m, ylab, titre, zero = NULL) {
  p <- ggplot(m, aes(x, med, colour = factor(n_rep), fill = factor(n_rep))) +
    geom_vline(data = vx, aes(xintercept = x), linetype = 2, colour = "grey55") +
    geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.15, colour = NA) +
    geom_line() + geom_point(size = 1.2) +
    facet_grid(effet ~ axe, scales = "free") +
    scale_colour_manual(values = col_rep, name = "repetitions") +
    scale_fill_manual(values = col_rep, name = "repetitions") +
    labs(title = titre, x = NULL, y = ylab)
  if (!is.null(zero)) p <- p + geom_hline(yintercept = zero, colour = "grey40", linewidth = 0.3)
  p
}
ggsave(file.path(sortie, "fig1_precision_noyau_faux.png"),
       fig_mis(mp, "correlation BLUP / valeur vraie (mediane et interquartile)",
               "Precision des valeurs genetiques quand un parametre du noyau est faux (tirets : valeur vraie)"),
       width = 10, height = 8.5, dpi = 150, bg = "white")
ggsave(file.path(sortie, "fig2_erreur_variance_noyau_faux.png"),
       fig_mis(me, "log10(V estimee x k / V vraie x k)",
               "Erreur sur les variances, echelle phenotypique, quand un parametre du noyau est faux", zero = 0),
       width = 10, height = 6.8, dpi = 150, bg = "white")

# --- 4. estimations de lambda et d ---------------------------------------------
est$fenetre <- factor(ifelse(est$grille == "r3", "r = 3 (vrai rayon)", "R = 6 (fenetre large)"),
                      levels = c("r = 3 (vrai rayon)", "R = 6 (fenetre large)"))
pe <- function(v, vr, nom) ggplot(est, aes(factor(n_rep), .data[[v]], fill = fenetre)) +
  geom_hline(yintercept = vr, linetype = 2, colour = "grey45") +
  geom_boxplot(outlier.size = 0.6, width = 0.6, position = position_dodge(0.7)) +
  scale_fill_manual(values = c("#41ab5d", "#fd8d3c"), name = NULL) +
  labs(x = "repetitions par genotype", y = nom)
ggsave(file.path(sortie, "fig3_estimation_lambda_d.png"),
       (pe("l_chapeau", vrai$lambda, "lambda estime") + ggtitle("a  lambda (tirets : vrai 1,5)")) +
       (pe("d_chapeau", vrai$d, "d estime") + ggtitle("b  d (tirets : vrai 0,5)")) +
       plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 10, height = 4.5, dpi = 150, bg = "white")

# --- 5. cartes lambda x d : logLik moyenne relative, precision et erreur IGE ----
r3 <- res[res$grille == "r3" & !is.na(res$logLik), ]
r3$dll <- r3$logLik - stats::ave(r3$logLik, r3$n_rep, r3$graine, FUN = max)
r3$prec_ige <- (r3$prec_intra_A + r3$prec_intra_B + r3$prec_inter_A + r3$prec_inter_B) / 4
r3$err_ige <- (r3$err_intra_A + r3$err_intra_B + r3$err_inter_A + r3$err_inter_B) / 4
carte <- stats::aggregate(cbind(dll, prec_ige, err_ige) ~ lambda + d + n_rep, r3, stats::median)
crete <- data.frame(lambda = seq(0.5, 3, 0.05)); crete$d <- vrai$d + (vrai$lambda - crete$lambda) / 2
crete <- crete[crete$d >= 0 & crete$d <= 1, ]
pc <- function(v, nom, pal, lim = NULL) ggplot(carte, aes(lambda, d, fill = .data[[v]])) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", .data[[v]])), size = 2.4) +
  geom_line(data = crete, aes(lambda, d), inherit.aes = FALSE, linetype = 3) +
  annotate("point", x = vrai$lambda, y = vrai$d, shape = 4, size = 3) +
  facet_wrap(~ factor(sprintf("%d repetitions", n_rep), levels = sprintf("%d repetitions", sort(unique(n_rep)))), nrow = 1) +
  (if (pal == "RdBu") scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = nom)
   else scale_fill_distiller(palette = pal, direction = 1, limits = lim, name = nom)) +
  labs(x = "lambda", y = "d") + theme(panel.grid = element_blank())
ggsave(file.path(sortie, "fig4_cartes_lambda_d.png"),
       (pc("dll", "delta logLik", "Greens", c(-15, 0)) + ggtitle("a  logLik relative au maximum du jeu (mediane)")) /
       (pc("prec_ige", "precision", "Blues") + ggtitle("b  precision moyenne des IGE intra et inter")) /
       (pc("err_ige", "log10 erreur", "RdBu") + ggtitle("c  erreur log10 moyenne des variances IGE x k")) +
       plot_annotation(caption = "croix : noyau vrai ; pointilles : d + lambda / 2 = constante (approximation du disque)"),
       width = 11, height = 10.5, dpi = 150, bg = "white")
close(journal)
cat("Sorties dans", normalizePath(sortie), "\n")
