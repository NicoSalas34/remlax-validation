# Bilan de test_noyau_pondere.R sur plusieurs graines.
#   Rscript test_noyau_pondere_bilan.R <dossier des verite_*_graine_*.csv>
dossier <- commandArgs(TRUE)[1]
f <- list.files(dossier, "^verite_._graine_.*\\.csv$", full.names = TRUE)
x <- do.call(rbind, lapply(f, read.csv))
seuil <- qchisq(0.95, 1) / 2
prof1 <- function(xv, ll) {
  sp <- splinefun(xv, ll, method = "natural")
  xx <- seq(min(xv), max(xv), length.out = 2001); yy <- sp(xx)
  ic <- range(xx[yy >= max(yy) - seuil])
  c(est = xx[which.max(yy)], bas = ic[1], haut = ic[2])
}
moy_sp <- function(z, e) (z[[paste0("prec_", e, "_A")]] + z[[paste0("prec_", e, "_B")]]) / 2
sink(file.path(dossier, "bilan.txt"), split = TRUE)
cat(length(f), "jeux,", nrow(x), "ajustements,", sum(is.na(x$logLik)), "echec(s),",
    sum(!x$conv, na.rm = TRUE), "non converge(s)\n")

# --- 1. profils de lambda et d par jeu ----------------------------------------
est <- do.call(rbind, lapply(split(x, list(x$verite, x$graine)), function(z) {
  p <- z[z$grille == "pond", ]
  M <- tapply(p$logLik, list(p$lambda, p$d), identity)
  pd <- prof1(as.numeric(colnames(M)), apply(M, 2, max))
  pl <- prof1(as.numeric(rownames(M)), apply(M, 1, max))
  nb <- z[z$grille == "nombre", ]; pn <- prof1(nb$d, nb$logLik)
  r8 <- p$logLik[p$lambda == 1 & p$d == 0.5]; fe <- z[z$grille == "fenetre", ]
  data.frame(verite = z$verite[1], graine = z$graine[1],
             l = pl[1], l_bas = pl[2], l_haut = pl[3], d = pd[1], d_bas = pd[2], d_haut = pd[3],
             dn = pn[1], dn_bas = pn[2], dn_haut = pn[3],
             dLL_R4 = fe$logLik[fe$r == 4] - r8, dLL_R6 = fe$logLik[fe$r == 6] - r8,
             dLL_R12 = fe$logLik[fe$r == 12] - r8)
}))
utils::write.csv(est, file.path(dossier, "estimations.csv"), row.names = FALSE)
cat("\n1. Profils (exponentiel + masse, R = 8). Verite A : lambda = 1, d = 0,5.\n")
for (v in c("A", "B")) {
  e <- est[est$verite == v, ]
  cat(sprintf("\nverite %s (%d jeux)\n", v, nrow(e)))
  cat(sprintf("  lambda : moyenne %.3f, RMSE %.3f, largeur IC %.3f, couverture de 1 %.2f\n",
              mean(e$l), sqrt(mean((e$l - 1)^2)), mean(e$l_haut - e$l_bas), mean(e$l_bas <= 1 & e$l_haut >= 1)))
  cat(sprintf("  d masse : moyenne %.3f, RMSE %.3f, largeur IC %.3f, couverture de 0,5 %.2f, d = 0 exclu %.2f, d au bord %.2f\n",
              mean(e$d), sqrt(mean((e$d - 0.5)^2)), mean(e$d_haut - e$d_bas),
              mean(e$d_bas <= 0.5 & e$d_haut >= 0.5), mean(e$d_bas > 0), mean(e$d %in% c(0, 1))))
  cat(sprintf("  d nombre (lambda = 1) : moyenne %.3f, largeur IC %.3f, d = 0 exclu %.2f, d au bord %.2f\n",
              mean(e$dn), mean(e$dn_haut - e$dn_bas), mean(e$dn_bas > 0), mean(e$dn %in% c(0, 1))))
  cat(sprintf("  fenetre, logLik - logLik(R = 8) : R = 4 med %.3f max|.| %.3f ; R = 6 med %.3f max|.| %.3f ; R = 12 med %.3f max|.| %.3f\n",
              median(e$dLL_R4), max(abs(e$dLL_R4)), median(e$dLL_R6), max(abs(e$dLL_R6)),
              median(e$dLL_R12), max(abs(e$dLL_R12))))
}

# --- 2. modeles compares --------------------------------------------------------
cmp <- do.call(rbind, lapply(split(x, list(x$verite, x$graine)), function(z) {
  p <- z[z$grille == "pond", ]
  vraie <- if (z$verite[1] == "A") p[p$lambda == 1 & p$d == 0.5, ] else z[z$grille == "puis", ]
  m <- list("noyau vrai" = vraie, "exp + masse, meilleur point" = p[which.max(p$logLik), ],
            "exp + masse, lambda 1 d 0,5" = p[p$lambda == 1 & p$d == 0.5, ],
            "puissance r 3, nombre" = z[z$grille == "puis", ],
            "chapitre 3" = z[z$grille == "chap3", ], "naif" = z[z$grille == "naif", ])
  do.call(rbind, lapply(names(m), function(k) data.frame(
    verite = z$verite[1], graine = z$graine[1], modele = k, dLL = m[[k]]$logLik - vraie$logLik,
    dge = moy_sp(m[[k]], "dge"), intra = moy_sp(m[[k]], "intra"), inter = moy_sp(m[[k]], "inter"),
    tbv = moy_sp(m[[k]], "tbv"),
    err_intra = (m[[k]]$err_intra_A + m[[k]]$err_intra_B) / 2,
    err_inter = (m[[k]]$err_inter_A + m[[k]]$err_inter_B) / 2)))
}))
utils::write.csv(cmp, file.path(dossier, "modeles.csv"), row.names = FALSE)
cat("\n2. Modeles compares (medianes ; precision = moyenne des deux especes ; dLL = logLik - logLik du noyau vrai)\n")
ag <- aggregate(cbind(dLL, dge, intra, inter, tbv, err_intra, err_inter) ~ verite + modele, cmp, median)
ag <- ag[order(ag$verite, -ag$dLL), ]
print(cbind(ag[, 1:2], round(ag[, -(1:2)], 3)), row.names = FALSE)
cat("\nPart des jeux ou le modele bat le noyau vrai en logLik, et ecart de TBV au noyau vrai (mediane, quartiles)\n")
for (v in c("A", "B")) for (k in setdiff(unique(cmp$modele), "noyau vrai")) {
  z <- cmp[cmp$verite == v & cmp$modele == k, ]; tv <- cmp[cmp$verite == v & cmp$modele == "noyau vrai", ]
  dt <- z$tbv - tv$tbv[match(z$graine, tv$graine)]
  cat(sprintf("  %s | %-28s bat le vrai %.2f ; dTBV %.3f [%.3f ; %.3f]\n", v, k, mean(z$dLL > 0),
              median(dt), quantile(dt, .25), quantile(dt, .75)))
}
sink()
