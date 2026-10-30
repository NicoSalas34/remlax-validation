# ==============================================================================
# Demonstration de remlax sur CPU : essai au champ simule avec tendance spatiale
# lissee par P-splines 2D, valeurs genetiques directes (DGE) et effets genetiques
# indirects (IGE) entre parcelles voisines, avec parente genomique.
#
# Dispositif : 150 genotypes, 3 repetitions, blocs incomplets de 5 parcelles
# (plan resolvable de type alpha), champ de 30 rangs x 15 colonnes.
# Chaque parcelle recoit l'IGE des genotypes de ses 4 voisines immediates
# (meme rang ou meme colonne), comme dans la competition entre parcelles
# adjacentes d'un essai varietal. L'incidence de voisinage est normalisee par
# ligne (norme L2), et la variance IGE est fixee par sa part de variance
# phenotypique, via l'exposition k de rx_exposure().
#
# Tous les modeles portent rep (fixe), bloc et la surface spl2d(col, rang).
# Ils different par la partie genetique :
#   A : DGE seul, vm(gid, K)
#   B : DGE + IGE independants, vm(gid, K) + mm(Zn, K)
#   C : DGE et IGE correles, str(~ gid + mm(Zn), struct = "us", K)
#
# Seconde partie (section 6) : le meme essai, avec 10 % de parcelles vides et
# un voisinage inconnu. Les poids valent delta^-lambda dans un rayon de r pas de
# grille, divises par n_i^d (dilution, n_i = nombre de voisines). La seconde
# partie retrouve r, lambda et d par profil de la vraisemblance REML (grille
# grossiere sur r, lambda, d, puis grille fine sur lambda, d au meilleur r), puis
# par maximisation directe de ce profil en lambda et d. Elle ajuste enfin le
# modele C avec le voisinage estime.
#
# Lancer depuis un dossier sans renv (remlax est dans la bibliotheque R de
# l'utilisateur) :
#   Rscript demo_remlax_spatial_cpu.R [graine] [dossier_sortie] [partie2 TRUE/FALSE]
# Sorties : quatre figures PNG et trois tableaux CSV dans le dossier de sortie.
# ==============================================================================

suppressPackageStartupMessages({
  library(remlax)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
graine <- if (length(args) >= 1) as.integer(args[1]) else 5L
sortie <- if (length(args) >= 2) args[2] else "sortie_demo_ige"
partie2 <- if (length(args) >= 3) as.logical(args[3]) else TRUE
dir.create(sortie, showWarnings = FALSE, recursive = TRUE)
set.seed(graine)

# ------------------------------------------------------------------------------
# 1. Valeurs vraies
# ------------------------------------------------------------------------------
# `ige` est la part de variance phenotypique portee par les IGE (sigma2_IGE x k),
# `r_dge_ige` la correlation genetique entre effet direct et effet indirect :
# negative, un genotype vigoureux deprime ses voisines.
vrai <- list(dge = 1.0, ige = 0.4, r_dge_ige = -0.5, bloc = 0.3, spatial = 1.2, residuelle = 0.8)
n_gen <- 150L; n_rep <- 3L; k_bloc <- 5L
n_col <- 15L; n_rang <- 30L                     # chaque repetition : 10 rangs
rang_rep <- n_rang / n_rep

# Marqueurs : deux sous-populations (derive de Balding-Nichols, Fst = 0,1),
# puis parente genomique de VanRaden.
n_snp <- 1500L
p0 <- runif(n_snp, 0.05, 0.5)
pop <- rep(1:2, length.out = n_gen)
P <- sapply(p0, function(p) rbeta(2, p * 9, (1 - p) * 9))
M <- matrix(rbinom(n_gen * n_snp, 2L, P[pop, ]), n_gen, n_snp)
noms_gen <- sprintf("G%03d", seq_len(n_gen))
rownames(M) <- noms_gen
K <- rx_grm(M, blend = 0.02, coding = "count")

# Tendance spatiale vraie : gradient de fertilite le long des rangs, une
# depression humide localisee et une ondulation entre colonnes. Elle est
# centree et mise a l'echelle sur les parcelles pour que sa variance vaille
# vrai$spatial ; la meme transformation s'applique sur la grille fine.
tendance_brute <- function(x, y)
  0.9 * (y - 15.5) / 14.5 -
  1.6 * exp(-((x - 11)^2 / 10 + (y - 9)^2 / 30)) +
  0.35 * sin(2 * pi * x / 10)

# ------------------------------------------------------------------------------
# 2. Dispositif et voisinage
# ------------------------------------------------------------------------------
d <- expand.grid(col = seq_len(n_col), rang = seq_len(n_rang))
d$rep <- factor(sprintf("R%d", (d$rang - 1L) %/% rang_rep + 1L))
d$bloc <- factor(sprintf("%s_r%02d_b%d", d$rep, d$rang, (d$col - 1L) %/% k_bloc + 1L))
d$parcelle <- sprintf("P%03d", seq_len(nrow(d)))
# Chaque genotype une fois par repetition, tire au hasard dans la repetition.
# Un vrai plan alpha optimise aussi les rencontres entre genotypes ; ce tirage
# garde la structure (repetitions completes, blocs incomplets) sans cette
# optimisation.
d$gid <- factor(NA, levels = noms_gen)
for (r in levels(d$rep)) d$gid[d$rep == r] <- sample(noms_gen)

# Voisines : les parcelles a un pas de grille en distance euclidienne, donc les
# 4 voisines de rang et de colonne (3 ou 2 en bordure). Les repetitions se
# touchent au champ : les voisines d'une autre repetition comptent aussi.
nb <- rx_neighbourhood(coord = cbind(d$rang, d$col), group = rep("A", nrow(d)), level = d$gid,
                       id = d$parcelle, rank = 1, window = "euclidean", kernel = "none",
                       normalise = TRUE)
Zn <- nb$level[["A<-A"]][, noms_gen]            # parcelles x genotypes voisins
Zd <- Matrix::sparseMatrix(i = seq_len(nrow(d)), j = as.integer(d$gid), x = 1,
                           dims = c(nrow(d), n_gen), dimnames = list(NULL, noms_gen))
# Exposition : la variance phenotypique apportee par un terme vaut sigma2 x k.
expo <- rx_exposure(Zd, Zn, K = K)
k_ige <- expo$k

# Valeurs genetiques vraies : (DGE, IGE) ~ N(0, Sigma (x) K)
s2_ige <- vrai$ige / k_ige
Sigma_vrai <- matrix(c(vrai$dge, rep(vrai$r_dge_ige * sqrt(vrai$dge * s2_ige), 2), s2_ige), 2)
U <- t(chol(K)) %*% matrix(rnorm(2 * n_gen), n_gen, 2) %*% chol(Sigma_vrai)
dge_vrai <- setNames(U[, 1], noms_gen)
ige_vrai <- setNames(U[, 2], noms_gen)
# Valeur totale d'un genotype en culture pure : son effet direct, plus l'effet
# qu'il exerce sur toutes les parcelles qui le voient. s = poids total emis par
# parcelle, en moyenne (environ 2 avec 4 voisines normalisees en L2).
s_emis <- sum(Zn) / nrow(d)
tbv_vrai <- dge_vrai + s_emis * ige_vrai

br <- tendance_brute(d$col, d$rang)
m_br <- mean(br); s_br <- sd(br)
tendance_vraie <- function(x, y) (tendance_brute(x, y) - m_br) / s_br * sqrt(vrai$spatial)

d$spatial_vrai <- tendance_vraie(d$col, d$rang)
d$g_vrai <- dge_vrai[as.character(d$gid)]
d$ige_recu <- as.numeric(Zn %*% ige_vrai)
b_vrai <- setNames(rnorm(nlevels(d$bloc), 0, sqrt(vrai$bloc)), levels(d$bloc))
d$bloc_vrai <- b_vrai[as.character(d$bloc)]
d$e_vrai <- rnorm(nrow(d), 0, sqrt(vrai$residuelle))
d$y <- 6 + 0.2 * (d$rep == "R2") - 0.1 * (d$rep == "R3") +
  d$spatial_vrai + d$g_vrai + d$ige_recu + d$bloc_vrai + d$e_vrai

cat(sprintf("Dispositif : %d parcelles, %d genotypes, %d repetitions, %d blocs de %d\n",
            nrow(d), n_gen, n_rep, nlevels(d$bloc), k_bloc))
cat(sprintf("Voisinage : %s voisines par parcelle ; exposition k = %.3f ; poids emis s = %.2f\n",
            paste(names(table(nb$n_neighbours[[1]])), collapse = "/"), k_ige, s_emis))

# ------------------------------------------------------------------------------
# 3. Ajustements
# ------------------------------------------------------------------------------
ajuster <- function(genetique) {
  t0 <- Sys.time()
  random <- stats::as.formula(paste("~", genetique, "+ bloc + spl2d(col, rang, nseg = c(7, 10))"))
  f <- rx_reml(y ~ rep, random = random, data = d, backend = "cpu", verbose = FALSE)
  f$duree <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  f
}
fits <- list(
  "A  DGE seul"            = ajuster("vm(gid, K = K)"),
  "B  DGE + IGE"           = ajuster("vm(gid, K = K) + mm(Zn, K = K, name = 'ige')"),
  "C  DGE + IGE correles"  = ajuster("str(~ gid + mm(Zn, name = 'ige'), struct = 'us', K = K, name = 'dge_ige')"))
noms_mod <- names(fits)
couleurs <- setNames(c("#c2a5cf", "#999999", "#1b7837"), noms_mod)
fC <- fits[[3]]
print(summary(fC))
# Les trois criteres derriere la ligne Convergence de summary()
for (m in noms_mod) with(fits[[m]], cat(sprintf(
  "%-24s criteres : decrement %s, gradient %s, hessien %s (%s valeur(s) propre(s) negative(s), %s direction(s) plate(s)) ; %d variance(s) a la borne\n",
  m, conv_decrement, conv_grad_rel, conv_hessien_ok, n_neg_eig, n_null_dir %||% NA, n_at_bound)))

# ------------------------------------------------------------------------------
# 4. Ce que l'inference retrouve
# ------------------------------------------------------------------------------
# BLUP par genotype : DGE et IGE (zero pour A, qui n'a pas d'IGE).
blups_de <- function(f) {
  b <- f$blups
  if (!is.null(b$dge_ige)) return(list(dge = b$dge_ige[noms_gen, 1], ige = b$dge_ige[noms_gen, 2]))
  list(dge = setNames(as.numeric(b$gid), rownames(b$gid))[noms_gen],
       ige = if (is.null(b$ige)) setNames(rep(0, n_gen), noms_gen)
             else setNames(as.numeric(b$ige), rownames(b$ige))[noms_gen])
}
g_est <- lapply(fits, blups_de)
top <- function(v, p = 0.1) names(sort(v, decreasing = TRUE))[seq_len(round(p * length(v)))]

# Composantes genetiques sur l'echelle phenotypique : sigma2_DGE, sigma2_IGE x k,
# et la correlation r(DGE, IGE).
compo_gen <- function(f) {
  v <- summary(f)$varcomp
  val <- function(nm) if (nm %in% rownames(v)) v[nm, "component"] else NA_real_
  et <- function(nm) if (nm %in% rownames(v)) v[nm, "std.error"] else NA_real_
  if ("dge_ige[1,1]" %in% rownames(v)) {
    s11 <- val("dge_ige[1,1]"); s21 <- val("dge_ige[2,1]"); s22 <- val("dge_ige[2,2]")
    return(data.frame(composante = c("DGE", "IGE x k"), estimation = c(s11, s22 * k_ige),
                      se = c(et("dge_ige[1,1]"), et("dge_ige[2,2]") * k_ige),
                      r = s21 / sqrt(s11 * s22)))
  }
  data.frame(composante = c("DGE", "IGE x k"), estimation = c(val("gid"), val("ige") * k_ige),
             se = c(et("gid"), et("ige") * k_ige), r = NA_real_)
}

# Environnement estime par C : surface spline + effets de repetition, centre sur
# les parcelles. Les repetitions absorbent une partie du gradient le long des
# rangs, d'ou leur ajout avant la comparaison avec la tendance vraie.
cf <- summary(fC)$coef.fixed
beta_rep <- c(R1 = 0, R2 = cf["repR2", "solution"], R3 = cf["repR3", "solution"])
tendance_C <- function(nx, ny) {
  tr <- rx_spatial_trend(fC, grid = c(nx, ny))
  tr$trend <- tr$trend + beta_rep[pmin(n_rep, floor((tr$rang - 0.5) / rang_rep) + 1L)]
  tr
}
tr_p <- tendance_C(n_col, n_rang)
moy_C <- mean(tr_p$trend)
d$spatial_est <- tr_p$trend[match(paste(d$col, d$rang), paste(tr_p$col, tr_p$rang))] - moy_C
for (i in seq_along(fits)) d[[paste0("res", LETTERS[i])]] <- residuals(fits[[i]])

tab <- do.call(rbind, lapply(noms_mod, function(m) {
  f <- fits[[m]]; g <- g_est[[m]]; s <- summary(f)
  tbv <- g$dge + s_emis * g$ige
  data.frame(modele = m, logLik = f$logLik, AIC = s$aic, converge = isTRUE(s$converged),
             secondes = round(f$duree, 1),
             r_dge = cor(g$dge, dge_vrai),
             r_ige = if (all(g$ige == 0)) NA_real_ else cor(g$ige, ige_vrai),
             r_tbv = cor(tbv, tbv_vrai),
             selection_10pc = length(intersect(top(tbv), top(tbv_vrai))) / length(top(tbv_vrai)))
}))
write.csv(tab, file.path(sortie, "comparaison_modeles.csv"), row.names = FALSE)

comp <- do.call(rbind, lapply(noms_mod, function(m) {
  v <- summary(fits[[m]])$varcomp
  data.frame(modele = m, composante = rownames(v), estimation = v$component, se = v$std.error)
}))
write.csv(comp, file.path(sortie, "composantes_variance.csv"), row.names = FALSE)
gen <- do.call(rbind, lapply(noms_mod, function(m) cbind(modele = m, compo_gen(fits[[m]]))))

cat("\nComparaison des modeles\n"); print(tab, digits = 4, row.names = FALSE)
cat("\nComposantes genetiques, echelle phenotypique (vraies : DGE", vrai$dge, ", IGE x k", vrai$ige,
    ", r", vrai$r_dge_ige, ")\n"); print(gen, digits = 3, row.names = FALSE)
cat(sprintf("\nCorrelation tendance estimee (C) / vraie sur les parcelles : %.3f\n",
            cor(d$spatial_est, d$spatial_vrai)))

# ------------------------------------------------------------------------------
# 5. Figures
# ------------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 10) +
          theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10),
                legend.key.height = unit(0.5, "cm")))
lim_sym <- function(...) { m <- max(abs(c(...))); c(-m, m) }
carte <- function(df, var, titre, lim = NULL) {
  ggplot(df, aes(col, rang)) + geom_tile(aes(fill = .data[[var]])) +
    scale_fill_distiller(palette = "RdBu", direction = 1, limits = lim, name = NULL) +
    coord_equal(expand = FALSE) + labs(title = titre, x = "colonne", y = "rang")
}
surface <- function(df, z, titre, lim) {
  ggplot(df, aes(col, rang)) + geom_raster(aes(fill = .data[[z]])) +
    geom_contour(aes(z = .data[[z]]), colour = "grey20", linewidth = 0.2, bins = 8) +
    scale_fill_distiller(palette = "RdBu", direction = 1, limits = lim, name = NULL) +
    coord_equal(expand = FALSE) + labs(title = titre, x = "colonne", y = "rang")
}

# Contours des blocs et des repetitions
bb <- aggregate(cbind(col, rang) ~ bloc, d, range)
blocs <- data.frame(xmin = bb$col[, 1] - 0.5, xmax = bb$col[, 2] + 0.5,
                    ymin = bb$rang[, 1] - 0.5, ymax = bb$rang[, 2] + 0.5)
reps <- data.frame(xmin = 0.5, xmax = n_col + 0.5, ymin = (seq_len(n_rep) - 1) * rang_rep + 0.5,
                   ymax = seq_len(n_rep) * rang_rep + 0.5, rep = levels(d$rep))
grille <- expand.grid(col = seq(1, n_col, length.out = 90), rang = seq(1, n_rang, length.out = 180))
grille$vrai <- tendance_vraie(grille$col, grille$rang)

# Figure 1 : donnees simulees. Le plan montre aussi une parcelle et ses voisines.
p0_ <- d[d$col == 8 & d$rang == 15, ]
vois <- d[abs(d$col - p0_$col) + abs(d$rang - p0_$rang) == 1, ]
p_plan <- ggplot() +
  geom_rect(data = reps, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = rep), alpha = 0.35) +
  geom_rect(data = blocs, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "grey30", linewidth = 0.3) +
  geom_tile(data = vois, aes(col, rang), fill = "#d6604d", alpha = 0.8) +
  geom_tile(data = p0_, aes(col, rang), fill = "black") +
  scale_fill_manual(values = c(R1 = "#8fbcd4", R2 = "#e3b778", R3 = "#b5d6a7"), name = NULL) +
  coord_equal(expand = FALSE) +
  labs(title = "a  Plan : 3 repetitions, blocs de 5 ;\n   une parcelle (noir) et ses 4 voisines (rouge)",
       x = "colonne", y = "rang")
p_tend <- surface(grille, "vrai", "b  Tendance spatiale vraie", lim_sym(grille$vrai))
p_g <- carte(d, "g_vrai", "c  DGE vrai de la parcelle", lim_sym(d$g_vrai))
p_i <- carte(d, "ige_recu", "d  IGE recu des voisines", lim_sym(d$ige_recu))
p_y <- ggplot(d, aes(col, rang, fill = y)) + geom_tile() +
  scale_fill_viridis_c(name = NULL) + coord_equal(expand = FALSE) +
  labs(title = "e  Phenotype observe", x = "colonne", y = "rang")
sources <- c("DGE", "IGE recu", "bloc", "spatial", "residuelle")
part <- data.frame(source = factor(sources, levels = sources),
                   variance = c(var(d$g_vrai), var(d$ige_recu), var(d$bloc_vrai), var(d$spatial_vrai),
                                var(d$e_vrai)))
p_part <- ggplot(part, aes(source, variance, fill = source)) + geom_col(width = 0.6) +
  scale_fill_manual(values = c("#1b7837", "#d6604d", "#999999", "#2166ac", "#bababa"), guide = "none") +
  labs(title = "f  Variances realisees sur les parcelles", x = NULL, y = "variance") +
  theme(panel.grid.major.y = element_line(colour = "grey90"))
fig1 <- (p_plan | p_tend | p_g) / (p_i | p_y | p_part) +
  plot_annotation(title = sprintf("Donnees simulees : %d parcelles, %d genotypes, r(DGE, IGE) = %.1f (graine %d)",
                                  nrow(d), n_gen, vrai$r_dge_ige, graine))
ggsave(file.path(sortie, "fig1_donnees_simulees.png"), fig1, width = 12, height = 10, dpi = 150, bg = "white")

# Figure 2 : la surface retrouvee par C
tr_fin <- tendance_C(90, 180)
tr_fin$trend <- tr_fin$trend - moy_C
L <- lim_sym(grille$vrai, tr_fin$trend)
p_v <- surface(grille, "vrai", "a  Tendance vraie", L)
p_e <- surface(tr_fin, "trend", "b  Tendance estimee par C (spl2d + rep)", L)
p_sc <- ggplot(d, aes(spatial_vrai, spatial_est)) +
  geom_abline(colour = "grey60", linetype = 2) + geom_point(size = 0.9, alpha = 0.6, colour = "#2166ac") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3,
           label = sprintf("r = %.3f", cor(d$spatial_vrai, d$spatial_est))) +
  coord_equal() + labs(title = "c  Par parcelle", x = "tendance vraie", y = "tendance estimee")
Lr <- lim_sym(d$resA, d$resC)
p_rA <- carte(d, "resA", "d  Residus de A (sans IGE)", Lr)
p_rC <- carte(d, "resC", "e  Residus de C (IGE correles)", Lr)
dens <- data.frame(modele = factor(rep(noms_mod, each = nrow(d)), levels = noms_mod),
                   res = c(d$resA, d$resB, d$resC))
p_dens <- ggplot(dens, aes(res, colour = modele)) + geom_density(linewidth = 0.7) +
  scale_colour_manual(values = couleurs, name = NULL) +
  labs(title = "f  Distribution des residus", x = "residu", y = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.85))
fig2 <- (p_v | p_e | p_sc) / (p_rA | p_rC | p_dens) +
  plot_annotation(title = "Tendance spatiale : vraie contre estimee")
ggsave(file.path(sortie, "fig2_surface.png"), fig2, width = 12, height = 10, dpi = 150, bg = "white")

# Figure 3 : la partie genetique
nuage <- function(vr, champ, mods, titre) {
  df <- do.call(rbind, lapply(mods, function(m)
    data.frame(vrai = vr, blup = g_est[[m]][[champ]], modele = m)))
  df$modele <- factor(df$modele, levels = noms_mod)
  ggplot(df, aes(vrai, blup, colour = modele)) +
    geom_abline(colour = "grey60", linetype = 2) + geom_point(size = 1, alpha = 0.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
    scale_colour_manual(values = couleurs, name = NULL) +
    labs(title = titre, x = "valeur vraie", y = "BLUP") +
    theme(legend.position = "inside", legend.position.inside = c(0.25, 0.88))
}
p_bd <- nuage(dge_vrai, "dge", noms_mod[c(1, 3)], "a  DGE : BLUP contre vrai")
p_bi <- nuage(ige_vrai, "ige", noms_mod[2:3], "b  IGE : BLUP contre vrai")
gen$composante <- factor(gen$composante, levels = c("DGE", "IGE x k"))
gen$modele <- factor(gen$modele, levels = noms_mod)
vv <- data.frame(composante = factor(c("DGE", "IGE x k"), levels = levels(gen$composante)),
                 vrai = c(vrai$dge, vrai$ige))
rC <- gen$r[gen$modele == noms_mod[3]][1]
p_vc <- ggplot(gen[!is.na(gen$estimation), ], aes(composante, estimation, colour = modele)) +
  geom_point(data = vv, aes(composante, vrai), inherit.aes = FALSE, shape = 95, size = 14) +
  geom_pointrange(aes(ymin = pmax(estimation - 2 * se, 0), ymax = estimation + 2 * se),
                  position = position_dodge(width = 0.6), size = 0.3) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.5, size = 3,
           label = sprintf("r(DGE, IGE) estime par C : %.2f (vrai %.2f)", rC, vrai$r_dge_ige)) +
  scale_colour_manual(values = couleurs, name = NULL) +
  labs(title = "c  Variances genetiques, echelle phenotypique", subtitle = "trait noir = valeur vraie",
       x = NULL, y = "variance (+/- 2 e.-t.)") +
  theme(legend.position = "none", panel.grid.major.y = element_line(colour = "grey90"))
ind <- rbind(data.frame(modele = tab$modele, critere = "r(DGE)", valeur = tab$r_dge),
             data.frame(modele = tab$modele, critere = "r(IGE)", valeur = tab$r_ige),
             data.frame(modele = tab$modele, critere = "r(valeur totale)", valeur = tab$r_tbv),
             data.frame(modele = tab$modele, critere = "15 meilleurs retrouves", valeur = tab$selection_10pc))
ind$modele <- factor(ind$modele, levels = noms_mod)
ind$critere <- factor(ind$critere, levels = unique(ind$critere))
p_ind <- ggplot(ind[!is.na(ind$valeur), ], aes(modele, valeur, fill = modele)) + geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", valeur)), vjust = -0.3, size = 2.8) +
  facet_wrap(~ critere, nrow = 1) + scale_fill_manual(values = couleurs, guide = "none") +
  scale_y_continuous(limits = c(0, 1.08)) + scale_x_discrete(labels = c("A", "B", "C")) +
  labs(title = "d  Precision et selection sur la valeur totale DGE + s IGE", x = NULL, y = NULL)
p_aic <- ggplot(tab, aes(factor(modele, levels = noms_mod), AIC - min(AIC), fill = modele)) +
  geom_col(width = 0.65) + scale_fill_manual(values = couleurs, guide = "none") +
  scale_x_discrete(labels = c("A", "B", "C")) +
  labs(title = "e  AIC relatif au meilleur modele", x = NULL, y = "delta AIC") +
  theme(panel.grid.major.y = element_line(colour = "grey90"))
fig3 <- (p_bd | p_bi | p_vc) / (p_ind + p_aic + plot_layout(widths = c(3, 1))) +
  plot_annotation(title = "Partie genetique : A = DGE seul, B = DGE + IGE, C = DGE et IGE correles")
ggsave(file.path(sortie, "fig3_genetique.png"), fig3, width = 14, height = 9, dpi = 150, bg = "white")

# ==============================================================================
# 6. Seconde partie : voisinage inconnu (r, lambda, d)
# ==============================================================================
# La normalisation L2 par ligne annulerait la dilution : elle divise chaque
# ligne par sa norme, apres la division par n_i^d. Les poids restent donc bruts
# ici, et d n'est identifiable que parce que n_i varie : bordures du champ et
# parcelles vides. Chaque point du profil est un ajustement REML complet du
# modele B (DGE et IGE independants), sans Hessien ni BLUP pour aller vite.
if (partie2) {
  cat("\n=== Partie 2 : voisinage a estimer ===\n")
  vrai2 <- list(r = 3L, lambda = 1.5, d = 0.5, dge = 1.0, ige = 0.4, r_dge_ige = -0.5,
                bloc = 0.3, residuelle = 0.8)
  set.seed(graine + 1000L)
  d2 <- d[sort(sample(nrow(d), round(0.9 * nrow(d)))), c("col", "rang", "rep", "bloc", "parcelle", "gid")]
  rownames(d2) <- NULL
  d2$bloc <- droplevels(d2$bloc)

  # Incidence parcelles x genotypes voisins pour (r, lambda, d) ; un genotype
  # absent du voisinage garde une colonne nulle.
  voisinage <- function(r, lambda, dd, normalise = FALSE) {
    # un scalaire nomme serait lu comme une liste nommee par paire
    r <- as.integer(unname(r)); lambda <- as.numeric(unname(lambda)); dd <- as.numeric(unname(dd))
    nb <- rx_neighbourhood(coord = cbind(d2$rang, d2$col), group = rep("A", nrow(d2)),
                           level = d2$gid, id = d2$parcelle, rank = r, reach = lambda,
                           dilution = dd, window = "euclidean", kernel = "power",
                           normalise = normalise)
    Z <- nb$level[["A<-A"]]
    P <- Matrix::sparseMatrix(i = match(colnames(Z), noms_gen), j = seq_len(ncol(Z)), x = 1,
                              dims = c(n_gen, ncol(Z)))
    Z <- Z %*% Matrix::t(P)
    dimnames(Z) <- list(NULL, noms_gen)
    methods::as(Z, "CsparseMatrix")
  }
  Zd2 <- Matrix::sparseMatrix(i = seq_len(nrow(d2)), j = as.integer(d2$gid), x = 1,
                              dims = c(nrow(d2), n_gen), dimnames = list(NULL, noms_gen))
  Zv <- voisinage(vrai2$r, vrai2$lambda, vrai2$d)
  k2 <- rx_exposure(Zd2, Zv, K = K)$k
  s2i <- vrai2$ige / k2
  Sig2 <- matrix(c(vrai2$dge, rep(vrai2$r_dge_ige * sqrt(vrai2$dge * s2i), 2), s2i), 2)
  U2 <- t(chol(K)) %*% matrix(rnorm(2 * n_gen), n_gen, 2) %*% chol(Sig2)
  dge2 <- setNames(U2[, 1], noms_gen); ige2 <- setNames(U2[, 2], noms_gen)
  s2_emis <- sum(Zv) / nrow(d2)
  tbv2 <- dge2 + s2_emis * ige2
  b2 <- setNames(rnorm(nlevels(d2$bloc), 0, sqrt(vrai2$bloc)), levels(d2$bloc))
  d2$ige_recu <- as.numeric(Zv %*% ige2)
  d2$y <- 6 + 0.2 * (d2$rep == "R2") - 0.1 * (d2$rep == "R3") + tendance_vraie(d2$col, d2$rang) +
    dge2[as.character(d2$gid)] + d2$ige_recu + b2[as.character(d2$bloc)] +
    rnorm(nrow(d2), 0, sqrt(vrai2$residuelle))
  cat(sprintf("%d parcelles semees ; voisinage vrai r = %d, lambda = %g, d = %g ; k = %.3f\n",
              nrow(d2), vrai2$r, vrai2$lambda, vrai2$d, k2))

  ajuster2 <- function(Z, correle = FALSE, leger = TRUE) {
    gen <- if (correle) "str(~ gid + mm(Z, name = 'ige'), struct = 'us', K = K, name = 'dge_ige')"
           else "vm(gid, K = K) + mm(Z, K = K, name = 'ige')"
    random <- stats::as.formula(paste("~", gen, "+ bloc + spl2d(col, rang, nseg = c(7, 10))"),
                                env = environment())
    t0 <- Sys.time()
    f <- if (leger) rx_reml(y ~ rep, random = random, data = d2, backend = "cpu", verbose = FALSE,
                            hessian = FALSE, blups = FALSE, pev = FALSE)
         else rx_reml(y ~ rep, random = random, data = d2, backend = "cpu", verbose = FALSE)
    f$duree <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    f
  }
  point <- function(r, lambda, dd) {
    f <- tryCatch(ajuster2(voisinage(r, lambda, dd)), error = function(e) conditionMessage(e))
    if (is.character(f)) return(data.frame(r = r, lambda = lambda, d = dd, logLik = NA, conv = NA,
                                           erreur = substr(gsub("[\r\n]+", " ", f), 1, 200)))
    data.frame(r = r, lambda = lambda, d = dd, logLik = f$logLik,
               conv = isTRUE(as.logical(f$conv_grad_rel)) || isTRUE(as.logical(f$conv_decrement)),
               erreur = "")
  }
  # Ajustements en parallele : un serveur remlax par processus, deux fils chacun.
  coeurs <- max(1L, min(8L, parallel::detectCores() %/% 2L))
  Sys.setenv(XLA_FLAGS = "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2")
  options(remlax.threads = 2L)
  en_parallele <- function(g) do.call(rbind, parallel::mclapply(seq_len(nrow(g)), function(i)
    point(g$r[i], g$lambda[i], g$d[i]), mc.cores = coeurs))

  # Etape 1 : grille grossiere. A r = 1, les 4 voisines sont a la meme distance
  # et lambda ne change rien.
  t0 <- Sys.time()
  g1 <- rbind(data.frame(r = 1L, lambda = 0, d = c(0, 0.5, 1)),
              expand.grid(r = 2:4, lambda = c(0, 1, 2, 3), d = c(0, 0.5, 1)))
  p1 <- en_parallele(g1)
  r_hat <- p1$r[which.max(p1$logLik)]
  cat(sprintf("Grille grossiere : %d ajustements en %.0f s ; meilleur r = %d\n", nrow(p1),
              as.numeric(difftime(Sys.time(), t0, units = "secs")), r_hat))

  # Etape 2 : grille fine en (lambda, d) au meilleur r.
  t0 <- Sys.time()
  g2 <- expand.grid(r = r_hat, lambda = if (r_hat == 1L) 0 else seq(0, 3, 0.25), d = seq(0, 1.25, 0.25))
  p2 <- en_parallele(g2)
  i2 <- which.max(p2$logLik)
  cat(sprintf("Grille fine : %d ajustements en %.0f s ; meilleur point lambda = %g, d = %g\n",
              nrow(p2), as.numeric(difftime(Sys.time(), t0, units = "secs")), p2$lambda[i2], p2$d[i2]))

  write.csv(rbind(cbind(etape = "grossiere", p1), cbind(etape = "fine", p2)),
            file.path(sortie, "profil_voisinage.csv"), row.names = FALSE)
  if (any(!is.na(c(p1$erreur, p2$erreur)) & nzchar(c(p1$erreur, p2$erreur))))
    print(unique(c(p1$erreur, p2$erreur)))

  # Etape 3 : maximisation directe du profil en (lambda, d), depuis le meilleur
  # point de la grille fine. Le gradient vient de differences avant, dont les
  # trois ajustements tournent en parallele.
  t0 <- Sys.time()
  memo <- new.env()
  eval_fg <- function(par) {
    cle <- paste(sprintf("%.6f", par), collapse = "_")
    if (!is.null(memo[[cle]])) return(memo[[cle]])
    h <- 0.02
    pts <- rbind(par, par + c(h, 0), par + c(0, h))
    res <- parallel::mclapply(1:3, function(i) point(r_hat, pts[i, 1], pts[i, 2]),
                              mc.cores = min(3L, coeurs))
    ll <- vapply(res, function(x) if (is.data.frame(x)) x$logLik else NA_real_, 0)
    if (any(!is.finite(ll))) {
      msg <- vapply(res, function(x) if (is.data.frame(x)) x$erreur else as.character(x), "")
      stop("maximisation directe : ajustement en echec a (", paste(round(par, 3), collapse = ", "),
           ") : ", paste(unique(msg), collapse = " | "))
    }
    memo[[cle]] <- list(f = -ll[1], g = -(ll[2:3] - ll[1]) / h)
    memo[[cle]]
  }
  if (r_hat == 1L) {
    o <- stats::optimize(function(x) -point(1L, 0, x)$logLik, c(0, 1.5))
    est <- c(lambda = 0, d = o$minimum); ll_est <- -o$objective; n_eval <- NA
  } else {
    o <- stats::optim(c(p2$lambda[i2], p2$d[i2]), function(p) eval_fg(p)$f, function(p) eval_fg(p)$g,
                      method = "L-BFGS-B", lower = c(0, 0), upper = c(4, 1.5),
                      control = list(maxit = 25, factr = 1e9))
    est <- c(lambda = o$par[1], d = o$par[2]); ll_est <- -o$value; n_eval <- length(ls(memo)) * 3L
  }
  cat(sprintf("Maximisation directe : lambda = %.2f, d = %.2f, logLik %.3f (%s ajustements, %.0f s)\n",
              est[["lambda"]], est[["d"]], ll_est, n_eval,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  options(remlax.threads = NULL)

  # Etape 4 : modele C complet sous trois voisinages. Le naif est celui de la
  # premiere partie (4 voisines, L2) ; l'oracle connait le voisinage vrai.
  Zs <- list("naif : r = 1, L2" = voisinage(1L, 0, 0, normalise = TRUE),
             "estime" = voisinage(r_hat, est[["lambda"]], est[["d"]]),
             "vrai (oracle)" = Zv)
  fin <- lapply(Zs, function(Z) ajuster2(Z, correle = TRUE, leger = FALSE))
  bil <- do.call(rbind, lapply(names(Zs), function(m) {
    f <- fin[[m]]; Z <- Zs[[m]]; v <- summary(f)$varcomp
    k_m <- rx_exposure(Zd2, Z, K = K)$k; s_m <- sum(Z) / nrow(d2)
    b <- f$blups$dge_ige[noms_gen, ]
    s11 <- v["dge_ige[1,1]", "component"]; s21 <- v["dge_ige[2,1]", "component"]
    s22 <- v["dge_ige[2,2]", "component"]
    data.frame(voisinage = m, logLik = f$logLik, dge = s11, ige_k = s22 * k_m,
               r_dge_ige = s21 / sqrt(s11 * s22), prec_dge = cor(b[, 1], dge2),
               prec_ige = cor(b[, 2], ige2), prec_tbv = cor(b[, 1] + s_m * b[, 2], tbv2),
               converge = isTRUE(summary(f)$converged), secondes = round(f$duree, 1))
  }))
  cat("\nModele C sous trois voisinages (vrai : DGE 1, IGE x k 0,4, r -0,5)\n")
  print(bil, digits = 3, row.names = FALSE)

  # Relation entre lambda et d a r fixe. La variance de l'IGE recu par la
  # parcelle i est proportionnelle a ||w_i||^2 = sum_j delta_ij^(-2 lambda) / n_i^(2 d).
  # L'echelle commune va dans sigma2_IGE : les donnees ne voient que la facon dont
  # ||w_i|| varie avec n_i. On ecrit sum_j delta_ij^(-2 lambda) = n_i M_i(lambda),
  # M_i moyenne de delta^(-2 lambda) sur les voisines de i, et
  #   log ||w_i|| = (1/2 - d) log n_i + (1/2) log M_i(lambda).
  # Si b(lambda) est la pente de log M_i contre log n_i, l'exposant vu par les
  # donnees vaut e = 1/2 - d + b(lambda)/2, et la courbe qui le garde constant
  # passe par le point vrai : d(lambda) = d_vrai + (b(lambda) - b(lambda_vrai)) / 2.
  # Approximation continue : si n_i varie comme l'aire d'un disque de rayon rho_i,
  # n_i ~ rho_i^2 et sum delta^(-2 lambda) ~ rho_i^(2 - 2 lambda), d'ou b = -lambda
  # et la courbe d + lambda / 2 = constante. Sur le dispositif, bordures et
  # parcelles vides coupent le voisinage comme un demi-plan : b est calcule.
  pente_b <- function(r, lambda) {
    nb <- rx_neighbourhood(cbind(d2$rang, d2$col), group = rep("A", nrow(d2)), id = d2$parcelle,
                           rank = r, reach = lambda, window = "euclidean", kernel = "power",
                           output = "unit")
    n_i <- nb$n_neighbours[[1]]; s2 <- Matrix::rowSums(nb$unit[[1]]^2); ok <- n_i > 0
    if (lambda == 0) return(0)
    unname(stats::coef(stats::lm(log(s2[ok] / n_i[ok]) ~ log(n_i[ok])))[2])
  }
  lam_g <- seq(0, 3, 0.05)
  b_g <- vapply(lam_g, function(l) pente_b(r_hat, l), 0)
  b_v <- pente_b(r_hat, vrai2$lambda); b_e <- pente_b(r_hat, est[["lambda"]])
  courbes <- rbind(
    data.frame(lambda = lam_g, d = vrai2$d + (b_g - b_v) / 2, quoi = "exposant constant, dispositif"),
    data.frame(lambda = lam_g, d = vrai2$d - (lam_g - vrai2$lambda) / 2, quoi = "d + lambda/2 constant (disque)"))
  courbes <- courbes[courbes$d >= 0 & courbes$d <= 1.25, ]
  # Crete du profil : meilleur d pour chaque lambda, par parabole sur les trois
  # points de grille autour du maximum (au bord de la grille, le point du bord).
  crete <- do.call(rbind, lapply(split(p2, p2$lambda), function(q) {
    q <- q[order(q$d), ]; j <- which.max(q$logLik)
    dd <- q$d[j]
    if (j > 1 && j < nrow(q)) {
      y <- q$logLik[(j - 1):(j + 1)]; h <- diff(q$d)[1]
      dd <- q$d[j] + h * (y[1] - y[3]) / (2 * (y[1] - 2 * y[2] + y[3]))
    }
    data.frame(lambda = q$lambda[1], d = dd)
  }))
  e_v <- 0.5 - vrai2$d + b_v / 2; e_e <- 0.5 - est[["d"]] + b_e / 2
  cat(sprintf(paste0("\nRelation lambda-d a r = %d : b(lambda) = %.3f au vrai, %.3f a l'estime ;\n",
                     "  exposant e = 1/2 - d + b/2 : %.3f au vrai, %.3f a l'estime ;\n",
                     "  d sur la courbe a lambda estime : %.3f (dispositif), %.3f (disque) ; d estime %.3f\n"),
              r_hat, b_v, b_e, e_v, e_e, vrai2$d + (b_e - b_v) / 2,
              vrai2$d - (est[["lambda"]] - vrai2$lambda) / 2, est[["d"]]))

  # Figure 4
  p1$dll <- pmax(p1$logLik - max(p1$logLik), -40)
  p_g <- ggplot(p1, aes(lambda, dll, colour = factor(d))) + geom_line(data = p1[p1$r > 1, ]) +
    geom_point(size = 1.5) +
    facet_wrap(~ paste("r =", r), nrow = 1) +
    scale_colour_manual(values = c("#d6604d", "#4393c3", "#1b7837"), name = "d") +
    labs(title = "a  Profil grossier : logLik relative au maximum (tronquee a -40)",
         x = "lambda", y = "delta logLik") +
    theme(panel.grid.major.y = element_line(colour = "grey90"))
  p2$dll <- pmax(p2$logLik - max(p2$logLik, ll_est), -15)
  seuil <- -stats::qchisq(0.95, 2) / 2
  pts <- data.frame(lambda = c(vrai2$lambda, est[["lambda"]]), d = c(vrai2$d, est[["d"]]),
                    quoi = c(sprintf("vrai (r = %d)", vrai2$r), "estime"))
  p_c <- ggplot(p2, aes(lambda, d)) + geom_raster(aes(fill = dll), interpolate = TRUE) +
    geom_contour(aes(z = dll), breaks = seuil, colour = "black", linewidth = 0.5) +
    geom_line(data = courbes, aes(linetype = quoi), colour = "white", linewidth = 0.8) +
    geom_point(data = crete, colour = "white", size = 1.2) +
    geom_point(data = pts, aes(shape = quoi), size = 3, colour = "black", fill = "white") +
    scale_shape_manual(values = c(21, 4), name = NULL) +
    scale_linetype_manual(values = c("solid", "22"), name = NULL,
                          guide = guide_legend(override.aes = list(colour = "black"))) +
    scale_fill_viridis_c(name = "delta\nlogLik") + coord_cartesian(expand = FALSE) +
    labs(title = sprintf("b  Profil fin a r = %d ; contour noir : region de confiance a 95 %%", r_hat),
         subtitle = "courbes blanches : exposant constant passant par le vrai ; points blancs : crete du profil",
         x = "lambda", y = "d") + theme(legend.position = "bottom", legend.box = "vertical")
  noyau <- function(r, lambda, dd, quoi) {
    o <- expand.grid(dx = -r:r, dy = -r:r); o$delta <- sqrt(o$dx^2 + o$dy^2)
    o <- o[o$delta > 0 & o$delta <= r, ]
    # poids relatif a la voisine la plus proche : n_i^-d est une echelle commune
    # aux parcelles interieures, absorbee par la variance IGE
    data.frame(delta = sort(unique(o$delta)), w = sort(unique(o$delta))^(-lambda), quoi = quoi)
  }
  ker <- rbind(noyau(vrai2$r, vrai2$lambda, vrai2$d, "vrai"),
               noyau(r_hat, est[["lambda"]], est[["d"]], "estime"))
  p_k <- ggplot(ker, aes(delta, w, colour = quoi)) + geom_line() + geom_point(size = 2) +
    scale_colour_manual(values = c(vrai = "black", estime = "#1b7837"), name = NULL) +
    labs(title = "c  Forme du noyau : poids relatif a la voisine la plus proche",
         x = "distance (pas de grille)", y = "poids relatif") + theme(panel.grid.major.y = element_line(colour = "grey90"))
  nu <- do.call(rbind, lapply(names(Zs), function(m)
    data.frame(vrai = ige2, blup = fin[[m]]$blups$dge_ige[noms_gen, 2], voisinage = m)))
  nu$voisinage <- factor(nu$voisinage, levels = names(Zs))
  p_n <- ggplot(nu, aes(vrai, blup)) + geom_point(size = 0.8, alpha = 0.6, colour = "#d6604d") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.5) +
    facet_wrap(~ voisinage, nrow = 1, scales = "free_y") +
    labs(title = "d  IGE : BLUP contre vrai selon le voisinage ajuste", x = "IGE vrai", y = "BLUP")
  ib <- rbind(data.frame(voisinage = bil$voisinage, critere = "r(IGE)", valeur = bil$prec_ige),
              data.frame(voisinage = bil$voisinage, critere = "r(valeur totale)", valeur = bil$prec_tbv),
              data.frame(voisinage = bil$voisinage, critere = "IGE x k (vrai 0,4)", valeur = bil$ige_k))
  ib$voisinage <- factor(ib$voisinage, levels = names(Zs))
  ib$critere <- factor(ib$critere, levels = unique(ib$critere))
  p_b <- ggplot(ib, aes(voisinage, valeur, fill = voisinage)) + geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%.2f", valeur)), vjust = -0.3, size = 2.8) +
    facet_wrap(~ critere, nrow = 1) +
    scale_fill_manual(values = c("#c2a5cf", "#1b7837", "#999999"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    scale_x_discrete(labels = c("naif", "estime", "oracle")) +
    labs(title = "e  Modele C : precision et variance IGE", x = NULL, y = NULL)
  fig4 <- p_g / (p_c | p_k) / (p_n | p_b) +
    plot_annotation(title = sprintf(
      "Voisinage inconnu : vrai r = %d, lambda = %g, d = %g ; estime r = %d, lambda = %.2f, d = %.2f",
      vrai2$r, vrai2$lambda, vrai2$d, r_hat, est[["lambda"]], est[["d"]]))
  ggsave(file.path(sortie, "fig4_voisinage_estime.png"), fig4, width = 14, height = 13, dpi = 150,
         bg = "white")
}

cat("\nFigures et tableaux ecrits dans", normalizePath(sortie), "\n")
