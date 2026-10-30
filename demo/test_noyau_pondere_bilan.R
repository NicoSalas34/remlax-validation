# Summary of test_noyau_pondere.R over several seeds.
#   Rscript test_noyau_pondere_bilan.R <folder of the truth_*_seed_*.csv files>
# Also reads the verite_*_graine_*.csv files of the earlier French version,
# whose columns and values are mapped to English on reading.
dir <- commandArgs(TRUE)[1]
f <- list.files(dir, "^truth_._seed_.*\\.csv$", full.names = TRUE)
if (!length(f)) f <- list.files(dir, "^verite_._graine_.*\\.csv$", full.names = TRUE)
fr_cols <- c(verite = "truth", graine = "seed", grille = "grid", noyau = "kernel", fenetre = "window",
             secondes = "seconds", erreur = "error")
fr_vals <- list(grid = c(pond = "weighted", fenetre = "window", nombre = "count", puis = "power",
                         chap3 = "ch3", naif = "naive"),
                dilution = c(masse = "mass", nombre = "count"))
to_english <- function(z) {
  nm <- ifelse(names(z) %in% names(fr_cols), fr_cols[names(z)], names(z))
  names(z) <- sub("^prec_", "acc_", nm)
  for (v in intersect(names(fr_vals), names(z)))
    z[[v]] <- ifelse(z[[v]] %in% names(fr_vals[[v]]), fr_vals[[v]][z[[v]]], z[[v]])
  z
}
x <- do.call(rbind, lapply(f, function(g) to_english(read.csv(g))))
threshold <- qchisq(0.95, 1) / 2
prof1 <- function(xv, ll) {
  sp <- splinefun(xv, ll, method = "natural")
  xx <- seq(min(xv), max(xv), length.out = 2001); yy <- sp(xx)
  ci <- range(xx[yy >= max(yy) - threshold])
  c(est = xx[which.max(yy)], low = ci[1], high = ci[2])
}
mean_sp <- function(z, e) (z[[paste0("acc_", e, "_A")]] + z[[paste0("acc_", e, "_B")]]) / 2
sink(file.path(dir, "summary.txt"), split = TRUE)
cat(length(f), "data sets,", nrow(x), "fits,", sum(is.na(x$logLik)), "failure(s),",
    sum(!x$conv, na.rm = TRUE), "not converged\n")

# --- 1. profiles of lambda and d per data set ----------------------------------
est <- do.call(rbind, lapply(split(x, list(x$truth, x$seed)), function(z) {
  p <- z[z$grid == "weighted", ]
  M <- tapply(p$logLik, list(p$lambda, p$d), identity)
  pd <- prof1(as.numeric(colnames(M)), apply(M, 2, max))
  pl <- prof1(as.numeric(rownames(M)), apply(M, 1, max))
  nb <- z[z$grid == "count", ]; pn <- prof1(nb$d, nb$logLik)
  r8 <- p$logLik[p$lambda == 1 & p$d == 0.5]; fe <- z[z$grid == "window", ]
  data.frame(truth = z$truth[1], seed = z$seed[1],
             l = pl[1], l_low = pl[2], l_high = pl[3], d = pd[1], d_low = pd[2], d_high = pd[3],
             dn = pn[1], dn_low = pn[2], dn_high = pn[3],
             dLL_R4 = fe$logLik[fe$r == 4] - r8, dLL_R6 = fe$logLik[fe$r == 6] - r8,
             dLL_R12 = fe$logLik[fe$r == 12] - r8)
}))
utils::write.csv(est, file.path(dir, "estimates.csv"), row.names = FALSE)
cat("\n1. Profiles (exponential + mass, R = 8). Truth A: lambda = 1, d = 0.5.\n")
for (v in c("A", "B")) {
  e <- est[est$truth == v, ]
  cat(sprintf("\ntruth %s (%d data sets)\n", v, nrow(e)))
  cat(sprintf("  lambda: mean %.3f, RMSE %.3f, CI width %.3f, coverage of 1 %.2f\n",
              mean(e$l), sqrt(mean((e$l - 1)^2)), mean(e$l_high - e$l_low), mean(e$l_low <= 1 & e$l_high >= 1)))
  cat(sprintf("  d mass: mean %.3f, RMSE %.3f, CI width %.3f, coverage of 0.5 %.2f, d = 0 excluded %.2f, d at the edge %.2f\n",
              mean(e$d), sqrt(mean((e$d - 0.5)^2)), mean(e$d_high - e$d_low),
              mean(e$d_low <= 0.5 & e$d_high >= 0.5), mean(e$d_low > 0), mean(e$d %in% c(0, 1))))
  cat(sprintf("  d count (lambda = 1): mean %.3f, CI width %.3f, d = 0 excluded %.2f, d at the edge %.2f\n",
              mean(e$dn), mean(e$dn_high - e$dn_low), mean(e$dn_low > 0), mean(e$dn %in% c(0, 1))))
  cat(sprintf("  window, logLik - logLik(R = 8): R = 4 med %.3f max|.| %.3f; R = 6 med %.3f max|.| %.3f; R = 12 med %.3f max|.| %.3f\n",
              median(e$dLL_R4), max(abs(e$dLL_R4)), median(e$dLL_R6), max(abs(e$dLL_R6)),
              median(e$dLL_R12), max(abs(e$dLL_R12))))
}

# --- 2. compared models ----------------------------------------------------------
cmp <- do.call(rbind, lapply(split(x, list(x$truth, x$seed)), function(z) {
  p <- z[z$grid == "weighted", ]
  true_k <- if (z$truth[1] == "A") p[p$lambda == 1 & p$d == 0.5, ] else z[z$grid == "power", ]
  m <- list("true kernel" = true_k, "exp + mass, best point" = p[which.max(p$logLik), ],
            "exp + mass, lambda 1 d 0.5" = p[p$lambda == 1 & p$d == 0.5, ],
            "power r 3, count" = z[z$grid == "power", ],
            "chapter 3" = z[z$grid == "ch3", ], "naive" = z[z$grid == "naive", ])
  do.call(rbind, lapply(names(m), function(k) data.frame(
    truth = z$truth[1], seed = z$seed[1], model = k, dLL = m[[k]]$logLik - true_k$logLik,
    dge = mean_sp(m[[k]], "dge"), intra = mean_sp(m[[k]], "intra"), inter = mean_sp(m[[k]], "inter"),
    tbv = mean_sp(m[[k]], "tbv"),
    err_intra = (m[[k]]$err_intra_A + m[[k]]$err_intra_B) / 2,
    err_inter = (m[[k]]$err_inter_A + m[[k]]$err_inter_B) / 2)))
}))
utils::write.csv(cmp, file.path(dir, "models.csv"), row.names = FALSE)
cat("\n2. Compared models (medians; accuracy = mean of the two species; dLL = logLik - logLik of the true kernel)\n")
ag <- aggregate(cbind(dLL, dge, intra, inter, tbv, err_intra, err_inter) ~ truth + model, cmp, median)
ag <- ag[order(ag$truth, -ag$dLL), ]
print(cbind(ag[, 1:2], round(ag[, -(1:2)], 3)), row.names = FALSE)
cat("\nShare of data sets where the model beats the true kernel in logLik, and TBV gap to the true kernel (median, quartiles)\n")
for (v in c("A", "B")) for (k in setdiff(unique(cmp$model), "true kernel")) {
  z <- cmp[cmp$truth == v & cmp$model == k, ]; tv <- cmp[cmp$truth == v & cmp$model == "true kernel", ]
  dt <- z$tbv - tv$tbv[match(z$seed, tv$seed)]
  cat(sprintf("  %s | %-28s beats the truth %.2f; dTBV %.3f [%.3f; %.3f]\n", v, k, mean(z$dLL > 0),
              median(dt), quantile(dt, .25), quantile(dt, .75)))
}
sink()
