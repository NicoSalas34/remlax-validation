# ==============================================================================
# 15_drops_multitrait.R -- multi-trait GWAS on the DROPS maize panel: grain
#                          yield in three trials as three traits, against asreml
# ------------------------------------------------------------------------------
# DATA. As in 07_drops_gwas.R: the DROPS panel (Millet et al. 2016, Plant
# Physiol 172:749-764; Millet et al. 2019, Sci Data 6:190011), shipped with
# statgenGWAS. Grain yield BLUEs of the hybrids in the trials Mur13W, Mur13R
# and Cra12R (genetic correlations from about 0.3 to 0.9), long format, one
# row per hybrid and trial; hybrids missing a trial keep their other records.
#
# NULL MODEL. y = trial + g + e, g ~ N(0, Sigma_g (x) K) with Sigma_g us
# between trials, e ~ N(0, diag(s2_e,trial)), K the VanRaden kinship of
# statgenGWAS plus a ridge of 1e-4 (as in 07).
#
# SCAN (rx_scan, fixed V). Tests "dir" (any effect, 3 df), "dir@common"
# (1 df), "dir@het" (2 df) and "dir[Cra12R]" (1 df).
#
# REFERENCES, same data and model, asreml:
#   (A) variance parameters FIXED at the remlax estimates (G.param and
#       R.param with constraint F), one fit per SNP with the SNP as a fixed
#       effect: conditional Wald statistics of
#           env:snp           in y ~ env + env:snp           (3 df)  = dir
#           snp               in y ~ env + snp               (1 df)  = dir@common
#           env:snp | snp     in y ~ env + snp + env:snp     (2 df)  = dir@het
#       on the 5 top SNPs of the scan and 5 random SNPs. At fixed V these are
#       the same statistics: tolerance 1e-4 relative.
#   (B) variance parameters re-estimated (REML refit with the SNP), against
#       rx_scan_exact() on the 3 top SNPs: tolerance 1e-3 relative.
#   The identity chi2(dir) = chi2(dir@common) + chi2(dir@het) is checked on
#   every SNP.
#
#   Rscript published_data/scripts/15_drops_multitrait.R
# ==============================================================================
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(f[1]), "..", "R", "common.R"))
})
suppressPackageStartupMessages({ library(statgenGWAS); library(asreml) })

start_checks("maize_drops_multitrait",
             "Millet et al. 2016 Plant Physiol; Millet et al. 2019 Sci Data")
TRIALS <- c("Mur13W", "Mur13R", "Cra12R")
RIDGE <- 1e-4; MAF <- 0.05

data(dropsMarkers); data(dropsMap); data(dropsPheno)
M <- as.matrix(dropsMarkers[, -1]); rownames(M) <- dropsMarkers$Ind
storage.mode(M) <- "double"
map <- data.frame(snp = dropsMap$SNP.names, chr = dropsMap$Chromosome, pos = dropsMap$Position)
map <- map[match(colnames(M), map$snp), ]
gd <- createGData(geno = dropsMarkers[, -1, drop = FALSE] |> `rownames<-`(dropsMarkers$Ind),
                  map = data.frame(chr = dropsMap$Chromosome, pos = dropsMap$Position,
                                   row.names = dropsMap$SNP.names))
K <- kinship(as.matrix(gd$markers), method = "vanRaden") + diag(RIDGE, nrow(M))
dimnames(K) <- list(rownames(M), rownames(M))

ph <- dropsPheno[dropsPheno$Experiment %in% TRIALS, c("Experiment", "Variety_ID", "grain.yield")]
ph <- ph[!is.na(ph$grain.yield) & ph$Variety_ID %in% rownames(M), ]
d <- data.frame(env = factor(ph$Experiment, levels = TRIALS),
                gen = factor(as.character(ph$Variety_ID), levels = rownames(M)),
                y = ph$grain.yield)
d <- d[order(d$env, d$gen), ]   # asreml dsum(~ units | env) needs rows grouped by env
d$unit <- factor(paste(d$gen, d$env))
cat(sprintf("hybrids %d, records %d (%s)\n", nlevels(droplevels(d$gen)), nrow(d),
            paste(table(d$env), collapse = "/")))

# ---- null model and scan -------------------------------------------------------
f <- rx_cpu(y ~ env, random = ~ us(gen, K = K), residual = ~ diag(env):units,
            data = d, trait = "env", unit = "unit")
print(f)
keep_snp <- { p <- colMeans(M) / 2; pmin(p, 1 - p) >= MAF }
Mk <- M[, keep_snp]; mapk <- map[keep_snp, ]
tests <- c("dir", "dir@common", "dir@het", "dir[Cra12R]")
t0 <- proc.time()[3]
res <- rx_scan(f, f$model, markers = Mk, map = mapk, incidences = c(dir = "gen"),
               tests = tests, maf = NULL, coding = "012", verbose = FALSE,
               backend = Sys.getenv("PD_BACKEND", "cpu"))
t_scan <- proc.time()[3] - t0
cat(sprintf("scan of %d SNPs in %.1f s\n", nrow(res), t_scan))
add_check("scan", "max |chi2(dir) - common - het| / chi2(dir)", "remlax",
          max(abs(res$chi2_dir - res[["chi2_dir@common"]] - res[["chi2_dir@het"]]) /
                pmax(res$chi2_dir, 1e-8), na.rm = TRUE), NA, tol = 1e-8, metric = "value")

# ---- (A) asreml at the remlax variance parameters --------------------------------
K_as <- as_asreml_G(K)
Sg <- f$sigmas$gen; Re <- diag(f$sigma_res)
sv_fixed <- function(call) {
  sv <- eval(as.call(c(as.list(call), start.values = TRUE)))$vparameters.table
  rownames(sv) <- sv$Component
  for (i in seq_along(TRIALS)) for (j in seq_len(i)) {
    nm <- sprintf("env:vm(gen, K_as)!env_%s:%s", TRIALS[i], TRIALS[j])
    if (nm %in% rownames(sv)) sv[nm, "Value"] <- Sg[i, j]
  }
  for (i in seq_along(TRIALS)) {
    nm <- sprintf("env_%s!R", TRIALS[i])
    if (nm %in% rownames(sv)) sv[nm, "Value"] <- Re[i]
  }
  sv$Constraint[sv$Constraint %in% c("P", "U")] <- "F"
  sv
}
# Wald statistic of the coefficients whose names match `rx`, in full precision
# from the coefficients and their covariance (asreml.options(Cfixed = TRUE)):
# the table of wald() rounds the statistic to four or five digits. Aliased
# coefficients (zero row of Cfixed) are left out, so that env:snp given snp
# keeps its T - 1 estimable contrasts.
asreml.options(Cfixed = TRUE)
wald_term <- function(a, rx) {
  C <- as.matrix(a$Cfixed); nm <- rownames(C); cf <- a$coefficients$fixed[nm, 1]
  i <- grep(rx, nm); i <- i[diag(C)[i] > 1e-12]
  b <- cf[i]
  drop(t(b) %*% solve(C[i, i, drop = FALSE], b))
}
ord <- order(res$p_dir)
set.seed(1)
pick <- unique(c(ord[1:5], sample(which(is.finite(res$p_dir)), 5)))
for (i in pick) {
  snp <- res$snp[i]
  d$snp <- Mk[as.character(d$gen), snp]
  calls <- list(
    dir = quote(asreml(y ~ env + env:snp, random = ~ us(env):vm(gen, K_as),
                       residual = ~ dsum(~ units | env), data = d, trace = FALSE)),
    common = quote(asreml(y ~ env + snp, random = ~ us(env):vm(gen, K_as),
                          residual = ~ dsum(~ units | env), data = d, trace = FALSE)),
    het = quote(asreml(y ~ env + snp + env:snp, random = ~ us(env):vm(gen, K_as),
                       residual = ~ dsum(~ units | env), data = d, trace = FALSE)))
  terms <- c(dir = ":snp$", common = "^snp$", het = ":snp$")
  cols <- c(dir = "chi2_dir", common = "chi2_dir@common", het = "chi2_dir@het")
  for (k in names(calls)) {
    sv <- sv_fixed(calls[[k]])
    a <- eval(as.call(c(as.list(calls[[k]]), G.param = quote(sv), R.param = quote(sv),
                        maxit = 1)))
    add_check(paste("A", snp), paste(cols[[k]], "(fixed V)"), "asreml",
              res[[cols[[k]]]][i], wald_term(a, terms[[k]]), tol = 1e-4)
  }
}

# ---- (B) re-estimated variance parameters: rx_scan_exact vs asreml refit ------
build <- function(Kx, X_sup = NULL) {
  Zl <- lapply(TRIALS, function(t) {
    Z <- Matrix::sparseMatrix(i = which(d$env == t), j = as.integer(d$gen[d$env == t]),
                              x = 1, dims = c(nrow(d), nlevels(d$gen)))
    Z })
  rx_model(d$y, cbind(stats::model.matrix(~ env, d), X_sup),
           list(rx_term("gen", Zl, K = Kx, struct = "us", levels = levels(d$gen))),
           rx_residual("diag", trait = d$env, unit = d$unit))
}
ex <- rx_scan_exact(res, "dir", build, function(W) K, Mk, mapk, c(dir = "gen"), top = 3L,
                    loco = FALSE, verbose = FALSE, fit_args = list(backend = "cpu"))
for (r in seq_len(nrow(ex))) {
  d$snp <- Mk[as.character(d$gen), ex$snp[r]]
  a <- stab(asreml(y ~ env + env:snp, random = ~ us(env):vm(gen, K_as),
                   residual = ~ dsum(~ units | env), data = d, trace = FALSE, maxit = 100))
  add_check(paste("B", ex$snp[r]), "chi2 dir (refit)", "asreml", ex$chi2_exact[r],
            wald_term(a, ":snp$"), tol = 1e-3)
  add_check(paste("B", ex$snp[r]), "p fixed V / p refit", "remlax",
            ex$p_fixed_V[r] / ex$p_exact[r], NA, kind = "context", metric = "value",
            note = "fixed V is conservative for large effects (Zhou & Stephens 2014)")
}
write_checks()
