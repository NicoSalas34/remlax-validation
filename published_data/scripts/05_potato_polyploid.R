# ==============================================================================
# 05_potato_polyploid.R -- tetraploid genomic BLUP on a potato panel
# ------------------------------------------------------------------------------
# Data: the autotetraploid potato panel distributed with sommer (now in the
# enhancer package as DT_polyploid, GT_polyploid, MP_polyploid): 221 clones
# genotyped at 3522 SNPs with allele calls such as "AAAG", 187 phenotyped for
# yield, tuber and processing traits. The panel is the one of the GWASpoly
# example (Rosyara et al. 2016, Plant Genome 9, doi:10.3835/plantgenome2015.08.0073).
#
# The SNP calls are turned into dosages 0..4 of the second allele (in
# alphabetical order) of each marker. Missing calls are imputed by the marker
# mean rounded to an integer dosage, and monomorphic markers dropped, once, so that every builder sees the
# same matrix.
#
# Checks:
#   * the tetraploid VanRaden matrix of rx_grm(ploidy = 4) against
#     AGHmatrix::Gmatrix(method = "VanRaden", ploidy = 4,
#     ploidy.correction = TRUE), element by element;
#   * GBLUP (y ~ 1 + G) on three traits against asreml and sommer::mmes
#     (components, log-likelihood, BLUPs).
#
# No variance component is printed for these models in the source papers, so
# there are no published checks here.
#
#   Rscript published_data/scripts/05_potato_polyploid.R
# ==============================================================================
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(.here, "..", "R", "common.R"))
suppressPackageStartupMessages({ library(asreml); library(sommer); library(enhancer); library(AGHmatrix) })

start_checks("potato_tetraploid", "Rosyara et al. 2016, doi:10.3835/plantgenome2015.08.0073")
data(DT_polyploid)

# --- dosages 0..4 --------------------------------------------------------------
calls <- as.matrix(GT_polyploid)
dose_of <- function(x) {
  al <- sort(unique(unlist(strsplit(x[!is.na(x)], ""))))
  if (length(al) < 2) return(rep(NA_real_, length(x)))   # monomorphic
  b <- al[2]
  vapply(x, function(s) if (is.na(s)) NA_real_ else sum(strsplit(s, "")[[1]] == b), 0)
}
D <- apply(calls, 2, dose_of)
dimnames(D) <- dimnames(calls)
poly <- colSums(!is.na(D)) > 0 & apply(D, 2, function(v) length(unique(v[!is.na(v)])) > 1)
D <- D[, poly]
na_rate <- mean(is.na(D))
# Imputed by the ROUNDED marker mean: AGHmatrix accepts integer dosages only,
# and both builders must receive the same matrix.
for (j in seq_len(ncol(D))) { v <- D[, j]; v[is.na(v)] <- round(mean(v, na.rm = TRUE)); D[, j] <- v }
poly2 <- apply(D, 2, function(v) length(unique(v)) > 1)
D <- D[, poly2]
cat(sprintf("clones %d, SNPs kept %d of %d, missing calls imputed %.2f%%\n",
            nrow(D), ncol(D), ncol(calls), 100 * na_rate))

# --- the two tetraploid G ------------------------------------------------------
G_rx <- rx_grm(D, ploidy = 4, coding = "count")
# rx_grm divides by sum(ploidy p (1 - p)), the polyploid VanRaden scaling.
# AGHmatrix uses it with ploidy.correction = TRUE; its default divides by the
# sum of the sample variances of the centred dosages instead, a different
# convention, reported as context only.
G_agh <- Gmatrix(D, method = "VanRaden", ploidy = 4, missingValue = NA, maf = 0,
                 ploidy.correction = TRUE)
G_agh <- G_agh[rownames(G_rx), colnames(G_rx)]
add_check("G matrix", "max |G_rx - G_AGHmatrix|", "AGHmatrix",
          max(abs(G_rx - G_agh)), NA, tol = 1e-10, metric = "value",
          note = "ploidy.correction = TRUE")
add_check("G matrix", "mean diagonal", "AGHmatrix", mean(diag(G_rx)), mean(diag(G_agh)),
          note = "ploidy.correction = TRUE")
G_agh0 <- Gmatrix(D, method = "VanRaden", ploidy = 4, missingValue = NA, maf = 0)
add_check("G matrix", "mean diagonal", "AGHmatrix default", mean(diag(G_rx)),
          mean(diag(G_agh0)), kind = "context",
          note = "default scaling: sum of sample variances of the centred dosages")

G <- G_rx + diag(1e-4, nrow(G_rx))
dimnames(G) <- dimnames(G_rx)
G_as <- as_asreml_G(G)
WS <- "2gb"

vcr <- function(vc, row) {
  if (!row %in% rownames(vc)) stop("no varcomp row ", row, " in: ", paste(rownames(vc), collapse = ", "))
  vc[row, "component"]
}
asr_blup <- function(a, prefix) {
  cr <- summary(a, coef = TRUE)$coef.random
  k <- startsWith(rownames(cr), prefix)
  setNames(cr[k, "solution"], sub(prefix, "", rownames(cr)[k], fixed = TRUE))
}
timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; attr(v, "secs") <- proc.time()[["elapsed"]] - t0; v }
TIMES <- list()

ph <- DT_polyploid
ph <- ph[ph$Name %in% rownames(G), ]
cat(sprintf("phenotyped clones with genotypes: %d\n", nrow(ph)))
for (tr in c("total_yield", "chip_color", "tuber_shape")) {
  d <- data.frame(id = factor(ph$Name, levels = rownames(G)), y = ph[[tr]])
  d <- d[!is.na(d$y), ]
  m <- paste0(tr, " gblup")
  cat(sprintf("\n=== %s: %d records ===\n", m, nrow(d)))
  f <- timed(rx_cpu(y ~ 1, random = ~ vm(id, K = G), data = d))
  print(f)
  a <- timed(stab(asreml(y ~ 1, random = ~ vm(id, G_as), data = d, trace = FALSE,
                         maxit = 100, workspace = WS)))
  vc <- summary(a)$varcomp; print(vc)
  so <- timed(mmes(y ~ 1, random = ~ vsm(ism(id), Gu = G), rcov = ~ units, data = d,
                   verbose = FALSE, dateWarning = FALSE, nIters = 200,
                   tolParConvLL = 1e-8, tolParConvNorm = 1e-8))
  TIMES[[m]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"), sommer = attr(so, "secs"))
  add_check(m, "sigma2 G", "asreml", f$sigmas$id[1, 1], vcr(vc, "vm(id, G_as)"))
  add_check(m, "sigma2 residual", "asreml", f$sigma_res[1, 1], vcr(vc, "units!R"))
  add_check(m, "sigma2 G", "sommer", f$sigmas$id[1, 1], so$theta[["vsm(ism(id), Gu = G)"]][1, 1])
  add_check(m, "sigma2 residual", "sommer", f$sigma_res[1, 1], so$theta[["units"]][1, 1])
  add_check(m, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik, tol = 1e-5, metric = "abs")
  u_rx <- setNames(f$blups$id[, 1], rownames(f$blups$id))
  add_blup_checks(m, "asreml", u_rx, asr_blup(a, "vm(id, G_as)_")[names(u_rx)])
  u_so <- so$uList[["vsm(ism(id), Gu = G)"]]
  add_blup_checks(m, "sommer", u_rx, setNames(u_so[, 1], rownames(u_so))[names(u_rx)])
}

cat("\nElapsed seconds (remlax includes process start-up):\n")
for (m in names(TIMES)) cat(sprintf("  %-26s %s\n", m,
  paste(sprintf("%s %.1f", names(TIMES[[m]]), TIMES[[m]]), collapse = " | ")))
write_checks()
