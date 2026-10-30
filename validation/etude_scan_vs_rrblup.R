# ==============================================================================
# Does the remlax scan reproduce rrBLUP SNP by SNP?
# ==============================================================================
# WHY AN EXTERNAL REFERENCE, AND WHY ON THIS CASE. `remlax.scan` accepts V
# matrices that GWAS packages cannot handle: two kernels, an AR1 residual,
# heterogeneous sections. There is therefore no external reference for the
# case of interest. There is one for the DEGENERATE case where both are meant
# to coincide:
#
#     one kernel, iid residual, one observation per genotype
#
# This is the EMMAX model. If both implementations agree there, the algebra of
# the scan is validated against proven third-party code. The rest is the same
# algebra with another V, already checked against the Wald test in
# tests/python/test_scan.py.
#
# WHY rrBLUP AND NOT statgenGWAS. Either would do. rrBLUP::GWAS with
# P3D = TRUE is the most cited EMMAX in plant breeding, and it is PURE R: it
# loads anywhere, without a compiled library. The script
# `etude_scan_vs_statgen.R`, next to this one, does the same with statgenGWAS
# for those who can load it.
#
# WHAT DIFFERS BETWEEN THE TWO, AND MUST THEREFORE BE RECONSTRUCTED. rrBLUP
# re-estimates the residual scale AT EACH SNP, on the residuals of the model
# that CONTAINS the SNP, and tests with F(1, n-2):
#
#     s2_j   = resid_j' Hinv resid_j / v2       v2 = n - 2
#     Fstat  = beta_j^2 / (s2_j Winv_pp)
#
# The scan keeps the scale of the NULL model and tests with a chi2 on 1 df.
# Section 3.5 of the design note documents this point: at fixed V, a SNP of
# large effect is tested under an overestimated residual variance, hence
# conservatively. The relation between the two is exact and is rebuilt without
# recomputing anything in dimension n:
#
#     Fstat_j = chi2_j * v2 / (y'Py - num_j^2/den_j)
#
# where num_j and den_j are the raw outputs of the scan and y'Py a scalar
# computed once. The comparison therefore covers THREE levels:
#   (A) beta_j           -> must agree to machine precision (GLS algebra)
#   (B) the rrBLUP score -> rebuilt from the scan, must agree as well
#                           (the convention gap is then fully explained)
#   (C) the size of the conservativeness, measured where it bites
#
# Run:  Rscript validation/etude_scan_vs_rrblup.R <output_dir>
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else file.path(tempdir(), "scan_vs_rrblup")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(rrBLUP))
cat("rrBLUP", as.character(packageVersion("rrBLUP")), "\n")

set.seed(20260907)
q <- 250L      # genotypes = observations: one per genotype, as in EMMAX
p <- 800L

# --- markers coded -1/0/1, the rrBLUP coding -----------------------------
frq <- runif(p, 0.10, 0.90)
G <- matrix(rbinom(q * p, 2, rep(frq, each = q)) - 1L, nrow = q)
G <- G[, apply(G, 2, stats::sd) > 0, drop = FALSE]
p <- ncol(G)
rownames(G) <- sprintf("g%03d", seq_len(q))
colnames(G) <- sprintf("snp%04d", seq_len(p))

K <- rrBLUP::A.mat(G)          # the rrBLUP kinship, not another one

# --- phenotype: polygenic background + ONE causal SNP of large effect -----
# The large effect is deliberate. The two conventions diverge most there, so
# the comparison is most informative there.
Lk <- t(chol(K + diag(1e-6, q)))
j0 <- sample.int(p, 1L)
y <- as.numeric(2.0 + 0.80 * scale(G[, j0])[, 1] +
                1.0 * (Lk %*% rnorm(q))[, 1] + rnorm(q, sd = 1.0))

geno <- data.frame(marker = colnames(G), chrom = 1L, pos = seq_len(p),
                   t(G), check.names = FALSE, stringsAsFactors = FALSE)
pheno <- data.frame(line = rownames(G), yield = y, stringsAsFactors = FALSE)

# --- rrBLUP: variance components once, then one test per SNP -----------
res <- rrBLUP::GWAS(pheno = pheno, geno = geno, K = K, n.PC = 0,
                    min.MAF = 0, P3D = TRUE, plot = FALSE)
# The components that P3D fixed: mixed.solve on the NULL model, exactly what
# GWAS() does internally before the loop.
ms <- rrBLUP::mixed.solve(y = y, K = K, X = matrix(1, q, 1), SE = FALSE)
cat(sprintf("components (mixed.solve): Vu = %.8g  Ve = %.8g\n", ms$Vu, ms$Ve))
cat(sprintf("causal SNP: %s  rrBLUP score = %.4f\n",
            colnames(G)[j0], res$yield[res$marker == colnames(G)[j0]]))

utils::write.csv(data.frame(snp = res$marker, score = res$yield),
                 file.path(out_dir, "rrblup_scores.csv"), row.names = FALSE)
utils::write.table(y, file.path(out_dir, "y.txt"), row.names = FALSE, col.names = FALSE)
utils::write.table(K, file.path(out_dir, "K.txt"), row.names = FALSE, col.names = FALSE)
utils::write.table(G, file.path(out_dir, "M.txt"), row.names = FALSE, col.names = FALSE)
writeLines(sprintf("%.17g", c(ms$Vu, ms$Ve, q, p, j0)),
           file.path(out_dir, "meta.txt"))   # Vu, Ve, q, p, index of the causal SNP
cat("written to", out_dir, "\n")
