# ==============================================================================
# Does the remlax scan reproduce statgenGWAS SNP by SNP?
# ==============================================================================
# WHY THIS TEST EXISTS, AND WHY IT IS THE ONLY ONE POSSIBLE. `remlax.scan`
# accepts V matrices that statgenGWAS cannot handle: two kernels, an AR1
# residual, heterogeneous sections. There is therefore NO external reference
# for the case of interest. There is one for the DEGENERATE case where both are
# meant to coincide exactly:
#
#     one kernel, iid residual, one observation per genotype
#
# This is the EMMAX model. If both implementations agree there, the algebra of
# the scan is validated against proven third-party code. The rest (two kernels,
# AR1) is the same algebra with another V, already checked against the Wald
# test in tests/python/test_scan.py.
#
# WHAT IS COMPARED, AND WHY NOT THE p-VALUES. The script compares `effect` and
# `effectSe`, the quantities both programs compute. The p-values depend on the
# REFERENCE DISTRIBUTION: statgenGWAS reports a Student test, remlax a chi2 with
# 1 df. Both are defensible at fixed V and differ only by the denominator
# degrees of freedom. Comparing p-values would mix an algebra disagreement
# (serious) with a distribution convention (harmless, and documented).
#
# THE VARIANCE COMPONENTS ARE FIXED. statgenGWAS estimates them by EMMA, remlax
# by REML. Nothing guarantees the same optimum to the 12th decimal, and a gap
# there would spread to every SNP. The test targets the algebra of the SCAN,
# not the estimation (tested elsewhere in remlax), so the components of
# statgenGWAS are imposed on the scan. The remaining disagreement can then only
# come from the scan.
#
# Run:  Rscript validation/etude_scan_vs_statgen.R <output_dir>
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else file.path(tempdir(), "scan_vs_statgen")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(statgenGWAS))
cat("statgenGWAS", as.character(packageVersion("statgenGWAS")), "\n")

set.seed(20260907)
q <- 250L      # genotypes = observations (one per genotype, as in EMMAX)
p <- 800L      # SNP

# --- genotypes as 0/1/2 dosages, varied frequencies ----------------------
frq <- runif(p, 0.10, 0.90)
G <- matrix(rbinom(q * p, 2, rep(frq, each = q)), nrow = q)
G <- G[, apply(G, 2, stats::sd) > 0, drop = FALSE]
p <- ncol(G)
rownames(G) <- sprintf("g%03d", seq_len(q))
colnames(G) <- sprintf("snp%04d", seq_len(p))

gmap <- data.frame(chr = rep(1L, p), pos = seq_len(p))
rownames(gmap) <- colnames(G)

# --- vanRaden kinship, the one statgenGWAS uses ---------------------------
K <- kinship(G, method = "vanRaden")

# --- phenotype: polygenic background + one causal SNP ---------------------
Lk <- t(chol(K + diag(1e-6, q)))
j0 <- sample.int(p, 1L)
y <- 2.0 + 0.45 * scale(G[, j0])[, 1] +
     1.0 * (Lk %*% rnorm(q))[, 1] + rnorm(q, sd = 1.0)
pheno <- data.frame(genotype = rownames(G), yield = y, stringsAsFactors = FALSE)

# --- statgenGWAS ---------------------------------------------------------
gd <- createGData(geno = G, map = gmap, pheno = pheno, kin = K)
res <- runSingleTraitGwas(gd, traits = "yield", kin = K,
                          GLSMethod = "single",     # no LOCO: a single V
                          remlAlgo = "EMMA",
                          MAF = 0,                  # no SNP dropped
                          thrType = "fixed", LODThr = 0)
gw <- as.data.frame(res$GWAResult$pheno)
vc <- res$GWASInfo$varComp
cat("variance components (EMMA):\n"); print(unlist(vc))

# --- what the scan must read back ---------------------------------------
# The dosages are taken AS statgenGWAS saw them (imputation and filters
# included). Comparing on other columns would compare nothing.
Gs <- gd$markers[rownames(G), gw$snp, drop = FALSE]
utils::write.csv(gw[, c("snp", "allFreq", "effect", "effectSe", "pValue", "LOD")],
                 file.path(out_dir, "statgen_results.csv"), row.names = FALSE)
utils::write.csv(data.frame(genotype = rownames(G), y = y),
                 file.path(out_dir, "pheno.csv"), row.names = FALSE)
utils::write.table(K, file.path(out_dir, "K.txt"), row.names = FALSE, col.names = FALSE)
utils::write.table(Gs, file.path(out_dir, "M.txt"), row.names = FALSE, col.names = FALSE)
writeLines(as.character(gw$snp), file.path(out_dir, "snp.txt"))
writeLines(jsonlite::toJSON(list(varComp = vc, q = q, p = length(gw$snp),
                                 causal_snp = colnames(G)[j0]),
                            auto_unbox = TRUE),
           file.path(out_dir, "meta.json"))
cat("written to", out_dir, ":", length(gw$snp), "SNPs,", q, "genotypes\n")
