# ==============================================================================
# 13_pig_cleveland2012.R -- genomic additive, dominance and epistatic variance
# on the PIC pig data
# ------------------------------------------------------------------------------
# Data: Cleveland, Hickey & Forni (2012, G3 2:429-435,
# doi:10.1534/g3.111.001453, CC BY 3.0), released by PIC for method
# comparison. 3534 genotyped pigs, 52,843 SNPs (imputed doses, some not
# integer), a pedigree of 6473 animals and five traits t1 to t5 already
# corrected for fixed effects. Downloaded from the Europe PMC supplementary
# bundle of PMC3337471 (file supp_2.4.429_FileS1.zip).
#
# Relationship matrices, built once and passed to every program:
#   G   VanRaden (2008) first method on the doses, SNPs with 0 < p < 1
#   D   genomic dominance of Vitezica et al. (2013, Genetics 195:1223), HWE
#       coding (-2q^2, 2pq, -2p^2 for 2, 1, 0 copies of the counted allele),
#       on the doses rounded to the nearest integer
#   AA  G # G (Hadamard), divided by its mean diagonal
#   A   pedigree relationship (nadiv), for the pedigree BLUP
#
# Models, for each trait, y = mu + random terms + e:
#   GBLUP     g ~ N(0, s2A G)
#   A+D       g + d, d ~ N(0, s2D D)
#   A+D+AA    g + d + aa, aa ~ N(0, s2AA AA)
#   ABLUP     pedigree a ~ N(0, s2A A), on the full pedigree (A^-1, MME engine)
#   3 traits  t3, t4, t5: us(trait) on G, us(trait) residual
#
# References, same model and data: asreml (all models), gaston lmm.aireml
# (GBLUP, A+D, A+D+AA), sommer::mmes (GBLUP, A+D), rrBLUP::mixed.solve
# (GBLUP). Published genomic h2 and d2 of Da et al. (2014, PLoS ONE
# 9:e87666, Table 1) come from other SNP filters and another dominance
# coding: they are context rows.
#
#   Rscript published_data/scripts/13_pig_cleveland2012.R
# ==============================================================================
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(.here, "..", "R", "common.R"))
suppressPackageStartupMessages({
  library(asreml); library(sommer); library(gaston); library(rrBLUP); library(nadiv)
})

start_checks("pig_cleveland2012", "Cleveland et al. 2012, doi:10.1534/g3.111.001453")
DA2014 <- "Da et al. 2014 Table 1 (GVCBLUP, other SNP filters and dominance coding)"

cache <- pd_cache("pic_cleveland2012")
zip <- file.path(cache, "supp_2.4.429_FileS1.zip")
if (!file.exists(file.path(cache, "FileS1", "phenotypes.txt"))) {
  if (!file.exists(zip)) {
    bundle <- file.path(cache, "epmc_supp.zip")
    utils::download.file("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3337471/supplementaryFiles",
                         bundle, mode = "wb", quiet = TRUE)
    utils::unzip(bundle, files = "supp_2.4.429_FileS1.zip", exdir = cache)
  }
  utils::unzip(zip, exdir = cache)
}
gfile <- file.path(cache, "genotypes.rds")
if (!file.exists(gfile)) {
  g <- data.table::fread(file.path(cache, "FileS1", "genotypes.txt"), header = TRUE, showProgress = FALSE)
  ids <- as.character(g$ID); g[, ID := NULL]
  M <- as.matrix(g); rm(g); rownames(M) <- ids
  saveRDS(M, gfile, compress = "xz")
}
M <- readRDS(gfile)
ph <- utils::read.csv(file.path(cache, "FileS1", "phenotypes.txt"), na.strings = ".")
ped <- utils::read.csv(file.path(cache, "FileS1", "pedigree.txt"))
ph$ID <- as.character(ph$ID)
stopifnot(all(rownames(M) %in% ph$ID) || all(ph$ID %in% rownames(M)))
cat(sprintf("genotyped %d, SNPs %d, pedigree %d\n", nrow(M), ncol(M), nrow(ped)))

# ---- relationship matrices --------------------------------------------------
p <- colMeans(M) / 2
keep <- p > 0 & p < 1
M <- M[, keep]; p <- p[keep]; q <- 1 - p
W <- sweep(M, 2, 2 * p)
G <- tcrossprod(W) / (2 * sum(p * q)); rm(W)
Mr <- round(M)
H <- matrix(0, nrow(Mr), ncol(Mr))
H[Mr == 2] <- (-2 * q^2)[col(Mr)[Mr == 2]]
H[Mr == 1] <- (2 * p * q)[col(Mr)[Mr == 1]]
H[Mr == 0] <- (-2 * p^2)[col(Mr)[Mr == 0]]
D <- tcrossprod(H) / sum((2 * p * q)^2); rm(H, Mr)
AA <- G * G; AA <- AA / mean(diag(AA))
ridge <- function(K) { K <- K + diag(1e-4, nrow(K)); dimnames(K) <- list(rownames(M), rownames(M)); K }
G <- ridge(G); D <- ridge(D); AA <- ridge(AA)
rm(M); invisible(gc())
cat(sprintf("SNPs kept %d; mean diag G %.3f, D %.3f, AA %.3f\n", sum(keep),
            mean(diag(G)), mean(diag(D)), mean(diag(AA))))

ped[ped == 0] <- NA
pp <- prepPed(ped)
Ai <- makeAinv(pp)$Ainv
dimnames(Ai) <- list(as.character(pp$ID), as.character(pp$ID))
ginv <- function(Mx) {
  T <- methods::as(methods::as(Matrix::forceSymmetric(Mx, "L"), "generalMatrix"), "TsparseMatrix")
  k <- T@i >= T@j
  g <- cbind(Row = T@i[k] + 1L, Column = T@j[k] + 1L, Value = T@x[k])
  g <- g[order(g[, "Row"], g[, "Column"]), ]
  attr(g, "rowNames") <- rownames(Mx); attr(g, "INVERSE") <- TRUE
  g
}
gA <- ginv(Ai)

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
da2014 <- list(t1 = c(0.03, 7.2e-7), t2 = c(0.27, 0.02), t3 = c(0.22, 0.07),
               t4 = c(0.35, 0.01), t5 = c(0.38, 0.05))

for (tr in paste0("t", 1:5)) {
  dt <- ph[!is.na(ph[[tr]]) & ph$ID %in% rownames(G), c("ID", tr)]
  names(dt)[2] <- "y"
  ids <- dt$ID
  dt$g <- factor(ids, levels = ids); dt$d <- dt$g; dt$aa <- dt$g
  Gs <- G[ids, ids]; Ds <- D[ids, ids]; AAs <- AA[ids, ids]
  Gs_as <- as_asreml_G(Gs); Ds_as <- as_asreml_G(Ds); AAs_as <- as_asreml_G(AAs)
  cat(sprintf("\n######## %s: %d genotyped pigs with a record\n", tr, nrow(dt)))

  mods <- list(
    GBLUP = list(rx = ~ vm(g, K = Gs), as = ~ vm(g, Gs_as), K = list(Gs),
                 so = ~ vsm(ism(g), Gu = Gs),
                 comp = c(A = "g"), asn = c(A = "vm(g, Gs_as)"),
                 son = c(A = "vsm(ism(g), Gu = Gs)")),
    "A+D" = list(rx = ~ vm(g, K = Gs) + vm(d, K = Ds), as = ~ vm(g, Gs_as) + vm(d, Ds_as),
                 K = list(Gs, Ds), so = ~ vsm(ism(g), Gu = Gs) + vsm(ism(d), Gu = Ds),
                 comp = c(A = "g", D = "d"), asn = c(A = "vm(g, Gs_as)", D = "vm(d, Ds_as)"),
                 son = c(A = "vsm(ism(g), Gu = Gs)", D = "vsm(ism(d), Gu = Ds)")),
    "A+D+AA" = list(rx = ~ vm(g, K = Gs) + vm(d, K = Ds) + vm(aa, K = AAs),
                    as = ~ vm(g, Gs_as) + vm(d, Ds_as) + vm(aa, AAs_as), K = list(Gs, Ds, AAs),
                    comp = c(A = "g", D = "d", AA = "aa"),
                    asn = c(A = "vm(g, Gs_as)", D = "vm(d, Ds_as)", AA = "vm(aa, AAs_as)")))
  for (m in names(mods)) {
    s <- mods[[m]]; lab <- paste(tr, m)
    cat(sprintf("\n=== %s ===\n", lab))
    f <- timed(rx_cpu(y ~ 1, random = s$rx, data = dt))
    print(f)
    a <- timed(stab(asreml(y ~ 1, random = s$as, data = dt, trace = FALSE, maxit = 100,
                           workspace = "4gb")))
    vc <- summary(a)$varcomp; print(vc)
    ga <- timed(lmm.aireml(dt$y, X = matrix(1, nrow(dt)), K = s$K, eps = 1e-10,
                           max_iter = 300, min_tau = 1e-10, verbose = FALSE))
    TIMES[[lab]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"), gaston = attr(ga, "secs"))
    for (k in seq_along(s$comp)) {
      cn <- names(s$comp)[k]; v <- f$sigmas[[s$comp[[k]]]][1, 1]
      add_check(lab, paste0("sigma2 ", cn), "asreml", v, vcr(vc, s$asn[[cn]]))
      add_check(lab, paste0("sigma2 ", cn), "gaston", v, ga$tau[k])
    }
    add_check(lab, "sigma2 residual", "asreml", f$sigma_res[1, 1], vcr(vc, "units!R"))
    add_check(lab, "sigma2 residual", "gaston", f$sigma_res[1, 1], ga$sigma2)
    add_check(lab, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
              tol = 1e-5, metric = "abs")
    add_check(lab, "logLik (asreml convention)", "gaston", f$logLik_asreml, ga$logL,
              tol = 1e-5, metric = "abs")
    u_rx <- setNames(f$blups$g[, 1], rownames(f$blups$g))
    add_blup_checks(paste(lab, "A"), "asreml", u_rx, asr_blup(a, "vm(g, Gs_as)_")[names(u_rx)])
    if (!is.null(s$so)) {
      so <- timed(mmes(y ~ 1, random = s$so, rcov = ~ units, data = dt, verbose = FALSE,
                       dateWarning = FALSE))
      TIMES[[lab]] <- c(TIMES[[lab]], sommer = attr(so, "secs"))
      for (k in seq_along(s$comp)) {
        cn <- names(s$comp)[k]; j <- match(s$son[[cn]], names(so$theta))
        add_check(lab, paste0("sigma2 ", cn), "sommer", f$sigmas[[s$comp[[k]]]][1, 1],
                  if (is.na(j)) NA else as.numeric(so$theta[[j]]))
      }
      add_check(lab, "sigma2 residual", "sommer", f$sigma_res[1, 1],
                as.numeric(so$theta[[match("units", names(so$theta))]]))
    }
    if (m == "GBLUP") {
      ms <- mixed.solve(dt$y, K = Gs, method = "REML")
      add_check(lab, "sigma2 A", "rrBLUP", f$sigmas$g[1, 1], ms$Vu)
      add_check(lab, "sigma2 residual", "rrBLUP", f$sigma_res[1, 1], ms$Ve)
    }
    if (m == "A+D") {
      tot <- f$sigmas$g[1, 1] + f$sigmas$d[1, 1] + f$sigma_res[1, 1]
      add_check(lab, "h2 (A share)", "Da 2014", f$sigmas$g[1, 1] / tot, da2014[[tr]][1],
                kind = "context", note = DA2014)
      add_check(lab, "d2 (D share)", "Da 2014", f$sigmas$d[1, 1] / tot, da2014[[tr]][2],
                kind = "context", note = DA2014)
    }
  }

  # Pedigree BLUP on the full pedigree, all animals with a record.
  lab <- paste(tr, "ABLUP")
  da <- ph[!is.na(ph[[tr]]), c("ID", tr)]; names(da)[2] <- "y"
  da$a <- factor(da$ID, levels = rownames(Ai))
  cat(sprintf("\n=== %s (%d records) ===\n", lab, nrow(da)))
  f <- timed(rx_cpu(y ~ 1, random = ~ vm(a, Kinv = Ai), data = da, engine = "mme"))
  print(f)
  a <- timed(stab(asreml(y ~ 1, random = ~ vm(a, gA), data = da, trace = FALSE, maxit = 100,
                         workspace = "2gb")))
  vc <- summary(a)$varcomp; print(vc)
  TIMES[[lab]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"))
  add_check(lab, "sigma2 A", "asreml", f$sigmas$a[1, 1], vcr(vc, "vm(a, gA)"))
  add_check(lab, "sigma2 residual", "asreml", f$sigma_res[1, 1], vcr(vc, "units!R"))
  add_check(lab, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
            tol = 1e-5, metric = "abs")
  u_rx <- setNames(f$blups$a[, 1], rownames(f$blups$a))
  add_blup_checks(paste(lab, "A"), "asreml", u_rx, asr_blup(a, "vm(a, gA)_")[names(u_rx)])
}

# ---- three traits, us(trait) on G and on the residual ------------------------
mname <- "t3+t4+t5 us G"
w <- ph[ph$ID %in% rownames(G) & stats::complete.cases(ph[, c("t3", "t4", "t5")]), ]
ids <- w$ID; Gs <- G[ids, ids]; Gs_as <- as_asreml_G(Gs)
w$g <- factor(ids, levels = ids)
dl <- do.call(rbind, lapply(c("t3", "t4", "t5"), function(t)
  data.frame(g = w$g, car = t, y = w[[t]])))
dl$car <- factor(dl$car, levels = c("t3", "t4", "t5"))
dl <- dl[order(dl$g, dl$car), ]
cat(sprintf("\n=== %s (%d pigs) ===\n", mname, length(ids)))
f <- timed(rx_cpu(y ~ car, random = ~ us(g, K = Gs), residual = ~ us(car):units,
                  data = dl, trait = "car", unit = "g"))
print(f)
a <- timed(stab(asreml(cbind(t3, t4, t5) ~ trait, random = ~ us(trait):vm(g, Gs_as),
                       residual = ~ units:us(trait), data = w, trace = FALSE, maxit = 100,
                       workspace = "4gb")))
vc <- summary(a)$varcomp; print(vc)
TIMES[[mname]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"))
tl <- c("t3", "t4", "t5")
get_us <- function(term) {
  S <- matrix(NA_real_, 3, 3)
  for (i in 1:3) for (j in 1:i) {
    k <- match(sprintf("%s!trait_%s:%s", term, tl[i], tl[j]), rownames(vc))
    S[i, j] <- S[j, i] <- vc[k[1], "component"]
  }
  S
}
for (blk in list(list("G", f$sigmas$g, get_us("trait:vm(g, Gs_as)")),
                 list("residual", f$sigma_res, get_us("units:trait")))) {
  Srx <- blk[[2]]; Sas <- blk[[3]]
  for (i in 1:3) add_check(mname, sprintf("%s var %s", blk[[1]], tl[i]), "asreml", Srx[i, i], Sas[i, i])
  for (i in 2:3) for (j in 1:(i - 1))
    add_check(mname, sprintf("%s cov %s-%s (/ sqrt(var var))", blk[[1]], tl[i], tl[j]), "asreml",
              abs(Srx[i, j] - Sas[i, j]) / sqrt(Sas[i, i] * Sas[j, j]), NA, tol = 1e-3,
              metric = "value")
}
add_check(mname, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
          tol = 1e-5, metric = "abs")

cat("\nRuntimes (s, elapsed):\n")
print(TIMES)
write_checks()
