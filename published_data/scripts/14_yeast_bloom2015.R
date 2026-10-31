# ==============================================================================
# 14_yeast_bloom2015.R -- additive and additive x additive (A # A) variance
# with strain repeatability, on the BYxRM yeast segregants
# ------------------------------------------------------------------------------
# Data: Bloom et al. (2015, Nat. Commun. 6:8712, doi:10.1038/ncomms9712,
# CC BY 4.0), Supplementary Data 5: 4390 haploid segregants of the BY x RM
# cross genotyped at 28,220 markers coded -1/+1, and 20 end-point growth
# traits with up to two replicates per segregant (pheno_raw).
#
# Model of the paper (Methods, "Variance component analysis"), fitted on the
# replicate records:
#   y = mu + Z a + Z i + Z p + e
#   a ~ N(0, A s2A)       A = rrBLUP::A.mat(markers) / 2 (fraction shared)
#   i ~ N(0, A#A s2AA)    Hadamard product of A with itself
#   p ~ N(0, I s2P)       strain effect beyond A and A#A (repeatability)
# The paper fits it by its own AI-REML code and reports, for each trait, the
# share of each component in the total variance (Supplementary Data 2,
# printed to 17 digits) for y = Zp + e, y = Za + Zp + e and
# y = Za + Zi + Zp + e. Those shares are compared here at 1e-3 absolute,
# the order of the stopping rule of an AI-REML iteration.
#
# Software reference, same model and data: asreml (two dense relationship
# matrices passed with INVERSE = FALSE).
#
# The four traits below are those quoted in the survey of candidate data
# (docs/candidate_datasets_animals_humans.md); the dense V of about 7000
# records makes each fit expensive on a CPU.
#
#   Rscript published_data/scripts/14_yeast_bloom2015.R [trait ...]
# ==============================================================================
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(.here, "..", "R", "common.R"))
suppressPackageStartupMessages({ library(asreml); library(rrBLUP) })

start_checks("yeast_bloom2015", "Bloom et al. 2015, doi:10.1038/ncomms9712")
PAPER <- "Bloom et al. 2015, Supplementary Data 2"

base <- "https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9712/MediaObjects/"
cache <- pd_cache("yeast_bloom2015")
zfile <- file.path(cache, "41467_2015_BFncomms9712_MOESM729_ESM.zip")
xfile <- file.path(cache, "41467_2015_BFncomms9712_MOESM726_ESM.xls")
if (!file.exists(xfile)) utils::download.file(paste0(base, basename(xfile)), xfile, mode = "wb", quiet = TRUE)
rd <- file.path(cache, "rdata")
if (!file.exists(file.path(rd, "phenotypes.RData"))) {
  if (!file.exists(zfile)) utils::download.file(paste0(base, basename(zfile)), zfile, mode = "wb", quiet = TRUE)
  tmp <- file.path(cache, "unz"); utils::unzip(zfile, exdir = tmp)
  for (z in list.files(tmp, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE))
    utils::unzip(z, exdir = rd)
}
load(file.path(rd, "genotypes.RData"))   # gdata: 4390 x 28220, -1/+1
load(file.path(rd, "phenotypes.RData"))  # pheno_raw: list of 20 traits
stopifnot(identical(rownames(gdata), names(pheno_raw[[1]])))

# Reference shares of Supplementary Data 2.
ref <- suppressMessages(readxl::read_excel(xfile, col_names = FALSE))
ref <- as.data.frame(ref)
row_of <- function(tr) which(ref[[1]] == tr)
num <- function(i, j) as.numeric(ref[i, j])

traits <- commandArgs(trailingOnly = TRUE)
if (!length(traits)) traits <- c("CobaltChloride", "CopperSulfate", "Hydroxyurea", "Zeocin")

# The paper's A is "the fraction of genome shared between pairs of
# segregants", and A # A "the fraction of pairs of markers shared". For
# haploid segregants coded -1/+1, rrBLUP::A.mat has a diagonal of 2: the
# fraction shared is A.mat / 2, and its Hadamard square has a unit diagonal.
# With these scalings the variance parameters are the shares of the paper.
A <- A.mat(gdata) / 2
AA <- A * A
cat(sprintf("segregants %d, markers %d; mean diag A %.4f, A#A %.4f\n", nrow(gdata), ncol(gdata),
            mean(diag(A)), mean(diag(AA))))
rm(gdata); invisible(gc())

vcr <- function(vc, row) {
  if (!row %in% rownames(vc)) stop("no varcomp row ", row, " in: ", paste(rownames(vc), collapse = ", "))
  vc[row, "component"]
}
timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; attr(v, "secs") <- proc.time()[["elapsed"]] - t0; v }
TIMES <- list()

for (tr in traits) {
  key <- grep(paste0("_", tr, "_"), names(pheno_raw), fixed = TRUE, value = TRUE)
  stopifnot(length(key) == 1L)
  pr <- pheno_raw[[key]]
  d <- do.call(rbind, lapply(names(pr), function(s) {
    v <- pr[[s]]; v <- v[!is.na(v)]
    if (length(v)) data.frame(strain = s, y = v) else NULL
  }))
  lev <- unique(d$strain)
  d$a <- factor(d$strain, levels = lev); d$i <- d$a; d$p <- d$a
  As <- A[lev, lev]; AAs <- AA[lev, lev]
  As_as <- as_asreml_G(As); AAs_as <- as_asreml_G(AAs)
  i0 <- row_of(tr)
  cat(sprintf("\n######## %s: %d records on %d segregants\n", tr, nrow(d), length(lev)))

  mods <- list(
    "P" = list(rx = ~ p, as = ~ p, comp = c(P = "p"), asn = c(P = "p"),
               ref = c(P = num(i0, 2), E = num(i0, 4))),
    "A+P" = list(rx = ~ vm(a, K = As) + p, as = ~ vm(a, As_as) + p,
                 comp = c(A = "a", P = "p"), asn = c(A = "vm(a, As_as)", P = "p"),
                 ref = c(A = num(i0, 6), P = num(i0, 8), E = num(i0, 10))),
    "A+AA+P" = list(rx = ~ vm(a, K = As) + vm(i, K = AAs) + p,
                    as = ~ vm(a, As_as) + vm(i, AAs_as) + p,
                    comp = c(A = "a", AA = "i", P = "p"),
                    asn = c(A = "vm(a, As_as)", AA = "vm(i, AAs_as)", P = "p"),
                    ref = c(A = num(i0, 12), AA = num(i0, 14), P = num(i0, 16), E = num(i0, 18))))
  for (m in names(mods)) {
    s <- mods[[m]]; lab <- paste(tr, m)
    cat(sprintf("\n=== %s ===\n", lab))
    f <- timed(rx_cpu(y ~ 1, random = s$rx, data = d))
    print(f)
    a <- timed(stab(asreml(y ~ 1, random = s$as, data = d, trace = FALSE, maxit = 100,
                           workspace = "5gb")))
    vc <- summary(a)$varcomp; print(vc)
    TIMES[[lab]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"))
    est <- c(sapply(s$comp, function(t) f$sigmas[[t]][1, 1]), E = f$sigma_res[1, 1])
    for (cn in names(s$comp))
      add_check(lab, paste0("sigma2 ", cn), "asreml", est[[cn]], vcr(vc, s$asn[[cn]]))
    add_check(lab, "sigma2 residual", "asreml", est[["E"]], vcr(vc, "units!R"))
    add_check(lab, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
              tol = 1e-5, metric = "abs")
    share <- est / sum(est)
    for (cn in names(s$ref))
      add_check(lab, paste0("share ", cn), "paper", share[[cn]], s$ref[[cn]], tol = 1e-3,
                metric = "abs", kind = "published", note = PAPER)
  }
}

cat("\nRuntimes (s, elapsed):\n")
print(round(do.call(rbind, TIMES), 1))
# One trait per job on the cluster: one file per trait, pooled by make_report.py.
write_checks(if (length(traits) == 1L)
  file.path(pd_root(), "results", sprintf("yeast_bloom2015_%s.csv", traits)) else NULL)
