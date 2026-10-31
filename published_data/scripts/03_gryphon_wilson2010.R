# ==============================================================================
# 03_gryphon_wilson2010.R -- pedigree animal model on the gryphon data
# ------------------------------------------------------------------------------
# Data: the simulated gryphon population of Wilson et al. (2010, J. Anim. Ecol.
# 79:13-26, doi:10.1111/j.1365-2656.2009.01639.x), shipped as DT_gryphon and
# A_gryphon (additive relationship matrix, 1309 individuals) in the enhancer
# package. 1084 records, 854 birth weights (BWT), 683 tarsus lengths (TARSUS).
#
# Published values: the ASReml-R tutorial that accompanies the paper, in its
# updated version (Martin & Videlier, "How to fit an animal model",
# https://wam-tutorial.github.io/en/02_01-univ_asreml.html for the univariate
# models, .../03_01-biv_asreml.html for the bivariate ones). The tutorial
# prints asreml components to 7 digits and log-likelihoods to 3 decimals.
# These are compared at that rounding (components 1e-3 relative, logLik 1e-3
# absolute: half a unit of the last printed decimal plus the stopping
# tolerance of asreml).
#
# The tutorial builds ainv from the pedigree file; here every program gets the
# SAME dense A_gryphon. The enhancer copy of the pedigree (P_gryphon) has no
# sire column and cannot be inverted by ainverse(), so A_gryphon is used. That
# the published components are reproduced shows that A_gryphon is the matrix
# of the tutorial.
#
# Software references, same model and data: asreml (components, logLik,
# BLUPs) and sommer::mmes (univariate models).
#
#   Rscript published_data/scripts/03_gryphon_wilson2010.R
# ==============================================================================
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(.here, "..", "R", "common.R"))
suppressPackageStartupMessages({ library(asreml); library(sommer); library(enhancer) })

start_checks("gryphon_wilson2010", "Wilson et al. 2010, doi:10.1111/j.1365-2656.2009.01639.x")
TUTO <- "wam-tutorial.github.io (updated Wilson et al. 2010 ASReml-R tutorial)"

data(DT_gryphon)
A <- A_gryphon
d <- DT_gryphon
d$ANIMAL <- factor(as.character(d$ANIMAL), levels = rownames(A))
d$unit <- factor(seq_len(nrow(d)))
d1 <- d[!is.na(d$BWT), ]
cat(sprintf("records %d, BWT %d, TARSUS %d, pedigree %d\n",
            nrow(d), nrow(d1), sum(!is.na(d$TARSUS)), nrow(A)))
A_as <- as_asreml_G(A)
WS <- "2gb"

# Component of an asreml varcomp table by exact row name.
vcr <- function(vc, row) {
  if (!row %in% rownames(vc)) stop("no varcomp row ", row, " in: ", paste(rownames(vc), collapse = ", "))
  vc[row, "component"]
}
# BLUPs of asreml for one term, named by level.
asr_blup <- function(a, prefix) {
  cr <- summary(a, coef = TRUE)$coef.random
  k <- startsWith(rownames(cr), prefix)
  setNames(cr[k, "solution"], sub(prefix, "", rownames(cr)[k], fixed = TRUE))
}
timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; attr(v, "secs") <- proc.time()[["elapsed"]] - t0; v }

# ------------------------------------------------------------------------------
# Univariate models 1 to 4 of the tutorial
# ------------------------------------------------------------------------------
uni <- list(
  model1 = list(fixed = BWT ~ 1,       rand_rx = ~ vm(ANIMAL, K = A),
                rand_as = ~ vm(ANIMAL, A_as),
                rand_so = ~ vsm(ism(ANIMAL), Gu = A),
                pub = c(ANIMAL = 3.395398, residual = 3.828602), loglik = NA,
                h2 = c(expr = "V1/(V1+V2)", pub = 0.4700163)),
  model2 = list(fixed = BWT ~ 1 + SEX, rand_rx = ~ vm(ANIMAL, K = A),
                rand_as = ~ vm(ANIMAL, A_as),
                rand_so = ~ vsm(ism(ANIMAL), Gu = A),
                pub = c(ANIMAL = 3.060441, residual = 2.938412), loglik = NA,
                h2 = c(expr = "V1/(V1+V2)", pub = 0.510171)),
  model3 = list(fixed = BWT ~ 1 + SEX, rand_rx = ~ vm(ANIMAL, K = A) + BYEAR,
                rand_as = ~ vm(ANIMAL, A_as) + BYEAR,
                rand_so = ~ vsm(ism(ANIMAL), Gu = A) + BYEAR,
                pub = c(ANIMAL = 2.7068665, BYEAR = 0.8862604, residual = 2.3092415),
                loglik = -1113.772,
                h2 = c(expr = "V1/(V1+V2+V3)", pub = 0.4586068)),
  model4 = list(fixed = BWT ~ 1 + SEX, rand_rx = ~ vm(ANIMAL, K = A) + BYEAR + MOTHER,
                rand_as = ~ vm(ANIMAL, A_as) + BYEAR + MOTHER,
                rand_so = ~ vsm(ism(ANIMAL), Gu = A) + BYEAR + MOTHER,
                pub = c(ANIMAL = 2.2985320, BYEAR = 0.8820313, MOTHER = 1.1184698,
                        residual = 1.6290034),
                loglik = -1094.815,
                # The tutorial prints h2 = 0.1487898 for model 4 with V1 in the
                # numerator, which is the BYEAR share in its component order.
                # The animal heritability is compared to its own components.
                h2 = c(expr = "V1/(V1+V2+V3+V4)", pub = NA)))
as_names <- c(ANIMAL = "vm(ANIMAL, A_as)", BYEAR = "BYEAR", MOTHER = "MOTHER",
              residual = "units!R")
# mmes rewrites a bare factor as an explicit vsm(ism()) term in its labels.
so_names <- c(ANIMAL = "vsm(ism(ANIMAL), Gu = A)",
              BYEAR = "sommer::vsm( sommer::ism( BYEAR ) )",
              MOTHER = "sommer::vsm( sommer::ism( MOTHER ) )", residual = "units")
TIMES <- list()

for (m in names(uni)) {
  s <- uni[[m]]
  cat(sprintf("\n=== %s: %s | %s ===\n", m, deparse(s$fixed), deparse(s$rand_rx)))
  f <- timed(rx_cpu(s$fixed, random = s$rand_rx, data = d1,
                    vpredict = c(h2 = s$h2[["expr"]])))
  print(f)
  a <- timed(stab(asreml(s$fixed, random = s$rand_as, data = d1, trace = FALSE,
                         maxit = 100, workspace = WS)))
  vc <- summary(a)$varcomp
  print(vc)
  so <- timed(mmes(s$fixed, random = s$rand_so, rcov = ~ units, data = d1,
                   verbose = FALSE, dateWarning = FALSE))
  TIMES[[m]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"), sommer = attr(so, "secs"))

  comps <- names(s$pub)
  rx_val <- function(cn) if (cn == "residual") f$sigma_res[1, 1] else f$sigmas[[cn]][1, 1]
  for (cn in comps) {
    add_check(m, paste0("sigma2 ", cn), "asreml", rx_val(cn), vcr(vc, as_names[[cn]]))
    j <- match(so_names[[cn]], names(so$theta))
    add_check(m, paste0("sigma2 ", cn), "sommer", rx_val(cn),
              if (is.na(j)) NA else as.numeric(so$theta[[j]]))
    add_check(m, paste0("sigma2 ", cn), "tutorial", rx_val(cn), s$pub[[cn]],
              kind = "published", note = TUTO)
  }
  add_check(m, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
            tol = 1e-5, metric = "abs")
  if (is.finite(s$loglik))
    add_check(m, "logLik (asreml convention)", "tutorial", f$logLik_asreml, s$loglik,
              tol = 1e-3, metric = "abs", kind = "published", note = TUTO)
  h2_rx <- f$vpredict$predictions$value[1]
  if (is.finite(s$h2[["pub"]]))
    add_check(m, "h2", "tutorial", h2_rx, as.numeric(s$h2[["pub"]]),
              kind = "published", note = TUTO)
  # BLUPs of the animals with a record (remlax reports those).
  u_rx <- setNames(f$blups$ANIMAL[, 1], rownames(f$blups$ANIMAL))
  u_as <- asr_blup(a, "vm(ANIMAL, A_as)_")[names(u_rx)]
  add_blup_checks(m, "asreml", u_rx, u_as)
  # mmes stacks the BLUPs of all terms in declaration order; the animal term
  # comes first and has one row per pedigree member.
  u_so <- setNames(so$u[seq_len(nrow(A)), 1], rownames(so$u)[seq_len(nrow(A))])
  stopifnot(identical(names(u_so), rownames(A)))
  add_blup_checks(m, "sommer", u_rx, u_so[names(u_rx)])
}

# ------------------------------------------------------------------------------
# Bivariate models of the tutorial (BWT, TARSUS), long format in remlax
# ------------------------------------------------------------------------------
long <- rbind(data.frame(d[, c("ANIMAL", "MOTHER", "BYEAR", "SEX", "unit")], trait = "BWT",    y = d$BWT),
              data.frame(d[, c("ANIMAL", "MOTHER", "BYEAR", "SEX", "unit")], trait = "TARSUS", y = d$TARSUS))
long$trait <- factor(long$trait, levels = c("BWT", "TARSUS"))
long <- long[!is.na(long$y), ]
long <- long[order(long$unit, long$trait), ]

pair <- function(M) c(M[1, 1], M[2, 1], M[2, 2])
asr_pair <- function(vc, term) {
  c(vcr(vc, paste0(term, "!trait_BWT:BWT")), vcr(vc, paste0(term, "!trait_TARSUS:BWT")),
    vcr(vc, paste0(term, "!trait_TARSUS:TARSUS")))
}
lab3 <- c("BWT", "cov", "TARSUS")
bi_check <- function(m, comp, rx3, ref3, software, kind = "software", note = "") {
  for (k in 1:3)
    add_check(m, paste0(comp, " ", lab3[k]), software, rx3[k], ref3[k], kind = kind, note = note)
}

# modela: animal us + residual us, mean per trait.
cat("\n=== modela: cbind(BWT, TARSUS) ~ trait, us animal, us residual ===\n")
fa_ <- timed(rx_cpu(y ~ trait, random = ~ us(ANIMAL, K = A), residual = ~ us(trait):units,
                    data = long, trait = "trait", unit = "unit", n_restarts = 3))
print(fa_)
aa <- timed(stab(asreml(cbind(BWT, TARSUS) ~ trait,
                        random = ~ us(trait, init = c(1, 0.1, 1)):vm(ANIMAL, A_as),
                        residual = ~ id(units):us(trait, init = c(1, 0.1, 1)),
                        data = d, trace = FALSE, maxit = 100, workspace = WS)))
vca <- summary(aa)$varcomp
print(vca)
TIMES$modela <- c(remlax = attr(fa_, "secs"), asreml = attr(aa, "secs"))
bi_check("modela", "animal", pair(fa_$sigmas$ANIMAL), asr_pair(vca, "trait:vm(ANIMAL, A_as)"), "asreml")
bi_check("modela", "residual", pair(fa_$sigma_res), asr_pair(vca, "units:trait"), "asreml")
bi_check("modela", "animal", pair(fa_$sigmas$ANIMAL), c(3.368397, 2.459809, 12.345792),
         "tutorial", "published", TUTO)
bi_check("modela", "residual", pair(fa_$sigma_res), c(3.849916, 3.313282, 17.646432),
         "tutorial", "published", TUTO)
add_check("modela", "logLik (asreml convention)", "asreml", fa_$logLik_asreml, aa$loglik,
          tol = 1e-5, metric = "abs")
add_check("modela", "logLik (asreml convention)", "tutorial", fa_$logLik_asreml, -2679.741,
          tol = 1e-3, metric = "abs", kind = "published",
          note = paste(TUTO, "(last line of the update() iteration log)"))

# modelb: trait-specific sex effect, us on animal, birth year and mother.
cat("\n=== modelb: + at(trait):sex, us animal + us byear + us mother ===\n")
fb <- timed(rx_cpu(y ~ trait + trait:SEX,
                   random = ~ us(ANIMAL, K = A) + us(BYEAR) + us(MOTHER),
                   residual = ~ us(trait):units,
                   data = long, trait = "trait", unit = "unit", n_restarts = 3))
print(fb)
ab <- timed(stab(asreml(cbind(BWT, TARSUS) ~ trait + at(trait):SEX,
                        random = ~ us(trait, init = c(1, 0.1, 1)):vm(ANIMAL, A_as) +
                          us(trait, init = c(1, 0.1, 1)):BYEAR +
                          us(trait, init = c(1, 0.1, 1)):MOTHER,
                        residual = ~ id(units):us(trait, init = c(1, 0.1, 1)),
                        data = d, trace = FALSE, maxit = 100, workspace = WS)))
vcb <- summary(ab)$varcomp
print(vcb)
TIMES$modelb <- c(remlax = attr(fb, "secs"), asreml = attr(ab, "secs"))
bi_check("modelb", "animal",   pair(fb$sigmas$ANIMAL), asr_pair(vcb, "trait:vm(ANIMAL, A_as)"), "asreml")
bi_check("modelb", "byear",    pair(fb$sigmas$BYEAR),  asr_pair(vcb, "trait:BYEAR"), "asreml")
bi_check("modelb", "mother",   pair(fb$sigmas$MOTHER), asr_pair(vcb, "trait:MOTHER"), "asreml")
bi_check("modelb", "residual", pair(fb$sigma_res),     asr_pair(vcb, "units:trait"), "asreml")
bi_check("modelb", "animal",   pair(fb$sigmas$ANIMAL), c(1.9893546, 3.3170404, 10.2294887), "tutorial", "published", TUTO)
bi_check("modelb", "byear",    pair(fb$sigmas$BYEAR),  c(0.9746385, 0.1624076, 3.7383721),  "tutorial", "published", TUTO)
bi_check("modelb", "mother",   pair(fb$sigmas$MOTHER), c(1.1445184, -1.5567306, 4.8206132), "tutorial", "published", TUTO)
bi_check("modelb", "residual", pair(fb$sigma_res),     c(1.8443110, 4.0142841, 12.4845955), "tutorial", "published", TUTO)
add_check("modelb", "logLik (asreml convention)", "asreml", fb$logLik_asreml, ab$loglik,
          tol = 1e-5, metric = "abs")
add_check("modelb", "logLik (asreml convention)", "tutorial", fb$logLik_asreml, -2380.246,
          tol = 1e-3, metric = "abs", kind = "published", note = TUTO)
rA <- fb$sigmas$ANIMAL[2, 1] / sqrt(fb$sigmas$ANIMAL[1, 1] * fb$sigmas$ANIMAL[2, 2])
add_check("modelb", "genetic correlation r_A", "tutorial", rA, 0.7353053,
          kind = "published", note = TUTO)
u_rx <- fb$blups$ANIMAL
cat("remlax BLUP matrix for the bivariate animal term:", dim(u_rx), "\n")

# modelc: maternal covariance fixed at zero (diag instead of us on mother).
cat("\n=== modelc: maternal covariance set to zero ===\n")
fc <- timed(rx_cpu(y ~ trait + trait:SEX,
                   random = ~ us(ANIMAL, K = A) + us(BYEAR) + diag(MOTHER),
                   residual = ~ us(trait):units,
                   data = long, trait = "trait", unit = "unit", n_restarts = 3))
print(fc)
ac <- timed(stab(asreml(cbind(BWT, TARSUS) ~ trait + at(trait):SEX,
                        random = ~ us(trait, init = c(1, 0.1, 1)):vm(ANIMAL, A_as) +
                          us(trait, init = c(1, 0.1, 1)):BYEAR +
                          diag(trait, init = c(1, 1)):MOTHER,
                        residual = ~ id(units):us(trait, init = c(1, 0.1, 1)),
                        data = d, trace = FALSE, maxit = 100, workspace = WS)))
print(summary(ac)$varcomp)
add_check("modelc", "logLik (asreml convention)", "asreml", fc$logLik_asreml, ac$loglik,
          tol = 1e-5, metric = "abs")
add_check("modelc", "logLik (asreml convention)", "tutorial", fc$logLik_asreml, -2386.045,
          tol = 1e-3, metric = "abs", kind = "published", note = TUTO)
add_check("modelb vs modelc", "LRT statistic", "tutorial",
          2 * (fb$logLik_asreml - fc$logLik_asreml), 11.59835,
          tol = 2e-3, metric = "abs", kind = "published",
          note = paste(TUTO, "; tolerance = twice the logLik rounding"))

cat("\nElapsed seconds (remlax includes process start-up):\n")
print(do.call(rbind, lapply(TIMES, function(x) x[c("remlax", "asreml", "sommer")])))
write_checks()
