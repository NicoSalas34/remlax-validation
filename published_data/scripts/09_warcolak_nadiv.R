# ==============================================================================
# 09_warcolak_nadiv.R -- dominance, epistasis and X-linked additive variance on
# a pedigree
# ------------------------------------------------------------------------------
# Data: the warcolak population shipped with nadiv (Wolak 2012, Methods Ecol.
# Evol. 3:792-796, doi:10.1111/j.2041-210X.2012.00213.x). 5400 individuals from
# 75 replicates of the three-generation design of Fairbairn & Roff (2006,
# Heredity 97:319-328), each with one record per trait. The data set is
# simulated, and its help page gives the code and the variances used:
#   trait1 = sex mean + a (0.4, A) + d (0.3, D) + e (0.3)
#   trait2 = sex mean + a (0.3, A) + s (0.1, S, X-linked, no global dosage
#            compensation) + d (0.3, D) + e (0.3)
# The simulated effects are columns of the data (t1_a, t1_d, ...); their
# realised variances are reported as context.
#
# Relationship inverses, built once by nadiv and passed to every program:
#   A^-1  makeAinv()          additive
#   D^-1  makeD()             dominance (Wolak's approximation for an outbred
#                             pedigree)
#   AA^-1 makeAA()            additive x additive epistasis, A # A
#   S^-1  makeS(ngdc)         X-linked additive, males half of females
#
# Models:
#   trait1  A          y ~ sex + a
#   trait1  A+D        y ~ sex + a + d             (the simulated model)
#   trait1  A+D+AA     y ~ sex + a + d + aa        (aa absent from the data)
#   trait2  A+S+D      y ~ sex + a + s + d         (the simulated model)
#   bivariate A+D      us(trait) on a and d, us(trait) residual
#
# remlax fits them with the mixed model equations engine and the sparse
# inverses (vm(id, Kinv = ...)). asreml gets the same inverses as ginv
# objects. sommer gets the dense matrices for the single-kernel model.
#
#   Rscript published_data/scripts/09_warcolak_nadiv.R
# ==============================================================================
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(.here, "..", "R", "common.R"))
suppressPackageStartupMessages({ library(asreml); library(sommer); library(nadiv) })

start_checks("warcolak_nadiv", "Wolak 2012, doi:10.1111/j.2041-210X.2012.00213.x (simulated data)")
cat(sprintf("  %-12s %s\n", "nadiv", as.character(utils::packageVersion("nadiv"))))
SIM <- "realised variance of the simulated effects (warcolak columns)"

data(warcolak)
ped <- warcolak[, 1:3]
ids <- as.character(warcolak$ID)
Ai  <- makeAinv(ped)$Ainv
Di  <- suppressWarnings(makeD(ped, invertD = TRUE)$Dinv)
AAi <- makeAA(ped)$AAinv
Si  <- suppressWarnings(makeS(warcolak[, 1:4], heterogametic = "M", DosageComp = "ngdc",
                              returnS = FALSE)$Sinv)
stopifnot(identical(rownames(makeA(ped)), ids))
for (nm in c("Ai", "Di", "AAi", "Si")) {
  M <- get(nm); dimnames(M) <- list(ids, ids); assign(nm, M)
}
cat(sprintf("individuals %d, nonzeros A^-1 %d, D^-1 %d, AA^-1 %d, S^-1 %d\n", length(ids),
            Matrix::nnzero(Ai), Matrix::nnzero(Di), Matrix::nnzero(AAi), Matrix::nnzero(Si)))

# asreml ginv: lower triangle in three columns, sorted, with rowNames.
ginv <- function(M) {
  T <- methods::as(methods::as(Matrix::forceSymmetric(M, "L"), "generalMatrix"), "TsparseMatrix")
  k <- T@i >= T@j
  g <- cbind(Row = T@i[k] + 1L, Column = T@j[k] + 1L, Value = T@x[k])
  g <- g[order(g[, "Row"], g[, "Column"]), ]
  attr(g, "rowNames") <- rownames(M); attr(g, "INVERSE") <- TRUE
  g
}
gA <- ginv(Ai); gD <- ginv(Di); gAA <- ginv(AAi); gS <- ginv(Si)

d <- data.frame(id = factor(ids, levels = ids), sex = warcolak$sex,
                y1 = warcolak$trait1, y2 = warcolak$trait2)
d$idD <- d$idAA <- d$idS <- d$id
WS <- "2gb"
rx_mme <- function(...) rx_cpu(..., engine = "mme")
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

# ------------------------------------------------------------------------------
# Univariate models
# ------------------------------------------------------------------------------
uni <- list(
  "trait1 A" = list(y = "y1", rx = ~ vm(id, Kinv = Ai), as = ~ vm(id, gA),
                    comp = c(A = "vm(id, gA)")),
  "trait1 A+D" = list(y = "y1", rx = ~ vm(id, Kinv = Ai) + vm(idD, Kinv = Di),
                      as = ~ vm(id, gA) + vm(idD, gD),
                      comp = c(A = "vm(id, gA)", D = "vm(idD, gD)"),
                      truth = c(A = "t1_a", D = "t1_d", residual = "t1_r")),
  "trait1 A+D+AA" = list(y = "y1",
                         rx = ~ vm(id, Kinv = Ai) + vm(idD, Kinv = Di) + vm(idAA, Kinv = AAi),
                         as = ~ vm(id, gA) + vm(idD, gD) + vm(idAA, gAA),
                         comp = c(A = "vm(id, gA)", D = "vm(idD, gD)", AA = "vm(idAA, gAA)")),
  "trait2 A+S+D" = list(y = "y2",
                        rx = ~ vm(id, Kinv = Ai) + vm(idS, Kinv = Si) + vm(idD, Kinv = Di),
                        as = ~ vm(id, gA) + vm(idS, gS) + vm(idD, gD),
                        comp = c(A = "vm(id, gA)", S = "vm(idS, gS)", D = "vm(idD, gD)"),
                        truth = c(A = "t2_a", S = "t2_s", D = "t2_d", residual = "t2_r")))
rx_term_of <- c(A = "id", D = "idD", AA = "idAA", S = "idS")

for (m in names(uni)) {
  s <- uni[[m]]
  fx <- stats::as.formula(paste(s$y, "~ sex"))
  cat(sprintf("\n=== %s: %s ===\n", m, deparse(s$rx)))
  f <- timed(rx_mme(fx, random = s$rx, data = d))
  print(f)
  a <- timed(stab(asreml(fx, random = s$as, data = d, trace = FALSE, maxit = 100,
                         workspace = WS)))
  vc <- summary(a)$varcomp
  print(vc)
  TIMES[[m]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"))
  for (cn in names(s$comp))
    add_check(m, paste0("sigma2 ", cn), "asreml", f$sigmas[[rx_term_of[[cn]]]][1, 1],
              vcr(vc, s$comp[[cn]]))
  add_check(m, "sigma2 residual", "asreml", f$sigma_res[1, 1], vcr(vc, "units!R"))
  add_check(m, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
            tol = 1e-5, metric = "abs")
  for (cn in names(s$comp)) {
    u_rx <- setNames(f$blups[[rx_term_of[[cn]]]][, 1], rownames(f$blups[[rx_term_of[[cn]]]]))
    if (stats::sd(u_rx) < 1e-6) next  # a component at the bound has null BLUPs
    u_as <- asr_blup(a, paste0(s$comp[[cn]], "_"))[names(u_rx)]
    add_blup_checks(paste(m, cn), "asreml", u_rx, u_as)
  }
  if (!is.null(s$truth))
    for (cn in names(s$truth)) {
      est <- if (cn == "residual") f$sigma_res[1, 1] else f$sigmas[[rx_term_of[[cn]]]][1, 1]
      add_check(m, paste0("sigma2 ", cn), "simulation", est,
                stats::var(warcolak[[s$truth[[cn]]]]), kind = "context", note = SIM)
    }
}

# ------------------------------------------------------------------------------
# sommer on the single-kernel model, dense A (mmes takes covariances)
# ------------------------------------------------------------------------------
A <- as.matrix(makeA(ped)); dimnames(A) <- list(ids, ids)
m <- "trait1 A"
so <- timed(mmes(y1 ~ sex, random = ~ vsm(ism(id), Gu = A), rcov = ~ units, data = d,
                 verbose = FALSE, dateWarning = FALSE))
f <- rx_mme(y1 ~ sex, random = ~ vm(id, Kinv = Ai), data = d)
add_check(m, "sigma2 A", "sommer", f$sigmas$id[1, 1],
          as.numeric(so$theta[[match("vsm(ism(id), Gu = A)", names(so$theta))]]))
add_check(m, "sigma2 residual", "sommer", f$sigma_res[1, 1],
          as.numeric(so$theta[[match("units", names(so$theta))]]))
TIMES[["trait1 A"]] <- c(TIMES[["trait1 A"]], sommer = attr(so, "secs"))

# ------------------------------------------------------------------------------
# Bivariate A+D, us(trait) on both genetic terms and on the residual
# ------------------------------------------------------------------------------
mname <- "bivariate A+D"
cat(sprintf("\n=== %s ===\n", mname))
dl <- rbind(data.frame(d[, c("id", "idD", "sex")], car = "t1", y = d$y1),
            data.frame(d[, c("id", "idD", "sex")], car = "t2", y = d$y2))
dl$car <- factor(dl$car, levels = c("t1", "t2"))
dl <- dl[order(dl$id, dl$car), ]
f <- timed(rx_mme(y ~ car + car:sex, random = ~ us(id, Kinv = Ai) + us(idD, Kinv = Di),
                  residual = ~ us(car):units, data = dl, trait = "car", unit = "id"))
print(f)
# asreml in wide format (a long column named `trait` is reserved in asreml).
a <- timed(stab(asreml(cbind(y1, y2) ~ trait + trait:sex,
                       random = ~ us(trait):vm(id, gA) + us(trait):vm(idD, gD),
                       residual = ~ units:us(trait), data = d, trace = FALSE,
                       maxit = 100, workspace = "4gb")))
vc <- summary(a)$varcomp
print(vc)
TIMES[[mname]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"))
tl <- c("y1", "y2")
get_us <- function(term) {
  S <- matrix(NA_real_, 2, 2)
  for (i in 1:2) for (j in 1:i) {
    k <- match(sprintf("%s!trait_%s:%s", term, tl[i], tl[j]), rownames(vc))
    S[i, j] <- S[j, i] <- vc[k[1], "component"]
  }
  S
}
for (blk in list(list("A", f$sigmas$id, get_us("trait:vm(id, gA)")),
                 list("D", f$sigmas$idD, get_us("trait:vm(idD, gD)")),
                 list("residual", f$sigma_res, get_us("units:trait")))) {
  Srx <- blk[[2]]; Sas <- blk[[3]]
  add_check(mname, paste(blk[[1]], "var trait1"), "asreml", Srx[1, 1], Sas[1, 1])
  add_check(mname, paste(blk[[1]], "var trait2"), "asreml", Srx[2, 2], Sas[2, 2])
  add_check(mname, paste(blk[[1]], "cov (/ sqrt(var var))"), "asreml",
            abs(Srx[1, 2] - Sas[1, 2]) / sqrt(Sas[1, 1] * Sas[2, 2]), NA,
            tol = 1e-3, metric = "value")
}
add_check(mname, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik,
          tol = 1e-5, metric = "abs")

cat("\nRuntimes (s, elapsed):\n")
print(TIMES)
write_checks()
