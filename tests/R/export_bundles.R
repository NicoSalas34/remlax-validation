# ==============================================================================
# export_bundles.R -- writes one serialised design per model, WITHOUT FITTING
# ANYTHING.
#
# WHY THIS SPLIT. On the cluster, R and JAX live in two separate containers and
# apptainer does not nest (chain of missing libraries: libsubid, libcrypt...).
# R therefore cannot call the solver. It writes the bundles, the JAX container
# fits them, and the comparison is made on files. This is also what makes the
# CPU/GPU parity verifiable, both backends reading STRICTLY the same input.
#
#   Rscript scripts/tests/export_bundles.R --out=/path/bundles [--root=<repo>]
#
# The battery covers the WHOLE catalogue, including the structures added on
# 01/09/2026 (lvr, ilv, anisotropic mtrn, own, dsum).
# ==============================================================================
suppressPackageStartupMessages({ library(Matrix) })
getarg <- function(k, d) { a <- commandArgs(TRUE); i <- grep(paste0("^--", k, "="), a)
  if (length(i)) sub(paste0("^--", k, "="), "", a[i[1]]) else d }
ROOT <- getarg("root", getarg("racine", "."))
OUT  <- getarg("out", file.path(ROOT, "output", "tests", "bundles"))
source(file.path(ROOT, "R", "charger_remlax.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

set.seed(20260901)
n_g <- 40; n_rep <- 4; n <- n_g * n_rep
gid <- factor(rep(sprintf("G%02d", 1:n_g), each = n_rep))
block <- factor(rep(1:n_rep, n_g))
# plausible GRM: product of centred markers
Mk <- matrix(rbinom(n_g * 300, 2, 0.3), n_g, 300)
Mc <- scale(Mk, scale = FALSE)
K <- tcrossprod(Mc) / mean(diag(tcrossprod(Mc)))
dimnames(K) <- list(levels(gid), levels(gid))
u  <- as.numeric(t(chol(K + diag(1e-6, n_g))) %*% rnorm(n_g))
d  <- data.frame(gid = gid, block = block, x = rnorm(n))
d$row <- factor(rep(1:10, length.out = n)); d$col <- factor(rep(1:16, each = 10))
d$row_n <- as.numeric(as.character(d$row)); d$col_n <- as.numeric(as.character(d$col))
d$site  <- factor(ifelse(seq_len(n) <= n / 2, "A", "B"))
d$y <- 3 + 0.5 * d$x + u[as.integer(gid)] + rnorm(n, 0, 0.8)
Zx <- model.matrix(~ 0 + gid, d) * d$x       # slope per genotype, SAME levels
colnames(Zx) <- levels(gid)
d$Zx <- Zx
Zx <- model.matrix(~ 0 + gid, d) * d$x       # slope per genotype, same levels
colnames(Zx) <- levels(gid)
d$Zx <- Zx                                   # reachable from the formula

# --- multi-trait data set (long format) --------------------------------------
tr <- c("t1", "t2", "t3")
dl <- do.call(rbind, lapply(tr, function(t_) transform(d, trait = t_, unit = seq_len(n))))
dl$trait <- factor(dl$trait)
Sg <- matrix(c(1, .6, .3, .6, 1.2, .4, .3, .4, .8), 3, 3)
Ug <- t(chol(Sg)) %*% matrix(rnorm(3 * n_g), 3, n_g)
dl$y <- 2 + Ug[cbind(as.integer(dl$trait), as.integer(dl$gid))] + rnorm(nrow(dl), 0, 0.7)

# --- spatial field, one observation per cell ---------------------------------
ds <- d[!duplicated(paste(d$row_n, d$col_n)), ]
# One observation per (section, column), as a structure between units within a
# section requires.
dsec <- d[!duplicated(paste(d$site, d$col)), ]

models <- list(
  iid        = list(f = y ~ 1 + x, r = ~ iid(gid), dat = d),
  vm_kin     = list(f = y ~ 1 + x, r = ~ vm(gid, K), dat = d),
  ar1        = list(f = y ~ 1, r = ~ ar1(col), dat = d),
  ar2        = list(f = y ~ 1, r = ~ ar2(col), dat = d),
  ar3        = list(f = y ~ 1, r = ~ ar3(col), dat = d),
  sar        = list(f = y ~ 1, r = ~ sar(col), dat = d),
  ma1        = list(f = y ~ 1, r = ~ ma1(col), dat = d),
  arma       = list(f = y ~ 1, r = ~ arma(col), dat = d),
  cor        = list(f = y ~ 1, r = ~ cor(col), dat = d),
  corb2      = list(f = y ~ 1, r = ~ corb(col, order = 2), dat = d),
  corg       = list(f = y ~ 1, r = ~ corg(block), dat = d),
  expo       = list(f = y ~ 1, r = ~ exp(col_n), dat = d),
  gau        = list(f = y ~ 1, r = ~ gau(col_n), dat = d),
  lvr        = list(f = y ~ 1, r = ~ lvr(col_n), dat = d),
  iexp       = list(f = y ~ 1, r = ~ iexp(row_n, col_n), dat = ds),
  ieuc       = list(f = y ~ 1, r = ~ ieuc(row_n, col_n), dat = ds),
  sph        = list(f = y ~ 1, r = ~ sph(row_n, col_n), dat = ds),
  cir        = list(f = y ~ 1, r = ~ cir(row_n, col_n), dat = ds),
  aexp       = list(f = y ~ 1, r = ~ aexp(row_n, col_n), dat = ds),
  mtrn_iso   = list(f = y ~ 1, r = ~ mtrn(row_n, col_n, phi = 3, nu = "0.8 F"), dat = ds),
  mtrn_aniso = list(f = y ~ 1, r = ~ mtrn(row_n, col_n, phi = 3, nu = "1.5 F",
                                          delta = 1.5, alpha = 0.3), dat = ds),
  own_exp    = list(f = y ~ 1, r = ~ own(col, expr = "exp(-lag*exp(p1))", n_par = 1), dat = d),
  own_2par   = list(f = y ~ 1, r = ~ own(col, expr = "exp(-(lag/exp(p1))^(1+tanh(p2)))",
                                         n_par = 2), dat = d),
  ar1ar1_res = list(f = y ~ 1, r = ~ iid(gid), res = ~ ar1(row):ar1(col), dat = ds),
  us_mv      = list(f = y ~ 1, r = ~ us(gid), dat = dl, trait = "trait", unit = "unit"),
  fa2_mv     = list(f = y ~ 1, r = ~ fa(gid, rank = 2), dat = dl, trait = "trait", unit = "unit"),
  chol2_mv   = list(f = y ~ 1, r = ~ chol(gid, rank = 2), dat = dl, trait = "trait", unit = "unit"),
  ante1_mv   = list(f = y ~ 1, r = ~ ante(gid, rank = 1), dat = dl, trait = "trait", unit = "unit"),
  rr2_mv     = list(f = y ~ 1, r = ~ rr(gid, rank = 2), dat = dl, trait = "trait", unit = "unit"),
  corh_mv    = list(f = y ~ 1, r = ~ corh(gid), dat = dl, trait = "trait", unit = "unit"),
  us_res_mv  = list(f = y ~ 1, r = ~ iid(gid), res = ~ us(trait):units, dat = dl,
                    trait = "trait", unit = "unit"),
  dsum       = list(f = y ~ 1, r = ~ iid(gid), res = ~ dsum(~ units | site), dat = d),
  # dsum + ar1: a single observation per (section, column). Otherwise two rows
  # share a unit, their correlation is 1 and R is singular.
  dsum_ar1   = list(f = y ~ 1, r = ~ iid(gid), res = ~ dsum(~ ar1(col) | site), dat = dsec),
  # str(): TWO incidences on the SAME levels (direct effect and slope on x),
  # with a free covariance between them. The columns of the second must carry
  # the same names as the levels, otherwise rx_term refuses, rightly.
  str_2      = list(f = y ~ 1, r = ~ str(~ gid + mm(Zx, name = "gid_x"), struct = "us"),
                    dat = d)
)

ok <- 0L; ko <- character(0)
for (nm in names(models)) {
  m <- models[[nm]]
  res <- tryCatch({
    dd <- m$dat
    tr_ <- if (!is.null(m$trait)) dd[[m$trait]] else NULL
    un_ <- if (!is.null(m$unit))  dd[[m$unit]]  else NULL
    X <- model.matrix(m$f, dd)
    asg <- attr(X, "assign")
    attr(X, "term_labels") <- c("(Intercept)", attr(terms(m$f), "term.labels"))[sort(unique(asg)) + 1L]
    terms_l <- .rx_parse_random(m$r, dd, trait = tr_)
    res_o <- .rx_parse_residual(m$res %||% "units", dd,
                                trait = if (!is.null(m$trait)) m$trait else NULL,
                                unit  = if (!is.null(m$unit))  m$unit  else NULL)
    mod <- rx_model(model.response(model.frame(m$f, dd)), X, terms_l, res_o, name = nm)
    rx_export(mod, file.path(OUT, nm))
    sprintf("%-12s n=%-5d p=%-2d %d parameter(s)", nm, mod$n, ncol(X), mod$n_par)
  }, error = function(e) structure(conditionMessage(e), class = "error_msg"))
  if (inherits(res, "error_msg")) { ko <- c(ko, sprintf("%-12s %s", nm, res)) }
  else { ok <- ok + 1L; cat("  ", res, "\n") }
}
cat(sprintf("\n%d/%d bundles written to %s\n", ok, length(models), OUT))
if (length(ko)) { cat("FAILURES :\n"); cat(paste0("  ", ko, collapse = "\n"), "\n"); quit(status = 1) }
