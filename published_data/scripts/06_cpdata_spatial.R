# ==============================================================================
# 06_cpdata_spatial.R -- genomic BLUP with field structure (sommer CP data)
# ------------------------------------------------------------------------------
# Data: DT_cpdata, GT_cpdata (enhancer package, formerly in sommer): 363 full
# sibs of an F1 cross between two heterozygous parents, 2889 markers coded
# -1/0/1, planted on a 13 x 36 grid (Row, Col; 363 of the 468 cells occupied).
# The help page of the data states that it was SIMULATED for fruit breeding
# applications. It is therefore no published experiment. It is kept because
# it is the example of the sommer paper (Covarrubias-Pazaran 2016, PLoS ONE
# 11:e0156744, doi:10.1371/journal.pone.0156744) for genomic models with field
# structure, and because it tests remlax on an incomplete grid.
#
# Models, identical in every program (G built once by VanRaden's first method):
#   gblup_rowcol   y ~ 1 + G + Row (iid) + Col (iid)
#   gblup_ar1res   y ~ 1 + G, residual AR1(Row) x AR1(Col) on the grid
#   gblup_ar1fld   y ~ 1 + G + random field AR1(Row) x AR1(Col) + iid residual
# for Yield and Firmness.
#
# Two differences between the programs surface here and are handled in the
# open. asreml fixes an AR1 correlation at its 0.98 bound; remlax is then
# refitted with the same value fixed. remlax drops the grid lines that carry
# no record (columns 28 and 32 to 36 for Firmness), so that the columns on
# either side of an empty interior column become adjacent; asreml is refitted
# on that collapsed grid as a diagnostic. The random field of asreml, fitted on
# the observed records only, collapses the same lines (the two programs agree
# there); its separable residual needs the full grid and keeps them.
#
# References: asreml (all models) and sommer::mmes (gblup_rowcol). For the
# separable residual, asreml needs the full grid: the 105 empty cells are added
# as rows with a missing response, sorted by Row then Col. remlax works on the
# observed cells. Both then maximise the same restricted likelihood of the
# observed records.
#
#   Rscript published_data/scripts/06_cpdata_spatial.R
# ==============================================================================
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(.here, "..", "R", "common.R"))
suppressPackageStartupMessages({ library(asreml); library(sommer); library(enhancer) })

start_checks("cpdata_simulated", "Covarrubias-Pazaran 2016, doi:10.1371/journal.pone.0156744 (simulated data)")
data(DT_cpdata)
stopifnot(identical(rownames(GT_cpdata), as.character(DT_cpdata$id)))
G <- vanraden(GT_cpdata + 1)
G_as <- as_asreml_G(G)
WS <- "2gb"
NR <- 13; NC <- 36

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

base <- data.frame(id = factor(DT_cpdata$id, levels = rownames(G)),
                   Row = factor(DT_cpdata$Row, levels = 1:NR),
                   Col = factor(DT_cpdata$Col, levels = 1:NC),
                   Yield = DT_cpdata$Yield, Firmness = DT_cpdata$Firmness)
# Full grid for asreml's separable residual: one row per cell, sorted.
grid <- expand.grid(Col = factor(1:NC, levels = 1:NC), Row = factor(1:NR, levels = 1:NR))[, c("Row", "Col")]
full <- merge(grid, base, by = c("Row", "Col"), all.x = TRUE, sort = FALSE)
full <- full[order(full$Row, full$Col), ]
full$id <- factor(as.character(full$id), levels = rownames(G))
cat(sprintf("records %d on a %d x %d grid (%d cells)\n", nrow(base), NR, NC, nrow(full)))

for (tr in c("Yield", "Firmness")) {
  d <- base[!is.na(base[[tr]]), ]; d$y <- d[[tr]]
  fu <- full; fu$y <- fu[[tr]]
  cat(sprintf("\n================ trait %s: %d records ================\n", tr, nrow(d)))

  # --- gblup_rowcol --------------------------------------------------------------
  m <- paste0(tr, " gblup_rowcol")
  cat("\n===", m, "===\n")
  f <- timed(rx_cpu(y ~ 1, random = ~ vm(id, K = G) + iid(Row) + iid(Col), data = d))
  print(f)
  a <- timed(stab(asreml(y ~ 1, random = ~ vm(id, G_as) + Row + Col, data = d,
                         trace = FALSE, maxit = 100, workspace = WS)))
  vc <- summary(a)$varcomp; print(vc)
  so <- timed(mmes(y ~ 1, random = ~ vsm(ism(id), Gu = G) + Row + Col, rcov = ~ units,
                   data = d, verbose = FALSE, dateWarning = FALSE, nIters = 200,
                   tolParConvLL = 1e-8, tolParConvNorm = 1e-8))
  TIMES[[m]] <- c(remlax = attr(f, "secs"), asreml = attr(a, "secs"), sommer = attr(so, "secs"))
  rx  <- c(G = f$sigmas$id[1, 1], Row = f$sigmas$Row[1, 1], Col = f$sigmas$Col[1, 1],
           residual = f$sigma_res[1, 1])
  asr <- c(G = vcr(vc, "vm(id, G_as)"), Row = vcr(vc, "Row"), Col = vcr(vc, "Col"),
           residual = vcr(vc, "units!R"))
  som <- c(G = so$theta[["vsm(ism(id), Gu = G)"]][1, 1],
           Row = so$theta[["sommer::vsm( sommer::ism( Row ) )"]][1, 1],
           Col = so$theta[["sommer::vsm( sommer::ism( Col ) )"]][1, 1],
           residual = so$theta[["units"]][1, 1])
  for (cn in names(rx)) {
    add_check(m, paste0("sigma2 ", cn), "asreml", rx[[cn]], asr[[cn]])
    add_check(m, paste0("sigma2 ", cn), "sommer", rx[[cn]], som[[cn]])
  }
  add_check(m, "logLik (asreml convention)", "asreml", f$logLik_asreml, a$loglik, tol = 1e-5, metric = "abs")
  u_rx <- setNames(f$blups$id[, 1], rownames(f$blups$id))
  add_blup_checks(m, "asreml", u_rx, asr_blup(a, "vm(id, G_as)_")[names(u_rx)])
  u_so <- so$uList[["vsm(ism(id), Gu = G)"]]
  add_blup_checks(m, "sommer", u_rx, setNames(u_so[, 1], rownames(u_so))[names(u_rx)])

  # --- gblup_ar1res: separable AR1 residual on the incomplete grid ---------------
  m <- paste0(tr, " gblup_ar1res")
  cat("\n===", m, "===\n")
  f2 <- timed(rx_cpu(y ~ 1, random = ~ vm(id, K = G), residual = ~ ar1(Row):ar1(Col),
                     data = d, n_restarts = 3))
  print(f2)
  a2 <- timed(stab(asreml(y ~ 1, random = ~ vm(id, G_as), residual = ~ ar1(Row):ar1(Col),
                          data = fu, trace = FALSE, maxit = 100, workspace = WS,
                          na.action = na.method(y = "include", x = "include"))))
  vc2 <- summary(a2)$varcomp; print(vc2)
  TIMES[[m]] <- c(remlax = attr(f2, "secs"), asreml = attr(a2, "secs"))
  empty_cols <- setdiff(levels(base$Col), as.character(unique(d$Col)))
  empty_rows <- setdiff(levels(base$Row), as.character(unique(d$Row)))
  cat(sprintf("grid lines without a record: rows {%s}, columns {%s}\n",
              paste(empty_rows, collapse = ","), paste(empty_cols, collapse = ",")))
  rho_res <- Filter(function(r) length(r) == 2 && all(is.finite(r)), f2$rho)[[1]]
  add_check(m, "sigma2 G", "asreml", f2$sigmas$id[1, 1], vcr(vc2, "vm(id, G_as)"))
  add_check(m, "sigma2 residual", "asreml", f2$sigma_res[1, 1], vcr(vc2, "Row:Col!R"))
  add_check(m, "rho Row", "asreml", rho_res[1], vcr(vc2, "Row:Col!Row!cor"))
  add_check(m, "rho Col", "asreml", rho_res[2], vcr(vc2, "Row:Col!Col!cor"))
  add_check(m, "logLik (asreml convention)", "asreml", f2$logLik_asreml, a2$loglik, tol = 1e-5, metric = "abs")
  u_rx <- setNames(f2$blups$id[, 1], rownames(f2$blups$id))
  add_blup_checks(m, "asreml", u_rx, asr_blup(a2, "vm(id, G_as)_")[names(u_rx)])

  # Diagnostic. remlax builds the grid from the levels that carry a record
  # (factor() in .rx_residual_ar1ar1 drops unused levels), so an empty column
  # vanishes and its two neighbours become adjacent. asreml refitted on the
  # grid with the empty columns removed tells whether this explains a gap.
  if (length(empty_cols) || length(empty_rows)) {
    kc <- setdiff(levels(base$Col), empty_cols); kr <- setdiff(levels(base$Row), empty_rows)
    fc <- fu[fu$Col %in% kc & fu$Row %in% kr, ]
    fc$Col <- factor(as.character(fc$Col), levels = kc); fc$Row <- factor(as.character(fc$Row), levels = kr)
    fc <- fc[order(fc$Row, fc$Col), ]
    a2c <- stab(asreml(y ~ 1, random = ~ vm(id, G_as), residual = ~ ar1(Row):ar1(Col),
                       data = fc, trace = FALSE, maxit = 100, workspace = WS,
                       na.action = na.method(y = "include", x = "include")))
    vc2c <- summary(a2c)$varcomp; print(vc2c)
    nt <- sprintf("asreml on the grid without the empty lines (%d x %d)", length(kr), length(kc))
    add_check(m, "logLik, grid without empty lines", "asreml", f2$logLik_asreml, a2c$loglik,
              tol = 1e-5, metric = "abs", kind = "context", note = nt)
    add_check(m, "rho Col, grid without empty lines", "asreml", rho_res[2],
              vcr(vc2c, "Row:Col!Col!cor"), kind = "context", note = nt)
    add_check(m, "sigma2 G, grid without empty lines", "asreml", f2$sigmas$id[1, 1],
              vcr(vc2c, "vm(id, G_as)"), kind = "context", note = nt)
  }

  # --- gblup_ar1fld: random AR1 x AR1 field plus nugget ---------------------------
  m <- paste0(tr, " gblup_ar1fld")
  cat("\n===", m, "===\n")
  f3 <- timed(rx_cpu(y ~ 1, random = ~ vm(id, K = G) + ar1(Row, Col), data = d, n_restarts = 5))
  print(f3)
  a3 <- timed(stab(asreml(y ~ 1, random = ~ vm(id, G_as) + ar1(Row):ar1(Col), data = d,
                          trace = FALSE, maxit = 100, workspace = WS)))
  vc3 <- summary(a3)$varcomp; print(vc3)
  TIMES[[m]] <- c(remlax = attr(f3, "secs"), asreml = attr(a3, "secs"))
  fld <- setdiff(names(f3$sigmas), "id")
  cor_rows <- c("Row:Col!Row!cor", "Row:Col!Col!cor")
  # asreml fixes an AR1 correlation that reaches 0.98 in absolute value
  # (bound "F"). The free fit is then not the same model: it is reported as
  # context, and remlax is refitted with that correlation fixed at the asreml
  # value, which is the comparison that counts.
  at_bound <- vc3[cor_rows, "bound"] == "F"
  free_kind <- if (any(at_bound)) "context" else "software"
  free_note <- if (any(at_bound)) paste("asreml fixed", paste(cor_rows[at_bound], collapse = ", "),
                                        "at", paste(signif(vc3[cor_rows[at_bound], "component"], 4), collapse = ", ")) else ""
  cmp_fld <- function(fit, label, kind, note) {
    add_check(label, "sigma2 G", "asreml", fit$sigmas$id[1, 1], vcr(vc3, "vm(id, G_as)"), kind = kind, note = note)
    add_check(label, "sigma2 field", "asreml", fit$sigmas[[fld]][1, 1], vcr(vc3, "Row:Col"), kind = kind, note = note)
    add_check(label, "rho Row", "asreml", fit$rho[[fld]][1], vcr(vc3, cor_rows[1]), kind = kind, note = note)
    add_check(label, "rho Col", "asreml", fit$rho[[fld]][2], vcr(vc3, cor_rows[2]), kind = kind, note = note)
    add_check(label, "sigma2 residual", "asreml", fit$sigma_res[1, 1], vcr(vc3, "units!R"), kind = kind, note = note)
    add_check(label, "logLik (asreml convention)", "asreml", fit$logLik_asreml, a3$loglik,
              tol = 1e-5, metric = "abs", kind = kind, note = note)
    if (kind == "software") {
      u_rx <- setNames(fit$blups$id[, 1], rownames(fit$blups$id))
      add_blup_checks(label, "asreml", u_rx, asr_blup(a3, "vm(id, G_as)_")[names(u_rx)], note = note)
    }
  }
  cmp_fld(f3, m, free_kind, free_note)
  if (any(at_bound)) {
    k <- match(paste0(fld, "!ar1ar1", which(at_bound)), f3$component_names)
    th <- f3$theta
    th[k] <- atanh(vc3[cor_rows[at_bound], "component"])
    f3c <- rx_cpu(y ~ 1, random = ~ vm(id, K = G) + ar1(Row, Col), data = d,
                  theta_init = th, fixed_theta = k)
    print(f3c)
    cmp_fld(f3c, paste0(m, " rho fixed"), "software",
            paste("remlax with", paste(cor_rows[at_bound], collapse = ", "), "fixed at the asreml value"))
  }
}

cat("\nElapsed seconds (remlax includes process start-up):\n")
for (m in names(TIMES)) cat(sprintf("  %-26s %s\n", m,
  paste(sprintf("%s %.1f", names(TIMES[[m]]), TIMES[[m]]), collapse = " | ")))
write_checks()
