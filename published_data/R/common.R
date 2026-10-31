# ==============================================================================
# common.R -- shared helpers for the validation of remlax on published data sets
# ------------------------------------------------------------------------------
# Each script of published_data/ fits remlax on a published data set and
# compares it with established software on the SAME model: asreml (REML, the
# reference of the field), sommer, gaston, rrBLUP, statgenGWAS. When the source
# paper reports estimates for a model that can be reproduced exactly, these are
# compared too, at the rounding of the paper.
#
# Every comparison becomes one row of a checks table, written by
# write_checks() to published_data/results/<dataset>.csv. The columns are the
# same for every data set, so that make_report.py can pool them.
#
#   kind = "software"   same model, same data, another program: pass/fail
#   kind = "published"  value printed in the source paper, same model: pass/fail
#                       at the rounding of the paper
#   kind = "context"    value from the paper under a model that differs (other
#                       software, other fixed effects, CV scheme): reported,
#                       never counted as a pass or a fail
#
# LOG-LIKELIHOOD CONVENTION. asreml omits the constant (n-p)/2 log(2 pi);
# remlax reports it in fit$logLik and gives fit$logLik_asreml without it.
#
# TOLERANCES (defaults, the same as the asreml suite of the campaign):
#   components 1e-3 relative, log-likelihood 1e-5 absolute, BLUPs 1e-3 relative
#   to their standard deviation. A component at the floor on both sides is
#   compared as two numerical zeros.
# ==============================================================================
suppressPackageStartupMessages({
  library(remlax)
  library(Matrix)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# Root of published_data/, found from the script being run.
pd_root <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) return(normalizePath(file.path(dirname(f[1]), ".."), mustWork = TRUE))
  normalizePath("published_data", mustWork = FALSE)
}

# Downloaded files live here and are not versioned (see .gitignore).
pd_cache <- function(...) {
  d <- file.path(pd_root(), "data_cache", ...)
  dir.create(dirname(d), recursive = TRUE, showWarnings = FALSE)
  d
}

.CHECKS <- new.env()
.CHECKS$rows <- list()
.CHECKS$dataset <- NA_character_
.CHECKS$source <- NA_character_

start_checks <- function(dataset, source) {
  .CHECKS$rows <- list()
  .CHECKS$dataset <- dataset
  .CHECKS$source <- source
  cat(sprintf("\n######## %s (%s)\n", dataset, source))
  cat(sprintf("remlax %s | R %s.%s\n", utils::packageVersion("remlax"),
              R.version$major, R.version$minor))
  for (p in c("asreml", "sommer", "gaston", "rrBLUP", "statgenGWAS", "BGLR", "enhancer"))
    if (requireNamespace(p, quietly = TRUE))
      cat(sprintf("  %-12s %s\n", p, as.character(utils::packageVersion(p))))
  invisible(NULL)
}

# Relative difference with a floor: two values below `floor_` are two zeros.
relc <- function(a, b, floor_ = 1e-6) {
  if (length(a) != 1 || length(b) != 1 || !is.finite(a) || !is.finite(b)) return(Inf)
  if (max(abs(a), abs(b)) < floor_) return(0)
  abs(a - b) / max(abs(b), floor_)
}

# One row of the checks table.
#   metric: "rel" (relc), "abs" (|a - b|), "value" (a is already the
#   discrepancy, b ignored), "ge" (pass if a >= tol, e.g. a correlation)
add_check <- function(model, quantity, software, remlax, reference = NA_real_,
                      tol = 1e-3, metric = c("rel", "abs", "value", "ge"),
                      kind = c("software", "published", "context"), note = "") {
  metric <- match.arg(metric); kind <- match.arg(kind)
  diff <- switch(metric,
                 rel   = relc(remlax, reference),
                 abs   = if (is.finite(remlax) && is.finite(reference)) abs(remlax - reference) else Inf,
                 value = remlax,
                 ge    = remlax)
  pass <- if (kind == "context") NA else
    if (metric == "ge") isTRUE(remlax >= tol) else isTRUE(diff <= tol)
  row <- data.frame(dataset = .CHECKS$dataset, source = .CHECKS$source,
                    model = model, quantity = quantity, software = software,
                    kind = kind, remlax = remlax, reference = reference,
                    metric = metric, diff = diff, tol = tol, pass = pass,
                    note = note, stringsAsFactors = FALSE)
  .CHECKS$rows[[length(.CHECKS$rows) + 1L]] <- row
  cat(sprintf("  %-10s %-26s %-34s %-10s %14.8g %14.8g  d=%9.2e  %s\n",
              kind, model, quantity, software, remlax, reference, diff,
              if (is.na(pass)) "ctx" else if (pass) "OK" else "FAIL"))
  invisible(row)
}

# BLUP agreement: correlation and largest difference relative to the SD of
# the reference BLUPs. Both vectors must be aligned on the same levels.
add_blup_checks <- function(model, software, u_rx, u_ref, tol = 1e-3, note = "") {
  stopifnot(length(u_rx) == length(u_ref))
  ok <- is.finite(u_rx) & is.finite(u_ref)
  r <- stats::cor(u_rx[ok], u_ref[ok])
  d <- max(abs(u_rx[ok] - u_ref[ok])) / stats::sd(u_ref[ok])
  add_check(model, "BLUP correlation", software, r, NA, tol = 1 - 1e-6,
            metric = "ge", note = note)
  add_check(model, "BLUP max |diff| / sd", software, d, NA, tol = tol,
            metric = "value", note = note)
}

write_checks <- function(file = NULL) {
  tab <- do.call(rbind, .CHECKS$rows)
  if (is.null(file)) {
    dir.create(file.path(pd_root(), "results"), showWarnings = FALSE, recursive = TRUE)
    file <- file.path(pd_root(), "results", paste0(.CHECKS$dataset, ".csv"))
  }
  utils::write.csv(tab, file, row.names = FALSE)
  counted <- tab[tab$kind != "context", ]
  cat(sprintf("\n%s: %d checks counted, %d pass, %d fail; %d context rows -> %s\n",
              .CHECKS$dataset, nrow(counted), sum(counted$pass), sum(!counted$pass),
              sum(tab$kind == "context"), file))
  invisible(tab)
}

# asreml stops on its own convergence rule; a few update() calls bring it to
# the optimum when the first call stops early.
stab <- function(a, k = 5) {
  for (i in seq_len(k)) a <- update(a, trace = FALSE)
  a
}

# VanRaden (2008) relationship matrix from a 0/1/2 dose matrix, with a small
# ridge so that asreml accepts it as a G-structure. The SAME matrix is passed
# to every program: the comparison then tests the solvers, not the builders.
vanraden <- function(M012, ridge = 1e-4) {
  p <- colMeans(M012, na.rm = TRUE) / 2
  keep <- p > 0 & p < 1
  M <- M012[, keep, drop = FALSE]
  p <- p[keep]
  W <- sweep(M, 2, 2 * p)
  W[is.na(W)] <- 0
  G <- tcrossprod(W) / (2 * sum(p * (1 - p)))
  G <- G + diag(ridge, nrow(G))
  dimnames(G) <- list(rownames(M012), rownames(M012))
  G
}

# asreml wants the inverse in sparse three-column form, or the matrix itself
# with attribute INVERSE = FALSE. The latter is the least error-prone.
as_asreml_G <- function(G) {
  attr(G, "INVERSE") <- FALSE
  G
}

# Shorthand used in every script.
rx_cpu <- function(...) rx_reml(..., backend = "cpu", verbose = FALSE)
