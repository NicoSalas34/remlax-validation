# ==============================================================================
# charger_remlax.R -- loads the installed remlax package for the scripts of this
#                     validation repository
# ------------------------------------------------------------------------------
# The validation scripts were born in the remlax repository, where they loaded
# the sources with source("R/remlax.R"), a file since split into several. They
# also test internal functions (.rx_parse_random, .rx_parse_residual, ...).
# This file therefore loads the INSTALLED package and makes its whole namespace
# visible from the global environment. The scripts see exactly what they saw
# with source(), but the code that answers is the tested package, at the
# installed version.
#
# The tested version is printed first. A validation log must state which
# version of remlax produced it.
# ==============================================================================
if (!requireNamespace("remlax", quietly = TRUE))
  stop("the remlax package is not installed: ",
       "remotes::install_github('NicoSalas34/remlax') or R CMD INSTALL <remlax repository>")
suppressPackageStartupMessages(library(remlax))
local({
  ns <- asNamespace("remlax")
  for (nm in ls(ns, all.names = TRUE))
    if (!exists(nm, envir = globalenv(), inherits = FALSE))
      assign(nm, get(nm, envir = ns), envir = globalenv())
})
cat("remlax", as.character(utils::packageVersion("remlax")),
    "installed in", dirname(system.file(package = "remlax")), "\n")
