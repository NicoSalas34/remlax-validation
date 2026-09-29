# ==============================================================================
# charger_remlax.R -- charge le paquet remlax installe pour les scripts de ce
#                     depot de validation
# ------------------------------------------------------------------------------
# Les scripts de validation sont nes dans le depot remlax, ou ils chargeaient les
# sources par source("R/remlax.R"). Ils testent aussi des fonctions internes
# (rx_python_cmd, rx_sparse_*, ...). On charge donc le paquet INSTALLE et l'on
# rend son espace de noms entier visible depuis l'environnement global : les
# scripts voient exactement ce qu'ils voyaient avec source(), mais c'est bien le
# code du paquet teste, a la version installee, qui repond.
#
# La version testee est imprimee en tete : un journal de validation doit dire
# contre quelle version de remlax il a ete produit.
# ==============================================================================
if (!requireNamespace("remlax", quietly = TRUE))
  stop("le paquet remlax n'est pas installe : ",
       "remotes::install_github('NicoSalas34/remlax') ou R CMD INSTALL <depot remlax>")
suppressPackageStartupMessages(library(remlax))
local({
  ns <- asNamespace("remlax")
  for (nm in ls(ns, all.names = TRUE))
    if (!exists(nm, envir = globalenv(), inherits = FALSE))
      assign(nm, get(nm, envir = ns), envir = globalenv())
})
cat("remlax", as.character(utils::packageVersion("remlax")),
    "installe dans", dirname(system.file(package = "remlax")), "\n")
