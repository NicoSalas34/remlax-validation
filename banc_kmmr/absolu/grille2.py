"""Voies entierement fixees pour construire la regle de choix du moteur.
Chaque voie fixe moteur, inversion et optimiseur ; une seule execution (froid),
la regle se lit sur le cout par evaluation et le nombre d'evaluations.
  CPU  : V_ai, V_lbfgsb (n <= 32000), mme_dense, mme_creux, mme_hybride
  A100 : V_ai (n <= 32000), mme_dense (les inversions creuse et hybride
         factorisent sur CPU, l'A100 n'y change rien)
"""
import importlib.util
spec = importlib.util.spec_from_file_location("g", "grille.py"); g = importlib.util.module_from_spec(spec)
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(g)
cas = sorted(g.cas, key=g.n_obs, reverse=True)
for cible, voies in (("cpu", ["V_ai", "V_lbfgsb", "mme_dense", "mme_creux", "mme_hybride"]), ("gpu", ["V_ai", "mme_dense"])):
    l = [f"{c} remlax {v} -" for c in cas for v in voies if not (v.startswith("V_") and g.n_obs(c) > 32000)]
    open(f"taches2_{cible}.txt", "w").write("\n".join(l) + "\n")
    print(cible, len(l))
