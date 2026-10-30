"""Fully fixed routes to build the engine choice rule.
Each route fixes engine, inversion and optimiser. One run only (cold): the
rule is read from the cost per evaluation and the number of evaluations.
  CPU  : V_ai, V_lbfgsb (n <= 32000), mme_dense, mme_sparse, mme_hybrid
  A100 : V_ai (n <= 32000), mme_dense (the sparse and hybrid inversions
         factorise on CPU, so the A100 changes nothing there)
The stored taches2_*.txt lists use the old route names mme_creux and
mme_hybride, which executer.R still accepts.
"""
import importlib.util
spec = importlib.util.spec_from_file_location("g", "grille.py"); g = importlib.util.module_from_spec(spec)
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(g)
cases = sorted(g.cases, key=g.n_obs, reverse=True)
for target, routes in (("cpu", ["V_ai", "V_lbfgsb", "mme_dense", "mme_sparse", "mme_hybrid"]), ("gpu", ["V_ai", "mme_dense"])):
    l = [f"{c} remlax {v} -" for c in cases for v in routes if not (v.startswith("V_") and g.n_obs(c) > 32000)]
    open(f"tasks2_{target}.txt", "w").write("\n".join(l) + "\n")
    print(target, len(l))
