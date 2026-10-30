"""Absolute benchmark with remlax v13 and three simulated data sets per point
(seeds g1 = original name, g2, g3: suffix _g2, _g3. sim.R draws the seed from
the name, and K stays the same for the same number of genotypes).

CPU  : kmmr (g2, g3; g1 already measured), remlax mme (auto), V_ai (n <= 32000),
       V_lbfgsb (n <= 8000: it only serves to compare the optimisers),
       mme_dense, mme_sparse and mme_hybrid (except us(8), us(12) and the
       10 % trial, which exceed 6 h on these routes)
A100 : kmmr (g2, g3), remlax mme, V_ai (n <= 32000), mme_dense
"""
import contextlib, importlib.util, io, random
spec = importlib.util.spec_from_file_location("g", "grille.py"); g = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(g)
SLOW = ("sim_q1000_t8_r1_sus_x0_o100", "sim_q1000_t12_r1_sus_x0_o100", "sim_q2000_t10_r1_sfa2_x0_o10")
sets = [(c, c) for c in g.cases] + [(c, c + s) for s in ("_g2", "_g3") for c in g.cases]
open("cases_g23.txt", "w").write("\n".join(n for c, n in sets if n != c) + "\n")
for target in ("cpu", "gpu"):
    l = []
    for c, name in sets:
        n = g.n_obs(c)
        if name != c:
            l.append(f"{name} kmmr - 0.000001")
        l.append(f"{name} remlax mme -")
        if n <= 32000:
            l.append(f"{name} remlax V_ai -")
        l.append(f"{name} remlax mme_dense -")
        if target == "cpu":
            if n <= 8000 and c not in SLOW:
                l.append(f"{name} remlax V_lbfgsb -")
            if c not in SLOW:
                l.append(f"{name} remlax mme_sparse -")
                l.append(f"{name} remlax mme_hybrid -")
    random.Random(13 + len(target)).shuffle(l)
    open(f"tasks3_{target}.txt", "w").write("\n".join(l) + "\n")
    print(target, len(l), "tasks")
