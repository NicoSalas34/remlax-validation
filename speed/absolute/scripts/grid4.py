"""_g4 and _g5 data sets for the multi-trait cases (t > 1), after the fix of
sim.R: the multi-trait _g2 and _g3 of grille3.py were copies of g1.
Same routes as grille3.py, kmmr included."""
import contextlib, importlib.util, io, random, re
spec = importlib.util.spec_from_file_location("g", "grille.py"); g = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(g)
SLOW = ("sim_q1000_t8_r1_sus_x0_o100", "sim_q1000_t12_r1_sus_x0_o100", "sim_q2000_t10_r1_sfa2_x0_o10")
multi = [c for c in g.cases if int(re.search(r"_t(\d+)_", c).group(1)) > 1]
sets = [(c, c + s) for s in ("_g4", "_g5") for c in multi]
open("cases_g45.txt", "w").write("\n".join(n for c, n in sets) + "\n")
for target in ("cpu", "gpu"):
    l = []
    for c, name in sets:
        n = g.n_obs(c)
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
    random.Random(17 + len(target)).shuffle(l)
    open(f"tasks4_{target}.txt", "w").write("\n".join(l) + "\n")
    print(target, len(l), "tasks")
