"""Jeux _g4 et _g5 pour les cas multicaracteres (t > 1), apres le correctif de
sim.R : les _g2 et _g3 multicaracteres de grille3.py etaient des copies de g1.
Memes voies que grille3.py, kmmr compris."""
import contextlib, importlib.util, io, random, re
spec = importlib.util.spec_from_file_location("g", "grille.py"); g = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(g)
DUR = ("sim_q1000_t8_r1_sus_x0_o100", "sim_q1000_t12_r1_sus_x0_o100", "sim_q2000_t10_r1_sfa2_x0_o10")
multi = [c for c in g.cas if int(re.search(r"_t(\d+)_", c).group(1)) > 1]
jeux = [(c, c + s) for s in ("_g4", "_g5") for c in multi]
open("cas_g45.txt", "w").write("\n".join(n for c, n in jeux) + "\n")
for cible in ("cpu", "gpu"):
    l = []
    for c, nom in jeux:
        n = g.n_obs(c)
        l.append(f"{nom} kmmr - 0.000001")
        l.append(f"{nom} remlax mme -")
        if n <= 32000:
            l.append(f"{nom} remlax V_ai -")
        l.append(f"{nom} remlax mme_dense -")
        if cible == "cpu":
            if n <= 8000 and c not in DUR:
                l.append(f"{nom} remlax V_lbfgsb -")
            if c not in DUR:
                l.append(f"{nom} remlax mme_creux -")
                l.append(f"{nom} remlax mme_hybride -")
    random.Random(17 + len(cible)).shuffle(l)
    open(f"taches4_{cible}.txt", "w").write("\n".join(l) + "\n")
    print(cible, len(l), "taches")
