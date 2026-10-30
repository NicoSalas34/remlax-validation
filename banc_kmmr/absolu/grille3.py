"""Banc absolu avec remlax v13 et trois jeux simules par point (graines g1 =
nom d'origine, g2, g3 : suffixe _g2, _g3 ; sim.R tire la graine du nom, K reste
la meme pour un meme nombre de genotypes).

CPU  : kmmr (g2, g3 ; g1 deja mesure), remlax mme (auto), V_ai (n <= 32000),
       V_lbfgsb (n <= 8000 : il ne sert qu'a comparer les optimiseurs),
       mme_dense, mme_creux et mme_hybride (sauf us(8), us(12) et l'essai a
       10 %, qui depassent 6 h sur ces voies)
A100 : kmmr (g2, g3), remlax mme, V_ai (n <= 32000), mme_dense
"""
import contextlib, importlib.util, io, random
spec = importlib.util.spec_from_file_location("g", "grille.py"); g = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(g)
DUR = ("sim_q1000_t8_r1_sus_x0_o100", "sim_q1000_t12_r1_sus_x0_o100", "sim_q2000_t10_r1_sfa2_x0_o10")
jeux = [(c, c) for c in g.cas] + [(c, c + s) for s in ("_g2", "_g3") for c in g.cas]
open("cas_g23.txt", "w").write("\n".join(n for c, n in jeux if n != c) + "\n")
for cible in ("cpu", "gpu"):
    l = []
    for c, nom in jeux:
        n = g.n_obs(c)
        if nom != c:
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
    random.Random(13 + len(cible)).shuffle(l)
    open(f"taches3_{cible}.txt", "w").write("\n".join(l) + "\n")
    print(cible, len(l), "taches")
