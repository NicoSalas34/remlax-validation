"""Grille du banc absolu : ecrit cas.txt (un cas par ligne) et taches_cpu.txt,
taches_gpu.txt (cas outil voie tol), une ligne par ajustement.

Axes
  n_fixeC  temps contre n a dim(C) fixe : GBLUP, q = 2000, r = 1 a 64
  n_q      temps contre n a q = n/4 : GBLUP, q = 500 a 16000
  p_t      temps contre t : us(t), q = 1000, une observation par couple
  p_struct temps contre p a donnees egales : t = 8, diag, fa1, fa2, fa3, us
  C_pred   dim(C) >> n par l'ensemble de prediction : q = 1000 phenotypes,
           1000 a 16000 genotypes dans la parente
  C_clair  dim(C) >> n par les essais clairsemes : t = 10, q = 2000, fa2,
           100 a 10 % des couples observes
La forme V (voie dense) n'est lancee que si V tient : n <= 16000 sur CPU,
n <= 32000 sur A100.
"""
def nom(q, t, r, s, x=0, o=100):
    return f"sim_q{q}_t{t}_r{r}_s{s}_x{x:g}_o{o:g}"

axes = {}
axes["n_fixeC"] = [nom(2000, 1, r, "id") for r in (1, 2, 4, 8, 16, 32, 64)]
axes["n_q"] = [nom(q, 1, 4, "id") for q in (500, 1000, 2000, 4000, 8000, 16000)]
axes["p_t"] = [nom(1000, t, 1, "us" if t > 1 else "id") for t in (1, 2, 3, 5, 8, 12)]
axes["p_struct"] = [nom(1000, 8, 1, s) for s in ("diag", "fa1", "fa2", "fa3", "us")]
axes["C_pred"] = [nom(1000, 1, 1, "id", x) for x in (0, 1, 3, 7, 15)]
axes["C_clair"] = [nom(2000, 10, 1, "fa2", 0, o) for o in (100, 50, 20, 10)]

cas = []
for a, l in axes.items():
    for c in l:
        if c not in cas:
            cas.append(c)

def n_obs(c):
    d = dict((k[0], k[1:]) for k in c.split("_")[1:])
    q, t, r, o = int(d["q"]), int(d["t"]), int(d["r"]), float(d["o"])
    return q * r * max(1, round(o * t / 100))

with open("cas.txt", "w") as f:
    f.write("\n".join(cas) + "\n")
with open("axes.csv", "w") as f:
    f.write("axe,cas\n")
    for a, l in axes.items():
        for c in l:
            f.write(f"{a},{c}\n")
for cible, nmax in (("cpu", 16000), ("gpu", 32000)):
    # les plus gros d'abord : ils fixent la duree de la campagne
    lignes = []
    for c in sorted(cas, key=n_obs, reverse=True):
        lignes.append(f"{c} kmmr - 0.000001")
        lignes.append(f"{c} remlax mme -")
        if n_obs(c) <= nmax:
            lignes.append(f"{c} remlax dense -")
    with open(f"taches_{cible}.txt", "w") as f:
        f.write("\n".join(lignes) + "\n")
    print(cible, len(lignes), "taches")
with open("taches_asreml.txt", "w") as f:
    f.write("\n".join(f"{c} asreml - -" for c in sorted(cas, key=n_obs)) + "\n")
print(len(cas), "cas")
for c in cas:
    print(c, n_obs(c))
