"""Grid of the absolute benchmark: writes cases.txt (one case per line) and
tasks_cpu.txt, tasks_gpu.txt (case tool route tol), one line per fit.

Axes
  n_fixeC  time against n at fixed dim(C): GBLUP, q = 2000, r = 1 to 64
  n_q      time against n at q = n/4: GBLUP, q = 500 to 16000
  p_t      time against t: us(t), q = 1000, one observation per pair
  p_struct time against p at equal data: t = 8, diag, fa1, fa2, fa3, us
  C_pred   dim(C) >> n through the prediction set: q = 1000 phenotyped,
           1000 to 16000 genotypes in the kinship
  C_clair  dim(C) >> n through sparse trials: t = 10, q = 2000, fa2,
           100 to 10 % of the pairs observed
The V form (dense route) is run only if V fits: n <= 16000 on CPU,
n <= 32000 on A100.

The axis names (n_fixeC, C_clair, ...) are kept as in the stored axes.csv.
The stored task lists of the 2026-10 campaign carry the old French file names
(cas.txt, taches_*.txt).
"""
def name(q, t, r, s, x=0, o=100):
    return f"sim_q{q}_t{t}_r{r}_s{s}_x{x:g}_o{o:g}"

axes = {}
axes["n_fixeC"] = [name(2000, 1, r, "id") for r in (1, 2, 4, 8, 16, 32, 64)]
axes["n_q"] = [name(q, 1, 4, "id") for q in (500, 1000, 2000, 4000, 8000, 16000)]
axes["p_t"] = [name(1000, t, 1, "us" if t > 1 else "id") for t in (1, 2, 3, 5, 8, 12)]
axes["p_struct"] = [name(1000, 8, 1, s) for s in ("diag", "fa1", "fa2", "fa3", "us")]
axes["C_pred"] = [name(1000, 1, 1, "id", x) for x in (0, 1, 3, 7, 15)]
axes["C_clair"] = [name(2000, 10, 1, "fa2", 0, o) for o in (100, 50, 20, 10)]

cases = []
for a, l in axes.items():
    for c in l:
        if c not in cases:
            cases.append(c)

def n_obs(c):
    d = dict((k[0], k[1:]) for k in c.split("_")[1:])
    q, t, r, o = int(d["q"]), int(d["t"]), int(d["r"]), float(d["o"])
    return q * r * max(1, round(o * t / 100))

with open("cases.txt", "w") as f:
    f.write("\n".join(cases) + "\n")
with open("axes.csv", "w") as f:
    f.write("axis,case\n")
    for a, l in axes.items():
        for c in l:
            f.write(f"{a},{c}\n")
for target, nmax in (("cpu", 16000), ("gpu", 32000)):
    # the largest first: they set the length of the campaign
    lines = []
    for c in sorted(cases, key=n_obs, reverse=True):
        lines.append(f"{c} kmmr - 0.000001")
        lines.append(f"{c} remlax mme -")
        if n_obs(c) <= nmax:
            lines.append(f"{c} remlax dense -")
    with open(f"tasks_{target}.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print(target, len(lines), "tasks")
with open("tasks_asreml.txt", "w") as f:
    f.write("\n".join(f"{c} asreml - -" for c in sorted(cases, key=n_obs)) + "\n")
print(len(cases), "cases")
for c in cases:
    print(c, n_obs(c))
