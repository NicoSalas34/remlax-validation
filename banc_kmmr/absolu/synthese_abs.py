"""Synthese du banc absolu : tableaux et figures (temps contre n, contre p,
contre dim(C)/n), lus dans les CSV de genotoul et de meso.

    python3 synthese_abs.py <dossier des CSV> <dossier de sortie>

Le dossier des CSV contient abs_cpu.csv, abs_gpu.csv (genotoul), abs_meso.csv
(meso) et les *_delais.csv, plus axes.csv (grille.py).

Temps retenus (secondes murales vues de R) :
  kmmr    appel kmmr + inversion de K et valeurs de depart faites dans R
  asreml  appel asreml + inversion de K faite dans R
  remlax  rx_fit complet, demarrage du serveur Python compris ; voie mme
          "chaud" (cache de compilation garni), "froid" garde a part
"""
import collections
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.expanduser("~/.claude/skills/figure-style"))

RES, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)


def lire(f):
    f = os.path.join(RES, f)
    if not os.path.exists(f):
        return []
    with open(f) as fh:
        return list(csv.DictReader(fh))


def num(x):
    try:
        v = float(x)
        return None if v != v else v
    except (TypeError, ValueError):
        return None


def nom_cas(c):
    d = dict((k[0], k[1:]) for k in c.split("_")[1:])
    q, t, r, x, o = int(d["q"]), int(d["t"]), int(d["r"]), float(d["x"]), float(d["o"])
    q_tot = int(round(q * (1 + x)))
    n = q * r * max(1, round(o * t / 100))
    return dict(q=q, t=t, r=r, s=d["s"], x=x, o=o, q_tot=q_tot, n=n, dimC=q_tot * t + t)


AXES = {}
for r in lire("axes.csv"):
    AXES.setdefault(r["axe"], []).append(r["cas"])

# --- mesures : une entree par (materiel, voie, cas) ---------------------------
# voies : asreml, kmmr, remlax (mme chaud), remlax_froid, remlax_V (forme V)
mes = {}
npar = {}
ll = {}
brut = collections.defaultdict(list)     # (materiel, voie, cas) -> une entree par execution


def voie_remlax(r):
    """voie de synthese d'une ligne remlax ; l'etiquette porte la repetition
    (r2_, r3_) et, pour la voie mme, le passage froid ou chaud."""
    e = r["etiquette"]
    if r["voie"] == "mme":
        return "remlax" if e.endswith("chaud") else "remlax_froid"
    return {"dense": "remlax_Veco", "V_ai": "remlax_V", "V_lbfgsb": "V_lbfgsb", "mme_dense": "mme_dense",
            "mme_creux": "mme_creux", "mme_hybride": "mme_hybride"}.get(r["voie"])


# asreml n'a de licence que sur meso : il y tourne sur un noeud RESERVE (abs_meso_x),
# avec kmmr et remlax en controle. Les premieres passes meso, sur noeuds partages
# surcharges, sont ecartees (abs_meso, abs_meso2).
FICHIERS = (("abs_cpu", "cpu"), ("abs_gpu", "a100"), ("abs_cpu_fixe", "cpu"), ("abs_gpu_fixe", "a100"),
            ("abs_cpu_rep", "cpu"), ("abs_gpu_rep", "a100"), ("abs_meso_x13", "meso"),
            ("abs_cpu_v13", "cpu"), ("abs_gpu_v13", "a100"), ("abs_cpu_v14", "cpu"), ("abs_gpu_v14", "a100"))
# remlax : une seule version par synthese (REMLAX_VERSION=v13 lit les seuls
# fichiers *_v13 ; vide = les passes v10). kmmr et asreml : toutes les passes.
# v14 ne modifie que la face du cone us. La passe v14 refait les cas us et
# complete les voies manquantes ; le reste est lu dans les fichiers v13.
VERSION_RX = os.environ.get("REMLAX_VERSION", "v14")
import re as _re
def version_de(fichier):
    if fichier.endswith("_v14"):
        return "v14"
    return "v13" if fichier.endswith(("_v13", "_x13")) else ""


def jeu_canonique(c):
    """trois jeux par point : g1, g2, g3 si t = 1 ; g1, g4, g5 si t > 1 (les g2 et
    g3 multicaracteres etaient des copies de g1, bug de graine du 2026-10-27)."""
    m = _re.search(r"_(g\d+)$", c)
    g = m.group(1) if m else "g1"
    t = int(_re.search(r"_t(\d+)_", c).group(1))
    return g in (("g1", "g2", "g3") if t == 1 else ("g1", "g4", "g5"))


def garder_rx(fichier, mat_g, r):
    """la ligne remlax appartient-elle a la version retenue ?"""
    v = version_de(fichier)
    if VERSION_RX != "v14":
        return v == VERSION_RX
    if not jeu_canonique(r["cas"]):
        return False
    if v == "v14":
        return True
    # un jeu refait en v14 pour une voie ne garde pas sa mesure v13
    return v == "v13" and (mat_g, r["voie"], r["cas"]) not in CLES_V14


def cas_de_base(c):
    """sim_..._g2 -> sim_... : les jeux g1, g2, g3 d'un meme point sont regroupes."""
    return _re.sub(r"_g\d+$", "", c)
CLES_V14 = set()
if VERSION_RX == "v14":
    for _f, _m in FICHIERS:
        if version_de(_f) == "v14":
            for _r in lire(_f + ".csv"):
                if _r.get("outil") == "remlax":
                    CLES_V14.add((_m, _r.get("voie"), _r.get("cas")))


for fichier, mat_g in FICHIERS:
    ch = lire(fichier + "_chrono.csv")
    t_opt, extra, t_iter = {}, collections.Counter(), {}
    for c in ch:
        cle = (c["cas"], c["outil"], c["voie"], c["etiquette"])
        if c["outil"] == "remlax" and c["etape"] == "optimisation":
            t_opt[cle] = num(c["secondes"])
        if (c["outil"] == "remlax" and c["etape"] in ("optimisation", "finition")) or \
           (c["outil"] == "kmmr" and c["etape"] in ("premiere_iteration", "iterations_suivantes")):
            t_iter[cle] = t_iter.get(cle, 0) + (num(c["secondes"]) or 0)
        if c["outil"] in ("kmmr", "asreml") and c["etape"] in ("inverse_K", "valeurs_depart"):
            extra[cle] += num(c["secondes"]) or 0
    for r in lire(fichier + ".csv"):
        # ligne tronquee par deux ecritures simultanees : ignoree
        if not (r.get("cas") or "").startswith("sim_") or r.get("outil") not in ("remlax", "kmmr", "asreml"):
            continue
        voie = voie_remlax(r) if r["outil"] == "remlax" else r["outil"]
        if voie is None:
            continue
        if r["outil"] == "remlax" and not garder_rx(fichier, mat_g, r):
            continue
        cle = (r["cas"], r["outil"], r["voie"], r["etiquette"])
        ok = (r.get("erreur") or "") == "" and num(r.get("logLik")) is not None and r.get("cas")
        s_ = num(r["secondes"])
        if ok and r["outil"] != "remlax":
            s_ += extra.get(cle, 0)
        ev = num(r.get("eval_s"))
        brut[(mat_g, voie, cas_de_base(r["cas"]))].append(dict(
            # vraisemblance comparee entre outils sur le premier jeu seulement :
            # g2 et g3 sont d'autres donnees
            s=s_ if ok else None, ll=num(r["logLik"]) if ok and cas_de_base(r["cas"]) == r["cas"] else None, it=r["n_iter"], err=(r.get("erreur") or "")[:120],
            moteur=r.get("moteur", ""), inv=r.get("inversion", ""), ev=ev if ok else None,
            nev=((t_opt[cle] / ev) if ok and ev and t_opt.get(cle) else None) if r["outil"] == "remlax"
                else (num(r["n_iter"]) if ok else None), noeud=r["noeud"],
            titer=t_iter.get(cle) if ok else None, niter=num(r["n_iter"]) if ok else None))
        if r["outil"] == "remlax" and num(r["n_par"]) is not None:
            npar[cas_de_base(r["cas"])] = int(num(r["n_par"]))
    for r in lire(fichier + "_delais.csv"):
        voie = voie_remlax(r) if r["outil"] == "remlax" else r["outil"]
        if voie is None:
            continue
        if r["outil"] == "remlax" and not garder_rx(fichier, mat_g, r):
            continue
        brut[(mat_g, voie, cas_de_base(r["cas"]))].append(dict(s=None, ll=None, it="", err="delai %s s" % r["delai"],
                                                   delai=float(r["delai"])))
        if voie == "remlax_froid":     # le froid a depasse : le chaud n'a pas tourne
            brut[(mat_g, "remlax", cas_de_base(r["cas"]))].append(dict(s=None, ll=None, it="", err="delai",
                                                           delai=float(r["delai"])))


def med(v):
    v = [x for x in v if x is not None]
    return float(np.median(v)) if v else None


for k, runs in brut.items():
    oks = [x for x in runs if x["s"] is not None]
    if oks:
        mes[k] = dict(s=med([x["s"] for x in oks]), smin=min(x["s"] for x in oks), smax=max(x["s"] for x in oks),
                      n_rep=len(oks), ll=med([x["ll"] for x in oks]), it=oks[0]["it"], err="",
                      moteur=oks[0]["moteur"], inv=oks[0]["inv"], ev=med([x["ev"] for x in oks]),
                      nev=med([x["nev"] for x in oks]), n_echecs=len(runs) - len(oks),
                      titer=med([x.get("titer") for x in oks]), niter=med([x.get("niter") for x in oks]))
        if mes[k]["ll"] is not None:
            ll.setdefault(k[2], []).append((k[0], k[1], mes[k]["ll"]))
    else:
        d = next((x for x in runs if x.get("delai")), None)
        mes[k] = dict(s=None, ll=None, it="", err=runs[0]["err"], n_rep=0, n_echecs=len(runs),
                      **({"delai": d["delai"]} if d else {}))

CAS = []
for l in AXES.values():
    for c in l:
        if c not in CAS:
            CAS.append(c)

# --- accord des vraisemblances ----------------------------------------------
# meilleure vraisemblance de chaque cas, et ecart de chaque voie a celle-ci
best = {c: max(v for _, _, v in ll[c]) for c in ll}


def ecart(mat, voie, c):
    m = mes.get((mat, voie, c))
    if not m or m["ll"] is None or c not in best:
        return None
    return m["ll"] - best[c]


# --- tableaux ---------------------------------------------------------------
VOIES = [("meso", "asreml", "asreml"), ("meso", "kmmr", "kmmr"), ("meso", "remlax", "remlax"),
         ("cpu", "kmmr", "kmmr"), ("cpu", "remlax", "remlax"), ("cpu", "remlax_V", "remlax V"),
         ("a100", "kmmr", "kmmr"), ("a100", "remlax", "remlax"), ("a100", "remlax_V", "remlax V")]
ENTETES = [("meso, 8 cœurs", 3), ("genotoul, 8 cœurs", 3), ("A100", 3)]


def f_s(x):
    if x is None:
        return ""
    return "%.0f" % x if x >= 100 else ("%.1f" % x if x >= 1 else "%.2f" % x)


def cellule(mat, voie, c):
    m = mes.get((mat, voie, c))
    if m is None:
        return '<td class="vide"></td>'
    if m["s"] is None:
        lab = "délai" if m.get("delai") else "échec"
        return '<td class="echec" title="%s">%s</td>' % (m["err"].replace('"', "'"), lab)
    e = ecart(mat, voie, c)
    cls = ' class="moins"' if e is not None and e < -0.05 else ""
    note = '<span class="ll">%+.2f</span>' % e if e is not None and abs(e) > 0.05 else ""
    return "<td%s>%s%s</td>" % (cls, f_s(m["s"]), note)


TITRES = {"n_fixeC": "Temps contre n, dim(C) fixée (GBLUP, q = 2000)",
          "n_q": "Temps contre n, q = n/4 (GBLUP)",
          "p_t": "Temps contre le nombre d'environnements (us(t), q = 1000, une observation par couple)",
          "p_struct": "Temps contre p à données égales (t = 8, q = 1000)",
          "C_pred": "dim(C) >> n : génotypes à prédire dans la parenté (q = 1000 phénotypés)",
          "C_clair": "dim(C) >> n : essais clairsemés (t = 10, q = 2000, fa2)"}


def descr(c):
    a = nom_cas(c)
    return a, "n = %d, q<sub>K</sub> = %d, t = %d, %s, dim C = %d, p = %s" % (
        a["n"], a["q_tot"], a["t"], a["s"], a["dimC"], npar.get(c, "?"))


with open(os.path.join(OUT, "tableaux.html"), "w") as f:
    for axe, liste in AXES.items():
        f.write('<h3>%s</h3>\n<div class="defile"><table class="banc"><thead><tr><th rowspan="2">cas</th>' % TITRES[axe])
        f.write("".join('<th colspan="%d">%s</th>' % (k, h) for h, k in ENTETES) + "</tr><tr>")
        f.write("".join("<th>%s</th>" % v[2] for v in VOIES) + "</tr></thead><tbody>\n")
        for c in liste:
            f.write('<tr><th scope="row">%s</th>%s</tr>\n' % (descr(c)[1], "".join(cellule(m, v, c) for m, v, _ in VOIES)))
        f.write("</tbody></table></div>\n")

# CSV long, une ligne par mesure
with open(os.path.join(OUT, "mesures.csv"), "w") as f:
    w = csv.writer(f)
    w.writerow(["cas", "materiel", "voie", "n", "q_K", "t", "struct", "dimC", "n_par", "secondes_mediane",
                "secondes_min", "secondes_max", "n_rep", "n_echecs", "eval_s", "n_eval",
                "logLik", "ecart_meilleur", "n_iter", "moteur", "inversion", "erreur"])
    for (mat, voie, c), m in sorted(mes.items(), key=lambda kv: (kv[0][2], kv[0][0], kv[0][1])):
        a = nom_cas(c)
        w.writerow([c, mat, voie, a["n"], a["q_tot"], a["t"], a["s"], a["dimC"], npar.get(c, ""),
                    m["s"], m.get("smin"), m.get("smax"), m.get("n_rep"), m.get("n_echecs"), m.get("ev"), m.get("nev"),
                    m["ll"], ecart(mat, voie, c), m["it"], m.get("moteur", ""), m.get("inv", ""), m["err"]])

# --- figures ----------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
from matplotlib.lines import Line2D
from figure_style import apply_figure_style

apply_figure_style(sizes=(8, 7, 6))


def p_theo(c):
    """parametres de variance du modele simule (genetique + residuelle)"""
    a = nom_cas(c); t = a["t"]
    if t == 1:
        return 2
    if a["s"] == "us":
        g = t * (t + 1) // 2
    elif a["s"] == "diag":
        g = t
    else:
        k_ = int(a["s"][2:]); g = t * k_ - k_ * (k_ - 1) // 2 + t
    return g + t


def fmt_s(v, _):
    if v >= 3600 and abs(v / 3600 - round(v / 3600)) < 1e-9:
        return "%d h" % round(v / 3600)
    return ("%g s" % v).replace(".", ",")


def axe_log_y(ax):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(mt.LogLocator(base=10, subs=(1, 3) if ax.get_ylim()[1] / ax.get_ylim()[0] < 300 else (1,)))
    ax.yaxis.set_major_formatter(mt.FuncFormatter(fmt_s))
    ax.yaxis.set_minor_formatter(mt.NullFormatter())


def ticks_x(ax, vals, labs, log=True):
    if log:
        ax.set_xscale("log")
    ax.set_xticks(vals); ax.set_xticklabels(labs)
    ax.xaxis.set_minor_locator(mt.NullLocator())


def k(x):
    return "%dk" % (x // 1000) if x >= 1000 else str(x)


# --- figures comparatives : une figure par axe, colonnes CPU | A100 ----------
# Couleur = logiciel ou moteur. Le materiel est donne par la colonne, jamais
# par le trait. Noms ecrits au bout des courbes.
COUL = {"remlax": "#1d6b58", "remlax_V": "#4f9fd9", "kmmr": "#d08a2e", "asreml": "#6f6f6f"}
NOMS = {"remlax": "remlax MME", "remlax_V": "remlax forme V", "kmmr": "kmmr", "asreml": "asreml (meso)"}
MAT = {"cpu": ["remlax", "remlax_V", "kmmr", ("meso", "asreml")], "a100": ["remlax", "remlax_V", "kmmr"]}
TITRE_COL = {"cpu": "CPU, 8 cœurs", "a100": "GPU A100"}


def serie(liste, xf, mat, voie):
    pts, fails, delais = [], [], []
    for c in liste:
        m = mes.get((mat, voie, c))
        if not m:
            continue
        if m["s"] is not None:
            pts.append((xf(c), m["s"], m.get("smin"), m.get("smax"), m.get("n_rep", 1)))
        elif m.get("delai"):
            delais.append(xf(c))
        else:
            fails.append(xf(c))
    return sorted(pts), fails, delais


def etiquettes_bout(ax, bouts):
    """Noms au bout des courbes, places une fois les echelles figees, decales
    verticalement (en echelle log) pour ne pas se chevaucher."""
    if not bouts:
        return
    import matplotlib.transforms as mtr
    y0, y1 = [np.log10(v) for v in ax.get_ylim()]
    gap = 0.07 * (y1 - y0)
    bouts = sorted(bouts, key=lambda b: b[1])
    ly = [np.log10(b[1]) for b in bouts]
    for i in range(1, len(ly)):
        ly[i] = max(ly[i], ly[i - 1] + gap)
    debord = ly[-1] - (y1 - 0.02 * (y1 - y0))
    if debord > 0:                       # tout redescendre si le dernier sort du cadre
        ly = [l - debord for l in ly]
    tr = mtr.offset_copy(ax.transData, fig=ax.figure, x=5, y=0, units="points")
    for (x, y, txt, col), l in zip(bouts, ly):
        ax.text(x, 10 ** l, txt, color=col, fontsize=7, va="center", ha="left", transform=tr, clip_on=False)


def panneau(ax, liste, xf, mat, xlog, ticks):
    bouts = []
    for v in MAT[mat]:
        m_, voie = (v if isinstance(v, tuple) else (mat, v))
        pts, fails, delais = serie(liste, xf, m_, voie)
        col = COUL[voie]; lw = 2.0 if voie == "remlax" else 1.3
        if pts:
            x = [p[0] for p in pts]; y = [p[1] for p in pts]
            ax.plot(x, y, "-o", color=col, lw=lw, ms=3.8, zorder=3 if voie == "remlax" else 2)
            for p in pts:
                if p[4] > 1 and p[3] > p[2] * 1.02:
                    ax.plot([p[0], p[0]], [p[2], p[3]], color=col, lw=1.0, alpha=0.55, zorder=1)
            bouts.append((x[-1], y[-1], NOMS[voie], col))
        for xx, mk in ((fails, "x"), (delais, "^")):
            if xx:
                ax.plot(xx, [0.955] * len(xx), ls="none", marker=mk, ms=6, mew=1.5, color=col, mfc="white",
                        transform=ax.get_xaxis_transform(), clip_on=False, zorder=4)
    vals, labs = ticks
    if xlog:
        ax.set_xscale("log")
    ax.set_xticks(vals); ax.set_xticklabels(labs); ax.xaxis.set_minor_locator(mt.NullLocator())
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(mt.FuncFormatter(fmt_s)); ax.yaxis.set_minor_formatter(mt.NullFormatter())
    ax.grid(axis="y", which="major", color="0.9", lw=0.6, zorder=0)
    if xlog:
        ax.set_xlim(min(vals) / 1.25, max(vals) * 2.6)      # place pour les noms a droite
    else:
        ax.set_xlim(min(vals) - 0.06 * (max(vals) - min(vals)), max(vals) + 0.45 * (max(vals) - min(vals)))
    return bouts


def figure_axe(nom, lignes, titre):
    """lignes : [(axe, xf, xlabel, titre de ligne, xlog, (vals, labs))]. Un bloc
    (sous-figure) par ligne, titre en tete du bloc ; colonnes CPU | A100."""
    fig = plt.figure(figsize=(7.2, 3.15 * len(lignes) + 0.5), layout="constrained")
    blocs = fig.subfigures(len(lignes), 1, hspace=0.04)
    blocs = np.atleast_1d(blocs)
    tous = []
    for i, ((axe, xf, xl, tl, xlog, ticks), bloc) in enumerate(zip(lignes, blocs)):
        bloc.suptitle(tl, x=0.01, ha="left", fontsize=8)
        axs = bloc.subplots(1, 2, sharey=True)
        liste = AXES.get(axe, [])
        for j, mat in enumerate(("cpu", "a100")):
            ax = axs[j]
            tous.append((ax, panneau(ax, liste, xf, mat, xlog, ticks)))
            ax.set_xlabel(xl)
            ax.set_title("%s   %s" % ("abcdefgh"[2 * i + j], TITRE_COL[mat]), loc="left", fontsize=8)
        axs[0].set_ylabel("temps de l'ajustement")
    fig.suptitle(titre, x=0.01, ha="left", fontsize=8.5, fontweight="bold")
    fig.get_layout_engine().set(w_pad=0.05, h_pad=0.06, rect=(0, 0.025, 1, 0.97))
    fig.canvas.draw()                    # echelles et mise en page figees
    for ax, bouts in tous:
        etiquettes_bout(ax, bouts)
    fig.text(0.01, 0.002, "× échec de l'ajustement   △ délai de 6 h dépassé   │ étendue des répétitions "
             "(médiane tracée). kmmr et remlax sur genotoul ; asreml sur un nœud meso réservé.",
             fontsize=6.5, color="0.35", ha="left", va="bottom")
    fig.savefig(os.path.join(OUT, nom + ".png"), dpi=200)
    fig.savefig(os.path.join(OUT, nom + ".pdf"))
    plt.close(fig)


A = AXES
fn = lambda c: nom_cas(c)["n"]
fr = lambda c: nom_cas(c)["dimC"] / nom_cas(c)["n"]
tk_n = lambda axe: (sorted({fn(c) for c in A.get(axe, [])}), [k(v) for v in sorted({fn(c) for c in A.get(axe, [])})])
tk_r = lambda axe: (sorted({fr(c) for c in A.get(axe, [])}), ["%g" % round(v, 1) for v in sorted({fr(c) for c in A.get(axe, [])})])
ps = sorted(A.get("p_struct", []), key=p_theo)
for nom in ("fig_comparatif",):
    f = os.path.join(OUT, nom + ".png")
    if os.path.exists(f):
        os.remove(f)
figure_axe("fig_temps_n", [
    ("n_fixeC", fn, "observations n", "dim(C) fixée : un seul caractère, q = 2000 génotypes, de 1 à 64 répétitions", True, tk_n("n_fixeC")),
    ("n_q", fn, "observations n", "dim(C) croît avec n : un seul caractère, 4 répétitions, q = n/4", True, tk_n("n_q"))],
    "Temps de l'ajustement contre le nombre d'observations (GBLUP, parenté dense)")
figure_axe("fig_temps_p", [
    ("p_t", p_theo, "paramètres de variance p", "us(t) pour t = 1 à 12 environnements ; n et dim(C) croissent avec t", True,
     (sorted({p_theo(c) for c in A.get("p_t", [])}), [str(v) for v in sorted({p_theo(c) for c in A.get("p_t", [])})])),
    ("p_struct", p_theo, "structure génétique (paramètres p)", "t = 8 environnements, mêmes données ; seule la structure change", False,
     ([p_theo(c) for c in ps], ["%s\n%d" % (nom_cas(c)["s"], p_theo(c)) for c in ps]))],
    "Temps de l'ajustement contre le nombre de paramètres de variance")
figure_axe("fig_temps_C", [
    ("C_pred", fr, "dim(C) / n", "génotypes à prédire : 1000 phénotypés, de 1000 à 16 000 dans la parenté", True, tk_r("C_pred")),
    ("C_clair", fr, "dim(C) / n", "essais clairsemés : t = 10, q = 2000, fa2, de 100 % à 10 % des couples observés", True, tk_r("C_clair"))],
    "Temps de l'ajustement quand C devient grande devant n")

# --- figure de la regle : voies entierement fixees de remlax -----------------
CV = {"remlax_V": "#4f9fd9", "mme_dense": "#1d6b58", "mme_hybride": "#8a6bb8", "mme_creux": "#c0504d"}
NV = {"remlax_V": "forme V (taille n)", "mme_dense": "MME, inversion dense", "mme_hybride": "MME, hybride",
      "mme_creux": "MME, creuse"}
fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.4), sharey="row")
for j, mat in enumerate(("cpu", "a100")):
    ax = axs[0, j]
    for voie in (("remlax_V", "mme_dense") if mat == "a100" else ("remlax_V", "mme_dense", "mme_hybride", "mme_creux")):
        pts = [((nom_cas(c)["n"] if voie == "remlax_V" else nom_cas(c)["dimC"]), mes[(mat, voie, c)]["ev"]) for c in CAS
               if mes.get((mat, voie, c), {}).get("ev")]
        if pts:
            x, y = zip(*pts)
            ax.plot(x, y, "o", ms=3.8, color=CV[voie], alpha=0.85, label=NV[voie])
    # pente 3 (factorisation dense) pour reperer
    xs = np.array([1e3, 3e4]); ax.plot(xs, 0.02 * (xs / 1e3) ** 3 * (0.05 if mat == "a100" else 1), color="0.75",
                                       lw=0.8, ls="--", zorder=0)
    ax.text(3e4, 0.02 * 30 ** 3 * (0.05 if mat == "a100" else 1), " pente 3", color="0.55", fontsize=6.5, va="center")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.yaxis.set_major_formatter(mt.FuncFormatter(fmt_s)); ax.yaxis.set_minor_formatter(mt.NullFormatter())
    ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: k(int(round(v))))); ax.xaxis.set_minor_formatter(mt.NullFormatter())
    ax.set_xlabel("taille de la matrice factorisée (n ou dim C)")
    ax.set_title(TITRE_COL[mat], fontweight="bold")
    ax.grid(color="0.92", lw=0.6, zorder=0)
    if j == 0:
        ax.set_ylabel("temps d'une évaluation")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    ax = axs[1, j]
    pts = []
    for c in CAS:
        v, m_ = mes.get((mat, "remlax_V", c)), mes.get((mat, "mme_dense", c))
        if v and m_ and v.get("ev") and m_.get("ev"):
            pts.append((nom_cas(c)["dimC"] / nom_cas(c)["n"], m_["ev"] / v["ev"]))
    if pts:
        x, y = zip(*pts)
        ax.plot(x, y, "o", ms=4, color="0.25")
    ax.set_xscale("log"); ax.set_yscale("log")
    lo, hi = 1e-4, 300
    ax.set_ylim(lo, hi)
    ax.axhspan(1, hi, color=CV["remlax_V"], alpha=0.10, zorder=0)
    ax.axhspan(lo, 1, color=CV["mme_dense"], alpha=0.08, zorder=0)
    ax.axhline(1, color="0.5", lw=0.8)
    ax.text(0.04, 0.95, "forme V moins chère", transform=ax.transAxes, color=CV["remlax_V"], va="top", fontsize=7)
    ax.text(0.04, 0.05, "MME moins cher", transform=ax.transAxes, color=CV["mme_dense"], va="bottom", fontsize=7)
    ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: ("%g" % v).replace(".", ",")))
    ax.yaxis.set_minor_formatter(mt.NullFormatter())
    ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: ("%g" % v).replace(".", ",")))
    ax.xaxis.set_minor_formatter(mt.NullFormatter())
    ax.set_xlabel("dim(C) / n")
    if j == 0:
        ax.set_ylabel("coût d'une évaluation\nMME dense / forme V")
for i, ax in enumerate(axs.flat):
    ax.text(-0.13 if i % 2 == 0 else -0.05, 1.04, "abcd"[i], transform=ax.transAxes, fontweight="bold", fontsize=9)
fig.text(0.01, 0.005, "Voies entièrement fixées de remlax, médiane des répétitions. Un point = un cas simulé.",
         fontsize=6.5, color="0.35")
fig.tight_layout(rect=(0, 0.03, 1, 1), h_pad=2.0)
fig.savefig(os.path.join(OUT, "fig_regle.png"), dpi=200)
fig.savefig(os.path.join(OUT, "fig_regle.pdf"))
plt.close(fig)

# --- figure decomposee : temps = fixe + n_eval x cout d'une evaluation -------
# remlax : n_eval = temps d'optimisation / cout d'une evaluation complete ;
# kmmr : une evaluation = une iteration AI. Le cout fixe est le reste du temps
# total (demarrage, preparation, compilation, inversion de K, sorties).
LIGNES_DEC = [
    ("n_fixeC", fn, "n", "n, dim(C) fixée (q = 2000)", True),
    ("n_q", fn, "n", "n, q = n/4", True),
    ("p_t", p_theo, "p", "us(t), t = 1 à 12", True),
    ("p_struct", p_theo, "structure", "t = 8, structures", False),
    ("C_pred", fr, "dim(C)/n", "génotypes à prédire", True),
    ("C_clair", fr, "dim(C)/n", "essais clairsemés", True),
]
COMPOS = [("cit", "coût d'une itération", True), ("niter", "nombre d'itérations", False),
          ("fixe", "coût fixe (hors itérations)", True)]


def valeur(m, cle):
    """Grandeurs mesurees : nombre d'iterations, temps passe dans les iterations
    (remlax : optimisation et finition ; kmmr : somme de ses iterations), et le
    reste du temps total."""
    if not m or m["s"] is None or not m.get("titer") or not m.get("niter"):
        return None
    if cle == "niter":
        return m["niter"]
    if cle == "cit":
        return m["titer"] / m["niter"]
    return max(m["s"] - m["titer"], 0.05)


def axe_log_lisible(ax):
    """Graduations 1-2-5 en echelle log, lisibles meme sur moins d'une decade."""
    ax.set_yscale("log")
    lo, hi = ax.get_ylim()
    subs = (1.0, 2.0, 5.0) if hi / lo < 100 else (1.0,)
    ax.yaxis.set_major_locator(mt.LogLocator(base=10, subs=subs))
    ax.yaxis.set_major_formatter(mt.FuncFormatter(fmt_s)); ax.yaxis.set_minor_formatter(mt.NullFormatter())


def figure_decomp(mat, nom):
    fig = plt.figure(figsize=(7.2, 1.95 * len(LIGNES_DEC) + 0.7), layout="constrained")
    axs = fig.subplots(len(LIGNES_DEC), 3)
    for i, (axe, xf, xl, tl, xlog) in enumerate(LIGNES_DEC):
        liste = AXES.get(axe, [])
        if axe == "p_struct":
            liste = sorted(liste, key=p_theo)
            vals = [p_theo(c) for c in liste]; labs = [nom_cas(c)["s"] for c in liste]
        else:
            vals = sorted({xf(c) for c in liste})
            labs = [k(v) if axe.startswith("n_") else ("%g" % round(v, 1)) for v in vals]
        for j, (cle, titre, ylog) in enumerate(COMPOS):
            ax = axs[i, j]
            bouts = []
            for voie in ("remlax", "kmmr"):
                pts, fails = [], []
                for c in liste:
                    m = mes.get((mat, voie, c))
                    v = valeur(m, cle)
                    if v is not None:
                        pts.append((xf(c), v))
                    elif m is not None and m["s"] is None:
                        fails.append(xf(c))
                pts.sort()
                if pts:
                    x, y = zip(*pts)
                    ax.plot(x, y, "-o", color=COUL[voie], lw=2.0 if voie == "remlax" else 1.3, ms=3.2)
                    bouts.append((x[-1], y[-1], NOMS[voie], COUL[voie]))
                if fails:
                    ax.plot(fails, [0.93] * len(fails), "x", ms=5, mew=1.3, color=COUL[voie],
                            transform=ax.get_xaxis_transform(), clip_on=False)
            if xlog:
                ax.set_xscale("log")
            ax.set_xticks(vals); ax.set_xticklabels(labs, fontsize=6); ax.xaxis.set_minor_locator(mt.NullLocator())
            if xlog:
                ax.set_xlim(min(vals) / 1.3, max(vals) * 1.3)
            if ylog:
                ax.set_yscale("log"); ax.autoscale_view(); axe_log_lisible(ax)
            else:
                ax.set_yscale("log"); ax.autoscale_view(); axe_log_lisible(ax)
                ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: "%g" % v))
            ax.grid(axis="y", color="0.92", lw=0.6)
            ax.tick_params(labelsize=6)
            if i == 0:
                ax.set_title(titre, fontsize=8, fontweight="bold")
            if j == 0:
                ax.set_ylabel(tl, fontsize=7)
            if i == len(LIGNES_DEC) - 1 or True:
                ax.set_xlabel(xl, fontsize=6.5, labelpad=1)
    fig.suptitle("Décomposition du temps de l'ajustement, %s : temps = coût fixe + nombre d'itérations × coût d'une itération"
                 % TITRE_COL[mat], x=0.01, ha="left", fontsize=8)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color=COUL[v], lw=2.0 if v == "remlax" else 1.3, marker="o", ms=3.2, label=NOMS[v])
                        for v in ("remlax", "kmmr")] + [Line2D([], [], color="0.3", ls="none", marker="x", label="échec")],
               loc="outside lower center", ncol=3, frameon=False, fontsize=7)
    fig.savefig(os.path.join(OUT, nom + ".png"), dpi=200)
    fig.savefig(os.path.join(OUT, nom + ".pdf"))
    plt.close(fig)


figure_decomp("cpu", "fig_decomposition_cpu")
figure_decomp("a100", "fig_decomposition_a100")

# --- CPU contre GPU ------------------------------------------------------------
# Gain du GPU = temps CPU 8 coeurs / temps A100, meme voie, meme cas, meme
# genotoul. Un point n'existe que si les deux materiels ont abouti.
def gain_gpu(liste, xf, voie):
    pts = []
    for c in liste:
        mc, mg = mes.get(("cpu", voie, c)), mes.get(("a100", voie, c))
        if mc and mg and mc["s"] and mg["s"]:
            pts.append((xf(c), mc["s"] / mg["s"]))
    return sorted(pts)


def fmt_gain(v, _):
    return ("×%g" % v).replace(".", ",")


tk_p = (sorted({p_theo(c) for c in A.get("p_t", [])}), [str(v) for v in sorted({p_theo(c) for c in A.get("p_t", [])})])
LIGNES_CG = [("n_fixeC", fn, "observations n", "dim(C) fixée, q = 2000", True, tk_n("n_fixeC")),
             ("n_q", fn, "observations n", "dim(C) croît avec n, q = n/4", True, tk_n("n_q")),
             ("p_t", p_theo, "paramètres de variance p", "us(t), t = 1 à 12", True, tk_p),
             ("p_struct", p_theo, "structure génétique (paramètres p)", "t = 8, seule la structure change", False,
              ([p_theo(c) for c in ps], ["%s\n%d" % (nom_cas(c)["s"], p_theo(c)) for c in ps])),
             ("C_pred", fr, "dim(C) / n", "génotypes à prédire", True, tk_r("C_pred")),
             ("C_clair", fr, "dim(C) / n", "essais clairsemés, fa2", True, tk_r("C_clair"))]
fig, axs = plt.subplots(3, 2, figsize=(7.2, 8.6), sharey=True, layout="constrained")
tous = []
for ax, (axe, xf, xl, tl, xlog, (vals, labs)), lettre in zip(axs.flat, LIGNES_CG, "abcdef"):
    bouts = []
    for voie in ("remlax", "remlax_V", "kmmr"):
        pts = gain_gpu(A.get(axe, []), xf, voie)
        if pts:
            x, y = zip(*pts)
            ax.plot(x, y, "-o", color=COUL[voie], lw=2.0 if voie == "remlax" else 1.3, ms=3.8,
                    zorder=3 if voie == "remlax" else 2)
            bouts.append((x[-1], y[-1], NOMS[voie], COUL[voie]))
    ax.axhline(1, color="0.45", lw=0.7, ls="-", zorder=1)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(mt.FuncFormatter(fmt_gain)); ax.yaxis.set_minor_formatter(mt.NullFormatter())
    if xlog:
        ax.set_xscale("log")
    ax.set_xticks(vals); ax.set_xticklabels(labs); ax.xaxis.set_minor_locator(mt.NullLocator())
    if xlog:
        ax.set_xlim(min(vals) / 1.25, max(vals) * 2.6)
    else:
        ax.set_xlim(min(vals) - 0.06 * (max(vals) - min(vals)), max(vals) + 0.45 * (max(vals) - min(vals)))
    ax.grid(axis="y", which="major", color="0.9", lw=0.6, zorder=0)
    ax.set_xlabel(xl)
    ax.set_title("%s   %s" % (lettre, tl), loc="left", fontsize=8)
    tous.append((ax, bouts))
for ax in axs[:, 0]:
    ax.set_ylabel("gain du GPU (temps CPU / temps A100)")
fig.suptitle("Gain de l'A100 sur 8 cœurs CPU, voie par voie", x=0.01, ha="left", fontsize=8.5, fontweight="bold")
fig.get_layout_engine().set(w_pad=0.05, h_pad=0.06, rect=(0, 0.03, 1, 0.955))
fig.canvas.draw()
for ax, bouts in tous:
    etiquettes_bout(ax, bouts)
fig.text(0.01, 0.002, "Au-dessus du trait gris, l'A100 est plus rapide. Un point manque quand l'un des deux matériels a "
         "échoué, a été arrêté ou n'a pas encore tourné.", fontsize=6.5, color="0.35", ha="left", va="bottom")
fig.savefig(os.path.join(OUT, "fig_cpu_gpu.png"), dpi=200)
fig.savefig(os.path.join(OUT, "fig_cpu_gpu.pdf"))
plt.close(fig)

# nuage : temps CPU contre temps A100, un point par cas et par voie
fig, ax = plt.subplots(figsize=(4.8, 4.6), layout="constrained")
tout_t = []
for voie in ("kmmr", "remlax_V", "remlax"):
    xy = [(mes[("cpu", voie, c)]["s"], mes[("a100", voie, c)]["s"]) for c in CAS
          if mes.get(("cpu", voie, c)) and mes.get(("a100", voie, c))
          and mes[("cpu", voie, c)]["s"] and mes[("a100", voie, c)]["s"]]
    if xy:
        x, y = zip(*xy); tout_t += list(x) + list(y)
        ax.scatter(x, y, s=16, color=COUL[voie], alpha=0.85, lw=0, label=NOMS[voie], zorder=3)
lo, hi = min(tout_t) / 1.5, max(tout_t) * 1.5
for g, ls in ((1, "-"), (10, "--"), (100, ":")):
    ax.plot([lo, hi], [lo / g, hi / g], color="0.5", lw=0.8, ls=ls, zorder=1)
    ax.text(hi, hi / g, "  ×%d" % g if g > 1 else "  égalité", fontsize=6.5, color="0.4", va="center", ha="left",
            clip_on=False)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
for a_ in (ax.xaxis, ax.yaxis):
    a_.set_major_formatter(mt.FuncFormatter(fmt_s)); a_.set_minor_formatter(mt.NullFormatter())
ax.set_xlabel("temps CPU, 8 cœurs"); ax.set_ylabel("temps GPU A100")
ax.grid(color="0.92", lw=0.6, zorder=0)
ax.legend(frameon=False, fontsize=7, loc="upper left")
ax.set_title("Temps CPU contre temps A100, tous les cas du banc", loc="left", fontsize=8.5, fontweight="bold")
fig.savefig(os.path.join(OUT, "fig_cpu_gpu_nuage.png"), dpi=200)
fig.savefig(os.path.join(OUT, "fig_cpu_gpu_nuage.pdf"))
plt.close(fig)

# --- remlax contre asreml, meme probleme ----------------------------------------
# asreml : noeud meso reserve. remlax MME v14 : CPU 8 coeurs de genotoul, puis A100.
FAM = {"n_fixeC": "GBLUP, un caractère", "n_q": "GBLUP, un caractère", "p_t": "multicaractère (us, fa, diag)",
       "p_struct": "multicaractère (us, fa, diag)", "C_pred": "dim(C) grande devant n", "C_clair": "dim(C) grande devant n"}
CFAM = {"GBLUP, un caractère": "#1d6b58", "multicaractère (us, fa, diag)": "#8a6bb8", "dim(C) grande devant n": "#d08a2e"}
MFAM = {"GBLUP, un caractère": "o", "multicaractère (us, fa, diag)": "s", "dim(C) grande devant n": "^"}
fam_de = {}
for axe, l in AXES.items():
    for c in l:
        fam_de.setdefault(c, FAM.get(axe))
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.9), sharex=True, sharey=True, layout="constrained")
paires = []
for c in CAS:
    ma = mes.get(("meso", "asreml", c))
    if ma and ma["s"]:
        paires.append((c, ma["s"]))
tt = [s_ for _, s_ in paires]
for c, _ in paires:
    for m_ in ("cpu", "a100"):
        mr = mes.get((m_, "remlax", c))
        if mr and mr["s"]:
            tt.append(mr["s"])
lo, hi = min(tt) / 1.6, max(tt) * 1.6
for ax, m_, lettre, tl in ((axs[0], "cpu", "a", "remlax MME sur CPU, 8 cœurs"), (axs[1], "a100", "b", "remlax MME sur GPU A100")):
    for g, ls in ((1, "-"), (10, "--"), (100, ":")):
        ax.plot([lo, hi], [lo / g, hi / g], color="0.55", lw=0.8, ls=ls, zorder=1)
        ax.text(hi, hi / g, " égalité" if g == 1 else " ×%d" % g, fontsize=6.5, color="0.4", va="center",
                ha="left", clip_on=False)
    vus = set()
    for c, ta in paires:
        mr = mes.get((m_, "remlax", c))
        if not (mr and mr["s"]):
            continue
        fa = fam_de.get(c)
        ax.scatter(ta, mr["s"], s=24, marker=MFAM[fa], color=CFAM[fa], lw=0, alpha=0.9, zorder=3,
                   label=fa if fa not in vus else None)
        vus.add(fa)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    for a_ in (ax.xaxis, ax.yaxis):
        a_.set_major_formatter(mt.FuncFormatter(fmt_s)); a_.set_minor_formatter(mt.NullFormatter())
    ax.set_aspect("equal")
    ax.grid(color="0.92", lw=0.6, zorder=0)
    ax.set_xlabel("temps asreml (CPU, 8 cœurs)")
    ax.set_title("%s   %s" % (lettre, tl), loc="left", fontsize=8)
axs[0].set_ylabel("temps remlax MME")
axs[0].legend(frameon=False, fontsize=6.5, loc="upper left")
fig.suptitle("remlax contre asreml sur les mêmes problèmes", x=0.01, ha="left", fontsize=8.5, fontweight="bold")
fig.get_layout_engine().set(w_pad=0.08, rect=(0, 0.06, 1, 1))
fig.text(0.01, 0.005, "Sous la diagonale, remlax est plus rapide. asreml tourne sur un nœud meso réservé, remlax v14 sur "
         "genotoul ; le nœud meso fait tourner remlax\n1,3 à 3 fois plus vite que genotoul, donc le panneau a "
         "désavantage remlax. %d problèmes, ceux qu'asreml a déjà traités." % len(paires),
         fontsize=6.5, color="0.35", ha="left", va="bottom")
fig.savefig(os.path.join(OUT, "fig_remlax_asreml.png"), dpi=200)
fig.savefig(os.path.join(OUT, "fig_remlax_asreml.pdf"))
plt.close(fig)

n_mes = sum(1 for m in mes.values() if m["s"] is not None)
print("mesures abouties :", n_mes, "sur", len(mes))
