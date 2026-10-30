"""Tableaux HTML du rapport de banc, lus dans les CSV de genotoul et d'asreml.

    python3 tableaux_html.py <dossier res genotoul> <banc_asreml.csv> > fragments.html
"""
import csv
import os
import sys

RES, ASR = sys.argv[1], sys.argv[2]
CAS = ["drops2_fa1", "drops2_fa3", "drops1_fa1", "drops1_us", "krause_gblup", "krause_fa2", "krause_us"]
NOM = {"drops2_fa1": "DROPS deux étapes, fa1", "drops2_fa3": "DROPS deux étapes, fa3",
       "drops1_fa1": "DROPS une étape, fa1", "drops1_us": "DROPS une étape, us(5)",
       "krause_gblup": "Krause, GBLUP", "krause_fa2": "Krause, fa2 × 5 traitements",
       "krause_us": "Krause, us(5) × 5 traitements"}


def lire(f):
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


def kmmr(cible):
    """kmmr a la tolerance 1e-6 : temps = appel + inversion de K et valeurs de
    depart faites dans R avant l'appel (premier fichier ou le cas figure)."""
    out = {}
    for v in ("_v3", "_v2", ""):
        f = os.path.join(RES, "banc_%s%s.csv" % (cible, v))
        ch = lire(f.replace(".csv", "_chrono.csv"))
        for r in lire(f):
            if r["outil"] != "kmmr" or r["voie"] != "ai_tol0.000001" or r["cas"] in out or num(r["logLik"]) is None:
                continue
            x = [c for c in ch if c["cas"] == r["cas"] and c["outil"] == "kmmr" and c["voie"] == "ai_tol0.000001"]
            n_runs = max(1, sum(1 for c in x if c["etape"] == "inverse_K"))
            extra = sum(num(c["secondes"]) or 0 for c in x if c["etape"] in ("inverse_K", "valeurs_depart")) / n_runs
            out[r["cas"]] = dict(ll=num(r["logLik"]), it=r["n_iter"], s=num(r["secondes"]) + extra)
    return out


def remlax(cible, versions, etiq="chaud"):
    """Derniere version disponible par cas. Le temps par evaluation de v8 est
    nul (appels servis par le cache) : on le prend alors dans la version
    precedente, dont l'evaluation est la meme."""
    out, ev = {}, {}
    for v in versions:
        lignes = lire(os.path.join(RES, "banc_%s_%s.csv" % (cible, v)))
        # a chaud si disponible, sinon a froid (meme version)
        cas_chaud = {r["cas"] for r in lignes if r["outil"] == "remlax" and r["voie"] == "mme" and r["etiquette"] == etiq}
        for r in lignes:
            if r["outil"] == "remlax" and r["voie"] == "mme" and (
                    r["etiquette"] == etiq or (r["etiquette"] == "froid" and r["cas"] not in cas_chaud)):
                if r["cas"] not in out:
                    out[r["cas"]] = dict(ll=num(r["logLik"]), it=r["n_iter"], s=num(r["secondes"]), v=v)
                e = num(r["eval_s"])
                if e is not None and e > 1e-3 and r["cas"] not in ev:
                    ev[r["cas"]] = e
    for c in out:
        out[c]["ev"] = ev.get(c)
    return out


def asreml():
    out = {}
    for r in lire(ASR):
        if num(r.get("logLik")) is not None:
            out[r["cas"]] = dict(ll=num(r["logLik"]), it=r["n_iter"], s=num(r["secondes"]))
    return out


def f_ll(x):
    return "" if x is None else "%.2f" % x


def f_s(x):
    if x is None:
        return ""
    return "%.0f" % x if x >= 100 else "%.1f" % x


def classe(ll, ref):
    if ll is None or ref is None:
        return ""
    if ll > ref + 0.01:
        return ' class="mieux"'
    if ll < ref - 0.01:
        return ' class="moins"'
    return ""


K = {c: kmmr(c) for c in ("cpu", "gpu")}
R = {c: remlax(c, ["v10", "v9", "v8", "v7", "v6"]) for c in ("cpu", "gpu")}
A = asreml()

print('<div class="defile"><table class="banc"><thead><tr><th rowspan="2">cas</th>'
      '<th colspan="3">log-vraisemblance REML à l’optimum</th>'
      '<th colspan="3">CPU, 8 cœurs (s)</th><th colspan="2">A100 (s)</th></tr>'
      '<tr><th>kmmr</th><th>asreml-R</th><th>remlax</th>'
      '<th>kmmr</th><th>remlax</th><th>asreml</th><th>kmmr</th><th>remlax</th></tr></thead><tbody>')
for c in CAS:
    kc, kg, rc, rg, a = K["cpu"].get(c, {}), K["gpu"].get(c, {}), R["cpu"].get(c, {}), R["gpu"].get(c, {}), A.get(c, {})
    llr = rc.get("ll") if rc.get("ll") is not None else rg.get("ll")
    lls = [x for x in (kc.get("ll"), a.get("ll"), llr) if x is not None]
    best = max(lls) if lls else None
    print("<tr><th scope=\"row\">%s</th><td%s>%s</td><td%s>%s</td><td%s>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"
          % (NOM[c], classe(kc.get("ll"), best), f_ll(kc.get("ll")), classe(a.get("ll"), best), f_ll(a.get("ll")),
             classe(llr, best), f_ll(llr), f_s(kc.get("s")), f_s(rc.get("s")), f_s(a.get("s")),
             f_s(kg.get("s")), f_s(rg.get("s"))))
print("</tbody></table></div>")

# iterations et cout par evaluation
print('<div class="defile"><table class="banc"><thead><tr><th>cas</th><th>it. kmmr</th><th>it. asreml</th>'
      '<th>pas AI essayés remlax</th><th>remlax s/éval CPU</th><th>remlax s/éval A100</th><th>version</th></tr></thead><tbody>')
for c in CAS:
    kc, rc, rg, a = K["cpu"].get(c, {}), R["cpu"].get(c, {}), R["gpu"].get(c, {}), A.get(c, {})
    print("<tr><th scope=\"row\">%s</th><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"
          % (NOM[c], kc.get("it", ""), a.get("it", ""), rc.get("it", rg.get("it", "")),
             "" if rc.get("ev") is None else "%.3f" % rc["ev"], "" if rg.get("ev") is None else "%.3f" % rg["ev"],
             "/".join(sorted({x for x in (rc.get("v"), rg.get("v")) if x}))))
print("</tbody></table></div>")

# progression des versions
print('<div class="defile"><table class="banc"><thead><tr><th>cas, matériel</th>'
      '<th>v1</th><th>v2</th><th>v3</th><th>v6</th><th>v10</th></tr></thead><tbody>')
for c, cible in (("drops2_fa3", "cpu"), ("drops1_fa1", "cpu"), ("krause_gblup", "cpu"), ("krause_us", "cpu"),
                 ("krause_fa2", "gpu"), ("krause_us", "gpu")):
    cells = []
    for v in ("", "v2", "v3", "v6", "v10"):
        f = os.path.join(RES, "banc_%s%s.csv" % (cible, "_" + v if v else ""))
        rr = [r for r in lire(f) if r["cas"] == c and r["outil"] == "remlax" and r["voie"] == "mme"]
        rr = [r for r in rr if r["etiquette"] == "chaud"] or rr
        r = rr[-1] if rr else None
        cells.append("" if r is None or num(r["secondes"]) is None else
                     "%s s <span class=\"ll\">%s</span>" % (f_s(num(r["secondes"])), f_ll(num(r["logLik"]))))
    print("<tr><th scope=\"row\">%s, %s</th>%s</tr>" % (NOM[c], "A100" if cible == "gpu" else "CPU",
                                                       "".join("<td>%s</td>" % x for x in cells)))
print("</tbody></table></div>")

# phases a chaud v8
for cible in ("cpu", "gpu"):
    ch = [r for r in lire(os.path.join(RES, "banc_%s_v10_chrono.csv" % cible))
          if r["outil"] == "remlax" and r["voie"] == "mme" and r["etiquette"] == "chaud"]
    if not ch:
        continue
    et = [("demarrage_serveur", "démarrage"), ("preparation", "préparation"), ("compilation", "1re évaluation"),
          ("optimisation", "optimisation"), ("finition", "finition"), ("effets", "sorties"),
          ("solveur_total", "total solveur")]
    print('<h4>%s</h4><div class="defile"><table class="banc phases"><thead><tr><th>cas</th>%s</tr></thead><tbody>'
          % ("CPU, 8 cœurs" if cible == "cpu" else "A100", "".join("<th>%s</th>" % e[1] for e in et)))
    for c in CAS:
        x = [r for r in ch if r["cas"] == c]
        if not x:
            continue
        vals = [sum(num(r["secondes"]) or 0 for r in x if r["etape"] == e[0]) for e in et]
        print("<tr><th scope=\"row\">%s</th>%s</tr>" % (NOM[c], "".join("<td>%s</td>" % f_s(v) for v in vals)))
    print("</tbody></table></div>")
