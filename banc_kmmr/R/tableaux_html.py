"""HTML tables of the benchmark report, read from the genotoul and asreml CSV files.

    python3 tableaux_html.py <genotoul res folder> <banc_asreml.csv> > fragments.html

The files are bench_<target>[_<version>].csv with bench_..._timings.csv. The
stored files of 2026-10-27 use the old French names (banc_*.csv,
banc_*_chrono.csv) and French columns and values: they are read too, and
mapped to English right after reading.
"""
import csv
import os
import sys

RES, ASR = sys.argv[1], sys.argv[2]
CASES = ["drops2_fa1", "drops2_fa3", "drops1_fa1", "drops1_us", "krause_gblup", "krause_fa2", "krause_us"]
NAME = {"drops2_fa1": "DROPS two stages, fa1", "drops2_fa3": "DROPS two stages, fa3",
        "drops1_fa1": "DROPS one stage, fa1", "drops1_us": "DROPS one stage, us(5)",
        "krause_gblup": "Krause, GBLUP", "krause_fa2": "Krause, fa2 × 5 treatments",
        "krause_us": "Krause, us(5) × 5 treatments"}


# old French columns and values -> English
FR_COLS = {"etiquette": "label", "cas": "case", "outil": "tool", "voie": "route", "noeud": "node",
           "fils": "threads", "converge": "converged", "secondes": "seconds", "moteur": "engine",
           "t_modele": "t_model", "erreur": "error", "etape": "step", "delai": "timeout"}
FR_ROUTES = {"blocs": "blocks", "mme_creux": "mme_sparse", "mme_hybride": "mme_hybrid"}
FR_LABELS = {"froid": "cold", "chaud": "warm", "unique": "single"}
# timing steps: kmmr and asreml steps of this benchmark, and remlax steps
# (old French names of the CSV files already written -> current names)
FR_STEPS = {"valeurs_depart": "start_values", "premiere_iteration": "first_iteration",
            "iterations_suivantes": "later_iterations", "hors_solveur": "outside_solver",
            "ajustement": "fit", "hors_iterations": "outside_iterations",
            "demarrage_serveur": "server_startup", "solveur_total": "solver_total", "lecture_R": "read_R",
             "origine": "origin", "demarrage": "startup", "peripherique": "device", "lecture": "read",
             "finition": "finish", "effets": "effects", "hessien": "hessian", "ecriture": "write",
             "transmission": "transfer"}


def english_label(x):
    """cold / warm / single, with an optional r2_ prefix."""
    head, sep, last = x.rpartition("_")
    return head + sep + FR_LABELS.get(last, last)


def to_english(r):
    r = {FR_COLS.get(k, k): v for k, v in r.items()}
    if "label" in r:
        r["label"] = english_label(r["label"])
    if "route" in r:
        r["route"] = FR_ROUTES.get(r["route"], r["route"])
    if "step" in r:
        r["step"] = FR_STEPS.get(r["step"], r["step"])
    return r


def read(f):
    if not os.path.exists(f):
        return []
    with open(f) as fh:
        return [to_english(r) for r in csv.DictReader(fh)]


def results_file(target, v):
    """English name first, then the old French one. v is "" or "_v3"."""
    f = os.path.join(RES, "bench_%s%s.csv" % (target, v))
    return f if os.path.exists(f) else os.path.join(RES, "banc_%s%s.csv" % (target, v))


def timings_file(f):
    g = f.replace(".csv", "_timings.csv")
    return g if os.path.exists(g) else f.replace(".csv", "_chrono.csv")


def num(x):
    try:
        v = float(x)
        return None if v != v else v
    except (TypeError, ValueError):
        return None


def kmmr(target):
    """kmmr at tolerance 1e-6: time = call + inversion of K and starting
    values done in R before the call (first file where the case appears)."""
    out = {}
    for v in ("_v3", "_v2", ""):
        f = results_file(target, v)
        ch = read(timings_file(f))
        for r in read(f):
            if r["tool"] != "kmmr" or r["route"] != "ai_tol0.000001" or r["case"] in out or num(r["logLik"]) is None:
                continue
            x = [c for c in ch if c["case"] == r["case"] and c["tool"] == "kmmr" and c["route"] == "ai_tol0.000001"]
            n_runs = max(1, sum(1 for c in x if c["step"] == "inverse_K"))
            extra = sum(num(c["seconds"]) or 0 for c in x if c["step"] in ("inverse_K", "start_values")) / n_runs
            out[r["case"]] = dict(ll=num(r["logLik"]), it=r["n_iter"], s=num(r["seconds"]) + extra)
    return out


def remlax(target, versions, label="warm"):
    """Latest version available per case. The time per evaluation of v8 is
    zero (calls served by the cache), so it is then taken from the previous
    version, whose evaluation is the same."""
    out, ev = {}, {}
    for v in versions:
        rows = read(results_file(target, "_" + v))
        # warm if available, otherwise cold (same version)
        warm_cases = {r["case"] for r in rows if r["tool"] == "remlax" and r["route"] == "mme" and r["label"] == label}
        for r in rows:
            if r["tool"] == "remlax" and r["route"] == "mme" and (
                    r["label"] == label or (r["label"] == "cold" and r["case"] not in warm_cases)):
                if r["case"] not in out:
                    out[r["case"]] = dict(ll=num(r["logLik"]), it=r["n_iter"], s=num(r["seconds"]), v=v)
                e = num(r["eval_s"])
                if e is not None and e > 1e-3 and r["case"] not in ev:
                    ev[r["case"]] = e
    for c in out:
        out[c]["ev"] = ev.get(c)
    return out


def asreml():
    out = {}
    for r in read(ASR):
        if num(r.get("logLik")) is not None:
            out[r["case"]] = dict(ll=num(r["logLik"]), it=r["n_iter"], s=num(r["seconds"]))
    return out


def f_ll(x):
    return "" if x is None else "%.2f" % x


def f_s(x):
    if x is None:
        return ""
    return "%.0f" % x if x >= 100 else "%.1f" % x


# the CSS classes (mieux, moins, defile, banc) belong to the report page
# stylesheet, so they keep their names
def css_class(ll, ref):
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

print('<div class="defile"><table class="banc"><thead><tr><th rowspan="2">case</th>'
      '<th colspan="3">REML log-likelihood at the optimum</th>'
      '<th colspan="3">CPU, 8 cores (s)</th><th colspan="2">A100 (s)</th></tr>'
      '<tr><th>kmmr</th><th>asreml-R</th><th>remlax</th>'
      '<th>kmmr</th><th>remlax</th><th>asreml</th><th>kmmr</th><th>remlax</th></tr></thead><tbody>')
for c in CASES:
    kc, kg, rc, rg, a = K["cpu"].get(c, {}), K["gpu"].get(c, {}), R["cpu"].get(c, {}), R["gpu"].get(c, {}), A.get(c, {})
    llr = rc.get("ll") if rc.get("ll") is not None else rg.get("ll")
    lls = [x for x in (kc.get("ll"), a.get("ll"), llr) if x is not None]
    best = max(lls) if lls else None
    print("<tr><th scope=\"row\">%s</th><td%s>%s</td><td%s>%s</td><td%s>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"
          % (NAME[c], css_class(kc.get("ll"), best), f_ll(kc.get("ll")), css_class(a.get("ll"), best), f_ll(a.get("ll")),
             css_class(llr, best), f_ll(llr), f_s(kc.get("s")), f_s(rc.get("s")), f_s(a.get("s")),
             f_s(kg.get("s")), f_s(rg.get("s"))))
print("</tbody></table></div>")

# iterations and cost per evaluation
print('<div class="defile"><table class="banc"><thead><tr><th>case</th><th>it. kmmr</th><th>it. asreml</th>'
      '<th>AI steps tried, remlax</th><th>remlax s/eval CPU</th><th>remlax s/eval A100</th><th>version</th></tr></thead><tbody>')
for c in CASES:
    kc, rc, rg, a = K["cpu"].get(c, {}), R["cpu"].get(c, {}), R["gpu"].get(c, {}), A.get(c, {})
    print("<tr><th scope=\"row\">%s</th><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"
          % (NAME[c], kc.get("it", ""), a.get("it", ""), rc.get("it", rg.get("it", "")),
             "" if rc.get("ev") is None else "%.3f" % rc["ev"], "" if rg.get("ev") is None else "%.3f" % rg["ev"],
             "/".join(sorted({x for x in (rc.get("v"), rg.get("v")) if x}))))
print("</tbody></table></div>")

# progress across versions
print('<div class="defile"><table class="banc"><thead><tr><th>case, hardware</th>'
      '<th>v1</th><th>v2</th><th>v3</th><th>v6</th><th>v10</th></tr></thead><tbody>')
for c, target in (("drops2_fa3", "cpu"), ("drops1_fa1", "cpu"), ("krause_gblup", "cpu"), ("krause_us", "cpu"),
                 ("krause_fa2", "gpu"), ("krause_us", "gpu")):
    cells = []
    for v in ("", "v2", "v3", "v6", "v10"):
        f = results_file(target, "_" + v if v else "")
        rr = [r for r in read(f) if r["case"] == c and r["tool"] == "remlax" and r["route"] == "mme"]
        rr = [r for r in rr if r["label"] == "warm"] or rr
        r = rr[-1] if rr else None
        cells.append("" if r is None or num(r["seconds"]) is None else
                     "%s s <span class=\"ll\">%s</span>" % (f_s(num(r["seconds"])), f_ll(num(r["logLik"]))))
    print("<tr><th scope=\"row\">%s, %s</th>%s</tr>" % (NAME[c], "A100" if target == "gpu" else "CPU",
                                                       "".join("<td>%s</td>" % x for x in cells)))
print("</tbody></table></div>")

# warm phases, v10
for target in ("cpu", "gpu"):
    ch = [r for r in read(timings_file(results_file(target, "_v10")))
          if r["tool"] == "remlax" and r["route"] == "mme" and r["label"] == "warm"]
    if not ch:
        continue
    steps = [("server_startup", "startup"), ("preparation", "preparation"), ("compilation", "1st evaluation"),
          ("optimisation", "optimisation"), ("finish", "finish"), ("effects", "outputs"),
          ("solver_total", "solver total")]
    print('<h4>%s</h4><div class="defile"><table class="banc phases"><thead><tr><th>case</th>%s</tr></thead><tbody>'
          % ("CPU, 8 cores" if target == "cpu" else "A100", "".join("<th>%s</th>" % e[1] for e in steps)))
    for c in CASES:
        x = [r for r in ch if r["case"] == c]
        if not x:
            continue
        vals = [sum(num(r["seconds"]) or 0 for r in x if r["step"] == e[0]) for e in steps]
        print("<tr><th scope=\"row\">%s</th>%s</tr>" % (NAME[c], "".join("<td>%s</td>" % f_s(v) for v in vals)))
    print("</tbody></table></div>")
