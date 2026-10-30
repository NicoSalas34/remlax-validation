"""Synthesis of the absolute benchmark: tables and figures (time against n,
against p, against dim(C)/n), read from the CSV files of genotoul and meso.

    python3 synthese_abs.py <CSV folder> <output folder>

The CSV folder holds abs_cpu.csv, abs_gpu.csv (genotoul), abs_meso.csv (meso),
their *_timeouts.csv and *_timings.csv companions, and axes.csv (grille.py).
Results written before the English pass carry French names (*_delais.csv,
*_chrono.csv) and French columns. Both are read: the English file first, then
the French one, and French columns and values are mapped to English on reading.

Times kept (wall-clock seconds seen from R):
  kmmr    kmmr call + inversion of K and starting values computed in R
  asreml  asreml call + inversion of K computed in R
  remlax  full rx_fit, Python server start included. For the mme route, the
          "warm" pass (compilation cache filled) is kept; "cold" is kept apart.
"""
import collections
import csv
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.expanduser("~/.claude/skills/figure-style"))

RES, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

# --- French results written before the English pass ---------------------------
FR_COLUMNS = {"etiquette": "label", "cas": "case", "outil": "tool", "voie": "route", "noeud": "node",
              "fils": "threads", "converge": "converged", "secondes": "seconds", "moteur": "engine",
              "t_modele": "t_model", "erreur": "error", "etape": "step", "delai": "timeout", "axe": "axis"}
FR_LABELS = {"froid": "cold", "chaud": "warm", "unique": "single"}
FR_ROUTES = {"mme_creux": "mme_sparse", "mme_hybride": "mme_hybrid"}
FR_STEPS = {"finition": "finish", "premiere_iteration": "first_iteration", "iterations_suivantes": "later_iterations",
            "valeurs_depart": "start_values", "ajustement": "fit", "hors_iterations": "outside_iterations",
            "hors_solveur": "outside_solver", "demarrage_serveur": "server_startup", "solveur_total": "solver_total",
            "lecture_R": "read_R", "origine": "origin", "demarrage": "startup", "peripherique": "device",
            "lecture": "read", "effets": "effects", "hessien": "hessian", "ecriture": "write",
            "transmission": "transfer"}
FR_FILES = {"_timings.csv": "_chrono.csv", "_timeouts.csv": "_delais.csv"}


def to_english(row):
    """French column names and values of an old row -> English."""
    r = {FR_COLUMNS.get(k, k): v for k, v in row.items()}
    lab = r.get("label")
    if lab:
        for fr, en in FR_LABELS.items():
            if lab == fr or lab.endswith("_" + fr):
                r["label"] = lab[:-len(fr)] + en
                break
    if "route" in r:
        r["route"] = FR_ROUTES.get(r["route"], r["route"])
    if "step" in r:
        r["step"] = FR_STEPS.get(r["step"], r["step"])
    return r


def read_rows(f):
    """Rows of RES/f in English; falls back on the old French file name."""
    names = [f] + [f[:-len(en)] + fr for en, fr in FR_FILES.items() if f.endswith(en)]
    for name in names:
        path = os.path.join(RES, name)
        if os.path.exists(path):
            with open(path) as fh:
                return [to_english(r) for r in csv.DictReader(fh)]
    return []


def num(x):
    try:
        v = float(x)
        return None if v != v else v
    except (TypeError, ValueError):
        return None


def parse_case(c):
    d = dict((k[0], k[1:]) for k in c.split("_")[1:])
    q, t, r, x, o = int(d["q"]), int(d["t"]), int(d["r"]), float(d["x"]), float(d["o"])
    q_tot = int(round(q * (1 + x)))
    n = q * r * max(1, round(o * t / 100))
    return dict(q=q, t=t, r=r, s=d["s"], x=x, o=o, q_tot=q_tot, n=n, dimC=q_tot * t + t)


AXES = {}
for r in read_rows("axes.csv"):
    AXES.setdefault(r["axis"], []).append(r["case"])

# --- measurements: one entry per (hardware, route, case) ----------------------
# routes: asreml, kmmr, remlax (warm mme), remlax_cold, remlax_V (V form)
meas = {}
npar = {}
ll = {}
raw = collections.defaultdict(list)     # (hardware, route, case) -> one entry per run


def remlax_route(r):
    """Synthesis route of a remlax row. The label carries the repetition
    (r2_, r3_) and, for the mme route, the cold or warm pass."""
    e = r["label"]
    if r["route"] == "mme":
        return "remlax" if e.endswith("warm") else "remlax_cold"
    return {"dense": "remlax_Veco", "V_ai": "remlax_V", "V_lbfgsb": "V_lbfgsb", "mme_dense": "mme_dense",
            "mme_sparse": "mme_sparse", "mme_hybrid": "mme_hybrid"}.get(r["route"])


# asreml has a licence on meso only: it runs there on a RESERVED node (abs_meso_x),
# with kmmr and remlax as controls. The first meso passes, on overloaded shared
# nodes, are discarded (abs_meso, abs_meso2).
FILES = (("abs_cpu", "cpu"), ("abs_gpu", "a100"), ("abs_cpu_fixe", "cpu"), ("abs_gpu_fixe", "a100"),
         ("abs_cpu_rep", "cpu"), ("abs_gpu_rep", "a100"), ("abs_meso_x13", "meso"),
         ("abs_cpu_v13", "cpu"), ("abs_gpu_v13", "a100"), ("abs_cpu_v14", "cpu"), ("abs_gpu_v14", "a100"))
# remlax: one version per synthesis (REMLAX_VERSION=v13 reads the *_v13 files
# only; empty = the v10 passes). kmmr and asreml: every pass.
# v14 changes only the face of the us cone. The v14 pass redoes the us cases
# and fills the missing routes. The rest is read from the v13 files.
RX_VERSION = os.environ.get("REMLAX_VERSION", "v14")


def version_of(file):
    if file.endswith("_v14"):
        return "v14"
    return "v13" if file.endswith(("_v13", "_x13")) else ""


def canonical_set(c):
    """Three data sets per point: g1, g2, g3 if t = 1; g1, g4, g5 if t > 1 (the
    multi-trait g2 and g3 were copies of g1, seed bug of 2026-10-27)."""
    m = re.search(r"_(g\d+)$", c)
    g = m.group(1) if m else "g1"
    t = int(re.search(r"_t(\d+)_", c).group(1))
    return g in (("g1", "g2", "g3") if t == 1 else ("g1", "g4", "g5"))


def keep_rx(file, hw, r):
    """Does this remlax row belong to the version kept?"""
    v = version_of(file)
    if RX_VERSION != "v14":
        return v == RX_VERSION
    if not canonical_set(r["case"]):
        return False
    if v == "v14":
        return True
    # a data set redone in v14 for a route drops its v13 measurement
    return v == "v13" and (hw, r["route"], r["case"]) not in KEYS_V14


def base_case(c):
    """sim_..._g2 -> sim_...: the g1, g2, g3 sets of one point are grouped."""
    return re.sub(r"_g\d+$", "", c)


KEYS_V14 = set()
if RX_VERSION == "v14":
    for _f, _m in FILES:
        if version_of(_f) == "v14":
            for _r in read_rows(_f + ".csv"):
                if _r.get("tool") == "remlax":
                    KEYS_V14.add((_m, _r.get("route"), _r.get("case")))


for file, hw in FILES:
    timings = read_rows(file + "_timings.csv")
    t_opt, extra, t_iter = {}, collections.Counter(), {}
    for c in timings:
        key = (c["case"], c["tool"], c["route"], c["label"])
        if c["tool"] == "remlax" and c["step"] == "optimisation":
            t_opt[key] = num(c["seconds"])
        if (c["tool"] == "remlax" and c["step"] in ("optimisation", "finish")) or \
           (c["tool"] == "kmmr" and c["step"] in ("first_iteration", "later_iterations")):
            t_iter[key] = t_iter.get(key, 0) + (num(c["seconds"]) or 0)
        if c["tool"] in ("kmmr", "asreml") and c["step"] in ("inverse_K", "start_values"):
            extra[key] += num(c["seconds"]) or 0
    for r in read_rows(file + ".csv"):
        # row truncated by two simultaneous writes: skipped
        if not (r.get("case") or "").startswith("sim_") or r.get("tool") not in ("remlax", "kmmr", "asreml"):
            continue
        route = remlax_route(r) if r["tool"] == "remlax" else r["tool"]
        if route is None:
            continue
        if r["tool"] == "remlax" and not keep_rx(file, hw, r):
            continue
        key = (r["case"], r["tool"], r["route"], r["label"])
        ok = (r.get("error") or "") == "" and num(r.get("logLik")) is not None and r.get("case")
        s_ = num(r["seconds"])
        if ok and r["tool"] != "remlax":
            s_ += extra.get(key, 0)
        ev = num(r.get("eval_s"))
        raw[(hw, route, base_case(r["case"]))].append(dict(
            # likelihood compared between tools on the first data set only:
            # g2 and g3 are other data
            s=s_ if ok else None, ll=num(r["logLik"]) if ok and base_case(r["case"]) == r["case"] else None, it=r["n_iter"], err=(r.get("error") or "")[:120],
            engine=r.get("engine", ""), inv=r.get("inversion", ""), ev=ev if ok else None,
            nev=((t_opt[key] / ev) if ok and ev and t_opt.get(key) else None) if r["tool"] == "remlax"
                else (num(r["n_iter"]) if ok else None), node=r["node"],
            titer=t_iter.get(key) if ok else None, niter=num(r["n_iter"]) if ok else None))
        if r["tool"] == "remlax" and num(r["n_par"]) is not None:
            npar[base_case(r["case"])] = int(num(r["n_par"]))
    for r in read_rows(file + "_timeouts.csv"):
        route = remlax_route(r) if r["tool"] == "remlax" else r["tool"]
        if route is None:
            continue
        if r["tool"] == "remlax" and not keep_rx(file, hw, r):
            continue
        raw[(hw, route, base_case(r["case"]))].append(dict(s=None, ll=None, it="", err="timeout %s s" % r["timeout"],
                                                       timeout=float(r["timeout"])))
        if route == "remlax_cold":     # the cold pass timed out: the warm pass did not run
            raw[(hw, "remlax", base_case(r["case"]))].append(dict(s=None, ll=None, it="", err="timeout",
                                                              timeout=float(r["timeout"])))


def med(v):
    v = [x for x in v if x is not None]
    return float(np.median(v)) if v else None


for k, runs in raw.items():
    oks = [x for x in runs if x["s"] is not None]
    if oks:
        meas[k] = dict(s=med([x["s"] for x in oks]), smin=min(x["s"] for x in oks), smax=max(x["s"] for x in oks),
                       n_rep=len(oks), ll=med([x["ll"] for x in oks]), it=oks[0]["it"], err="",
                       engine=oks[0]["engine"], inv=oks[0]["inv"], ev=med([x["ev"] for x in oks]),
                       nev=med([x["nev"] for x in oks]), n_failures=len(runs) - len(oks),
                       titer=med([x.get("titer") for x in oks]), niter=med([x.get("niter") for x in oks]))
        if meas[k]["ll"] is not None:
            ll.setdefault(k[2], []).append((k[0], k[1], meas[k]["ll"]))
    else:
        d = next((x for x in runs if x.get("timeout")), None)
        meas[k] = dict(s=None, ll=None, it="", err=runs[0]["err"], n_rep=0, n_failures=len(runs),
                       **({"timeout": d["timeout"]} if d else {}))

CASES = []
for l in AXES.values():
    for c in l:
        if c not in CASES:
            CASES.append(c)

# --- agreement of the likelihoods ---------------------------------------------
# best likelihood of each case, and difference of each route from it
best = {c: max(v for _, _, v in ll[c]) for c in ll}


def diff_best(hw, route, c):
    m = meas.get((hw, route, c))
    if not m or m["ll"] is None or c not in best:
        return None
    return m["ll"] - best[c]


# --- tables -------------------------------------------------------------------
ROUTES = [("meso", "asreml", "asreml"), ("meso", "kmmr", "kmmr"), ("meso", "remlax", "remlax"),
          ("cpu", "kmmr", "kmmr"), ("cpu", "remlax", "remlax"), ("cpu", "remlax_V", "remlax V"),
          ("a100", "kmmr", "kmmr"), ("a100", "remlax", "remlax"), ("a100", "remlax_V", "remlax V")]
HEADERS = [("meso, 8 cores", 3), ("genotoul, 8 cores", 3), ("A100", 3)]


def f_s(x):
    if x is None:
        return ""
    return "%.0f" % x if x >= 100 else ("%.1f" % x if x >= 1 else "%.2f" % x)


def cell(hw, route, c):
    m = meas.get((hw, route, c))
    if m is None:
        return '<td class="vide"></td>'
    if m["s"] is None:
        lab = "timeout" if m.get("timeout") else "failed"
        return '<td class="echec" title="%s">%s</td>' % (m["err"].replace('"', "'"), lab)
    e = diff_best(hw, route, c)
    cls = ' class="moins"' if e is not None and e < -0.05 else ""
    note = '<span class="ll">%+.2f</span>' % e if e is not None and abs(e) > 0.05 else ""
    return "<td%s>%s%s</td>" % (cls, f_s(m["s"]), note)


TITLES = {"n_fixeC": "Time against n, dim(C) fixed (GBLUP, q = 2000)",
          "n_q": "Time against n, q = n/4 (GBLUP)",
          "p_t": "Time against the number of environments (us(t), q = 1000, one observation per pair)",
          "p_struct": "Time against p with the same data (t = 8, q = 1000)",
          "C_pred": "dim(C) >> n: genotypes to predict in the relationship matrix (q = 1000 phenotyped)",
          "C_clair": "dim(C) >> n: sparse trials (t = 10, q = 2000, fa2)"}


def describe(c):
    a = parse_case(c)
    return a, "n = %d, q<sub>K</sub> = %d, t = %d, %s, dim C = %d, p = %s" % (
        a["n"], a["q_tot"], a["t"], a["s"], a["dimC"], npar.get(c, "?"))


# the CSS classes (vide, echec, moins, defile, banc) belong to the report page
# that includes this fragment, and keep their names
with open(os.path.join(OUT, "tables.html"), "w") as f:
    for axis, cases in AXES.items():
        f.write('<h3>%s</h3>\n<div class="defile"><table class="banc"><thead><tr><th rowspan="2">case</th>' % TITLES[axis])
        f.write("".join('<th colspan="%d">%s</th>' % (k, h) for h, k in HEADERS) + "</tr><tr>")
        f.write("".join("<th>%s</th>" % v[2] for v in ROUTES) + "</tr></thead><tbody>\n")
        for c in cases:
            f.write('<tr><th scope="row">%s</th>%s</tr>\n' % (describe(c)[1], "".join(cell(m, v, c) for m, v, _ in ROUTES)))
        f.write("</tbody></table></div>\n")

# long CSV, one row per measurement
with open(os.path.join(OUT, "measurements.csv"), "w") as f:
    w = csv.writer(f)
    w.writerow(["case", "hardware", "route", "n", "q_K", "t", "struct", "dimC", "n_par", "seconds_median",
                "seconds_min", "seconds_max", "n_rep", "n_failures", "eval_s", "n_eval",
                "logLik", "diff_best", "n_iter", "engine", "inversion", "error"])
    for (hw, route, c), m in sorted(meas.items(), key=lambda kv: (kv[0][2], kv[0][0], kv[0][1])):
        a = parse_case(c)
        w.writerow([c, hw, route, a["n"], a["q_tot"], a["t"], a["s"], a["dimC"], npar.get(c, ""),
                    m["s"], m.get("smin"), m.get("smax"), m.get("n_rep"), m.get("n_failures"), m.get("ev"), m.get("nev"),
                    m["ll"], diff_best(hw, route, c), m["it"], m.get("engine", ""), m.get("inv", ""), m["err"]])

# --- figures ------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
from matplotlib.lines import Line2D
from figure_style import apply_figure_style

apply_figure_style(sizes=(8, 7, 6))


def p_theo(c):
    """variance parameters of the simulated model (genetic + residual)"""
    a = parse_case(c); t = a["t"]
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
    return ("%g s" % v)


def log_y_axis(ax):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(mt.LogLocator(base=10, subs=(1, 3) if ax.get_ylim()[1] / ax.get_ylim()[0] < 300 else (1,)))
    ax.yaxis.set_major_formatter(mt.FuncFormatter(fmt_s))
    ax.yaxis.set_minor_formatter(mt.NullFormatter())


def x_ticks(ax, vals, labs, log=True):
    if log:
        ax.set_xscale("log")
    ax.set_xticks(vals); ax.set_xticklabels(labs)
    ax.xaxis.set_minor_locator(mt.NullLocator())


def k(x):
    return "%dk" % (x // 1000) if x >= 1000 else str(x)


# --- comparison figures: one figure per axis, columns CPU | A100 ---------------
# Colour = software or engine. The hardware is given by the column, never by
# the line. Names written at the end of the curves.
COLOURS = {"remlax": "#1d6b58", "remlax_V": "#4f9fd9", "kmmr": "#d08a2e", "asreml": "#6f6f6f"}
NAMES = {"remlax": "remlax MME", "remlax_V": "remlax V form", "kmmr": "kmmr", "asreml": "ASReml-R (meso)"}
HW_ROUTES = {"cpu": ["remlax", "remlax_V", "kmmr", ("meso", "asreml")], "a100": ["remlax", "remlax_V", "kmmr"]}
COL_TITLE = {"cpu": "CPU, 8 cores", "a100": "GPU A100"}


def series(cases, xf, hw, route):
    pts, fails, timeouts = [], [], []
    for c in cases:
        m = meas.get((hw, route, c))
        if not m:
            continue
        if m["s"] is not None:
            pts.append((xf(c), m["s"], m.get("smin"), m.get("smax"), m.get("n_rep", 1)))
        elif m.get("timeout"):
            timeouts.append(xf(c))
        else:
            fails.append(xf(c))
    return sorted(pts), fails, timeouts


def end_labels(ax, ends):
    """Names at the end of the curves, placed once the scales are frozen, and
    shifted vertically (on the log scale) so that they do not overlap."""
    if not ends:
        return
    import matplotlib.transforms as mtr
    y0, y1 = [np.log10(v) for v in ax.get_ylim()]
    gap = 0.07 * (y1 - y0)
    ends = sorted(ends, key=lambda b: b[1])
    ly = [np.log10(b[1]) for b in ends]
    for i in range(1, len(ly)):
        ly[i] = max(ly[i], ly[i - 1] + gap)
    overflow = ly[-1] - (y1 - 0.02 * (y1 - y0))
    if overflow > 0:                       # move everything down if the last one leaves the frame
        ly = [l - overflow for l in ly]
    tr = mtr.offset_copy(ax.transData, fig=ax.figure, x=5, y=0, units="points")
    for (x, y, txt, col), l in zip(ends, ly):
        ax.text(x, 10 ** l, txt, color=col, fontsize=7, va="center", ha="left", transform=tr, clip_on=False)


def panel(ax, cases, xf, hw, xlog, ticks):
    ends = []
    for v in HW_ROUTES[hw]:
        m_, route = (v if isinstance(v, tuple) else (hw, v))
        pts, fails, timeouts = series(cases, xf, m_, route)
        col = COLOURS[route]; lw = 2.0 if route == "remlax" else 1.3
        if pts:
            x = [p[0] for p in pts]; y = [p[1] for p in pts]
            ax.plot(x, y, "-o", color=col, lw=lw, ms=3.8, zorder=3 if route == "remlax" else 2)
            for p in pts:
                if p[4] > 1 and p[3] > p[2] * 1.02:
                    ax.plot([p[0], p[0]], [p[2], p[3]], color=col, lw=1.0, alpha=0.55, zorder=1)
            ends.append((x[-1], y[-1], NAMES[route], col))
        for xx, mk in ((fails, "x"), (timeouts, "^")):
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
        ax.set_xlim(min(vals) / 1.25, max(vals) * 2.6)      # room for the names on the right
    else:
        ax.set_xlim(min(vals) - 0.06 * (max(vals) - min(vals)), max(vals) + 0.45 * (max(vals) - min(vals)))
    return ends


def axis_figure(name, rows, title):
    """rows: [(axis, xf, xlabel, row title, xlog, (vals, labs))]. One block
    (subfigure) per row, with its title on top; columns CPU | A100."""
    fig = plt.figure(figsize=(7.2, 3.15 * len(rows) + 0.5), layout="constrained")
    blocks = fig.subfigures(len(rows), 1, hspace=0.04)
    blocks = np.atleast_1d(blocks)
    every = []
    for i, ((axis, xf, xl, tl, xlog, ticks), block) in enumerate(zip(rows, blocks)):
        block.suptitle(tl, x=0.01, ha="left", fontsize=8)
        axs = block.subplots(1, 2, sharey=True)
        cases = AXES.get(axis, [])
        for j, hw in enumerate(("cpu", "a100")):
            ax = axs[j]
            every.append((ax, panel(ax, cases, xf, hw, xlog, ticks)))
            ax.set_xlabel(xl)
            ax.set_title("%s   %s" % ("abcdefgh"[2 * i + j], COL_TITLE[hw]), loc="left", fontsize=8)
        axs[0].set_ylabel("time of the fit")
    fig.suptitle(title, x=0.01, ha="left", fontsize=8.5, fontweight="bold")
    fig.get_layout_engine().set(w_pad=0.05, h_pad=0.06, rect=(0, 0.025, 1, 0.97))
    fig.canvas.draw()                    # scales and layout frozen
    for ax, ends in every:
        end_labels(ax, ends)
    fig.text(0.01, 0.002, "× failed fit   △ 6 h limit exceeded   │ range of the replicates "
             "(median plotted). kmmr and remlax on genotoul; ASReml-R on a reserved meso node.",
             fontsize=6.5, color="0.35", ha="left", va="bottom")
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=200)
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    plt.close(fig)


A = AXES
fn = lambda c: parse_case(c)["n"]
fr = lambda c: parse_case(c)["dimC"] / parse_case(c)["n"]
tk_n = lambda axis: (sorted({fn(c) for c in A.get(axis, [])}), [k(v) for v in sorted({fn(c) for c in A.get(axis, [])})])
tk_r = lambda axis: (sorted({fr(c) for c in A.get(axis, [])}), ["%g" % round(v, 1) for v in sorted({fr(c) for c in A.get(axis, [])})])
ps = sorted(A.get("p_struct", []), key=p_theo)
for name in ("fig_comparatif",):          # figure of an earlier version, removed from the output folder
    f = os.path.join(OUT, name + ".png")
    if os.path.exists(f):
        os.remove(f)
axis_figure("fig_time_n", [
    ("n_fixeC", fn, "observations n", "dim(C) fixed: one trait, q = 2000 genotypes, 1 to 64 replicates", True, tk_n("n_fixeC")),
    ("n_q", fn, "observations n", "dim(C) grows with n: one trait, 4 replicates, q = n/4", True, tk_n("n_q"))],
    "Time of the fit against the number of observations (GBLUP, dense relationship matrix)")
axis_figure("fig_time_p", [
    ("p_t", p_theo, "variance parameters p", "us(t) for t = 1 to 12 environments; n and dim(C) grow with t", True,
     (sorted({p_theo(c) for c in A.get("p_t", [])}), [str(v) for v in sorted({p_theo(c) for c in A.get("p_t", [])})])),
    ("p_struct", p_theo, "genetic structure (parameters p)", "t = 8 environments, same data; only the structure changes", False,
     ([p_theo(c) for c in ps], ["%s\n%d" % (parse_case(c)["s"], p_theo(c)) for c in ps]))],
    "Time of the fit against the number of variance parameters")
axis_figure("fig_time_C", [
    ("C_pred", fr, "dim(C) / n", "genotypes to predict: 1000 phenotyped, 1000 to 16 000 in the relationship matrix", True, tk_r("C_pred")),
    ("C_clair", fr, "dim(C) / n", "sparse trials: t = 10, q = 2000, fa2, 100 % to 10 % of the pairs observed", True, tk_r("C_clair"))],
    "Time of the fit when C grows large against n")

# --- figure of the rule: fully fixed remlax routes -----------------------------
CV = {"remlax_V": "#4f9fd9", "mme_dense": "#1d6b58", "mme_hybrid": "#8a6bb8", "mme_sparse": "#c0504d"}
NV = {"remlax_V": "V form (size n)", "mme_dense": "MME, dense inversion", "mme_hybrid": "MME, hybrid",
      "mme_sparse": "MME, sparse"}
fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.4), sharey="row")
for j, hw in enumerate(("cpu", "a100")):
    ax = axs[0, j]
    for route in (("remlax_V", "mme_dense") if hw == "a100" else ("remlax_V", "mme_dense", "mme_hybrid", "mme_sparse")):
        pts = [((parse_case(c)["n"] if route == "remlax_V" else parse_case(c)["dimC"]), meas[(hw, route, c)]["ev"]) for c in CASES
               if meas.get((hw, route, c), {}).get("ev")]
        if pts:
            x, y = zip(*pts)
            ax.plot(x, y, "o", ms=3.8, color=CV[route], alpha=0.85, label=NV[route])
    # slope 3 (dense factorisation) as a guide
    xs = np.array([1e3, 3e4]); ax.plot(xs, 0.02 * (xs / 1e3) ** 3 * (0.05 if hw == "a100" else 1), color="0.75",
                                       lw=0.8, ls="--", zorder=0)
    ax.text(3e4, 0.02 * 30 ** 3 * (0.05 if hw == "a100" else 1), " slope 3", color="0.55", fontsize=6.5, va="center")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.yaxis.set_major_formatter(mt.FuncFormatter(fmt_s)); ax.yaxis.set_minor_formatter(mt.NullFormatter())
    ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: k(int(round(v))))); ax.xaxis.set_minor_formatter(mt.NullFormatter())
    ax.set_xlabel("size of the factorised matrix (n or dim C)")
    ax.set_title(COL_TITLE[hw], fontweight="bold")
    ax.grid(color="0.92", lw=0.6, zorder=0)
    if j == 0:
        ax.set_ylabel("time of one evaluation")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    ax = axs[1, j]
    pts = []
    for c in CASES:
        v, m_ = meas.get((hw, "remlax_V", c)), meas.get((hw, "mme_dense", c))
        if v and m_ and v.get("ev") and m_.get("ev"):
            pts.append((parse_case(c)["dimC"] / parse_case(c)["n"], m_["ev"] / v["ev"]))
    if pts:
        x, y = zip(*pts)
        ax.plot(x, y, "o", ms=4, color="0.25")
    ax.set_xscale("log"); ax.set_yscale("log")
    lo, hi = 1e-4, 300
    ax.set_ylim(lo, hi)
    ax.axhspan(1, hi, color=CV["remlax_V"], alpha=0.10, zorder=0)
    ax.axhspan(lo, 1, color=CV["mme_dense"], alpha=0.08, zorder=0)
    ax.axhline(1, color="0.5", lw=0.8)
    ax.text(0.04, 0.95, "V form cheaper", transform=ax.transAxes, color=CV["remlax_V"], va="top", fontsize=7)
    ax.text(0.04, 0.05, "MME cheaper", transform=ax.transAxes, color=CV["mme_dense"], va="bottom", fontsize=7)
    ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: ("%g" % v)))
    ax.yaxis.set_minor_formatter(mt.NullFormatter())
    ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: ("%g" % v)))
    ax.xaxis.set_minor_formatter(mt.NullFormatter())
    ax.set_xlabel("dim(C) / n")
    if j == 0:
        ax.set_ylabel("cost of one evaluation\ndense MME / V form")
for i, ax in enumerate(axs.flat):
    ax.text(-0.13 if i % 2 == 0 else -0.05, 1.04, "abcd"[i], transform=ax.transAxes, fontweight="bold", fontsize=9)
fig.text(0.01, 0.005, "Fully fixed remlax routes, median of the replicates. One point = one simulated case.",
         fontsize=6.5, color="0.35")
fig.tight_layout(rect=(0, 0.03, 1, 1), h_pad=2.0)
fig.savefig(os.path.join(OUT, "fig_rule.png"), dpi=200)
fig.savefig(os.path.join(OUT, "fig_rule.pdf"))
plt.close(fig)

# --- decomposed figure: time = fixed + n_eval x cost of one evaluation ---------
# remlax: n_eval = optimisation time / cost of one full evaluation.
# kmmr: one evaluation = one AI iteration. The fixed cost is the rest of the
# total time (start-up, preparation, compilation, inversion of K, outputs).
DEC_ROWS = [
    ("n_fixeC", fn, "n", "n, dim(C) fixed (q = 2000)", True),
    ("n_q", fn, "n", "n, q = n/4", True),
    ("p_t", p_theo, "p", "us(t), t = 1 to 12", True),
    ("p_struct", p_theo, "structure", "t = 8, structures", False),
    ("C_pred", fr, "dim(C)/n", "genotypes to predict", True),
    ("C_clair", fr, "dim(C)/n", "sparse trials", True),
]
PARTS = [("cit", "cost of one iteration", True), ("niter", "number of iterations", False),
         ("fixed", "fixed cost (outside the iterations)", True)]


def value(m, key):
    """Measured quantities: number of iterations, time spent in the iterations
    (remlax: optimisation and finish; kmmr: sum of its iterations), and the
    rest of the total time."""
    if not m or m["s"] is None or not m.get("titer") or not m.get("niter"):
        return None
    if key == "niter":
        return m["niter"]
    if key == "cit":
        return m["titer"] / m["niter"]
    return max(m["s"] - m["titer"], 0.05)


def readable_log_axis(ax):
    """1-2-5 ticks on a log scale, readable even over less than a decade."""
    ax.set_yscale("log")
    lo, hi = ax.get_ylim()
    subs = (1.0, 2.0, 5.0) if hi / lo < 100 else (1.0,)
    ax.yaxis.set_major_locator(mt.LogLocator(base=10, subs=subs))
    ax.yaxis.set_major_formatter(mt.FuncFormatter(fmt_s)); ax.yaxis.set_minor_formatter(mt.NullFormatter())


def decomposition_figure(hw, name):
    fig = plt.figure(figsize=(7.2, 1.95 * len(DEC_ROWS) + 0.7), layout="constrained")
    axs = fig.subplots(len(DEC_ROWS), 3)
    for i, (axis, xf, xl, tl, xlog) in enumerate(DEC_ROWS):
        cases = AXES.get(axis, [])
        if axis == "p_struct":
            cases = sorted(cases, key=p_theo)
            vals = [p_theo(c) for c in cases]; labs = [parse_case(c)["s"] for c in cases]
        else:
            vals = sorted({xf(c) for c in cases})
            labs = [k(v) if axis.startswith("n_") else ("%g" % round(v, 1)) for v in vals]
        for j, (key, title, ylog) in enumerate(PARTS):
            ax = axs[i, j]
            ends = []
            for route in ("remlax", "kmmr"):
                pts, fails = [], []
                for c in cases:
                    m = meas.get((hw, route, c))
                    v = value(m, key)
                    if v is not None:
                        pts.append((xf(c), v))
                    elif m is not None and m["s"] is None:
                        fails.append(xf(c))
                pts.sort()
                if pts:
                    x, y = zip(*pts)
                    ax.plot(x, y, "-o", color=COLOURS[route], lw=2.0 if route == "remlax" else 1.3, ms=3.2)
                    ends.append((x[-1], y[-1], NAMES[route], COLOURS[route]))
                if fails:
                    ax.plot(fails, [0.93] * len(fails), "x", ms=5, mew=1.3, color=COLOURS[route],
                            transform=ax.get_xaxis_transform(), clip_on=False)
            if xlog:
                ax.set_xscale("log")
            ax.set_xticks(vals); ax.set_xticklabels(labs, fontsize=6); ax.xaxis.set_minor_locator(mt.NullLocator())
            if xlog:
                ax.set_xlim(min(vals) / 1.3, max(vals) * 1.3)
            if ylog:
                ax.set_yscale("log"); ax.autoscale_view(); readable_log_axis(ax)
            else:
                ax.set_yscale("log"); ax.autoscale_view(); readable_log_axis(ax)
                ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: "%g" % v))
            ax.grid(axis="y", color="0.92", lw=0.6)
            ax.tick_params(labelsize=6)
            if i == 0:
                ax.set_title(title, fontsize=8, fontweight="bold")
            if j == 0:
                ax.set_ylabel(tl, fontsize=7)
            if i == len(DEC_ROWS) - 1 or True:
                ax.set_xlabel(xl, fontsize=6.5, labelpad=1)
    fig.suptitle("Decomposition of the time of the fit, %s: time = fixed cost + number of iterations × cost of one iteration"
                 % COL_TITLE[hw], x=0.01, ha="left", fontsize=8)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color=COLOURS[v], lw=2.0 if v == "remlax" else 1.3, marker="o", ms=3.2, label=NAMES[v])
                        for v in ("remlax", "kmmr")] + [Line2D([], [], color="0.3", ls="none", marker="x", label="failed")],
               loc="outside lower center", ncol=3, frameon=False, fontsize=7)
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=200)
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    plt.close(fig)


decomposition_figure("cpu", "fig_decomposition_cpu")
decomposition_figure("a100", "fig_decomposition_a100")

# --- CPU against GPU -------------------------------------------------------------
# GPU speed-up = CPU time on 8 cores / A100 time, same route, same case, same
# cluster (genotoul). A point exists only if both devices succeeded.
def gpu_speedup(cases, xf, route):
    pts = []
    for c in cases:
        mc, mg = meas.get(("cpu", route, c)), meas.get(("a100", route, c))
        if mc and mg and mc["s"] and mg["s"]:
            pts.append((xf(c), mc["s"] / mg["s"]))
    return sorted(pts)


def fmt_speedup(v, _):
    return ("×%g" % v)


tk_p = (sorted({p_theo(c) for c in A.get("p_t", [])}), [str(v) for v in sorted({p_theo(c) for c in A.get("p_t", [])})])
CG_ROWS = [("n_fixeC", fn, "observations n", "dim(C) fixed, q = 2000", True, tk_n("n_fixeC")),
           ("n_q", fn, "observations n", "dim(C) grows with n, q = n/4", True, tk_n("n_q")),
           ("p_t", p_theo, "variance parameters p", "us(t), t = 1 to 12", True, tk_p),
           ("p_struct", p_theo, "genetic structure (parameters p)", "t = 8, only the structure changes", False,
            ([p_theo(c) for c in ps], ["%s\n%d" % (parse_case(c)["s"], p_theo(c)) for c in ps])),
           ("C_pred", fr, "dim(C) / n", "genotypes to predict", True, tk_r("C_pred")),
           ("C_clair", fr, "dim(C) / n", "sparse trials, fa2", True, tk_r("C_clair"))]
fig, axs = plt.subplots(3, 2, figsize=(7.2, 8.6), sharey=True, layout="constrained")
every = []
for ax, (axis, xf, xl, tl, xlog, (vals, labs)), letter in zip(axs.flat, CG_ROWS, "abcdef"):
    ends = []
    for route in ("remlax", "remlax_V", "kmmr"):
        pts = gpu_speedup(A.get(axis, []), xf, route)
        if pts:
            x, y = zip(*pts)
            ax.plot(x, y, "-o", color=COLOURS[route], lw=2.0 if route == "remlax" else 1.3, ms=3.8,
                    zorder=3 if route == "remlax" else 2)
            ends.append((x[-1], y[-1], NAMES[route], COLOURS[route]))
    ax.axhline(1, color="0.45", lw=0.7, ls="-", zorder=1)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(mt.FuncFormatter(fmt_speedup)); ax.yaxis.set_minor_formatter(mt.NullFormatter())
    if xlog:
        ax.set_xscale("log")
    ax.set_xticks(vals); ax.set_xticklabels(labs); ax.xaxis.set_minor_locator(mt.NullLocator())
    if xlog:
        ax.set_xlim(min(vals) / 1.25, max(vals) * 2.6)
    else:
        ax.set_xlim(min(vals) - 0.06 * (max(vals) - min(vals)), max(vals) + 0.45 * (max(vals) - min(vals)))
    ax.grid(axis="y", which="major", color="0.9", lw=0.6, zorder=0)
    ax.set_xlabel(xl)
    ax.set_title("%s   %s" % (letter, tl), loc="left", fontsize=8)
    every.append((ax, ends))
for ax in axs[:, 0]:
    ax.set_ylabel("GPU speed-up (CPU time / A100 time)")
fig.suptitle("Speed-up of the A100 over 8 CPU cores, route by route", x=0.01, ha="left", fontsize=8.5, fontweight="bold")
fig.get_layout_engine().set(w_pad=0.05, h_pad=0.06, rect=(0, 0.03, 1, 0.955))
fig.canvas.draw()
for ax, ends in every:
    end_labels(ax, ends)
fig.text(0.01, 0.002, "Above the grey line the A100 is faster. A point is missing when either device "
         "failed, was stopped or has not run yet.", fontsize=6.5, color="0.35", ha="left", va="bottom")
fig.savefig(os.path.join(OUT, "fig_cpu_gpu.png"), dpi=200)
fig.savefig(os.path.join(OUT, "fig_cpu_gpu.pdf"))
plt.close(fig)

# scatter: CPU time against A100 time, one point per case and per route
fig, ax = plt.subplots(figsize=(4.8, 4.6), layout="constrained")
all_t = []
for route in ("kmmr", "remlax_V", "remlax"):
    xy = [(meas[("cpu", route, c)]["s"], meas[("a100", route, c)]["s"]) for c in CASES
          if meas.get(("cpu", route, c)) and meas.get(("a100", route, c))
          and meas[("cpu", route, c)]["s"] and meas[("a100", route, c)]["s"]]
    if xy:
        x, y = zip(*xy); all_t += list(x) + list(y)
        ax.scatter(x, y, s=16, color=COLOURS[route], alpha=0.85, lw=0, label=NAMES[route], zorder=3)
lo, hi = min(all_t) / 1.5, max(all_t) * 1.5
for g, ls in ((1, "-"), (10, "--"), (100, ":")):
    ax.plot([lo, hi], [lo / g, hi / g], color="0.5", lw=0.8, ls=ls, zorder=1)
    ax.text(hi, hi / g, "  ×%d" % g if g > 1 else "  equality", fontsize=6.5, color="0.4", va="center", ha="left",
            clip_on=False)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
for a_ in (ax.xaxis, ax.yaxis):
    a_.set_major_formatter(mt.FuncFormatter(fmt_s)); a_.set_minor_formatter(mt.NullFormatter())
ax.set_xlabel("CPU time, 8 cores"); ax.set_ylabel("A100 GPU time")
ax.grid(color="0.92", lw=0.6, zorder=0)
ax.legend(frameon=False, fontsize=7, loc="upper left")
ax.set_title("CPU time against A100 time, every case of the benchmark", loc="left", fontsize=8.5, fontweight="bold")
fig.savefig(os.path.join(OUT, "fig_cpu_gpu_scatter.png"), dpi=200)
fig.savefig(os.path.join(OUT, "fig_cpu_gpu_scatter.pdf"))
plt.close(fig)

# --- remlax against asreml, same problem -------------------------------------------
# asreml: reserved meso node. remlax MME v14: 8 CPU cores of genotoul, then A100.
FAM = {"n_fixeC": "GBLUP, one trait", "n_q": "GBLUP, one trait", "p_t": "multi-trait (us, fa, diag)",
       "p_struct": "multi-trait (us, fa, diag)", "C_pred": "dim(C) large against n", "C_clair": "dim(C) large against n"}
CFAM = {"GBLUP, one trait": "#1d6b58", "multi-trait (us, fa, diag)": "#8a6bb8", "dim(C) large against n": "#d08a2e"}
MFAM = {"GBLUP, one trait": "o", "multi-trait (us, fa, diag)": "s", "dim(C) large against n": "^"}
family_of = {}
for axis, l in AXES.items():
    for c in l:
        family_of.setdefault(c, FAM.get(axis))
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.9), sharex=True, sharey=True, layout="constrained")
pairs = []
for c in CASES:
    ma = meas.get(("meso", "asreml", c))
    if ma and ma["s"]:
        pairs.append((c, ma["s"]))
tt = [s_ for _, s_ in pairs]
for c, _ in pairs:
    for m_ in ("cpu", "a100"):
        mr = meas.get((m_, "remlax", c))
        if mr and mr["s"]:
            tt.append(mr["s"])
lo, hi = min(tt) / 1.6, max(tt) * 1.6
for ax, m_, letter, tl in ((axs[0], "cpu", "a", "remlax MME on CPU, 8 cores"), (axs[1], "a100", "b", "remlax MME on A100 GPU")):
    for g, ls in ((1, "-"), (10, "--"), (100, ":")):
        ax.plot([lo, hi], [lo / g, hi / g], color="0.55", lw=0.8, ls=ls, zorder=1)
        ax.text(hi, hi / g, " equality" if g == 1 else " ×%d" % g, fontsize=6.5, color="0.4", va="center",
                ha="left", clip_on=False)
    seen = set()
    for c, ta in pairs:
        mr = meas.get((m_, "remlax", c))
        if not (mr and mr["s"]):
            continue
        fa = family_of.get(c)
        ax.scatter(ta, mr["s"], s=24, marker=MFAM[fa], color=CFAM[fa], lw=0, alpha=0.9, zorder=3,
                   label=fa if fa not in seen else None)
        seen.add(fa)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    for a_ in (ax.xaxis, ax.yaxis):
        a_.set_major_formatter(mt.FuncFormatter(fmt_s)); a_.set_minor_formatter(mt.NullFormatter())
    ax.set_aspect("equal")
    ax.grid(color="0.92", lw=0.6, zorder=0)
    ax.set_xlabel("ASReml-R time (CPU, 8 cores)")
    ax.set_title("%s   %s" % (letter, tl), loc="left", fontsize=8)
axs[0].set_ylabel("remlax MME time")
axs[0].legend(frameon=False, fontsize=6.5, loc="upper left")
fig.suptitle("remlax against ASReml-R on the same problems", x=0.01, ha="left", fontsize=8.5, fontweight="bold")
fig.get_layout_engine().set(w_pad=0.08, rect=(0, 0.06, 1, 1))
fig.text(0.01, 0.005, "Below the diagonal remlax is faster. ASReml-R ran on a reserved meso node, remlax v14 on genotoul; "
         "the meso node runs remlax\n1.3 to 3 times faster than genotoul, so the panel is to the "
         "disadvantage of remlax. %d problems, those ASReml-R completed." % len(pairs),
         fontsize=6.5, color="0.35", ha="left", va="bottom")
fig.savefig(os.path.join(OUT, "fig_remlax_asreml.png"), dpi=200)
fig.savefig(os.path.join(OUT, "fig_remlax_asreml.pdf"))
plt.close(fig)

n_ok = sum(1 for m in meas.values() if m["s"] is not None)
print("successful measurements:", n_ok, "of", len(meas))
