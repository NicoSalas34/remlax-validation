# =============================================================================
# summarise_benchc.py: tables and figures of the "complexity" benchmark
#                      (bench_complexite.R)
# -----------------------------------------------------------------------------
# Reads benchmarks/results/benchc_<date>_*.csv (asreml, remlax cpu, remlax gpu)
# and writes, in --out (default: the --results folder):
#   benchc_<date>_table.csv   one row per (case, n, software)
#   benchc_<date>_wide.md     wide table per (case, n): asreml wall /
#        iterations / s per iteration; remlax CPU and GPU: startup,
#        compilation, solver, iterations, evaluations, s per iteration,
#        s per evaluation, wall; logLik difference to asreml
#   fig_benchc_wall.png       total time against n, per case
#   fig_benchc_iter.png       time per iteration against n_par
#   fig_benchc_compile.png    XLA compilation against n and case
#   fig_benchc_breakdown.png  startup / compilation / solver
#   fig_benchc_density.png    wall time against q at fixed n
#
#   python benchmarks/summarise_benchc.py [--date 2026-09-28] [--results benchmarks/results] [--out <dir>]
# In the session Python kernel, apply_figure_style() (figure-style skill) is
# applied when available.
# =============================================================================
import argparse, glob, os, sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import textwrap

ap = argparse.ArgumentParser()
ap.add_argument("--date", default=None)
ap.add_argument("--results", default=os.path.join(os.path.dirname(__file__) or ".", "results"))
ap.add_argument("--out", default=None, help="output folder (default: the --results folder)")
A = ap.parse_args([] if "ipykernel" in sys.modules or not sys.argv[1:] else None)
RES = A.results
OUT = A.out or RES
os.makedirs(OUT, exist_ok=True)
files = sorted(glob.glob(os.path.join(RES, "benchc_*.csv")))
files = [f for f in files if not f.endswith("_table.csv")]
if A.date:
    files = [f for f in files if A.date in os.path.basename(f)]
assert files, "no benchc_*.csv in %s" % RES
df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
if "reps" not in df.columns:
    df["reps"] = np.nan
if "q" not in df.columns:
    df["q"] = np.nan
# The CSVs of the size axis (first version of the script) have no reps column:
# 4 replicates per genotype for iid, grm, ar1ar1, ige; 2 for the multi-trait
# cases (q = n1 %/% 2). The density axis carries them explicitly.
missing = df["reps"].isna()
df.loc[missing & df["case"].astype(str).str.startswith("us"), "reps"] = 2
df["reps"] = df["reps"].fillna(4).astype(int)
# q is absent from the CSVs of the size axis: q = n / reps, or (n / t) / reps for the multi-trait cases
def _q(r):
    if np.isfinite(r["q"]) if isinstance(r["q"], float) else r["q"] is not None:
        return r["q"]
    # first version of the script: q = n %/% 4, or (n %/% t) %/% 2 for the multi-trait cases
    case = str(r["case"]); n = int(r["n"])
    if case.startswith("us"):
        t = int(case.replace("usK", "").replace("us", ""))
        return (n // t) // 2
    if case == "ige":
        return round(n / 4)
    return n // 4
df["q"] = df.apply(_q, axis=1)
DATE = A.date or df["date"].max()

# ---- software x backend -> series ----------------------------------------------
def series_of(r):
    if r["software"] == "asreml":
        return "asreml"
    return "remlax GPU (A100)" if r["backend"] == "gpu" else "remlax CPU (4 cores)"
df["series"] = df.apply(series_of, axis=1)
# allocated CPU cores, read from the run tag (cpu4 / cpu8); 4 by default
df["cores"] = df["tag"].astype(str).str.extract(r"cpu(\d+)")[0].fillna("4").astype(int)
df.loc[df["backend"] == "gpu", "cores"] = 8
CASE_ORDER = ["iid", "grm", "ar1ar1", "us3", "us6", "us9", "us12", "usK3", "usK6", "ige"]
CASE_LABEL = {"iid": "one random factor", "grm": "GRM (dense kinship)", "ar1ar1": "AR1 x AR1 field",
              "us3": "3 traits, us + us", "us6": "6 traits, us + us", "us9": "9 traits, us + us",
              "us12": "12 traits, us + us", "usK3": "3 traits x GRM", "usK6": "6 traits x GRM",
              "ige": "Salas et al. (2026)"}
df["case"] = pd.Categorical(df["case"], [c for c in CASE_ORDER if c in set(df["case"])], ordered=True)
df = df.sort_values(["case", "n", "series"]).reset_index(drop=True)
# nominal size: the multi-trait cases lose a few observations to rounding
df["n_nom"] = df["n"].apply(lambda n: int(min([500, 2000, 8000, 16000, 32000], key=lambda s: abs(s - n))))
# n_par: the same for the three programs; the remlax value is taken when asreml lacks it
npar = df[df["status"] == "ok"].groupby(["case", "n_nom", "reps"], observed=True)["n_par"].max()
df = df.merge(npar.rename("n_par_ref").reset_index(), on=["case", "n_nom", "reps"], how="left")

# ---- logLik difference to asreml ---------------------------------------------------
# one asreml row per (case, n): after a rerun, the last one wins
# one row per (case, n, reps, series): a successful rerun replaces a failure
prio = df["status"].map({"ok": 0, "timeout": 1, "error": 2}).fillna(3)
df = df.assign(_prio=prio).sort_values(["case", "n_nom", "reps", "series", "_prio"], kind="stable")
dup = df.duplicated(["case", "n_nom", "reps", "series"], keep="first")
if dup.any():
    print("duplicate rows removed:", int(dup.sum()))
df = df[~dup].drop(columns="_prio").sort_values(["case", "n", "series"]).reset_index(drop=True)
ll_as = {(str(r["case"]), int(r["n_nom"]), int(r["reps"])): float(r["logLik"])
         for _, r in df[(df["series"] == "asreml") & (df["status"] == "ok")].iterrows()}
def dll(r):
    if r["status"] != "ok" or r["series"] == "asreml":
        return np.nan
    v = ll_as.get((str(r["case"]), int(r["n_nom"]), int(r["reps"])))
    return np.nan if v is None else float(r["logLik"]) - v
df["dlogLik_vs_asreml"] = df.apply(dll, axis=1)

long_cols = ["case", "n_nom", "n", "reps", "q", "n_par_ref", "series", "cores", "status", "wall_s", "startup_s", "compile_s", "solver_s",
             "n_iter", "n_eval", "s_per_iter", "s_per_eval", "eval_s", "converged", "logLik", "dlogLik_vs_asreml", "error"]
df[long_cols].to_csv(os.path.join(OUT, "benchc_%s_table.csv" % DATE), index=False)

# ---- wide markdown table ----------------------------------------------------------
def fmt(v, d=1):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ""
    return ("%%.%df" % d) % v
def cell_time(r):
    if r is None:
        return "not run"
    if r["status"] == "timeout":
        return "> %s (stopped)" % fmt(r["wall_s"], 0)
    if r["status"] != "ok":
        return "failed"
    return fmt(r["wall_s"], 1)
lines = ["# Benchmark at equal model: asreml, remlax CPU, remlax GPU (%s)" % DATE, "",
         "Times in seconds. `wall` is the complete call from R. For remlax, `startup` is the Python process, "
         "JAX import and bundle transfer; `compile` is XLA compilation, paid once; `solver` is the optimisation "
         "(L-BFGS-B and Newton polishing). `s/iter` is solver time per L-BFGS-B iteration (remlax) or per AI "
         "iteration (asreml); `s/eval` is per evaluation of -2 logL and its gradient. `dlogLik` is remlax minus "
         "asreml at the optimum (asreml convention); positive means remlax reached a higher REML likelihood. "
         "A fit stopped by the per-fit time limit is shown as `> limit`. CPU runs use 4 cores, except the density axis at n = 2000 (8 cores, `cores` column of the CSV); asreml is single-threaded for most of its work and the CPU scaling of remlax between 4 and 8 cores is below 20 %.", ""]
df_all = df
# size axis = the plan of the first campaign: 4 replicates (2 for the multi-trait cases)
is_us = df_all["case"].astype(str).str.startswith("us")
df = df_all[((df_all["reps"] == 4) & ~is_us) | ((df_all["reps"] == 2) & is_us)]
lines += ["## Compact view, size axis (4 replicates per genotype, 2 for the multi-trait cases): wall time (s) and solver time per iteration (s)", "",
          "| case | n | n_par | asreml wall | remlax CPU wall | remlax GPU wall | asreml s/iter | CPU s/iter | GPU s/iter | GPU compile | GPU dlogLik |", "|" + "---|" * 11]
for case_ in df["case"].cat.categories:
    sub = df[df["case"] == case_]
    for n in sorted(sub["n_nom"].unique()):
        s = sub[sub["n_nom"] == n]
        def get(series):
            x = s[s["series"] == series]
            return None if x.empty else x.iloc[0]
        a, cpu, gpu = get("asreml"), get("remlax CPU (4 cores)"), get("remlax GPU (A100)")
        def spi(r): return fmt(r["s_per_iter"], 3) if r is not None and r["status"] == "ok" else ""
        lines.append("| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            case_, n, fmt(s["n_par_ref"].max(), 0), cell_time(a), cell_time(cpu), cell_time(gpu), spi(a), spi(cpu), spi(gpu),
            fmt(gpu["compile_s"], 1) if gpu is not None and gpu["status"] == "ok" else "",
            ("%+.1e" % gpu["dlogLik_vs_asreml"]) if gpu is not None and gpu["status"] == "ok" and np.isfinite(gpu["dlogLik_vs_asreml"]) else ""))
lines.append("")
for case_ in df["case"].cat.categories:
    sub = df[df["case"] == case_]
    lines += ["## %s (%s)" % (case_, CASE_LABEL.get(case_, case_)), "",
              "| n | n_par | asreml wall | asreml iters | asreml s/iter | remlax CPU wall | CPU startup | CPU compile | CPU solver | CPU iters | CPU evals | CPU s/iter | CPU s/eval | CPU dlogLik | remlax GPU wall | GPU startup | GPU compile | GPU solver | GPU iters | GPU evals | GPU s/iter | GPU s/eval | GPU dlogLik |",
              "|" + "---|" * 23]
    for n in sorted(sub["n_nom"].unique()):
        s = sub[sub["n_nom"] == n]
        def get(series):
            x = s[s["series"] == series]
            return None if x.empty else x.iloc[0]
        a, cpu, gpu = get("asreml"), get("remlax CPU (4 cores)"), get("remlax GPU (A100)")
        row = [str(n), fmt(s["n_par_ref"].max(), 0), cell_time(a),
               fmt(a["n_iter"], 0) if a is not None and a["status"] == "ok" else "", fmt(a["s_per_iter"], 3) if a is not None and a["status"] == "ok" else ""]
        for r in (cpu, gpu):
            if r is None or r["status"] != "ok":
                row += [cell_time(r), "", "", "", "", "", "", "", ""]
            else:
                row += [fmt(r["wall_s"], 1), fmt(r["startup_s"], 1), fmt(r["compile_s"], 1), fmt(r["solver_s"], 1),
                        fmt(r["n_iter"], 0), fmt(r["n_eval"], 0), fmt(r["s_per_iter"], 4), fmt(r["s_per_eval"], 4),
                        ("%+.1e" % r["dlogLik_vs_asreml"]) if np.isfinite(r["dlogLik_vs_asreml"]) else ""]
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
# ---- density axis: fixed n, varying replicates per genotype ------------------------
dens = df_all[df_all["case"].isin(["grm", "usK3", "usK6", "ige"]) & (df_all["n_nom"].isin([2000, 8000]))]
if dens["reps"].nunique() > 1:
    lines += ["## Density axis: n fixed, replicates per genotype from 40 to 1 (q = genotypes)", "",
              "| case | n | reps | q | n_par | asreml wall | remlax CPU wall | remlax GPU wall | asreml s/iter | CPU s/iter | GPU s/iter | GPU compile | GPU solver | GPU dlogLik |", "|" + "---|" * 14]
    for case_ in [c_ for c_ in CASE_ORDER if c_ in set(dens["case"].astype(str))]:
        for n in sorted(dens["n_nom"].unique()):
            for reps in sorted(dens["reps"].unique(), reverse=True):
                s = dens[(dens["case"] == case_) & (dens["n_nom"] == n) & (dens["reps"] == reps)]
                if s.empty:
                    continue
                def get(series):
                    x = s[s["series"] == series]
                    return None if x.empty else x.iloc[0]
                a, cpu, gpu = get("asreml"), get("remlax CPU (4 cores)"), get("remlax GPU (A100)")
                def spi(r): return fmt(r["s_per_iter"], 3) if r is not None and r["status"] == "ok" else ""
                qv = s["q"].dropna()
                lines.append("| %s | %d | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                    case_, n, reps, fmt(qv.iloc[0], 0) if len(qv) else "", fmt(s["n_par_ref"].max(), 0), cell_time(a), cell_time(cpu), cell_time(gpu),
                    spi(a), spi(cpu), spi(gpu),
                    fmt(gpu["compile_s"], 1) if gpu is not None and gpu["status"] == "ok" else "",
                    fmt(gpu["solver_s"], 1) if gpu is not None and gpu["status"] == "ok" else "",
                    ("%+.1e" % gpu["dlogLik_vs_asreml"]) if gpu is not None and gpu["status"] == "ok" and np.isfinite(gpu["dlogLik_vs_asreml"]) else ""))
    lines.append("")
open(os.path.join(OUT, "benchc_%s_wide.md" % DATE), "w").write("\n".join(lines))

# ---- figures -------------------------------------------------------------------------
try:
    apply_figure_style()          # noqa: F821  (figure-style skill, session kernel)
except NameError:
    pass
COL = {"asreml": "#7f7f7f", "remlax CPU (4 cores)": "#1f77b4", "remlax GPU (A100)": "#d62728"}
MK = {"asreml": "s", "remlax CPU (4 cores)": "o", "remlax GPU (A100)": "^"}
cases = list(df["case"].cat.categories)
ok = df[df["status"] == "ok"]
to = df[df["status"] == "timeout"]


def _sup(fig, text, **kw):
    """Figure title wrapped to the figure width (about 1 character per 0.09 inch at 10 pt)."""
    width = max(40, int(fig.get_figwidth() / 0.105))
    kw.setdefault("x", 0.01); kw.setdefault("ha", "left")
    fig.suptitle("\n".join(textwrap.wrap(text, width)), **kw)

def _ticks_log(ax, axis="y"):
    f = mpl.ticker.FuncFormatter(lambda v, _: ("%g" % v) if v < 1000 else ("%gk" % (v / 1000)))
    (ax.yaxis if axis == "y" else ax.xaxis).set_major_formatter(f)

# A. total time against n, one panel per case
ncol = 5; nrow = int(np.ceil(len(cases) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(13, 2.9 * nrow), sharex=True, sharey=True, squeeze=False)
for i, case_ in enumerate(cases):
    ax = axes[i // ncol][i % ncol]
    for series, col in COL.items():
        s = ok[(ok["case"] == case_) & (ok["series"] == series)].sort_values("n_nom")
        if not s.empty:
            ax.plot(s["n_nom"], s["wall_s"], marker=MK[series], color=col, lw=1.4, ms=5, label=series)
        t = to[(to["case"] == case_) & (to["series"] == series)]
        if not t.empty:
            ax.scatter(t["n_nom"], t["wall_s"], marker="^" if MK[series] != "^" else "v", facecolors="none",
                       edgecolors=col, s=40, zorder=3)
            for _, r in t.iterrows():
                ax.annotate("", xy=(r["n_nom"], r["wall_s"] * 1.8), xytext=(r["n_nom"], r["wall_s"]),
                            arrowprops=dict(arrowstyle="->", color=col, lw=0.8))
    ax.set_xscale("log"); ax.set_yscale("log")
    _ticks_log(ax, "x"); _ticks_log(ax, "y")
    ax.set_title("%s: %s" % (case_, CASE_LABEL.get(case_, case_)), loc="left")
    ax.grid(True, which="major", lw=0.3, alpha=0.5)
    ax.margins(0.08)
for j in range(len(cases), nrow * ncol):
    axes[j // ncol][j % ncol].axis("off")
for r_ in range(nrow):
    axes[r_][0].set_ylabel("wall time of the fit (s)")
for c_ in range(ncol):
    axes[nrow - 1][c_].set_xlabel("observations n")
axes[0][0].legend(loc="upper left", frameon=False)
_sup(fig, "Total fit time at equal model and data: asreml on 4 CPU cores, remlax on 4 CPU cores, remlax on one A100. "
             "Open markers with arrows: stopped at the per-fit time limit.", x=0.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(os.path.join(OUT, "fig_benchc_wall.png"), dpi=200)

# B. time per iteration per case (ordered by n_par), one panel per n: points, no lines
ns = sorted(ok["n_nom"].unique())
fig2, axes2 = plt.subplots(1, len(ns), figsize=(3.3 * len(ns), 3.6), sharey=True, squeeze=False)
xpos = {case_: i for i, case_ in enumerate(cases)}
for k, n in enumerate(ns):
    ax = axes2[0][k]
    for series, col in COL.items():
        s = ok[(ok["n_nom"] == n) & (ok["series"] == series)]
        if s.empty:
            continue
        ax.scatter([xpos[c_] for c_ in s["case"].astype(str)], s["s_per_iter"], marker=MK[series], color=col, s=28, label=series, zorder=3)
    ax.set_yscale("log"); _ticks_log(ax, "y")
    ax.set_xticks(range(len(cases))); ax.set_xticklabels(cases, rotation=60, ha="right")
    ax.set_title("n = %s" % ("%gk" % (n / 1000) if n >= 1000 else n), loc="left")
    ax.grid(True, axis="y", lw=0.3, alpha=0.5); ax.margins(x=0.06, y=0.1)
axes2[0][0].set_ylabel("solver time per iteration (s)")
h, l = axes2[0][0].get_legend_handles_labels()
fig2.legend(h, l, loc="upper left", frameon=False, bbox_to_anchor=(0.01, 0.995), ncol=3)
_sup(fig2, "Time per iteration: asreml (average information) against remlax (L-BFGS-B); cases ordered by number of variance parameters",
              x=0.01, y=0.90, ha="left")
fig2.tight_layout(rect=(0, 0, 1, 0.86))
fig2.savefig(os.path.join(OUT, "fig_benchc_iter.png"), dpi=200)

# C. XLA compilation against n, per case (remlax only)
rx = ok[(ok["series"] != "asreml") & (ok["compile_s"] > 0.05)]
fig3, ax3 = plt.subplots(figsize=(6.5, 4))
cmap = plt.get_cmap("tab10")
for i, case_ in enumerate(cases):
    for series, ls in (("remlax CPU (4 cores)", "--"), ("remlax GPU (A100)", "-")):
        s = rx[(rx["case"] == case_) & (rx["series"] == series)].sort_values("n_nom")
        if s.empty:
            continue
        ax3.plot(s["n_nom"], s["compile_s"], ls=ls, marker=MK[series], color=cmap(i % 10), lw=1.0, ms=4,
                 label=case_ if series.endswith("(A100)") else None)
ax3.set_xscale("log"); ax3.set_yscale("log"); _ticks_log(ax3, "x"); _ticks_log(ax3, "y")
ax3.set_xlabel("observations n"); ax3.set_ylabel("XLA compilation time (s), paid once per fit")
ax3.set_title("XLA compilation, paid once per fit\nsolid: A100; dashed: 4 CPU cores", loc="left")
ax3.grid(True, lw=0.3, alpha=0.5); ax3.margins(0.08)
ax3.legend(frameon=False, ncol=2, title="case")
fig3.tight_layout()
fig3.savefig(os.path.join(OUT, "fig_benchc_compile.png"), dpi=200)

# D. startup / compilation / solver breakdown on A100, per case, at each n
gpu_ok = ok[ok["series"] == "remlax GPU (A100)"]
fig4, axes4 = plt.subplots(1, len(ns), figsize=(3.4 * len(ns), 3.9), sharey=False, squeeze=False)
parts = [("startup_s", "startup (Python, JAX, CUDA, bundle)", "#bdbdbd"), ("compile_s", "XLA compilation (once per fit)", "#fdae6b"), ("solver_s", "solver (L-BFGS-B + Newton)", "#d62728")]
for k, n in enumerate(ns):
    ax = axes4[0][k]
    x = np.arange(len(cases)); base = np.zeros(len(cases))
    for col, lab, colr in parts:
        vals = np.array([gpu_ok[(gpu_ok["case"] == case_) & (gpu_ok["n_nom"] == n)][col].sum() for case_ in cases])
        ax.bar(x, vals, 0.7, bottom=base, color=colr, edgecolor="white", lw=0.4, label=lab if k == 0 else None)
        base += vals
    for xi, tot in zip(x, base):
        if tot > 0:
            ax.text(xi, tot, "%.0f" % tot if tot >= 10 else "%.1f" % tot, ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x); ax.set_xticklabels(cases, rotation=60, ha="right")
    ax.set_title("n = %s" % ("%gk" % (n / 1000) if n >= 1000 else n), loc="left")
    ax.grid(True, axis="y", lw=0.3, alpha=0.5); ax.margins(y=0.12)
axes4[0][0].set_ylabel("wall time on one A100 (s)")
h, l = axes4[0][0].get_legend_handles_labels()
fig4.legend(h, l, loc="upper left", frameon=False, ncol=3, bbox_to_anchor=(0.01, 0.995))
_sup(fig4, "Where the A100 time goes: below n = 2k the fit is startup and compilation; the solver dominates from n = 8k",
              x=0.01, y=0.90, ha="left")
fig4.tight_layout(rect=(0, 0, 1, 0.86))
fig4.savefig(os.path.join(OUT, "fig_benchc_breakdown.png"), dpi=200)

# E. density axis: wall time against q at fixed n
if dens["reps"].nunique() > 1:
    dcases = [c_ for c_ in CASE_ORDER if c_ in set(dens["case"].astype(str))]
    dns = sorted(dens["n_nom"].unique())
    fig5, axes5 = plt.subplots(len(dns), len(dcases), figsize=(3.1 * len(dcases), 3.0 * len(dns)), sharex="col", sharey=True, squeeze=False)
    okd = dens[dens["status"] == "ok"]; tod = dens[dens["status"] == "timeout"]
    for i, n in enumerate(dns):
        for j, case_ in enumerate(dcases):
            ax = axes5[i][j]
            for series, col in COL.items():
                s_ = okd[(okd["case"] == case_) & (okd["n_nom"] == n) & (okd["series"] == series)].sort_values("q")
                if not s_.empty:
                    ax.plot(s_["q"], s_["wall_s"], marker=MK[series], color=col, lw=1.4, ms=5, label=series)
                t_ = tod[(tod["case"] == case_) & (tod["n_nom"] == n) & (tod["series"] == series)]
                if not t_.empty:
                    ax.scatter(t_["q"], t_["wall_s"], marker="v", facecolors="none", edgecolors=col, s=40, zorder=3)
            ax.set_xscale("log"); ax.set_yscale("log"); _ticks_log(ax, "x"); _ticks_log(ax, "y")
            ax.grid(True, lw=0.3, alpha=0.5); ax.margins(0.1)
            if i == 0:
                ax.set_title("%s: %s" % (case_, CASE_LABEL.get(case_, case_)), loc="left")
            if j == 0:
                ax.set_ylabel("n = %s\nwall time (s)" % ("%gk" % (n / 1000) if n >= 1000 else n))
            if i == len(dns) - 1:
                ax.set_xlabel("genotypes q (n / replicates)")
    axes5[0][0].legend(loc="upper left", frameon=False)
    _sup(fig5, "Density of the problem at fixed n: many genotypes and few replicates make asreml's equations large and dense; "
                  "remlax's solver works on the n x n matrix and its cost depends only weakly on q. "
                  "Open markers: stopped at the time limit.", x=0.01, ha="left")
    fig5.tight_layout(rect=(0, 0, 1, 0.90))
    fig5.savefig(os.path.join(OUT, "fig_benchc_density.png"), dpi=200)
print("written:", os.path.join(OUT, "benchc_%s_wide.md" % DATE), "+ figures;", len(df), "rows,",
      (df["status"] == "ok").sum(), "ok,", (df["status"] == "timeout").sum(), "timeout,", (df["status"] == "error").sum(), "error")