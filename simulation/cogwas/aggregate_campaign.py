"""Aggregate the co-GWAS simulation campaign (run_cascade.R) into tables and figures.

    python3 simulation/cogwas/aggregate_campaign.py <campaign_dir> <out_dir>

Tables (CSV, and a Markdown summary):
  variance_recovery.csv   bias, relative bias, RMSE and 95 % CI coverage per
                          scenario x null model x parameter
  gwas_power.csv          power per effect type (found / simulated), binomial 95 % CI
  gwas_calibration.csv    lambda and false peaks per test
  cogwas.csv              pair power, lambda, significant pairs outside the causal windows
  exact_refit.csv         gain of the exact refit on the top SNPs
Figures: fig_variances, fig_gwas, fig_cogwas (PDF and PNG).
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.expanduser("~/.claude/skills/figure-style"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from figure_style import apply_figure_style, panel_letter  # noqa: E402

CAMP, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
SCEN = ["S0", "S1", "S2"]
SCEN_LABEL = {"S0": "additive only", "S1": "+ major interactions", "S2": "+ interaction background"}
VAR = ["no_kernel", "kernels"]
VAR_LABEL = {"no_kernel": "null without interaction kernels", "kernels": "null with interaction kernels"}
COL = {"no_kernel": "#8c8c8c", "kernels": "#1f5fa8"}


def read(kind):
    files = sorted(glob.glob(os.path.join(CAMP, "*_%s.csv" % kind)))
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else pd.DataFrame()


vc, vcomp, gw, gsum, pr, ex, truth = (read(k) for k in ("vc", "vcomp", "gw", "gsum", "pr", "ex", "truth"))
n_rep = truth.groupby("scenario")["rep"].nunique().to_dict()
print("replicates per scenario:", n_rep)

# ---------------------------------------------------------------- variances
NAME = {"var_DGE": "gen[1,1]", "var_IGE_intra": "gen[2,2]", "var_IGE_inter": "ige_inter",
        "var_bin": "bin", "var_row": "bin_row", "var_col": "bin_col", "var_residual": "residual",
        "var_k_inter": "k_inter", "var_k_intra": "k_intra", "var_k_sim": "k_sim"}
se = vcomp.set_index(["tag", "variant", "name"])["se"]
vc["se"] = [se.get((t, v, NAME[p]), np.nan) if p in NAME else np.nan
            for t, v, p in zip(vc["tag"], vc["variant"], vc["parameter"])]
cov = vcomp[vcomp["name"] == "gen[2,1]"][["tag", "variant", "value", "se"]].rename(columns={"value": "estimated"})
tcov = truth.assign(simulated=truth["r_g"] * np.sqrt(truth["var_dge"] * truth["var_ige"]))[["tag", "simulated"]]
cov = cov.merge(tcov, on="tag").merge(vc[["tag", "variant", "scenario", "rep", "converged"]].drop_duplicates(),
                                      on=["tag", "variant"])
cov["parameter"] = "cov_DGE_IGE"
vc = pd.concat([vc, cov[vc.columns]], ignore_index=True)
vc["error"] = vc["estimated"] - vc["simulated"]
vc["rel_error"] = np.where(vc["simulated"].abs() > 1e-12, vc["error"] / vc["simulated"].abs(), np.nan)
vc["covered"] = np.where(vc["se"].notna(), (vc["error"].abs() <= 1.96 * vc["se"]), np.nan)
rows = []
for (s, v, p), g in vc.groupby(["scenario", "variant", "parameter"]):
    g = g[g["converged"]]
    rows.append(dict(scenario=s, null_model=v, parameter=p, n=len(g), simulated_mean=g["simulated"].mean(),
                     estimated_mean=g["estimated"].mean(), bias=g["error"].mean(),
                     bias_se=g["error"].std(ddof=1) / np.sqrt(len(g)) if len(g) > 1 else np.nan,
                     relative_bias=g["rel_error"].mean(), rmse=np.sqrt((g["error"] ** 2).mean()),
                     coverage_95=g["covered"].mean() if g["covered"].notna().any() else np.nan))
var_tab = pd.DataFrame(rows)
var_tab.to_csv(os.path.join(OUT, "variance_recovery.csv"), index=False)


# ---------------------------------------------------------------- GWAS
def binom_ci(k, n):
    if n == 0:
        return (np.nan, np.nan)
    lo = stats.beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return lo, hi


pw = []
for (s, v, e), g in gw.groupby(["scenario", "variant", "effect"]):
    k, n = int(g["found"].sum()), len(g)
    lo, hi = binom_ci(k, n)
    pw.append(dict(scenario=s, null_model=v, effect=e, found=k, simulated=n, power=k / n, ci_low=lo, ci_high=hi))
gwas_power = pd.DataFrame(pw)
gwas_power.to_csv(os.path.join(OUT, "gwas_power.csv"), index=False)
gcal = gsum.groupby(["scenario", "variant", "test"]).agg(
    lambda_median=("lambda", "median"), lambda_q05=("lambda", lambda x: x.quantile(0.05)),
    lambda_q95=("lambda", lambda x: x.quantile(0.95)), reps=("rep", "nunique"),
    reps_with_false_peak=("false_peaks", lambda x: int((x > 0).sum())),
    false_peaks_total=("false_peaks", "sum")).reset_index()
gcal.to_csv(os.path.join(OUT, "gwas_calibration.csv"), index=False)

# ---------------------------------------------------------------- co-GWAS
cg = []
for (s, v, e), g in pr.groupby(["scenario", "variant", "effect"]):
    sim = g[g["simulated"]]
    k, n = int(sim["found"].sum()), len(sim)
    lo, hi = binom_ci(k, n)
    cg.append(dict(scenario=s, null_model=v, effect=e, power=k / n if n else np.nan, ci_low=lo, ci_high=hi,
                   simulated=n, lambda_median=g["lambda"].median(),
                   reps_with_sig_outside=int((g["sig_outside"] > 0).sum()), reps=len(g),
                   sig_outside_total=int(g["sig_outside"].sum())))
cog_tab = pd.DataFrame(cg)
cog_tab.to_csv(os.path.join(OUT, "cogwas.csv"), index=False)

ex_tab = pd.DataFrame()
if len(ex):
    ex["log10_gain"] = np.log10(ex["p_fixed_V"]) - np.log10(ex["p_exact"])
    ex_tab = ex.groupby(["scenario", "test"]).agg(n=("snp", "size"), median_gain_log10=("log10_gain", "median"),
                                                  converged=("converged", "mean")).reset_index()
    ex_tab.to_csv(os.path.join(OUT, "exact_refit.csv"), index=False)

# ---------------------------------------------------------------- figures
apply_figure_style(sizes=(8, 7, 6))
rng = np.random.default_rng(1)
OFFSET = {"no_kernel": -0.17, "kernels": 0.17}
MARK = {"S0": "o", "S1": "s", "S2": "^"}


def strip(ax, data, value, title, ylabel, zero=0.0):
    for i, s in enumerate(SCEN):
        for v in VAR:
            y = data[(data["scenario"] == s) & (data["variant"] == v)][value].dropna().to_numpy()
            if not len(y):
                continue
            x = i + OFFSET[v] + rng.uniform(-0.07, 0.07, len(y))
            ax.scatter(x, y, s=8, color=COL[v], alpha=0.55, linewidths=0, zorder=2, marker=MARK[s])
            ax.hlines(np.median(y), i + OFFSET[v] - 0.12, i + OFFSET[v] + 0.12, color="black", lw=1.4, zorder=3)
    ax.axhline(zero, color="#444444", lw=0.6, ls="--", zorder=1)
    ax.set_xticks(range(len(SCEN)))
    ax.set_xticklabels(SCEN)
    ax.set_xlim(-0.5, len(SCEN) - 0.5)
    ax.set_title(title, loc="left")
    ax.set_ylabel(ylabel)
    ax.margins(y=0.06)


MOFF = {("S0", "no_kernel"): -0.3, ("S0", "kernels"): -0.2, ("S1", "no_kernel"): -0.05,
        ("S1", "kernels"): 0.05, ("S2", "no_kernel"): 0.2, ("S2", "kernels"): 0.3}


def power_points(ax, tab, effects, labels, scen):
    """Power with its binomial 95 % CI: shape = scenario, colour = null model."""
    for i, e in enumerate(effects):
        for sc in scen:
            for v in VAR:
                r = tab[(tab["effect"] == e) & (tab["scenario"] == sc) & (tab["null_model"] == v)]
                if not len(r) or not np.isfinite(r["power"].iloc[0]):
                    continue
                x = i + MOFF[(sc, v)]
                ax.vlines(x, r["ci_low"].iloc[0], r["ci_high"].iloc[0], color=COL[v], lw=0.8)
                ax.plot(x, r["power"].iloc[0], MARK[sc], color=COL[v], ms=4, mec="white", mew=0.4)
    ax.set_xticks(range(len(effects)))
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.5, len(effects) - 0.5)
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("power (found / simulated)")


def scenario_legend(ax, scen, loc="upper right"):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker=MARK[sc], ls="", color="#555555", label=SCEN_LABEL[sc].replace("\n", " "))
         for sc in scen]
    ax.legend(handles=h, frameon=False, loc=loc)


def legend_models(fig, y=0.995):
    """One legend for the figure: colour = null model, shape = scenario."""
    from matplotlib.lines import Line2D
    top = [Line2D([], [], marker="o", ls="", color=COL[v], label=VAR_LABEL[v]) for v in VAR]
    top.append(Line2D([], [], color="black", lw=1.4, label="median"))
    bottom = [Line2D([], [], marker=MARK[sc], ls="", color="#555555", label="%s  %s" % (sc, SCEN_LABEL[sc]))
              for sc in SCEN]
    h = [x for pair in zip(top, bottom) for x in pair]     # legend fills by column: one row each
    fig.legend(handles=h, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, y),
               columnspacing=1.2, handletextpad=0.4)


# Figure 1: variances
d = vc[vc["converged"]].copy()
d["rel_pct"] = 100 * d["rel_error"]
panels = [("var_DGE", "rel_pct", "DGE variance", "relative error (%)"),
          ("var_IGE_intra", "rel_pct", "IGE variance (wheat)", "relative error (%)"),
          ("cov_DGE_IGE", "error", "DGE-IGE covariance", "estimate - simulated"),
          ("var_IGE_inter", "rel_pct", "IGE variance (alfalfa)", "relative error (%)"),
          ("var_residual", "rel_pct", "Residual variance", "relative error (%)"),
          ("var_k_inter", "error", "Kernel wheat x alfalfa", "estimate - simulated"),
          ("var_k_intra", "error", "Kernel wheat x wheat", "estimate - simulated"),
          ("var_k_sim", "error", "Similarity kernel", "estimate - simulated")]
fig, axes = plt.subplots(2, 4, figsize=(7.2, 4.6))
for ax, (p, val, title, yl), let in zip(axes.flat, panels, "abcdefgh"):
    strip(ax, d[d["parameter"] == p], val, title, yl)
    panel_letter(ax, let)
legend_models(fig)
fig.tight_layout(rect=(0, 0, 1, 0.91))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, "fig_variances.%s" % ext), dpi=300)
plt.close(fig)

# Figure 2: GWAS
EFF = ["major direct", "major indirect intra", "major indirect inter", "similarity"]
EFF_LABEL = {"major direct": "direct", "major indirect intra": "indirect\n(wheat)",
             "major indirect inter": "indirect\n(alfalfa)", "similarity": "similarity\n(same locus)"}
fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=[1.4, 1.2, 1]))
ax = axes[0]
power_points(ax, gwas_power, EFF, [EFF_LABEL[e] for e in EFF], SCEN)
ax.set_title("Major QTL found by the GWAS", loc="left")
panel_letter(ax, "a")
ax = axes[1]
gs = gsum.copy()
TESTS = ["dir", "ind", "sim|dir+ind", "inter ind"]
for i, t in enumerate(TESTS):
    for v in VAR:
        y = gs[(gs["test"] == t) & (gs["variant"] == v)]["lambda"].to_numpy()
        x = i + OFFSET[v] + rng.uniform(-0.07, 0.07, len(y))
        ax.scatter(x, y, s=6, color=COL[v], alpha=0.5, linewidths=0)
        if len(y):
            ax.hlines(np.median(y), i + OFFSET[v] - 0.12, i + OFFSET[v] + 0.12, color="black", lw=1.4)
ax.axhline(1, color="#444444", lw=0.6, ls="--")
ax.set_xticks(range(len(TESTS)))
ax.set_xticklabels(["direct", "indirect\n(wheat)", "similarity\n| dir+ind", "indirect\n(alfalfa)"])
ax.set_ylabel("genomic control lambda")
ax.set_title("Calibration of the GWAS tests", loc="left")
ax.margins(y=0.08)
panel_letter(ax, "b")
ax = axes[2]
if len(ex):
    for t, c in (("dir", "#1f5fa8"), ("ind", "#b5651d")):
        e = ex[ex["test"] == t]
        ax.scatter(-np.log10(e["p_fixed_V"]), -np.log10(e["p_exact"]), s=7, color=c, alpha=0.6, linewidths=0,
                   label="direct" if t == "dir" else "indirect (wheat)")
    lim = [0, max(-np.log10(ex["p_exact"]).max(), -np.log10(ex["p_fixed_V"]).max()) * 1.05]
    ax.plot(lim, lim, color="#444444", lw=0.6, ls="--")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.legend(frameon=False, loc="upper left")
ax.set_xlabel("-log10 p, fixed V")
ax.set_ylabel("-log10 p, exact refit")
ax.set_title("Exact refit of the top SNPs", loc="left")
panel_letter(ax, "c")
legend_models(fig, y=1.0)
fig.tight_layout(rect=(0, 0, 1, 0.83))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, "fig_gwas.%s" % ext), dpi=300)
plt.close(fig)

# Figure 3: co-GWAS
fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0))
ax = axes[0]
K = ["pair intra", "pair inter"]
power_points(ax, cog_tab, K, ["wheat x wheat\npair", "wheat x alfalfa\npair"], ["S1", "S2"])
ax.set_title("Simulated pairs found", loc="left")
panel_letter(ax, "a")
for ax, kind, let in ((axes[1], "pair intra", "b"), (axes[2], "pair inter", "c")):
    strip(ax, pr[pr["effect"] == kind], "lambda",
          "Calibration, %s" % ("wheat x wheat" if kind == "pair intra" else "wheat x alfalfa"),
          "genomic control lambda", zero=1.0)
    panel_letter(ax, let)
legend_models(fig, y=1.0)
fig.tight_layout(rect=(0, 0, 1, 0.83))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, "fig_cogwas.%s" % ext), dpi=300)
plt.close(fig)

# ---------------------------------------------------------------- summary
def md(tab):
    """Markdown table without the tabulate dependency."""
    cols = list(tab.columns)
    fmt = lambda x: ("%.4g" % x) if isinstance(x, (float, np.floating)) else str(x)  # noqa: E731
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in tab.itertuples(index=False)]
    return "\n".join(lines)


with open(os.path.join(OUT, "summary.md"), "w") as f:
    f.write("# Co-GWAS simulation campaign\n\nReplicates per scenario: %s\n\n" % n_rep)
    for name, tab in (("Variance recovery", var_tab), ("GWAS power", gwas_power), ("GWAS calibration", gcal),
                      ("Co-GWAS", cog_tab), ("Exact refit", ex_tab)):
        f.write("## %s\n\n%s\n\n" % (name, md(tab) if len(tab) else "none"))
print("written to", OUT)
