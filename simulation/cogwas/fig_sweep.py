"""Figures of the design x effect-size sweep (run_cascade.R, interaction null model).

    python3 simulation/cogwas/fig_sweep.py <sweep_dir> <out_dir>

The number of wheat plants is fixed (1,920); fewer genotypes means more
plants per genotype and more repeats of each pair of genotypes.
fig_sweep_recovery  BLUP accuracy and log relative error of the variances
fig_sweep_power     detection rate of each simulated QTL and pair
Also writes sweep_summary.csv.
"""
import glob
import os
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.expanduser("~/.claude/skills/figure-style"))
from figure_style import apply_figure_style, panel_letter  # noqa: E402

N_PLANTS = 1920
EFFECTS = [("DGE", "DGE", "var_DGE"), ("IGE intra", "IGE wheat", "var_IGE_intra"),
           ("IGE inter", "IGE alfalfa", "var_IGE_inter"), ("GxG intra", "DGE x IGE wheat", "var_GxG_intra"),
           ("GxG inter", "DGE x IGE alfalfa", "var_GxG_inter")]
DETECT = [("gw", "major direct", "Major direct QTL"), ("gw", "major indirect intra", "Major IGE QTL, wheat"),
          ("gw", "major indirect inter", "Major IGE QTL, alfalfa"), ("gw", "similarity", "Similarity QTL"),
          ("pr", "pair intra", "Interaction pair, wheat"), ("pr", "pair inter", "Interaction pair, alfalfa")]


def load(run, key):
    files = sorted(glob.glob(os.path.join(run, f"*_{key}.csv")))
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else pd.DataFrame()


def palette(levels):
    cm = mpl.colormaps["Blues"]
    return {g: cm(0.4 + 0.55 * k / max(1, len(levels) - 1)) for k, g in enumerate(levels)}


def setup_x(ax, lines):
    ax.set_xticks(range(len(lines)))
    ax.set_xticklabels([f"{n}\n{N_PLANTS / n:.0f}" for n in lines])
    ax.set_xlim(-0.4, len(lines) - 0.6)


def series(ax, df, value, lines, cols, floor=None):
    """Replicates as small points, medians joined by a line, one series per effect level."""
    levels = sorted(cols)
    for k, g in enumerate(levels):
        d = df[df.var_gxg == g]
        off = (k - (len(levels) - 1) / 2) * 0.12
        med = []
        for i, n in enumerate(lines):
            v = d.loc[d.n_lines == n, value].dropna().to_numpy()
            if floor is not None:
                low = v < floor
                if low.any():
                    ax.scatter(np.full(low.sum(), i + off), np.full(low.sum(), floor), marker="v", s=12,
                               facecolors="none", edgecolors=cols[g], lw=0.7, clip_on=False, zorder=3)
                v = np.maximum(v, floor)
                ax.scatter(np.full((~low).sum(), i + off), v[~low], s=5, color=cols[g], alpha=0.6, lw=0, zorder=3)
            else:
                ax.scatter(np.full(len(v), i + off), v, s=5, color=cols[g], alpha=0.6, lw=0, zorder=3)
            med.append(np.median(v) if len(v) else np.nan)
        ax.plot(np.arange(len(lines)) + off, med, color=cols[g], lw=1.3, marker="o", ms=3, zorder=4)
    setup_x(ax, lines)


def fig_recovery(run, out):
    acc = load(run, "acc"); vc = load(run, "vc")
    lines = sorted(acc.n_lines.unique(), reverse=True)
    cols = palette(sorted(acc.var_gxg.unique()))
    vc = vc.copy(); vc["lre"] = np.log(vc.estimated / vc.simulated)
    fig, axes = plt.subplots(2, 5, figsize=(7.2, 4.4), sharey="row")
    for j, (eff, lab, par) in enumerate(EFFECTS):
        series(axes[0, j], acc[acc.effect == eff], "accuracy", lines, cols)
        axes[0, j].set_title(lab, loc="left")
        axes[0, j].axhline(0, color="0.75", lw=0.5, zorder=0)
        series(axes[1, j], vc[vc.parameter == par], "lre", lines, cols, floor=-4)
        axes[1, j].axhline(0, color="0.3", lw=0.7, ls="--", zorder=0)
        axes[1, j].set_ylim(-4.3, 3)
    axes[0, 0].set_ylabel("BLUP accuracy")
    axes[1, 0].set_ylabel("log(estimated / simulated)\nvariance")
    for ax in axes[1]:
        ax.set_xlabel("lines\nplants per line", fontsize=6)
    for ax, l in zip(axes[:, 0], "ab"):
        panel_letter(ax, l, dx=-0.55)
    handles = [mpl.lines.Line2D([], [], color=cols[g], marker="o", ms=3, lw=1.3, label=f"{g:.2f}")
               for g in sorted(cols)]
    fig.legend(handles=handles, title="DGE x IGE variance per neighbour (x var DGE)", ncol=3,
               loc="upper center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    nrep = acc.groupby(["n_lines", "var_gxg"]).rep.nunique().min()
    fig.text(0.01, 0.003, f"{N_PLANTS} wheat plants in every design; {nrep} replicates per cell; "
             "open triangles: variance estimated at zero", fontsize=6, color="0.35")
    fig.tight_layout(rect=(0, 0.02, 1, 0.9))
    save(fig, out, "fig_sweep_recovery")


def fig_power(run, out):
    gw = load(run, "gw"); pr = load(run, "pr")
    src = {"gw": gw, "pr": pr}
    lines = sorted(gw.n_lines.unique(), reverse=True)
    cols = palette(sorted(gw.var_gxg.unique()))
    rows = []
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.2), sharey=True)
    for ax, (k, eff, lab), l in zip(axes.flat, DETECT, "abcdef"):
        d = src[k]; d = d[d.effect == eff]
        for j, g in enumerate(sorted(cols)):
            off = (j - 1) * 0.12
            rate, n_ = [], []
            for n in lines:
                f = d.loc[(d.n_lines == n) & (d.var_gxg == g), "found"].astype(float)
                rate.append(f.mean() if len(f) else np.nan); n_.append(len(f))
                rows.append(dict(effect=lab, n_lines=n, var_gxg=g, found=f.sum(), simulated=len(f)))
            ax.plot(np.arange(len(lines)) + off, rate, color=cols[g], marker="o", ms=3, lw=1.3)
        setup_x(ax, lines)
        ax.set_ylim(-0.04, 1.04)
        ax.set_title(lab, loc="left")
        panel_letter(ax, l, dx=-0.3)
    for ax in axes[:, 0]:
        ax.set_ylabel("detection rate")
    for ax in axes[1]:
        ax.set_xlabel("lines\nplants per line", fontsize=6)
    handles = [mpl.lines.Line2D([], [], color=cols[g], marker="o", ms=3, lw=1.3, label=f"{g:.2f}")
               for g in sorted(cols)]
    fig.legend(handles=handles, title="DGE x IGE variance per neighbour (x var DGE)", ncol=3,
               loc="upper center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    fig.text(0.01, 0.003, "Found: a chip SNP (or SNP pair) within 5 cM of the simulated locus passes Bonferroni 5 %. "
             "Additive major QTL: 10 % of their component.", fontsize=6, color="0.35")
    fig.tight_layout(rect=(0, 0.02, 1, 0.9))
    save(fig, out, "fig_sweep_power")
    pd.DataFrame(rows).to_csv(os.path.join(out, "sweep_power.csv"), index=False)


def summary(run, out):
    vc = load(run, "vc"); acc = load(run, "acc"); gs = load(run, "gsum"); pr = load(run, "pr")
    vc = vc.assign(lre=np.log(vc.estimated / vc.simulated), zero=vc.estimated < 1e-4 * vc.simulated.abs())
    a = vc.groupby(["n_lines", "var_gxg", "parameter"]).agg(simulated=("simulated", "mean"),
                                                             estimated=("estimated", "mean"),
                                                             median_lre=("lre", "median"),
                                                             at_zero=("zero", "sum"),
                                                             converged=("converged", "mean")).reset_index()
    a.to_csv(os.path.join(out, "sweep_variances.csv"), index=False)
    acc.groupby(["n_lines", "var_gxg", "effect"]).accuracy.median().unstack().to_csv(
        os.path.join(out, "sweep_accuracy.csv"))
    lam = pd.concat([gs.rename(columns={"test": "what"})[["n_lines", "var_gxg", "what", "lambda", "false_peaks"]],
                     pr[pr.effect.str.startswith("outside")].rename(columns={"effect": "what", "sig_outside": "false_peaks"})
                     [["n_lines", "var_gxg", "what", "lambda", "false_peaks"]]])
    lam.groupby(["what", "n_lines", "var_gxg"]).agg(lambda_median=("lambda", "median"),
                                                   false_peaks=("false_peaks", "sum")).to_csv(
        os.path.join(out, "sweep_calibration.csv"))


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    run, out = sys.argv[1], sys.argv[2]
    apply_figure_style()
    fig_recovery(run, out)
    fig_power(run, out)
    summary(run, out)
    print("written to", out)
