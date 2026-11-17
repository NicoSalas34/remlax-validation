"""Repetitions x QTL share grid at 181 lines (diag_qtl_power.R --bins --n_major 2 --null additive).

    python3 studies/cogwas/fig_qtl_grid.py <grid_dir> <out_dir>

fig_qtl_grid   detection rate (best chip SNP, one-step at fixed V) against plants per line,
               one series per QTL share, with the rule chi2 - 1 = N s R r2 (lines)
Table: qtl_grid_summary.csv.
"""
import glob
import os
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.expanduser("~/.claude/skills/figure-style"))
from figure_style import apply_figure_style, panel_letter  # noqa: E402

N_LINES, N_FAM, NOISE = 181, 106, 3.4
TARGETS = [("DGE", "Direct effect (DGE)"), ("IGE_wheat", "IGE of wheat neighbours"),
           ("IGE_alfalfa", "IGE of alfalfa neighbours")]


def load(d, key):
    files = sorted(glob.glob(os.path.join(d, f"*_{key}.csv")))
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else pd.DataFrame()


def reliability(target, bins, v_alf):
    """Own reliability of a genotype effect: sigma2 / (sigma2 + noise / observations per genotype)."""
    n_w = 160 * bins
    if target == "DGE":
        return 1.0 / (1.0 + NOISE / (n_w / N_LINES))
    if target == "IGE_wheat":
        return 0.1 / (0.1 + NOISE / (2 * n_w / N_LINES))
    return v_alf / (v_alf + NOISE / (6 * n_w / N_FAM))


def predicted(row, thr, v_alf):
    n = N_LINES if row.target != "IGE_alfalfa" else N_FAM
    ncp = n * row.share * reliability(row.target, row.bins, v_alf) * (row.r2_tag if np.isfinite(row.r2_tag) else 1)
    return stats.ncx2.sf(thr, 1, ncp)


def main(d, out):
    q, c, tr = load(d, "qtl"), load(d, "cal"), load(d, "truth")
    q = q[q.method == "p3d"].copy()
    q.loc[q.target == "IGE_alfalfa", "share"] *= 4      # family dose 0..2 against plant dose 0..4
    v_alf = tr.var_inter_fam.median()
    thr = {t: stats.chi2.isf(0.05 / c[(c.target == t)].n_chip.median(), 1) for t, _ in TARGETS}
    q["pred"] = [predicted(r, thr[r.target], v_alf) for r in q.itertuples()]
    q["plants_per_line"] = 160 * q.bins / N_LINES
    s = q.groupby(["target", "bins", "plants_per_line", "pve_major"]).agg(
        n_qtl=("found_tag", "size"), detected=("found_tag", "mean"), detected_at_causal=("found_causal", "mean"),
        predicted=("pred", "mean"), r2_tag=("r2_tag", "median")).reset_index()
    cal = c[c.method == "p3d"].groupby(["target", "bins", "pve_major"]).agg(
        lambda_median=("lambda", "median"), false_peaks=("false_peaks", "sum"), reps=("lambda", "size")).reset_index()
    s = s.merge(cal, on=["target", "bins", "pve_major"])
    os.makedirs(out, exist_ok=True)
    s.to_csv(os.path.join(out, "qtl_grid_summary.csv"), index=False)

    shares = sorted(s.pve_major.unique())
    cm = mpl.colormaps["Blues"]
    col = {p: cm(0.35 + 0.6 * k / max(1, len(shares) - 1)) for k, p in enumerate(shares)}
    ppl = sorted(s.plants_per_line.unique())
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), sharey=True)
    for ax, (t, lab), l in zip(axes, TARGETS, "abc"):
        for p in shares:
            d_ = s[(s.target == t) & (s.pve_major == p)].sort_values("plants_per_line")
            ax.plot(d_.plants_per_line, d_.predicted, color=col[p], lw=1, ls="--", zorder=2)
            lo, hi = zip(*[stats.beta.ppf([0.025, 0.975], k + 0.5, n - k + 0.5)
                           for k, n in zip((d_.detected * d_.n_qtl).round(), d_.n_qtl)])
            ax.errorbar(d_.plants_per_line, d_.detected, yerr=[np.clip(d_.detected - np.array(lo), 0, None), np.clip(np.array(hi) - d_.detected, 0, None)],
                        color=col[p], marker="o", ms=3.5, lw=1.3, capsize=0, zorder=3)
        ax.axhline(0.8, color="0.6", lw=0.6, ls=":", zorder=0)
        ax.set_xscale("log", base=2)
        ax.set_xticks(ppl)
        ax.set_xticklabels([f"{v:.1f}\n{int(round(v * N_LINES)):,}".replace(",", " ") for v in ppl])
        ax.set_xlabel("plants per line\nwheat plants", fontsize=6)
        ax.set_title(lab, loc="left")
        ax.set_ylim(-0.04, 1.04)
        panel_letter(ax, l, dx=-0.3)
    axes[0].set_ylabel("detection rate (best chip SNP)")
    h = [mpl.lines.Line2D([], [], color=col[p], marker="o", ms=3.5, lw=1.3, label=f"{p:.2f}") for p in shares]
    h += [mpl.lines.Line2D([], [], color="0.4", lw=1, ls="--", label="rule N s R r²")]
    axes[-1].legend(handles=h, title="QTL share", frameon=False, fontsize=6, title_fontsize=6,
                    loc="upper left", bbox_to_anchor=(1.0, 1.0))
    nrep = int(s.reps.min())
    fig.text(0.01, 0.003, f"181 lines, 106 families; 2 major QTL per component, r_g = 0; one-step test at fixed V; "
             f"{nrep} or more replicates per cell; bars: 95 % interval; dotted line: 80 %.", fontsize=6, color="0.35")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"fig_qtl_grid.{ext}"), dpi=300)
    plt.close(fig)
    print(s.to_string())


if __name__ == "__main__":
    apply_figure_style()
    main(sys.argv[1], sys.argv[2])
