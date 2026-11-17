"""Figures and tables of diag_qtl_power.R: one-step against two-step GWAS of the major QTL.

    python3 studies/cogwas/fig_qtl_diag.py <diag_dir> <out_dir>

fig_qtl_chi2     test statistic at the causal SNP and at its best chip tag, per method
fig_qtl_power    detection rate against the share of the major QTL, per method
Tables: qtl_summary.csv (detection, statistics), qtl_calibration.csv (lambda, false peaks).
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

METHODS = [("p3d", "one-step, fixed V", "#9ecae1"), ("loco", "one-step, LOCO", "#4292c6"),
           ("exact", "one-step, exact refit", "#08519c"),
           ("blup_k", "two-step BLUP, with K", "#fdae6b"), ("blup_ols", "two-step BLUP, no K (inflated)", "#e6550d"),
           ("blue_k", "two-step BLUE, with K", "#bcbddc"), ("blue_ols", "two-step BLUE, no K (inflated)", "#6a51a3")]
LAB = {m: l for m, l, _ in METHODS}
COL = {m: c for m, _, c in METHODS}
TARGETS = [("DGE", "Direct effect (DGE)"), ("IGE_wheat", "IGE of wheat neighbours"),
           ("IGE_alfalfa", "IGE of alfalfa neighbours")]


def load(d, key):
    files = sorted(glob.glob(os.path.join(d, f"*_{key}.csv")))
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else pd.DataFrame()


def chi2(p):
    return stats.chi2.isf(np.clip(p, 1e-300, 1), 1)


def tables(q, c, out):
    q = q.assign(chi2_tag=chi2(q.p_tag), chi2_causal=chi2(q.p_causal))
    # diag_qtl_power.R took the family dose (0..2) for the alfalfa share, the plant dose is 0..4:
    # the family effect is beta x 2 x family dose, so the share is 4 times larger
    q.loc[q.target == "IGE_alfalfa", "share"] *= 4
    keys = ["n_lines", "pve_major", "r_g_target", "method", "target"]
    s = q.groupby(keys).agg(n_qtl=("found_tag", "size"), detected=("found_tag", "mean"),
                            detected_gc=("found_tag_gc", "mean"), detected_at_causal=("found_causal", "mean"),
                            chi2_causal_mean=("chi2_causal", "mean"), chi2_tag_mean=("chi2_tag", "mean"),
                            r2_tag_median=("r2_tag", "median"), share_median=("share", "median")).reset_index()
    cal = c.groupby(keys).agg(lambda_median=("lambda", "median"), false_peaks=("false_peaks", "sum"),
                              reps=("lambda", "size"), n_chip=("n_chip", "median")).reset_index()
    s = s.merge(cal, on=keys)
    s.to_csv(os.path.join(out, "qtl_summary.csv"), index=False)
    cal.to_csv(os.path.join(out, "qtl_calibration.csv"), index=False)
    return q, s


def fig_chi2(q, c, out, pve=0.10, rg=-0.6):
    d = q[(q.pve_major == pve) & (q.r_g_target == rg) & (q.n_lines == 181)]
    cc = c[(c.pve_major == pve) & (c.r_g_target == rg) & (c.n_lines == 181)]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.6), sharey="row")
    for j, (t, lab) in enumerate(TARGETS):
        dt = d[d.target == t]
        thr = chi2(0.05 / cc[cc.target == t].n_chip.median()) if len(cc[cc.target == t]) else np.nan
        ax = axes[0, j]
        for i, (m, _, col) in enumerate(METHODS):
            v = dt.loc[dt.method == m, ["chi2_causal", "chi2_tag"]].dropna()
            x = np.full(len(v), i) + np.random.default_rng(i).uniform(-0.18, 0.18, len(v))
            ax.scatter(x, v.chi2_causal, s=6, color=col, lw=0, alpha=0.8, zorder=3)
            ax.scatter(x, v.chi2_tag, s=8, facecolors="none", edgecolors=col, lw=0.5, zorder=3)
            if len(v):
                ax.plot([i - 0.3, i + 0.3], [v.chi2_causal.median()] * 2, color="0.2", lw=1, zorder=4)
        ax.axhline(thr, color="0.3", ls="--", lw=0.7, zorder=0)
        ax.set_ylim(bottom=0)
        ax.set_xticks(range(len(METHODS)))
        ax.set_xticklabels([])
        ax.set_title(lab, loc="left")
        ax.set_xlim(-0.6, len(METHODS) - 0.4)
        # lambda per replicate
        ax = axes[1, j]
        ct = cc[cc.target == t]
        for i, (m, _, col) in enumerate(METHODS):
            v = ct.loc[ct.method == m, "lambda"].dropna()
            x = np.full(len(v), i) + np.random.default_rng(i + 7).uniform(-0.18, 0.18, len(v))
            ax.scatter(x, v, s=6, color=col, lw=0, zorder=3)
            if len(v):
                ax.plot([i - 0.3, i + 0.3], [v.median()] * 2, color="0.2", lw=1, zorder=4)
        ax.axhline(1, color="0.3", ls="--", lw=0.7, zorder=0)
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
        ax.set_xticks(range(len(METHODS)))
        ax.set_xticklabels([LAB[m] for m, _, _ in METHODS], rotation=55, ha="right", fontsize=6)
        ax.set_xlim(-0.6, len(METHODS) - 0.4)
    axes[0, 0].set_ylabel("chi-square of the QTL")
    axes[1, 0].set_ylabel("genomic inflation lambda")
    panel_letter(axes[0, 0], "a", dx=-0.5)
    panel_letter(axes[1, 0], "b", dx=-0.5)
    h = [mpl.lines.Line2D([], [], ls="", marker="o", ms=3, color="0.4", label="at the causal SNP"),
         mpl.lines.Line2D([], [], ls="", marker="o", ms=3.5, mfc="none", mec="0.4", label="best chip SNP within 5 cM"),
         mpl.lines.Line2D([], [], ls="--", color="0.3", lw=0.7, label="Bonferroni 5 %")]
    fig.legend(handles=h, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    nrep = cc.groupby(["method", "target"]).size().min() if len(cc) else 0
    fig.text(0.01, 0.003, f"Major QTL at {pve:.2f} of their genetic component, r_g = {rg}; {nrep} replicates; "
             "181 lines, 106 families.", fontsize=6, color="0.35")
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    save(fig, out, "fig_qtl_chi2")


def fig_power(q, out, rg=0.0):
    d = q[(q.r_g_target == rg) & (q.n_lines == 181)]
    pves = sorted(d.pve_major.unique())
    if len(pves) < 2:
        return
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), sharey=True)
    for ax, (t, lab), l in zip(axes, TARGETS, "abc"):
        dt = d[d.target == t]
        for m, name, col in METHODS:
            r = dt[dt.method == m].groupby("pve_major").found_tag.mean().reindex(pves)
            ax.plot(pves, r.values, color=col, marker="o", ms=3, lw=1.3, label=name)
        ax.set_title(lab, loc="left")
        ax.set_xticks(pves)
        ax.set_xlabel("QTL share of its genetic variance", fontsize=6)
        ax.set_ylim(-0.04, 1.04)
        panel_letter(ax, l, dx=-0.3)
    axes[0].set_ylabel("detection rate (best chip SNP)")
    axes[-1].legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.text(0.01, 0.003, f"r_g = {rg}; 181 lines, 106 families; found: chip SNP within 5 cM below Bonferroni 5 %.",
             fontsize=6, color="0.35")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, out, "fig_qtl_power")


def fig_lines(q, out, pve=0.10, rg=-0.6):
    d = q[(q.pve_major == pve) & (q.r_g_target == rg)]
    nl = sorted(d.n_lines.unique())
    if len(nl) < 2:
        return
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), sharey=True)
    for ax, (t, lab), l in zip(axes, TARGETS, "abc"):
        dt = d[d.target == t]
        for m, name, col in METHODS:
            r = dt[dt.method == m].groupby("n_lines").found_tag.mean().reindex(nl)
            if m.startswith("blue"):
                r[r.index >= 724] = np.nan  # stage 1 saturated: 1,872 fixed effects for 1,920 plants
            ax.plot(range(len(nl)), r.values, color=col, marker="o", ms=3, lw=1.3, label=name)
        ax.set_xticks(range(len(nl)))
        ax.set_xticklabels([f"{n}\n{1920 / n:.1f}" for n in nl])
        ax.set_xlabel("lines\nplants per line", fontsize=6)
        ax.set_title(lab, loc="left")
        ax.set_ylim(-0.04, 1.04)
        panel_letter(ax, l, dx=-0.3)
    axes[0].set_ylabel("detection rate (best chip SNP)")
    axes[-1].legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.text(0.01, 0.003, f"1920 wheat plants; major QTL at {pve:.2f} of their genetic component, r_g = {rg}. "
             "BLUE not shown at 724 lines (stage 1 saturated).",
             fontsize=6, color="0.35")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, out, "fig_qtl_lines")


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    d, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    apply_figure_style()
    q, c = load(d, "qtl"), load(d, "cal")
    q, s = tables(q, c, out)
    fig_chi2(q, c, out)
    fig_power(q, out)
    fig_lines(q, out)
    print(s.to_string())
