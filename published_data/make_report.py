"""Pool published_data/results/*.csv into docs/published-data.md and a figure.

Each script of published_data/scripts writes one checks table with the columns
of R/common.R (dataset, source, model, quantity, software, kind, remlax,
reference, metric, diff, tol, pass, note). This script reads them all, writes
the summary page and draws one panel per data set: every counted check as its
discrepancy divided by its tolerance, so that 1 is the pass line whatever the
quantity.

    python3 published_data/make_report.py
"""
from __future__ import annotations

import glob
import math
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(ROOT)
OUT_MD = os.path.join(REPO, "docs", "published-data.md")
OUT_FIG = os.path.join(ROOT, "results", "fig_published_data.png")


def load() -> pd.DataFrame:
    files = sorted(glob.glob(os.path.join(ROOT, "results", "*.csv")))
    files = [f for f in files if not f.endswith("_timing.csv")]
    if not files:
        sys.exit("no results/*.csv")
    tab = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    tab["pass"] = tab["pass"].map({True: True, False: False, "TRUE": True, "FALSE": False})
    return tab


def ratio(row) -> float:
    """Discrepancy over tolerance; for 'ge' checks (correlations) the gap to 1
    over the gap allowed."""
    if row["metric"] == "ge":
        # Only correlations (tolerance just below 1) have a discrepancy scale;
        # the other "ge" checks are signs or booleans and stay off the figure.
        if not (0.99 < row["tol"] < 1):
            return math.nan
        return (1.0 - row["remlax"]) / (1.0 - row["tol"])
    if row["tol"] and row["tol"] > 0:
        return row["diff"] / row["tol"]
    return math.nan


def fmt(x) -> str:
    if pd.isna(x):
        return ""
    if isinstance(x, (int,)) or (isinstance(x, float) and x.is_integer() and abs(x) < 1e6):
        return f"{x:.0f}"
    return f"{x:.4g}"


def summary_md(tab: pd.DataFrame) -> str:
    counted = tab[tab["kind"] != "context"]
    lines = []
    lines.append("| data set | source | checks | pass | fail | context rows |")
    lines.append("|---|---|---|---|---|---|")
    for ds, g in tab.groupby("dataset", sort=False):
        c = g[g["kind"] != "context"]
        lines.append(f"| {ds} | {g['source'].iloc[0]} | {len(c)} | {int(c['pass'].sum())} | "
                     f"{int((~c['pass'].astype(bool)).sum())} | {int((g['kind'] == 'context').sum())} |")
    lines.append(f"| total | | {len(counted)} | {int(counted['pass'].sum())} | "
                 f"{int((~counted['pass'].astype(bool)).sum())} | {int((tab['kind'] == 'context').sum())} |")
    return "\n".join(lines)


def by_software_md(tab: pd.DataFrame) -> str:
    counted = tab[tab["kind"] != "context"]
    lines = ["| reference | checks | pass | largest discrepancy / tolerance |", "|---|---|---|---|"]
    for sw, g in counted.groupby("software"):
        r = g.apply(ratio, axis=1)
        lines.append(f"| {sw} | {len(g)} | {int(g['pass'].sum())} | {fmt(r.max())} |")
    return "\n".join(lines)


def rows_md(g: pd.DataFrame) -> str:
    lines = ["| data set | model | quantity | reference | remlax | reference value | discrepancy | note |",
             "|---|---|---|---|---|---|---|---|"]
    for _, r in g.iterrows():
        note = "" if pd.isna(r["note"]) else str(r["note"]).replace("|", "/")
        disc = 1.0 - r["remlax"] if r["metric"] == "ge" and 0.99 < r["tol"] < 1 else r["diff"]
        val = f"{r['remlax']:.7f}" if r["metric"] == "ge" and 0.99 < r["tol"] < 1 else fmt(r["remlax"])
        lines.append(f"| {r['dataset']} | {r['model']} | {r['quantity']} | {r['software']} | "
                     f"{val} | {fmt(r['reference'])} | {fmt(disc)} | {note} |")
    return "\n".join(lines)


def figure(tab: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    counted = tab[tab["kind"] != "context"].copy()
    counted["ratio"] = counted.apply(ratio, axis=1)
    counted = counted[counted["ratio"].notna()]
    counted["lr"] = np.log10(np.clip(counted["ratio"].astype(float), 1e-12, 1e6))
    sets = list(dict.fromkeys(counted["dataset"]))
    softs = sorted(counted["software"].unique())
    cmap = plt.get_cmap("tab10")
    col = {s: cmap(i % 10) for i, s in enumerate(softs)}
    fig, ax = plt.subplots(figsize=(8.5, 0.55 * len(sets) + 1.6))
    rng = np.random.default_rng(1)
    for i, ds in enumerate(sets):
        g = counted[counted["dataset"] == ds]
        y = i + rng.uniform(-0.25, 0.25, len(g))
        ax.scatter(g["lr"], y, s=14, c=[col[s] for s in g["software"]], alpha=0.8,
                   edgecolors="none")
    ax.axvline(0, color="0.3", lw=1, ls="--")
    ax.set_yticks(range(len(sets)))
    ax.set_yticklabels(sets)
    ax.invert_yaxis()
    ax.set_xlabel("log10(discrepancy / tolerance), pass left of the dashed line; "
                  "values below 1e-12 drawn at -12")
    for s in softs:
        ax.scatter([], [], c=[col[s]], s=14, label=s)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=200)
    print("figure:", OUT_FIG)


def main() -> None:
    tab = load()
    counted = tab[tab["kind"] != "context"]
    fails = counted[~counted["pass"].astype(bool)]
    ctx = tab[tab["kind"] == "context"]
    pub = tab[tab["kind"] == "published"]
    md = []
    md.append("# Validation on published data sets\n")
    md.append("<!-- Generated by published_data/make_report.py from published_data/results/*.csv. "
              "Do not edit by hand. -->\n")
    md.append("Each script of `published_data/scripts/` fits remlax on a published data set and "
              "fits the same model with established programs. A check passes when the "
              "discrepancy is within the tolerance of `published_data/R/common.R`. Rows of "
              "kind `context` compare with values from the source paper obtained under a "
              "different model or program. They are reported and never counted.\n")
    md.append("## Summary\n")
    md.append(summary_md(tab) + "\n")
    md.append("## By reference program\n")
    md.append(by_software_md(tab) + "\n")
    md.append("## Failed checks\n")
    md.append(rows_md(fails) + "\n" if len(fails) else "None.\n")
    md.append("## Values printed in the source papers\n")
    md.append(rows_md(pub) + "\n" if len(pub) else "None.\n")
    md.append("## Context rows\n")
    md.append(rows_md(ctx) + "\n" if len(ctx) else "None.\n")
    md.append("## Rerun\n")
    md.append("```\nfor s in published_data/scripts/0[1-9]_*.R; do Rscript \"$s\"; done\n"
              "python3 published_data/make_report.py\n```\n")
    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(md))
    print("report:", OUT_MD)
    figure(tab)


if __name__ == "__main__":
    main()
