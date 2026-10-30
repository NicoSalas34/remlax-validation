#!/usr/bin/env python3
"""Validation figure: remlax estimates against the reference, all checks.

    python3 validation/make_validation_fig.py --csv validation/results/checks_<date>.csv \
            --out validation/results/fig_validation.png [--pairs-out validation/results/pairs_<date>.csv]

The (remlax, reference) pairs are read from the `detail` column of the CSV,
which quotes the log line as it is. Three notations coexist there and are
recognised:
    "a vs b"                       one pair, remlax first
    "remlax a b c | <ref> d e f"   two lists of the same length
    "s2g a/b | s2e c/d"            remlax/reference fractions
A pair is kept only if the label does not explicitly say otherwise
("asreml a <= remlax b" is an ordering, not a pair; "a <= b" is ignored).
The right panel shows the -2logL (or logLik) gaps for the checks whose label
contains "logL".

Every pair kept is written to --pairs-out, so that the figure can be checked
row by row.

Check labels are matched in French (stored logs before the English pass of
the suites) and in English.
"""
import argparse
import csv
import os
import re
from collections import OrderedDict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"
RX_VS = re.compile(r"(?<![\w.])(%s)\s+vs\.?\s+(%s)(?![\w.])" % (NUM, NUM))
RX_FRAC = re.compile(r"(?<![\w./])(%s)/(%s)(?![\w./])" % (NUM, NUM))
RX_LISTS = re.compile(r"remlax(?: theta| phi)?\s+((?:%s\s*)+)\|\s*(?:asreml|nlme|lme4|sommer)(?: phi| theta)?\s+((?:%s\s*)+)" % (NUM, NUM))
RX_BAR = re.compile(r"(?:remlkit|remlax)\s+(%s)\s*\|\s*(?:R|asreml)\s+(%s)" % (NUM, NUM))

SOFTWARE = OrderedDict([
    ("asreml", ("asreml_structures", "asreml_catalogue", "asreml_rest", "asreml_blup", "ige_model")),
    ("lme4", ("lme4_sommer", "lme4_extra")),
    ("sommer", ("sommer_extra",)),
    ("nlme", ("nlme",)),
    ("closed form / internal", ("algebra", "sparse_parity", "sparse_vs_dense", "sparse_prec", "pytest")),
])
COLOUR = {"asreml": "#1f77b4", "lme4": "#2ca02c", "sommer": "#d62728", "nlme": "#9467bd",
           "closed form / internal": "#7f7f7f"}


def software_of(suite, check):
    for sw, suites in SOFTWARE.items():
        if suite in suites:
            if suite == "lme4_sommer":
                return "sommer" if "sommer" in check.lower() else "lme4"
            return sw
    return "other"


# Labels whose comparison is not remlax against a reference. French keywords
# match the stored logs, English ones the logs of the translated suites.
EXCLUDE = ("ameliore", "signe de", "emboitement", "domine", "moins bon", "constante",
           "a bouge", "decroit", "repliqu", "==", "identique a lme4",
           "voie explicite", "backend", "annonce", "ddl du terme", "degeneree", "disponible",
           "entre dans", "separees", "sections ameliorent", "GONFLE", "n'atteint pas",
           "improve", "sign of", "nesting", "dominate", "worse", "constant",
           "moved", "decreas", "replicat", "identical to lme4",
           "explicit route", "announced", "df of the term", "degenerate", "available",
           "enters", "separate", "INFLATE", "does not reach")
# Checks that compare sorted sets of values (French "ensemble", English "set").
SET_KEYS = ("ensemble", "(set)")


def pairs_in(detail, check):
    """(remlax, reference) pairs read from `detail`, or [] if the label says the
    comparison is not remlax against a reference (two remlax models, a
    simulated truth, an ordering)."""
    d = detail.replace("\\|", "|")
    out = []
    if "<=" in d or ">=" in d or any(k.lower() in check.lower() for k in EXCLUDE):
        return out
    m = RX_LISTS.search(d)
    if m:
        a = [float(x) for x in m.group(1).split()]
        b = [float(x) for x in m.group(2).split()]
        if any(k in check for k in SET_KEYS):  # the check compares sorted sets
            a, b = sorted(a), sorted(b)
        if len(a) == len(b):
            return list(zip(a, b))
    m = RX_BAR.search(d)
    if m:
        return [(float(m.group(1)), float(m.group(2)))]
    for m in RX_VS.finditer(d):
        out.append((float(m.group(1)), float(m.group(2))))
    if out:
        return out
    for m in RX_FRAC.finditer(d):
        out.append((float(m.group(1)), float(m.group(2))))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pairs-out", default=None)
    a = ap.parse_args(argv)
    rows = list(csv.DictReader(open(a.csv, newline="")))

    pts, gaps = [], []
    for r in rows:
        if r["verdict"] != "pass" and r["verdict"] != "FAIL":
            continue
        if r["suite"] in ("stress", "pytest", "cpu_gpu_parity"):
            continue                         # counters, not remlax / reference pairs
        sw = software_of(r["suite"], r["check"])
        is_ll = "logl" in r["check"].lower() or "loglik" in r["detail"].lower()
        for x, y in pairs_in(r["detail"], r["check"]):
            rec = dict(suite=r["suite"], check=r["check"], software=sw, remlax=x, reference=y,
                       kind="logLik" if is_ll else "estimate")
            if is_ll:
                gaps.append(rec)
            else:
                pts.append(rec)
    if a.pairs_out:
        with open(a.pairs_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["kind", "software", "suite", "check", "remlax", "reference"])
            w.writeheader()
            for rec in pts + gaps:
                w.writerow({k: rec[k] for k in w.fieldnames})

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [1.15, 1]})
    for sw in list(SOFTWARE) + ["other"]:
        P = [p for p in pts if p["software"] == sw]
        if not P:
            continue
        x = np.array([p["reference"] for p in P]); y = np.array([p["remlax"] for p in P])
        ax1.scatter(x, y, s=18, alpha=0.75, color=COLOUR.get(sw, "k"), label="%s (%d)" % (sw, len(P)), zorder=3)
    allv = np.array([p["reference"] for p in pts] + [p["remlax"] for p in pts])
    lo, hi = np.nanmin(allv), np.nanmax(allv)
    ax1.plot([lo, hi], [lo, hi], color="k", lw=0.8, zorder=1)
    ax1.set_xscale("symlog", linthresh=1e-2); ax1.set_yscale("symlog", linthresh=1e-2)
    ax1.set_xlabel("reference estimate"); ax1.set_ylabel("remlax estimate")
    ax1.set_title("Variance components, correlations, fixed effects, SE, BLUP\n(%d pairs read from the check log)" % len(pts))
    ax1.legend(frameon=False, fontsize=8)

    offset = 0
    for sw in list(SOFTWARE) + ["other"]:
        G = [g for g in gaps if g["software"] == sw]
        if not G:
            continue
        gap = np.abs(np.array([g["remlax"] - g["reference"] for g in G]))
        gap = np.where(gap == 0, 1e-13, gap)
        ax2.scatter(np.arange(len(G)) + offset, gap, s=18, color=COLOUR.get(sw, "k"),
                    label="%s (%d)" % (sw, len(G)), zorder=3)
        offset += len(G)
    ax2.axhline(1e-6, color="k", lw=0.8, ls="--"); ax2.text(0.5, 1.3e-6, "tolerance 1e-6 (lme4, sommer, nlme, closed form)", fontsize=7)
    ax2.axhline(1e-5, color="k", lw=0.8, ls=":"); ax2.text(0.5, 1.3e-5, "tolerance 1e-5 (asreml, own stopping rule)", fontsize=7)
    ax2.set_yscale("log"); ax2.set_ylim(5e-14, 10.0)
    ax2.set_xlabel("check (in table order)"); ax2.set_ylabel("|logLik gap| (remlax vs reference)")
    ax2.set_title("Log-likelihood at the optimum\n(%d checks; 0 drawn at 1e-13)" % len(gaps))
    ax2.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(a.out, dpi=160)
    print("figure:", a.out, "| estimate pairs", len(pts), "| logLik gaps", len(gaps))
    if gaps:
        gg = np.abs(np.array([g["remlax"] - g["reference"] for g in gaps]))
        print("max logLik gap:", gg.max())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
