"""Which covariate detects which mechanism? A study of the three scan families.

THE QUESTION. The model of Sato et al. (2021, 2024) carries, per SNP, the
effect of the focal genotype and the focal x neighbour PRODUCT. The product is
an allelic similarity, hence a DGE x IGE. The indirect genetic effects model of
Muir (2005) and Bijma (2011) carries the effect of the focal genotype and the
additive ALLELE DOSAGE OF THE NEIGHBOUR. These are two distinct biological
hypotheses, and neither model contains the other. If the real mechanism is
additive, does the Sato product term see it? If the mechanism is a
similarity, does the additive term see it? And does the model that carries
all three families lose anything by spending degrees of freedom?

WHAT IS SIMULATED. A randomised block design, a polygenic DGE + IGE background
with a covariance (`us` 2x2 structure) and an AR1 residual within blocks. One
causal SNP whose mechanism varies. The other SNPs are null. The null model of
the scan is THE TRUE model. The study concerns the fixed-effect covariate, not
the specification of the variance components. That second point is settled
elsewhere (see docs/scan.md).

THE EFFECT IS CALIBRATED IN VARIANCE, NOT IN COEFFICIENT. The three covariates
do not have the same spread. The neighbour dosage is a mean over J neighbours,
so it is less spread than the focal dosage, and the product even less.
Comparing at equal coefficient would compare signals of different amplitude.
The standard deviation of the CONTRIBUTION of the causal term is therefore
fixed, identical for the four mechanisms.

Run (from the repository root):
    python validation/etude_scan_familles.py
"""
import importlib.util
import os
import sys
import time

import numpy as np
from scipy.stats import chi2 as chi2_dist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import moteur_remlax  # noqa: E402
sys.path.insert(0, os.path.join(moteur_remlax.ENGINE_DIR, "..", "..", "tests-python"))  # outils_dispositifs, if the remlax repository sits alongside


def _load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


_scan = _load("remlax_scan", os.path.join(moteur_remlax.ENGINE_DIR, "remlax", "scan.py"))
_t = _load("t_scan", os.path.join(ROOT, "tests", "python", "test_scan.py"))
project, scan, lambda_gc = _scan.project, _scan.scan, _scan.lambda_gc

RNG = np.random.default_rng(90210)

# The real mechanisms at the causal SNP.
MECHANISMS = ("direct", "additive neighbour", "similarity", "neighbour + similarity")

# The candidate models, as `scan` specifications.
TESTS = [
    "dir",              # standard GWAS
    "ind",              # additive IGE alone (neighbour dosage)
    "sim",              # similarity alone
    "dir+sim",          # the Sato model: focal + product
    "dir+ind",          # the additive IGE model: focal + neighbour dosage
    "dir+ind+sim",      # the three families, 3 df
    "ind|dir+sim",      # is there an ADDITIVE neighbour effect beyond Sato?
    "sim|dir+ind",      # is there an INTERACTION beyond the additive effects?
]

N_REP = 300
SIGNAL = 0.40       # standard deviation of the causal contribution, in units of y


def contributions(d, j):
    """The three covariates of SNP j, at observation level."""
    u = d["Zdir"] @ d["Mdir"][:, j]
    v = d["Zind"] @ d["Mind"][:, j]
    return {"direct": u, "additive neighbour": v, "similarity": u * v}


def main():
    t0 = time.time()
    d = _t.dispositif(q=90, n_bloc=4, nrow=9, ncol=10, p=300)
    V, _, _ = _t.ajuste(d)
    n, p, X = d["n"], d["p"], d["X"]
    Vi = np.linalg.solve(V, np.eye(n))
    XtVi = X.T @ Vi
    Am = XtVi @ X
    vbeta = np.linalg.pinv(Am)
    Lv = np.linalg.cholesky(V + 1e-9 * np.eye(n))
    b_fixed = np.array([1.0, 1.3, 0.7, -0.4])

    def fit_of(y):
        beta = np.linalg.solve(Am, XtVi @ y)
        resid = y - X @ beta
        Py = Vi @ resid - Vi @ X @ (vbeta @ (XtVi @ resid))
        return dict(Vi=Vi, vbeta=vbeta, Py=Py)

    def run_scan(y):
        proj = project(fit_of(y), X, {"dir": d["Zdir"], "ind": d["Zind"]})
        return scan(proj, {"dir": d["Mdir"], "ind": d["Mind"]}, tests=TESTS)

    # ---- (I) calibration: no causal SNP ----------------------------------
    lam = {t: [] for t in TESTS}
    for _ in range(N_REP):
        y = X @ b_fixed + Lv @ RNG.normal(size=n)
        r = run_scan(y)
        for t in TESTS:
            lam[t].append(lambda_gc(r["chi2_" + t], int(r["df_" + t][0])))
    print("calibration (lambda, median over %d replicates, no causal SNP)" % N_REP)
    for t in TESTS:
        print("    %-14s %.3f" % (t, np.median(lam[t])))

    # ---- (II) power by mechanism ------------------------------------------
    det = {m: {t: [] for t in TESTS} for m in MECHANISMS}
    for mech in MECHANISMS:
        for _ in range(N_REP):
            j0 = int(RNG.integers(p))
            c = contributions(d, j0)
            if mech == "neighbour + similarity":
                # both components with equal weight, each standardised first
                sig = (c["additive neighbour"] / c["additive neighbour"].std()
                       + c["similarity"] / c["similarity"].std())
            else:
                sig = c[mech]
            # CALIBRATION IN VARIANCE: the causal contribution always has
            # standard deviation SIGNAL, whatever the mechanism.
            y = X @ b_fixed + SIGNAL * (sig / sig.std()) + Lv @ RNG.normal(size=n)
            r = run_scan(y)
            for t in TESTS:
                df = int(r["df_" + t][0])
                s = chi2_dist.isf(0.05 / p, df)
                det[mech][t].append(bool(r["chi2_" + t][j0] > s))

    print("\npower: detection of the causal SNP at the Bonferroni threshold "
          "(%d replicates, signal standard deviation %.2f)" % (N_REP, SIGNAL))
    header = "  %-14s" % "tested model" + "".join("%24s" % m for m in MECHANISMS)
    print(header)
    print("  " + "-" * (len(header) - 2))
    tab = {}
    for t in TESTS:
        line = "  %-14s" % t
        for m in MECHANISMS:
            v = 100 * float(np.mean(det[m][t]))
            tab[(t, m)] = v
            line += "%23.0f%%" % v
        print(line)

    np.savez(os.path.join(ROOT, "validation", "results", "scan_families.npz"),
             tests=np.array(TESTS), mechanisms=np.array(MECHANISMS),
             lam=np.array([np.median(lam[t]) for t in TESTS]),
             power=np.array([[tab[(t, m)] for m in MECHANISMS] for t in TESTS]),
             n_rep=N_REP, signal=SIGNAL, n=n, p=p, q=d["q"])
    print("\n%.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
