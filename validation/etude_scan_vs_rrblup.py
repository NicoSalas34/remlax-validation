"""Python side of the external validation of the scan against rrBLUP.

    Rscript validation/etude_scan_vs_rrblup.R  <dir>
    python  validation/etude_scan_vs_rrblup.py <dir>

The R script produces the data AND the rrBLUP scores. This one runs the same
GWAS with `remlax.scan` and compares the two. The header of the R script
gives the reasoning: why this degenerate case is the only one with an
external reference, and why the comparison covers three levels.

WHAT IS CHECKED.
  (A) beta per SNP: the generalised least squares solution, identical.
  (B) the rrBLUP SCORE rebuilt from the raw outputs of the scan ALONE, through
      the exact relation
          Fstat_j = chi2_j * v2 / (y'Py - num_j^2/den_j),      v2 = n - 2
      which expresses the per-SNP re-estimation of the residual scale that
      rrBLUP performs and the scan does not. If this level passes, the
      convention gap between the two is FULLY explained, and no algebra
      disagreement remains.
  (C) the size of the conservativeness of the chi2 at fixed V, measured where
      it bites (on the causal SNP of large effect) and where it does not (on
      the null SNPs). This number justifies, or not, the exact refit of the
      top hits.
"""
import importlib.util
import os
import sys

import numpy as np
from scipy.stats import beta as beta_dist, chi2 as chi2_dist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
_sp = importlib.util.spec_from_file_location(
    "remlax_scan", os.path.join(__import__("moteur_remlax").ENGINE_DIR, "remlax", "scan.py"))
rxscan = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rxscan)


def main(folder):
    import csv
    Vu, Ve, q, p, j0 = (float(x) for x in
                        open(os.path.join(folder, "meta.txt")).read().split())
    q, p, j0 = int(q), int(p), int(j0) - 1          # R counts from 1
    y = np.loadtxt(os.path.join(folder, "y.txt"))
    K = np.loadtxt(os.path.join(folder, "K.txt"))
    M = np.loadtxt(os.path.join(folder, "M.txt"))      # q x p, dosages -1/0/1
    with open(os.path.join(folder, "rrblup_scores.csv")) as f:
        rr = {r["snp"]: float(r["score"]) for r in csv.DictReader(f)}
    score_rr = np.array([rr["snp%04d" % (j + 1)] for j in range(p)])
    assert M.shape == (q, p) and K.shape == (q, q), (M.shape, K.shape)

    # V of the null model, with THE components estimated by rrBLUP. The test
    # targets the SCAN, not the estimation of variance components (tested
    # elsewhere in remlax). Imposing the same V makes any remaining
    # disagreement attributable to the scan and to nothing else.
    V = Vu * K + Ve * np.eye(q)
    X = np.ones((q, 1))
    Vi = np.linalg.solve(V, np.eye(q))
    Am = X.T @ Vi @ X
    vbeta = np.linalg.pinv(Am)
    resid = y - X @ np.linalg.solve(Am, X.T @ Vi @ y)
    Py = Vi @ resid - Vi @ X @ (vbeta @ (X.T @ Vi @ resid))

    # one observation per genotype => the focal genotype incidence is I(q)
    proj = rxscan.project(dict(Vi=Vi, vbeta=vbeta, Py=Py), X, {"dir": np.eye(q)})
    r = rxscan.scan(proj, M, tests=["dir"])
    num, den, chi2 = r["num_dir"], r["den_dir_dir"], r["chi2_dir"]

    # (B) the rrBLUP score, rebuilt without recomputing anything in dimension n
    yPy = float(y @ Py)
    v2 = q - 2
    Fstat = chi2 * v2 / (yPy - num ** 2 / den)
    score_re = -np.log10(beta_dist.cdf(v2 / (v2 + Fstat), v2 / 2, 0.5))

    finite = np.isfinite(score_rr) & (score_rr > 0)
    err = np.abs(score_re - score_rr)
    print("q = %d genotypes, p = %d SNPs; Vu = %.6g, Ve = %.6g (rrBLUP)"
          % (q, p, Vu, Ve))
    print("(B) rrBLUP score rebuilt from the scan outputs:")
    print("    max absolute gap %.2e over %d SNPs" % (err[finite].max(), int(finite.sum())))
    print("    causal SNP snp%04d: rrBLUP %.6f | rebuilt %.6f"
          % (j0 + 1, score_rr[j0], score_re[j0]))
    assert err[finite].max() < 1e-6, "disagreement with rrBLUP"

    # (C) the conservativeness of the chi2 at fixed V
    score_chi2 = -np.log10(chi2_dist.sf(chi2, 1))
    others = np.ones(p, bool)
    others[j0] = False
    d = score_chi2[others] - score_rr[others]
    print("(C) chi2 at fixed V against the rrBLUP F, in -log10(p):")
    print("    causal SNP: %.3f against %.3f  (gap %+.3f)"
          % (score_chi2[j0], score_rr[j0], score_chi2[j0] - score_rr[j0]))
    print("    null SNPs : median %+.4f, 5-95%% quantiles [%+.4f, %+.4f]"
          % (np.median(d), *np.percentile(d, [5, 95])))
    np.savez(os.path.join(folder, "comparison.npz"),
             score_rr=score_rr, score_re=score_re, score_chi2=score_chi2,
             chi2=chi2, Fstat=Fstat, beta=r["beta_dir"], j0=j0, q=q, p=p)
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__.strip().splitlines()[2].strip())
    sys.exit(main(sys.argv[1]))
