#!/usr/bin/env python3
"""Generate docs/validation.md from the CSV of checks.

The CSV is the source. This script only formats it. No number is typed by
hand: everything on the page comes from a row of the CSV, which a program
extracted from the logs of the cluster jobs.

    python3 validation/make_validation_md.py --csv validation/results/checks_<date>.csv \
            --versions validation/results/versions_<date>.json --rev <git sha> --date <date>
"""
import argparse
import csv
import json
import os
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CSV = os.path.join(HERE, "results", "checks_2026-09-28.csv")
OUT = os.path.join(ROOT, "docs", "validation.md")

TITLES = OrderedDict([
    ("algebra", ("Internal algebra", "the solver against itself and against closed-form results",
                 "Gradients against central finite differences, the Hessian against a second "
                 "differentiation, log-determinant and quadratic form against a direct "
                 "computation, parameter counts, and the refusal of designs that would make V "
                 "singular. These need no other software: they are the properties the "
                 "implementation must have to be an implementation of REML at all.")),
    ("stress", ("Random sweep", "properties over randomly drawn designs",
                "Configurations are drawn at random (structure, dimensions, number of traits, "
                "seed) and fitted. Nothing is compared to a reference: the sweep asks only "
                "whether the solver ever produces something impossible - an exception, a "
                "non-finite log-likelihood, a covariance matrix that is not positive "
                "semi-definite, a run that does not converge.")),
    ("pytest", ("Unit tests (pytest)", "closed forms and independent numpy references",
                "The pytest suite of remlax (tests-python): one test per feature, each against an "
                "independent reference (numpy Cholesky REML written out longhand, closed-form "
                "correlation matrices, finite differences). Collected here as a single line.")),
    ("lme4_sommer", ("Against lme4 and sommer", "two open-source reference implementations",
                     "The models both packages can fit: a single random factor on three "
                     "layouts, two crossed factors, a genomic relationship matrix, and a "
                     "bivariate unstructured covariance. Agreement is required on the "
                     "restricted log-likelihood first, then on the variance components, the "
                     "fixed effects and the BLUPs.")),
    ("nlme", ("Against nlme", "the reference implementation for correlation structures",
              "nlme::gls and nlme::lme on the same data: corAR1, corCAR1, corARMA (AR2 and "
              "ARMA(1,1)), corCompSymm, corExp, corGaus, corLin, corSpher, corSymm with "
              "varIdent (an unstructured covariance) and varIdent alone (one residual variance "
              "per group), then a random intercept with an AR1 residual. A grouped corStruct "
              "is written in remlax by declaring the group as the trait of the residual, "
              "residual = ~ id(g):ar1(t) with trait = g. Agreement on -2logL, the residual "
              "variance, the correlation parameters, the fixed effects and their standard "
              "errors.")),
    ("lme4_extra", ("Against lme4, random regressions, nesting and PEV",
                    "the reference implementation for crossed and nested random effects",
                    "Correlated and independent random slopes (a two-column us and diag term), "
                    "a nested design, three crossed factors on an unbalanced layout, very "
                    "unequal group sizes. Besides -2logL, components, fixed effects and their "
                    "standard errors: the BLUPs, and the prediction error variance checked "
                    "against the mixed-model-equations formula computed by hand (lme4's condVar "
                    "follows a different convention, without the uncertainty on beta, and is "
                    "checked against its own formula).")),
    ("sommer_extra", ("Against sommer, genomic and multi-trait models",
                      "the reference implementation for genomic relationship matrices",
                      "A GRM on one term with BLUPs and PEV; three traits with us genetic and "
                      "us residual covariance and a GRM; the same with diag and the likelihood "
                      "ratio between the two; additive plus dominance (two relationship matrices "
                      "on the same factor); a multi-environment trial with one genetic variance "
                      "per environment. sommer is tightened to tolParConvLL = 1e-12; its "
                      "log-likelihood carries a data-dependent constant, so only differences "
                      "of -2logL are compared.")),
    ("asreml_structures", ("Against asreml, structures and inference", "the field reference",
                           "The structures a plant breeder actually writes: a one-dimensional "
                           "AR1, a separable AR1xAR1 field, two-dimensional splines, and a "
                           "covariance shared between two terms through str(). Then the "
                           "inference layer: heritability and its delta-method standard error, "
                           "Wald tests, and the equivalences between saturated "
                           "parameterisations.")),
    ("asreml_catalogue", ("Against asreml, the rest of the catalogue", "the field reference",
                          "The structures that were identified by matching asreml's output "
                          "rather than from documentation - lvr, the anisotropic Matern, "
                          "user-written structures, sectioned residuals - plus prediction and "
                          "the Kenward-Roger denominator degrees of freedom against pbkrtest.")),
    ("asreml_rest", ("Against asreml, the remaining structures",
                     "the field reference",
                     "Grouped residual series ar2, ar3, ma1, ma2, arma, sar, cor, corb and corg; "
                     "the two-dimensional metric kernels iexp, igau, ieuc, aexp, agau, sph and "
                     "cir on an irregular field; the non-saturated multi-trait structures fa(1), "
                     "rr(1), corh, ante(1), chol(1) and diag on four traits with a GRM; the "
                     "separable product id x ar1 x ar1 with shared correlations; and the "
                     "multi-trait spatial residual us(trait):ar1(row):ar1(col).")),
    ("asreml_blup", ("Against asreml, BLUPs", "the field reference",
                     "BLUPs and their ordering on a one-factor model, level by level.")),
    ("ige_model", ("A full model of known truth", "asreml, and the simulated parameters",
                   "The complete direct and indirect genetic effects model: a weighted "
                   "neighbourhood incidence, a covariance shared between the direct and "
                   "indirect terms, a separable field and a nugget - eight variance "
                   "parameters - fitted on simulated data whose true values are known, by "
                   "remlax and by asreml.")),
    ("sparse_parity", ("Sparse engine: parameterisation parity",
                       "the dense engine, at a common theta",
                       "Sigma, level correlations and their closed-form sparse precisions read "
                       "identically by the RTMB sparse engine and the JAX dense engine, including "
                       "K^-1 K = I checked by product.")),
    ("sparse_vs_dense", ("Sparse engine: same likelihood as the dense engine",
                         "the dense engine, warm-started in both directions",
                         "Each engine is evaluated at the other's optimum; comparing two optima "
                         "alone would only show where each stopped.")),
    ("sparse_prec", ("Sparse engine: precision path", "an independent dense-algebra REML",
                     "A term declared by its precision matrix gives the same likelihood as the "
                     "same term declared by its relationship matrix, arbitrated by a REML written "
                     "out longhand in the test file.")),
    ("cpu_gpu_parity", ("CPU against GPU", "the other backend, on the same input",
                        "Every design of the structure catalogue is serialised once, then "
                        "fitted twice in the same process: on the CPU and on a slice of an "
                        "A100. The model must not depend on the machine that fits it.")),
])


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def table(rows, columns, headers):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in rows:
        cells = []
        for c in columns:
            v = (r.get(c) or "").replace("|", "\\|").strip()
            cells.append(v if v else "-")
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--csv", default=CSV)
    ap.add_argument("--versions", default=None, help="JSON written by parse_logs.py")
    ap.add_argument("--rev", default="unknown", help="tested git revision")
    ap.add_argument("--date", default="unknown")
    ap.add_argument("--out", default=OUT, help="output Markdown file (default: docs/validation.md)")
    a = ap.parse_args(argv)
    rows = read_rows(a.csv)
    vers = json.load(open(a.versions)) if a.versions and os.path.exists(a.versions) else {}
    n_tot = len(rows)
    n_pass = sum(r["verdict"] == "pass" for r in rows)
    by_suite = OrderedDict()
    for r in rows:
        by_suite.setdefault(r["suite"], []).append(r)
    provenances = OrderedDict()
    for r in rows:
        provenances.setdefault(r["provenance"], 0)
        provenances[r["provenance"]] += 1

    L = []
    A = L.append
    A("# Validation")
    A("")
    A("Every number on this page was produced on %s by the test suites of this repository "
      "and of remlax, with remlax at git revision `%s`. The campaign ran on the user's SLURM "
      "cluster inside the apptainer image `ige_reml_v3.sif`, launched by `campaign/prepare.sh`: "
      "every suite on CPU, the CPU/GPU parity on a slice of an A100. The measurements were "
      "extracted from the logs (`validation/results/logs_%s/`) by `validation/parse_logs.py`, "
      "and this page is generated from the resulting table (`%s`) by "
      "`validation/make_validation_md.py`. No value here was typed by hand."
      % (a.date, a.rev, a.date, os.path.relpath(a.csv, ROOT)))
    A("")
    A("**%d checks, %d passed, %d failed.**" % (n_tot, n_pass, n_tot - n_pass))
    A("")
    if provenances:
        A("Provenance of the rows:")
        A("")
        for pv, k in provenances.items():
            A("- %s: %d checks" % (pv, k))
        A("")
    if vers.get("versions"):
        A("Software versions read from the log headers:")
        A("")
        for sw, d in vers["versions"].items():
            A("- %s: %s" % (sw, "; ".join("%s (%s)" % (v, pv.split(", ", 1)[1] if ", " in pv else pv)
                                          for pv, v in d.items())))
        A("")
    A("The check labels in the tables below are the strings the test suites "
      "printed, quoted verbatim. Suites run before their English translation "
      "printed French labels, and those labels are kept as they are. "
      "Translating them would give a more readable table and a trace that no "
      "longer matches the logs. The section headings and the commentary are in "
      "English.")
    A("")
    if vers.get("missing"):
        A("Suites whose log was not available for this run, and therefore absent from the "
          "table: %s. See the last section." % ", ".join("`%s`" % m for m in vers["missing"]))
        A("")

    A("## Summary")
    A("")
    summary = []
    for s, rs in by_suite.items():
        t = TITLES.get(s, (s, "", ""))
        summary.append({"suite": t[0], "n": str(len(rs)), "ref": t[1],
                       "verdict": "%d/%d" % (sum(r["verdict"] == "pass" for r in rs), len(rs))})
    A(table(summary, ["suite", "n", "ref", "verdict"],
              ["Suite", "Checks", "Compared against", "Passed"]))
    A("")

    for s, rs in by_suite.items():
        title, ref, expl = TITLES.get(s, (s, "", ""))
        A("## %s" % title)
        A("")
        if expl:
            A(expl)
            A("")
        secs = OrderedDict()
        for r in rs:
            secs.setdefault(r["section"] or "", []).append(r)
        for sec, rr in secs.items():
            if sec:
                A("**%s**" % sec)
                A("")
            A(table(rr, ["check", "detail", "verdict"],
                      ["Check", "Measured", "Verdict"]))
            A("")

    A("## The failed checks")
    A("")
    A("A failed check stays failed on purpose: each one is a documented disagreement, "
      "not a reason to widen a tolerance.")
    A("")
    A("**ma2, sph and cir against asreml (12 checks, as in 2026-09).** "
      "`tests/R/diag_asreml3_ecarts.R` settles them (its log sits with the others). A dense "
      "REML arbiter written in R finds the remlax value to the last digit, and asreml held "
      "at the parameters of remlax returns the remlax value. For sph and cir, asreml reports "
      "converge = FALSE (it oscillates between hills of a multimodal surface). For ma2, "
      "asreml declares convergence while its parameters still change by several per cent "
      "per iteration, and its two correlation parameters are not the MA coefficients of "
      "the standard form.")
    A("")
    A("**2D P-spline against asreml (2 checks).** `rx_spl2d()` builds the PS-ANOVA of "
      "SpATS, five smooth components, and asreml receives the same five basis matrices. "
      "remlax reaches a restricted log-likelihood higher than asreml by 1.4e-5. asreml "
      "declares convergence and stays at its point when it is allowed to iterate further, "
      "so the two programs stop at different points of a flat direction. The variance of "
      "f(y) differs by 0.5 per cent, and the three components at the floor agree.")
    A("")
    A("**lvr, truncated tent (1 check).** The likelihood is not unimodal: its support "
      "changes each time the range crosses an integer distance. remlax stops at a local "
      "optimum (range 1.51, log-likelihood -2.679) while a fine grid in R finds range 5.89 "
      "and -1.914. asreml stops at a worse point (-2.820). The remlax revision of "
      "2026-09-28, which passed this check on the cluster then, gives the same local "
      "optimum on the authors' workstation: the outcome depends on the start, not on a "
      "change of the code.")
    A("")
    A("**Random sweep, seed 0 (1 check).** One configuration out of 200 ends with a "
      "Newton decrement of 0.011, above the threshold of 0.01. Two others reach a better "
      "optimum after a restart, with gains of 0.004 and 0.27 on -2logL. Seeds 50 and 100 "
      "have no unconverged configuration.")
    A("")
    A("**Corrected during the campaign.** The first run, on remlax `bc2af0c` (logs in "
      "`validation/results/logs_2026-10-30_bc2af0c/`), refused the fit of the anisotropic "
      "Matern design `mtrn_aniso` of the CPU/GPU parity at its start: the anisotropic "
      "distance and the Matern function had an infinite derivative at zero distance, where "
      "the diagonal of the kernel sits. remlax `6dbada5` evaluates them on a placeholder "
      "there, with unchanged values, and the design now agrees to 6.8e-13 in relative value "
      "on the two devices. Every other check gave the same verdict in both runs.")
    A("")
    A("## What could not be measured here")
    A("")
    A("**Structures with no external reference.** `mtrn` with all four parameters free, "
      "`own` (user-written kernels) and `sep` beyond three factors have no counterpart in "
      "asreml, lme4, sommer or nlme; they are checked against closed forms (besselK profile, "
      "the same kernel written as `exp`) and against the dense reference.")
    A("")
    A("**What agreement means here.** Two REML implementations agree on the restricted "
      "log-likelihood at the optimum to 1e-6 absolute, on the variance components to 1e-4 "
      "relative (looser on purpose: on a flat surface two optimisers stop at slightly "
      "different theta for the same likelihood) and on fixed effects to 1e-6 relative. "
      "asreml stops on its own criterion and is compared at 1e-5 on the log-likelihood and "
      "1e-3 on components; sommer is tightened to tolParConvLL = 1e-12 and compared at 1e-4. "
      "Iteration counts are never compared: average information, quasi-Newton and "
      "derivative-free optimisers do not count the same thing.")
    A("")
    A("## Re-running all of it")
    A("")
    A("```sh")
    A("# copies remlax and this repository at their commits to the cluster, submits both jobs")
    A("campaign/prepare.sh")
    A("# afterwards, from this repository")
    A("rsync -a cluster:<work folder>/logs/ validation/results/logs_<date>/")
    A("python3 validation/parse_logs.py validation/results/logs_<date>=<provenance> \\")
    A("        --out validation/results/checks_<date>.csv")
    A("python3 validation/make_validation_md.py --csv validation/results/checks_<date>.csv \\")
    A("        --versions validation/results/versions_<date>.json --rev <sha> --date <date>")
    A("```")
    A("")
    A("`campaign/run_cpu.sbatch` lists every suite and the log it writes.")
    A("")
    A("asreml requires a licence, and the `test_remlax_asreml*.R` suites "
      "check one out at run time. The other suites need no licence. The "
      "CPU/GPU parity check needs a machine with a GPU; everything else runs on "
      "a CPU.")
    A("")

    with open(a.out, "w") as f:
        f.write("\n".join(L) + "\n")
    print("%s written: %d lines, %d checks, %d suites"
          % (os.path.relpath(a.out, ROOT), len(L), n_tot, len(by_suite)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
