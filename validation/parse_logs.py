#!/usr/bin/env python3
"""Extract the table of checks from the logs of the validation suites.

    python3 validation/parse_logs.py <logs>[=<provenance>] [<logs2>[=<provenance2>] ...] \
            --out validation/results/checks_<date>.csv

WHY A PARSER RATHER THAN TYPING. The test suites print one line per check,
with a verdict. Copying these lines into a table by hand invites a typo in a
number that a paper will later cite. The parser cannot get a number wrong. It
takes the line as it is.

WHAT IT READS. The formats the suites produce:
  - a verdict line: two spaces, a label, at least two spaces, then OK or FAIL
    (ECHEC in older logs), then a numeric detail (main R suites);
  - the variant `  [OK] label   detail` of the sparse engine suites;
  - the variant `  OK    label  detail` of test_remlax_asreml_blup.R;
  - a section title:   === 3. 2D spline ===
  - for the random sweep, the counters of its summary;
  - for pytest, the final line `N passed, M skipped`;
  - for the CPU/GPU parity, the table of parite_gpu.py.

SEVERAL DIRECTORIES. The same log may exist in several copies (cluster and
local). The directories are read in the order given and the FIRST copy found
of a log is authoritative, so give the cluster first. Each directory carries
its provenance after `=`.

Labels are copied verbatim from the logs. Logs written before the English pass
of the suites (such as the 2026-09 campaign) carry French labels, and the
table keeps them. Translating them here would give a more readable table and a
false trace.
"""
import argparse
import csv
import json
import os
import re
from collections import Counter, OrderedDict

# (file, suite, reference): this order sets the order of the final table.
SOURCES = [
    ("01_core.log", "algebra", "internal algebra"),
    ("00_pytest.log", "pytest", "internal algebra and closed forms (pytest)"),
    ("02_stress_s0.log", "stress", "properties (no reference implementation)"),
    ("02_stress_s50.log", "stress", "properties (no reference implementation)"),
    ("02_stress_s100.log", "stress", "properties (no reference implementation)"),
    ("03_lme4_sommer.log", "lme4_sommer", "lme4 / sommer"),
    ("07_nlme.log", "nlme", "nlme"),
    ("08_lme4_extra.log", "lme4_extra", "lme4"),
    ("09_sommer_extra.log", "sommer_extra", "sommer"),
    ("04_asreml.log", "asreml_structures", "asreml"),
    ("05_asreml2.log", "asreml_catalogue", "asreml / pbkrtest / R"),
    ("15_asreml3.log", "asreml_rest", "asreml"),
    ("06b_asreml_blup.log", "asreml_blup", "asreml"),
    ("06_ige.log", "ige_model", "asreml + known truth"),
    # The three tmb logs exist only in the 2026-09 campaign. The RTMB engine
    # was removed from remlax afterwards, together with its suites.
    ("10_tmb_parity.log", "sparse_parity", "dense engine (same theta)"),
    ("11_tmb_vs_dense.log", "sparse_vs_dense", "dense engine (warm start both ways)"),
    ("12_tmb_prec.log", "sparse_prec", "independent dense-algebra REML"),
    ("21_parite.log", "cpu_gpu_parity", "the other backend (same input bundle)"),
]

# FAIL is the verdict of the English suites, ECHEC that of older French logs.
VERDICT = re.compile(r"^  (?P<name>.{2,70}?)\s{2,}(?P<verdict>OK|FAIL|ECHEC)\s*(?P<detail>.*)$")
VERDICT_BRACKETS = re.compile(r"^\s*\[(?P<verdict>OK|FAIL|ECHEC)\]\s+(?P<name>.+?)(?:\s{2,}(?P<detail>.*))?$")
VERDICT_HEAD = re.compile(r"^  (?P<verdict>OK|FAIL|ECHEC)\s+(?P<name>.+?)(?:\s{2,}(?P<detail>.*))?$")
SECTION = re.compile(r"^=+\s*(?P<t>[A-Z]?\d*[a-z]?\..*?)\s*=+\s*$")
PARITY = re.compile(r"^(?P<m>[a-z0-9_.]+)\s+(?P<c>-?\d+\.\d+)\s+(?P<g>-?\d+\.\d+)\s+"
                    r"(?P<rel>[\d.]+e[+-]\d+)\s+(?P<dth>[\d.]+e[+-]\d+)\s+"
                    r"(?P<sc>[\d.]+)\s+(?P<sg>[\d.]+)", re.M)
PYTEST = re.compile(r"^(?:=+ )?(?P<p>\d+) passed(?:, (?P<s>\d+) skipped)?(?:, (?P<f>\d+) failed)?"
                    r"(?:, (?P<e>\d+) errors?)? in (?P<t>[\d.]+)s", re.M)
# Counters of the stress summary, each with its English and older French label.
COUNTERS = OrderedDict([
    ("exception", ["exception"]),
    ("non-finite logLik", ["non-finite logLik", "logLik non finie"]),
    ("not converged", ["not converged", "non convergé"]),
    ("Sigma not PSD", ["Sigma not PSD", "Sigma non PSD"]),
    ("unstable optimum", ["unstable optimum", "optimum instable"]),
])
# Software versions, read from the headers the suites print. They are read in
# the first EIGHT lines of each log (the suites print their versions first) and
# in the whole of 00_versions.log, a file the job writes with the line
# [R] ... | asreml ... and the line [jax] ... .
VERSIONS = [
    ("R", re.compile(r"R version (\d+\.\d+\.\d+)")),
    ("R", re.compile(r"\| R (\d+) (\d+\.\d+)")),
    ("asreml", re.compile(r"asreml (\d+\.\d+\.\d+[\d.]*)")),
    ("lme4", re.compile(r"lme4 (\d+\.\d+[\d.-]*)")),
    ("sommer", re.compile(r"sommer (\d+\.\d+\.\d+)")),
    ("nlme", re.compile(r"nlme (\d+\.\d+[-.]\d+)")),
    ("pbkrtest", re.compile(r"pbkrtest (\d+\.\d+\.\d+)")),
    ("jax", re.compile(r"\[jax\] (\d+\.\d+\.\d+)")),
    # Only used by the 2026-09 logs (the RTMB engine has since left remlax).
    ("RTMB", re.compile(r"RTMB (\d+\.\d+[\d.]*)")),
]
N_VERSION_LINES = 8


def find_log(root, name):
    """The same log may appear at several depths, depending on the job."""
    for d, _, fs in os.walk(root):
        if name in fs:
            return os.path.join(d, name)
    return None


def verdict_rows(txt, suite, ref, prov, src):
    section, out = "", []
    for line in txt.splitlines():
        m = SECTION.match(line.strip())
        if m:
            section = m.group("t").strip()
            continue
        m = VERDICT.match(line.rstrip()) or VERDICT_BRACKETS.match(line.rstrip()) \
            or VERDICT_HEAD.match(line.rstrip())
        if m:
            out.append(dict(suite=suite, section=section, check=m.group("name").strip(),
                            reference=ref,
                            verdict="pass" if m.group("verdict") == "OK" else "FAIL",
                            detail=(m.group("detail") or "").strip(),
                            provenance=prov, source_log=src))
    return out


def stress_rows(txt, suite, ref, prov, src, seed):
    cnt = {}
    for k, labels in COUNTERS.items():
        alt = "|".join(re.escape(x) for x in labels)
        m = re.search(r"^  (?:%s)\s+(\d+)\s*$" % alt, txt, re.M)
        cnt[k] = int(m.group(1)) if m else None
    m = re.search(r"^(\d+)/(\d+) (?:without any problem|sans aucun probleme)", txt, re.M)
    ok, tot = (int(m.group(1)), int(m.group(2))) if m else (0, 0)
    trip = re.findall(r"^\s+\((\d+), ([\d.]+), ([\d.]+)\)\s*$", txt, re.M)
    impossible = ["exception", "non-finite logLik", "not converged", "Sigma not PSD"]
    detail = "clean %d/%d; unstable optima %d" % (ok, tot, cnt.get("unstable optimum") or 0)
    if trip:
        detail += " (best-of-restart gains %s in -2logL)" % ", ".join(
            "%.4f" % (float(a) - float(b)) for _, a, b in trip)
    return [dict(suite=suite, section="random configuration sweep, seed %d" % seed,
                 check="%d configurations: no exception, finite logLik, PSD Sigma, converged" % tot,
                 reference=ref,
                 verdict="pass" if all((cnt.get(k) or 0) == 0 for k in impossible) else "FAIL",
                 detail=detail, provenance=prov, source_log=src)]


def pytest_rows(txt, suite, ref, prov, src):
    m = None
    for m in PYTEST.finditer(txt):
        pass
    if m is None:
        return []
    p, s, f, e = (int(m.group(k) or 0) for k in ("p", "s", "f", "e"))
    return [dict(suite=suite, section="pytest, tests/python",
                 check="%d tests collected by pytest" % (p + s + f + e), reference=ref,
                 verdict="pass" if f == 0 and e == 0 else "FAIL",
                 detail="%d passed, %d skipped, %d failed, %d errors in %s s" % (p, s, f, e, m.group("t")),
                 provenance=prov, source_log=src)]


def parity_rows(txt, suite, ref, prov, src, tol=1e-8):
    out = []
    for d in (m.groupdict() for m in PARITY.finditer(txt)):
        out.append(dict(
            suite=suite, section="structure catalogue, 34 serialised designs",
            check="design '%s': same -2logL on CPU and on an A100 MIG slice" % d["m"],
            reference=ref, verdict="pass" if float(d["rel"]) <= tol else "FAIL",
            detail="-2logL %s vs %s | rel gap %s | max|dtheta| %s | %s s CPU vs %s s MIG"
                   % (d["c"], d["g"], d["rel"], d["dth"], d["sc"], d["sg"]),
            provenance=prov, source_log=src))
    return out


def versions_in(txt):
    v = {}
    for name, rx in VERSIONS:
        m = rx.search(txt)
        if m and name not in v:
            v[name] = ".".join(g for g in m.groups() if g)
    return v


COLS = ["suite", "section", "check", "reference", "verdict", "detail",
        "provenance", "source_log"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("logs", nargs="+",
                    help="log directory or directories, each written path[=provenance]")
    ap.add_argument("--out", required=True)
    ap.add_argument("--provenance", default="re-measured, provenance not given")
    ap.add_argument("--versions-out", default=None,
                    help="JSON of the software versions read in the logs (default: next to --out)")
    a = ap.parse_args(argv)

    roots = []
    for spec in a.logs:
        path, _, prov = spec.partition("=")
        roots.append((path, prov or a.provenance))

    rows, missing, versions, origin = [], [], OrderedDict(), {}
    for name, suite, ref in SOURCES:
        p, prov = None, None
        for root, prov_r in roots:
            p = find_log(root, name)
            if p is not None:
                prov = prov_r
                break
        if p is None:
            missing.append(name)
            continue
        origin[name] = p
        txt = open(p, errors="replace").read()
        for k, v in versions_in("\n".join(txt.splitlines()[:N_VERSION_LINES])).items():
            versions.setdefault(k, OrderedDict()).setdefault(prov, v)
        if suite == "stress":
            seed = int(re.search(r"_s(\d+)\.log$", name).group(1))
            rows += stress_rows(txt, suite, ref, prov, name, seed)
        elif suite == "cpu_gpu_parity":
            rows += parity_rows(txt, suite, ref, prov, name)
        elif suite == "pytest":
            rows += pytest_rows(txt, suite, ref, prov, name)
        else:
            rows += verdict_rows(txt, suite, ref, prov, name)

    for root, prov_r in roots:
        pv = find_log(root, "00_versions.log")
        if pv:
            for k, v in versions_in(open(pv, errors="replace").read()).items():
                versions.setdefault(k, OrderedDict())[prov_r] = v
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLS})
    vout = a.versions_out or re.sub(r"checks_", "versions_", a.out).replace(".csv", ".json")
    with open(vout, "w") as f:
        json.dump({"versions": versions, "logs": origin, "missing": missing}, f, indent=1)

    print("%d checks written to %s" % (len(rows), a.out))
    print("by suite :", dict(Counter(r["suite"] for r in rows)))
    print("verdicts :", dict(Counter(r["verdict"] for r in rows)))
    print("versions :", json.dumps(versions))
    if missing:
        print("MISSING logs (suites not covered):", ", ".join(missing))
    failed = [r["check"] for r in rows if r["verdict"] == "FAIL"]
    if failed:
        print("FAILURES:", failed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
