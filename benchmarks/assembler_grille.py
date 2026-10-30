#!/usr/bin/env python3
"""Assemble the grid CSVs and separate what must be separated.

WHY THIS SCRIPT RATHER THAN A CROSS TABLE. The numbers of iterations differ
between backends on the same cell (419 on the card against 317 on CPU,
measured). A total time therefore mixes TWO effects. One is the unit cost of an
evaluation, which is what the engines are compared on. The other is the length
of the optimisation path, which depends on the conditioning and on the order of
the reductions. Both are reported, together with the time per iteration that
separates them.

It runs on PARTIAL CSVs: the jobs write cell by cell, so a grid in progress
can be assembled and read without waiting.

The CSVs written before the English pass of the scripts (French column names
and values) are read as well: they are mapped to English right after reading.

    python3 benchmarks/assembler_grille.py [csv or glob ...] [--out grid_assembled.csv]
"""
import csv, glob, sys

# Old French column names and values of grille.R, mapped to English.
FR_COLUMNS = {"cellule": "cell", "cas": "case", "moteur": "engine",
              "p_attendu": "p_expected", "fit_paroi_s": "fit_wall_s",
              "fit_interne_s": "fit_internal_s", "tronque": "truncated",
              "statut": "status", "erreur": "error"}
FR_VALUES = {"creux": "sparse", "echec": "fail"}


def to_english(r):
    r = {FR_COLUMNS.get(k, k): v for k, v in r.items()}
    for k in ("case", "engine", "status"):
        if k in r:
            r[k] = FR_VALUES.get(r[k], r[k])
    if "cell" in r and r["cell"]:
        r["cell"] = "-".join(FR_VALUES.get(x, x) for x in r["cell"].split("-"))
    return r


def read(patterns):
    rows = []
    for m in patterns:
        for f in sorted(glob.glob(m)):
            for r in csv.DictReader(open(f)):
                r = to_english(r)
                if r.get("status") == "ok":
                    for k in ("n", "q", "t", "p", "n_iter"):
                        r[k] = int(float(r[k])) if r[k] not in ("", "NA") else None
                    for k in ("fit_wall_s", "fit_internal_s", "logLik"):
                        r[k] = float(r[k]) if r[k] not in ("", "NA") else float("nan")
                    rows.append(r)
    return rows


def main(patterns, out="grid_assembled.csv"):
    rows = read(patterns)
    if not rows:
        print("no rows: the jobs may not have written anything yet"); return 1
    for r in rows:
        # The time per iteration is the quantity COMPARABLE between engines: it
        # removes the path length, which is a property of the surface and of
        # the order of the reductions, not of the engine.
        r["s_per_iter"] = (r["fit_internal_s"] / r["n_iter"]) if r["n_iter"] else float("nan")

    fields = ["cell", "case", "engine", "backend", "n", "q", "t", "p",
              "n_iter", "fit_internal_s", "fit_wall_s", "s_per_iter", "logLik", "tag"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(sorted(rows, key=lambda r: (r["n"], r["q"], r["p"], r["cell"])))

    print("%d cells in %s\n" % (len(rows), out))
    print("=== number of iterations against p, per engine ===")
    by = {}
    for r in rows:
        by.setdefault((r["engine"], r["p"]), []).append(r["n_iter"])
    for k in sorted(by):
        v = by[k]
        print("  engine %-6s p=%-4d : %3d cells, iterations %4d to %4d (median %4d)"
              % (k[0], k[1], len(v), min(v), max(v), sorted(v)[len(v) // 2]))

    print("\n=== cost per iteration, converged cells only ===")
    print("  %-13s %6s %6s %4s %6s %8s %12s" % ("cell", "n", "q", "p", "iter", "total s", "s/iter"))
    for r in sorted(rows, key=lambda r: (r["cell"], r["n"], r["q"], r["p"])):
        if r.get("truncated") == "TRUE":
            continue                      # truncated: not a fit time
        print("  %-13s %6d %6d %4d %6d %8.1f %12.4f"
              % (r["cell"], r["n"], r["q"], r["p"], r["n_iter"],
                 r["fit_internal_s"], r["s_per_iter"]))

    print("\n=== the price of the wrong engine: dense case forced on the sparse engine ===")
    ref = {(r["n"], r["q"], r["t"]): r for r in rows if r["cell"] == "sparse-sparse"}
    for r in sorted(rows, key=lambda r: (r["n"], r["q"], r["p"])):
        if r["cell"] != "dense-sparse":
            continue
        b = ref.get((r["n"], r["q"], r["t"]))
        if b:
            print("  n=%-6d q=%-5d p=%-3d : dense case %9.1f s against sparse case %7.2f s  -> x%.0f"
                  % (r["n"], r["q"], r["p"], r["fit_internal_s"], b["fit_internal_s"],
                     r["fit_internal_s"] / max(b["fit_internal_s"], 1e-9)))

    print("\n=== whole card against CPU, same case and SAME dense engine ===")
    idx = {}
    for r in rows:
        if r["engine"] == "dense":
            idx.setdefault((r["case"], r["n"], r["q"], r["t"]), {})[r["backend"]] = r
    for k in sorted(idx):
        d = idx[k]
        if "cpu" in d and "gpu" in d:
            c, g = d["cpu"], d["gpu"]
            print("  case %-6s n=%-6d q=%-5d p=%-3d | CPU %9.1f s (%3d it) | card %8.1f s (%3d it)"
                  "  -> total x%.1f, per iteration x%.1f"
                  % (k[0], k[1], k[2], k[3], c["fit_internal_s"], c["n_iter"],
                     g["fit_internal_s"], g["n_iter"],
                     c["fit_internal_s"] / max(g["fit_internal_s"], 1e-9),
                     c["s_per_iter"] / max(g["s_per_iter"], 1e-9)))
    return 0


if __name__ == "__main__":
    argv = sys.argv[1:]
    out = "grid_assembled.csv"
    if "--out" in argv:
        i = argv.index("--out")
        out = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]
    m = argv or ["hpc/*/remlax/out/grid*.csv", "hpc/*/grid*.csv",
                 "hpc/*/remlax/out/grille*.csv", "hpc/*/grille*.csv"]
    sys.exit(main(m, out))
