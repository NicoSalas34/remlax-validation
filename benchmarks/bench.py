#!/usr/bin/env python3
"""remlax benchmarks: time, compilation and memory, by backend.

    python3 benchmarks/bench.py --backend cpu --suite all --out results/cpu.csv

NO dependency on the host: this script runs as is on a compute node, in a
container, without network. It writes one CSV and one line per measurement.

WHAT IS MEASURED, AND WHY THIS SPLIT
------------------------------------
V is formed DENSE (n x n, float64) and factorised. The cost is O(n^3) per
evaluation and the footprint O(n^2). The suites separate the axes that decide
between CPU and GPU.

  scaling      n grows, model fixed. The dominant axis.
  params       number of variance parameters grows, n fixed. The number of
               L-BFGS-B iterations depends on it, the cost of one evaluation
               does not.
  structures   one catalogue structure per row, fixed size. Isolates the cost
               of BUILDING V (kernels, Bessel quadrature) from the cost of its
               factorisation.
  compile      first call of the objective (XLA compilation included) against
               the later calls. On GPU, compilation can dominate the compute
               for small models. This measurement tells a user whether the GPU
               is worth it.

PROTOCOL. Seed fixed per case. One warm-up outside the measurement, then --reps
repetitions, reported as median and range. The time reported by `fit_*` is the
wall time of a COMPLETE fit (initialisation, L-BFGS-B, Newton polishing,
Hessian if requested), which is what the user pays.
"""
import argparse
import csv
import gc
import json
import os
import platform
import resource
import statistics
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import moteur_remlax  # noqa: E402,F401  (engine of the installed R package)


# ==============================================================================
# Designs
# ==============================================================================
def design_iid(n, q=None, seed=0):
    """One random factor + iid residual: the simplest model that exercises the
    whole path (assembly of V, Cholesky, gradient, polishing)."""
    q = q or max(5, n // 5)
    rng = np.random.default_rng(seed)
    lev = rng.integers(0, q, n)
    u = rng.normal(size=q)
    y = 2.0 + u[lev] + rng.normal(size=n, scale=0.8)
    term = dict(name="g", struct="iid", t=1, rank=0, q=q,
                zi=np.arange(n), zj=lev.astype(np.int64), zx=np.ones(n), LK=None)
    res = dict(struct="iid", t=1, rank=0,
               trait=np.zeros(n, np.int64), unit=np.arange(n))
    return [term], res, y, np.ones((n, 1))


def design_traits(n_unit, t, seed=0, struct="us"):
    """t traits in long format: n = n_unit * t observations, free t x t Sigma
    on the genetic term and a diagonal residual between traits."""
    q = max(6, n_unit // 4)
    rng = np.random.default_rng(seed)
    n = n_unit * t
    unit = np.tile(np.arange(n_unit), t)
    trait = np.repeat(np.arange(t), n_unit)
    lev = rng.integers(0, q, n_unit)
    lev_long = np.tile(lev, t)
    zi = np.arange(n)
    zj = trait * q + lev_long                     # LEVEL varies fastest
    U = rng.normal(size=(q, t))
    y = 1.0 + U[lev_long, trait] + rng.normal(size=n, scale=0.7)
    term = dict(name="g", struct=struct, t=t, rank=(2 if struct == "fa" else 0),
                q=q, zi=zi, zj=zj.astype(np.int64), zx=np.ones(n), LK=None)
    res = dict(struct="diag", t=t, rank=0, trait=trait.astype(np.int64), unit=unit)
    X = np.zeros((n, t))
    X[np.arange(n), trait] = 1.0                  # one mean per trait
    return [term], res, y, X


def design_genomic(n, q, t=1, seed=0, n_marq=None, ige=False):
    """Dense GENOMIC kinship: the case that really dominates compute time.

    WHY THIS SUITE EXISTS. With K = I, term_factor takes the `LK is None`
    branch and forms NO product: (L_Sigma (x) I) acts block by block. A sweep
    at K = I therefore measures the cheapest possible assembly of V, and
    underestimates the cost of the case a breeder cares about. With a dense K,
    one evaluation also pays

        ZL = Z (I (x) L_K)      n t q^2
        B  = ZL (L_Sigma (x) I) n t m q
        V += B B'               n^2 m q

    At n = q = 2000 this is about 16 Gflop against 2.7 Gflop for the Cholesky,
    so assembly dominates. These are dense matrix products, where a GPU is most
    favoured. This axis is therefore measured separately.

    ige=True gives the shape of the model with direct and indirect genetic
    effects: ONE term with t = 2 and a free 2x2 covariance, the first trait
    column carrying the direct incidence and the second a WEIGHTED neighbour
    incidence. A covariance between two genetic terms is expressed this way
    here, and it multiplies the cost of the product above by t^2.
    """
    rng = np.random.default_rng(seed)
    n_marq = n_marq or max(200, q // 2)
    # K = W W' / m on biallelic markers, centred and scaled, plus a tiny ridge.
    # This is the VanRaden construction, the one a breeder actually supplies,
    # and it is FULL.
    p = rng.uniform(0.05, 0.95, n_marq)
    W = rng.binomial(2, p, size=(q, n_marq)).astype(np.float64)
    W -= 2.0 * p
    W /= np.sqrt(2.0 * np.sum(p * (1.0 - p)))
    K = W @ W.T
    K[np.diag_indices(q)] += 1e-6 * np.trace(K) / q
    LK = np.linalg.cholesky(K)

    lev = rng.integers(0, q, n)                      # replicates per genotype
    struct = "us" if t > 1 else "iid"
    if ige and t == 2:
        # trait 0: direct effect of the plant genotype
        # trait 1: indirect effect, weighted sum of the neighbours
        zi = np.concatenate([np.arange(n), np.arange(n)])
        neigh = rng.integers(0, q, n)
        zj = np.concatenate([lev, q + neigh]).astype(np.int64)
        w = rng.uniform(0.3, 1.0, n)                 # neighbour weights
        zx = np.concatenate([np.ones(n), w])
    else:
        zi = np.concatenate([np.arange(n)] * t)
        zj = np.concatenate([k * q + lev for k in range(t)]).astype(np.int64)
        zx = np.ones(n * t)
    U = rng.normal(size=(q, t))
    y = 1.0 + U[lev, 0] + rng.normal(size=n, scale=0.8)
    term = dict(name="g", struct=struct, t=t, rank=0, q=q,
                zi=zi, zj=zj, zx=zx, LK=LK)
    res = dict(struct="iid", t=1, rank=0,
               trait=np.zeros(n, np.int64), unit=np.arange(n))
    return [term], res, y, np.ones((n, 1))


def design_level(n_row, n_col, kind, seed=0):
    """Structured field between LEVELS on an n_row x n_col grid: ar1, separable
    product, or metric kernel on coordinates. One level per cell."""
    rng = np.random.default_rng(seed)
    q = n_row * n_col
    n = q
    rr, cc = np.meshgrid(np.arange(n_row), np.arange(n_col), indexing="ij")
    coord = np.column_stack([rr.ravel().astype(float), cc.ravel().astype(float)])
    y = rng.normal(size=n)
    term = dict(name="f", struct="iid", t=1, rank=0, q=q,
                zi=np.arange(n), zj=np.arange(n), zx=np.ones(n), LK=None)
    term["lvl"] = kind
    if kind == "ar1ar1":
        term["dims"] = (n_row, n_col)
    elif kind == "ar1":
        pass                                       # 1D on the order of the levels
    if kind in ("iexp", "igau", "ieuc", "sph", "cir", "aexp", "agau", "mtrn", "lvr"):
        term["coord"] = coord
    if kind == "exp" or kind == "gau":
        term["coord"] = coord[:, :1]
    if kind == "mtrn":
        term["lvl_opts"] = {"est_phi": 1.0, "est_nu": 1.0}
    if kind == "mtrn_aniso":
        term["lvl"] = "mtrn"
        term["coord"] = coord
        term["lvl_opts"] = {"est_phi": 1.0, "est_nu": 1.0,
                            "est_delta": 1.0, "est_alpha": 1.0}
    res = dict(struct="iid", t=1, rank=0,
               trait=np.zeros(n, np.int64), unit=np.arange(n))
    return [term], res, y, np.ones((n, 1))


# ==============================================================================
# Measurement
# ==============================================================================
def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def gpu_mem_mb():
    """Peak device memory use, when JAX exposes it."""
    try:
        import jax
        st = jax.local_devices()[0].memory_stats() or {}
        peak = st.get("peak_bytes_in_use") or st.get("bytes_in_use")
        return (peak / 1e6) if peak else float("nan")
    except Exception:
        return float("nan")


def timed_fit(terms, res, y, X, reps, hessian, **kw):
    """Repeated fits. Returns the median, min, max and the last result."""
    from remlax.fit import fit_reml
    ts = []
    out = None
    for _ in range(reps):
        gc.collect()
        t0 = time.perf_counter()
        out = fit_reml(terms, res, y, X, verbose=False, hessian=hessian,
                       blups=False, **kw)
        ts.append(time.perf_counter() - t0)
    return statistics.median(ts), min(ts), max(ts), out


def timed_objective(terms, res, y, X, reps):
    """Separates compilation from execution.

    make_objective returns a jitted function. Its FIRST call pays the trace and
    the XLA compilation, the later calls do not. The difference is exactly what
    makes a small model faster on CPU than on GPU.
    """
    import jax
    from remlax.model import make_objective
    from remlax.start import initial_theta
    fun_sc, fun, sc, Zs = make_objective(terms, res, y, X)
    th = np.asarray(initial_theta(terms, res, y))
    t0 = time.perf_counter()
    v, g = fun(th)
    jax.block_until_ready(v)
    t_first = time.perf_counter() - t0
    ts = []
    for k in range(reps):
        thk = th + 1e-3 * (k + 1)
        t0 = time.perf_counter()
        v, g = fun(thk)
        jax.block_until_ready(v)
        ts.append(time.perf_counter() - t0)
    return t_first, statistics.median(ts), float(v)


# ==============================================================================
# Suites
# ==============================================================================
def suite_scaling(rows, args, meta):
    ns = [int(v) for v in args.ns.split(",")] if args.ns else \
         [250, 500, 1000, 2000, 3000, 4000, 6000, 8000, 12000, 16000]
    ns = [n for n in ns if n <= args.nmax]
    for n in ns:
        terms, res, y, X = design_iid(n, seed=1)
        # A complete fit costs O(n^3) per evaluation and several tens of
        # evaluations. Repeating beyond n = 4000 buys a range at the price of
        # hours of compute, so repetitions are kept where they are informative.
        reps = args.reps if n <= args.reps_upto else 1
        try:
            med, lo, hi, out = timed_fit(terms, res, y, X, reps, hessian=False)
        except Exception as e:                      # OOM: record the limit
            rows.append(dict(meta, suite="scaling", case="iid", n=n, n_par=2,
                             status="failed", error=type(e).__name__ + ": " + str(e)[:200]))
            print("[scaling] n=%d FAIL %s" % (n, type(e).__name__), flush=True)
            break
        rows.append(dict(meta, suite="scaling", case="iid", n=n, q=terms[0]["q"],
                         n_par=int(out["n_par"]), status="ok",
                         fit_s_median=med, fit_s_min=lo, fit_s_max=hi,
                         # WITHOUT n_iter, a fit time difference between two
                         # backends cannot be used: it may come from the cost
                         # of one evaluation or from their number.
                         n_iter=int(out.get("n_iter") or -1),
                         logLik=float(out["logLik"]),
                         newton_decrement=float(out.get("newton_decrement") or float("nan")),
                         rss_mb=rss_mb(), gpu_mb=gpu_mem_mb()))
        print("[scaling] n=%-6d %8.2f s  logLik %.6f" % (n, med, out["logLik"]), flush=True)


def suite_compile(rows, args, meta):
    ns = [int(v) for v in args.ns.split(",")] if args.ns else \
         [250, 500, 1000, 2000, 4000, 8000, 12000, 16000]
    ns = [n for n in ns if n <= args.nmax]
    for n in ns:
        terms, res, y, X = design_iid(n, seed=2)
        try:
            t_first, t_run, v = timed_objective(terms, res, y, X, args.reps)
        except Exception as e:
            rows.append(dict(meta, suite="compile", case="iid", n=n,
                             status="failed", error=type(e).__name__ + ": " + str(e)[:200]))
            break
        rows.append(dict(meta, suite="compile", case="iid", n=n, status="ok",
                         first_call_s=t_first, run_call_s=t_run,
                         compile_s=max(t_first - t_run, 0.0),
                         compile_share=(t_first - t_run) / t_first if t_first > 0 else float("nan"),
                         neg2logl=v, rss_mb=rss_mb(), gpu_mb=gpu_mem_mb()))
        print("[compile] n=%-6d first %7.3f s | execution %7.4f s | compilation share %5.1f %%"
              % (n, t_first, t_run, 100 * (t_first - t_run) / max(t_first, 1e-12)), flush=True)


def suite_params(rows, args, meta):
    n_unit = args.n_unit
    for t in [int(v) for v in (args.ts or "1,2,3,4,5,6,8").split(",")]:
        for struct in ("us", "fa") if t >= 3 else ("us",):
            terms, res, y, X = design_traits(n_unit, t, seed=3, struct=struct)
            n = len(y)
            if n > args.nmax:
                continue
            try:
                med, lo, hi, out = timed_fit(terms, res, y, X, args.reps, hessian=True)
            except Exception as e:
                rows.append(dict(meta, suite="params", case="%s(t=%d)" % (struct, t),
                                 n=n, status="failed",
                                 error=type(e).__name__ + ": " + str(e)[:200]))
                continue
            rows.append(dict(meta, suite="params", case="%s(t=%d)" % (struct, t),
                             n=n, n_traits=t, n_par=int(out["n_par"]), status="ok",
                             fit_s_median=med, fit_s_min=lo, fit_s_max=hi,
                             n_iter=int(out.get("n_iter") or -1),
                             logLik=float(out["logLik"]),
                             rss_mb=rss_mb(), gpu_mb=gpu_mem_mb()))
            print("[params] %-8s t=%d n=%-6d p=%-3d %8.2f s"
                  % (struct, t, n, out["n_par"], med), flush=True)


def suite_structures(rows, args, meta):
    nr, nc = args.grid_rows, args.grid_cols
    cases = [("iid", None), ("us(3)", None), ("fa(4,2)", None),
           ("ar1", "ar1"), ("ar1ar1", "ar1ar1"),
           ("iexp", "iexp"), ("igau", "igau"), ("sph", "sph"),
           ("mtrn", "mtrn"), ("mtrn_aniso", "mtrn_aniso")]
    for name, kind in cases:
        if kind is None:
            if name == "iid":
                terms, res, y, X = design_iid(nr * nc, seed=4)
            elif name == "us(3)":
                terms, res, y, X = design_traits((nr * nc) // 3, 3, seed=4, struct="us")
            else:
                terms, res, y, X = design_traits((nr * nc) // 4, 4, seed=4, struct="fa")
        else:
            terms, res, y, X = design_level(nr, nc, kind, seed=4)
        n = len(y)
        try:
            t_first, t_run, v = timed_objective(terms, res, y, X, args.reps)
            if args.no_fit:
                med, out = float("nan"), {"n_par": -1, "logLik": float("nan")}
            else:
                med, lo, hi, out = timed_fit(terms, res, y, X, 1, hessian=False)
            ok, err = "ok", ""
        except Exception as e:
            t_first = t_run = v = med = float("nan")
            out = {"n_par": -1, "logLik": float("nan")}
            ok, err = "failed", type(e).__name__ + ": " + str(e)[:200]
        rows.append(dict(meta, suite="structures", case=name, n=n, status=ok, error=err,
                         first_call_s=t_first, run_call_s=t_run,
                         compile_s=(t_first - t_run) if t_first == t_first else float("nan"),
                         fit_s_median=med, n_par=int(out["n_par"]),
                         logLik=float(out["logLik"]),
                         rss_mb=rss_mb(), gpu_mb=gpu_mem_mb()))
        print("[struct] %-12s n=%-5d compilation %8.3f s | execution %8.4f s | fit %8.2f s %s"
              % (name, n, t_first - t_run if t_first == t_first else float("nan"),
                 t_run, med, err), flush=True)


def suite_genomic(rows, args, meta):
    """Dense kinship: q grows at fixed n, then n grows at fixed q/n."""
    qs = [int(v) for v in (args.qs or "250,500,1000,1500,2000,3000,4000").split(",")]
    for q in qs:
        n = args.n_gen
        if q > n:
            continue
        terms, res, y, X = design_genomic(n, q, t=1, seed=7)
        try:
            t_first, t_run, v = timed_objective(terms, res, y, X, args.reps)
            if args.no_fit:
                # On dense kinship a complete fit needs HUNDREDS of
                # evaluations: 435 s measured for 0.247 s per evaluation at
                # q = 2000. The evaluation is therefore timed, since it carries
                # the comparison between backends. The number of iterations is
                # measured by the suites that can measure it.
                med, out = float('nan'), {'n_par': -1, 'logLik': float('nan')}
            else:
                med, lo, hi, out = timed_fit(terms, res, y, X, 1, hessian=False)
            ok, err, ni = "ok", "", int(out.get("n_iter") or -1)
        except Exception as e:
            t_first = t_run = med = float("nan")
            out = {"n_par": -1, "logLik": float("nan")}
            ok, err, ni = "failed", type(e).__name__ + ": " + str(e)[:200], -1
        # Expected assembly Gflop per evaluation, to relate time to the amount
        # of arithmetic and not to size alone.
        gflop = (n * q * q + n * n * q) / 1e9
        rows.append(dict(meta, suite="genomic", case="iid+K(q=%d)" % q, n=n, q=q,
                         status=ok, error=err, first_call_s=t_first, run_call_s=t_run,
                         compile_s=(t_first - t_run) if t_first == t_first else float("nan"),
                         fit_s_median=med, n_iter=ni, n_par=int(out["n_par"]),
                         logLik=float(out["logLik"]), assembly_gflop=gflop,
                         rss_mb=rss_mb(), gpu_mb=gpu_mem_mb()))
        print("[genomic] q=%-5d n=%-5d assembly %7.1f Gflop | evaluation %8.4f s "
              "| fit %8.2f s %s" % (q, n, gflop, t_run, med, err), flush=True)


def suite_crossterm(rows, args, meta):
    """Covariance between genetic terms, with dense kinship.

    t = 1: a single genetic term
    t = 2: two correlated terms, shape of the direct/indirect model (weighted
           neighbour incidence on the second)
    t = 3: three correlated terms, free 3x3 covariance
    """
    n, q = args.n_gen, args.q_gen
    for t in [int(v) for v in (args.ts or "1,2,3").split(",")]:
        for ige in ([False, True] if t == 2 else [False]):
            terms, res, y, X = design_genomic(n, q, t=t, seed=8, ige=ige)
            name = "us(%d)+K%s" % (t, " [DGE/IGE]" if ige else "")
            try:
                t_first, t_run, v = timed_objective(terms, res, y, X, args.reps)
                if args.no_fit:
                    # On dense kinship a complete fit needs HUNDREDS of
                    # evaluations: 435 s measured for 0.247 s per evaluation at
                    # q = 2000. The evaluation is therefore timed, since it
                    # carries the comparison between backends. The number of
                    # iterations is measured by the suites that can measure it.
                    med, out = float('nan'), {'n_par': -1, 'logLik': float('nan')}
                else:
                    med, lo, hi, out = timed_fit(terms, res, y, X, 1, hessian=True)
                ok, err, ni = "ok", "", int(out.get("n_iter") or -1)
            except Exception as e:
                t_first = t_run = med = float("nan")
                out = {"n_par": -1, "logLik": float("nan")}
                ok, err, ni = "failed", type(e).__name__ + ": " + str(e)[:200], -1
            gflop = (n * t * q * q + n * n * t * q) / 1e9
            rows.append(dict(meta, suite="crossterm", case=name, n=n, q=q, n_traits=t,
                             status=ok, error=err, first_call_s=t_first, run_call_s=t_run,
                             compile_s=(t_first - t_run) if t_first == t_first else float("nan"),
                             fit_s_median=med, n_iter=ni, n_par=int(out["n_par"]),
                             logLik=float(out["logLik"]), assembly_gflop=gflop,
                             rss_mb=rss_mb(), gpu_mb=gpu_mem_mb()))
            print("[crossterm] %-18s p=%-3d assembly %7.1f Gflop | evaluation %8.4f s "
                  "| fit %8.2f s %s"
                  % (name, out["n_par"], gflop, t_run, med, err), flush=True)


# ==============================================================================
# Where the density lives: in G or in R
# ==============================================================================
def design_gr(n_unit, t, seed=0, g_dense=False, r_dense=False, q_gen=None,
              n_marq=400):
    """The same model, the same number of parameters, the density MOVED.

    THE QUESTION. The mixed model equations contain Z' R^-1 Z + G^-1. A density
    in G fills only its own block, from which it spreads through elimination. A
    density in R fills Z' R^-1 Z for ALL terms at once, including those with a
    diagonal precision. The prediction is that the LOCATION of the density
    matters more than its amount, and that a us residual rules out the sparse
    route when a dense genomic kinship does not.

    WHY TWO PATHS ARE NEEDED TO RAISE THE PARAMETER COUNT. The number of
    parameters is not a free setting: a us on t traits imposes t(t+1)/2.
    Raising p THROUGH THE RESIDUAL confounds p with the structure of R, and any
    conclusion on p would be a conclusion on R. Both paths are offered at equal
    p:

        g_dense=True, r_dense=False   us on the GENETIC term, diagonal R
        g_dense=False, r_dense=True   diagonal G, us on the RESIDUAL

    At equal t, these two cases have the SAME number of covariance parameters
    and the same problem size. They differ only by the side where the density
    sits, which is exactly the comparison sought.
    """
    rng = np.random.default_rng(seed)
    n = n_unit * t
    unit = np.tile(np.arange(n_unit), t)
    trait = np.repeat(np.arange(t), n_unit)
    q = q_gen or max(6, n_unit // 4)

    # --- G side: us if dense, DIAG otherwise ------------------------------------
    # `diag` is required, not `iid`: iid carries ONE parameter whatever t, while
    # a diag residual carries t. The number of parameters would then not match
    # between the two sides, and the time difference would be confounded with a
    # difference in p. That would destroy the comparison this suite exists for.
    # With diag on both sides, "us on G" and "us on R" have exactly
    # t + t(t+1)/2 parameters.
    struct_g = "us" if g_dense else "diag"
    LK = None
    if g_dense:
        # dense genomic kinship, built after VanRaden on markers
        M = rng.binomial(2, 0.3, size=(q, n_marq)).astype(np.float64)
        pf = M.mean(0) / 2
        W = M - 2 * pf
        K = W @ W.T / (2 * np.sum(pf * (1 - pf)))
        K += np.eye(q) * 1e-6 * float(np.mean(np.diag(K)))
        LK = np.linalg.cholesky(K)

    zj = np.tile(rng.integers(0, q, size=n_unit), t)
    terms = [dict(name="g", struct=struct_g, t=t, rank=0, q=q,
                  zi=np.arange(n), zj=zj.astype(np.int64),
                  zx=np.ones(n), LK=LK)]

    # --- R side: us between traits if dense, diag otherwise --------------------
    res = dict(struct=("us" if r_dense else "diag"), t=t, rank=0,
               trait=trait.astype(np.int64), unit=unit)

    y = 2.0 + rng.normal(size=n)
    X = np.ones((n, 1))
    return terms, res, y, X


def suite_gr(rows, args, meta):
    """The G x R cells, at matched p, q and n."""
    ts = [int(v) for v in (args.ts or "2,4,8").split(",")]
    for t in ts:
        n_unit = max(50, args.n_gen // t)
        for g_dense, r_dense in ((False, False), (True, False),
                                 (False, True), (True, True)):
            name = "G%s/R%s" % ("dense" if g_dense else "sparse",
                                "dense" if r_dense else "sparse")
            terms, res, y, X = design_gr(n_unit, t, seed=11, g_dense=g_dense,
                                         r_dense=r_dense, q_gen=args.q_gen)
            n = len(y)
            if n > args.nmax:
                continue
            try:
                t_first, t_run, v = timed_objective(terms, res, y, X, args.reps)
                if args.no_fit:
                    med, out = float("nan"), {"n_par": -1, "logLik": float("nan"),
                                              "n_iter": -1}
                else:
                    med, lo, hi, out = timed_fit(terms, res, y, X, 1, hessian=False)
                ok, err = "ok", ""
            except Exception as e:
                t_first = t_run = v = med = float("nan")
                out = {"n_par": -1, "logLik": float("nan"), "n_iter": -1}
                ok, err = "failed", type(e).__name__ + ": " + str(e)[:200]
            rows.append(dict(meta, suite="gr", case=name, n=n, q=terms[0]["q"],
                             n_traits=t, g_dense=int(g_dense), r_dense=int(r_dense),
                             status=ok, error=err,
                             first_call_s=t_first, run_call_s=t_run,
                             compile_s=(t_first - t_run) if t_first == t_first else float("nan"),
                             fit_s_median=med, n_par=int(out["n_par"]),
                             n_iter=int(out.get("n_iter") or -1),
                             logLik=float(out["logLik"]),
                             rss_mb=rss_mb(), gpu_mb=gpu_mem_mb()))
            print("[gr] %-14s t=%-2d n=%-6d q=%-5d p=%-4s eval %8.4f s | fit %9.2f s %s"
                  % (name, t, n, terms[0]["q"], out["n_par"], t_run, med, err), flush=True)


SUITES = {"scaling": suite_scaling, "compile": suite_compile,
          "params": suite_params, "structures": suite_structures, "gr": suite_gr,
          "genomic": suite_genomic, "crossterm": suite_crossterm}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--backend", default="cpu", choices=["cpu", "gpu"])
    ap.add_argument("--suite", default="all",
                    help="all, or a list: scaling,compile,params,structures")
    ap.add_argument("--out", required=True, help="path of the output CSV")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--reps-upto", dest="reps_upto", type=int, default=4000,
                    help="beyond this n, a single repetition (cost in n^3)")
    ap.add_argument("--nmax", type=int, default=16000,
                    help="upper bound on n (device memory)")
    ap.add_argument("--ns", default=None, help="explicit list of n")
    ap.add_argument("--ts", default=None, help="list of numbers of traits")
    ap.add_argument("--qs", default=None, help="list of kinship sizes (genomic suite)")
    ap.add_argument("--n-gen", dest="n_gen", type=int, default=4000,
                    help="fixed n of the genomic and crossterm suites")
    ap.add_argument("--q-gen", dest="q_gen", type=int, default=2000,
                    help="fixed q of the crossterm suite")
    ap.add_argument("--no-fit", dest="no_fit", action="store_true",
                    help="time only the EVALUATIONS, not the fits. "
                         "On dense kinship a fit costs 100 to 900 s while one "
                         "evaluation costs 0.05 to 1.3 s. The evaluation carries "
                         "the cost comparison between backends. The fit mostly adds "
                         "the number of iterations, measured separately.")
    ap.add_argument("--n-unit", dest="n_unit", type=int, default=600)
    ap.add_argument("--grid-rows", dest="grid_rows", type=int, default=40)
    ap.add_argument("--grid-cols", dest="grid_cols", type=int, default=30)
    ap.add_argument("--tag", default="", help="free label (e.g. a100-mig-1g10gb)")
    a = ap.parse_args(argv)

    import jax
    from remlax.device import pick_device, device_report
    dev, plat = pick_device(a.backend)
    ctx = jax.default_device(dev)
    ctx.__enter__()

    meta = dict(backend=plat, device=str(dev), device_kind=getattr(dev, "device_kind", ""),
                tag=a.tag, jax=jax.__version__, numpy=np.__version__,
                python=platform.python_version(), machine=platform.machine(),
                n_cpu=os.cpu_count(),
                xla_flags=os.environ.get("XLA_FLAGS", ""),
                omp=os.environ.get("OMP_NUM_THREADS", ""),
                # os.cpu_count() returns the cores of the MACHINE, not those the
                # cgroup allocates. On a 192-core node allocated 8 cores, JAX
                # would size its pool on 192 and oversubscribe by a factor of
                # 24. Both counts are therefore logged, and the measured
                # throughput tells which one was used.
                n_cpu_affinity=(len(os.sched_getaffinity(0))
                                if hasattr(os, "sched_getaffinity") else -1),
                mkl=os.environ.get("MKL_NUM_THREADS", ""),
                openblas=os.environ.get("OPENBLAS_NUM_THREADS", ""),
                utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    print("[bench] %s" % json.dumps(meta), flush=True)
    print("[bench] devices: %s" % device_report(), flush=True)

    # Warm-up outside the measurement: the first compilation of the process also
    # pays the backend initialisation, which has nothing to do with the model.
    w = design_iid(200, seed=99)
    timed_fit(*w, reps=1, hessian=False)

    rows = []
    names = list(SUITES) if a.suite == "all" else a.suite.split(",")
    for nm in names:
        if nm not in SUITES:
            raise SystemExit("unknown suite: %s (choose from: %s)" % (nm, ", ".join(SUITES)))
        print("\n=== suite %s ===" % nm, flush=True)
        SUITES[nm](rows, a, meta)

    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("\n[bench] %d measurements written to %s" % (len(rows), a.out), flush=True)
    ctx.__exit__(None, None, None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
