"""Where GPU and CPU differ, and why: the paired benchmark.

WHY THIS FILE REPLACES THE EARLIER GRID
=======================================
The grid compared FIT TIMES between backends. That measurement is biased, and
a reviewer would see it: the two backends do not follow the same optimisation
path. Measured on one cell, 419 iterations on the card against 317 on CPU, for
the same optimum within 8.6e-10 in relative terms. Floating-point addition is
not associative, and the reductions run in a different order. The gradient
therefore differs in the last bit, the line search accepts a different step,
and the paths diverge. The number of iterations is decided at the end of the
run by rounding noise. Comparing two total times thus compares a hardware cost
MULTIPLIED by a random path length.

THE CLEAN UNIT OF MEASUREMENT: THE EVALUATION AT AN IMPOSED THETA
=================================================================
An evaluation of the likelihood at the SAME theta does exactly the same
arithmetic on both sides. No optimiser, no path: the ratio is purely a
hardware one, deterministic and repeatable. The common theta is obtained ONCE
by a reference fit, then imposed on both backends.

THE DECOMPOSITION INTO THREE PIECES, MEASURED SEPARATELY
========================================================
Instead of a global regression, three quantities are measured, and their
product must rebuild the total time. Each one can be checked alone, and the
rebuild is a test this benchmark can FAIL.

  1. compilation   paid once, a property of the device (XLA compiles for the
                   target); more expensive on GPU
  2. evaluation    steady state; depends on the structure but NOT on the number
                   of free parameters, since the Cholesky of V dominates and
                   ignores p
  3. evaluations   n_iter + 2p: the finite-difference Hessian costs exactly 2p
                   gradients, so this part is PREDICTABLE

HYPOTHESIS UNDER TEST, AND THE EXPECTED RESULT
==============================================
The GPU pays a more expensive compilation and wins on each evaluation. At
small p there are few evaluations, compilation dominates, and the GPU loses.
At large p the number of evaluations grows, through the iterations AND through
the 2p of the Hessian. Compilation is then amortised and the GPU wins. The
CPU/GPU crossover therefore happens in p, and not only in n. This can be
refuted.

WHAT MAKES THE DESIGN DEFENSIBLE
================================
  - paired       same data, same seed, same theta; only the backend changes
  - fixed cores  16 everywhere, since they drive the compilation time
  - whole card   asserted by the calling script, which exits with an error on
                 a slice: a slice gives plausible numbers that are silently a
                 pessimistic bound
  - repetitions  median and range, since times vary between repetitions
  - empty cell   DECLARED: the sparse engine has no GPU route, CHOLMOD being
                 CPU only. The grid has six cells and not eight. This absence
                 is a result: it makes 'sparse' and 'GPU' exclusive.

HOW p IS VARIED, AND THE BIAS THAT REMAINS
==========================================
On STRICTLY IDENTICAL DATA: same n, q, t, same Z, same y, same seed. Only the
covariance STRUCTURE on the genetic term changes, which is also what a real
user does. With t traits:

    iid   1 parameter        us       t(t+1)/2 parameters
    diag  t parameters       fa(t,k)  in between

The remaining bias is stated here. Changing the structure also changes the
sparsity PATTERN of G^-1 and the geometry of the surface. p is not a physical
dimension of the problem but a property of the structure, so no single knob
turns it alone. There are two safeguards.
  - The cost per evaluation is reported apart from the number of evaluations.
    This separates a pattern effect (unit cost) from a dimension effect
    (number of iterations).
  - The family is given on each row, because reduced rank does not compare
    with full rank. In an earlier measurement, fa at t=8 needed 258 iterations
    against 26 for us, WITH FEWER parameters (31 against 44).

Usage:
    python3 bench_backend.py --backend cpu --out cpu.csv --tag cpu16
    python3 bench_backend.py --backend gpu --out gpu.csv --tag card
"""
import argparse, csv, json, os, platform, statistics, sys, time
import numpy as np


# ==============================================================================
# Design: ONE data set per (n, q, t), reused for every p
# ==============================================================================
def design(n_unit, t, q, seed=7, LK=None):
    """Multi-trait design with q genotypes. The data do NOT depend on the fitted
    structure, which allows p to vary on identical data.

    BEWARE OF A TRAP ALREADY MET: if the number of units falls below q, some
    levels are never observed and the system becomes nearly singular. This
    then reads as a cost of p when the design is degenerate. It is forbidden
    here rather than detected afterwards.
    """
    if n_unit < 2 * q:
        raise ValueError("n_unit = %d < 2q = %d: insufficient replication, "
                         "some levels would be unobserved or unreplicated"
                         % (n_unit, 2 * q))
    rng = np.random.default_rng(seed)
    n = n_unit * t
    unit = np.repeat(np.arange(n_unit), t)
    trait = np.tile(np.arange(t), n_unit)
    lev = unit % q
    # genetic effects correlated between traits, full structure
    A = rng.normal(0, 1.0, (t, t)) / np.sqrt(t)
    G = A @ A.T + 0.5 * np.eye(t)
    u = rng.multivariate_normal(np.zeros(t), G, size=q)
    # The kinship must be used FOR SIMULATION as well as for fitting. Fitting a K
    # on data that do not carry it does not change the COST, since the
    # arithmetic is the same, but it makes the likelihoods and the number of
    # iterations unrepresentative. The two other factories already did this;
    # this one did not.
    if LK is not None:
        u = np.asarray(LK) @ u
    y = 1.5 + u[lev, trait] + rng.normal(0, 0.8, n)
    X = np.zeros((n, t))
    X[np.arange(n), trait] = 1.0                       # one mean per trait
    # COLUMN INDEX: dense_Z builds Z with t*q columns, the random effect being
    # the vector vec(Sigma (x) K) ordered BY TRAIT. The BLUP code computes
    # kron(Ls, Kq) then reshape(t, q), so the column is
    #     trait * q + level
    # and not the level alone. With zj = lev, every observation loaded the same
    # block of q columns. The product Z (Sigma (x) K) Z' then involved only
    # Sigma[0,0], so iid, diag and us returned EXACTLY the same likelihood
    # despite different p. The pairing check of the benchmark caught this
    # defect, not a review.
    zj = (trait * q + lev).astype(np.int64)
    return dict(n=n, n_unit=n_unit, t=t, q=q, y=y, X=X,
                zi=np.arange(n), zj=zj, zx=np.ones(n),
                trait=trait.astype(np.int64), unit=unit)



# ==============================================================================
# KINSHIP: the same design, with or without a dense K
# ==============================================================================
def kinship_chol(q, n_marq=1500, seed=11):
    """Cholesky factor of a DENSE VanRaden genomic kinship.

    This is the structural axis. The same model with K = I is structurally
    sparse (diagonal precision); with a genomic K it is structurally dense
    (full K^-1). The term receives it through the LK field, which term_factor
    treats as a 'fixed' level. BEWARE, as the code of model.py says: a K given
    WITHOUT an 'lvl' key is treated as fixed, and falling back to 'id' would
    silently ignore the kinship. This failure already happened on the GPU port
    of the IGE model, where two runs labelled 'kinship' were exact duplicates of
    the runs without it.
    """
    rng = np.random.default_rng(seed)
    fr = rng.uniform(0.05, 0.95, n_marq)
    M = rng.binomial(2, fr, size=(q, n_marq)).astype(np.float64)
    Zc = M - 2.0 * fr
    K = Zc @ Zc.T / (2.0 * np.sum(fr * (1.0 - fr)))
    K += np.eye(q) * 1e-6                      # minimal regularisation
    return np.linalg.cholesky(K)


# ==============================================================================
# Spatial design: field trial in rows x columns, ar1 x ar1 residual
# ==============================================================================
def design_spatial(nr, nc, q, seed=13, LK=None):
    """The classic field trial: one plot per cell of a grid.

    The residual carries ar1(row) x ar1(col), THE structure of agronomic trials
    since Gilmour and Cullis. It is separable, so its two correlation
    parameters add to the variance parameters.
    """
    rng = np.random.default_rng(seed)
    n = nr * nc
    if n < 2 * q:
        raise ValueError("grid %dx%d = %d plots < 2q = %d" % (nr, nc, n, 2 * q))
    row = np.repeat(np.arange(nr), nc)
    col = np.tile(np.arange(nc), nr)
    lev = rng.permutation(np.repeat(np.arange(q), int(np.ceil(n / q)))[:n])
    u = rng.normal(0, 1.0, q) if LK is None else LK @ rng.normal(0, 1.0, q)
    y = 2.0 + u[lev] + rng.normal(0, 0.7, n)
    return dict(n=n, n_unit=n, t=1, q=q, y=y, X=np.ones((n, 1)),
                zi=np.arange(n), zj=lev.astype(np.int64), zx=np.ones(n),
                trait=np.zeros(n, np.int64), unit=np.arange(n),
                row=row, col=col, dims=(nr, nc), LK=LK)


# ==============================================================================
# Direct and indirect genetic effects: two incidences, one us(2)
# ==============================================================================
def design_ige(n_unit, q, neighbours=4, seed=17, LK=None):
    """DGE and IGE: the effect of a genotype on itself and on its neighbours.

    ONE term with t = 2, whose two 'traits' are the direct and the indirect
    effect, with a us(2) between them. This is the classic parametrisation,
    and the correlation between the two is the quantity of interest. The two
    columns of Z differ: identity of the plot genotype for the direct effect,
    sum of the neighbours for the indirect one. The indirect incidence is DENSE
    in columns even when the precision is diagonal, which fills the block of
    cross products. This is the second fill-in mechanism identified in the
    session, distinct from the density of K.
    """
    rng = np.random.default_rng(seed)
    if n_unit < 2 * q:
        raise ValueError("n_unit = %d < 2q = %d" % (n_unit, 2 * q))
    gen = rng.permutation(np.repeat(np.arange(q), int(np.ceil(n_unit / q)))[:n_unit])
    # neighbourhood: `neighbours` plots drawn at random, weight 1/neighbours
    vz_i, vz_j, vz_x = [], [], []
    for i in range(n_unit):
        vs = rng.choice(n_unit, size=neighbours, replace=False)
        for v in vs:
            vz_i.append(i); vz_j.append(q + gen[v]); vz_x.append(1.0 / neighbours)
    zi = np.concatenate([np.arange(n_unit), np.asarray(vz_i)])
    zj = np.concatenate([gen, np.asarray(vz_j)]).astype(np.int64)   # direct: 0..q-1
    zx = np.concatenate([np.ones(n_unit), np.asarray(vz_x)])        # indirect: q..2q-1
    ud = rng.normal(0, 1.0, q); ui = 0.6 * ud + 0.8 * rng.normal(0, 0.6, q)
    if LK is not None:
        ud, ui = LK @ ud, LK @ ui
    y = np.full(n_unit, 2.0)
    np.add.at(y, np.arange(n_unit), ud[gen])
    for i, j, x in zip(vz_i, vz_j, vz_x):
        y[i] += x * ui[j - q]
    y += rng.normal(0, 0.7, n_unit)
    return dict(n=n_unit, n_unit=n_unit, t=2, q=q, y=y, X=np.ones((n_unit, 1)),
                zi=zi, zj=zj, zx=zx,
                trait=np.zeros(n_unit, np.int64), unit=np.arange(n_unit), LK=LK)


def model(d, struct, rank=0, res_struct="diag", res_lvl=None):
    """LK is passed on FROM the design. This is the only place where the
    'with or without kinship' axis is decided, so that it cannot be lost on
    the way. A K given without an 'lvl' key is treated as a fixed level by
    term_factor, which is the intended behaviour.
    """
    terms = [dict(name="g", struct=struct, t=d["t"], rank=rank, q=d["q"],
                  zi=d["zi"], zj=d["zj"], zx=d["zx"], LK=d.get("LK"))]
    res = dict(struct=res_struct, t=d["t"], rank=0,
               trait=d["trait"], unit=d["unit"])
    if res_lvl == "ar1ar1":
        # separable spatial residual: C_unit = AR1(rows) x AR1(columns)
        res.update(lvl="ar1ar1", lvl_parts=("ar1", "ar1"), dims=d["dims"])
    return terms, res


# ==============================================================================
# Measurements
# ==============================================================================
def time_eval(fit_reml, terms, res, d, theta, reps):
    """Evaluation AT AN IMPOSED THETA: maxiter = 0 bypasses the optimiser, so
    this evaluates and does not move. It is the paired measurement between
    backends. The first pass pays the compilation, which the solver reports.
    """
    ts, out = [], None
    for _ in range(reps):
        t0 = time.perf_counter()
        out = fit_reml(terms, res, d["y"], d["X"], theta_init=theta,
                       maxiter=0, polish=0, hessian=False, blups=False, verbose=False)
        ts.append(time.perf_counter() - t0)
    return out, ts


def time_fit(fit_reml, terms, res, d, maxiter, hessian):
    t0 = time.perf_counter()
    out = fit_reml(terms, res, d["y"], d["X"], maxiter=maxiter, polish=2,
                   hessian=hessian, blups=False, verbose=False)
    return out, time.perf_counter() - t0


def structures(t, families):
    """(name, struct, rank) for each requested structure, by increasing p.

    SATURATION OF THE REDUCED RANK, not to be read as a defect. fa(t, k) is a
    restriction of us(t) only if the rank is strictly insufficient. At t = 2,
    a rank-1 factor plus the specific variances represents ANY 2x2 matrix.
    fa(2,1) and us(2) are therefore the SAME model and return the same
    likelihood, fa even carrying MORE parameters (6 against 5, measured). The
    nesting check must allow it. The case stays in the benchmark because it
    shows exactly what the session measured elsewhere: reduced rank does not
    always save parameters, and it costs iterations (23 against 14 here at
    t=2, for the same optimum).
    """
    cat = []
    for f in families:
        if f == "iid":
            cat.append(("iid", "iid", 0))
        elif f == "diag":
            cat.append(("diag", "diag", 0))
        elif f == "us":
            cat.append(("us", "us", 0))
        elif f.startswith("fa"):
            k = int(f[2:])
            if k < t:
                cat.append(("fa%d" % k, "fa", k))
    return cat


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--backend", default="cpu", choices=["cpu", "gpu"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="")
    # DRIVEN BY THE TOTAL, NOT BY THE NUMBER OF UNITS. The dense engine forms V
    # explicitly: at n = 32000 it weighs 8 GB in float64, and three must be
    # alive. Fixing n_unit and letting n = n_unit * t grow with t made large t
    # infeasible. With n fixed, the p axis is covered at CONSTANT problem size,
    # which is also what this benchmark measures.
    ap.add_argument("--n", type=int, default=4000,
                    help="TOTAL observations; n_unit = n / t")
    # q must stay below n_unit/2 for the largest requested t, otherwise some
    # levels would be unreplicated. That defect distorted the earlier grid.
    ap.add_argument("--q", type=int, default=200)
    ap.add_argument("--ts", default="2,4,8")
    ap.add_argument("--families", default="iid,diag,fa1,fa2,fa4,us",
                    help="covariance structures for the multi-trait family")
    ap.add_argument("--models", default="multi,spatial,ige",
                    help="MODEL families: multi, spatial, ige")
    ap.add_argument("--kinship", default="both", choices=["both", "yes", "no"],
                    help="with a dense genomic kinship, without, or both")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--maxiter", type=int, default=3000)
    ap.add_argument("--no-fit", action="store_true",
                    help="evaluations only; the complete fit is the expensive item")
    a = ap.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import moteur_remlax  # noqa: F401
    import jax
    from remlax.device import pick_device, device_report
    from remlax.fit import fit_reml
    dev, plat = pick_device(a.backend)
    ctx = jax.default_device(dev)
    ctx.__enter__()

    # The number of VISIBLE cores drives the compilation time. A different
    # value between two jobs makes their times incomparable, so it is recorded.
    try:
        ncpu = len(os.sched_getaffinity(0))
    except AttributeError:
        ncpu = os.cpu_count()
    meta = dict(backend=plat, device=str(dev), tag=a.tag,
                device_kind=getattr(dev, "device_kind", ""),
                n_cpu_affinity=ncpu, omp=os.environ.get("OMP_NUM_THREADS", ""),
                jax=jax.__version__, python=platform.python_version(),
                utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    print("[bench] %s" % json.dumps(meta), flush=True)
    print("[bench] devices: %s" % device_report(), flush=True)

    # Warm-up outside the measurement: the very first compilation of the
    # process also pays the backend initialisation, unrelated to the model.
    dw = design(200, 2, 40, seed=99)
    tw, rw = model(dw, "diag")
    fit_reml(tw, rw, dw["y"], dw["X"], maxiter=3, polish=0, hessian=False,
             blups=False, verbose=False)

    rows = []
    # ---- catalogue of families, each one WITH and WITHOUT kinship ----------
    # These are the models actually fitted in plant breeding, not synthetic
    # designs: basic GBLUP, multi-trait and multi-environment, reduced-rank
    # factors, field trial, indirect effects.
    def cells(a):
        ts = [int(v) for v in a.ts.split(",")]
        fam = a.models.split(",")
        for kin in ([False, True] if a.kinship == "both" else [a.kinship == "yes"]):
            LK = kinship_chol(a.q) if kin else None
            lab = "K_genomic" if kin else "K_identity"
            if "multi" in fam:
                for t in ts:
                    d = design(a.n // t, t, a.q, LK=LK)
                    for name, st, rk in structures(t, a.families.split(",")):
                        yield ("multi/%s(t=%d)" % (name, t), lab, d,
                               dict(struct=st, rank=rk))
            if "spatial" in fam:
                nr = int(round(a.n ** 0.5 * 1.25)); nc = max(1, a.n // nr)
                d = design_spatial(nr, nc, a.q, LK=LK)
                yield ("spatial/ar1xar1", lab, d,
                       dict(struct="iid", rank=0, res_struct="iid", res_lvl="ar1ar1"))
            if "ige" in fam:
                d = design_ige(a.n, a.q, LK=LK)
                yield ("ige/us(2)", lab, d, dict(struct="us", rank=0, res_struct="iid"))

    ts = [int(v) for v in a.ts.split(",")]
    n_unit_min = a.n // max(ts)
    if n_unit_min < 2 * a.q:
        raise SystemExit("q = %d too large: at the largest t = %d only %d units remain, "
                         "at least 2q = %d are needed. Lower --q or raise --n."
                         % (a.q, max(ts), n_unit_min, 2 * a.q))
    for name, lab, d, kw in cells(a):
        t = d["t"]
        if True:
            terms, res = model(d, **kw)
            base = dict(meta, n=d["n"], n_unit=d["n_unit"], t=t, q=a.q,
                        structure=name, kinship=lab,
                        family=kw["struct"], rank=kw.get("rank", 0))
            # --- 1. reference fit: gives the common theta AND the total
            try:
                out, wall = time_fit(fit_reml, terms, res, d, a.maxiter,
                                     hessian=not a.no_fit)
            except Exception as e:
                rows.append(dict(base, phase="fit", status="fail",
                                 error=type(e).__name__ + ": " + str(e)[:180]))
                print("[fit ] %-22s %-12s FAIL %s" % (name, lab, type(e).__name__), flush=True)
                continue
            p = int(out["n_par"])
            rows.append(dict(base, phase="fit", status="ok", p=p,
                             n_iter=int(out.get("n_iter", -1)),
                             fit_wall_s=wall, fit_internal_s=float(out["seconds"]),
                             compile_s=float(out.get("compile_s", float("nan"))),
                             eval_s=float(out.get("eval_s", float("nan"))),
                             logLik=float(out["logLik"]),
                             newton_decrement=float(out.get("newton_decrement", float("nan"))),
                             grad_rel=float(out.get("grad_rel", float("nan"))),
                             conv_decrement=out.get("conv_decrement"),
                             truncated=bool(out.get("n_iter", 0) >= a.maxiter),
                             n_eval_predicted=int(out.get("n_iter", 0)) + 2 * p * (0 if a.no_fit else 1)))
            print("[fit ] %-22s %-12s p=%-3d iter=%-4d internal %8.2f s (compile %6.2f s)"
                  % (name, lab, p, out.get("n_iter", -1), out["seconds"],
                     out.get("compile_s", float("nan"))), flush=True)

            # --- 2. evaluations AT AN IMPOSED THETA: the paired measurement
            th = np.asarray(out["theta"], dtype=np.float64)
            try:
                oe, ts = time_eval(fit_reml, terms, res, d, th, a.reps)
            except Exception as e:
                rows.append(dict(base, phase="eval", status="fail",
                                 error=type(e).__name__ + ": " + str(e)[:180]))
                continue
            rows.append(dict(base, phase="eval", status="ok", p=p, reps=a.reps,
                             eval_wall_median_s=statistics.median(ts),
                             eval_wall_min_s=min(ts), eval_wall_max_s=max(ts),
                             compile_s=float(oe.get("compile_s", float("nan"))),
                             eval_s=float(oe.get("eval_s", float("nan"))),
                             logLik=float(oe["logLik"]),
                             # pairing check: the likelihood at an imposed theta
                             # must be IDENTICAL between backends. If it differs,
                             # it is not the same computation and the time
                             # ratio means nothing.
                             theta_imposed=" ".join("%.10g" % v for v in th)))
            print("[eval] %-22s %-12s p=%-3d compile %7.3f s | evaluation %8.5f s | logLik %.9f"
                  % (name, lab, p, oe.get("compile_s", float("nan")),
                     oe.get("eval_s", float("nan")), oe["logLik"]), flush=True)

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
    print("\n[bench] %d rows written to %s" % (len(rows), a.out), flush=True)
    ctx.__exit__(None, None, None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
