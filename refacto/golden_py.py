"""Python safety net of the remlax cleanup: engine values at fixed theta.

    python3 golden_py.py capture <out.pkl> <remlax_repo>
    python3 golden_py.py compare <ref.pkl> <new.pkl>

For each case of cases_py.py and two values of theta, it records:
  - dense engine: -2logL and gradient (model.make_objective), and value,
    gradient and AI matrix (ai.make_fga);
  - low-memory engine: value, gradient and AI (econome.make_fga_econome);
  - MME engine, for each inversion (dense, sparse, hybrid with the separator
    of the case): value, gradient and AI, the value-only then full evaluation
    path at the same theta, and the final outputs (beta, vbeta, BLUPs, PEV);
  - blocks engine when the model splits: value and gradient;
  - whole fits (AI on the MME and dense engines, L-BFGS-B dense): theta,
    -2logL, iterations, evaluations, message, beta.
compare reports every entry that is not identical to the bit.

Runs pinned to one CPU core, without the compilation cache, like golden.R.
"""
import os
import sys

# One core: XLA on CPU is not deterministic across runs when the process may
# use several cores, even with the thread flags below (about one run in six
# differed at the 1e-13 level on the dense AI matrix). Pinned to one core,
# repeated captures are identical. The first allowed core is used.
if hasattr(os, "sched_setaffinity"):
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
os.environ["JAX_PLATFORMS"] = "cpu"
os.environ["XLA_FLAGS"] = "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1"
os.environ["REMLAX_JAX_CACHE"] = "0"
# Folder where a development checkout compiles the selected inversion kernel
# (src/selinv_c.c) when REMLAX_SELINV_LIB is not set.
os.environ.setdefault("REMLAX_CACHE", os.path.join(os.environ.get("TMPDIR", "/tmp"), "golden_py_cache"))

import pickle  # noqa: E402

import numpy as np  # noqa: E402

VOLATILE = {"secondes", "compile_s", "eval_s", "first_call_s", "secondes_optim", "chrono",
            "durees", "temps"}


def _try(fn):
    try:
        return fn()
    except Exception as e:                       # the failure itself is recorded
        return {"erreur": "%s: %s" % (type(e).__name__, e)}


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items() if k not in VOLATILE}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if hasattr(o, "shape") and not isinstance(o, np.ndarray):
        return np.asarray(o)
    return o


def capture(out, repo):
    sys.path.insert(0, os.path.join(repo, "inst", "python"))
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from cases_py import CASES
    import remlax  # noqa: F401
    from remlax.model import make_objective
    from remlax.ai import make_fga
    from remlax.mme import MME
    from remlax.fit import fit_reml
    try:
        from remlax.start import initial_theta
    except ImportError:
        from remlax.fit import initial_theta
    from remlax import econome
    from remlax import blocs

    res_all = {}
    for name, (build, sep) in CASES.items():
        terms, res, y, X = build()
        th0 = initial_theta(terms, res, y)
        rng = np.random.default_rng(13)
        thetas = [th0, th0 + rng.normal(0, 0.3, len(th0))]
        r = {"thetas": thetas}
        kinv_only = any(tm.get("Kinv") is not None and tm.get("LK") is None for tm in terms)
        if not kinv_only:
            def dense():
                fsc, fj, sc, _ = make_objective(terms, res, y, X)
                fga = make_fga(terms, res, y, X)
                return [dict(fj=fj(th), fga=fga(th)) for th in thetas]
            r["dense"] = _try(dense)

            def low():
                Plan = getattr(econome, "LowMemoryPlan", None) or econome.PlanEconome
                mk = getattr(econome, "make_fga_low_memory", None) or econome.make_fga_econome
                plan = Plan(terms, res, y, X)
                fga = mk(plan)
                return [fga(th) for th in thetas]
            r["econome"] = _try(low)

            def blk():
                oos = getattr(blocs, "out_of_scope", None) or blocs.hors_perimetre
                Plan = getattr(blocs, "BlockPlan", None) or blocs.PlanBlocs
                mk = getattr(blocs, "make_objective_blocks", None) or blocs.make_objective_blocs
                if oos(terms, res, None) is not None:
                    return "hors perimetre"
                plan = Plan(terms, res, y, X, verbose=False)
                nb = plan.n_blocks if hasattr(plan, "n_blocks") else plan.n_blocs
                if nb < 2:
                    return "un seul bloc"
                fsc, fj, sc = mk(plan)
                return [fj(th) for th in thetas]
            r["blocs"] = _try(blk)
        for inv in ("dense", "creuse", "hybride"):
            def mme(inv=inv):
                kw = dict(separateur=sep) if inv == "hybride" and sep else {}
                M = MME(terms, res, y, X, inversion=inv, verbose=False, **kw)
                ev_fn = getattr(M, "evaluate", None) or M.evaluer
                out_fn = getattr(M, "outputs", None) or M.sorties
                ev = []
                for th in thetas:
                    e = ev_fn(th, True)
                    ev.append({k: e[k] for k in ("f", "g", "AI") if k in e})
                M2 = MME(terms, res, y, X, inversion=inv, verbose=False, **kw)
                ev2 = getattr(M2, "evaluate", None) or M2.evaluer
                a = ev2(thetas[1], False)
                b = ev2(thetas[1], True)
                return dict(inversion=M.inversion, evals=ev,
                            two_step=dict(f=a.get("f"), g=b.get("g"), AI=b.get("AI")),
                            sorties=out_fn(thetas[0], pev=True))
            r["mme_" + inv] = _try(mme)
        fits = {}
        kw = dict(verbose=False, hessian=False, pev=False)
        fits["ai_mme"] = _try(lambda: fit_reml(terms, res, y, X, optimizer="ai", moteur="mme", **kw))
        if not kinv_only:
            fits["ai_dense"] = _try(lambda: fit_reml(terms, res, y, X, optimizer="ai", **kw))
            fits["lbfgsb_dense"] = _try(lambda: fit_reml(terms, res, y, X, maxiter=200, **kw))
        keep = ("theta", "neg2_reml", "n_iter", "n_eval", "n_polish", "scipy_message",
                "scipy_success", "beta", "moteur", "inversion", "grad_rel", "max_grad")
        r["fits"] = {k: ({kk: v[kk] for kk in keep if kk in v} if "erreur" not in v else v)
                     for k, v in fits.items()}
        res_all[name] = _clean(r)
        print("captured", name, flush=True)
    with open(out, "wb") as f:
        pickle.dump(res_all, f)


def _diff(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                out.append("%s/%s: present on one side only" % (path, k))
            else:
                _diff(a[k], b[k], "%s/%s" % (path, k), out)
        return
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b):
            out.append("%s: lengths %d and %d" % (path, len(a), len(b)))
            return
        for i, (x, y) in enumerate(zip(a, b)):
            _diff(x, y, "%s[%d]" % (path, i), out)
        return
    if isinstance(a, (np.ndarray, float, int, np.floating, np.integer)) and \
            isinstance(b, (np.ndarray, float, int, np.floating, np.integer)):
        A, B = np.asarray(a), np.asarray(b)
        if A.shape != B.shape or not np.array_equal(A, B, equal_nan=True):
            gap = (float(np.nanmax(np.abs(A - B))) if A.shape == B.shape and A.size else float("nan"))
            out.append("%s: differs (max gap %.3g)" % (path, gap))
        return
    if a != b:
        out.append("%s: %r != %r" % (path, a, b))


def compare(ref, new):
    with open(ref, "rb") as f:
        A = pickle.load(f)
    with open(new, "rb") as f:
        B = pickle.load(f)
    out = []
    _diff(A, B, "", out)
    if out:
        print("%d DIFFERENCE(S)" % len(out))
        for line in out[:200]:
            print("  " + line)
        sys.exit(1)
    n = sum(1 for _ in A)
    print("%d cases compared: NO DIFFERENCE, identical to the bit." % n)


if __name__ == "__main__":
    if sys.argv[1] == "capture":
        capture(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "compare":
        compare(sys.argv[2], sys.argv[3])
    else:
        sys.exit("usage: golden_py.py capture <out.pkl> <repo> | compare <ref.pkl> <new.pkl>")
