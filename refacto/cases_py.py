"""Model cases of the Python safety net (golden_py.py).

Copied from tests-python (test_econome, test_blocs, test_face_cone) at commit
bbc7340 of remlax, so that later changes to the tests cannot move the
reference. Three cases are added for MME routes the tests reach less often: a
pedigree precision Kinv, a parametrised dense kernel (exp) and a residual on a
grid with missing cells.
"""
import numpy as np


def terme_facteur(name, codes, q, t=1, struct="iid", rank=0, LK=None, **extra):
    codes = np.asarray(codes, dtype=np.int64)
    n = len(codes)
    zi = np.tile(np.arange(n, dtype=np.int64), t)
    zj = np.concatenate([a * q + codes for a in range(t)]).astype(np.int64)
    d = dict(name=name, struct=struct, t=int(t), rank=int(rank), q=int(q),
             zi=zi, zj=zj, zx=np.ones(n * t), LK=LK)
    d.update(extra)
    return d


def terme_long(name, codes, trait, q, t, struct="us", rank=0, LK=None, **extra):
    codes = np.asarray(codes, dtype=np.int64)
    trait = np.asarray(trait, dtype=np.int64)
    n = len(codes)
    d = dict(name=name, struct=struct, t=int(t), rank=int(rank), q=int(q),
             zi=np.arange(n, dtype=np.int64), zj=(trait * q + codes).astype(np.int64),
             zx=np.ones(n), LK=LK)
    d.update(extra)
    return d


def residuelle(struct, trait, unit, t, rank=0, lvl="id", **extra):
    d = dict(struct=struct, t=int(t), rank=int(rank),
             trait=np.asarray(trait, np.int64), unit=np.asarray(unit, np.int64),
             n_unit=int(np.max(unit)) + 1, lvl=lvl)
    d.update(extra)
    return d


def _chol_spd(q, rng):
    A = rng.normal(size=(q, q))
    K = A @ A.T / q + 0.5 * np.eye(q)
    return np.linalg.cholesky(K)


def cas_ige(graine=1, n_bacs=3, n_ble=40, n_luz=12, t_b=3, q_b=15, q_l=6):
    rng = np.random.default_rng(graine)
    trait, unit, esp, bac, geno = [], [], [], [], []
    u = 0
    for b in range(n_bacs):
        for _ in range(n_ble):
            g = rng.integers(q_b)
            for a in range(t_b):
                if rng.random() < 0.85:
                    trait.append(a); unit.append(u); esp.append(0); bac.append(b); geno.append(g)
            u += 1
        for _ in range(n_luz):
            trait.append(0); unit.append(u); esp.append(1); bac.append(b)
            geno.append(rng.integers(q_l)); u += 1
    trait, unit, esp, bac, geno = map(np.asarray, (trait, unit, esp, bac, geno))
    n = len(trait)
    perm = rng.permutation(n)
    trait, unit, esp, bac, geno = trait[perm], unit[perm], esp[perm], bac[perm], geno[perm]
    rb, rl = np.flatnonzero(esp == 0), np.flatnonzero(esp == 1)
    terms = []
    terms.append(dict(name="gen_ble", struct="us", t=t_b, rank=0, q=q_b,
                      zi=rb, zj=trait[rb] * q_b + geno[rb], zx=np.ones(len(rb)),
                      LK=_chol_spd(q_b, rng), lvl="fixed"))
    terms.append(dict(name="gen_luz", struct="iid", t=1, rank=0, q=q_l,
                      zi=rl, zj=geno[rl], zx=np.ones(len(rl)), LK=_chol_spd(q_l, rng),
                      lvl="fixed"))
    terms.append(dict(name="plante", struct="iid", t=1, rank=0, q=int(unit.max()) + 1,
                      zi=np.arange(n), zj=unit, zx=rng.random(n) + 0.5, LK=None))
    for b in range(n_bacs):
        rr = np.flatnonzero(bac == b)
        terms.append(dict(name="spat_%d" % b, struct="iid", t=1, rank=0, q=6, zi=rr,
                          zj=rng.integers(0, 6, len(rr)), zx=np.ones(len(rr)), LK=None))
        k = 4
        zi = np.repeat(rr, k)
        terms.append(dict(name="iee_%d" % b, struct="iid", t=1, rank=0, q=n_ble + n_luz,
                          zi=zi, zj=rng.integers(0, n_ble + n_luz, len(zi)),
                          zx=rng.random(len(zi)), LK=None))
    secs = [dict(name="ble", struct="us", t=t_b, rank=0, trait=trait[rb], unit=unit[rb],
                 n_unit=int(unit.max()) + 1, lvl="id", rows=rb),
            dict(name="luz", struct="iid", t=1, rank=0, trait=np.zeros(len(rl), np.int64),
                 unit=unit[rl], n_unit=int(unit.max()) + 1, lvl="id", rows=rl)]
    res = dict(struct="iid", t=1, rank=0, trait=np.zeros(n, np.int64), unit=unit,
               n_unit=int(unit.max()) + 1, sections=secs)
    X = np.zeros((n, t_b + 1))
    X[rb, trait[rb]] = 1.0
    X[rl, t_b] = 1.0
    y = rng.normal(size=n) + X @ rng.normal(size=t_b + 1)
    return terms, res, y, X


def cas_ar1(graine=2):
    rng = np.random.default_rng(graine)
    nr, nc, t = 6, 7, 2
    cell = np.arange(nr * nc)
    unit = np.tile(cell, t)
    trait = np.repeat(np.arange(t), nr * nc)
    n = len(unit)
    res = residuelle("us", trait, unit, t, lvl="ar1ar1", dims=(nr, nc))
    q = 10
    g = rng.integers(0, q, nr * nc)[unit]
    terms = [terme_long("g", g, trait, q, t, struct="fa", rank=1),
             terme_facteur("col", cell[unit] % nc, nc, lvl="ar1")]
    X = np.zeros((n, t)); X[np.arange(n), trait] = 1.0
    y = rng.normal(size=n)
    return terms, res, y, X


def cas_sections_bloc(graine=3):
    rng = np.random.default_rng(graine)
    nA, nB = 60, 40
    n = nA + nB
    perm = rng.permutation(n)
    rowsA, rowsB = np.sort(perm[:nA]), np.sort(perm[nA:])
    secA = dict(struct="iid", t=1, rank=0, trait=np.zeros(nA, np.int64),
                unit=np.arange(nA), n_unit=nA, lvl="id", rows=rowsA, name="A")
    secB = dict(struct="iid", t=1, rank=0, trait=np.zeros(nB, np.int64),
                unit=np.arange(nB), n_unit=nB, lvl="ar1", rows=rowsB, name="B")
    res = dict(struct="iid", t=1, rank=0, trait=np.zeros(n, np.int64), unit=np.arange(n),
               n_unit=n, sections=[secA, secB])
    terms = [terme_facteur("g", rng.integers(0, 12, n), 12)]
    X = np.ones((n, 1))
    y = rng.normal(size=n)
    return terms, res, y, X


def cas_essais(graine=4, n_essais=4, n_geno=25, rep=3):
    rng = np.random.default_rng(graine)
    lignes = []
    u = 0
    for e in range(n_essais):
        for g in range(n_geno):
            for _ in range(rep):
                for a in range(2):
                    if rng.random() < 0.9:
                        lignes.append((e, g, u, a))
                u += 1
    E, Gn, U, A = map(np.asarray, zip(*lignes))
    n = len(E)
    perm = rng.permutation(n)
    E, Gn, U, A = E[perm], Gn[perm], U[perm], A[perm]
    terms = [
        dict(name="geno", struct="us", t=2, rank=0, q=n_geno, zi=np.arange(n),
             zj=A * n_geno + Gn, zx=np.ones(n), LK=_chol_spd(n_geno, rng), lvl="fixed"),
        terme_facteur("geno_ar1", Gn, n_geno, lvl="ar1"),
    ]
    for e in range(n_essais):
        rr = np.flatnonzero(E == e)
        blocs_e = rng.integers(0, 5, len(rr))
        terms.append(dict(name="bloc_%d" % e, struct="us", t=2, rank=0, q=5, zi=rr,
                          zj=A[rr] * 5 + blocs_e, zx=np.ones(len(rr)), LK=None))
    r0 = np.flatnonzero(E == 0)
    terms.append(dict(name="rang_e0", struct="iid", t=1, rank=0, q=8, zi=r0,
                      zj=rng.integers(0, 8, len(r0)), zx=np.ones(len(r0)), LK=None, lvl="ar1"))
    res = residuelle("us", A, U, 2)
    X = np.zeros((n, 2 * n_essais))
    X[np.arange(n), E * 2 + A] = 1.0
    y = rng.normal(size=n) + X @ rng.normal(size=X.shape[1])
    return terms, res, y, X


def cas_rang_2(graine=4, q=80, t=4, rep=2):
    rng = np.random.default_rng(graine)
    L = rng.normal(size=(t, 2))
    u = rng.normal(size=(q, 2)) @ L.T
    codes, trait = [], []
    for g in range(q):
        for a in range(t):
            for _ in range(rep):
                codes.append(g); trait.append(a)
    codes, trait = np.asarray(codes), np.asarray(trait)
    n = len(codes)
    X = np.zeros((n, t)); X[np.arange(n), trait] = 1.0
    y = X @ rng.normal(size=t) + u[codes, trait] + rng.normal(size=n)
    terms = [terme_long("gen", codes, trait, q, t, struct="us")]
    res = residuelle("diag", trait, np.arange(n), t)
    return terms, res, y, X


def _pedigree_A(n_found=8, n_off=32, rng=None):
    """Additive relationship matrix of a random pedigree."""
    q = n_found + n_off
    sire = -np.ones(q, dtype=int); dam = -np.ones(q, dtype=int)
    for i in range(n_found, q):
        sire[i], dam[i] = rng.choice(i, 2, replace=False)
    A = np.zeros((q, q))
    for i in range(q):
        s, d = sire[i], dam[i]
        A[i, i] = 1.0 + (0.5 * A[s, d] if s >= 0 and d >= 0 else 0.0)
        for j in range(i):
            v = 0.5 * ((A[j, s] if s >= 0 else 0.0) + (A[j, d] if d >= 0 else 0.0))
            A[i, j] = A[j, i] = v
    return A


def cas_pedigree(graine=5, rep=3):
    """A term given only by its precision Kinv (MME route only)."""
    rng = np.random.default_rng(graine)
    A = _pedigree_A(rng=rng)
    q = A.shape[0]
    codes = np.repeat(np.arange(q), rep)
    n = len(codes)
    Ai = np.linalg.inv(A)
    tm = terme_facteur("anim", codes, q)
    tm.update(Kinv=Ai, Kinv_logdet=float(-np.linalg.slogdet(A)[1]), LK=None)
    terms = [tm, terme_facteur("bloc", rng.integers(0, 4, n), 4)]
    X = np.ones((n, 1))
    y = rng.normal(size=n) + 0.8 * np.linalg.cholesky(A)[codes] @ rng.normal(size=q)
    res = residuelle("iid", np.zeros(n), np.arange(n), 1)
    return terms, res, y, X


def cas_noyau_exp(graine=6, q=30, rep=3):
    """Parametrised dense kernel between levels (exp on 1D coordinates)."""
    rng = np.random.default_rng(graine)
    coord = np.sort(rng.uniform(0, 20, q))
    codes = np.repeat(np.arange(q), rep)
    n = len(codes)
    terms = [terme_facteur("champ", codes, q, lvl="exp", coord=coord)]
    X = np.ones((n, 1))
    y = rng.normal(size=n)
    res = residuelle("iid", np.zeros(n), np.arange(n), 1)
    return terms, res, y, X


def cas_grille(graine=7, nr=8, nc=9, q=12, p_vide=0.15):
    """Residual iid x ar1 x ar1 on a grid with missing cells (MME grid route)."""
    rng = np.random.default_rng(graine)
    cells = np.flatnonzero(rng.random(nr * nc) >= p_vide)
    n = len(cells)
    terms = [terme_facteur("g", rng.integers(0, q, n), q)]
    res = residuelle("iid", np.zeros(n), cells, 1, lvl="ar1ar1", dims=(nr, nc))
    res["n_unit"] = nr * nc
    X = np.ones((n, 1))
    y = rng.normal(size=n)
    return terms, res, y, X


# name -> (builder, separator for the hybrid inversion or None)
CASES = {
    "ige": (cas_ige, ["gen_ble", "gen_luz"]),
    "ar1": (cas_ar1, None),
    "sections": (cas_sections_bloc, None),
    "essais": (cas_essais, ["geno", "geno_ar1", "rang_e0"]),
    "rang_2": (cas_rang_2, None),
    "pedigree": (cas_pedigree, None),
    "noyau_exp": (cas_noyau_exp, None),
    "grille": (cas_grille, None),
}
