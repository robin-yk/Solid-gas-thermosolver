"""Bulk vacancy aggregates of any size, by transition-matrix Monte Carlo.

Model (assumptions agreed for this pilot)
    Energy is pairwise additive: every vacancy pair within the cutoff adds the
    ZHA2017 bulk pair energy E(r) = A/r - B/r^2 - C/r^6. Clusters of three or
    more therefore have energies (sums of pairs), with no extra many-body term.
    Every vacancy holds two Ti3+ among its three Ti neighbours, one electron per
    Ti (SHI2012: double occupation unstable). This caps the local vacancy
    fraction at 1/4 (Ti2O3 composition), the same cap as the LOCAL closure.

Output: ln Q(N) for a periodic box of n rutile O sites, N = 0..n/4 vacancies.
The box becomes one finite-state domain of the equilibrium engine: state N has
free energy -kT ln Q(N) and N vacancies. Its grand partition function is
therefore exact for the box, clusters of every size included. Boxes do not
interact (finite-size approximation, declared).
"""
import ctypes
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from . import particle as pt

HERE = Path(__file__).resolve().parent
CACHE = HERE.parent / 'outputs' / 'mc'
BOX = (5, 5, 7)                    # a, a, c repeats: 700 O, 350 Ti (half box >= 1 nm)


def _lib():
    src = HERE / 'mc_bulk.c'
    so = HERE / '_build' / ('mc_bulk_' + hashlib.sha1(src.read_bytes()).hexdigest()[:12] + '.so')
    if not so.exists():
        so.parent.mkdir(exist_ok=True)
        subprocess.run(['gcc', '-O2', '-shared', '-fPIC', '-o', str(so), str(src), '-lm'], check=True)
    lib = ctypes.CDLL(str(so))
    P = np.ctypeslib.ndpointer
    common = [ctypes.c_int, ctypes.c_int, P(np.int32), P(np.int32), P(np.float64),
              P(np.int32), P(np.int32), P(np.int32), P(np.int32), P(np.float64), ctypes.c_double]
    lib.canon_run.argtypes = common + [ctypes.c_long, ctypes.c_int, ctypes.c_uint64,
                                       ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_long)]
    lib.grow.argtypes = common + [ctypes.c_uint64]
    return lib


def lattice(box=BOX, cutoff_nm=1.0):
    """O and Ti positions of a periodic rutile box; O-O pairs within the cutoff
    with their pair energy; the three Ti of every O."""
    from .build import zha2017
    nx, ny, nz = box
    L = np.array([nx * pt.A_NM, ny * pt.A_NM, nz * pt.C_NM])
    if cutoff_nm > 0.5 * L.min():
        raise ValueError('box too small for the cutoff')
    u = pt.U
    o_frac = np.array([[u, u, 0], [-u, -u, 0], [.5 + u, .5 - u, .5], [.5 - u, .5 + u, .5]])
    ti_frac = np.array([[0, 0, 0], [.5, .5, .5]])
    cells = np.array([[i, j, k] for i in range(nx) for j in range(ny) for k in range(nz)])
    cellv = np.array([pt.A_NM, pt.A_NM, pt.C_NM])
    O = ((cells[:, None, :] + o_frac[None]) * cellv).reshape(-1, 3) % L
    Ti = ((cells[:, None, :] + ti_frac[None]) * cellv).reshape(-1, 3) % L

    def mind(d):
        return d - L * np.round(d / L)
    nO = len(O)
    start, idx, J = [0], [], []
    for i in range(nO):
        d = np.sqrt((mind(O - O[i]) ** 2).sum(1))
        sel = np.where((d > 1e-9) & (d <= cutoff_nm + 1e-12))[0]
        idx += sel.tolist(); J += [zha2017(r) for r in d[sel]]
        start.append(len(idx))
    oti = []
    for i in range(nO):
        d = np.sqrt((mind(Ti - O[i]) ** 2).sum(1))
        near = np.where(d < 0.21)[0]
        if len(near) != 3:
            raise RuntimeError('each O must have three Ti neighbours')
        oti += near.tolist()
    return dict(nO=nO, nTi=len(Ti), start=np.array(start, np.int32), idx=np.array(idx, np.int32),
                J=np.array(J, float), oti=np.array(oti, np.int32))


def ln_q(T, cutoff_nm=1.0, box=BOX, moves_per_n=20000, sample_every=50, seed=12345):
    """ln Q(N), N = 0..n/4, for the box at temperature T (cached on disk).

    Q(N+1)/Q(N) = <W>_N / (N+1), W = sum over every empty site that can take a
    vacancy of exp(-h_i/kT) (Widom insertion, summed exactly over all sites,
    averaged over canonical samples at N). The chain then grows by one vacancy
    drawn with those same weights, and continues at N+1.
    """
    key = dict(T=T, cutoff=cutoff_nm, box=list(box), moves=moves_per_n, every=sample_every, seed=seed,
               src=hashlib.sha1((HERE / 'mc_bulk.c').read_bytes()).hexdigest()[:12])
    tag = hashlib.sha1(json.dumps(key, sort_keys=True).encode()).hexdigest()[:16]
    f = CACHE / f'lnq_{tag}.json'
    if f.exists():
        return np.array(json.loads(f.read_text())['lnQ'])
    lat = lattice(box, cutoff_nm)
    lib = _lib()
    nO, nTi = lat['nO'], lat['nTi']
    beta = 1.0 / (pt.KB * T)
    occ = np.zeros(nO, np.int32); owner = -np.ones(nTi, np.int32); slot = -np.ones(2 * nO, np.int32)
    h = np.zeros(nO)
    args = (nO, nTi, lat['start'], lat['idx'], lat['J'], lat['oti'], occ, owner, slot, h, beta)
    lnQ = [0.0]
    for n in range(nO // 4):
        acc, ns = ctypes.c_double(0.0), ctypes.c_long(0)
        lib.canon_run(*args, moves_per_n if n > 0 else sample_every, sample_every, seed + n,
                      ctypes.byref(acc), ctypes.byref(ns))
        if ns.value == 0 or not np.isfinite(acc.value):
            break
        lnQ.append(lnQ[-1] + acc.value - np.log(ns.value) - np.log(n + 1))
        if lib.grow(*args, seed + 7919 * (n + 1)) != n + 1:
            break
    lnQ = np.array(lnQ)
    CACHE.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(dict(key=key, lnQ=lnQ.tolist(), n_sites=nO)))
    return lnQ
