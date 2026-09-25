"""Vacancy transport on the same free energy: how fast the 600 C state is reached.

Every O domain exchanges vacancies with its neighbours. For one link i-j with
N_ij site pairs (umol/g) and barrier B above the higher endpoint energy,

    F = N Gamma0 exp(-(B + max(eps_i, eps_j) - eps_i)/kT) theta_i (1 - theta_j)
    R = the same with i and j exchanged
    J = Lambda(F, R) * (mu_i - mu_j)/kT,   Lambda = logarithmic mean.

Gamma0 = kT/h (the v12 transport-gate prefactor). mu is the equilibrium
chemical potential of a vacancy on that site, including the Ti3+ term for
LOCAL. For NEUTRAL, (mu_i - mu_j)/kT = ln(F/R), so J = F - R: exact
site-exclusion hopping. For LOCAL the same flux carries the extra
electron-entropy force, which gives the ambipolar factor 3 in the dilute
limit. J (mu_i - mu_j) >= 0 on every link, so the free energy falls
monotonically and the only rest point is the equilibrium solver's answer.

Links
    Bulk and continuum shells: finite-volume faces, N = (z lambda^2 / 6)
    A n_O / dr per particle, so that dilute flat-energy flux is Fick's law with
    D = z lambda^2 Gamma / 6 (Iddir step C geometry, particle.bulk_hop).
    Explicit discrete sites: the WU (Wu 2023) or JUG (Jug 2005) near-surface
    edges, with atomic O layer n -> trilayer (n-1)//3 + 1, site BRI/IPL/SBR by
    (n-1) % 3. Consecutive sites without a literature edge use the bulk
    barrier. One partner per site, N = min(C_i, C_j) (declared).
"""
import csv
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

from . import particle as pt
from .build import A_LI, HERE

H_PLANCK = 4.135667696e-15          # eV s
SITES = ('BRI', 'IPL', 'SBR')
BULK_SETS = ('IDDIR_2007', 'TSETSERIS_2011', 'BAKULIN_2021', 'KNAUP_2012', 'UBERUAGA_2011')


def _edges():
    """specification/transport_edges.csv: near-surface networks and bulk barriers.

    Atomic O layer n is trilayer (n-1)//3 + 1, site BRI/IPL/SBR by (n-1) % 3.
    A pair given in both directions keeps the lower barrier, the one from the
    higher-energy end in the paper's own energetics.
    """
    surface, bulk = {}, {}
    with open(HERE / 'specification' / 'transport_edges.csv') as f:
        for r in csv.DictReader(f):
            a, b, B = int(r['from_atomic_O_layer']), int(r['to_atomic_O_layer']), float(r['barrier_eV'])
            if r['network'] in BULK_SETS:
                bulk[r['network']] = B
            elif a != b and a > 0:
                site = lambda n: ((n - 1) // 3 + 1, SITES[(n - 1) % 3])   # noqa: E731
                net = r['network'].split('_')[0]
                key = tuple(sorted((site(a), site(b))))
                edges = surface.setdefault(net, {})
                edges[key] = min(B, edges.get(key, B))
    return surface, bulk


_SURFACE, _BULK = _edges()
SURFACE_SETS = {net: [(a, b, B) for (a, b), B in e.items()] for net, e in _SURFACE.items()}
BULK_BARRIERS = {min(_BULK, key=_BULK.get): min(_BULK.values()),
                 max(_BULK, key=_BULK.get): max(_BULK.values())}


def _face_contact(R, z_face, dr, mass):
    lam, z = pt.bulk_hop()
    n_o = 2 * pt.RHO / pt.M * pt.NA * 1e-21
    area = 4 * math.pi * (R - z_face) ** 2
    return z * lam ** 2 / 6 * area * n_o / dr / (mass * pt.NA) * 1e6


class Network:
    def __init__(self, model, lay, bulk_barrier, surface_set='WU'):
        self.model, self.lay = model, lay
        o = sorted(lay.o, key=lambda r: r['z'])
        self.idx = np.array([r['idx'] for r in o])
        self.C = np.array([r['C'] for r in o])
        self.eps = np.array([r['eps'] for r in o])
        self.z = np.array([r['z'] for r in o])
        self.local = lay.closure == 'LOCAL'
        if lay.closure not in ('NEUTRAL', 'LOCAL') or lay.pairs:
            raise ValueError('transport is defined for NEUTRAL/LOCAL without pairs')
        n = len(o)
        if self.local:
            regs = sorted({r['region'] for r in o})
            pos = {r: k for k, r in enumerate(regs)}
            self.reg = np.array([pos[r['region']] for r in o])
            self.cti = np.array([lay.ti[r][1] for r in regs])
        mass = pt.particle_mass(lay.d_nm)
        B = bulk_barrier
        links = {}
        if lay.kind == 'discrete':
            tags = [model.tags[i] for i in self.idx]
            where = {t: k for k, t in enumerate(tags) if isinstance(t, tuple)}
            for a, b, barrier in SURFACE_SETS[surface_set]:
                i, j = sorted((where[a], where[b]))
                links[(i, j)] = (barrier, min(self.C[i], self.C[j]))
            n_expl = len(where)
            for i in range(n_expl - 1):
                links.setdefault((i, i + 1), (B, min(self.C[i], self.C[i + 1])))
            top = lay.edges[0]
            links[(n_expl - 1, n_expl)] = (B, _face_contact(lay.R, top, self.z[n_expl] - self.z[n_expl - 1], mass))
            for i in range(n_expl, n - 1):
                links[(i, i + 1)] = (B, _face_contact(lay.R, lay.edges[i - n_expl + 1],
                                                      self.z[i + 1] - self.z[i], mass))
        else:
            for i in range(n - 1):
                links[(i, i + 1)] = (B, _face_contact(lay.R, lay.edges[i + 1], self.z[i + 1] - self.z[i], mass))
        pairs = sorted(links)
        self.li = np.array([p[0] for p in pairs]); self.lj = np.array([p[1] for p in pairs])
        self.lB = np.array([links[p][0] for p in pairs]); self.lN = np.array([links[p][1] for p in pairs])
        S = lil_matrix((n, n))
        groups = [np.where(self.reg == self.reg[k])[0] if self.local else [k] for k in range(n)]
        for i, j in pairs:
            for a in (i, j):
                for b in list(groups[i]) + list(groups[j]):
                    S[a, b] = 1
        for k in range(n):
            S[k, k] = 1
        self.sparsity = S.tocsr()

    def beta_mu(self, th, beta):
        th = np.clip(th, 1e-300, 1 - 1e-16)
        bm = beta * self.eps + np.log(th) - np.log1p(-th)
        if self.local:
            y = 2 * np.bincount(self.reg, self.C * th, len(self.cti)) / self.cti
            y = np.clip(y, 1e-300, 1 - 1e-16)
            bm = bm + 2 * (np.log(y) - np.log1p(-y))[self.reg]
        return bm

    def rhs(self, t, th, beta, g0):
        th = np.clip(th, 1e-300, 1 - 1e-16)
        bm = self.beta_mu(th, beta)
        i, j = self.li, self.lj
        emax = np.maximum(self.eps[i], self.eps[j])
        k = self.lN * g0 * np.exp(-beta * self.lB)
        F = k * np.exp(-beta * (emax - self.eps[i])) * th[i] * (1 - th[j])
        Rv = k * np.exp(-beta * (emax - self.eps[j])) * th[j] * (1 - th[i])
        ell = np.log(F) - np.log(Rv)
        lam = Rv * np.where(np.abs(ell) < 1e-8, 1 + ell / 2, np.expm1(ell) / np.where(ell == 0, 1, ell))
        J = lam * (bm[i] - bm[j])
        out = np.zeros_like(th)
        np.add.at(out, i, -J); np.add.at(out, j, J)
        return out / self.C

    def theta_of(self, sol):
        return np.array([sol.occupancy(i) for i in self.idx])

    def surface(self, th, beta):
        """Bridging coverage from the node state (same definition as build.surface_coverage)."""
        if self.lay.kind == 'discrete':
            return float(th[list(self.idx).index(self.lay.bri)])
        bm0 = self.beta_mu(th, beta)[0]                       # outermost shell
        G0 = -A_LI
        if not self.local:
            return 1.0 / (1.0 + math.exp(beta * G0 - bm0))
        lo, hi = 0.0, 0.25
        for _ in range(200):
            x = 0.5 * (lo + hi)
            f = beta * G0 + math.log(x / (1 - x)) + 2 * math.log(4 * x / (1 - 4 * x)) - bm0
            lo, hi = (lo, x) if f > 0 else (x, hi)
        return 0.5 * (lo + hi)


def relax(model, lay, N, T_start, T_run, bulk_barrier, surface_set='WU',
          t_obs=(1.0, 10.0, 60.0, 600.0), rtol=1e-7, theta0=None):
    """Start from equilibrium at T_start, hold at T_run. Returns surface coverage
    at each observation time, the equilibrium value at T_run, and the time after
    which the coverage stays within 1 % of it."""
    net = Network(model, lay, bulk_barrier, surface_set)
    beta = 1.0 / (pt.KB * T_run)
    g0 = pt.KB * T_run / H_PLANCK
    th0 = net.theta_of(model.solve(N, T_start)) if theta0 is None else np.asarray(theta0, float)
    th_eq = net.theta_of(model.solve(N, T_run))
    s_eq = net.surface(th_eq, beta)
    t_end = max(t_obs)
    grid = np.unique(np.r_[np.geomspace(1e-9, t_end, 400), t_obs])
    res = solve_ivp(net.rhs, (0.0, t_end), th0, method='BDF', t_eval=grid, args=(beta, g0),
                    jac_sparsity=net.sparsity, rtol=rtol, atol=1e-14 * np.ones_like(th0))
    if not res.success:
        raise RuntimeError(res.message)
    cov = np.array([net.surface(res.y[:, k], beta) for k in range(len(res.t))])
    off = np.abs(cov / s_eq - 1) > 0.01
    t_1pct = float(res.t[np.where(off)[0][-1] + 1]) if off.any() and not off[-1] else (0.0 if not off.any() else math.inf)
    drift = float(abs(net.C @ res.y[:, -1] - net.C @ th0) / N)
    at = {t: float(cov[np.searchsorted(res.t, t)]) for t in t_obs}
    return dict(coverage=at, coverage_eq=s_eq, coverage_start=float(cov[0]), t_1pct=t_1pct,
                mass_drift=drift, theta_end=res.y[:, -1], theta_eq=th_eq, t=res.t, theta=res.y,
                weights=net.C)
