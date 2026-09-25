"""Constrained free-energy minimum over a finite state table.

A domain d holds C_d umol of cells per gram; each cell is in one state a with
energy E_a, degeneracy g_a, v_a oxygen vacancies and charge q_a = 2 v_a - e_a
(e_a = Ti3+ electrons). Populations x_a = C_d p_a minimise

    F = sum x_a E_a + kT sum x_a ln(x_a / (C_d g_a))

subject to sum v_a x_a = N (measured inventory) and, for every charge region r,
sum_{a in r} q_a x_a = 0 (local neutrality). Domains outside every region
carry implicit electrons (neutral closure) and q = 0. The minimiser is Gibbs form,

    p_a = g_a exp(-E_a/kT + mu v_a + chi_r q_a) / Z_d,

so the unknowns are only mu and one chi per region. Each region's charge is
increasing in its own chi (slope = capacity-weighted variance of q), and the
inventory is increasing in mu once every region is neutral (slope = the Schur
complement of that covariance). Both are therefore bracketed scalar roots: an
inner Newton-bisection for all chi at once, an outer one for mu. Neither can
fail to converge, however far the start is from the answer.
"""
import numpy as np

from .particle import KB


class Model:
    """domains: list of dicts with C, E, v and optional e, g, region, tag."""

    def __init__(self, domains):
        S = max(len(d['E']) for d in domains)
        D = len(domains)
        self.C = np.array([d['C'] for d in domains], float)
        self.E = np.zeros((D, S)); self.v = np.zeros((D, S)); self.q = np.zeros((D, S))
        self.logg = np.full((D, S), -np.inf)
        for i, d in enumerate(domains):
            n = len(d['E'])
            v = np.asarray(d['v'], float); e = np.asarray(d.get('e', [0] * n), float)
            self.E[i, :n] = d['E']; self.v[i, :n] = v; self.q[i, :n] = 2 * v - e
            self.logg[i, :n] = np.log(np.asarray(d.get('g', [1] * n), float))
        if np.any(self.C <= 0):
            raise ValueError('every domain needs positive capacity')
        # Outside any region the electrons are implicit (neutral closure).
        region = [d.get('region', -1) for d in domains]
        self.q[np.array(region) < 0] = 0
        ids = {r: k for k, r in enumerate(sorted({r for r in region if r >= 0}))}
        self.region = np.array([ids.get(r, -1) for r in region])
        self.nreg = len(ids)
        self.tags = [d.get('tag') for d in domains]

    def _p(self, y_mu, chi, beta):
        L = self.logg - beta * self.E + y_mu * self.v
        if self.nreg:
            L = L + np.r_[0.0, chi][self.region + 1][:, None] * self.q
        w = np.exp(L - L.max(axis=1)[:, None])
        return w / w.sum(axis=1)[:, None]

    def _neutralise(self, y_mu, chi, beta, tol, N):
        """chi that makes every region neutral at this mu, with d(inventory)/d(mu)."""
        R, inr = self.nreg, self.region >= 0
        ri, C = self.region[inr], self.C[inr]
        lo, hi = np.full(R, -BOUND), np.full(R, BOUND)
        for _ in range(MAX_ITER):
            p = self._p(y_mu, chi, beta)[inr]
            mq = (p * self.q[inr]).sum(1)
            h = np.bincount(ri, C * mq, R)
            dh = np.bincount(ri, C * ((p * self.q[inr] ** 2).sum(1) - mq ** 2), R)
            if np.all((np.abs(h) <= tol * N) | (hi - lo <= COLLAPSE * np.maximum(1, np.abs(chi)))):
                break
            hi = np.where(h > 0, chi, hi); lo = np.where(h > 0, lo, chi)
            step = chi - h / np.maximum(dh, 1e-300)
            chi = np.where((step > lo) & (step < hi), step, 0.5 * (lo + hi))
        else:
            raise RuntimeError('local neutrality did not converge')
        mv = (p * self.v[inr]).sum(1)
        b = np.bincount(ri, C * ((p * self.v[inr] * self.q[inr]).sum(1) - mv * mq), R)
        return chi, float(np.sum(b * b / dh))

    def solve(self, N, T, tol=1e-12):
        beta = 1.0 / (KB * T)
        if not 0 < N < float(np.sum(self.C * self.v.max(axis=1))):
            raise ValueError('inventory outside the represented capacity')
        y, chi, lo, hi = 0.0, np.zeros(self.nreg), -BOUND, BOUND
        for it in range(1, MAX_ITER + 1):
            schur = 0.0
            if self.nreg:
                chi, schur = self._neutralise(y, chi, beta, tol, N)
            p = self._p(y, chi, beta)
            mv = (p * self.v).sum(1)
            g = float(self.C @ mv) - N
            # Stop at the tolerance, or when the bracket has shrunk to rounding.
            if abs(g) <= tol * N or hi - lo <= COLLAPSE * max(1, abs(y)):
                break
            dg = float(self.C @ ((p * self.v ** 2).sum(1) - mv ** 2)) - schur
            hi, lo = (y, lo) if g > 0 else (hi, y)
            step = y - g / max(dg, 1e-300)
            y = step if lo < step < hi else 0.5 * (lo + hi)
        else:
            raise RuntimeError('equilibrium did not converge')
        x = self.C[:, None] * p
        inr = self.region >= 0
        chg = np.bincount(self.region[inr], (x * self.q).sum(1)[inr], self.nreg)
        return Solution(self, x, y * KB * T, chi, it, float((x * self.v).sum() - N),
                        float(np.abs(chg).max(initial=0.0)))


BOUND, MAX_ITER, COLLAPSE = 300.0, 400, 1e-14


class Solution:
    def __init__(self, model, x, mu_eV, chi, iterations, inv_res, charge_res):
        self.model, self.x, self.mu_eV, self.chi = model, x, mu_eV, chi
        self.iterations, self.inventory_residual, self.charge_residual = iterations, inv_res, charge_res

    def vacancies(self, tag):
        """Vacancy population (umol/g) summed over domains with this tag."""
        m = self.model
        return float(sum((self.x[i] * m.v[i]).sum() for i, t in enumerate(m.tags) if t == tag))
