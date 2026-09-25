"""Fixed-inventory neutral vacancies with local Ti3+ compensation (Note 2b)."""
import numpy as np
from scipy.optimize import brentq

MOLAR = 79.866
KB = 8.617333262145e-5

def solve(inventory=94., radius=450., temperature=873.15, A=1.31, xi=.5):
    if not all(np.isfinite(v) for v in [inventory, radius, temperature, A, xi]):
        raise ValueError('Inputs must be finite')
    target = inventory * MOLAR * 1e-6 / 2
    if not (0 < target < .25 and radius > 0 and temperature > 0 and A >= 0 and xi > 0):
        raise ValueError('Require 0 < vacancy fraction < 0.25 and positive radius, temperature, decay length')
    edge = min(radius, max(20 * xi, 2.))
    edges = np.unique(np.r_[np.linspace(0, edge, 2001), np.linspace(edge, radius, 201), min(2., radius)])
    depth = (edges[1:] + edges[:-1]) / 2
    weights = ((radius - edges[:-1])**3 - (radius - edges[1:])**3) / radius**3
    energy = -A * np.exp(-depth / xi)
    def occupancy(mu, E):
        lo = np.zeros_like(E, dtype=float); hi = np.full_like(E, .25, dtype=float)
        for _ in range(56):
            x = (lo + hi) / 2
            with np.errstate(divide='ignore'):
                h = np.log(x / (1-x)) + 2*np.log(4*x / (1-4*x))
            mask = E + KB*temperature*h < mu
            lo = np.where(mask, x, lo); hi = np.where(mask, hi, x)
        return (lo+hi)/2
    mu = brentq(lambda m: weights @ occupancy(m, energy)-target, -A-100*KB*temperature, 100*KB*temperature, xtol=1e-14)
    x = occupancy(mu, energy)
    return dict(surface=float(occupancy(mu, np.array(-A))), interior=float(occupancy(mu, np.array(-A*np.exp(-radius/xi)))), mean=target,
                shell=float(weights[depth < min(2., radius)] @ x[depth < min(2., radius)]/target),
                closure=float(weights @ x-target), mu=mu)
