"""Constrained Gibbs-energy minimisation for the C-H-O-Ti system.

The unknown is the mole vector `n` over SPECIES. Total Gibbs energy is
minimised subject to elemental balance, with gases treated as ideal and
mixing, and the condensed phases as pure (unit activity, no mixing term).
"""

import numpy as np
from scipy.optimize import minimize

from .shomate import mu0, R_KJ, R_ATM  # noqa: F401  (R_ATM re-exported for scripts)

SPECIES = ['CO2', 'H2', 'CO', 'H2O', 'CH4', 'TiO2', 'Ti4O7']

ELEMENTS = ['C', 'H', 'O', 'Ti']

#                C  H  O  Ti
ELEM = np.array([[1, 0, 2, 0],   # CO2
                 [0, 2, 0, 0],   # H2
                 [1, 0, 1, 0],   # CO
                 [0, 2, 1, 0],   # H2O
                 [1, 4, 0, 0],   # CH4
                 [0, 0, 2, 1],   # TiO2  -> 2 O, 1 Ti
                 [0, 0, 7, 4]],  # Ti4O7 -> 7 O, 4 Ti
                dtype=float)

SOLIDS = {'TiO2', 'Ti4O7'}
gas_idx = [i for i, s in enumerate(SPECIES) if s not in SOLIDS]
sol_idx = [i for i, s in enumerate(SPECIES) if s in SOLIDS]

# Ti4O7 averages Ti(3.5+): 2 of its 4 Ti are Ti(3+), the other 2 Ti(4+).
TI3_PER_TI4O7 = 2.0


def G_total(n, T, P_tot=1.0):
    """Total Gibbs energy in kJ for a mole vector `n` at temperature T (K)."""
    ng = sum(n[i] for i in gas_idx)
    if ng <= 0:
        return 1e10
    return (sum(n[i] * (mu0(SPECIES[i], T) + R_KJ * T * np.log(n[i] / ng * P_tot))
                for i in gas_idx if n[i] > 1e-30)
            + sum(n[i] * mu0(SPECIES[i], T) for i in sol_idx))


def solve(b, T, inits, P_tot=1.0, tol=1e-8):
    """Minimise G subject to ELEM.T @ n == b.

    `inits` is a list of starting mole vectors; SLSQP is run from each and the
    lowest feasible optimum wins. Multiple starts matter here because the
    solid reduction branch and the all-gas branch are separate local minima.
    Returns the mole vector, or None if no start converged feasibly.
    """
    cons = [{'type': 'eq', 'fun': lambda n, j=j: ELEM[:, j].dot(n) - b[j]}
            for j in range(ELEM.shape[1])]
    bnds = [(1e-30, None)] * len(SPECIES)
    best, bestG = None, 1e10
    for n_init in inits:
        try:
            s = minimize(lambda n: G_total(n, T, P_tot), n_init, method='SLSQP',
                         bounds=bnds, constraints=cons,
                         options={'ftol': 1e-13, 'maxiter': 1000, 'disp': False})
            if s.fun < bestG and residual(s.x, b) < tol:
                bestG, best = s.fun, s
        except Exception:
            pass
    return np.maximum(best.x, 0) if best is not None else None


def residual(n, b):
    """Largest absolute elemental-balance violation."""
    return max(abs(ELEM[:, j].dot(n) - b[j]) for j in range(ELEM.shape[1]))


def gas_fractions(n, as_percent=True):
    """Gas-phase mole fractions, indexed like SPECIES (solids come back as 0)."""
    ng = sum(n[i] for i in gas_idx)
    scale = 100.0 if as_percent else 1.0
    return [n[i] / ng * scale if i in gas_idx else 0.0 for i in range(len(SPECIES))]


def ti3_percent(n, n_Ti_total):
    """Percentage of total Ti present as Ti(3+), carried by Ti4O7."""
    if n_Ti_total <= 0:
        return 0.0
    return min(TI3_PER_TI4O7 * n[SPECIES.index('Ti4O7')] / n_Ti_total, 1.0) * 100.0


def moles_of_gas(P_tot, V_cm3, T):
    """Ideal-gas moles filling V_cm3 at P_tot atm and T kelvin."""
    return P_tot * V_cm3 / (R_ATM * T)
