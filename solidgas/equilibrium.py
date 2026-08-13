"""Constrained Gibbs-energy minimisation for the C-H-O-Ti system.

The unknown is the mole vector `n` over SPECIES. Total Gibbs energy is
minimised subject to elemental balance, with gases treated as ideal and
mixing, and the condensed phases as pure (unit activity, no mixing term).
"""

import numpy as np
from scipy.optimize import minimize

from .shomate import mu0, R_KJ, R_ATM  # noqa: F401  (R_ATM re-exported for scripts)

GASES = ['CO2', 'H2', 'CO', 'H2O', 'CH4']

# Titanium phases, name -> (Ti atoms, Ti(3+) atoms) per formula unit.
#
# Every Magneli member Ti_n O_(2n-1) carries exactly two Ti(3+), for any n:
#   oxygen charge  = -2(2n-1) = -4n+2, so the Ti must supply +4n-2
#   3x + 4(n - x)  = 4n - x = 4n - 2   ->   x = 2
# so Ti2O3, Ti3O5, Ti4O7 and Ti5O9 all sit at two Ti(3+) per formula.
TI_PHASES = {
    'TiO2':  (1, 0),   # rutile, all Ti(4+)
    'Ti5O9': (5, 2),
    'Ti4O7': (4, 2),
    'Ti3O5': (3, 2),
    'Ti2O3': (2, 2),
}

# Which titanium phases take part. Extend this list (and add the matching
# Shomate entry) to bring another member of the series into the calculation.
ACTIVE_TI_PHASES = ['TiO2', 'Ti4O7', 'Ti3O5', 'Ti2O3']

FORMULAS = {
    'CO2': {'C': 1, 'O': 2},
    'H2': {'H': 2},
    'CO': {'C': 1, 'O': 1},
    'H2O': {'H': 2, 'O': 1},
    'CH4': {'C': 1, 'H': 4},
    'TiO2': {'O': 2, 'Ti': 1},
    'Ti5O9': {'O': 9, 'Ti': 5},
    'Ti4O7': {'O': 7, 'Ti': 4},
    'Ti3O5': {'O': 5, 'Ti': 3},
    'Ti2O3': {'O': 3, 'Ti': 2},
}

SPECIES = GASES + ACTIVE_TI_PHASES
ELEMENTS = ['C', 'H', 'O', 'Ti']

# Built from FORMULAS rather than typed out, so adding a phase cannot put the
# element matrix out of step with the species list.
ELEM = np.array([[FORMULAS[s].get(e, 0) for e in ELEMENTS] for s in SPECIES],
                dtype=float)

SOLIDS = set(ACTIVE_TI_PHASES)
gas_idx = [i for i, s in enumerate(SPECIES) if s not in SOLIDS]
sol_idx = [i for i, s in enumerate(SPECIES) if s in SOLIDS]


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


def phase_seeds(n_Ti_total, extents=(0.3, 0.9)):
    """Solid-phase starting vectors, one family per titanium phase.

    The reduced phases sit in separate basins of the Gibbs surface, so a sweep
    seeded only near rutile will never discover them. This walks each phase
    from lightly to nearly fully converted.
    """
    out = []
    for i, p in enumerate(ACTIVE_TI_PHASES):
        for fr in ((0.0,) if p == 'TiO2' else extents):
            v = np.zeros(len(ACTIVE_TI_PHASES))
            v[i] = n_Ti_total * fr / TI_PHASES[p][0]
            v[0] += n_Ti_total * (1 - fr)
            out.append(np.maximum(v, 1e-30))
    return out


def residual(n, b):
    """Largest absolute elemental-balance violation."""
    return max(abs(ELEM[:, j].dot(n) - b[j]) for j in range(ELEM.shape[1]))


def gas_fractions(n, as_percent=True):
    """Gas-phase mole fractions, indexed like SPECIES (solids come back as 0)."""
    ng = sum(n[i] for i in gas_idx)
    scale = 100.0 if as_percent else 1.0
    return [n[i] / ng * scale if i in gas_idx else 0.0 for i in range(len(SPECIES))]


def ti3_percent(n, n_Ti_total):
    """Percentage of total Ti present as Ti(3+), summed over all phases."""
    if n_Ti_total <= 0:
        return 0.0
    ti3 = sum(n[SPECIES.index(p)] * TI_PHASES[p][1] for p in ACTIVE_TI_PHASES)
    return min(ti3 / n_Ti_total, 1.0) * 100.0


def phase_split(n, n_Ti_total):
    """Fraction of total Ti sitting in each titanium phase, as percentages."""
    if n_Ti_total <= 0:
        return {p: 0.0 for p in ACTIVE_TI_PHASES}
    return {p: n[SPECIES.index(p)] * TI_PHASES[p][0] / n_Ti_total * 100.0
            for p in ACTIVE_TI_PHASES}


def mean_valence(n, n_Ti_total):
    """Average Ti oxidation state across the solid."""
    if n_Ti_total <= 0:
        return float('nan')
    ti3 = sum(n[SPECIES.index(p)] * TI_PHASES[p][1] for p in ACTIVE_TI_PHASES)
    return 4.0 - ti3 / n_Ti_total


def moles_of_gas(P_tot, V_cm3, T):
    """Ideal-gas moles filling V_cm3 at P_tot atm and T kelvin."""
    return P_tot * V_cm3 / (R_ATM * T)
