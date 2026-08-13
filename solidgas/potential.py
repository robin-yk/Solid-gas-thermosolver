"""Oxygen-potential route to the same answer, without the solids as unknowns.

Titanium does not enter the gas below about 2000 K (TiO(g) sits at
dHf = +54.4 kJ/mol), so the solid never exchanges anything with the gas but
oxygen. That splits the problem in two:

  1. equilibrate the C-H-O gas on its own, and read off its oxygen potential
  2. compare that against each condensed phase's critical oxygen potential,
     which is a plain function of temperature computed straight from the
     Shomate data - no optimisation at all

The line compounds are then present or absent rather than continuous unknowns,
which is what they physically are (unit activity). That removes the whole
reason the full minimisation needed multiple starts: the amount of an absent
phase is a variable whose optimum is on its own lower bound, and SLSQP has to
crawl there through a log singularity that never had to exist.

mu_O is quoted per mole of atomic oxygen, referenced to 1/2 O2, in kJ.
"""

import numpy as np
from scipy.optimize import minimize

from .shomate import mu0, R_KJ
from .equilibrium import ACTIVE_TI_PHASES, FORMULAS, GASES

# Redox couples that can report the gas oxygen potential. Both must agree at
# equilibrium; the disagreement is a useful check on the gas solution.
COUPLES = [('H2', 'H2O'), ('CO', 'CO2')]


def mu_O_from_gas(y, T):
    """Oxygen potential of a gas mixture, from whichever couples are present.

    `y` maps species name to mole fraction. Returns (mu_O, spread) where
    spread is the disagreement between couples, zero at a converged gas.
    """
    vals = []
    for red, ox in COUPLES:
        a, b = y.get(red, 0.0), y.get(ox, 0.0)
        if a > 1e-25 and b > 1e-25:
            vals.append(mu0(ox, T) - mu0(red, T) + R_KJ * T * np.log(b / a))
    if not vals:
        return float('nan'), float('nan')
    return float(np.mean(vals)), float(max(vals) - min(vals))


def mu_O_critical(lower, higher, T):
    """Critical oxygen potential for `higher` -> `lower`, per mole of O removed.

    `higher` is the more oxidised phase. Below this potential the reduced phase
    is the stable one. Both are line compounds, so this is a pure function of T.
    """
    fh, fl = FORMULAS[higher], FORMULAS[lower]
    # scale so both sides carry the same titanium
    n_hi = fl['Ti'] / fh['Ti']
    dO = n_hi * fh['O'] - fl['O']
    if dO <= 0:
        raise ValueError(f'{lower} is not more reduced than {higher}')
    return (n_hi * mu0(higher, T) - mu0(lower, T)) / dO


def phase_ladder(T, phases=None):
    """Critical potentials down the series, most oxidised first."""
    ph = list(phases or ACTIVE_TI_PHASES)
    ph.sort(key=lambda p: -FORMULAS[p]['O'] / FORMULAS[p]['Ti'])
    return [(ph[i], ph[i + 1], mu_O_critical(ph[i + 1], ph[i], T))
            for i in range(len(ph) - 1)]


def stable_phase(mu_O, T, phases=None):
    """Which condensed phase a given oxygen potential holds, and by how much.

    Returns (phase, margin). Margin is how far mu_O sits above the potential
    that would take it to the next phase down, in kJ per mole of O: positive
    means that step does not happen, and the size says by how much.
    """
    ladder = phase_ladder(T, phases)
    for higher, lower, crit in ladder:
        if mu_O > crit:
            return higher, mu_O - crit
    return ladder[-1][1], float('nan')


def gas_equilibrium(b_CHO, T, P_tot=1.0, species=None):
    """Equilibrate a C-H-O gas on its own. Small, convex, converges cold.

    `b_CHO` is the (C, H, O) element vector. Returns mole numbers keyed by
    species name.
    """
    sp = list(species or GASES)
    E = np.array([[FORMULAS[s].get(e, 0) for e in ('C', 'H', 'O')] for s in sp],
                 dtype=float)

    # Solve normalised. Mole fractions do not depend on how much gas there is
    # at fixed pressure, and a charge of ~1e-4 mol leaves SLSQP working on
    # equality constraints of that magnitude, where it stops about three digits
    # early. Scaling the basis to order one and multiplying back afterwards is
    # exact and buys roughly three orders of magnitude on the quotients.
    scale = float(np.sum(b_CHO))
    if scale <= 0:
        return dict.fromkeys(sp, 0.0)
    b = np.asarray(b_CHO, dtype=float) / scale

    def G(n):
        ng = n.sum()
        if ng <= 0:
            return 1e10
        return sum(n[i] * (mu0(sp[i], T) + R_KJ * T * np.log(max(n[i], 1e-300) / ng * P_tot))
                   for i in range(len(sp)) if n[i] > 1e-30)

    cons = [{'type': 'eq', 'fun': lambda n, j=j: E[:, j].dot(n) - b[j]}
            for j in range(3)]
    start = np.full(len(sp), b.sum() / len(sp))
    r = minimize(G, start, method='SLSQP', bounds=[(1e-30, None)] * len(sp),
                 constraints=cons, options={'ftol': 1e-14, 'maxiter': 500})
    n = np.maximum(r.x, 0.0) * scale
    return dict(zip(sp, n))


def survey(b_CHO, T, P_tot=1.0, phases=None):
    """Full answer for one temperature: gas, its oxygen potential, the phase."""
    n = gas_equilibrium(b_CHO, T, P_tot)
    tot = sum(n.values())
    y = {k: v / tot for k, v in n.items()}
    mu_O, spread = mu_O_from_gas(y, T)
    phase, margin = stable_phase(mu_O, T, phases)
    return {'y': y, 'mu_O': mu_O, 'couple_spread': spread,
            'phase': phase, 'margin': margin}
