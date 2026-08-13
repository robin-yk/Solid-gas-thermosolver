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

from .shomate import mu0, R_KJ, R_ATM
from .equilibrium import ACTIVE_TI_PHASES, FORMULAS, GASES
from . import waldner

# Which dataset supplies the titanium-bearing condensed phases. NIST-JANAF has
# no Magneli member between rutile and Ti4O7, and those are the ones that
# decide the hot end, so 'waldner' is the default. The two must never be mixed
# within one comparison: they put rutile 6.0 kJ/mol apart and Ti4O7 12.0
# kJ/mol apart, which is the same size as the margins being measured.
SOURCES = {
    'waldner': waldner.mu0,
    'nist': mu0,
}
DEFAULT_SOURCE = 'waldner'


def solid_mu0(name, T, source=DEFAULT_SOURCE):
    return SOURCES[source](name, T)


def ti_phases(source=DEFAULT_SOURCE):
    """Titanium phases the given source actually carries, most oxidised first.

    NIST-JANAF stops at Ti4O7, so asking it for the shallow Magneli members
    returns a shorter list rather than a KeyError. That difference is the whole
    reason the Waldner assessment is here.
    """
    if source == 'waldner':
        names = list(waldner.WALDNER)
    else:
        from .shomate import SHOMATE
        names = [p for p in ACTIVE_TI_PHASES if p in SHOMATE]
    return sorted(names, key=lambda p: -_ratio(p))


def _ratio(name):
    if name in FORMULAS:
        f = FORMULAS[name]
        return f['O'] / f['Ti']
    n, o = waldner._stoich(name)
    return o / n

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


def _stoich(name):
    if name in FORMULAS:
        return FORMULAS[name]['Ti'], FORMULAS[name]['O']
    return waldner._stoich(name)


def mu_O_to_pO2(mu_O, T):
    """Oxygen partial pressure, in atm, that carries a given oxygen potential.

    mu_O = 0.5 * [mu0(O2, T) + R T ln pO2]. The mu0(O2) term is the one that is
    easy to drop, and dropping it moves the answer by fourteen orders of
    magnitude in pressure - though never by anything at all in a margin, since
    margins are differences of two potentials and it cancels.
    """
    return float(np.exp((2.0 * mu_O - mu0('O2', T)) / (R_KJ * T)))


def pO2_to_mu_O(pO2, T):
    """Inverse of mu_O_to_pO2."""
    return 0.5 * (mu0('O2', T) + R_KJ * T * np.log(pO2))


def inert_atmosphere(T, V_cm3, n_Ti_total, pO2_initial=0.0, P_tot=1.0,
                     source=DEFAULT_SOURCE, phases=None):
    """Equilibrium under a gas that carries no oxygen of its own, such as N2.

    With nothing to accept oxygen but the gas phase itself, the oxide gives up
    oxygen as O2 until the potential it sets matches the phase boundary. The
    amount is set by the vessel: a bigger headspace holds more O2 at the same
    pressure and so pulls the solid further.

    `pO2_initial` is the oxygen already in the carrier gas. It matters more
    than the temperature does, which is the real answer to 'what does pure N2
    do' - the number comes off the cylinder's purity spec, not a table.
    """
    ph = [p for p in (phases or ti_phases(source)) if _ratio(p) < 2.0]
    best, best_crit = None, -np.inf
    for p in ph:
        c = mu_O_critical(p, 'TiO2', T, source)
        if c > best_crit:
            best, best_crit = p, c

    pO2_eq = mu_O_to_pO2(best_crit, T)
    n_O2_released = max(pO2_eq - pO2_initial, 0.0) * V_cm3 / (R_ATM * T)

    n_ti, n_o = _stoich(best)
    dO = 2 * n_ti - n_o                      # oxygen atoms freed per formula
    n_phase = 2.0 * n_O2_released / dO       # O2 carries two oxygen atoms
    capped = n_phase * n_ti > n_Ti_total
    if capped:
        n_phase = n_Ti_total / n_ti

    return {'phase': best, 'mu_O_boundary': best_crit, 'pO2_eq': pO2_eq,
            'n_O2': n_O2_released, 'n_phase': n_phase,
            'ti3_percent': min(2.0 * n_phase / n_Ti_total, 1.0) * 100.0,
            'O_removed_fraction': n_phase * dO / (2.0 * n_Ti_total),
            'solid_limited': capped,
            'reducing': pO2_initial < pO2_eq}


def mu_O_critical(lower, higher, T, source=DEFAULT_SOURCE):
    """Critical oxygen potential for `higher` -> `lower`, per mole of O removed.

    `higher` is the more oxidised phase. Below this potential the reduced phase
    is the stable one. Both are line compounds, so this is a pure function of T.

    Note this is an absolute oxygen potential on the elemental reference, not a
    pressure. Converting to pO2 needs mu0(O2) = -T*S(O2), which is not carried
    here; differences of two potentials - which is all any margin is - do not
    need it, since the O2 term cancels.
    """
    th, oh = _stoich(higher)
    tl, ol = _stoich(lower)
    n_hi = tl / th                      # scale so both sides carry the same Ti
    dO = n_hi * oh - ol
    if dO <= 0:
        raise ValueError(f'{lower} is not more reduced than {higher}')
    return (n_hi * solid_mu0(higher, T, source) - solid_mu0(lower, T, source)) / dO


def phase_ladder(T, phases=None, source=DEFAULT_SOURCE):
    """Critical potentials down the series, most oxidised first."""
    ph = list(phases) if phases else ti_phases(source)
    ph.sort(key=lambda p: -_ratio(p))
    return [(ph[i], ph[i + 1], mu_O_critical(ph[i + 1], ph[i], T, source))
            for i in range(len(ph) - 1)]


def first_phase_margin(mu_O, T, source=DEFAULT_SOURCE, phases=None):
    """Margin against whichever reduced phase forms most easily from rutile.

    Every reduced phase is compared straight to rutile rather than down a
    chain, because the question is which one appears first out of untouched
    rutile. Returns (phase, margin) with margin in kJ per mole of oxygen:
    positive means nothing reduces.
    """
    ph = [p for p in (phases or ti_phases(source)) if _ratio(p) < 2.0]
    best, best_margin = None, float('inf')
    for p in ph:
        m = mu_O - mu_O_critical(p, 'TiO2', T, source)
        if m < best_margin:
            best, best_margin = p, m
    return best, best_margin


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
    phase, margin = stable_phase(mu_O, T, phases or ACTIVE_TI_PHASES)
    first, first_m = first_phase_margin(mu_O, T)
    return {'y': y, 'mu_O': mu_O, 'couple_spread': spread,
            'phase': phase, 'margin': margin,
            'first_phase': first, 'first_margin': first_m}
