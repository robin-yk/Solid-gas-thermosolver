"""Species, phases, and the element matrix for the C-H-O-Ti system.

The one place the system's vocabulary is written down: which gases exist,
which titanium and cerium condensed phases the registry carries, how many
Ti(3+) each Magneli member holds, and the element-count matrix that every
balance closes against. No solver lives here - `activeset.py` reads these
tables and `potential.py` reads them too, and neither may disagree with
the other about what a phase is made of.
"""

import numpy as np

from .shomate import R_ATM

GASES = ['CO2', 'H2', 'CO', 'H2O', 'CH4']

# Titanium phases, name -> (Ti atoms, Ti(3+) atoms) per formula unit.
#
# Every Magneli member Ti_n O_(2n-1) carries exactly two Ti(3+), for any n:
#   oxygen charge  = -2(2n-1) = -4n+2, so the Ti must supply +4n-2
#   3x + 4(n - x)  = 4n - x = 4n - 2   ->   x = 2
# so Ti2O3, Ti3O5, Ti4O7 and Ti5O9 all sit at two Ti(3+) per formula.
TI_PHASES = {
    'TiO2':   (1, 0),  # rutile, all Ti(4+)
    'Ti20O39': (20, 2),
    'Ti10O19': (10, 2),
    'Ti9O17': (9, 2),
    'Ti8O15': (8, 2),
    'Ti7O13': (7, 2),
    'Ti6O11': (6, 2),
    'Ti5O9': (5, 2),
    'Ti4O7': (4, 2),
    'Ti3O5': (3, 2),
    'Ti2O3': (2, 2),
}

# Ceria, the second oxide system the same engine runs. Ce2O3 carries two Ce(3+)
# per formula, the same way every Magneli member carries two Ti(3+).
CE_PHASES = {'CeO2': (1, 0), 'Ce2O3': (2, 2)}

# Every condensed phase, keyed name -> (metal atoms, reduced-metal atoms) per
# formula unit. The readout is the same count for either oxide.
METAL_PHASES = {**TI_PHASES, **CE_PHASES}

# Which titanium phases take part. Extend this list (and add the matching
# Shomate entry) to bring another member of the series into the calculation.
ACTIVE_TI_PHASES = ['TiO2', 'Ti10O19', 'Ti8O15', 'Ti6O11', 'Ti5O9',
                    'Ti4O7', 'Ti3O5', 'Ti2O3']

FORMULAS = {
    'CO2': {'C': 1, 'O': 2},
    'H2': {'H': 2},
    'CO': {'C': 1, 'O': 1},
    'H2O': {'H': 2, 'O': 1},
    'CH4': {'C': 1, 'H': 4},
    # Carried for the inert case only. Under N2 nothing but O2 can accept
    # the oxygen a line compound gives up, so it is the whole gas phase
    # story there; in a C-H-O feed it is negligible and left out, because
    # its equilibrium amount is a difference of two much larger extents.
    'O2': {'O': 2},
    'N2': {'N': 2},
    'CeO2': {'O': 2, 'Ce': 1},
    'Ce2O3': {'O': 3, 'Ce': 2},
    'TiO2': {'O': 2, 'Ti': 1},
    'Ti20O39': {'O': 39, 'Ti': 20},
    'Ti10O19': {'O': 19, 'Ti': 10},
    'Ti9O17': {'O': 17, 'Ti': 9},
    'Ti8O15': {'O': 15, 'Ti': 8},
    'Ti7O13': {'O': 13, 'Ti': 7},
    'Ti6O11': {'O': 11, 'Ti': 6},
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


def moles_of_gas(P_tot, V_cm3, T):
    """Ideal-gas moles filling V_cm3 at P_tot atm and T kelvin."""
    return P_tot * V_cm3 / (R_ATM * T)


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
