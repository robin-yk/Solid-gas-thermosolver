"""Solid-gas equilibrium thermodynamics for the C-H-O-Ti system.

Gibbs-energy minimisation over CO2/H2/CO/H2O/CH4 gases in contact with
rutile TiO2 and the Magneli phase Ti4O7, using NIST-JANAF Shomate data.
"""

from .shomate import (SHOMATE, VALID_RANGE, R_KJ, R_ATM, mu0, enthalpy, entropy,
                      heat_capacity, Kp)
from . import waldner
from .potential import (
    mu_O_from_gas, mu_O_critical, phase_ladder, stable_phase,
    first_phase_margin, ti_phases, solid_mu0, gas_equilibrium, survey,
    SOURCES, DEFAULT_SOURCE,
)
from .equilibrium import (
    GASES, TI_PHASES, ACTIVE_TI_PHASES, FORMULAS,
    SPECIES, ELEMENTS, ELEM, SOLIDS, gas_idx, sol_idx,
    G_total, solve, residual, phase_seeds, gas_fractions, ti3_percent, phase_split,
    mean_valence, moles_of_gas,
)

__all__ = [
    'SHOMATE', 'VALID_RANGE', 'R_KJ', 'R_ATM', 'mu0', 'enthalpy', 'entropy',
    'heat_capacity', 'Kp',
    'GASES', 'TI_PHASES', 'ACTIVE_TI_PHASES', 'FORMULAS',
    'SPECIES', 'ELEMENTS', 'ELEM', 'SOLIDS', 'gas_idx', 'sol_idx',
    'G_total', 'solve', 'residual', 'phase_seeds', 'gas_fractions', 'ti3_percent',
    'phase_split', 'mean_valence', 'moles_of_gas',
    'waldner', 'mu_O_from_gas', 'mu_O_critical', 'phase_ladder',
    'stable_phase', 'first_phase_margin', 'ti_phases', 'solid_mu0',
    'gas_equilibrium', 'survey', 'SOURCES', 'DEFAULT_SOURCE',
]
