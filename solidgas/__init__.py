"""Solid-gas equilibrium thermodynamics for the C-H-O-Ti system.

Active-set Gibbs-energy minimisation over a CO2/H2/CO/H2O/N2/O2 gas in
contact with rutile TiO2 and the Magneli ladder down to Ti2O3, on
NIST-JANAF gas data and the Waldner-Eriksson condensed-phase assessment.
"""

from .shomate import (SHOMATE, VALID_RANGE, R_KJ, R_ATM, mu0, enthalpy, entropy,
                      heat_capacity, Kp)
from . import waldner
from .potential import (
    mu_O_from_gas, mu_O_critical, phase_ladder, stable_phase,
    first_phase_margin, ti_phases, solid_mu0, gas_equilibrium, survey,
    mu_O_to_pO2, pO2_to_mu_O, inert_atmosphere,
    SOURCES, DEFAULT_SOURCE,
)
from .species import (
    GASES, TI_PHASES, CE_PHASES, METAL_PHASES, ACTIVE_TI_PHASES, FORMULAS,
    SPECIES, ELEMENTS, ELEM, SOLIDS, gas_idx, sol_idx, moles_of_gas,
    ti3_percent, phase_split, mean_valence,
)

__all__ = [
    'SHOMATE', 'VALID_RANGE', 'R_KJ', 'R_ATM', 'mu0', 'enthalpy', 'entropy',
    'heat_capacity', 'Kp',
    'GASES', 'TI_PHASES', 'CE_PHASES', 'METAL_PHASES', 'ACTIVE_TI_PHASES',
    'FORMULAS', 'SPECIES', 'ELEMENTS', 'ELEM', 'SOLIDS', 'gas_idx', 'sol_idx',
    'moles_of_gas', 'ti3_percent', 'phase_split', 'mean_valence',
    'waldner', 'mu_O_from_gas', 'mu_O_critical', 'phase_ladder',
    'stable_phase', 'first_phase_margin', 'ti_phases', 'solid_mu0',
    'gas_equilibrium', 'survey', 'SOURCES', 'DEFAULT_SOURCE',
    'mu_O_to_pO2', 'pO2_to_mu_O', 'inert_atmosphere',
]
