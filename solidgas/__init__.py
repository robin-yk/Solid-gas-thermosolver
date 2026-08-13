"""Solid-gas equilibrium thermodynamics for the C-H-O-Ti system.

Gibbs-energy minimisation over CO2/H2/CO/H2O/CH4 gases in contact with
rutile TiO2 and the Magneli phase Ti4O7, using NIST-JANAF Shomate data.
"""

from .shomate import (SHOMATE, VALID_RANGE, R_KJ, R_ATM, mu0, enthalpy, entropy,
                      heat_capacity, Kp)
from .equilibrium import (
    SPECIES, ELEMENTS, ELEM, SOLIDS, gas_idx, sol_idx,
    G_total, solve, residual, gas_fractions, ti3_percent, moles_of_gas,
)

__all__ = [
    'SHOMATE', 'VALID_RANGE', 'R_KJ', 'R_ATM', 'mu0', 'enthalpy', 'entropy',
    'heat_capacity', 'Kp',
    'SPECIES', 'ELEMENTS', 'ELEM', 'SOLIDS', 'gas_idx', 'sol_idx',
    'G_total', 'solve', 'residual', 'gas_fractions', 'ti3_percent', 'moles_of_gas',
]
