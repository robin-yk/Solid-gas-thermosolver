"""Two calculations for the reduced-rutile manuscript.

`activeset` answers what a gas charge does at equilibrium: the stable
Ti-O phase assemblage, the gas composition, the CO2 conversion and the
margin to the nearest reduced phase, by Gibbs-energy minimisation with
the condensed phases carried as an active set so an absent phase is
exactly zero rather than a numerical floor.

`vacancy_inventory` answers what a titration can support on its own:
how the measured oxygen removal compares with the oxygen the rutile
(110) surface could hold, from lattice geometry alone. It returns a
bound, not a partition. Where the below-surface oxygen sits - subsurface
point defects, bulk defects, extended defects - is not assigned here,
because no measurement in this work distinguishes them.

`species`, `shomate` and `waldner` are the data layer both read.
"""

from .shomate import (SHOMATE, VALID_RANGE, R_KJ, R_ATM, mu0, enthalpy, entropy,
                      heat_capacity, Kp)
from . import waldner
from . import activeset
from . import vacancy_inventory
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
    'waldner', 'activeset', 'vacancy_inventory',
]
