"""Two calculations for the reduced-rutile manuscript.

`activeset` answers what a gas charge does at equilibrium: the stable
Ti-O phase assemblage, the gas composition, the CO2 conversion and the
margin to the nearest reduced phase, by Gibbs-energy minimisation with
the condensed phases carried as an active set so an absent phase is
exactly zero rather than a numerical floor.

`vacancy_population` answers why the CO rate peaks after 600 C while
the measured oxygen removal keeps rising: only isolated surface
vacancies activate CO2, and they associate away by a thermally
activated second-order loss. It is fitted to the reduction series and
returns the effective partition of each sample's inventory.

`species`, `shomate` and `waldner` are the data layer both read.
"""

from .shomate import (SHOMATE, VALID_RANGE, R_KJ, R_ATM, mu0, enthalpy, entropy,
                      heat_capacity, Kp)
from . import waldner
from . import activeset
from . import vacancy_population
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
    'waldner', 'activeset', 'vacancy_population',
]
