"""Dump everything the browser solver needs, straight from the package.

One deliberate transformation happens here and nowhere else. The WGS page's
`shomate(c, T)` reads seven coefficients and computes H as A t + ... - E/t + F,
which lands on the absolute enthalpy because NIST's F already carries the
formation term. This package instead keeps eight coefficients and adds dHf
explicitly. Rather than touch that function - it is copied into the page
unchanged, and generalising it is how the gas path gets broken - the constant
is folded into the coefficient here:

    F' = F - H + dHf

so the seven-coefficient evaluator reproduces this package's mu0 exactly.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
sys.path.insert(0, str(ROOT))

from solidgas import waldner
from solidgas.shomate import SHOMATE, BREAKPOINT, DEFAULT_BREAKPOINT, R_KJ, R_ATM
from solidgas.equilibrium import FORMULAS, METAL_PHASES, ACTIVE_TI_PHASES

GASES = ['CO2', 'H2', 'CO', 'H2O', 'CH4', 'O2']
CE_SOLIDS = ['CeO2', 'Ce2O3']


def shomate_entry(name):
    lo, hi, dHf = SHOMATE[name]

    def fold(c):
        A, B, C, D, E, F, G, H = c
        return [A, B, C, D, E, F - H + dHf, G]

    return {'lo': fold(lo), 'hi': fold(hi), 'dHf': dHf,
            'brk': BREAKPOINT.get(name, DEFAULT_BREAKPOINT)}


def main():
    data = {
        'R_KJ': R_KJ, 'R_ATM': R_ATM,
        'shomate': {s: shomate_entry(s) for s in GASES + CE_SOLIDS},
        'calphad': {name: {'dHf298': waldner.dHf(name),
                           'S298': waldner.WALDNER[name]['S298'],
                           'ranges': [{'T0': r0, 'T1': r1, 'form': form,
                                       'c': list(c), 'dHt': dHt}
                                      for (r0, r1, form, c, dHt)
                                      in waldner.WALDNER[name]['ranges']]}
                    for name in ACTIVE_TI_PHASES},
        'formulas': {s: FORMULAS[s]
                     for s in GASES + CE_SOLIDS + ACTIVE_TI_PHASES},
        'metal_phases': {p: list(METAL_PHASES[p])
                         for p in CE_SOLIDS + ACTIVE_TI_PHASES},
        'gases': GASES,
        'systems': {
            'Ti-O': {'metal': 'Ti', 'host': 'TiO2', 'mw': 79.87,
                     'solids': ACTIVE_TI_PHASES, 'thermo': 'calphad',
                     'source': 'Waldner & Eriksson, Calphad 23 (1999) 189-218'},
            'Ce-O': {'metal': 'Ce', 'host': 'CeO2', 'mw': 172.12,
                     'solids': CE_SOLIDS, 'thermo': 'shomate',
                     'source': 'NIST-JANAF; one coefficient set across the range'},
        },
    }
    out = DATA / 'oxide_data.json'
    with out.open('w') as fh:
        json.dump(data, fh, separators=(',', ':'), sort_keys=True)
    print(f"data/oxide_data.json written ({out.stat().st_size / 1024:.1f} kB), "
          f"{len(data['shomate'])} Shomate + {len(data['calphad'])} CALPHAD entries")


if __name__ == '__main__':
    main()
