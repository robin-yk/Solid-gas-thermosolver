"""Dump the active-set solver's data layer as JSON, for the browser engine.

One source of truth: the tables in solidgas/shomate.py and solidgas/waldner.py.
The JavaScript engine evaluates the same coefficients with the same formulas;
tests/test_activeset_port.py holds the two engines together, and a freshness
test compares this file (and the copy inlined in the page) against the live
package so the exports cannot drift.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
sys.path.insert(0, str(ROOT))

from solidgas.activeset import (GASES, SOLIDS, GAS_FORMULA, STOICH, TI3,
                                ATOMIC_WEIGHT, MODE_NOTE)
from solidgas.shomate import SHOMATE, BREAKPOINT, DEFAULT_BREAKPOINT, R_KJ, R_ATM
from solidgas.waldner import WALDNER


def build():
    return {
        'R_KJ': R_KJ, 'R_ATM': R_ATM,
        'gases': GASES, 'solids': SOLIDS,
        'gas_formula': GAS_FORMULA,
        'stoich': {s: list(STOICH[s]) for s in SOLIDS},
        'ti3': TI3,
        'atomic_weight': ATOMIC_WEIGHT,
        'shomate': {g: {'lo': SHOMATE[g][0], 'hi': SHOMATE[g][1],
                        'dHf': SHOMATE[g][2],
                        'brk': BREAKPOINT.get(g, DEFAULT_BREAKPOINT)}
                    for g in GASES},
        'waldner': {s: {'dHf298': WALDNER[s]['dHf298'],
                        'S298': WALDNER[s]['S298'],
                        'ranges': [{'T0': r0, 'T1': r1, 'form': form,
                                    'c': list(c), 'dHt': dHt}
                                   for (r0, r1, form, c, dHt)
                                   in WALDNER[s]['ranges']]}
                    for s in SOLIDS},
        'mode_note': MODE_NOTE,
        'provenance': {
            'gases': 'NIST-JANAF 4th ed. Shomate (solidgas/shomate.py)',
            'solids': ('Waldner & Eriksson, Calphad 23 (1999) 189-218, '
                       'Appendix pp. 216-218 (solidgas/waldner.py)'),
        },
    }


def main():
    data = build()
    with (DATA / 'activeset_data.json').open('w') as fh:
        json.dump(data, fh, indent=1, sort_keys=True)
    print(f"activeset_data.json written "
          f"({len(json.dumps(data)) / 1024:.0f} kB)")


if __name__ == '__main__':
    main()
