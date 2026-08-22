"""Dump the thermodynamic tables as JSON for the browser solver.

The page recomputes everything client-side rather than displaying a
precomputed sweep, so it needs the same coefficients the Python package uses.
Exporting them instead of retyping them is the point: a transcription slip
here would be invisible and would move every number on the page.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
sys.path.insert(0, str(ROOT))

from solidgas import shomate as S
from solidgas import waldner as W
from solidgas.equilibrium import ACTIVE_TI_PHASES, FORMULAS, GASES, TI_PHASES

def main():
    gases = GASES + ['O2', 'N2']
    payload = {
        'R_KJ': S.R_KJ,
        'R_ATM': S.R_ATM,
        'gases': gases,
        'shomate': {k: {'lo': v[0], 'hi': v[1], 'dHf': v[2],
                        'brk': S.BREAKPOINT.get(k, S.DEFAULT_BREAKPOINT),
                        'range': list(S.VALID_RANGE.get(k, (298.0, 6000.0)))}
                    for k, v in S.SHOMATE.items()},
        'waldner': {k: {'dHf': W.dHf(k), 'S298': v['S298'],
                        'ranges': [{'lo': r[0], 'hi': r[1], 'form': r[2],
                                    'c': list(r[3]), 'dHt': r[4]}
                                   for r in v['ranges']]}
                    for k, v in W.WALDNER.items()},
        'formulas': FORMULAS,
        'ti_phases': {k: list(v) for k, v in TI_PHASES.items()},
        'active': ACTIVE_TI_PHASES,
        'waldner_order': sorted(W.WALDNER,
                                key=lambda p: -W._stoich(p)[1] / W._stoich(p)[0]),
        'ti7o13_printed': W.TI7O13_PRINTED,
    }
    out = DATA / 'thermo_data.json'
    with out.open('w') as fh:
        json.dump(payload, fh, separators=(',', ':'))
    print(f'data/thermo_data.json written ({out.stat().st_size / 1024:.1f} kB), '
          f'{len(payload["shomate"])} NIST entries, {len(payload["waldner"])} Waldner entries')


if __name__ == '__main__':
    main()
