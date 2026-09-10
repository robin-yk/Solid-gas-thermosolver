"""Independent 80-digit reference for the fresh-feed RWGS gas equilibrium.

This oracle is intentionally smaller than the Ti-O active-set reference.  It
evaluates the NIST Shomate polynomials with mpmath and solves the one-reaction
extent at high precision.  The production Python and browser implementations
are both compared with the committed result.
"""

import json
import sys
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from solidgas.shomate import (BREAKPOINT, DEFAULT_BREAKPOINT, R_KJ,
                              SHOMATE)


mp.mp.dps = 80
OUT = ROOT / 'data' / 'reference_rwgs_gas_only.json'
TEMPERATURES_C = [400, 600, 800, 900, 1000, 1200, 1500]
RATIOS = ['0.1', '0.5', '1', '3', '10', '100', '10000', '10000000']


def m(value):
    return mp.mpf(str(value))


def mu0(name, temperature):
    lo, hi, dhf = SHOMATE[name]
    row = lo if temperature <= m(BREAKPOINT.get(name, DEFAULT_BREAKPOINT)) else hi
    A, B, C, D, E, F, G, H = map(m, row)
    t = temperature / 1000
    enthalpy = m(dhf) + (A*t + B*t*t/2 + C*t**3/3 + D*t**4/4
                         - E/t + F - H)
    entropy = A*mp.log(t) + B*t + C*t*t/2 + D*t**3/3 - E/(2*t*t) + G
    return enthalpy - temperature*entropy/1000


def solve(temperature_c, ratio_text):
    temperature = m(temperature_c) + m('273.15')
    ratio = m(ratio_text)
    co2 = 1 / (1 + ratio)
    h2 = ratio / (1 + ratio)
    dg = (mu0('CO', temperature) + mu0('H2O', temperature)
          - mu0('CO2', temperature) - mu0('H2', temperature))
    kp = mp.exp(-dg / (m(R_KJ) * temperature))
    low, high = mp.mpf(0), min(co2, h2)
    for _ in range(500):
        extent = (low + high) / 2
        residual = extent**2 - kp*(co2 - extent)*(h2 - extent)
        if residual < 0:
            low = extent
        else:
            high = extent
    extent = (low + high) / 2
    return {
        'T_C': temperature_c,
        'H2_CO2_ratio': ratio_text,
        'Kp': mp.nstr(kp, 70),
        'extent_per_mol_feed': mp.nstr(extent, 70),
        'CO2_conversion_pct': mp.nstr(100*extent/co2, 70),
        'H2_utilization_pct': mp.nstr(100*extent/h2, 70),
        'CO_yield_per_mol_feed_pct': mp.nstr(100*extent, 70),
    }


def build():
    return {
        'title': 'High-precision reference: fresh-feed RWGS gas equilibrium',
        'generator': 'scripts/oracle_rwgs.py',
        'precision': 'mpmath, 80 working digits, 70 reported digits',
        'method': 'one-reaction ideal-gas extent; independent bisection',
        'scope': 'CO2 + H2 <=> CO + H2O; one mol initial gas; no solid',
        'rows': [solve(t, r) for t in TEMPERATURES_C for r in RATIOS],
    }


if __name__ == '__main__':
    OUT.write_text(json.dumps(build(), indent=2) + '\n')
    print(f'{OUT.relative_to(ROOT)} written')
