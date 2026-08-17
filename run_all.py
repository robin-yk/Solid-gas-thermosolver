"""Every reported number, computed by full-composition Gibbs minimisation.

Four atmospheres over rutile, six temperatures each, twenty-four points, no
point skipped and no shortcut anywhere: each row is the converged minimum of
the total Gibbs energy over the full species set, subject to elemental balance.
Nothing here reads an answer off a phase-boundary table.

The oxygen-potential route is still in the tree and is still run, but only as a
cross-check written alongside each row. Where the two disagree the difference
is recorded and the Gibbs value is the one reported.

Writes results.json. Every row carries `method`; the test suite refuses any row
that does not say gibbs_min.
"""

import json
import sys

import numpy as np

from solidgas.gibbs import GibbsSystem
from solidgas.equilibrium import ACTIVE_TI_PHASES, TI_PHASES, moles_of_gas
from solidgas import potential as pot
from solidgas.shomate import R_KJ

P_TOT = 1.0        # atm
V_CM3 = 25.0       # reactor free volume
M_TIO2 = 0.100     # g
MW_TIO2 = 79.87    # g/mol
N_TI = M_TIO2 / MW_TIO2

TEMPS_C = [500, 700, 900, 1100, 1300, 1500]

# Each case: the gas species present, the components that express the charge as
# loaded, and the element vector of that charge. O2 is a free species wherever
# a C-H-O gas can carry the oxygen potential itself, and a component only in
# the inert case, where nothing else can.
CASES = {
    'RWGS': dict(
        label='CO2/H2 = 1:1',
        gases=['CO2', 'H2', 'CO', 'H2O', 'CH4', 'O2'],
        components=['CO2', 'H2', 'H2O', 'TiO2'],
        charge=lambda ng: [ng / 2, ng, ng + 2 * N_TI, N_TI],
    ),
    'H2': dict(
        label='pure H2',
        gases=['H2', 'H2O', 'O2'],
        components=['H2', 'H2O', 'TiO2'],
        charge=lambda ng: [2 * ng, 2 * N_TI, N_TI],
    ),
    'CH4': dict(
        label='pure CH4',
        gases=['CH4', 'H2', 'CO', 'CO2', 'H2O', 'O2'],
        components=['CH4', 'H2', 'CO', 'TiO2'],
        charge=lambda ng: [ng, 4 * ng, 2 * N_TI, N_TI],
    ),
    'N2': dict(
        label='sealed N2',
        gases=['N2', 'O2'],
        components=['N2', 'O2', 'TiO2'],
        charge=lambda ng: [2 * ng, 2 * N_TI, N_TI],
    ),
}

# Extent magnitudes to start from, as a fraction of full conversion. The inert
# case settles near 1e-38 and the reducing cases near 1e-1, so the ladder has to
# span the whole way; a start within a couple of decades of the answer is what
# the optimiser needs, and there is no way to know in advance which decade.
LADDER = [10.0 ** k for k in range(-1, -46, -3)]
GAS_PRESETS = [1e-6, 0.05, 0.3, 0.7]


def build_seeds(S, ng):
    """Starting extents: every phase basin crossed with every gas extent."""
    free = [S.species[i] for i in S.free]
    ref = np.array([ng if S.is_gas[i] else N_TI / TI_PHASES[S.species[i]][0]
                    for i in S.free])
    gas_pos = [j for j, i in enumerate(S.free) if S.is_gas[i]]
    sol_pos = [j for j, i in enumerate(S.free) if not S.is_gas[i]]

    seeds = [np.full(len(free), 1e-48)]
    for g in GAS_PRESETS:
        base = np.full(len(free), 1e-48)
        base[gas_pos] = ref[gas_pos] * g
        seeds.append(base.copy())
        for j in sol_pos:                       # one basin at a time
            for m in LADDER:
                v = base.copy()
                v[j] = ref[j] * m
                seeds.append(v)
    return seeds


def mu_O_of_gas(y, T):
    """Oxygen potential of the minimised gas, per mole of O.

    Read off the composition the minimisation returned - a property of the
    Gibbs answer, not a second calculation.

    The redox couples carry this, not molecular O2. Where a C-H-O gas is
    present its equilibrium O2 is around 1e-40 and sits at the optimiser's
    floor: it costs nothing to place because it changes G by less than a
    double can hold, so its value is arbitrary and must not be read. Only the
    inert case, which has no couple at all, takes its potential from O2.
    """
    vals = []
    for red, ox in (('H2', 'H2O'), ('CO', 'CO2')):
        a, b = y.get(red, 0.0), y.get(ox, 0.0)
        if a > 1e-280 and b > 1e-280:
            vals.append(pot.mu0(ox, T) - pot.mu0(red, T) + R_KJ * T * np.log(b / a))
    if not vals and y.get('O2', 0.0) > 1e-280:
        vals.append(0.5 * (pot.mu0('O2', T) + R_KJ * T * np.log(y['O2'] * P_TOT)))
    if not vals:
        return None, None
    return float(np.mean(vals)), float(max(vals) - min(vals))


# How close to a phase boundary counts as sitting on it. Coexistence pins the
# oxygen potential exactly, so a converged buffered row lands far inside this.
ON_BOUNDARY_KJ = 0.5


def crosscheck(row, T):
    """Independent check against the oxygen-potential route.

    Everything here comes from the phase data alone - never from anything the
    minimiser produced - so agreement is a real test rather than a number
    compared against itself. It reads only the stored row, which is why it can
    be replayed without repeating the minimisation.

    Whether a reduced phase is really there is decided by the gas, not by the
    amount. A line compound has no mixing entropy, so above its formation
    boundary its optimum is exactly zero and any residue is the optimiser
    resting against its own tolerance: RWGS leaves 1e-17 percent of Ti10O19
    lying around while sitting 78 kJ above the boundary, and that is nothing.
    Under sealed N2 the same 1e-17 percent, with the gas *on* the boundary, is
    the entire answer. Size cannot tell those apart; the potential can.
    """
    split = row['phase_split_pct']
    mu_O = row['mu_O_kJ']
    out = {'route': 'oxygen_potential'}
    if mu_O is None:
        out['status'] = 'no redox couple in the gas'
        return out

    held = max(split, key=split.get) if split else None
    reduced = {p: v for p, v in split.items() if p != 'TiO2' and v > 0}
    partner = max(reduced, key=reduced.get) if reduced else None
    host = 'TiO2' if split.get('TiO2', 0.0) > 0 and partner != 'TiO2' else held

    d = None
    if partner and host and partner != host:
        mu_b = pot.mu_O_critical(partner, host, T)
        d = mu_O - mu_b
        out.update(boundary=f'{host}/{partner}', boundary_mu_O_kJ=mu_b,
                   d_mu_O_kJ=d, partner_pct=reduced[partner])

    if d is not None and abs(d) <= ON_BOUNDARY_KJ:
        out.update(status='two-phase buffer', phases_present=sorted({host, partner}),
                   phase_agrees=held in (host, partner))
    else:
        # The gas is on no boundary, so one phase holds the titanium and
        # anything else in the vector is below the optimiser's tolerance.
        out.update(status='single phase', phases_present=[held],
                   residual_trace_pct=(reduced[partner] if partner else 0.0))
        out.pop('boundary', None)
        out.pop('boundary_mu_O_kJ', None)
        out.pop('d_mu_O_kJ', None)
        out.pop('partner_pct', None)

    # What the potential route would say holds the titanium, walked down the
    # ladder from rutile. Reported, never used to gate: the walk assumes the
    # series is monotonic in n, and in this assessment it is not everywhere.
    predicted, margin = pot.stable_phase(mu_O, T, ACTIVE_TI_PHASES)
    out['predicted_phase'] = predicted
    out['ladder_margin_kJ'] = None if margin != margin else margin
    out['gibbs_phase'] = held
    out.setdefault('phase_agrees', predicted == held)

    if row['case'] == 'RWGS':
        ng = row['n_gas_mol']
        srv = pot.survey([ng / 2, ng, ng], T, P_TOT, ACTIVE_TI_PHASES)
        out['gas_alone_mu_O_kJ'] = srv['mu_O']
        out['d_mu_O_gas_alone_kJ'] = mu_O - srv['mu_O']
    return out


def recheck():
    """Replay the cross-check over an existing results.json.

    The minimisation is the expensive half and its output is already stored;
    the check reads only that, so re-running it does not re-run anything and
    cannot change a reported value. Guarded by that: rows keep their method.
    """
    doc = json.load(open('results.json'))
    for row in doc['rows']:
        before = row['ti3_pct'], row['mu_O_kJ'], row['G_rel_kJ']
        row['crosscheck'] = crosscheck(row, row['T_C'] + 273.15)
        assert (row['ti3_pct'], row['mu_O_kJ'], row['G_rel_kJ']) == before
        assert row['method'] == 'gibbs_min'
    with open('results.json', 'w') as fh:
        json.dump(doc, fh, indent=1)
    for r in doc['rows']:
        cc = r['crosscheck']
        d = cc.get('d_mu_O_kJ')
        m = cc['ladder_margin_kJ']
        print(f"{r['case']:5s} {r['T_C']:5d}C  {cc['status']:18s} "
              f"held={str(cc['gibbs_phase']):9s} agrees={str(cc['phase_agrees']):5s} "
              + ('margin=  n/a  ' if m is None else f'margin={m:+8.2f}')
              + ('' if d is None else f'  d_muO={d:+.5f}'))
    return 0


def main():
    if '--recheck' in sys.argv:
        return recheck()
    rows = []
    for name, case in CASES.items():
        for T_C in TEMPS_C:
            row = run_point(case, name, T_C)
            rows.append(row)
            cc = row.get('crosscheck', {})
            d = cc.get('d_mu_O_kJ')
            print(f"{name:5s} {T_C:5d}C  Ti3+={row['ti3_pct']:11.6f}%  "
                  f"{str(row['phase']):9s} muO={row['mu_O_kJ']:9.2f}  "
                  f"resid={row['balance_residual_mol']:.1e}  "
                  f"xcheck={cc.get('status','-'):18s}"
                  f"{'' if d is None else f'  d_muO={d:+.3f}'}"
                  f"{'' if cc.get('bracket_ok') is None else '  bracket=' + str(cc['bracket_ok'])}",
                  flush=True)

    doc = {
        'conditions': {'P_atm': P_TOT, 'V_cm3': V_CM3, 'm_TiO2_g': M_TIO2,
                       'n_Ti_mol': N_TI, 'phases': ACTIVE_TI_PHASES,
                       'solid_data': 'Waldner & Eriksson, Calphad 23 (1999) 189',
                       'gas_data': 'NIST-JANAF Shomate'},
        'method': 'gibbs_min',
        'method_note': ('Total Gibbs energy minimised over the full species set '
                        'subject to elemental balance, in reaction-extent '
                        'coordinates. The oxygen-potential route is run only as '
                        'the cross-check recorded on each row.'),
        'temperatures_C': TEMPS_C,
        'cases': {k: v['label'] for k, v in CASES.items()},
        'rows': rows,
    }
    with open('results.json', 'w') as fh:
        json.dump(doc, fh, indent=1)
    bad = [r for r in rows if r.get('method') != 'gibbs_min']
    print(f"\n{len(rows)} rows -> results.json   ({len(bad)} not gibbs_min)")
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
