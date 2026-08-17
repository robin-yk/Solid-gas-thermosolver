"""The production active-set solver against the 80-digit oracle.

The oracle (oracle_tio.py) shares nothing with solidgas/activeset.py but the
coefficient tables; it solves a different formulation at 80 digits. Agreement
here is therefore a real measurement of the production solver, not a
comparison of a number with itself - which is what the old oxygen-potential
"cross-check" turned out to be.
"""

import ast
import json
import math
import pathlib

import pytest

from solidgas import activeset as A

ROOT = pathlib.Path(__file__).resolve().parent.parent
REF = ROOT / 'reference_results_high_precision.json'

CASES = {
    'rwgs_1_1': {'CO2': 1, 'H2': 1},
    'h2_rich': {'CO2': 1, 'H2': 3},
    'product_mix': {'CO2': 3, 'H2': 3, 'CO': 2, 'H2O': 2},
    'pure_h2': {'H2': 1},
    'pure_n2': {'N2': 1},
    'n2_10ppm_o2': {'N2': 999990, 'O2': 10},
}
TEMPS = [500, 700, 900, 1100, 1300, 1500]

Y_TOL = 1e-9            # measured worst 2e-13 against the oracle
R_TOL_KJ = 1e-8         # measured worst 8e-12 kJ
LOG_TOL = 1e-6          # trace amounts, in log10; measured worst 2e-13


def _num(v):
    return float(v) if isinstance(v, str) else v


@pytest.fixture(scope='module')
def oracle():
    assert REF.exists(), 'run python3 oracle_tio.py first'
    doc = json.loads(REF.read_text())
    return {(r['case'], r['T_C']): r for r in doc['rows']}, doc


@pytest.fixture(scope='module')
def solved():
    out = {}
    for case, feed in CASES.items():
        for t in TEMPS:
            out[(case, t)] = A.solve(feed, t + 273.15)
    return out


def test_active_sets_match_the_oracle(oracle, solved):
    rows, _ = oracle
    for k, r in solved.items():
        assert r['active_condensed_phases'] == \
            rows[k]['active_condensed_phases'], k


def test_gas_fractions_match_the_oracle(oracle, solved):
    rows, _ = oracle
    for k, r in solved.items():
        for g, v in rows[k]['gas_fractions'].items():
            ov = _num(v)
            if ov < 1e-25:
                continue
            d = abs(r['gas_fractions'][g] - ov) / ov
            assert d < Y_TOL, (k, g, d)


def test_reduced_costs_match_the_oracle(oracle, solved):
    rows, _ = oracle
    for k, r in solved.items():
        for s, rc in rows[k]['inactive_phase_reduced_costs'].items():
            ov = rc['per_mol_O_kJ']
            mine = r['inactive_phase_reduced_costs'][s]['per_mol_O_kJ']
            if ov == '-inf':
                assert math.isinf(mine) and mine < 0, (k, s)
                continue
            if ov is None:
                assert mine is None, (k, s)
                continue
            assert abs(mine - _num(ov)) < R_TOL_KJ, (k, s, mine, ov)


def test_trace_amounts_match_the_oracle_in_log10(oracle, solved):
    """The sealed-N2 ladder from 1e-41 mol up: resolvable, and resolved."""
    rows, _ = oracle
    for k, r in solved.items():
        for s, ov in rows[k]['log10_trace_amounts'].items():
            assert s in r['log10_trace_amounts'], (k, s)
            assert abs(r['log10_trace_amounts'][s] - ov) < LOG_TOL, (k, s)


def test_zero_reduction_is_exactly_zero(solved):
    """What used to be 1e-17 percent of residue is now the number it always
    was: zero, backed by a positive reduced cost instead of a threshold."""
    for (case, t), r in solved.items():
        if case in ('rwgs_1_1', 'h2_rich', 'product_mix', 'n2_10ppm_o2'):
            assert r['reduced_pct'] == 0.0, (case, t)
            for s, v in r['solid_moles'].items():
                if s != 'TiO2':
                    assert v == 0.0, (case, t, s)
            rmin = min(_finite_r(r))
            assert rmin > 0, (case, t, rmin)


def _finite_r(r):
    out = []
    for s, rc in r['inactive_phase_reduced_costs'].items():
        v = rc['per_mol_O_kJ']
        if v is not None and not math.isinf(v):
            out.append(v)
    return out


def test_excluded_phases_are_exactly_zero_everywhere(solved):
    for k, r in solved.items():
        for s, v in r['solid_moles'].items():
            if s in r['active_condensed_phases']:
                assert v > 0, (k, s)
            else:
                assert v == 0.0, (k, s)


def test_inactive_reduced_costs_are_nonnegative_at_the_winner(solved):
    for k, r in solved.items():
        for s, rc in r['inactive_phase_reduced_costs'].items():
            v = rc['per_formula_kJ']
            if v is None or math.isinf(v):
                continue
            assert v >= -A.TOL_KKT_KJ, (k, s, v)


def test_balance_closes_per_element(solved):
    for k, r in solved.items():
        for e, v in r['balance_by_element'].items():
            assert abs(v) < 1e-15, (k, e, v)


def test_method_and_mode_are_stamped(solved):
    for k, r in solved.items():
        assert r['method'] == 'gibbs_min', k
        assert 'CH4' in r['mode'] and 'excluded' in r['mode'], k
        assert r['converged'] is True, k


def test_result_carries_every_contract_field(solved):
    required = [
        'method', 'active_condensed_phases', 'inactive_phase_reduced_costs',
        'reduced_cost_basis', 'common_reference_G_kJ', 'G_total_kJ',
        'balance_by_element', 'solver_tolerance', 'reporting_detection_limit',
        'boundary_or_degenerate', 'log10_trace_amounts', 'delta_n_gas',
        'feed_mole_fractions', 'x_in_TiO2x', 'oxygen_to_gas_mol_O',
    ]
    r = next(iter(solved.values()))
    for field in required:
        assert field in r, field
    for s, rc in r['inactive_phase_reduced_costs'].items():
        assert set(rc) == {'per_formula_kJ', 'per_mol_O_kJ',
                           'oxygen_removed_per_extent'}, s


def test_phase_order_permutation_invariance():
    shuffled = ['Ti4O7', 'TiO2', 'Ti2O3', 'Ti10O19', 'Ti6O11',
                'Ti8O15', 'Ti3O5', 'Ti5O9']
    a = A.solve({'N2': 1}, 1173.15)
    b = A.solve({'N2': 1}, 1173.15, solids=shuffled)
    assert a['active_condensed_phases'] == b['active_condensed_phases']
    assert a['solid_moles'] == b['solid_moles']
    assert a['G_total_kJ'] == b['G_total_kJ']
    assert a['inactive_phase_reduced_costs'] == \
        b['inactive_phase_reduced_costs']


def test_all_candidates_share_one_gibbs_reference():
    """The reference is the loaded charge, so it cannot depend on which
    phases are allowed - solve with different registries and it must not
    move by even a bit."""
    full = A.solve({'CO2': 1, 'H2': 1}, 1173.15)
    narrow = A.solve({'CO2': 1, 'H2': 1}, 1173.15,
                     solids=['TiO2', 'Ti4O7'])
    assert full['common_reference_G_kJ'] == narrow['common_reference_G_kJ']
    assert full['G_total_kJ'] - full['G_rel_kJ'] == \
        pytest.approx(full['common_reference_G_kJ'], rel=1e-14)


def test_stable_set_flips_exactly_at_the_boundary(oracle):
    """Above/below gate: composition boundary from the oracle's data block."""
    _, doc = oracle
    ratio = _num(doc['boundary_validation']['h2_h2o_boundary']
                 ['H2O_over_H2_critical_ratio'])
    below = A.solve({'H2': 1, 'H2O': ratio * (1 - 1e-3)}, 1173.15)
    above = A.solve({'H2': 1, 'H2O': ratio * (1 + 1e-3)}, 1173.15)
    assert below['active_condensed_phases'] == ['TiO2', 'Ti10O19']
    assert above['active_condensed_phases'] == ['TiO2']
    assert above['reduced_pct'] == 0.0
    assert below['reduced_pct'] > 0


def test_initial_solid_other_than_rutile_re_oxidises():
    """Loading Ti4O7 under an oxidising product gas must walk back up the
    ladder, not sit still: the host is one more phase, not a special case."""
    r = A.solve({'CO2': 3, 'H2': 3, 'CO': 2, 'H2O': 2}, 1173.15,
                initial_solid='Ti4O7')
    assert 'Ti4O7' not in r['active_condensed_phases'] or \
        len(r['active_condensed_phases']) == 2
    assert r['oxygen_to_gas_mol_O'] < 0        # gas oxygen moved into the solid
    assert r['x_in_TiO2x'] < 0.25              # more oxidised than Ti4O7


def test_gas_charge_temperature_is_an_input():
    hot = A.solve({'CO2': 1, 'H2': 1}, 1173.15)
    cold = A.solve({'CO2': 1, 'H2': 1}, 1173.15, T_charge_K=298.15)
    assert cold['n_gas_charge_mol'] > hot['n_gas_charge_mol'] * 3.9
    assert cold['method'] == 'gibbs_min'


def test_legacy_solvers_are_not_imported_by_production():
    """Isolation, enforced: the legacy minimisers stay for regression but the
    production module must not touch them."""
    src = (ROOT / 'solidgas' / 'activeset.py').read_text()
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            assert node.module not in ('gibbs', 'equilibrium'), node.module
            for alias in node.names:
                assert alias.name not in ('minimise', 'solve', 'GibbsSystem')
        if isinstance(node, ast.Import):
            for alias in node.names:
                assert 'gibbs' not in alias.name
                assert 'equilibrium' not in alias.name


def test_legacy_minimiser_agrees_where_doubles_can_see():
    """The old reaction-extent SLSQP as a control group, on a mid-range point
    with no underflow anywhere: same gas, same (absent) reduction."""
    import numpy as np
    from solidgas.gibbs import GibbsSystem
    from solidgas.equilibrium import ACTIVE_TI_PHASES, moles_of_gas

    T = 1173.15
    ng = moles_of_gas(1.0, 25.0, T)
    n_metal = 0.100 / (A.ATOMIC_WEIGHT['Ti'] + 2 * A.ATOMIC_WEIGHT['O'])
    S = GibbsSystem(['CO2', 'H2', 'CO', 'H2O', 'O2'], ACTIVE_TI_PHASES,
                    ['CO2', 'H2', 'H2O', 'TiO2'])
    b = S.b_from(CO2=ng / 2, H2=ng / 2, TiO2=n_metal)
    seeds = []
    for frac in (0.2, 0.5, 0.8):
        v = np.full(len(S.free), 1e-48)
        for j, i in enumerate(S.free):
            if S.is_gas[i]:
                v[j] = ng * frac * 0.5
        seeds.append(v)
    old = S.minimise(b, T, seeds, 1.0)
    assert old is not None
    y_old = S.gas_fractions(old['n'])

    new = A.solve({'CO2': 1, 'H2': 1}, T)
    for g in ('CO2', 'H2', 'CO', 'H2O'):
        assert y_old[g] == pytest.approx(new['gas_fractions'][g], rel=2e-5), g
