"""Invariants of the 80-digit reference file itself.

reference_results_high_precision.json is the anchor everything else is tested
against, so its own claims get checked before anything is compared to it:
coverage, method stamps, KKT uniqueness, element closure at oracle precision,
the boundary flip pair, and the data-against-data legs (NIST vs Waldner
offsets and heat capacities) that stand in for an external phase-diagram
digitisation.
"""

import json
import math
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
REF = ROOT / 'data' / 'reference_results_high_precision.json'

CASES = ['rwgs_1_1', 'h2_rich', 'product_mix', 'pure_h2', 'pure_n2',
         'n2_10ppm_o2']
TEMPS = [500, 700, 900, 1100, 1300, 1500]


@pytest.fixture(scope='module')
def doc():
    assert REF.exists(), 'run python3 oracle_tio.py first'
    return json.loads(REF.read_text())


def _num(v):
    return float(v) if isinstance(v, str) else v


def test_every_case_and_temperature_is_present(doc):
    keys = {(r['case'], r['T_C']) for r in doc['rows']}
    for c in CASES:
        for t in TEMPS:
            assert (c, t) in keys, (c, t)
    assert ('h2_h2o_below_boundary', 900) in keys
    assert ('h2_h2o_above_boundary', 900) in keys


def test_method_is_gibbs_min_on_every_row(doc):
    assert doc['method'] == 'gibbs_min'
    for r in doc['rows']:
        assert r['method'] == 'gibbs_min', (r['case'], r['T_C'])
        assert r['converged'] is True


def test_kkt_winner_is_unique_up_to_degeneracy(doc):
    """Convex program: feasible + KKT-clean must single out one active set."""
    for r in doc['rows']:
        passing = [c for c in r['candidates']
                   if c['feasible'] and c['kkt_ok']]
        assert len(passing) >= 1, (r['case'], r['T_C'])
        if len(passing) > 1:
            assert r['boundary_or_degenerate'], (r['case'], r['T_C'])


def test_element_residuals_close_at_oracle_precision(doc):
    for r in doc['rows']:
        for e, v in r['element_residuals_mol'].items():
            assert abs(float(v)) < 1e-55, (r['case'], r['T_C'], e, v)


def test_excluded_phases_are_exactly_zero(doc):
    for r in doc['rows']:
        for s, v in r['solid_moles'].items():
            if s not in r['active_condensed_phases']:
                assert _num(v) == 0.0, (r['case'], r['T_C'], s, v)
            else:
                assert _num(v) > 0.0, (r['case'], r['T_C'], s, v)


def test_inactive_reduced_costs_are_nonnegative(doc):
    for r in doc['rows']:
        for s, rc in r['inactive_phase_reduced_costs'].items():
            v = rc['per_formula_kJ']
            if v in (None, '-inf'):
                continue
            assert _num(v) > -1e-25, (r['case'], r['T_C'], s, v)


def test_dry_feeds_sit_on_the_buffer_and_wet_feeds_do_not(doc):
    """The physics the -inf rule encodes: bone-dry H2 and N2 always reduce
    infinitesimally; 10 ppm O2 in the same N2 stops it entirely."""
    rows = {(r['case'], r['T_C']): r for r in doc['rows']}
    for t in TEMPS:
        assert rows[('pure_h2', t)]['active_condensed_phases'] == \
            ['TiO2', 'Ti10O19'], t
        assert rows[('pure_n2', t)]['active_condensed_phases'] == \
            ['TiO2', 'Ti10O19'], t
        assert rows[('n2_10ppm_o2', t)]['active_condensed_phases'] == \
            ['TiO2'], t
        assert _num(rows[('n2_10ppm_o2', t)]['reduced_pct']) == 0.0


def test_oxidised_feeds_report_reduction_as_exact_zero(doc):
    """The number that used to print as 1e-17 percent."""
    rows = {(r['case'], r['T_C']): r for r in doc['rows']}
    for c in ('rwgs_1_1', 'h2_rich', 'product_mix', 'n2_10ppm_o2'):
        for t in TEMPS:
            r = rows[(c, t)]
            assert r['active_condensed_phases'] == ['TiO2'], (c, t)
            assert _num(r['reduced_pct']) == 0.0, (c, t)
            assert _num(r['solid_moles']['Ti10O19']) == 0.0, (c, t)


def test_pure_n2_trace_is_present_not_a_floor(doc):
    """Under sealed N2 the trace IS the answer; it must be positive, tiny,
    and grow monotonically with temperature."""
    rows = {(r['case'], r['T_C']): r for r in doc['rows']}
    last = 0.0
    for t in TEMPS:
        r = rows[('pure_n2', t)]
        v = _num(r['reduced_pct'])
        assert 0 < v < 1e-6, (t, v)
        assert v > last, t
        last = v
        assert 'Ti10O19' in r['log10_trace_amounts'], t


def test_boundary_pair_flips_across_the_h2o_h2_ratio(doc):
    rows = {r['case']: r for r in doc['rows'] if r['T_C'] == 900}
    below = rows['h2_h2o_below_boundary']
    above = rows['h2_h2o_above_boundary']
    assert below['active_condensed_phases'] == ['TiO2', 'Ti10O19']
    assert above['active_condensed_phases'] == ['TiO2']
    r = above['inactive_phase_reduced_costs']['Ti10O19']['per_mol_O_kJ']
    assert 0 < _num(r) < 0.05, r     # 1e-3 in ratio ~ RT*1e-3 ~ 0.01 kJ


def test_nist_waldner_offsets_match_the_pinned_values(doc):
    """Data-against-data: same offsets the package pins in test_waldner."""
    off = doc['boundary_validation']['dHf_offsets_waldner_minus_nist_kJ']
    for name, expect in (('TiO2', -6.03), ('Ti4O7', -12.00), ('Ti2O3', -5.29)):
        assert off[name] == pytest.approx(expect, abs=0.05), name


def test_heat_capacities_agree_between_assessments(doc):
    cp = doc['boundary_validation']['cp_relative_difference']
    for t, phases in cp.items():
        for name, rel in phases.items():
            assert rel < 0.01, (t, name, rel)


def test_boundary_oxygen_pressures_rise_with_temperature(doc):
    """Ellingham direction, stated in the invariant form: the equilibrium O2
    pressure over a reducible oxide rises monotonically with temperature -
    here by 28 decades across the range - and at fixed temperature a deeper
    phase always needs a lower pressure. (The elemental-reference mu_O_crit
    itself falls with T, because -S*T dominates both sides; that slope is not
    the invariant, the pressure is.)"""
    po = doc['boundary_validation']['equilibrium_pO2_atm']
    for phase in ('Ti10O19', 'Ti4O7', 'Ti2O3'):
        vals = [float(po[str(t)][phase]) for t in TEMPS]
        assert all(b > a for a, b in zip(vals, vals[1:])), (phase, vals)
    for t in TEMPS:
        row = {p: float(po[str(t)][p]) for p in ('Ti10O19', 'Ti4O7', 'Ti2O3')}
        assert row['Ti10O19'] > row['Ti4O7'] > row['Ti2O3'], (t, row)


def test_reduced_cost_bases_differ_only_by_the_oxygen_count(doc):
    for r in doc['rows']:
        for s, rc in r['inactive_phase_reduced_costs'].items():
            pf, po = rc['per_formula_kJ'], rc['per_mol_O_kJ']
            dO = rc['oxygen_removed_per_extent']
            if pf in (None, '-inf') or po in (None, '-inf') or not dO:
                continue
            assert _num(pf) == pytest.approx(_num(po) * dO, rel=1e-12)
