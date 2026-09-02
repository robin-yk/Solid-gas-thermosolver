"""The isolated-surface-vacancy population model.

Two things are tested separately, because they fail for different reasons:
reproducing Table S1 at the parameters the supplement reports, and finding
those parameters again from a cold start. The second one turns on the
residual definition, which the supplement does not record, so every fitted
number here is computed under all three.
"""

import math
import pathlib

import pytest

from solidgas import vacancy_population as VP

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = VP.REPORTED
NU = 10.0 ** R['log10_nu']

# The populations the supplement's own run gives at those parameters.
EXPECTED = {
    'A600':  (8.950, 0.0035, 4.49),
    'R500':  (13.950, 0.0008, 22.90),
    'R600':  (15.240, 0.0567, 78.70),
    'R800':  (7.810, 7.5200, 177.77),
    'R1000': (0.177, 15.1600, 765.87),
}


@pytest.fixture(scope='module')
def series():
    return VP.load_series()


@pytest.fixture(scope='module')
def part(series):
    return {q['sample']: q for q in VP.partition(series, R['ns'], R['E_loss'], NU)}


# ------------------------------------------------------------------ the data

def test_the_five_samples_are_the_ones_the_table_lists(series):
    assert [r['sample'] for r in series] == list(EXPECTED)


def test_the_kelvin_temperatures_are_the_ones_specified(series):
    want = {'A600': 873.15, 'R500': 773.15, 'R600': 873.15,
            'R800': 1073.15, 'R1000': 1273.15}
    for r in series:
        assert r['T_red_K'] == pytest.approx(want[r['sample']])


# ----------------------------------------------- fixed parameters, Table S1

def test_the_closure_is_exact_and_not_merely_close(part):
    """n_below is closed form, not integrated, so this is an identity."""
    for s, q in part.items():
        assert abs(q['closure_umol_g']) <= VP.CLOSURE_TOL, s
        assert q['closure_umol_g'] == 0.0, s


def test_the_populations_reproduce_the_supplement(part):
    for s, (iso, assoc, below) in EXPECTED.items():
        q = part[s]
        assert q['n_iso_umol_g'] == pytest.approx(iso, abs=5e-3), s
        assert q['n_assoc_umol_g'] == pytest.approx(assoc, abs=5e-3), s
        assert q['n_below_umol_g'] == pytest.approx(below, abs=5e-3), s


def test_no_population_is_negative(part):
    for s, q in part.items():
        for k in ('n_iso_umol_g', 'n_assoc_umol_g', 'n_below_umol_g'):
            assert q[k] >= -1e-12, (s, k)


def test_the_normalising_sample_is_reproduced_exactly(part):
    assert part['R600']['r_CO_fitted'] == pytest.approx(
        part['R600']['r_CO_measured'], rel=1e-12)


def test_the_active_region_is_never_over_filled(part):
    for s, q in part.items():
        assert q['n_iso_umol_g'] + q['n_assoc_umol_g'] <= R['ns'] + 1e-9, s


# ------------------------------------------------------------ the integrator

@pytest.mark.parametrize('E_loss,nu_eff,T_K', [
    (0.5, 1e13, 773.15), (0.5, 1e10, 773.15), (1.0, 1e13, 1273.15),
    (3.76, 1e13, 1273.15), (8.0, 1e13, 773.15),
])
def test_it_stays_valid_where_an_explicit_step_would_not(E_loss, nu_eff, T_K):
    """The search visits k_loss of order 1e9 per second, where an explicit
    step of any practical size overflows or goes negative. The reduction
    plus an implicit solve has no such limit."""
    iso, assoc, below = VP.integrate(781.2, 3600, T_K, 15.3, E_loss, nu_eff)
    assert all(math.isfinite(v) for v in (iso, assoc, below))
    assert iso >= -1e-9 and assoc >= -1e-9 and below >= -1e-9
    assert abs(iso + assoc + below - 781.2) <= VP.CLOSURE_TOL
    assert iso <= 15.3 + 1e-9


def test_bdf_and_radau_land_on_the_same_answer(series):
    for r in series:
        a = VP.integrate(r['nV_total_umol_g'], r['t_red_s'], r['T_red_K'],
                         R['ns'], R['E_loss'], NU, method='Radau')
        b = VP.integrate(r['nV_total_umol_g'], r['t_red_s'], r['T_red_K'],
                         R['ns'], R['E_loss'], NU, method='BDF')
        assert a[0] == pytest.approx(b[0], rel=1e-6), r['sample']


def test_it_refuses_impossible_inputs():
    for bad in ((0.0, 1800, 873.15), (94.0, 0.0, 873.15)):
        with pytest.raises(ValueError):
            VP.integrate(bad[0], bad[1], bad[2], 15.3, 2.18, 5.5e6)
    with pytest.raises(ValueError):
        VP.integrate(94.0, 1800, 873.15, 0.0, 2.18, 5.5e6)


# ---------------------------------------------------------------- residuals

def test_the_three_residual_definitions_are_genuinely_different(series, part):
    """Raw least squares barely sees the slow samples; the rates span
    86-fold. The definitions agree about which sample misses worst and
    disagree about everything else, including where the optimum is - which
    is why the supplement has to say which it used."""
    rows = list(part.values())
    rates = [q['r_CO_measured'] for q in rows]
    assert max(rates) / min(rates) > 80

    def share(kind, sample):
        v = VP.residuals(rows, kind)
        tot = sum(x * x for x in v)
        i = [q['sample'] for q in rows].index(sample)
        return v[i] * v[i] / tot

    # the slowest sample carries three times the weight under log
    assert share('log', 'A600') > 2.0 * share('raw', 'A600')
    # and the optima are far apart
    a = VP.fit(series, kind='log', iters=250)
    b = VP.fit(series, kind='raw', iters=250)
    assert abs(a['log10_nu_eff'] - b['log10_nu_eff']) > 2.0


def test_an_unknown_residual_is_refused(part):
    with pytest.raises(ValueError):
        VP.residuals(list(part.values()), 'huber')


# ------------------------------------------------------------------ the fit

def test_the_log_residual_recovers_the_reported_parameters(series):
    """The whole disagreement was the objective. Under the log residual the
    cold-start fit lands on the supplement's own three numbers."""
    f = VP.fit(series, kind='log', iters=250)
    assert f['ns_umol_g'] == pytest.approx(R['ns'], rel=2e-3)
    assert f['E_loss_eV'] == pytest.approx(R['E_loss'], rel=2e-3)
    assert f['log10_nu_eff'] == pytest.approx(R['log10_nu'], rel=2e-3)
    assert not f['ns_at_bound']


def test_the_relative_residual_agrees_with_the_log_one(series):
    f = VP.fit(series, kind='relative', iters=250)
    assert f['ns_umol_g'] == pytest.approx(R['ns'], rel=1e-2)
    assert f['E_loss_eV'] == pytest.approx(R['E_loss'], rel=1e-2)


def test_raw_least_squares_lands_somewhere_else_entirely(series):
    """Not a bug in either - a consequence of an 86-fold spread in the data.
    Reported because a reader who picks the obvious objective will not get
    the supplement's numbers and needs to know why."""
    f = VP.fit(series, kind='raw', iters=250)
    assert f['E_loss_eV'] < 0.9 * R['E_loss']
    assert f['log10_nu_eff'] < R['log10_nu'] - 2.0
    assert f['ns_umol_g'] == pytest.approx(R['ns'], rel=0.02)


def test_the_capacity_is_the_bridging_row_and_not_the_subsurface_layer(series):
    for kind in VP.RESIDUALS:
        f = VP.fit(series, kind=kind, iters=250)
        assert 13.6 <= f['ns_umol_g'] <= 17.0, (kind, f['ns_umol_g'])


# -------------------------------------------------------- identifiability

def test_the_barrier_and_the_log_prefactor_are_one_parameter(series):
    c = VP.correlation(series, R['ns'], R['E_loss'], NU, kind='log')
    assert c['rho_E_lognu'] == pytest.approx(0.997, abs=0.002)


def test_the_capacity_is_nearly_independent_of_the_kinetic_pair(series):
    c = VP.correlation(series, R['ns'], R['E_loss'], NU, kind='log')
    assert abs(c['rho_ns_E']) < 0.1, c['rho_ns_E']
    assert abs(c['rho_ns_lognu']) < 0.1, c['rho_ns_lognu']


def test_pinning_the_prefactor_at_an_attempt_frequency_breaks_the_ratio(series):
    """The supplement's probe, recomputed under its own residual: nu = 1e13
    refits to 3.76 eV and then predicts R800/R600 = 0.94 against a measured
    0.51. Under raw least squares it refits to 3.52 eV instead and the
    ratio survives, which is the same warning the residual choice carries
    everywhere else in this file."""
    r = VP.fit_at_fixed_nu(series, 1e13, kind='log', iters=200)
    assert r['E_loss_eV'] == pytest.approx(3.76, abs=0.03)
    ratio = VP.rate_ratio(series, r['ns_umol_g'], r['E_loss_eV'], 1e13)
    assert ratio == pytest.approx(0.94, abs=0.02)
    assert VP.measured_ratio(series) == pytest.approx(0.51, abs=0.01)

    raw = VP.fit_at_fixed_nu(series, 1e13, kind='raw', iters=200)
    assert raw['E_loss_eV'] == pytest.approx(3.52, abs=0.05)
    assert VP.rate_ratio(series, raw['ns_umol_g'], raw['E_loss_eV'],
                         1e13) == pytest.approx(0.51, abs=0.02)


def test_the_fitted_point_reproduces_the_measured_ratio(series):
    assert VP.rate_ratio(series, R['ns'], R['E_loss'], NU) == pytest.approx(
        VP.measured_ratio(series), abs=0.01)


# -------------------------------------------------------------- the shape

def test_the_isolated_population_peaks_while_the_inventory_climbs(part):
    order = list(EXPECTED)
    iso = [part[s]['n_iso_umol_g'] for s in order]
    tot = [part[s]['nV_total_umol_g'] for s in order]
    assert tot == sorted(tot)
    assert max(iso) == iso[2], 'the isolated population peaks at R600'
    assert iso[-1] < 0.02 * max(iso), 'and collapses by R1000'


def test_the_narrative_numbers_hold(part):
    iso = {s: part[s]['n_iso_umol_g'] for s in EXPECTED}
    top = max(iso.values())
    assert iso['R500'] / top > 0.9, 'saturates by 500 C'
    assert iso['R800'] / top == pytest.approx(0.5, abs=0.03), 'halved at 800 C'
    assert iso['R1000'] / top == pytest.approx(0.012, abs=0.004), '1% at 1000 C'
    assert (part['R1000']['nV_total_umol_g']
            / part['R600']['nV_total_umol_g']) == pytest.approx(8.3, abs=0.2)


def test_most_of_the_inventory_is_below_the_region_in_every_reduced_sample(part):
    """The supplement says 'in every sample'. The annealed one is 33 percent
    below, not most: its total is 13.45 umol-O/g and the active region can
    nearly hold that alone. The four reduced samples are 62 to 98 percent
    below, so the sentence wants 'in every reduced sample'."""
    below = {s: part[s]['below_fraction'] for s in EXPECTED}
    assert below['A600'] == pytest.approx(0.33, abs=0.02)
    for s in ('R500', 'R600', 'R800', 'R1000'):
        assert below[s] > 0.5, (s, below[s])


# ------------------------------------------------------------ the step rule

def test_the_linear_partition_is_needed_and_a_step_rule_fails(series):
    """Filling the active region before anything goes below it over-predicts
    the annealed sample by half, under every residual definition."""
    lin = {q['sample']: q for q in VP.partition(series, R['ns'], R['E_loss'], NU)}
    stp = {q['sample']: q for q in VP.partition_step(series, R['ns'],
                                                     R['E_loss'], NU)}
    assert stp['A600']['r_CO_fitted'] > 1.4 * lin['A600']['r_CO_measured']
    for kind in VP.RESIDUALS:
        a = sum(v * v for v in VP.residuals(list(lin.values()), kind))
        b = sum(v * v for v in VP.residuals(list(stp.values()), kind))
        assert b > 20 * a, (kind, a, b)


# ------------------------------------------------------ the standing caveat

def test_the_module_says_what_it_must_not_be_used_for():
    src = ' '.join((ROOT / 'solidgas' / 'vacancy_population.py')
                   .read_text().split())
    for phrase in ('turnover-frequency denominator', 'does not predict',
                   'not an elementary barrier',
                   'say nothing that separates n_assoc from n_below'):
        assert phrase in src, phrase
