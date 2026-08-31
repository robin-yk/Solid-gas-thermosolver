"""Gates for the particle engine.

The claims worth defending are not "the code runs".  They are: the cap
on reduction is a consequence and not a clamp, the dilute isotherm has
the exponent of a compensated vacancy, the depth profile has no fitted
shape parameter left in it, and the inventory that goes in comes back
out.  Each of those has a test here that fails if it stops being true.
"""

import json
import math
import os

import pytest

from solidgas import particle as P
from solidgas.particle import freeenergy as F
from solidgas.particle import material as M
from solidgas.particle import profile as PR
from solidgas.particle import spacecharge as SC
from solidgas.particle import transport as T

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, 'data', 'reference_particle.json')


@pytest.fixture(scope='module')
def p():
    return P.load()


@pytest.fixture(scope='module')
def part(p):
    return P.particle(p)


@pytest.fixture(scope='module')
def flagship(part):
    return P.equilibrate(95.0, 600.0, part)


@pytest.fixture(scope='module')
def ref():
    with open(REF) as fh:
        return json.load(fh)


# ------------------------------------------------------------ free energy

def test_the_cap_is_a_consequence_of_the_entropy_not_a_clamp():
    """Nothing in the code compares theta to 0.25 and stops it.

    The cation sublattice runs out of Ti(3+) sites there, so the
    configurational term diverges on its own.  Push the chemical
    potential absurdly high and the occupancy still cannot pass the
    quarter, and it approaches from below rather than sitting on it."""
    kt = M.kt_eV(600.0)
    last = 0.0
    for mu in (0.5, 1.0, 2.0, 3.0, 4.0):
        th = F.theta_of_mu(mu, 0.0, kt)
        assert th < F.THETA_CAP
        assert th > last
        last = th
    # past here the gap to the cap is smaller than a double can hold
    # next to a quarter, which is arithmetic running out, not a clamp
    assert F.THETA_CAP - last < 1e-11

    # the signature that distinguishes an entropic divergence from a
    # clamp: the gap to the cap closes exponentially, at exactly the rate
    # the second logarithm sets.  A clamp would close it abruptly.
    a, b = 3.0, 4.0
    ga = F.THETA_CAP - F.theta_of_mu(a, 0.0, kt)
    gb = F.THETA_CAP - F.theta_of_mu(b, 0.0, kt)
    slope = (math.log(gb) - math.log(ga)) / (b - a)
    assert slope == pytest.approx(-1.0 / (2.0 * kt), rel=1e-4)

    # and it really does diverge: halving the gap adds 2 ln 2 to phi,
    # without bound
    for gap in (1e-6, 1e-9, 1e-12):
        assert (F.phi(F.THETA_CAP - gap)
                - F.phi(F.THETA_CAP - 2 * gap)) == pytest.approx(
                    2.0 * math.log(2.0), rel=1e-5)
    assert F.phi(F.THETA_CAP - 1e-12) > 50.0
    assert math.isinf(F.phi(F.THETA_CAP))


def test_the_uncompensated_form_has_no_cap_at_all():
    """The counterexample the compensated form exists to rule out.

    Left on the oxygen sublattice alone, the same conditions drive the
    site fraction towards one - a completely deoxygenated lattice, which
    is not a phase of titanium oxide."""
    kt = M.kt_eV(600.0)
    th = F.theta_of_mu(2.0, 0.0, kt, compensated=False)
    assert th > 0.999
    assert th > 3.9 * F.THETA_CAP


def test_the_dilute_isotherm_is_the_minus_one_sixth_law():
    """theta ~ p(O2)^(-1/6), the doubly ionised compensated vacancy.

    Two logarithms add to 3 kT ln theta in the dilute limit, so
    d ln theta / d mu_V is one third of a kT rather than one.  With
    mu_V = -mu_O = -(1/2) kT ln p that makes the pressure exponent
    -1/6.  Left on the oxygen sublattice alone the slope is 1/kT and the
    exponent is -1/2, the neutral vacancy - and an enrichment between
    two sites that is the CUBE of the right one.

    The law is asymptotic, so it is tested where the material is
    genuinely dilute; the exact asymptote itself is checked separately
    against the full inversion."""
    kt = M.kt_eV(600.0)
    e = 1.31
    a, b = -2.0, -1.8            # theta of order 1e-7 at both ends

    def slope(comp):
        ta = F.theta_of_mu(a, e, kt, compensated=comp)
        tb = F.theta_of_mu(b, e, kt, compensated=comp)
        return (math.log(tb) - math.log(ta)) / (b - a)

    assert F.theta_of_mu(b, e, kt) < 1e-6
    assert slope(True) == pytest.approx(1.0 / (3.0 * kt), rel=1e-5)
    assert slope(False) == pytest.approx(1.0 / kt, rel=1e-5)

    # the same statement as a pressure exponent
    assert -0.5 * kt * slope(True) == pytest.approx(-1.0 / 6.0, abs=1e-5)
    assert -0.5 * kt * slope(False) == pytest.approx(-0.5, abs=1e-5)


def test_the_chemical_potential_is_the_derivative_of_the_free_energy():
    kt = M.kt_eV(600.0)
    for th in (1e-4, 1e-2, 0.05, 0.15, 0.24):
        h = th * 1e-6
        num = (F.free_energy_per_site(th + h, 5.7, kt)
               - F.free_energy_per_site(th - h, 5.7, kt)) / (2 * h)
        assert num == pytest.approx(F.mu_of_theta(th, 5.7, kt), rel=1e-6)


def test_the_inverse_is_exact_across_the_whole_domain_and_any_ratio():
    """Including within a rounding of both ends, and for ratios other
    than four - the surface bridging row is the reason that matters."""
    for ratio in (1.0, 2.0, 4.0, 6.0):
        top = F.cap(ratio)
        probes = [top * 10 ** (-k) for k in range(1, 26)]
        probes += [top * (1 - 10 ** (-k)) for k in range(1, 14)]
        probes += [top * f for f in (0.1, 0.3, 0.5, 0.7, 0.9)]
        for th in probes:
            back = F.theta_of_phi(F.phi(th, ratio), ratio=ratio)
            assert back == pytest.approx(th, rel=1e-12)


def test_the_dilute_asymptote_is_the_real_limit():
    kt = M.kt_eV(600.0)
    for mu in (-1.0, 0.0, 0.3):
        exact = F.theta_of_mu(mu, 5.7, kt)
        assert F.dilute_theta(mu, 5.7, kt) == pytest.approx(exact, rel=1e-12)


# -------------------------------------------------------------- material

def test_the_decay_length_reproduces_the_dft_triple_with_nothing_fitted():
    """E(z) = E_bulk - dE exp(-z/xi) has two parameters and the DFT set
    has three energies with the bulk as the asymptote, so both are
    determined.  The test is that the profile passes through the surface
    and first-subsurface points exactly, not approximately."""
    p = P.load()
    for key in p['vacancy_energetics']['presets']:
        q = p['vacancy_energetics']['presets'][key]
        s = M.segregation(p, key)
        assert M.formation_energy(0.0, s) == pytest.approx(
            q['surface_eV'], abs=1e-12)
        assert M.formation_energy(s['d110_nm'], s) == pytest.approx(
            q['subsurface_eV'], abs=1e-12)
        assert M.formation_energy(1e4, s) == pytest.approx(
            q['bulk_eV'], abs=1e-12)
        assert 1.0 < s['xi_in_layers'] < 4.0


def test_the_decay_length_is_not_a_borrowed_number():
    """The ceria study assumed 0.7 nm as 'two or three atomic layers'.
    Ours comes out of our own DFT and lands next to it without being
    told to, which is a check on both."""
    p = P.load()
    assert M.segregation(p, 'sx')['xi_nm'] == pytest.approx(0.5427, abs=5e-4)
    assert M.segregation(p, 'gga')['xi_nm'] == pytest.approx(0.7250, abs=5e-4)


def test_a_non_monotone_triple_is_refused_rather_than_fitted():
    p = P.load()
    p = json.loads(json.dumps(p))
    p['vacancy_energetics']['presets']['sx']['subsurface_eV'] = 9.0
    with pytest.raises(ValueError, match='monotonic'):
        M.segregation(p, 'sx')


def test_the_two_reduction_conventions_do_not_get_mixed(part):
    """theta counts oxygen SITES; TiO(2-x) has x = 2 theta.  The ceria
    report calls this out as the easiest factor of two to lose."""
    th = M.theta_from_loading(95.0, part)
    assert th == pytest.approx(95.0 / part['N_oxygen_umol_g'], rel=1e-14)
    assert M.loading_from_theta(th, part) == pytest.approx(95.0, rel=1e-14)
    assert 2.0 * th == pytest.approx(0.0075866, abs=1e-6)


def test_the_bridging_row_is_not_smeared_into_the_bulk_density(part):
    """A full d110 oxygen layer holds four oxygens per (1x1) cell and the
    bridging row holds one, so the row is a quarter of a layer.  Giving
    it a bulk density would quadruple the cheapest sites in the
    particle."""
    p = P.load()
    a, c = p['lattice_constants_A']['a'], p['lattice_constants_A']['c']
    full_layer_cm2 = 4e16 / (c * a * math.sqrt(2.0))
    assert M.bridging_density_cm2(p) * 4.0 == pytest.approx(
        full_layer_cm2, rel=1e-12)
    assert part['N_surface_umol_g'] == pytest.approx(13.557, abs=1e-2)

    # a d110 slab of the interior, for comparison: four times the row
    seg = part['segregation']
    r_in = part['radius_nm'] - 0.5 * seg['d110_nm']
    slab_frac = 1.0 - ((r_in - seg['d110_nm']) / r_in) ** 3
    assert part['N_interior_umol_g'] * slab_frac == pytest.approx(
        4.0 * part['N_surface_umol_g'], rel=0.02)


# --------------------------------------------------------------- profile

def test_the_inventory_that_goes_in_comes_back_out(part):
    for vo in (1.0, 20.0, 95.0, 400.0):
        r = P.equilibrate(vo, 600.0, part)
        assert (r['held_umol_g'] + r['umol_extended']) == pytest.approx(
            vo, rel=1e-11)


def test_the_layer_sum_is_what_runs_and_the_smear_costs_what_it_costs(
        part, flagship):
    """The interior is a layer stack, not a continuum.  The smeared form
    is kept only so the size of the approximation can be quoted; if it
    ever became the default this test would notice."""
    mu, kt = flagship['mu_eV'], flagship['kt_eV']
    exact = PR.interior_average(mu, kt, part)
    smeared = PR.interior_average_smeared(mu, kt, part)
    assert flagship['theta_interior_mean'] == pytest.approx(exact, rel=1e-14)
    assert smeared != pytest.approx(exact, rel=1e-6)
    assert abs(smeared - exact) / exact == pytest.approx(0.00517, abs=5e-4)


def test_the_deep_truncation_is_exact_not_a_cutoff(part, flagship):
    """interior_average stops summing once layers stop differing from the
    bulk.  Forcing it to keep going must not move the answer."""
    mu, kt = flagship['mu_eV'], flagship['kt_eV']
    loose = PR.interior_average(mu, kt, part, rel_tol=1e-15)
    tight = PR.interior_average(mu, kt, part, rel_tol=0.0)
    assert loose == pytest.approx(tight, rel=1e-14)


def test_the_near_surface_shell_holds_far_more_than_its_volume(
        part, flagship):
    """The finding the whole construction exists to produce: a shell too
    thin to see holds a disproportionate share of the oxygen."""
    rows = PR.shell_split(flagship, part)
    one_nm = rows[0]
    assert one_nm['depth_nm'] == 1.0
    assert one_nm['volume_fraction'] < 0.01
    assert one_nm['inventory_fraction'] > 0.10
    assert (one_nm['inventory_fraction'] / one_nm['volume_fraction']) > 10.0
    for a, b in zip(rows[:-1], rows[1:]):
        assert b['inventory_fraction'] >= a['inventory_fraction']


def test_the_surface_saturates_long_before_the_bulk_does(flagship):
    """The bridging row costs 1.31 eV less than the bulk, so at any
    loading that registers on a titration it is already reconstructed.
    It is a small reservoir, not the main one."""
    assert flagship['surface']['phase'] == '1x2'
    assert flagship['umol_surface'] < 0.1 * flagship['loading_umol_g']
    assert flagship['theta_first_layer'] > 10.0 * flagship[
        'theta_bulk_plateau']


def test_the_reconstruction_pins_the_chemical_potential(part):
    """Inside the two-phase window the added-row phase is a line
    compound, so mu cannot move: the extra oxygen converts surface
    instead of raising the occupancy anywhere."""
    ph = part['surface_phases']
    kt = M.kt_eV(600.0)
    mu_t = F.mu_of_theta(ph['theta_transition'],
                         part['segregation']['e_surface_eV'], kt)
    span = part['N_surface_umol_g'] * (ph['theta_reconstructed_eff']
                                       - ph['theta_transition'])
    base = PR.inventory(mu_t, kt, part)
    lo = P.equilibrate(base + 0.2 * span, 600.0, part)
    hi = P.equilibrate(base + 0.8 * span, 600.0, part)
    assert lo['mu_eV'] == pytest.approx(hi['mu_eV'], abs=1e-9)
    assert lo['limited_by'] == 'surface reconstruction'
    assert hi['surface']['reconstructed_fraction'] > lo['surface'][
        'reconstructed_fraction']


def test_shear_planes_take_the_inventory_the_lattice_cannot_dissolve(part):
    """Rutile does not dissolve point defects all the way to Ti2O3; past
    x = 0.008 it throws crystallographic shear planes.  This is the one
    piece of the construction the ceria system has no analogue for, and
    dropping it would let the model report a bulk site fraction the
    crystallography does not permit."""
    r = P.equilibrate(400.0, 600.0, part)
    assert r['limited_by'] == 'crystallographic shear'
    assert r['umol_extended'] > 0.0
    assert r['theta_bulk_plateau'] == pytest.approx(
        part['theta_solubility'], rel=1e-9)
    free = P.equilibrate(400.0, 600.0, part, shear_limit=False)
    assert free['theta_bulk_plateau'] > part['theta_solubility']
    assert free['umol_extended'] == 0.0

    # the ceiling is close enough to the shipped loading to matter: the
    # particle cannot hold much more than it already does as point
    # defects, whatever the gas does
    ceiling = r['point_defect_ceiling_umol_g']
    assert ceiling == pytest.approx(112.4, abs=0.5)
    assert P.equilibrate(95.0, 600.0, part)['limited_by'] is None
    assert 95.0 / ceiling > 0.8


def test_a_loading_past_the_cation_sublattice_is_refused(part):
    with pytest.raises(ValueError, match='Ti2O3'):
        P.equilibrate(part['N_oxygen_umol_g'] * 0.3, 600.0, part)


def test_the_polaron_ratio_is_an_assumption_with_a_measured_size(part):
    """At the bridging row the local cation reservoir is one Ti(3+) site
    per vacancy, not four; ratio four says the polarons drain into the
    subsurface, which is what DFT reports.  This bounds what that choice
    is worth rather than leaving it undeclared."""
    kt = M.kt_eV(600.0)
    four = P.equilibrate(95.0, 600.0, part)
    one = P.equilibrate(95.0, 600.0, part, ratio=1.0)

    # in the dilute limit the ratio enters phi only as an additive
    # constant, so it shifts the chemical potential and leaves the
    # distribution where it was: the choice is close to a gauge
    assert one['mu_eV'] - four['mu_eV'] == pytest.approx(
        -2.0 * kt * math.log(4.0), rel=0.05)
    assert one['theta_bulk_plateau'] == pytest.approx(
        four['theta_bulk_plateau'], rel=0.02)
    assert one['umol_interior'] == pytest.approx(
        four['umol_interior'], rel=0.01)

    # where it is not a gauge is near the cap, and the first subsurface
    # layer is near enough to the cap to feel it
    e4 = four['theta_first_layer'] / four['theta_bulk_plateau']
    e1 = one['theta_first_layer'] / one['theta_bulk_plateau']
    assert e1 > 1.1 * e4


def test_a_hotter_particle_needs_less_driving_force(part):
    mus = [P.equilibrate(95.0, t, part)['mu_eV'] for t in (400, 600, 800)]
    assert mus[0] > mus[1] > mus[2]


# -------------------------------------------------------------- transport

def test_the_crank_series_reproduces_the_reference_fourier_number():
    """The ceria study solves for the Fourier number giving twenty
    percent uptake and gets 0.003912.  Landing on the same number from
    our own series is a cross-check on both."""
    assert T.tau_for_uptake(0.2) == pytest.approx(0.003912, abs=2e-6)
    assert T.uptake(T.tau_for_uptake(0.2)) == pytest.approx(0.2, abs=1e-12)


def test_transport_is_not_the_limit_at_the_shipped_condition(part, p):
    v = T.verdict(600.0, part, exposure_s=300.0, p=p)
    assert v['equilibrium_reachable']
    assert v['t_homogenise_s'] < 1e-2
    assert v['ratio'] < 1e-4


def test_the_margin_on_that_claim_is_reported_not_hidden(part, p):
    """Saying transport is fast is only useful with the distance to it
    not being fast.  The barrier that would close the gap is about an eV
    above the literature one."""
    v = T.verdict(600.0, part, exposure_s=300.0, p=p)
    assert v['barrier_margin_eV'] > 0.5
    assert v['E_m_that_would_break_it_eV'] == pytest.approx(1.645, abs=0.01)
    slow = T.verdict(600.0, part, exposure_s=300.0, p=p,
                     e_m_eV=v['E_m_that_would_break_it_eV'] + 0.2)
    assert not slow['equilibrium_reachable']


# ----------------------------------------------------------- space charge

def test_local_electroneutrality_is_checked_rather_than_assumed(
        part, p, flagship):
    """The ceria study found screening much shorter than segregation and
    said so was a consequence of its concentration.  Here they are the
    same length, so the assumption is on the edge, and the assessment
    has to say that out loud."""
    a = SC.assess(flagship, part, p)
    assert 0.7 < a['ratio_debye_over_xi'] < 1.5
    assert 'cannot be separated' in a['verdict']
    for eps in (50.0, 170.0):
        r = SC.assess(flagship, part, p, eps_r=eps)['ratio_debye_over_xi']
        assert 0.5 < r < 2.0


def test_how_far_outside_the_ceria_regime_this_sits(part, p, flagship):
    c = SC.ceria_comparison(flagship, part, p)
    assert c['dilution_factor'] == pytest.approx(13.2, abs=0.3)
    assert c['screening_lengthened_by'] == pytest.approx(3.6, abs=0.1)


# --------------------------------------------------------- oracle parity

def test_the_oracle_agrees_on_the_free_energy(ref):
    for row in ref['phi']:
        got = F.phi(float(row['theta']), float(row['ratio']))
        assert got == pytest.approx(float(row['phi']), rel=1e-13)
    for row in ref['phi_inverse']:
        got = F.theta_of_phi(float(row['x']), ratio=float(row['ratio']))
        assert got == pytest.approx(float(row['theta']), rel=1e-12)


def test_the_oracle_agrees_on_the_power_law(ref):
    kt = float(ref['power_law']['kt_eV'])
    e = float(ref['power_law']['e_form_eV'])
    for row in ref['power_law']['rows']:
        mu = float(row['mu_eV'])
        assert F.theta_of_mu(mu, e, kt) == pytest.approx(
            float(row['theta']), rel=1e-12)
        assert F.theta_of_mu(mu, e, kt, compensated=False) == pytest.approx(
            float(row['theta_oxygen_only']), rel=1e-12)


def test_the_oracle_agrees_on_the_decay_length(ref, p):
    for key, want in ref['segregation'].items():
        s = M.segregation(p, key)
        assert s['xi_nm'] == pytest.approx(float(want['xi_nm']), rel=1e-13)
        assert s['dE_seg_eV'] == pytest.approx(
            float(want['dE_seg_eV']), rel=1e-13)
        for n, e in enumerate(want['E_at_layers_eV']):
            assert M.formation_energy(n * s['d110_nm'], s) == pytest.approx(
                float(e), rel=1e-13, abs=1e-15)


def test_the_oracle_agrees_on_the_site_budget(ref, part):
    q = ref['particle']
    assert part['radius_nm'] == pytest.approx(float(q['radius_nm']), rel=1e-13)
    assert part['N_surface_umol_g'] == pytest.approx(
        float(q['N_surface_umol_g']), rel=1e-13)
    assert part['N_interior_umol_g'] == pytest.approx(
        float(q['N_interior_umol_g']), rel=1e-13)


def test_the_oracle_agrees_on_every_equilibrium_row(ref, part):
    for row in ref['equilibrium']:
        r = P.equilibrate(row['VO_umol_g'], row['T_C'], part)
        assert r['mu_eV'] == pytest.approx(float(row['mu_eV']), abs=1e-11)
        for k, want in (('theta_interior_mean', 'theta_interior_mean'),
                        ('theta_bulk_plateau', 'theta_bulk_plateau'),
                        ('theta_first_layer', 'theta_first_layer'),
                        ('umol_surface', 'umol_surface'),
                        ('umol_interior', 'umol_interior')):
            assert r[k] == pytest.approx(float(row[want]), rel=1e-9), (
                '%s at %g C %g umol' % (k, row['T_C'], row['VO_umol_g']))
        assert r['limited_by'] == row['limited_by']


def test_the_oracle_agrees_on_the_crank_series(ref):
    for f, tau in ref['crank']['tau_for_uptake'].items():
        assert T.tau_for_uptake(float(f)) == pytest.approx(
            float(tau), rel=1e-9)
    for tau, u in ref['crank']['uptake_at_tau'].items():
        assert T.uptake(float(tau)) == pytest.approx(float(u), rel=1e-13)


# ------------------------------------------------------- the prose is bound

def _doc():
    with open(os.path.join(ROOT, 'docs', 'particle-model.md')) as fh:
        return fh.read()


def test_every_number_the_documentation_quotes_is_the_engine_s(part, p,
                                                               flagship):
    """The documentation states results, not intentions, so it can go
    stale silently.  Each number below is recomputed and matched against
    the string in the file, so editing one without the other fails."""
    doc = _doc()
    seg = part['segregation']
    a = SC.assess(flagship, part, p)
    v = T.verdict(600.0, part, exposure_s=300.0, p=p)
    shell = PR.shell_split(flagship, part)[0]
    smear = abs(PR.interior_average_smeared(flagship['mu_eV'],
                                            flagship['kt_eV'], part)
                - flagship['theta_interior_mean']
                ) / flagship['theta_interior_mean']

    claims = [
        ('1.31 eV', seg['dE_seg_eV'], '%.2f eV'),
        ('0.543 nm', seg['xi_nm'], '%.3f nm'),
        ('1.67', seg['xi_in_layers'], '%.2f'),
        ('0.2336 eV', flagship['mu_eV'], '%.4f eV'),
        ('0.0654', flagship['theta_first_layer'], '%.4f'),
        ('0.00334', flagship['theta_bulk_plateau'], '%.5f'),
        ('6.8 μmol-O/g', flagship['umol_surface'], '%.1f μmol-O/g'),
        ('88.2 μmol-O/g', flagship['umol_interior'], '%.1f μmol-O/g'),
        ('112.4 μmol-O/g', flagship['point_defect_ceiling_umol_g'],
         '%.1f μmol-O/g'),
        ('0.57 nm', a['debye_nm_at_bulk'], '%.2f nm'),
        ('0.54 nm', a['xi_nm'], '%.2f nm'),
        ('0.74', SC.assess(flagship, part, p, eps_r=50.0)[
            'ratio_debye_over_xi'], '%.2f'),
        ('1.37', SC.assess(flagship, part, p, eps_r=170.0)[
            'ratio_debye_over_xi'], '%.2f'),
        ('1.645 eV', v['E_m_that_would_break_it_eV'], '%.3f eV'),
        ('0.54 ms', v['t_homogenise_s'] * 1e3, '%.2f ms'),
        ('12.4 % of the oxygen', 100 * shell['inventory_fraction'],
         '%.1f %% of the oxygen'),
        ('0.7 % of the volume', 100 * shell['volume_fraction'],
         '%.1f %% of the volume'),
        ('0.5 %', 100 * smear, '%.1f %%'),
    ]
    for text, value, fmt in claims:
        assert text in doc, 'documentation no longer says %r' % text
        assert (fmt % value) == text, (
            'documentation says %r; the engine says %r' % (text, fmt % value))

    # the GGA decay length, and the two branch numbers, in the table
    assert '0.725 nm' in doc
    assert '%.3f nm' % M.segregation(p, 'gga')['xi_nm'] == '0.725 nm'
    lo = PR.inventory(flagship['mu_transition_eV'], flagship['kt_eV'], part,
                      surf_fraction=part['surface_phases'][
                          'theta_transition'])
    hi = PR.inventory(flagship['mu_transition_eV'], flagship['kt_eV'], part,
                      surf_fraction=part['surface_phases'][
                          'theta_reconstructed_eff'])
    assert '33.2 → 37.7 μmol-O/g' in doc
    assert '%.1f → %.1f μmol-O/g' % (lo, hi) == '33.2 → 37.7 μmol-O/g'


def test_the_command_line_report_still_runs():
    import subprocess
    for args in ([], ['--ladder']):
        out = subprocess.run(
            ['python3', os.path.join(ROOT, 'scripts', 'run_particle.py')]
            + args, capture_output=True, text=True, cwd=ROOT)
        assert out.returncode == 0, out.stderr
        assert 'umol-O/g' in out.stdout
