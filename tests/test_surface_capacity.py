"""The surface-capacity bound: is it arithmetic, and is it a bound?

The claim this module supports is an inequality, so the gates check the
inequality and the geometry it is built from, not a fitted number.
"""

import json
import math
import pathlib

import pytest

from solidgas import vacancy_inventory as V

ROOT = pathlib.Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def g():
    return V.load()


# ------------------------------------------------------------- the geometry

def test_the_cell_volume_reproduces_the_density_of_rutile(g):
    """4.25 g/cm3 is a tabulated number this file never stores."""
    assert V.crystal_density_g_cm3(g) == pytest.approx(4.25, abs=0.01)


def test_the_110_cell_is_the_one_the_literature_quotes(g):
    """c by a*sqrt(2) is 19.22 square angstrom."""
    assert V.surface_cell_area_cm2(g) * 1e16 == pytest.approx(19.22, abs=0.01)


def test_one_bridging_oxygen_per_cell_sets_the_areal_density(g):
    area = V.surface_cell_area_cm2(g)
    assert V.bridging_areal_density_cm2(g) == pytest.approx(1.0 / area, rel=1e-12)


def test_the_oxygen_site_density_is_two_per_formula_unit(g):
    assert V.oxygen_site_density_umol_g(g) == pytest.approx(
        2.0 / 79.866 * 1e6, rel=1e-9)


# --------------------------------------------------------------- the area

def test_a_smaller_particle_exposes_more_area(g):
    assert (V.specific_surface_area(g, d_um=0.9)
            > V.specific_surface_area(g, d_um=1.6))


def test_area_scales_as_one_over_diameter(g):
    a, b = V.specific_surface_area(g, d_um=0.9), V.specific_surface_area(g, d_um=1.8)
    assert a / b == pytest.approx(2.0, rel=1e-12)


def test_a_measured_bet_area_overrides_the_sphere_estimate(g):
    assert V.specific_surface_area(g, d_um=0.9, bet_m2_g=5.0) == 5.0


def test_it_refuses_to_guess_an_area(g):
    with pytest.raises(ValueError):
        V.specific_surface_area(g)


# ------------------------------------------------------------- the capacity

def test_the_reference_capacities_are_what_the_manuscript_quotes(g):
    assert V.bridging_oxygen_capacity(g, d_um=0.9) == pytest.approx(13.56, abs=0.01)
    assert V.reconstruction_deficiency(g, d_um=0.9) == pytest.approx(6.78, abs=0.01)
    assert V.bridging_oxygen_capacity(g, d_um=1.6) == pytest.approx(7.63, abs=0.01)


def test_the_reconstructed_surface_holds_exactly_half_the_row(g):
    for d in (0.9, 1.6):
        assert (V.reconstruction_deficiency(g, d_um=d)
                == pytest.approx(0.5 * V.bridging_oxygen_capacity(g, d_um=d),
                                 rel=1e-12))


def test_the_surface_is_a_vanishing_share_of_the_solid(g):
    """A micron particle has of order a ten-thousandth of its oxygen on top.

    This is the whole reason the bound bites, and it is why a slab ratio
    must never be used in its place.
    """
    frac = (V.bridging_oxygen_capacity(g, d_um=0.9)
            / V.oxygen_site_density_umol_g(g))
    assert frac < 1e-3


# ---------------------------------------------------------------- the bound

def test_the_two_fractions_are_complements(g):
    for inv in (30.0, 95.0, 220.0):
        b = V.surface_fraction_bound(g, inv, d_um=0.9)
        assert (b['surface_fraction_max'] + b['below_surface_fraction_min']
                == pytest.approx(1.0, abs=1e-12))
        assert (b['surface_fraction_max_full_row']
                + b['below_surface_fraction_min_full_row']
                == pytest.approx(1.0, abs=1e-12))


def test_the_bound_falls_as_more_oxygen_comes_out(g):
    got = [V.surface_fraction_bound(g, v, d_um=0.9)['surface_fraction_max']
           for v in (30.0, 60.0, 95.0, 130.0, 190.0, 220.0)]
    assert all(b < a for a, b in zip(got, got[1:])), got


def test_the_full_row_bound_is_looser_than_the_reconstructed_one(g):
    """Stripping the whole row is the weakest claim, so it must bound less."""
    b = V.surface_fraction_bound(g, 95.0, d_um=0.9)
    assert b['surface_fraction_max_full_row'] > b['surface_fraction_max']


def test_the_reference_bound_is_what_the_manuscript_would_quote(g):
    b = V.surface_fraction_bound(g, 95.0, d_um=0.9)
    assert b['inventory_over_capacity'] == pytest.approx(7.01, abs=0.01)
    assert b['surface_fraction_max'] * 100 == pytest.approx(7.14, abs=0.01)
    assert b['below_surface_fraction_min'] * 100 == pytest.approx(92.86, abs=0.01)


def test_a_bigger_particle_makes_the_bound_tighter(g):
    small = V.surface_fraction_bound(g, 95.0, d_um=0.9)['surface_fraction_max']
    big = V.surface_fraction_bound(g, 95.0, d_um=1.6)['surface_fraction_max']
    assert big < small


def test_the_fractions_never_exceed_one(g):
    """An inventory below the surface capacity is allowed and must not
    report more than all of itself on the surface."""
    b = V.surface_fraction_bound(g, 1.0, d_um=0.9)
    assert b['surface_fraction_max'] == 1.0
    assert b['below_surface_fraction_min'] == 0.0


def test_it_refuses_a_nonpositive_inventory(g):
    for bad in (0.0, -5.0):
        with pytest.raises(ValueError):
            V.surface_fraction_bound(g, bad, d_um=0.9)


# --------------------------------------------------------------- conversions

def test_x_and_the_ti3_fraction_agree_with_each_other(g):
    """One vacancy leaves two Ti(3+), so the Ti(3+) fraction is 2x."""
    for inv in (30.0, 95.0, 220.0):
        x = V.inventory_to_x(g, inv)
        assert V.inventory_to_ti3_fraction(g, inv) == pytest.approx(2 * x, rel=1e-12)


def test_the_reference_composition_is_near_the_shear_range(g):
    """95 umol-O/g is x = 0.0076, which is where rutile stops dissolving
    point defects. The number is quoted; nothing here acts on it."""
    assert V.inventory_to_x(g, 95.0) == pytest.approx(0.00759, abs=1e-5)


# ------------------------------------------------------------- what it isn't

def test_no_defect_energy_reaches_this_calculation():
    """The bound must stay arithmetic. If a DFT number ever gets wired in,
    it stops being a bound and becomes a model, and this catches that."""
    src = (ROOT / 'solidgas' / 'vacancy_inventory.py').read_text()
    for word in ('literature_dft', 'formation_energy', 'eV', 'chemical_potential'):
        assert word not in src, word
    geo = json.loads((ROOT / 'data' / 'rutile_geometry.json').read_text())
    assert 'eV' not in json.dumps(geo)


def test_the_only_structural_input_is_declared(g):
    """The 0.5 ML deficiency is a literature value and the file says so."""
    assert g['reconstruction']['areal_O_deficiency_ML'] == 0.5
    assert 'Onishi' in g['reconstruction']['provenance']
