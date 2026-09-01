"""How much oxygen the rutile surface can hold, against how much was removed.

A titration measures one number: total oxygen removed, in umol-O per gram.
It does not say where that oxygen came from. This module answers the one
question the measurement can support on its own - whether the surface
could have held it - and it answers it with lattice geometry and nothing
else. No defect energies, no chemical potential, no partition.

The answer is a BOUND, not an estimate:

    f_surface <= (surface oxygen capacity) / (measured inventory)

and its complement is a floor on how much of the inventory must lie below
the surface. A bound cannot be wrong about the direction of the inequality
the way a point estimate can be wrong about a number. Assigning the
below-surface inventory among subsurface point defects, bulk defects and
extended defects needs measurements this work does not have, and is not
attempted here or anywhere else in the package.

Two capacities are worth quoting, and the honest presentation is both:

  the full bridging row, every two-coordinate oxygen on (110) removed,
  which no experiment has ever reached and which needs no assumption at
  all beyond the lattice;

  the (1x2) added-row deficiency, half a bridging monolayer, which is
  what the reconstructed surface actually holds. That 0.5 is a structural
  literature value and it is the only non-geometric input here.

Three-coordinate in-plane oxygen is not counted. It is not the oxygen the
gas takes, and counting it would inflate the capacity and weaken the
bound in the direction that flatters the conclusion.
"""

import json
import os

N_AVO = 6.02214076e23
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(_ROOT, 'data', 'rutile_geometry.json')


def load(path=DATA_PATH):
    with open(path) as fh:
        return json.load(fh)


# ------------------------------------------------------------------ geometry

def surface_cell_area_cm2(g):
    """Area of the (110) 1x1 surface cell: c by a*sqrt(2)."""
    a = g['lattice_constants_A']['a'] * 1e-8
    c = g['lattice_constants_A']['c'] * 1e-8
    return a * (2.0 ** 0.5) * c


def bridging_areal_density_cm2(g):
    """Two-coordinate bridging oxygen per square centimetre of (110)."""
    return g['surface_cell']['bridging_O_per_1x1_cell'] / surface_cell_area_cm2(g)


def crystal_density_g_cm3(g):
    """Rutile density from the unit cell, not from a table."""
    a = g['lattice_constants_A']['a'] * 1e-8
    c = g['lattice_constants_A']['c'] * 1e-8
    volume = a * a * c
    return g['formula_units_per_cell'] * g['molar_mass_g_mol'] / (N_AVO * volume)


def specific_surface_area(g, d_um=None, bet_m2_g=None):
    """Area per gram, in m2/g. A measured BET area wins when both are given.

    From a diameter it is the smooth-sphere value 6/(rho d), which is a
    lower bound on the real area of a rough particle - so it is also the
    conservative choice for the capacity bound below.
    """
    if bet_m2_g is not None:
        return float(bet_m2_g)
    if d_um is None:
        raise ValueError('give a BET area or a particle diameter')
    d_cm = d_um * 1e-4
    return 6.0 / (crystal_density_g_cm3(g) * d_cm) * 1e-4


# ------------------------------------------------------------------ capacity

def bridging_oxygen_capacity(g, d_um=None, bet_m2_g=None):
    """Every bridging oxygen on the exposed (110), in umol-O per gram."""
    area_cm2_g = specific_surface_area(g, d_um, bet_m2_g) * 1e4
    return bridging_areal_density_cm2(g) * area_cm2_g / N_AVO * 1e6


def reconstruction_deficiency(g, d_um=None, bet_m2_g=None):
    """Oxygen held by the fully (1x2)-reconstructed surface, umol-O per gram.

    Half a bridging monolayer, the added-row Ti2O3 stoichiometry.
    """
    ml = g['reconstruction']['areal_O_deficiency_ML']
    return ml * bridging_oxygen_capacity(g, d_um, bet_m2_g)


# ----------------------------------------------------------------- inventory

def oxygen_site_density_umol_g(g):
    """All oxygen in the solid, umol-O per gram. 2 / M for a dioxide."""
    return g['oxygen_per_formula_unit'] / g['molar_mass_g_mol'] * 1e6


def inventory_to_x(g, umol_o_g):
    """x in TiO(2-x) for a measured oxygen removal."""
    return umol_o_g / oxygen_site_density_umol_g(g) * g['oxygen_per_formula_unit']


def inventory_to_ti3_fraction(g, umol_o_g):
    """Fraction of titanium reduced to Ti(3+). One vacancy leaves two."""
    n_ti = 1.0 / g['molar_mass_g_mol'] * 1e6
    return min(2.0 * umol_o_g / n_ti, 1.0)


def surface_fraction_bound(g, umol_o_g, d_um=None, bet_m2_g=None):
    """The bound, with both capacities and everything they are read from.

    `surface_max` uses the reconstructed-surface deficiency, which is the
    physically reachable ceiling; `surface_max_full_row` strips the whole
    bridging row, which nothing has been observed to do and which is the
    weakest claim that needs no structural input at all.
    """
    if not umol_o_g > 0:
        raise ValueError('the measured inventory must be positive')
    full = bridging_oxygen_capacity(g, d_um, bet_m2_g)
    recon = reconstruction_deficiency(g, d_um, bet_m2_g)
    return {
        'measured_umol_o_g': float(umol_o_g),
        'area_m2_g': specific_surface_area(g, d_um, bet_m2_g),
        'density_g_cm3': crystal_density_g_cm3(g),
        'bridging_capacity_umol_o_g': full,
        'reconstruction_deficiency_umol_o_g': recon,
        'inventory_over_capacity': umol_o_g / full,
        'surface_fraction_max': min(recon / umol_o_g, 1.0),
        'below_surface_fraction_min': max(1.0 - recon / umol_o_g, 0.0),
        'surface_fraction_max_full_row': min(full / umol_o_g, 1.0),
        'below_surface_fraction_min_full_row': max(1.0 - full / umol_o_g, 0.0),
        'x_in_TiO2mx': inventory_to_x(g, umol_o_g),
        'Ti3_fraction': inventory_to_ti3_fraction(g, umol_o_g),
        'basis': ('lattice geometry; the (1x2) 0.5 ML deficiency is the only '
                  'structural literature input'),
    }
