"""The material card: what the particle is made of, in one place.

Everything here is either read from data/rutile_dft.json or derived from
it by closed formula.  Nothing is fitted and nothing is invented.

The one derivation worth reading is the segregation decay length.  The
DFT set gives three vacancy formation energies - bridging row, first
subsurface oxygen layer, bulk - and the profile

    E(z) = E_bulk - dE * exp(-z / xi)

has exactly two free parameters once E_bulk is pinned as its asymptote.
Two parameters, two remaining data points: dE = E_bulk - E_surface from
z = 0, and then xi from the subsurface point at z = d110,

    xi = -d110 / ln( (E_bulk - E_sub) / dE ).

So the depth profile has no shape assumption left in it.  The previous
engine placed intermediate layers at a declared fraction of the
subsurface-to-bulk interval, and that fraction was labelled "assumed" in
the dataset; the exponential reproduces it to within a few hundredths
without being told to.
"""

import json
import math
import os

KB_EV = 8.617333262e-5          # eV / K
N_AVO = 6.02214076e23

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
DATA_PATH = os.path.join(_ROOT, 'data', 'rutile_dft.json')


def load(path=DATA_PATH):
    with open(path) as fh:
        return json.load(fh)


def segregation(p, preset=None):
    """Depth profile of the vacancy formation energy, from the DFT triple.

    Returns the three input energies, the surface-to-bulk gap, and the
    decay length that carries one into the other through the subsurface
    point.  Raises if the triple is not monotone towards the bulk, in
    which case a single decaying exponential is the wrong ansatz and
    saying so is better than fitting something else."""
    ve = p['vacancy_energetics']
    key = preset or ve['default_preset']
    q = ve['presets'][key]
    e_s, e_ss, e_b = q['surface_eV'], q['subsurface_eV'], q['bulk_eV']
    d110_nm = p['lattice_constants_A']['a'] / math.sqrt(2.0) / 10.0

    de = e_b - e_s
    rem = e_b - e_ss
    if not (de > 0.0 and 0.0 < rem < de):
        raise ValueError(
            'the %s triple (%.3f / %.3f / %.3f eV) does not decay '
            'monotonically towards the bulk, so E(z) = E_bulk - dE '
            'exp(-z/xi) cannot pass through all three' % (key, e_s, e_ss, e_b))
    xi_nm = -d110_nm / math.log(rem / de)

    return {'preset': key, 'label': q.get('label', key),
            'e_surface_eV': e_s, 'e_subsurface_eV': e_ss, 'e_bulk_eV': e_b,
            'dE_seg_eV': de, 'xi_nm': xi_nm, 'd110_nm': d110_nm,
            'xi_in_layers': xi_nm / d110_nm,
            'provenance': ve.get('provenance', '')}


def formation_energy(z_nm, seg):
    """E(z) in eV at depth z below the bridging plane.  z >= 0."""
    return seg['e_bulk_eV'] - seg['dE_seg_eV'] * math.exp(
        -max(0.0, z_nm) / seg['xi_nm'])


def density_g_cm3(p):
    a = p['lattice_constants_A']['a']
    c = p['lattice_constants_A']['c']
    return 2.0 * p['molar_mass_g_mol'] / (N_AVO * a * a * c * 1e-24)


def bridging_density_cm2(p):
    """Areal density of bridging oxygens on (110): one per 1x1 cell."""
    a = p['lattice_constants_A']['a']
    c = p['lattice_constants_A']['c']
    return 1e16 / (c * a * math.sqrt(2.0))


def particle(p, d_um=None, bet_m2_g=None, preset=None):
    """Geometry and site budget for one spherical particle, per gram.

    Oxygen is counted in umol-O/g throughout, which is the unit the
    titration measures in.  Two site pools:

      surface  - the bridging row, one oxygen per (1x1) cell, an areal
                 density rather than a volume one, and the site the
                 4.39 eV number belongs to;
      interior - everything below the bridging plane, at the bulk oxygen
                 density, carrying a continuous depth profile.

    Splitting them is not decoration.  A full d110 oxygen slab holds
    four oxygens per cell and the bridging row holds one, so smearing
    the row into a uniform continuum would give it four times the sites
    it has, and it is the cheapest site in the particle."""
    if bet_m2_g is not None:
        area_cm2_g = bet_m2_g * 1e4
    else:
        if d_um is None:
            d_um = p['defaults']['particle_diameter_um']
        area_cm2_g = 6.0 / (density_g_cm3(p) * d_um * 1e-4)

    n_o = 2.0 / p['molar_mass_g_mol'] * 1e6              # umol-O/g, all of it
    n_surface = area_cm2_g * bridging_density_cm2(p) / N_AVO * 1e6
    seg = segregation(p, preset)

    # radius that reproduces this specific surface-to-volume ratio, so
    # the sphere the profile is integrated over is the same object the
    # site budget came from
    r_nm = 3.0 / (density_g_cm3(p) * area_cm2_g) * 1e7

    sat = p.get('saturation', {})
    x_shear = sat.get('x_max_shear')
    ph = p.get('surface_phases', {})
    return {'area_m2_g': area_cm2_g * 1e-4,
            'radius_nm': r_nm,
            'diameter_um': 2.0 * r_nm * 1e-3,
            'N_surface_umol_g': n_surface,
            'N_interior_umol_g': n_o - n_surface,
            'N_oxygen_umol_g': n_o,
            'segregation': seg,
            'x_max_shear': x_shear,
            'theta_solubility': None if x_shear is None else 0.5 * x_shear,
            'surface_phases': {
                'theta_transition': ph.get('theta_transition'),
                'theta_reconstructed_eff': ph.get('theta_reconstructed_eff'),
                'transition_interval': ph.get('transition_interval')}}


def theta_from_loading(vo_umol_g, part):
    """Particle-average oxygen site fraction from a titration number.

    theta counts vacant oxygen SITES, so TiO(2-x) has x = 2 theta: two
    oxygen sites per formula unit.  Mixing the two conventions is the
    single easiest way to be wrong by a factor of two here."""
    return vo_umol_g / part['N_oxygen_umol_g']


def loading_from_theta(theta, part):
    return theta * part['N_oxygen_umol_g']


def kt_eV(t_c):
    return KB_EV * (t_c + 273.15)
