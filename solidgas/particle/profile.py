"""Equilibrium vacancy profile in one spherical particle.

The particle is a sphere of radius R.  Oxygen sits in two pools:

  the bridging row, an areal density on the surface, the cheapest site
  in the particle and the one the surface DFT number belongs to;

  the interior, at the bulk oxygen density, with a formation energy that
  relaxes towards the bulk value over the decay length derived in
  material.py.

One chemical potential mu is shared by both - that is what equilibrium
means - and mu is fixed by the one thing an experiment actually measures,
the total oxygen removed.  So the calculation is: find the mu whose
pools together hold the titrated inventory.

Layers, not a smear.  The interior is a stack of oxygen layers d110
apart, layer k centred at z = k*d110, and the continuum integral used
here is the limit of that sum.  It starts at z = d110/2 so that layer
one occupies [d110/2, 3 d110/2] and gets the subsurface energy the DFT
set assigns it.  Both forms are implemented and the gates hold them
against each other; the continuum is the one that runs fast, the sum is
the one that is obviously right.

Two things happen in rutile that do not happen in the ceria system this
construction is taken from, and both cut the profile off before the
cation sublattice fills:

  the (110) face reconstructs.  Past a bridging-row coverage of about
  0.17 the (1x1) vacancy lattice gas gives way to the added-row Ti2O3
  (1x2) phase, a first-order transition, so mu pins while the surface
  converts;

  the bulk throws crystallographic shear planes.  Past x = 0.008 the
  point-defect solution is not what forms - Magneli phases are - so mu
  pins again and further inventory leaves the point-defect count as
  extended defects.
"""

import math

from . import freeenergy as F
from . import material as M

# Gauss-Legendre, 10 points on [-1, 1].  The near-surface integrand is a
# smooth exponential over each graded panel, so this is exact to machine
# precision and there is no grid to tune.
_GL_X = [-0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
         -0.4333953941292472, -0.1488743389816312, 0.1488743389816312,
         0.4333953941292472, 0.6794095682990244, 0.8650633666889845,
         0.9739065285171717]
_GL_W = [0.0666713443086881, 0.1494513491505806, 0.2190863625159820,
         0.2692667193099963, 0.2955242247147529, 0.2955242247147529,
         0.2692667193099963, 0.2190863625159820, 0.1494513491505806,
         0.0666713443086881]

# Panel edges in units of xi.  Graded, because everything happens in the
# first few xi and nothing happens after twenty.
_PANELS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]


def _theta(mu, e, kt, compensated, ratio):
    return F.theta_of_mu(mu, e, kt, compensated, ratio)


def _layer_geometry(part):
    """Radial shells of the interior, one per d110 oxygen layer.

    Layer k is centred at depth k*d110 below the bridging plane and
    occupies the shell between depths (k-1/2) d110 and (k+1/2) d110.
    The shells tile the interior sphere exactly, so their volumes need
    no normalising."""
    seg = part['segregation']
    d = seg['d110_nm']
    r_in = part['radius_nm'] - 0.5 * d
    return d, r_in, r_in ** 3 / 3.0


def interior_average(mu, kt, part, compensated=True,
                     ratio=F.POLARON_RATIO, rel_tol=1e-15):
    """Volume-averaged interior site fraction at chemical potential mu.

    Exact, not quadrature.  The interior is a finite stack of oxygen
    layers, and past a few tens of decay lengths every layer is at the
    bulk energy to the last bit of a double, so the sum splits into an
    explicit near-surface part and a closed-form remainder:

        <theta> = theta_bulk + sum_k V_k (theta_k - theta_bulk) / V.

    The sum stops when a layer's contribution falls below rel_tol of the
    running total, which for the rutile numbers is about seventy layers
    - twenty-three nanometres of a four-hundred-and-fifty nanometre
    particle."""
    seg = part['segregation']
    d, r_in, vol = _layer_geometry(part)
    th_bulk = _theta(mu, seg['e_bulk_eV'], kt, compensated, ratio)

    excess, k = 0.0, 1
    n_max = int(math.floor(r_in / d))
    while k <= n_max:
        ro = r_in - (k - 1) * d
        ri = max(0.0, r_in - k * d)
        v = (ro ** 3 - ri ** 3) / 3.0
        e = M.formation_energy(k * d, seg)
        term = v * (_theta(mu, e, kt, compensated, ratio) - th_bulk)
        excess += term
        if abs(term) <= rel_tol * max(abs(excess), 1e-300) and k > 4:
            break
        k += 1
    return th_bulk + excess / vol


def interior_average_smeared(mu, kt, part, compensated=True,
                             ratio=F.POLARON_RATIO):
    """The same average with the layers smeared into a continuum.

    This is the form the ceria study integrates, and it is not used to
    run: it is kept so the cost of the smear can be quoted rather than
    assumed away.  A layer stack whose energy changes by a large factor
    from one layer to the next is not well represented by its own
    average, and at the rutile decay length the two differ by about half
    a percent in the interior mean."""
    seg = part['segregation']
    xi, z0 = seg['xi_nm'], 0.5 * seg['d110_nm']
    d, r_in, vol = _layer_geometry(part)
    th_bulk = _theta(mu, seg['e_bulk_eV'], kt, compensated, ratio)

    excess = 0.0
    for a, b in zip(_PANELS[:-1], _PANELS[1:]):
        za, zb = a * xi, b * xi
        if za >= r_in:
            break
        zb = min(zb, r_in)
        half, mid = 0.5 * (zb - za), 0.5 * (zb + za)
        for x, w in zip(_GL_X, _GL_W):
            z = mid + half * x
            e = M.formation_energy(z + z0, seg)
            dth = _theta(mu, e, kt, compensated, ratio) - th_bulk
            excess += w * half * dth * (r_in - z) ** 2
    return th_bulk + excess / vol


def surface_state(mu, kt, part, compensated=True, ratio=F.POLARON_RATIO,
                  reconstruction=True):
    """Bridging-row coverage at mu, and how much of it is reconstructed.

    Below the transition the row is a (1x1) lattice gas.  At mu_t the
    added-row (1x2) phase appears; it is a line compound with a fixed
    areal oxygen deficiency of 0.5 in bridging-monolayer units, so it
    has no configurational entropy of its own and coexistence happens at
    a single mu.  Above mu_t the row is fully reconstructed and its
    coverage stops responding."""
    seg = part['segregation']
    e_s = seg['e_surface_eV']
    th_11 = _theta(mu, e_s, kt, compensated, ratio)
    if not reconstruction:
        return {'theta': th_11, 'phase': '1x1', 'reconstructed_fraction': 0.0,
                'mu_transition_eV': None}

    ph = part.get('surface_phases') or {}
    th_t = ph.get('theta_transition')
    th_r = ph.get('theta_reconstructed_eff')
    if th_t is None or th_r is None:
        return {'theta': th_11, 'phase': '1x1', 'reconstructed_fraction': 0.0,
                'mu_transition_eV': None}
    mu_t = F.mu_of_theta(th_t, e_s, kt, compensated, ratio)
    if mu < mu_t:
        return {'theta': th_11, 'phase': '1x1',
                'reconstructed_fraction': 0.0, 'mu_transition_eV': mu_t}
    if mu > mu_t:
        return {'theta': th_r, 'phase': '1x2',
                'reconstructed_fraction': 1.0, 'mu_transition_eV': mu_t}
    return {'theta': th_t, 'phase': 'coexistence',
            'reconstructed_fraction': 0.0, 'mu_transition_eV': mu_t}


def inventory(mu, kt, part, compensated=True, ratio=F.POLARON_RATIO,
              reconstruction=True, surf_fraction=None):
    """Total oxygen held as point defects at mu, in umol-O/g.

    surf_fraction overrides the surface coverage during two-phase
    coexistence, where mu alone does not determine it."""
    s = surface_state(mu, kt, part, compensated, ratio, reconstruction)
    th_s = s['theta'] if surf_fraction is None else surf_fraction
    th_i = interior_average(mu, kt, part, compensated, ratio)
    return (part['N_surface_umol_g'] * th_s
            + part['N_interior_umol_g'] * th_i)


def _bracket(target, kt, part, compensated, ratio, reconstruction):
    """A mu interval that contains the answer, from the dilute estimate."""
    seg = part['segregation']
    th = target / part['N_oxygen_umol_g']
    th = min(max(th, 1e-300), 0.5 / ratio)
    guess = seg['e_bulk_eV'] + kt * (3.0 * math.log(th)
                                     + 2.0 * math.log(ratio))
    lo, hi = guess - 1.0, guess + 1.0
    for _ in range(60):
        if inventory(lo, kt, part, compensated, ratio, reconstruction) <= target:
            break
        lo -= (hi - lo)
    for _ in range(60):
        if inventory(hi, kt, part, compensated, ratio, reconstruction) >= target:
            break
        hi += (hi - lo)
    return lo, hi


def equilibrate(vo_umol_g, t_c, part, compensated=True,
                ratio=F.POLARON_RATIO, reconstruction=True,
                shear_limit=True):
    """The whole calculation: titrated inventory in, profile out.

    Returns mu, where the oxygen sits, and - separately - whatever the
    point-defect phases cannot hold.  Inventory past the shear limit is
    reported as extended defects rather than being pushed into a site
    fraction that would not survive the crystallography."""
    kt = M.kt_eV(t_c)
    seg = part['segregation']
    cap_total = part['N_oxygen_umol_g'] / ratio
    if vo_umol_g > cap_total:
        raise ValueError(
            'a loading of %.4g umol-O/g is past the cation sublattice '
            'itself (%.4g umol-O/g is theta = 1/%g everywhere, which is '
            'Ti2O3, not reduced rutile)' % (vo_umol_g, cap_total, ratio))

    mu_cs = None
    if shear_limit and part.get('theta_solubility'):
        mu_cs = F.mu_of_theta(part['theta_solubility'], seg['e_bulk_eV'],
                              kt, compensated, ratio)

    ph = part.get('surface_phases') or {}
    th_t = ph.get('theta_transition') if reconstruction else None
    th_r = ph.get('theta_reconstructed_eff') if reconstruction else None
    mu_t = None if th_t is None else F.mu_of_theta(
        th_t, seg['e_surface_eV'], kt, compensated, ratio)

    mu, limited_by, extended, th_s = None, None, 0.0, None

    # The added-row phase is a line compound, so the surface converts at
    # one chemical potential.  Inside that window mu does not respond to
    # inventory at all and bisecting on it would land on either edge, so
    # the window is tested for first and the lever rule does the rest.
    if mu_t is not None and (mu_cs is None or mu_t <= mu_cs):
        lo_inv = inventory(mu_t, kt, part, compensated, ratio,
                           reconstruction, surf_fraction=th_t)
        hi_inv = inventory(mu_t, kt, part, compensated, ratio,
                           reconstruction, surf_fraction=th_r)
        if lo_inv <= vo_umol_g <= hi_inv:
            mu = mu_t
            frac = ((vo_umol_g - lo_inv) / (hi_inv - lo_inv)
                    if hi_inv > lo_inv else 0.0)
            th_s = th_t + frac * (th_r - th_t)
            limited_by = 'surface reconstruction'

    if mu is None:
        lo, hi = _bracket(vo_umol_g, kt, part, compensated, ratio,
                          reconstruction)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if inventory(mid, kt, part, compensated, ratio,
                         reconstruction) < vo_umol_g:
                lo = mid
            else:
                hi = mid
        mu = 0.5 * (lo + hi)

    if mu_cs is not None and mu > mu_cs:
        mu, th_s, limited_by = mu_cs, None, 'crystallographic shear'
        extended = vo_umol_g - inventory(mu, kt, part, compensated, ratio,
                                         reconstruction)

    surf = surface_state(mu, kt, part, compensated, ratio, reconstruction)
    if th_s is not None:
        surf = dict(surf, theta=th_s, phase='coexistence',
                    reconstructed_fraction=frac)
    th_i = interior_average(mu, kt, part, compensated, ratio)
    held_s = part['N_surface_umol_g'] * surf['theta']
    held_i = part['N_interior_umol_g'] * th_i

    return {
        'mu_eV': mu, 'kt_eV': kt, 'T_C': t_c,
        'compensated': compensated, 'polaron_ratio': ratio,
        'loading_umol_g': vo_umol_g,
        'theta_particle': vo_umol_g / part['N_oxygen_umol_g'],
        'surface': surf,
        'theta_interior_mean': th_i,
        'theta_bulk_plateau': _theta(mu, seg['e_bulk_eV'], kt,
                                     compensated, ratio),
        'theta_first_layer': _theta(
            mu, M.formation_energy(seg['d110_nm'], seg), kt,
            compensated, ratio),
        'umol_surface': held_s, 'umol_interior': held_i,
        'umol_extended': extended,
        'held_umol_g': held_s + held_i,
        'limited_by': limited_by,
        'mu_shear_eV': mu_cs, 'mu_transition_eV': mu_t,
        'point_defect_ceiling_umol_g': (
            None if mu_cs is None else inventory(
                mu_cs, kt, part, compensated, ratio, reconstruction)),
    }


def depth_profile(res, part, n=160, z_max_nm=None):
    """theta against depth, for drawing.  Log-graded, because the whole
    structure lives inside the first few nanometres of a 450 nm sphere."""
    seg = part['segregation']
    kt = res['kt_eV']
    zmax = z_max_nm if z_max_nm is not None else min(
        part['radius_nm'], 40.0 * seg['xi_nm'])
    out = [{'z_nm': 0.0, 'theta': res['surface']['theta'],
            'pool': 'bridging row'}]
    for i in range(n):
        z = 0.5 * seg['d110_nm'] + (zmax - 0.5 * seg['d110_nm']) * (
            (i / (n - 1.0)) ** 3)
        out.append({'z_nm': z,
                    'theta': _theta(res['mu_eV'],
                                    M.formation_energy(z, seg), kt,
                                    res['compensated'],
                                    res['polaron_ratio']),
                    'pool': 'interior'})
    return out


def shell_split(res, part, depths_nm=(1.0, 2.0, 5.0, 10.0)):
    """How much of the inventory sits within each depth of the surface.

    The point of the whole exercise: a shell a few nanometres deep on a
    micron particle is a rounding error in volume and can still hold a
    large share of the oxygen.  Counted over the same layer stack the
    solver uses, so the shares add up to the solved inventory."""
    seg = part['segregation']
    kt, mu = res['kt_eV'], res['mu_eV']
    d, r_in, vol = _layer_geometry(part)
    total = res['held_umol_g']
    rows = []
    for dep in depths_nm:
        if dep >= r_in:
            continue
        acc, k = 0.0, 1
        while k <= int(math.floor(r_in / d)):
            z_out = (k - 1) * d
            if z_out >= dep:
                break
            ro = r_in - z_out
            ri = max(r_in - min(k * d, dep), 0.0)
            v = (ro ** 3 - ri ** 3) / 3.0
            acc += v * _theta(mu, M.formation_energy(k * d, seg), kt,
                              res['compensated'], res['polaron_ratio'])
            k += 1
        umol = (res['umol_surface']
                + part['N_interior_umol_g'] * acc / vol)
        rows.append({'depth_nm': dep,
                     'volume_fraction': 1.0 - ((r_in - dep) / r_in) ** 3,
                     'umol_g': umol,
                     'inventory_fraction': umol / total if total else 0.0})
    return rows
