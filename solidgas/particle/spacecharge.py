"""Is local electroneutrality allowed here?

The radial profile assumes each vacancy drags its two polarons with it,
pointwise: theta_polaron(r) = 4 theta_V(r) at every radius.  That is the
ceria study's Model 3.  Its Model 4 drops the assumption, lets the two
species equilibrate independently in a self-consistent electrostatic
potential, and asks what changes.  For the concentrated ceria case the
answer was: almost nothing, over more than ninety-five percent of the
radius, because the carrier density was so high that screening was
essentially instantaneous.  The report is explicit that this collapse is
a consequence of the concentration and that a more dilute system would
show a wider space-charge layer.

That is the direction this system lies in.  The rutile particle at the
shipped loading runs an order of magnitude more dilute than the
twenty-percent-reduced ceria, so the comparison has to be computed
rather than inherited.

The computation that settles it is one line.  Screening happens over the
Debye length

    lambda_D = sqrt( eps_r eps_0 kT / sum_i z_i^2 e^2 n_i ),

chemical segregation happens over xi, and the question is which is
shorter.  If lambda_D << xi, charge follows chemistry and pointwise
electroneutrality is safe.  If they are comparable, the two cannot be
separated and the profile's outer few nanometres carry an unquantified
error in the direction of more surface enrichment, not less.
"""

import math

from .material import KB_EV, N_AVO

E_CHARGE = 1.602176634e-19       # C
EPS0_F_PER_CM = 8.8541878128e-14  # F/cm

# Static dielectric constant of rutile is large and strongly anisotropic:
# about 86 perpendicular to c and 170 along it at room temperature, and
# it falls with temperature.  A hundred is the round number in the middle
# and the length below goes as its square root, so the choice moves the
# answer by tens of percent, not by orders of magnitude.
EPS_R_RUTILE = 100.0


def oxygen_density_cm3(p):
    """Oxygen atoms per cm^3 in perfect rutile."""
    a = p['lattice_constants_A']['a']
    c = p['lattice_constants_A']['c']
    return 4.0 / (a * a * c * 1e-24)


def debye_length_nm(theta_v, p, t_c, eps_r=EPS_R_RUTILE,
                    ratio=4.0):
    """Screening length at a local vacancy site fraction.

    Both charged species count, each weighted by the square of its
    charge: the vacancy carries +2 and there are theta_v of them per
    oxygen site, the polaron carries -1 and there are ratio*theta_v of
    them per cation site, with half as many cation sites as oxygen
    sites.  Their contributions are 4*theta_v and ratio*theta_v/2 times
    the oxygen density."""
    n_o = oxygen_density_cm3(p)
    kt_J = KB_EV * (t_c + 273.15) * E_CHARGE
    z2n = (4.0 * theta_v + 0.5 * ratio * theta_v) * n_o
    if z2n <= 0.0:
        return float('inf')
    lam_cm = math.sqrt(eps_r * EPS0_F_PER_CM * kt_J
                       / (E_CHARGE * E_CHARGE * z2n))
    return lam_cm * 1e7


def assess(res, part, p, eps_r=EPS_R_RUTILE):
    """Compare screening against segregation, and say what it licenses.

    Reported at the bulk plateau, which is the most dilute part of the
    particle and therefore where screening is weakest and the assumption
    is under the most strain."""
    seg = part['segregation']
    xi = seg['xi_nm']
    th = res['theta_bulk_plateau']
    lam = debye_length_nm(th, p, res['T_C'], eps_r, res['polaron_ratio'])
    lam_surf = debye_length_nm(max(res['theta_first_layer'], 1e-30), p,
                               res['T_C'], eps_r, res['polaron_ratio'])
    r = lam / xi
    if r < 0.2:
        verdict = ('screening is much shorter than the segregation depth, '
                   'so charge follows chemistry and pointwise '
                   'electroneutrality is safe')
    elif r < 2.0:
        verdict = ('screening and segregation happen over the same '
                   'distance, so they cannot be separated; the outer few '
                   'nanometres of this profile carry an error that this '
                   'model does not quantify, and the ceria study finds '
                   'that error runs towards MORE surface enrichment')
    else:
        verdict = ('screening is longer than the segregation depth, so a '
                   'genuine space-charge layer extends past the '
                   'chemically perturbed shell and pointwise '
                   'electroneutrality is not a safe assumption here')
    return {'eps_r': eps_r,
            'theta_bulk_plateau': th,
            'debye_nm_at_bulk': lam,
            'debye_nm_at_first_layer': lam_surf,
            'xi_nm': xi,
            'ratio_debye_over_xi': r,
            'radius_nm': part['radius_nm'],
            'fraction_of_radius': lam / part['radius_nm'],
            'verdict': verdict}


def ceria_comparison(res, part, p, ceria_theta=0.05, eps_r_ceria=25.0):
    """How dilute this is next to the case the construction came from.

    The ceria study ran a 100 nm particle at a volume-averaged vacancy
    site fraction of 0.05 and found local electroneutrality held almost
    everywhere.  Screening length goes as the inverse square root of
    carrier density, so a system this much more dilute screens over a
    correspondingly longer distance, and that is the whole reason the
    conclusion does not transfer for free."""
    th = res['theta_particle']
    return {'theta_here': th, 'theta_ceria': ceria_theta,
            'dilution_factor': ceria_theta / th if th else float('inf'),
            'screening_lengthened_by': math.sqrt(ceria_theta / th) if th
            else float('inf'),
            'note': ('the ceria conclusion was reported as concentration '
                     'dependent by its own authors; this factor is how far '
                     'outside that regime the present case sits')}
