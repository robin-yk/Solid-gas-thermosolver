"""Fifty-digit reference for the particle engine, written independently.

Nothing here imports solidgas.  The dataset is read straight from JSON
and every quantity is recomputed from the algebra, with mpmath at fifty
decimal places and deliberately dumb numerics: one uniform substitution
and plain bisection, where the engine uses two analytic branches.  If a
number agrees between the two it is not because they share a mistake.

    python scripts/oracle_particle.py            regenerate the reference
    python scripts/oracle_particle.py --check    compare without writing
"""

import json
import os
import sys

from mpmath import mp, mpf, exp, log, sqrt, pi, cosh

mp.dps = 50

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data', 'rutile_dft.json')
OUT = os.path.join(ROOT, 'data', 'reference_particle.json')

KB = mpf('8.617333262e-5')
NAVO = mpf('6.02214076e23')


def load():
    with open(DATA) as fh:
        return json.load(fh)


# ------------------------------------------------------------ free energy

def phi(theta, ratio):
    """ln(t/(1-t)) + 2 ln(rt/(1-rt)).  Monotone on (0, 1/r)."""
    p = ratio * theta
    return log(theta / (1 - theta)) + 2 * log(p / (1 - p))


def _phi_y(y, cap):
    """The same quantity in the logistic coordinate theta = cap*sigmoid(y).

    Worth the change of variable: r*theta is exactly sigmoid(y), so the
    polaron logit is exactly y and the whole cation term collapses to
    2y.  What is left,

        phi = ln( cap / (1 + e^-y - cap) ) + 2y,

    has no subtraction of nearly equal numbers anywhere, at either end,
    which the theta form does have - it loses every digit of 1 - r*theta
    once theta is within a rounding of the cap."""
    return log(cap / (1 + exp(-y) - cap)) + 2 * y


def _theta_y(y, cap):
    if y < 0:
        e = exp(y)
        return cap * e / (1 + e)
    return cap / (1 + exp(-y))


def theta_of_phi(x, ratio, span=mpf(1000), steps=200):
    """Invert phi by bisecting the single logistic coordinate.

    theta = cap*sigmoid(y) maps the real line onto (0, cap) monotonically,
    so one bracket covers the dilute and the saturated end alike and
    there is no branch to get wrong.  Relative accuracy in theta equals
    absolute accuracy in y at both ends, and two hundred halvings of a
    two-thousand-wide bracket is well past fifty digits."""
    cap = 1 / mpf(ratio)
    lo, hi = -span, span
    for _ in range(steps):
        mid = (lo + hi) / 2
        if _phi_y(mid, cap) < x:
            lo = mid
        else:
            hi = mid
    return _theta_y((lo + hi) / 2, cap)


def theta_of_mu(mu, e, kt, ratio):
    return theta_of_phi((mu - e) / kt, ratio)


# ------------------------------------------------------------- material

def segregation(p, key):
    q = p['vacancy_energetics']['presets'][key]
    e_s, e_ss, e_b = (mpf(str(q['surface_eV'])), mpf(str(q['subsurface_eV'])),
                      mpf(str(q['bulk_eV'])))
    a = mpf(str(p['lattice_constants_A']['a']))
    d110 = a / sqrt(mpf(2)) / 10
    de = e_b - e_s
    xi = -d110 / log((e_b - e_ss) / de)
    return {'e_s': e_s, 'e_ss': e_ss, 'e_b': e_b, 'de': de,
            'xi': xi, 'd110': d110}


def e_of_z(z, seg):
    return seg['e_b'] - seg['de'] * exp(-z / seg['xi'])


def particle(p, d_um, key):
    a = mpf(str(p['lattice_constants_A']['a']))
    c = mpf(str(p['lattice_constants_A']['c']))
    mm = mpf(str(p['molar_mass_g_mol']))
    rho = 2 * mm / (NAVO * a * a * c * mpf('1e-24'))
    area = 6 / (rho * mpf(str(d_um)) * mpf('1e-4'))          # cm^2/g
    sigma = mpf('1e16') / (c * a * sqrt(mpf(2)))             # bridging /cm^2
    n_o = 2 / mm * mpf('1e6')
    n_s = area * sigma / NAVO * mpf('1e6')
    r_nm = 3 / (rho * area) * mpf('1e7')
    return {'rho': rho, 'area': area, 'N_s': n_s, 'N_i': n_o - n_s,
            'N_o': n_o, 'r_nm': r_nm, 'seg': segregation(p, key)}


# ------------------------------------------------------------ the profile
def layer_stack(part):
    """Every oxygen layer of the interior: its volume and its energy.

    Built once.  Layers deep enough that E(z) equals the bulk value to
    the working precision share one entry - exact arithmetic at fifty
    digits, not a truncation, and what makes summing a fourteen-hundred
    layer particle affordable."""
    seg = part['seg']
    d, r_in = seg['d110'], part['r_nm'] - seg['d110'] / 2
    n = int(r_in / d)
    near, deep_v = [], mpf(0)
    floor = mpf(10) ** (-(mp.dps + 5))
    for k in range(1, n + 1):
        ro = r_in - (k - 1) * d
        ri = r_in - k * d
        if ri < 0:
            ri = mpf(0)
        v = (ro ** 3 - ri ** 3) / 3
        e = e_of_z(k * d, seg)
        if seg['e_b'] - e < floor:
            deep_v += v
        else:
            near.append((v, e))
    core = r_in - n * d
    if core > 0:
        deep_v += core ** 3 / 3
    return {'near': near, 'deep_v': deep_v, 'vol': r_in ** 3 / 3,
            'e_b': seg['e_b'], 'n_layers': n, 'n_distinct': len(near)}


def interior_mean(mu, kt, part, ratio):
    """Volume-weighted mean over the layer stack, summed to the centre."""
    st = part.setdefault('_stack', layer_stack(part))
    total = st['deep_v'] * theta_of_mu(mu, st['e_b'], kt, ratio)
    for v, e in st['near']:
        total += v * theta_of_mu(mu, e, kt, ratio)
    return total / st['vol']

def surface_theta(mu, kt, part, ratio, ph):
    seg = part['seg']
    th_t = mpf(str(ph['theta_transition']))
    th_r = mpf(str(ph['theta_reconstructed_eff']))
    mu_t = seg['e_s'] + kt * phi(th_t, ratio)
    if mu < mu_t:
        return theta_of_mu(mu, seg['e_s'], kt, ratio), mu_t
    if mu > mu_t:
        return th_r, mu_t
    return th_t, mu_t


def held(mu, kt, part, ratio, ph):
    th_s, _ = surface_theta(mu, kt, part, ratio, ph)
    return (part['N_s'] * th_s
            + part['N_i'] * interior_mean(mu, kt, part, ratio))


def equilibrate(vo, t_c, part, ratio, ph, theta_sol):
    kt = KB * (mpf(str(t_c)) + mpf('273.15'))
    seg = part['seg']
    mu_cs = seg['e_b'] + kt * phi(mpf(str(theta_sol)), ratio)

    lo, hi = mpf(-40), mpf(40)
    for _ in range(140):
        mid = (lo + hi) / 2
        if held(mid, kt, part, ratio, ph) < mpf(str(vo)):
            lo = mid
        else:
            hi = mid
    mu = (lo + hi) / 2
    limited = None
    extended = mpf(0)
    if mu > mu_cs:
        mu = mu_cs
        extended = mpf(str(vo)) - held(mu, kt, part, ratio, ph)
        limited = 'crystallographic shear'

    th_s, mu_t = surface_theta(mu, kt, part, ratio, ph)
    th_i = interior_mean(mu, kt, part, ratio)

    # If the bisection landed on the transition potential the surface is
    # mid-conversion and mu cannot say how far: the line compound has no
    # width.  Detected here from the converged mu, and closed by putting
    # whatever the interior did not take onto the row.
    if limited is None and abs(mu - mu_t) < mpf('1e-25'):
        th_s = (mpf(str(vo)) - part['N_i'] * th_i) / part['N_s']
        limited = 'surface reconstruction'
    hs, hi_ = part['N_s'] * th_s, part['N_i'] * th_i
    return {'mu_eV': mu, 'kt_eV': kt,
            'theta_surface': th_s, 'mu_transition_eV': mu_t,
            'theta_interior_mean': th_i,
            'theta_bulk_plateau': theta_of_mu(mu, seg['e_b'], kt, ratio),
            'theta_first_layer': theta_of_mu(
                mu, e_of_z(seg['d110'], seg), kt, ratio),
            'umol_surface': hs, 'umol_interior': hi_,
            'umol_extended': extended, 'mu_shear_eV': mu_cs,
            'limited_by': limited}


# -------------------------------------------------------------- transport

def uptake(tau, terms=400):
    s = mpf(0)
    for n in range(1, terms + 1):
        s += exp(-mpf(n * n) * pi ** 2 * tau) / mpf(n * n)
    return 1 - 6 / pi ** 2 * s


def tau_for_uptake(f, steps=250):
    lo, hi = mpf('1e-14'), mpf(20)
    for _ in range(steps):
        mid = sqrt(lo * hi)
        if uptake(mid) < mpf(str(f)):
            lo = mid
        else:
            hi = mid
    return sqrt(lo * hi)


# --------------------------------------------------------------- assembly

def build():
    p = load()
    ratio = mpf(4)
    ph = p['surface_phases']
    theta_sol = mpf(str(p['saturation']['x_max_shear'])) / 2
    key = p['vacancy_energetics']['default_preset']
    part = particle(p, p['defaults']['particle_diameter_um'], key)

    out = {'dps': mp.dps, 'polaron_ratio': 4,
           'note': ('Independent fifty-digit reference for the particle '
                    'engine. Recomputed from data/rutile_dft.json with a '
                    'single logistic substitution and plain bisection; no '
                    'engine code is imported.')}

    out['phi'] = [{'theta': str(t), 'ratio': int(r),
                   'phi': str(phi(mpf(t), mpf(r)))}
                  for r in (1, 2, 4, 6)
                  for t in ('1e-12', '1e-6', '0.001', '0.05')
                  if mpf(t) < 1 / mpf(r)]
    out['phi_inverse'] = [{'x': str(x), 'ratio': 4,
                           'theta': str(theta_of_phi(mpf(x), ratio))}
                          for x in ('-120', '-60', '-30', '-1.9459101490553132',
                                    '0', '5', '20', '60')]

    # the -1/6 law: theta against ln p(O2) with mu_V = -(1/2) kT ln p
    kt = KB * mpf('873.15')
    rows = []
    for lnp in ('-40', '-30', '-20', '-10'):
        mu = -kt * mpf(lnp) / 2
        rows.append({'ln_pO2': lnp, 'mu_eV': str(mu),
                     'theta': str(theta_of_mu(mu, mpf('1.31'), kt, ratio)),
                     'theta_oxygen_only': str(
                         1 / (1 + exp((mpf('1.31') - mu) / kt)))})
    out['power_law'] = {'kt_eV': str(kt), 'e_form_eV': '1.31', 'rows': rows}

    out['segregation'] = {}
    for k in p['vacancy_energetics']['presets']:
        s = segregation(p, k)
        out['segregation'][k] = {
            'dE_seg_eV': str(s['de']), 'xi_nm': str(s['xi']),
            'd110_nm': str(s['d110']),
            'E_at_layers_eV': [str(e_of_z(n * s['d110'], s))
                               for n in range(4)]}

    out['particle'] = {'diameter_um': p['defaults']['particle_diameter_um'],
                       'radius_nm': str(part['r_nm']),
                       'area_m2_g': str(part['area'] / mpf('1e4')),
                       'N_surface_umol_g': str(part['N_s']),
                       'N_interior_umol_g': str(part['N_i']),
                       'N_oxygen_umol_g': str(part['N_o'])}

    out['equilibrium'] = []
    for t_c, vo in ((600.0, 95.0), (600.0, 5.0), (600.0, 35.0),
                    (600.0, 400.0), (400.0, 95.0), (800.0, 95.0)):
        r = equilibrate(vo, t_c, part, ratio, ph, theta_sol)
        out['equilibrium'].append({
            'T_C': t_c, 'VO_umol_g': vo,
            'mu_eV': str(r['mu_eV']),
            'theta_surface': str(r['theta_surface']),
            'theta_interior_mean': str(r['theta_interior_mean']),
            'theta_bulk_plateau': str(r['theta_bulk_plateau']),
            'theta_first_layer': str(r['theta_first_layer']),
            'umol_surface': str(r['umol_surface']),
            'umol_interior': str(r['umol_interior']),
            'umol_extended': str(r['umol_extended']),
            'mu_transition_eV': str(r['mu_transition_eV']),
            'mu_shear_eV': str(r['mu_shear_eV']),
            'limited_by': r['limited_by']})
        sys.stderr.write('  equilibrium %g C %g umol -> mu %s\n'
                         % (t_c, vo, mp.nstr(r['mu_eV'], 12)))

    out['crank'] = {'tau_for_uptake': {f: str(tau_for_uptake(f))
                                       for f in ('0.2', '0.5', '0.9', '0.99')},
                    'uptake_at_tau': {t: str(uptake(mpf(t)))
                                      for t in ('0.001', '0.01', '0.1', '0.5')}}
    return out


def main():
    check = '--check' in sys.argv
    out = build()
    if check:
        with open(OUT) as fh:
            old = json.load(fh)
        print('identical' if old == out else 'DIFFERS')
        return 0 if old == out else 1
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
        fh.write('\n')
    print('wrote %s' % OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
