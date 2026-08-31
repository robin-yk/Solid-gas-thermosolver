"""Can the particle reach the profile the equilibrium model draws?

An equilibrium profile is only the answer if the oxygen had time to
arrange itself.  The ceria study this engine follows makes that check
its first two models: a Crank sphere-diffusion solution for the
reduction step, and then a closed-system relaxation showing how long the
resulting front takes to flatten once the gas is switched off.  The same
two questions decide whether the rutile profile here is a prediction or
a decoration.

The transport numbers:

    D = a^2 nu0 exp(-E_m / kT)

with the oxygen jump distance a, an attempt frequency nu0, and the
migration barrier E_m.  For a sphere the relevant clock is the Fourier
number tau = D t / R^2; the profile is flat to a percent by tau of
roughly one half, and the classical Crank series gives the fractional
uptake at any tau.

The conclusion for rutile at the shipped conditions is short: with the
literature bulk barrier, homogenisation takes under a millisecond and
transport is not the limit.  That is a strong statement and it comes
with a margin, which is why the module also inverts the question and
reports what barrier WOULD make transport matter.  The group's own
recovery measurement implies a higher effective barrier than the
literature one - charged migration, trapping, aggregation - so the
margin is the number worth quoting, not the time.
"""

import math

from .material import KB_EV

# Crank's series for the fractional uptake of a sphere with a fixed
# surface concentration, M(t)/M(inf) = 1 - (6/pi^2) sum n^-2 exp(-n^2 pi^2 tau).
_CRANK_TERMS = 200


def diffusivity(t_c, e_m_eV, nu0_per_s, a_nm):
    """Tracer diffusivity in nm^2/s from a jump distance and a barrier."""
    kt = KB_EV * (t_c + 273.15)
    return a_nm * a_nm * nu0_per_s * math.exp(-e_m_eV / kt)


def uptake(tau):
    """Fractional approach to equilibrium of a sphere at Fourier number tau."""
    if tau <= 0.0:
        return 0.0
    s = 0.0
    for n in range(1, _CRANK_TERMS + 1):
        x = n * n * math.pi * math.pi * tau
        if x > 700.0:
            break
        s += math.exp(-x) / (n * n)
    return max(0.0, min(1.0, 1.0 - 6.0 / (math.pi ** 2) * s))


def tau_for_uptake(f):
    """Fourier number at which the sphere is a fraction f of the way there."""
    if not 0.0 < f < 1.0:
        raise ValueError('uptake fraction must be strictly between 0 and 1')
    lo, hi = 1e-12, 10.0
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if uptake(mid) < f:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def timescales(t_c, part, e_m_eV=None, nu0_per_s=None, a_nm=None,
               p=None, fractions=(0.5, 0.9, 0.99)):
    """How long the particle needs, against how long the experiment gives it.

    a defaults to the rutile c axis, which is the jump the literature
    barrier belongs to."""
    if p is not None:
        kin = p.get('kinetics_cgmc', {})
        e_m_eV = e_m_eV if e_m_eV is not None else kin.get('E_m_bulk_eV')
        nu0_per_s = nu0_per_s if nu0_per_s is not None else kin.get('nu0_per_s')
        a_nm = a_nm if a_nm is not None else (
            p['lattice_constants_A']['c'] / 10.0)
    if None in (e_m_eV, nu0_per_s, a_nm):
        raise ValueError('need E_m, nu0 and a, or the dataset to read them from')

    r_nm = part['radius_nm']
    d = diffusivity(t_c, e_m_eV, nu0_per_s, a_nm)
    scale = r_nm * r_nm / d
    return {'T_C': t_c, 'E_m_eV': e_m_eV, 'nu0_per_s': nu0_per_s,
            'a_nm': a_nm, 'radius_nm': r_nm,
            'D_nm2_per_s': d, 'D_cm2_per_s': d * 1e-14,
            'tau_unit_s': scale,
            'times_s': {f: tau_for_uptake(f) * scale for f in fractions}}


def barrier_for_time(t_c, part, t_target_s, nu0_per_s, a_nm,
                     fraction=0.99):
    """The migration barrier that would make homogenisation take t_target_s.

    Inverting the clock is the honest way to state the margin: rather
    than claiming transport is fast, say how much slower it would have
    to be before it stopped being fast, and let the reader judge whether
    the gap between that and the literature barrier is comfortable."""
    kt = KB_EV * (t_c + 273.15)
    r_nm = part['radius_nm']
    d_needed = tau_for_uptake(fraction) * r_nm * r_nm / t_target_s
    return kt * math.log(a_nm * a_nm * nu0_per_s / d_needed)


def verdict(t_c, part, exposure_s, p=None, fraction=0.99, **kw):
    """Is the equilibrium profile reachable inside the experiment's window?

    Returns the two times, their ratio, and the barrier that would close
    the gap.  A ratio far below one is what licenses using an
    equilibrium profile at all; the barrier is what says how far from
    that licence the system is."""
    ts = timescales(t_c, part, p=p, fractions=(fraction,), **kw)
    t_hom = ts['times_s'][fraction]
    e_break = barrier_for_time(t_c, part, exposure_s, ts['nu0_per_s'],
                               ts['a_nm'], fraction)
    ts.update({
        'fraction': fraction,
        'exposure_s': exposure_s,
        't_homogenise_s': t_hom,
        'ratio': t_hom / exposure_s if exposure_s else float('inf'),
        'equilibrium_reachable': t_hom < exposure_s,
        'E_m_that_would_break_it_eV': e_break,
        'barrier_margin_eV': e_break - ts['E_m_eV'],
    })
    return ts
