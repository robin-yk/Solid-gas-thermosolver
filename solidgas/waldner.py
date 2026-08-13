"""Ti-O condensed phases from the Waldner & Eriksson CALPHAD assessment.

    P. Waldner and G. Eriksson, "Thermodynamic modelling of the system
    titanium-oxygen", Calphad 23 (1999) 189-218. Appendix, pp. 216-218.

Why a second source at all: NIST-JANAF has no Magneli member between rutile
and Ti4O7. Those are exactly the phases that decide whether the reverse
water-gas-shift conclusion survives at the hot end, because they form more
easily than Ti4O7 does. This assessment carries the whole series, Ti4O7
through Ti10O19 plus Ti20O39.

The two sources must not be mixed for the condensed phases. Waldner puts
rutile at dHf298 = -944750 J/mol against NIST's -938722, a 6.0 kJ/mol offset
that is comparable to the margins being measured. Everything titanium-bearing
comes from here; the gases stay on NIST.

Gibbs energy is built the CALPHAD way,

    G(T) = dHf298 + int(Cp dT) - T * [S298 + int(Cp/T dT)]

integrated piecewise from 298.15 K. Where a phase changes polymorph inside its
range the paper prints a second enthalpy on that range's row - 11757 J for
Ti3O5 at 450 K, 1138 J for Ti2O3 at 470 K. Those are transition enthalpies,
not formation enthalpies, and they must be added at the join along with
dS = dH/T_trans. Enthalpy and entropy both step there; Gibbs energy does not,
which is the check. Dropping them looks harmless because G stays continuous
either way, but it is not: the error grows as dH*(1 - T/T_trans), which for
Ti3O5 reaches -35 kJ/mol by 1500 C and puts the phase in the wrong place
entirely.
"""

import numpy as np

# Cp forms used in the appendix.
#   'poly'    Cp = a + b*T + c*T^2 + d*T^-2       (rutile vacancy end-member,
#                                                  Ti2O3, Ti3O5, TiO, Ti3O2)
#   'magneli' Cp = k0 + k1*T^-0.5 + k2*T^-2 + k3*T^-3   (rutile, Magneli series)

# phase -> dict(dHf298, S298, ranges=[(T_lo, T_hi, form, coeffs, dH_trans), ...])
# dH_trans is the enthalpy added on entering that range, in J/mol; zero for the
# first range and for any join that is only a change of fit.
WALDNER = {
    # (Ti4+)1(O2-)2 / RUTILE end-member, p.217. k0 here is the same 77.83762
    # that each successive Magneli member adds, which is the series' internal
    # check: every step up in n adds exactly one TiO2 unit.
    'TiO2': dict(dHf298=-944750.0, S298=50.46, ranges=[
        (298.15, 6000.0, 'magneli', (77.83762, 0.0, -3367841.0, 402940672.0), 0.0)]),

    'Ti2O3': dict(dHf298=-1526171.0, S298=77.463, ranges=[
        (298.15, 470.0, 'poly', (30.3934, 0.199918, 0.355968e-04, 235600.0), 0.0),
        (470.0, 2115.0, 'poly', (147.67400, 0.34742300e-02, 0.92030000e-09, -4790850.0), 1138.0)]),

    'Ti3O5': dict(dHf298=-2469900.0, S298=127.75, ranges=[
        (298.15, 450.0, 'poly', (-23.9073, 0.84031, -0.8084413e-03, 0.0), 0.0),
        (450.0, 5000.0, 'poly', (158.99208, 0.50207950e-01, 0.0, 0.0), 11757.0)]),

    # Magneli series Ti_n O_(2n-1), p.218. k1 is the same for every member --
    # the crystallographic shear contribution does not depend on n.
    'Ti4O7': dict(dHf298=-3416525.0, S298=195.39, ranges=[
        (298.15, 5000.0, 'magneli', (364.36711, -2110.8923, -2519517.0, -0.1445482e+09), 0.0)]),
    'Ti5O9': dict(dHf298=-4360918.5, S298=248.95446, ranges=[
        (298.15, 5000.0, 'magneli', (442.20473, -2110.8923, -5887358.0, 0.2583924e+09), 0.0)]),
    'Ti6O11': dict(dHf298=-5306479.2, S298=300.50153, ranges=[
        (298.15, 5000.0, 'magneli', (520.04235, -2110.8923, -9255199.0, 0.6613331e+09), 0.0)]),
    # dHf298 is misprinted in the paper as -3415327.0, which is Ti4O7's value
    # to four figures and breaks a series that is otherwise linear in n to
    # under 600 J. See TI7O13_PRINTED and the test that pins this.
    'Ti7O13': dict(dHf298=None, S298=351.91071, ranges=[
        (298.15, 5000.0, 'magneli', (597.87997, -2110.8923, -12623040.0, 0.10642738e+10), 0.0)]),
    'Ti8O15': dict(dHf298=-7196482.4, S298=403.17028, ranges=[
        (298.15, 5000.0, 'magneli', (675.71759, -2110.8923, -15990881.0, 0.14672144e+10), 0.0)]),
    'Ti9O17': dict(dHf298=-8141241.7, S298=454.36515, ranges=[
        (298.15, 5000.0, 'magneli', (753.55521, -2110.8923, -19358722.0, 0.18701551e+10), 0.0)]),
    'Ti10O19': dict(dHf298=-9085912.9, S298=505.51533, ranges=[
        (298.15, 5000.0, 'magneli', (831.39283, -2110.8923, -22726564.0, 0.22730958e+10), 0.0)]),
    'Ti20O39': dict(dHf298=-18533039.0, S298=1015.3802, ranges=[
        (298.15, 5000.0, 'magneli', (1609.7690, -2110.8923, -56404974.0, 0.63025026e+10), 0.0)]),
}

TI7O13_PRINTED = -3415327.0


def _fit_ti7o13():
    """Recover Ti7O13's formation enthalpy from the rest of the series.

    dHf298 is linear in n across the Magneli members to better than 600 J, so
    the misprint is repaired by least squares on the neighbours rather than by
    guessing.
    """
    ns, hs = [], []
    for name, d in WALDNER.items():
        if name.startswith('Ti') and d['dHf298'] is not None and 'O' in name:
            n_ti = _stoich(name)[0]
            if 4 <= n_ti <= 10 and name != 'TiO2':
                ns.append(n_ti)
                hs.append(d['dHf298'])
    slope, intercept = np.polyfit(ns, hs, 1)
    return float(slope * 7 + intercept), float(slope), float(intercept)


def _stoich(name):
    """(n_Ti, n_O) from a formula string like 'Ti10O19'."""
    import re
    m = re.fullmatch(r'Ti(\d*)O(\d*)', name)
    if not m:
        raise ValueError(name)
    return (int(m.group(1) or 1), int(m.group(2) or 1))


def _cp(form, c, T):
    if form == 'poly':
        a, b, cc, d = c
        return a + b * T + cc * T**2 + d / T**2
    k0, k1, k2, k3 = c
    return k0 + k1 * T**-0.5 + k2 / T**2 + k3 / T**3


def _int_cp(form, c, T0, T1):
    """int Cp dT from T0 to T1."""
    if form == 'poly':
        a, b, cc, d = c
        F = lambda T: a * T + b * T**2 / 2 + cc * T**3 / 3 - d / T
    else:
        k0, k1, k2, k3 = c
        F = lambda T: k0 * T + 2 * k1 * np.sqrt(T) - k2 / T - k3 / (2 * T**2)
    return F(T1) - F(T0)


def _int_cp_over_t(form, c, T0, T1):
    """int Cp/T dT from T0 to T1."""
    if form == 'poly':
        a, b, cc, d = c
        F = lambda T: a * np.log(T) + b * T + cc * T**2 / 2 - d / (2 * T**2)
    else:
        k0, k1, k2, k3 = c
        F = lambda T: k0 * np.log(T) - 2 * k1 / np.sqrt(T) - k2 / (2 * T**2) - k3 / (3 * T**3)
    return F(T1) - F(T0)


def _segments(name, T):
    """Range slices from 298.15 K up to T, with the transition entered at each."""
    out, lo = [], 298.15
    for (r0, r1, form, c, dHt) in WALDNER[name]['ranges']:
        if T <= r0:
            break
        hi = min(T, r1)
        if hi > max(lo, r0):
            out.append((max(lo, r0), hi, form, c, dHt, r0))
            lo = hi
        if T <= r1:
            break
    if out and out[-1][1] < T:      # extrapolate the last range if asked
        r0, r1, form, c, dHt = WALDNER[name]['ranges'][-1]
        out.append((out[-1][1], T, form, c, 0.0, r0))
    return out


def dHf(name):
    d = WALDNER[name]
    return _fit_ti7o13()[0] if d['dHf298'] is None else d['dHf298']


def heat_capacity(name, T):
    """Cp in J/mol/K."""
    for (r0, r1, form, c, _dHt) in WALDNER[name]['ranges']:
        if T <= r1:
            return _cp(form, c, T)
    r0, r1, form, c, _dHt = WALDNER[name]['ranges'][-1]
    return _cp(form, c, T)


def enthalpy(name, T):
    """H in J/mol, referenced to the elements. Steps at each transition."""
    h = dHf(name)
    for (t0, t1, form, c, dHt, T_tr) in _segments(name, T):
        h += dHt + _int_cp(form, c, t0, t1)
    return h


def entropy(name, T):
    """S in J/mol/K. Steps by dH/T_trans at each transition."""
    s = WALDNER[name]['S298']
    for (t0, t1, form, c, dHt, T_tr) in _segments(name, T):
        s += (dHt / T_tr if dHt else 0.0) + _int_cp_over_t(form, c, t0, t1)
    return s


def mu0(name, T):
    """G = H - T*S in kJ/mol, to match the rest of the package."""
    return (enthalpy(name, T) - T * entropy(name, T)) / 1000.0
