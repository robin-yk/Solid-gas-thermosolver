"""Shomate polynomials for the gas and solid species used in this repo.

Every coefficient set is NIST-JANAF (4th ed.). Each entry is

    (low-T coefficients, high-T coefficients, dHf298 in kJ/mol)

and `mu0` switches between the two sets at 1000 K. The condensed phases carry
one set repeated twice, so the switch is a no-op for them.
"""

import numpy as np

# Gas constant in kJ/mol/K, and in cm3*atm/mol/K for the ideal-gas charge.
R_KJ = 8.314472e-3
R_ATM = 82.06

SHOMATE = {
    'CH4':  ([-0.703029, 108.4773, -42.52157, 5.862788, 0.678565, -76.84376, 158.7163, -74.87310],
             [85.81217, 11.26467, -2.114146, 0.138190, -26.42221, -153.5327, 224.4143, -74.87310], -74.873),
    'CO2':  ([24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431, -393.5224],
             [58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, 263.6125, -393.5224], -393.510),
    'CO':   ([25.56759, 6.096130, 4.054656, -2.671301, 0.131021, -118.0089, 227.3665, -110.5271],
             [35.15070, 1.300095, -0.205921, 0.013550, -3.282780, -127.8375, 231.7120, -110.5271], -110.527),
    'H2':   ([33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.7080, 0.0],
             [18.563083, 12.257357, -2.859786, 0.268238, 1.977990, -1.147438, 156.2881, 0.0], 0.0),
    'H2O':  ([30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, 223.3967, -241.8264],
             [30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, 223.3967, -241.8264], -241.826),
    # O2 splits at 700 K, not 1000 - see BREAKPOINT.
    'O2':   ([31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0.0],
             [30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0], 0.0),
    'N2':   ([28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 226.4168, 0.0],
             [19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 212.3900, 0.0], 0.0),
    # rutile TiO2 - NIST-JANAF 4th ed (298-2000 K)
    'TiO2': ([67.29830, 18.70940, -11.57900, 2.449561, -1.485471, -964.5140, 117.8630, -938.7220],
             [67.29830, 18.70940, -11.57900, 2.449561, -1.485471, -964.5140, 117.8630, -938.7220], -938.722),
    # Ti4O7 (Magneli) - NIST-JANAF 4th ed (298-1950 K)
    'Ti4O7': ([248.9442, 81.06458, -36.39967, 6.641765, -5.509993, -3500.652, 446.0018, -3404.525],
              [248.9442, 81.06458, -36.39967, 6.641765, -5.509993, -3500.652, 446.0018, -3404.525], -3404.525),
    # Ti3O5 beta - NIST-JANAF 4th ed (298-2050 K). This is the polymorph that
    # matters here: beta overtakes alpha at 450 K, far below the 773 K floor.
    'Ti3O5': ([158.9933, 50.20633, 0.001071, -0.000225, -0.000037, -2495.827, 335.0547, -2446.192],
              [158.9933, 50.20633, 0.001071, -0.000225, -0.000037, -2495.827, 335.0547, -2446.192], -2446.192),
    # Ti2O3 solid - NIST-JANAF 4th ed. The 470-2115 K fit is used across the
    # board: the whole working range sits inside it, and G is continuous with
    # the 298-470 K fit at the join (Cp, S and H each jump there, across the
    # metal-insulator transition, but the free energy does not).
    # Note the liquid entry in the same NIST page is a different phase
    # (dHf = -1418.46) and must not be used here.
    'Ti2O3': ([147.5509, 3.649536, -0.081337, 0.011812, -4.774029, -1579.849, 230.2715, -1520.884],
              [147.5509, 3.649536, -0.081337, 0.011812, -4.774029, -1579.849, 230.2715, -1520.884], -1520.884),
    # Ti3O5 alpha (298-1500 K), kept for reference and for the polymorph test.
    # Not a member of SPECIES - it is never the stable form in this range.
    'Ti3O5_alpha': ([233.4149, -59.35757, 68.20422, -20.59934, -5.901490, -2546.516, 393.3579, -2459.150],
                    [233.4149, -59.35757, 68.20422, -20.59934, -5.901490, -2546.516, 393.3579, -2459.150], -2459.150),
    # Ceria, carried so the same engine can run a second oxide system. Both
    # entries have one coefficient set standing in for the whole range rather
    # than the usual low/high pair, which is how they arrive from the notebook
    # this work started in; the page says so where it shows them.
    'CeO2': ([69.3600, 9.4400, -7.4800, 2.1400, -1.0940, -1133.920, 133.900, -1088.700],
             [69.3600, 9.4400, -7.4800, 2.1400, -1.0940, -1133.920, 133.900, -1088.700],
             -1088.700),
    'Ce2O3': ([117.400, 17.800, -12.200, 3.200, -2.200, -1869.200, 246.600, -1796.200],
              [117.400, 17.800, -12.200, 3.200, -2.200, -1869.200, 246.600, -1796.200],
              -1796.200),
}


# Where each entry switches from its low- to its high-temperature coefficients.
# NIST does not use one breakpoint for everything; 1000 K is the common case.
BREAKPOINT = {'O2': 700.0, 'N2': 500.0}
DEFAULT_BREAKPOINT = 1000.0


def _coeffs(name, T):
    lo, hi, _dHf = SHOMATE[name]
    return lo if T <= BREAKPOINT.get(name, DEFAULT_BREAKPOINT) else hi


# Temperature window each entry is fitted over, in K. Ti4O7 is the binding one:
# it stops at 1950 K, i.e. 1677 C.
VALID_RANGE = {
    'CH4': (298.0, 6000.0), 'CO2': (298.0, 6000.0), 'CO': (298.0, 6000.0),
    'H2': (298.0, 2500.0), 'H2O': (500.0, 1700.0),
    'O2': (100.0, 2000.0), 'N2': (100.0, 2000.0),
    'TiO2': (298.0, 2000.0), 'Ti4O7': (298.0, 1950.0),
    'Ti3O5': (298.0, 2050.0), 'Ti3O5_alpha': (298.0, 1500.0),
    'Ti2O3': (470.0, 2115.0),
    # No fit range was published with the ceria coefficients - they arrive from
    # the notebook as a single set with no range attached. This records the range
    # the tools actually evaluate them over, which is not the same as a range the
    # fit is known to be good on. See test_ceria_fit_is_not_self_consistent.
    'CeO2': (298.15, 2000.0),
    'Ce2O3': (298.15, 2000.0),
}


def enthalpy(name, T):
    """Standard molar enthalpy in kJ/mol, referenced to the elements."""
    _lo, _hi, dHf = SHOMATE[name]
    A, B, C, D, E, F, _G, H = _coeffs(name, T)
    t = T / 1000.0
    return dHf + (A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F - H)


def heat_capacity(name, T):
    """Standard molar heat capacity in J/mol/K."""
    A, B, C, D, E, _F, _G, _H = _coeffs(name, T)
    t = T / 1000.0
    return A + B * t + C * t**2 + D * t**3 + E / t**2


def entropy(name, T):
    """Standard molar entropy in J/mol/K."""
    A, B, C, D, E, _F, G, _H = _coeffs(name, T)
    t = T / 1000.0
    return A * np.log(t) + B * t + C * t**2 / 2 + D * t**3 / 3 - E / (2 * t**2) + G


def mu0(name, T):
    """Standard chemical potential G = H - T*S in kJ/mol."""
    _lo, _hi, dHf = SHOMATE[name]
    A, B, C, D, E, F, G, H = _coeffs(name, T)
    t = T / 1000.0
    Ht = dHf + (A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F - H)
    St = A * np.log(t) + B * t + C * t**2 / 2 + D * t**3 / 3 - E / (2 * t**2) + G
    return Ht - T * St / 1000.0


def Kp(reaction, T):
    """Equilibrium constant from standard potentials.

    `reaction` maps species name to stoichiometric coefficient, negative for
    reactants. RWGS is {'CO2': -1, 'H2': -1, 'CO': 1, 'H2O': 1}.
    """
    dG = sum(nu * mu0(name, T) for name, nu in reaction.items())
    return np.exp(-dG / (R_KJ * T))
