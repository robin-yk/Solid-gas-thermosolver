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
    # rutile TiO2 - NIST-JANAF 4th ed (298-2000 K)
    'TiO2': ([67.29830, 18.70940, -11.57900, 2.449561, -1.485471, -964.5140, 117.8630, -938.7220],
             [67.29830, 18.70940, -11.57900, 2.449561, -1.485471, -964.5140, 117.8630, -938.7220], -938.722),
    # Ti4O7 (Magneli) - NIST-JANAF 4th ed (298-1950 K)
    'Ti4O7': ([248.9442, 81.06458, -36.39967, 6.641765, -5.509993, -3500.652, 446.0018, -3404.525],
              [248.9442, 81.06458, -36.39967, 6.641765, -5.509993, -3500.652, 446.0018, -3404.525], -3404.525),
}

# Temperature window each entry is fitted over, in K. Ti4O7 is the binding one:
# it stops at 1950 K, i.e. 1677 C.
VALID_RANGE = {
    'CH4': (298.0, 6000.0), 'CO2': (298.0, 6000.0), 'CO': (298.0, 6000.0),
    'H2': (298.0, 2500.0), 'H2O': (500.0, 1700.0),
    'TiO2': (298.0, 2000.0), 'Ti4O7': (298.0, 1950.0),
}


def enthalpy(name, T):
    """Standard molar enthalpy in kJ/mol, referenced to the elements."""
    lo, hi, dHf = SHOMATE[name]
    A, B, C, D, E, F, _G, H = lo if T <= 1000 else hi
    t = T / 1000.0
    return dHf + (A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F - H)


def heat_capacity(name, T):
    """Standard molar heat capacity in J/mol/K."""
    lo, hi, _dHf = SHOMATE[name]
    A, B, C, D, E, _F, _G, _H = lo if T <= 1000 else hi
    t = T / 1000.0
    return A + B * t + C * t**2 + D * t**3 + E / t**2


def entropy(name, T):
    """Standard molar entropy in J/mol/K."""
    lo, hi, _dHf = SHOMATE[name]
    A, B, C, D, E, _F, G, _H = lo if T <= 1000 else hi
    t = T / 1000.0
    return A * np.log(t) + B * t + C * t**2 / 2 + D * t**3 / 3 - E / (2 * t**2) + G


def mu0(name, T):
    """Standard chemical potential G = H - T*S in kJ/mol."""
    lo, hi, dHf = SHOMATE[name]
    c = lo if T <= 1000 else hi
    A, B, C, D, E, F, G, H = c
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
