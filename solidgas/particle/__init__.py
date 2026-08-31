"""Vacancy distribution in one oxide particle, from a titrated inventory.

The engine is built on the equilibrium construction used for reduced
ceria nanoparticles - a lattice gas on the oxygen sublattice whose
electrons are paid for on the cation sublattice, distributed over a
continuous radial formation-energy profile and pinned by conservation
of the measured inventory - and the rutile physics is put on top of it.

    from solidgas import particle
    p    = particle.load()
    part = particle.particle(p, d_um=0.9)
    res  = particle.equilibrate(95.0, 600.0, part)

What rutile adds to that base: a distinguished bridging-oxygen row with
a quarter of a full layer's site density, a (110) surface that
reconstructs before the cation sublattice fills, and a bulk that throws
crystallographic shear planes at x = 0.008 rather than dissolving point
defects all the way to Ti2O3.
"""

from .freeenergy import (POLARON_RATIO, THETA_CAP, cap, phi,
                         phi_oxygen_only, mu_of_theta, theta_of_mu,
                         theta_of_phi, dilute_theta, free_energy_per_site)
from .material import (KB_EV, load, segregation, formation_energy,
                       particle, density_g_cm3, bridging_density_cm2,
                       theta_from_loading, loading_from_theta, kt_eV)
from .profile import (equilibrate, inventory, interior_average,
                      interior_average_smeared, surface_state,
                      depth_profile, shell_split)

__all__ = [
    'POLARON_RATIO', 'THETA_CAP', 'cap', 'phi', 'phi_oxygen_only',
    'mu_of_theta', 'theta_of_mu', 'theta_of_phi', 'dilute_theta',
    'free_energy_per_site', 'KB_EV', 'load', 'segregation',
    'formation_energy', 'particle', 'density_g_cm3',
    'bridging_density_cm2', 'theta_from_loading', 'loading_from_theta',
    'kt_eV', 'equilibrate', 'inventory', 'interior_average',
    'interior_average_smeared', 'surface_state', 'depth_profile',
    'shell_split',
]
