"""Build one equilibrium model per (energy map, closure, particle diameter).

Energy maps
  PAB, HAM   Pabisiak 2007 GGA (layers 1-4) and Hameeuw 2006 LDA (layers 1-3),
             per site class BRI / IPL / SBR, relative to bulk.
  LI_SBR1    Li 2015 sX surface -1.31 on layer-1 BRI, subsurface -0.72 on
             layer-1 SBR. Other explicit sites 0 (assumed).
  LI_L2      Li 2015 sX surface -1.31 on layer-1 BRI, subsurface -0.72 on all
             four layer-2 O sites. Other explicit sites 0 (assumed).
  LI_CONT    Manuscript Note 2b: G(z) = -A exp(-z/xi), A = 1.31 eV, continuum
             O sites; bridging coverage taken as x(z = 0). xi = 0.25, 0.5, 1 nm.
Closures
  NEUTRAL    electrons implicit; one constraint (inventory).
  LOCAL      every trilayer (or continuum shell) neutral: 2 N_V = N_Ti3+, with
             ideal Ti3+ configurational entropy (manuscript Note 2b form).
"""
import csv
import math
from pathlib import Path

import numpy as np

from . import particle as pt
from .engine import Model

HERE = Path(__file__).resolve().parent.parent
A_LI, CONT_XI = 1.31, (0.25, 0.5, 1.0)
DIAMETERS = (900.0, 600.0, 300.0)
CLOSURES = ('NEUTRAL', 'LOCAL')


def energy_sets():
    sets = {}
    with open(HERE / 'specification' / 'energy_sets.csv') as f:
        for r in csv.DictReader(f):
            sets.setdefault(r['scenario'], {})[(r['layer'], r['site'])] = float(r['relative_E_eV'])
    return sets


def discrete_maps():
    s = energy_sets()
    pab = {(int(k), site): e for (k, site), e in s['PABISIAK_2007_GGA'].items()}
    ham = {(int(k), site): e for (k, site), e in s['HAMEEUW_2006_LDA'].items()}
    li = s['LI_2015_SX']
    surf, sub = li[('SURFACE', 'GENERIC')], li[('SUBSURFACE', 'GENERIC')]
    return {
        'PAB': (4, pab),
        'HAM': (3, ham),
        'LI_SBR1': (1, {(1, 'BRI'): surf, (1, 'SBR'): sub}),
        'LI_L2': (2, {(1, 'BRI'): surf, **{(2, x): sub for x in ('BRI', 'IPL', 'SBR')}}),
    }


def _two_state(C, E, tag, region, electron=False):
    d = dict(C=C, E=[0.0, E], tag=tag, region=region)
    d.update(v=[0, 0], e=[0, 1]) if electron else d.update(v=[0, 1])
    return d


def build_discrete(n_layers, emap, closure, d_nm):
    local = closure == 'LOCAL'
    o, ti = pt.layer_sites(d_nm, n_layers)
    doms = []
    for k, site, C, _z in o:
        doms.append(_two_state(C, emap.get((k, site), 0.0), (k, site), k if local else -1))
    bulk_o = pt.O_TOTAL - sum(c for _, _, c, _ in o)
    doms.append(_two_state(bulk_o, 0.0, 'bulk', 0 if local else -1))
    if local:
        for k, C in enumerate(ti, 1):
            doms.append(_two_state(C, 0.0, 'Ti', k, electron=True))
        doms.append(_two_state(pt.TI_TOTAL - sum(ti), 0.0, 'Ti', 0, electron=True))
    return Model(doms)


def continuum_edges(R, n=3000, z0=1e-4):
    return np.concatenate([[0.0], np.geomspace(z0, R, n)])


def build_continuum(xi, closure, d_nm, n=3000):
    """Note 2b as a finite-volume model: each spherical shell is one O domain
    (and, for LOCAL, one Ti domain with half its sites) and one region."""
    local = closure == 'LOCAL'
    R = d_nm / 2
    z = continuum_edges(R, n)
    doms = []
    for j in range(n):
        C = pt.O_TOTAL * pt.shell_fraction(R, z[j], z[j + 1])
        E = -A_LI * math.exp(-0.5 * (z[j] + z[j + 1]) / xi)
        doms.append(_two_state(C, E, 'shell', j if local else -1))
        if local:
            doms.append(_two_state(C / 2, 0.0, 'Ti', j, electron=True))
    return Model(doms)


def continuum_surface_x(mu_eV, T, local):
    """Site fraction at z = 0 from the stationarity condition
    -A + kT[ln x/(1-x) + 2 ln 4x/(1-4x)] = mu (LOCAL) or the Fermi form."""
    kT = pt.KB * T
    if not local:
        return 1.0 / (1.0 + math.exp((-A_LI - mu_eV) / kT))
    lo, hi = 0.0, 0.25
    for _ in range(200):
        x = 0.5 * (lo + hi)
        f = -A_LI + kT * (math.log(x / (1 - x)) + 2 * math.log(4 * x / (1 - 4 * x))) - mu_eV
        lo, hi = (lo, x) if f > 0 else (x, hi)
    return 0.5 * (lo + hi)


def scenarios():
    """Yield (energy_map, closure, diameter, builder, theta_bri(sol, T))."""
    for d_nm in DIAMETERS:
        for closure in CLOSURES:
            for name, (nl, emap) in discrete_maps().items():
                def theta(sol, T, d_nm=d_nm):
                    C = pt.layer_sites(d_nm, 1)[0][0][2]
                    return sol.vacancies((1, 'BRI')) / C
                yield name, closure, d_nm, (lambda nl=nl, e=emap, c=closure, d=d_nm: build_discrete(nl, e, c, d)), theta
            for xi in CONT_XI:
                def theta(sol, T, c=closure):
                    return continuum_surface_x(sol.mu_eV, T, c == 'LOCAL')
                yield f'LI_CONT_xi{xi:g}', closure, d_nm, (lambda xi=xi, c=closure, d=d_nm: build_continuum(xi, c, d)), theta
