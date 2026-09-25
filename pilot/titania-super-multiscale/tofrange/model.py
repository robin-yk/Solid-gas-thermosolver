"""The model: one rutile powder, every element on at once.

A sample is an equal-mass mixture of spheres at 8 diameters over 900-1600 nm
(same inventory per gram in each). Each sphere is one state table:

  surface   four explicit (110) trilayers on a share f110 of the surface, with
            bridging (BRI), in-plane (IPL, basal) and sub-bridging (SBR) O and
            two Ti per 1x1 cell, at their crystallographic depths. Energies
            from one map; sites a map does not name are at 0 (bulk). The Li
            maps give no basal value: basal = bridging + 0.11 eV (Matsunaga).
            The rest of the surface is bulk-like O with no reactive sites.
  (1x2)     every pair of bridging sites and the four layer-1 Ti under them
            form one cell: two bridging sites (00, one, both vacant) with
            0-4 Ti3+ on the four Ti, or the reconstructed Ti2O3 row (one
            vacancy equivalent, energy eps_BRI + dG) that keeps two Ti3+ on
            two of its Ti, with 0-2 more on the other two.
  bulk      below the trilayers (1.30 nm), O sits in periodic boxes of 700
            sites. A box with N vacancies has free energy -kT ln Q(N) from
            Monte Carlo with the ZHA2017 pair energy summed over every pair
            within a cutoff; clusters of any size that fit are included.
  charge    every vacancy leaves two electrons on Ti (Ti3+, one per Ti).
            Only the particle is neutral; every atomic plane is a charged
            shell in a medium of the rutile permittivity (Parker 1961, 873 K).
            Layer-1 Ti carry the RET2018 penalty, 0.2 eV.
  history   600 C equilibrium: the treatment state relaxes to it within
            seconds (checked with hop kinetics in an earlier revision).

Inputs without a source are the model's parameters: energy map, aggregate
cutoff, reconstruction dG, f110, permittivity axis. The measured inventory is
fixed; the CO rate enters only at the end, TOF = rate / reactive sites.
"""
import csv
import math
from pathlib import Path

import numpy as np

from . import particle as pt
from .engine import Model

HERE = Path(__file__).resolve().parent.parent
T_EQ = 873.15
K_EXPL = 4
N_BULK, Z0_BULK = 240, 0.02
S0_PENALTY = 0.2
MATSUNAGA = 0.11
ZHA = (0.753, 11.83, -383.45)
PLANES = {'BRI': 0, 'IPL': 1, 'SBR': 2}

# Parameters (no source fixes them): the model is solved over all of them.
MAPS = ('PAB', 'HAM', 'LI_SBR1', 'LI_L2')
CUTOFFS = (0.28, 0.34, 0.40, 0.45, 0.6, 1.0)
DG = (-0.4, -0.2, 0.0, 0.2, 0.4)
F110 = (1.0, 0.75, 0.5)
EPS = {'a_axis': 64.0, 'c_axis': 107.0}
DIAMETERS = tuple(np.linspace(900.0, 1600.0, 4))   # equal mass each; 8 would change the size factor by 1 %

REACTIVE = {'BRI': None, 'ISO_z2': 2, 'ISO_z4': 4, 'ISO_z8': 8, 'BRI+BASAL': None}
POOLS = ('bridging', 'reconstructed_row', 'basal_L1', 'L1_subbridging', 'subsurface_L2_4', 'bulk')


def energy_sets():
    sets = {}
    with open(HERE / 'specification' / 'energy_sets.csv') as f:
        for r in csv.DictReader(f):
            sets.setdefault(r['scenario'], {})[(r['layer'], r['site'])] = float(r['relative_E_eV'])
    return sets


def energy_maps():
    s = energy_sets()
    pab = {(int(k), site): e for (k, site), e in s['PABISIAK_2007_GGA'].items()}
    ham = {(int(k), site): e for (k, site), e in s['HAMEEUW_2006_LDA'].items()}
    li = s['LI_2015_SX']
    surf, sub = li[('SURFACE', 'GENERIC')], li[('SUBSURFACE', 'GENERIC')]
    basal = {(1, 'IPL'): surf + MATSUNAGA}
    return {
        'PAB': pab,
        'HAM': ham,
        'LI_SBR1': {(1, 'BRI'): surf, (1, 'SBR'): sub, **basal},
        'LI_L2': {(1, 'BRI'): surf, **{(2, x): sub for x in ('BRI', 'IPL', 'SBR')}, **basal},
    }


def zha2017(r_nm):
    """ZHA2017 bulk vacancy pair energy, eV (r converted to angstrom)."""
    A, B, C = ZHA
    r = 10.0 * r_nm
    return A / r - B / r ** 2 - C / r ** 6


def _two_state(C, E, tag, region, shell, electron=False):
    d = dict(C=C, E=[0.0, E], tag=tag, region=region, shell=shell)
    d.update(v=[0, 0], e=[0, 1]) if electron else d.update(v=[0, 1])
    return d


def recon_cell(E, dG, eT):
    """States of one (1x2) cell with its four layer-1 Ti.
    Returns the state table and the O-state of each state (0: 00, 1: one
    vacant, 2: both vacant, 3: reconstructed row)."""
    o_states = [(0.0, 0, 0, 1, 4), (E, 1, 0, 2, 4), (2 * E, 2, 0, 1, 4), (E + dG, 1, 2, 1, 2)]
    Es, vs, es, gs, os_ = [], [], [], [], []
    for o, (Eo, v, e_int, g, n_free) in enumerate(o_states):
        for n in range(n_free + 1):
            Es.append(Eo + n * eT); vs.append(v); es.append(e_int + n)
            gs.append(g * math.comb(n_free, n)); os_.append(o)
    return dict(E=Es, v=vs, e=es, g=gs), os_


class Layout:
    def __init__(self, **kw):
        self.__dict__.update(kw)


def particle(energy_map, d_nm, *, cutoff, dG, f110, eps_r, closure='GLOBAL', s0=None, n_bulk=N_BULK):
    """One sphere. closure='LOCAL' (each trilayer and bulk shell neutral) is
    the eps -> 0 limit; the solver starts from it, and tests use it."""
    from .aggregates import BOX, ln_q
    if closure not in ('GLOBAL', 'LOCAL'):
        raise ValueError(closure)
    s0 = (S0_PENALTY if closure == 'GLOBAL' else 0.0) if s0 is None else s0
    R = d_nm / 2
    kT = pt.KB * T_EQ
    lnq = ln_q(T_EQ, cutoff)
    box_sites = BOX[0] * BOX[1] * BOX[2] * 4
    doms, radii = [], {}

    def add(dom, depth):
        doms.append(dom)
        radii[dom['shell']] = R - depth
        return len(doms) - 1

    emap = energy_maps()[energy_map]
    o_sites, ti_layers = pt.layer_sites(d_nm, K_EXPL)
    o_sites = [(k, site, C * f110, z) for k, site, C, z in o_sites]
    ti_layers = [C * f110 for C in ti_layers]
    idx = {}
    cells = ti_layers[0] / 4               # the Ti plane sets the cell count
    for k, site, C, z in o_sites:
        E = emap.get((k, site), 0.0)
        if (k, site) == (1, 'BRI'):
            states, cell_o = recon_cell(E, dG, s0)
            idx['cell'] = add(dict(C=cells, tag=(1, 'BRI'), region=1, shell=(1, 0), shell2=(1, 1), **states), z)
            # The Ti plane lies deeper, so its area is 3e-4 smaller: the bridging
            # O left over stay plain bridging sites.
            idx['bri_rest'] = add(_two_state(C - 2 * cells, E, (1, 'BRI'), 1, (1, 0)), z)
            continue
        idx[(k, site)] = add(_two_state(C, E, (k, site), k, (k, PLANES[site])), z)
    for k, C in enumerate(ti_layers[1:], 2):
        add(_two_state(C, 0.0, 'Ti', k, (k, 1), electron=True), (k - 1) * pt.D110 + pt.H_BRI)
    top = K_EXPL * pt.D110
    edges = top + np.concatenate([[0.0], np.geomspace(Z0_BULK, R - top, n_bulk)])
    edges[-1] = R
    o_rem = pt.O_TOTAL - sum(c for _, _, c, _ in o_sites)
    ti_rem = pt.TI_TOTAL - sum(ti_layers)
    whole = pt.shell_fraction(R, top, R)
    n = np.arange(len(lnq))
    for j in range(n_bulk):
        lo, hi = edges[j], edges[j + 1]
        f = pt.shell_fraction(R, lo, hi) / whole
        zc, reg, sh = 0.5 * (lo + hi), K_EXPL + 1 + j, ('bulk', j)
        add(dict(C=o_rem * f / box_sites, E=(-kT * lnq).tolist(), v=n.tolist(), tag='bulk',
                 region=reg, shell=sh, box=box_sites), zc)
        add(_two_state(ti_rem * f, 0.0, 'Ti', reg, sh, electron=True), zc)
    es = dict(radii=radii, eps=eps_r, mass_g=pt.particle_mass(d_nm)) if closure == 'GLOBAL' else None
    model = Model(doms, electrostatics=es)
    return model, Layout(energy_map=energy_map, d_nm=d_nm, cutoff=cutoff, dG=dG, f110=f110, eps_r=eps_r,
                         closure=closure, c_bri=pt.bridging_capacity(d_nm) * f110, cell_o=np.array(cell_o),
                         idx=idx, edges=edges)


def surface(sol, lay):
    """(theta on unreconstructed bridging sites, their capacity, reconstructed fraction)."""
    m, o = sol.model, lay.cell_o
    i, r = lay.idx['cell'], lay.idx['bri_rest']
    x = np.bincount(o, sol.x[i][:len(o)], 4)
    vac = x[1] + 2 * x[2] + float(sol.x[r][1])
    total = 2 * m.C[i] + m.C[r]
    sites = total - 2 * x[3]
    f_rec = 2 * x[3] / total
    return (vac / sites if sites > 0 else 0.0), lay.c_bri * (1 - f_rec), f_rec


def reactive_sites(sol, lay):
    """Reactive sites (umol/g) for each definition, plus theta and f_rec."""
    th, c, f_rec = surface(sol, lay)
    basal = float(sol.x[lay.idx[(1, 'IPL')]] @ sol.model.v[lay.idx[(1, 'IPL')]])
    out = {}
    for name, z in REACTIVE.items():
        n = c * th * (1 - th) ** z if z else c * th
        out[name] = n + basal if name == 'BRI+BASAL' else n
    return out, th, f_rec


def populations(sol, lay):
    """Vacancies (umol/g) in each pool."""
    m, x, o = sol.model, sol.x, lay.cell_o
    i = lay.idx['cell']
    per = x[i][:len(o)] * m.v[i][:len(o)]
    out = dict.fromkeys(POOLS, 0.0)
    out['bridging'] = float(per[o != 3].sum())
    out['reconstructed_row'] = float(per[o == 3].sum())
    for k, tag in enumerate(m.tags):
        if k == i or tag == 'Ti':
            continue
        n = float(x[k] @ m.v[k])
        if tag == (1, 'BRI'):
            out['bridging'] += n
        elif tag == (1, 'IPL'):
            out['basal_L1'] += n
        elif tag == (1, 'SBR'):
            out['L1_subbridging'] += n
        elif tag == 'bulk':
            out['bulk'] += n
        else:
            out['subsurface_L2_4'] += n
    return out
