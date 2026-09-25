"""Particle models: one state table per (energy map, closure, diameter, options).

Discrete maps (PAB, HAM, LI_SBR1, LI_L2)
    Four explicit (110) trilayers. Each holds BRI (1), IPL (2) and SBR (1) O
    sites and 2 Ti per 1x1 cell, at their crystallographic depths. Sites a map
    does not name have energy 0 (bulk). Below them, radial bulk shells take the
    remaining O and Ti.
Continuum map (LI_CONT, manuscript Note 2b)
    Spherical shells from the surface to the centre, G(z) = -A exp(-z/xi).

Closures
    NEUTRAL  O sites only, electrons implicit.
    LOCAL    each trilayer (discrete) or shell (continuum) is a neutral region
             with its own Ti3+ pool (ideal Ti configurational entropy).
    GLOBAL   the same regions, plus electrostatics. Every atomic plane is a
             charged shell (vacancies at O planes, Ti3+ at Ti planes). Ti in
             the surface trilayer carries the RET2018 S0 penalty, 0.2 eV above
             the subsurface optimum. All other Ti are equal (declared).

Bulk aggregates of any size (aggregates=cutoff_nm)
    Bulk O below the four explicit trilayers is grouped into periodic boxes of
    aggregates.BOX sites. A box with N vacancies has free energy -kT ln Q(N)
    from transition-matrix Monte Carlo (pairwise-additive ZHA2017 energy
    within the cutoff; two Ti3+ per vacancy on its own Ti neighbours).
Reconstruction
    ('fixed', f): a fraction f of the surface is Ti2O3-(1x2). It holds 0.5
        vacancy equivalents per 1x1 cell (one per Ti2O3 row unit, i.e. per
        1x2 cell) with its own Ti3+, and no reactive bridging sites.
    ('state', dG): each 1x2 surface cell is either two bridging sites
        (4 occupancy states) or reconstructed, with one vacancy equivalent,
        internal Ti3+ and energy eps_BRI + dG. The reconstructed fraction is
        then an equilibrium output. Discrete maps only.

Bulk vacancy pairs (pairs=True)
    Bulk O below the four explicit trilayers (1.30 nm) can hold vacancy pairs.
    The pair energy is the ZHA2017 bulk potential E(r) = A/r - B/r^2 - C/r^6,
    taken at every rutile O-O separation up to its 1 nm range. Each distance
    enters with Mayer weight (n_k/2)(exp(-E_k/kT) - 1), which is exact to
    second order in vacancy fugacity. The pair states are therefore built for
    one temperature.
"""
import csv
import math
from pathlib import Path

import numpy as np

from . import particle as pt
from .engine import Model

HERE = Path(__file__).resolve().parent.parent
K_EXPL = 4
N_BULK, Z0_BULK = 240, 0.02
N_CONT, Z0_CONT = 1500, 1e-4
A_LI = 1.31
S0_PENALTY = 0.2
ZHA = (0.753, 11.83, -383.45)
PAIR_RMAX = 1.0
PAIR_DEPTH = K_EXPL * pt.D110
PLANES = {'BRI': 0, 'IPL': 1, 'SBR': 2}
RECON_ROW_VACANCIES = 0.5          # per 1x1 cell: one O short per Ti2O3 row unit (1x2 cell)


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
        'PAB': pab,
        'HAM': ham,
        'LI_SBR1': {(1, 'BRI'): surf, (1, 'SBR'): sub},
        'LI_L2': {(1, 'BRI'): surf, **{(2, x): sub for x in ('BRI', 'IPL', 'SBR')}},
    }


def zha2017(r_nm):
    """ZHA2017 bulk vacancy pair energy, eV (r converted to angstrom)."""
    A, B, C = ZHA
    r = 10.0 * r_nm
    return A / r - B / r ** 2 - C / r ** 6


def pair_states(eps_site, T):
    """Pair-cell states per O site: none, or a pair at each O-O distance."""
    kT = pt.KB * T
    E, g = [0.0], [1.0]
    for r, n in pt.oxygen_shells(PAIR_RMAX):
        Ek = zha2017(r)
        if Ek >= 0:
            raise ValueError('ZHA2017 is attractive at every rutile O-O distance within 1 nm')
        E.append(2 * eps_site + Ek)
        g.append(0.5 * n * -math.expm1(Ek / kT))
    return dict(E=E, v=[0] + [2] * (len(E) - 1), g=g)


def _recon_cell(C_bri, E, dG, ti_cap, eT):
    """States of one 1x2 surface cell: two bridging sites (00, one, both
    vacant) or the reconstructed row (one vacancy equivalent, energy E + dG).

    LOCAL (ti_cap given): the cell also owns its layer-1 Ti, four per cell,
    each empty or holding one Ti3+ electron (energy eT). The reconstructed row
    keeps two of them for its own two electrons, so only two stay free. The
    number of cells is set by the Ti plane (ti_cap / 4); the bridging O left
    over (the Ti plane lies slightly deeper, so its area is 3e-4 smaller) stay
    plain bridging sites. Returns (cells, state dict, O-state of each state).
    """
    o_states = [(0.0, 0, 0, 1, 4), (E, 1, 0, 2, 4), (2 * E, 2, 0, 1, 4), (E + dG, 1, 2, 1, 2)]
    if ti_cap is None:
        return C_bri / 2, dict(E=[s[0] for s in o_states], v=[s[1] for s in o_states],
                               e=[s[2] for s in o_states], g=[s[3] for s in o_states]), [0, 1, 2, 3]
    cells = min(C_bri / 2, ti_cap / 4)
    Es, vs, es, gs, os_ = [], [], [], [], []
    for o, (Eo, v, e_int, g, n_free) in enumerate(o_states):
        for n in range(n_free + 1):
            Es.append(Eo + n * eT); vs.append(v); es.append(e_int + n)
            gs.append(g * math.comb(n_free, n)); os_.append(o)
    return cells, dict(E=Es, v=vs, e=es, g=gs), os_


def _two_state(C, E, tag, region, shell, electron=False):
    d = dict(C=C, E=[0.0, E], tag=tag, region=region, shell=shell)
    d.update(v=[0, 0], e=[0, 1]) if electron else d.update(v=[0, 1])
    return d


class Layout:
    """What the scenario and kinetics code need to read back from a model."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def build(name, closure, d_nm, *, xi=0.5, eps_r=None, pairs=False, T=873.15,
          n_cont=N_CONT, n_bulk=N_BULK, z0_cont=Z0_CONT, s0=None, amp=A_LI,
          aggregates=None, recon=None, f110=1.0, basal_shift=None):
    """aggregates: MC cutoff in nm (bulk boxes replace bulk O sites), or None.
    recon: None, ('fixed', f) or ('state', dG_eV); see the module docstring.
    f110: (110) share of the surface. The rest carries no explicit surface
        sites or surface energy (other facets have no sourced energies).
    basal_shift: layer-1 in-plane (basal) O energy = bridging energy + shift
        (Matsunaga 2014: +0.11 eV); used for the Li maps, which give no basal value."""
    if closure not in ('NEUTRAL', 'LOCAL', 'GLOBAL'):
        raise ValueError(closure)
    if (closure == 'GLOBAL') != (eps_r is not None):
        raise ValueError('GLOBAL needs eps_r, and only GLOBAL takes it')
    charged = closure != 'NEUTRAL'
    s0 = (S0_PENALTY if closure == 'GLOBAL' else 0.0) if s0 is None else s0
    R = d_nm / 2
    doms, radii, o_list, ti_cap, pair_idx = [], {}, [], {}, []

    def add(dom, depth=None, eps=None):
        doms.append(dom)
        if not charged:
            dom['region'] = -1
            dom.pop('shell')
        elif dom.get('shell') is not None:
            radii[dom['shell']] = R - depth
        idx = len(doms) - 1
        if eps is not None:
            o_list.append(dict(idx=idx, z=depth, eps=eps, C=dom['C'] * dom.get('box', 1),
                               region=dom['region'], box=dom.get('box', 1)))
        return idx

    lnq = None
    if aggregates is not None:
        from .aggregates import BOX, ln_q
        lnq = ln_q(T, aggregates)
        box_sites = BOX[0] * BOX[1] * BOX[2] * 4
        kT = pt.KB * T

    def o_domain(C, E, tag, region, shell, lo):
        """Two-state O sites, or aggregate boxes below the explicit trilayers."""
        if lnq is None or lo < PAIR_DEPTH - 1e-12:
            return _two_state(C, E, tag, region, shell)
        n = np.arange(len(lnq))
        return dict(C=C / box_sites, E=(-kT * lnq + n * E).tolist(), v=n.tolist(),
                    tag=tag, region=region, shell=shell, box=box_sites)

    def add_pairs(C, eps_site, region, shell, zc, lo):
        # Pairs only in shells that lie wholly below the explicit trilayers.
        if lo >= PAIR_DEPTH - 1e-12:
            o_list[-1]['eligible'] = True
            if pairs:
                pair_idx.append(add(dict(C=C, tag='pair', region=region, shell=shell,
                                         **pair_states(eps_site, T)), zc))

    if recon and recon[0] == 'state' and name == 'LI_CONT':
        raise ValueError('reconstruction states need explicit bridging sites (discrete maps)')
    if name == 'LI_CONT':
        edges = np.concatenate([[0.0], np.geomspace(z0_cont, R, n_cont)])
        for j in range(n_cont):
            lo, hi = edges[j], edges[j + 1]
            zc = 0.5 * (lo + hi)
            C = pt.O_TOTAL * pt.shell_fraction(R, lo, hi)
            E = -amp * math.exp(-zc / xi)
            add(o_domain(C * f110, E, 'shell', j, j, lo), zc, eps=E)
            if f110 < 1.0:
                add(o_domain(C * (1 - f110), 0.0, 'other_facet', j, j, lo), zc, eps=0.0)
            if charged:
                eT = s0 if zc < pt.D110 else 0.0
                i = add(_two_state(C / 2, eT, 'Ti', j, j, electron=True), zc)
                ti_cap[j] = (i, C / 2)
            add_pairs(C, E, j, j, zc, lo)
        lay = dict(kind='continuum', xi=xi, edges=edges)
    else:
        emap = dict(discrete_maps()[name])
        if basal_shift is not None:
            emap[(1, 'IPL')] = emap[(1, 'BRI')] + basal_shift
        o_sites, ti_layers = pt.layer_sites(d_nm, K_EXPL)
        o_sites = [(k, site, C * f110, z) for k, site, C, z in o_sites]
        ti_layers = [C * f110 for C in ti_layers]
        ti_rec, ti_explicit = ti_layers[0], sum(ti_layers)
        bri = basal = None
        cell_o, bri_rest = None, None
        if recon and recon[0] == 'state' and closure == 'GLOBAL':
            raise ValueError('reconstruction states are built for NEUTRAL and LOCAL only')
        if recon and recon[0] == 'fixed' and charged:
            # The Ti2O3 rows hold their own Ti3+ (two per row vacancy): those Ti
            # leave the free layer-1 Ti pool.
            ti_layers[0] -= 2 * RECON_ROW_VACANCIES * recon[1] * pt.bridging_capacity(d_nm) * f110
        for k, site, C, z in o_sites:
            E = emap.get((k, site), 0.0)
            if (k, site) == (1, 'BRI') and recon and recon[0] == 'state':
                cells, states, cell_o = _recon_cell(C, E, recon[1], ti_layers[0] if charged else None,
                                                    s0 if charged else 0.0)
                bri = add(dict(C=cells, tag=(k, site), region=k, shell=(k, PLANES[site]), **states),
                          z, eps=E)
                if C - 2 * cells > 0:
                    bri_rest = add(_two_state(C - 2 * cells, E, (k, site), k, (k, PLANES[site])), z, eps=E)
                if charged:
                    ti_layers[0] = 0.0          # every layer-1 Ti now sits in a cell
                continue
            if (k, site) == (1, 'BRI') and recon and recon[0] == 'fixed':
                C = C * (1 - recon[1])
            i = add(_two_state(C, E, (k, site), k, (k, PLANES[site])), z, eps=E)
            if (k, site) == (1, 'BRI'):
                bri = i
            if (k, site) == (1, 'IPL'):
                basal = i
        if charged:
            for k, C in enumerate(ti_layers, 1):
                if C <= 0:
                    ti_cap[k] = (bri, ti_rec)
                    continue
                eT = s0 if k == 1 else 0.0
                i = add(_two_state(C, eT, 'Ti', k, (k, 1), electron=True),
                        (k - 1) * pt.D110 + pt.H_BRI)
                ti_cap[k] = (i, C)
        top = K_EXPL * pt.D110
        edges = top + np.concatenate([[0.0], np.geomspace(Z0_BULK, R - top, n_bulk)])
        edges[-1] = R
        o_rem = pt.O_TOTAL - sum(c for _, _, c, _ in o_sites)
        ti_rem = pt.TI_TOTAL - ti_explicit     # layer Ti moved into cells or reserved stay in layer 1
        whole = pt.shell_fraction(R, top, R)
        for j in range(n_bulk):
            lo, hi = edges[j], edges[j + 1]
            f = pt.shell_fraction(R, lo, hi) / whole
            zc, reg, sh = 0.5 * (lo + hi), K_EXPL + 1 + j, ('bulk', j)
            add(o_domain(o_rem * f, 0.0, 'bulk', reg, sh, lo), zc, eps=0.0)
            if charged:
                i = add(_two_state(ti_rem * f, 0.0, 'Ti', reg, sh, electron=True), zc)
                ti_cap[reg] = (i, ti_rem * f)
            add_pairs(o_rem * f, 0.0, reg, sh, zc, lo)
        lay = dict(kind='discrete', bri=bri, basal=basal, edges=edges, cell_o=cell_o, bri_rest=bri_rest)
    es = None
    if closure == 'GLOBAL':
        es = dict(radii=radii, eps=eps_r, mass_g=pt.particle_mass(d_nm))
    model = Model(doms, electrostatics=es)
    c_bri = pt.bridging_capacity(d_nm) * f110
    f_fix = recon[1] if recon and recon[0] == 'fixed' else 0.0
    return model, Layout(name=name, closure=closure, d_nm=d_nm, R=R, T=T, pairs=pairs,
                         c_bri=c_bri, o=o_list, ti=ti_cap, pair_idx=pair_idx, recon=recon,
                         aggregates=aggregates, f_fixed=f_fix,
                         fixed=dict(v=RECON_ROW_VACANCIES * f_fix * c_bri) if f_fix else {}, **lay)


def surface_coverage(sol, lay):
    """Bridging-vacancy coverage of the outer surface.

    Discrete maps: the layer-1 BRI site. Continuum: the site fraction at z = 0
    (manuscript Note 2b's theta = x(R)), from the stationarity condition there.
    """
    if lay.kind == 'discrete':
        return sol.occupancy(lay.bri)
    kT = pt.KB * sol.T
    G0 = -A_LI
    if lay.closure == 'NEUTRAL':
        return 1.0 / (1.0 + math.exp((G0 - sol.mu_eV) / kT))
    if lay.closure == 'GLOBAL':
        # The potential is continuous at the surface; shell 0 is the outermost.
        return 1.0 / (1.0 + math.exp((G0 + 2 * sol.phi[0] - sol.mu_eV) / kT))
    lo, hi = 0.0, 0.25
    for _ in range(200):
        x = 0.5 * (lo + hi)
        f = G0 + kT * (math.log(x / (1 - x)) + 2 * math.log(4 * x / (1 - 4 * x))) - sol.mu_eV
        lo, hi = (lo, x) if f > 0 else (x, hi)
    return 0.5 * (lo + hi)


def surface_state(sol, lay):
    """(theta on unreconstructed bridging sites, their capacity, reconstructed fraction)."""
    if lay.recon and lay.recon[0] == 'state':
        o = np.asarray(lay.cell_o)
        x = np.bincount(o, sol.x[lay.bri][:len(o)], 4)
        cells = sol.model.C[lay.bri]
        vac, sites = x[1] + 2 * x[2], 2 * (cells - x[3])
        if lay.bri_rest is not None:
            vac += float(sol.x[lay.bri_rest][1]); sites += sol.model.C[lay.bri_rest]
        f_rec = 2 * x[3] / (2 * cells + (sol.model.C[lay.bri_rest] if lay.bri_rest is not None else 0.0))
        theta = vac / sites if sites > 0 else 0.0
        return theta, lay.c_bri * (1 - f_rec), f_rec
    return surface_coverage(sol, lay), lay.c_bri * (1 - lay.f_fixed), lay.f_fixed


def pair_fraction(sol, lay):
    """Share of the bulk-pair-eligible vacancies that sit in pairs."""
    if not lay.pair_idx:
        return 0.0
    paired = sum(float(sol.x[i] @ sol.model.v[i]) for i in lay.pair_idx)
    mono = sum(float(sol.x[o['idx']] @ sol.model.v[o['idx']]) for o in lay.o
               if o.get('eligible'))
    return paired / (paired + mono) if paired + mono > 0 else 0.0
